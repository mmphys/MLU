/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: MultiFit.cpp
 
 Copyright (C) 2019
 
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 
 See the full license in the file "LICENSE" in the top level distribution directory
 *************************************************************************************/
/*  END LEGAL */

//#include "CommonLatAn.hpp"
#include "Common.hpp"
#include "MultiFit.hpp"

// Uncomment the next line if your cmath doesn't define M_PI etc by default
//#define _USE_MATH_DEFINES
#include <cmath>
//#include <iomanip>
#include <ios>
#include <iostream>
//#include <sys/stat.h>
//#include <LatAnalyze/Core/OptParser.hpp>
//#include <LatAnalyze/Eigen/Dense>
//#include <LatAnalyze/Functional/CompiledModel.hpp>
//#include <LatAnalyze/Io/Io.hpp>
//#include <LatAnalyze/Statistics/MatSample.hpp>
//#include <LatAnalyze/Core/Math.hpp>
//#include <LatAnalyze/Numerical/MinuitMinimizer.hpp>
//#include <LatAnalyze/Numerical/NloptMinimizer.hpp>
//#include <LatAnalyze/Core/Plot.hpp>
//#include <LatAnalyze/Statistics/XYSampleData.hpp>

#include <Minuit2/FCNBase.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/VariableMetricMinimizer.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MnPrint.h>

// Indices for operators in correlator names
static constexpr int idxSrc{ 0 };
static constexpr int idxSnk{ 1 };
static const char * pSrcSnk[] = { "src", "snk" };

using scalar = double;
using Fold = Common::Fold<scalar>;
using Model = Common::Model<scalar>;
using Matrix = Eigen::MatrixXd; // dynamic sized matrix of complex double
//using SampleD = Common::SampleD; // Bootstrap sample double
//using SampleC = Common::SampleC; // Bootstrap sample std::complex<double>

// Minimisation for multi-exponential fit

class MultiExpModel : public ROOT::Minuit2::FCNBase {
public:
  static constexpr double pi{ M_PI };
  explicit MultiExpModel( std::vector<Fold> &&Corr, const std::vector<std::string> &OpNames,
                         bool bFactor, const std::string &sOpNameConcat,
                         int NumExponents, int Verbosity, const std::string &OutputBaseName,
                         Common::SeedType Seed, int NSamples );
  virtual ~MultiExpModel() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const std::vector<double> & ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  // These are definitely not part of the FCNBase interface
  std::vector<Common::ValWithEr> PerformFit( bool bCorrelated, int tMin, int tMax, int Skip,
        bool bSaveCorr, int MaxIt, double Tolerance, double RelEnergySep, double &ChiSq, int &dof );
  const int NumOps;
  const std::vector<std::string> &OpNames;
  bool bFactor;
  const std::string &sOpNameConcat;
  const int Verbosity;
  const std::string OutputBaseName;
  const Common::SeedType Seed;
  enum ModelType { sinh = -1, exp, cosh };
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int NumExponents;
  const int NumParams;
  const std::vector<Fold> Corr;
  const std::vector<std::string> ParamNames;
  private:
  // These variables only used during a fit
  int tMin = -1;
  int tMax = -1;
  int NtCorr;
  int Extent;
  Matrix Covar;
  Matrix CovarInv;
  std::vector<double> VarianceInv;
  bool bCorrelated;
  std::string OutputRunBase;
  int idx;
  // Helper functions
  inline int AlphaIndex( int snk, int src) const { return snk * NumOps + src; };
  inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; };
  void MakeCovar( void );
  void DumpParameters( const std::string &msg, const ROOT::Minuit2::MnUserParameters &par );
  std::vector<std::string> MakeParamNames();
};

MultiExpModel::MultiExpModel( std::vector<Fold> &&Corr_, const std::vector<std::string> &opNames_,
                             bool bfactor_, const std::string &sopNameConcat_,
                             int numExponents_, int verbosity_, const std::string &outputBaseName_, Common::SeedType seed_, int nSamples_ )
  : NumOps{ static_cast<int>( opNames_.size() ) },
    OpNames{ opNames_ },
    bFactor{ bfactor_ },
    sOpNameConcat{ sopNameConcat_ },
    Verbosity{ verbosity_ },
    OutputBaseName{ outputBaseName_ },
    Seed{ seed_ },
    NumFiles{ static_cast<int>( Corr_.size() ) },
    NSamples{ ( nSamples_ <= 0 || nSamples_ > Corr_[0].NumSamples() ) ? Corr_[0].NumSamples() : nSamples_ },
    Nt      { Corr_[0].Nt() },
    NumExponents{ numExponents_ },
    NumParams{ NumExponents * ( 1 + NumOps ) },
    Corr{ std::move( Corr_ ) },
    ParamNames( MakeParamNames() )
{
}

std::vector<std::string> MultiExpModel::MakeParamNames()
{
  assert( NumExponents > 0 && "Number of exponents in the fit should be positive" );
  std::vector<std::string> Names( NumParams );
  for( int e = 0; e < NumExponents; e++ ) {
    Names[e] = "E" + std::to_string( e );
    for( int op = 0; op < NumOps; op++ )
      Names[NumExponents + MELIndex( op, e )] = "A" + OpNames[op] + std::to_string( e );
  }
  return Names;
}

// Make the covariance matrix, for only the timeslices we are interested in
void MultiExpModel::MakeCovar( void )
{
  assert( tMin >= 0 && tMin < Nt && "Error: tMin invalid" );
  assert( tMax >= 0 && tMax < Nt && "Error: tMax invalid" );
  assert( NtCorr > 1 && "Error: tMin >= tMax" );
  assert( Extent == NumFiles * NtCorr && "Error: Extent != NumFiles * NtCorr" );
  if( bCorrelated )
  {
    std::cout << "Creating covariance matrix ...\n";
    Covar.resize( Extent, Extent );
    CovarInv.resize( Extent, Extent );
    if( Verbosity )
      std::cout << "Covariance matrix is " << Covar.rows() << " x " << Covar.cols() << "\n";
  }
  // Make covariance
  VarianceInv.resize( Extent );
  {
    const scalar * CentralX = nullptr;
    const scalar * CentralY = nullptr;
    for( int x = 0; x < Extent; x++ )
    {
      const int t1{ x % NtCorr };
      if( !t1 )
        CentralX = Corr[x / NtCorr][Fold::idxCentral] + tMin;
      for( int y = bCorrelated ? 0 : x; y <= x; y++ )
      {
        const int t2{ y % NtCorr };
        if( !t2 )
          CentralY = Corr[y / NtCorr][Fold::idxCentral] + tMin;
        double z = 0;
        const scalar * DataX{ CentralX };
        const scalar * DataY{ CentralY };
        for( int i = 0; i < NSamples; i++ )
        {
          DataX += Nt;
          DataY += Nt;
          z += ( DataX[t1] - CentralX[t1] ) * ( DataY[t2] - CentralY[t2] );
        }
        z /= NSamples;
        if( bCorrelated )
        {
          Covar( x, y ) = z;
          if( x != y )
            Covar( y, x ) = z;
        }
        if( y == x )
          VarianceInv[x] = 1. / z;
      }
    }
  }
  if( bCorrelated )
  {
    if( !Common::IsFinite( Covar ) )
      throw std::runtime_error( "Covariance matrix isn't finite" );
    CovarInv = Covar.inverse();
    if( !Common::IsFinite( CovarInv ) )
      throw std::runtime_error( "Inverse covariance matrix isn't finite" );
    const double CondNumber{ Covar.norm() * CovarInv.norm() };
    const int CondDigits{ static_cast<int>( std::log10( CondNumber ) + 0.5 ) };
    if( CondDigits >= 12 )
      std::cout << "WARNING see https://en.wikipedia.org/wiki/Condition_number\n";
    std::cout << "Covariance matrix condition number=" << CondNumber
              << ", i.e. potential loss of " << CondDigits << " digits.\n";
  }
  if( !Common::IsFinite( VarianceInv ) )
    throw std::runtime_error( "Variance vector isn't finite" );
}

// Compute chi-squared given parameters for multi-exponential fit
  // Energy Levels        Parameters 0 ... NumExponents - 1
  //  NB: second and subsequent energy levels are deltas
  // Overlap Coefficients Parameters NumExponents ... NumExponents ( NumOps + 1 ) - 1

double MultiExpModel::operator()( const std::vector<double> & par ) const
{
  const double * const Energy{ par.data() };
  const double * const Coeff{ Energy + NumExponents };
  // Calculate the theory errors for these model parameters
  const int HalfNt{ Nt - 1 }; // NB: This is half of the ORIGINAL correlator time dimension
  double SinhCoshAdjust[NumExponents];
  for( int e = 0; e < NumExponents; ++e )
    SinhCoshAdjust[e] = std::exp( - Energy[e] * HalfNt );
  double ModelError[Extent]; // Should happily fit on stack
  for( int f = 0; f < NumFiles; f++ )
  {
    const int snk{ Corr[f].Name_.op[idxSnk] };
    const int src{ Corr[f].Name_.op[idxSrc] };
    const scalar * CorrData = Corr[f][idx] + tMin;
    for( int t = 0; t < NtCorr; t++ )
    {
      double z = 0;
      for( int e = 0; e < NumExponents; ++e )
      {
        double d;
        switch( Corr[f].parity )
        {
          case Common::Parity::Even:
            d = SinhCoshAdjust[e] * std::cosh( - Energy[e] * ( tMin + t - HalfNt ) );
            break;
          case Common::Parity::Odd:
            d = SinhCoshAdjust[e] * std::sinh( - Energy[e] * ( tMin + t - HalfNt ) );
            break;
          default:
            d = std::exp( - Energy[e] * ( tMin + t ) );
            break;
        }
        z += d * Coeff[MELIndex( src, e )] * Coeff[MELIndex( snk, e )];
      }
      z -= *CorrData++;
      if( !std::isfinite( z ) )
      {
        if( Verbosity )
        {
          static int iOOB{ 0 };
          std::cout << "index " << idx << ", " << ++iOOB << "th overflow\n";
        }
        return std::numeric_limits<double>::max();
      }
      const int iWrite{ f * NtCorr + t };
      ModelError[iWrite] = z;
    }
  }
  // The inverse of a symmetric matrix is also symmetric, so only calculate half the matrix
  double chi2 = 0;
  for( int i = 0; i < Extent; ++i )
  {
    if( !bCorrelated )
      chi2 += ModelError[i] * VarianceInv[i] * ModelError[i];
    else
    {
      for( int j = 0; j <= i; ++j )
      {
        double z = ModelError[i] * CovarInv(i, j) * ModelError[j];
        chi2 += z;
        if( i != j )
          chi2 += z;
      }
    }
  }
  if( Verbosity > 2 )
    std::cout << "index " << idx << ", chi^2=" << chi2 << ", E0=" << Energy[0] << ", E1=" << Energy[1] << "\n";
  if( chi2 < 0 )
  {
    //static int count{ 0 };
    //std::cout << "Chi^2 < 0" << ++count << "\n";
    throw std::runtime_error( "Chi^2 < 0" );
  }
  return chi2;
}

void MultiExpModel::DumpParameters( const std::string &msg, const ROOT::Minuit2::MnUserParameters &par )
{
  static const char NewLine[] = "\n";
  static const char Tab[] = "\t";
  if( !msg.empty() )
    std::cout << msg << NewLine;
  for( int p = 0; p < NumParams; p++ )
    std::cout << Tab << ParamNames[p] << Tab << par.Value( p ) << NewLine;
}

// Perform a fit
std::vector<Common::ValWithEr> MultiExpModel::PerformFit( bool Bcorrelated_, int TMin_, int TMax_,
  int Skip, bool bSaveCorr, int MaxIt, double Tolerance, double RelEnergySep, double &ChiSq, int &dof )
{
  std::cout << "=========================================\nPerforming "
            << ( Bcorrelated_ ? "correlated" : "uncorrelated" )
            << " fit on timeslices " << TMin_ << " to " << TMax_ << "\n";

  // Drop the covariance & variance if the fit timeslices change
  const bool bCorrelatedChanged{ bCorrelated != Bcorrelated_ };
  const bool bTimesChanged{ TMin_ != tMin || TMax_ != tMax };
  bCorrelated = Bcorrelated_;
  tMin = TMin_;
  tMax = TMax_;
  if( bTimesChanged || bCorrelatedChanged )
  {
    OutputRunBase = OutputBaseName;
    OutputRunBase.append( 1, '.' );
    OutputRunBase.append( bCorrelated ? "corr" : "uncorr" );
    OutputRunBase.append( 1, '_' );
    OutputRunBase.append( std::to_string( tMin ) );
    OutputRunBase.append( 1, '_' );
    OutputRunBase.append( std::to_string( tMax ) );
  }
  if( bTimesChanged )
  {
    NtCorr = TMax_ - TMin_ + 1;
    Extent = NumFiles * NtCorr;
    VarianceInv.resize( 0 );
    CovarInv.resize( 0, 0 );
    Covar.resize( 0, 0 );
  }
  dof = Extent;
  dof -= NumParams;
  // Make the covariance matrix / vector if need be
  if( ( bCorrelated && Covar.rows() != Extent ) || ( !bCorrelated && VarianceInv.size() != Extent ) )
    MakeCovar();

  // Make initial guesses for the parameters
  // For each Exponent, I need the delta_E + a constant for each operator
  std::vector<double> par( NumParams );
  {
    // Take a starting guess for the parameters - same as LatAnalyze
    double * const Energy{ par.data() };
    double * const Coeff{ Energy + NumExponents };
    const int tGuess{ Nt / 4 };
    Energy[0] = 0;
    for( int f = 0; f < NumFiles; f++ )
    {
      //const int i{ AlphaIndex( snk, src ) * Nt + tGuess };
      const double * const c{ Corr[f][Fold::idxCentral] };
      double E0 = std::log( c[tGuess] / c[tGuess + 1] );
      Energy[0] += E0;
      double MELGuess = std::sqrt( std::abs( c[tGuess] / ( std::exp( - E0 * tGuess ) ) ) );
      Coeff[MELIndex( Corr[f].Name_.op[idxSnk], 0 )] = MELGuess;
      Coeff[MELIndex( Corr[f].Name_.op[idxSrc], 0 )] = MELGuess;
    }
    Energy[0] /= NumFiles;
    // Now guess Higher exponents - same as LatAnalyze
    //static const double MELFactor{ std::sqrt( 0.5 ) };
    static const double MELFactor{ std::sqrt( 2 ) };
    for( int e = 1; e < NumExponents; ++e )
    {
      Energy[e] = Energy[e - 1] + ( Energy[e - 1] - ( e > 1 ? Energy[e - 2] : 0 ) ) * 0.5;
      for( int o = 0; o < NumOps; ++o )
        Coeff[ MELIndex( o, e ) ] = Coeff[ MELIndex( o, e - 1 ) ] * MELFactor;
    }
  }

  // Create minimisation parameters
  std::vector<double> Error( NumParams, 0.1 );
  ROOT::Minuit2::MnUserParameters parInit( par, Error );
  //for( int e = 1; e < NumExponents; ++e )
    //for( int o = 0; o < NumOps; ++o )
      //parInit.SetLowerLimit( NumExponents + MELIndex(o, e), 0 );
  //parInit.SetValue(0,0.184);
  //parInit.Fix(0);
  DumpParameters( "    Initial guess:", parInit );
  static const int StrategyLevel{ 1 }; // for parameter ERRORS (see MnStrategy) 0=low, 1=medium, 2=high
  ROOT::Minuit2::MnMigrad minimizer( *this, parInit, StrategyLevel );

  // Make somewhere to store the results of the fit for each bootstrap sample
  Model ModelParams( OpNames, NumExponents, NumFiles, tMin, tMax, dof, bFactor, NSamples, NumParams+1 );
  std::vector<Common::Fold<double>> ModelCorr( NumFiles ); // correlators resulting from the fit params
  for( int f = 0; f < NumFiles; f++ )
    ModelCorr[f].resize( NSamples, Nt );

  // Make somewhere to sort the results of each fit by energy level
  std::vector<std::vector<double>> SortingHat{ static_cast<size_t>( NumExponents ) };
  for( int e = 0; e < NumExponents; ++e )
    SortingHat[e].resize( NumOps + 1 );

  const int HalfNt{ Nt - 1 }; // NB: This is half of the ORIGINAL correlator time dimension
  for( int LoopIdx = Fold::idxCentral - Skip; LoopIdx < NSamples; LoopIdx++ )
  {
    idx = ( LoopIdx >= 0 && LoopIdx < NSamples ) ? LoopIdx : Fold::idxCentral;
    ROOT::Minuit2::FunctionMinimum min = minimizer( MaxIt, Tolerance );
    ROOT::Minuit2::MnUserParameterState state = min.UserState();
    if( !state.IsValid() )
      throw std::runtime_error( "Fit " + std::to_string( LoopIdx ) + " did not converge" );
    // Throw away the first few fit results - just use them as a seed for the next fit
    if( LoopIdx >= Fold::idxCentral )
    {
      const ROOT::Minuit2::MnUserParameters &upar{ state.Parameters() };
      if( idx == Fold::idxCentral || Verbosity > 2 )
      {
        const std::string FitResult{ "    Fit result:" };
        if( Verbosity )
          std::cout << FitResult << state << "\n";
        DumpParameters( FitResult, upar );
      }
      const double ThisChiSq{ state.Fval() };
      if( idx == Fold::idxCentral )
      {
        ChiSq = ThisChiSq;
        std::cout << "    Chi^2=" << ChiSq << ", dof=" << dof << ", chi^2/dof=" << ChiSq / dof
                  << "\n\t computing statistics\n";
      }
      // Save the fit parameters for this replica, sorted by E_0
      for( int e = 0; e < NumExponents; ++e )
      {
        SortingHat[e][0] = upar.Value( e );
        for( int o = 0; o < NumOps; ++o )
          SortingHat[e][o + 1] = upar.Value( MELIndex(o, e) + NumExponents );
      }
      std::sort( SortingHat.begin(), SortingHat.end(),
                []( const std::vector<double> &l, const std::vector<double> &r )
                { return l[0] < r[0]; } );
      double * const FitParams{ ModelParams[idx] };
      for( int e = 0; e < NumExponents; ++e )
      {
        FitParams[e] = SortingHat[e][0];
        for( int o = 0; o < NumOps; ++o )
          FitParams[MELIndex(o, e) + NumExponents] = SortingHat[e][o + 1];
      }
      FitParams[MELIndex(0, NumExponents) + NumExponents] = ThisChiSq / dof;
      if( idx == Fold::idxCentral )
      {
        // Check whether energy levels are separated by minimum separation
        for( int e = 1; e < NumExponents; e++ )
        {
          double RelSep = FitParams[e] / FitParams[e - 1];
          if( ( RelSep - 1 ) < RelEnergySep || RelSep * RelEnergySep > 1 )
            throw std::runtime_error( "Fit failed energy separation criteria: "
                                     + std::to_string( 1 + RelEnergySep ) + " < E_n / E_{n-1} < "
                                     + std::to_string( 1 / RelEnergySep ) );
        }
      }
      // Save the reconstructed correlator values for this replica
      {
        double SinhCoshAdjust[NumExponents];
        for( int e = 0; e < NumExponents; ++e )
          SinhCoshAdjust[e] = std::exp( - FitParams[e] * HalfNt );
        for( int f = 0; f < NumFiles; f++ )
        {
          const Common::Parity parity{ Corr[f].parity };
          const int snk{ Corr[f].Name_.op[idxSnk] };
          const int src{ Corr[f].Name_.op[idxSrc] };
          scalar * mc{ ModelCorr[f][idx] };
          for( int t = 0; t < Nt; t++ )
          {
            double z = 0;
            for( int e = 0; e < NumExponents; ++e )
            {
              double d;
              switch( parity )
              {
                case Common::Parity::Even:
                  d = SinhCoshAdjust[e] * std::cosh( - FitParams[e] * ( t - HalfNt ) );
                  break;
                case Common::Parity::Odd:
                  d = SinhCoshAdjust[e] * std::sinh( - FitParams[e] * ( t - HalfNt ) );
                  break;
                default:
                  d = std::exp( - FitParams[e] * t );
                  break;
              }
              d *= FitParams[MELIndex(src,e)+NumExponents] * FitParams[MELIndex(snk,e)+NumExponents];
              z += d;
            }
            *mc++ = z;
          }
        }
      }
    }
  }
  std::string sModelBase{ OutputRunBase };
  sModelBase.append( 1, '.' );
  sModelBase.append( sOpNameConcat );
  ModelParams.Write( Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ) );
  //Common::SummariseBootstrap(ModelParams, sModelBase, Seed, "params" );
  for( int f = 0; f < NumFiles; f++ )
  {
    const int snk{ Corr[f].Name_.op[idxSnk] };
    const int src{ Corr[f].Name_.op[idxSrc] };
    std::string sSink{ OpNames[snk] };
    std::size_t pos = sSink.find_last_of( '_' );
    if( pos != std::string::npos )
      sSink.resize( pos );
    std::string sSrc{ OpNames[src] };
    pos = sSrc.find_last_of( '_' );
    if( pos != std::string::npos )
      sSrc.resize( pos );
    const std::string SummaryBase{ OutputRunBase + '.' + sSink + '_' + sSrc };
    if( bSaveCorr )
      ModelCorr[f].Write( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, DEF_FMT ) );
    ModelCorr[f].WriteSummary( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, TEXT_EXT ));
  }
  // Return the statistics on the fit results
  std::vector<Common::ValWithEr> Results( NumParams );
  std::vector<double> data( NSamples );
  for( int p = 0; p < NumParams; p++ ) {
    double Central = ModelParams[Fold::idxCentral][p];
    std::size_t Count{ 0 };
    for( int j = 0; j < NSamples; ++j ) {
      double d = ModelParams[j][p];
      if( std::isfinite( d ) )
        data[Count++] = d;
    }
    Results[p].Get( Central, data, Count );
  }
  return Results;
}

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"ti", CL::SwitchType::Single, nullptr},
      {"tf", CL::SwitchType::Single, nullptr},
      {"dti", CL::SwitchType::Single, "1"},
      {"dtf", CL::SwitchType::Single, "1"},
      {"sep", CL::SwitchType::Single, "0.2"},
      {"delta", CL::SwitchType::Single, "6"},
      {"skip", CL::SwitchType::Single, "10"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "0.0001"},
      {"v", CL::SwitchType::Single, "0"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"e", CL::SwitchType::Single, "1"},
      {"n", CL::SwitchType::Single, "0"},
      {"f", CL::SwitchType::Flag, nullptr},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      const int ti_start{ cl.SwitchValue<int>("ti") };
      const int tf_start{ cl.SwitchValue<int>("tf") };
      const int dti_max{ cl.SwitchValue<int>("dti") };
      const int dtf_max{ cl.SwitchValue<int>("dtf") };
      const double RelEnergySep{ cl.SwitchValue<double>("sep") };
      const int delta{ cl.SwitchValue<int>("delta") };
      const int Skip{ cl.SwitchValue<int>("skip") };
      const int MaxIterations{ cl.SwitchValue<int>("iter") }; // Max iteration count, 0=unlimited
      const double Tolerance{ cl.SwitchValue<double>("tol") }; // Actual tolerance is 10^{-3} * this
      const int Verbosity{ cl.SwitchValue<int>("v") };
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      const int NumExponents{ cl.SwitchValue<int>( "e" ) };
      const int NSamples{ cl.SwitchValue<int>("n") };
      //const std::string model{ opt.optionValue("m") };
      const bool bFactor{ cl.GotSwitch( "f" ) };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };

      if( Skip < 0 )
        throw std::invalid_argument( "Skip must be >= 0" );
      if( MaxIterations < 0 )
        throw std::invalid_argument( "MaxIterations must be >= 0" );

      bShowUsage = false;
      std::size_t i = 0;
      Common::SeedType Seed = 0;
      std::vector<std::string> OpNames;
      std::vector<Fold> Corr;
      std::cout << "Loading folded correlators\n";
      for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
      {
        std::string GroupName;
        Corr.emplace_back();
        Corr[i].Read( sFileName, GroupName, "  ", &OpNames );
        if( i == 0 )
        {
          outBaseFileName.append( Corr[0].Name_.Base );
          Seed = Corr[0].Name_.Seed;
        }
        else
        {
          static const std::string sFile{ "File " + std::to_string( i ) + " " };
          static const std::string sBad{ " doesn't match " };
          if( !Common::EqualIgnoreCase( Corr[i].Name_.Base, Corr[0].Name_.Base ) )
            throw std::runtime_error( sFile+"base "+Corr[i].Name_.Base+sBad+Corr[0].Name_.Base );
          if( !Common::EqualIgnoreCase( Corr[i].Name_.Type, Corr[0].Name_.Type ) )
            throw std::runtime_error( sFile+"type "+Corr[i].Name_.Type+sBad+Corr[0].Name_.Type );
          if( Corr[i].Name_.Seed != Corr[0].Name_.Seed )
            throw std::runtime_error(sFile+"seed "+Corr[i].Name_.SeedString+sBad+Corr[0].Name_.SeedString);
          if( !Common::EqualIgnoreCase( Corr[i].Name_.Ext, Corr[0].Name_.Ext ) )
            throw std::runtime_error( sFile+"extension "+Corr[i].Name_.Ext+sBad+Corr[0].Name_.Ext );
        }
        i++;
      }
      std::vector<std::string> OpNameOrig{ OpNames };
      if( !bFactor )
      {
        OpNames.clear();
        for( Fold & f : Corr )
        {
          using It = const typename std::vector<std::string>::iterator;
          for( int i = 0; i < 2; i++ )
          {
            const std::string Op{ OpNameOrig[f.Name_.op[i]] + "_" + pSrcSnk[i] };
            It it{ std::find( OpNames.begin(), OpNames.end(), Op ) };
            if( it == OpNames.end() )
            {
              f.Name_.op[i] = static_cast<int>( OpNames.size() );
              OpNames.emplace_back( std::move( Op ) );
            }
            else
              f.Name_.op[i] = static_cast<int>( it - OpNames.begin() );
          }
        }
      }
      std::sort( OpNameOrig.begin(), OpNameOrig.end() );
      std::string sOpNameConcat{ OpNameOrig[0] };
      for( std::size_t i = 1; i < OpNameOrig.size(); i++ )
      {
        sOpNameConcat.append( 1, '_' );
        sOpNameConcat.append( OpNameOrig[i] );
      }
      MultiExpModel m ( std::move( Corr ), OpNames, bFactor, sOpNameConcat, NumExponents, Verbosity, outBaseFileName, Seed, NSamples );
      std::string sSummaryBase{ outBaseFileName };
      sSummaryBase.append( 1, '.' );
      if( doCorr )
        sSummaryBase.append( "corr" );
      else
        sSummaryBase.append( "uncorr" );
      sSummaryBase.append( 1, '.' );
      sSummaryBase.append( sOpNameConcat );
      static const char Sep[] = " ";
      const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
      if( Common::FileExists( sFitFilename ) )
        throw std::runtime_error( "Output file " + sFitFilename + " already exists" );
      std::ofstream s;
      for( int tf = tf_start; tf - tf_start < dtf_max; tf++ )
      {
        bool bNeedHeader = true;
        for( int ti = ti_start; ti - ti_start < dti_max; ti++ )
        {
          if( tf - ti + 1 >= delta )
          {
            try
            {
              double ChiSq;
              int dof;
              auto params = m.PerformFit( doCorr, ti, tf, Skip, bSaveCorr, MaxIterations, Tolerance,
                                          RelEnergySep, ChiSq, dof );
              if( bNeedHeader )
              {
                bNeedHeader = false;
                if( !s.is_open() )
                {
                  s.open( sFitFilename );
                  s << "# Fit parameters\n# " << sSummaryBase << "\n# Seed " << Seed
                    << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
                }
                else
                {
                  // two blank lines at start of new data block
                  s << "\n" << std::endl;
                }
                // Name the data series
                s << "# [tf=" << tf << "]" << std::endl;
                // Column names, with the series value embedded in the column header (best I can do atm)
                s << "tf=" << tf << Sep << "ti";
                for( int p = 0; p < m.NumParams; p++ )
                  s << Sep << m.ParamNames[p] << Sep << m.ParamNames[p] << "_low" << Sep << m.ParamNames[p] << "_high" << Sep << m.ParamNames[p] << "_check";
                s << " ChiSq Dof ChiSqPerDof" << std::endl;
              }
              s << tf << Sep << ti;
              for( int p = 0; p < m.NumParams; p++ )
                s << Sep << params[p];
              s << Sep << ChiSq << Sep << dof << Sep << ( ChiSq / dof ) << std::endl;
            }
            catch(const std::exception &e)
            {
              std::cout << "Error: " << e.what() << "\n";
            }
          }
        }
      }
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> Bootstrap1 [Bootstrap2 ...]\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, where <options> are:\n"
    "--ti   Initial fit time\n"
    "--tf   Final   fit time\n"
    "--dti  Number of initial fit times (default 1)\n"
    "--dtf  Number of final   fit times (default 1)\n"
    "--sep  Minimum relative separation between energy levels (default 0.2)\n"
    "--delta Minimum number of timeslices in fit range (default 6)\n"
    "--skip Number of fits to throw away prior to central fit (default 10)\n"
    "--iter Max iteration count, 0 (default) = unlimited\n"
    "--tol  Tolerance of required fits * 10^3 (default 0.0001)\n"
    "-v     Verbosity, 0 = normal (default), 1=debug, 2=save covar, 3=Every iteration\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "Flags:\n"
    "-f         Factorising operators (default non-factorising)\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--help     This message\n";
  }
  return iReturn;
}
