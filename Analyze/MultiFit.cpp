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

#ifndef DEBUG_NUM
#define DEBUG_NUM 6
#endif

using Matrix = Eigen::MatrixXd; // dynamic sized matrix of complex double
using SampleD = Common::SampleD; // Bootstrap sample double
using SampleC = Common::SampleC; // Bootstrap sample std::complex<double>

// Minimisation for multi-exponential fit

class MultiExpModel : public ROOT::Minuit2::FCNBase {
public:
  static constexpr double pi{ M_PI };
  using vd_t = std::vector<double>;
  explicit MultiExpModel( SampleD &&Corr, int NumOps, const std::vector<std::string> &OpNames,
                         int NumExponents, int Verbosity, const std::string &OutputBaseName,
                         Common::SeedType Seed, int Fold, int NSamples );
  virtual ~MultiExpModel() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const vd_t & ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  // These are definitely not part of the FCNBase interface
  std::vector<Common::ValWithEr> PerformFit( bool bCorrelated, int tMin, int tMax,
                int Skip, bool bSaveCorr, int MaxIt, double Tolerance, double &ChiSq, int &dof );
  const int NumOps;
  const std::vector<std::string> &OpNames;
  const int Verbosity;
  const std::string OutputBaseName;
  const Common::SeedType Seed;
  enum ModelType { sinh = -1, exp, cosh };
  const ModelType Model;
  bool bAlternating;
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int NumExponents;
  const int NumParams;
  const SampleD Corr;
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

MultiExpModel::MultiExpModel( SampleD &&corr_, int numOps_, const std::vector<std::string> &opNames_, int numExponents_, int verbosity_, const std::string &outputBaseName_, Common::SeedType seed_, int fold_, int nSamples_ )
  : NumOps{ numOps_ },
    OpNames{ opNames_ },
    Verbosity{ verbosity_ },
    OutputBaseName{ outputBaseName_ },
    Seed{ seed_ },
    Model{ fold_ < 0 ? ModelType::sinh : fold_ > 0 ? ModelType::cosh : ModelType::exp },
    bAlternating{ abs( fold_ ) > 1 },
    NumFiles{ numOps_ * numOps_ },
    NSamples{ ( nSamples_ <= 0 || nSamples_ > corr_.NumSamples() ) ? corr_.NumSamples() : nSamples_ },
    Nt      { corr_.Nt() / NumFiles },
    NumExponents{ numExponents_ },
    NumParams{ NumExponents * ( 1 + NumOps ) },
    Corr{ std::move( corr_ ) },
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
  assert( Corr.Nt() == NumFiles * Nt && "Bug" );
  assert( tMin >= 0 && tMin < Nt && "Error: tMin invalid" );
  assert( tMax >= 0 && tMax < Nt && "Error: tMax invalid" );
  assert( NtCorr > 1 && "Error: tMin >= tMax" );
  if( bCorrelated ) {
    std::cout << "Creating covariance matrix ...\n";
    Covar.resize( Extent, Extent );
    CovarInv.resize( Extent, Extent );
    if( Verbosity )
      std::cout << "Covariance matrix is " << Covar.rows() << " x " << Covar.cols() << "\n";
  }

  // Always make variance (since it's so fast)
  if( VarianceInv.size() != Extent )
    VarianceInv.resize( Extent );

  // Make covariance
  for( int x = 0; x < Extent; ++x ) {
    const int f1{ x / NtCorr };
    const int t1{ x % NtCorr };
    const int ReadX{ Nt * f1 + t1 + tMin };
    const double ExpectedX{ Corr[SampleD::idxCentral][ReadX] };
    const int yMin{ bCorrelated ? 0 : x };
    for( int y = yMin; y <= x; ++y ) {
      const int f2{ y / NtCorr };
      const int t2{ y % NtCorr };
      const int ReadY{ Nt * f2 + t2 + tMin };
      const double ExpectedY{ Corr[SampleD::idxCentral][ReadY] };
      double z = 0;
      for( int i = 0; i < NSamples; ++i )
        z += ( Corr[i][ReadX] - ExpectedX ) * ( Corr[i][ReadY] - ExpectedY );
      z /= NSamples;
      if( bCorrelated ) {
        Covar(x, y) = z;
        if( x != y )
          Covar(y, x) = z;
      }
      if( y == x )
        VarianceInv[x] = 1. / z;
    }
  }
  if( bCorrelated ) {
    assert( Common::IsFinite( Covar ) && "Error: covariance matrix isn't finite" );
    CovarInv = Covar.inverse();
    assert( Common::IsFinite( CovarInv ) && "Error: inverse covariance matrix isn't finite" );
#ifdef DEBUG
    std::cout << "Covariance matrix is " << Covar.rows() << " x " << Covar.cols() << "\n";
    std::cout << "Inverse covariance matrix is " << CovarInv.rows() << " x " << CovarInv.cols() << "\n";
    if( Verbosity > 1 ) {
      std::cout << "Correlation matrix:\n";
      for( int i = 0; i < ( NumFiles > 1 ? 1 : 0 ) * NtCorr + DEBUG_NUM; ++i ) {
        std::cout << " " << i << ":";
        for( int j = 0; j < DEBUG_NUM; ++j )
          std::cout << " " << Covar(i, j) / ( sqrt( Covar(i, i) ) * sqrt( Covar(j, j) ) );
        std::cout << "\n";
      }
    }
#endif
    if( Verbosity > 1 ) {
      //Latan::Io::save(Covar, OutputRunBase + "Covar.h5");
      //Latan::Io::save(CovarInv, OutputRunBase + "CovarInv.h5");
      //Latan::DMat CovarProd{ Covar * CovarInv };
      //Latan::Io::save(CovarProd, OutputRunBase + "CovarProd.h5");
    }
    const double CondNumber{ Covar.norm() * CovarInv.norm() };
    const int CondDigits{ static_cast<int>( std::log10( CondNumber ) + 0.5 ) };
    if( CondDigits >= 12 )
      std::cout << "ERROR see https://en.wikipedia.org/wiki/Condition_number\n";
    std::cout << "Covariance matrix condition number=" << CondNumber
              << ", i.e. potential loss of " << CondDigits << " digits.\n";
  }
  assert( Common::IsFinite( VarianceInv ) && "Error: variance vector isn't finite" );
}

// Compute chi-squared given parameters for multi-exponential fit
  // Energy Levels        Parameters 0 ... NumExponents - 1
  //  NB: second and subsequent energy levels are deltas
  // Overlap Coefficients Parameters NumExponents ... NumExponents ( NumOps + 1 ) - 1

double MultiExpModel::operator()( const vd_t & par ) const
{
  const double * const Energy{ par.data() };
  const double * const Coeff{ Energy + NumExponents };
  // Calculate the theory errors for these model parameters
  double ModelError[Extent]; // Should happily fit on stack
  const int HalfNt{ Nt - 1 }; // NB: This is half of the ORIGINAL correlator time dimension
  for( int snk = 0; snk < NumOps; ++snk )
    for( int src = 0; src < NumOps; ++src )
    {
      const int alpha{ AlphaIndex( snk, src ) };
      // Check which fit model to use based on operator product parity
      ModelType ThisModel{ Model };
      if( bAlternating && ( snk & 1 ) != ( src & 1 ) )
      {
        if( ThisModel == cosh )
          ThisModel = sinh;
        else if( ThisModel == sinh )
          ThisModel = cosh;
      }
      for( int t = 0; t < NtCorr; ++t )
      {
        double z = 0;
        for( int e = 0; e < NumExponents; ++e )
        {
          double dThis;
          switch( ThisModel )
          {
            case cosh:
              dThis = std::exp( - Energy[e] * HalfNt ) * std::cosh( - Energy[e] * ( t + tMin - HalfNt ) );
              break;
            case sinh:
              dThis = std::exp( - Energy[e] * HalfNt ) * std::sinh( - Energy[e] * ( t + tMin - HalfNt ) );
              break;
            default:
              dThis = std::exp( - Energy[e] * ( t + tMin ) );
              break;
          }
          z += dThis * Coeff[MELIndex( src, e )] * Coeff[MELIndex( snk, e )];
        }
        const int iRead{ alpha * Nt + t + tMin };
        z -= Corr[idx][iRead];
        if( !std::isfinite( z ) )
        {
          if( Verbosity )
          {
            static int iOOB{ 0 };
            std::cout << "index " << idx << ", " << ++iOOB << "th overflow\n";
          }
          return std::numeric_limits<double>::max();
        }
        const int iWrite{ alpha * NtCorr + t };
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
                  int Skip, bool bSaveCorr, int MaxIt, double Tolerance, double &ChiSq, int &dof )
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
  if( bTimesChanged || bCorrelatedChanged ) {
    OutputRunBase = OutputBaseName;
    OutputRunBase.append( 1, '.' );
    OutputRunBase.append( bCorrelated ? "corr" : "uncorr" );
    OutputRunBase.append( 1, '_' );
    OutputRunBase.append( std::to_string( tMin ) );
    OutputRunBase.append( 1, '_' );
    OutputRunBase.append( std::to_string( tMax ) );
  }
  if( bTimesChanged ) {
    NtCorr = TMax_ - TMin_ + 1;
    Extent = NumFiles * NtCorr;
    VarianceInv.resize( 0 );
    CovarInv.resize( 0, 0 );
    Covar.resize( 0, 0 );
  }
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
    const double * const m{ Corr[SampleD::idxCentral] };
    const int tGuess{Nt / 4};
    Energy[0] = 0;
    for( int snk = 0; snk < NumOps; ++snk )
      for( int src = 0; src < NumOps; ++src ) {
        const int i{ AlphaIndex( snk, src ) * Nt + tGuess };
        double E0 = std::log( m[i] / m[i + 1] );
        Energy[0] += E0;
        if( snk == src )
          Coeff[MELIndex( snk, 0 )] = std::sqrt( std::abs( m[i] / ( std::exp( - E0 * tGuess ) ) ) );
      }
    Energy[0] /= NumFiles;
    // Now guess Higher exponents - same as LatAnalyze
    static const double MELFactor{ std::sqrt( 0.5 ) };
    for( int e = 1; e < NumExponents; ++e ) {
      Energy[e] = Energy[e - 1] + ( Energy[e - 1] - ( e > 1 ? Energy[e - 2] : 0 ) ) * 0.5;
      for( int o = 0; o < NumOps; ++o )
        Coeff[ MELIndex( o, e ) ] = Coeff[ MELIndex( o, e - 1 ) ] * MELFactor;
    }
  }

  // Create minimisation parameters
  std::vector<double> Error( NumParams, 0.1 );
  ROOT::Minuit2::MnUserParameters upar( par, Error );
  //upar.SetLowerLimit( 0, 0.1 );
  //upar.SetValue(0,0.184);
  //upar.Fix(0);
  DumpParameters( "Initial guess:", upar );
  static const int StrategyLevel{ 1 }; // for parameter ERRORS (see MnStrategy) 0=low, 1=medium, 2=high
  ROOT::Minuit2::MnMigrad minimizer( *this, upar, StrategyLevel );

  // Make somewhere to store the results of the fit for each bootstrap sample
  SampleD ModelParams( NSamples, NumParams ); // Model parameters determined by fitting function
  std::vector<std::vector<SampleC>> ModelCorr( NumOps ); // correlators resulting from the fit params
  for( int iSnk = 0; iSnk < NumOps; iSnk++ )
  {
    ModelCorr[iSnk].resize( NumOps );
    for( int iSrc = 0; iSrc < NumOps; iSrc++ )
      ModelCorr[iSnk][iSrc].resize( NSamples, Nt ); // Complex component will be zero
  }

  // Make somewhere to sort the results of each fit by energy level
  std::vector<std::vector<double>> SortingHat{ static_cast<size_t>( NumExponents ) };
  for( int e = 0; e < NumExponents; ++e )
    SortingHat[e].resize( NumOps + 1 );

  const int HalfNt{ Nt - 1 }; // NB: This is half of the ORIGINAL correlator time dimension
  for( int LoopIdx = SampleD::idxCentral - Skip; LoopIdx < NSamples; LoopIdx++ )
  {
    idx = ( LoopIdx >= 0 && LoopIdx < NSamples ) ? LoopIdx : SampleD::idxCentral;
    ROOT::Minuit2::FunctionMinimum min = minimizer( MaxIt, Tolerance );
    ROOT::Minuit2::MnUserParameterState state = min.UserState();
    if( !state.IsValid() )
      throw std::runtime_error( "Fit " + std::to_string( LoopIdx ) + " did not converge" );
    // Throw away the first few fit results - just use them as a seed for the next fit
    if( LoopIdx >= SampleD::idxCentral )
    {
      upar = state.Parameters();
      if( idx == SampleD::idxCentral || Verbosity > 2 )
      {
        const std::string FitResult{ "Fit result:" };
        if( Verbosity )
          std::cout << FitResult << state << "\n";
        DumpParameters( FitResult, upar );
      }
      if( idx == SampleD::idxCentral )
      {
        ChiSq = state.Fval();
        dof = Extent;
        //if( bCorrelated ) dof *= dof;
        dof -= NumParams;
        std::cout << "Chi^2=" << ChiSq << ", dof=" << dof << ", chi^2/dof=" << ChiSq / dof
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
      // Save the reconstructed correlator values for this replica
      for( int snk = 0; snk < NumOps; ++snk )
        for( int src = 0; src < NumOps; ++src )
        {
          std::complex<double> * const mc{ ModelCorr[snk][src][idx] };
          // Check which fit model to use based on operator product parity
          ModelType ThisModel{ Model };
          if( bAlternating && ( snk & 1 ) != ( src & 1 ) )
          {
            if( ThisModel == cosh )
              ThisModel = sinh;
            else if( ThisModel == sinh )
              ThisModel = cosh;
          }
          for( int t = 0; t < Nt; ++t )
          {
            double z = 0;
            for( int e = 0; e < NumExponents; ++e )
            {
              switch( ThisModel )
              {
                case exp:
                  z += FitParams[MELIndex( src, e ) + NumExponents]
                     * FitParams[MELIndex( snk, e ) + NumExponents]
                     * std::exp( - FitParams[e] * t );
                  break;
                case cosh:
                  z += FitParams[MELIndex( src, e ) + NumExponents]
                     * FitParams[MELIndex( snk, e ) + NumExponents]
                     * std::exp( - FitParams[e] * HalfNt )
                     * std::cosh(- FitParams[e] * ( t - HalfNt ) );
                  break;
                case sinh:
                  z += FitParams[MELIndex( src, e ) + NumExponents]
                     * FitParams[MELIndex( snk, e ) + NumExponents]
                     * std::exp( - FitParams[e] * HalfNt )
                     * std::sinh(- FitParams[e] * ( t - HalfNt ) );
                  break;
              }
            }
            mc[t] = z;
          }
        }
    }
  }
  std::string sModelBase{ OutputRunBase };
  sModelBase.append( 1, '.' );
  sModelBase.append( OpNames[0] );
  for( std::size_t i = 1; i < OpNames.size(); i++ ) {
    sModelBase.append( 1, '_' );
    sModelBase.append( OpNames[i] );
  }
  ModelParams.Write( Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ),
                     Common::sModel.c_str() );
  //Common::SummariseBootstrap(ModelParams, sModelBase, Seed, "params" );
  for( int snk = 0; snk < NumOps; ++snk )
    for( int src = 0; src < NumOps; ++src )
    {
      const std::string SummaryBase{ OutputRunBase + '.' + OpNames[snk] + '_' + OpNames[src] };
      Common::SampleC &c{ ModelCorr[snk][src] };
      if( bSaveCorr )
        c.Write( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, DEF_FMT ) );
      Common::SummariseBootstrapCorr( c, SummaryBase, Seed );
    }
  // Return the statistics on the fit results
  std::vector<Common::ValWithEr> Results( NumParams );
  std::vector<double> data( NSamples );
  for( int p = 0; p < NumParams; p++ ) {
    double Central = ModelParams[SampleD::idxCentral][p];
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
      {"skip", CL::SwitchType::Single, "10"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "0.1"},
      {"s", CL::SwitchType::Single, "0"},
      {"v", CL::SwitchType::Single, "0"},
      {"o", CL::SwitchType::Single, "" },
      {"f", CL::SwitchType::Single, "0"},
      {"e", CL::SwitchType::Single, "0"},
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      const int ti{ cl.SwitchValue<int>("ti") };
      const int tf{ cl.SwitchValue<int>("tf") };
      const int dti_max{ cl.SwitchValue<int>("dti") };
      const int dtf_max{ cl.SwitchValue<int>("dtf") };
      const int Skip{ cl.SwitchValue<int>("skip") };
      const int MaxIterations{ cl.SwitchValue<int>("iter") }; // Max iteration count, 0=unlimited
      const double Tolerance{ cl.SwitchValue<double>("tol") }; // Actual tolerance is 10^{-3} * this
      const int shift{ cl.SwitchValue<int>("s") };
      const int fold{ cl.SwitchValue<int>("f") };
      int NumExponents{ cl.SwitchValue<int>( "e" ) };
      const int NSamples{ cl.SwitchValue<int>("n") };
      //const std::string model{ opt.optionValue("m") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      const int Verbosity{ cl.SwitchValue<int>("v") };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };

      int NumOps = static_cast<int>( std::sqrt( static_cast<double>( NumFiles ) ) + 0.5 );
      if( NumOps * NumOps != NumFiles )
        throw std::invalid_argument( "Number of files should be a perfect square" );
      if( Skip < 0 )
        throw std::invalid_argument( "Skip must be >= 0" );
      if( MaxIterations < 0 )
        throw std::invalid_argument( "MaxIterations must be >= 0" );

      std::size_t i = 0;
      Common::SeedType Seed = 0;
      std::vector<std::string> OpNames;
      std::vector<Common::FileNameAtt> FileNames;
      for( const std::string &sFileName : cl.Args )
      {
        FileNames.emplace_back( sFileName, OpNames );
        if( i == 0 )
        {
          outBaseFileName.append( FileNames[0].Base );
          Seed = FileNames[0].Seed;
        }
        else
        {
          static const std::string sFile{ "File " };
          static const std::string sBad{ " doesn't match " };
          if( !Common::EqualIgnoreCase( FileNames[i].Base, FileNames[0].Base ) )
            throw std::runtime_error( sFile + std::to_string( i ) + " base " + FileNames[i].Base + sBad + FileNames[0].Base );
          if( !Common::EqualIgnoreCase( FileNames[i].Type, FileNames[0].Type ) )
            throw std::runtime_error( sFile + std::to_string( i ) + " type " + FileNames[i].Type + sBad + FileNames[0].Type );
          if( FileNames[i].Seed != FileNames[0].Seed )
            throw std::runtime_error( sFile + std::to_string( i ) + " seed " + std::to_string( FileNames[i].Seed ) + sBad + std::to_string( FileNames[0].Seed ) );
          if( !Common::EqualIgnoreCase( FileNames[i].Ext, FileNames[0].Ext ) )
            throw std::runtime_error( sFile + std::to_string( i ) + " extension " + FileNames[i].Ext + sBad + FileNames[0].Ext );
        }
        i++;
      }
      if( OpNames.size() != NumOps )
        throw std::runtime_error( std::to_string( OpNames.size() ) + " operators provided, but " + std::to_string( NumOps ) + " expected for " + std::to_string( NumFiles ) + " files" );
      for( int snk = 0; snk < NumOps; ++snk )
        for( int src = 0; src < NumOps; ++src )
        {
          Common::FileNameAtt &f{ FileNames[snk * NumOps + src] };
          if( f.op[0] != src || f.op[1] != snk )
            throw std::runtime_error( "Warning: Operator order should be sink-major, source minor" );
        }

      bShowUsage = false;
      if( NumExponents == 0 )
        NumExponents = NumOps;
      SampleD sd;
      {
        SampleC sc;
        sc.Read( cl.Args, fold, shift, NumOps, "  " );
        const int NumSamples{ sc.NumSamples() };
        const int Nt{ sc.Nt() };
        sd.resize( NumSamples, Nt );
        const std::complex<double> * pSrc = sc[SampleC::idxAux];
        double * pDst = sd[SampleD::idxAux];
        for( int i = SampleC::idxAux; i != NumSamples; i++ )
          for( int t = 0; t != Nt; t++ )
            *pDst++ = (pSrc++)->real();
      }
      MultiExpModel m ( std::move( sd ), NumOps, OpNames, NumExponents, Verbosity, outBaseFileName, Seed, fold, NSamples );
      std::string sSummaryBase{ outBaseFileName };
      sSummaryBase.append( 1, '.' );
      if( doCorr )
        sSummaryBase.append( "corr" );
      else
        sSummaryBase.append( "uncorr" );
      sSummaryBase.append( 1, '.' );
      sSummaryBase.append( OpNames[0] );
      for( int op = 1; op < NumOps; op++ ) {
        sSummaryBase.append( 1, '_' );
        sSummaryBase.append( OpNames[op] );
      }
      static const char Sep[] = " ";
      const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
      if( Common::FileExists( sFitFilename ) )
        throw std::runtime_error( "Output file " + sFitFilename + " already exists" );
      std::ofstream s( sFitFilename );
      s << "# Fit parameters\n# " << sSummaryBase << "\n# Seed " << Seed
        << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
      for( int dtf = 0; dtf < dtf_max; dtf++ ) {
        // two blank lines at start of new data block
        if( dtf )
          s << "\n" << std::endl;
        // Name the data series
        s << "# [tf=" << ( tf + dtf ) << "]" << std::endl;
        // Column names, with the series value embedded in the column header (best I can do atm)
        s << "tf=" << ( tf + dtf ) << Sep << "ti";
        for( int p = 0; p < m.NumParams; p++ )
          s << Sep << m.ParamNames[p] << Sep << m.ParamNames[p] << "ErLow" << Sep << m.ParamNames[p] << "ErHigh" << Sep << m.ParamNames[p] << "Check";
        s << " ChiSq Dof ChiSqPerDof" << std::endl;
        for( int dti = 0; dti < dti_max; dti++ )
        {
          if( ti + dti < tf + dtf )
          {
            double ChiSq;
            int dof;
            auto params = m.PerformFit( doCorr, ti + dti, tf + dtf, Skip, bSaveCorr, MaxIterations, Tolerance, ChiSq, dof );
            s << ( tf + dtf ) << Sep << ( ti + dti );
            for( int p = 0; p < m.NumParams; p++ )
              s << Sep << params[p];
            s << Sep << ChiSq << Sep << dof << Sep << ( ChiSq / dof ) << std::endl;
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
    "--dti  Number of final   fit times (default 1)\n"
    "--skip Number of fits to throw away prior to central fit (default 10)\n"
    "--iter Max iteration count, 0 (default) = unlimited\n"
    "--tol  Tolerance of required fits * 10^3 (default 0.1)\n"
    "-s     time Shift (default 0)\n"
    "-v     Verbosity, 0 = normal (default), 1=debug, 2=save covar, 3=Every iteration\n"
    "-o     Output base filename\n"
    "-f     Fold the correlator, 0=don't fold (default), 1=+parity, -1=-parity\n"
    "-e     number of Exponents (default same as number of operators)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--help     This message\n";
  }
  return iReturn;
}
