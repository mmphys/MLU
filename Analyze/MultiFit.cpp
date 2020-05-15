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

// Forward declaration of fitter for multi-exponential fit
class Fitter;

// Several of these will be running at the same time on different threads during a fit
class MultiExpModel : public ROOT::Minuit2::FCNBase
{
public:
  static constexpr int StrategyLevel{ 1 }; // for parameter ERRORS (see MnStrategy) 0=low, 1=medium, 2=high
  const Fitter &parent;
  // These variables change on each iteration
  int idx;
  Matrix Covar;
  Matrix CovarInv;
  std::vector<double> VarianceInv;
  std::vector<std::vector<double>> SortingHat;
protected:
  const bool bFreezeCovar;
  const bool bCorrelated;
  std::vector<Common::Fold<double>> &ModelCorr; // Fill this with data for model correlator
  const ROOT::Minuit2::MnUserParameters &parInitialGuess;
  // Helper functions
  void MakeCovar( void );
public:
  MultiExpModel(const Fitter &fitter, const ROOT::Minuit2::MnUserParameters &parInitialGuess,
                std::vector<Common::Fold<double>> &ModelCorr, bool bFreezeCovar, bool bCorrelated);
  void UpdateGuess( ROOT::Minuit2::MnUserParameters &parGuess, int MaxIt, double Tolerance );
  void FitOne( const int idx, Model &ModelParams, const std::string &SaveCorMatFileName,
               int MaxIt, double Tolerance, double RelEnergySep, int dof, double &ChiSq );
// These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const std::vector<double> &ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
};

class Fitter
{
public:
  //static constexpr double pi{ M_PI };
  const int NumOps;
  const std::vector<std::string> &OpNames;
  const bool bFactor;
  const std::string &sOpNameConcat;
  const int Verbosity;
  const std::string OutputBaseName;
  const Common::SeedType Seed;
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int NumExponents;
  const int NumParams;
  const std::vector<Fold> Corr;
  const std::vector<std::string> ParamNames;

  // These variables set once at the start of each fit
  int tMin = -1;
  int tMax = -1;
  int NtCorr;
  int Extent;
  bool bCorrelated;

  // Helper functions
  //inline int AlphaIndex( int snk, int src) const { return snk * NumOps + src; };
  std::vector<std::string> MakeParamNames() const;
  inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; };
  void DumpParameters( const std::string &msg, const ROOT::Minuit2::MnUserParameters &par ) const;

public:
  explicit Fitter( std::vector<Fold> &&Corr, const std::vector<std::string> &OpNames,
                         bool bFactor, const std::string &sOpNameConcat,
                         int NumExponents, int Verbosity, const std::string &OutputBaseName,
                         Common::SeedType Seed, int NSamples );
  virtual ~Fitter() {}
  std::vector<Common::ValWithEr<scalar>>
    PerformFit(bool bCorrelated, bool bFreezeCovar, int tMin, int tMax, int Skip, bool bSaveCorr,
          bool bSaveCMat, int MaxIt, double Tolerance, double RelEnergySep, double &ChiSq, int &dof );
};

MultiExpModel::MultiExpModel(const Fitter &fitter_, const ROOT::Minuit2::MnUserParameters &parInitialGuess_,
                             std::vector<Common::Fold<double>> &ModelCorr_, bool bFreezeCovar_, bool bCorrelated_)
: parent{ fitter_ },
  idx{ Fold::idxCentral },
  SortingHat{ static_cast<size_t>( fitter_.NumExponents ) },
  bFreezeCovar{bFreezeCovar_},
  bCorrelated{bCorrelated_},
  ModelCorr{ModelCorr_},
  parInitialGuess{ parInitialGuess_ }
{
  VarianceInv.resize( parent.Extent );
  if( bCorrelated )
  {
    CovarInv.resize( parent.Extent, parent.Extent );
    Covar.resize( parent.Extent, parent.Extent );
  }
  if( bFreezeCovar )
    MakeCovar();
  // Make somewhere to sort the results of each fit by energy level
  for( int e = 0; e < parent.NumExponents; ++e )
    SortingHat[e].resize( parent.NumOps + 1 );
}

// Make the covariance matrix, for only the timeslices we are interested in
void MultiExpModel::MakeCovar( void )
{
  const int NumBoot{ parent.Corr[0].NumSamples() };
  // Make covariance
  const scalar * CentralX = nullptr;
  const scalar * CentralY = nullptr;
  const scalar * ReplicaX = nullptr;
  const scalar * ReplicaY = nullptr;
  for( int x = 0; x < parent.Extent; x++ )
  {
    const int t1{ x % parent.NtCorr };
    if( !t1 )
    {
      CentralX = parent.Corr[x / parent.NtCorr][idx] + parent.tMin;
      ReplicaX = parent.Corr[x / parent.NtCorr][ 0 ] + parent.tMin;
    }
    for( int y = bCorrelated ? 0 : x; y <= x; y++ )
    {
      const int t2{ y % parent.NtCorr };
      if( !t2 )
      {
        CentralY = parent.Corr[y / parent.NtCorr][idx] + parent.tMin;
        ReplicaY = parent.Corr[y / parent.NtCorr][ 0 ] + parent.tMin;
      }
      double z = 0;
      const scalar * DataX{ ReplicaX };
      const scalar * DataY{ ReplicaY };
      for( int i = 0; i < NumBoot; i++ )
      {
        z += ( DataX[t1] - CentralX[t1] ) * ( DataY[t2] - CentralY[t2] );
        DataX += parent.Nt;
        DataY += parent.Nt;
      }
      z /= NumBoot;
      if( bCorrelated )
      {
        Covar( x, y ) = z;
        if( x != y )
          Covar( y, x ) = z;
      }
      if( y == x )
        VarianceInv[x] = 1. / sqrt( z );
    }
  }
  if( bCorrelated )
  {
    // Turn covariance matrix into correlation matrix
    for( int x = 0; x < parent.Extent; x++ )
      for( int y = 0; y <= x; y++ )
      {
        if( y == x )
          Covar( x, y ) = 1.;
        else
        {
          double z = Covar( x, y ) * VarianceInv[x] * VarianceInv[y];
          Covar( x, y ) = z;
          Covar( y, x ) = z;
        }
      }
    if( !Common::IsFinite( Covar ) )
      throw std::runtime_error( "Covariance matrix isn't finite" );
    CovarInv = Covar.inverse();
    if( !Common::IsFinite( CovarInv ) )
      throw std::runtime_error( "Inverse covariance matrix isn't finite" );
    if( idx == Fold::idxCentral )
    {
      const double CondNumber{ Covar.norm() * CovarInv.norm() };
      const int CondDigits{ static_cast<int>( std::log10( CondNumber ) + 0.5 ) };
      if( CondDigits >= 12 )
        std::cout << "    WARNING see https://en.wikipedia.org/wiki/Condition_number\n";
      std::cout << "    Covariance matrix condition number " << CondNumber
                << ", potential loss of " << CondDigits << " digits\n";
    }
  }
  if( !Common::IsFinite( VarianceInv ) )
    throw std::runtime_error( "Variance vector isn't finite" );
}

// Perform an uncorrelated fit on the central timeslice to update parameters guessed
void MultiExpModel::UpdateGuess( ROOT::Minuit2::MnUserParameters &parGuess, int MaxIt, double Tolerance )
{
  if( !bFreezeCovar )
    MakeCovar();
  ROOT::Minuit2::MnMigrad minimizer( *this, parInitialGuess, StrategyLevel );
  ROOT::Minuit2::FunctionMinimum min = minimizer( MaxIt, Tolerance );
  ROOT::Minuit2::MnUserParameterState state = min.UserState();
  if( !state.IsValid() )
    throw std::runtime_error( "Uncorrelated fit on central replica did not converge" );
  // Save the parameters and errors as the updated guess
  const ROOT::Minuit2::MnUserParameters &upar{ state.Parameters() };
  for( int i = 0; i < parent.NumParams; ++i )
  {
    parGuess.SetValue( i, upar.Value( i ) );
    parGuess.SetError( i, upar.Error( i ) );
  }
}

// Perform a fit. Only central fit updates ChiSq
void MultiExpModel::FitOne(const int idx_, Model &ModelParams, const std::string &SaveCorMatFileName,
                           int MaxIt, double Tolerance, double RelEnergySep, int dof, double &ChiSq)
{
  idx = idx_;
  // Re-use correlation matrix from uncorrelated fit, or from central replica if frozen
  if( !bFreezeCovar )
    MakeCovar();

  // Save correlation matrix for central replica if requested
  // NB: Do this here so that only one thread writes the file
  if( bCorrelated && !SaveCorMatFileName.empty() )
  {
    const std::string Sep{ " " };
    const std::string NewLine{ "\n" };
    std::ofstream s{ SaveCorMatFileName };
    s << "# Correlation matrix\n# Files: " << parent.NumFiles << "\n# NtCorr: " << parent.NtCorr
      << "\n# gnuplot: plot '" << SaveCorMatFileName
      << "' matrix columnheaders rowheaders with image pixels" << NewLine << parent.Extent;
    for( int f = 0; f < parent.NumFiles; f++ )
      for( int t = 0; t < parent.NtCorr; t++ )
        s << Sep << parent.OpNames[f] << t + parent.tMin;
    s << NewLine;
    for( int f = 0; f < parent.NumFiles; f++ )
      for( int t = 0; t < parent.NtCorr; t++ )
      {
        int i = f * parent.NtCorr + t;
        s << parent.OpNames[f] << t + parent.tMin;
        for( int j = 0; j < parent.Extent; j++ )
          s << Sep << Covar( i, j );
        s << NewLine;
      }
  }
  ROOT::Minuit2::MnMigrad minimizer( *this, parInitialGuess, StrategyLevel );
  ROOT::Minuit2::FunctionMinimum min = minimizer( MaxIt, Tolerance );
  ROOT::Minuit2::MnUserParameterState state = min.UserState();
  if( !state.IsValid() )
    throw std::runtime_error( "Fit " + std::to_string( idx ) + " did not converge" );
  const ROOT::Minuit2::MnUserParameters &upar{ state.Parameters() };
  if( idx == Fold::idxCentral || parent.Verbosity > 2 )
  {
    const std::string FitResult{ std::string(4, ' ') + (bCorrelated ? "C" : "Unc" ) + "orrelated fit result:" };
    if( parent.Verbosity )
      std::cout << FitResult << state << "\n";
    parent.DumpParameters( FitResult, upar );
#ifdef DEBUG_COMPARE_WITH_PREVIOUS
    std::vector<double> vd;
    Common::ReadArray( vd,
                      "/Users/s1786208/data/Study1/C1/GFWall/fit/test_thaw/"
                      "h1_l_p_0_0_0.corr_6_15.g5_gT5.model.4147798751.h5",
                      "data_C" );
    std::cout << "    Ratios (new/old):\n";
    for( int i = 0; i < 10; i++ )
      std::cout << "\t" << i << ": " << ( upar.Value( i ) / vd[i] ) << "\n";
    std::cout << "Set breakpoint here\n";
#endif
  }
  const double ThisChiSq{ state.Fval() };
  if( idx == Fold::idxCentral )
  {
    ChiSq = ThisChiSq;
    std::cout << "    Chi^2=" << ThisChiSq << ", dof=" << dof << ", chi^2/dof=" << ThisChiSq / dof << "\n";
  }
  // Save the fit parameters for this replica, sorted by E_0
  for( int e = 0; e < parent.NumExponents; ++e )
  {
    SortingHat[e][0] = upar.Value( e );
    for( int o = 0; o < parent.NumOps; ++o )
      SortingHat[e][o + 1] = upar.Value( parent.MELIndex(o, e) + parent.NumExponents );
  }
  std::sort( SortingHat.begin(), SortingHat.end(),
            []( const std::vector<double> &l, const std::vector<double> &r )
            { return l[0] < r[0]; } );
  double * const FitParams{ ModelParams[idx] };
  for( int e = 0; e < parent.NumExponents; ++e )
  {
    FitParams[e] = SortingHat[e][0];
    for( int o = 0; o < parent.NumOps; ++o )
      FitParams[parent.MELIndex(o, e) + parent.NumExponents] = SortingHat[e][o + 1];
  }
  FitParams[parent.MELIndex(0, parent.NumExponents) + parent.NumExponents] = ThisChiSq / dof;
  if( idx == Fold::idxCentral )
  {
    // Check whether energy levels are separated by minimum separation
    for( int e = 1; e < parent.NumExponents; e++ )
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
    const int HalfNt{ parent.Nt - 1 }; // NB: This is half of the ORIGINAL correlator time dimension
    double SinhCoshAdjust[parent.NumExponents];
    for( int e = 0; e < parent.NumExponents; ++e )
      SinhCoshAdjust[e] = std::exp( - FitParams[e] * HalfNt );
    for( int f = 0; f < parent.NumFiles; f++ )
    {
      const Common::Parity parity{ parent.Corr[f].parity };
      const int snk{ parent.Corr[f].Name_.op[idxSnk] };
      const int src{ parent.Corr[f].Name_.op[idxSrc] };
      scalar * mc{ ModelCorr[f][idx] };
      for( int t = 0; t < parent.Nt; t++ )
      {
        double z = 0;
        for( int e = 0; e < parent.NumExponents; ++e )
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
          d *=   FitParams[ parent.MELIndex( src, e ) + parent.NumExponents ]
               * FitParams[ parent.MELIndex( snk, e ) + parent.NumExponents ];
          z += d;
        }
        *mc++ = z;
      }
    }
  }
}

// Compute chi-squared given parameters for multi-exponential fit
  // Energy Levels        Parameters 0 ... NumExponents - 1
  //  NB: second and subsequent energy levels are deltas
  // Overlap Coefficients Parameters NumExponents ... NumExponents ( NumOps + 1 ) - 1

double MultiExpModel::operator()( const std::vector<double> & par ) const
{
  const double * const Energy{ par.data() };
  const double * const Coeff{ Energy + parent.NumExponents };
  // Calculate the theory errors for these model parameters
  const int HalfNt{ parent.Nt - 1 }; // NB: This is half of the ORIGINAL correlator time dimension
  double SinhCoshAdjust[parent.NumExponents];
  for( int e = 0; e < parent.NumExponents; ++e )
    SinhCoshAdjust[e] = std::exp( - Energy[e] * HalfNt );
  double ModelError[parent.Extent]; // Should happily fit on stack
  for( int f = 0; f < parent.NumFiles; f++ )
  {
    const int snk{ parent.Corr[f].Name_.op[idxSnk] };
    const int src{ parent.Corr[f].Name_.op[idxSrc] };
    const scalar * CorrData = parent.Corr[f][idx] + parent.tMin;
    for( int t = 0; t < parent.NtCorr; t++ )
    {
      double z = 0;
      for( int e = 0; e < parent.NumExponents; ++e )
      {
        double d;
        switch( parent.Corr[f].parity )
        {
          case Common::Parity::Even:
            d = SinhCoshAdjust[e] * std::cosh( - Energy[e] * ( parent.tMin + t - HalfNt ) );
            break;
          case Common::Parity::Odd:
            d = SinhCoshAdjust[e] * std::sinh( - Energy[e] * ( parent.tMin + t - HalfNt ) );
            break;
          default:
            d = std::exp( - Energy[e] * ( parent.tMin + t ) );
            break;
        }
        z += d * Coeff[parent.MELIndex( src, e )] * Coeff[parent.MELIndex( snk, e )];
      }
      z -= *CorrData++;
      if( !std::isfinite( z ) )
      {
        if( parent.Verbosity )
        {
          std::cout << "Bootstrap replica " << idx << ", timeslice " << t << " overflow\n";
        }
        return std::numeric_limits<double>::max();
      }
      const int iWrite{ f * parent.NtCorr + t };
      ModelError[iWrite] = z * VarianceInv[iWrite];
    }
  }
  // The inverse of a symmetric matrix is also symmetric, so only calculate half the matrix
  double chi2 = 0;
  for( int i = 0; i < parent.Extent; ++i )
  {
    if( !bCorrelated )
      chi2 += ModelError[i] * ModelError[i];
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
  if( parent.Verbosity > 2 )
    std::cout << "index " << idx << ", chi^2=" << chi2 << ", E0=" << Energy[0] << ", E1=" << Energy[1] << "\n";
  if( chi2 < 0 )
  {
    //static int count{ 0 };
    //std::cout << "Chi^2 < 0" << ++count << "\n";
    throw std::runtime_error( "Chi^2 < 0 on replica " + std::to_string( idx ) );
  }
  return chi2;
}

Fitter::Fitter( std::vector<Fold> &&Corr_, const std::vector<std::string> &opNames_,
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

std::vector<std::string> Fitter::MakeParamNames() const
{
  assert( NumExponents > 0 && "Number of exponents in the fit should be positive" );
  std::vector<std::string> Names( NumParams );
  for( int e = 0; e < NumExponents; e++ ) {
    Names[e] = "E" + std::to_string( e );
    for( int op = 0; op < NumOps; op++ )
      Names[NumExponents + MELIndex( op, e )] = OpNames[op] + std::to_string( e );
  }
  return Names;
}

void Fitter::DumpParameters( const std::string &msg, const ROOT::Minuit2::MnUserParameters &par ) const
{
  static const char NewLine[] = "\n";
  static const char Tab[] = "\t";
  if( !msg.empty() )
    std::cout << msg << NewLine;
  std::size_t FieldLen{ 0 };
  for( int p = 0; p < NumParams; p++ )
  {
    std::size_t ThisLen = ParamNames[p].length();
    if( FieldLen < ThisLen )
      FieldLen = ThisLen;
  }
  FieldLen += 4;
  for( int p = 0; p < NumParams; p++ )
  {
    std::string sep( FieldLen - ParamNames[p].length(), ' ' );
    std::cout << Tab << ParamNames[p] << sep << par.Value( p ) << " +/- " << par.Error( p ) << NewLine;
  }
}

// Perform a fit
std::vector<Common::ValWithEr<scalar>> Fitter::PerformFit( bool Bcorrelated_, bool bFreezeCovar,
                    int TMin_, int TMax_, int Skip, bool bSaveCorr, bool bSaveCMat, int MaxIt,
                    double Tolerance, double RelEnergySep, double &ChiSq, int &dof )
{
  {
    std::stringstream ss;
    ss << "Performing " << ( Bcorrelated_ ? "" : "un" )
       << "correlated fit on timeslices " << TMin_ << " to " << TMax_
       << " using " << omp_get_max_threads() << " Open MP threads";
    const std::string &sMsg{ ss.str() };
    std::cout << std::string( sMsg.length(), '=' ) << "\n" << sMsg << "\n";
  }
  bCorrelated = Bcorrelated_;
  tMin = TMin_;
  tMax = TMax_;
  NtCorr = TMax_ - TMin_ + 1;
  Extent = NumFiles * NtCorr;
  dof = Extent;
  // If there's only one file, then there's no info for the operators to factorise
  dof -= NumFiles == 1 ? NumExponents * 2 : NumParams;

  // Make somewhere to store the results of the fit for each bootstrap sample
  Model ModelParams( OpNames, NumExponents, NumFiles, tMin, tMax, dof, bFactor, bFreezeCovar, NSamples, NumParams+1 );
  ModelParams.Seed_ = Corr[0].Seed_;
  ModelParams.SeedMachine_ = Corr[0].SeedMachine_;
  {
    std::vector<std::string> ColNames{ ParamNames };
    ColNames.push_back( "ChiSqPerDof" );
    ModelParams.SetColumnNames( ColNames );
  }
  for( const Fold &f : Corr )
    ModelParams.FileList.emplace_back( f.Name_.Filename );

  // See whether this fit already exists
  bool bPerformFit{ true };
  std::string OutputRunBase{ OutputBaseName };
  OutputRunBase.append( 1, '.' );
  OutputRunBase.append( bCorrelated ? "corr" : "uncorr" );
  OutputRunBase.append( 1, '_' );
  OutputRunBase.append( std::to_string( tMin ) );
  OutputRunBase.append( 1, '_' );
  OutputRunBase.append( std::to_string( tMax ) );
  std::string sModelBase{ OutputRunBase };
  sModelBase.append( 1, '.' );
  sModelBase.append( sOpNameConcat );
  const std::string ModelFileName{ Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ) };
  if( Common::FileExists( ModelFileName ) )
  {
    Model PreBuilt;
    PreBuilt.Read( ModelFileName, "Pre-built: " );
    if( !PreBuilt.Compatible( NumExponents, dof, tMin, tMax, OpNames ) )
      throw std::runtime_error( "Pre-existing fits not compatible with parameters from this run" );
    bPerformFit = PreBuilt.NewParamsMorePrecise( bFreezeCovar, NSamples );
    if( !bPerformFit )
    {
      ChiSq = dof * PreBuilt.getSummaryData()[NumParams].Central; // Last summary has chi squared per dof
      ModelParams = std::move( PreBuilt );
    }
    else
      std::cout << "Overwriting\n";
  }

  if( bPerformFit )
  {
    // If we're saving the correlation matrix, it should have a similar name to the model
    const std::string SaveCorrMatrixFileName{ bCorrelated && bSaveCMat
                            ? Common::MakeFilename( sModelBase, Common::sCormat, Seed, TEXT_EXT ) : "" };
    // Make somewhere to hold the correlators corresponding to the fitted model
    std::vector<Common::Fold<double>> ModelCorr( NumFiles ); // correlators resulting from the fit params
    for( int f = 0; f < NumFiles; f++ )
    {
      ModelCorr[f].resize( NSamples, Nt );
      ModelCorr[f].Seed_ = ModelParams.Seed_;
      ModelCorr[f].SeedMachine_ = ModelParams.SeedMachine_;
    }
    // Make initial guesses for the parameters
    // For each Exponent, I need the delta_E + a constant for each operator
    std::vector<double> par( NumParams );
    std::vector<double> Error( NumParams ); // I actually have no idea what the parameter errors are
    static constexpr double ErrorFactor = 0.1;
    {
      // Take a starting guess for the parameters - same as LatAnalyze
      double * const Energy{ par.data() };
      double * const Coeff{ Energy + NumExponents };
      double * const EnergyErr{ Error.data() };
      double * const CoeffErr{ EnergyErr + NumExponents };
      const int tGuess{ Nt / 4 };
      Energy[0] = 0;
      for( int f = 0; f < NumFiles; f++ )
      {
        //const int i{ AlphaIndex( snk, src ) * Nt + tGuess };
        const double * const c{ Corr[f][Fold::idxCentral] };
        double E0 = std::log( c[tGuess] / c[tGuess + 1] );
        Energy[0] += E0;
        double MELGuess = std::sqrt( std::abs( c[tGuess] / ( std::exp( - E0 * tGuess ) ) ) );
        int i = MELIndex( Corr[f].Name_.op[idxSnk], 0 );
        Coeff[i] = MELGuess;
        CoeffErr[i] = MELGuess * ErrorFactor;
        i = MELIndex( Corr[f].Name_.op[idxSrc], 0 );
        Coeff[i] = MELGuess;
        CoeffErr[i] = MELGuess * ErrorFactor;
      }
      Energy[0] /= NumFiles;
      EnergyErr[0] = Energy[0] * ErrorFactor;
      // Now guess Higher exponents - same as LatAnalyze
      //static const double MELFactor{ std::sqrt( 0.5 ) };
      static const double MELFactor{ std::sqrt( 2 ) };
      for( int e = 1; e < NumExponents; ++e )
      {
        Energy[e] = Energy[e - 1] + ( Energy[e - 1] - ( e > 1 ? Energy[e - 2] : 0 ) ) * 0.5;
        EnergyErr[e] = Energy[e] * ErrorFactor;
        for( int o = 0; o < NumOps; ++o )
        {
          int i = MELIndex( o, e );
          Coeff[i] = Coeff[ MELIndex( o, e - 1 ) ] * MELFactor;
          CoeffErr[i] = Coeff[i] * ErrorFactor;
        }
      }
    }
    
    // Create minimisation parameters
    ROOT::Minuit2::MnUserParameters parInit( par, Error );
    Error.clear();
    //for( int e = 1; e < NumExponents; ++e )
    //for( int o = 0; o < NumOps; ++o )
    //parInit.SetLowerLimit( NumExponents + MELIndex(o, e), 0 );
    //parInit.SetValue(0,0.184);
    //parInit.Fix(0);
    DumpParameters( "    Initial guess:", parInit );
    if( bCorrelated )
    {
      // Perform an uncorrelated fit on the central replica, and use that as the guess for every replica
      MultiExpModel fitModel( *this, parInit, ModelCorr, bFreezeCovar, false );
      fitModel.UpdateGuess( parInit, MaxIt, Tolerance );
      DumpParameters( "    Uncorrelated fit on central replica:", parInit );
    }
    // Protected by CancelCritSec
    {
      int Abort{0};
      std::string sError;
      #pragma omp parallel default(shared)
      {
        try
        {
          MultiExpModel fitModel( *this, parInit, ModelCorr, bFreezeCovar, bCorrelated );
          #pragma omp for schedule(dynamic)
          for( int idx = Fold::idxCentral; idx < NSamples; ++idx )
          {
            try
            {
              // Use uncorrelated fit as guess for correlated fit
              if( !Abort )
                fitModel.FitOne( idx, ModelParams, idx == Fold::idxCentral ? SaveCorrMatrixFileName : "",
                                 MaxIt, Tolerance, RelEnergySep, dof, ChiSq );
            }
            catch(const std::exception &e)
            {
              #pragma omp critical(CancelCritSec)
              {
                if( !Abort ) { Abort = 1; sError = e.what(); }
              }
            }
          }
        }
        catch(const std::exception &e)
        {
          #pragma omp critical(CancelCritSec)
          {
            if( !Abort ) { Abort = 1; sError = e.what(); }
          }
        }
      }
      if( Abort )
        throw std::runtime_error( sError );
    }
    ModelParams.MakeCorrSummary( "Params" );
    ModelParams.Write( ModelFileName );
    //ModelParams.WriteSummary( Common::MakeFilename( sModelBase, Common::sModel, Seed, TEXT_EXT ) );
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
      ModelCorr[f].MakeCorrSummary( nullptr );
      ModelCorr[f].WriteSummary( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, TEXT_EXT ));
    }
  }
  // Return the statistics on the fit results
  const int NumSummaries{ ModelParams.NumSamples() }; // because we might have read back an old fit
  std::vector<Common::ValWithEr<scalar>> Results( NumParams );
  std::vector<double> data( NumSummaries );
  for( int p = 0; p < NumParams; p++ ) {
    double Central = ModelParams[Fold::idxCentral][p];
    std::size_t Count{ 0 };
    for( int j = 0; j < NumSummaries; ++j ) {
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
  // Can't do this because of Minuit2     std::ios_base::sync_with_stdio( false );
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
      {"delta", CL::SwitchType::Single, "3"},
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
      {"freeze", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"savecmat", CL::SwitchType::Flag, nullptr},
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
      const bool bFreezeCovar{ cl.GotSwitch( "freeze" ) };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };
      const bool bSaveCMat{ cl.GotSwitch("savecmat") };

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
        Corr.emplace_back();
        Corr[i].Read( sFileName, "  ", &OpNames );
        if( i == 0 )
        {
          outBaseFileName.append( Corr[0].Name_.Base );
          Seed = Corr[0].Name_.Seed;
        }
        else
        {
          Corr[0].IsCompatible( Corr[i], nullptr, true );
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
      Fitter m ( std::move( Corr ), OpNames, bFactor, sOpNameConcat, NumExponents, Verbosity, outBaseFileName, Seed, NSamples );
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
      //if( Common::FileExists( sFitFilename ) )
        //throw std::runtime_error( "Output file " + sFitFilename + " already exists" );
      std::ofstream s;
      for( int tf = tf_start; tf - tf_start < dtf_max; tf++ )
      {
        bool bNeedHeader = true;
        for( int ti = ti_start; ti - ti_start < dti_max; ti++ )
        {
          if( tf - ti + 1 >= delta )
          {
            // Log what file we're processing and when we started
            const std::chrono::system_clock::time_point start{ std::chrono::system_clock::now() };
            std::time_t now = std::chrono::system_clock::to_time_t( start );
            try
            {
              double ChiSq;
              int dof;
              auto params = m.PerformFit(doCorr, bFreezeCovar, ti, tf, Skip, bSaveCorr, bSaveCMat,
                                         MaxIterations, Tolerance, RelEnergySep, ChiSq, dof);
              if( bNeedHeader )
              {
                bNeedHeader = false;
                if( !s.is_open() )
                {
                  s.open( sFitFilename );
                  Common::SummaryHeader<scalar>( s, sSummaryBase );
                  s << "# Seed " << Seed << std::endl;
                }
                else
                {
                  // two blank lines at start of new data block
                  s << "\n" << std::endl;
                }
                // Name the data series
                s << "# [tf=" << tf << "]" << std::endl;
                // Column names, with the series value embedded in the column header (best I can do atm)
                s << "tf=" << tf << Sep << "ti" << Sep << "dof";
                for( int p = 0; p < m.NumParams; p++ )
                  s << Sep << m.ParamNames[p] << Sep << m.ParamNames[p] << "_low" << Sep << m.ParamNames[p] << "_high" << Sep << m.ParamNames[p] << "_check";
                s << " ChiSqPerDof" << std::endl;
              }
              s << tf << Sep << ti << Sep << dof;
              for( int p = 0; p < m.NumParams; p++ )
                s << Sep << params[p];
              s << Sep << ( ChiSq / dof ) << std::endl;
            }
            catch(const std::exception &e)
            {
              std::cout << "Error: " << e.what() << "\n";
            }
            // Mention that we're finished, what the time is and how long it took
            const std::chrono::system_clock::time_point stop{ std::chrono::system_clock::now() };
            now = std::chrono::system_clock::to_time_t( stop );
            const std::chrono::duration<double> duration_seconds = stop - start;
            //const double hours{ ( duration_seconds.count() + 0.5 ) / 3600 };
            std::string sNow{ std::ctime( &now ) };
            while( sNow.length() && sNow[sNow.length() - 1] == '\n' )
              sNow.resize( sNow.length() - 1 );
            std::stringstream ss;
            ss << sNow << ". Total duration " << std::fixed << std::setprecision(1)
                      << duration_seconds.count() << " seconds.\n";
            std::cout << ss.str();
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
    "--delta Minimum number of timeslices in fit range (default 3)\n"
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
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--savecmat Save correlation matrix\n"
    "--help     This message\n";
  }
  return iReturn;
}
