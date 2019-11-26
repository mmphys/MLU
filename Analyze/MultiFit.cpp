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

#include "CommonLatAn.hpp"
#include "MultiFit.hpp"

// Uncomment the next line if your cmath doesn't define M_PI etc by default
//#define _USE_MATH_DEFINES
#include <cmath>
//#include <iomanip>
#include <ios>
#include <iostream>
//#include <sys/stat.h>
#include <LatAnalyze/Core/OptParser.hpp>
//#include <LatAnalyze/Eigen/Dense>
//#include <LatAnalyze/Functional/CompiledModel.hpp>
#include <LatAnalyze/Io/Io.hpp>
//#include <LatAnalyze/Statistics/MatSample.hpp>
//#include <LatAnalyze/Core/Math.hpp>
#include <LatAnalyze/Numerical/MinuitMinimizer.hpp>
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

// Minimisation for multi-exponential fit

class MultiExpModel : public ROOT::Minuit2::FCNBase {
public:
  static constexpr double pi{ M_PI };
  using vd_t = std::vector<double>;
  explicit MultiExpModel( const Latan::DMatSample Corr, int NumOps, const std::vector<std::string> &OpNames, int NumExponents, int Verbosity, const std::string &OutputBaseName, Common::SeedType Seed );
  virtual ~MultiExpModel() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const vd_t & ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  // These are definitely not part of the FCNBase interface
  std::vector<Common::ValWithEr> PerformFit( bool bCorrelated, int tMin, int tMax, double &ChiSq, int &dof, bool bSaveCorr );
  const int NumOps;
  const std::vector<std::string> &OpNames;
  const int Verbosity;
  const std::string OutputBaseName;
  const Common::SeedType Seed;
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int NumExponents;
  const int NumParams;
  const Latan::DMatSample Corr;
  const std::vector<std::string> ParamNames;
  private:
  // These variables only used during a fit
  int tMin = -1;
  int tMax = -1;
  int NtCorr;
  int Extent;
  Latan::DMat Covar;
  Latan::DMat CovarInv;
  std::vector<double> VarianceInv;
  bool bCorrelated;
  std::string OutputRunBase;
  Latan::Index idx;
  // Helper functions
  inline int AlphaIndex( int snk, int src) const { return snk * NumOps + src; };
  inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; };
  void MakeCovar( void );
  void DumpParameters( const std::string &msg, const ROOT::Minuit2::MnUserParameters &par );
  std::vector<std::string> MakeParamNames();
};

MultiExpModel::MultiExpModel( const Latan::DMatSample corr_, int numOps_, const std::vector<std::string> &opNames_, int numExponents_, int verbosity_, const std::string &outputBaseName_, Common::SeedType seed_ )
  : NumOps{ numOps_ },
    OpNames{ opNames_ },
    Verbosity{ verbosity_ },
    OutputBaseName{ outputBaseName_ },
    Seed{ seed_ },
    NumFiles{ numOps_ * numOps_ },
    NSamples{ static_cast<int>( corr_.size() ) },
    Nt      { static_cast<int>( corr_[Latan::central].rows() / NumFiles ) },
    NumExponents{ numExponents_ },
    NumParams{ NumExponents * ( 1 + NumOps ) },
    Corr{ corr_ },
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
  assert( Corr[Latan::central].rows() == NumFiles * Nt && "Bug" );
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
    const double ExpectedX{ Corr[Latan::central]( ReadX,  0 ) };
    const int yMin{ bCorrelated ? 0 : x };
    for( int y = yMin; y <= x; ++y ) {
      const int f2{ y / NtCorr };
      const int t2{ y % NtCorr };
      const int ReadY{ Nt * f2 + t2 + tMin };
      const double ExpectedY{ Corr[Latan::central]( ReadY, 0 ) };
      double z = 0;
      for( int i = 0; i < NSamples; ++i )
        z += ( Corr[i]( ReadX, 0 ) - ExpectedX ) * ( Corr[i]( ReadY, 0 ) - ExpectedY );
      z /= NSamples;
      if( bCorrelated ) {
        Covar(x, y) = z;
        if( x != y )
          Covar(y, x) = z;
      }
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
      Latan::Io::save(Covar, OutputRunBase + "Covar.h5");
      Latan::Io::save(CovarInv, OutputRunBase + "CovarInv.h5");
      Latan::DMat CovarProd{ Covar * CovarInv };
      Latan::Io::save(CovarProd, OutputRunBase + "CovarProd.h5");
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

double MultiExpModel::operator()( const vd_t & par ) const {
  const double * const Energy{ par.data() };
  const double * const Coeff{ Energy + NumExponents };
  // Calculate the theory errors for these model parameters
  double ModelError[Extent]; // Should happily fit on stack
  for( int snk = 0; snk < NumOps; ++snk )
    for( int src = 0; src < NumOps; ++src )
      for( int t = 0; t < NtCorr; ++t ) {
        double z = 0;
        double CumulativeEnergy = 0;
        for( int e = 0; e < NumExponents; ++e ) {
          CumulativeEnergy -= Energy[e];
          z += Coeff[MELIndex( src, e )] * Coeff[MELIndex( snk, e )] * std::exp( CumulativeEnergy * ( t + tMin ) );
        }
        const int iRead{ AlphaIndex( snk, src ) * Nt + t + tMin };
        const int iWrite{ AlphaIndex( snk, src ) * NtCorr + t };
        ModelError[iWrite] = Corr[idx]( iRead, 0 ) - z;
      }
  static int iCall{ 0 };
  iCall++;
  if( !Common::IsFinite( ModelError, Extent ) ) {
    if( Verbosity ) {
      static int iOOB{ 0 };
      std::cout << "Call " << iCall << ", " << ++iOOB << "th overflow\n";
    }
    return std::numeric_limits<double>::max();
  }

  // The inverse of a symmetric matrix is also symmetric, so only calculate half the matrix
  double chi2 = 0;
  for( int i = 0; i < Extent; ++i ) {
    if( !bCorrelated )
      chi2 += ModelError[i] * VarianceInv[i] * ModelError[i];
    else {
      for( int j = 0; j <= i; ++j ) {
        double z = ModelError[i] * CovarInv(i, j) * ModelError[j];
        chi2 += z;
        if( i != j )
          chi2 += z;
      }
    }
  }
  if( Verbosity > 2 )
    std::cout << "Call " << iCall << ", chi^2=" << chi2 << ", E0=" << Energy[0] << ", E1=" << Energy[1] << "\n";
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
std::vector<Common::ValWithEr> MultiExpModel::PerformFit( bool Bcorrelated_, int TMin_, int TMax_, double &ChiSq, int &dof, bool bSaveCorr )
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
    VarianceInv.empty();
    CovarInv.resize( 0, 0 );
    Covar.resize( 0, 0 );
  }
  // Make the covariance matrix / vector if need be
  if( ( bCorrelated && Covar.rows() != Extent ) || ( !bCorrelated && VarianceInv.size() != Extent ) )
    MakeCovar();

  // Results for each bootstrap sample
  Latan::DMatSample ModelParams( NSamples, NumParams, 2 );      // model parameters (from the fit), NB Absolute energies
  std::vector<Latan::DMatSample> ModelCorr( NumOps * NumOps );  // correlators resulting from the fit parameters
  for( int i = 0; i < NumOps * NumOps; i++ ) {
    ModelCorr[i].resize( NSamples );
    ModelCorr[i].resizeMat( Nt, 2 ); // Complex component will be zero
  }

  // Make initial guesses for the parameters
  // For each Exponent, I need the delta_E + a constant for each operator
  std::vector<double> par( NumParams );
  {
    // Take a starting guess for the parameters - same as LatAnalyze
    double * const Energy{ par.data() };
    double * const Coeff{ Energy + NumExponents };
    const Latan::DMat & m{ Corr[Latan::central] };
    const int tGuess{Nt / 4};
    Energy[0] = 0;
    for( int snk = 0; snk < NumOps; ++snk )
      for( int src = 0; src < NumOps; ++src ) {
        const int i{ AlphaIndex( snk, src ) * Nt + tGuess };
        double E0 = std::log( m( i, 0 ) / m( i + 1, 0 ) );
        Energy[0] += E0;
        if( snk == src )
          Coeff[snk] = std::sqrt( std::abs( m( i, 0 ) / ( std::exp( - E0 * tGuess ) ) ) );
      }
    Energy[0] /= NumFiles;
    // Now guess Higher exponents - same as LatAnalyze
    const double MELFactor{ std::sqrt( 0.5 ) };
    for( int e = 1; e < NumExponents; ++e ) {
      Energy[e] = Energy[e - 1] * ( 1 << ( e - 1 ) ); // Cumulative
      for( int o = 0; o < NumOps; ++o )
        Coeff[ MELIndex( o, e ) ] = Coeff[ MELIndex( o, e - 1 ) ] * MELFactor;
    }
  }

  // Create minimisation parameters
  std::vector<double> Error( NumParams, 0.1 );
  ROOT::Minuit2::MnUserParameters upar( par, Error );
#ifdef _DEBUG
  //upar.SetLowerLimit( 0, 0.1 );
  //upar.SetValue(0,0.184);
  //upar.Fix(0);
#endif
  if( Verbosity )
    DumpParameters( "Initial guess:", upar );

  //ROOT::Minuit2::Minuit2Minimizer minimizer( ROOT::Minuit2::kMigrad ); // even though this is default
  //minimizer.SetPrintLevel( Verbosity );
  //minimizer.SetFunction();
  //VariableIterator i;
  //ROOT::Minuit2::VariableMetricMinimizer minimizer;
  ROOT::Minuit2::MnMigrad minimizer( *this, upar );
  //FOR_STAT_ARRAY( Fit, <#i#> ) {
  //for (idx = -ModelParams.offset; idx < ModelParams.size(); idx++) {
  for (int LoopIdx = -1; LoopIdx < NSamples; LoopIdx++) {
    idx = ( LoopIdx >= 0 && LoopIdx < NSamples ) ? LoopIdx : Latan::central;
    ROOT::Minuit2::FunctionMinimum min = minimizer();
    ROOT::Minuit2::MnUserParameterState state = min.UserState();
    //assert( minimizer.SetVariables( par.begin(), par.end() ) == NumParams && "Error: could not set variables" );
    //minimizer.SetFunction( *this );
    //ROOT::Minuit2::FunctionMinimum min = minimizer.Minimize( *this, par, Error );
    //ROOT::Minuit2::MnUserParameterState state = min.UserState();
    //bool bOK{ minimizer.Minimize() };
    if( !state.IsValid() )
      throw std::runtime_error( "ERROR: Fit " + std::to_string( idx ) + " did not converge on an answer\n" );
    // Throw away the first fit result - just use it as a seed for the next fit
    if( LoopIdx >= -1 )
    {
      //ROOT::Minuit2::MnUserParameters params = min.UserParameters();
      //std::cout << "Fit result:" << min << std::endl;
      upar = state.Parameters();
      if( idx == Latan::central || Verbosity > 2 ) {
        const std::string FitResult{ "Fit result:" };
        if( Verbosity )
          std::cout << FitResult << state << "\n";
        DumpParameters( FitResult, upar );
      }
      if( idx == Latan::central ) {
        ChiSq = state.Fval();
        dof = Extent;
        if( bCorrelated )
          dof *= dof;
        dof -= NumParams;
        std::cout << "Chi^2=" << ChiSq << ", dof=" << dof << ", chi^2/dof=" << ChiSq / dof
                  << "\n\t computing statistics\n";
      }
      // Save the fit parameters for this replica
      double Cumulative;
      Latan::DMat & FitParams{ ModelParams[idx] };
      for( int j = 0; j < NumParams; ++j ) {
        double z = upar.Value( j ); // Need this so Energy and Coeff are up-to-date
        if( j > 0 && j < NumExponents )
          Cumulative += z;
        else
          Cumulative  = z;
        FitParams(j, 0) = Cumulative;
        FitParams(j, 1) = 0;
      }
      // Save the reconstructed correlator values for this replica
      for( int snk = 0; snk < NumOps; ++snk )
        for( int src = 0; src < NumOps; ++src ) {
          Latan::DMat &mc{ ModelCorr[AlphaIndex( snk, src )][idx] };
          for( int t = 0; t < Nt; ++t ) {
            double z = 0;
            for( int e = 0; e < NumExponents; ++e )
              z +=   FitParams( MELIndex( src, e ) + NumExponents, 0)
                   * FitParams( MELIndex( snk, e ) + NumExponents, 0)
                   * std::exp( - FitParams( e, 0 ) * t );
            mc(t,0) = z;
            mc(t,1) = 0;
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
  Latan::Io::save(ModelParams, Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ));
  //Common::SummariseBootstrap(ModelParams, sModelBase, Seed, "params" );
  for( int snk = 0; snk < NumOps; ++snk )
    for( int src = 0; src < NumOps; ++src ) {
      const int i{ AlphaIndex( snk, src ) };
      const std::string SummaryBase{ OutputRunBase + '.' + OpNames[snk] + '_' + OpNames[src] };
      if( bSaveCorr )
        Latan::Io::save(ModelCorr[i], Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, DEF_FMT ) );
      Common::SummariseBootstrapCorr(ModelCorr[i], SummaryBase, Seed );
    }
  /*if( NumOps == 2) {  // I've only coded for two operators
    Latan::DMatSample OptimalCorr( NSamples, Nt, 2 );
    for( int degrees = -90; degrees <= 90; degrees ++ )
    {
      //if( degrees % 3 && degrees % 5 ) continue;
      const double theta{ pi * degrees / 180 };
      const double costheta{ cos( theta ) };
      const double sintheta{ sin( theta ) };
      const double cos_sq_theta{ costheta * costheta };
      const double sin_sq_theta{ sintheta * sintheta };
      const double cos_sin_theta{ costheta * sintheta };
      for (Latan::Index i = Latan::central; i < NSamples; i++) {
        const double A_P1{ ModelParams[i](3, 0) };
        const double A_A1{ ModelParams[i](5, 0) };
        const double Op_PP{ cos_sq_theta / ( A_P1 * A_P1 ) };
        const double Op_AP{ cos_sin_theta / ( A_P1 * A_A1 ) };
        const double Op_AA{ sin_sq_theta / ( A_A1 * A_A1 ) };
        for( int t = 0; t < Nt; t++ ) {
          OptimalCorr[i](t,0) = Op_PP * ModelCorr[0][i](t,0) + Op_AA * ModelCorr[3][i](t,0)
                              + Op_AP * ( ModelCorr[1][i](t,0) + ModelCorr[2][i](t,0) );
          OptimalCorr[i](t,1) = 0;
        }
      }
      std::string sOptimusPrime{ OutputRunBase };
      sOptimusPrime.append( ".theta_" );
      sOptimusPrime.append( std::to_string(degrees) );
      if( bSaveCorr )
        Latan::Io::save( OptimalCorr, Common::MakeFilename( sOptimusPrime, Common::sBootstrap, Seed, DEF_FMT ) );
      Common::SummariseBootstrapCorr( OptimalCorr, sOptimusPrime, Seed );
    }
  }*/
  // Return the statistics on the fit results
  std::vector<Common::ValWithEr> Results( NumParams );
  std::vector<double> data( NSamples );
  for( int p = 0; p < NumParams; p++ ) {
    double Central = ModelParams[Latan::central](p,0);
    std::size_t Count{ 0 };
    for( int j = 0; j < NSamples; ++j ) {
      double d = ModelParams[j](p, 0);
      if( std::isfinite( d ) )
        data[Count++] = d;
    }
    Results[p].Get( Central, data, Count );
  }
  return Results;
}

int main(int argc, char *argv[])
{
  using namespace Latan;
  std::ios_base::sync_with_stdio( false );
  // parse arguments /////////////////////////////////////////////////////////
  OptParser opt;
  opt.addOption("" , "ti"       , OptParser::OptType::value  , false,
                "initial fit time");
  opt.addOption("" , "tf"       , OptParser::OptType::value  , false,
                "final fit time");
  opt.addOption("" , "dti"       , OptParser::OptType::value , true,
                "number of initial fits", "1");
  opt.addOption("" , "dtf"       , OptParser::OptType::value , true,
                "number of final fits", "1");
  opt.addOption("t" , "thinning", OptParser::OptType::value  , true,
                "thinning of the time interval", "1");
  opt.addOption("s", "shift"    , OptParser::OptType::value  , true,
                "time variable shift", "0");
  opt.addOption("m", "model"    , OptParser::OptType::value  , true,
                "fit model (exp|exp2|exp3|cosh|cosh2|cosh3|explin|<interpreter code>)", "cosh");
  opt.addOption("" , "nPar"     , OptParser::OptType::value  , true,
                "number of model parameters for custom models "
                "(-1 if irrelevant)", "-1");
  opt.addOption("" , "svd"      , OptParser::OptType::value  , true,
                "singular value elimination threshold", "0.");
  opt.addOption("v", "verbosity", OptParser::OptType::value  , true,
                "minimizer verbosity level (0|1|2)", "0");
  opt.addOption("o", "output",    OptParser::OptType::value  , true,
                "output base filename", "");
  opt.addOption("" , "uncorr"   , OptParser::OptType::trigger, true,
                "only do the uncorrelated fit");
  opt.addOption("f", "fold"     , OptParser::OptType::value, true,
                "fold the correlator (0=don't fold, 1=+parity, -1=-parity)", "0");
  opt.addOption("p", "plot"     , OptParser::OptType::trigger, true,
                "show the fit plot");
  opt.addOption("r", "range"    , OptParser::OptType::value  , true,
                "vertical range multiplier in plots", "20.");
  opt.addOption("h", "heatmap"  , OptParser::OptType::trigger, true,
                "show the fit correlation heatmap");
  opt.addOption("", "help"      , OptParser::OptType::trigger, true,
                "show this help message and exit");
  opt.addOption("e", "exponents", OptParser::OptType::value  , true,
                "number of exponents", "0");
  opt.addOption("" , "uncorr"   , OptParser::OptType::trigger, true,
                "perform uncorrelated fit (default=correlated");
  opt.addOption("" , "opnames"   , OptParser::OptType::trigger, true,
                "operator names (default=op_n");
  opt.addOption("" , "nosave"   , OptParser::OptType::trigger, true,
                "don't save bootstrap replicas of correlators or model");
  opt.addOption("" , "savemodel"   , OptParser::OptType::trigger, true,
                "save bootstrap of model fit results");
  if (!opt.parse(argc, argv) or (opt.getArgs().size() < 1) or opt.gotOption("help"))
  {
    std::cerr << "usage: " << argv[0] << " <options> <correlator file>" << std::endl;
    std::cerr << std::endl << "Possible options:" << std::endl << opt << std::endl;
    
    return EXIT_FAILURE;
  }
  const int ti{ opt.optionValue<int>("ti") };
  const int tf{ opt.optionValue<int>("tf") };
  const int dti_max{ opt.optionValue<int>("dti") };
  const int dtf_max{ opt.optionValue<int>("dtf") };
  //int thinning            = opt.optionValue<int>("t");
  const int shift{ opt.optionValue<int>("s") };
  //Index nPar              = opt.optionValue<Index>("nPar");
  //double svdTol           = opt.optionValue<double>("svd");
  //bool doCorr             = !opt.gotOption("uncorr");
  const int fold{ opt.optionValue<int>("fold") };
  //bool doPlot             = opt.gotOption("p");
  //double plotrange        = opt.optionValue<Index>("range");
  //bool doHeatmap          = opt.gotOption("h");
  const std::string model{ opt.optionValue("m") };
  std::string outBaseFileName{ opt.optionValue("o") };
  const bool doCorr{ !opt.gotOption("uncorr") };
  const int Verbosity{ opt.optionValue<int>("v") }; // 0 = normal, 1=debug, 2=save covar, 3=Every iteration
  const bool bSaveCorr{ opt.gotOption("savecorr") };
  /*Minimizer::Verbosity verbosity;
  switch ( Verbosity )
  {
    case 0:
      verbosity = Minimizer::Verbosity::Silent;
      break;
    case 1:
      verbosity = Minimizer::Verbosity::Normal;
      break;
    case 2:
      verbosity = Minimizer::Verbosity::Debug;
      break;
    default:
      std::cerr << "error: wrong verbosity level" << std::endl;
      return EXIT_FAILURE;
  }*/
  
  // load correlators /////////////////////////////////////////////////////////
  const int NumFiles{ static_cast<int>( opt.getArgs().size() ) };
  int NumOps = 1; // Number of operators. Should equal square root of number of files
  while( NumOps * NumOps < NumFiles )
    NumOps++;
  if( NumOps * NumOps != NumFiles ) {
    std::cerr << "Number of files should be a perfect square" << std::endl;
    return EXIT_FAILURE;
  }
  std::vector<std::string> OpNames;
  Common::SeedType Seed = 0;
  try{
    {
      std::vector<Common::FileNameAtt> FileNames;
      std::size_t i = 0;
      for( const std::string &sFileName : opt.getArgs() ) {
        FileNames.emplace_back( sFileName, OpNames );
        if( i == 0 ) {
          outBaseFileName.append( FileNames[0].Base );
          Seed = FileNames[0].Seed;
        } else {
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
        for( int src = 0; src < NumOps; ++src ) {
          Common::FileNameAtt &f{ FileNames[snk * NumOps + src] };
          if( f.op[0] != src || f.op[1] != snk )
            throw std::runtime_error( "Warning: Operator order should be sink-major, source minor" );
        }
    }
    int NumExponents{ opt.optionValue<int>( "exponents" ) };
    if( NumExponents == 0 )
      NumExponents = NumOps;
    MultiExpModel m ( Common::ReadBootstrapCorrs( opt.getArgs(), fold, shift, NumOps ), NumOps, OpNames, NumExponents, Verbosity, outBaseFileName, Seed );
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
        double ChiSq;
        int dof;
        auto params = m.PerformFit( doCorr, ti + dti, tf + dtf, ChiSq, dof, bSaveCorr );
        s << ( tf + dtf ) << Sep << ( ti + dti );
        for( int p = 0; p < m.NumParams; p++ )
          s << Sep << params[p];
        s << Sep << ChiSq << Sep << dof << Sep << ( ChiSq / dof ) << std::endl;
      }
    }
  }
  catch(const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
