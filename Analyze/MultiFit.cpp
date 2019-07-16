#include "Common.hpp"
#include "MultiFit.hpp"

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
  using vd_t = std::vector<double>;
  explicit MultiExpModel( const Latan::DMatSample Corr, int tMin, int tMax, int NumOps, int NumExponents, int Verbosity );
  virtual ~MultiExpModel() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const vd_t & ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  // These are definitely not part of the FCNBase interface
  void PerformFit( bool bCorrelated );
private:
  const int Verbosity;
  const int NumExponents;
  const int NumOps;
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int tMin;
  const int tMax;
  const int NtCorr;
  const Latan::DMatSample Corr;
  bool bGotCovar = false;
  Latan::DMat Covar;
  Latan::DMat CovarInv;
  // Helper functions
  void MakeCovar( bool bCorrelated );
  inline int AlphaIndex( int src, int snk) const { return snk * NumOps + src; };
  inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; };
};

MultiExpModel::MultiExpModel( const Latan::DMatSample corr_, int TMin_, int TMax_, int numOps_, int numExponents_, int verbosity_ )
  : Verbosity{ verbosity_ },
    NumExponents{ numExponents_ },
    NumOps{ numOps_ },
    NumFiles{ numOps_ * numOps_ },
    NSamples{ static_cast<int>( corr_.size() ) },
    Nt      { static_cast<int>( corr_[Latan::central].rows() / NumFiles ) },
    tMin{ TMin_ },
    tMax{ TMax_ },
    NtCorr{ TMax_ - TMin_ + 1 },
    Corr{ corr_ },
    Covar{ NumFiles * NtCorr, NumFiles * NtCorr },
    CovarInv{ NumFiles * NtCorr, NumFiles * NtCorr }
{
  assert( numExponents_ > 0 && "Number of exponents in the fit should be positive" );
}

// Make the covariance matrix, for only the timeslices we are interested in
void MultiExpModel::MakeCovar( bool bCorrelated )
{
  if( bGotCovar )
    return;
  if( bCorrelated )
    std::cout << "Creating covariance matrix ...\n";
  else
    std::cout << "Creating variance vector ...\n";
  assert( Corr[Latan::central].rows() == NumOps * NumOps * Nt && "Bug" );
  assert( tMin >= 0 && tMin < Nt && "Error: tMin invalid" );
  assert( tMax >= 0 && tMax < Nt && "Error: tMax invalid" );
  assert( NtCorr > 1 && "Error: tMin >= tMax" );
  const int Extent{ NumFiles * NtCorr };
  for( int x = 0; x < Extent; ++x ) {
    const int f1{ x / NtCorr };
    const int t1{ x % NtCorr };
    const int ReadX{ Nt * f1 + t1 + tMin };
    const double ExpectedX{ Corr[Latan::central]( ReadX,  0 ) };
    for( int y = 0; y <= x; ++y ) {
      if( !bCorrelated && x != y ) {
        Covar(x, y) = 0;
        Covar(y, x) = 0;
      } else {
        const int f2{ y / NtCorr };
        const int t2{ y % NtCorr };
        const int ReadY{ Nt * f2 + t2 + tMin };
        const double ExpectedY{ Corr[Latan::central]( ReadY, 0 ) };
        double z = 0;
        for( int i = 0; i < NSamples; ++i )
          z += ( Corr[i]( ReadX, 0 ) - ExpectedX ) * ( Corr[i]( ReadY, 0 ) - ExpectedY );
        z /= NSamples;
        Covar(x, y) = z;
        if( x != y )
          Covar(y, x) = z;
      }
    }
  }
  assert( Common::IsFinite( Covar, !bCorrelated ) && "Error: covariance matrix isn't finite" );
  if( !bCorrelated )
    for( int i = 0; i < Extent; i++ )
      for( int j = 0; j < Extent; j++ )
        if( i == j )
          CovarInv( i, i ) = 1. / Covar( i, i );
        else
          CovarInv( i, j ) = 0;
  else {
    CovarInv = Covar.inverse();
    bGotCovar = true;
  }
  assert( Common::IsFinite( CovarInv, !bCorrelated ) && "Error: inverse covariance matrix isn't finite" );
#ifdef DEBUG
  if( bCorrelated ) {
    std::cout << "Covariance matrix is " << Covar.rows() << " x " << Covar.cols() << "\n";
    std::cout << "Inverse covariance matrix is " << CovarInv.rows() << " x " << CovarInv.cols() << "\n";
    if( Verbosity ) {
      std::cout << "Correlation matrix:\n";
      for( int i = 0; i < ( NumFiles > 1 ? 1 : 0 ) * NtCorr + DEBUG_NUM; ++i ) {
        std::cout << " " << i << ":";
        for( int j = 0; j < DEBUG_NUM; ++j )
          std::cout << " " << Covar(i, j) / ( sqrt( Covar(i, i) ) * sqrt( Covar(j, j) ) );
        std::cout << "\n";
      }
    }
    Latan::Io::save(Covar, "Debug.Covar.h5");
    Latan::Io::save(CovarInv, "Debug.CovarInv.h5");
    Latan::DMat CovarProd{ Covar * CovarInv };
    Latan::Io::save(CovarProd, "Debug.CovarProd.h5");
    const double CondNumber{ Covar.norm() * CovarInv.norm() };
    const int CondDigits{ static_cast<int>( std::log10( CondNumber ) + 0.5 ) };
    if( CondDigits >= 12 )
      std::cout << "ERROR: "; // see https://en.wikipedia.org/wiki/Condition_number
    std::cout << "Covariance matrix condition number=" << CondNumber
      << ", i.e. potential loss of " << CondDigits << " digits.\n";
  }
#endif
}

// Compute chi-squared given parameters for multi-exponential fit
  // Energy Levels        Parameters 0 ... NumExponents - 1
  //  NB: second and subsequent energy levels are deltas
  // Overlap Coefficients Parameters NumExponents ... NumExponents ( NumOps + 1 ) - 1

double MultiExpModel::operator()( const vd_t & par ) const {
  const int NeedSize{ NumExponents * ( 1 + NumOps ) };
  assert( par.size() == NeedSize && "Error: number of parameters" );
  const double * const Energy{ par.data() };
  const double * const Coeff{ Energy + NumExponents };
  // Calculate the theory errors for these model parameters
  const int Extent{ NumFiles * NtCorr };
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
        const int iRead{ AlphaIndex( src, snk ) * Nt + t + tMin };
        const int iWrite{ AlphaIndex( src, snk ) * NtCorr + t };
        ModelError[iWrite] = Corr[Latan::central]( iRead, 0 ) - z;
      }
  //assert( Common::IsFinite( ModelError, Extent ) && "Error: Model errors aren't finite" );
  static int iCall{ 0 };
  iCall++;
  if( !Common::IsFinite( ModelError, Extent ) ) {
    static int iOOB{ 0 };
    std::cout << "Call " << iCall << ", " << ++iOOB << "th overflow\n";
    return std::numeric_limits<double>::max();
  }

  // The inverse of a symmetric matrix is also symmetric, so only calculate half the matrix
  double chi2 = 0;
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j < Extent; ++j ) {
      double z = ModelError[i] * CovarInv(i, j) * ModelError[j];
      chi2 += z;
      //if( i != j )
        //chi2 += z;
    }
  if( Verbosity )
    std::cout << "Call " << iCall << ", chi^2=" << chi2 << ", E0=" << Energy[0] << ", E1=" << Energy[1] << "\n";
  return chi2;
}

// Perform a fit
  void MultiExpModel::PerformFit( bool bCorrelated )
{
  std::cout << "=========================================\nPerforming ";
  if( bCorrelated )
    std::cout << "correlated";
  else
    std::cout << "uncorrelated";
  std::cout << " fit\n";
  MakeCovar( bCorrelated );
  //ROOT::Minuit2::Minuit2Minimizer minimizer( ROOT::Minuit2::kMigrad ); // even though this is default
  //minimizer.SetPrintLevel( Verbosity );
  //minimizer.SetFunction();
  //VariableIterator i;
  const int NeedSize{ NumExponents * ( 1 + NumOps ) };
  std::vector<double> par( NeedSize );
  double * const Energy{ par.data() };
  double * const Coeff{ Energy + NumExponents };
  // Take a starting guess for the parameters - same as LatAnalyze
  const Latan::DMat & m{ Corr[Latan::central] };
  const int tGuess{Nt / 4};
  Energy[0] = 0;
  for( int snk = 0; snk < NumOps; ++snk )
    for( int src = 0; src < NumOps; ++src ) {
      const int i{ AlphaIndex( src, snk ) * Nt + tGuess };
      double E0 = std::log( m( i, 0 ) / m( i + 1, 0 ) );
      Energy[0] += E0;
      if( snk == src )
        Coeff[snk] = std::sqrt( std::abs( m( i, 0 ) / ( std::exp( - E0 * tGuess ) ) ) );
    }
  Energy[0] /= NumFiles;
  // Now guess Higher exponents - same as LatAnalyze
  const double MELFactor{ std::sqrt( 0.5 ) };
  for( int e = 1; e < NumExponents; ++e ) {
    //Energy[e] = 2 * Energy[e - 1]; // Absolute
    Energy[e] = Energy[e - 1] * ( 1 << ( e - 1 ) ); // Cumulative
    for( int o = 0; o < NumOps; ++o )
      Coeff[ MELIndex( o, e ) ] = Coeff[ MELIndex( o, e - 1 ) ] * MELFactor;
  }
  std::vector<double> Error( NeedSize, 0. );
  std::cout << "Initial guess:\n";
  for( int i = 0; i < par.size(); ++i )
    std::cout << "\tParameter[" << i << "]=" << par[i] << ", +/- " << Error[i] << "\n";

  // Create minimisation parameters
  ROOT::Minuit2::MnUserParameters upar( par, Error );
  //upar.SetPrecision(0.001);
  //upar.SetLowerLimit( 0, 0.1 );
  //for( int e = 1; e < NumExponents; ++e )
    //upar.SetLowerLimit( e, 0. );
#ifdef _DEBUG
  //upar.SetValue(0,0.184);
  upar.Fix(0);
#endif

  //ROOT::Minuit2::VariableMetricMinimizer minimizer;
  ROOT::Minuit2::MnMigrad minimizer( *this, upar );
  ROOT::Minuit2::FunctionMinimum min = minimizer();
  ROOT::Minuit2::MnUserParameterState state = min.UserState();
  //assert( minimizer.SetVariables( par.begin(), par.end() ) == NeedSize && "Error: could not set variables" );
  //minimizer.SetFunction( *this );
  //ROOT::Minuit2::FunctionMinimum min = minimizer.Minimize( *this, par, Error );
  //ROOT::Minuit2::MnUserParameterState state = min.UserState();
  bool bOK{ state.IsValid() };
  //bool bOK{ minimizer.Minimize() };
  std::cout << "state.IsValid() = " << bOK << "\n";
  ROOT::Minuit2::MnUserParameters params = min.UserParameters();
  std::cout << "Fit result:" << state << std::endl;
  //std::cout << "Fit result:" << min << std::endl;
}

/*static const char DefaultOutputStem[] = "@corr@.@type@";

static void Usage( const char * pExecutableName )
{
  std::string cmd{ pExecutableName };
  auto pos = cmd.find_last_of('/');
  if( pos != std::string::npos && pos < cmd.length() - 1 )
    cmd = cmd.substr( pos + 1 );
  std::cerr << "usage: " << cmd << " <options> ContractionFile1 [ContractionFile2 ...]"
  << "\nOutput stem (if present) must contain '@corr@' and '@type', e.g. '" << DefaultOutputStem << "'\n";
}*/

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
                "output file", "");
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
  if (!opt.parse(argc, argv) or (opt.getArgs().size() < 1) or opt.gotOption("help"))
  {
    std::cerr << "usage: " << argv[0] << " <options> <correlator file>" << std::endl;
    std::cerr << std::endl << "Possible options:" << std::endl << opt << std::endl;
    
    return EXIT_FAILURE;
  }
  const int ti{ opt.optionValue<int>("ti") };
  const int tf{ opt.optionValue<int>("tf") };
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
  const std::string outFileName{ opt.optionValue("o") };
  const bool doCorr{ !opt.gotOption("uncorr") };
  const int Verbosity{ opt.optionValue<int>("v") };
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
  
  // load correlator /////////////////////////////////////////////////////////
  const std::vector<std::string> & FileName( opt.getArgs() );

  const int NumFiles{ static_cast<int>( opt.getArgs().size() ) };
  int NumOps = 1; // Number of operators. Should equal square root of number of files
  while( NumOps * NumOps < NumFiles )
    NumOps++;
  if( NumOps * NumOps != NumFiles ) {
    std::cerr << "Number of files should be a perfect square" << std::endl;
    return EXIT_FAILURE;
  }
  int NumExponents{ opt.optionValue<int>( "exponents" ) };
  if( NumExponents == 0 )
    NumExponents = NumOps;
  MultiExpModel m ( Common::ReadBootstrapCorrs( FileName, fold, shift, NumOps ), ti, tf, NumOps, NumExponents, Verbosity );
  m.PerformFit( false );
  if( doCorr )
    m.PerformFit( true );

  return EXIT_SUCCESS;
}