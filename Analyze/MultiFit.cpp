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
//#include <LatAnalyze/Io/Io.hpp>
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
  explicit MultiExpModel( const Latan::DMatSample Corr, int tMin, int tMax, int NumOps, int NumExponents = 2 );
  virtual ~MultiExpModel() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const vd_t & ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  // These are definitely not part of the FCNBase interface
  void PerformFit( void );
private:
  const int NumExponents;
  const int NumOps;
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int tMin;
  const int tMax;
  const int NtCorr;
  const Latan::DMatSample Corr;
  const Latan::DMat Covar;
  const Latan::DMat CovarInv;
  // Helper functions
  Latan::DMat MakeCovar( void );
  inline int AlphaIndex( int src, int snk) const { return src * NumOps + snk; };
  inline int MELIndex( int op, int EnergyLevel) const { return op * NumExponents + EnergyLevel; };
};

MultiExpModel::MultiExpModel( const Latan::DMatSample corr_, int TMin_, int TMax_, int numOps_, int numExponents_ )
  : NumExponents{ numExponents_ },
    NumOps{ numOps_ },
    NumFiles{ numOps_ * numOps_ },
    NSamples{ static_cast<int>( corr_.size() ) },
    Nt      { static_cast<int>( corr_[Latan::central].rows() / NumFiles ) },
    tMin{ TMin_ },
    tMax{ TMax_ },
    NtCorr{ TMax_ - TMin_ + 1 },
    Corr{ corr_ },
    Covar{ MakeCovar() },
    CovarInv{ Covar.inverse() }
{
  assert( numExponents_ > 0 && "Number of exponents in the fit should be positive" );
  std::cout << "Covariance matrix is " << Covar.rows() << " x " << Covar.cols() << "\n";
  std::cout << "Inverse covariance matrix is " << CovarInv.rows() << " x " << CovarInv.cols() << "\n";
  std::cout << "Correlation matrix:\n";
  for( int i = 0; i < ( NumFiles > 1 ? 1 : 0 ) * NtCorr + DEBUG_NUM; ++i ) {
    std::cout << " " << i << ":";
    for( int j = 0; j < DEBUG_NUM; ++j ) {
      std::cout << " " << Covar(i, j) / ( sqrt( Covar(i, i) ) * sqrt( Covar(j, j) ) );
    }
    std::cout << "\n";
  }
  assert( Common::IsFinite( CovarInv ) && "Error: inverse covariance matrix isn't finite" );
}

// Make the covariance matrix, for only the timeslices we are interested in
Latan::DMat MultiExpModel::MakeCovar( void )
{
  std::cout << "Creating covariance matrix ... ";
  assert( Corr[Latan::central].rows() == NumOps * NumOps * Nt && "Bug" );
  assert( tMin >= 0 && tMin < Nt && "Error: tMin invalid" );
  assert( tMax >= 0 && tMax < Nt && "Error: tMax invalid" );
  assert( NtCorr > 1 && "Error: tMin >= tMax" );
  const int Extent{ NumFiles * NtCorr };
  Latan::DMat covar( Extent, Extent );
  for( int x = 0; x < Extent; ++x ) {
    const int f1{ x / NtCorr };
    const int t1{ x % NtCorr };
    const int ReadX{ Nt * f1 + t1 + tMin };
    const double ExpectedX{ Corr[Latan::central]( ReadX,  0 ) };
    for( int y = 0; y <= x; ++y ) {
      const int f2{ y / NtCorr };
      const int t2{ y % NtCorr };
      const int ReadY{ Nt * f2 + t2 + tMin };
      const double ExpectedY{ Corr[Latan::central]( ReadY, 0 ) };
      double z = 0;
      for( int i = 0; i < NSamples; ++i )
        z += ( Corr[i]( ReadX, 0 ) - ExpectedX ) * ( Corr[i]( ReadY, 0 ) - ExpectedY );
      z /= NSamples;
      covar(x, y) = z;
      if( x != y )
        covar(y, x) = z;
    }
  }
  assert( Common::IsFinite( covar ) && "Error: covariance matrix isn't finite" );
  std::cout << "done.\n";
  return covar;
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
  for( int src = 0; src < NumOps; ++src )
    for( int snk = 0; snk < NumOps; ++snk )
      for( int t = 0; t < NtCorr; ++t ) {
        double z = 0;
        double CumulativeEnergy = 0;
        for( int e = 0; e < NumExponents; ++e ) {
          CumulativeEnergy -= Energy[e];
          z += Coeff[MELIndex( src, e )] * Coeff[MELIndex( snk, e )] * std::exp( CumulativeEnergy * t );
        }
        const int iRead{ AlphaIndex( src, snk ) * Nt + t + tMin };
        const int iWrite{ AlphaIndex( src, snk ) * NtCorr + t };
        ModelError[iWrite] = Corr[Latan::central]( iRead, 0 ) - z;
      }
  assert( Common::IsFinite( ModelError, Extent ) && "Error: Model errors aren't finite" );

  double chi2 = 0;
  // TODO: I think this should be symmetric - only calculate half the non-diagonals if so
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j < Extent; ++j )
      chi2 += ModelError[i] * CovarInv(i, j) * ModelError[j];
  return chi2;
}

// Perform a fit
void MultiExpModel::PerformFit( void )
{
  const int NeedSize{ NumExponents * ( 1 + NumOps ) };
  std::vector<double> par( NeedSize );
  double * const Energy{ par.data() };
  double * const Coeff{ Energy + NumExponents };
  double CoeffGuess;
  // Take a starting guess for the parameters - same as LatAnalyze
  {
    const int tGuess{Nt / 4};
    const Latan::DMat & m{ Corr[Latan::central] };
    Energy[0] = std::log( m( tGuess, 0 ) / m( tGuess + 1, 0 ) );
    CoeffGuess = m( tGuess, 0 ) / ( std::exp( - Energy[0] * tGuess ) );
    for( int o = 0; o < NumOps; ++o )
      Coeff[MELIndex(o, 0)] = std::sqrt( CoeffGuess / NumOps );
  }
  for( int e = 1; e < NumExponents; ++e ) {
    //Energy[e] = 2 * Energy[e - 1]; // Absolute
    Energy[e] = Energy[e - 1] * ( 1 << ( e - 1 ) ); // Cumulative
    CoeffGuess *= 0.5;
    for( int o = 0; o < NumOps; ++o )
      Coeff[MELIndex(o, e)] = std::sqrt( CoeffGuess / NumOps );
  }
  std::vector<double> Error( NeedSize, 0. );
  std::cout << "Initial guess:\n";
  for( int i = 0; i < par.size(); ++i )
    std::cout << "\tParameter[" << i << "]=" << par[i] << ", +/- " << Error[i] << "\n";

  // Create minimisation parameters
  ROOT::Minuit2::MnUserParameters upar( par, Error );
  //upar.SetPrecision(0.001);
  //upar.SetLowerLimit( 0, 0. );
  for( int e = 0; e < NumExponents; ++e )
    upar.SetLowerLimit( e, 0. );

  //ROOT::Minuit2::VariableMetricMinimizer minimizer2;
  //ROOT::Minuit2::FunctionMinimum min = minimizer2.Minimize( *this, par, Error );
  ROOT::Minuit2::MnMigrad minimizer( *this, upar );
  ROOT::Minuit2::FunctionMinimum min = minimizer();
  ROOT::Minuit2::MnUserParameterState state = min.UserState();
  std::cout << "state.IsValid() = " << state.IsValid() << "\n";
  ROOT::Minuit2::MnUserParameters params = min.UserParameters();
  std::cout << "Fit result:" << state << std::endl;
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
                "minimizer verbosity level (0|1|2)",
#ifdef _DEBUG
                "2");
#else
                "0");
#endif
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
  if (!opt.parse(argc, argv) or (opt.getArgs().size() < 1) or opt.gotOption("help"))
  {
    std::cerr << "usage: " << argv[0] << " <options> <correlator file>" << std::endl;
    std::cerr << std::endl << "Possible options:" << std::endl << opt << std::endl;
    
    return EXIT_FAILURE;
  }
  int ti                  = opt.optionValue<int>("ti");
  int tf                  = opt.optionValue<int>("tf");
  //int thinning            = opt.optionValue<int>("t");
  int shift               = opt.optionValue<int>("s");
  std::string model       = opt.optionValue("m");
  //Index nPar              = opt.optionValue<Index>("nPar");
  //double svdTol           = opt.optionValue<double>("svd");
  std::string outFileName = opt.optionValue("o");
  //bool doCorr             = !opt.gotOption("uncorr");
  int fold                = opt.optionValue<int>("fold");
  //bool doPlot             = opt.gotOption("p");
  //double plotrange        = opt.optionValue<Index>("range");
  //bool doHeatmap          = opt.gotOption("h");
  Minimizer::Verbosity verbosity;
  switch (opt.optionValue<unsigned int>("v"))
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
  }
  
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
  MultiExpModel m ( Common::ReadBootstrapCorrs( FileName, fold, shift, NumOps ), ti, tf, NumOps, NumExponents );
  m.PerformFit();

  return EXIT_SUCCESS;
}
