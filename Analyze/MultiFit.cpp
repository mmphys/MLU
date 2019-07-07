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
#include <Minuit2/FunctionMinimum.h>

// Minimisation for multi-exponential fit

class MultiExpModel : public ROOT::Minuit2::FCNBase {
public:
  using vd_t = std::vector<double>;
  using cdouble = std::complex<double>;
  explicit MultiExpModel( const Latan::DMatSample Corr, int NumOps, int tMin, int tMax, int NumExponents = 2 );
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
  const int NSamples;
  const int Nt;
  const int tMin, tMax;
  const Latan::DMatSample Corr;
  const Latan::DMat Covar;
  const Latan::DMat CovarInv;
};

MultiExpModel::MultiExpModel( const Latan::DMatSample corr_, int numOps_, int TMin_, int TMax_, int numExponents_ )
  : NumExponents{ numExponents_ },
    NumOps{ numOps_ }, tMin{ TMin_ }, tMax{ TMax_ },
    NSamples{ static_cast<int>( corr_.size() ) },
    Nt      { static_cast<int>( corr_[Latan::central].rows() / ( NumOps * NumOps ) ) },
    Corr{ corr_ },
    Covar{ Corr.covarianceMatrix( Corr ) }, // TODO: only construct this for the timeslices we need
    CovarInv{ Covar.inverse() }
{
  assert( numExponents_ > 0 && "Number of exponents in the fit should be positive" );
  assert( corr_[Latan::central].rows() == NumOps * NumOps * Nt && "Bug" );
  assert( tMin >= 0 && tMin < Nt && "Error: tMin invalid" );
  assert( tMax >= 0 && tMax < Nt && "Error: tMax invalid" );
  assert( tMin < tMax && "Error: tMin >= tMax" );
  std::cout << "Covariance matrix is " << Covar.rows() << " x " << Covar.cols() << "\n";
  std::cout << "Inverse covariance matrix is " << CovarInv.rows() << " x " << CovarInv.cols() << "\n";
}

// Compute chi-squared given parameters for multi-exponential fit
  // Energy Levels        Parameters 0 ... NumExponents - 1
  // Overlap Coefficients Parameters NumExponents ... NumExponents ( NumExponents + 1 ) - 1

double MultiExpModel::operator()( const vd_t & par ) const {
  const int NeedSize{ NumExponents + 2 * NumExponents * NumOps };
  assert( par.size() == NeedSize && "Error: number of parameters" );
  const double * const Energy{ par.data() };
  const cdouble * const Coeff{ reinterpret_cast<const cdouble *> ( Energy + NumExponents ) };
  // Calculate the theory errors for these model parameters
  const int Extent{ NumOps * NumOps * Nt };
  double ModelError[Extent]; // Should happily fit on stack
  for( int src = 0; src < NumOps; ++src )
    for( int snk = 0; snk < NumOps; ++snk )
      for( int t = tMin; t < tMax; ++t ) {
        const int i{ ( src * NumOps + snk ) * NumOps + t };
        cdouble z = 0;
        for( int e = 0; e < NumExponents; ++e )
          z += Coeff[src * NumExponents + e] * std::conj( Coeff[snk * NumExponents + e] ) * std::exp(0 - Energy[e] * t );
        ModelError[i] = Corr[Latan::central]( i, 0 ) - z.real();
      }

  double chi2 = 0;
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j < Extent; ++j )
      chi2 += ModelError[i] * CovarInv(i, j) * ModelError[j];
  return chi2;
}

// Perform a fit
void MultiExpModel::PerformFit( void )
{
  const int NeedSize{ NumExponents + 2 * NumExponents * NumOps };
  std::vector<double> par( NeedSize );
  double * const Energy{ par.data() };
  cdouble * const Coeff{ reinterpret_cast<cdouble *> ( Energy + NumExponents ) };
  // Take a starting guess for the parameters - same as LatAnalyze
  {
    const int tGuess{Nt / 4};
    const Latan::DMat & m{ Corr[Latan::central] };
    Energy[0] = std::log( m( tGuess, 0 ) / m( tGuess + 1, 0 ) );
    Coeff[0] = m( tGuess, 0 ) / ( std::exp( -Energy[0] * tGuess ) );
  }
  for( int i = 1; i < NumExponents; ++i ) {
    Energy[i] = 2 * Energy[i - 1];
    Coeff[i] = Coeff[i - 1] * 0.5;
  }
  std::vector<double> Error( NeedSize, 0.1 );
  ROOT::Minuit2::VariableMetricMinimizer minimizer;
  ROOT::Minuit2::FunctionMinimum min = minimizer.Minimize( *this, par, Error );
  //std::cout << "Minimum: " << min << std::endl;
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
  OptParser            opt;
  bool                 parsed, doPlot, doHeatmap, doCorr;
  std::string          model, outFileName, outFmt;
  Index                ti, tf, nPar, thinning;
  double               svdTol,plotrange;
  Minimizer::Verbosity verbosity;
  
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
  parsed = opt.parse(argc, argv);
  if (!parsed or (opt.getArgs().size() < 1) or opt.gotOption("help"))
  {
    std::cerr << "usage: " << argv[0] << " <options> <correlator file>" << std::endl;
    std::cerr << std::endl << "Possible options:" << std::endl << opt << std::endl;
    
    return EXIT_FAILURE;
  }
  ti           = opt.optionValue<Index>("ti");
  tf           = opt.optionValue<Index>("tf");
  thinning     = opt.optionValue<Index>("t");
  int shift    = opt.optionValue<int>("s");
  model        = opt.optionValue("m");
  nPar         = opt.optionValue<Index>("nPar");
  svdTol       = opt.optionValue<double>("svd");
  outFileName  = opt.optionValue<std::string>("o");
  doCorr       = !opt.gotOption("uncorr");
  int fold     = opt.optionValue<int>("fold");
  doPlot       = opt.gotOption("p");
  plotrange    = opt.optionValue<Index>("range");
  doHeatmap    = opt.gotOption("h");
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
  MultiExpModel m ( Common::ReadBootstrapCorrs( FileName, fold, shift ), NumOps, ti, tf, NumExponents );
  m.PerformFit();

  return EXIT_SUCCESS;
}
