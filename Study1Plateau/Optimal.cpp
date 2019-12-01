//
//  Optimal.cpp
//  Utility
//
//  Created by Michael Marshall on 25/11/2019.
//  Copyright © 2019 sopa. All rights reserved.
//

#include "Optimal.hpp"

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
  opt.addOption("" , "save"   , OptParser::OptType::trigger, true,
                "save bootstrap replicas of correlators or model");
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
  const bool bSave{ opt.gotOption("save") };
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
    Latan::DMatSample Corr=Common::ReadBootstrapCorrs( opt.getArgs(), fold, shift, NumOps );
    std::string sSummaryBase{ outBaseFileName };
    sSummaryBase.append( 1, '.' );
    if( doCorr )
      sSummaryBase.append( "corr" );
    else
      sSummaryBase.append( "uncorr" );
    std::string OutputRunBase{sSummaryBase};
    sSummaryBase.append( 1, '.' );
    sSummaryBase.append( OpNames[0] );
    for( int op = 1; op < NumOps; op++ ) {
      sSummaryBase.append( 1, '_' );
      sSummaryBase.append( OpNames[op] );
    }
    OutputRunBase.append( 1, '_' );
    OutputRunBase.append( std::to_string( ti ) );
    OutputRunBase.append( 1, '_' );
    OutputRunBase.append( std::to_string( tf ) );
    static const char Sep[] = " ";
    const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
    if( !Common::FileExists( sFitFilename ) )
      throw std::runtime_error( "Output file " + sFitFilename + " already exists" );
    const int NumFiles{ NumOps * NumOps };
    const int NSamples{static_cast<int>(Corr.size())};
    const int Nt      { static_cast<int>( Corr[Latan::central].rows() / NumFiles ) };
    const int NumParams{ NumExponents * ( 1 + NumOps ) };
    std::string sModelBase{ OutputRunBase };
    sModelBase.append( 1, '.' );
    sModelBase.append( OpNames[0] );
    for( std::size_t i = 1; i < OpNames.size(); i++ ) {
      sModelBase.append( 1, '_' );
      sModelBase.append( OpNames[i] );
    }
    Latan::DMatSample ModelParams( NSamples, NumParams, 2 );
    ModelParams = Latan::Io::load<Latan::DMatSample>(Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ));
    if( NumOps == 2) {  // I've only coded for two operators
      Latan::DMatSample OptimalCorr( NSamples, Nt, 2 );
      for( int degrees = -90; degrees <= 90; degrees ++ )
      {
        double costheta, sintheta;
        switch( ( degrees + 360 ) % 360 ) {
          case 0:
            costheta = 1;
            sintheta = 0;
            break;
          case 90:
            costheta = 0;
            sintheta = 1;
            break;
          case 180:
            costheta = -1;
            sintheta = 0;
            break;
          case 270:
            costheta = 0;
            sintheta = -1;
            break;
          default:
          {
            const double theta{ M_PI * degrees / 180 };
            costheta = cos( theta );
            sintheta = sin( theta );
          }
        }
        const double cos_sq_theta{ costheta * costheta };
        const double sin_sq_theta{ sintheta * sintheta };
        const double cos_sin_theta{ costheta * sintheta };
        const double A_P1{ ModelParams[Latan::central](4, 0) };
        const double A_A1{ ModelParams[Latan::central](5, 0) };
        for (Latan::Index i = Latan::central; i < NSamples; i++) {
          const double Op_PP{ cos_sq_theta / ( A_P1 * A_P1 ) };
          const double Op_AP{ cos_sin_theta / ( A_P1 * A_A1 ) };
          const double Op_AA{ sin_sq_theta / ( A_A1 * A_A1 ) };
          for( int t = 0; t < Nt; t++ ) {
            double z = Op_PP * Corr[i](0 * Nt + t,0) + Op_AA * Corr[i](3 * Nt + t,0)
                                + Op_AP * ( Corr[i](1 * Nt + t,0) + Corr[i](2 * Nt + t,0) );
            OptimalCorr[i](t,0) = z;
            OptimalCorr[i](t,1) = 0;
          }
        }
        std::string sOptimusPrime{ OutputRunBase };
        sOptimusPrime.append( ".theta_" );
        sOptimusPrime.append( std::to_string(degrees) );
        if( bSave )
          Latan::Io::save( OptimalCorr, Common::MakeFilename( sOptimusPrime, Common::sBootstrap, Seed, DEF_FMT ) );
        Common::SummariseBootstrapCorr( OptimalCorr, sOptimusPrime, Seed );
      }
    }
  }
  catch(const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
