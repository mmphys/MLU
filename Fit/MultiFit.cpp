/*************************************************************************************
 
 Fast (using OpenMP) multi-model fits to lattice QCD correlators

 Source file: MultiFit.cpp
 
 Copyright (C) 2019-2022
 
 Author: Michael Marshall <Mike@lqcd.me>

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

#include "Fitter.hpp"
#include <MLU/DebugInfo.hpp>
#include <chrono>

// Uncomment the next line to disable Minuit2 (for testing)
//#undef HAVE_MINUIT2

// Indices for operators in correlator names
const char * pSrcSnk[2] = { "src", "snk" };

// This should be the only place which knows about different fitters

Fitter * MakeFitterGSL( const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                        std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                        CovarParams &&cp );
#ifdef HAVE_MINUIT2
Fitter * MakeFitterMinuit2(const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                           std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                           CovarParams &&cp );
#endif

Fitter * MakeFitter( const Common::CommandLine &cl, const DataSet &ds,
                     std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                     CovarParams &&cp )
{
  Fitter * f;
  std::string FitterArgs{ cl.SwitchValue<std::string>( "fitter" ) };
  std::string FitterType{ Common::ExtractToSeparator( FitterArgs ) };
  if( Common::EqualIgnoreCase( FitterType, "GSL" ) )
    f = MakeFitterGSL( FitterArgs, cl, ds, ModelArgs, opNames, std::move( cp ) );
#ifdef HAVE_MINUIT2
  else if( Common::EqualIgnoreCase( FitterType, "Minuit2" ) )
    f = MakeFitterMinuit2( FitterArgs, cl, ds, ModelArgs, opNames, std::move( cp ) );
#endif
  else
    throw std::runtime_error( "Unrecognised fitter: " + FitterType );
  if( !f )
    throw std::runtime_error( "Unrecognised " + FitterType + " fitter options: " + FitterArgs );
  return f;
}

int main(int argc, const char *argv[])
{
#ifndef HAVE_MINUIT2
  std::ios_base::sync_with_stdio( false );
#endif
  static const char DefaultEnergySep[] = "0"; // Default used to be 0.2 until 10 Jul 2021
  static const char DefaultHotelling[] = "0.05";
  static const char DefaultMinDP[] = "3";
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"analytic", CL::SwitchType::Flag, nullptr},
      {"testrun", CL::SwitchType::Flag, nullptr},
      {"overwrite", CL::SwitchType::Flag, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"chisqdof", CL::SwitchType::Single, "0"},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"mindof", CL::SwitchType::Single, "1"},
      {"mindp", CL::SwitchType::Single, DefaultMinDP},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"summary", CL::SwitchType::Single, "1"},
      {"savecmat", CL::SwitchType::Flag, nullptr},
      {"v", CL::SwitchType::Single, "0"},
      {"guess", CL::SwitchType::Single, nullptr},
      // ModelDefaultParams
      {"e", CL::SwitchType::Single, "1"},
      {Common::sOverlapAltNorm.c_str(), CL::SwitchType::Flag, nullptr},
      // Covariance parameters
      {"covsrc", CL::SwitchType::Single, "Bootstrap"},
      {"covboot", CL::SwitchType::Single, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      // Other params
      {"fitter", CL::SwitchType::Single, "GSL"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"opnames", CL::SwitchType::Flag, nullptr},
      {"debug-signals", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      if( cl.GotSwitch( "debug-signals" ) )
        Common::Grid_debug_handler_init();
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      Common::MakeAncestorDirs( outBaseFileName );
      const int NSamples{ cl.SwitchValue<int>("n") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bOpSort{ !cl.GotSwitch("opnames") };

      // Walk the list of parameters on the command-line, loading correlators and making models
      bShowUsage = false;
      std::vector<std::string> OpName;
      std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ ) << "Loading folded correlators\n";
      // Split each argument at the first comma (so the first part can be treated as a filename to glob
      const std::size_t NumArgs{ cl.Args.size() };
      DataSet ds( NSamples );
      std::vector<Model::Args> ModelArgs;
      std::vector<int> ModelFitRange;
      std::vector<std::string> FitRangeSpec;
      for( std::size_t ArgNum = 0; ArgNum < NumArgs; ++ArgNum )
      {
        // First parameter (up to comma) is the filename we're looking for
        std::string FileToGlob{ Common::ExtractToSeparator( cl.Args[ArgNum] ) };
        bool bGlobEmpty{ true };
        // Anything after the comma is a list of arguments
        Model::Args vArgs;
        vArgs.FromString( cl.Args[ArgNum], true );
        for( const std::string &sFileName : Common::glob( &FileToGlob, &FileToGlob + 1, inBase.c_str() ) )
        {
          bGlobEmpty = false;
          Common::FileNameAtt Att( sFileName, &OpName );
          const bool bIsCorr{ Common::EqualIgnoreCase( Att.Type, Common::sFold ) };
          if( bIsCorr )
          {
            // This is a correlator - load it
            std::string PrintPrefix( 2, ' ' );
            if( !cl.Args[ArgNum].empty() )
            {
              PrintPrefix.append( cl.Args[ArgNum] );
              PrintPrefix.append( 1, ' ' );
            }
            ds.LoadCorrelator( std::move( Att ), Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT,
                               PrintPrefix.c_str() );
            // Get the fit range (if present)
            int FitRange;
            bool bGotFitRange;
            std::string ThisFitRange = vArgs.Remove( "t", &bGotFitRange );
            if( bGotFitRange )
            {
              std::istringstream is( ThisFitRange );
              if( is >> FitRange && ( is.eof() || ( is >> std::ws && is.eof() ) ) )
                ; // Use the single number fit range we just decoded
              else
              {
                // Treat this as a new fit range spec
                FitRange = static_cast<int>( FitRangeSpec.size() );
                FitRangeSpec.emplace_back( std::move( ThisFitRange ) );
              }
            }
            else
            {
              if( ModelFitRange.empty() )
                throw std::runtime_error( "The first correlator must specify a fit range" );
              FitRange = ModelFitRange.back(); // Use same fit range as previous correlator
            }
            ModelFitRange.push_back( FitRange );
            ModelArgs.emplace_back( vArgs );
          }
          else
          {
            ds.LoadModel( std::move( Att ), cl.Args[ArgNum] );
          }
        }
        if( bGlobEmpty )
          std::cout << "Warning: No files matched " << FileToGlob << Common::NewLine;
      }
      if( ds.corr.empty() )
        throw std::runtime_error( "At least one correlator must be loaded to perform a fit" );
      // Make sure models refer to all fit ranges
      const int MinDP{ cl.SwitchValue<int>( "mindp" ) };
      Common::FitRanges fitRanges( FitRangeSpec, MinDP );
      {
        std::vector<bool> ModelFitRangeUsed( fitRanges.size(), false );
        for( std::size_t i = 0; i < ds.corr.size(); ++i )
        {
          if( ModelFitRange[i] < 0 || ModelFitRange[i] >= ModelFitRangeUsed.size() )
          {
            std::ostringstream es;
            es << "Model " << i << " refers to non-existent fit range " << ModelFitRange[i];
            throw std::runtime_error( es.str().c_str() );
          }
          ModelFitRangeUsed[ModelFitRange[i]] = true;
        }
        for( std::size_t i = 0; i < ModelFitRangeUsed.size(); ++i )
          if( !ModelFitRangeUsed[i] )
            throw std::runtime_error( "Fit range " + std::to_string( i ) + " not referred to by any model" );
      }
      // Check whether each fit range is valid for each model using it
      for( std::size_t i = 0; i < ds.corr.size(); ++i )
      {
        if( !fitRanges[ModelFitRange[i]].Validate( ds.corr[i].Nt() ) )
        {
          std::stringstream oss;
          oss << "Fit range " << fitRanges[ModelFitRange[i]]
              << " not valid for correlator " << ds.corr[i].Name_.Filename;
          throw std::runtime_error( oss.str().c_str() );
        }
      }
      // Describe the number of replicas
      ds.SortOpNames( OpName );
      std::cout << "Using ";
      if( ds.NSamples == ds.MaxSamples )
        std::cout << "all ";
      else
        std::cout << "first " << ds.NSamples << " of ";
      std::cout << ds.MaxSamples << " bootstrap replicas";
      if( ds.NSamples != ds.MaxSamples )
        std::cout << " (all " << ds.MaxSamples << " for var/covar)";
      std::cout << Common::NewLine;

      // Work out how covariance matrix should be constructed and tell user
      CovarParams cp{ cl, ds };
      std::cout << cp << Common::NewLine;

      // I'll need a sorted, concatenated list of the operators in the fit for filenames
      std::string sOpNameConcat;
      if( bOpSort )
      {
        sOpNameConcat = OpName[0];
        for( std::size_t i = 1; i < OpName.size(); i++ )
        {
          sOpNameConcat.append( 1, '_' );
          sOpNameConcat.append( OpName[i] );
        }
      }
      else
      {
        for( std::size_t i = 0; i < ds.corr.size(); ++i )
        {
          const std::size_t ThisSize{ ds.corr[i].Name_.op.size() };
          for( std::size_t j = ThisSize; j; )
          {
            if( i || j != ThisSize )
              sOpNameConcat.append( Common::Underscore );
            sOpNameConcat.append( OpName[ds.corr[i].Name_.op[--j]] );
          }
        }
      }

      // Now make base filenames for output
      // If the output base is given, but doesn't end in '/', then it already contains the correct name
      if( outBaseFileName.empty() || outBaseFileName.back() == '/' )
        outBaseFileName.append( ds.corr[0].Name_.Base ); // Very simplistic, but better is hard to automate
      outBaseFileName.append( 1, '.' );
      outBaseFileName.append( doCorr ? "corr" : "uncorr" );
      const std::string sSummaryBase{ outBaseFileName + Common::Period + sOpNameConcat };
      outBaseFileName.append( 1, '_' );
      const std::size_t outBaseFileNameLen{ outBaseFileName.size() };
      const Common::SeedType Seed{ ds.corr[0].Name_.Seed };

      // All the models are loaded
      std::unique_ptr<Fitter> m{ MakeFitter( cl, ds, ModelArgs, OpName, std::move( cp ) ) };
      std::size_t CountTotal{ 0 };
      std::size_t CountOK{ 0 };
      const auto Start{ std::chrono::steady_clock::now() };
      for( Common::FitRangesIterator it = fitRanges.begin(); !it.PastEnd(); ++it )
      {
        // Log what file we're processing and when we started
        const auto then{ std::chrono::steady_clock::now() };
        ++CountTotal;
        // Set fit ranges
        std::vector<std::vector<int>> fitTimes( ds.corr.size() );
        for( int i = 0; i < ds.corr.size(); ++i )
        {
          const Common::FitTime &ft{ it[ModelFitRange[i]] };
          const int Extent{ ft.tf - ft.ti + 1 };
          fitTimes[i].resize( Extent );
          for( int j = 0; j < Extent; ++j )
            fitTimes[i][j] = ft.ti + j;
        }
        ds.SetFitTimes( fitTimes );
        try
        {
          {
            std::stringstream ss;
            ss << ( doCorr ? "C" : "unc" ) << "orrelated " << m->Type() << " fit on timeslices " << it.to_string( "-", ", " );
            Fitter::SayNumThreads( ss );
            const std::string &sMsg{ ss.str() };
            if( !m->bTestRun )
              std::cout << std::string( sMsg.length(), '=' ) << Common::NewLine;
            std::cout << sMsg << Common::NewLine;
          }
          double ChiSq;
          int dof;
          outBaseFileName.resize( outBaseFileNameLen );
          outBaseFileName.append( it.to_string( Common::Underscore ) );
          outBaseFileName.append( 1, '.' );
          m->PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat, Seed );
          ++CountOK;
        }
        catch(const std::exception &e)
        {
          std::cout << "Error: " << e.what() << "\n";
        }
        // Mention that we're finished, what the time is and how long it took
        if( !m->bTestRun )
        {
          const auto now{ std::chrono::steady_clock::now() };
          const auto mS{ std::chrono::duration_cast<std::chrono::milliseconds>( now - then ) };
          std::cout << "Fit duration ";
          if( mS.count() < 1100 )
            std::cout << mS.count() << " milliseconds.\n";
          else
          {
            const auto S{ std::chrono::duration_cast<std::chrono::duration<double>>( now - then ) };
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(1) << S.count() << " seconds.\n";
            std::cout << ss.str();
          }
        }
      }
      const auto now{ std::chrono::steady_clock::now() };
      auto S{ std::chrono::duration_cast<std::chrono::seconds>( now - Start ).count() };
      std::cout << std::string( 50, '=' ) << Common::NewLine
                << CountOK << " fits succeeded of " << CountTotal << " attempted\n"
                << "Total duration ";
      constexpr unsigned int DaySeconds{ 24 * 60 * 60 };
      if( S >= DaySeconds )
      {
        std::cout << S / DaySeconds << '-';
        S %= DaySeconds;
      }
      const int Hours{ static_cast<int>( S / 3600 ) };
      S %= 3600;
      const int Minutes{ static_cast<int>( S / 60 ) };
      S %= 60;
      std::cout << std::setfill( '0' ) << std::setw(2) << Hours
         << ':' << std::setfill( '0' ) << std::setw(2) << Minutes
         << ':' << std::setfill( '0' ) << std::setw(2) << S << Common::NewLine;
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
    " [Options] BootstrapOrModel[,params]...\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, \n"
    "reading in parameters from previous fit results (Models).\n"
    "Options:\n"
    "--Hotelling Minimum Hotelling Q-value on central replica (default " << DefaultHotelling << ")\n"
    "--chisqdof  Maximum chi^2 / dof on central replica\n"
    "--sep    Minimum relative separation between energy levels (default " << DefaultEnergySep << ")\n"
    "--retry  Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter   Max iteration count, 0 (default) = unlimited\n"
    "--tol    Tolerance of required fits (default 1e-7)\n"
    "--mindof Minimum degrees of freedom (default 1)\n"
    "--mindp  Minimum number of data points per fit range (default " << DefaultMinDP << ")\n"
    "--fitter (GSL|Minuit2)[,options] fitter (default GSL,Levenberg-Marquardt)\n"
    "         GSL options: lm, lmaccel, dogleg, ddogleg, subspace2D\n"
    "--covsrc source[,options[,...]] build (co)variance from source, i.e. one of\n"
    "         Binned    Use the already binned data\n"
    "         Raw       Use the raw (unbinned) data\n"
    "         Rebin     Rebin the raw data using bin size(s) specified\n"
    "         Bootstrap Use bootstrap replicas (default)\n"
    "         H5,f[,g],d Load INVERSE covariance from .h5 file f, group g, dataset d\n"
    "--covboot How many bootstrap replicas in covariance (-1=no bootstrap)\n"
    "--guess   List of specific values to use for inital guess\n"
    "--summary 0 no summaries; 1 model_td.seq.txt only; 2 model_td and model.seq.txt\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecmat Save correlation matrix\n"
    "--analytic Analytic derivatives for GSL (default: numeric)\n"
    "--opnames  Disable sorting and deduplicating operator name list\n"
    "--testrun  Don't perform fits - just say which fits would be attempted\n"
    "--overwrite Overwite always. Default: only overwrite smaller Nboot\n"
    "--debug-signals Trap signals (code courtesy of Grid)\n"
    "--" << Common::sOverlapAltNorm << "  Alternate normalisation for overlap factors. DEPRECATED\n"
    "--help     This message\n"
    "Parameters accepted by all models:\n"
    " e         Number of exponentials\n"
    " model     Model type {Exp, Cosh, Sinh, ThreePoint, R3, Const}\n"
    " t         Fit range: [R[n:]start:stop[:numstart=1[:numstop=numstart]])\n"
    "Parameters accepted by models with overlap coefficients:\n"
    " SrcSnk    Force source and sink to be different (by appending 'src' and 'snk')\n"
    "Parameters accepted by models for single objects (e.g. 2pt-functions):\n"
    " ObjectID  Object identifier (defaults to base of filename)\n"
    "Parameters accepted by models for dual objects (e.g. 3pt-functions):\n"
    " Src       ObjectID for source\n"
    " Snk       ObjectID for sink\n"
    "Parameters accepted by R3 model:\n"
    " C2Model   Which model to use for 2pt: Exp (default); Cosh; Sinh\n";
  }
  Common::Grid_exit_handler_disable = true;
  return iReturn;
}
