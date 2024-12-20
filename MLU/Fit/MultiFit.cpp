/**
 
 Fast (using OpenMP) multi-model fits to lattice QCD correlators

 Source file: MultiFit.cpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "Fitter.hpp"
#include "ModelCommon.hpp"
#include <MLU/DebugInfo.hpp>
#include <chrono>

// Indices for operators in correlator names
const char * pSrcSnk[2] = { "src", "snk" };

int main(int argc, const char *argv[])
{
#ifndef HAVE_MINUIT2
  std::ios_base::sync_with_stdio( false );
#endif
  static const char DefaultEnergySep[] = "0"; // Default used to be 0.2 until 10 Jul 2021
  static const char DefaultHotelling[] = "0.05";
  static const char DefaultMinDP[] = "1";
  static const char DefaultErrDig[] = "2";
  static const char szFileWildcard[] = "*?[";
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"analytic", CL::SwitchType::Flag, nullptr},
      {"testrun", CL::SwitchType::Flag, nullptr},
      {"central", CL::SwitchType::Flag, nullptr},
      {"overwrite", CL::SwitchType::Flag, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"chisqdof", CL::SwitchType::Single, "0"},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"shrink", CL::SwitchType::Single, "0"},
      {"mindof", CL::SwitchType::Single, "1"},
      {"mindp", CL::SwitchType::Single, DefaultMinDP},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"summary", CL::SwitchType::Single, "1"},
      {"nopolap", CL::SwitchType::Single, ""},
      {"savecmat", CL::SwitchType::Flag, nullptr},
      {"v", CL::SwitchType::Single, "0"},
      {"guess", CL::SwitchType::Single, nullptr},
      {"strict",  CL::SwitchType::Single, "0"},
      {"maxE",  CL::SwitchType::Single, "10"},
      {"errdig", CL::SwitchType::Single, DefaultErrDig},
      // ModelDefaultParams
      {"e", CL::SwitchType::Single, "1"},
      {"N", CL::SwitchType::Single, "0"},
      {"dispersion", CL::SwitchType::Single, nullptr},
      {MLU::sOverlapAltNorm.c_str(), CL::SwitchType::Flag, nullptr},
      // Covariance parameters
      {"covsrc", CL::SwitchType::Single, nullptr},
      {"covboot", CL::SwitchType::Single, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      // Other params
      {"fitter", CL::SwitchType::Single, "GSL"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"extra", CL::SwitchType::Single, ""},
      {"showname", CL::SwitchType::Flag, nullptr},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"opnames", CL::SwitchType::Flag, nullptr},
      {"debug-signals", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      if( cl.GotSwitch( "debug-signals" ) )
        MLU::Grid_debug_handler_init();
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      const bool inBaseWildcard{ inBase.find_first_of( szFileWildcard ) != std::string::npos };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      // If output base is given, but doesn't end in '/', then it already contains the correct name
      bool bAppendCorr0Name{ outBaseFileName.empty() || outBaseFileName.back() == '/' };
      MLU::MakeAncestorDirs( outBaseFileName );
      const int NSamples{ cl.SwitchValue<int>("n") };
      const std::string NameExtra{ cl.SwitchValue<std::string>("extra") };
      const bool bShowName{ cl.GotSwitch("showname") }; // Show base output name - don't run
      const bool bTestRun{ cl.GotSwitch("testrun") }; // Load files and say which fits attempted
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bOpSort{ !cl.GotSwitch("opnames") };

      if( bShowName && !bOpSort )
        throw std::runtime_error( "--showname is incompatible with --opnames" );
      if( bShowName && bTestRun )
        throw std::runtime_error( "--showname is incompatible with --testrun" );
      // Walk the list of parameters on the command-line, loading correlators and making models
      bShowUsage = false;
      std::vector<std::string> OpName;
      if( !bShowName )
        std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ )
                  << "Loading folded correlators\n";
      // Split each argument at the first comma (so the first part can be treated as a filename to glob
      const std::size_t NumArgs{ cl.Args.size() };
      DataSet ds( NSamples );
      std::vector<Model::Args> ModelArgs;
      std::vector<int> ModelFitRange;
      std::vector<std::string> FitRangeSpec;
      bool bGotCorrelator{ false };
      for( std::size_t ArgNum = 0; ArgNum < NumArgs; ++ArgNum )
      {
        // First parameter (up to comma) is the filename we're looking for
        std::string FileToGlob{ MLU::ExtractToSeparator( cl.Args[ArgNum] ) };
        std::vector<std::string> Filenames;
        if( !inBaseWildcard && FileToGlob.find_first_of( szFileWildcard ) == std::string::npos )
          Filenames.push_back( MLU::PreferSeed( inBase + FileToGlob ) );
        else
          Filenames = MLU::glob( &FileToGlob, &FileToGlob + 1, inBase.c_str() );
        bool bGlobEmpty{ true };
        // Anything after the comma is a list of arguments
        Model::Args vArgs;
        vArgs.FromString( cl.Args[ArgNum], true );
        for( const std::string &sFileName : Filenames )
        {
          bGlobEmpty = false;
          MLU::FileNameAtt Att( sFileName, &OpName );
          const bool bIsCorr{ MLU::EqualIgnoreCase( Att.Type, MLU::sFold ) };
          if( bIsCorr )
          {
            bGotCorrelator = true;
            if( bAppendCorr0Name )
            {
              bAppendCorr0Name = false;
              outBaseFileName.append( Att.GetBaseExtra() ); // Simplistic - better hard to automate
            }
            if( !bShowName ) // No need to load if all I need is output name
            {
            // This is a correlator - load it
            std::string PrintPrefix( 2, ' ' );
            if( !cl.Args[ArgNum].empty() )
            {
              PrintPrefix.append( cl.Args[ArgNum] );
              PrintPrefix.append( 1, ' ' );
            }
            ds.LoadCorrelator( std::move( Att ), MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT,
                              PrintPrefix.c_str() );
            }
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
          else if( !bShowName ) // No need to load if all I need is output name
          {
            ds.LoadModel( std::move( Att ), cl.Args[ArgNum] );
          }
        }
        if( bGlobEmpty )
          throw std::runtime_error( "No files matched " + FileToGlob );
      }
      if( !bGotCorrelator )
        throw std::runtime_error( "At least one correlator must be loaded to perform a fit" );
      // Make sure models refer to all fit ranges
      const int MinDP{ cl.SwitchValue<int>( "mindp" ) };
      MLU::FitRanges fitRanges( FitRangeSpec, MinDP );
      {
        std::vector<bool> ModelFitRangeUsed( fitRanges.size(), false );
        for( std::size_t i = 0; i < ModelFitRange.size(); ++i )
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
      // DataSet will be empty if we're just getting output name - so won't trigger checks/errors
      for( std::size_t i = 0; i < ds.corr.size(); ++i )
      {
        if( !fitRanges[ModelFitRange[i]].Validate( ds.corr[i]->Nt() ) )
        {
          std::stringstream oss;
          oss << "Fit range " << fitRanges[ModelFitRange[i]]
              << " not valid for correlator " << ds.corr[i]->Name_.Filename;
          throw std::runtime_error( oss.str().c_str() );
        }
      }
      if( bShowName )
      {
        std::sort( OpName.begin(), OpName.end(), MLU::LessThanIgnoreCase );
      }
      else
      {
      // Describe the number of replicas
      ds.SortOpNames( OpName );
      std::cout << "Using ";
      if( ds.NSamples == ds.MaxSamples )
        std::cout << "all ";
      else
        std::cout << "first " << ds.NSamples << " of ";
      std::cout << ds.MaxSamples << MLU::Space << ds.corr[0]->getData().SeedTypeString()
                << " replicas";
      if( ds.NSamples != ds.MaxSamples )
        std::cout << " (all " << ds.MaxSamples << " for var/covar)";
      std::cout << MLU::NewLine;
      }

      // Work out how covariance matrix should be constructed and tell user
      CovarParams cp{ cl, ds, bShowName };
      if( !bShowName )
        std::cout << cp << MLU::NewLine;

      // I'll need a sorted, concatenated list of the operators in the fit for filenames
      std::string sOpNameConcat;
      if( bOpSort )
      {
        for( std::size_t i = 0; i < OpName.size(); i++ )
        {
          if( i )
            sOpNameConcat.append( 1, '_' );
          sOpNameConcat.append( OpName[i] );
        }
      }
      else
      {
        for( std::size_t i = 0; i < ds.corr.size(); ++i )
        {
          const std::size_t ThisSize{ ds.corr[i]->Name_.op.size() };
          for( std::size_t j = ThisSize; j; )
          {
            if( i || j != ThisSize )
              sOpNameConcat.append( MLU::Underscore );
            sOpNameConcat.append( OpName[ds.corr[i]->Name_.op[--j]] );
          }
        }
      }

      // Now make base filenames for output
      if( !NameExtra.empty() )
      {
        if( NameExtra[0] != '_' )
          outBaseFileName.append( 1, '.' );
        outBaseFileName.append( NameExtra );
      }
      outBaseFileName.append( 1, '.' );
      outBaseFileName.append( doCorr ? "corr" : "uncorr" );
      outBaseFileName.append( 1, '_' );
      const std::size_t outBaseFileNameLen{ outBaseFileName.size() };
      if( bShowName )
      {
        // Just show the resulting filenames
        sOpNameConcat.append( 1, '.' );
        sOpNameConcat.append( MLU::sModel );
        for( MLU::FitRangesIterator it = fitRanges.begin(); !it.PastEnd(); ++it )
        {
          outBaseFileName.resize( outBaseFileNameLen );
          outBaseFileName.append( it.AbbrevString() );
          outBaseFileName.append( 1, '.' );
          outBaseFileName.append( sOpNameConcat );
          std::cout << outBaseFileName << MLU::NewLine;
        }
      }
      else
      {
      // All the models are loaded
      std::unique_ptr<Fitter> m{ Fitter::Make( MultiFitCreateParams{ OpName, cl }, ds,
                                               std::move( ModelArgs ), std::move( cp ), true ) };
      std::size_t CountTotal{ 0 };
      std::size_t CountOK{ 0 };
      bool bAllParamsResolved{ true };
      const auto Start{ std::chrono::steady_clock::now() };
      for( MLU::FitRangesIterator it = fitRanges.begin(); !it.PastEnd(); ++it )
      {
        // Log what file we're processing and when we started
        const auto then{ std::chrono::steady_clock::now() };
        ++CountTotal;
        // Set fit ranges
        std::vector<std::vector<int>> fitTimes( ds.corr.size() );
        for( int i = 0; i < ds.corr.size(); ++i )
          fitTimes[i] = it[ModelFitRange[i]].GetFitTimes();
        ds.SetFitTimes( fitTimes );
        try
        {
          {
            std::stringstream ss;
            if( m->bTestRun )
              ss << "Test run of ";
            ss << (doCorr ? "C" : "Unc") << "orrelated " << m->Type() << " fit on timeslices " << it;
            Fitter::SayNumThreads( ss );
            const std::string &sMsg{ ss.str() };
            if( !m->bTestRun )
              std::cout << std::string( sMsg.length(), '=' ) << MLU::NewLine;
            std::cout << sMsg << MLU::NewLine;
          }
          double ChiSq;
          int dof;
          outBaseFileName.resize( outBaseFileNameLen );
          outBaseFileName.append( it.AbbrevString() );
          outBaseFileName.append( 1, '.' );
          bAllParamsResolved =
          m->PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat );
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
            std::cout << mS.count() << " milliseconds";
          else
          {
            const auto S{ std::chrono::duration_cast<std::chrono::duration<double>>( now - then ) };
            std::ostringstream ss;
            ss << std::fixed << std::setprecision(1) << S.count() << " seconds";
            std::cout << ss.str();
          }
          std::time_t ttNow;
          std::time( &ttNow );
          std::string sNow{ std::ctime( &ttNow ) };
          while( sNow.length() && sNow.back() == '\n' )
            sNow.resize( sNow.length() - 1 );
          std::cout << ". " << sNow << MLU::NewLine;
        }
      }
      const auto now{ std::chrono::steady_clock::now() };
      auto S{ std::chrono::duration_cast<std::chrono::seconds>( now - Start ).count() };
      std::cout << std::string( 50, '=' ) << MLU::NewLine
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
         << ':' << std::setfill( '0' ) << std::setw(2) << S << MLU::NewLine;
      // Error if only one fit requested and it failed ... or all params not resolved
      if( CountTotal == 1 )
      {
        if( CountOK == 0 )
          iReturn = 2;
        else if( !bAllParamsResolved )
          iReturn = 3;
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
    " [Options] BootstrapOrModel[,params]...\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, \n"
    "reading in parameters from previous fit results (Models).\n"
    "Options:\n"
    "--Hotelling Minimum Hotelling Q-value on central replica (default " << DefaultHotelling << ")\n"
    "--chisqdof  Maximum chi^2 / dof on central replica\n"
    "--sep    Minimum relative separation between energy levels (default " << DefaultEnergySep << ")\n"
    "--shrink Ledoit and Wolf shrinkage (default: 0)\n"
    "--retry  Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter   Max iteration count, 0 (default) = unlimited\n"
    "--tol    Tolerance of required fits (default 1e-7)\n"
    "--mindof Minimum degrees of freedom (default 1)\n"
    "--mindp  Minimum number of data points per fit range (default " << DefaultMinDP << ")\n"
    "--fitter (GSL|Minuit2)[,options] fitter (default GSL,Levenberg-Marquardt)\n"
    "         GSL options: lm, lmaccel, dogleg, ddogleg, subspace2D\n"
    "--strict Mask. Check field shows whether parameters are non-zero and unique\n"
    "         Strictness bits: Off=+/- 1 sigma; On=every replica (default 0)\n"
    "         Bit 0: Difference from 0\n"
    "         Bit 1: Difference from other parameters in same series\n"
    "--maxE   Maximum energy (default 10 - decays so fast effectively undetermined)\n"
    "--errdig Number of significant figures in error (default " << DefaultErrDig << ")\n"
    "--covsrc source[,options[,...]] build (co)variance from source, i.e. one of\n"
    "         Binned    Use the already binned data\n"
    "         Bootstrap Use bootstrap replicas\n"
    "         Raw       Use the raw (unbinned) data\n"
    "         Rebin n[,n[...] Rebin raw data using bin size(s) specified, 0=auto\n"
    "         Reboot[NumBoot] Re-bootstrap binned data, optional NumBoot samples\n"
    "                         Can be used with rebin\n"
    "         H5,f[,g],d Load INVERSE covariance from .h5 file f, group g, dataset d\n"
    "         Default: Binned if available, otherwise Bootstrap\n"
    "--covboot Build covariance using secondary bootstrap & this num replicas, 0=all\n"
    "--guess   List of specific values to use for inital guess\n"
    "--summary 0 no summaries; 1 model_td.seq.txt only; 2 model_td and model.seq.txt\n"
    "--nopolap List of overlap coefficients which don't depend on momenta\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "--extra Extra text to include after the base output name\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-N     Use lattice dispersion relation for boosted energies with N=L/a\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecmat Save correlation matrix\n"
    "--analytic Analytic derivatives for GSL (default: numeric)\n"
    "--opnames  Disable sorting and deduplicating operator name list\n"
    "--testrun  Don't perform fits - load files, say which fits would be attempted\n"
    "--showname Don't perform fits - get the output file base name up to '.model'\n"
    "--dispersion Which dispersion relation to use (default: LatFreeScalar)\n"
    "--central  Don't use the central replica as guess for each replica\n"
    "--overwrite Overwite always. Default: only overwrite smaller Nboot\n"
    "--debug-signals Trap signals (code courtesy of Grid)\n"
    "--" << MLU::sOverlapAltNorm << "  Alternate normalisation for overlap factors. DEPRECATED\n"
    "--help     This message\n"
    "Parameters accepted by all models:\n"
    " e         Number of exponentials\n"
    " model     Model type {Exp, Cosh, Sinh, ThreePoint, R3, Const}\n"
    " t         Fit range: [R[n:]start:stop[:numstart=1[:numstop=numstart]])\n"
    "Parameters accepted by models with overlap coefficients:\n"
    " SrcSnk    Force source and sink to be different (by appending 'src' and 'snk')\n"
    "Parameters accepted by models for single objects (e.g. 2pt-functions):\n"
    " ObjectID  Object identifier (defaults to base of filename)\n"
    " Energy    Energy parameter name (defaults to " << MLU::ModelBase::EnergyPrefix << ")\n"
    "Parameters accepted by models for dual objects (e.g. 3pt-functions):\n"
    " Src       ObjectID for source\n"
    " Snk       ObjectID for sink\n"
    " ESrc      Energy parameter name at source (defaults to "
        << MLU::ModelBase::EnergyPrefix << ")\n"
    " ESnk      Energy parameter name at sink   (defaults to "
        << MLU::ModelBase::EnergyPrefix << ")\n"
    " " << MLU::ModelBase::EDiffPrefix
        << "     Name of derived parameter for energy difference (default: "
        << MLU::ModelBase::EDiffPrefix << ")\n"
    "Parameters accepted by R3 model:\n"
    " C2e       Number of exponents for 2pt (default: Num needed for C3)\n"
    " C2Model   Which model to use for 2pt: Exp (default); Cosh; Sinh\n"
    " Raw       R3 ratio was constructed as a raw ratio - no overlap coefficients\n"
    "Parameters accepted by constant model:\n"
    " const     Constant name to use in output (defaults to '"
                << MLU::ModelBase::ConstantPrefix << "')\n";
  }
  MLU::Grid_exit_handler_disable = true;
  return iReturn;
}
