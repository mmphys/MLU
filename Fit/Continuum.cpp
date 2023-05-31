/*************************************************************************************
 
 Chiral continuum fit
 
 Source file: Continuum.cpp

 Copyright (C) 2023
 
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

#include "Continuum.hpp"
#include "Fitter.hpp"
#include "ModelContinuum.hpp"
#include <MLU/DebugInfo.hpp>
#include <chrono>

// Indices for operators in correlator names
const char * pSrcSnk[2] = { "src", "snk" };

struct EnsembleInfo
{
  ModelFilePtr m;
  Model::Args Args;
};

struct EnsembleLess
{
  bool operator()( const EnsembleInfo &lhs, const EnsembleInfo &rhs )
  {
    int i{ Common::CompareIgnoreCase( lhs.m->Ensemble, rhs.m->Ensemble ) };
    if( i )
      return i < 0;
    const Common::Momentum &lp{ lhs.m->Name_.GetFirstNonZeroMomentum().second };
    const Common::Momentum &rp{ rhs.m->Name_.GetFirstNonZeroMomentum().second };
    return lp < rp;
  }
};

void SortModels( DataSet &ds, std::vector<Model::Args> &ModelArgs )
{
  assert( ds.constFile.size() == ModelArgs.size() && "Should be one ModelArgs per ds.constFile" );
  if( ds.constFile.empty() )
    throw std::runtime_error( "Nothing to fit" );
  // Grab the models and corresponding arguments
  std::vector<EnsembleInfo> ei( ds.constFile.size() );
  for( std::size_t i = 0; i < ds.constFile.size(); ++i )
  {
    ei[i].m.reset( ds.constFile[i].release() );
    ei[i].Args = std::move( ModelArgs[i] );
  }
  // Sort
  std::sort( ei.begin(), ei.end(), EnsembleLess() );
  // Put them back
  for( std::size_t i = 0; i < ds.constFile.size(); ++i )
  {
    ds.constFile[i].reset( ei[i].m.release() );
    ModelArgs[i] = std::move( ei[i].Args );
  }
  // Debug
  //for( const ModelFilePtr &m : ds.constFile )
    //std::cout << "  " << m->Name_.Filename << Common::NewLine;
}

int main(int argc, const char *argv[])
{
#ifndef HAVE_MINUIT2
  std::ios_base::sync_with_stdio( false );
#endif
  static const char DefaultFormFactor[] = "f0";
  static const char DefaultEnergySep[] = "0";
  static const char DefaultHotelling[] = "0.05";
  static const char DefaultErrDig[] = "2";
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"f", CL::SwitchType::Single, DefaultFormFactor},
      {"overwrite", CL::SwitchType::Flag, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"chisqdof", CL::SwitchType::Single, "0"},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"mindof", CL::SwitchType::Single, "1"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"summary", CL::SwitchType::Single, "1"},
      {"v", CL::SwitchType::Single, "0"},
      {"strict",  CL::SwitchType::Single, "0"},
      {"maxE",  CL::SwitchType::Single, "10"},
      {"errdig", CL::SwitchType::Single, DefaultErrDig},
      // ModelDefaultParams
      {"e", CL::SwitchType::Single, "1"},
      {"N", CL::SwitchType::Single, "0"},
      // Covariance parameters
      {"covsrc", CL::SwitchType::Single, nullptr},
      {"covboot", CL::SwitchType::Single, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      // Other params
      {"fitter", CL::SwitchType::Single, "GSL"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"debug-signals", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      if( cl.GotSwitch( "debug-signals" ) )
        Common::Grid_debug_handler_init();
      const std::string sFFValue{ cl.SwitchValue<std::string>("f") };
      [[maybe_unused]]
      const Common::FormFactor ff{ Common::FromString<Common::FormFactor>( sFFValue ) };
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      Common::MakeAncestorDirs( outBaseFileName );
      const int NSamples{ cl.SwitchValue<int>("n") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };

      // Walk the list of parameters on the command-line, loading correlators and making models
      bShowUsage = false;
      std::vector<std::string> OpName;
      std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ ) << "Loading folded correlators\n";
      // Split each argument at the first comma (so the first part can be treated as a filename to glob
      const std::size_t NumArgs{ cl.Args.size() };
      DataSet ds( NSamples );
      std::vector<Model::Args> ModelArgs;
      for( std::size_t ArgNum = 0; ArgNum < NumArgs; ++ArgNum )
      {
        // First parameter (up to comma) is the filename we're looking for
        std::string FileToGlob{ Common::ExtractToSeparator( cl.Args[ArgNum] ) };
        std::vector<std::string> Filenames{ Common::glob( &FileToGlob, &FileToGlob + 1,
                                                          inBase.c_str() ) };
        bool bGlobEmpty{ true };
        // Anything after the comma is a list of arguments
        Model::Args vArgs;
        vArgs.FromString( cl.Args[ArgNum], true );
        std::string Ensemble{ vArgs.Remove( Common::sEnsemble ) };
        if( vArgs.find( sFF ) == vArgs.end() )
          vArgs.emplace( sFF, sFFValue );
        for( const std::string &sFileName : Filenames )
        {
          bGlobEmpty = false;
          // This is a correlator - load it
          std::string PrintPrefix( 2, ' ' );
          if( !cl.Args[ArgNum].empty() )
          {
            PrintPrefix.append( cl.Args[ArgNum] );
            PrintPrefix.append( 1, ' ' );
          }
          ds.LoadModel( Common::FileNameAtt( sFileName, &OpName ), PrintPrefix.c_str(),
                        Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT
                        | Common::COMPAT_DISABLE_ENSEMBLE );
          ModelArgs.emplace_back( vArgs );
          if( !Ensemble.empty() )
            ds.constFile.back()->Ensemble = Ensemble;
        }
        if( bGlobEmpty )
          throw std::runtime_error( "No files matched " + FileToGlob );
      }
      // Sort the models
      SortModels( ds, ModelArgs );
      // Describe the number of replicas
      ds.SortOpNames( OpName );
      std::cout << "Using ";
      if( ds.NSamples == ds.MaxSamples )
        std::cout << "all ";
      else
        std::cout << "first " << ds.NSamples << " of ";
      std::cout << ds.MaxSamples << Common::Space << ds.front().getData().SeedTypeString()
                << " replicas";
      if( ds.NSamples != ds.MaxSamples )
        std::cout << " (all " << ds.MaxSamples << " for var/covar)";
      std::cout << Common::NewLine;

      // Work out how covariance matrix should be constructed and tell user
      CovarParams cp{ cl, ds };
      std::cout << cp << Common::NewLine;

      // I'll need a sorted, concatenated list of the operators in the fit for filenames
      std::string sOpNameConcat;
      sOpNameConcat = OpName[0];
      for( std::size_t i = 1; i < OpName.size(); i++ )
      {
        sOpNameConcat.append( 1, '_' );
        sOpNameConcat.append( OpName[i] );
      }

      // Now make base filenames for output
      // If the output base is given, but doesn't end in '/', then it already contains the correct name
      if( outBaseFileName.empty() || outBaseFileName.back() == '/' )
      {
        // Simplistic, but better is hard to automate
        const Common::FileNameAtt &f{ ds.front().Name_ };
        if( f.BaseShortParts.size() < 3 || f.Spectator.empty() )
          throw std::runtime_error( "Please include output filename in -o" );
        outBaseFileName.append( f.BaseShortParts[0] );
        outBaseFileName.append( 1, '_' );
        outBaseFileName.append( Common::MesonName( f.BaseShortParts[1], f.Spectator ) );
        outBaseFileName.append( 1, '_' );
        outBaseFileName.append( Common::MesonName( f.BaseShortParts[2], f.Spectator ) );
        for( std::size_t i = 3; i < f.BaseShortParts.size(); ++i )
        {
          outBaseFileName.append( 1, '_' );
          outBaseFileName.append( f.BaseShortParts[i] );
        }
      }
      outBaseFileName.append( 1, '.' );
      outBaseFileName.append( doCorr ? "corr" : "uncorr" );

      // All the models are loaded
      std::unique_ptr<Fitter> f{ Fitter::Make( cl, ds, std::move( ModelArgs ), OpName,
                                               std::move( cp ) ) };
      bool bAllParamsResolved{ true };
      const auto Start{ std::chrono::steady_clock::now() };
      // Log what file we're processing and when we started
      const auto then{ std::chrono::steady_clock::now() };
      // Set fit ranges
      std::vector<std::vector<int>> fitTimes( ds.corr.size() );
      for( int i = 0; i < ds.corr.size(); ++i )
        fitTimes[i] = { 0 };
      ds.SetFitTimes( fitTimes );
      bool bFitOK{ true };
      try
      {
        {
          std::stringstream ss;
          if( f->bTestRun )
            ss << "Test run of ";
          ss << (doCorr ? "C" : "Unc") << "orrelated " << f->Type() << " chiral continuum fit";
          Fitter::SayNumThreads( ss );
          const std::string &sMsg{ ss.str() };
          if( !f->bTestRun )
            std::cout << std::string( sMsg.length(), '=' ) << Common::NewLine;
          std::cout << sMsg << Common::NewLine;
        }
        double ChiSq;
        int dof;
        bAllParamsResolved =
        f->PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat );
      }
      catch(const std::exception &e)
      {
        bFitOK = false;
        std::cout << "Error: " << e.what() << "\n";
      }
      // Mention that we're finished, what the time is and how long it took
      if( !f->bTestRun )
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
        std::cout << ". " << sNow << Common::NewLine;
      }
      const auto now{ std::chrono::steady_clock::now() };
      auto S{ std::chrono::duration_cast<std::chrono::seconds>( now - Start ).count() };
      std::cout << std::string( 50, '=' ) << "\nFit ";
      if( bFitOK )
        std::cout << "succeeded";
      else
        std::cout << "failed";
      std::cout << ". Total duration ";
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
      // Error if only one fit requested and it failed ... or all params not resolved
      if( !bFitOK )
        iReturn = 2;
      else if( !bAllParamsResolved )
        iReturn = 3;
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
    " [Options] Model[,params]...\n"
    "Perform a chiral continuum fit of the per Ensemble data\n"
    "Options:\n"
    "-f       Form factor: f0, fplus, fpar or fperp (default: " << DefaultFormFactor << ")\n"
    "--Hotelling Minimum Hotelling Q-value on central replica (default " << DefaultHotelling << ")\n"
    "--chisqdof  Maximum chi^2 / dof on central replica\n"
    "--sep    Minimum relative separation between energies (default " << DefaultEnergySep << ")\n"
    "--retry  Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter   Max iteration count, 0 (default) = unlimited\n"
    "--tol    Tolerance of required fits (default 1e-7)\n"
    "--mindof Minimum degrees of freedom (default 1)\n"
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
    "--summary 0 no summaries; 1 model_td.seq.txt only; 2 model_td and model.seq.txt\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-N     Use lattice dispersion relation for boosted energies with N=L/a\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--debug-signals Trap signals (code courtesy of Grid)\n"
    "--help     This message\n";
  }
  Common::Grid_exit_handler_disable = true;
  return iReturn;
}
