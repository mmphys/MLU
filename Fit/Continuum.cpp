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
#include "ModelContinuum.hpp"
#include <MLU/DebugInfo.hpp>
#include <chrono>

// Indices for operators in correlator names
const char * pSrcSnk[2] = { "src", "snk" };

struct ModelWithArgs
{
  ModelFilePtr m;
  Model::Args Args;
};

struct ModelWithArgsLess
{
  bool operator()( const ModelWithArgs &lhs, const ModelWithArgs &rhs )
  {
    int i{ Common::CompareIgnoreCase( lhs.m->Ensemble, rhs.m->Ensemble ) };
    if( i )
      return i < 0;
    const Common::Momentum &lp{ lhs.m->Name_.GetFirstNonZeroMomentum().second };
    const Common::Momentum &rp{ rhs.m->Name_.GetFirstNonZeroMomentum().second };
    return lp < rp;
  }
};

CreateParams::CreateParams( const std::vector<std::string> &O, const Common::CommandLine &C )
: Model::CreateParams( O, C ),
  EnsembleMap{
    {"C1",  {24, 64}},
    {"C2",  {24, 64}},
    {"F1M", {48, 96}},
    {"M1",  {32, 64}},
    {"M2",  {32, 64}},
    {"M3",  {32, 64}} }
{
}

ContinuumFit::ContinuumFit( Common::CommandLine &cl_ )
: cl{cl_},
  NumSamples{cl.SwitchValue<int>("n")},
  doCorr{ !cl.GotSwitch( "uncorr" ) },
  CovarBlock{ !cl.GotSwitch( "noblock" ) },
  sFFValue{ cl.SwitchValue<std::string>("f") },
  ff{ Common::FromString<Common::FormFactor>( sFFValue ) },
  inBase{ cl.SwitchValue<std::string>("i") },
  outBaseFileName{ cl.SwitchValue<std::string>("o") },
  ds{NumSamples}
{
  Common::MakeAncestorDirs( outBaseFileName );
}

void ContinuumFit::LoadModels()
{
  std::cout << std::setprecision( std::numeric_limits<double>::max_digits10+2 )
            << "Performing chiral continuum fit from ensemble results" << std::endl;
  // Walk the list of parameters on the command-line, loading correlators and making models
  for( std::size_t ArgNum = 0; ArgNum < cl.Args.size(); ++ArgNum )
  {
    // First parameter is the filename we're looking for, optionally followed by comma and arguments
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
      // This is a correlator - load it
      std::string PrintPrefix( 2, ' ' );
      if( !cl.Args[ArgNum].empty() )
      {
        PrintPrefix.append( cl.Args[ArgNum] );
        PrintPrefix.append( 1, ' ' );
      }
      Common::FileNameAtt fna{ sFileName, &OpName };
      if( ( ff == Common::FormFactor::fplus || ff == Common::FormFactor::fperp )
         && !fna.HasNonZeroMomentum() )
        std::cout << PrintPrefix << "Ignoring zero momentum file " << sFileName << Common::NewLine;
      else
      {
        bGlobEmpty = false;
        ds.LoadModel( std::move( fna ), PrintPrefix.c_str(),
                     Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT
                     | Common::COMPAT_DISABLE_ENSEMBLE );
        ModelArgs.emplace_back( vArgs );
        if( !Ensemble.empty() )
          ds.constFile.back()->Ensemble = Ensemble;
      }
    }
    if( bGlobEmpty )
      throw std::runtime_error( "No files matched " + FileToGlob );
  }
}

void ContinuumFit::SortModels()
{
  assert( ds.constFile.size() == ModelArgs.size() && "Should be one ModelArgs per ds.constFile" );
  if( ds.constFile.empty() )
    throw std::runtime_error( "Nothing to fit" );
  // Sort the models and corresponding arguments by ensemble
  std::vector<ModelWithArgs> ei( ds.constFile.size() );
  for( std::size_t i = 0; i < ds.constFile.size(); ++i )
  {
    ei[i].m.reset( ds.constFile[i].release() );
    ei[i].Args = std::move( ModelArgs[i] );
  }
  // Sort
  std::sort( ei.begin(), ei.end(), ModelWithArgsLess() );
  // Put them back
  for( std::size_t i = 0; i < ds.constFile.size(); ++i )
  {
    ds.constFile[i].reset( ei[i].m.release() );
    ModelArgs[i] = std::move( ei[i].Args );
  }

  // Describe the number of replicas
  std::cout << "Using ";
  if( ds.NSamples == ds.MaxSamples )
    std::cout << "all ";
  else
    std::cout << "first " << ds.NSamples << " of ";
  std::cout << ds.MaxSamples << Common::Space << ds.front().getData().SeedTypeString()
            << " replicas";
  if( ds.NSamples != ds.MaxSamples )
    std::cout << " (all " << ds.MaxSamples << " for var/covar)";
  std::cout << std::endl;
  // Debug
  //for( const ModelFilePtr &m : ds.constFile )
    //std::cout << "  " << m->Name_.Filename << std::endl;
}

void ContinuumFit::LoadExtra()
{
  const std::vector<std::string> &vExtra{ cl.SwitchStrings( "model" ) };
  const std::string PrintPrefix( 2, ' ' );
  const unsigned int CompareFlags{ Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT
                                 | Common::COMPAT_DISABLE_ENSEMBLE };
  for( const std::string &sFileName : vExtra )
  {
    Common::FileNameAtt fna;
    fna.Filename = sFileName;
    fna.Type = "Model constants";
    if( ds.empty() )
    {
      fna.Ext = "h5";
      fna.bSeedNum = true;
      fna.Seed = 0;
    }
    else
    {
      const Common::FileNameAtt &S{ ds.front().Name_ };
      fna.Ext = S.Ext;
      fna.bSeedNum = S.bSeedNum;
      fna.Seed = S.Seed;
    }
    ds.LoadModel( std::move( fna ), "", PrintPrefix.c_str(), CompareFlags );
  }
}

void ContinuumFit::MakeOutputFilename()
{
  // Sort operator names (also renumbers references to these names)
  ds.SortOpNames( OpName );
  // Make sorted, concatenated list of operators in the fit for filenames
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
    const Common::FileNameAtt &fna{ ds.front().Name_ };
    if( fna.BaseShortParts.size() < 3 || fna.Spectator.empty() )
      throw std::runtime_error( "Please include output filename in -o" );
    outBaseFileName.append( fna.BaseShortParts[0] );
    outBaseFileName.append( 1, '_' );
    outBaseFileName.append( Common::MesonName( fna.BaseShortParts[1], fna.Spectator ) );
    outBaseFileName.append( 1, '_' );
    outBaseFileName.append( Common::MesonName( fna.BaseShortParts[2], fna.Spectator ) );
    for( std::size_t i = 3; i < fna.BaseShortParts.size(); ++i )
    {
      outBaseFileName.append( 1, '_' );
      outBaseFileName.append( fna.BaseShortParts[i] );
    }
  }
  outBaseFileName.append( 1, '.' );
  outBaseFileName.append( doCorr ? "corr" : "uncorr" );
  outBaseFileName.append( 1, '_' );
  outBaseFileName.append( sFFValue );
  outBaseFileName.append( 1, '.' );
}

void ContinuumFit::MakeCovarBlock()
{
  for( std::size_t i = 0; i < ds.mCovar.size1; ++i )
    for( std::size_t j = 0; j < i; ++j )
      if( i != j && !Common::EqualIgnoreCase( ds.constFile[i]->Ensemble, ds.constFile[j]->Ensemble ) )
      {
        ds.mCovar(i,j) = 0;
        ds.mCovar(j,i) = 0;
      }
}

int ContinuumFit::DoFit()
{
  int iReturn{ EXIT_SUCCESS };

  // Set fit ranges
  std::vector<std::vector<int>> fitTimes( ModelArgs.size() );
  for( int i = 0; i < f->ModelArgs.size(); ++i )
    fitTimes[i] = { static_cast<int>( f->model[i]->GetFitColumn() ) };
  ds.SetFitTimes( fitTimes, false );
  if( CovarBlock )
    MakeCovarBlock();

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
        std::cout << std::string( sMsg.length(), '=' ) << std::endl;
      std::cout << sMsg << std::endl;
    }
    double ChiSq;
    int dof;
    if( !f->PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat ) )
      iReturn = 3; // Not all parameters resolved
  }
  catch(const std::exception &e)
  {
    iReturn = 2;
    std::cout << "Fitter error: " << e.what() << std::endl;
  }
  return iReturn;
}

int ContinuumFit::Run()
{
  LoadModels();
  SortModels();
  MakeOutputFilename();
  LoadExtra();
  // Work out how covariance matrix should be constructed and tell user
  // TODO: No choice - we have to make it from the resampled data, block diagonal per Ensemble
  CovarParams cp{ cl, ds };
  std::cout << cp << std::endl;
  std::cout << "Covariance matrix from resampled data, block diagonal per Ensemble" << std::endl;

  // Make the fitter
  f.reset( Fitter::Make( CreateParams{OpName,cl}, ds, std::move(ModelArgs), std::move(cp), false ) );
  // Now do the fit
  return DoFit();
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
      {"model", CL::SwitchType::Multiple, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"chisqdof", CL::SwitchType::Single, "0"},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"mindof", CL::SwitchType::Single, "1"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"noblock", CL::SwitchType::Flag, nullptr},
      {"summary", CL::SwitchType::Single, "1"},
      {"nopolap", CL::SwitchType::Single, ""},
      {"v", CL::SwitchType::Single, "0"},
      {"strict",  CL::SwitchType::Single, "0"},
      {"maxE",  CL::SwitchType::Single, "10"},
      {"errdig", CL::SwitchType::Single, DefaultErrDig},
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
      bShowUsage = false;
      if( cl.GotSwitch( "debug-signals" ) )
        Common::Grid_debug_handler_init();
      const auto Start{ std::chrono::steady_clock::now() };
      ContinuumFit cf{ cl };
      iReturn = cf.Run();
      const auto Stop{ std::chrono::steady_clock::now() };
      // Say whether this worked, what the time is and how long it took
      switch( iReturn )
      {
        case 0:
          std::cout << "OK";
          break;
        case 3:
          std::cout << "Warning: not all parameters resolved";
          break;
        default:
          std::cout << "Error " << iReturn;
          break;
      }
      const auto mS{ std::chrono::duration_cast<std::chrono::milliseconds>( Stop - Start ) };
      std::cout << ". Fit duration ";
      if( mS.count() < 1100 )
        std::cout << mS.count() << " milliseconds";
      else
      {
        const auto S{ std::chrono::duration_cast<std::chrono::duration<double>>( Stop - Start ) };
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << S.count() << " seconds";
        std::cout << ss.str();
      }
      std::time_t ttNow;
      std::time( &ttNow );
      std::string sNow{ std::ctime( &ttNow ) };
      while( sNow.length() && sNow.back() == '\n' )
        sNow.resize( sNow.length() - 1 );
      std::cout << ". " << sNow << std::endl;
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
    "--model  An additional model file, e.g. to load lattice spacings\n"
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
    "--noblock Use full covariance matrix (rather than block diagonal per ensemble)\n"
    "--summary 0 no summaries; 1 model_td.seq.txt only; 2 model_td and model.seq.txt\n"
    "--nopolap List of overlap coefficients which don't depend on momenta\n"
    "          TODO: This is ignored for now. Should be moved into MultiFit only\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--debug-signals Trap signals (code courtesy of Grid)\n"
    "--help     This message\n";
  }
  Common::Grid_exit_handler_disable = true;
  return iReturn;
}
