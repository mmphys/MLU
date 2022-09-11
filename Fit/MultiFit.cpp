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
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"analytic", CL::SwitchType::Flag, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"mindof", CL::SwitchType::Single, "1"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"savecmat", CL::SwitchType::Flag, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      {"v", CL::SwitchType::Single, "0"},
      {"srcsnk", CL::SwitchType::Flag, nullptr},
      {"guess", CL::SwitchType::Single, nullptr},
      // ModelDefaultParams
      {"e", CL::SwitchType::Single, "1"},
      // Other params
      {"fitter", CL::SwitchType::Single, "GSL"},
      {"delta", CL::SwitchType::Single, "3"},
      {"t", CL::SwitchType::Single, nullptr},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"opnames", CL::SwitchType::Flag, nullptr},
      {"covsrc", CL::SwitchType::Single, "Bootstrap"},
      {"covboot", CL::SwitchType::Single, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      const int delta{ cl.SwitchValue<int>("delta") };
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      Common::MakeAncestorDirs( outBaseFileName );
      const int NSamples{ cl.SwitchValue<int>("n") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bOpSort{ !cl.GotSwitch("opnames") };

      if( !cl.GotSwitch( "t") )
        throw std::runtime_error( "No fit ranges specified" );
      Common::FitRanges fitRanges( Common::ArrayFromString( cl.SwitchValue<std::string>( "t" ) ) );

      // Walk the list of parameters on the command-line, loading correlators and making models
      bShowUsage = false;
      std::vector<std::string> OpName;
      std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ ) << "Loading folded correlators\n";
      // Split each argument at the first comma (so the first part can be treated as a filename to glob
      const std::size_t NumArgs{ cl.Args.size() };
      DataSet ds( NSamples );
      std::vector<Model::Args> ModelArgs;
      std::vector<int> ModelFitRange;
      std::vector<int> ModelFitRangeCount( fitRanges.size() );
      for( std::size_t ArgNum = 0; ArgNum < NumArgs; ++ArgNum )
      {
        // First parameter (up to comma) is the filename we're looking for
        std::string FileToGlob{ Common::ExtractToSeparator( cl.Args[ArgNum] ) };
        // Anything after the comma is a list of arguments
        Model::Args vArgs;
        vArgs.FromString( cl.Args[ArgNum], true );
        // Get the fit range (if present)
        const int ThisFitRange{ vArgs.Remove<int>( "Range", 0 ) };
        if( ThisFitRange < 0 || ThisFitRange >= fitRanges.size() )
          throw std::runtime_error( "Fit range " + std::to_string( ThisFitRange ) + " invalid" );
        for( const std::string &sFileName : Common::glob( &FileToGlob, &FileToGlob + 1, inBase.c_str() ) )
        {
          Common::FileNameAtt Att( sFileName );
          const bool bIsCorr{ Common::EqualIgnoreCase( Att.Type, Common::sFold ) };
          if( bIsCorr )
          {
            // This is a correlator - load it
            Att.ParseOpNames( OpName );
            std::string PrintPrefix( 2, ' ' );
            if( !cl.Args[ArgNum].empty() )
            {
              PrintPrefix.append( cl.Args[ArgNum] );
              PrintPrefix.append( 1, ' ' );
            }
            ds.LoadCorrelator( std::move( Att ), Common::COMPAT_DEFAULT, PrintPrefix.c_str() );
            // Now see whether the first parameter is a fit range (i.e. integer index of a defined fit range)
            if( !fitRanges[ThisFitRange].Validate( ds.corr.back().Nt() ) )
            {
              std::stringstream oss;
              oss << "Fit range " << fitRanges[ThisFitRange] << " not valid";
              throw std::runtime_error( oss.str() );
            }
            ModelArgs.emplace_back( vArgs );
            ModelFitRange.push_back( ThisFitRange );
            ++ModelFitRangeCount[ThisFitRange];
          }
          else
          {
            ds.LoadModel( std::move( Att ), cl.Args[ArgNum] );
          }
        }
      }
      if( ds.corr.empty() )
        throw std::runtime_error( "At least one correlator must be loaded to perform a fit" );
      for( int i = 0; i < ModelFitRangeCount.size(); ++i )
        if( ModelFitRangeCount[i] == 0 )
          throw std::runtime_error( "Models don't refer to all fit ranges" );
      ds.SortOpNames( OpName );
      // Describe the number of replicas
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
      outBaseFileName.append( ds.corr[0].Name_.Base );
      outBaseFileName.append( 1, '.' );
      outBaseFileName.append( doCorr ? "corr" : "uncorr" );
      const std::string sSummaryBase{ outBaseFileName + Common::Period + sOpNameConcat };
      outBaseFileName.append( 1, '_' );
      const std::size_t outBaseFileNameLen{ outBaseFileName.size() };
      const Common::SeedType Seed{ ds.corr[0].Name_.Seed };

      // All the models are loaded
      std::unique_ptr<Fitter> m{ MakeFitter( cl, ds, ModelArgs, OpName, std::move( cp ) ) };
      const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
      std::ofstream s;
      for( Common::FitRangesIterator it = fitRanges.begin(); !it.PastEnd(); ++it )
      {
        // Set fit ranges
        int MinExtent = delta;
        std::vector<std::vector<int>> fitTimes( ds.corr.size() );
        for( int i = 0; MinExtent >= delta && i < ds.corr.size(); ++i )
        {
          const Common::FitTime &ft{ it[ModelFitRange[i]] };
          const int Extent{ ft.tf - ft.ti + 1 };
          if( i == 0 || MinExtent > Extent )
            MinExtent = Extent;
          if( Extent >= delta )
          {
            fitTimes[i].resize( Extent );
            for( int j = 0; j < Extent; ++j )
              fitTimes[i][j] = ft.ti + j;
          }
        }
        if( MinExtent >= delta )
        {
          ds.SetFitTimes( fitTimes );
          // Log what file we're processing and when we started
          std::time_t then;
          std::time( &then );
          try
          {
            {
              std::stringstream ss;
              ss << ( doCorr ? "C" : "unc" ) << "orrelated " << m->Type() << " fit on timeslices " << it.to_string( "-", ", " );
              Fitter::SayNumThreads( ss );
              const std::string &sMsg{ ss.str() };
              std::cout << std::string( sMsg.length(), '=' ) << Common::NewLine << sMsg << Common::NewLine;
            }
            double ChiSq;
            int dof;
            outBaseFileName.resize( outBaseFileNameLen );
            outBaseFileName.append( it.to_string( Common::Underscore ) );
            outBaseFileName.append( 1, '.' );
            auto params = m->PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat, Seed );
          }
          catch(const std::exception &e)
          {
            std::cout << "Error: " << e.what() << "\n";
          }
          // Mention that we're finished, what the time is and how long it took
          std::time_t now;
          std::time( &now );
          double dNumSecs = std::difftime( now, then );
          std::string sNow{ std::ctime( &now ) };
          while( sNow.length() && sNow[sNow.length() - 1] == '\n' )
            sNow.resize( sNow.length() - 1 );
          std::stringstream ss;
          ss << sNow << ". Total duration " << std::fixed << std::setprecision(1)
                     << dNumSecs << " seconds.\n";
          std::cout << ss.str();
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
    " <options> Bootstrap1[[,FitRange1],model1[,params1]] [Bootstrap2[[,FitRange1],model2[,params2]] ...]\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, where:\n"
    " FitRange chooses the n'th fit range (default=0)\n"
    " model    is one of {Exp, Cosh, Sinh, ThreePoint, Const}\n"
    " params   Depend on the model:\n"
    " Exp/Cosh/Sinh\n"
    "   Param1 Snk_Src, sink and source names\n"
    "   Param2 'n' to normalise by energy\n"
    "<options> are:\n"
    "--Hotelling Minimum Hotelling Q-value on central replica (default " << DefaultHotelling << ")\n"
    "--sep    Minimum relative separation between energy levels (default " << DefaultEnergySep << ")\n"
    "--delta  Minimum number of timeslices in fit range (default 3)\n"
    "--retry  Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter   Max iteration count, 0 (default) = unlimited\n"
    "--tol    Tolerance of required fits (default 1e-7)\n"
    "--mindof Minimum degrees of freedom (default 1)\n"
    "--fitter (GSL|Minuit2)[,options] fitter (default GSL,Levenberg-Marquardt)\n"
    "         GSL options: lm, lmaccel, dogleg, ddogleg, subspace2D\n"
    "--covsrc source[,options[,...]] build (co)variance from source, i.e. one of\n"
    "         Binned    Use the already binned data\n"
    "         Raw       Use the raw (unbinned) data\n"
    "         Rebin     Rebin the raw data using bin size(s) specified\n"
    "         Bootstrap Use bootstrap replicas (default)\n"
    "         H5,f[,g],d  Load INVERSE covariance from hdf5 file f, group g, dataset d"
    "--covboot How many bootstrap replicas in covariance (-1=no bootstrap)\n"
    "--guess  List of specific values to use for inital guess"
    "-t     Fit range1[,range2[,...]] (start:stop[:numstart=1[:numstop=numstart]])\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--savecmat Save correlation matrix\n"
    "--analytic Analytic derivatives for GSL (default: numeric)\n"
    "--srcsnk   Append _src and _snk to overlap coefficients (ie force different)\n"
    "--opnames  Disable sorting and deduplicating operator name list\n"
    "--help     This message\n";
  }
  return iReturn;
}
