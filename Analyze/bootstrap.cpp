/*
 * bootstrap.cpp - based on ...
 * resample.cpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2016 Antonin Portelli
 *   Modified           2019 Michael Marshall
 *
 * LatAnalyze 3 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatAnalyze 3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatAnalyze 3.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "CommonGrid.hpp"
#include "CommonLatAn.hpp"

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Io/Hdf5File.hpp>
#include <H5CompType.h> // HDF5 Library

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif

// Show me the averages for each timeslice
void ShowTimeSliceAvg(const Latan::Dataset<Latan::DMat> &data) {
  const int nFile{static_cast<int>(data.size())};
  if( !nFile )
    std::cout << "No timeslices to average over!" << std::endl;
  else {
    const int nt{static_cast<int>(data[0].rows())};
    std::cout << "Averages of " << nFile << " files over " << nt << " timeslices:" << std::endl;
    std::vector<std::complex<double>> Avg(nt);
    for (unsigned int t = 0; t < nt; ++t) {
      Avg[t] = 0;
      for (unsigned int j = 0; j < nFile; ++j) {
        std::complex<double> d( data[j](t,0), data[j](t,1) );
        //std::cout << "d[" << j << "]=" << d << std::endl;
        Avg[t] += d;
        //std::cout << "a[" << j << "]=" << a << std::endl;
      }
      Avg[t] /= nFile;
      std::cout << "C(" << t << ")=" << Avg[t] << std::endl;
    }
  }
}

struct BootstrapParams
{
  std::string outStem;
  std::string OutFileExt;
  Latan::SeedType seed;
  Latan::Index nSample;
  Latan::Index binSize;
  bool bShowAverages;
  bool bSaveSummaries;
  bool bSaveBootstrap;
  bool dumpBoot;
};

bool PerformBootstrap( const BootstrapParams &par, Latan::Dataset<Latan::DMat> &in, const std::string &Contraction )//, int momentum_squared )
{
  const bool bPerformBootstrap{ par.bSaveSummaries || par.bSaveBootstrap };
  //if( !bPerformBootstrap && !par.dumpBoot )
    //throw std::runtime_error( "No bootstrap to perform" );
  // Skip bootstrap if output exists
  const std::string sOutFileBase{ Common::tokenReplaceCopy( par.outStem, "corr", Contraction ) };
  std::string sOutFileName{ Common::MakeFilename( sOutFileBase, Common::sBootstrap, par.seed, par.OutFileExt ) };
  if( Common::FileExists( sOutFileName ) || (!par.bSaveBootstrap && Common::FileExists(
      Common::MakeFilename(sOutFileBase, Common::SummaryNames[0], par.seed, TEXT_EXT))))
  {
    std::cout << "Skipping " << Contraction << " because output already exists" << std::endl;
    return false;
  }
  // Make sure there are at least 5x more bootstrap samples than data points
  std::cout << Contraction << std::endl;
  if( par.bShowAverages )
    ShowTimeSliceAvg( in );
  if( bPerformBootstrap ) {
    if( ( par.nSample / 5 ) < in.size() )
      throw std::runtime_error( "nSample=" + std::to_string( par.nSample ) + " too small for dataset with " + std::to_string( in.size() ) + " elements" );
    std::cout << "-- resampling (" << par.nSample << " samples)..." << std::endl;
    if( par.binSize != 1 )
      in.bin( par.binSize );
    Latan::DMatSample out = in.bootstrapMean( par.nSample, par.seed );
    if( par.bSaveBootstrap ) {
      std::cout << "Saving sample to " << sOutFileName << std::endl;
      Latan::Io::save<Latan::DMatSample>(out, sOutFileName, Latan::File::Mode::write, Common::sBootstrap);
    }
    if( par.bSaveSummaries )
      Common::SummariseBootstrapCorr( out, sOutFileBase, par.seed );//, momentum_squared );
  }
  if( par.dumpBoot ) { // Dump the bootstrap sequences
    std::string SeqFileName{ Common::MakeFilename( sOutFileBase, "bootseq", par.seed, TEXT_EXT ) };
    std::cout << "Saving bootstrap sequence to '" << SeqFileName << "' ...";
    std::ofstream file(SeqFileName);
    file << "# bootstrap sequences" << std::endl;
    file << "#      bin size: " << par.binSize << std::endl;
    file << "#          seed: " << std::to_string( par.seed ) << std::endl;
    in.dumpBootstrapSeq(file, par.nSample, par.seed);
    std::cout << std::endl;
  }
  return true;
}

/*struct Momentum {
  int x;
  int y;
  int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
};

struct SubStrings {
  std::string              Name;
  std::vector<std::string> Strings;
};

const std::string & MesonSuffix( const std::string MesonTypeName )
{
  static const std::string axial{ "ax" };
  static const std::string axial_name{ "_ax" };
  static const std::string other_name{};
  return axial.compare( MesonTypeName ) == 0 ? axial_name : other_name;
}*/

// Make a default manifest
// Heavy-light meson decays. 2pt function with current inserted and momentum on light meson

void Z2Bootstrap( const BootstrapParams &bsParams, const Common::Momentum &mom, const std::string &Heavy, const std::string &Light )
{
  static constexpr unsigned int Nt{64};
  static constexpr unsigned int CStart{3000};
  static constexpr unsigned int CSkip{40};
  static constexpr unsigned int CCount{10};
  static constexpr unsigned int CEnd{CStart + CSkip * CCount};
  static const std::string sPath{
#ifdef DEBUG
    "/Users/s1786208"               // Laptop
#else
    "/home/dp008/dp008/dc-mars3"  // Tesseract
#endif
    "/data/201910Plat/mesons/C1/Z2/" };
  static const std::string Dot{ "." };    // used inside filenames
  static const std::string H5{ Dot + "h5" };    // used inside filenames
  static const std::string Sep{ "_" };    // used inside filenames
  static const std::string Space{ " " };  // whitespace as field separator / human readable info
  static const std::string Sink{ Sep + "sink" };    // used inside filenames
  static const std::string sProp{ "prop" + Sep };    // used inside filenames
  static const std::string sMom0{ Sep + "p" + Sep + Common::Momentum(0,0,0).to_string( Sep ) };
  const std::string sMom{ Sep + "p" + Sep + mom.to_string( Sep )};
  const std::string sMomNeg{ Sep + "p" + Sep + mom.to_string( Sep, true )};
  using Algebra = Grid::Gamma::Algebra;
  static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
  static const std::vector<std::string> algNames = { "g5", "gT5" };
  assert( alg.size() == algNames.size() && "Abbreviations should match gamma algebra" );
  static const int NumAlgebra{ static_cast<int>( alg.size() ) };
  static const int NumCorr{ NumAlgebra * NumAlgebra };
  const std::string CorrPrefix{ Heavy + Sep + Light + sMom };
  std::vector<std::string> CorrSuffixes( NumCorr );
  for( int iSrc = 0; iSrc < NumAlgebra; iSrc++ )
    for( int iSnk = 0; iSnk < NumAlgebra; iSnk++ ) {
      const int iCorr{ iSnk * NumAlgebra + iSrc };
      std::stringstream ss;
      ss << Sep << algNames[iSnk] << Sep << algNames[iSrc];
      CorrSuffixes[iCorr] = ss.str();
    }
  const unsigned int NumEnds{ Common::EqualIgnoreCase( Heavy, Light ) ? 1u : 2u };
  const std::size_t NumSamplesT{ NumEnds * CCount };
  const std::size_t NumSamples{ NumSamplesT * Nt };
  std::vector<Latan::Dataset<Latan::DMat>> bsData( NumCorr );
  std::vector<Latan::Dataset<Latan::DMat>> bsDataT( NumCorr );
  std::vector<Common::Correlator> buffer( NumCorr );
  for( int iCorr = 0 ; iCorr < NumCorr; iCorr++ ) {
    bsData[ iCorr ].resize( NumSamples );
    bsDataT[ iCorr ].resize( NumSamplesT );
    buffer[iCorr].resize( Nt );
  }
  int CorrIndex{ 0 };
  BootstrapParams par{ bsParams };
  par.bSaveSummaries = true;
  par.bSaveBootstrap = false;
  for( int t = 0; t < Nt; ++t ) {
    const std::string tPrefix{ Sep + "t" + std::to_string( t ) };
    int CorrIndexT{ 0 };
    for( unsigned int iEnd = 0; iEnd < NumEnds; iEnd++ ) {
      const std::string & Left{ iEnd == 0 ? Heavy : Light };
      const std::string & Right{ iEnd == 0 ? Light : Heavy };
      for( int iConfig = CStart; iConfig < CEnd; iConfig += CSkip ) {
        const std::string sFileName{ sPath + sProp + Left + tPrefix + sMom + Sep + sProp + Right + tPrefix + sMom0 + Sink + sMomNeg + Dot + std::to_string( iConfig ) + H5};
        Common::ReadComplexArray(buffer, alg, sFileName, 0);
        for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
        {
          Common::CopyCorrelator( bsData[iCorr][CorrIndex], buffer[iCorr], t );
          Common::CopyCorrelator( bsDataT[iCorr][CorrIndexT], buffer[iCorr] );
        }
        CorrIndex++;
        CorrIndexT++;
      }
    }
    // Now perform a bootstrap of this timeslice
    for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
      PerformBootstrap( par, bsDataT[iCorr], CorrPrefix + tPrefix + CorrSuffixes[iCorr] );
  }
  // Now perform a bootstrap overall
  par.bSaveBootstrap = true;
  for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
    PerformBootstrap( par, bsData[iCorr], CorrPrefix + CorrSuffixes[iCorr] );
}

/*****************************************************************

 Perform a bootstrap of all the files specified on the command line.
 Organise all the files by contraction, and sort them by trajectory.
 Ensure trajectories are processed in the same order so the
 bootstrap replicas for each contraction can be combined later during analysis

*****************************************************************/

static void Usage( const char * pExecutableName, Latan::OptParser &opt )
{
  std::string cmd{ pExecutableName };
  auto pos = cmd.find_last_of('/');
  if( pos != std::string::npos && pos < cmd.length() - 1 )
    cmd = cmd.substr( pos + 1 );
  std::cerr << "usage: " << pExecutableName << " <options> ContractionFile1 [ContractionFile2 ...]"
  << "\nPossible options:\n" << opt << std::endl;
}

/*int DebugTest1()
{
  static const std::string Filename{ "prop_h1_t0_p_0_0_0_prop_h1_t0_p_0_0_0_sink_p_0_0_0.3000.h5" };
  MesonFile f{ Filename, 0, false, 0, Grid::Gamma::Algebra::GammaZ,
    Grid::Gamma::Algebra::GammaT };
  std::vector<std::complex<double>> buf( 64 );
  f.GetCorellator( buf );
  return 1;
}*/

int main(int argc, char *argv[])
{
  /*int i = DebugTest1();
  if( i )
    return i;*/
  // parse command line/options ////////////////////////////////////////////////////////
  Latan::OptParser       opt;
  opt.addOption("n", "nsample"   , Latan::OptParser::OptType::value,   true,
                "number of samples", DEF_NSAMPLE);
  opt.addOption("b", "bin"       , Latan::OptParser::OptType::value,   true,
                "bin size", "1");
  opt.addOption("r", "seed"      , Latan::OptParser::OptType::value,   true,
                "random generator seed (default: random)");
  opt.addOption("o", "output"    , Latan::OptParser::OptType::value,   true,
                "output file", "@corr@");
  opt.addOption("f", "format"    , Latan::OptParser::OptType::value,   true,
                "output file format", DEF_FMT);
  opt.addOption("d", "dump-boot" , Latan::OptParser::OptType::trigger, true,
                "dump bootstrap sequence");
  opt.addOption("t", "txt"       , Latan::OptParser::OptType::trigger, true,
                "don't save text versions of correlators (for GnuPlot)");
  opt.addOption("a", "average"   , Latan::OptParser::OptType::trigger, true,
                "show correlator averages");
  opt.addOption("x", "exclude"   , Latan::OptParser::OptType::value,   true,
                "exclude files (colon separated list)", "");
  opt.addOption("" , "help"      , Latan::OptParser::OptType::trigger, true,
                "show this help message and exit");
  bool parsed = opt.parse(argc, argv) and !opt.gotOption("help") /*and opt.getArgs().size() > 0*/;
  // If the output stem was specified, ensure it includes "@corr@"
  BootstrapParams par;
  par.bSaveBootstrap = true;
  if( parsed )
  {
    par.outStem = opt.optionValue("o");
    if( par.outStem.find( "@corr@" ) == std::string::npos )
    {
      parsed=false;
      std::cout << "Output stem " << par.outStem << " invalid" << std::endl;
    }
  }
  if (!parsed) {
    Usage( argv[0], opt );
    return EXIT_FAILURE;
  }

  par.nSample = opt.optionValue<Latan::Index>("n");
  par.binSize = opt.optionValue<Latan::Index>("b");
  if (opt.gotOption("r"))
    par.seed = opt.optionValue<Latan::SeedType>("r");
  else {
    std::random_device rd;
    par.seed = rd();
  }
  par.OutFileExt = opt.optionValue("f");
  par.dumpBoot = opt.gotOption("d");
  par.bShowAverages = opt.gotOption("a");
  par.bSaveSummaries = !opt.gotOption("t");
  std::cout << "Creating bootstrap output " << par.outStem << ", seed=" << par.seed << std::endl;

  int iReturn = EXIT_SUCCESS;
  try
  {
    ///////////////////////////
    // If there are files specified on the command line,
    // parse the input file names, grouping by correlator, indexed by trajectory.
    ///////////////////////////
    if( opt.getArgs().size() )
    {
      Common::Manifest Manifest{ opt.getArgs(), opt.optionValue("x") };
      // Walk the list of contractions, performing a separate bootstrap for each
      int BootstrapCount = 0;
      Latan::DMatSample out(par.nSample);
      for( auto itc = Manifest.begin(); itc != Manifest.end(); itc++ )
      {
        const std::string &Contraction{itc->first};
        const Common::TrajList &l{itc->second};
        const unsigned int nFile{ static_cast<unsigned int>(l.TrajFile.size())};
        std::cout << "Loading " << nFile << " files for " << Contraction << std::endl;
        unsigned int j = 0; // which trajectory are we processing
        Common::Correlator buf;
        Latan::Dataset<Latan::DMat> data(nFile);
        for( auto it = l.TrajFile.begin(); it != l.TrajFile.end(); it++, j++ )
        {
          const int &traj{it->first}; // Trajectory number
          const Common::TrajFile &tf{*it->second};
          const std::string &Filename{tf.Filename};
          std::cout << "\t" << traj << "\t" << Filename << std::endl;
          Common::ReadComplexArray( buf, Filename );
          Common::CopyCorrelator( data[j], buf );
          PerformBootstrap( par, data, Contraction );
          BootstrapCount++;
        }
      }
      std::cout << "Bootstraps written for " << BootstrapCount
      << " / " << Manifest.size() << " correlators" << std::endl;
    } else {
      ///////////////////////////
      // Nothing specified on the command-line
      // Perform default bootstraps
      // Heavy-light semi-leptonics
      ///////////////////////////
      const std::string qHeavy{ "h1" };
      const std::string qLight{ "l" };
      std::cout << "No files specified. Making default manifest for "
      << qHeavy << " and " << qLight << std::endl;
      static const Common::Momentum Z2Momenta[] = {
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 1, 1, 1 },
        { 2, 0, 0 },
      };
      for( const Common::Momentum &m : Z2Momenta )
      {
        Z2Bootstrap( par, m, qHeavy, qLight );
        Z2Bootstrap( par, m, qHeavy, qHeavy );
        Z2Bootstrap( par, m, qLight, qLight );
      }
    }
  }
  catch(const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
    } catch( ... ) {
      std::cerr << "Error: Unknown exception" << std::endl;
    }
  return iReturn;
}
