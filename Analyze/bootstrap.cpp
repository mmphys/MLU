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

using Algebra = Grid::Gamma::Algebra;
static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
static const std::vector<std::string> algNames = { "g5", "gT5" };
static const int NumAlgebra{ static_cast<int>( alg.size() ) };
static const int NumCorr{ NumAlgebra * NumAlgebra };

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
  Common::SeedType seed;
  int nSample;
  int binSize;
  bool bSaveSummaries;
  bool bSaveBootstrap;
  bool PerformBootstrap( Latan::Dataset<Latan::DMat> &in, const std::string &Contraction ) const;
};

bool BootstrapParams::PerformBootstrap( Latan::Dataset<Latan::DMat> &in, const std::string &Contraction ) const//, int momentum_squared )
{
  const bool bPerformBootstrap{ bSaveSummaries || bSaveBootstrap };
  //if( !bPerformBootstrap && !dumpBoot )
    //throw std::runtime_error( "No bootstrap to perform" );
  // Skip bootstrap if output exists
  const std::string sOutFileBase{ Common::tokenReplaceCopy( outStem, "corr", Contraction ) };
  std::string sOutFileName{ Common::MakeFilename( sOutFileBase, Common::sBootstrap, seed, "h5" ) };
  if( Common::FileExists( sOutFileName ) || (!bSaveBootstrap && Common::FileExists(
      Common::MakeFilename(sOutFileBase, Common::SummaryNames[0], seed, TEXT_EXT))))
  {
    std::cout << "Skipping " << Contraction << " because output already exists" << std::endl;
    return false;
  }
  // Make sure there are at least 5x more bootstrap samples than data points
  std::cout << Contraction << std::endl;
  //if( bShowAverages )
    //ShowTimeSliceAvg( in );
  if( bPerformBootstrap ) {
    if( ( nSample / 5 ) < in.size() )
      throw std::runtime_error( "nSample=" + std::to_string( nSample ) + " too small for dataset with " + std::to_string( in.size() ) + " elements" );
    std::cout << "-- resampling (" << nSample << " samples)..." << std::endl;
    if( binSize != 1 )
      in.bin( binSize );
    Latan::DMatSample out = in.bootstrapMean( nSample, seed );
    if( bSaveBootstrap ) {
      std::cout << "Saving sample to " << sOutFileName << std::endl;
      Latan::Io::save<Latan::DMatSample>(out, sOutFileName, Latan::File::Mode::write, Common::sBootstrap);
    }
    if( bSaveSummaries )
      Common::SummariseBootstrapCorr( out, sOutFileBase, seed );//, momentum_squared );
  }
  // Dump the bootstrap sequences
  /*if( dumpBoot )
  {
    std::string SeqFileName{ Common::MakeFilename( sOutFileBase, "bootseq", seed, TEXT_EXT ) };
    std::cout << "Saving bootstrap sequence to '" << SeqFileName << "' ...";
    std::ofstream file(SeqFileName);
    file << "# bootstrap sequences" << std::endl;
    file << "#      bin size: " << binSize << std::endl;
    file << "#          seed: " << std::to_string( seed ) << std::endl;
    in.dumpBootstrapSeq(file, nSample, seed);
    std::cout << std::endl;
  }*/
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

enum StudySubject { Z2, GFWall };

void Study1Bootstrap( StudySubject Study, const BootstrapParams &bsParams, const Common::Momentum &mom, const std::string &Heavy, const std::string &Light )
{
  static constexpr unsigned int Nt{64};
  static constexpr unsigned int CStart{3000};
  static constexpr unsigned int CSkip{40};
  static constexpr unsigned int CCount{3};//10};
  static constexpr unsigned int CEnd{CStart + CSkip * CCount};
  static const std::string sPath{
#ifdef DEBUG
    "/Users/s1786208"               // Laptop
#else
    "/home/dp008/dp008/dc-mars3"  // Tesseract
#endif
    "/data/201910Plat/mesons/C1/" };
  const std::string StudyPath{ sPath + ( Study == Z2 ? "Z2" : "GFWall" ) + "/" };
  static const std::string Dot{ "." };    // used inside filenames
  static const std::string H5{ Dot + "h5" };    // used inside filenames
  static const std::string Sep{ "_" };    // used inside filenames
  static const std::string Space{ " " };  // whitespace as field separator / human readable info
  static const std::string Sink{ Sep + "sink" };    // used inside filenames
  static const std::string sProp{ "prop" + Sep };    // used inside filenames
  static const std::string sMom0{ Sep + "p" + Sep + Common::Momentum(0,0,0).to_string( Sep ) };
  const std::string sMom{ Sep + "p" + Sep + mom.to_string( Sep )};
  const std::string sMomNeg{ Sep + "p" + Sep + mom.to_string( Sep, true )};
  assert( alg.size() == algNames.size() && "Abbreviations should match gamma algebra" );
  const std::string CorrPrefix{ Heavy + Sep + Light + sMom };
  std::vector<std::string> CorrSuffixes( NumCorr );
  for( int iSrc = 0; iSrc < NumAlgebra; iSrc++ )
    for( int iSnk = 0; iSnk < NumAlgebra; iSnk++ ) {
      const int iCorr{ iSnk * NumAlgebra + iSrc };
      std::stringstream ss;
      ss << Sep << algNames[iSnk] << Sep << algNames[iSrc];
      CorrSuffixes[iCorr] = ss.str();
    }
  const int TimeSliceInc{ 4 };//Study == GFWall ? 4 : 1 };
  const unsigned int NumEnds{ Common::EqualIgnoreCase( Heavy, Light ) ? 1u : 2u };
  const std::size_t NumSamplesT{ NumEnds * CCount };
  const std::size_t NumSamples{ NumSamplesT * ( Nt / TimeSliceInc ) };
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
  for( int t = 0; t < Nt; t += TimeSliceInc ) {
    const std::string tPrefix{ Sep + "t" + ( Study == Z2 ? "" : Sep ) + std::to_string( t ) };
    int CorrIndexT{ 0 };
    for( unsigned int iEnd = 0; iEnd < NumEnds; iEnd++ ) {
      const std::string & Left{ iEnd == 0 ? Heavy : Light };
      const std::string & Right{ iEnd == 0 ? Light : Heavy };
      for( int iConfig = CStart; iConfig < CEnd; iConfig += CSkip )
      {
        std::string sFileName{ StudyPath };
        switch( Study )
        {
          case Z2:
            sFileName.append( sProp + Left + tPrefix + sMom + Sep + sProp + Right + tPrefix + sMom0 + Sink + sMomNeg );
            break;
          default:
            sFileName.append( Left + Sep + Right + sMom + tPrefix );
            break;
        }
        sFileName.append( Dot + std::to_string( iConfig ) + H5 );
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
    // If there's more than one configuration, perform a bootstrap of this timeslice
    for( int iCorr = 0; CCount > 1 && iCorr < NumCorr; ++iCorr )
      par.PerformBootstrap( bsDataT[iCorr], CorrPrefix + tPrefix + CorrSuffixes[iCorr] );
  }
  // Now perform a bootstrap overall
  par.bSaveBootstrap = true;
  for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
    par.PerformBootstrap( bsData[iCorr], CorrPrefix + CorrSuffixes[iCorr] );
}

/*****************************************************************

 Perform a bootstrap of all the files specified on the command line.
 Organise all the files by contraction, and sort them by trajectory.
 Ensure trajectories are processed in the same order so the
 bootstrap replicas for each contraction can be combined later during analysis

*****************************************************************/

#ifdef DEBUG
int DebugTest1()
{
  std::string File1{ "../../mesons/C1/Z2/prop_h1_t0_p_0_0_0_prop_h1_t0_p_0_0_0_sink_p_0_0_0.3000.h5" };
  std::string File2{ "../../../201907HL/contract/c2p0_h1_l_p2_0/sink_0_l_0.0.source_0_0.3200.h5" };
  Common::CorrelatorFile<std::complex<double>> f;
  std::string GroupName;
  std::vector<Common::Gamma::Algebra> Alg{ Common::Gamma::Algebra::Gamma5, Common::Gamma::Algebra::GammaTGamma5};
  f.Read( File1, GroupName, Alg );
  Alg.clear();
  GroupName.clear();
  f.Read( File1, GroupName, Alg );
  Alg.clear();
  GroupName.clear();
  f.Read( File2, GroupName, Alg );
  return 1;
}
#endif

int main(const int argc, const char *argv[])
{
  // int i = DebugTest1(); if( i ) return i;
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ false };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"n", CL::SwitchType::Single, DEF_NSAMPLE},
      {"b", CL::SwitchType::Single, "1"},
      {"r", CL::SwitchType::Single, nullptr},
      {"o", CL::SwitchType::Single, "@corr@" },
      {"t", CL::SwitchType::Flag, nullptr},
      {"x", CL::SwitchType::Multiple, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( cl.GotSwitch( "help" ) )//|| cl.Args.size() == 0 )
      bShowUsage = true;
    else
    {
      BootstrapParams par;
      par.bSaveBootstrap = true;
      par.outStem = cl.SwitchValue<std::string>( "o" );
      if( par.outStem.find( "@corr@" ) == std::string::npos )
      {
        std::string Error{ "Output stem \"" };
        Error.append( par.outStem );
        Error.append( "\" invalid" );
        throw std::invalid_argument( Error );
      }
      par.nSample = cl.SwitchValue<int>( "n" );
      par.binSize = cl.SwitchValue<int>( "b" );
      if( par.binSize != 1 )
        throw std::invalid_argument( "Binning not implemented yet" );
      if( cl.GotSwitch( "r" ) )
        par.seed = cl.SwitchValue<Common::SeedType>( "r" );
      else
      {
        std::random_device rd;
        par.seed = rd();
      }
      par.bSaveSummaries = !cl.GotSwitch( "t" );
      std::cout << "Creating bootstrap output " << par.outStem << ", seed=" << par.seed << std::endl;
      
      // If there are files specified on the command line,
      // parse the input file names, grouping by correlator, indexed by trajectory.
      if( cl.Args.size() )
      {
        Common::Manifest Manifest{ cl.Args, cl.SwitchStrings( "x" ) };
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
          //std::vector<double> buf;
          Latan::Dataset<Latan::DMat> data(nFile);
          for( auto it = l.TrajFile.begin(); it != l.TrajFile.end(); it++, j++ )
          {
            const int &traj{it->first}; // Trajectory number
            const Common::TrajFile &tf{*it->second};
            const std::string &Filename{tf.Filename};
            std::cout << "\t" << traj << "\t" << Filename << std::endl;
            std::string GroupName;
            Common::ReadComplexArray( buf, Filename, GroupName );
            //Common::ReadRealArray( buf, Filename, "/", "bootstrap" );
            Common::CopyCorrelator( data[j], buf );
            //Common::CopyRealCorrelator( data[j], buf );
          }
          par.PerformBootstrap( data, Contraction );
          BootstrapCount++;
        }
        std::cout << "Bootstraps written for " << BootstrapCount
        << " / " << Manifest.size() << " correlators" << std::endl;
      }
      else
      {
        // Nothing specified on the command-line
        // Perform default bootstraps
        // Heavy-light semi-leptonics
        const std::string qHeavy{ "h1" };
        const std::string qLight{ "l" };
        std::cout << "No files specified. Making default manifest for "
        << qHeavy << " and " << qLight << std::endl;
        static const Common::Momentum StudyMomenta[] = {
          { 0, 0, 0 },
          { 1, 0, 0 },/*
                       { 1, 1, 0 },
                       { 1, 1, 1 },
                       { 2, 0, 0 },*/
        };
        const StudySubject Study{ GFWall };
        for( const Common::Momentum &m : StudyMomenta )
        {
          Study1Bootstrap( Study, par, m, qHeavy, qLight );
          Study1Bootstrap( Study, par, m, qHeavy, qHeavy );
          Study1Bootstrap( Study, par, m, qLight, qLight );
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
  }
  if( bShowUsage || iReturn != EXIT_SUCCESS )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> ContractionFile1 [ContractionFile2 ...]\n"
    "Perform a bootstrap of the specified files, where <options> are:\n"
    "-n     Number of samples\n"
    "-b     Bin size (not implemented yet)\n"
    "-r     Random number seed (unspecified=random seed)\n"
    "-o     Output prefix\n"
    "-t     don't save .Txt correlators for GNUplot\n"
    "-x     eXclude file (may be repeated)\n"
    "--help This message\n";
  }
  return iReturn;
}
