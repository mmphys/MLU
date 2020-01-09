/**

 Mike's lattice QCD utilities
 
 Source file: bootstrap.cpp
 
 Copyright (C) 2019
 
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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
**/

// Perform a bootstrap

#include "Common.hpp"
using namespace Common;

//#include "CommonGrid.hpp"
//#include "CommonLatAn.hpp"

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <set>
//#include <LatAnalyze/Core/OptParser.hpp>
//#include <LatAnalyze/Statistics/Dataset.hpp>
//#include <LatAnalyze/Io/Io.hpp>
//#include <LatAnalyze/Io/Hdf5File.hpp>

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif

// Show me the averages for each timeslice
/*void ShowTimeSliceAvg(const Latan::Dataset<Latan::DMat> &data) {
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
}*/

struct TrajFile
{
  bool bHasTimeslice = false;
  int  Timeslice;
  TrajFile( bool bHasTimeslice_, int  Timeslice_ ) : bHasTimeslice{ bHasTimeslice_ }, Timeslice{ Timeslice_ } {}
};

// This describes one contraction and each of its trajectory files
struct TrajList
{
  std::string Name;                     // name of the contraction
  std::map<std::string, TrajFile> FileInfo; // Filenames with corresponding timeslice info
  TrajList( const std::string &Name_ ) : Name{ Name_ } {}
};

// This is a list of all the contractions we've been asked to process
class Manifest : public std::map<std::string, TrajList> {
  void Construct( const std::vector<std::string> &Files, const std::vector<std::string> &Ignore );
public:
  // Process list of files on the command-line, breaking them up into individual trajectories
  Manifest( const std::vector<std::string> &Files, const std::string &sIgnore );
  Manifest( const std::vector<std::string> &Files, const std::vector<std::string> &Ignore )
  { Construct( Files, Ignore ); }
};

/*enum ExtractFilenameReturn {Good, Bad, No_trajectory};

// Extract the contraction name and trajectory number from filename
static ExtractFilenameReturn ExtractFilenameParts(const std::string &Filename, std::string &Contraction, int &traj)
{
  ExtractFilenameReturn r{Bad};
  Contraction.clear();
  traj = 0;
  auto pos = Filename.find_last_of('/');
  if( pos == std::string::npos )
    pos = 0;
  else
    pos++;
  auto delim = Filename.find_first_of('.', pos);
  if( delim != std::string::npos )
  {
    Contraction = Filename.substr(pos, delim - pos);
    pos = delim + 1;
    delim = Filename.find_first_of('.', pos);
    if( delim != std::string::npos && delim != pos )
    {
      r = Good;
      while( r == Good && pos < delim )
      {
        auto c = Filename[pos++] - '0';
        if( c >=0 && c < 10 )
          traj = traj * 10 + c;
        else
          r = No_trajectory;
      }
    }
  }
  return r;
}*/

void Manifest::Construct(const std::vector<std::string> &Args, const std::vector<std::string> &Ignore)
{
  // Now walk the list of arguments.
  // Any file that's not in the ignore list gets added to the manifest
  if( Args.size() == 0 )
    return;
  bool parsed = true;
  std::map<std::string, TrajList> & Contractions = (* this);
  for( const std::string &Filename : Args )
  {
    // See whether this file is in the ignore list
    std::size_t iIgnore = 0;
    while( iIgnore < Ignore.size() && Ignore[iIgnore].compare(Filename) )
      iIgnore++;
    if( iIgnore < Ignore.size() )
      std::cout << "Ignoring " << Filename << std::endl;
    else if( !FileExists(Filename))
    {
      parsed = false;
      std::cout << "Error: " << Filename << " doesn't exist" << std::endl;
    }
    else
    {
      // Parse the name. Not expecting a type, so if present, put it back on the end of Base
      FileNameAtt Name_{ Filename };
      if( !Name_.Type.empty() )
      {
        // Not sure whether I should bother doing this?
        Name_.Base.append( 1, '.' );
        Name_.Base.append( Name_.Type );
        Name_.Type.clear();
      }
      std::string Contraction{ Name_.Base };
      bool bHasTimeslice{ false };
      int Timeslice_{ 0 };
      ExtractTimeslice( Contraction, bHasTimeslice, Timeslice_ );
      // Replace momentum with momentum squared
      bool bGotMomentum{ false };
      int p2{ 0 };
      ExtractP2( Contraction, bGotMomentum, p2 );
      // Sort all the separate, underscore delimited component parts
      {
        for( char &c : Contraction )
          if( c == '_' )
            c = ' ';
        std::vector<std::string> vs{ ArrayFromString<std::string>( Contraction ) };
        std::sort(vs.begin(), vs.end());
        bool bFirst = true;
        for( const std::string &s : vs )
        {
          if( bFirst )
          {
            bFirst = false;
            Contraction = s;
          }
          else
          {
            Contraction.append(1, '_');
            Contraction.append( s );
          }
        }
      }
      // Add the momentum into the correlator if it was present
      if( bGotMomentum )
      {
        Contraction.append( "_p2_" );
        Contraction.append( std::to_string( p2 ) );
      }
      // Look for the contraction list this file belongs to
      auto itc = Contractions.find( Contraction );
      if( itc == Contractions.end() )
        itc = Contractions.emplace( Contraction, TrajList( Contraction ) ).first;
      TrajList & cl{ itc->second };
      auto it = cl.FileInfo.find( Name_.Filename );
      if( it == cl.FileInfo.end() )
        cl.FileInfo.emplace( Name_.Filename, TrajFile( bHasTimeslice, Timeslice_ ) );
      else
        std::cout << "Ignoring repetition of " << Name_.Filename << std::endl;
    }
  }
  if( !parsed )
    throw std::runtime_error( "Error parsing command-line arguments" );
}

Manifest::Manifest(const std::vector<std::string> &Args, const std::string &sIgnore)
{
  // Get the list of files to ignore
  std::vector<std::string> Ignore;
  for( std::size_t start = 0; start < sIgnore.length() ; ) {
    auto end = sIgnore.find(':', start);
    std::size_t SepLen;
    if( end == std::string::npos ) {
      SepLen = 0;
      end = sIgnore.length();
    }
    else
      SepLen = 1;
    if( end > start )
      Ignore.push_back( sIgnore.substr( start, end - start ) );
    start = end + SepLen;
  }
  Construct( Args, Ignore );
}

class BootstrapParams
{
public:
  std::string outStem;
  Common::SeedType seed;
  int nSample;
  int binSize;
  int TimesliceDetail;
  int PerformBootstrap( std::vector<Common::CorrelatorFileC> &f, const std::string Prefix ) const;
private:
  template <class Iter> int PerformBootstrap( const Iter &first, const Iter &last, const std::string &Prefix, bool bAlignTimeslices, bool bSaveBootstrap, bool bSaveSummaries ) const;
};

template <class Iter> int BootstrapParams::PerformBootstrap( const Iter &first, const Iter &last, const std::string &Prefix, bool bAlignTimeslices, bool bSaveBootstrap, bool bSaveSummaries ) const
{
  if( !bSaveSummaries && !bSaveBootstrap )
    return 0;
  const int NumFiles{ static_cast<int>( std::distance( first, last ) ) };
  int iCount{ 0 };
  const int Nt{ first->Nt() };
  const int NumOps{ first->NumOps() };
  const std::vector<Common::Gamma::Algebra> &Alg{ first->Alg() };
  Common::SampleC in( NumFiles, Nt );
  for( int Snk = 0; Snk < NumOps; Snk++ )
  {
    static const char pszSep[] = "_";
    std::string sSnk{ Common::Gamma::NameShort( Alg[Snk], pszSep ) };
    for( int Src = 0; Src < NumOps; Src++ )
    {
      std::string sSrc{ Common::Gamma::NameShort( Alg[Src], pszSep ) };
      // Skip bootstrap if output exists
      const std::string sOutBase{ outStem + Prefix + sSnk + sSrc };
      const std::string sOutFile{ Common::MakeFilename( sOutBase, Common::sBootstrap, seed, DEF_FMT ) };
      if( Common::FileExists( sOutFile ) || (!bSaveBootstrap && Common::FileExists( Common::MakeFilename(sOutBase, CorrSumm::SummaryNames[0], seed, TEXT_EXT))))
      {
        std::cout << sOutBase << " skipped - output already exists" << std::endl;
      }
      else
      {
        //if( bShowAverages )
        //ShowTimeSliceAvg( in );
        // std::cout << sOutBase << " " << nSample << " samples" << std::endl;
        assert( binSize == 1 && "Binning still needs to be implemented" );
        // Copy the correlators into my sample, aligning timeslices if need be
        std::complex<double> * pDst = in[0];
        for( Iter i = first; i != last; i++ )
        {
          const int CorrelatorTimeslice{ i->Timeslice() };
          int TOffset{ bAlignTimeslices ? CorrelatorTimeslice : 0 };
          // Only need to say which correlators contribute for first gamma structure
          if( Src == 0 && Snk == 0 )
          {
            std::cout << "  t=";
            if( TOffset )
              std::cout << TOffset << "->0";
            else
              std::cout << CorrelatorTimeslice;
            std::cout << '\t' << i->Name_.Base << "." << i->Name_.SeedString << std::endl;
          }
          std::complex<double> * pSrc = (*i)( Alg[Snk], Alg[Src] ) + TOffset;
          for( int t = 0; t < Nt; t++ )
          {
            *pDst++ = *pSrc++;
            if( ++TOffset == Nt )
              pSrc -= Nt;
          }
        }
        Common::SampleC out = in.Bootstrap( nSample, seed );
        if( bSaveBootstrap )
        {
          std::cout << "  " << nSample << " samples to " << sOutFile << std::endl;
          out.Write( sOutFile );
        }
        if( bSaveSummaries )
          Common::SummariseBootstrapCorr( out, sOutBase, seed );//, momentum_squared );
        iCount++;
      }
    }
  }
  return iCount;
}

int BootstrapParams::PerformBootstrap( std::vector<Common::CorrelatorFileC> &f, const std::string Prefix ) const
{
  auto NumFilesRaw{ f.size() };
  const int NumFiles{ static_cast<int>( NumFilesRaw ) };
  if( NumFiles != NumFilesRaw || NumFiles < 1 )
    throw std::invalid_argument( "Can't perform a bootstrap without correlators" );
  // Make sure there are at least 5x more bootstrap samples than data points
  if( ( nSample / 5 ) < NumFiles )
    throw std::runtime_error( "nSample=" + std::to_string( nSample ) + " too small for dataset with " + std::to_string( NumFiles ) + " elements" );
  int iCount{ 0 };
  // sort all of the files by timeslice info
  using File = Common::CorrelatorFileC;
  std::sort( f.begin(), f.end(), [&](const File &l, const File &r)
  {
    // Sort first by timeslice
    if( l.bHasTimeslice != r.bHasTimeslice )
      return r.bHasTimeslice;
    if( l.bHasTimeslice && r.bHasTimeslice && l.Timeslice_ != r.Timeslice_ )
      return l.Timeslice_ < r.Timeslice_;
    // Sort by trajectory
    if( l.Name_.Seed != r.Name_.Seed )
      return l.Name_.Seed < r.Name_.Seed;
    // Now sort by base filename
    return l.Name_.Base.compare( r.Name_.Base ) < 0;
  } );
  // Now perform bootstraps for any individual timeslices
  using Iter = typename std::vector<Common::CorrelatorFileC>::iterator;
  Iter first = f.begin();
  Iter last  = f.end();
  if( TimesliceDetail > 0 )
  {
    // If more than one timeslice, save summary info for individual timeslices
    int Timeslice{ f[0].Timeslice() };
    Iter i = std::find_if( first + 1, last, [Timeslice]( const File &cf ) { return cf.Timeslice() != Timeslice; } );
    if( i != last )
    {
      while( first != last )
      {
        std::string PrefixT{ Prefix };
        PrefixT.append( "_t_" );
        PrefixT.append( std::to_string( Timeslice ) );
        iCount += PerformBootstrap( first, i, PrefixT, false, TimesliceDetail > 1, true );
        first = i;
        if( first != last )
        {
          Timeslice = first->Timeslice();
          i = std::find_if( first + 1, last, [Timeslice]( const File &cf ) { return cf.Timeslice() != Timeslice; } );
        }
      }
      first = f.begin();
    }
  }
  // Now perform a single bootstrap, combining all the separate timeslices
  iCount += PerformBootstrap( first, last, Prefix, true, true, true );
  return iCount;
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

enum StudySubject { Z2=1, GFWall };

void Study1Bootstrap( StudySubject Study, const BootstrapParams &bsParams, const Common::Momentum &mom,
          std::vector<Common::Gamma::Algebra> Alg, const std::string &Heavy, const std::string &Light )
{
  /*using Algebra = Common::Gamma::Algebra;
  static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
  static const std::vector<std::string> algNames = { "g5", "gT5" };
  static const int NumAlgebra{ static_cast<int>( alg.size() ) };
  static const int NumCorr{ NumAlgebra * NumAlgebra };*/

  static constexpr unsigned int Nt{64};
  static constexpr unsigned int CStart{3000};
  static constexpr unsigned int CSkip{40};
  static constexpr unsigned int CCount{1};//10};
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
  static const std::string H5{ Dot + DEF_FMT };    // used inside filenames
  static const std::string Sep{ "_" };    // used inside filenames
  static const std::string Space{ " " };  // whitespace as field separator / human readable info
  static const std::string Sink{ Sep + "sink" };    // used inside filenames
  static const std::string sProp{ "prop" + Sep };    // used inside filenames
  static const std::string sMom0{ Sep + "p" + Sep + Common::Momentum(0,0,0).to_string( Sep ) };
  const std::string sMom{ Sep + "p" + Sep + mom.to_string( Sep )};
  const std::string sMomNeg{ Sep + "p" + Sep + mom.to_string( Sep, true )};
  //assert( alg.size() == algNames.size() && "Abbreviations should match gamma algebra" );
  const std::string CorrPrefix{ Heavy + Sep + Light + sMom };
  /*std::vector<std::string> CorrSuffixes( NumCorr );
  for( int iSrc = 0; iSrc < NumAlgebra; iSrc++ )
    for( int iSnk = 0; iSnk < NumAlgebra; iSnk++ ) {
      const int iCorr{ iSnk * NumAlgebra + iSrc };
      std::stringstream ss;
      ss << Sep << algNames[iSnk] << Sep << algNames[iSrc];
      CorrSuffixes[iCorr] = ss.str();
    }*/
  const int TimeSliceInc{ Study == GFWall ? 4 : 1 };
  const unsigned int NumEnds{ Common::EqualIgnoreCase( Heavy, Light ) ? 1u : 2u };
  /*const std::size_t NumSamplesT{ NumEnds * CCount };
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
  par.bSaveBootstrap = false;*/
  int CorrIndex{ 0 };
  std::vector<Common::CorrelatorFileC> InFiles( NumEnds * CCount * ( Nt / TimeSliceInc ) );
  for( int t = 0; t < Nt; t += TimeSliceInc ) {
    const std::string tPrefix{ Sep + "t" + ( Study == Z2 ? "" : Sep ) + std::to_string( t ) };
    //int CorrIndexT{ 0 };
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
        /*Common::ReadComplexArray(buffer, alg, sFileName, 0);
        for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
        {
          Common::CopyCorrelator( bsData[iCorr][CorrIndex], buffer[iCorr], t );
          Common::CopyCorrelator( bsDataT[iCorr][CorrIndexT], buffer[iCorr] );
        }
        CorrIndex++;
        CorrIndexT++;*/
        std::string GroupName;
        InFiles[CorrIndex++].Read( sFileName, GroupName, Alg, &t );
      }
    }
    // If there's more than one configuration, perform a bootstrap of this timeslice
    //for( int iCorr = 0; CCount > 1 && iCorr < NumCorr; ++iCorr )
      //par.PerformBootstrap( bsDataT[iCorr], CorrPrefix + tPrefix + CorrSuffixes[iCorr] );
  }
  // Now perform a bootstrap overall
  /*par.bSaveBootstrap = true;
  for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
    par.PerformBootstrap( bsData[iCorr], CorrPrefix + CorrSuffixes[iCorr] );*/
  bsParams.PerformBootstrap( InFiles, Heavy + Sep + Light );
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
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"n", CL::SwitchType::Single, DEF_NSAMPLE},
      {"b", CL::SwitchType::Single, "1"},
      {"r", CL::SwitchType::Single, nullptr},
      {"o", CL::SwitchType::Single, "" },
      {"a", CL::SwitchType::Single, ""},
      {"s", CL::SwitchType::Single, nullptr},
      {"t", CL::SwitchType::Single, "1"},
      {"x", CL::SwitchType::Multiple, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) )
    {
      std::vector<Common::Gamma::Algebra> Alg = Common::ArrayFromString<Common::Gamma::Algebra>( cl.SwitchValue<std::string>( "a" ) );
      BootstrapParams par;
      par.TimesliceDetail = cl.SwitchValue<int>( "t" );
      if( par.TimesliceDetail < 0 || par.TimesliceDetail > 2 )
        throw std::invalid_argument( "Timeslice detail " + std::to_string( par.TimesliceDetail ) + " invalid" );
      par.outStem = cl.SwitchValue<std::string>( "o" );
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

      // If there are files specified on the command line,
      // parse the input file names, grouping by correlator, indexed by trajectory.
      const bool bGotStudy{ cl.GotSwitch( "s" ) };
      if( cl.Args.size() )
      {
        if( bGotStudy )
          throw std::invalid_argument( "Can't specify study " + cl.SwitchValue<std::string>( "s" )
                                      + " with command-line arguments" );
        bShowUsage = false;
        Manifest Manifest{ cl.Args, cl.SwitchStrings( "x" ) };
        // Walk the list of contractions, performing a separate bootstrap for each
        int BootstrapCount = 0;
        for( auto itc = Manifest.begin(); itc != Manifest.end(); itc++ )
        {
          const std::string &Contraction{itc->first};
          const TrajList &l{itc->second};
          const unsigned int nFile{ static_cast<unsigned int>( l.FileInfo.size() ) };
          std::cout << "Loading " << nFile << " files for " << Contraction << std::endl;
          unsigned int j = 0; // which trajectory are we processing
          std::vector<Common::CorrelatorFileC> InFiles( nFile );
          for( auto it = l.FileInfo.begin(); it != l.FileInfo.end(); it++, j++ )
          {
            const std::string &Filename{ it->first }; // Trajectory number
            const TrajFile &tf{ it->second };
            std::string GroupName;
            InFiles[j].Read( Filename, GroupName, Alg, tf.bHasTimeslice ? &tf.Timeslice : nullptr );
          }
          par.PerformBootstrap( InFiles, Contraction );
          BootstrapCount++;
        }
        std::cout << "Bootstraps written for " << BootstrapCount
        << " / " << Manifest.size() << " correlators" << std::endl;
      }
      else if( bGotStudy )
      {
        StudySubject Study{ static_cast<StudySubject>( cl.SwitchValue<int>( "s" ) ) };
        switch( Study )
        {
          case Z2:
          case GFWall:
          {
            // Nothing specified on the command-line
            // Perform default bootstraps
            // Heavy-light semi-leptonics
            const std::string qHeavy{ "h1" };
            const std::string qLight{ "l" };
            std::cout << "No files specified. Making default manifest for "
            << qHeavy << " and " << qLight << std::endl;
            static const Common::Momentum StudyMomenta[] = {
              { 0, 0, 0 },/*
              { 1, 0, 0 },
              { 1, 1, 0 },
              { 1, 1, 1 },
              { 2, 0, 0 },*/
            };
            for( const Common::Momentum &m : StudyMomenta )
            {
              Study1Bootstrap( Study, par, m, Alg, qHeavy, qLight );
              Study1Bootstrap( Study, par, m, Alg, qHeavy, qHeavy );
              Study1Bootstrap( Study, par, m, Alg, qLight, qLight );
            }
          }
            bShowUsage = false;
            break;
          default:
            std::cout << "Study " << Study << " undefined\n";
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
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> ContractionFile1 [ContractionFile2 ...]\n"
    "Perform a bootstrap of the specified files, where <options> are:\n"
    "-n     Number of samples (" DEF_NSAMPLE ")\n"
    "-b     Bin size (not implemented yet)\n"
    "-r     Random number seed (unspecified=random)\n"
    "-o     Output prefix\n"
    "-a     list of gamma algebras we're interested in\n"
    "-s     Perform bootstrap for specified study number\n"
    "-t     timeslice detail 0 (none), 1 (.txt=default) or 2 (.txt+.h5)\n"
    "-x     eXclude file (may be repeated)\n"
    "--help This message\n";
  }
  return iReturn;
}
