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

enum StudySubject{ Z2=1, GFWW, GFPW };
enum GroupMomenta{ None, Squared, Abs };

enum BinOrder{ Auto, Old, VeryOld };

std::istream& operator>>(std::istream& is, BinOrder &binOrder )
{
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, "Old" ) )
      binOrder = BinOrder::Old;
    else if( EqualIgnoreCase( s, "VeryOld" ) )
      binOrder = BinOrder::VeryOld;
    else if( EqualIgnoreCase( s, "Auto" ) )
      binOrder = BinOrder::Auto;
    else
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

struct OperatorAttributes
{
  Common::Gamma::Algebra alg;
  Common::RPS rps;
};

const OperatorAttributes DefaultOperatorAttributes[]{
  {Common::Gamma::Algebra::Gamma5, { Common::Reality::Real, Common::Parity::Even, Common::Sign::Negative} },
  {Common::Gamma::Algebra::GammaTGamma5, { Common::Reality::Real, Common::Parity::Odd,  Common::Sign::Positive} },
  {Common::Gamma::Algebra::GammaXGamma5, { Common::Reality::Imag, Common::Parity::Even, Common::Sign::Positive} },
  {Common::Gamma::Algebra::GammaYGamma5, { Common::Reality::Imag, Common::Parity::Even, Common::Sign::Positive} },
  {Common::Gamma::Algebra::GammaZGamma5, { Common::Reality::Imag, Common::Parity::Even, Common::Sign::Positive} },
};
static constexpr int NumDefaultOperatorAttributes{ sizeof( DefaultOperatorAttributes ) / sizeof( DefaultOperatorAttributes[0] )  };

// Get the attributes of the specified product of operators
const Common::RPS &DefaultOperatorAttribute( Common::Gamma::Algebra Alg )
{
  for( int i = 0; i < NumDefaultOperatorAttributes; i++ )
    if( DefaultOperatorAttributes[i].alg == Alg )
      return DefaultOperatorAttributes[i].rps;
  std::stringstream ss;
  ss << "I don't know how to fold " << Alg;
  throw std::runtime_error( ss.str() );
}

#ifdef DEBUG_TEST
bool Debug()
{
  for( int i = 0; i < NumDefaultOperatorAttributes; i++ )
    for( int j = 0; j < NumDefaultOperatorAttributes; j++ )
    {
      Common::RPS rps{ DefaultOperatorAttributes[i].rps * DefaultOperatorAttributes[j].rps };
      std::cout << DefaultOperatorAttributes[i].alg << "-" << DefaultOperatorAttributes[j].alg
                << ": " << rps << "\n";
    }
  return true;
}
#endif

struct TrajFile
{
  bool bHasTimeslice = false;
  int  Timeslice;
  bool bGotMomentum = false;
  Common::Momentum p;
  int Config;
  TrajFile( bool bHasTimeslice_, int Timeslice_, bool bGotMomentum_, Common::Momentum p_, int Config_ )
  : bHasTimeslice{bHasTimeslice_}, Timeslice{Timeslice_}, bGotMomentum{bGotMomentum_}, p{p_}, Config{Config_} {}
};

// This describes one contraction and each of its trajectory files
struct TrajList
{
  std::string Name;                     // name of the contraction
  std::string sShortPrefix;
  std::string sShortSuffix;
  std::vector<Common::Gamma::Algebra> Alg;
  bool bHasDeltaT;
  int DeltaT;
  std::map<std::string, TrajFile> FileInfo; // Filenames with corresponding timeslice info
  TrajList(const std::string &Name_, const std::string &sShortPrefix_, const std::string &sShortSuffix_,
           std::vector<Common::Gamma::Algebra> Alg_, bool bHasDeltaT_, int DeltaT_)
  : Name{Name_}, sShortPrefix{sShortPrefix_}, sShortSuffix{ sShortSuffix_}, Alg{Alg_},
    bHasDeltaT{bHasDeltaT_}, DeltaT{DeltaT_} {}
};

// This is a list of all the contractions we've been asked to process
struct Manifest : public std::map<std::string, TrajList>
{
  // Process list of files on the command-line, breaking them up into individual trajectories
  Manifest(const std::vector<std::string> &Files, const std::vector<std::string> &Ignore,
           bool bSwapQuarks, const GroupMomenta GroupP);
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

Manifest::Manifest(const std::vector<std::string> &Args, const std::vector<std::string> &Ignore,
                   bool bSwapQuarks, const GroupMomenta GroupP)
{
  static const std::string Sep{ "_" };
  // Now walk the list of arguments.
  // Any file that's not in the ignore list gets added to the manifest
  if( Args.empty() )
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
      if( !Name_.bSeedNum )
        throw std::runtime_error( "Contraction files must contain a configuration number" );
      if( !Name_.Type.empty() )
      {
        // Not sure whether I should bother doing this?
        Name_.Base.append( 1, '.' );
        Name_.Base.append( Name_.Type );
        Name_.Type.clear();
      }
      // Extract attributes from contraction name
      std::string Contraction{ Name_.Base };
      Momentum mom;
      bool bGotMomentum{ mom.Extract( Contraction ) };
      bool bHasTimeslice{ false };
      int Timeslice_{ 0 };
      ExtractTimeslice( Contraction, bHasTimeslice, Timeslice_ );
      bool bHasDeltaT{ false };
      int DeltaT{ 0 };
      ExtractDeltaT( Contraction, bHasDeltaT, DeltaT );
      std::vector<Common::Gamma::Algebra> vAlg{ ExtractGamma( Contraction ) };
      if( vAlg.size() > 1 )
        throw std::runtime_error( "Multiple gamma insertions unsupported" );
      bool b3pt{ bHasDeltaT || !vAlg.empty() };
      // Sort all the separate, underscore delimited component parts
      if( !b3pt && bSwapQuarks )
      {
        std::vector<std::string> vs{ ArrayFromString( Contraction, Common::Underscore ) };
        Contraction.clear();
        std::sort(vs.begin(), vs.end());
        for( int i = 0; i < vs.size(); i++ )
        {
          if( i )
            Contraction.append(1, '_');
          Contraction.append( vs[i] );
        }
      }
      // Add attributes back into contraction string
      const std::string sShortPrefix{ Contraction };
      for( Common::Gamma::Algebra a : vAlg )
        Contraction.append( Gamma::NameShort( a, Sep.c_str() ) );
      if( bHasDeltaT )
      {
        Contraction.append( "_dt_" );
        Contraction.append( std::to_string( DeltaT ) );
      }
      // Add the momentum into the correlator if it was present
      std::string sShortSuffix;
      if( bGotMomentum )
      {
        sShortSuffix = "_p_";
        switch( GroupP )
        {
          case GroupMomenta::Squared:
            sShortSuffix = mom.p2_string( Sep );
            break;
          case GroupMomenta::Abs:
            sShortSuffix.append( mom.abs().to_string( Sep ) );
            break;
          default:
            sShortSuffix.append( mom.to_string( Sep ) );
            break;
        }
        Contraction.append( sShortSuffix );
      }
      // Look for the contraction list this file belongs to
      auto itc = Contractions.find( Contraction );
      if( itc == Contractions.end() )
        itc = Contractions.emplace( Contraction, TrajList( Contraction, sShortPrefix, sShortSuffix, std::move(vAlg), bHasDeltaT, DeltaT ) ).first;
      TrajList & cl{ itc->second };
      auto it = cl.FileInfo.find( Name_.Filename );
      if( it == cl.FileInfo.end() )
        cl.FileInfo.emplace(Name_.Filename,TrajFile( bHasTimeslice, Timeslice_, bGotMomentum, mom, Name_.Seed));
      else
        std::cout << "Ignoring repetition of " << Name_.Filename << std::endl;
    }
  }
  if( !parsed )
    throw std::runtime_error( "Remove non-existent files (or use '-x filename' to eXclude)" );
}

template <class Iter>
std::vector<Common::ConfigCount> CountConfigs( const Iter &first, const Iter &last )
{
  std::vector<Common::ConfigCount> cc;
  using CCMap = std::map<int, int>;
  CCMap ConfigCount;
  for( Iter it = first; it != last; ++it )
  {
    const int ConfigNum{ static_cast<int>( it->Name_.Seed ) };
    auto p = ConfigCount.find( ConfigNum );
    if( p == ConfigCount.end() )
      ConfigCount.insert( { ConfigNum, 1 } );
    else
      p->second++;
  }
  for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
    cc.emplace_back( it->first, it->second );
  return cc;
}

class BootstrapParams
{
public:
  std::string outStem;
  std::string MachineName;
  Common::SeedType seed;
  Common::SampleC RandomSample;
  int nSample;
  bool binAuto;
  BinOrder binOrder;
  int binSize;
  int TimesliceDetail;
  bool bFactorised;
  bool bVerboseSummaries;
  int PerformBootstrap(std::vector<Common::CorrelatorFileC> &f, const TrajList &Traj ) const;
  void Study1Bootstrap(StudySubject Study, const std::string &StudyPath, const Common::Momentum &mom,
                       std::vector<Common::Gamma::Algebra> Alg,
                       const std::string &Heavy, const std::string &Light, bool Factorised) const;

private:
  template <class Iter>
  int PerformBootstrap(const Iter &first, const Iter &last, const TrajList &Traj, const std::string &Suffix,
                       bool bAlignTimeslices, bool bSaveBootstrap, bool bSaveSummaries) const;
};

template <class Iter>
int BootstrapParams::PerformBootstrap(const Iter &first, const Iter &last, const TrajList &Traj,
                  const std::string &Suffix, bool bAlignTimeslices, bool bSaveBootstrap, bool bSaveSummaries) const
{
  using File = Common::CorrelatorFileC;
  static const std::string sBlanks{ "  " };
  const int NumFiles{ static_cast<int>( std::distance( first, last ) ) };
  if( NumFiles == 0 || ( !bSaveSummaries && !bSaveBootstrap ) )
    return 0;
  // Count how many files there are for each config - needs to work regardless of sort
  assert( sizeof( int ) == sizeof( Common::SeedType ) );
  using CCMap = std::map<int, int>;
  const std::vector<Common::ConfigCount> ConfigCount{ CountConfigs( first, last ) };
  const int NumConfigs{ static_cast<int>( ConfigCount.size() ) };
  // Sort ... unless we want the old sort order (e.g. to check old results can be replicated)
  if( binOrder == BinOrder::Auto )
  {
    // sort all of the files by config, then timeslice
    std::sort( first, last, [&](const File &l, const File &r)
    {
      // Sort by trajectory
      if( l.Name_.Seed != r.Name_.Seed )
        return l.Name_.Seed < r.Name_.Seed;
      // Sort by timeslice
      if( l.bHasTimeslice != r.bHasTimeslice )
        return r.bHasTimeslice;
      if( l.bHasTimeslice && r.bHasTimeslice && l.Timeslice_ != r.Timeslice_ )
        return l.Timeslice_ < r.Timeslice_;
      // Sort by base filename
      return l.Name_.Base.compare( r.Name_.Base ) < 0;
    } );
  }
  if( !binAuto && NumFiles % binSize )
    std::cout << "Warning: last bin partially filled (" << ( NumFiles % binSize ) << " of " << binSize << ")\n";
  const bool b3pt{ !Traj.Alg.empty() };
  std::string Prefix;
  if( b3pt )
  {
    if( Traj.Alg.size() > 1 )
      throw std::runtime_error( "I don't know how to handle more than one current insertion at a time" );
    Prefix = Traj.sShortPrefix;
  }
  else
  {
    Prefix = Traj.Name + Suffix;
  }
  int iCount{ 0 };
  const int Nt{ first->Nt() };
  const std::vector<Common::Gamma::Algebra> &AlgSrc{ first->AlgSrc() };
  const std::vector<Common::Gamma::Algebra> &AlgSnk{ first->AlgSnk() };
  for( int Snk = 0; Snk < first->NumSnk(); Snk++ )
  {
    static const char pszSep[] = "_";
    std::string sSnk{ Common::Gamma::NameShort( AlgSnk[Snk], pszSep ) };
    if( b3pt )
    {
      if( Traj.bHasDeltaT )
      {
        sSnk.append( "_dt_" );
        sSnk.append( std::to_string( Traj.DeltaT ) );
      }
      sSnk.append( Traj.sShortSuffix );
      sSnk.append( Suffix );
      sSnk.append( Common::Gamma::NameShort( Traj.Alg[0], pszSep ) );
    }
    for( int Src = 0; Src < ( ( !b3pt && bFactorised ) ? Snk + 1 : first->NumSrc() ); Src++ )
    {
      std::string sSrc{ Common::Gamma::NameShort( AlgSrc[Src], pszSep ) };
      // Skip bootstrap if output exists
      const std::string sOutBase{ outStem + Prefix + sSnk + sSrc };
      const std::string sOutFile{ Common::MakeFilename( sOutBase, sBootstrap, seed, DEF_FMT ) };
      const std::string sSummary{ Common::MakeFilename( sOutBase, sBootstrap, seed, TEXT_EXT)};
      if( Common::FileExists( sOutFile ) || ( !bSaveBootstrap && Common::FileExists( sSummary ) ) )
      {
        std::cout << sBlanks << sOutBase << " skipped - output already exists" << std::endl;
      }
      else
      {
        // Copy the correlators into my sample, aligning timeslices if need be
        const bool bDifferent{ Snk != Src };
        const int OpFactor{ ( !b3pt && bFactorised && bDifferent ) ? 2 : 1 };
        // If we are determining the bin size automatically, it depends on whether we have more than one config
        int binSize{ ( binAuto ? 1 : this->binSize ) * OpFactor };
        int NumBinnedSamples{ ( NumFiles + binSize - 1 ) / binSize };
        if( binAuto && NumConfigs > 1 ) // && bAlignTimeslices )
        {
          binSize = ConfigCount[0].Count * OpFactor;
          NumBinnedSamples = static_cast<int>( ConfigCount.size() );
        }
        Common::SampleC out( nSample, Nt );
        out.ConfigCount = ConfigCount;
        out.FileList.reserve( NumFiles );
        out.binSize = binSize;
        Common::SampleC in( NumBinnedSamples, Nt );
        int Bin{0};
        int WhichBin{0};
        std::complex<double> * pDst = in[0];
        for( Iter it = first; it != last; )
        {
          const Common::CorrelatorFileC &file{ *it++ };
          const int CorrelatorTimeslice{ file.Timeslice() };
          const int TOffset{ bAlignTimeslices ? CorrelatorTimeslice : 0 };
          out.FileList.emplace_back( file.Name_.Filename );
          for( int o = 0; o < OpFactor; o++ )
          {
            const std::complex<double> * const pSrc = file( AlgSnk[o ? Src : Snk], AlgSrc[o ? Snk : Src] );
            for( int t = 0; t < Nt; t++ )
            {
              const std::ptrdiff_t idx{ ( t + TOffset ) % Nt };
              std::complex<double> z{ pSrc[ idx ] };
              if( Bin )
                pDst[t] += z;
              else
                pDst[t]  = z;
            }
            // If we are binning, divide by number of items in this bin. Make sure to do the last bin
            if( ++Bin == binSize || ( it == last && o == OpFactor - 1 ) )
            {
              if( Bin == 1 )
                pDst += Nt;
              else
                for( int t = 0; t < Nt; t++ )
                  *pDst++ /= Bin;
              Bin = 0;
              if( it != last && binAuto && NumConfigs > 1 )
              {
                const int newBinSize{ ConfigCount[++WhichBin].Count * OpFactor };
                if( out.binSize && out.binSize != newBinSize )
                  out.binSize = 0; // Indicates the bin size varies per ConfigCount
                binSize = newBinSize;
              }
            }
          }
        }
        std::cout << sBlanks << nSample << " samples to " << sOutFile << std::endl;
        if( RandomSample.RandNum() )
          in.Bootstrap( out, RandomSample );
        else
          in.Bootstrap( out, seed, &MachineName );
        // Now save the audit data for the bootstrap
        out.MakeCorrSummary( nullptr );
        if( bSaveBootstrap )
          out.Write( sOutFile );
        if( bSaveSummaries )
          out.WriteSummary( sSummary, bVerboseSummaries );
        iCount++;
      }
    }
  }
  return iCount;
}

int BootstrapParams::PerformBootstrap( std::vector<Common::CorrelatorFileC> &f, const TrajList &Traj ) const
{
  //const std::string& Prefix{ Traj.Name };
  auto NumFilesRaw{ f.size() };
  const int NumFiles{ static_cast<int>( NumFilesRaw ) };
  if( NumFiles != NumFilesRaw || NumFiles < 1 )
    throw std::invalid_argument( "Can't perform a bootstrap without correlators" );
  // Make sure there are at least 5x more bootstrap samples than data points
  //if( ( nSample / 5 ) < NumFiles )
    //throw std::runtime_error( "nSample=" + std::to_string( nSample ) + " too small for dataset with " + std::to_string( NumFiles ) + " elements" );
  int iCount{ 0 };
  using File = Common::CorrelatorFileC;
  // Now perform bootstraps for any individual timeslices
  using Iter = typename std::vector<Common::CorrelatorFileC>::iterator;
  Iter first = f.begin();
  Iter last  = f.end();
  // sort files by timeslice if we are doing timeslice summaries ... or old sort order asked for
  if( binOrder != BinOrder::Auto || TimesliceDetail > 0 )
  {
    std::sort( f.begin(), f.end(), [&](const File &l, const File &r)
    {
      // Sort first by timeslice
      if( l.bHasTimeslice != r.bHasTimeslice )
        return r.bHasTimeslice;
      if( l.bHasTimeslice && r.bHasTimeslice && l.Timeslice_ != r.Timeslice_ )
        return l.Timeslice_ < r.Timeslice_;
      // Very old sort order was config then filename
      if( binOrder == BinOrder::VeryOld )
      {
        if( l.Name_.Seed != r.Name_.Seed )
          return l.Name_.Seed < r.Name_.Seed;
        return l.Name_.Base.compare( r.Name_.Base ) < 0;
      }
      // More recent sort order is filename then config
      int iCompare = l.Name_.Base.compare( r.Name_.Base );
      if( iCompare )
        return iCompare < 0;
      return l.Name_.Seed < r.Name_.Seed;
    } );
  }
  if( TimesliceDetail > 0 )
    {
    // If more than one timeslice, save summary info for individual timeslices
    int Timeslice{ f[0].Timeslice() };
    Iter i = std::find_if( first + 1, last, [Timeslice]( const File &cf ) { return cf.Timeslice() != Timeslice; } );
    if( i != last )
    {
      while( first != last )
      {
        std::string PrefixT{ "_t_" };
        PrefixT.append( std::to_string( Timeslice ) );
        iCount += PerformBootstrap( first, i, Traj, PrefixT, false, TimesliceDetail > 1, true );
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
  iCount += PerformBootstrap( first, last, Traj, "", true, true, true );
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

void BootstrapParams::Study1Bootstrap(StudySubject Study, const std::string &StudyPath,
                      const Common::Momentum &mom, std::vector<Common::Gamma::Algebra> Alg,
                      const std::string &Heavy, const std::string &Light, bool Factorised ) const
{
  /*using Algebra = Common::Gamma::Algebra;
  static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
  static const std::vector<std::string> algNames = { "g5", "gT5" };
  static const int NumAlgebra{ static_cast<int>( alg.size() ) };
  static const int NumCorr{ NumAlgebra * NumAlgebra };*/

  static constexpr unsigned int Nt{64};
  static constexpr unsigned int CStart{3000};
  static constexpr unsigned int CSkip{40};
  const unsigned int CCount{ Study == Z2 ? 10u : 1u };//10};
  const unsigned int CEnd{CStart + CSkip * CCount};
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
  const int TimeSliceInc{ Study == Z2 ? 1 : 4 };
  const unsigned int NumEnds{ Factorised && !Common::EqualIgnoreCase( Heavy, Light ) ? 2u : 1u };
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
        InFiles[CorrIndex++].Read( sFileName, Alg, Alg, &t );
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
  const std::string sShortPrefix{ Heavy + Sep + Light };
  const std::string sShortSuffix{ mom.p2_string( Sep ) };
  PerformBootstrap( InFiles, TrajList( sShortPrefix + sShortSuffix, sShortPrefix, sShortSuffix, {}, false, 0 ) );
}

/*****************************************************************

 Perform a bootstrap of all the files specified on the command line.
 Organise all the files by contraction, and sort them by trajectory.
 Ensure trajectories are processed in the same order so the
 bootstrap replicas for each contraction can be combined later during analysis

*****************************************************************/

int main(const int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  const std::string MachineName{ Common::GetHostName() };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
#ifdef DEBUG_TEST
    if( Debug() )
      return iReturn;
#endif
    const std::initializer_list<CL::SwitchDef> list = {
      {"n", CL::SwitchType::Single, DEF_NSAMPLE},
      {"b", CL::SwitchType::Single, "0"},
      {"border", CL::SwitchType::Single, "Auto"},
      {"r", CL::SwitchType::Single, nullptr},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"a", CL::SwitchType::Single, ""},
      {"c", CL::SwitchType::Single, nullptr},
      {"g", CL::SwitchType::Single, "" },
      {"d", CL::SwitchType::Single, "" },
      {"s", CL::SwitchType::Single, nullptr},
      {"t", CL::SwitchType::Single, "0"},
      {"m", CL::SwitchType::Single, nullptr},
      {"x", CL::SwitchType::Multiple, nullptr},
      {"f", CL::SwitchType::Flag, nullptr},
      {"p2",CL::SwitchType::Flag, nullptr},
      {"pa",CL::SwitchType::Flag, nullptr},
      {"show", CL::SwitchType::Flag, nullptr},
      {"sort", CL::SwitchType::Flag, nullptr},
      {"terse", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) )
    {
      std::vector<Common::Gamma::Algebra> Alg = Common::ArrayFromString<Common::Gamma::Algebra>( cl.SwitchValue<std::string>( "a" ) );
      std::vector<bool> Alg3ptNeg;
      std::vector<Common::Gamma::Algebra> Alg3pt;
      if( cl.GotSwitch( "c" ) )
      {
        Alg3pt = Common::ArrayFromString<Common::Gamma::Algebra>(cl.SwitchValue<std::string>( "c" ), &Alg3ptNeg);
        if( Alg3pt.empty() )
          throw std::invalid_argument( "At least one current must be specified for three-point function" );
      }
      const bool b3pt{ !Alg3pt.empty() };
      BootstrapParams par;
      par.TimesliceDetail = cl.SwitchValue<int>( "t" );
      if( par.TimesliceDetail < 0 || par.TimesliceDetail > 2 )
        throw std::invalid_argument( "Timeslice detail " + std::to_string( par.TimesliceDetail ) + " invalid" );
      const std::string InStem{ cl.SwitchValue<std::string>( "i" ) };
      par.outStem = cl.SwitchValue<std::string>( "o" );
      par.nSample = cl.SwitchValue<int>( "n" );
      // Binning
      par.binSize = cl.SwitchValue<int>( "b" );
      par.binAuto = par.binSize == 0;
      if( par.binSize < 0 )
        throw std::runtime_error( "Bin size must be positive if specified" );
      par.binOrder = cl.SwitchValue<BinOrder>( "border" );
      if( par.binAuto && par.binOrder != BinOrder::Auto )
        throw std::runtime_error( "Auto binning only works with Auto order" );
      par.bFactorised = cl.GotSwitch( "f" );
      par.bVerboseSummaries = !cl.GotSwitch( "terse" );
      const std::string &DefaultGroup{ cl.SwitchValue<std::string>( "g" ) };
      const std::string &DefaultDataSet{ cl.SwitchValue<std::string>( "d" ) };
      bool bSwapQuarks{ cl.GotSwitch( "sort" ) };
      bool bShowOnly{ cl.GotSwitch( "show" ) };
      if( !b3pt )
        bSwapQuarks = !bSwapQuarks;
      GroupMomenta GroupP{ GroupMomenta::None };
      {
        const bool p2{ cl.GotSwitch( "p2" ) };
        if( cl.GotSwitch( "pa" ) )
        {
          if( p2 )
            throw std::invalid_argument( "Can't group momenta by both p^2 and abs( p )" );
          GroupP = GroupMomenta::Abs;
        }
        else if( p2 )
          GroupP = GroupMomenta::Squared;
      }
      if( cl.GotSwitch( "m" ) )
        par.MachineName = cl.SwitchValue<std::string>( "m" );
      else
        par.MachineName = MachineName;
      if( par.MachineName.empty() )
        throw std::invalid_argument( "Machine name can't be empty" );
      if( cl.GotSwitch( "r" ) )
      {
        bool bGotSeed{ true };
        try
        {
          par.seed = cl.SwitchValue<Common::SeedType>( "r" );
        }
        catch(const std::exception &e)
        {
          bGotSeed = false;
        }
        if( !bGotSeed )
        {
          par.RandomSample.Read( cl.SwitchValue<std::string>( "r" ) );
          if( !par.RandomSample.RandNum() )
            throw std::runtime_error( "No random numbers in " + cl.SwitchValue<std::string>( "r" ) );
          par.seed = par.RandomSample.Seed_;
        }
      }
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
        Manifest Manifest{ glob( cl.Args.begin(), cl.Args.end(), InStem.c_str() ),
                           cl.SwitchStrings( "x" ), bSwapQuarks, GroupP };
        // Walk the list of contractions, performing a separate bootstrap for each
        int BootstrapCount = 0;
        if( bShowOnly )
          std::cout << "Contraction, Files, Configs, Config..., Count..." << Common::NewLine;
        for( auto itc = Manifest.begin(); itc != Manifest.end(); itc++ )
        {
          const std::string &Contraction{itc->first};
          const TrajList &l{itc->second};
          const unsigned int nFile{ static_cast<unsigned int>( l.FileInfo.size() ) };
          if( bShowOnly )
          {
            // Count each configuration
            using CCMap = std::map<int, int>;
            CCMap ConfigCount;
            for( auto it = l.FileInfo.begin(); it != l.FileInfo.end(); ++it )
            {
              const int ConfigNum{ static_cast<int>( it->second.Config ) };
              auto p = ConfigCount.find( ConfigNum );
              if( p == ConfigCount.end() )
                ConfigCount.insert( { ConfigNum, 1 } );
              else
                p->second++;
            }
            //Print a summary of each configuration
            static const std::string SepTab{ ",\t" };
            std::cout << Contraction << SepTab << nFile << SepTab << ConfigCount.size();
            for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
              std::cout << Common::Comma << it->first;
            for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
              std::cout << Common::Comma << it->second;
            std::cout << Common::NewLine;
          }
          else
          {
            std::cout << "Loading " << nFile << " files for " << Contraction << Common::NewLine;
            const bool b3pt{ !Alg3pt.empty() };
            std::vector<Common::CorrelatorFileC> InFiles( nFile );
            typename std::map<std::string, TrajFile>::const_iterator it = l.FileInfo.begin();
            for( unsigned int j = 0; it != l.FileInfo.end(); ++j, ++it )
            {
              const std::string &Filename{ it->first }; // Trajectory number
              const TrajFile &tf{ it->second };
              std::cout << "  t=" << tf.Timeslice << ( ( tf.bHasTimeslice && tf.Timeslice ) ? "->0" : "   " )
                        << '\t' << Filename << std::endl;
              std::string GroupName{ DefaultGroup };
              InFiles[j].Read( Filename, b3pt ? Alg3pt : Alg, Alg,
                               tf.bHasTimeslice ? &tf.Timeslice : nullptr, nullptr, &GroupName,
                               b3pt && tf.bGotMomentum && tf.p.IsNeg() ? &Alg3ptNeg : nullptr, DefaultDataSet.c_str() );
            }
            try
            {
              par.PerformBootstrap( InFiles, l );
              BootstrapCount++;
            }
            catch(const std::exception &e)
            {
              std::cerr << "Error: " << e.what() << std::endl;
            }
          }
        }
        std::cout << "Bootstraps written for " << BootstrapCount
        << " / " << Manifest.size() << " correlators" << std::endl;
      }
      else if( bGotStudy )
      {
        const int iStudy{ cl.SwitchValue<int>( "s" ) };
        std::cout << "Study " << iStudy << ": ";
        const StudySubject Study{ static_cast<StudySubject>( iStudy ) };
        switch( Study )
        {
          case Z2:
          case GFWW:
          case GFPW:
          {
            // Nothing specified on the command-line
            // Perform default bootstraps
            // Heavy-light semi-leptonics
            const std::string qHeavy{ "h1" };
            const std::string qLight{ "l" };
            std::cout << "Making manifest for " << qHeavy << " and " << qLight << std::endl;
            static const Common::Momentum StudyMomenta[] = {
              { 0, 0, 0 },
              { 1, 0, 0 },
              { 1, 1, 0 },
              { 1, 1, 1 },
              { 2, 0, 0 },
            };
            const bool Factorised{ Study != GFPW };
            for( const Common::Momentum &m : StudyMomenta )
            {
              par.Study1Bootstrap( Study, InStem, m, Alg, qHeavy, qLight, Factorised );
              if( !Factorised )
                par.Study1Bootstrap( Study, InStem, m, Alg, qLight, qHeavy, false );
              par.Study1Bootstrap( Study, InStem, m, Alg, qHeavy, qHeavy, false );
              par.Study1Bootstrap( Study, InStem, m, Alg, qLight, qLight, false );
            }
          }
            bShowUsage = false;
            break;
          default:
            std::cout << " undefined" << std::endl;
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
    " <options> ContractionFile1 [ContractionFile2 ...]\n"
    "Perform a bootstrap of the specified files, where <options> are:\n"
    "-n     Number of samples (" DEF_NSAMPLE ")\n"
    "-b     Bin size, or 0 (default)=auto (1 config=no binning, else 1 bin/config)\n"
    "--border Bin Order: `Auto' (default)=config then timeslice then filename\n"
    "        `Old'=timeslice/filename/config, `VeryOld'=timeslice/config/filename\n"
    "-r     Random number seed (unspecified=random)\n"
    "-i     Input  prefix\n"
    "-o     Output prefix\n"
    "-a     list of gamma Algebras we're interested in at source (and sink for 2pt)\n"
    "-c     list of gamma algebras for current insertion         (Enable 3-pt mode)\n"
    "-g     Group name to read correlators from\n"
    "-d     DataSet name to read correlators from\n"
    "-s     Perform bootstrap for specified study number\n"
    "-t     timeslice detail 0 (none=default), 1 (.txt) or 2 (.txt+.h5)\n"
    "-m     Machine name (default: " << MachineName << ")\n"
    "-x     eXclude file (may be repeated)\n"
    "Flags:\n"
    "-f     Factorising operators (e.g. g5-gT5 same as gT5-g5)\n"
    "--p2   group momenta by P^2\n"
    "--pa   group momenta by Abs( p )\n"
    "--show Show how files would be bootstrapped, but don't execute\n"
    "--sort Disable(enable) sort correlator before group in 2pt(3pt) mode\n"
    "--terse Terse summaries (no file list)\n"
    "--help This message\n";
  }
  return iReturn;
}
