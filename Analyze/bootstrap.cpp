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

enum StudySubject { Z2=1, GFWW, GFPW };

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
  TrajFile( bool bHasTimeslice_, int  Timeslice_ ) : bHasTimeslice{ bHasTimeslice_ }, Timeslice{ Timeslice_ } {}
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
           bool bSwapQuarks, bool GroupLikeMomenta);
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
                   bool bSwapQuarks, bool GroupLikeMomenta)
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
        for( char &c : Contraction )
          if( c == '_' )
            c = ' ';
        std::vector<std::string> vs{ ArrayFromString<std::string>( Contraction ) };
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
        if( GroupLikeMomenta )
          sShortSuffix = mom.p2_string( Sep );
        else
        {
          sShortSuffix = "_p_";
          sShortSuffix.append( mom.to_string( Sep ) );
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
        cl.FileInfo.emplace( Name_.Filename, TrajFile( bHasTimeslice, Timeslice_ ) );
      else
        std::cout << "Ignoring repetition of " << Name_.Filename << std::endl;
    }
  }
  if( !parsed )
    throw std::runtime_error( "Remove non-existent files (or use '-x filename' to eXclude)" );
}

class BootstrapParams
{
public:
  std::string outStem;
  std::string MachineName;
  Common::SeedType seed;
  int nSample;
  int binSize;
  int TimesliceDetail;
  bool bFactorised;
  bool bFold;
  bool bT0Abs;
  bool bBackwardsWave;
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
  static const std::string sBlanks{ "  " };
  if( !bSaveSummaries && !bSaveBootstrap )
    return 0;
  const int NumFiles{ static_cast<int>( std::distance( first, last ) ) };
  if( NumFiles % binSize )
  {
    //throw std::runtime_error(  std::to_string(NumFiles) + " files not divisible by bin size " + std::to_string(binSize) );
    std::cout << "Warning: last bin partially filled (" << ( NumFiles % binSize ) << " correlators)\n";
  }
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
  const int NtHalf{ Nt / 2 + 1 };
  const int NtDest{ bFold ? NtHalf : Nt };
  const int FoldFactor{ bFold ? 2 : 1 };
  const std::string &sType{ bFold ? Common::sFold : Common::sBootstrap };
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
      const std::string sOutFile{ Common::MakeFilename( sOutBase, sType, seed, DEF_FMT ) };
      const std::string sSummary{ Common::MakeFilename( sOutBase, sType, seed, TEXT_EXT)};
      if( Common::FileExists( sOutFile ) || ( !bSaveBootstrap && Common::FileExists( sSummary ) ) )
      {
        std::cout << sBlanks << sOutBase << " skipped - output already exists" << std::endl;
      }
      else
      {
        // Copy the correlators into my sample, aligning timeslices if need be
        const bool bDifferent{ Snk != Src };
        const int OpFactor{ ( !b3pt && bFactorised && bDifferent ) ? 2 : 1 };
        Common::SampleC in( ( NumFiles * OpFactor * FoldFactor + binSize - 1 ) / binSize, NtDest );
        using Fold = Common::Fold<double>;
        Fold f;
        if( bFold )
        {
          // I haven't updated the folding code for 3 point functions ... folding should perhaps be deleted
          const Common::RPS att{DefaultOperatorAttribute(AlgSrc[Src])*DefaultOperatorAttribute(AlgSnk[Snk])};
          f.resize( nSample, NtDest );
          f.NtUnfolded = Nt;
          f.parity = att.parity;
          f.reality = att.reality;
          f.sign = att.sign;
          f.Conjugated = bFactorised;
          f.t0Negated = false;
        }
        int Bin{0};
        std::complex<double> * pDst = in[0];
        std::map<int, int> ConfigCount;
        std::vector<std::string> FileList;
        FileList.reserve( last - first );
        for( Iter it = first; it != last; )
        {
          const Common::CorrelatorFileC &file{ *it++ };
          const int CorrelatorTimeslice{ file.Timeslice() };
          const int TOffset{ bAlignTimeslices ? ( bBackwardsWave ? Nt - CorrelatorTimeslice : CorrelatorTimeslice ) : 0 };
          {
            // increment count
            assert( sizeof( int ) == sizeof( file.Name_.Seed ) );
            const int ConfigNum{ static_cast<int>( file.Name_.Seed ) };
            auto p = ConfigCount.find( ConfigNum );
            if( p == ConfigCount.end() )
              ConfigCount.insert( { ConfigNum, 1 } );
            else
              p->second++;
            // Save info about this file
            FileList.emplace_back( file.Name_.Filename );
          }
          // Only need to say which correlators contribute for first gamma structure
          if( Src == 0 && Snk == 0 )
          {
            std::cout << "  t=";
            if( TOffset )
              std::cout << TOffset << "->0";
            else
              std::cout << CorrelatorTimeslice;
            std::cout << '\t' << file.Name_.Base << "." << file.Name_.SeedString << std::endl;
          }
          for( int o = 0; o < OpFactor; o++ )
          {
            const std::complex<double> * const pSrc = file( AlgSnk[o ? Src : Snk], AlgSrc[o ? Snk : Src] );
            for( int h = 0; h < FoldFactor; h++ )
            {
              for( int t = 0; t < NtDest; t++ )
              {
                const std::ptrdiff_t idx{ ( ( h == 0 ? t : Nt - t ) + TOffset ) % Nt };
                std::complex<double> z{ pSrc[ idx ] };
                if( bFold && t != 0 )
                {
                  if( h && f.parity == Common::Parity::Odd )
                    z = -z;
                  if( f.sign == Common::Sign::Negative )
                    z = -z;
                }
                if( Bin )
                  pDst[t] += z;
                else
                  pDst[t]  = z;
              }
              // If we are binning, divide by number of items in this bin. Make sure to do the last bin
              if( ++Bin == binSize || ( o == OpFactor - 1 && h == FoldFactor - 1 && it == last ) )
              {
                if( binSize == 1 )
                  pDst += NtDest;
                else
                  for( int t = 0; t < NtDest; t++ )
                    *pDst++ /= binSize;
                Bin = 0;
              }
            }
          }
        }
        if( bSaveBootstrap )
          std::cout << sBlanks << nSample << " samples to " << sOutFile << std::endl;
        Common::SampleC out = in.Bootstrap( nSample, seed, &MachineName );
        // Now save the audit data for the bootstrap
        out.ConfigCount.reserve( ConfigCount.size() );
        for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
          out.ConfigCount.emplace_back( it->first, it->second );
        out.SampleSize = static_cast<int>( ( FileList.size() + binSize - 1 ) / binSize );
        out.binSize = binSize;
        out.FileList = std::move( FileList );
        if( bFold )
        {
          double * dst{ f[Fold::idxCentral] };
          const double * src{ SampleC::Traits::ScalarPtr( out[SampleC::idxCentral] )};
          if( f.reality == Common::Reality::Imag )
            src++;
          for( int s = SampleC::idxCentral; s < nSample; s++ )
          {
            // Timeslice 0 is a little special
            double d = *src++;
            src++;
            if( s == SampleC::idxCentral && bT0Abs && d < 0 )
              f.t0Negated = true;
            if( f.t0Negated )
              d = std::abs( d );
            *dst++ = d;
            // Now do all the other timeslices
            for( int t = 1; t < NtHalf; t++ )
            {
              *dst++ = *src++;
              src++;
            }
          }
          f.Seed_ = out.Seed_;
          f.SeedMachine_ = out.SeedMachine_;
          f.MakeCorrSummary( nullptr );
          f.SampleSize = out.SampleSize;
          f.binSize = out.binSize;
          f.FileList = std::move( out.FileList );
          f.ConfigCount = std::move( out.ConfigCount );
          if( bSaveBootstrap )
            f.Write( sOutFile );
          if( bSaveSummaries )
            f.WriteSummary( sSummary );
        }
        else
        {
          out.MakeCorrSummary( nullptr );
          if( bSaveBootstrap )
            out.Write( sOutFile );
          if( bSaveSummaries )
            out.WriteSummary( sSummary );
        }
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
  // sort all of the files by timeslice info
  using File = Common::CorrelatorFileC;
  std::sort( f.begin(), f.end(), [&](const File &l, const File &r)
  {
    // Sort first by timeslice
    if( l.bHasTimeslice != r.bHasTimeslice )
      return r.bHasTimeslice;
    if( l.bHasTimeslice && r.bHasTimeslice && l.Timeslice_ != r.Timeslice_ )
      return l.Timeslice_ < r.Timeslice_;
    // Sort by base filename
    {
      int iCompare = l.Name_.Base.compare( r.Name_.Base );
      if( iCompare )
        return iCompare < 0;
    }
    // Sort by trajectory
    //if( l.Name_.Seed != r.Name_.Seed )
      return l.Name_.Seed < r.Name_.Seed;
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
      {"b", CL::SwitchType::Single, "1"},
      {"r", CL::SwitchType::Single, nullptr},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"a", CL::SwitchType::Single, ""},
      {"c", CL::SwitchType::Single, ""},
      {"g", CL::SwitchType::Single, "" },
      {"s", CL::SwitchType::Single, nullptr},
      {"t", CL::SwitchType::Single, "0"},
      {"m", CL::SwitchType::Single, nullptr},
      {"x", CL::SwitchType::Multiple, nullptr},
      {"v", CL::SwitchType::Flag, nullptr},
      {"p", CL::SwitchType::Flag, nullptr},
      {"f", CL::SwitchType::Flag, nullptr},
      {"h", CL::SwitchType::Flag, nullptr},
      {"z", CL::SwitchType::Flag, nullptr},
      {"q", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) )
    {
      std::vector<Common::Gamma::Algebra> Alg = Common::ArrayFromString<Common::Gamma::Algebra>( cl.SwitchValue<std::string>( "a" ) );
      std::vector<Common::Gamma::Algebra> Alg3pt = Common::ArrayFromString<Common::Gamma::Algebra>( cl.SwitchValue<std::string>( "c" ) );
      BootstrapParams par;
      par.TimesliceDetail = cl.SwitchValue<int>( "t" );
      if( par.TimesliceDetail < 0 || par.TimesliceDetail > 2 )
        throw std::invalid_argument( "Timeslice detail " + std::to_string( par.TimesliceDetail ) + " invalid" );
      const std::string InStem{ cl.SwitchValue<std::string>( "i" ) };
      par.outStem = cl.SwitchValue<std::string>( "o" );
      par.nSample = cl.SwitchValue<int>( "n" );
      par.binSize = cl.SwitchValue<int>( "b" );
      par.bFactorised = cl.GotSwitch( "f" );
      par.bFold = cl.GotSwitch( "h" );
      par.bT0Abs = !cl.GotSwitch( "z" ); // Add the switch to turn this off
      par.bBackwardsWave = cl.GotSwitch( "b" );
      const bool bSwapQuarks{ !cl.GotSwitch( "q" ) }; // Add the switch to turn this off
      if( cl.GotSwitch( "m" ) )
        par.MachineName = cl.SwitchValue<std::string>( "m" );
      else
        par.MachineName = MachineName;
      if( par.MachineName.empty() )
        throw std::invalid_argument( "Machine name can't be empty" );
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
        Manifest Manifest{ glob( cl.Args.begin(), cl.Args.end(), InStem.c_str() ),
                           cl.SwitchStrings( "x" ), bSwapQuarks, cl.GotSwitch( "p" ) };
        // Walk the list of contractions, performing a separate bootstrap for each
        int BootstrapCount = 0;
        for( auto itc = Manifest.begin(); itc != Manifest.end(); itc++ )
        {
          const std::string &Contraction{itc->first};
          const TrajList &l{itc->second};
          const unsigned int nFile{ static_cast<unsigned int>( l.FileInfo.size() ) };
          const bool b3pt{ !l.Alg.empty() };
          std::cout << "Loading " << nFile << " files for " << Contraction << std::endl;
          unsigned int j = 0; // which trajectory are we processing
          std::vector<Common::CorrelatorFileC> InFiles( nFile );
          for( auto it = l.FileInfo.begin(); it != l.FileInfo.end(); it++, j++ )
          {
            const std::string &Filename{ it->first }; // Trajectory number
            const TrajFile &tf{ it->second };
            std::string GroupName{ cl.SwitchValue<std::string>( "g" ) };
            InFiles[j].Read( Filename, b3pt ? Alg3pt : Alg, Alg, tf.bHasTimeslice ? &tf.Timeslice : nullptr, nullptr, &GroupName );
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
    "-b     Bin size (default 1, i.e. no binning)\n"
    "-r     Random number seed (unspecified=random)\n"
    "-i     Input  prefix\n"
    "-o     Output prefix\n"
    "-a     list of gamma algebras we're interested in at source (and sink for 2pt)\n"
    "-c     list of gamma algebras for current insertion in 3pt functions\n"
    "-g     Group name to read correlators from\n"
    "-s     Perform bootstrap for specified study number\n"
    "-t     timeslice detail 0 (none=default), 1 (.txt) or 2 (.txt+.h5)\n"
    "-m     Machine name (default: " << MachineName << ")\n"
    "-x     eXclude file (may be repeated)\n"
    "Flags:"
    "-v     Backwards propagating wave (adjust timeslice by -t)"
    "-p     group momenta by P^2\n"
    "-f     Factorising operators (e.g. g5-gT5 same as gT5-g5)\n"
    "-h     Take forward/backward Half-waves as separate bootstrap inputs (i.e. fold)\n"
    "-z     Disable taking the absolute value of C(0) when folding\n"
    "-q     Disable Quark factorisation (default: sort correlator components & group)\n"
    "--help This message\n";
  }
  return iReturn;
}
