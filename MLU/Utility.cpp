/**

 Mike's lattice QCD utilities
 
 Source file: Utility.cpp
 
 Copyright (C) 2019 - 2023
 
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

// Common utilities (no dependencies other than c++ stdlib)

#include "Common.hpp"

#include <cstddef>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <sys/stat.h>

// This should be the only use of the AutoConf header
// Any useful values should be exposed as functions

#include <MLUconfig.h>
#include <MLU/JackBoot.hpp>

extern "C" const char * MLUVersionInfoHuman()
{
  static const char PackageString[] = MLU_PACKAGE_STRING ", " MLU_GIT_SUMMARY;
  return PackageString;
}

BEGIN_COMMON_NAMESPACE

// Text required for summaries of correlators
namespace CorrSumm {
  const char sep[] = " ";
  const char Comment[] = "# ";
  const char * FieldNames[NumFields] = { "corr", "exp", "cosh" };
};

const std::string sBootstrap{ "bootstrap" };
const std::string sFold{ "fold" };
const std::string sModel{ "model" };
const std::string sParams{ "params" };
const std::string sCovmat{ "covmat" };
const std::string sCovmatIn{ "covmat_in" };
const std::string sCovmatInv{ "covmat_inv" };
const std::string sCovmatInvCholesky{ sCovmatInv + "chol" };
const std::string sCormat{ "cormat" };
const std::string sCormatCholesky{ sCormat + "_chol" };
const std::string sCormatInvCholesky{ sCormat + "_invchol" };
const std::string sNtUnfolded{ "NtUnfolded" };
const std::string st0Negated{ "t0Negated" };
const std::string sConjugated{ "Conjugated" };
const std::string sTI{ "TI" };
const std::string sTF{ "TF" };
const std::string sDoF{ "DoF" };
const std::string sNumExponents{ "NumExponents" };
const std::string sNumFiles{ "NumFiles" };
const std::string sCovarFrozen{ "CovarFrozen" };
const std::string sStdErrorMean{ "StdErrorMean" };
const std::string sErrorScaled{ "ErrorScaled" };
const std::string sCovariance{ "Covariance" };
const std::string sCovarianceIn{ sCovariance + "In" };
const std::string sCovarianceInv{ sCovariance + "Inv" };
const std::string sCovarianceInvCholesky{ sCovarianceInv + "Cholesky" };
const std::string sCorrelation{ "Correlation" };
const std::string sCorrelationInv{ sCorrelation + "Inv" };
const std::string sCorrelationCholesky{ sCorrelation + "Cholesky" };
const std::string sCorrelationInvCholesky{ sCorrelation + "InvCholesky" };
const std::string sFitInput{ "FitInput" };
const std::string sModelPrediction{ "ModelPrediction" };
const std::string sOperators{ "Operators" };
const std::string sSummaryDSName{ "Summary" };
const std::string sSummaryNames{ "SummaryNames" };
const std::string sColumnNames{ "ColumnNames" };
const std::string sSampleSize{ "SampleSize" };
const std::string sEnsemble{ "Ensemble" };
const std::string sConfigCount{ "ConfigCount" };
const std::string sFileList{ "FileList" };
const std::string sBootstrapList{ "BootstrapList" };
const std::string sBinSize{ "BinSize" };
const std::string sCovarSource{ "CovarSource" };
const std::string sCovarRebin{ "CovarRebin" };
const std::string sCovarSampleSize{ "CovarSampleSize" };
const std::string sCovarNumBoot{ "CovarNumBoot" };
const std::string sGuess{ "Guess" };
const std::string sParam{ "Param" };
const std::string sFitTime{ "FitTime" };
const std::string sModelTypeList{ "ModelType" };
const std::string sModelArgsList{ "ModelArgs" };
const std::string sNE{ " != " };
const std::vector<std::string> sCorrSummaryNames{ "corr", "bias", "exp", "cosh" };
const std::string sChiSqPerDof{ "ChiSqPerDof" };
const std::string sPValue{ "pvalue" };
const std::string sPValueH{ "pvalueH" };

// String containing switch name to enable alternate overlap normalisation
// (see definition of ModelOverlap in Fit/ModelCommon.hpp)
const std::string sOverlapAltNorm{ "AltOver" };

const std::vector<std::string> DefaultModelStats{ Common::sChiSqPerDof, Common::sPValue, Common::sPValueH };

void MomentumMap::Parse( std::string &BaseShort )
{
  // Now get the remaining momenta and save them
  std::smatch match;
  for( int pLoop = 0; pLoop < 2; ++pLoop )
  {
    const std::regex Pattern{ Momentum::MakeRegex( pLoop ) };
    while( std::regex_search( BaseShort, match, Pattern ) )
    {
      // Extract the momentum
      const std::string sMom{ match[1] };
      std::vector<Momentum> np;
      np.reserve( 1 );
      if( pLoop == 0 )
        np.emplace_back( std::stoi( match[2] ), std::stoi( match[3] ), std::stoi( match[4] ) );
      else
        np.emplace_back( std::stoi( match[2] ) );
      const std::string sSuffix{ match.suffix() };
      BaseShort = match.prefix();
      BaseShort.append( sSuffix );
      // Now see whether we already have it
      auto it = find( sMom );
      if( it == end() )
        emplace( std::make_pair( sMom, np[0] ) );
      else if( it->second != np[0] )
      {
        std::stringstream ss;
        ss << "Repeated momentum " << sMom << CommaSpace << it->second << " != " << np[0];
        throw std::runtime_error( ss.str().c_str() );
      }
    }
  }
}

// These are the attributes I like to use in my filenames
void FileNameAtt::Parse( const std::string &Filename_, std::vector<std::string> * pOpNames,
                         const std::vector<std::string> * pIgnoreMomenta,
                         const std::vector<std::string> * pIgnoreRegEx,
                         bool bPreBootstrap )
{
  clear();
  Filename = Filename_;
  std::size_t pos = Filename.find_last_of( '/' );
  if( pos == std::string::npos )
    pos = 0;
  else
  {
    // Work out whether the last subdirectory is a Spectator directory
    if( pos )
    {
      std::size_t Prev = Filename.find_last_of( '/', pos - 1 );
      if( Prev == std::string::npos )
        Prev = 0;
      else
        ++Prev;
      // Spectator directories look like "3*_{spec}[p2]" ... with only one underscore
      if( Filename[Prev] == '3' )
      {
        std::string SpecDir{ Filename.substr( Prev, pos - Prev ) };
        Prev = SpecDir.find_first_of( '_', 1 );
        if( Prev != std::string::npos && Prev < SpecDir.length() - 1
           && SpecDir.find_first_of( '_', Prev + 1 ) == std::string::npos )
        {
          this->SpecDir = SpecDir;
          Spectator = SpecDir.substr( Prev + 1 );
          std::size_t SpecLen{ Spectator.length() };
          std::size_t SuffixLen{ Common::Momentum::DefaultPrefixSquared.length() };
          if( SpecLen > SuffixLen
              && Common::EqualIgnoreCase( Spectator.substr( SpecLen - SuffixLen ),
                                          Common::Momentum::DefaultPrefixSquared ) )
          {
            Spectator.resize( SpecLen - SuffixLen );
            bSpectatorGotSuffix = true;
          }
        }
      }
    }
    // Save the directory
    Dir = Filename.substr( 0, ++pos );
  }
  Base = Filename.substr( pos );
  NameNoExt = Base;
  int i = 0;
  while( i < ( bPreBootstrap ? 2 : 3 ) && ( pos = Base.find_last_of( '.' ) ) != std::string::npos )
  {
    switch( i )
    {
      case 0:
        Ext = Base.substr( pos + 1 );
        break;
      case 1:
        SeedString = Base.substr( pos + 1 );
        try {
          Seed = RandomCache::Seed( SeedString );
          bSeedNum = true;
        } catch(...) {
          std::cout << "Ignoring invalid seed in " << Filename << std::endl;
        }
        break;
      case 2:
        Type = Base.substr( pos + 1 );
        break;
    }
    Base.resize(pos);
    if( i == 0 )
      NameNoExt = Base;
    i++;
  }
  // If there are extra segments, the last contains operator names. Otherwise extract 2 op names
  if( !bPreBootstrap )
  {
    pos = Base.find_last_of( '.' );
    opNames = ParseOpNames( pos == std::string::npos ? 2 : INT_MAX );
    op.resize( opNames.size() );
    std::vector<std::string> &OpNames{ pOpNames ? * pOpNames : opNames };
    for( int i = 0; i < opNames.size(); ++i )
    {
      // We want to be able to compare the indices to see whether source and sink same
      // ... even when a global operator list not used
      op[i] = 0;
      while( op[i] < OpNames.size() && !EqualIgnoreCase( opNames[i], OpNames[op[i]] ) )
        op[i]++;
      if( op[i] == OpNames.size() )
        OpNames.emplace_back( opNames[i] );
    }
  }
  // Save any extra segments after '.'. I parse backward, so [0] is last, [1] is second last, etc
  while( ( pos = Base.find_last_of( '.' ) ) != std::string::npos )
  {
    Extra.push_back( Base.substr( pos + 1 ) );
    Base.resize( pos );
  }

  // Remove momenta we are going to ignore in the filename
  if( pIgnoreMomenta )
  {
    Momentum pIgnore;
    for( const std::string &s : *pIgnoreMomenta )
    {
      while( pIgnore.Extract( Base, s ) ) // There might be more than one copy
        ;
    }
  }
  if( pIgnoreRegEx )
  {
    for( const std::string &s : *pIgnoreRegEx )
    {
      // Remove strings from the filename. No need to expose these strings (for now)
      while( ExtractToken( Base, s ) )
        ;
    }
  }
  std::string BaseShort = Base;
  p.Parse( BaseShort );
  // Extract other attributes from filename
  ExtractTimeslice( BaseShort, bGotTimeslice, Timeslice );
  ExtractDeltaT( BaseShort, bGotDeltaT, DeltaT );
  Gamma.clear();
  Gamma = ExtractGamma( BaseShort );
  BaseShortParts = ArrayFromString( BaseShort, Underscore );
  if( !Spectator.empty() && BaseShortParts.size() >= 3 )
  {
    Quark.reserve( 2 );
    Quark.push_back( BaseShortParts[2] );
    Quark.push_back( BaseShortParts[1] );
    Meson.emplace_back( MakeMesonName( Quark[0] ) );
    Meson.emplace_back( MakeMesonName( Quark[1] ) );
    auto pp = p.find( Momentum::DefaultPrefix );
    if( pp != p.end() )
    {
      auto pps = p.find( Momentum::SinkPrefix );
      if( pps != p.end() )
      {
        // We have two momenta
        MesonP.emplace_back( *pp );
        MesonP.emplace_back( *pps );
      }
      else
      {
        MesonP.resize( 2 ); // They will both default to zero momentum
        // We only have one momentum - this goes with the lightest meson
        const bool bSourceHeavier{ QuarkWeight(BaseShortParts[2]) >= QuarkWeight(BaseShortParts[1]) };
        MesonP[bSourceHeavier ? 1 : 0].p = pp->second;
        MesonP[bSourceHeavier ? 0 : 1].p.bp2 = pp->second.bp2;
      }
    }
  }
  else if( Spectator.empty() && BaseShortParts.size() >= 2 )
  {
    Quark.reserve( 2 );
    Quark.push_back( BaseShortParts[1] ); // Source
    Quark.push_back( BaseShortParts[0] ); // Sink
    Meson.emplace_back( ::Common::ConcatHeavyFirst( Quark[1], Quark[0] ) );
    typename MomentumMap::iterator it{ p.find( Momentum::DefaultPrefix ) };
    if( it != p.end() )
      MesonP.emplace_back( *it );
  }
  // If we know momenta, create a version of the meson names with momenta and DefaultPrefix
  for( std::size_t i = 0; i < MesonP.size(); ++i )
    MesonMom.push_back( Meson[i] + MesonP[i].p.FileString( Momentum::DefaultPrefix ) );
}

// Parse operator names from the end of Base, building up a list of all the operator names
// NB: Because I parse from the end, op[0] is the last operator, op[1] second last, etc
std::vector<std::string> FileNameAtt::ParseOpNames( int NumOps )
{
  static const char Sep[] = "_.";
  constexpr std::size_t NumSeps{ sizeof( Sep ) / sizeof( Sep[0] ) - 1 };
  std::vector<std::string> o;
  std::size_t LastPos = Base.length();
  int i = 0;
  for( ; LastPos && i < NumOps; ++i )
  {
    // Search for the separators - on the last one, check for a period as well as underscore
    LastPos--; // Skip past the last separator
    std::size_t pos = Base.find_last_of( Sep, LastPos, NumSeps );
    if( pos == std::string::npos )
    {
      o.push_back( Base.substr( 0, LastPos + 1 ) );
      LastPos = 0;
    }
    else
    {
      o.push_back( Base.substr( pos + 1, LastPos - pos ) );
      LastPos = pos;
      if( Base[pos] == '.' )
        i = NumOps - 1;
    }
  }
  if( i != NumOps )
    throw std::runtime_error( "Can't extract " + std::to_string( NumOps ) + " operators from " + Base );
  if( i )
  {
    Base.resize( LastPos ); // Shorten the base
  }
  return o;
}

// Append the extra info to the string
void FileNameAtt::AppendExtra( std::string &s, int Last, int First ) const
{
  if( Extra.size() )
  {
    if( Last < 0 )
      Last = 0;
    if( First < 0 )
      First = static_cast<int>( Extra.size() - 1 + First );
    if( First < 0 || First >= Extra.size() )
      First = static_cast<int>( Extra.size() - 1 );
    for( int i = First; i >= Last; --i )
    {
      s.append( 1, '.' );
      s.append( Extra[i] );
    }
  }
}

void FileNameAtt::AppendOps( std::string &s, const std::string &FirstSeparator,
                             std::vector<std::string> * pOpNames ) const
{
  const std::vector<std::string> &OpNames{ pOpNames ? * pOpNames : opNames };
  bool bFirst{ true };
  for( std::size_t i = op.size(); i--; )
  {
    s.append( bFirst ? FirstSeparator : Underscore );
    bFirst = false;
    s.append( OpNames[op[i]] );
  }
}

std::string FileNameAtt::GetBaseShort( int First, int Last ) const
{
  std::string s;
  const int iSize{ static_cast<int>( BaseShortParts.size() ) };
  if( iSize )
  {
    // Negative inputs are counts from the end
    if( First < 0 )
      First += iSize;
    if( Last < 0 )
      Last += iSize;
    // Now we range check
    if( First < 0 )
      First = 0;
    if( Last >= iSize )
      Last = iSize - 1;
    while( First <= Last )
    {
      if( !s.empty() )
        s.append( 1, '_' );
      s.append( BaseShortParts[First++] );
    }
  }
  return s;
}

std::string FileNameAtt::GetBaseShortExtra( int First, int Last ) const
{
  std::string s{ GetBaseShort( First, Last ) };
  AppendExtra( s );
  return s;
}

std::string FileNameAtt::GetBaseExtra( int Last, int First ) const
{
  std::string s{ Base };
  AppendExtra( s, Last, First );
  return s;
}

// Make a new name based on this one, overriding specified elements
std::string FileNameAtt::DerivedName( const std::string &Suffix, const std::string &Snk, const std::string &Src,
                                      const std::string &Ext ) const
{
  std::string s{ Base };
  s.append( 1, '_' );
  s.append( Snk );
  s.append( 1, '_' );
  s.append( Src );
  AppendExtra( s );
  s.append( 1, '.' );
  s.append( Type );
  s.append( 1, '.' );
  s.append( SeedString );
  s.append( 1, '.' );
  s.append( Ext );
  return s;
}

const MomentumPair &FileNameAtt::GetMomentum( const std::string &Name ) const
{
  auto it = p.find( Name );
  if( it == p.end() )
    it = p.find( Name + Momentum::SquaredSuffix );
  if( it == p.end() )
    throw std::runtime_error( "Momentum " + Name + " not available" );
  return *it;
}

bool FileNameAtt::HasNonZeroMomentum() const
{
  for( auto it = p.begin(); it != p.end(); ++it )
    if( it->second )
      return true;
  return false;
}

const MomentumPair &FileNameAtt::GetFirstNonZeroMomentum() const
{
  if( p.empty() )
    throw std::runtime_error( "No momenta available" );
  for( auto it = p.begin(); it != p.end(); ++it )
    if( it->second )
      return *it;
  return *p.begin();
}

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext)
{
  const char Sep = '.';
  std::string s{ Base };
  s.append( 1, Sep );
  s.append( Type );
  s.append( 1, Sep );
  s.append( RandomCache::SeedString( Seed ) );
  s.append( 1, Sep );
  s.append( Ext );
  return s;
}

std::string MakeFilename( const std::string &Base, const std::string &Type, const std::string &Ext )
{
  return MakeFilename( Base, Type, RandomCache::DefaultSeed(), Ext );
}

std::string ExistingFilename( const std::string &Base, const std::string &Type, SeedType Seed,
                              const std::string &Ext )
{
  // Check the file exists
  std::string Filename{ MakeFilename( Base, Type, Seed, Ext ) };
  if( !FileExists( Filename ) )
  {
    // If it doesn't, see whether there's another name we can use
    std::string SearchName{ Base };
    SearchName.append( 1, '.' );
    SearchName.append( Type );
    SearchName.append( ".*." );
    SearchName.append( Ext );
    const std::string *Begin{ &SearchName };
    std::vector<std::string> v{ glob( Begin, Begin + 1 ) };
    if( !v.empty() )
      Filename = v[0];
  }
  return Filename;
}

/// Find existing file with specified Seed, but fall back to another seed if missing
std::string ExistingAnySeed( const std::string &sFilename )
{
  if( !FileExists( sFilename ) )
  {
    std::string New{ sFilename };
    bool bBad{};
    std::array<std::string, 2> Ext;
    for( std::size_t i = 0; !bBad && i < Ext.size(); ++i )
    {
      std::size_t pos{ New.find_last_of( "/." ) };
      if( pos == std::string::npos || New[pos] != '.' )
        bBad = true;
      else
      {
        Ext[i] = New.substr( pos + 1 );
        New.resize( pos );
      }
    }
    if( !bBad )
    {
      New.append( ".*." );
      New.append( Ext[0] );
      const std::string *Begin{ &New };
      std::vector<std::string> v{ glob( Begin, Begin + 1 ) };
      if( !v.empty() )
        return v[0];
    }
  }
  return sFilename;
}

/// Override Seed in filename to specified seed (if it exists)
std::string PreferSeed( const std::string &sFilename, SeedType NewSeed )
{
  std::string New{ sFilename };
  std::array<std::string, 2> Ext;
  for( std::size_t i = 0; i < Ext.size(); ++i )
  {
    std::size_t pos{ New.find_last_of( "/." ) };
    if( pos == std::string::npos || New[pos] != '.' )
      return sFilename;
    Ext[i] = New.substr( pos + 1 );
    New.resize( pos );
  }
  New.append( 1, '.' );
  New.append( RandomCache::SeedString( NewSeed ) );
  New.append( 1, '.' );
  New.append( Ext[0] );
  if( FileExists( New ) )
    return New;
  return sFilename;
}

std::string PreferSeed( const std::string &sFilename )
{
  return PreferSeed( sFilename, RandomCache::DefaultSeed() );
}

// If present, remove Token from a string. Return true if removed
bool ExtractToken( std::string &Prefix, const std::string &Token )
{
  bool bExtracted{ false };
  std::smatch match;
  const std::regex pattern{ "(^|_)" + Token + "(_|$)" };
  while( std::regex_search( Prefix, match, pattern  ) )
  {
    if( bExtracted )
      throw std::runtime_error( "Multiple " + Token + " tokens in " + Prefix );
    bExtracted = true;
    std::string s{ match.prefix() };
    if( match[1].length() && match[2].length() )
      s.append( match[2] );
    s.append( match.suffix() );
    Prefix = s;
  }
  return bExtracted;
}

// If present, remove integer preceded by Token from a string
void ExtractInteger( std::string &Prefix, bool &bHasValue, int &Value, const std::string Token )
{
  std::smatch match;
  const std::regex pattern{ "_" + Token + "_?([0-9]+)" };
  while( std::regex_search( Prefix, match, pattern  ) )
  {
    int ThisValue = std::stoi( match[1] );
    if( !bHasValue )
    {
      Value = ThisValue;
      bHasValue = true;
    }
    else if( ThisValue != Value )
      throw std::runtime_error( "Multiple regex " + Token + ": " + std::to_string(Value) + " and " + std::to_string(ThisValue) );
    const std::string sSuffix{ match.suffix() };
    Prefix = match.prefix();
    Prefix.append( sSuffix );
  }
}

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &Prefix, bool &bHasTimeslice, int &Timeslice )
{
  ExtractInteger( Prefix, bHasTimeslice, Timeslice, "[tT]" );
}

// Strip out DeltaT from a string if present
void ExtractDeltaT( std::string &Prefix, bool &bHasDeltaT, int &DeltaT )
{
  ExtractInteger( Prefix, bHasDeltaT, DeltaT, "[dD][tT]" );
}

// Append DeltaT to string
void AppendDeltaT( std::string &s, int DeltaT )
{
  s.append( "_dt_" );
  s.append( std::to_string( DeltaT ) );
}

// Remove any gammas from Prefix
std::vector<Gamma::Algebra> ExtractGamma( std::string &Prefix )
{
  static const std::regex GammaPattern{ "_(g[^_]+)" };

  std::vector<Gamma::Algebra> v;
  std::smatch match;
  std::string Search{ std::move( Prefix ) };
  Prefix.clear();
  while( std::regex_search( Search, match, GammaPattern ) )
  {
    Prefix.append( match.prefix() );
    Gamma::Algebra a;
    std::stringstream ss( match[1] );
    if( ss >> a && ( ss.eof() || ( ss >> std::ws && ss.eof() ) ) )
      v.push_back( a );
    else
    {
      Prefix.append( 1, '_' );
      Prefix.append( match[1] );
    }
    Search = match.suffix();
  }
  Prefix.append( Search );
  return v;
}

// Append Gamma to string
void AppendGamma( std::string &s, Gamma::Algebra g, const std::string &Sep )
{
  std::ostringstream os;
  os << Sep << g;
  s.append( os.str() );
}

// Remove any gammas from Prefix
void ReplaceGamma( std::string &Prefix, Gamma::Algebra gFrom, Gamma::Algebra gTo )
{
  std::regex GammaPattern;
  {
    std::ostringstream ss;
    ss << Common::Underscore << gFrom << "(_|$)";
    GammaPattern = ss.str();
  }
  std::string GammaReplace;
  {
    std::ostringstream ss;
    ss << Common::Underscore << gTo;
    GammaReplace = ss.str();
  }
  std::smatch match;
  std::string Search{ std::move( Prefix ) };
  Prefix.clear();
  while( std::regex_search( Search, match, GammaPattern ) )
  {
    Prefix.append( match.prefix() );
    Prefix.append( GammaReplace );
    Prefix.append( match[1] );
    Search = match.suffix();
  }
  Prefix.append( Search );
}

const std::string sUnknown{ "unknown" };

const std::string sReality{ "Reality" };
const std::string Reality_TextEquiv_Real{ "real" };
const std::string Reality_TextEquiv_Imaginary{ "imaginary" };

std::ostream& operator<<(std::ostream& os, const Reality &reality)
{
  switch( reality )
  {
    case Reality::Real:
      os << Reality_TextEquiv_Real;
      break;
    case Reality::Imag:
      os << Reality_TextEquiv_Imaginary;
      break;
    default:
      os << sUnknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Reality &reality)
{
  reality = Reality::Unknown;
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Reality_TextEquiv_Real ) )
      reality = Reality::Real;
    else if( EqualIgnoreCase( s, Reality_TextEquiv_Imaginary ) )
      reality = Reality::Imag;
    else if( !EqualIgnoreCase( s, sUnknown ) )
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

const std::string sParity{ "Parity" };
const std::string Parity_TextEquiv_Even{ "even" };
const std::string Parity_TextEquiv_Odd{ "odd" };

std::ostream& operator<<(std::ostream& os, const Parity &parity)
{
  switch( parity )
  {
    case Parity::Even:
      os << Parity_TextEquiv_Even;
      break;
    case Parity::Odd:
      os << Parity_TextEquiv_Odd;
      break;
    default:
      os << sUnknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Parity &parity)
{
  parity = Parity::Unknown;
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Parity_TextEquiv_Even ) )
      parity = Parity::Even;
    else if( EqualIgnoreCase( s, Parity_TextEquiv_Odd ) )
      parity = Parity::Odd;
    else if( !EqualIgnoreCase( s, sUnknown ) )
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

const std::string sSign{ "Sign" };
const std::string Sign_TextEquiv_Positive{ "positive" };
const std::string Sign_TextEquiv_Negative{ "negative" };

std::ostream& operator<<(std::ostream& os, const Sign &sign)
{
  switch( sign )
  {
    case Sign::Positive:
      os << Sign_TextEquiv_Positive;
      break;
    case Sign::Negative:
      os << Sign_TextEquiv_Negative;
      break;
    default:
      os << sUnknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Sign &sign)
{
  sign = Sign::Unknown;
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Sign_TextEquiv_Positive ) )
      sign = Sign::Positive;
    else if( EqualIgnoreCase( s, Sign_TextEquiv_Negative ) )
      sign = Sign::Negative;
    else if( !EqualIgnoreCase( s, sUnknown ) )
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

/// TODO: Far from the best implementation
template <typename T>
void SummaryHeader( std::ostream &s, const std::string & sOutFileName,
                    const char * pszHeaderComments )
{
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  if( !s )
    throw std::runtime_error( "Unable to create " + sOutFileName );
  s << "# Summary: " << sOutFileName << NewLine
    << std::boolalpha << std::setprecision(std::numeric_limits<scalar_type>::digits10+2);
  const char * psz = nullptr;
  if( Traits::is_complex )
  {
    switch( sizeof( scalar_type ) )
    {
      case sizeof( float ):
        psz = "std::complex<float>";
        break;
      case sizeof( double ):
        psz = "std::complex<double>";
        break;
    }
  }
  else
  {
    switch( sizeof( scalar_type ) )
    {
      case sizeof( float ):
        psz = "float";
        break;
      case sizeof( double ):
        psz = "double";
        break;
    }
  }
  if( psz ) s << "# Representation: " << psz << NewLine;
  if( pszHeaderComments ) s << pszHeaderComments;
}

template void SummaryHeader<float>( std::ostream &s,
                            const std::string &sOutFileName, const char * pszHeaderComments );
template void SummaryHeader<double>( std::ostream &s,
                            const std::string &sOutFileName, const char * pszHeaderComments );
template void SummaryHeader<std::complex<float>>( std::ostream &s,
                            const std::string &sOutFileName, const char * pszHeaderComments );
template void SummaryHeader<std::complex<double>>( std::ostream &s,
                            const std::string &sOutFileName, const char * pszHeaderComments );

template <typename T>
void SummaryHelper( const std::string & sOutFileName, const T * pData, const int nt,
                    const int nSample, const char * pszHeaderComments, bool bFolded )
{
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  using namespace CorrSumm;
  std::ofstream s( sOutFileName );
  SummaryHeader<T>( s, sOutFileName, pszHeaderComments );
  // Now write all the field names
  s << "t";
  for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
  {
    for(int f = 0; f < NumFields; f++) // each field
    {
      std::string n{ FieldNames[f] };
      if( i )
        n.append( "_im" );
      s << sep;
      Common::ValWithEr<T>::Header( n, s, sep );
    }
  }
  s << NewLine;
  // Now perform summaries
  const int tMid{ bFolded ? nt : nt / 2 }; // Summaries may want to change sign for -ve wave
  scalar_type Central = -777; // Just to silence compiler warning
  std::vector<scalar_type> Data( nSample );
  for(int t = 0; t < nt; t++ )
  {
    // Index of previous and next timeslice relative to current timeslice
    const int Next{ ( t == nt - 1 ? 1 - nt :  1 ) * Traits::scalar_count };
    const int Prev{ ( t == 0      ? nt - 1 : -1 ) * Traits::scalar_count };
    s << t;
    for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
    {
      for(int f = 0; f < NumFields; f++) // each field
      {
        const scalar_type * p = Traits::ScalarPtr( pData ) + t * Traits::scalar_count + i;
        int Count = 0;
        for(int n = -1; n < nSample; n++, p += nt * Traits::scalar_count )
        {
          scalar_type d;
          switch( f )
          {
            case 0:
              d = * p;
              break;
            case 1: // mass
              if( t == 0 )
                d = NaN;
              else if( t <= tMid )
                d = std::log( p[Prev] / *p );
              else
                d = -std::log( *p / p[Next] );
              break;
            case 2: // cosh mass
              d = std::acosh( ( p[Prev] + p[Next] ) / ( *p * 2 ) );
              break;
            default:
              d = 0;
          }
          if( n == -1 )
            Central = d;
          else if( std::isfinite( d ) )
            Data[Count++] = d;
        }
        Common::ValWithEr<scalar_type> v( Central, Data, Count );
        s << sep << v;
      }
    }
    s << NewLine;
  }
}

template void SummaryHelper<float>
 ( const std::string &sOutFileName, const float *pData,
   const int nt, const int nSample, const char * pszHeaderComments, bool bFolded );
template void SummaryHelper<double>
 ( const std::string &sOutFileName, const double *pData,
   const int nt, const int nSample, const char * pszHeaderComments, bool bFolded );
template void SummaryHelper<std::complex<float>>
 ( const std::string &sOutFileName, const std::complex<float> *pData,
   const int nt, const int nSample, const char * pszHeaderComments, bool bFolded );
template void SummaryHelper<std::complex<double>>
 ( const std::string &sOutFileName, const std::complex<double> *pData,
   const int nt, const int nSample, const char * pszHeaderComments, bool bFolded );

const std::array<std::string, 3> aSampleSource{ "Binned", "Raw", "Bootstrap" };

std::ostream& operator<<( std::ostream& os, SampleSource sampleSource )
{
  const int enumIdx{ static_cast<int>( sampleSource ) };
  if( enumIdx >= 0 && enumIdx < aSampleSource.size() )
    return os << aSampleSource[enumIdx];
  return os << "Unknown" << std::to_string( enumIdx );
}

std::istream& operator>>( std::istream& is, SampleSource &sampleSource )
{
  std::string sEnum;
  if( is >> sEnum )
  {
    const int enumIdx{ IndexIgnoreCase( aSampleSource, sEnum ) };
    if( enumIdx != aSampleSource.size() )
    {
      sampleSource = static_cast<SampleSource>( enumIdx );
      return is;
    }
  }
  throw std::runtime_error( "SampleSource \"" + sEnum + "\" unrecognised" );
}

// Read an array (real or complex) from an HDF5 file
#define template_ReadArrayArgs( T ) \
std::vector<T> &buffer, const std::string &FileName, const std::string &ObjectName, \
const char *PrintPrefix, std::string * pGroupName

template<typename T> void ReadArray( template_ReadArrayArgs( T ) )
{
  ::H5::H5File f;
  ::H5::Group  g;
  H5::OpenFileGroup( f, g, FileName, PrintPrefix, pGroupName );
  ::H5::DataSet ds = g.openDataSet(ObjectName);
  ::H5::DataSpace dsp = ds.getSpace();
  const int nDims{dsp.getSimpleExtentNdims()};
  if( nDims != 1 )
    throw std::runtime_error("Object " + ObjectName + " in " + FileName + " has " + std::to_string( nDims ) + " dimensions" );
  hsize_t Nt;
  dsp.getSimpleExtentDims( &Nt );
  hsize_t BufferSize{ static_cast<hsize_t>( buffer.size() ) };
  if( BufferSize == 0 )
    buffer.resize( Nt );
  else if( BufferSize != Nt )
    throw std::runtime_error("Object " + ObjectName + " in " + FileName + " has " + std::to_string( Nt ) + " entries, doesn't match Nt=" + std::to_string( BufferSize ) );
  ds.read( &buffer[0], H5::Equiv<T>::Type );
  if( !IsFinite( buffer ) )
     throw std::runtime_error( "Error: Infinite/NaN values in " + FileName );
}

// Make sure all the specialisations I support are instantiated
template void ReadArray<float>( template_ReadArrayArgs( float ) );
template void ReadArray<double>( template_ReadArrayArgs( double ) );
template void ReadArray<std::complex<float>>( template_ReadArrayArgs( std::complex<float> ) );
template void ReadArray<std::complex<double>>( template_ReadArrayArgs( std::complex<double> ) );

std::ostream& operator<<( std::ostream& os, const CommandLine &cl)
{
  os << "Command-line has " << cl.Args.size() << " arguments and " << cl.Switches.size() << " switches";
  for( int i = 0; i < cl.Args.size(); i++ )
    os << "\nArg[" << i << "]=\"" << cl.Args[i] << "\"";
  int i = 0;
  //for( const CommandLine::SwitchPair &p : cl.Switches ) {
  for( CommandLine::SwitchMap::const_iterator p = cl.Switches.begin(); p != cl.Switches.end(); ++p ) {
    const std::string &Switch{ p->first };
    const std::vector<std::string> &Values{ p->second };
    os << "\nSwitch[" << i << "]=\"" << Switch << "\", was specified";
    const std::size_t sz{ Values.size() };
    if( sz ) {
      os << " " << Values.size() << " time";
      if( sz > 1 )
        os << "s";
    }
    for( int j = 0; j < Values.size(); j++ )
      os << "\n      [" << i << "][" << j << "]=\"" << Values[j] << "\"";
    ++i;
  }
  return os;
}

bool CommandLine::IsValuePresent( const char * & p )
{
  const char * const pOriginal{ p };
  while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
    p++;
  if( *p == '=' ) {
    ++p;
    while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
      p++;
  }
  // There's a value present if there are characters in the string
  // ... or if we skipped past whitespace and/or equal sign to get to end of string
  return *p || pOriginal != p;
}

void CommandLine::Parse( int argc, const char *argv[], const std::vector<SwitchDef> &defs )
{
  // Save the name of the executable without any path
  Name = argv[0];
  auto pos = Name.find_last_of('/');
  if( pos != std::string::npos )
    Name = Name.substr( pos + 1 );
  // Now parse the command-line
  bool bError{ false };
  int SwitchNo = -1; // Not waiting for switch value
  bool bInMultiChar{ false };
  for( int i = 1; i < argc; i++ ) {
    const char *p = argv[i];
    if( *p == '-' ) { // Switch
      if( SwitchNo != -1 ) { // We were waiting on a switch that requires a value
        Switches[defs[SwitchNo].Switch].push_back( "" );
        SwitchNo = -1;
      }
      std::string SwitchName;
      if( *++p == '-' ) {
        // multi-char switch.
        const char * pEnd = ++p;
        if( *p == 0 ) {
          // This is the special switch "--"
          p -= 2;
          bInMultiChar = true;
        } else {
          while( *pEnd && *pEnd != '=' && *pEnd != ' ' && *pEnd != '\t' && *pEnd != '\r' && *pEnd != '\n' )
            pEnd++;
        }
        SwitchName = std::string( p, pEnd - p );
        p = pEnd;
      } else {
        // Single character switch
        if( *p == 0 ) {
          std::cerr << "Unrecognised switch \"-\"" << std::endl;
          bError = true;
          continue;
        }
        SwitchName.resize( 1 );
        SwitchName[0] = *p++;
      }
      // I've just decoded a short or long switch name - see whether it's valid
      SwitchNo = 0;
      while( SwitchNo < defs.size() && SwitchName.compare( defs[SwitchNo].Switch ) )
        SwitchNo++;
      if( SwitchNo >= defs.size() ) {
        std::cerr << "Unrecognised switch \"" << argv[i] << "\"" << std::endl;
        bError = true;
        SwitchNo = -1;
        continue;
      }
      if( defs[SwitchNo].Type == SwitchType::Flag ) {
        // This is just a switch - it should not have a value
        // Swallow any trailing whitespace
        while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
          p++;
        if( *p ) {
          std::cerr << "Switch \"" << SwitchName << "\" should not have a value \"" << p << "\"" << std::endl;
          bError = true;
        } else if( GotSwitch( SwitchName ) ) {
          std::cerr << "Switch \"" << SwitchName << "\" should not be repeated" << std::endl;
          bError = true;
        } else
          Switches[SwitchName]; // First time we've seen this flag
        SwitchNo = -1;
        continue;
      } else if( defs[SwitchNo].Type == SwitchType::Single && GotSwitch( SwitchName ) ) {
        std::cerr << "Switch \"" << SwitchName << "\" should not be repeated" << std::endl;
        bError = true;
        SwitchNo = -1;
        continue;
      }
      // Use the remainder of this switch as a value ... or wait for next param if empty
      if( !IsValuePresent( p ) )
        continue;
    }
    if( SwitchNo != -1 ) {
      // Get the value of the current switch
      Switches[defs[SwitchNo].Switch].push_back( p );
      if( !bInMultiChar )
        SwitchNo = -1;
    }
    else // Argument
      Args.push_back( p );
  }
  if( SwitchNo != -1 && !bInMultiChar ) {
    std::cerr << "Still processing switch \"" << defs[SwitchNo].Switch << "\"" << std::endl;
    bError = true;
  }
  if( bError ) {
    // std::cerr << (*this) << std::endl;
    throw std::runtime_error( "Invalid command-line" );
  }
  // Now put defaults in for any missing switches
  for( int i = 0; i < defs.size(); i++ ) {
    if( defs[i].Type == Single && defs[i].Default && !GotSwitch( defs[i].Switch ) )
      Switches[defs[i].Switch].push_back( defs[i].Default );
  }
}
END_COMMON_NAMESPACE
