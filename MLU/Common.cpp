/**

 Mike's lattice QCD utilities
 
 Source file: Common.cpp
 
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
    Meson.emplace_back( MakeMesonName( BaseShortParts[2] ) );
    Meson.emplace_back( MakeMesonName( BaseShortParts[1] ) );
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
      MesonMom = Meson;
      for( std::size_t i = 0; i < MesonMom.size(); ++i )
        MesonMom[i].append( MesonP[i].p.FileString( Momentum::DefaultPrefix ) );
    }
  }
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

template <typename T>
void CorrelatorFile<T>::resize( int NumSnk, int NumSrc, int Nt )
{
  const std::size_t OldMemSize{ static_cast<std::size_t>( NumSnk_ ) * NumSrc_ * Nt_ };
  const std::size_t NewMemSize{ static_cast<std::size_t>( NumSnk  ) * NumSrc  * Nt  };
  if( NumSnk_ != NumSnk || NumSrc_ != NumSrc )
  {
    AlgSnk_.clear();
    AlgSrc_.clear();
    NumSnk_ = NumSnk;
    NumSrc_ = NumSrc;
  }
  Nt_ = Nt;
  if( OldMemSize != NewMemSize )
  {
    if( NewMemSize == 0 )
      m_pData.reset( nullptr );
    else
      m_pData.reset( new T[ NewMemSize ] );
  }
}

template <typename T>
void CorrelatorFile<T>::swap( CorrelatorFile &o )
{
  using std::swap;
  swap( NumSnk_, o.NumSnk_ );
  swap( NumSrc_, o.NumSrc_ );
  swap( Nt_, o.Nt_ );
  swap( m_pData, o.m_pData );
  swap( AlgSnk_, o.AlgSnk_ );
  swap( AlgSrc_, o.AlgSrc_ );
  Name_.swap( o.Name_ );
  swap( bHasTimeslice, o.bHasTimeslice );
  swap( Timeslice_, o.Timeslice_ );
}

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void CorrelatorFile<T>::Read(const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk,
                             std::vector<Gamma::Algebra> &AlgSrc, const int * pTimeslice,
                             const char * PrintPrefix, std::string *pGroupName,
                             const std::vector<NegateStar> *pAlgSnkNeg, const char * pDSName )
{
  const bool bSameAlgebras{&AlgSnk == &AlgSrc};
  if( pAlgSnkNeg )
  {
    assert( pAlgSnkNeg->size() == AlgSnk.size() );
    assert( pAlgSnkNeg->size() && "Algebras unknown. Negating makes no sense" );
    assert( !bSameAlgebras && "Work out whether two negatives make a positive, then implement" );
  }
  // Parse the name. Not expecting a type, so if present, put it back on the end of Base
  Name_.Parse( FileName );
  if( !Name_.Type.empty() )
  {
    // Not sure whether I should bother doing this?
    Name_.Base.append( 1, '.' );
    Name_.Base.append( Name_.Type );
    Name_.Type.clear();
  }
  if( !Name_.bSeedNum )
    throw std::runtime_error( "Configuration number missing from " + FileName );
  // Strip out timeslice info if present
  bHasTimeslice = pTimeslice;
  if( bHasTimeslice )
    Timeslice_ = * pTimeslice;
  else
  {
    std::string sCopy{ Name_.Base };
    ExtractTimeslice( sCopy, bHasTimeslice, Timeslice_ );
  }
  // Now load the file
  ::H5::H5File f;
  ::H5::Group  g;
  std::string localGroupName;
  if( !pGroupName )
    pGroupName = &localGroupName;
  H5::OpenFileGroup( f, g, FileName, PrintPrefix, pGroupName );
  bool bOK = false;
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  try // to load a single correlator
  {
    if( ( AlgSnk.empty() || ( AlgSnk.size() == 1 && AlgSnk[0] == Gamma::Algebra::Unknown ) )
      && ( AlgSrc.empty() || ( AlgSrc.size() == 1 && AlgSrc[0] == Gamma::Algebra::Unknown ) ) )
    {
      ::H5::DataSet ds = g.openDataSet( ( pDSName && *pDSName ) ? pDSName : "correlator" );
      ::H5::DataSpace dsp = ds.getSpace();
      if( dsp.getSimpleExtentNdims() == 1 )
      {
        hsize_t Dim[1];
        dsp.getSimpleExtentDims( Dim );
        if( Dim[0] <= std::numeric_limits<int>::max() )
        {
          resize( 1, 1, static_cast<int>( Dim[0] ) );
          AlgSnk_.resize( 1 );
          AlgSnk_[0] = Gamma::Algebra::Unknown;
          AlgSrc_.resize( 1 );
          AlgSrc_[0] = Gamma::Algebra::Unknown;
          if( AlgSnk.empty() )
            AlgSnk.push_back( Gamma::Algebra::Unknown );
          if( !bSameAlgebras && AlgSrc.empty() )
            AlgSrc.push_back( Gamma::Algebra::Unknown );
          T * const pData{ (*this)[0] };
          ds.read( pData, H5::Equiv<T>::Type );
          // Negate this correlator if requested
          if( pAlgSnkNeg )
            for( int i = 0; i < Nt_; i++ )
              ApplyNegateStar( pData[i], (*pAlgSnkNeg)[0] );
          bOK = true;
        }
      }
    }
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
    ::H5::Exception::clearErrorStack();
  }
  if( !bOK )
  {
    try // to load from array of correlators indexed by gamma matrix
    {
      unsigned short NumVec;
      ::H5::Attribute a;
      a = g.openAttribute("_Grid_vector_size");
      a.read( ::H5::PredType::NATIVE_USHORT, &NumVec );
      a.close();
      // Must be a perfect square and have at least as many as entries as requested
      const unsigned short NumFileOps{static_cast<unsigned short>( std::sqrt( NumVec ) + 0.5 )};
      if( NumFileOps * NumFileOps == NumVec && NumFileOps >= AlgSnk.size() && NumFileOps >= AlgSrc.size() )
      {
        bOK = true;
        std::vector<int> count;
        for( unsigned short i = 0; bOK && i < NumVec; i++ ) // Loop through all vectors in file
        {
          bOK = false;
          ::H5::Group gi = g.openGroup( *pGroupName + "_" + std::to_string( i ) );
          ::H5::DataSet ds = gi.openDataSet( "corr" );
          ::H5::DataSpace dsp = ds.getSpace();
          if( dsp.getSimpleExtentNdims() == 1 )
          {
            hsize_t Dim[1];
            dsp.getSimpleExtentDims( Dim );
            if( Dim[0] <= std::numeric_limits<int>::max() )
            {
              const int ThisNt{ static_cast<int>( Dim[0] ) };
              if( i == 0 )
              {
                // First correlator - resize and save operators if known
                resize(AlgSnk.size() ? static_cast<int>( AlgSnk.size() ) : NumFileOps,
                       AlgSrc.size() ? static_cast<int>( AlgSrc.size() ) : NumFileOps,
                       ThisNt);
                count.resize( NumSnk_ * NumSrc_, 0 ); // I want to check each operator combination appears once only
                if( AlgSnk.size() )
                {
                  AlgSnk_.resize( AlgSnk.size() );
                  std::copy( AlgSnk.cbegin(), AlgSnk.cend(), AlgSnk_.begin() );
                }
                else
                {
                  AlgSnk_.clear();
                  AlgSnk_.reserve( NumFileOps );
                }
                if( AlgSrc.size() )
                {
                  AlgSrc_.resize( AlgSrc.size() );
                  std::copy( AlgSrc.cbegin(), AlgSrc.cend(), AlgSrc_.begin() );
                }
                else
                {
                  AlgSrc_.clear();
                  AlgSrc_.reserve( NumFileOps );
                }
              }
              else if( ThisNt != Nt_ )
              {
                break;
              }
              // Read the gamma algebra strings and make sure they are valid
              const Gamma::Algebra gSnk{ H5::ReadGammaAttribute( gi, "gamma_snk" ) };
              int idxSnk;
              for( idxSnk = 0; idxSnk < AlgSnk_.size() && AlgSnk_[idxSnk] != gSnk; idxSnk++ )
                ;
              if( idxSnk == AlgSnk_.size() && AlgSnk_.size() < NumSnk_ )
                AlgSnk_.push_back( gSnk );
              bOK = true; // We can safely ignore gamma structures we're not interested in
              if( idxSnk < AlgSnk_.size() )
              {
                const Gamma::Algebra gSrc{ H5::ReadGammaAttribute( gi, "gamma_src" ) };
                int idxSrc;
                for( idxSrc = 0; idxSrc < AlgSrc_.size() && AlgSrc_[idxSrc] != gSrc; idxSrc++ )
                  ;
                if( idxSrc == AlgSrc_.size() && AlgSrc_.size() < NumSrc_ )
                  AlgSrc_.push_back( gSrc );
                if( idxSrc < AlgSrc_.size() )
                {
                  const int idx{ idxSnk * NumSrc_ + idxSrc };
                  T * const pData{ (*this)[idx] };
                  ds.read( pData, H5::Equiv<T>::Type );
                  count[idx]++;
                  if( pAlgSnkNeg )
                    for( int i = 0; i < Nt_; i++ )
                      ApplyNegateStar( pData[i], (*pAlgSnkNeg)[idxSnk] );
                }
              }
            }
          }
        }
        // Make sure that everything we wanted was loaded once and only once
        for( int i = 0; bOK && i < NumSrc_ * NumSnk_; i++ )
          if( count[i] != 1 )
            bOK = false;
      }
    }
    catch(const ::H5::Exception &)
    {
      bOK = false;
      ::H5::Exception::clearErrorStack();
    }
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  if( !bOK )
    throw std::runtime_error( "Unable to read sample from " + FileName );
  if( !IsFinite() )
    throw std::runtime_error( "Values read are not all finite" );
  // If I'm discovering which operators are in the file, copy them back to caller
  // Bear in mind that caller may have passed in the same array for each gamma algebra
  bool bCopyBackSrc{ AlgSrc.empty() };
  bool bCopyBackSnk{ !bSameAlgebras && AlgSnk.empty() };
  if( bCopyBackSrc )
  {
    AlgSrc.resize( AlgSrc_.size() );
    std::copy( AlgSrc_.cbegin(), AlgSrc_.cend(), AlgSrc.begin() );
  }
  if( bCopyBackSnk )
  {
    AlgSnk.resize( AlgSnk_.size() );
    std::copy( AlgSnk_.cbegin(), AlgSnk_.cend(), AlgSnk.begin() );
  }
}

template <typename T>
void CorrelatorFile<T>::WriteSummary( const std::string &Prefix, const std::vector<Gamma::Algebra> &AlgSnk, const std::vector<Gamma::Algebra> &AlgSrc )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  const int nt{ Nt() };
  const std::vector<Gamma::Algebra> &MySnk{ AlgSnk.size() ? AlgSnk : AlgSnk_ };
  const std::vector<Gamma::Algebra> &MySrc{ AlgSrc.size() ? AlgSrc : AlgSrc_ };
  const int NumSnk{ static_cast<int>( MySnk.size() ) };
  const int NumSrc{ static_cast<int>( MySrc.size() ) };
  std::string sOutFileName{ Prefix };
  sOutFileName.append( Name_.Base );
  std::size_t Len{ sOutFileName.length() };
  std::string sSuffix( 1, '.' );
  sSuffix.append( Name_.SeedString );
  sSuffix.append( 1, '.' );
  sSuffix.append( TEXT_EXT );
  for( int Snk = 0; Snk < NumSnk; Snk++ )
  {
    static const char pszSep[] = "_";
    sOutFileName.resize( Len );
    sOutFileName.append( Common::Gamma::NameShort( AlgSnk[Snk], pszSep ) );
    std::size_t Len2{ sOutFileName.length() };
    for( int Src = 0; Src < NumSrc; Src++ )
    {
      sOutFileName.resize( Len2 );
      sOutFileName.append( Common::Gamma::NameShort( AlgSrc[Src], pszSep ) );
      sOutFileName.append( sSuffix );
      SummaryHelper( sOutFileName, (*this)( AlgSnk[Snk], AlgSrc[Src] ), nt );
    }
  }
}

template class CorrelatorFile<double>;
template class CorrelatorFile<float>;
template class CorrelatorFile<std::complex<double>>;
template class CorrelatorFile<std::complex<float>>;

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

/*template class Sample<double>;
template class Sample<float>;
template class Sample<std::complex<double>>;
template class Sample<std::complex<float>>;*/

const std::string ModelBase::EnergyPrefix{ "E" };
const std::string ModelBase::EDiffPrefix{ "EDiff" };
const std::string ModelBase::ConstantPrefix{ "K" };

template <typename T>
JackBootColumn<T> Model<T>::Column( const Param::Key &k, std::size_t Index )
{
  Params::const_iterator it = params.FindPromiscuous( k );
  if( it == params.cend() )
  {
    std::ostringstream os;
    os << k.FullName( Index, std::numeric_limits<std::size_t>::max() ) << " not found";
    throw std::runtime_error( os.str().c_str() );
  }
  const Param &p{ it->second };
  if( Index >= p.size )
  {
    std::ostringstream os;
    os << k.FullName( Index, p.size ) << " not found - only " << p.size << " elements";
    throw std::runtime_error( os.str().c_str() );
  }
  return Base::Column( p.GetOffset( Index, Param::Type::All ) );
}

template <typename T>
JackBootColumn<T> Model<T>::ColumnFrozen( const Param::Key &k, std::size_t Index )
{
  Params::const_iterator it = params.FindPromiscuous( k );
  if( it == params.cend() )
  {
    std::ostringstream os;
    os << k.FullName( Index, std::numeric_limits<std::size_t>::max() ) << " not found";
    throw std::runtime_error( os.str().c_str() );
  }
  const Param &p{ it->second };
  if( Index >= p.size )
  {
    std::ostringstream os;
    os << k.FullName( Index, p.size ) << " not found - only " << p.size << " elements";
    throw std::runtime_error( os.str().c_str() );
  }
  return Base::ColumnFrozen( p.GetOffset( Index, Param::Type::All ) );
}

template <typename T>
void Sample<T>::resize( int NumReplicas, int Nt )
{
  const bool NtChanged{ Nt_ != Nt };
  if( NtChanged || NumSamples() != NumReplicas )
  {
    if( NtChanged )
    {
      Nt_ = Nt;
      AllocSummaryBuffer();
      ColumnNames.clear();
      Raw.clear();
      Binned[0].clear();
      Data[0].clear();
    }
    Data[0].resize( NumReplicas, Nt );
    for( int i = 1; i < NumJackBoot; ++i )
    {
      Binned[i].clear();
      Data[i].clear();
    }
  }
}

template <typename T>
void Sample<T>::clear( bool bClearName )
{
  Nt_ = 0;
  Raw.clear();
  for( int i = 0; i < NumJackBoot; ++i )
  {
    Binned[i].clear();
    Data[i].clear();
  }
  SummaryNames.clear();
  m_SummaryData.clear();
  ColumnNames.clear();
  SeedMachine_.clear();
  binSize = 1;
  SampleSize = 0;
  ConfigCount.clear();
  FileList.clear();
  if( bClearName )
    Name_.clear();
}

template <typename T>
void Sample<T>::WriteColumnNames( std::ostream &s ) const
{
  const std::string Sep{ " " };
  if( Nt_ == 0 )
    throw std::runtime_error( "Can't write header - Nt_ = 0" );
  std::string Buffer{ "t" };
  const std::size_t BufferLen{ Buffer.size() };
  const std::string * ps{ &Buffer };
  for( std::size_t t = 0; t < Nt_; t++ )
  {
    if( ColumnNames.empty() )
    {
      Buffer.resize( BufferLen );
      Buffer.append( std::to_string( t ) );
    }
    else
      ps = &ColumnNames[t];
    if( t )
      s << Sep;
    Common::ValWithEr<T>::Header( *ps, s, Sep );
  }
}

// Initialise *pNumSamples either to 0, or to the size of the first Sample before first call
template <typename T> template <typename U>
void Sample<T>::IsCompatible( const Sample<U> &o, int * pNumSamples, unsigned int CompareFlags,
                              bool bSuffix ) const
{
  static const std::string sPrefix{ "Incompatible " + sBootstrap + " samples - " };
  const std::string sSuffix{ bSuffix ? ":\n  " + Name_.Filename + "\nv " + o.Name_.Filename : "" };
  if( !( CompareFlags & COMPAT_DISABLE_BASE ) && !Common::EqualIgnoreCase( o.Name_.Base, Name_.Base ) )
    throw std::runtime_error( sPrefix + "base " + o.Name_.Base + sNE + Name_.Base + sSuffix );
  if( !( CompareFlags & COMPAT_DISABLE_TYPE ) && !Common::EqualIgnoreCase( o.Name_.Type, Name_.Type ) )
    throw std::runtime_error( sPrefix + "type " + o.Name_.Type + sNE + Name_.Type + sSuffix );
  if( !Common::EqualIgnoreCase( o.Name_.Ext, Name_.Ext ) )
    throw std::runtime_error( sPrefix + "extension " + o.Name_.Ext + sNE + Name_.Ext + sSuffix );
  if( pNumSamples )
  {
    if( *pNumSamples == 0 || *pNumSamples > NumSamples() )
      *pNumSamples = NumSamples();
    if( *pNumSamples > o.NumSamples() )
      *pNumSamples = o.NumSamples();
  }
  else if( o.NumSamples() != NumSamples() )
    throw std::runtime_error( sPrefix + "NumSamples " + std::to_string(o.NumSamples()) +
                             sNE + std::to_string(NumSamples()) + sSuffix );
  if( !( CompareFlags & COMPAT_DISABLE_NT ) && o.Nt() != Nt_ )
    throw std::runtime_error( sPrefix + "Nt " + std::to_string(o.Nt()) +
                             sNE + std::to_string(Nt_) + sSuffix );
  if( o.SampleSize && SampleSize && o.SampleSize != SampleSize )
    throw std::runtime_error( sPrefix + "SampleSize " + std::to_string(o.SampleSize) +
                             sNE + std::to_string(SampleSize) + sSuffix );
  // NB: If the random numbers were different, the random number cache would have caught this
  if( o.Seed() != Seed() )
    throw std::runtime_error( "Seed " + o.SeedString() + sNE + SeedString() );
  const std::size_t CSize {   ConfigCount.size() };
  const std::size_t CSizeO{ o.ConfigCount.size() };
  if( CSize && CSizeO )
  {
    if( CSizeO != CSize )
      throw std::runtime_error( sPrefix + "Number of configs " +
                          std::to_string(CSizeO) + sNE + std::to_string(CSize) + sSuffix );
    for( std::size_t i = 0; i < CSize; i++ )
    {
      const Common::ConfigCount &l{   ConfigCount[i] };
      const Common::ConfigCount &r{ o.ConfigCount[i] };
      if( r.Config != l.Config )
        throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                 sNE + std::to_string(l.Config) + sSuffix );
      /*if( !( CompareFlags & COMPAT_DISABLE_CONFIG_COUNT ) && r.Count != l.Count )
        throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                 ", NumTimeslices " + std::to_string(r.Count) + sNE +
                                 std::to_string(l.Count) + sSuffix );*/
    }
  }
  if( !Ensemble.empty() && !o.Ensemble.empty() && !EqualIgnoreCase( Ensemble, o.Ensemble ) )
    throw std::runtime_error( "Ensemble " + Ensemble + " != " + o.Ensemble );
}

// Take ownership of the FileList
template <typename T>
template <typename U> void Sample<T>::CopyAttributes( const Sample<U> &in )
{
  binSize = in.binSize;
  SampleSize = in.SampleSize;
  ConfigCount = in.ConfigCount;
  Ensemble = in.Ensemble;
  if( FileList.empty() )
    FileList = in.FileList;
  // Copy the random numbers for this sample
  SeedMachine_ = in.SeedMachine_;
  SetSeed( in.Seed() );
}

template void Sample<float>::CopyAttributes( const Sample<float> &in );
template void Sample<float>::CopyAttributes( const Sample<double> &in );
template void Sample<float>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<float>::CopyAttributes( const Sample<std::complex<double>> &in );
template void Sample<double>::CopyAttributes( const Sample<float> &in );
template void Sample<double>::CopyAttributes( const Sample<double> &in );
template void Sample<double>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<double>::CopyAttributes( const Sample<std::complex<double>> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<float> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<double> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<std::complex<double>> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<float> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<double> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<std::complex<double>> &in );

// Auto bin: all measurements on each config binned into a single measurement
template <typename T> void Sample<T>::BinAuto( int JackBootDest )
{
  if( !NumSamplesRaw() )
    throw std::runtime_error( "No raw data to bin" );
  // Work out how many bins there are and whether they are all the same size
  int NumSamplesRaw{ 0 };
  int NewBinSize{ ConfigCount[0].Count };
  for( const Common::ConfigCount &cc : ConfigCount )
  {
    if( cc.Count < 1 )
      throw std::runtime_error( "ConfigCount invalid" );
    NumSamplesRaw += cc.Count;
    if( NewBinSize && NewBinSize != cc.Count )
      NewBinSize = 0; // Indicates the bin size varies per ConfigCount
  }
  if( NumSamplesRaw != this->NumSamplesRaw() )
    throw std::runtime_error( "ConfigCount doesn't match raw data" );
  // Rebin
  Matrix<T> &mSrc = getRaw();
  Matrix<T> &mDst = resizeBinned( static_cast<int>( ConfigCount.size() ), JackBootDest );
  if( JackBootDest == 0 )
    binSize = NewBinSize; // File attribute is based on the primary sample
  Vector<T> vSrc, vDst;
  std::size_t rowSrc = 0;
  for( std::size_t Bin = 0; Bin < ConfigCount.size(); ++Bin )
  {
    for( std::size_t i = 0; i < ConfigCount[Bin].Count; ++i )
    {
      for( std::size_t t = 0; t < Nt_; ++t )
      {
        if( i )
          mDst( Bin, t ) += mSrc( rowSrc, t );
        else
          mDst( Bin, t ) = mSrc( rowSrc, t );
      }
      ++rowSrc;
    }
    for( std::size_t t = 0; t < Nt_; ++t )
      mDst( Bin, t ) /= ConfigCount[Bin].Count;
  }
}

template <typename T> void Sample<T>::BinFixed( int binSize_, int JackBootDest )
{
  if( binSize_ < 1 )
    throw std::runtime_error( "BinSize " + std::to_string( binSize_ ) + " invalid" );
  if( !NumSamplesRaw() )
    throw std::runtime_error( "No raw data to bin" );
  // Work out how many samples there are and check they match ConfigCount
  int NumSamplesRaw{ 0 };
  for( const Common::ConfigCount &cc : ConfigCount )
  {
    if( cc.Count != ConfigCount[0].Count )
      throw std::runtime_error( "Can't manually bin - raw sample counts differ per config" );
    NumSamplesRaw += cc.Count;
  }
  if( NumSamplesRaw != this->NumSamplesRaw() )
    throw std::runtime_error( "ConfigCount doesn't match raw data" );
  const bool PartialLastBin = NumSamplesRaw % binSize_;
  const int NewNumSamples{ static_cast<int>( NumSamplesRaw / binSize_ + ( PartialLastBin ? 1 : 0 ) ) };
  // Rebin
  Matrix<T> &mSrc = getRaw();
  Matrix<T> &mDst = resizeBinned( NewNumSamples, JackBootDest );
  if( JackBootDest == 0 )
    binSize = binSize_; // File attribute is based on the primary sample
  Vector<T> vSrc, vDst;
  std::size_t rowSrc = 0;
  int ThisBinSize{ binSize_ };
  for( std::size_t Bin = 0; Bin < NewNumSamples; ++Bin )
  {
    if( PartialLastBin && Bin == NewNumSamples - 1 )
      ThisBinSize = NumSamplesRaw % binSize_;
    for( std::size_t i = 0; i < ThisBinSize; ++i )
    {
      for( std::size_t t = 0; t < Nt_; ++t )
      {
        if( i )
          mDst( Bin, t ) += mSrc( rowSrc, t );
        else
          mDst( Bin, t ) = mSrc( rowSrc, t );
      }
      ++rowSrc;
    }
    for( std::size_t t = 0; t < Nt_; ++t )
      mDst( Bin, t ) /= ThisBinSize;
  }
}

/// Perform bootstrap or jackknife, depending on Seed
template <typename T>
void Sample<T>::Resample( int JackBootDest, int NumReplicas, int BinnedSource, SeedType Seed )
{
  ValidateJackBoot( JackBootDest );
  ValidateJackBoot( BinnedSource );
  Data[JackBootDest].Seed = Seed;
  Data[JackBootDest].Resample( Binned[BinnedSource], NumReplicas );
  if( JackBootDest == 0 )
  {
    SampleSize = NumSamplesBinned( BinnedSource );
    SeedMachine_ = GetHostName();
  }
}

template<typename T> inline typename std::enable_if<!SampleTraits<T>::is_complex>::type
CopyOldFormat( T &Dest, const double &Real, const double &Imag )
{
  if( Imag )
    throw std::runtime_error( "Complex sample has imaginary component" );
  using Scalar = typename SampleTraits<T>::scalar_type;
  Dest = static_cast<Scalar>( Real );
}

template<typename T> inline typename std::enable_if<SampleTraits<T>::is_complex>::type
CopyOldFormat( T &Dest, const double &Real, const double &Imag )
{
  using Scalar = typename SampleTraits<T>::scalar_type;
  Dest.real( static_cast<Scalar>( Real ) );
  Dest.imag( static_cast<Scalar>( Imag ) );
}

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void Sample<T>::Read( const char *PrintPrefix, std::string *pGroupName )
{
  if( !Name_.bSeedNum )
    throw std::runtime_error( "Seed missing from " + Name_.Filename );
  if( Name_.Type.empty() )
    throw std::runtime_error( "Type missing from " + Name_.Filename );
  const SeedType NewSeed{ RandomCache::DefaultSeed() };
  ::H5::H5File f;
  ::H5::Group  g;
  H5::OpenFileGroup( f, g, Name_.Filename, PrintPrefix, pGroupName );
  bool bOK{ false };
  bool bNoReturn{ false };
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  clear( false );
  try // to load from LatAnalyze format
  {
    SetSeed( Name_.Seed ); // Seed comes from filename if not loaded from hdf5
    unsigned short att_Type;
    ::H5::Attribute a;
    a = g.openAttribute("type");
    a.read( ::H5::PredType::NATIVE_USHORT, &att_Type );
    a.close();
    if( att_Type == 2 && Name_.Seed == NewSeed )
    {
      unsigned long  att_nSample;
      a = g.openAttribute("nSample");
      a.read( ::H5::PredType::NATIVE_ULONG, &att_nSample );
      a.close();
      std::vector<double> buffer;
      for( int i = idxCentral; i < NumSamples(); i++ )
      {
        ::H5::DataSet ds = g.openDataSet( "data_" + ( i == idxCentral ? "C" : "S_" + std::to_string( i ) ) );
        ::H5::DataSpace dsp = ds.getSpace();
        if( dsp.getSimpleExtentNdims() == 2 )
        {
          hsize_t Dim[2];
          dsp.getSimpleExtentDims( Dim );
          if( Dim[1] != 2 || Dim[0] * att_nSample > std::numeric_limits<int>::max() )
            break;
          if( i == idxCentral )
          {
            resize( static_cast<int>( att_nSample ), static_cast<int>( Dim[0] ) );
            buffer.resize( 2 * Nt_ );
          }
          else if( Dim[0] != static_cast<unsigned int>( Nt_ ) )
            break;
          ds.read( buffer.data(), ::H5::PredType::NATIVE_DOUBLE );
          if( !Common::IsFinite( buffer ) )
            throw std::runtime_error( "NANs at row " + std::to_string(i) + "in " + Name_.Filename );
          for( int t = 0; t < Nt_; t++ )
            CopyOldFormat( Data[0]( i, t ), buffer[t], buffer[t + Nt_] );
          // Check whether this is the end
          if( i == 0 )
            bNoReturn = true;
          if( i == NumSamples() - 1 )
            bOK = true;
        }
      }
    }
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
    ::H5::Exception::clearErrorStack();
  }
  catch(...)
  {
    ::H5::Exception::setAutoPrint(h5at, f5at_p);
    clear( false );
    throw;
  }
  if( !bOK && !bNoReturn )
  {
    clear( false );
    try // to load from my format
    {
      unsigned short att_nAux;
      ::H5::Attribute a;
      // Auxilliary rows removed 30 Apr 2023. Was never used
      a = g.openAttribute("nAux");
      a.read( ::H5::PredType::NATIVE_USHORT, &att_nAux );
      a.close();
      if( att_nAux )
        throw std::runtime_error( std::to_string( att_nAux ) + " unexpected auxilliary rows in "
                                 + Name_.Filename );
      unsigned int att_nSample;
      a = g.openAttribute("nSample");
      a.read( ::H5::PredType::NATIVE_UINT, &att_nSample );
      a.close();
      if( att_nSample > std::numeric_limits<int>::max() )
        throw std::runtime_error( "nSample " + std::to_string( att_nSample ) + " invalid" );
      try
      {
        SeedType tmp;
        a = g.openAttribute( RandomCache::sSeed );
        a.read( H5::Equiv<SeedType>::Type, &tmp );
        a.close();
        SetSeed( tmp );
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        SetSeed( 0 ); // In this format, not present means 0
      }
      try
      {
        a = g.openAttribute( sEnsemble );
        a.read( a.getStrType(), Ensemble );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        Ensemble.clear();
      }
      try
      {
        a = g.openAttribute( RandomCache::sSeedMachine );
        a.read( a.getStrType(), SeedMachine_ );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        SeedMachine_.clear();
      }
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sSummaryNames);
        SummaryNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sColumnNames);
        ColumnNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        int tmp;
        a = g.openAttribute( sSampleSize );
        a.read( ::H5::PredType::NATIVE_INT, &tmp );
        a.close();
        SampleSize = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        int tmp;
        a = g.openAttribute( sBinSize );
        a.read( ::H5::PredType::NATIVE_INT, &tmp );
        a.close();
        binSize = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      // Try to load the FileList from dataSet first, falling back to attribute
      bool bGotFileList{ false };
      try
      {
        ::H5::DataSet ds = g.openDataSet( sFileList );
        FileList = H5::ReadStrings( ds );
        ds.close();
        bGotFileList = true;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      if( !bGotFileList )
      {
        try
        {
          a = g.openAttribute(sFileList);
          FileList = H5::ReadStrings( a );
          a.close();
          bGotFileList = true;
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
        }
      }
      try
      {
        a = g.openAttribute(sConfigCount);
        ::H5::DataSpace dsp = a.getSpace();
        const int rank{ dsp.getSimpleExtentNdims() };
        if( rank != 1 )
          throw std::runtime_error( sConfigCount + " number of dimensions=" + std::to_string( rank ) + ", expecting 1" );
        hsize_t NumConfig;
        dsp.getSimpleExtentDims( &NumConfig );
        if( NumConfig > std::numeric_limits<int>::max() )
          throw std::runtime_error( sConfigCount + " too many items in ConfigCount: " + std::to_string( NumConfig ) );
        ConfigCount.resize( NumConfig );
        a.read( H5::Equiv<Common::ConfigCount>::Type, &ConfigCount[0] );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        ReadAttributes( g );
      }
      catch(const std::exception &e)
      {
        throw;
      }
      catch(...)
      {
        throw std::runtime_error( "Extended attributes for " + DefaultGroupName() + " missing" );
      }
      // Load Binned data
      try
      {
        H5::ReadMatrix( g, "samples_Binned", Binned[0] );
        if( Binned[0].size1 != SampleSize )
          throw std::runtime_error( "SampleSize = " + std::to_string( SampleSize ) + " but "
                      + std::to_string( Binned[0].size1 ) + " binned replicas read" );
        if( Nt_ == 0 )
          Nt_ = static_cast<int>( Binned[0].size2 );
        else if( Binned[0].size2 != Nt_ )
          throw std::runtime_error( "nT = " + std::to_string( Nt_ ) + " but "
                      + std::to_string( Binned[0].size2 ) + " binned columns read" );
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        Binned[0].clear();
      }
      // Load raw data
      try
      {
        H5::ReadMatrix( g, "samples_Raw", Raw );
        if( Nt_ == 0 )
          Nt_ = static_cast<int>( Raw.size2 );
        else if( Raw.size2 != Nt_ )
          throw std::runtime_error( "nT = " + std::to_string( Nt_ ) + " but "
                      + std::to_string( Raw.size2 ) + " raw columns read" );
        if( FileList.size() && Raw.size1 != FileList.size() )
          std::cout << "* Warning: FileList has " << FileList.size() << " entries but "
                    << Raw.size1 << " raw replicas read\n";
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        Raw.clear();
      }
      // If random numbers present, put them in cache (which checks compatibility)
      if( Seed() == NewSeed )
      {
        try
        {
          Matrix<fint> Random;
          H5::ReadMatrix( g, RandomCache::sRandom, Random );
          if( Random.size2 != SampleSize )
            throw std::runtime_error( "SampleSize = " + std::to_string( SampleSize ) + " but "
                                     + std::to_string( Random.size2 ) + " random columns read" );
          if( Random.size1 != att_nSample )
            throw std::runtime_error( "nSample = " + std::to_string( att_nSample ) + " but "
                                     + std::to_string( Random.size1 ) + " random replicas read" );
          // Put random numbers in cache - will throw if they don't match previous entries
          RandomCache::Global.Put( Seed(), Random );
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
        }
      }
      const int NewNumReplicas{ NewSeed == SeedWildcard ? SampleSize // Jackknife
                                      : static_cast<int>( RandomCache::DefaultNumReplicas() ) };
      // If I don't need to resample, load optional items
      const bool bNeedResample{ Seed() != NewSeed || att_nSample < NewNumReplicas };
      if( bNeedResample && !CanResample() )
      {
        static const char szTo[]{ " -> " };
        std::ostringstream os;
        os << "Can't resample";
        if( att_nSample < NewNumReplicas )
          os << Space << att_nSample << szTo << NewNumReplicas << " replicas";
        if( Seed() != Name_.Seed)
          os << " seed " << SeedString() << szTo
             << RandomCache::SeedString( Name_.Seed );
        throw std::runtime_error( os.str().c_str() );
      }
      if( !bNeedResample )
      {
        try
        {
          // Load resampled data
          Data[0].Read( g, "data" );
          if( Data[0].NumReplicas() != att_nSample )
            throw std::runtime_error( "nSample = " + std::to_string( att_nSample ) + " but "
                                     + std::to_string( Data[0].NumReplicas() ) + " replicas read" );
          if( Nt_ == 0 )
            Nt_ = static_cast<int>( Data[0].extent() );
          else if( Data[0].extent() != Nt_ )
            throw std::runtime_error( "nT = " + std::to_string( Nt_ ) + " but "
                        + std::to_string( Data[0].extent() ) + " columns read" );
          if( !SummaryNames.empty() )
          {
            AllocSummaryBuffer();
            try // to load the summaries
            {
              ::H5::DataSet ds = g.openDataSet( sSummaryDSName );
              ::H5::DataSpace dsp = ds.getSpace();
              if( dsp.getSimpleExtentNdims() != 2 )
                throw std::runtime_error( "Dimension error reading " + sSummaryDSName );
              hsize_t Dim[2];
              dsp.getSimpleExtentDims( Dim );
              if( Dim[0] != SummaryNames.size() || Dim[1] != static_cast<unsigned int>( Nt_ ) )
                throw std::runtime_error( "Bad size reading " + sSummaryDSName );
              const ::H5::DataType ErrType{ ds.getDataType() };
              if( ErrType == H5::Equiv<ValWithErOldV1<scalar_type>>::Type )
              {
                const std::size_t Count{ SummaryNames.size() * Nt_ };
                std::vector<ValWithErOldV1<scalar_type>> tmp( Count );
                ds.read( &tmp[0], ErrType );
                for( std::size_t i = 0; i < Count; ++i )
                {
                  m_SummaryData[i].Central = tmp[i].Central;
                  m_SummaryData[i].Low     = tmp[i].Low;
                  m_SummaryData[i].High    = tmp[i].High;
                  m_SummaryData[i].Check   = tmp[i].Check;
                  m_SummaryData[i].Min     = 0;
                  m_SummaryData[i].Max     = 0;
                }
              }
              else
                ds.read( &m_SummaryData[0], H5::Equiv<ValWithEr<scalar_type>>::Type );
            }
            catch(const ::H5::Exception &)
            {
              ::H5::Exception::clearErrorStack();
              MakeCorrSummary(); // Rebuild summaries if I can't load them
            }
          }
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
        }
      }
      // Either resampled or binned data must be available
      if( Data[0].NumReplicas() == 0 && Binned[0].Empty() )
        throw std::runtime_error( "Neither binned nor resampled data available" );
      if( bNeedResample && CanResample() )
      {
        Data[0].Seed = NewSeed;
        Data[0].Resample( Binned[0], NewNumReplicas );
        MakeCorrSummary();
      }
      // If we get here then we've loaded OK
      ValidateAttributes();
      bOK = true;
    }
    catch(const ::H5::Exception &)
    {
      bOK = false;
      ::H5::Exception::clearErrorStack();
    }
    catch(...)
    {
      ::H5::Exception::setAutoPrint(h5at, f5at_p);
      clear( false );
      throw;
    }
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  if( !bOK )
  {
    clear( false );
    throw std::runtime_error( "Unable to read sample from " + Name_.Filename );
  }
}

template <typename T>
void Sample<T>::Write( const std::string &FileName, const char * pszGroupName )
{
  SetName( FileName );
  const std::string GroupName{ pszGroupName == nullptr || *pszGroupName == 0 ? DefaultGroupName() : pszGroupName };
  bool bOK = false;
  try // to write in my format
  {
    ::H5::H5File f( FileName, H5F_ACC_TRUNC );
    hsize_t Dims[2];
    Dims[0] = 1;
    ::H5::DataSpace ds1( 1, Dims );
    ::H5::Group g = f.createGroup( GroupName );
    ::H5::Attribute a = g.createAttribute( "nAux", ::H5::PredType::STD_U16LE, ds1 );
    int tmp{ NumExtraSamples - 1 };
    a.write( ::H5::PredType::NATIVE_INT, &tmp );
    a.close();
    a = g.createAttribute( "nSample", ::H5::PredType::STD_U32LE, ds1 );
    tmp = NumSamples();
    a.write( ::H5::PredType::NATIVE_INT, &tmp );
    a.close();
    int NumAttributes = 2;
    if( Seed() )
    {
      const SeedType tmp{ Seed() };
      a = g.createAttribute( RandomCache::sSeed, ::H5::PredType::STD_U32LE, ds1 );
      a.write( H5::Equiv<SeedType>::Type, &tmp );
      a.close();
      NumAttributes++;
    }
    if( !Ensemble.empty() )
    {
      a = g.createAttribute( sEnsemble, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, Ensemble );
      a.close();
      NumAttributes++;
    }
    if( !SeedMachine_.empty() )
    {
      a = g.createAttribute( RandomCache::sSeedMachine, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, SeedMachine_ );
      a.close();
      NumAttributes++;
    }
    if( SummaryNames.size() )
    {
      H5::WriteAttribute( g, sSummaryNames, SummaryNames );
      NumAttributes++;
    }
    if( ColumnNames.size() )
    {
      H5::WriteAttribute( g, sColumnNames, ColumnNames );
      NumAttributes++;
    }
    if( SampleSize )
    {
      a = g.createAttribute( sSampleSize, ::H5::PredType::STD_U32LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &SampleSize );
      a.close();
      NumAttributes++;
    }
    if( binSize != 1 )
    {
      a = g.createAttribute( sBinSize, ::H5::PredType::STD_I32LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &binSize );
      a.close();
      NumAttributes++;
    }
    ::H5::DataSpace dsp;
    if( ConfigCount.size() )
    {
      Dims[0] = ConfigCount.size();
      dsp = ::H5::DataSpace( 1, Dims );
      a = g.createAttribute( sConfigCount, H5::Equiv<Common::ConfigCount>::Type, dsp );
      a.write( H5::Equiv<Common::ConfigCount>::Type, &ConfigCount[0] );
      a.close();
      dsp.close();
      NumAttributes++;
    }
    NumAttributes += WriteAttributes( g );
    if( FileList.size() )
      H5::WriteStringData( g, sFileList, FileList );
    // If I don't have binned data, I must save the bootstrap (because it can't be recreated)
    if( Binned[0].Empty() || RandomCache::SaveFatFiles() )
      Data[0].Write( g, "data" );
    // Only save random numbers if asked, and then only for bootstrap (jackknife don't need random)
    if( RandomCache::SaveFatFiles() && Seed() != SeedWildcard && SampleSize )
    {
      const fint &cSeed{ Seed() };
      Matrix<fint> mRnd{ RandomCache::Global.Get( cSeed, SampleSize, Data[0].extent() ) };
      Dims[0] = mRnd.size1;
      Dims[1] = mRnd.size2;
      dsp = ::H5::DataSpace( 2, Dims );
      // The values in this table go from 0 ... NumSamples_ - 1, so choose a space-minimising size in the file
      ::H5::DataSet ds = g.createDataSet( RandomCache::sRandom, SampleSize<=static_cast<int>(std::numeric_limits<std::uint8_t>::max())
        ? ::H5::PredType::STD_U8LE : SampleSize<=static_cast<int>(std::numeric_limits<std::uint16_t>::max())
        ? ::H5::PredType::STD_U16LE: ::H5::PredType::STD_U32LE, dsp );
      ds.write( mRnd.Data(), H5::Equiv<std::uint_fast32_t>::Type );
      ds.close();
      dsp.close();
    }
    if( !Raw.Empty() )
      H5::WriteMatrix( g, "samples_Raw", Raw );
    if( !Binned[0].Empty() )
      H5::WriteMatrix( g, "samples_Binned", Binned[0] );
    if( SummaryNames.size() )
    {
      Dims[0] = SummaryNames.size();
      Dims[1] = Nt_;
      dsp = ::H5::DataSpace( 2, Dims );
      ::H5::DataSet ds = g.createDataSet( sSummaryDSName, H5::Equiv<ValWithEr<scalar_type>>::Type, dsp );
      ds.write( &SummaryData(0,0), H5::Equiv<ValWithEr<scalar_type>>::Type );
      ds.close();
      dsp.close();
    }
    g = f.openGroup( "/" );
    a = g.createAttribute( "_Grid_dataset_threshold", ::H5::PredType::STD_U32LE, ds1);
    a.write( ::H5::PredType::NATIVE_INT, &NumAttributes );
    a.close();
    ds1.close();
    g.close();
    bOK = true;
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
  }
  if( !bOK )
    throw std::runtime_error( "Unable to write sample to " + FileName + ", group " + GroupName );
}

// Make a summary of the data
// If a name is specified, simply save the central replica
// Otherwise, calculate exp and cosh mass as well

template <typename T>
void Sample<T>::MakeCorrSummary()
{
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  // Now perform summaries
  AllocSummaryBuffer();
  const int tMid{ bFolded() ? Nt_ : Nt_ / 2 };
  JackBoot<scalar_type> Buffer( NumSamples(), Nt_ );
  std::vector<ValEr> veBuf( Nt_ );
  const int NumRealFields{ static_cast<int>( SummaryNames.size() ) / scalar_count };
  for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
  {
    for( int f = 0; f < std::min( NumRealFields, 4 ); ++f ) // each field
    {
      const int SummaryRow{ f + i * NumRealFields };
      // Fill in the JackBoot sample for this field
      if( f == 1 )
      {
        // Same as previous row - except average is taken from average of all samples
        JackBoot<scalar_type>::MakeMean( Buffer.GetCentral(), Buffer.Replica );
      }
      else
      {
        for( int n = idxCentral; n < NumSamples(); ++n )
        {
          for( int t = 0; t < Nt_; ++t )
          {
            // Index of previous and next timeslice relative to current timeslice
            const int tNext{ ( t == Nt_ - 1 ? 0       : t + 1 ) };
            const int tPrev{ ( t == 0       ? Nt_ - 1 : t - 1 ) };
            const scalar_type z{ Traits::RealImag( (*this)(n,t), i ) };
            const scalar_type zNext{ Traits::RealImag( (*this)(n,tNext), i ) };
            const scalar_type zPrev{ Traits::RealImag( (*this)(n,tPrev), i ) };
            scalar_type d;
            switch( f )
            {
              case 0: // central value
                d = z;
                break;
              case 2: // exponential mass
                if( t == 0 )
                  d = NaN;
                else if( t <= tMid )
                  d = std::log( zPrev / z );
                else
                  d = std::log( zNext / z );
                break;
              case 3: // cosh mass
                d = std::acosh( ( zPrev + zNext ) / ( z * 2 ) );
                break;
            }
            Buffer(n,t) = d;
          }
        }
      }
      // Now get the averages
      Buffer.MakeStatistics( veBuf );
      for( int t = 0; t < Nt_; ++t )
        SummaryData( SummaryRow, t ) = veBuf[t];
    }
  }
}

template <typename T>
void Sample<T>::WriteSummary( const std::string &sOutFileName, bool bVerboseSummary )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  if( SummaryNames.empty() )
    throw std::runtime_error( "Summaries can't be written because they've not been created" );
  std::ofstream s( sOutFileName );
  SummaryHeader<T>( s, sOutFileName );
  SummaryComments( s, bVerboseSummary );
  SummaryColumnNames( s );
  s << NewLine;
  SummaryContents( s );
  s << NewLine;
}

template <typename T>
void Sample<T>::SummaryComments( std::ostream & s, bool bVerboseSummary ) const
{
  s << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha
    << "# Seed: " << SeedString() << NewLine;
  if( !SeedMachine_.empty() ) s << "# Seed machine: " << SeedMachine_ << NewLine;
  s << "# Bootstrap: " << NumSamples() << NewLine << "# SampleSize: " << SampleSize << NewLine;
  if( binSize != 1 ) s << "# BinSize: " << binSize << NewLine;
  if( !ConfigCount.empty() )
  {
    s << "# Configs: " << ConfigCount.size();
    for( const Common::ConfigCount &cc : ConfigCount )
      s << CommaSpace << cc;
    s << NewLine;
  }
  s << "# NumColumns: " << Nt();
  if( ColumnNames.empty() )
    s << " (rows t0 ... t" << ( Nt() - 1 ) << ")";
  else
  {
    s << NewLine << "# " << sColumnNames << ": ";
    for( std::size_t i = 0; i < ColumnNames.size(); ++i )
    {
      if( i )
        s << CommaSpace;
      s << ColumnNames[i];
    }
  }
  s << NewLine;
  if( bVerboseSummary )
  {
    s << "# FileCount: " << FileList.size() << NewLine;
    std::size_t i = 0;
    for( const std::string &f : FileList )
      s << "# File " << i++ << ": " << f << NewLine;
  }
}

template <typename T>
void Sample<T>::SummaryColumnNames( std::ostream &os ) const
{
  if( ColumnNames.empty() )
  {
    // First column is timeslice, then each of the summary columns per timeslice
    os << "t";
    for( const std::string &Name : SummaryNames )
    {
      os << Space;
      ValWithEr<T>::Header( Name, os, Space );
    }
  }
  else
  {
    // Show the bootstrap central replica for every column on one row
    bool bFirst{ true };
    for( const std::string &Name : ColumnNames )
    {
      if( bFirst )
        bFirst = false;
      else
        os << Space;
      ValWithEr<T>::Header( Name, os, Space );
    }
  }
}

template <typename T>
void Sample<T>::SummaryContents( std::ostream &os ) const
{
  os << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
  if( ColumnNames.empty() )
  {
    // One row per timeslice, showing each of the summary values
    for( int t = 0; t < Nt_; ++t )
    {
      os << t;
      for( int f = 0; f < SummaryNames.size(); f++ )
        os << Space << SummaryData( f, t );
      if( t != Nt_ - 1 )
        os << NewLine;
    }
  }
  else
  {
    // A single row, showing the first summary field for each column
    for( int t = 0; t < Nt_; t++ )
    {
      if( t )
        os << Space;
      os << SummaryData( t );
    }
  }
}

template class Sample<double>;
template class Sample<float>;
template class Sample<std::complex<double>>;
template class Sample<std::complex<float>>;

template <typename T>
const std::string &Fold<T>::DefaultGroupName()
{
  return sFold;
}

template <typename T>
bool Fold<T>::bFolded()
{
  return true;
}

template <typename T>
void Fold<T>::SummaryComments( std::ostream &s, bool bVerboseSummary ) const
{
  Base::SummaryComments( s, bVerboseSummary );
  if( NtUnfolded ) s << "# NtUnfolded: " << NtUnfolded << NewLine;
  if( reality != Reality::Unknown ) s << "# Reality: " << reality << NewLine;
  if( parity != Parity::Unknown ) s << "# Parity: " << parity << NewLine;
  if( sign != Sign::Unknown ) s << "# Sign: " << sign << NewLine;
  if( t0Negated ) s << "# timeslice 0 negated: true" << NewLine;
  if( Conjugated ) s << "# conjugate operator average: true" << NewLine;
}

template <typename T>
void Fold<T>::ReadAttributes( ::H5::Group &g )
{
  Base::ReadAttributes( g );
  ::H5::Attribute a;
  NtUnfolded = 0;
  try
  {
    a = g.openAttribute(sNtUnfolded);
    a.read( ::H5::PredType::NATIVE_INT, &NtUnfolded );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  reality = Reality::Unknown;
  try
  {
    a = g.openAttribute(sReality);
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    if( EqualIgnoreCase( s, Reality_TextEquiv_Real ) )
      reality = Reality::Real;
    else if( EqualIgnoreCase( s, Reality_TextEquiv_Imaginary ) )
      reality = Reality::Imag;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  parity = Parity::Unknown;
  try
  {
    a = g.openAttribute(sParity);
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    if( EqualIgnoreCase( s, Parity_TextEquiv_Even ) )
      parity = Parity::Even;
    else if( EqualIgnoreCase( s, Parity_TextEquiv_Odd ) )
      parity = Parity::Odd;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  sign = Sign::Unknown;
  try
  {
    a = g.openAttribute(sSign);
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    if( EqualIgnoreCase( s, Sign_TextEquiv_Positive ) )
      sign = Sign::Positive;
    else if( EqualIgnoreCase( s, Sign_TextEquiv_Negative ) )
      sign = Sign::Negative;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  t0Negated = false;
  std::int8_t i8;
  try
  {
    a = g.openAttribute(st0Negated);
    a.read( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    if( i8 )
      t0Negated = true;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  Conjugated = false;
  try
  {
    a = g.openAttribute(sConjugated);
    a.read( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    if( i8 )
      Conjugated = true;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  BootstrapList.clear();
  try
  {
    a = g.openAttribute(sBootstrapList);
    BootstrapList = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
}

template <typename T>
int Fold<T>::WriteAttributes( ::H5::Group &g ) const
{
  int iReturn = Sample<T>::WriteAttributes( g );
  const hsize_t OneDimension{ 1 };
  ::H5::DataSpace ds1( 1, &OneDimension );
  ::H5::Attribute a;
  if( NtUnfolded )
  {
    a = g.createAttribute( sNtUnfolded, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &NtUnfolded );
    a.close();
    iReturn++;
  }
  if( reality == Reality::Real || reality == Reality::Imag )
  {
    const std::string &s{reality==Reality::Real?Reality_TextEquiv_Real:Reality_TextEquiv_Imaginary};
    a = g.createAttribute( sReality, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, s );
    a.close();
    iReturn++;
  }
  if( parity == Parity::Even || parity == Parity::Odd )
  {
    const std::string &s{ parity == Parity::Even ? Parity_TextEquiv_Even : Parity_TextEquiv_Odd };
    a = g.createAttribute( sParity, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, s );
    a.close();
    iReturn++;
  }
  if( sign == Sign::Positive || sign == Sign::Negative )
  {
    const std::string &s{ sign == Sign::Positive ? Sign_TextEquiv_Positive : Sign_TextEquiv_Negative};
    a = g.createAttribute( sSign, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, s );
    a.close();
    iReturn++;
  }
  const std::int8_t i8{ 1 };
  if( t0Negated )
  {
    a = g.createAttribute( st0Negated, ::H5::PredType::STD_U8LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    iReturn++;
  }
  if( Conjugated )
  {
    a = g.createAttribute( sConjugated, ::H5::PredType::STD_U8LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    iReturn++;
  }
  if( BootstrapList.size() )
  {
    H5::WriteAttribute( g, sBootstrapList, BootstrapList );
    iReturn++;
  }
  return iReturn;
}

template class Fold<double>;
template class Fold<float>;
template class Fold<std::complex<double>>;
template class Fold<std::complex<float>>;

template <typename T>
UniqueNameSet Model<T>::GetStatColumnNames() const
{
  UniqueNameSet Names;
  const std::size_t NumParams{ params.NumScalars( Common::Param::Type::All ) };
  const std::size_t NumColumns{ Base::ColumnNames.size() };
  for( std::size_t i = NumParams; i < NumColumns; ++i )
    Names.insert( Base::ColumnNames[i] );
  return Names;
}

template <typename T>
Model<T>::Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns )
: Base::Sample( NumSamples, static_cast<int>( Params_.NumScalars( Param::Type::All ) + ExtraColumns.size() ) ),
  params{Params_}
{
  CommonConstruct( ExtraColumns );
}

template <typename T>
Model<T>::Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns,
                 int CovarSampleSize_, bool CovarFrozen_, SampleSource CovarSource_, std::vector<int> CovarRebin_, int CovarNumBoot_ )
: Base::Sample( NumSamples, static_cast<int>( Params_.NumScalars( Param::Type::All ) + ExtraColumns.size() ) ),
  params{Params_}, CovarSampleSize{CovarSampleSize_}, CovarFrozen{CovarFrozen_},
  CovarSource{CovarSource_}, CovarRebin{CovarRebin_}, CovarNumBoot{CovarNumBoot_}
{
  CommonConstruct( ExtraColumns );
}

template <typename T>
void Model<T>::ReadAttributes( ::H5::Group &g )
{
  Base::ReadAttributes( g );
  ::H5::Attribute a;
  a = g.openAttribute(sDoF);
  a.read( ::H5::PredType::NATIVE_INT, &dof );
  a.close();
  std::int8_t i8;
  a = g.openAttribute(sCovarFrozen);
  a.read( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  CovarFrozen = ( i8 != 0 );
  try
  {
    a = g.openAttribute( sNumExponents );
    a.read( ::H5::PredType::NATIVE_INT, &NumExponents );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadVector( g, sGuess+s_C, Guess );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceIn+s_C, CovarIn );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovariance+s_C, Covar );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelation+s_C, Correl );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationCholesky+s_C, CorrelCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceInv+s_C, CovarInv );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationInv+s_C, CorrelInv );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationInvCholesky+s_C, CorrelInvCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceInvCholesky+s_C, CovarInvCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    StdErrorMean.Read( g, sStdErrorMean );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    FitInput.Read( g, sFitInput );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    ModelPrediction.Read( g, sModelPrediction );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    ErrorScaled.Read( g, sErrorScaled );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute( sCovarSource );
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    std::stringstream ss( s );
    ss >> CovarSource;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
    CovarSource = SS::Binned; // Default if missing
  }
  CovarRebin.clear();
  try
  {
    a = g.openAttribute( sCovarRebin );
    ::H5::DataSpace dsp = a.getSpace();
    const int rank{ dsp.getSimpleExtentNdims() };
    if( rank != 1 )
      throw std::runtime_error( sCovarRebin + " dimensions " + std::to_string( rank ) + ", expecting 1" );
    hsize_t Num;
    dsp.getSimpleExtentDims( &Num );
    if( Num > std::numeric_limits<int>::max() )
      throw std::runtime_error( sCovarRebin + " too many items " + std::to_string( Num ) );
    std::vector<int> Buffer( Num );
    a.read( ::H5::PredType::NATIVE_INT, &Buffer[0] );
    a.close();
    CovarRebin = std::move( Buffer );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute(sCovarSampleSize);
    a.read( ::H5::PredType::NATIVE_INT, &CovarSampleSize );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    CovarSampleSize = this->SampleSize;
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute(sCovarNumBoot);
    a.read( ::H5::PredType::NATIVE_INT, &CovarNumBoot );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    CovarNumBoot = 0;
    ::H5::Exception::clearErrorStack();
  }
  ModelType.clear();
  try
  {
    a = g.openAttribute( sModelTypeList );
    ModelType = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  ModelArgs.clear();
  try
  {
    a = g.openAttribute( sModelArgsList );
    ModelArgs = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  ::H5::Group gParam;
  if( H5::OpenOptional( gParam, g, sParam ) )
  {
    OldFormatNumExponents = std::numeric_limits<int>::min();
    // Read new format containing parameters
    gParam.close();
    params.ReadH5( g, sParam );
    // Get the fit times
    std::size_t NumFitTimes;
    std::string sAttName{ sFitTime };
    const std::size_t Len{ sAttName.length() };
    a = g.openAttribute( sAttName );
    a.read( H5::Equiv<std::size_t>::Type, &NumFitTimes );
    a.close();
    FitTimes.resize( NumFitTimes );
    for( std::size_t i = 0; i < NumFitTimes; ++i )
    {
      sAttName.resize( Len );
      sAttName.append( std::to_string( i ) );
      a = g.openAttribute( sAttName );
      ::H5::DataSpace dsp = a.getSpace();
      const int nDims{ dsp.getSimpleExtentNdims() };
      bool bError{ true };
      if( nDims == 1 )
      {
        hsize_t Dim;
        dsp.getSimpleExtentDims( &Dim );
        if( Dim > 0 && Dim <= std::numeric_limits<std::size_t>::max() )
        {
          FitTimes[i].resize( static_cast<std::size_t>( Dim ) );
          a.read( ::H5::PredType::NATIVE_INT, &FitTimes[i][0] );
          bError = false;
        }
      }
      if( bError )
        throw std::runtime_error( "Error reading attribute " + sAttName  );
      a.close();
    }
  }
  else
  {
    // Read old format
    a = g.openAttribute(sNumExponents);
    a.read( ::H5::PredType::NATIVE_INT, &OldFormatNumExponents );
    a.close();
    int NumFiles;
    a = g.openAttribute(sNumFiles);
    a.read( ::H5::PredType::NATIVE_INT, &NumFiles );
    a.close();
    int ti;
    a = g.openAttribute(sTI);
    a.read( ::H5::PredType::NATIVE_INT, &ti );
    a.close();
    int tf;
    a = g.openAttribute(sTF);
    a.read( ::H5::PredType::NATIVE_INT, &tf );
    a.close();
    a = g.openAttribute(sOperators);
    OldFormatOpNames = H5::ReadStrings( a );
    a.close();
    // Recreate FitTimes
    const int Extent{ tf - ti + 1 };
    if( Extent <= 0 )
      throw std::runtime_error( "TI " + std::to_string( ti ) + " > TF " + std::to_string( tf ) );
    FitTimes.resize( NumFiles );
    for( std::vector<int> &v : FitTimes )
    {
      v.resize( Extent );
      for( int i = 0; i < Extent; ++i )
        v[i] = ti + i;
    }
  }
}

template <typename T>
void Model<T>::ValidateAttributes()
{
  if( OldFormatNumExponents != std::numeric_limits<int>::min() )
  {
    // Validate old format
    const int NumOps{ static_cast<int>( OldFormatOpNames.size() ) + 1 };
    const int NumExpected{ OldFormatNumExponents * NumOps + 1 };
    if( Base::Nt_ != NumExpected )
    {
      std::ostringstream s;
      s << "Have " << Base::Nt_ << " parameters, but expected " << NumExpected << ", i.e. " << OldFormatNumExponents
      << " exponents * ( " << ( NumOps - 1 ) << " operators + 1 energy ) + chi squared per degree of freedom";
      throw std::runtime_error( s.str().c_str() );
    }
    // Now build parameters
    params.clear();
    Param::Key k;
    k.Object.push_back( this->Name_.Base );
    std::size_t Len{ k.Object[0].find_first_of( '.' ) };
    if( Len != std::string::npos )
      k.Object[0].resize( Len );
    k.Name = EnergyPrefix;
    params.Add( k, OldFormatNumExponents, true, Param::Type::Variable );
    for( const std::string &s : OldFormatOpNames )
    {
      k.Name = s;
      params.Add( k, OldFormatNumExponents, false, Param::Type::Variable );
    }
    params.AssignOffsets();
    OldFormatOpNames.clear();
    // If there's more than one exponent, data need reordering
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Data[0].GetCentral() );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Data[0].GetReplicaMean() );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Data[0].Replica );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Binned[0] );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Raw );
  }
  Base::ValidateAttributes();
  if( Base::Nt_ < NumParams() )
  {
    std::ostringstream s;
    s << "Model has " << NumParams() << " parameters, but " << Base::Nt_ << " columns";
    throw std::runtime_error( s.str().c_str() );
  }
}

template <typename T>
int Model<T>::WriteAttributes( ::H5::Group &g ) const
{
  int iReturn = Base::WriteAttributes( g ) + 4;
  const hsize_t OneDimension{ 1 };
  ::H5::DataSpace ds1( 1, &OneDimension );
  ::H5::Attribute a;
  a = g.createAttribute( sNumExponents, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &NumExponents );
  a.close();
  a = g.createAttribute( sDoF, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &dof );
  a.close();
  std::int8_t i8{ static_cast<std::int8_t>( CovarFrozen ? 1 : 0 ) };
  a = g.createAttribute( sCovarFrozen, ::H5::PredType::STD_U8LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  H5::WriteVector( g, sGuess+s_C, Guess );
  H5::WriteMatrix( g, sCovarianceIn+s_C, CovarIn );
  H5::WriteMatrix( g, sCovariance+s_C, Covar );
  H5::WriteMatrix( g, sCorrelation+s_C, Correl );
  H5::WriteMatrix( g, sCorrelationCholesky+s_C, CorrelCholesky );
  H5::WriteMatrix( g, sCovarianceInv+s_C, CovarInv );
  H5::WriteMatrix( g, sCorrelationInv+s_C, CorrelInv );
  H5::WriteMatrix( g, sCorrelationInvCholesky+s_C, CorrelInvCholesky );
  H5::WriteMatrix( g, sCovarianceInvCholesky+s_C, CovarInvCholesky );
  StdErrorMean.Write( g, sStdErrorMean );
  FitInput.Write( g, sFitInput );
  ModelPrediction.Write( g, sModelPrediction );
  ErrorScaled.Write( g, sErrorScaled );
  {
    std::ostringstream os;
    os << CovarSource;
    a = g.createAttribute( sCovarSource, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, os.str() );
    a.close();
    iReturn++;
  }
  if( !CovarRebin.empty() )
  {
    hsize_t Dims[1] = { CovarRebin.size() };
    ::H5::DataSpace dsp( 1, Dims );
    a = g.createAttribute( sCovarRebin, ::H5::PredType::NATIVE_INT, dsp );
    a.write( ::H5::PredType::NATIVE_INT, &CovarRebin[0] );
    a.close();
    dsp.close();
    iReturn++;
  }
  a = g.createAttribute( sCovarSampleSize, ::H5::PredType::STD_U32LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &CovarSampleSize );
  a.close();
  if( CovarNumBoot )
  {
    a = g.createAttribute( sCovarNumBoot, ::H5::PredType::STD_U32LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &CovarNumBoot );
    a.close();
    iReturn++;
  }
  if( ModelType.size() )
    H5::WriteAttribute( g, sModelTypeList, ModelType );
  if( ModelArgs.size() )
    H5::WriteAttribute( g, sModelArgsList, ModelArgs );
  params.WriteH5( g, sParam );
  // Write the fit times
  std::size_t NumFitTimes{ FitTimes.size() };
  std::string sAttName{ sFitTime };
  const std::size_t Len{ sAttName.length() };
  a = g.createAttribute( sAttName, ::H5::PredType::STD_U32LE, ds1 );
  a.write( H5::Equiv<std::size_t>::Type, &NumFitTimes );
  a.close();
  for( std::size_t i = 0; i < NumFitTimes; ++i )
  {
    sAttName.resize( Len );
    sAttName.append( std::to_string( i ) );
    hsize_t Dims[1] = { FitTimes[i].size() };
    ::H5::DataSpace dsp( 1, Dims );
    a = g.createAttribute( sAttName, ::H5::PredType::STD_U16LE, dsp );
    a.write( ::H5::PredType::NATIVE_INT, &FitTimes[i][0] );
    a.close();
    dsp.close();
    iReturn++;
  }
  return iReturn;
}

template <typename T>
bool Model<T>::CheckParameters( int Strictness, scalar_type MonotonicUpperLimit )
{
  bool bResult{ !this->SummaryNames.empty() };
  if( bResult && ( Strictness >= 0 || MonotonicUpperLimit < std::numeric_limits<scalar_type>::max() ) )
  {
    const bool bVeryStrictNonZero{ ( Strictness & 1 ) != 0 };
    const bool bVeryStrictDifferent{ ( Strictness & 2 ) != 0 };
    using TypeVE = Common::ValWithEr<scalar_type>;
    scalar_type TypeVE::* const pLowZ { bVeryStrictNonZero ? &TypeVE::Min : &TypeVE::Low };
    scalar_type TypeVE::* const pHighZ{ bVeryStrictNonZero ? &TypeVE::Max : &TypeVE::High };
    scalar_type TypeVE::* const pLowD { bVeryStrictDifferent ? &TypeVE::Min : &TypeVE::Low };
    scalar_type TypeVE::* const pHighD{ bVeryStrictDifferent ? &TypeVE::Max : &TypeVE::High };
    for( const Params::value_type &it : params )
    {
      //const Param::Key &pk{ it.first };
      const Param &p{ it.second };
      if( p.type == Param::Type::Variable )
      {
        const std::size_t o{ p.GetOffset( 0, Param::Type::All ) };
        for( std::size_t i = 0; i < p.size; ++i )
        {
          TypeVE &VEoi{ Base::SummaryData( static_cast<int>( o + i ) ) };
          // All elements must be statistically different from zero
          bool bOK{  ( VEoi.*pLowZ > 0 && VEoi.*pHighZ > 0 )
                  || ( VEoi.*pLowZ < 0 && VEoi.*pHighZ < 0 ) };
          // Monotonic parameters (energies) must be < limit
          if( bOK && p.bMonotonic && VEoi.Central > MonotonicUpperLimit )
            bOK = false;
          // All elements must be different from each other
          for( std::size_t j = 0; bOK && j < i; ++j )
          {
            TypeVE &VEoj{ Base::SummaryData( static_cast<int>( o + j ) ) };
            bOK = VEoi.*pLowD > VEoj.*pHighD || VEoi.*pHighD < VEoj.*pLowD;
          }
          // Save results
          VEoi.Check = bOK ? 1 : 0;
          if( !bOK )
            bResult = false;
        }
      }
    }
  }
  return bResult;
}

template <typename T>
void Model<T>::SummaryComments( std::ostream & s, bool bVerboseSummary ) const
{
  Base::SummaryComments( s, true );
  for( std::size_t i = 0; i < FitTimes.size(); ++i )
  {
    s << "# Fit Times " << i << Colon;
    for( const int t : FitTimes[i] )
      s << Space << t;
    s << NewLine;
  }
  for( std::size_t i = 0; i < ModelType.size(); ++i )
    s << "# ModelType " << i << ": " << ModelType[i] << NewLine;
  for( std::size_t i = 0; i < ModelArgs.size(); ++i )
    s << "# ModelArgs " << i << ": " << ModelArgs[i] << NewLine;
  s << "# NumExponents: " << NumExponents << NewLine;
  if( !params.empty() )
  {
    s << "# Params: ";
    bool bFirst{ true };
    for( const Params::value_type &p : params )
    {
      if( bFirst )
        bFirst = false;
      else
        s << CommaSpace;
      s << p;
    }
    s << NewLine;
  }
  s << "# Frozen covariance: " << CovarFrozen << NewLine
    << "# Covariance source: " << CovarSource << NewLine;
  if( !CovarRebin.empty() )
  {
    s << "# Covariance rebin:";
    for( auto i : CovarRebin )
      s << Space << i;
    s << NewLine;
  }
  // Now write info about statistics
  if( dof < 1 )
    s << "# Extrapolation : All p-values=1\n";
  // Chi^2 per dof
  int iCol{ this->GetColumnIndexNoThrow( sChiSqPerDof ) };
  if( iCol >= 0 )
    s << "# " << sChiSqPerDof << " : " << (*this)(Model<T>::idxCentral,iCol)
      << " (dof=" << dof << ")" << NewLine;
  // Hotelling pValue
  s << "# " << sPValueH << " : Hotelling (p=" << dof << ", m-1=" << ( CovarSampleSize - 1 ) << ")";
  bool bHotellingUsable{ Common::HotellingDist::Usable( dof, CovarSampleSize - 1 ) };
  if( !bHotellingUsable )
    s << ". Not usable - set to ChiSquared p-value instead";
  else
  {
    iCol = this->GetColumnIndexNoThrow( sPValueH );
    if( iCol >= 0 )
      s << " = " << (*this)(Model<T>::idxCentral,iCol);
  }
  s << NewLine;
  // Chi^2 pValue
  iCol = this->GetColumnIndexNoThrow( sPValue );
  if( iCol >= 0 )
    s << "# " << sPValue << " : Chi^2 = " << (*this)(Model<T>::idxCentral,iCol) << NewLine;
}

const char ModelBase::SummaryColumnPrefix[] = "ti tf tiLabel tfLabel NumDataPoints dof SampleSize CovarSampleSize";

template <typename T>
void Model<T>::SummaryColumnNames( std::ostream &os ) const
{
  os << SummaryColumnPrefix << Common::Space;
  Base::SummaryColumnNames( os );
  for( std::size_t i = 1; i < FitTimes.size(); ++i )
    os << Space << "ti" << i << Space << "tf" << i;
}

template <typename T>
void Model<T>::SummaryContentsPrefix( std::ostream &os ) const
{
  if( FitTimes.size() )
  {
    os << FitTimes[0][0] << Space << FitTimes[0].back();
    // Write tiLabel
    os << Space << FitTimes[0][0];
    for( std::size_t i = 1; i < FitTimes.size(); ++i )
      os << Underscore << FitTimes[i][0] << Underscore << FitTimes[i].back();
    // Write tfLabel
    os << Space << FitTimes[0].back();
    for( std::size_t i = 1; i < FitTimes.size(); ++i )
      os << Underscore << FitTimes[i][0] << Underscore << FitTimes[i].back();
  }
  else
    os << "0 0 0 0";
  os << Space << GetExtent();
  os << Space << dof;
  os << Space << this->SampleSize;
  os << Space << CovarSampleSize;
}

template <typename T>
void Model<T>::SummaryContentsSuffix( std::ostream &os ) const
{
  // Write all the fit times
  for( std::size_t i = 1; i < FitTimes.size(); ++i )
    os << Space << FitTimes[i][0] << Space << FitTimes[i].back();
  os << NewLine;
}

template <typename T>
void Model<T>::SummaryContents( std::ostream &os ) const
{
  SummaryContentsPrefix( os );
  os << Space;
  Base::SummaryContents( os );
  SummaryContentsSuffix( os );
}

template <typename T> std::vector<ValWithEr<typename Model<T>::scalar_type>>
Model<T>::GetValWithEr( const Params &ParamNames, const UniqueNameSet &StatNames ) const
{
  using Scalar = typename Model<T>::scalar_type;
  const ValWithEr<Scalar> Zero(0,0,0,0,0,0);
  const std::size_t NumScalars{ ParamNames.NumScalars( Param::Type::All ) };
  const std::size_t NumStats{ StatNames.size() };
  std::vector<ValWithEr<Scalar>> v;
  v.reserve( NumScalars + NumStats );
  for( const Params::value_type &it : ParamNames )
  {
    const Param::Key &k{ it.first };
    const Param &p{ it.second };
    const std::size_t NumToWrite{ p.size };
    Params::const_iterator itMe{ params.find( k ) };
    const std::size_t NumHave{ itMe == params.cend() ? 0 : itMe->second.size };
    // Write out the parameters I have
    if( NumHave )
    {
      const Param &pMe{ itMe->second };
      std::size_t MyOffset{ pMe.GetOffset( 0, Param::Type::All ) };
      for( std::size_t i = 0; i < std::min( NumHave, NumToWrite ); ++i )
        v.emplace_back( Base::SummaryData( static_cast<int>( MyOffset + i ) ) );
    }
    // Now write out dummy values for those I don't have
    for( std::size_t i = NumHave; i < NumToWrite; ++i )
      v.emplace_back( Zero );
  }
  // Now write out all the statistics columns
  for( const std::string &StatName : StatNames )
  {
    const int idx{ Base::GetColumnIndexNoThrow( StatName ) };
    v.emplace_back( idx < 0 ? Zero : Base::SummaryData( idx ) );
  }
  return v;
}

template <typename T>
void Model<T>::WriteSummaryTD( const DataSet<T> &ds, const std::string &sOutFileName,
                               bool bVerboseSummary )
{
  const ValWithEr<T> veZero( 0, 0, 0, 0, 0, 0 );
  using Scalar = typename is_complex<T>::Scalar;
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  std::ofstream ofs( sOutFileName );
  SummaryHeader<T>( ofs, sOutFileName );
  SummaryComments( ofs, bVerboseSummary );
  // Write column names
  static constexpr int idxData{ 0 };
  //static constexpr int idxTheory{ 1 };
  ofs << "# Correlator data used in fit\n";
  ofs << "field seq fitseq model t tfit ";
  ValWithEr<T>::Header( "data", ofs );
  ofs << Space;
  ValWithEr<T>::Header( "theory", ofs );
  ofs << NewLine;
  // Grab the model data and theory with bootstrapped errors
  const int Extent{ GetExtent() };
  if( !Extent )
    return;
  std::array<VectorView<T>, 2> vv; // View into data and theory replica data. Used to populate Value
  std::vector<T> Buffer; // Scratch buffer for ValWithEr<T>
  std::array<std::vector<ValWithEr<T>>,2> Value; // Data and theory values with errors
  for( std::size_t i = 0; i < Value.size(); ++i )
  {
    JackBoot<T> &ThD{ i == idxData ? FitInput : ModelPrediction };
    ThD.MakeStatistics( Value[i] );
  }
  // Write theory and data CORRELATOR values, with each data point in fit on a separate line
  int Seq = 0;
  int FitSeq = 0;
  std::vector<std::vector<int>> SortedTimes{ FitTimes };
  std::ostringstream osBuffer;
  osBuffer << std::boolalpha << std::setprecision(std::numeric_limits<scalar_type>::digits10+2);
  for( std::size_t m = 0; m < FitTimes.size(); ++m )
  {
    std::sort( SortedTimes[m].begin(), SortedTimes[m].end() );
    for( int t = 0, idxSort = 0; t < ds.corr[m].Nt(); ++t )
    {
      const bool bInFit{ t == SortedTimes[m][idxSort] };
      std::ostream &os{ bInFit ? dynamic_cast<std::ostream &>( ofs ) : osBuffer };
      os << "corr " << Seq << Space << ( bInFit ? FitSeq : -1 ) << Space << m
         << Space << t << Space << ( bInFit ? t : -1 );
      if( bInFit )
      {
        for( std::size_t i = 0; i < Value.size(); ++i )
          os << Space << Value[i][FitSeq];
        ++idxSort;
        ++FitSeq;
      }
      else
        os << Space << ds.corr[m].SummaryData( t ) << Space << veZero;
      os << NewLine;
      ++Seq;
    }
  }
  if( !osBuffer.str().empty() )
  {
    ofs << "\n\n# Correlator data not used in fit\n" << osBuffer.str();
    osBuffer.str( "" );
  }
  // Write theory and data EFFECTIVE MASS, with each data point in fit on a separate line
  ofs << "\n\n# Effective masses used in fit\n";
  Seq = 0;
  FitSeq = 0;
  ValWithEr<T> v;
  for( std::size_t m = 0; m < FitTimes.size(); ++m )
  {
    ++Seq;    // Skip past t=0
    ++FitSeq; // Skip past t=0
    for( int t = 1, idxSort = 1; t < ds.corr[m].Nt(); ++t )
    {
      // Effective masses are for the point half-way between the measurements
      // Keeping track of 2 * the effective timeslice allows us to use integer arithmetic
      const int TwoT{ 2 * t - 1 }; // Average of this and previous timeslice
      const int TwoTFit{ idxSort >= FitTimes[m].size() ? std::numeric_limits<int>::max()
                                              : FitTimes[m][idxSort] + FitTimes[m][idxSort - 1] };
      const bool bShowFit{ TwoTFit <= TwoT };
      if( bShowFit )
      {
        // Show the fit
        const double dblT{ 0.5 * TwoTFit };
        const int DeltaT{ FitTimes[m][idxSort] - FitTimes[m][idxSort - 1] };
        const Scalar DeltaTInv{ static_cast<Scalar>( DeltaT == 1 ? 1. : ( 1. / DeltaT ) ) };
        ofs << "log " << Seq << Space << FitSeq << Space << m << Space << dblT << Space << dblT;
        for( std::size_t i = 0; i < Value.size(); ++i )
        {
          JackBoot<T> &ThD{ i == idxData ? FitInput : ModelPrediction };
          Buffer.resize( ThD.NumReplicas() ); // Shouldn't vary, but must be right size for check
          int bootCount{ 0 };
          for( int bootrep = 0; bootrep < ThD.NumReplicas(); ++bootrep )
          {
            Buffer[bootCount] = std::log( ThD(bootrep,FitSeq-1) / ThD(bootrep,FitSeq) ) * DeltaTInv;
            if( IsFinite( Buffer[bootCount] ) )
              ++bootCount;
          }
          v.Get( std::log( ThD(JackBoot<T>::idxCentral,FitSeq - 1)
                          / ThD(JackBoot<T>::idxCentral,FitSeq) ) * DeltaTInv, Buffer, bootCount,
                ThD.Seed == SeedWildcard );
          ofs << Space << v;
        }
        ofs << NewLine;
        ++idxSort;
        ++FitSeq;
      }
      if( TwoT != TwoTFit )
      {
        // Show the raw data
        const double dblT{ 0.5 * TwoT };
        osBuffer << "log " << Seq;
        if( bShowFit )
          osBuffer << ".5"; // Fit and data points out by 0.5 - make sure seq is different
        osBuffer << Space << -1 << Space << m << Space << dblT << Space << -1
                 << Space << ds.corr[m].SummaryData( 2, t ) << Space << veZero << NewLine;
      }
      ++Seq;
    }
  }
  if( !osBuffer.str().empty() )
  {
    ofs << "\n\n# Effective masses not used in fit\n" << osBuffer.str();
    osBuffer.str( "" );
  }
}

template <typename T>
void Model<T>::ReorderOldFormat( int NumOps, int NumExponents, Vector<T> &v )
{
  if( NumExponents > 1 && NumOps > 1 )
  {
    std::vector<T> Buffer( NumOps * NumExponents );
    if( v.size < Buffer.size() )
      throw std::runtime_error( "Model<T>::ReorderOldFormat( Vector<T> & )" );
    for( int e = 0; e < NumExponents; ++e )
      for( int o = 0; o < NumOps; ++o )
        Buffer[o * NumExponents + e] = v[e * NumOps + o];
    for( std::size_t j = 0; j < Buffer.size(); ++j )
      v[j] = Buffer[j];
  }
}

template <typename T>
void Model<T>::ReorderOldFormat( int NumOps, int NumExponents, Matrix<T> &m )
{
  if( NumExponents > 1 && NumOps > 1 )
  {
    std::vector<T> Buffer( NumOps * NumExponents );
    if( m.size2 < Buffer.size() )
      throw std::runtime_error( "Model<T>::ReorderOldFormat( Matrix<T> & )" );
    for( std::size_t i = 0; i < m.size1; ++i )
    {
      for( int e = 0; e < NumExponents; ++e )
        for( int o = 0; o < NumOps; ++o )
          Buffer[o * NumExponents + e] = m(i, e * NumOps + o);
      for( std::size_t j = 0; j < Buffer.size(); ++j )
        m(i,j) = Buffer[j];
    }
  }
}

template <typename T> void Model<T>::CommonConstruct( const std::vector<std::string> &ExtraColumns )
{
  std::vector<std::string> Cols{ params.GetNames( Param::Type::All, false ) };
  Cols.reserve( Cols.size() + ExtraColumns.size() );
  for( const std::string &s : ExtraColumns )
    Cols.push_back( s );
  this->SetColumnNames( Cols );
}

template class Model<double>;
template class Model<float>;
template class Model<std::complex<double>>;
template class Model<std::complex<float>>;

template <typename T>
void DataSet<T>::clear()
{
  NSamples = 0;
  Extent = 0;
  corr.clear();
  FitTimes.clear();
  constFile.clear();
  constMap.clear();
}

// Specify which times I'm fitting to, as a list of timeslices for each correlator
template <typename T>
void DataSet<T>::SetFitTimes( const std::vector<std::vector<int>> &fitTimes_ )
{
  if( fitTimes_.size() != corr.size() )
    throw std::runtime_error( std::to_string( fitTimes_.size() ) + " FitTimes but "
                             + std::to_string( corr.size() ) + " correlators" );
  std::vector<std::vector<int>> ft{ fitTimes_ };
  for( int i = 0; i < ft.size(); ++i )
  {
    if( !ft[i].empty() )
    {
      std::sort( ft[i].begin(), ft[i].end() );
      if( ft[i][0] < 0 || ft[i].back() >= corr[i].Nt()
         || std::adjacent_find( ft[i].begin(), ft[i].end() ) != ft[i].end() )
        throw std::runtime_error( "FitTimes[" + std::to_string( i ) + "]=[" + std::to_string( ft[i][0] )
                                 + "..." + std::to_string( ft[i].back() ) + "] invalid" );
    }
  }
  SetValidatedFitTimes( std::move( ft ) );
}

template <typename T>
void DataSet<T>::SetFitTimes( int tMin, int tMax )
{
  const int extent_{ tMax - tMin + 1 };
  std::vector<std::vector<int>> ft{ corr.size() };
  for( int i = 0; i < corr.size(); ++i )
  {
    if( tMin < 0 || extent_ < 1 || tMax >= corr[i].Nt() )
      throw std::runtime_error("Fit range ["+std::to_string(tMin)+", "+std::to_string(tMax)+"] invalid");
    ft[i].reserve( extent_ );
    for( int j = tMin; j <= tMax; ++j )
      ft[i].push_back( j );
  }
  SetValidatedFitTimes( std::move( ft ) );
}

// Cache the data for this data set
template <typename T>
void DataSet<T>::SetValidatedFitTimes( std::vector<std::vector<int>> &&FitTimes_ )
{
  std::size_t Extent_{ GetExtent( FitTimes_ ) };
  if( Extent_ == 0 )
    throw std::runtime_error( "Fit range empty" );
  if( Extent_ > std::numeric_limits<int>::max() )
    throw std::runtime_error( "Fit range stupidly big" );
  Extent = static_cast<int>( Extent_ );
  FitTimes = std::move( FitTimes_ );
  CacheRawData();
}

template <typename T>
void DataSet<T>::CacheRawData()
{
  SampleSource ss{ CovarSource };
  int idxJackBoot{ idxJackBootCovarSource };
  bool bFrozen{ bFrozenCovarSource };
  // How many data rows are available?
  const int NumJackBoot{ NumSamples( SS::Bootstrap, 0 ) };
  if( NumJackBoot < 1 )
    throw std::runtime_error( "DataSet<T>::CacheRawData JackBoot data unavailable" );
  // Is the requested source available
  const int NumRequestedSource{ NumSamples( ss, idxJackBoot ) };
  if( NumRequestedSource < 1 )
  {
    std::ostringstream os;
    os << "DataSet<T>::CacheRawData " << ss << '[' << idxJackBoot << "] unavailable";
    throw std::runtime_error( os.str().c_str() );
  }
  // Binned data must be available if unfrozen
  const int NumBinned{ NumSamples( SS::Binned, 0 ) };
  if( !bFrozen && NumBinned < 1 )
  {
    std::ostringstream os;
    os << "DataSet<T>::CacheRawData " << ss << '[' << idxJackBoot << "] unfrozen impossible (no binned data)";
    throw std::runtime_error( os.str().c_str() );
  }
  // Work out what to copy
  const bool FullyUnfrozen{ idxJackBoot == 0 && ss == SS::Binned && !bFrozen };
  const bool bCovarSrcData{ idxJackBoot == 0 && ss == SS::Bootstrap };
  const bool CopyCovarSrc{ !FullyUnfrozen && !bCovarSrcData };
  mFitData.resize( NumJackBoot, Extent );
  if( FullyUnfrozen )
  {
    mCovar.clear();
    mCovarCentralMean.clear();
  }
  else
  {
    mCovar.resize( Extent, Extent );
    mCovarCentralMean.resize( Extent );
  }
  if( bFrozen )
    mBinned.clear();
  else
    mBinned.resize( NumBinned, Extent );
  Matrix<T> CovarBuffer;
  if( CopyCovarSrc )
    CovarBuffer.resize( NumRequestedSource, Extent );
  // Copy the resampled data
  for( int idx = static_cast<int>( JackBootBase::idxReplicaMean ); idx < NumJackBoot; ++idx )
  {
    int dst{ 0 };
    for( int f = 0; f < corr.size(); ++f )
    {
      const JackBoot<T> &mSrc{ corr[f].getData() };
      for( int t : FitTimes[f] )
        mFitData( idx, dst++ ) = mSrc( idx, t );
    }
  }
  // Copy the binned data
  if( !bFrozen )
  {
    for( int idx = 0; idx < NumBinned; ++idx )
    {
      int dst{ 0 };
      for( int f = 0; f < corr.size(); ++f )
      {
        const Matrix<T> &mSrc{ corr[f].getBinned( 0 ) };
        for( int t : FitTimes[f] )
          mBinned( idx, dst++ ) = mSrc( idx, t );
      }
    }
  }
  // Copy the covariance source
  if( CopyCovarSrc )
  {
    for( int idx = 0; idx < NumRequestedSource; ++idx )
    {
      int dst{ 0 };
      for( int f = 0; f < corr.size(); ++f )
      {
        const Matrix<T> &mSrc{ corr[f].get( ss, idxJackBoot ) };
        for( int t : FitTimes[f] )
          CovarBuffer( idx, dst++ ) = mSrc( idx, t );
      }
    }
    // Make correlation matrix from this data
    if( ss == SampleSource::Bootstrap )
    {
      // For Bootstrap, we check environment for whether to use mean of source or mean of resamples
      const std::size_t idxMean{ JackBootBase::getCovarMeanIdx() };
      int dst{ 0 };
      for( int f = 0; f < corr.size(); ++f )
      {
        const JackBoot<T> &mSrc{ corr[f].getData( idxJackBoot ) };
        for( int t : FitTimes[f] )
          mCovarCentralMean[dst++] = mSrc( idxMean, t );
      }
    }
    else
    {
      // For binned/raw data we create the mean from the data
      JackBoot<T>::MakeMean( mCovarCentralMean, CovarBuffer );
    }
    JackBoot<T>::MakeCovar( mCovar, mCovarCentralMean, CovarBuffer, corr[0].Norm( ss ) );
    CovarBuffer.clear();
  }
  else if( bCovarSrcData )
  {
    // Make correlation matrix from this data
    mCovarCentralMean = mFitData.GetCovarMean();
    JackBoot<T>::MakeCovar( mCovar, mCovarCentralMean, mFitData.Replica, corr[0].Norm( ss ) );
  }
}

// Get the constants from the appropriate timeslice
template <typename T>
void DataSet<T>::GetFixed( int idx, Vector<T> &vResult, const std::vector<FixedParam> &Params ) const
{
  // Now copy the data into vResult
  for( const FixedParam &p : Params )
  {
    const Param &param{ GetConstantParam( p.src ) };
    std::size_t SrcIdx{ param.GetOffset( 0, Param::Type::All ) };
    for( std::size_t i = 0; i < p.Count; ++i )
      vResult[p.idx + i] = constFile[p.src.File](idx, SrcIdx + i);
  }
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// From what is already a bootstrap replica, and therefore only the central replica
/*template <typename T>
void DataSet<T>::MakeCovarFromBootstrap( SS ss, Matrix<T> &Covar ) const
{
  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, replica );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - vCentral[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += z[i] * z[j];
  }
  const T Norm{ static_cast<T>( 1 ) / static_cast<T>( m.size1 ) };
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * Norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

template <typename T>
void DataSet<T>::MakeVarFromBootstrap( SS ss, Vector<T> &Var ) const
{
  Var.resize( Extent );
  for( int i = 0; i < Extent; ++i )
    Var[i] = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, replica );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - vCentral[i];
    for( int i = 0; i < Extent; ++i )
      Var[i] += Squared( z[i] );
  }
  const T Norm{ static_cast<T>( 1 ) / static_cast<T>( m.size1 ) };
  for( int i = 0; i < Extent; ++i )
    Var[i] *= Norm;
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// From the underlying raw/(re)binned data (without bootstrapping)
template <typename T>
void DataSet<T>::MakeCovarFromNonBootstrap( int idx, SS ss, Matrix<T> &Covar ) const
{
  VectorView<T> Mean;
  VectorView<fint> Replace;
  if( idx == Fold<T>::idxCentral )
  {
    // Central replica - no replacement
    Mean = vCentral;
    Replace = RandomCentral( ss );
  }
  else
  {
    // bootstrap replica
    Mean.MapRow( mBoot, idx );
    Replace.MapRow( RandomNumbers( ss ), idx );
  }

  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, Replace[replica] );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - Mean[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += z[i] * z[j];
  }
  const T Norm{ static_cast<T>( 1 ) / ( static_cast<T>( m.size1 - 1 ) * static_cast<T>( m.size1 ) ) };
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * Norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

template <typename T>
void DataSet<T>::MakeVarFromNonBootstrap( int idx, SS ss, Vector<T> &Var ) const
{
  VectorView<T> Mean;
  VectorView<fint> Replace;
  if( idx == Fold<T>::idxCentral )
  {
    // Central replica - no replacement
    Mean = vCentral;
    Replace = RandomCentral( ss );
  }
  else
  {
    // bootstrap replica
    Mean.MapRow( mBoot, idx );
    Replace.MapRow( RandomNumbers( ss ), idx );
  }

  Var.resize( Extent );
  for( int i = 0; i < Extent; ++i )
    Var[i] = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, Replace[replica] );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - Mean[i];
    for( int i = 0; i < Extent; ++i )
      Var[i] += Squared( z[i] );
  }
  const T Norm{ static_cast<T>( 1 ) / ( static_cast<T>( m.size1 - 1 ) * static_cast<T>( m.size1 ) ) };
  for( int i = 0; i < Extent; ++i )
    Var[i] *= Norm;
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// Call the appropriate one of the previous two functions, depending on what the data are
template <typename T>
void DataSet<T>::MakeCovariance( int idx, SS ss, Matrix<T> &Covar ) const
{
  if( ss == SampleSource::Bootstrap || ( ss == SampleSource::Raw && corr[0].bRawBootstrap ) )
  {
    if( idx != Fold<T>::idxCentral )
      throw std::runtime_error( "Can only make covariance from bootstrap on central replica" );
    MakeCovarFromBootstrap( ss, Covar );
  }
  else
    MakeCovarFromNonBootstrap( idx, ss, Covar );
}

template <typename T>
void DataSet<T>::MakeVariance( int idx, SS ss, Vector<T> &Var ) const
{
  if( ss == SampleSource::Bootstrap || ( ss == SampleSource::Raw && corr[0].bRawBootstrap ) )
  {
    if( idx != Fold<T>::idxCentral )
      throw std::runtime_error( "Can only make variance from bootstrap on central replica" );
    MakeVarFromBootstrap( ss, Var );
  }
  else
    MakeVarFromNonBootstrap( idx, ss, Var );
}

// Make covariance matrix using a secondary bootstrap - Random numbers must be provided by caller
// I.e. unfrozen covariance matrix
template <typename T> void DataSet<T>::
MakeCovariance( int idx, SS ss, Matrix<T> &Covar, const MatrixView<fint> &Random ) const
{
  if( ss == SampleSource::Bootstrap || ( ss == SampleSource::Raw && corr[0].bRawBootstrap ) )
    throw std::runtime_error( "Can't rebootstrap a bootstrap" );
  const Matrix<T> &m{ Cache( ss ) };
  if( m.size1 != Random.size2() )
  {
    std::ostringstream os;
    os << ss << " data have " << m.size1 << " samples, but random numbers expect "
       << Random.size2() << " samples";
    throw std::runtime_error( os.str().c_str() );
  }
  VectorView<T> Mean;
  VectorView<fint> Replace;
  if( idx == Fold<T>::idxCentral )
  {
    // Central replica - no replacement
    Mean = vCentral;
    Replace = RandomCentral( ss );
  }
  else
  {
    // bootstrap replica
    Mean.MapRow( mBoot, idx );
    Replace.MapRow( RandomNumbers( ss ), idx );
  }

  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const T InvBootLen{ static_cast<T>( 1 ) / static_cast<T>( m.size1 ) };
  for( int replica = 0; replica < Random.size1(); ++replica )
  {
    // Perform the covariance (inner) bootstrap
    for( int i = 0; i < Extent; ++i )
      z[i] = 0;
    for( int inner = 0; inner < m.size1; ++inner )
    {
      Data.MapRow( m, Replace[ Random( replica, inner ) ] );
      for( int i = 0; i < Extent; ++i )
        z[i] += 0;
    }
    // Now contribute this to the covariance
    for( int i = 0; i < Extent; ++i )
      z[i] = z[i] * InvBootLen - Mean[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += z[i] * z[j];
  }
  const T Norm{ static_cast<T>( 1 ) / static_cast<T>( Random.size1() ) };
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * Norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}*/

// Write covariance matrix to file
template <typename T>
void DataSet<T>::SaveMatrixFile( const Matrix<T> &m, const std::string &Type, const std::string &Filename,
                                 std::vector<std::string> &Abbreviations,
                                 const std::vector<std::string> *FileComments, const char *pGnuplotExtra ) const
{
  // For now, all the matrices I write are square and match dataset, but this restriction can be safely relaxed
  if( m.size1 != Extent || m.size1 != m.size2 )
    throw std::runtime_error( Type + " matrix dimensions (" + std::to_string( m.size1 ) + Common::CommaSpace
                             + std::to_string( m.size2 ) + ") don't match dataset" );
  // Header describing what the covariance matrix is for and how to plot it with gnuplot
  std::ofstream s{ Filename };
  s << "# Matrix: " << Filename << "\n# Type: " << Type << "\n# Files: " << corr.size() << Common::NewLine;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    s << "# File" << f << ": " << corr[f].Name_.NameNoExt << "\n# Times" << f << ":";
    for( int t : FitTimes[f] )
      s << Common::Space << t;
    s << Common::NewLine;
    // Say which operators are in each file
    if( FileComments && !(*FileComments)[f].empty() )
      s << (*FileComments)[f]; // These are expected to have a trailing NewLine
    // Make default abbreviation - a letter for each file
    if( Abbreviations.size() < f )
      Abbreviations.emplace_back( 1, 'A' + f );
  }
  // Save a command which can be used to plot this file
  s << "# gnuplot: set xtics rotate noenhanced font 'Arial,4'; set ytics noenhanced font 'Arial,4';"
    << " set title noenhanced '" << Filename << "'; unset key; ";
  if( pGnuplotExtra && *pGnuplotExtra )
    s << pGnuplotExtra << "; ";
  s << "plot '" << Filename << "' matrix columnheaders rowheaders with image pixels\n";
  // Now save the column names. First column is matrix extent
  s << "# The next row contains column headers, starting with matrix extent\n" << Extent;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    for( int t : FitTimes[f] )
      s << Common::Space << Abbreviations[f] << t;
  }
  s << Common::NewLine << std::setprecision( std::numeric_limits<T>::max_digits10 );
  // Now print the actual covariance matrix
  int i{ 0 };
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    for( int t : FitTimes[f] )
    {
      s << Abbreviations[f] << t;
      for( int j = 0; j < Extent; ++j )
        s << Common::Space << m( i, j );
      s << Common::NewLine;
      ++i;
    }
  }
}

// Add a constant to my list of known constants - make sure it isn't already there
template <typename T>
void DataSet<T>::AddConstant( const Param::Key &Key, std::size_t File )
{
  static const char Invalid[] = " invalid";
  std::ostringstream os;
  os << "DataSet::AddConstant " << Key << Space;
  if( constMap.find( Key ) != constMap.end() )
  {
    os << "loaded from multiple model files";
    throw std::runtime_error( os.str().c_str() );
  }
  if( File >= constFile.size() )
  {
    os << "file " << File << Invalid;
    throw std::runtime_error( os.str().c_str() );
  }
  constMap.insert( { Key, ConstantSource( File, Key ) } );
}

template <typename T>
const Param &DataSet<T>::GetConstantParam( const ConstantSource &cs ) const
{
  Params::const_iterator it{ constFile[cs.File].params.find( cs.pKey ) };
  if( it == constFile[cs.File].params.cend() )
  {
    std::ostringstream os;
    os << "Constant not found: File " << cs.File << CommaSpace << cs.pKey;
    throw std::runtime_error( os.str().c_str() );
  }
  return it->second;
}

template <typename T>
int DataSet<T>::LoadCorrelator( Common::FileNameAtt &&FileAtt, unsigned int CompareFlags,
                                const char * PrintPrefix )
{
  if( corr.size() >= std::numeric_limits<int>::max() )
    throw std::runtime_error( "More than an integer's worth of correlators loaded!!??" );
  const int i{ static_cast<int>( corr.size() ) };
  corr.emplace_back();
  corr[i].SetName( std::move( FileAtt ) );
  corr[i].Read( PrintPrefix );
  // See whether this correlator is compatible with prior correlators
  if( i )
  {
    corr[0].IsCompatible( corr[i], &NSamples, CompareFlags );
    if( MaxSamples > corr[i].NumSamples() )
      MaxSamples = corr[i].NumSamples();
  }
  else
  {
    MaxSamples = corr[i].NumSamples();
    if( NSamples == 0 || NSamples > MaxSamples )
      NSamples = MaxSamples;
  }
  return i;
}

// Load a model file
template <typename T>
void DataSet<T>::LoadModel( Common::FileNameAtt &&FileAtt, const std::string &Args )
{
  // This is a pre-built model (i.e. the result of a previous fit)
  const std::size_t i{ constFile.size() };
  constFile.emplace_back();
  constFile[i].SetName( std::move( FileAtt ) );
  constFile[i].Read( "  " );
  // Keep track of minimum number of replicas across all files
  if( NSamples == 0 )
    NSamples = constFile[i].NumSamples();
  else if( NSamples > constFile[i].NumSamples() )
    NSamples = constFile[i].NumSamples();
  // Now see which parameters we want to read from the model
  if( Args.find_first_not_of( Common::WhiteSpace ) == std::string::npos )
  {
    // Load every constant in this file
    for( const Params::value_type &it : constFile[i].params )
      AddConstant( it.first, i );
  }
  else
  {
    // Load only those constants specifically asked for
    using ReadMap = std::map<Param::Key, Param::Key, Param::Key::Less>;
    using KVR = KeyValReader<Param::Key, Param::Key, Param::Key::Less>;
    ReadMap vThisArg{ KVR::Read( Common::ArrayFromString( Args ), &EqualSign, true ) };
    for( const ReadMap::value_type &it : vThisArg )
    {
      const Param::Key OldKey{ it.first };
      const Param::Key NewKey{ it.second };
      Params::iterator itP{ constFile[i].params.find( OldKey ) };
      if( itP == constFile[i].params.end() )
      {
        std::ostringstream es;
        es << "Parameter " << OldKey << " not found";
        throw std::runtime_error( es.str().c_str() );
      }
      AddConstant( NewKey.empty() ? OldKey : NewKey, i );
    }
  }
}

template <typename T>
std::vector<std::string> DataSet<T>::GetModelFilenames() const
{
  std::vector<std::string> v;
  for( const Model<T> &m : constFile )
    v.emplace_back( m.Name_.Filename );
  return v;
}

// Sort the operator names and renumber all the loaded correlators referring to them
template <typename T>
void DataSet<T>::SortOpNames( std::vector<std::string> &OpNames )
{
  int NumOps{ static_cast<int>( OpNames.size() ) };
  if( OpNames.size() > 1 )
  {
    // Sort the names
    UniqueNames OpSorted;
    for( int i = 0; i < NumOps; ++i )
      OpSorted.emplace( std::move( OpNames[i] ), i );
    // Extract the sorted names and indices (to renumber operators in correlator names)
    std::vector<std::string> SortedNames;
    std::vector<int> SortIndex( NumOps );
    SortedNames.reserve( NumOps );
    int idx{ 0 };
    for( UniqueNames::iterator it = OpSorted.begin(); it != OpSorted.end(); ++it )
    {
      SortedNames.emplace_back( it->first );
      SortIndex[it->second] = idx++;
    }
    // Renumber the operators and save the sorted operator names
    for( auto &f : corr )
      for( int i = 0; i < f.Name_.op.size(); ++i )
        f.Name_.op[i] = SortIndex[f.Name_.op[i]];
    OpNames = SortedNames;
  }
}

template <typename T>
int DataSet<T>::NumSamples( SampleSource ss, int idxJackBoot ) const
{
  int NSB{ corr.empty() ? 0 : corr[0].NumSamples( ss, idxJackBoot ) };
  for( std::size_t i = 1; NSB && i < corr.size(); ++i )
  {
    const int ThisNSB{ corr[i].NumSamples( ss, idxJackBoot ) };
    if( NSB > ThisNSB )
      NSB = ThisNSB;
  }
  return NSB;
}

template <typename T>
void DataSet<T>::Rebin( const std::vector<int> &NewSize )
{
  constexpr int DestJackBoot{ 1 };
  RebinSize.clear();
  RebinSize.reserve( corr.size() );
  for( std::size_t i = 0; i < corr.size(); ++i )
  {
    if( !corr[i].NumSamplesRaw() )
    {
      std::ostringstream os;
      os << "Raw samples unavailable rebinning corr " << i << CommaSpace << corr[i].Name_.Filename;
      throw std::runtime_error( os.str().c_str() );
    }
    // Rebin to the size specified (last bin size repeats to end).
    // bin size 0 => auto (i.e. all measurements on same config binned together)
    RebinSize.push_back( NewSize.empty() ? 0 : NewSize[i < NewSize.size() ? i : NewSize.size() - 1] );
    if( RebinSize.back() )
      corr[i].BinFixed( RebinSize.back(), DestJackBoot );
    else
      corr[i].BinAuto( DestJackBoot );
    if( corr[i].NumSamplesBinned( DestJackBoot ) != corr[0].NumSamplesBinned( DestJackBoot ) )
    {
      std::ostringstream os;
      os << "Rebinned corr " << i << " has " << corr[i].NumSamplesBinned( DestJackBoot )
         << " samples, others have " << corr[0].NumSamplesBinned( DestJackBoot );
      throw std::runtime_error( os.str().c_str() );
    }
  }
}

template class DataSet<double>;
template class DataSet<float>;
template class DataSet<std::complex<double>>;
template class DataSet<std::complex<float>>;

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
