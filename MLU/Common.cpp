/**

 Mike's lattice QCD utilities
 
 Source file: Common.cpp
 
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
const std::string sRawBootstrap{ "RawBootstrap" };
const std::string sTI{ "TI" };
const std::string sTF{ "TF" };
const std::string sDoF{ "DoF" };
const std::string sNumExponents{ "NumExponents" };
const std::string sNumFiles{ "NumFiles" };
const std::string sCovarFrozen{ "CovarFrozen" };
const std::string s_C{ "_C" };
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
const std::string sAuxNames{ "AuxNames" };
const std::string sSummaryDSName{ "Summary" };
const std::string sSummaryNames{ "SummaryNames" };
const std::string sColumnNames{ "ColumnNames" };
const std::string sSeed{ "Seed" };
const std::string sSeedMachine{ "SeedMachine" };
const std::string sRandom{ "Random" };
const std::string sSampleSize{ "SampleSize" };
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

const double NaN{ std::nan( "" ) };

const std::vector<std::string> DefaultModelStats{ Common::sChiSqPerDof, Common::sPValue, Common::sPValueH };

// Does the specified file exist?
bool FileExists( const std::string& Filename )
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

// Make the ancestor directories leading up to last element
// NB: Don't test whether this worked, so that it if this is done simultaneously
// by many threads, it will still work
void MakeAncestorDirs( const std::string& Filename )
{
  bool WasSlash{ true };
  std::size_t SegmentStart = 0;
  for( std::size_t i = 0; i < Filename.length(); ++i )
  {
    const bool IsSlash{ Filename[i] == '/' };
    if( IsSlash )
    {
      if( !WasSlash )
      {
        // First slash after non-slash - try to make this bit
        WasSlash = true;
        const std::size_t SegmentLen{ i - SegmentStart };
        // Don't bother with '.' or '..'
        if( SegmentLen > 2 || Filename[i - SegmentLen] != '.'
           || ( SegmentLen == 2 && Filename[i - SegmentLen + 1] != '.' ) )
        {
          // Make this part of the name
          mkdir( Filename.substr( 0, i ).c_str(), 0777 );
        }
      }
    }
    else
    {
      if( WasSlash )
      {
        // First non-slash after slash - remember where this segment starts
        WasSlash = false;
        SegmentStart = i;
      }
    }
  }
}

// Wrapper for posix gethostname()
std::string GetHostName()
{
  char Buffer[256];
  const int BufLen{ sizeof( Buffer ) - 1 };
  if( gethostname( Buffer, BufLen ) )
    throw std::runtime_error( "gethostname() returned error " + std::to_string( errno ) );
  Buffer[BufLen] = 0;
  return std::string( Buffer );
}

// Read a bootstrap replica from an HDF5 group
template <typename T>
void BootRep<T>::Read( ::H5::Group &g, const std::string &Name )
{
  H5::ReadVector( g, Name + s_C, Central );
  H5::ReadMatrix( g, Name, Replica );
}

// Write a bootstrap replica to an HDF5 group
template <typename T>
void BootRep<T>::Write( ::H5::Group &g, const std::string &Name ) const
{
  H5::WriteVector( g, Name + s_C, Central );
  H5::WriteMatrix( g, Name, Replica );
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
          Seed = FromString<unsigned int>( SeedString );
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
  Momentum pIgnore;
  if( pIgnoreMomenta )
  {
    for( const std::string &s : *pIgnoreMomenta )
    {
      while( pIgnore.Extract( Base, s ) ) // There might be more than one copy
        ; //momIgnore.emplace_back( s, pIgnore );
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
  // Now get the remaining momenta and save them
  std::smatch match;
  std::regex Pattern{ Momentum::MakeRegex( "([pP][[:alnum:]]*)" ) };
  for( int pLoop = 0; pLoop < 2; ++pLoop )
  {
    while( std::regex_search( BaseShort, match, Pattern ) )
    {
      // Extract the momentum
      const std::string sMom{ match[1] };
      std::vector<FileNameMomentum> fnp;
      fnp.reserve( 1 );
      if( pLoop == 0 )
      {
        fnp.emplace_back( sMom, std::stoi( match[2] ), std::stoi( match[3] ), std::stoi( match[4] ) );
      }
      else
      {
        fnp.emplace_back( sMom, std::stoi( match[2] ) );
      }
      const std::string sSuffix{ match.suffix() };
      BaseShort = match.prefix();
      BaseShort.append( sSuffix );
      // Now see whether we already have it
      auto it = p.find( sMom );
      if( it == p.end() )
        p.emplace( std::make_pair( sMom, fnp[0] ) );
      else if( ! ( it->second == fnp[0] ) )
      {
        std::stringstream ss;
        ss << "Repeated momentum " << sMom << CommaSpace << it->second << " != " << pIgnore;
        throw std::runtime_error( ss.str() );
      }
    }
    // Second time through the loop, we're looking for p^2
    if( pLoop == 0 )
    {
      const std::string sPattern{ "_([pP][[:alnum:]]*)" + Momentum::SquaredSuffix + "_(-?[0-9]+)" };
      Pattern = sPattern;
    }
  }
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
      MesonMom = Meson;
      auto pps = p.find( Momentum::SinkPrefix );
      if( pps != p.end() )
      {
        // We have two momenta
        MesonMom[0].append( pp->second.FileString() );
        Common::FileNameMomentum fnp{ pps->second };
        fnp.Name = Momentum::DefaultPrefix;
        MesonMom[1].append( fnp.FileString() );
      }
      else
      {
        // We only have one momentum - this goes with the lightest meson
        const bool bSourceHeavier{ QuarkWeight(BaseShortParts[2]) >= QuarkWeight(BaseShortParts[1]) };
        MesonMom[bSourceHeavier ? 1 : 0].append( pp->second.FileString() );
        Common::FileNameMomentum fnp( Momentum::DefaultPrefix, pp->second.bp2 );
        MesonMom[bSourceHeavier ? 0 : 1].append( fnp.FileString() );
      }
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
  o.reserve( NumOps );
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
  std::string s{ GetBaseShort() };
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

const FileNameMomentum &FileNameAtt::GetMomentum( const std::string &Name ) const
{
  auto it = p.find( Name );
  if( it == p.end() )
    it = p.find( Name + Momentum::SquaredSuffix );
  if( it == p.end() )
    throw std::runtime_error( "Momentum " + Name + " not available" );
  return it->second;
}

bool FileNameAtt::HasNonZeroMomentum() const
{
  for( auto it = p.begin(); it != p.end(); ++it )
    if( it->second )
      return true;
  return false;
}

const FileNameMomentum &FileNameAtt::GetFirstNonZeroMomentum() const
{
  if( p.empty() )
    throw std::runtime_error( "No momenta available" );
  for( auto it = p.begin(); it != p.end(); ++it )
    if( it->second )
      return it->second;
  return p.begin()->second;
}

void FileNameAtt::AppendMomentum( std::string &s, const FileNameMomentum &fnp, const std::string &Name ) const
{
  if( fnp.bp2 )
    s.append( fnp.p2_string( Underscore, Name ) );
  else
  {
    s.append( Underscore );
    s.append( Name );
    s.append( fnp.to_string( Underscore ) );
  }
}

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext)
{
  const char Sep = '.';
  std::string s{ Base };
  s.append( 1, Sep );
  s.append( Type );
  s.append( 1, Sep );
  s.append( std::to_string( Seed ) );
  s.append( 1, Sep );
  s.append( Ext );
  return s;
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

void GenerateRandom( std::vector<fint> &Random, SeedType Seed, std::size_t NumBoot, std::size_t NumSamples )
{
  const std::size_t RandomLen{ NumSamples * NumBoot };
  if( RandomLen / NumSamples != NumBoot || NumSamples > std::numeric_limits<fint>::max() )
    throw std::runtime_error( "Too many bootstrap replicas " + std::to_string( NumBoot )
                             + " (with " + std::to_string( NumSamples ) + " per replica)" );
  Random.resize( RandomLen );
  std::mt19937                        engine( Seed );
  std::uniform_int_distribution<fint> random( 0, static_cast<fint>( NumSamples - 1 ) );
  for( std::size_t i = 0; i < RandomLen; ++i )
    Random[i] = random( engine );
}

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

template <typename T> const std::string Model<T>::EnergyPrefix{ "E" };

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
  if( !OldFormatOpNames.empty() )
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
    ReorderOldFormat( NumOps, OldFormatNumExponents, this->m_pData, this->NumSamples_ + this->NumExtraSamples );
    ReorderOldFormat( NumOps, OldFormatNumExponents, this->m_pDataRaw, this->NumSamplesRaw_ );
    ReorderOldFormat( NumOps, OldFormatNumExponents, this->m_pDataBinned, this->NumSamplesBinned_ );
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
int Model<T>::WriteAttributes( ::H5::Group &g )
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
    TypeVE *VE = this->getSummaryData();
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
          // All elements must be statistically different from zero
          bool bOK{  ( VE[o+i].*pLowZ > 0 && VE[o+i].*pHighZ > 0 )
                  || ( VE[o+i].*pLowZ < 0 && VE[o+i].*pHighZ < 0 ) };
          // Monotonic parameters (energies) must be < limit
          if( bOK && p.bMonotonic && VE[o+i].Central > MonotonicUpperLimit )
            bOK = false;
          // All elements must be different from each other
          for( std::size_t j = 0; bOK && j < i; ++j )
            bOK = VE[o+i].*pLowD > VE[o+j].*pHighD || VE[o+i].*pHighD < VE[o+j].*pLowD;
          // Save results
          VE[o+i].Check = bOK ? 1 : 0;
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
    s << "# " << sChiSqPerDof << " : " << (*this)[Model<T>::idxCentral][iCol]
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
      s << " = " << (*this)[Model<T>::idxCentral][iCol];
  }
  s << NewLine;
  // Chi^2 pValue
  iCol = this->GetColumnIndexNoThrow( sPValue );
  if( iCol >= 0 )
    s << "# " << sPValue << " : Chi^2 = " << (*this)[Model<T>::idxCentral][iCol] << NewLine;
}

template <typename T>
void Model<T>::SummaryColumnNames( std::ostream &os ) const
{
  for( std::size_t i = 0; i < FitTimes.size(); ++i )
  {
    if( i == 0 )
      os << "ti tf";
    else
      os << Space << "ti" << i << Space << "tf" << i;
  }
  os << " tfLabel NumDataPoints SampleSize dof CovarSampleSize ";
  Base::SummaryColumnNames( os );
}

template <typename T>
void Model<T>::SummaryContents( std::ostream &os ) const
{
  std::size_t NumDataPoints{ 0 };
  std::ostringstream tfl;
  for( std::size_t i = 0; i < FitTimes.size(); ++i )
  {
    NumDataPoints += FitTimes[i].size();
    if( i )
    {
      os << Space;
      tfl << Underscore << FitTimes[i][0] << Underscore;
    }
    os << FitTimes[i][0] << Space << FitTimes[i].back();
    tfl << FitTimes[i].back();
  }
  os << Space << tfl.str();
  os << Space << NumDataPoints;
  os << Space << this->SampleSize;
  os << Space << dof;
  os << Space << CovarSampleSize;
  os << Space;
  Base::SummaryContents( os );
}

template <typename T>
void Model<T>::WriteSummaryTD( const std::string &sOutFileName, bool bVerboseSummary )
{
  using Scalar = typename is_complex<T>::Scalar;
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  std::ofstream os( sOutFileName );
  SummaryHeader<T>( os, sOutFileName );
  SummaryComments( os, bVerboseSummary );
  // Write column names
  static constexpr int idxData{ 0 };
  //static constexpr int idxTheory{ 1 };
  static const char * pFieldNames[] = { "data", "theory" };
  os << "field seq model t ";
  ValWithEr<T>::Header( pFieldNames[0], os );
  os << Space;
  ValWithEr<T>::Header( pFieldNames[1], os );
  os << NewLine;
  // Grab the model data and theory with bootstrapped errors
  const int Extent{ GetExtent() };
  if( !Extent )
    return;
  std::array<VectorView<T>, 2> vv;
  std::vector<T> Buffer;
  std::array<std::vector<ValWithEr<T>>,2> Value;
  for( std::size_t i = 0; i < Value.size(); ++i )
  {
    BootRep<T> &ThD{ i == idxData ? FitInput : ModelPrediction };
    Value[i].resize( Extent );
    for( int j = 0; j < Extent; ++j )
    {
      ThD.MapColumn( vv[0], j );
      Value[i][j].Get( ThD.Central[j], vv[0], Buffer );
    }
  }
  // Write theory and data values
  int idx = 0;
  for( std::size_t m = 0; m < FitTimes.size(); ++m )
    for( const int t : FitTimes[m] )
    {
      os << "corr " << idx << Space << m << Space << t;
      for( std::size_t i = 0; i < Value.size(); ++i )
        os << Space << Value[i][idx];
      os << NewLine;
      ++idx;
    }
  // Write effective masses
  Buffer.resize( FitInput.size() );
  idx = 0;
  ValWithEr<T> v;
  for( std::size_t m = 0; m < FitTimes.size(); ++m )
  {
    ++idx;
    for( std::size_t tidx = 1; tidx < FitTimes[m].size(); ++tidx )
    {
      const double t = 0.5 * ( FitTimes[m][tidx] + FitTimes[m][tidx - 1] );
      const int DeltaT{ FitTimes[m][tidx] - FitTimes[m][tidx - 1] };
      const Scalar DeltaTInv{ static_cast<Scalar>( DeltaT == 1 ? 1. : ( 1. / DeltaT ) ) };
      os << "log " << idx << Space << m << Space << t;
      for( std::size_t i = 0; i < Value.size(); ++i )
      {
        BootRep<T> &ThD{ i == idxData ? FitInput : ModelPrediction };
        ThD.MapColumn( vv[0], idx - 1 );
        ThD.MapColumn( vv[1], idx );
        int bootCount{ 0 };
        for( int bootrep = 0; bootrep < ThD.size(); ++bootrep )
        {
          Buffer[bootCount] = std::log( vv[0][bootrep] / vv[1][bootrep] ) * DeltaTInv;
          if( IsFinite( Buffer[bootCount] ) )
            ++bootCount;
        }
        v.Get( std::log( ThD.Central[idx - 1] / ThD.Central[idx] ) * DeltaTInv, Buffer, bootCount );
        os << Space << v;
      }
      os << NewLine;
      ++idx;
    }
  }
}

template <typename T>
void Model<T>::ReorderOldFormat( int NumOps, int NumExponents, std::unique_ptr<T[]> &pData, int Num )
{
  T * p{ pData.get() };
  if( p && NumExponents > 1 && NumOps > 1 )
  {
    std::vector<T> Buffer( NumOps * NumExponents );
    while( Num-- )
    {
      for( int e = 0; e < NumExponents; ++e )
        for( int o = 0; o < NumOps; ++o )
          Buffer[o * NumExponents + e] = p[e * NumOps + o];
      for( int j = 0; j < Buffer.size(); ++j )
        p[j] = Buffer[j];
      p += this->Nt();
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

// Make non-random numbers 0 ... m.size1 - 1 to simplify central replica code
// If the sample is already bootstrapped, then we won't need these numbers
template <typename T>
void DataSet<T>::MakeCentralReplicaNonRandom()
{
  for( int i = 0; i < 2; ++i )
  {
    std::vector<fint> &c{ RandomCentralBuf[i] };
    const SS ss{ i == 0 ? SS::Raw : SS::Binned };
    const int Num{ corr.empty() || ( ss == SS::Raw && corr[0].bRawBootstrap ) ? 0 : corr[0].NumSamples( ss ) };
    if( !Num )
      c.clear();
    else
    {
      c.resize( Num );
      for( fint i = 0; i < c.size(); ++i )
        c[i] = i;
    }
  }
}

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
  // Get the central replica
  vCentral.resize( Extent );
  {
    int dst{ 0 };
    for( int f = 0; f < corr.size(); ++f )
    {
      const T * pSrc{ corr[f][Fold<T>::idxCentral] };
      for( int t : FitTimes[f] )
        vCentral[dst++] = pSrc[t];
    }
  }
  // Cache the raw data
  for( int i = 0; i < 3; ++i )
  {
    SampleSource ss{ i == 0 ? SampleSource::Raw : i == 1 ? SampleSource::Binned : SampleSource::Bootstrap };
    Matrix<T> &m{ i == 0 ? mRaw : i == 1 ? mBinned : mBoot };
    const int iCount{ corr[0].NumSamples( ss ) };
    bool bGotMatchingCount{ iCount != 0 };
    for( int f = 1; bGotMatchingCount && f < corr.size(); ++f )
      bGotMatchingCount = ( corr[f].NumSamples( ss ) == iCount );
    if( !bGotMatchingCount )
      m.clear();
    else
    {
      m.resize( iCount, Extent );
      for( int idx = 0; idx < iCount; ++idx )
      {
        int dst{ 0 };
        for( int f = 0; f < corr.size(); ++f )
        {
          const T * pSrc{ corr[f].get( ss, idx ) };
          for( int t : FitTimes[f] )
            m( idx, dst++ ) = pSrc[t];
        }
      }
    }
  }
}

// Get the constants from the appropriate timeslice
template <typename T>
void DataSet<T>::GetFixed( int idx, Vector<T> &vResult, const std::vector<FixedParam> &Params ) const
{
  // Get a pointer to the raw data in each model for replica idx
  std::vector<const T *> Src( constFile.size() );
  for( int f = 0; f < constFile.size(); ++f )
    Src[f] = constFile[f][idx];
  // Now copy the data into vResult
  for( const FixedParam &p : Params )
  {
    std::size_t SrcIdx{ p.src.param.GetOffset( 0, Param::Type::All ) };
    for( std::size_t i = 0; i < p.Count; ++i )
      vResult[p.idx + i] = Src[p.src.File][SrcIdx + i];
  }
}

// Initialise random numbers I will need if I'm to get (co)variance on non-central replicas
template <typename T>
bool DataSet<T>::InitRandomNumbers( SS ss )
{
  bool bMade{ false };
  if( ss == SS::Binned || ss == SS::Raw )
  {
    const int i{ ss == SS::Raw ? 0 : 1 };
    MatrixView<fint> &r{ RandomViews[i] };
    std::vector<fint> &Buffer{ RandomBuffer[i] };
    const int Count{corr.empty() || ( ss == SS::Raw && corr[0].bRawBootstrap ) ? 0 : corr[0].NumSamples( ss )};
    if( !Count )
    {
      r.clear();
      Buffer.clear();
    }
    else
    {
      // Make Random numbers
      if(   corr[0].RandNum() // Do I have original bootstrap random numbers
         && corr[0].SampleSize == Count ) // Are the random numbers over same range (0...SampleSize-1)
      {
        // Re-use existing random numbers
        // I assume enough replicas are available because this is checked on creation and load
        Buffer.clear();
        r.Map( corr[0].RandNum(), MaxSamples, Count );
      }
      else
      {
        // Generate random numbers
        GenerateRandom( Buffer, corr[0].Name_.Seed, MaxSamples, Count );
        r.Map( Buffer.data(), MaxSamples, Count );
        bMade = true;
      }
    }
  }
  return bMade;
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// From what is already a bootstrap replica, and therefore only the central replica
template <typename T>
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
}

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
void DataSet<T>::AddConstant( const Param::Key &Key, std::size_t File, const Param &param )
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
  constMap.insert( { Key, ConstantSource( File, param ) } );
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
    OriginalBinSize = corr[i].binSize; // Because this will be overwritten
    MakeCentralReplicaNonRandom();
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
      AddConstant( it.first, i, it.second );
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
      AddConstant( NewKey.empty() ? OldKey : NewKey, i, itP->second );
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
int DataSet<T>::NumSamplesBinned() const
{
  int NSB{ corr.empty() ? 0 : corr[0].NumSamplesBinned() };
  for( std::size_t i = 1; NSB && i < corr.size(); ++i )
  {
    const int ThisNSB{ corr[i].NumSamplesBinned() };
    if( NSB > ThisNSB )
      NSB = ThisNSB;
  }
  return NSB;
}

template <typename T>
void DataSet<T>::Rebin( const std::vector<int> &NewSize )
{
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
      corr[i].Bin( RebinSize.back(), SS::Raw );
    else
      corr[i].Bin( SS::Raw );
    if( corr[i].NumSamplesRaw() != corr[0].NumSamplesRaw() )
    {
      std::ostringstream os;
      os << "Rebinned corr " << i << " has " << corr[i].NumSamplesRaw()
         << " samples, others have " << corr[i].NumSamplesRaw();
      throw std::runtime_error( os.str().c_str() );
    }
  }
  MakeCentralReplicaNonRandom();
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
