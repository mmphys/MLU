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
#include <H5CompType.h>

// This should be the only use of the AutoConf header
// Any useful values should be exposed as functions

#include <MLUconfig.h>

extern "C" const char * MLUVersionInfoHuman()
{
  static const char PackageString[] = MLU_PACKAGE_STRING ", " MLU_GIT_SUMMARY;
  return PackageString;
}

BEGIN_COMMON_NAMESPACE

// I will use the GNU Scientific Library, but I will have GSL errors throw exceptions
class GSLLibraryGlobal
{
private:
  gsl_error_handler_t * pOldErrorHandler;
  static void GSLErrorHandler(const char * reason, const char * file, int line, int gsl_errno );

public:
  GSLLibraryGlobal();
  ~GSLLibraryGlobal();
};

GSLLibraryGlobal gslLibraryGlobal;

GSLLibraryGlobal::GSLLibraryGlobal()
{
  pOldErrorHandler = gsl_set_error_handler( GSLErrorHandler );
}

GSLLibraryGlobal::~GSLLibraryGlobal()
{
  gsl_set_error_handler( pOldErrorHandler );
}

void GSLLibraryGlobal::GSLErrorHandler(const char * reason, const char * file, int line, int gsl_errno )
{
  std::stringstream ss;
  ss << "gsl: " << file << ":" << line << ": errno " << gsl_errno << ": " << reason;
  throw std::runtime_error( ss.str() );
}

// Text required for summaries of correlators
namespace CorrSumm {
  const char sep[] = " ";
  const char Comment[] = "# ";
  const char * FieldNames[NumFields] = { "corr", "exp", "cosh" };
};

extern const std::string Empty{ "" };
extern const std::string Space{ " " };
extern const std::string WhiteSpace{ " \t\n\r\f\v" };
extern const std::string Underscore{ "_" };
extern const std::string Period{ "." };
extern const std::string Colon{ ":" };
extern const std::string NewLine{ "\n" };
extern const std::string Comma{ "," };
extern const std::string CommaSpace{ ", " };
extern const std::string Hash{ "#" };
const std::string sBootstrap{ "bootstrap" };
const std::string sFold{ "fold" };
const std::string sModel{ "model" };
const std::string sParams{ "params" };
const std::string sCormat{ "cormat" };
const std::string sNtUnfolded{ "NtUnfolded" };
const std::string st0Negated{ "t0Negated" };
const std::string sConjugated{ "Conjugated" };
const std::string sTI{ "TI" };
const std::string sTF{ "TF" };
const std::string sDoF{ "DoF" };
const std::string sNumExponents{ "NumExponents" };
const std::string sNumFiles{ "NumFiles" };
const std::string sFactorised{ "Factorised" };
const std::string sCovarFrozen{ "CovarFrozen" };
const std::string sOperators{ "Operators" };
const std::string sAuxNames{ "AuxNames" };
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
const std::string sNE{ " != " };
const std::vector<std::string> sCorrSummaryNames{ "corr", "bias", "exp", "cosh" };
const double NaN{ std::nan( "" ) };

namespace Gamma
{
  const std::array<const std::string, nGamma> name{ // Long name, per Grid
    "MinusGamma5",
    "Gamma5",
    "MinusGammaT",
    "GammaT",
    "MinusGammaTGamma5",
    "GammaTGamma5",
    "MinusGammaX",
    "GammaX",
    "MinusGammaXGamma5",
    "GammaXGamma5",
    "MinusGammaY",
    "GammaY",
    "MinusGammaYGamma5",
    "GammaYGamma5",
    "MinusGammaZ",
    "GammaZ",
    "MinusGammaZGamma5",
    "GammaZGamma5",
    "MinusIdentity",
    "Identity",
    "MinusSigmaXT",
    "SigmaXT",
    "MinusSigmaXY",
    "SigmaXY",
    "MinusSigmaXZ",
    "SigmaXZ",
    "MinusSigmaYT",
    "SigmaYT",
    "MinusSigmaYZ",
    "SigmaYZ",
    "MinusSigmaZT",
    "SigmaZT",
  };
  const std::array<const std::string, nGamma> nameShort{ // My abbreviations
    "g-5",
    "g5",
    "g-T",
    "gT",
    "g-T5",
    "gT5",
    "g-X",
    "gX",
    "g-X5",
    "gX5",
    "g-Y",
    "gY",
    "g-Y5",
    "gY5",
    "g-Z",
    "gZ",
    "g-Z5",
    "gZ5",
    "g-I",
    "gI",
    "g-sXT",
    "gsXT",
    "g-sXY",
    "gsXY",
    "g-sXZ",
    "gsXZ",
    "g-sYT",
    "gsYT",
    "g-sYZ",
    "gsYZ",
    "g-sZT",
    "gsZT",
  };
  const std::string sGammaSpatial{ "gXYZ" };

  std::string NameShort( Algebra alg, const char * pszPrefix, const char * pszOpSuffix )
  {
    std::string sName;
    const bool bGotSuffix{ pszOpSuffix && *pszOpSuffix };
    if( alg != Algebra::Unknown || bGotSuffix )
    {
      if( pszPrefix && *pszPrefix )
        sName.append( pszPrefix );
      if( alg == Algebra::Spatial )
        sName.append( sGammaSpatial );
      else if( alg != Algebra::Unknown )
      {
        int idx = static_cast<int>( alg );
        sName.append( ( idx >= 0 && idx < nGamma ) ? nameShort[idx] : "?" );
      }
      if( bGotSuffix )
        sName.append( pszOpSuffix );
    }
    return sName;
  }

  std::ostream& operator<<(std::ostream& os, const Algebra &a)
  {
    return os << NameShort( a );
  }

  std::istream& operator>>(std::istream& is, Algebra &a)
  {
    std::string s;
    if( is >> s )
    {
      int i;
      for( i = 0; i < nGamma && !EqualIgnoreCase( s, nameShort[i] ); i++ )
        ;
      if( i < nGamma )
        a = static_cast<Algebra>( i );
      else if( EqualIgnoreCase( s, sGammaSpatial ) )
        a = Algebra::Spatial;
      else
        is.setstate( std::ios_base::failbit );
    }
    return is;
  }
};

// Default delimeters for the next couple of functions
extern const char szDefaultDelimeters[] = " \t,";

// Remove anything past the last delimeter from string, returning the removed part in suffix
// Return success / fail
bool ExtractSuffix( std::string &String, std::string &Suffix, const char * pszDelimeters )
{
  if( !pszDelimeters || !*pszDelimeters )
    pszDelimeters = szDefaultDelimeters;
  std::size_t NumDelims{ 0 };
  while( pszDelimeters[NumDelims] )
    NumDelims++;
  std::size_t Len{ String.length() };
  bool bFoundDelim{ false };
  while( !bFoundDelim && Len )
  {
    const char c{ String[--Len] };
    for( std::size_t i = 0; !bFoundDelim && i < NumDelims; i++ )
    {
      bFoundDelim = ( pszDelimeters[i] == c );
      if( bFoundDelim )
      {
        Suffix = String.substr( Len + 1 );
        Trim( Suffix );
        // Skip past multiple delimeters if they are all whitespace
        if( c == ' ' || c == '\t' || c == '\r' || c == '\n' )
        {
          while( Len && ( String[Len - 1] == ' ' || String[Len - 1] == '\t'
                         || String[Len - 1] == '\r' || String[Len - 1] == '\n' ) )
            --Len;
        }
        String.resize( Len );
      }
    }
  }
  return bFoundDelim;
}

// Remove the directory from the start of FileName (leave the trailing '/' in place)
std::string ExtractDirPrefix( std::string &FileName )
{
  std::string Dir;
  std::size_t pos{ FileName.find_last_of( '/' ) };
  if( pos != std::string::npos )
  {
    Dir = FileName.substr( 0, pos + 1 );
    FileName.erase( 0, pos + 1 );
  }
  return Dir;
}

// Split String into an array using specified delimeters

std::vector<std::string> Split( const std::string &String, const char * pszDelimeters )
{
  std::vector<std::string> a;
  if( !pszDelimeters || !*pszDelimeters )
    pszDelimeters = szDefaultDelimeters;
  std::size_t NumDelims{ 0 };
  while( pszDelimeters && pszDelimeters[NumDelims] )
    NumDelims++;
  const std::size_t Len{ String.length() };
  std::size_t Start{ 0 };
  while( Start < Len )
  {
    // Look for the next delimeter
    std::size_t Pos{ Start };
    bool bFoundDelim{ false };
    char c = 0;
    while( Pos < Len && !bFoundDelim )
    {
      c = String[Pos];
      for( std::size_t i = 0; !bFoundDelim && i < NumDelims; i++ )
        bFoundDelim = ( pszDelimeters[i] == c );
      if( !bFoundDelim )
        Pos++;
    }
    // Append this substring to list of items to return
    a.push_back( String.substr( Start, Pos - Start ) );
    // Skip past this delimeter
    Start = Pos + 1;
    // Skip past multiple delimeters if they are all whitespace
    if( c == ' ' || c == '\t' || c == '\r' || c == '\n' )
    {
      while( Start < Len && ( String[Start] == ' ' || String[Start] == '\t'
                             || String[Start] == '\r' || String[Start] == '\n' ) )
        ++Start;
    }
  }
  return a;
}

// Extract suffix, then split strings. Default delimeters '.' and '_' respectively
bool ExtractSuffixSplit( std::string &String, std::vector<std::string> &Suffii,
                        const char * pszStringDelim, const char * pszSuffixDelim )
{
  if( !pszStringDelim || !*pszStringDelim )
    pszStringDelim = Period.c_str();
  if( !pszSuffixDelim || !*pszSuffixDelim )
    pszSuffixDelim = Underscore.c_str();
  std::string Suffix;
  const bool bExtracted{ ExtractSuffix( String, Suffix, pszStringDelim ) };
  if( bExtracted )
    Suffii = Split( Suffix, pszSuffixDelim );
  return bExtracted;
}

// Zipper merge v1 and v2 if same size (otherwise just append)
std::vector<std::string> ZipperMerge( const std::vector<std::string> &v1, const std::vector<std::string> &v2 )
{
  std::vector<std::string> myFileList;
  // Merge the file lists
  if( v1.size() == v2.size() )
  {
    // Same length - merge the lists like a zipper
    myFileList.reserve( v1.size() + v2.size() );
    auto p2 = v2.begin();
    for( auto p1 = v1.begin(); p1 != v1.end(); ++p1, ++p2 )
    {
      myFileList.emplace_back( *p1 );
      myFileList.emplace_back( *p2 );
    }
  }
  else
  {
    myFileList = v1;
    myFileList.insert( myFileList.end(), v2.begin(), v2.end() );
  }
  return myFileList;
}

// Dump the environment to stdout, prefixed by optional message
void DumpEnv(int argc, const char * const *argv, const char * pStr )
{
  static const char sIndent2[]{ "    " };
  static const char * sIndent1{ sIndent2 + 2 };
  static const char sQuote[]{ "\"" };
  if( pStr )
    std::cout << pStr << std::endl;
  std::cout << sIndent1 << argc << " arguments:" << std::endl;
  for( int i = 0; i < argc; i++ )
  {
    std::cout << sIndent2 << "argv[" << i << "] = ";
    if( argv[i] )
      std::cout << sQuote << argv[i] << sQuote;
    else
      std::cout << "nullptr";
    std::cout << std::endl;
  }
}

// Does the specified file exist?
bool FileExists( const std::string& Filename )
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
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

// Sort the list of values, then extract the lower and upper 68th percentile error bars
template<typename T>
void ValWithEr<T>::Get( T Central_, std::vector<T> &Data, std::size_t Count )
{
  assert( Data.size() >= Count && "ValWithErr<T>::Get() data too small" );
  Central = Central_;
  if( Count == 0 )
  {
    if( Data.size() == 0 )
    {
      // I wasn't trying to take an estimate over many samples
      Min = Central;
      Max = Central;
      Low = Central;
      High  = Central;
      Check = 1;
    }
    else
    {
      // Tried to take an estimate over many samples, but all values NaN
      Min = NaN;
      Max = NaN;
      Low = NaN;
      High  = NaN;
      Check = 0;
    }
    return;
  }
  const typename std::vector<T>::iterator itStart{ Data.begin() };
  const typename std::vector<T>::iterator itEnd  { itStart + Count };
  std::sort( itStart, itEnd );
  std::size_t Index = 0.16 * Count + 0.5;
  if( Index >= Count )
    Index = Count - 1;
  Min  = Data[0];
  Max  = Data[Count - 1];
  Low  = Data[Index];
  High = Data[Count - 1 - Index];
  Check = static_cast<T>( Count ) / Data.size();
}

template class ValWithEr<float>;
template class ValWithEr<double>;
template class ValWithEr<long double>;

const Momentum p0(0,0,0);
const std::string Momentum::DefaultPrefix{ "p" };
const std::string Momentum::SinkPrefix{ DefaultPrefix + "s" };
const std::string Momentum::SquaredSuffix{ "2" };

std::string Momentum::p2_string( const std::string &separator, const std::string Prefix ) const
{
  std::string s{ separator };
  s.append( Prefix );
  s.append( SquaredSuffix );
  s.append( separator );
  s.append( std::to_string( p2() ) );
  return s;
}

std::string Momentum::to_string( const std::string &separator, bool bNegative ) const
{
  return std::to_string( bNegative ? -x : x ) + separator
  + std::to_string( bNegative ? -y : y ) + separator
  + std::to_string( bNegative ? -z : z );
}

std::string Momentum::to_string4d( const std::string &separator, bool bNegative ) const
{
  return to_string( separator, bNegative ) + separator + "0";
}

std::regex Momentum::MakeRegex( const std::string &MomName )
{
  std::string sPattern( Common::Underscore );
  sPattern.append( MomName );
  sPattern.append( "_(-?[0-9]+)_(-?[0-9]+)_(-?[0-9]+)" );
  return std::regex( sPattern );
}

// Strip out momentum info from string if present
bool Momentum::Extract( std::string &Prefix, const std::string &MomName, bool IgnoreSubsequentZeroNeg )
{
  bool bGotMomentum = false;
  std::smatch match;
  const std::regex Pattern{ MakeRegex( MomName ) };
  while( std::regex_search( Prefix, match, Pattern ) )
  {
    int px{ std::stoi( match[1] ) };
    int py{ std::stoi( match[2] ) };
    int pz{ std::stoi( match[3] ) };
    if( !bGotMomentum )
    {
      bGotMomentum = true;
      x = px;
      y = py;
      z = pz;
    }
    else if( x != px || y != py || z != pz )
    {
      if( IgnoreSubsequentZeroNeg && !(*this) ) // i.e. previous momentum was zero
      {
        x = px;
        y = py;
        z = pz;
      }
      else if(IgnoreSubsequentZeroNeg && ((px==0 && py==0 && pz==0) || (px==-x && py==-y && pz==-z)))
        ;
      else
      {
        static const std::string Sep{ "," };
        static const std::string m1{ to_string( Sep ) };
        x = px;
        y = py;
        z = pz;
        throw std::runtime_error( "Multiple momenta: " + m1 + " != " + to_string( Sep ) );
      }
    }
    const std::string sSuffix{ match.suffix() };
    Prefix = match.prefix();
    Prefix.append( sSuffix );
  }
  return bGotMomentum;
}

void Momentum::Replace( std::string &s, const std::string &MomName, bool bNegative ) const
{
  std::string Replace( Common::Underscore );
  Replace.append( MomName );
  Replace.append( Common::Underscore );
  Replace.append( to_string( Common::Underscore, bNegative ) );
  std::smatch match;
  const std::regex Pattern{ MakeRegex( MomName ) };
  bool bFound{ false };
  std::string Search{ std::move( s ) };
  while( std::regex_search( Search, match, Pattern ) )
  {
    s.append( match.prefix() );
    s.append( Replace );
    Search = match.suffix();
    bFound = true;
  }
  if( !bFound )
  {
    std::ostringstream ss;
    ss << "momentum " << MomName << " not found in " << s;
    throw std::runtime_error( ss.str() );
  }
  s.append( Search );
}

// Take a squared momentum and make an approximate 3d momentum from it
Momentum FileNameMomentum::FromSquared( const int p2Original )
{
  int p2{ p2Original };
  Momentum p;
  const bool bNegative{ p2 < 0 };
  if( bNegative )
    p2 = - p2;
  for( int i = 0; i < 3 && p2; ++i )
  {
    int root = 1;
    while( ( root + 1 ) * ( root + 1 ) <= p2 )
      root++;
    p2 -= root * root;
    switch( i )
    {
      case 0:
        p.x = root;
        break;
      case 1:
        p.y = root;
        break;
      case 2:
        p.z = root;
        break;
    }
  }
  if( p2 )
    throw std::runtime_error( std::to_string( p2Original ) + " can't be expressed as a three-momentum" );
  if( bNegative )
    p.x = -p.x;
  return p;
}

std::ostream& operator<<( std::ostream& os, const Momentum &p )
{
  return os << p.to_string( Underscore );
}

std::istream& operator>>( std::istream& is, Momentum &p )
{
  char c = 0;
  if(!(is>>p.x && is.get(c) && c==Underscore[0] && is>>p.y && is.get(c) && c==Underscore[0] && is>>p.z))
    is.setstate( std::ios_base::failbit );
  return is;
}

// These are the attributes I like to use in my filenames
void FileNameAtt::Parse( const std::string &Filename_, std::vector<std::string> * pOpNames,
                         const std::vector<std::string> * pIgnoreMomenta )
{
  clear();
  Filename = Filename_;
  std::size_t pos = Filename.find_last_of( '/' );
  if( pos == std::string::npos )
    pos = 0;
  else
    Dir = Filename.substr( 0, ++pos );
  Base = Filename.substr( pos );
  NameNoExt = Base;
  int i = 0;
  while( i < 3 && ( pos = Base.find_last_of( '.' ) ) != std::string::npos )
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
  // Now see whether we can extract operator names
  if( pOpNames )
    ParseOpNames( *pOpNames, 2, pIgnoreMomenta );
  else
    ParseShort( pIgnoreMomenta );
}

void FileNameAtt::ParseShort( const std::vector<std::string> * pIgnoreMomenta )
{
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
  BaseShort = Base;
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
  {
    // Remove hit info from the filename, but no need to expose it
    bool bHasHit = false;
    int Hit = 0;
    ExtractInteger( BaseShort, bHasHit, Hit, "[hH][iI][tT]" );
  }
  Gamma.clear();
  Gamma = ExtractGamma( BaseShort );
}

// Parse operator names from the end of Base, building up a list of all the operator names
// NB: Because I parse from the end, op[0] is the last operator, op[1] second last, etc
std::vector<std::string> FileNameAtt::ParseOpNames( int NumOps, const std::vector<std::string> * pIgnoreMomenta )
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
    else if( i != NumOps - 1 && Base[pos] == '.' )
    {
      break;
    }
    else
    {
      o.push_back( Base.substr( pos + 1, LastPos - pos ) );
      LastPos = pos;
    }
  }
  if( i != NumOps )
    throw std::runtime_error( "Can't extract " + std::to_string( NumOps ) + " operators from " + Base );
  if( i )
  {
    Base.resize( LastPos ); // Shorten the base
  }
  ParseShort( pIgnoreMomenta );
  return o;
}

void FileNameAtt::ParseOpNames( std::vector<std::string> &OpNames, int NumOps,
                                const std::vector<std::string> * pIgnoreMomenta )
{
  std::vector<std::string> Ops{ ParseOpNames( NumOps, pIgnoreMomenta ) };
  op.clear();
  op.reserve( NumOps );
  for( std::string &s : Ops )
  {
    int iOp = 0;
    while( iOp < OpNames.size() && !EqualIgnoreCase( s, OpNames[iOp] ) )
      iOp++;
    if( iOp == OpNames.size() )
      OpNames.emplace_back( std::move( s ) );
    op.push_back( iOp );
  }
}

void FileNameAtt::ParseExtra( unsigned int MaxElements, const std::vector<std::string> * pIgnoreMomenta )
{
  std::size_t pos;
  while( MaxElements-- && ( pos = Base.find_last_of( '.' ) ) != std::string::npos )
  {
    Extra.push_back( Base.substr( pos + 1 ) );
    Base.resize( pos );
  }
  ParseShort( pIgnoreMomenta );
}

// Append the extra info to the string
void FileNameAtt::AppendExtra( std::string &s, int Last, int First ) const
{
  int Num{ static_cast<int>( Extra.size() ) };
  if( First >= 0 && First < Num )
    Num = First;
  while( Num-- > Last )
  {
    s.append( 1, '.' );
    s.append( Extra[Num] );
  }
}

std::string FileNameAtt::GetBaseExtra( int Last, int First ) const
{
  std::string s{ Base };
  AppendExtra( s, Last, First );
  return s;
}

/*std::string FileNameAtt::GetBaseShortExtra( int Last, int First ) const
{
  std::string s{ BaseShort };
  AppendExtra( s, Last, First );
  return s;
}*/

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
    if( match[1].length() )
      s.append( match[1] );
    else if( match[2].length() )
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
  while( std::regex_search( Search, match, GammaPattern ) )
  {
    Prefix.append( match.prefix() );
    Prefix.append( GammaReplace );
    Prefix.append( match[1] );
    Search = match.suffix();
  }
  Prefix.append( Search );
}

// Make the same HDF5 complex type Grid uses
template<typename T> static ::H5::CompType MakeComplex()
{
  ::H5::CompType myComplex( sizeof( std::complex<T> ) );
  myComplex.insertMember("re", 0 * sizeof(T), H5::Equiv<T>::Type);
  myComplex.insertMember("im", 1 * sizeof(T), H5::Equiv<T>::Type);
  return myComplex;
}

// Make an HDF5 type representing a value with an error
template<typename T> static ::H5::CompType MakeValWithErOldV1()
{
  ::H5::CompType myType( sizeof( ValWithErOldV1<T> ) );
  myType.insertMember("Central", offsetof(ValWithErOldV1<T>, Central), H5::Equiv<T>::Type);
  myType.insertMember("Low",     offsetof(ValWithErOldV1<T>, Low    ), H5::Equiv<T>::Type);
  myType.insertMember("High",    offsetof(ValWithErOldV1<T>, High   ), H5::Equiv<T>::Type);
  myType.insertMember("Check",   offsetof(ValWithErOldV1<T>, Check  ), H5::Equiv<T>::Type);
  return myType;
}

// Make an HDF5 type representing a value with an error
template<typename T> static ::H5::CompType MakeValWithEr()
{
  ::H5::CompType myType( sizeof( ValWithEr<T> ) );
  myType.insertMember("Min",     offsetof(ValWithEr<T>, Min    ), H5::Equiv<T>::Type);
  myType.insertMember("Low",     offsetof(ValWithEr<T>, Low    ), H5::Equiv<T>::Type);
  myType.insertMember("Central", offsetof(ValWithEr<T>, Central), H5::Equiv<T>::Type);
  myType.insertMember("High",    offsetof(ValWithEr<T>, High   ), H5::Equiv<T>::Type);
  myType.insertMember("Max",     offsetof(ValWithEr<T>, Max    ), H5::Equiv<T>::Type);
  myType.insertMember("Check",   offsetof(ValWithEr<T>, Check  ), H5::Equiv<T>::Type);
  return myType;
}

// Make an HDF5 type representing a config number and count
static ::H5::CompType MakeConfigCount()
{
  ::H5::CompType myType( sizeof( ConfigCount ) );
  myType.insertMember("Config", offsetof(ConfigCount, Config), ::H5::PredType::NATIVE_UINT32);
  myType.insertMember("Count",  offsetof(ConfigCount, Count ), ::H5::PredType::NATIVE_UINT32);
  return myType;
}

const ::H5::PredType& H5::Equiv<float>      ::Type{ ::H5::PredType::NATIVE_FLOAT };
const ::H5::PredType& H5::Equiv<double>     ::Type{ ::H5::PredType::NATIVE_DOUBLE };
const ::H5::PredType& H5::Equiv<long double>::Type{ ::H5::PredType::NATIVE_LDOUBLE };
const ::H5::StrType   H5::Equiv<std::string>::Type{ ::H5::PredType::C_S1, H5T_VARIABLE };
const ::H5::StrType&  H5::Equiv<char *>     ::Type{   H5::Equiv<std::string>::Type };
const ::H5::PredType& H5::Equiv<std::uint_fast32_t>::Type{ sizeof( std::uint_fast32_t ) == 4
                                                    ? ::H5::PredType::STD_U32LE : ::H5::PredType::STD_U64LE };
const ::H5::CompType  H5::Equiv<std::complex<float>>      ::Type{ MakeComplex<float>() };
const ::H5::CompType  H5::Equiv<std::complex<double>>     ::Type{ MakeComplex<double>() };
const ::H5::CompType  H5::Equiv<std::complex<long double>>::Type{ MakeComplex<long double>() };
const ::H5::CompType  H5::Equiv<ValWithEr<float>>      ::Type{ MakeValWithEr<float>() };
const ::H5::CompType  H5::Equiv<ValWithEr<double>>     ::Type{ MakeValWithEr<double>() };
const ::H5::CompType  H5::Equiv<ValWithEr<long double>>::Type{ MakeValWithEr<long double>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<float>>      ::Type{ MakeValWithErOldV1<float>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<double>>     ::Type{ MakeValWithErOldV1<double>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<long double>>::Type{ MakeValWithErOldV1<long double>() };
const ::H5::CompType  H5::Equiv<ConfigCount> ::Type{ MakeConfigCount() };

// Open the specified HDF5File and group
void H5::OpenFileGroup(::H5::H5File &f, ::H5::Group &g, const std::string &FileName, const char *PrintPrefix,
                       std::string * pGroupName, unsigned int flags)
{
  const bool bFindGroupName{ !pGroupName || pGroupName->empty() };
  std::string localGroupName;
  if( !pGroupName )
    pGroupName = &localGroupName;
  f.openFile( FileName, flags );
  g = f.openGroup( bFindGroupName ? std::string("/") : *pGroupName );
  if( bFindGroupName ) {
    *pGroupName = GetFirstGroupName( g );
    g = g.openGroup( *pGroupName );
  }
  if( PrintPrefix )
    std::cout << PrintPrefix << FileName << " (" << *pGroupName << ")\n";
}

// Get first groupname from specified group
std::string H5::GetFirstGroupName( ::H5::Group & g )
{
  hsize_t n = g.getNumObjs();
  for( hsize_t i = 0; i < n; ++i ) {
    H5G_obj_t t = g.getObjTypeByIdx( i );
    if( t == H5G_GROUP )
      return g.getObjnameByIdx( i );
  }
  return std::string();
}

// Read the gamma algebra attribute string and make sure it's valid
Gamma::Algebra H5::ReadGammaAttribute( ::H5::Group &g, const char * pAttName )
{
  std::string sGamma;
  ::H5::Attribute a = g.openAttribute( pAttName );
  ::H5::StrType s = a.getStrType();
  a.read( s, sGamma );
  a.close();
  for( int idxGamma = 0; idxGamma < Gamma::nGamma; idxGamma++ )
    if( EqualIgnoreCase( sGamma, Gamma::name[idxGamma] ) )
      return static_cast<Gamma::Algebra>( idxGamma );
  throw ::H5::Exception( "Common::ReadGammaAttribute", "Invalid gamma algebra string" );
}

template<> void H5::ReadStringsHelper<::H5::Attribute>( const ::H5::Attribute &a, const ::H5::StrType &aType, char * * MDString )
{
  a.read( aType, ( void * ) MDString );
}
template<> void H5::ReadStringsHelper<::H5::DataSet>( const ::H5::DataSet &ds, const ::H5::StrType &aType, char * * MDString )
{
  ds.read( ( void * ) MDString, aType );
}

// Make a multi-dimensional string attribute
void H5::WriteAttribute( ::H5::Group &g, const std::string &AttName, const std::vector<std::string> &vs)
{
  std::unique_ptr<char *[]> RawArray( new char * [vs.size()] );
  const char * * const MDString{ const_cast<const char * *>( RawArray.get() ) };
  for( std::size_t i = 0; i < vs.size(); i++ )
    MDString[i] = vs[i].c_str();
  const hsize_t NDimension{ vs.size() };
  ::H5::DataSpace dsN( 1, &NDimension );
  ::H5::Attribute a{ g.createAttribute( AttName, Equiv<std::string>::Type, dsN ) };
  a.write( Equiv<std::string>::Type, MDString );
  a.close();
}

// Make a multi-dimensional string attribute
void H5::WriteStringData( ::H5::Group &g, const std::string &DSName, const std::vector<std::string> &vs)
{
  std::unique_ptr<char *[]> RawArray( new char * [vs.size()] );
  const char * * const MDString{ const_cast<const char * *>( RawArray.get() ) };
  for( std::size_t i = 0; i < vs.size(); i++ )
    MDString[i] = vs[i].c_str();
  const hsize_t NDimension{ vs.size() };
  ::H5::DataSpace dsN( 1, &NDimension );
  ::H5::DataSet ds{ g.createDataSet( DSName, Equiv<std::string>::Type, dsN ) };
  ds.write( MDString, Equiv<std::string>::Type );
  ds.close();
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
template void ReadArray<long double>( template_ReadArrayArgs( long double ) );
template void ReadArray<std::complex<float>>( template_ReadArrayArgs( std::complex<float> ) );
template void ReadArray<std::complex<double>>( template_ReadArrayArgs( std::complex<double> ) );
template void ReadArray<std::complex<long double>>( template_ReadArrayArgs(std::complex<long double>));

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

void CommandLine::SkipPastSep( const char * & p )
{
  while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
    p++;
  if( *p == '=' ) {
    ++p;
    while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
      p++;
  }
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
      // Swallow any trailing whitespace
      while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
        p++;
      if( defs[SwitchNo].Type == SwitchType::Flag ) {
        // This is just a switch - it should not have a value
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
      SkipPastSep( p );
      if( *p == 0 )
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
