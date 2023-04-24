/*****************************************************************************
 
 String manipulation

 Source file: String.hpp

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
 ****************************************************************************/
/*  END LEGAL */

#ifndef MLU_String_hpp
#define MLU_String_hpp namespace Common {
#define MLU_String_hpp_end };

#include <algorithm>
#include <array>
#include <istream>
#include <ios>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

MLU_String_hpp

extern const std::string Empty;
extern const std::string Space;
extern const std::string WhiteSpace;
extern const std::string Underscore;
extern const std::string Period;
extern const std::string Colon;
extern const std::string NewLine;
extern const std::string Comma;
extern const std::string CommaSpace;
extern const std::string Hash;
extern const std::string EqualSign;

// Default delimeters for ... functions
extern const char szDefaultDelimeters[];

// Compare characters, by default ignoring case
template <class _CharT>
inline int Compare( _CharT a, _CharT b, bool bIgnoreCase = true )
{
  if( bIgnoreCase )
  {
    a = std::toupper( a );
    b = std::toupper( b );
  }
  if( a == b )
    return 0;
  if( a < b )
    return -1;
  return 1;
}

// Check that a stream is empty ... or only contains white space to end of stream
template <class _CharT, class _Traits>
inline bool StreamEmpty( std::basic_istream<_CharT, _Traits> & s )
{
  return s.eof() || ( s >> std::ws && s.eof() );
}

// Return true if the next, non-whitespace character in the stream is the character specified (and if so, consume it)
template <class _CharT, class _Traits>
inline bool NextCharIs( std::basic_istream<_CharT, _Traits> & s, _CharT c, bool bOptional = false, bool bIgnoreCase = true )
{
  if( StreamEmpty( s ) )
    return bOptional;
  _CharT d = s.peek();
  const bool bIsEqual{ !Compare( c, d, bIgnoreCase ) };
  if( bIsEqual )
    s.get();
  return bIsEqual;
}

enum class NegateStar { None, Negate, Star, NegateStar };

std::ostream & operator<<( std::ostream &os, const NegateStar &ns );
std::istream & operator>>( std::istream& is, NegateStar &ns );

// Compare two strings, case insensitive
inline int CompareIgnoreCase(const std::string &s1, const std::string &s2)
{
  const std::size_t Len1{ s1.size() };
  const std::size_t Len2{ s2.size() };
  const std::size_t Len{ std::min( Len1, Len2 ) };
  for( std::size_t i = 0; i < Len; ++i )
  {
    int compare = Compare( s1[i], s2[i] );
    if( compare )
      return compare;
  }
  if( Len1 == Len2 )
    return 0;
  if( Len1 < Len2 )
    return -1;
  return 1;
}

// Compare two strings, case insensitive
inline bool EqualIgnoreCase(const std::string &s1, const std::string &s2)
{
  if( s1.size() != s2.size() )
    return false;
  return !CompareIgnoreCase( s1, s2 );
}

// Find the index of the specified string in the array, or array size if not found
template <typename C = std::vector<std::string>>
inline int IndexIgnoreCase( const C &v, const std::string &s )
{
  assert( v.size() < std::numeric_limits<int>::max() && "Terribly inefficient search!" );
  int idx = 0;
  for( ; idx < v.size() && !Common::EqualIgnoreCase( v[idx], s ); ++idx )
    ;
  return idx;
}

// Ensure that the vector has the required length and no duplicates

template <typename T> void NoDuplicates( const std::vector<T> &v, const std::string &sErrorPrefix, std::size_t MinSize )
{
  if( v.size() < MinSize )
  {
    std::stringstream ss;
    ss << sErrorPrefix << " contains " << v.size() << " entries, but " << MinSize
       << (MinSize == 1 ? " is" : " are") << " required";
    throw std::runtime_error( ss.str() );
  }
  if( v.size() > 1 )
  {
    std::vector<T> vc{ v };
    std::sort( vc.begin(), vc.end() );
    const auto dup = std::adjacent_find( vc.begin(), vc.end() );
    if( dup != vc.end() )
    {
      std::stringstream ss;
      ss << sErrorPrefix << " contains duplicates, e.g. " << *dup;
      throw std::runtime_error( ss.str() );
    }
  }
}

struct LessCaseInsensitive
{
  bool operator()( const std::string &lhs, const std::string &rhs ) const { return CompareIgnoreCase( lhs, rhs ) < 0;}
};

// A list of case insensitive, but unique names, each mapped to an int
using UniqueNames = std::map<std::string, int, Common::LessCaseInsensitive>;
using UniqueNameSet = std::set<std::string, Common::LessCaseInsensitive>;

inline std::ostream &operator<<( std::ostream &os, const UniqueNameSet &uns )
{
  bool bFirst{ true };
  for( const std::string &s : uns )
  {
    if( bFirst )
      bFirst = false;
    else
      os << Space;
    os << s;
  }
  return os;
}

inline bool operator==( const UniqueNameSet &lhs, const UniqueNameSet &rhs )
{
  if( lhs.size() != rhs.size() )
    return false;
  UniqueNameSet::const_iterator l{ lhs.cbegin() };
  UniqueNameSet::const_iterator r{ rhs.cbegin() };
  for( ; l != lhs.cend(); ++l, ++r )
    if( !EqualIgnoreCase( *l, *r ) )
      return false;
  return true;
}

inline bool operator!=( const UniqueNameSet &lhs, const UniqueNameSet &rhs )
{
  return !( lhs == rhs );
}

// Remove leading and trailing whitespace from string
// Returns true if the string contains something after the trim
inline bool Trim( std::string &s )
{
  bool bHasString( !s.empty() );
  if( bHasString )
  {
    std::size_t first{ s.find_first_not_of( WhiteSpace ) };
    if( first == std::string::npos )
    {
      s.clear();
      bHasString = false;
    }
    else
    {
      std::size_t OnePastLast{ s.find_last_not_of( WhiteSpace ) + 1 };
      if( first )
        s = s.substr( first, OnePastLast - first );
      else if( OnePastLast != s.length() )
        s.resize( OnePastLast );
    }
  }
  return bHasString;
}

// Generic conversion from a string to any type
template<typename T> inline T FromString( const std::string &String )
{
  T t;
  std::istringstream iss( String );
  if( std::is_same<T, bool>::value )
    iss >> std::boolalpha;
  if( !( iss >> t && StreamEmpty( iss ) ) )
    throw std::invalid_argument( "Argument \"" + String + "\" is not type " + typeid(T).name() );
  return t;
}

template<typename T> inline T FromString( std::string &&String )
{
  T t;
  std::istringstream iss( std::move( String ) );
  if( std::is_same<T, bool>::value )
    iss >> std::boolalpha;
  if( !( iss >> t && StreamEmpty( iss ) ) )
    throw std::invalid_argument( "Argument \"" + String + "\" is not type " + typeid(T).name() );
  return t;
}

// Converting a string to a string makes a copy
template<> inline std::string FromString<std::string>( const std::string &String )
{
  std::string s{ String };
  Trim( s );
  return s;
}

// Converting a string to a string can also move
template<> inline std::string FromString<std::string>( std::string &&String )
{
  std::string s{ std::move( String ) };
  Trim( s );
  return s;
}

// Generic conversion from a string to any type - with a default if missing
template<typename T> inline T FromString( const std::string &String, T Default )
{
  std::istringstream iss( String );
  if( StreamEmpty( iss ) )
    return Default;
  T t;
  if( !( iss >> t && StreamEmpty( iss ) ) )
    throw std::invalid_argument( "Argument \"" + String + "\" is not type " + typeid(T).name() );
  return t;
}

// Converting a string to a string makes a copy - unless the string is empty, in which case we get the default
template<> inline std::string FromString<std::string>( const std::string &String, std::string Default )
{
  std::string s{ String };
  if( Trim( s ) )
    return s;
  return Default;
}

// Generic conversion from a string to an array of any type (comma or space separated)
template<typename T> inline typename std::enable_if<!std::is_same<T, std::string>::value, std::vector<T>>::type
  ArrayFromString( const std::string &String, std::vector<NegateStar> *pNeg = nullptr )
{
  std::string s{ String };
  for( std::size_t pos = 0; ( pos = s.find( ',', pos ) ) != std::string::npos; )
    s[pos] = ' ';
  std::istringstream iss( s );
  std::vector<T> v;
  if( pNeg )
    pNeg->clear();
  for( T t; !StreamEmpty( iss ); )
  {
    NegateStar NS = NegateStar::None;
    if( pNeg )
      iss >> NS;
    if( !( iss >> t ) )
      throw std::invalid_argument( "ArrayFromString: \"" + String + "\" is not type " + typeid(T).name() );
    v.push_back( t );
    if( pNeg )
      pNeg->push_back( NS );
  }
  return v;
}

// Split string into an array of strings using separator(s). Trim leading and trailing white space
inline std::vector<std::string> ArrayFromString( const std::string &String, const std::string &Separators = Comma )
{
  std::vector<std::string> v;
  std::string s{ String };
  std::size_t pos = 0;
  do
  {
    // Skip past leading white space
    while( pos < s.size() && std::isspace( s[pos] ) )
      pos++;
    if( pos < s.size() )
    {
      // There is definitely a string here
      std::size_t end = s.find_first_of( Separators, pos );
      if( end == std::string::npos )
        end = s.size();
      const std::size_t LastSep{ end };
      // Now remove trailing white space
      while( end > pos && std::isspace( s[end - 1] ) )
        --end;
      v.emplace_back( &s[pos], end - pos );
      // Start looking for next sub-string beyond the separator
      pos = LastSep;
      if( pos < s.size() )
        ++pos;
    }
  }
  while( pos < s.size() );
  return v;
}

// Extract and return up to first separator. Trim leading and trailing white space
inline std::string ExtractToSeparator( std::string &String, const std::string &Separators = Comma )
{
  std::string s;
  if( !String.empty() )
  {
    const std::size_t first{ String.find_first_not_of( WhiteSpace ) };
    if( first == std::string::npos )
      String.clear();
    else
    {
      std::size_t pos = String.find_first_of( Separators, first );
      if( pos == std::string::npos )
        s = std::move( String );
      else
      {
        // Return the string from first non-whitespace up to the separator
        s = String.substr( first, pos - first );
        // Now trim the trailing string
        const std::size_t part2_start{ String.find_first_not_of( WhiteSpace, pos + 1 ) };
        if( part2_start == std::string::npos )
          String.clear();
        else
        {
          const std::size_t part2_stop{ String.find_last_not_of( WhiteSpace ) };
          String = String.substr( part2_start, part2_stop - part2_start + 1 );
        }
      }
      Trim( s );
    }
  }
  return s;
}

// Remove the number from the end of a string
inline bool ExtractTrailing( std::string &s, int &i )
{
  std::size_t pos{ s.find_last_not_of( WhiteSpace ) };
  if( pos != std::string::npos && std::isdigit( s[pos] ) )
  {
    while( pos && std::isdigit( s[pos - 1] ) )
      --pos;
    if( pos && s[pos - 1] == '-' )
      --pos;
    std::stringstream ss( s.substr( pos ) );
    if( ss >> i )
    {
      s.resize( pos );
      return true;
    }
  }
  return false;
}

inline std::string AppendSlash( const std::string &String )
{
  std::string s{ String };
  std::size_t Len = s.length();
  if( Len && s[Len - 1] != '/' )
    s.append( 1, '/' );
  return s;
}

inline void TrimTrailingSlash( std::string &String )
{
  std::size_t Len = String.length();
  if( Len > 1 && String[Len - 1] != '/' )
    String.resize( Len - 1 );
}

inline std::string &Append( std::string &Dest, const std::string &Source,
                            const std::string &Separator = Underscore )
{
  if( Source.length() )
  {
    if( Dest.length() )
      Dest.append( Separator );
    Dest.append( Source );
  }
  return Dest;
}

// Remove anything past the last delimeter from string, returning the removed part in suffix
// Return success / fail
bool ExtractSuffix( std::string &String, std::string &Suffix, const char * pszDelimeters = nullptr );

/// Replace the extension in FileName with NewExtension (appending if missing)
std::string ReplaceExtension( const std::string &FileName, const std::string &NewExtension );
void        ReplaceExtension(       std::string &FileName, const std::string &NewExtension );

// Remove the directory from the start of FileName (leave the trailing '/' in place)
std::string ExtractDirPrefix( std::string &FileName );
std::string GetDirPrefix( const std::string &FileName );

// Split String into an array using specified delimeters
std::vector<std::string> Split( const std::string &String, const char * pszDelimeters = nullptr );

// Extract suffix, then split strings. Default delimeters '.' and '_' respectively
bool ExtractSuffixSplit( std::string &String, std::vector<std::string> &Suffii,
                        const char * pszStringDelim = nullptr, const char * pszSuffixDelim = nullptr );

// Zipper merge v1 and v2 if same size (otherwise just append)
std::vector<std::string> ZipperMerge( const std::vector<std::string> &v1, const std::vector<std::string> &v2 );

// Dump the environment to stdout, prefixed by optional message
void DumpEnv(int argc, const char * const *argv, const char * pStr = nullptr );

MLU_String_hpp_end
#endif
