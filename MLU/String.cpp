/*****************************************************************************
 
 String manipulation

 Source file: String.cpp
 
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

#include "String.hpp"
#include <iostream>

MLU_String_hpp

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
extern const std::string EqualSign{ "=" };

// Default delimeters for the next couple of functions
extern const char szDefaultDelimeters[] = " \t,";

std::ostream & operator<<( std::ostream &os, const NegateStar &ns )
{
  switch( ns )
  {
    case NegateStar::None:
      break;
    case NegateStar::Negate:
      os << '-';
      break;
    case NegateStar::Star:
      os << '*';
      break;
    case NegateStar::NegateStar:
      os << "-*";
      break;
    default:
      os << "(NegateStarUnknown=" << static_cast<std::underlying_type_t<NegateStar>>( ns ) << ')';
      break;
  }
  return os;
}

std::istream & operator>>( std::istream& is, NegateStar &ns )
{
  char Buffer[2];
  int i = 0;
  while( i < 2 && !StreamEmpty( is ) )
  {
    Buffer[i] = is.peek();
    if( ( Buffer[i] == '-' || Buffer[i] == '*' ) && ( i == 0 || Buffer[i] != Buffer[i-1] ) )
    {
      is.get();
      ++i;
    }
    else
      break;
  }
  if( i == 0 )
    ns = NegateStar::None;
  else if( i == 2 )
    ns = NegateStar::NegateStar;
  else if( Buffer[0] == '-' )
    ns = NegateStar::Negate;
  else
    ns = NegateStar::Star;
  return is;
}

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

void ReplaceExtension( std::string &FileName, const std::string &NewExtension )
{
  std::size_t pos{ FileName.find_last_of( "/." ) };
  if( pos != std::string::npos && FileName[pos] == '.' )
    FileName.resize( pos + 1 );
  else
    FileName.append( 1, '.' );
  FileName.append( NewExtension );
}

std::string ReplaceExtension( const std::string &FileName, const std::string &NewExtension )
{
  std::string NewName;
  std::size_t pos{ FileName.find_last_of( "/." ) };
  if( pos != std::string::npos && FileName[pos] == '.' )
    NewName = FileName.substr( 0, pos + 1 );
  else
  {
    NewName = FileName;
    NewName.append( 1, '.' );
  }
  NewName.append( NewExtension );
  return NewName;
}

// Remove the directory from the start of FileName (leave the trailing '/' in place)
std::string ExtractDirPrefix( std::string &FileName )
{
  std::string Dir = GetDirPrefix( FileName );
  if( Dir.size() )
    FileName.erase( 0, Dir.size() );
  return Dir;
}

std::string GetDirPrefix( const std::string &FileName )
{
  std::string Dir;
  std::size_t pos{ FileName.find_last_of( '/' ) };
  if( pos != std::string::npos )
    Dir = FileName.substr( 0, pos + 1 );
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

MLU_String_hpp_end
