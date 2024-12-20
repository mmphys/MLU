/**
 
 Posix wrappers

 Source file: Posix.hpp

 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#ifndef MLU_Posix_hpp
#define MLU_Posix_hpp

#include <MLU/MLUFirst.hpp>

#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

// posix
#include <glob.h>
#include <sys/stat.h>
#include <unistd.h>

BEGIN_MLU_NAMESPACE

// Does the specified file exist?
bool FileExists( const std::string& Filename );

/// Wrapper for posix gethostname()
const std::string &GetHostName();

/// When not empty, specifies the value to be returned by GetHostName()
void SetHostName( const std::string &NewHostName );

/**
 Make the ancestor directories leading up to last element
 
 NB: Don't test whether this worked, so that it if this is done simultaneously by many threads, it will still work
 */
void MakeAncestorDirs( const std::string& Filename );

// Wrapper for posix glob
template<typename Iter>
std::vector<std::string> glob( const Iter &first, const Iter &last, const char * pszPrefix = nullptr )
{
  // Perform the glob
  std::string NameBuffer;
  if( pszPrefix && *pszPrefix )
    NameBuffer = pszPrefix;
  const std::size_t PrefixLen{ NameBuffer.length() };
  glob_t globBuf;
  memset( &globBuf, 0, sizeof( globBuf ) );
  int iFlags = GLOB_BRACE | GLOB_TILDE | GLOB_NOSORT; // | GLOB_NOCHECK | GLOB_NOMAGIC;
  for( Iter i = first; i != last; i++ )
  {
    NameBuffer.resize( PrefixLen );
    NameBuffer.append( *i );
    const int globResult{ glob( NameBuffer.c_str(), iFlags, NULL, &globBuf ) };
    if( !globResult )
      iFlags |= GLOB_APPEND;
    else if( globResult != GLOB_NOMATCH )
    {
      globfree( &globBuf );
      throw std::runtime_error( "glob() returned error " + std::to_string( globResult ) );
    }
  }

  // collect all the filenames into a std::list<std::string>
  std::vector<std::string> Filenames;
  Filenames.reserve( globBuf.gl_pathc );
  for( size_t i = 0; i < globBuf.gl_pathc; ++i )
    Filenames.push_back( std::string( globBuf.gl_pathv[i] ) );
  globfree(&globBuf);
  return Filenames;
}

END_MLU_NAMESPACE
#endif
