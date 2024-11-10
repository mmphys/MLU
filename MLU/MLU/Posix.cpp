/*****************************************************************************
 
 Posix wrappers

 Source file: Posix.cpp
 
 Copyright (C) 2019-2023

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

#include <MLUconfig.h>
#include "Posix.hpp"

BEGIN_MLU_NAMESPACE

/*****************************************************************************

 Does the specified file exist?

*****************************************************************************/

bool FileExists( const std::string& Filename )
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

/*****************************************************************************

 Allow host name to be get and set

*****************************************************************************/

static std::string MyHostName;

// Wrapper for posix gethostname()
void SetHostName( const std::string &NewHostName )
{
  MyHostName = NewHostName;
}

// Wrapper for posix gethostname()
const std::string &GetHostName()
{
  if( MyHostName.empty() )
  {
    // Get hostname from environment variable
    const char * pDefaultString{ std::getenv( "MLUHost" ) };
    if( pDefaultString && * pDefaultString )
      MyHostName = pDefaultString;
    else
    {
      // Otherwise get the real hostname
      char Buffer[256];
      const int BufLen{ sizeof( Buffer ) - 1 };
      if( gethostname( Buffer, BufLen ) )
        throw std::runtime_error( "gethostname() returned error " + std::to_string( errno ) );
      Buffer[BufLen] = 0;
      MyHostName = std::string( Buffer );
      // Trim the hostname if it ends in ".local" or ".home"
      std::size_t pos = MyHostName.find_last_of( "." );
      if( pos != std::string::npos && pos && pos != MyHostName.length() - 1
         && ( MyHostName.compare( pos + 1, std::string::npos, "local" ) == 0
             || MyHostName.compare( pos + 1, std::string::npos, "home" ) == 0 ) )
        MyHostName.resize( pos - 1 );
    }
  }
  return MyHostName;
}

/*****************************************************************************

 Make the ancestor directories leading up to last element
 NB: Don't test whether this worked, so that it if this is done simultaneously by many threads, it will still work

*****************************************************************************/

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

END_MLU_NAMESPACE
