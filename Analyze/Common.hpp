/*************************************************************************************
 
 Common utilities (no dependencies other than c++ stdlib)
 
 Source file: Common.hpp
 
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
 *************************************************************************************/
/*  END LEGAL */

#ifndef Common_hpp
#define Common_hpp

#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <H5Cpp.h>

#define BEGIN_COMMON_NAMESPACE namespace Common {
#define END_COMMON_NAMESPACE   };

BEGIN_COMMON_NAMESPACE

class H5File : public H5::H5File {
public:
  using H5::H5File::H5File; // Inherit the base class's constructors
  // Return the same HDF5 complex type Grid uses
  static H5::CompType & ComplexType( void );
};

// Get first groupname from specified group
std::string GetFirstGroupName( H5::Group & g );

// Does the specified file exist?
bool FileExists(const std::string& Filename);

// Are all the floating point numbers pointed to finite
template <typename T, typename I> inline bool IsFinite( const T * d, I n )
{
  while( n-- )
    if( !std::isfinite( *d++ ) )
      return false;
  return true;
}

// Are all the floating point numbers in this vector finite
template <typename T> inline bool IsFinite( const std::vector<T> & v )
{
  for( const T n : v )
    if( !std::isfinite( n ) )
      return false;
  return true;
}

// Read a complex array from an HDF5 file

void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      const std::string &GroupName  = std::string(),
                      const std::string &ObjectName = std::string( "correlator" ) );

// This describes one contraction and each of its trajectory files
struct TrajList
{
  std::string Name;                     // name of the contraction
  std::map<int, std::string> Filename;  // list of all the files, sorted by trajectory number
};

// This is a list of all the contractions we've been asked to process
class Manifest : public std::map<std::string, TrajList> {
public:
  // Process list of files on the command-line, breaking them up into individual trajectories
  explicit Manifest(const std::vector<std::string> &Files, const std::string &sIgnore = "");
};

// Generic conversion from a string to any type
template<typename T> inline T FromString( const std::string &String ) {
  T t;
  std::istringstream iss( String );
  if( !( iss >> t ) || ( iss >> std::ws && !iss.eof() ) ) {
    std::cerr << "Error converting \"" << String << "\"" << std::endl;
    exit( EXIT_FAILURE );
  }
  return t;
}
template<> inline std::string FromString<std::string>( const std::string &String ) {
  return String; // Yes, this makes a copy
}

struct CommandLine;
std::ostream& operator<<( std::ostream& os, const CommandLine &cl);

struct CommandLine {
  using SwitchMap = std::map<std::string, std::vector<std::string>>;
  using SwitchPair = std::pair<std::string, std::vector<std::string>>;
  
  enum SwitchType { Flag, Single, Multiple };
  struct SwitchDef {
    const char * Switch;
    SwitchType   Type;
    const char * Default;
    SwitchDef( const char * switch_, SwitchType type_ = Single, const char * default_ = nullptr )
    : Switch{ switch_ }, Type{ type_ }, Default{ default_ } {}
  };
  
  std::vector<std::string> Args;
  SwitchMap                Switches;
  
private:
  static void SkipPastSep( const char * & p );
  
public:
  CommandLine( int argc, char *argv[], const std::vector<SwitchDef> &defs );

  inline bool GotSwitch( const std::string &SwitchName ) {
    return Switches.find( SwitchName ) != Switches.end(); }

  template<typename T> inline T SwitchValue( const std::string &Switch, int Subscript = 0 ) {
    SwitchMap::const_iterator it = Switches.find( Switch );
    if( it == Switches.end() ) {
      std::cerr << "Switch \"" << Switch << "\" not found" << std::endl;
      exit( EXIT_FAILURE );
    }
    const std::vector<std::string> &v{ it->second };
    if( static_cast<std::size_t>( Subscript ) >= v.size() ) {
      std::cerr << "Switch \"" << Switch << "\", Subscript " << Subscript << " not found" << std::endl;
      exit( EXIT_FAILURE );
    }
    return FromString<T>( v[Subscript] );
  }
};

END_COMMON_NAMESPACE
#endif // Common_hpp
