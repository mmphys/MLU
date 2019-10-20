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

#include <cassert>
#include <cctype>
#include <climits>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <H5Cpp.h>

#include <LatAnalyze/Statistics/Dataset.hpp>

// Default output file extension for binary data
#ifndef DEF_FMT
#define DEF_FMT "h5"
#endif
// Default output file extension for text based summaries of binary data
#ifndef TEXT_EXT
#define TEXT_EXT "txt"
#endif

#define BEGIN_COMMON_NAMESPACE namespace Common {
#define END_COMMON_NAMESPACE   };

BEGIN_COMMON_NAMESPACE

extern const std::string sBootstrap;
extern const std::string sModel;

using SeedType = unsigned int;

// Does the specified file exist?
inline bool FileExists(const std::string& Filename)
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

extern const double NaN;

class ValWithEr
{
public:
  double Central;
  double ErLow;
  double ErHigh;
  double Check;
  void Get( double dCentral, std::vector<double> &Data, std::size_t Count );
  ValWithEr() = default;
  ValWithEr( double dCentral, std::vector<double> &Data, std::size_t Count )
  { Get( dCentral, Data, Count ); }
};

inline std::ostream & operator<<( std::ostream &os, const ValWithEr &v )
{
  return os << v.Central << " " << v.ErLow << " " << v.ErHigh << " " << v.Check;
}

// Attributes for filenames in form base.type.seq.ext
struct FileNameAtt
{
  std::string Filename; // Full name
  std::string Base;
  std::string Type;
  SeedType    Seed;
  std::string Ext;
  std::vector<int> op;
  explicit FileNameAtt( const std::string &Filename );
  FileNameAtt( const std::string &Filename, std::vector<std::string> &OpNames );
};

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext);

// For now, this is just an alias
using Correlator = std::vector<std::complex<double>>;

inline void CopyCorrelator( Correlator &dst, const Correlator &src, int iOffset = 0, bool bSwapRealImag = false )
{
  const std::size_t Nt{ src.size() };
  if( Nt == 0 )
    throw new std::runtime_error( "Can't copy an uninitialised Correlator" );
  if( Nt > INT_MAX )
    throw new std::runtime_error( "Correlator size > INT_MAX" );
  if( dst.size() == 0 )
    dst.resize( Nt );
  else if( dst.size() != Nt )
    throw new std::runtime_error( "Can't assign correlator of length " + std::to_string( Nt ) + " to correlator of length " + std::to_string( dst.size() ) );
  const std::size_t dt{ ( iOffset < 0 ) ? Nt - ( -iOffset % Nt ) : iOffset % Nt };
  for( std::size_t t = 0; t < Nt; t++ )
  {
    const std::complex<double> & z{ src[( t + dt ) % Nt] };
    if( bSwapRealImag ) {
      dst[t].real( z.imag() );
      dst[t].imag( z.real() );
    }
    else
      dst[t] = z;
  }
}

// A vector of correlators
class vCorrelator : public std::vector<Correlator>
{
public:
  using base_type = typename std::vector<Correlator>;
  using size_type = base_type::size_type;
  void ResizeCorrelators( int Nt )
  {
    for( Correlator &c : * this )
      c.resize( Nt );
  }
  // This next probably isn't needed - the assign will set the size of each Correlator
  //vCorrelator( size_type Size, int Nt ) : base_type( Size )
  //{ ResizeCorrelators( Nt ); }
};

// Compare two strings, case insensitive
inline bool EqualIgnoreCase(const std::string & s1, const std::string & s2)
{
  const std::size_t Len{ s1.size() };
  bool bEqual = ( s2.size() == Len );
  for( std::size_t i = 0; bEqual && i < Len; i++ )
    bEqual = ( s1[ i ] == s2[ i ] ) || ( std::toupper( s1[ i ] ) == std::toupper( s2[ i ] ) );
  return bEqual;
}

// Generic conversion from a string to any type
template<typename T> inline T FromString( const std::string &String ) {
  T t;
  std::istringstream iss( String );
  if( !( iss >> t ) || ( iss >> std::ws && !iss.eof() ) )
    throw new std::invalid_argument( String );
  return t;
}
// Converting a string to a string makes a copy
template<> inline std::string FromString<std::string>( const std::string &String )
  { return String; }

// Generic representation of momentum
struct Momentum
{
  const int x;
  const int y;
  const int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
  inline bool isZero() const { return x==0 && y==0 && z==0; }
  std::string to_string( const std::string &separator, bool bNegative = false ) const {
    return std::to_string( bNegative ? -x : x ) + separator
    + std::to_string( bNegative ? -y : y ) + separator
    + std::to_string( bNegative ? -z : z );
  }
};

// My implementation of H5File - adds a definition of complex type
class H5File : public H5::H5File {
public:
  using H5::H5File::H5File; // Inherit the base class's constructors
  // Return the same HDF5 complex type Grid uses
  static H5::CompType & ComplexType( void );
};

// Get first groupname from specified group
std::string GetFirstGroupName( H5::Group & g );

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

struct TrajFile
{
  const std::string Filename;  // Name of the file
  explicit TrajFile( const std::string & Filename_ ) : Filename{ Filename_ } {}
};

// This describes one contraction and each of its trajectory files
struct TrajList
{
  std::string Name;                     // name of the contraction
  std::map<int, std::unique_ptr<TrajFile>> TrajFile;  // list of all the files, sorted by trajectory number
  TrajList( const std::string &Name_ ) : Name{ Name_ } {}
};

// This is a list of all the contractions we've been asked to process
class Manifest : public std::map<std::string, TrajList> {
public:
  // Process list of files on the command-line, breaking them up into individual trajectories
  Manifest( const std::vector<std::string> &Files, const std::string &sIgnore );
};

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
