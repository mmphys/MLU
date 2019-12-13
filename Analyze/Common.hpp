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

// test whether a type is complex
template<typename T> struct is_complex                  : public std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : public std::true_type {};

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
  inline explicit operator bool() const { return x!=0 || y!=0 || z!=0; }
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

// Open the specified HDF5File and group
void OpenHdf5FileGroup(H5::H5File &f, H5::Group &g, const std::string &FileName, std::string &GroupName, unsigned int flags = H5F_ACC_RDONLY );

// Traits for samples
template<typename ST, typename V = void> struct SampleTraits : public std::false_type{};
template<typename ST> struct SampleTraits<ST, typename std::enable_if<std::is_floating_point<ST>::value>::type> : public std::true_type
{
  static constexpr bool is_complex = false;
  using scalar_type = ST;
};
template<typename ST> struct SampleTraits<std::complex<ST>, typename std::enable_if<SampleTraits<ST>::value>::type> : public std::true_type
{
  static constexpr bool is_complex = true;
  using scalar_type = typename SampleTraits<ST>::scalar_type;
};

// This is for a sample of correlators. There is a central replica, specified aux replicas and variable number of samples

template <typename T, int NumAuxSamples_ = 0>
class Sample
{
  int NumSamples_;
  int Nt_;
  std::unique_ptr<T[]> m_pData;
public:
  static constexpr int NumAuxSamples{ NumAuxSamples_ };
  static constexpr int NumExtraSamples{ NumAuxSamples + 1 }; // Auxiliary samples + central replica
  static constexpr int idxCentral{ -1 };
  static constexpr int idxAux{ idxCentral - NumAuxSamples };
  inline int NumSamples() const { return NumSamples_; }
  inline int Nt() const { return Nt_; }
  void resize( int NumSamples, int Nt )
  {
    NumSamples_ = NumSamples;
    Nt_ = Nt;
    m_pData.reset( new T[ static_cast<std::size_t>( NumSamples + NumExtraSamples ) * Nt ] );
  }
  void Read ( const std::string &FileName, std::string &GroupName );
  void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  Sample( int NumSamples, int Nt ) { resize( NumSamples, Nt ); }
  Sample( const std::string &FileName, std::string &GroupName ) { Read( FileName, GroupName ); }
  T * operator[]( int Sample )
  {
    if( Sample < idxAux || Sample > NumSamples_ )
      throw std::out_of_range( "Sample " + std::to_string( Sample ) );
    return & m_pData[static_cast<std::size_t>( Sample + NumExtraSamples ) * Nt_];
  }
  const T * operator[]( int Sample ) const
  {
    if( Sample < idxAux || Sample > NumSamples_ )
      throw std::out_of_range( "Sample " + std::to_string( Sample ) );
    return & m_pData[static_cast<std::size_t>( Sample + NumExtraSamples ) * Nt_];
  }
};

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T, int NumAuxSamples_>
void Sample<T, NumAuxSamples_>::Read( const std::string &FileName, std::string &GroupName )
{
  H5::H5File f;
  H5::Group  g;
  OpenHdf5FileGroup( f, g, FileName, GroupName );
  bool bOK = false;
  H5E_auto2_t h5at;
  void      * f5at_p;
  H5::Exception::getAutoPrint(h5at, &f5at_p);
  H5::Exception::dontPrint();
  try // to load from LatAnalyze format
  {
    unsigned short att_Type;
    H5::Attribute a;
    a = g.openAttribute("type");
    a.read( H5::PredType::NATIVE_USHORT, &att_Type );
    a.close();
    if( att_Type == 2 )
    {
      unsigned long  att_nSample;
      a = g.openAttribute("nSample");
      a.read( H5::PredType::NATIVE_ULONG, &att_nSample );
      a.close();
      std::vector<double> buffer;
      for( int i = idxCentral; i < NumSamples_; i++ )
      {
        H5::DataSet ds = g.openDataSet( "data_" + ( i == idxCentral ? "C" : "S_" + std::to_string( i ) ) );
        H5::DataSpace dsp = ds.getSpace();
        int nDims{ dsp.getSimpleExtentNdims() };
        if( nDims == 2 )
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
          ds.read( buffer.data(), H5::PredType::NATIVE_DOUBLE );
          if( !IsFinite( buffer ) )
            break;
          T * p { (*this)[i] };
          for( int t = 0; t < Nt_ ; t++ )
          {
            p->real( buffer[t] );
            p->imag( buffer[t + Nt_] );
            p++;
          }
          // Check whether this is the end
          if( i == NumSamples_ - 1 )
          {
            bOK = true;
            // This format has no auxiliary rows - set them to zero
            for( int j = idxAux; j < idxCentral; j++ )
            {
              p = (*this)[j];
              for( int t = 0; t < Nt_ ; t++ )
                p = 0;
            }
          }
        }
      }
    }
  }
  catch(const H5::Exception &)
  {
    bOK = false;
    H5::Exception::clearErrorStack();
  }
  if( !bOK )
  {
    try // to load from my format
    {
      unsigned short att_nAux;
      H5::Attribute a;
      a = g.openAttribute("nAux");
      a.read( H5::PredType::NATIVE_USHORT, &att_nAux );
      a.close();
      if( att_nAux == NumAuxSamples )
      {
        unsigned int att_nSample;
        a = g.openAttribute("nSample");
        a.read( H5::PredType::NATIVE_UINT, &att_nSample );
        a.close();
        H5::DataSet ds = g.openDataSet( "data_C" );
        H5::DataSpace dsp = ds.getSpace();
        int nDims{ dsp.getSimpleExtentNdims() };
        if( nDims == 1 )
        {
          hsize_t Dim[2];
          dsp.getSimpleExtentDims( Dim );
          if( Dim[0] * att_nSample <= std::numeric_limits<int>::max() )
          {
            resize( static_cast<int>( att_nSample ), static_cast<int>( Dim[0] ) );
            ds.read( (*this)[idxCentral], H5File::ComplexType() );
            dsp.close();
            ds.close();
            ds = g.openDataSet( "data_S" );
            dsp = ds.getSpace();
            nDims = dsp.getSimpleExtentNdims();
            if( nDims == 2 )
            {
              dsp.getSimpleExtentDims( Dim );
              if( Dim[0] == att_nSample && Dim[1] == static_cast<unsigned int>( Nt_ ) )
              {
                ds.read( (*this)[0], H5File::ComplexType() );
                if( NumAuxSamples == 0 )
                  bOK = true;
                else
                {
                  dsp.close();
                  ds.close();
                  ds = g.openDataSet( "data_Aux" );
                  dsp = ds.getSpace();
                  nDims = dsp.getSimpleExtentNdims();
                  if( nDims == 2 )
                  {
                    dsp.getSimpleExtentDims( Dim );
                    if( Dim[0] == static_cast<unsigned int>( NumAuxSamples ) && Dim[1] == static_cast<unsigned int>( Nt_ ) )
                    {
                      ds.read( (*this)[idxAux], H5File::ComplexType() );
                      bOK = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    catch(const H5::Exception &)
    {
      bOK = false;
      H5::Exception::clearErrorStack();
    }
  }
  H5::Exception::setAutoPrint(h5at, f5at_p);
  if( !bOK )
    throw std::runtime_error( "Unable to read sample from " + FileName );
}

template <typename T, int NumAuxSamples_>
void Sample<T, NumAuxSamples_>::Write( const std::string &FileName, const char * pszGroupName )
{
  const std::string GroupName{ pszGroupName == nullptr || *pszGroupName == 0 ? "/bootstrap" : pszGroupName };
  bool bOK = false;
  try // to write in my format
  {
    H5::H5File f( FileName, H5F_ACC_TRUNC );
    H5::Group  g = f.openGroup( "/" );
    hsize_t Dims[2];
    Dims[0] = 1;
    H5::DataSpace ds1( 1, Dims );
    H5::Attribute a = g.createAttribute( "_Grid_dataset_threshold", H5::PredType::STD_U32LE, ds1);
    int tmp{ 6 };
    a.write( H5::PredType::NATIVE_INT, &tmp );
    a.close();
    g.close();
    g = f.createGroup( GroupName );
    a = g.createAttribute( "nAux", H5::PredType::STD_U16LE, ds1);
    tmp = NumAuxSamples;
    a.write( H5::PredType::NATIVE_INT, &tmp );
    a.close();
    a = g.createAttribute( "nSample", H5::PredType::STD_U32LE, ds1);
    a.write( H5::PredType::NATIVE_INT, &NumSamples_ );
    a.close();
    ds1.close();
    Dims[0] = Nt_;
    ds1 = H5::DataSpace( 1, Dims );
    H5::DataSet ds = g.createDataSet( "data_C", H5File::ComplexType(), ds1 );
    ds.write( (*this)[idxCentral], H5File::ComplexType() );
    ds.close();
    ds1.close();
    Dims[0] = NumSamples_;
    Dims[1] = Nt_;
    ds1 = H5::DataSpace( 2, Dims );
    ds = g.createDataSet( "data_S", H5File::ComplexType(), ds1 );
    ds.write( (*this)[0], H5File::ComplexType() );
    ds.close();
    ds1.close();
    if( NumAuxSamples )
    {
      Dims[0] = NumAuxSamples;
      ds1 = H5::DataSpace( 2, Dims );
      ds = g.createDataSet( "data_Aux", H5File::ComplexType(), ds1 );
      ds.write( (*this)[idxAux], H5File::ComplexType() );
      ds.close();
      ds1.close();
    }
    bOK = true;
  }
  catch(const H5::Exception &)
  {
    bOK = false;
  }
  if( !bOK )
    throw std::runtime_error( "Unable to write sample to " + FileName + ", group " + GroupName );
}

// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      std::string &GroupName,
                      const std::string &ObjectName = std::string( "correlator" ) );

void ReadRealArray(std::vector<double> &buffer, const std::string &FileName,
                   std::string &GroupName,
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

  inline int NumValues( const std::string &Switch ) {
    int iNumValues{ 0 };
    SwitchMap::const_iterator it = Switches.find( Switch );
    if( it != Switches.end() ) {
      iNumValues = static_cast<int>( it->second.size() );
    }
    return iNumValues;
  }
  
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

std::ostream& operator<<( std::ostream& os, const CommandLine &cl);

END_COMMON_NAMESPACE
#endif // Common_hpp
