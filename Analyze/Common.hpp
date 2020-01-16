/**

 Mike's lattice QCD utilities
 
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
**/

// Common utilities (no dependencies other than c++ stdlib)

#ifndef Common_hpp
#define Common_hpp

#include <array>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <regex>
#include <string>
#include <vector>

// HDF5 Library
#include <H5Cpp.h>
#include <H5CompType.h>

// Eigen dense matrices
#include <Eigen/Dense>

//#include <LatAnalyze/Statistics/Dataset.hpp>

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

// Text required for summaries of correlators
namespace CorrSumm {
  extern const char sep[];
  extern const char Comment[];
  extern const char NewLine[];

  static constexpr int NumSummaries{ 4 };
  extern const char * SummaryNames[NumSummaries];

  extern const char FieldNames[];
  extern const char FieldNames2[];
  extern const char * SummaryHeader[NumSummaries];
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

namespace Gamma
{
  enum class Algebra
  {
    Unknown = -1,
    MinusGamma5,        // 0
    Gamma5,             // 1
    MinusGammaT,        // 2
    GammaT,             // 3
    MinusGammaTGamma5,  // 4
    GammaTGamma5,       // 5
    MinusGammaX,        // 6
    GammaX,             // 7
    MinusGammaXGamma5,  // 8
    GammaXGamma5,       // 9
    MinusGammaY,        // 10
    GammaY,             // 11
    MinusGammaYGamma5,  // 12
    GammaYGamma5,       // 13
    MinusGammaZ,        // 14
    GammaZ,             // 15
    MinusGammaZGamma5,  // 16
    GammaZGamma5,       // 17
    MinusIdentity,      // 18
    Identity,           // 19
    MinusSigmaXT,       // 20
    SigmaXT,            // 21
    MinusSigmaXY,       // 22
    SigmaXY,            // 23
    MinusSigmaXZ,       // 24
    SigmaXZ,            // 25
    MinusSigmaYT,       // 26
    SigmaYT,            // 27
    MinusSigmaYZ,       // 28
    SigmaYZ,            // 29
    MinusSigmaZT,       // 30
    SigmaZT,            // 31
  };
  static constexpr unsigned int                nGamma = 32;
  extern const std::array<std::string, nGamma> name;      // Long name, per Grid
  extern const std::array<std::string, nGamma> nameShort; // My abbreviations
  std::string NameShort( Algebra alg, const char * pszPrefix = nullptr );
  std::ostream& operator<<(std::ostream& os, const Algebra &a);
  std::istream& operator>>(std::istream& is, Algebra &a);
};

extern const std::string sBootstrap;
extern const std::string sModel;

using SeedType = unsigned int;

// Does the specified file exist?
bool FileExists( const std::string& Filename );

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

// Are all the floating point numbers in this matrix finite
template <typename T> inline bool IsFinite( const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & m, bool bDiagonalsOnly = false )
{
  for( Eigen::Index row = 0; row < m.rows(); ++row )
    for( Eigen::Index col = 0; col < m.cols(); ++col )
      if( ( !bDiagonalsOnly || row == col ) && !std::isfinite( m( row, col ) ) )
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

// Generic representation of momentum
struct Momentum
{
  int x;
  int y;
  int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
  inline explicit operator bool() const { return x!=0 || y!=0 || z!=0; }
  inline int p2() const { return x * x + y * y + z * z; }
  inline std::string p2_string( const std::string &separator ) const {
    std::string s{ separator };
    s.append( "p2" );
    s.append( separator );
    s.append( std::to_string( p2() ) );
    return s;
  }
  std::string to_string( const std::string &separator, bool bNegative = false ) const {
    return std::to_string( bNegative ? -x : x ) + separator
    + std::to_string( bNegative ? -y : y ) + separator
    + std::to_string( bNegative ? -z : z );
  }
  std::string to_string4d( const std::string &separator, bool bNegative = false ) const {
    return to_string( separator, bNegative ) + separator + "0";
  }
};

// Attributes for filenames in form base.type.seed.ext
struct FileNameAtt
{
  std::string Filename; // Full (and unadulterated) original filename
  std::string Dir;      // Directory part of the filename (with trailing '/')
  std::string Base;     // Base of the filename
  std::string Type;
  std::string SeedString;
  bool        bSeedNum = false;
  SeedType    Seed = 0; // Numeric value of SeedString, but only if bSeedNum == true
  std::string Ext;
  std::vector<int> op;
  // reset contents
  void clear()
  {
    Filename.clear();
    Dir.clear();
    Base.clear();
    Type.clear();
    SeedString.clear();
    bSeedNum = false;
    Seed = 0;
    Ext.clear();
    op.clear();
  }

  void Parse( const std::string &Filename_ );
  FileNameAtt() = default;
  explicit FileNameAtt( const std::string &Filename ) { Parse( Filename ); }
  FileNameAtt( const std::string &Filename, std::vector<std::string> &OpNames );
};

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext);

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &s, bool &bHasTimeslice, int & Timeslice );

// Strip out momentum info from string if present
void ExtractP2( std::string &s, bool &bGotMomentum, int &P2, bool IgnoreSubsequentZero = true );

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

// Generic conversion from a string to any type
template<typename T> inline T FromString( const std::string &String ) {
  T t;
  std::istringstream iss( String );
  if( !( iss >> t ) || ( iss >> std::ws && !iss.eof() ) )
    throw std::invalid_argument( "Argument \"" + String + "\" is not type " + typeid(T).name() );
  return t;
}

// Converting a string to a string makes a copy
template<> inline std::string FromString<std::string>( const std::string &String )
  { return String; }

// Generic conversion from a string to an array of any type (comma or space separated)
template<typename T> inline std::vector<T> ArrayFromString( const std::string &String ) {
  std::string s{ String };
  for( std::size_t pos = 0; ( pos = s.find( ',', pos ) ) != std::string::npos; )
    s[pos] = ' ';
  std::istringstream iss( s );
  T t;
  std::vector<T> v;
  while( iss >> std::ws && !iss.eof() )
  {
    if( !( iss >> t ) )
      throw std::invalid_argument( "ArrayFromString: \"" + String + "\" is not type " + typeid(T).name() );
    v.push_back( t );
  }
  return v;
}

// My implementation of H5File - adds a definition of complex type
class H5File : public H5::H5File {
public:
  using H5::H5File::H5File; // Inherit the base class's constructors
  // Return the same HDF5 complex type Grid uses
  static H5::CompType & ComplexType( int fpsize = sizeof( double ) );
};

// Get first groupname from specified group
std::string GetFirstGroupName( H5::Group & g );

// Open the specified HDF5File and group
void OpenHdf5FileGroup(H5::H5File &f, H5::Group &g, const std::string &FileName, std::string &GroupName, unsigned int flags = H5F_ACC_RDONLY );

// Read the gamma algebra attribute string and make sure it's valid
Gamma::Algebra ReadGammaAttribute( H5::Group &g, const char * pAttName );

// Enumeration describing whether signal is in complex or real part
enum class Signal
{
  Unknown = 0,
  Real,
  Imag
};

std::ostream& operator<<(std::ostream& os, const Signal &sig);
std::istream& operator>>(std::istream& is, Signal &sig);

// Traits for samples
// SampleTraits<ST>::value is true for supported types (floating point and complex)
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

// Correlator file. Could be either single correlator, or multiple gammas

template <typename T>
class CorrelatorFile
{
  // Data members
public:
  using scalar_type = typename SampleTraits<T>::scalar_type;
  static constexpr int scalar_size { sizeof( scalar_type ) };
  static constexpr bool is_complex { SampleTraits<T>::is_complex };
private:
  int NumOps_ = 0;
  int Nt_ = 0;
  std::unique_ptr<T[]> m_pData;
  std::vector<Gamma::Algebra> Alg_;
public:
  FileNameAtt Name_;
  //std::string Prefix; // this is the processed prefix name
  bool bHasTimeslice = false;
  int  Timeslice_ = 0;
  // Swap function so that this type is sortable
  /*friend void swap( CorrelatorFile &l, CorrelatorFile &r )
  {
    int i;
    i = l.NumOps_;
    l.NumOps_ = r.NumOps_;
    r.NumOps_ = i;
    i = l.Nt_;
    l.Nt_ = r.Nt_;
    r.Nt_ = i;
    std::swap( l.m_pData, r.m_pData );
    std::swap( l.Alg_, r.Alg_ );
    //swap( l.Name_, r.Name_ );
    std::swap( l.Prefix, r.Prefix );
    i = l.Timeslice_;
    l.Timeslice_ = r.Timeslice_;
    r.Timeslice_ = i;
    bool b = l.bHasTimeslice;
    l.bHasTimeslice = r.bHasTimeslice;
    r.bHasTimeslice = b;
  }*/
  // Member functions
private:
  inline void RangeCheck( int Sample ) const
  {
    if( Sample < 0 || Sample >= NumOps_ * NumOps_ )
      throw std::out_of_range( "Sample " + std::to_string( Sample ) );
  }
  inline int GammaIndex( Gamma::Algebra g ) const
  {
    int idx;
    for( idx = 0; idx < Alg_.size() && Alg_[idx] != g; idx++ )
      ;
    return idx;
  }
public:
  template <typename DontSpecialise=T>
  inline typename std::enable_if<!SampleTraits<DontSpecialise>::is_complex, Signal>::type
    getSignal( int Sample ) const { return ( NumOps_==0 || Nt_==0 ) ? Signal::Unknown : Signal::Real; }
  template <typename DontSpecialise=T>
  inline typename std::enable_if<SampleTraits<DontSpecialise>::is_complex, Signal>::type
    getSignal( int Sample ) const
  {
    if( NumOps_ == 0 || Nt_ == 0 )
      return Signal::Unknown;
    const T * p = (*this)[Sample];
    int cntReal{ 0 };
    int cntImag{ 0 };
    for( int i = 0; i < Nt_; i++ )
    {
      if( std::abs( p->real() ) >= std::abs( p->imag() ) )
        cntReal++;
      else
        cntImag++;
      p++;
    }
    const int Quorum{ Nt_ * 2 / 3 };
    if( cntReal > Quorum )
      return Signal::Real;
    if( cntImag > Quorum )
      return Signal::Imag;
    return Signal::Unknown;
  }
  inline Signal getSignal( Gamma::Algebra gSink, Gamma::Algebra gSource ) const
  {
    return getSignal( GammaIndex( gSink ) * NumOps_ + GammaIndex( gSource ) );
  }
  inline int NumOps() const { return NumOps_; }
  inline int Nt() const { return Nt_; }
  inline int Timeslice() const { return bHasTimeslice ? Timeslice_ : 0; }
  inline const std::vector<Gamma::Algebra> &Alg() { return Alg_; }
  void resize( int NumOps, int Nt )
  {
    if( NumOps_ != NumOps )
      Alg_.clear();
    if( NumOps_ != NumOps || Nt_ != Nt )
    {
      NumOps_ = NumOps;
      Nt_ = Nt;
      if( NumOps_ == 0 || Nt_ == 0 )
        m_pData.reset( nullptr );
      else
        m_pData.reset( new T[ static_cast<std::size_t>( NumOps * NumOps ) * Nt ] );
    }
  }
  void Read ( const std::string &FileName, std::string &GroupName, std::vector<Gamma::Algebra> &Alg,
             const int * pTimeslice = nullptr );
  //void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void WriteSummaries ( const std::string &Prefix, const std::vector<Common::Gamma::Algebra> &AlgSpecific );
  T * operator[]( int Sample )
  {
    RangeCheck( Sample );
    return & m_pData[static_cast<std::size_t>( Sample ) * Nt_];
  }
  const T * operator[]( int Sample ) const
  {
    RangeCheck( Sample );
    return & m_pData[static_cast<std::size_t>( Sample ) * Nt_];
  }
  T * operator()( Gamma::Algebra gSink, Gamma::Algebra gSource )
  { return (*this)[ GammaIndex( gSink ) * NumOps_ + GammaIndex( gSource ) ]; }
  const T * operator()( Gamma::Algebra gSink, Gamma::Algebra gSource ) const
  { return (*this)[ GammaIndex( gSink ) * NumOps_ + GammaIndex( gSource ) ]; }
  bool IsFinite() { return Common::IsFinite( reinterpret_cast<scalar_type *>( m_pData.get() ),
      static_cast<size_t>( NumOps_ * NumOps_ * ( SampleTraits<T>::is_complex ? 2 : 1 ) ) * Nt_ ); }
  // Constructors (copy operations missing for now - add them if they become needed)
  CorrelatorFile() {}
  CorrelatorFile( CorrelatorFile && ) = default; // Move constructor
  CorrelatorFile(const std::string &FileName, std::string &GroupName, std::vector<Gamma::Algebra> &Alg,
                 int * pTimeslice = nullptr )
  { Read( FileName, GroupName, Alg, pTimeslice ); }
  // Operators
  CorrelatorFile& operator=(CorrelatorFile && r) = default; // Move assignment
};

using CorrelatorFileC = CorrelatorFile<std::complex<double>>;
using CorrelatorFileD = CorrelatorFile<double>;

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void CorrelatorFile<T>::Read( const std::string &FileName, std::string &GroupName,
                              std::vector<Gamma::Algebra> &Alg, const int * pTimeslice )
{
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
  H5::H5File f;
  H5::Group  g;
  OpenHdf5FileGroup( f, g, FileName, GroupName );
  bool bOK = false;
  H5E_auto2_t h5at;
  void      * f5at_p;
  H5::Exception::getAutoPrint(h5at, &f5at_p);
  H5::Exception::dontPrint();
  try // to load a single correlator
  {
    if( Alg.size() == 0 || ( Alg.size() == 1 && Alg[0] == Gamma::Algebra::Unknown ) )
    {
      H5::DataSet ds = g.openDataSet( "correlator" );
      H5::DataSpace dsp = ds.getSpace();
      int nDims{ dsp.getSimpleExtentNdims() };
      if( nDims == 1 )
      {
        hsize_t Dim[1];
        dsp.getSimpleExtentDims( Dim );
        if( Dim[0] <= std::numeric_limits<int>::max() )
        {
          resize( 1, static_cast<int>( Dim[0] ) );
          Alg_.resize( 1 );
          Alg_[0] = Gamma::Algebra::Unknown;
          ds.read( (*this)[0], H5File::ComplexType( scalar_size ) );
          bOK = true;
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
    try // to load from array of correlators indexed by gamma matrix
    {
      unsigned short NumVec;
      H5::Attribute a;
      a = g.openAttribute("_Grid_vector_size");
      a.read( H5::PredType::NATIVE_USHORT, &NumVec );
      a.close();
      // Must be a perfect square and have at least as many as entries as requested
      const unsigned short NumFileOps{static_cast<unsigned short>( std::sqrt( NumVec ) + 0.5 )};
      if( NumFileOps * NumFileOps == NumVec && NumFileOps >= Alg.size() )
      {
        bOK = true;
        std::vector<int> count;
        for( unsigned short i = 0; bOK && i < NumVec; i++ )
        {
          bOK = false;
          H5::Group gi = g.openGroup( GroupName + "_" + std::to_string( i ) );
          H5::DataSet ds = gi.openDataSet( "corr" );
          H5::DataSpace dsp = ds.getSpace();
          int nDims{ dsp.getSimpleExtentNdims() };
          if( nDims == 1 )
          {
            hsize_t Dim[1];
            dsp.getSimpleExtentDims( Dim );
            if( Dim[0] <= std::numeric_limits<int>::max() )
            {
              const int ThisNt{ static_cast<int>( Dim[0] ) };
              if( i == 0 )
              {
                // First correlator - resize and save operators if known
                resize( Alg.size() ? Alg.size() : NumFileOps, ThisNt );
                count.resize( NumOps_ * NumOps_, 0 );
                if( Alg.size() )
                {
                  Alg_.resize( Alg.size() );
                  std::copy( Alg.cbegin(), Alg.cend(), Alg_.begin() );
                }
                else
                {
                  Alg_.clear();
                  Alg_.reserve( Alg.size() );
                }
              }
              else if( ThisNt != Nt_ )
              {
                break;
              }
              // Read the gamma algebra strings and make sure they are valid
              const Gamma::Algebra gSnk{ ReadGammaAttribute( gi, "gamma_snk" ) };
              int idxSnk;
              for( idxSnk = 0; idxSnk < Alg_.size() && Alg_[idxSnk] != gSnk; idxSnk++ )
                ;
              if( idxSnk == Alg_.size() && Alg_.size() < NumOps_ )
                Alg_.push_back( gSnk );
              bOK = true; // We can safely ignore gamma structures we're not interested in
              if( idxSnk < Alg_.size() )
              {
                const Gamma::Algebra gSrc{ ReadGammaAttribute( gi, "gamma_src" ) };
                int idxSrc;
                for( idxSrc = 0; idxSrc < Alg_.size() && Alg_[idxSrc] != gSrc; idxSrc++ )
                  ;
                if( idxSrc == Alg_.size() && Alg_.size() < NumOps_ )
                  Alg_.push_back( gSrc );
                if( idxSrc < Alg_.size() )
                {
                  const int idx{ idxSnk * NumOps_ + idxSrc };
                  ds.read( (*this)[idx], H5File::ComplexType( scalar_size ) );
                  count[idx]++;
                }
              }
            }
          }
        }
        // Make sure that everything we wanted was loaded once and only once
        for( int i = 0; bOK && i < NumOps_ * NumOps_; i++ )
          if( count[i] != 1 )
            bOK = false;
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
  if( !IsFinite() )
    throw std::runtime_error( "Values read are not all finite" );
  // If I'm discovering which operators are in the file, copy them out
  if( Alg.empty() )
  {
    Alg.resize( Alg_.size() );
    std::copy( Alg_.cbegin(), Alg_.cend(), Alg.begin() );
  }
}

template <typename T>
void CorrelatorFile<T>::WriteSummaries ( const std::string &Prefix, const std::vector<Gamma::Algebra> &AlgSpecific )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  const int nt{ Nt() };
  const std::vector<Gamma::Algebra> &Alg{ AlgSpecific.size() ? AlgSpecific : Alg_ };
  const int NumOps{ static_cast<int>( Alg.size() ) };
  for( int Snk = 0; Snk < NumOps; Snk++ )
  {
    static const char pszSep[] = "_";
    std::string sSnk{ Common::Gamma::NameShort( Alg[Snk], pszSep ) };
    for( int Src = 0; Src < NumOps; Src++ )
    {
      std::string sSrc{ Common::Gamma::NameShort( Alg[Src], pszSep ) };
      const std::string sOutBaseShort{ Name_.Base + sSnk + sSrc };
      const std::string sOutBase{ Prefix + sOutBaseShort };
      for(int f = 0; f < NumSummaries; f++)
      {
        std::string sOutFileName{MakeFilename(sOutBase, SummaryNames[f], Name_.Seed, TEXT_EXT)};
        std::ofstream s( sOutFileName );
        s << Comment << SummaryHeader[f] << NewLine << Comment << sOutBaseShort
          << NewLine << Comment << "Config " << Name_.Seed
          << NewLine << Comment << "Signal " << this->getSignal( Alg[Snk], Alg[Src] )
          << NewLine << FieldNames << ( ( f == 0 ) ? FieldNames2 : "" )
          << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
        const T * const pc = (*this)( Alg[Snk], Alg[Src] );
        for( int t = 0; t < nt; t++ )
        {
          if( f == 0 )
          {
            double re = pc[t].real();
            double im = pc[t].imag();
            s << t << sep << re << sep << 0 << sep << 0 << sep << 1
                   << sep << im << sep << 0 << sep << 0 << sep << 1 << std::endl;
          }
          else
          {
            double DThis;
            switch(f)
            {
              case 1: // mass
                DThis = std::log( abs( pc[t].real() / pc[(t + 1 + nt) % nt].real() ) );
                break;
              case 2: // cosh mass
                DThis = std::acosh((pc[(t - 1 + nt) % nt].real() + pc[(t + 1) % nt].real()) / (2 * pc[t].real()));
                break;
                case 3: // sinh mass
                DThis = std::asinh((pc[(t - 1 + nt) % nt].real() - pc[(t + 1) % nt].real()) / (2 * pc[t].real()));
                break;
              default:
                DThis = 0;
            }
            s << t << sep << DThis << sep << 0 << sep << 0 << sep << 1 << std::endl;
          }
        }
      }
    }
  }
}

// This is for a sample of correlators. There is a central replica, specified aux replicas and variable number of samples

template <typename T, int NumAuxSamples_ = 0>
class Sample
{
  int NumSamples_;
  int Nt_;
  std::unique_ptr<T[]> m_pData;
public:
  Signal Signal_ = Signal::Unknown;
  using scalar_type = typename SampleTraits<T>::scalar_type;
  static constexpr int scalar_size { sizeof( scalar_type ) };
  static constexpr int NumAuxSamples{ NumAuxSamples_ };
  static constexpr int NumExtraSamples{ NumAuxSamples + 1 }; // Auxiliary samples + central replica
  static constexpr int idxCentral{ -1 };
  static constexpr int idxAux{ idxCentral - NumAuxSamples };
  inline int NumSamples() const { return NumSamples_; }
  inline int Nt() const { return Nt_; }
  void resize( int NumSamples, int Nt )
  {
    if( NumSamples_ != NumSamples || Nt_ != Nt )
    {
      NumSamples_ = NumSamples;
      Nt_ = Nt;
      if( Nt == 0 )
        m_pData.reset( nullptr );
      else
        m_pData.reset( new T[ static_cast<std::size_t>( NumSamples + NumExtraSamples ) * Nt ] );
    }
  }
  bool IsFinite() { return Common::IsFinite( reinterpret_cast<scalar_type *>( m_pData.get() ),
      static_cast<size_t>( ( NumSamples_ + NumExtraSamples ) * ( SampleTraits<T>::is_complex ? 2 : 1 ) ) * Nt_ ); }
  Sample<T, NumAuxSamples_> Bootstrap( int NumBootSamples, SeedType Seed, Signal sig );
  void Read ( const std::string &FileName, std::string &GroupName );
  void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void Read ( const std::vector<std::string> & FileName, int Fold, int Shift, int NumOps );
  explicit Sample( int NumSamples = 0, int Nt = 0 ) : NumSamples_{0}, Nt_{0} { resize( NumSamples, Nt ); }
  Sample( const std::string &FileName, std::string &GroupName ) : NumSamples_{0}, Nt_{0}
    { Read( FileName, GroupName ); }
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
  inline void ZeroSlot( int idxSlot = idxCentral )
  {
    if( Nt_ )
    {
      T * p{ (*this)[idxSlot] };
      for( int t = 0; t < Nt_; t++ )
        *p++ = 0;
    }
  }
  void MakeMean( int idxSlot = idxCentral )
  {
    ZeroSlot( idxSlot );
    if( NumSamples_ )
    {
      T * const dst{ (*this)[idxSlot] };
      const T * src{ (*this)[0] };
      for( int i = 0; i < NumSamples_; i++ )
        for( int t = 0; t < Nt_; t++ )
          dst[t] += *src++;
      for( int t = 0; t < Nt_; t++ )
        dst[t] /= NumSamples_;
    }
  }
};

using SampleC = Sample<std::complex<double>>;
using SampleD = Sample<double>;

// Perform bootstrap
template <typename T, int NumAuxSamples_>
Sample<T, NumAuxSamples_> Sample<T, NumAuxSamples_>::Bootstrap( int NumBootSamples, SeedType Seed, Signal sig )
{
  using fint = std::uint_fast32_t;
  std::mt19937                        engine( Seed );
  std::uniform_int_distribution<fint> random( 0, NumSamples_ - 1 );
  Sample<T, NumAuxSamples_>           boot( NumBootSamples, Nt_ );
  boot.Signal_ = sig;

  // Compute the mean, then copy all the extra info to the bootstrap I will return
  MakeMean();
  std::copy( (*this)[idxAux], (*this)[0], boot[idxAux] );

  // Now make the bootstrap replicas
  T * dst{ boot[0] };
  for( int i = 0; i < NumBootSamples; i++, dst += Nt_ )
  {
    // Initialise this sum for this bootstrap sample to zero
    for( int t = 0; t < Nt_; t++)
      dst[t] = 0;
    for( int s = 0; s < NumSamples_; s++ )
    {
      const T * src{ (*this)[ random( engine ) ] };
      for( int t = 0; t < Nt_; t++ )
        dst[t] += *src++;
    }
    // Turn the sum into an average
    for( int t = 0; t < Nt_; t++)
      dst[t] /= NumSamples_;
  }
  return boot;
}

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
          if( !IsFinite() )
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
        Signal Sig{ Signal::Unknown };
        try
        {
          a = g.openAttribute("Signal");
          std::string s{};
          H5::StrType sType = a.getStrType();
          a.read( sType, s );
          a.close();
          if( EqualIgnoreCase( s, "real" ) )
            Sig = Signal::Real;
          else if( EqualIgnoreCase( s, "imag" ) )
            Sig = Signal::Imag;
        }
        catch(const H5::Exception &)
        {
          H5::Exception::clearErrorStack();
        }
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
            ds.read( (*this)[idxCentral], H5File::ComplexType( scalar_size ) );
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
                ds.read( (*this)[0], H5File::ComplexType( scalar_size ) );
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
                      ds.read( (*this)[idxAux], H5File::ComplexType( scalar_size ) );
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
    a = g.createAttribute( "nAux", H5::PredType::STD_U16LE, ds1 );
    tmp = NumAuxSamples;
    a.write( H5::PredType::NATIVE_INT, &tmp );
    a.close();
    a = g.createAttribute( "nSample", H5::PredType::STD_U32LE, ds1 );
    a.write( H5::PredType::NATIVE_INT, &NumSamples_ );
    a.close();
    if( Signal_ == Signal::Real || Signal_ == Signal::Imag )
    {
      std::stringstream ss;
      ss << Signal_;
      const std::string s{ ss.str() };
      H5::StrType sType( H5::PredType::C_S1, s.length() );
      a = g.createAttribute( "Signal", sType, ds1 );
      a.write( sType, s );
      a.close();
    }
    ds1.close();
    Dims[0] = Nt_;
    ds1 = H5::DataSpace( 1, Dims );
    H5::DataSet ds = g.createDataSet( "data_C", H5File::ComplexType( scalar_size ), ds1 );
    ds.write( (*this)[idxCentral], H5File::ComplexType( scalar_size ) );
    ds.close();
    ds1.close();
    Dims[0] = NumSamples_;
    Dims[1] = Nt_;
    ds1 = H5::DataSpace( 2, Dims );
    ds = g.createDataSet( "data_S", H5File::ComplexType( scalar_size ), ds1 );
    ds.write( (*this)[0], H5File::ComplexType( scalar_size ) );
    ds.close();
    ds1.close();
    if( NumAuxSamples )
    {
      Dims[0] = NumAuxSamples;
      ds1 = H5::DataSpace( 2, Dims );
      ds = g.createDataSet( "data_Aux", H5File::ComplexType( scalar_size ), ds1 );
      ds.write( (*this)[idxAux], H5File::ComplexType( scalar_size ) );
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

template <typename T, int NumAuxSamples_>
void Sample<T, NumAuxSamples_>::Read( const std::vector<std::string> & FileName, int Fold, int Shift, int NumOps )
{
  if( abs( Fold ) > 2 )
    throw std::runtime_error( "Error: Fold=" + std::to_string( Fold ) + " invalid" );
  const bool bAlternateFold{ abs( Fold ) > 1 };
  if( bAlternateFold )
  {
    if( Fold > 0 )
      --Fold;
    else
      ++Fold;
  }
  const int FirstFold{ Fold };
  const int NumFiles{ static_cast<int>( FileName.size() ) };
  // These will be initialised when we read in the first correlator
  int NtSource = -999; // Number of timeslices in the source file(s)
  int NtDest   = -888; // Number of timeslices after folding
  int NSamples = -777; // Number of bootstrap samples per correlator
  std::vector<double> ReadBuffer;
  SampleC f;
  for( int i = 0; i < NumFiles; ++i )
  {
    std::string GroupName;
    f.Read( FileName[i], GroupName );
    std::cout << FileName[i] << " (" << GroupName << ")\n";
    if( i == 0 )
    {
      // Initialise defaults first time around
      NSamples = f.NumSamples();
      NtSource = f.Nt();
      NtDest   = Fold ? NtSource / 2 + 1 : NtSource;
      std::cout << "\tNt=" << NtSource << " x " << NSamples << " samples";
      if( Fold )
        std::cout << " folded into " << NtDest << " timeslices with "
                  << (bAlternateFold ? "alternating, " : "")
                  << (Fold == 1 ? "positive" : "negative") << " parity";
      std::cout << "\n";
      resize( NSamples, NtDest * NumFiles );
    }
    else
    {
      if( NSamples != f.NumSamples() )
        throw std::runtime_error( "Error: inconsistent number of samples" );
      if( NtSource != f.Nt() )
        throw std::runtime_error( "Error: inconsistent number of timeslices" );
      if( bAlternateFold )
      {
        Fold = FirstFold;
        const int row{ i / NumOps };
        const int col{ i % NumOps };
        if( row & 1 )
          Fold *= -1;
        if( col & 1 )
          Fold *= -1;
      }
    }
    // Now copy the data from this bootstrap into the combined bootstrap
    T * pDest = (*this)[idxAux] + NtDest * i;
    for( int i = idxAux; i < NSamples; i++, pDest += NtDest * ( NumFiles - 1 ) )
    {
      T * ReadBuffer = (*this)[i];
      for( int k = 0; k < NtDest; ++k )
      {
        T Source1 = ReadBuffer[ ( NtSource + Shift + k ) % NtSource ];
        if( Fold != 0 )
        {
          T Source2 = ReadBuffer[ ( 2 * NtSource + Shift - k ) % NtSource ];
          if( Fold > 0 )
            Source1 = ( Source1 + Source2 ) * 0.5;
          else
            Source1 = ( std::abs( Source2 ) + std::abs( Source1 ) ) * 0.5;
            //m( Dest, 0 ) = ( Source2 - Source1 ) * 0.5;
        }
        *pDest++ = Source1;
      }
    }
  }
}

// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      std::string &GroupName,
                      const std::string &ObjectName = std::string( "correlator" ) );

void ReadRealArray(std::vector<double> &buffer, const std::string &FileName,
                   std::string &GroupName,
                   const std::string &ObjectName = std::string( "correlator" ) );

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
  
  std::string              Name;
  std::vector<std::string> Args;
  SwitchMap                Switches;
  
private:
  static void SkipPastSep( const char * & p );
  
public:
  void Parse( int argc, const char *argv[], const std::vector<SwitchDef> &defs );

  CommandLine() = default;
  CommandLine( int argc, const char *argv[], const std::vector<SwitchDef> &defs )
  { Parse( argc, argv, defs ); }

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
  
  inline const std::vector<std::string> & SwitchStrings( const std::string &Switch )
  {
    SwitchMap::const_iterator it = Switches.find( Switch );
    if( it == Switches.end() )
    {
      static const std::vector<std::string> v;
      return v;
    }
    const std::vector<std::string> &v{ it->second };
    return v;
  }

  template<typename T> inline T SwitchValue( const std::string &Switch, int Subscript = 0 )
  {
    const std::vector<std::string> &v{ SwitchStrings( Switch ) };
    if( static_cast<std::size_t>( Subscript ) >= v.size() )
      throw std::invalid_argument( "Switch " + Switch + "[" + std::to_string( Subscript ) + "] not found" );
    return FromString<T>( v[Subscript] );
  }
};

std::ostream& operator<<( std::ostream& os, const CommandLine &cl);

// Make summary files of a bootstrap of a correlator
void SummariseBootstrapCorr(const SampleC &out, const std::string & sOutFileBase, SeedType Seed );//, int momentum_squared);

END_COMMON_NAMESPACE
#endif // Common_hpp
