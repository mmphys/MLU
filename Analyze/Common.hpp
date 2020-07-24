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

// std c++
#include <array>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <regex>
#include <string>
#include <vector>

// posix
#include <glob.h>
#include <unistd.h>

// HDF5 Library
#include <H5Cpp.h>
#include <H5CompType.h>
#include <H5public.h>

// Eigen dense matrices
#include <Eigen/Dense>

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

// Check that a stream is empty ... or only contains white space to end of stream
template <class _CharT, class _Traits>
inline bool StreamEmpty( std::basic_istream<_CharT, _Traits> & s )
{
  return s.eof() || ( s >> std::ws && s.eof() );
}

// Text required for summaries of correlators
namespace CorrSumm {
  extern const char sep[];
  extern const char Comment[];
  static constexpr int NumFields{ 3 };
  extern const char * FieldNames[NumFields];
};

extern const std::string Space;
extern const std::string Underscore;
extern const std::string Period;
extern const std::string NewLine;
extern const std::string Comma;

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
  extern const std::array<const std::string, nGamma> name;      // Long name, per Grid
  extern const std::array<const std::string, nGamma> nameShort; // My abbreviations
  std::string NameShort( Algebra alg, const char * pszPrefix = nullptr );
  std::ostream& operator<<(std::ostream& os, const Algebra &a);
  std::istream& operator>>(std::istream& is, Algebra &a);
};

extern const std::string sBootstrap;
extern const std::string sFold;
extern const std::string sModel;
extern const std::string sParams;
extern const std::string sCormat;
extern const std::string sNtUnfolded;
extern const std::string st0Negated;
extern const std::string sConjugated;
extern const std::string sTI;
extern const std::string sTF;
extern const std::string sDoF;
extern const std::string sNumExponents;
extern const std::string sNumFiles;
extern const std::string sFactorised;
extern const std::string sCovarFrozen;
extern const std::string sOperators;
extern const std::string sAuxNames;
extern const std::string sSummaryNames;
extern const std::string sColumnNames;
extern const std::string sSeed;
extern const std::string sSeedMachine;
extern const std::string sRandom;
extern const std::string sSampleSize;
extern const std::string sConfigCount;
extern const std::string sFileList;
extern const std::string sBootstrapList;
extern const std::string sBinSize;
extern const std::string sNE;

extern const std::vector<std::string> sCorrSummaryNames;

using SeedType = unsigned int;

// Generic conversion from a string to any type
template<typename T> inline T FromString( const std::string &String ) {
  T t;
  std::istringstream iss( String );
  if( !( iss >> t && StreamEmpty( iss ) ) )
    throw std::invalid_argument( "Argument \"" + String + "\" is not type " + typeid(T).name() );
  return t;
}

// Converting a string to a string makes a copy
template<> inline std::string FromString<std::string>( const std::string &String )
  { return String; }

// Generic conversion from a string to an array of any type (comma or space separated)
template<typename T>
inline std::vector<T> ArrayFromString( const std::string &String, std::vector<bool> *pNeg = nullptr )
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
    bool bSign{ pNeg && iss.peek() == '-' };
    if( bSign )
      assert( iss.get() == '-' );
    if( !( iss >> t ) )
      throw std::invalid_argument( "ArrayFromString: \"" + String + "\" is not type " + typeid(T).name() );
    v.push_back( t );
    if( pNeg )
      pNeg->push_back( bSign );
  }
  return v;
}

inline std::string AppendSlash( const std::string & String )
{
  std::string s{ String };
  std::size_t Len = s.length();
  if( Len && s[Len - 1] != '/' )
    s.append( 1, '/' );
  return s;
}

// Default delimeters for the next couple of functions
extern const char szDefaultDelimeters[];

// Remove anything past the last delimeter from string, returning the removed part in suffix
// Return success / fail
bool ExtractSuffix( std::string &String, std::string &Suffix, const char * pszDelimeters = nullptr );

// Split String into an array using specified delimeters
std::vector<std::string> Split( const std::string &String, const char * pszDelimeters = nullptr );

// Extract suffix, then split strings. Default delimeters '.' and '_' respectively
bool ExtractSuffixSplit( std::string &String, std::vector<std::string> &Suffii,
                        const char * pszStringDelim = nullptr, const char * pszSuffixDelim = nullptr );

// Dump the environment to stdout, prefixed by optional message
void DumpEnv(int argc, const char * const *argv, const char * pStr = nullptr );

// Does the specified file exist?
bool FileExists( const std::string& Filename );

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
  int iFlags = GLOB_BRACE | GLOB_TILDE | GLOB_NOSORT | GLOB_NOCHECK;// | GLOB_NOMAGIC
  for( Iter i = first; i != last; i++ )
  {
    NameBuffer.resize( PrefixLen );
    NameBuffer.append( *i );
    const int globResult{ glob( NameBuffer.c_str(), iFlags, NULL, &globBuf ) };
    if( globResult )
    {
      if( globResult == GLOB_NOMATCH )
        continue;
      globfree( &globBuf );
      throw std::runtime_error( "glob() returned error " + std::to_string( globResult ) );
    }
    iFlags |= GLOB_APPEND;
  }

  // collect all the filenames into a std::list<std::string>
  std::vector<std::string> Filenames;
  Filenames.reserve( globBuf.gl_pathc );
  for( size_t i = 0; i < globBuf.gl_pathc; ++i )
    Filenames.push_back( std::string( globBuf.gl_pathv[i] ) );
  globfree(&globBuf);
  return Filenames;
}

// Wrapper for posix gethostname()
std::string GetHostName();

extern const double NaN;

// test whether a type is complex
template<typename T> struct is_complex                  : public std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : public std::true_type {};

// Component-wise absolute value for complex types
template<typename T> typename std::enable_if<is_complex<T>::value, T>::type
ComponentAbs( T c ) { return { std::abs( c.real() ), std::abs( c.imag() ) }; }

// Component-wise absolute value for scalars
template<typename T> typename std::enable_if<!(is_complex<T>::value), T>::type
ComponentAbs( T r ) { return std::abs( r ); }

// IsFinite() for floats and complex types
template<typename T> inline typename std::enable_if<is_complex<T>::value, bool>::type
IsFinite( const T &c ) { return std::isfinite( c.real() ) && std::isfinite( c.imag() ); }
template<typename T> inline typename std::enable_if<std::is_floating_point<T>::value, bool>::type
IsFinite( const T &c ) { return std::isfinite( c ); }

// Are all the floating point numbers pointed to finite
template <typename T, typename I> inline bool IsFinite( const T * d, I n )
{
  while( n-- )
    if( !IsFinite( *d++ ) )
      return false;
  return true;
}

// Are all the floating point numbers in this vector finite
template <typename T> inline bool IsFinite( const std::vector<T> & v )
{
  for( const T &n : v )
    if( !IsFinite( n ) )
      return false;
  return true;
}

// Are all the floating point numbers in this matrix finite
template <typename T> inline bool IsFinite( const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & m, bool bDiagonalsOnly = false )
{
  for( Eigen::Index row = 0; row < m.rows(); ++row )
    for( Eigen::Index col = 0; col < m.cols(); ++col )
      if( ( !bDiagonalsOnly || row == col ) && !IsFinite( m( row, col ) ) )
        return false;
  return true;
}

template <typename T = double>
class ValWithEr
{
public:
  T Central;
  T Low;
  T High;
  T Check;
  void Get( T dCentral, std::vector<T> &Data, std::size_t Count );
  ValWithEr() = default;
  ValWithEr( T dCentral, std::vector<T> &Data, std::size_t Count )
  { Get( dCentral, Data, Count ); }
  ValWithEr( T Central, T Low, T High, T Check = 1 );
};

template <typename T> ValWithEr<T>::ValWithEr( T central_, T low_, T high_, T check_ )
: Central{central_}, Low{low_}, High{high_}, Check{check_} {}

template <typename T>
inline std::ostream & operator<<( std::ostream &os, const ValWithEr<T> &v )
{
  return os << v.Central << " " << v.Low << " " << v.High << " " << v.Check;
}

// Generic representation of momentum
struct Momentum
{
  int x;
  int y;
  int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
  Momentum() : Momentum(0,0,0) {}
  inline explicit operator bool() const { return x!=0 || y!=0 || z!=0; }
  inline Momentum operator-() const { return Momentum(-x, -y, -z); }
  inline Momentum abs() const { return Momentum( x<0?-x:x, y<0?-y:y, z<0?-z:z ); }
  inline bool operator==(const Momentum &m) const { return x==m.x || y==m.y || z==m.z; }
  inline bool IsNeg() const { return x<0 || ( x==0 && ( y<0 || ( y==0 && z < 0 ))); }
  inline bool EqualsNeg(const Momentum &m) const { return x==-m.x || y==-m.y || z==-m.z; }
  inline int p2() const { return x * x + y * y + z * z; }
  std::string p2_string  ( const std::string &separator ) const;
  std::string to_string  ( const std::string &separator, bool bNegative = false ) const;
  std::string to_string4d( const std::string &separator, bool bNegative = false ) const;
  // Strip out momentum info from string if present
  bool Extract( std::string &s, bool IgnoreSubsequentZeroNeg = true );
};

std::ostream& operator<<( std::ostream& os, const Momentum &p );
std::istream& operator>>( std::istream& is, Momentum &p );

// Attributes for filenames in form base.type.seed.ext
struct FileNameAtt
{
  std::string Filename; // Full (and unadulterated) original filename
  std::string NameNoExt;// Original filename, without path and without extension
  std::string Dir;      // Directory part of the filename (with trailing '/')
  std::string Base;     // Base of the filename
  std::string Type;
  std::string SeedString;
  bool        bSeedNum = false;
  SeedType    Seed = 0; // Numeric value of SeedString, but only if bSeedNum == true
  std::string Ext;
  std::vector<int> op;
  std::vector<std::string> Extra;
  // reset contents
  void clear()
  {
    Filename.clear();
    NameNoExt.clear();
    Dir.clear();
    Base.clear();
    Type.clear();
    SeedString.clear();
    bSeedNum = false;
    Seed = 0;
    Ext.clear();
    op.clear();
    Extra.clear();
  }

  void Parse( const std::string &Filename_, std::vector<std::string> * pOpNames = nullptr );
  FileNameAtt() = default;
  explicit FileNameAtt( const std::string &Filename, std::vector<std::string> * pOpNames = nullptr )
    { Parse( Filename, pOpNames ); }
  void ParseExtra( unsigned int MaxElements = UINT_MAX );
};

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext);

// If present, remove integer preceded by Token from a string
void ExtractInteger( std::string &Prefix, bool &bHasValue, int &Value, const std::string Token );

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &s, bool &bHasTimeslice, int & Timeslice );

// Strip out DeltaT from a string if present
void ExtractDeltaT( std::string &Prefix, bool &bHasDeltaT, int &DeltaT );

// Remove any gammas from Prefix
std::vector<Gamma::Algebra> ExtractGamma( std::string &Prefix );

struct ConfigCount {
  int         Config;
  int         Count = 0;
  ConfigCount( int Config_, int Count_ ) : Config{ Config_ }, Count{ Count_ } {}
  ConfigCount() = default;
  ConfigCount( const ConfigCount & ) = default;
  ConfigCount( ConfigCount && ) = default;
  ConfigCount& operator=( const ConfigCount & ) = default;
  ConfigCount& operator=( ConfigCount && ) = default;
};

inline std::ostream & operator<<( std::ostream &os, const ConfigCount &cc )
{
  return os << "{" << cc.Config << "," << cc.Count << "}";
}

// My implementation of H5File - adds a definition of complex type
namespace H5 {
  template <typename T> struct Equiv;
  template<> struct Equiv<float>       { static const ::H5::PredType& Type; };
  template<> struct Equiv<double>      { static const ::H5::PredType& Type; };
  template<> struct Equiv<long double> { static const ::H5::PredType& Type; };
  template<> struct Equiv<std::string> { static const ::H5::StrType Type; };
  template<> struct Equiv<char *>      { static const ::H5::StrType& Type; };
  template<> struct Equiv<std::uint_fast32_t>{ static const ::H5::PredType& Type; };
  template<> struct Equiv<std::complex<float>>      { static const ::H5::CompType Type; };
  template<> struct Equiv<std::complex<double>>     { static const ::H5::CompType Type; };
  template<> struct Equiv<std::complex<long double>>{ static const ::H5::CompType Type; };
  template<> struct Equiv<ValWithEr<float>>         { static const ::H5::CompType Type; };
  template<> struct Equiv<ValWithEr<double>>        { static const ::H5::CompType Type; };
  template<> struct Equiv<ValWithEr<long double>>   { static const ::H5::CompType Type; };
  template<> struct Equiv<ConfigCount> { static const ::H5::CompType Type; };

  /**
   Open the specified HDF5File
   Read from the specified group name (if pGroupName is non-null and *pGroupName not empty)
   Otherwise read from the first group in the file, returning the group name if pGroupName not null
   */
  void OpenFileGroup(::H5::H5File &f, ::H5::Group &g, const std::string &FileName, const char *PrintPrefix = nullptr,
                     std::string * pGroupName = nullptr, unsigned int flags = H5F_ACC_RDONLY);

  // Get first groupname from specified group
  std::string GetFirstGroupName( ::H5::Group & g );

  // Read the gamma algebra attribute string and make sure it's valid
  Gamma::Algebra ReadGammaAttribute( ::H5::Group &g, const char * pAttName );

  template<typename AttributeOrDataSet>
  std::vector<std::string> ReadStrings( const AttributeOrDataSet &a )
  {
    ::H5::DataSpace dsp = a.getSpace();
    const ::H5::StrType aType = a.getStrType();
    const int rank{ dsp.getSimpleExtentNdims() };
    if( rank != 1 )
      throw std::runtime_error( sOperators + " number of dimensions=" + std::to_string( rank ) + ", expecting 1" );
    hsize_t NumOps;
    dsp.getSimpleExtentDims( &NumOps );
    std::vector<std::string> vs( NumOps );
    char * * MDString = new char *[NumOps];
    for( std::size_t i = 0; i < NumOps; i++ )
      MDString[i] = nullptr;
    try
    {
      a.read( aType, (void *)MDString );
    }
    catch(...)
    {
      for( std::size_t i = 0; i < NumOps; i++ )
        if( MDString[i] )
          delete [] MDString[i];
      delete [] MDString;
      throw;
    }
    for( std::size_t i = 0; i < NumOps; i++ )
    {
      vs[i] = MDString[i];
      delete [] MDString[i];
    }
    delete [] MDString;
    return vs;
  }

  // Make a multi-dimensional string attribute
  void WriteAttribute( ::H5::Group &g, const std::string &AttName, const std::vector<std::string> &vs );
};

// Enumeration describing whether signal is in complex or real part
enum class Reality
{
  Imag = -1,
  Unknown = 0,
  Real
};
inline Reality operator*( const Reality l, const Reality r )
{
  if( ( l != Reality::Imag && l != Reality::Real ) || ( r != Reality::Imag && r != Reality::Real ) )
    return Reality::Unknown;
  if( l == r )
    return Reality::Real;
  return Reality::Imag;
}

extern const std::string sUnknown;
extern const std::string sReality;
extern const std::string Reality_TextEquiv_Real;
extern const std::string Reality_TextEquiv_Imaginary;

std::ostream& operator<<(std::ostream& os, const Reality &reality);
std::istream& operator>>(std::istream& is, Reality &reality);

// Enumeration describing whether signal is in complex or real part
enum class Parity
{
  Odd = -1,
  Unknown,
  Even
};
inline Parity operator*( const Parity l, const Parity r )
{
  if( ( l != Parity::Even && l != Parity::Odd ) || ( r != Parity::Even && r != Parity::Odd ) )
    return Parity::Unknown;
  if( l == r )
    return Parity::Even;
  return Parity::Odd;
}
extern const std::string sParity;
extern const std::string Parity_TextEquiv_Even;
extern const std::string Parity_TextEquiv_Odd;

std::ostream& operator<<(std::ostream& os, const Parity &parity);
std::istream& operator>>(std::istream& is, Parity &parity);

// Enumeration describing whether signal is in complex or real part
enum class Sign
{
  Negative = -1,
  Unknown = 0,
  Positive
};
inline Sign operator*( const Sign l, const Sign r )
{
  if( ( l != Sign::Positive && l != Sign::Negative ) || ( r != Sign::Positive && r != Sign::Negative ) )
    return Sign::Unknown;
  if( l == r )
    return Sign::Positive;
  return Sign::Negative;
}
extern const std::string sSign;
extern const std::string Sign_TextEquiv_Positive;
extern const std::string Sign_TextEquiv_Negative;

std::ostream& operator<<(std::ostream& os, const Sign &sign);
std::istream& operator>>(std::istream& is, Sign &sign);

struct RPS
{
  Reality reality = Reality::Unknown;
  Parity  parity = Parity::Unknown;
  Sign    sign = Sign::Unknown;
  RPS(){};
  RPS( Reality reality_, Parity parity_, Sign sign_ ) : reality{reality_}, parity{parity_}, sign{sign_} {}
};

inline RPS operator*( const RPS &l, const RPS &r )
{
  RPS rps;
  rps.reality = l.reality * r.reality;
  rps.parity = l.parity * r.parity;
  rps.sign = l.sign * r.sign;
  if( l.reality == Common::Reality::Imag && r.reality == Common::Reality::Imag )
    rps.sign = rps.sign * Common::Sign::Negative;
  return rps;
}
inline std::ostream& operator<<(std::ostream& os, const RPS &rps)
{
  os << rps.reality << ", " << rps.parity << ", " << rps.sign;
  return os;
}

struct FoldProp
{
  RPS rps;
  bool t0Abs = true;
  bool Conjugate = false;
  bool Parse( const char * const pc, std::size_t const Len )
  {
    bool bOK = true;
    bool bRealImagOnly = true;
    static constexpr int NumOptions{ 8 };
    static const char Options[] = "RIEOPN0C";
    int iCount[NumOptions];
    for( int i = 0; i < NumOptions; i++ )
      iCount[i] = 0;
    for( int i = 0; bOK && i < Len; i++ )
    {
      char c = pc[i];
      if( c >= 'a' && c <= 'z' )
        c -= 'a' - 'A';
      int idx = 0;
      while( idx < NumOptions && c != Options[idx] )
        idx++;
      if( idx < NumOptions )
        iCount[idx]++;
      else
        bOK = false;
    }
    if( bOK )
    {
      if( iCount[0] == 1 && iCount[1] == 0 )
        rps.reality = Common::Reality::Real;
      else if( iCount[0] == 0 && iCount[1] == 1 )
        rps.reality = Common::Reality::Imag;
      else if( iCount[0] != 0 || iCount[1] != 0 )
        bOK = false;
      if( iCount[2] == 1 && iCount[3] == 0 )
        rps.parity = Common::Parity::Even;
      else if( iCount[2] == 0 && iCount[3] == 1 )
        rps.parity = Common::Parity::Odd;
      else if( iCount[2] != 0 || iCount[3] != 0 )
        bOK = false;
      if( iCount[4] == 1 && iCount[5] == 0 )
        rps.sign = Common::Sign::Positive;
      else if( iCount[4] == 0 && iCount[5] == 1 )
        rps.sign = Common::Sign::Negative;
      else if( iCount[4] != 0 || iCount[5] != 0 )
        bOK = false;
      if( iCount[6] == 1 )
        t0Abs = false;
      else if( iCount[6] != 0 )
        bOK = false;
      if( iCount[7] == 1 )
        Conjugate = true;
      else if( iCount[7] != 0 )
        bOK = false;
    }
    if( !bOK )
    {
      const std::string s{ pc, Len };
      throw std::runtime_error( "Option string '" + s + "' invalid" );
    }
    for( int i = 2; bRealImagOnly && i < NumOptions; ++i )
      if( iCount[i] )
        bRealImagOnly = false;
    return bRealImagOnly;
  }
  FoldProp() = default;
  explicit FoldProp( const char * pc, std::size_t Len ) { Parse( pc, Len ); }
};

inline std::ostream & operator<<( std::ostream &os, const FoldProp &f )
{
  return os << f.rps << (f.t0Abs ? ", abs( C(0) )" : "" )
            << ( f.Conjugate ? " and conjugate operators" : "" );
}

// Traits for samples
// SampleTraits<ST>::value is true for supported types (floating point and complex)
template<typename ST, typename V = void> struct SampleTraits : public std::false_type{};
template<typename ST> struct SampleTraits<ST, typename std::enable_if<std::is_floating_point<ST>::value>::type> : public std::true_type
{
  using scalar_type = ST;
  static constexpr bool is_complex = false;
  static constexpr int scalar_count = 1;
  static inline scalar_type * ScalarPtr( ST * p ) { return p; };
  static inline const scalar_type * ScalarPtr( const ST * p ) { return p; };
};
template<typename ST> struct SampleTraits<std::complex<ST>, typename std::enable_if<SampleTraits<ST>::value>::type> : public std::true_type
{
  using scalar_type = typename SampleTraits<ST>::scalar_type;
  static constexpr bool is_complex = true;
  static constexpr int scalar_count = 2;
  static inline scalar_type * ScalarPtr( std::complex<ST> * p )
      { return reinterpret_cast<scalar_type *>( p ); };
  static inline const scalar_type * ScalarPtr( const std::complex<ST> * p )
      { return reinterpret_cast<const scalar_type *>( p ); };
};

/*
template<typename TSrc, typename TDst>
inline typename std::enable_if<!SampleTraits<TSrc>::is_complex || SampleTraits<TDst>::is_complex>::type
CopyBuffer( const TSrc * Begin, const TSrc * End, TDst * Dst, Signal s = Signal::Unknown )
{
  while( Begin != End )
    *Dst++ = *Begin++;
}

template<typename TSrc, typename TDst>
inline typename std::enable_if<SampleTraits<TSrc>::is_complex && !SampleTraits<TDst>::is_complex>::type
CopyBuffer( const TSrc * Begin, const TSrc * End, TDst * Dst, Signal s = Signal::Unknown )
{
  while( Begin != End )
  {
    *Dst++ = s == Signal::Imag ? Begin->Imag() : Begin->Real();
    Begin++;
  }
}
*/

template <typename T>
void SummaryHeader(std::ostream &s, const std::string & sOutFileName,
                   const char * pszHeaderComments = nullptr)
{
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  if( !s )
    throw std::runtime_error( "Unable to create " + sOutFileName );
  s << "# Summary: " << sOutFileName << NewLine
    << std::boolalpha << std::setprecision(std::numeric_limits<scalar_type>::digits10+2);
  const char * psz = nullptr;
  if( Traits::is_complex )
  {
    switch( sizeof( scalar_type ) )
    {
      case sizeof( float ):
        psz = "std::complex<float>";
        break;
      case sizeof( double ):
        psz = "std::complex<double>";
        break;
      case sizeof( long double ):
        psz = "std::complex<long double>";
        break;
    }
  }
  else
  {
    switch( sizeof( scalar_type ) )
    {
      case sizeof( float ):
        psz = "float";
        break;
      case sizeof( double ):
        psz = "double";
        break;
      case sizeof( long double ):
        psz = "long double";
        break;
    }
  }
  if( psz ) s << "# Representation: " << psz << NewLine;
  if( pszHeaderComments ) s << pszHeaderComments;
}

// Make a summary of the data
// Assume there is a central replica followed by nSample copies (i.e. nSample + 1 total)

template <typename T>
void SummaryHelper( const std::string & sOutFileName, const T * pData, const int nt,
                    const int nSample = 0, const char * pszHeaderComments = nullptr,
                    bool bFolded = false )
{
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  using namespace CorrSumm;
  std::ofstream s( sOutFileName );
  SummaryHeader<T>( s, sOutFileName, pszHeaderComments );
  // Now write all the field names
  s << "t";
  for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
  {
    for(int f = 0; f < NumFields; f++) // each field
    {
      std::string n{ FieldNames[f] };
      if( i )
        n.append( "_im" );
      s << sep << n << sep << n << "_low " << n << "_high " << n << "_check";
    }
  }
  s << NewLine;
  // Now perform summaries
  const int tMid{ bFolded ? nt : nt / 2 }; // Summaries may want to change sign for -ve wave
  scalar_type Central = -777; // Just to silence compiler warning
  std::vector<scalar_type> Data( nSample );
  for(int t = 0; t < nt; t++ )
  {
    // Index of previous and next timeslice relative to current timeslice
    const int Next{ ( t == nt - 1 ? 1 - nt :  1 ) * Traits::scalar_count };
    const int Prev{ ( t == 0      ? nt - 1 : -1 ) * Traits::scalar_count };
    s << t;
    for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
    {
      for(int f = 0; f < NumFields; f++) // each field
      {
        const scalar_type * p = Traits::ScalarPtr( pData ) + t * Traits::scalar_count + i;
        int Count = 0;
        for(int n = -1; n < nSample; n++, p += nt * Traits::scalar_count )
        {
          scalar_type d;
          switch( f )
          {
            case 0:
              d = * p;
              break;
            case 1: // mass
              if( t == 0 )
                d = NaN;
              else if( t <= tMid )
                d = std::log( p[Prev] / *p );
              else
                d = -std::log( *p / p[Next] );
              break;
            case 2: // cosh mass
              d = std::acosh( ( p[Prev] + p[Next] ) / ( *p * 2 ) );
              break;
            default:
              d = 0;
          }
          if( n == -1 )
            Central = d;
          else if( std::isfinite( d ) )
            Data[Count++] = d;
        }
        /*if( i==0 && f==1 && t==20 )
        {
          std::ofstream o( "BeforeSort.txt" );
          o << -1 << sep << Central << NewLine;
          for( int i = 0; i < Count; i++ )
            o << i << sep << Data[i] << NewLine;
        }*/
        Common::ValWithEr<scalar_type> v( Central, Data, Count );
        s << sep << v;
        /*if( i==0 && f==1 && t==20 )
        {
          std::ofstream o( "AfterSort.txt" );
          o << -1 << sep << Central << NewLine;
          for( int i = 0; i < Count; i++ )
            o << i << sep << Data[i] << NewLine;
        }*/
      }
    }
    s << NewLine;
  }
}

// Correlator file. Could be either single correlator, or multiple gammas

template <typename T>
class CorrelatorFile
{
  // Data members
public:
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  static constexpr bool is_complex { Traits::is_complex };
private:
  int NumSnk_ = 0;
  int NumSrc_ = 0;
  int Nt_ = 0;
  std::unique_ptr<T[]> m_pData;
  std::vector<Gamma::Algebra> AlgSnk_;
  std::vector<Gamma::Algebra> AlgSrc_;
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
    if( Sample < 0 || Sample >= NumSnk_ * NumSrc_ )
      throw std::out_of_range( "Sample " + std::to_string( Sample ) );
  }
  inline int SinkIndex( Gamma::Algebra g ) const
  {
    int idx;
    for( idx = 0; idx < AlgSnk_.size() && AlgSnk_[idx] != g; idx++ )
      ;
    return idx;
  }
  inline int SourceIndex( Gamma::Algebra g ) const
  {
    int idx;
    for( idx = 0; idx < AlgSrc_.size() && AlgSrc_[idx] != g; idx++ )
      ;
    return idx;
  }
public:
  inline int NumSnk() const { return NumSnk_; }
  inline int NumSrc() const { return NumSrc_; }
  inline int NumOps() const { return NumSnk() * NumSrc(); }
  inline int Nt() const { return Nt_; }
  inline int Timeslice() const { return bHasTimeslice ? Timeslice_ : 0; }
  inline const std::vector<Gamma::Algebra> &AlgSnk() const { return AlgSnk_; }
  inline const std::vector<Gamma::Algebra> &AlgSrc() const { return AlgSrc_; }
  void resize( int NumSnk, int NumSrc, int Nt )
  {
    const std::size_t OldMemSize{ static_cast<std::size_t>( NumSnk_ ) * NumSrc_ * Nt_ };
    const std::size_t NewMemSize{ static_cast<std::size_t>( NumSnk  ) * NumSrc  * Nt  };
    if( NumSnk_ != NumSnk || NumSrc_ != NumSrc )
    {
      AlgSnk_.clear();
      AlgSrc_.clear();
      NumSnk_ = NumSnk;
      NumSrc_ = NumSrc;
    }
    Nt_ = Nt;
    if( OldMemSize != NewMemSize )
    {
      if( NewMemSize == 0 )
        m_pData.reset( nullptr );
      else
        m_pData.reset( new T[ NewMemSize ] );
    }
  }
  void Read (const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk, std::vector<Gamma::Algebra> &AlgSrc,
             const int * pTimeslice = nullptr, const char * PrintPrefix = nullptr,
             std::string *pGroupName = nullptr, const std::vector<bool> *pAlgSnkNeg = nullptr);
  //void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void WriteSummary(const std::string &Prefix, const std::vector<Gamma::Algebra> &AlgSnk,
                    const std::vector<Gamma::Algebra> &AlgSrc);
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
  { return (*this)[ SinkIndex( gSink ) * NumSrc_ + SourceIndex( gSource ) ]; }
  const T * operator()( Gamma::Algebra gSink, Gamma::Algebra gSource ) const
  { return (*this)[ SinkIndex( gSink ) * NumSrc_ + SourceIndex( gSource ) ]; }
  bool IsFinite() { return Common::IsFinite( reinterpret_cast<scalar_type *>( m_pData.get() ),
      static_cast<size_t>( NumSnk_ * NumSrc_ * ( SampleTraits<T>::is_complex ? 2 : 1 ) ) * Nt_ ); }
  // Constructors (copy operations missing for now - add them if they become needed)
  CorrelatorFile() {}
  CorrelatorFile( CorrelatorFile && ) = default; // Move constructor
  CorrelatorFile(const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk, std::vector<Gamma::Algebra> &AlgSrc,
                 const int * pTimeslice = nullptr, const char * PrintPrefix = nullptr, std::string *pGroupName = nullptr)
  { Read( FileName, AlgSnk, AlgSrc, pTimeslice, PrintPrefix, pGroupName ); }
  // Operators
  CorrelatorFile& operator=(CorrelatorFile && r) = default; // Move assignment
};

using CorrelatorFileC = CorrelatorFile<std::complex<double>>;
using CorrelatorFileD = CorrelatorFile<double>;

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void CorrelatorFile<T>::Read(const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk,
                             std::vector<Gamma::Algebra> &AlgSrc, const int * pTimeslice,
                             const char * PrintPrefix, std::string *pGroupName,
                             const std::vector<bool> *pAlgSnkNeg)
{
  const bool bSameAlgebras{&AlgSnk == &AlgSrc};
  if( pAlgSnkNeg )
  {
    assert( pAlgSnkNeg->size() == AlgSnk.size() );
    assert( pAlgSnkNeg->size() && "Algebras unknown. Negating makes no sense" );
    assert( !bSameAlgebras && "Work out whether two negatives make a positive, then implement" );
  }
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
  ::H5::H5File f;
  ::H5::Group  g;
  std::string localGroupName;
  if( !pGroupName )
    pGroupName = &localGroupName;
  H5::OpenFileGroup( f, g, FileName, PrintPrefix, pGroupName );
  bool bOK = false;
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  try // to load a single correlator
  {
    if( ( AlgSnk.empty() || ( AlgSnk.size() == 1 && AlgSnk[0] == Gamma::Algebra::Unknown ) )
      && ( AlgSrc.empty() || ( AlgSrc.size() == 1 && AlgSrc[0] == Gamma::Algebra::Unknown ) ) )
    {
      ::H5::DataSet ds = g.openDataSet( "correlator" );
      ::H5::DataSpace dsp = ds.getSpace();
      int nDims{ dsp.getSimpleExtentNdims() };
      if( nDims == 1 )
      {
        hsize_t Dim[1];
        dsp.getSimpleExtentDims( Dim );
        if( Dim[0] <= std::numeric_limits<int>::max() )
        {
          resize( 1, 1, static_cast<int>( Dim[0] ) );
          AlgSnk_.resize( 1 );
          AlgSnk_[0] = Gamma::Algebra::Unknown;
          AlgSrc_.resize( 1 );
          AlgSrc_[0] = Gamma::Algebra::Unknown;
          if( AlgSnk.empty() )
            AlgSnk.push_back( Gamma::Algebra::Unknown );
          if( !bSameAlgebras && AlgSrc.empty() )
            AlgSrc.push_back( Gamma::Algebra::Unknown );
          T * const pData{ (*this)[0] };
          ds.read( pData, H5::Equiv<T>::Type );
          // Negate this correlator if requested
          if( pAlgSnkNeg && (*pAlgSnkNeg)[0] )
            for( int i = 0; i < Nt_; i++ )
              pData[i] = -pData[i];
          bOK = true;
        }
      }
    }
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
    ::H5::Exception::clearErrorStack();
  }
  if( !bOK )
  {
    try // to load from array of correlators indexed by gamma matrix
    {
      unsigned short NumVec;
      ::H5::Attribute a;
      a = g.openAttribute("_Grid_vector_size");
      a.read( ::H5::PredType::NATIVE_USHORT, &NumVec );
      a.close();
      // Must be a perfect square and have at least as many as entries as requested
      const unsigned short NumFileOps{static_cast<unsigned short>( std::sqrt( NumVec ) + 0.5 )};
      if( NumFileOps * NumFileOps == NumVec && NumFileOps >= AlgSnk.size() && NumFileOps >= AlgSrc.size() )
      {
        bOK = true;
        std::vector<int> count;
        for( unsigned short i = 0; bOK && i < NumVec; i++ ) // Loop through all vectors in file
        {
          bOK = false;
          ::H5::Group gi = g.openGroup( *pGroupName + "_" + std::to_string( i ) );
          ::H5::DataSet ds = gi.openDataSet( "corr" );
          ::H5::DataSpace dsp = ds.getSpace();
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
                resize(AlgSnk.size() ? static_cast<int>( AlgSnk.size() ) : NumFileOps,
                       AlgSrc.size() ? static_cast<int>( AlgSrc.size() ) : NumFileOps,
                       ThisNt);
                count.resize( NumSnk_ * NumSrc_, 0 ); // I want to check each operator combination appears once only
                if( AlgSnk.size() )
                {
                  AlgSnk_.resize( AlgSnk.size() );
                  std::copy( AlgSnk.cbegin(), AlgSnk.cend(), AlgSnk_.begin() );
                }
                else
                {
                  AlgSnk_.clear();
                  AlgSnk_.reserve( NumFileOps );
                }
                if( AlgSrc.size() )
                {
                  AlgSrc_.resize( AlgSrc.size() );
                  std::copy( AlgSrc.cbegin(), AlgSrc.cend(), AlgSrc_.begin() );
                }
                else
                {
                  AlgSrc_.clear();
                  AlgSrc_.reserve( NumFileOps );
                }
              }
              else if( ThisNt != Nt_ )
              {
                break;
              }
              // Read the gamma algebra strings and make sure they are valid
              const Gamma::Algebra gSnk{ H5::ReadGammaAttribute( gi, "gamma_snk" ) };
              int idxSnk;
              for( idxSnk = 0; idxSnk < AlgSnk_.size() && AlgSnk_[idxSnk] != gSnk; idxSnk++ )
                ;
              if( idxSnk == AlgSnk_.size() && AlgSnk_.size() < NumSnk_ )
                AlgSnk_.push_back( gSnk );
              bOK = true; // We can safely ignore gamma structures we're not interested in
              if( idxSnk < AlgSnk_.size() )
              {
                const Gamma::Algebra gSrc{ H5::ReadGammaAttribute( gi, "gamma_src" ) };
                int idxSrc;
                for( idxSrc = 0; idxSrc < AlgSrc_.size() && AlgSrc_[idxSrc] != gSrc; idxSrc++ )
                  ;
                if( idxSrc == AlgSrc_.size() && AlgSrc_.size() < NumSrc_ )
                  AlgSrc_.push_back( gSrc );
                if( idxSrc < AlgSrc_.size() )
                {
                  const int idx{ idxSnk * NumSrc_ + idxSrc };
                  T * const pData{ (*this)[idx] };
                  ds.read( pData, H5::Equiv<T>::Type );
                  count[idx]++;
                  if( pAlgSnkNeg && (*pAlgSnkNeg)[idxSnk] )
                    for( int i = 0; i < Nt_; i++ )
                      pData[i] = -pData[i];
                }
              }
            }
          }
        }
        // Make sure that everything we wanted was loaded once and only once
        for( int i = 0; bOK && i < NumSrc_ * NumSnk_; i++ )
          if( count[i] != 1 )
            bOK = false;
      }
    }
    catch(const ::H5::Exception &)
    {
      bOK = false;
      ::H5::Exception::clearErrorStack();
    }
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  if( !bOK )
    throw std::runtime_error( "Unable to read sample from " + FileName );
  if( !IsFinite() )
    throw std::runtime_error( "Values read are not all finite" );
  // If I'm discovering which operators are in the file, copy them back to caller
  // Bear in mind that caller may have passed in the same array for each gamma algebra
  bool bCopyBackSrc{ AlgSrc.empty() };
  bool bCopyBackSnk{ !bSameAlgebras && AlgSnk.empty() };
  if( bCopyBackSrc )
  {
    AlgSrc.resize( AlgSrc_.size() );
    std::copy( AlgSrc_.cbegin(), AlgSrc_.cend(), AlgSrc.begin() );
  }
  if( bCopyBackSnk )
  {
    AlgSnk.resize( AlgSnk_.size() );
    std::copy( AlgSnk_.cbegin(), AlgSnk_.cend(), AlgSnk.begin() );
  }
}

template <typename T>
void CorrelatorFile<T>::WriteSummary( const std::string &Prefix, const std::vector<Gamma::Algebra> &AlgSnk, const std::vector<Gamma::Algebra> &AlgSrc )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  const int nt{ Nt() };
  const std::vector<Gamma::Algebra> &MySnk{ AlgSnk.size() ? AlgSnk : AlgSnk_ };
  const std::vector<Gamma::Algebra> &MySrc{ AlgSrc.size() ? AlgSrc : AlgSrc_ };
  const int NumSnk{ static_cast<int>( MySnk.size() ) };
  const int NumSrc{ static_cast<int>( MySrc.size() ) };
  std::string sOutFileName{ Prefix };
  sOutFileName.append( Name_.Base );
  std::size_t Len{ sOutFileName.length() };
  std::string sSuffix( 1, '.' );
  sSuffix.append( Name_.SeedString );
  sSuffix.append( 1, '.' );
  sSuffix.append( TEXT_EXT );
  for( int Snk = 0; Snk < NumSnk; Snk++ )
  {
    static const char pszSep[] = "_";
    sOutFileName.resize( Len );
    sOutFileName.append( Common::Gamma::NameShort( AlgSnk[Snk], pszSep ) );
    std::size_t Len2{ sOutFileName.length() };
    for( int Src = 0; Src < NumSrc; Src++ )
    {
      sOutFileName.resize( Len2 );
      sOutFileName.append( Common::Gamma::NameShort( AlgSrc[Src], pszSep ) );
      sOutFileName.append( sSuffix );
      SummaryHelper( sOutFileName, (*this)( AlgSnk[Snk], AlgSrc[Src] ), nt );
    }
  }
}

/**
 This is for a sample of anything, but usually a sample of correlators.
 There is always a central replica (specified by user, not calculated as average) in addition to NumSamples copies
 If AuxNames are specified, then there are 3 * AuxNames.count() auxiliarry records (value, low and high)
 Fields are in memory in reverse order, so that index -1 is always the central value
*/

template <typename T>
class Sample
{
public:
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  using fint = std::uint_fast32_t;
  static constexpr bool is_complex { Traits::is_complex };
  static constexpr int scalar_count { Traits::scalar_count };
  static constexpr int idxCentral{ -1 };
protected:
  int NumSamples_;
  int Nt_;
  std::unique_ptr<T[]> m_pData;
  std::unique_ptr<fint[]> m_pRandNum; // Random numbers used to generate the bootstrap sample
  std::vector<std::string> AuxNames;
  int NumExtraSamples; // Same length as AuxNames, but minimum of 1 for central replica
  std::vector<std::string> SummaryNames;
  std::unique_ptr<ValWithEr<scalar_type>[]> m_pSummaryData;
  std::vector<std::string> ColumnNames;
  inline void AllocSummaryBuffer( std::size_t NewSize )
  {
    if( is_complex )
      NewSize *= 2;
    const std::size_t OldSize{ SummaryNames.size() * Nt_ * ( is_complex ? 2 : 1 ) };
    if( NewSize != OldSize )
    {
      if( NewSize )
        m_pSummaryData.reset( new ValWithEr<scalar_type>[ NewSize ] );
      else
        m_pSummaryData.reset( nullptr );
    }
  }
public:
  FileNameAtt Name_;
  SeedType Seed_ = 0; // This is only ever set by the bootstrap code
  std::string SeedMachine_; // name of the machine that ran the bootstrap
  int binSize = 1;
  int SampleSize = 0; // Number of samples (after binning) used to create bootstrap
  std::vector<Common::ConfigCount> ConfigCount; // Info on every config in the bootstrap in order
  std::vector<std::string> FileList; // Info on every config in the bootstrap in order
  inline int NumSamples() const { return NumSamples_; }
  inline int Nt() const { return Nt_; }
  inline const fint * RandNum() const { return m_pRandNum.get(); };
  inline ValWithEr<scalar_type> * getSummaryData( int idx = 0 )
  {
    const int iSize{ static_cast<int>( SummaryNames.size() ) };
    if( Nt_ == 0 || idx < 0 || idx >= iSize )
      throw "Summary " + std::to_string(idx) + "/" + std::to_string(iSize) + " doesn't exist";
    return m_pSummaryData.get() + idx * Nt_;
  }
  inline const ValWithEr<scalar_type> * getSummaryData( int idx = 0 ) const
  {
    const int iSize{ static_cast<int>( SummaryNames.size() ) };
    if( Nt_ == 0 || idx < 0 || idx >= iSize )
      throw "Summary " + std::to_string(idx) + "/" + std::to_string(iSize) + " doesn't exist";
    return m_pSummaryData.get() + idx * Nt_;
  }
  void WriteSummaryData( std::ostream &s, int idx = 0 ) const
  {
    s << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
    const ValWithEr<scalar_type> * p{ getSummaryData( idx ) };
    for( int t = 0; t < Nt_; t++ )
      s << ( t == 0 ? "" : " " ) << *p++;
  }
  void SetSummaryNames( const std::vector<std::string> &summaryNames_ )
  {
    AllocSummaryBuffer( summaryNames_.size() * Nt_ );
    SummaryNames = summaryNames_;
    if( is_complex )
      for( const std::string & s : summaryNames_ )
        SummaryNames.push_back( s + "_im" );
  }
  void SetColumnNames( const std::vector<std::string> &columnNames_ )
  {
    if( columnNames_.size() != Nt_ )
      throw std::runtime_error( "There should be names for " + std::to_string( Nt_ ) + " columns" );
    ColumnNames = columnNames_;
  }
  void WriteColumnNames( std::ostream &s ) const
  {
    const std::string Sep{ " " };
    if( Nt_ == 0 )
      throw std::runtime_error( "Can't write header - Nt_ = 0" );
    std::string Buffer{ "t" };
    const std::size_t BufferLen{ Buffer.size() };
    const std::string * ps{ &Buffer };
    for( std::size_t t = 0; t < Nt_; t++ )
    {
      if( ColumnNames.empty() )
      {
        Buffer.resize( BufferLen );
        Buffer.append( std::to_string( t ) );
      }
      else
        ps = &ColumnNames[t];
      s << ( t == 0 ? "" : Sep ) << *ps << Sep << *ps << "_low" << Sep << *ps << "_high" << Sep << *ps << "_check";
    }
  }
  const std::vector<std::string> & GetColumnNames() const { return ColumnNames; }
  int GetColumnIndex( const std::string & ColumnName ) const
  {
    for( int i = 0; i < ColumnNames.size(); i++ )
      if( EqualIgnoreCase( ColumnNames[i], ColumnName ) )
        return i;
    throw std::runtime_error( "Column " + ColumnName + " not found" );
  }
  void resize( int NumSamples, int Nt, std::vector<std::string> * pAuxNames = nullptr )
  {
    const int NewNumExtraSamples{ static_cast<int>( pAuxNames && pAuxNames->size() ? pAuxNames->size() : 1 ) };
    if( Nt_ != Nt )
    {
      AllocSummaryBuffer( SummaryNames.size() * Nt );
      ColumnNames.clear();
    }
    if( NumSamples_ != NumSamples || Nt_ != Nt || NumExtraSamples != NewNumExtraSamples )
    {
      NumSamples_ = NumSamples;
      Nt_ = Nt;
      NumExtraSamples = NewNumExtraSamples;
      if( Nt == 0 )
        m_pData.reset( nullptr );
      else
        m_pData.reset( new T[ static_cast<std::size_t>( NumSamples + NumExtraSamples ) * Nt ] );
      m_pRandNum.reset( nullptr );
    }
    if( pAuxNames )
      AuxNames = * pAuxNames;
    else
      AuxNames.clear();
  }
  bool IsFinite() { return Common::IsFinite( reinterpret_cast<scalar_type *>( m_pData.get() ),
      static_cast<size_t>( ( NumSamples_ + NumExtraSamples ) * ( SampleTraits<T>::is_complex ? 2 : 1 ) ) * Nt_ ); }
  template <typename U>
  void IsCompatible( const Sample<U> &o, int * pNumSamples = nullptr/*, bool bCompareBase = true*/ ) const;
  template <typename U> void CopyAttributes( const Sample<U> &o );
  Sample<T> Bootstrap(int NumBootSamples, SeedType Seed, const std::string * pMachineName = nullptr,
                      std::vector<std::string> * pAuxNames = nullptr);
  void Read (const std::string &FileName, const char *PrintPrefix = nullptr, std::vector<std::string> * pOpNames = nullptr,
             std::string * pGroupName = nullptr );
  void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void MakeCorrSummary( const char * pAvgName );
  void WriteSummary( const std::string &sOutFileName, bool bVerboseSummary = false );
  explicit Sample(int NumSamples = 0, int Nt = 0, std::vector<std::string> * pAuxNames = nullptr,
                  std::vector<std::string> * pSummaryNames = nullptr)
    : NumSamples_{0}, Nt_{0}, NumExtraSamples{0}
  {
    resize( NumSamples, Nt, pAuxNames );
    if( pSummaryNames )
      SetSummaryNames( *pSummaryNames );
  }
  Sample(const std::string &FileName, const char *PrintPrefix = nullptr, std::vector<std::string> * pOpNames = nullptr,
         std::string * pGroupName = nullptr)
      : NumSamples_{0}, Nt_{0}, NumExtraSamples{0} { Read( FileName, PrintPrefix, pOpNames, pGroupName ); }
  inline T * operator[]( int Sample )
  {
    if( Sample < -NumExtraSamples || Sample > NumSamples_ )
      throw std::out_of_range( "Sample " + std::to_string( Sample ) );
    return & m_pData[static_cast<std::size_t>( Sample + NumExtraSamples ) * Nt_];
  }
  inline const T * operator[]( int Sample ) const
  {
    if( Sample < -NumExtraSamples || Sample > NumSamples_ )
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
protected:
  template<typename ST> inline typename std::enable_if<SampleTraits<ST>::is_complex>::type
  CopyOldFormat( ST * Dest, const std::vector<double> & Src )
  {
    for( int t = 0; t < Nt_; t++ )
    {
      Dest->real( Src[t] );
      Dest->imag( Src[t + Nt_] );
      Dest++;
    }
  }
  template<typename ST> inline typename std::enable_if<!SampleTraits<ST>::is_complex>::type
  CopyOldFormat( ST * Dest, const std::vector<double> & Src )
  {
    for( int t = 0; t < Nt_; t++ )
    {
      if( Src[t + Nt_] )
        throw new std::runtime_error( "Complex sample has imaginary component" );
      *Dest++ = Src[t];
    }
  }
public: // Override these for specialisations
  virtual const std::string & DefaultGroupName() { return sBootstrap; }
  virtual bool bFolded() { return false; }
  // Descendants should call base first
  virtual void SummaryComments( std::ostream & s, bool bVerboseSummary = false ) const
  {
    static const char SeedPrefix[] = "# Seed: ";
    s << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
    if( Seed_ ) s << SeedPrefix << Seed_ << NewLine;
    else if( !Name_.SeedString.empty() ) s << SeedPrefix << Name_.SeedString << NewLine;
    if( !SeedMachine_.empty() ) s << "# Seed machine: " << SeedMachine_ << NewLine;
    s << "# Bootstrap: " << NumSamples_ << NewLine << "# SampleSize: " << SampleSize << NewLine;
    if( binSize != 1 ) s << "# BinSize: " << binSize << NewLine;
    if( !ConfigCount.empty() )
    {
      s << "# Configs: " << ConfigCount.size();
      for( const Common::ConfigCount &cc : ConfigCount )
        s << ", " << cc;
      s << NewLine;
    }
    if( bVerboseSummary )
    {
      s << "# FileCount: " << FileList.size() << NewLine;
      std::size_t i = 0;
      for( const std::string &f : FileList )
        s << "# File " << i++ << ": " << f << NewLine;
    }
  }
  virtual void ReadAttributes( ::H5::Group &g ) {}
  virtual void ValidateAttributes() {} // Called once data read to validate attributes against data
  virtual int WriteAttributes( ::H5::Group &g ) { return 0; }
};

using SampleC = Sample<std::complex<double>>;
using SampleD = Sample<double>;

// Initialise *pNumSamples either to 0, or to the size of the first Sample before first call
template <typename T> template <typename U>
void Sample<T>::IsCompatible( const Sample<U> &o, int * pNumSamples/*, bool bCompareBase*/ ) const
{
  static const std::string sPrefix{ "Incompatible " + sBootstrap + " samples - " };
  std::string sSuffix{ ":\n  " + Name_.Filename + "\n+ " + o.Name_.Filename };
  if( /*bCompareBase &&*/ !Common::EqualIgnoreCase( o.Name_.Base, Name_.Base ) )
    throw std::runtime_error( sPrefix + "base " + o.Name_.Base + sNE + Name_.Base + sSuffix );
  if( !Common::EqualIgnoreCase( o.Name_.Type, Name_.Type ) )
    throw std::runtime_error( sPrefix + "type " + o.Name_.Type + sNE + Name_.Type + sSuffix );
  if( !Common::EqualIgnoreCase( o.Name_.Ext, Name_.Ext ) )
    throw std::runtime_error( sPrefix + "extension " + o.Name_.Ext + sNE + Name_.Ext + sSuffix );
  if( pNumSamples )
  {
    if( *pNumSamples == 0 || *pNumSamples > NumSamples_)
      *pNumSamples = NumSamples_;
    if( *pNumSamples > o.NumSamples() )
      *pNumSamples = o.NumSamples();
  }
  else if( o.NumSamples() != NumSamples_ )
    throw std::runtime_error( sPrefix + "NumSamples " + std::to_string(o.NumSamples()) +
                             sNE + std::to_string(NumSamples_) + sSuffix );
  if( o.Nt() != Nt_ )
    throw std::runtime_error( sPrefix + "Nt " + std::to_string(o.Nt()) +
                             sNE + std::to_string(Nt_) + sSuffix );
  if( o.SampleSize && SampleSize && o.SampleSize != SampleSize )
    throw std::runtime_error( sPrefix + "SampleSize " + std::to_string(o.SampleSize) +
                             sNE + std::to_string(SampleSize) + sSuffix );
  if( m_pRandNum && o.m_pRandNum )
  {
    // Compare the actual random numbers
    const std::size_t Len{ static_cast<std::size_t>( NumSamples_ ) * SampleSize };
    const fint * l{ m_pRandNum.get() };
    const fint * r{ o.m_pRandNum.get() };
    std::size_t i = 0;
    while( i < Len && l[i] == r[i] )
      ++i;
    if( i != Len )
      throw std::runtime_error( "Random number seeds don't match" );
  }
  else
  {
    if( o.Seed_ != Seed_ )
      throw std::runtime_error( "Seed " + std::to_string( o.Seed_ ) + sNE + std::to_string( Seed_ ) );
    if( !Common::EqualIgnoreCase( o.SeedMachine_, SeedMachine_ ) )
      throw std::runtime_error( "Machine " + o.SeedMachine_ + sNE + SeedMachine_ );
  }
  const std::size_t CSize {   ConfigCount.size() };
  const std::size_t CSizeO{ o.ConfigCount.size() };
  if( CSize && CSizeO )
  {
    if( CSizeO != CSize )
      throw std::runtime_error( sPrefix + "Number of configs " +
                          std::to_string(CSizeO) + sNE + std::to_string(CSize) + sSuffix );
    for( std::size_t i = 0; i < CSize; i++ )
    {
      const Common::ConfigCount &l{   ConfigCount[i] };
      const Common::ConfigCount &r{ o.ConfigCount[i] };
      if( r.Config != l.Config )
        throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                 sNE + std::to_string(l.Config) + sSuffix );
      if( r.Count != l.Count )
        throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                 ", NumTimeslices " + std::to_string(r.Count) + sNE +
                                 std::to_string(l.Count) + sSuffix );
    }
  }
}

// Take ownership of the FileList
template <typename T>
template <typename U> void Sample<T>::CopyAttributes( const Sample<U> &in )
{
  binSize = in.binSize;
  SampleSize = in.SampleSize;
  ConfigCount = in.ConfigCount;
  if( FileList.empty() )
    FileList = in.FileList;
  // Copy the random numbers for this sample
  Seed_ = in.Seed_;
  SeedMachine_ = in.SeedMachine_;
  const fint * pSrc{ in.RandNum() };
  if( pSrc )
  {
    if( NumSamples_ > in.NumSamples() )
      throw std::runtime_error( "Can't copy random numbers from a smaller sample" );
    const std::size_t Len{ static_cast<std::size_t>( NumSamples_ ) * SampleSize };
    fint * pDst{ new fint[ Len ] };
    m_pRandNum.reset( pDst );
    const fint * pSrc{ in.RandNum() };
    for( std::size_t i = 0; i < Len; ++i )
      pDst[i] = pSrc[i];
  }
  else
    m_pRandNum.reset( nullptr );
}

/**
 Perform bootstrap
 Compute the mean on this sample, saving this as the central value
 Returned sample has specified extra fields ... but only the central value is copied from original
 */
template <typename T>
Sample<T> Sample<T>::Bootstrap(int NumBootSamples, SeedType Seed, const std::string * pMachineName,
                               std::vector<std::string> * pAuxNames)
{
  std::mt19937                        engine( Seed );
  std::uniform_int_distribution<fint> random( 0, NumSamples_ - 1 );
  Sample<T> boot( NumBootSamples, Nt_, pAuxNames ? pAuxNames : &AuxNames );

  // Compute the mean, then copy all the extra info to the bootstrap I will return
  MakeMean();
  std::copy( (*this)[idxCentral], (*this)[0], boot[idxCentral] );

  // Master thread can choose random samples while other threads bootstrap
  T * dst{ boot[0] };
  fint * rnd{ new fint[ static_cast<std::size_t>( NumBootSamples ) * NumSamples_ ] };
  boot.m_pRandNum.reset( rnd );
  auto start = std::chrono::steady_clock::now();
  int nthreads;
  #pragma omp parallel
  #pragma omp single
  {
    nthreads = omp_get_num_threads();
    for( int i = 0; i < NumBootSamples; i++, dst += Nt_, rnd += NumSamples_ )
    {
      // Generate the random numbers for this sample
      for( int s = 0; s < NumSamples_; s++ )
        rnd[s] = random( engine );
      // Schedule a task to perform the bootstrap for this replica
      #pragma omp task firstprivate( dst, rnd )
      {
        for( int s = 0; s < NumSamples_; s++ )
        {
          const T * src{ (*this)[ rnd[s] ] };
          for( int t = 0; t < Nt_; t++ )
          {
            if( !s )
              dst[t]  = *src++;
            else
              dst[t] += *src++;
          }
        }
        // Turn the sum into an average
        if( NumSamples_ > 1 )
          for( int t = 0; t < Nt_; t++)
            dst[t] /= NumSamples_;
      }
    }
  }
  auto diff = std::chrono::steady_clock::now() - start;
  std::cout << "  bootstrapping with " << nthreads << " threads took " << std::chrono::duration <double, std::milli> (diff).count() << " ms\n";
  boot.Seed_ = Seed;
  if( pMachineName && ! pMachineName->empty() )
    boot.SeedMachine_ = *pMachineName;
  else
    boot.SeedMachine_ = GetHostName();
  return boot;
}

template <typename T>
void Sample<T>::WriteSummary( const std::string &sOutFileName, bool bVerboseSummary )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  if( SummaryNames.empty() )
    throw std::runtime_error( "Summaries can't be written because they've not been created" );
  std::ofstream s( sOutFileName );
  SummaryHeader<T>( s, sOutFileName );
  SummaryComments( s, bVerboseSummary );
  // Now write all the field names
  const std::string Sep{ " " };
  s << "t";
  for( const std::string &n : SummaryNames )
    s << Sep << n << Sep << n << "_low " << n << "_high " << n << "_check";
  s << NewLine;
  const ValWithEr<scalar_type> * p{ m_pSummaryData.get() };
  for( int t = 0; t < Nt_; t++, p++ )
  {
    s << t;
    for( std::size_t f = 0; f < SummaryNames.size(); f++ )
      s << Sep << p[f * Nt_];
    s << NewLine;
  }
}

// Make a summary of the data
// If a name is specified, simply save the central replica
// Otherwise, calculate exp and cosh mass as well

template <typename T>
void Sample<T>::MakeCorrSummary( const char * pAvgName )
{
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  // Now perform summaries
  const int NumFields{ pAvgName ? 2 : static_cast<int>( sCorrSummaryNames.size() ) };
  if( pAvgName )
  {
    std::string s1{ pAvgName };
    std::string s2{ "bias" };
    std::vector<std::string> v{ s1, s2 };
    SetSummaryNames( v );
  }
  else
    SetSummaryNames( sCorrSummaryNames );
  const int tMid{ bFolded() ? Nt_ : Nt_ / 2 };
  std::vector<scalar_type> Data( NumSamples_ );
  ValWithEr<scalar_type> * pDest = m_pSummaryData.get();
  for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
  {
    for(int f = 0; f < NumFields; f++) // each field
    {
      for(int t = 0; t < Nt_; t++, pDest++ )
      {
        if( f != 1 )
        {
          // Index of previous and next timeslice relative to current timeslice
          const int Next{ ( t == Nt_ - 1 ? 1 - Nt_ :  1 ) * scalar_count };
          const int Prev{ ( t == 0       ? Nt_ - 1 : -1 ) * scalar_count };
          const scalar_type * p = Traits::ScalarPtr( (*this)[idxCentral] ) + t * scalar_count + i;
          int Count = 0;
          for(int n = idxCentral; n < NumSamples_; n++, p += Nt_ * scalar_count )
          {
            scalar_type d;
            switch( f )
            {
              case 0: // central value
                d = * p;
                break;
              case 2: // exponential mass
                if( t == 0 )
                  d = NaN;
                else if( t <= tMid )
                  d = std::log( p[Prev] / *p );
                else
                  d = -std::log( *p / p[Next] );
                break;
              case 3: // cosh mass
                d = std::acosh( ( p[Prev] + p[Next] ) / ( *p * 2 ) );
                break;
              default:
                d = 0;
            }
            if( n == idxCentral )
              pDest->Central = d;
            else if( std::isfinite( d ) )
              Data[Count++] = d;
          }
          pDest->Get( pDest->Central, Data, Count );
          if( f == 0 )
          {
            pDest[Nt_] = *pDest;
            scalar_type d = 0;
            for( int i = 0; i < Count; i++ )
              d += Data[i];
            d /= Count;
            pDest[Nt_].Central = d;
          }
        }
      }
    }
  }
}

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void Sample<T>::Read(const std::string &FileName, const char *PrintPrefix, std::vector<std::string> * pOpNames, std::string *pGroupName)
{
  Name_.Parse( FileName, pOpNames );
  if( !Name_.bSeedNum )
    throw std::runtime_error( "Seed missing from " + FileName );
  if( Name_.Type.empty() )
    throw std::runtime_error( "Type missing from " + FileName );
  ::H5::H5File f;
  ::H5::Group  g;
  H5::OpenFileGroup( f, g, FileName, PrintPrefix, pGroupName );
  bool bOK{ false };
  bool bFiniteError{ false };
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  try // to load from LatAnalyze format
  {
    unsigned short att_Type;
    ::H5::Attribute a;
    a = g.openAttribute("type");
    a.read( ::H5::PredType::NATIVE_USHORT, &att_Type );
    a.close();
    if( att_Type == 2 )
    {
      unsigned long  att_nSample;
      a = g.openAttribute("nSample");
      a.read( ::H5::PredType::NATIVE_ULONG, &att_nSample );
      a.close();
      std::vector<double> buffer;
      for( int i = idxCentral; i < NumSamples_; i++ )
      {
        ::H5::DataSet ds = g.openDataSet( "data_" + ( i == idxCentral ? "C" : "S_" + std::to_string( i ) ) );
        ::H5::DataSpace dsp = ds.getSpace();
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
          ds.read( buffer.data(), ::H5::PredType::NATIVE_DOUBLE );
          if( !Common::IsFinite( buffer ) )
          {
            bFiniteError = true;
            break;
          }
          T * p { (*this)[i] };
          CopyOldFormat( p, buffer );
          // Check whether this is the end
          if( i == NumSamples_ - 1 )
          {
            bOK = true;
            p = (*this)[-NumExtraSamples];
            // This format has no auxiliary rows - set them to zero
            for( int j = -NumExtraSamples; j < idxCentral; j++ )
              for( int t = 0; t < Nt_ ; t++ )
                *p++ = 0;
          }
        }
      }
    }
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
    ::H5::Exception::clearErrorStack();
  }
  if( !bOK && !bFiniteError )
  {
    try // to load from my format
    {
      unsigned short att_nAux;
      ::H5::Attribute a;
      a = g.openAttribute("nAux");
      a.read( ::H5::PredType::NATIVE_USHORT, &att_nAux );
      a.close();
      unsigned int att_nSample;
      a = g.openAttribute("nSample");
      a.read( ::H5::PredType::NATIVE_UINT, &att_nSample );
      a.close();
      Seed_ = 0;
      try
      {
        unsigned int tmp;
        a = g.openAttribute( sSeed );
        a.read( ::H5::PredType::NATIVE_UINT, &tmp );
        a.close();
        Seed_ = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      SeedMachine_.clear();
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute( sSeedMachine );
        a.read( a.getStrType(), SeedMachine_ );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      std::vector<std::string> myAuxNames;
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sAuxNames);
        myAuxNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      std::vector<std::string> mySummaryNames;
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sSummaryNames);
        mySummaryNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      std::vector<std::string> myColumnNames;
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sColumnNames);
        myColumnNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      SampleSize = 0;
      try
      {
        int tmp;
        a = g.openAttribute( sSampleSize );
        a.read( ::H5::PredType::NATIVE_INT, &tmp );
        a.close();
        SampleSize = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      binSize = 1;
      try
      {
        int tmp;
        a = g.openAttribute( sBinSize );
        a.read( ::H5::PredType::NATIVE_INT, &tmp );
        a.close();
        binSize = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      std::vector<std::string> myFileList;
      try
      {
        a = g.openAttribute(sFileList);
        myFileList = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      std::vector<Common::ConfigCount> myConfigCount;
      try
      {
        a = g.openAttribute(sConfigCount);
        ::H5::DataSpace dsp = a.getSpace();
        const int rank{ dsp.getSimpleExtentNdims() };
        if( rank != 1 )
          throw std::runtime_error( sConfigCount + " number of dimensions=" + std::to_string( rank ) + ", expecting 1" );
        hsize_t NumConfig;
        dsp.getSimpleExtentDims( &NumConfig );
        if( NumConfig > std::numeric_limits<int>::max() )
          throw std::runtime_error( sConfigCount + " too many items in ConfigCount: " + std::to_string( NumConfig ) );
        myConfigCount.resize( NumConfig );
        a.read( H5::Equiv<Common::ConfigCount>::Type, &myConfigCount[0] );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      const unsigned short ExpectedNAux{ static_cast<unsigned short>( myAuxNames.size() ? myAuxNames.size() - 1 : 0 ) };
      if( att_nAux == ExpectedNAux )
      {
        try
        {
          ReadAttributes( g );
        }
        catch(const std::exception &e)
        {
          throw;
        }
        catch(...)
        {
          throw std::runtime_error( "Extended attributes for " + DefaultGroupName() + " missing" );
        }
        ::H5::DataSet ds = g.openDataSet( "data_C" );
        ::H5::DataSpace dsp = ds.getSpace();
        int nDims{ dsp.getSimpleExtentNdims() };
        if( ( att_nAux == 0 && nDims == 1 ) || ( att_nAux != 0 && nDims == 2 ) )
        {
          hsize_t Dim[2];
          dsp.getSimpleExtentDims( Dim );
          if( ( att_nAux == 0 && Dim[0] * ( att_nSample + 1 ) <= std::numeric_limits<int>::max() )
           || ( att_nAux != 0 && Dim[0] == att_nAux + 1 && Dim[1] * ( att_nAux + 1 + att_nSample ) <= std::numeric_limits<int>::max() ) )
          {
            SummaryNames = std::move( mySummaryNames );
            resize( static_cast<int>( att_nSample ), static_cast<int>( Dim[ att_nAux == 0 ? 0 : 1 ] ), &myAuxNames );
            ColumnNames = std::move( myColumnNames );
            FileList = std::move( myFileList );
            ConfigCount = std::move( myConfigCount );
            ds.read( (*this)[-NumExtraSamples], H5::Equiv<T>::Type );
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
                ds.read( (*this)[0], H5::Equiv<T>::Type );
                if( !IsFinite() )
                  bFiniteError = true;
                else
                {
                  try
                  {
                    dsp.close();
                    ds.close();
                    ds = g.openDataSet( sRandom );
                    dsp = ds.getSpace();
                    if( dsp.getSimpleExtentNdims() == 2 )
                    {
                      dsp.getSimpleExtentDims( Dim );
                      if( Dim[0] == att_nSample && Dim[1] == SampleSize )
                      {
                        fint * rnd{ new fint[ static_cast<std::size_t>( att_nSample ) * SampleSize ] };
                        m_pRandNum.reset( rnd );
                        ds.read( m_pRandNum.get(), H5::Equiv<std::uint_fast32_t>::Type );
                      }
                    }
                  }
                  catch(const ::H5::Exception &)
                  {
                    ::H5::Exception::clearErrorStack();
                    m_pRandNum.reset( nullptr );
                  }
                  if( SummaryNames.empty() )
                    bOK = true;
                  else
                  {
                    dsp.close();
                    ds.close();
                    ds = g.openDataSet( "Summary" );
                    dsp = ds.getSpace();
                    nDims = dsp.getSimpleExtentNdims();
                    if( nDims == 2 )
                    {
                      dsp.getSimpleExtentDims( Dim );
                      if( Dim[0] == SummaryNames.size() && Dim[1] == static_cast<unsigned int>(Nt_) )
                      {
                        ds.read( m_pSummaryData.get(), H5::Equiv<ValWithEr<scalar_type>>::Type );
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
      if( bOK )
        ValidateAttributes();
    }
    catch(const ::H5::Exception &)
    {
      bOK = false;
      ::H5::Exception::clearErrorStack();
    }
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  if( bFiniteError )
    throw std::runtime_error( "NANs in " + FileName );
  if( !bOK )
    throw std::runtime_error( "Unable to read sample from " + FileName );
}

template <typename T>
void Sample<T>::Write( const std::string &FileName, const char * pszGroupName )
{
  const std::string GroupName{ pszGroupName == nullptr || *pszGroupName == 0 ? DefaultGroupName() : pszGroupName };
  bool bOK = false;
  try // to write in my format
  {
    ::H5::H5File f( FileName, H5F_ACC_TRUNC );
    hsize_t Dims[2];
    Dims[0] = 1;
    ::H5::DataSpace ds1( 1, Dims );
    ::H5::Group g = f.createGroup( GroupName );
    ::H5::Attribute a = g.createAttribute( "nAux", ::H5::PredType::STD_U16LE, ds1 );
    int tmp{ NumExtraSamples - 1 };
    a.write( ::H5::PredType::NATIVE_INT, &tmp );
    a.close();
    a = g.createAttribute( "nSample", ::H5::PredType::STD_U32LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &NumSamples_ );
    a.close();
    int NumAttributes = 2;
    if( Seed_ )
    {
      const unsigned int tmp{ Seed_ };
      a = g.createAttribute( sSeed, sizeof( Seed_ ) == 4 ? ::H5::PredType::STD_U32LE : ::H5::PredType::STD_U64LE, ds1 );
      a.write( ::H5::PredType::NATIVE_UINT, &tmp );
      a.close();
      NumAttributes++;
    }
    if( !SeedMachine_.empty() )
    {
      a = g.createAttribute( sSeedMachine, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, SeedMachine_ );
      a.close();
      NumAttributes++;
    }
    if( AuxNames.size() )
    {
      H5::WriteAttribute( g, sAuxNames, AuxNames );
      NumAttributes++;
    }
    if( SummaryNames.size() )
    {
      H5::WriteAttribute( g, sSummaryNames, SummaryNames );
      NumAttributes++;
    }
    if( ColumnNames.size() )
    {
      H5::WriteAttribute( g, sColumnNames, ColumnNames );
      NumAttributes++;
    }
    if( SampleSize )
    {
      a = g.createAttribute( sSampleSize, ::H5::PredType::STD_U32LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &SampleSize );
      a.close();
      NumAttributes++;
    }
    if( binSize != 1 )
    {
      a = g.createAttribute( sBinSize, ::H5::PredType::STD_I32LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &binSize );
      a.close();
      NumAttributes++;
    }
    ::H5::DataSpace dsp;
    if( ConfigCount.size() )
    {
      Dims[0] = ConfigCount.size();
      dsp = ::H5::DataSpace( 1, Dims );
      a = g.createAttribute( sConfigCount, H5::Equiv<Common::ConfigCount>::Type, dsp );
      a.write( H5::Equiv<Common::ConfigCount>::Type, &ConfigCount[0] );
      a.close();
      dsp.close();
      NumAttributes++;
    }
    if( FileList.size() )
    {
      H5::WriteAttribute( g, sFileList, FileList );
      NumAttributes++;
    }
    NumAttributes += WriteAttributes( g );
    if( NumExtraSamples == 1 )
    {
      Dims[0] = Nt_;
      dsp = ::H5::DataSpace( 1, Dims );
    }
    else
    {
      Dims[0] = NumExtraSamples;
      Dims[1] = Nt_;
      dsp = ::H5::DataSpace( 2, Dims );
    }
    ::H5::DataSet ds = g.createDataSet( "data_C", H5::Equiv<T>::Type, dsp );
    ds.write( (*this)[-NumExtraSamples], H5::Equiv<T>::Type );
    ds.close();
    dsp.close();
    Dims[0] = NumSamples_;
    Dims[1] = Nt_;
    dsp = ::H5::DataSpace( 2, Dims );
    ds = g.createDataSet( "data_S", H5::Equiv<T>::Type, dsp );
    ds.write( (*this)[0], H5::Equiv<T>::Type );
    ds.close();
    dsp.close();
    if( SummaryNames.size() )
    {
      Dims[0] = SummaryNames.size();
      Dims[1] = Nt_;
      dsp = ::H5::DataSpace( 2, Dims );
      ds = g.createDataSet( "Summary", H5::Equiv<ValWithEr<scalar_type>>::Type, dsp );
      ds.write( m_pSummaryData.get(), H5::Equiv<ValWithEr<scalar_type>>::Type );
      ds.close();
      dsp.close();
    }
    if( m_pRandNum )
    {
      Dims[0] = NumSamples_;
      Dims[1] = SampleSize;
      dsp = ::H5::DataSpace( 2, Dims );
      // The values in this table go from 0 ... NumSamples_ + 1, so choose a space-minimising sie in the file
      ds = g.createDataSet( sRandom, SampleSize<=static_cast<int>(std::numeric_limits<std::uint8_t>::max())+1
        ? ::H5::PredType::STD_U8LE : SampleSize<=static_cast<int>(std::numeric_limits<std::uint16_t>::max())+1
        ? ::H5::PredType::STD_U16LE: ::H5::PredType::STD_U32LE, dsp );
      ds.write( m_pRandNum.get(), H5::Equiv<std::uint_fast32_t>::Type );
      ds.close();
      dsp.close();
    }
    g = f.openGroup( "/" );
    a = g.createAttribute( "_Grid_dataset_threshold", ::H5::PredType::STD_U32LE, ds1);
    a.write( ::H5::PredType::NATIVE_INT, &NumAttributes );
    a.close();
    ds1.close();
    g.close();
    bOK = true;
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
  }
  if( !bOK )
    throw std::runtime_error( "Unable to write sample to " + FileName + ", group " + GroupName );
}

template <typename T>
class Fold : public Sample<T>
{
public:
  int NtUnfolded = 0;
  Reality reality = Reality::Unknown;
  Parity parity = Parity::Unknown;
  Sign sign = Sign::Unknown;
  bool t0Negated = false;
  bool Conjugated = false;
  std::vector<std::string> BootstrapList;
  using Base = Sample<T>;
  using Base::Base;
  /*explicit Fold( int NumSamples = 0, int Nt = 0 ) : Base::Sample( NumSamples,  Nt ) {}
  Fold( const std::string &FileName, std::string &GroupName, const char *PrintPrefix = nullptr,
        std::vector<std::string> * pOpNames = nullptr )
  : Base::Sample( FileName, GroupName, PrintPrefix, pOpNames ) {}*/
  virtual const std::string & DefaultGroupName() { return sFold; }
  virtual bool bFolded() { return true; }
  virtual void SummaryComments( std::ostream & s ) const
  {
    Base::SummaryComments( s );
    if( NtUnfolded ) s << "# NtUnfolded: " << NtUnfolded << NewLine;
    if( reality != Reality::Unknown ) s << "# Reality: " << reality << NewLine;
    if( parity != Parity::Unknown ) s << "# Parity: " << parity << NewLine;
    if( sign != Sign::Unknown ) s << "# Sign: " << sign << NewLine;
    if( t0Negated ) s << "# timeslice 0 negated: true" << NewLine;
    if( Conjugated ) s << "# conjugate operator average: true" << NewLine;
  }
  virtual void ReadAttributes( ::H5::Group &g )
  {
    Base::ReadAttributes( g );
    ::H5::Attribute a;
    NtUnfolded = 0;
    try
    {
      a = g.openAttribute(sNtUnfolded);
      a.read( ::H5::PredType::NATIVE_INT, &NtUnfolded );
      a.close();
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
    reality = Reality::Unknown;
    try
    {
      a = g.openAttribute(sReality);
      std::string s{};
      a.read( a.getStrType(), s );
      a.close();
      if( EqualIgnoreCase( s, Reality_TextEquiv_Real ) )
        reality = Reality::Real;
      else if( EqualIgnoreCase( s, Reality_TextEquiv_Imaginary ) )
        reality = Reality::Imag;
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
    parity = Parity::Unknown;
    try
    {
      a = g.openAttribute(sParity);
      std::string s{};
      a.read( a.getStrType(), s );
      a.close();
      if( EqualIgnoreCase( s, Parity_TextEquiv_Even ) )
        parity = Parity::Even;
      else if( EqualIgnoreCase( s, Parity_TextEquiv_Odd ) )
        parity = Parity::Odd;
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
    sign = Sign::Unknown;
    try
    {
      a = g.openAttribute(sSign);
      std::string s{};
      a.read( a.getStrType(), s );
      a.close();
      if( EqualIgnoreCase( s, Sign_TextEquiv_Positive ) )
        sign = Sign::Positive;
      else if( EqualIgnoreCase( s, Sign_TextEquiv_Negative ) )
        sign = Sign::Negative;
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
    t0Negated = false;
    std::int8_t i8;
    try
    {
      a = g.openAttribute(st0Negated);
      a.read( ::H5::PredType::NATIVE_INT8, &i8 );
      a.close();
      if( i8 )
        t0Negated = true;
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
    Conjugated = false;
    try
    {
      a = g.openAttribute(sConjugated);
      a.read( ::H5::PredType::NATIVE_INT8, &i8 );
      a.close();
      if( i8 )
        Conjugated = true;
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
    BootstrapList.clear();
    try
    {
      a = g.openAttribute(sBootstrapList);
      BootstrapList = H5::ReadStrings( a );
      a.close();
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
    }
  }
  virtual int WriteAttributes( ::H5::Group &g )
  {
    int iReturn = Sample<T>::WriteAttributes( g );
    const hsize_t OneDimension{ 1 };
    ::H5::DataSpace ds1( 1, &OneDimension );
    ::H5::Attribute a;
    if( NtUnfolded )
    {
      a = g.createAttribute( sNtUnfolded, ::H5::PredType::STD_U16LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &NtUnfolded );
      a.close();
      iReturn++;
    }
    if( reality == Reality::Real || reality == Reality::Imag )
    {
      const std::string &s{reality==Reality::Real?Reality_TextEquiv_Real:Reality_TextEquiv_Imaginary};
      a = g.createAttribute( sReality, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, s );
      a.close();
      iReturn++;
    }
    if( parity == Parity::Even || parity == Parity::Odd )
    {
      const std::string &s{ parity == Parity::Even ? Parity_TextEquiv_Even : Parity_TextEquiv_Odd };
      a = g.createAttribute( sParity, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, s );
      a.close();
      iReturn++;
    }
    if( sign == Sign::Positive || sign == Sign::Negative )
    {
      const std::string &s{ sign == Sign::Positive ? Sign_TextEquiv_Positive : Sign_TextEquiv_Negative};
      a = g.createAttribute( sSign, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, s );
      a.close();
      iReturn++;
    }
    const std::int8_t i8{ 1 };
    if( t0Negated )
    {
      a = g.createAttribute( st0Negated, ::H5::PredType::STD_U8LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT8, &i8 );
      a.close();
      iReturn++;
    }
    if( Conjugated )
    {
      a = g.createAttribute( sConjugated, ::H5::PredType::STD_U8LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT8, &i8 );
      a.close();
      iReturn++;
    }
    if( BootstrapList.size() )
    {
      H5::WriteAttribute( g, sBootstrapList, BootstrapList );
      iReturn++;
    }
    return iReturn;
  }
};

template <typename T>
class Model : public Sample<T>
{
public:
  int NumExponents = 0;
  int NumFiles = 0; // TODO: this duplicates FileList ... but FileList optional and added later ...
  int ti = 0;
  int tf = 0;
  int dof = 0;
  bool Factorised = false;
  bool CovarFrozen = false;
  std::vector<std::string> OpNames;
  using Base = Sample<T>;
  Model() : Base::Sample{} {}
  Model( std::vector<std::string> OpNames_, int NumExponents_, int NumFiles_, int ti_, int tf_,
         int dof_, bool Factorised_, bool CovarFrozen_, int NumSamples, int Nt )
  : Base::Sample{ NumSamples, Nt }, OpNames{ OpNames_ }, NumExponents{ NumExponents_ },
    NumFiles{ NumFiles_}, ti{ ti_ }, tf{ tf_ }, dof{ dof_ }, Factorised{ Factorised_ }, CovarFrozen{ CovarFrozen_ } {}
  virtual const std::string & DefaultGroupName() { return sModel; }
  // Are the new parameters compatible with the pre-performed fit?
  inline bool Compatible(int numExponents_, int Dof_, int Ti_, int Tf_,
                         const std::vector<std::string> &opNames_) const
  {
    return std::tie( NumExponents , dof , ti , tf , OpNames  )
        == std::tie( numExponents_, Dof_, Ti_, Tf_, opNames_ );
  }
  inline bool NewParamsMorePrecise( bool covarFrozen_, int NumSamples ) const
  {
    // Unfreezing the covariance matrix makes a big difference!
    if( ( CovarFrozen && !covarFrozen_ ) || ( !CovarFrozen && covarFrozen_ ) )
      return CovarFrozen;
    return NumSamples > Base::NumSamples_;
  }
  void Read ( const std::string &FileName, const char *PrintPrefix = nullptr, std::string *pGroupName =nullptr )
  {
    Base::Read( FileName, PrintPrefix, nullptr, pGroupName );
    Base::Name_.ParseExtra();
  }
  virtual void ReadAttributes( ::H5::Group &g )
  {
    Base::ReadAttributes( g );
    ::H5::Attribute a;
    a = g.openAttribute(sNumExponents);
    a.read( ::H5::PredType::NATIVE_INT, &NumExponents );
    a.close();
    a = g.openAttribute(sNumFiles);
    a.read( ::H5::PredType::NATIVE_INT, &NumFiles );
    a.close();
    a = g.openAttribute(sTI);
    a.read( ::H5::PredType::NATIVE_INT, &ti );
    a.close();
    a = g.openAttribute(sTF);
    a.read( ::H5::PredType::NATIVE_INT, &tf );
    a.close();
    a = g.openAttribute(sDoF);
    a.read( ::H5::PredType::NATIVE_INT, &dof );
    a.close();
    a = g.openAttribute(sOperators);
    OpNames = H5::ReadStrings( a );
    a.close();
    std::int8_t i8;
    a = g.openAttribute(sFactorised);
    a.read( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    Factorised = ( i8 != 0 );
    a = g.openAttribute(sCovarFrozen);
    a.read( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    CovarFrozen = ( i8 != 0 );
  }
  virtual void ValidateAttributes()
  {
    Base::ValidateAttributes();
    const int NumOps{ static_cast<int>( OpNames.size() ) };
    const int NumExpected{ NumExponents * ( NumOps + 1 ) + 1 };
    if( Base::Nt_ != NumExpected )
    {
      std::ostringstream s;
      s << "Have " << Base::Nt_ << " parameters, but expected " << NumExpected << ", i.e. "
        << NumExponents << " exponents * ( " << NumOps << " operators + 1 energy ) + chi squared per degree of freedom";
      throw std::runtime_error( s.str().c_str() );
    }
  }
  virtual int WriteAttributes( ::H5::Group &g )
  {
    int iReturn = Base::WriteAttributes( g ) + 8;
    const hsize_t OneDimension{ 1 };
    ::H5::DataSpace ds1( 1, &OneDimension );
    ::H5::Attribute a;
    a = g.createAttribute( sNumExponents, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &NumExponents );
    a.close();
    a = g.createAttribute( sNumFiles, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &NumFiles );
    a.close();
    a = g.createAttribute( sTI, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &ti );
    a.close();
    a = g.createAttribute( sTF, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &tf );
    a.close();
    a = g.createAttribute( sDoF, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &dof );
    a.close();
    std::int8_t i8{ static_cast<std::int8_t>( Factorised ? 1 : 0 ) };
    a = g.createAttribute( sFactorised, ::H5::PredType::STD_U8LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    i8 = static_cast<std::int8_t>( CovarFrozen ? 1 : 0 );
    a = g.createAttribute( sCovarFrozen, ::H5::PredType::STD_U8LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    H5::WriteAttribute( g, sOperators, OpNames );
    return iReturn;
  }
  virtual void SummaryComments( std::ostream & s ) const
  {
    Base::SummaryComments( s );
    s << "# NumExponents: " << NumExponents << NewLine
      << "# NumFiles: " << NumFiles << NewLine
      << "# Factorised: " << Factorised << NewLine
      << "# Frozen covariance: " << CovarFrozen << NewLine;
  }
};

// Read a complex array from an HDF5 file
template<typename T>
void ReadArray(std::vector<T> &buffer, const std::string &FileName, const std::string &ObjectName = std::string( "correlator" ),
               const char *PrintPrefix = nullptr, std::string * pGroupName = nullptr );

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
      throw std::invalid_argument("Switch " + Switch + ( Subscript ? "[" + std::to_string( Subscript ) + "]" : "" )
                                  + " not found");
    return FromString<T>( v[Subscript] );
  }
};

std::ostream& operator<<( std::ostream& os, const CommandLine &cl);

END_COMMON_NAMESPACE
#endif // Common_hpp
