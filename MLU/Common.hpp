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

#include <MLUconfig.h>
#include <MLU/FitRange.hpp>
#include <MLU/GSLVecMat.hpp>
#include <MLU/Param.hpp>

// std c++
#include <array>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <ios>
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

// Default output file extension for text based summaries of binary data
#ifndef TEXT_EXT
#define TEXT_EXT "txt"
#endif

#define BEGIN_COMMON_NAMESPACE namespace Common {
#define END_COMMON_NAMESPACE   };

extern "C" const char * MLUVersionInfoHuman();

BEGIN_COMMON_NAMESPACE

// Compatibility flags for filename comparisons

static constexpr unsigned int COMPAT_DEFAULT{ 0 };
static constexpr unsigned int COMPAT_DISABLE_BASE{ 1 };
static constexpr unsigned int COMPAT_DISABLE_NT{ 2 };
/*static constexpr unsigned int COMPAT_DISABLE_CONFIG_COUNT{ 4 };*/
static constexpr unsigned int COMPAT_DISABLE_TYPE{ 8 };

template <typename T> int sgn( T x )
{
  return (T(0) < x) - (x < T(0));
}

// Text required for summaries of correlators
namespace CorrSumm {
  extern const char sep[];
  extern const char Comment[];
  static constexpr int NumFields{ 3 };
  extern const char * FieldNames[NumFields];
};

extern const std::string sBootstrap;
extern const std::string sFold;
extern const std::string sModel;
extern const std::string sParams;
extern const std::string sCovmat;
extern const std::string sCovmatIn;
extern const std::string sCovmatInv;
extern const std::string sCovmatInvCholesky;
extern const std::string sCormat;
extern const std::string sCormatCholesky;
extern const std::string sCormatInvCholesky;
extern const std::string sNtUnfolded;
extern const std::string st0Negated;
extern const std::string sConjugated;
extern const std::string sRawBootstrap;
extern const std::string sTI;
extern const std::string sTF;
extern const std::string sDoF;
extern const std::string sNumExponents;
extern const std::string sNumFiles;
extern const std::string sCovarFrozen;
extern const std::string s_C;
extern const std::string sStdErrorMean;
extern const std::string sErrorScaled;
extern const std::string sCovariance;
extern const std::string sCovarianceIn;
extern const std::string sCovarianceInv;
extern const std::string sCovarianceInvCholesky;
extern const std::string sCorrelation;
extern const std::string sCorrelationCholesky;
extern const std::string sCorrelationInvCholesky;
extern const std::string sFitInput;
extern const std::string sModelPrediction;
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
extern const std::string sCovarSource;
extern const std::string sCovarRebin;
extern const std::string sCovarSampleSize;
extern const std::string sCovarNumBoot;
extern const std::string sGuess;
extern const std::string sParam;
extern const std::string sFitTime;
extern const std::string sNE;
extern const std::vector<std::string> sCorrSummaryNames;
extern const std::string sChiSqPerDof;
extern const std::string sPValue;
extern const std::string sPValueH;

// String containing switch name to enable alternate overlap normalisation
// (see definition of ModelOverlap in Fit/ModelCommon.hpp)
extern const std::string sOverlapAltNorm;

extern const std::vector<std::string> DefaultModelStats;

using SeedType = unsigned int;

template<typename T = int> struct StartStopStepIterator;

template<typename T = int> struct StartStopStep
{
  using Scalar = T;
  static constexpr Scalar Default = 1; // This is the default step if none specified
  using Iterator = StartStopStepIterator<T>;
  Scalar Start;
  Scalar Stop;
  Scalar Step;
  std::vector<Scalar> AlwaysIterate; // A list of special numbers that will always be iterated
public:
  Scalar StepValue( Scalar step_ ) const { return Start == Stop ? 0 : ( Stop < Start ? -1 : 1 ) * ( step_ == 0 ? Default : std::abs( step_ ) ); }
  void   SetStep  ( Scalar step_ ) { Step = StepValue( step_ ); }
  void   SetSingle( Scalar Single ) { Start = Single; Stop = Single; SetStep( 0 ); }
  StartStopStep( Scalar start_, Scalar stop_, Scalar step_ = Default ) : Start{start_}, Stop{stop_}, Step{ StepValue( step_ ) } {}
  StartStopStep() : StartStopStep( 0, 0 ) {}
  inline bool InRange( Scalar x ) const
  {
    Scalar Min;
    Scalar Max;
    if( Stop < Start )
    {
      Min = Stop;
      Max = Start;
    }
    else
    {
      Min = Start;
      Max = Stop;
    }
    return x >= Min && x <= Max;
  }
  inline std::string to_string( const std::string &Sep = Colon ) const
    { return std::to_string(Start) + Sep + std::to_string(Stop) + Sep + std::to_string(Step); }
  Iterator begin() const;
  Iterator end() const;
};
template<typename T = int> std::ostream & operator<<( std::ostream &os, const StartStopStep<T> &sss )
{
  return os << sss.to_string();
}

// Parse a Start:Stop[:Step] combination
template<typename T = int> std::istream & operator>>( std::istream &is,       StartStopStep<T> &sss )
{
  const char Sep = ':';
  T Start, Stop, Step = 0;
  bool bOK = false;
  if( is >> Start && NextCharIs( is, Sep ) && is >> Stop )
  {
    bool bGotStep{ NextCharIs( is, Sep ) };
    if( !bGotStep || is >> Step )
    {
      bOK = true;
      sss.Start = Start;
      sss.Stop  = Stop;
      sss.SetStep( Step );
    }
  }
  if( !bOK )
    is.setstate( std::ios_base::failbit );
  return is;
}

template<typename T> struct StartStopStepIterator : public StartStopStep<T>
{
  friend class StartStopStep<T>;
public:
  using Scalar = T;
  using Iterator = StartStopStepIterator<T>;
  using Base = StartStopStep<T>;
  static constexpr Scalar Default = Base::Default;
  //using Base::Base;
protected:
  std::vector<bool> IterateAtEnd; // true if we still need to iterate through corresponding value at end of iteration
  bool   bAtEnd = true;
  Scalar Position;
  std::size_t IdxIterated = std::numeric_limits<std::size_t>::max();
public:
  StartStopStepIterator() : Base() {}
  StartStopStepIterator( const Base &sss ) : Base( sss ), IterateAtEnd( Base::AlwaysIterate.size(), true )
  {
    bAtEnd = false;
    Position = Base::Start;
  }
  inline bool AtEnd() const { return bAtEnd && IdxIterated >= IterateAtEnd.size(); }
  inline Scalar operator*() const
  {
    if( AtEnd() )
      throw std::runtime_error( "Cannot dereference StartStopStepIterator at end" );
    return bAtEnd ? Base::AlwaysIterate[IdxIterated] : Position;
  }
  Iterator &operator++()
  {
    if( bAtEnd )
    {
      if( IdxIterated < IterateAtEnd.size() )
        ++IdxIterated;
    }
    else
    {
      // Tick items I've already iterated off my list
      for( IdxIterated = 0; IdxIterated < Base::AlwaysIterate.size(); ++IdxIterated )
        if( Position == Base::AlwaysIterate[IdxIterated] )
        {
          IterateAtEnd[IdxIterated] = false;
          break;
        }
      // Move on to the next item
      if( Base::Step )
      {
        Position += Base::Step;
        if( ( Base::Step > 0 && Position <= Base::Stop ) || ( Base::Step < 0 && Position >= Base::Stop ) )
          return *this;
      }
      // Falling through - I'm now going to see whether there are any special numbers I need to also iterate
      bAtEnd = true;
      IdxIterated = 0;
    }
    // Skip forward to the next special number that I haven't ticked off my list
    while( IdxIterated < IterateAtEnd.size() && !IterateAtEnd[IdxIterated] )
      ++IdxIterated;
    return *this;
  }
  inline Iterator operator++(int) { Iterator it(*this); return ++it; }
  inline bool operator==( const Iterator & o ) const
  {
    bool bMe{ AtEnd() };
    bool bOther{ o.AtEnd() };
    bool bEq{ bMe == bOther };
    if( bEq && !bMe )
      bEq = **this == *o;
    return bEq;
  }
  inline bool operator!=( const Iterator & o ) const { return !( *this == o ) ; }
  inline std::string to_string() const { return bAtEnd ? "end" : std::to_string( Position ); }
};

template<typename T> StartStopStepIterator<T> StartStopStep<T>::begin() const
{
  Iterator it( *this );
  it.bAtEnd = false;
  it.Position = it.Start;
  return it;
};

template<typename T> StartStopStepIterator<T> StartStopStep<T>::end() const
{
  Iterator it( *this );
  it.bAtEnd = true;
  return it;
};

// Does the specified file exist?
bool FileExists( const std::string& Filename );

// This is a key/value file, and we ignore anything after a comment character
// Use a std::MultiMap if the keys aren't unique
template<class Key, class T, class Compare = std::less<Key>, class Map = std::map<Key, T, Compare>>
class KeyValReader
{
protected:
  template<class V = T> static typename std::enable_if< std::is_same<V, std::string>::value, bool>::type
  GetValue( std::istringstream &iss, V &v )
  {
    bool bOK = static_cast<bool>( std::getline( iss, v ) );
    if( bOK )
      Common::Trim( v );
    return bOK;
  }
  template<class V = T> static typename std::enable_if<!std::is_same<V, std::string>::value, bool>::type
  GetValue( std::istringstream &iss, V &v ) { return iss >> v && Common::StreamEmpty( iss ); }
  template<class V = T> static void MyEmplace( std::multimap<Key, V, Compare> &m, Key &&k, V &&t )
  {
    // Multimap, so won't fail
    m.emplace( std::make_pair( std::move( k ), std::forward<V>( t ) ) );
  }
  template<class V = T> static void MyEmplace( std::map<Key, V, Compare> &m, Key &&k, V &&t )
  {
    // Map - might be already present
    using iterator = typename std::map<Key, V, Compare>::iterator;
    using pair = std::pair<iterator, bool>;
    pair a{ m.emplace( std::make_pair( std::move( k ), std::forward<V>( t ) ) ) };
    if( !a.second )
    {
      std::ostringstream es;
      es << "KeyValReader::ReadLine() duplicate key " << a.first->first;
      throw std::runtime_error( es.str().c_str() );
    }
  }
  static void ReadLine( Map &m, const std::string &s, bool bOptional = false )
  {
    std::istringstream iss( s );
    if( !Common::StreamEmpty( iss ) )
    {
      Key k;
      if( !( iss >> k ) )
        throw std::runtime_error( "KeyValReader::ReadLine() bad key " + s );
      T t{};
      if( !bOptional || !Common::StreamEmpty( iss ) )
      {
        if( !GetValue( iss, t ) )
          throw std::runtime_error( "KeyValReader::ReadLine() bad value " + s );
      }
      MyEmplace( m, std::move( k ), std::move( t ) );
    }
  }
public:
  static Map Read( const std::vector<std::string> &v, const std::string *Separator = &EqualSign,
                   bool bOptional = false )
  {
    Map m;
    std::string Tmp;
    for( const std::string &s : v )
    {
      const std::size_t pos = Separator ? s.find_first_of( *Separator ) : std::string::npos;
      if( pos != std::string::npos )
      {
        Tmp = s;
        Tmp[pos] = ' ';
      }
      ReadLine( m, ( pos == std::string::npos ) ? s : Tmp, bOptional );
    }
    return m;
  }
  static Map Read( std::istream &is, const std::string &sComment = Hash )
  {
    if( is.bad() )
      throw std::runtime_error( "KeyValReader::Read() input stream bad" );
    Map m;
    while( !StreamEmpty( is ) )
    {
      std::string s;
      if( std::getline( is, s ) )
      {
        // Get rid of anything past the comment
        std::size_t pos = s.find_first_of( sComment );
        if( pos != std::string::npos )
          s.resize( pos );
        ReadLine( m, s );
      }
    }
    return m;
  }
  static Map Read( const std::string &Filename, const char *ErrorMsgType = nullptr,
                   const std::string &sComment = Hash )
  {
    if( FileExists( Filename ) )
    {
      std::ifstream s( Filename );
      if( !s.bad() )
        return Read( s, sComment );
    }
    std::ostringstream os;
    os << ( ErrorMsgType ? ErrorMsgType : "Key/value file" ) << " \"" << Filename << "\" not found";
    throw std::runtime_error( os.str().c_str() );
  }
};

// Make the ancestor directories leading up to last element
// NB: Don't test whether this worked, so that it if this is done simultaneously
// by many threads, it will still work
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

// Wrapper for posix gethostname()
std::string GetHostName();

template<typename T> struct BootRep
{
  static constexpr int idxCentral{ -1 };
  Vector<T> Central;
  Matrix<T> Replica;
  int extent() const { return Central.size; }
  int size() const { return Replica.size1; }
  void resize( int Extent, int NumBoot ) { Central.resize( Extent ); Replica.resize( NumBoot, Extent ); }
  BootRep() = default;
  BootRep( int Extent, int NumBoot ) { resize( Extent, NumBoot ); }
  void Read( ::H5::Group &g, const std::string &Name );
  void Write( ::H5::Group &g, const std::string &Name ) const;
  inline T * operator[]( int idx )
  {
    if( Central.size == 0 )
      throw std::runtime_error( "Can't access empty BootRep" );
    if( idx == -1 )
      return Central.data;
    if( Replica.size2 != Central.size )
      throw std::runtime_error( "BoootRep extent mismatch" );
    if( idx >= Replica.size1 )
      throw std::runtime_error( "BoootRep index " + std::to_string( idx ) + " >= NumBoot "
                               + std::to_string( Replica.size1 ) );
    return Replica.data + Replica.tda * idx;
  }
  inline void MapView( Vector<T> &v, int idx )
  {
    if( idx == idxCentral )
      v.MapView( Central );
    else
      v.MapRow( Replica, idx );
  }
  inline void MapView( VectorView<T> &v, int idx ) const
  {
    if( idx == idxCentral )
      v = Central;
    else
      v.MapRow( Replica, idx );
  }
  inline void MapColumn( VectorView<T> &v, int Column ) const
  {
    v.MapColumn( Replica, Column );
  }
};

// Attributes for filenames in form base.type.seed.ext
struct FileNameAtt
{
  std::string Filename; // Full (and unadulterated) original filename
  std::string NameNoExt;// Original filename, without path and without extension
  std::string Dir;      // Directory part of the filename (with trailing '/')
  std::string Base;     // Base of the filename (containing momentum, current, delta T and timeslice)
  std::vector<std::string> BaseShortParts;// BaseShort split into parts at separator '_'
  std::string Type;
  std::string SeedString;
  bool        bSeedNum = false;
  SeedType    Seed = 0; // Numeric value of SeedString, but only if bSeedNum == true
  std::string Ext;
  std::vector<int> op; // These indices only match opNames if no global opName vector used
  std::vector<std::string> opNames; // These will always be present
  std::vector<std::string> Extra;
  bool bGotTimeslice = false;
  int         Timeslice = 0;
  using MomentumMap = std::map<std::string, FileNameMomentum, LessCaseInsensitive>;
  MomentumMap p;
  bool bGotDeltaT = false;
  int         DeltaT;
  std::vector<Gamma::Algebra> Gamma; // Gamma algebras extracted from name
  // SpecDir and Spectator are either both set or neither is
  std::string SpecDir;
  std::string Spectator;
  bool bSpectatorGotSuffix = false;
  std::vector<std::string> Meson; // 0 is source, 1 is sink
  std::vector<std::string> MesonMom; // With momenta
  // reset contents
  void clear()
  {
    Filename.clear();
    NameNoExt.clear();
    Dir.clear();
    Base.clear();
    BaseShortParts.clear();
    Type.clear();
    SeedString.clear();
    bSeedNum = false;
    Seed = 0;
    Ext.clear();
    op.clear();
    opNames.clear();
    Extra.clear();
    bGotTimeslice = false;
    p.clear();
    bGotDeltaT = false;
    Gamma.clear();
    SpecDir.clear();
    Spectator.clear();
    bSpectatorGotSuffix = false;
    Meson.clear();
    MesonMom.clear();
  }
  void swap( FileNameAtt &o )
  {
    std::swap( Filename, o.Filename );
    std::swap( NameNoExt, o.NameNoExt );
    std::swap( Dir, o.Dir );
    std::swap( Base, o.Base );
    std::swap( BaseShortParts, o.BaseShortParts );
    std::swap( Type, o.Type );
    std::swap( SeedString, o.SeedString );
    std::swap( bSeedNum, o.bSeedNum );
    std::swap( Seed, o.Seed );
    std::swap( Ext, o.Ext );
    std::swap( op, o.op );
    std::swap( opNames, o.opNames );
    std::swap( Extra, o.Extra );
    std::swap( bGotTimeslice, o.bGotTimeslice );
    std::swap( Timeslice, o.Timeslice );
    std::swap( p, o.p );
    std::swap( bGotDeltaT, o.bGotDeltaT );
    std::swap( DeltaT, o.DeltaT );
    std::swap( Gamma, o.Gamma );
    std::swap( SpecDir, o.SpecDir );
    std::swap( Spectator, o.Spectator );
    std::swap( bSpectatorGotSuffix, o.bSpectatorGotSuffix );
    std::swap( Meson, o.Meson );
    std::swap( MesonMom, o.MesonMom );
  }
  friend void swap( FileNameAtt &l, FileNameAtt &r )
  {
    l.swap( r );
  }
  void Parse( const std::string &Filename, std::vector<std::string> * pOpNames = nullptr,
              const std::vector<std::string> * pIgnoreMomenta = nullptr,
              const std::vector<std::string> * pIgnoreRegEx = nullptr,
              bool bPreBootstrap = false );
  FileNameAtt() = default;
  explicit FileNameAtt( const std::string &Filename, std::vector<std::string> * pOpNames = nullptr,
                        const std::vector<std::string> * pIgnoreMomenta = nullptr,
                        const std::vector<std::string> * pIgnoreRegEx = nullptr,
                        bool bPreBootstrap = false )
    { Parse( Filename, pOpNames, pIgnoreMomenta, pIgnoreRegEx, bPreBootstrap ); }
  // Append the extra info to the string
  void AppendExtra( std::string &s, int Last = 0, int First = INT_MAX ) const;
  void AppendOps( std::string &s, const std::string &FirstSeparator,
                  std::vector<std::string> * pOpNames = nullptr ) const;
  void AppendOps( std::string &s, std::vector<std::string> * pOpNames = nullptr ) const
  { AppendOps( s, Extra.size() ? Period : Underscore, pOpNames ); }
  std::string GetBaseShort( int First = 0, int Last = UINT_MAX ) const;
  std::string GetBaseShortExtra( int First = 0, int Last = UINT_MAX ) const;
  std::string GetBaseExtra( int Last = 0, int First = -1 ) const;
  // Make a new name based on this one, overriding specified elements
  std::string DerivedName( const std::string &Suffix, const std::string &Snk, const std::string &Src,
                           const std::string &Ext ) const;
  const FileNameMomentum &GetMomentum( const std::string &Name = Momentum::DefaultPrefix ) const;
  bool HasNonZeroMomentum() const;
  const FileNameMomentum &GetFirstNonZeroMomentum() const;
  void AppendMomentum( std::string &s, const FileNameMomentum &fnp, const std::string &Name ) const;
  void AppendMomentum( std::string &s, const FileNameMomentum &fnp )const{AppendMomentum( s, fnp, fnp.Name );}
  std::string MakeMesonName( const std::string &Quark ) const
  { return ::Common::MakeMesonName( Quark, Spectator ); }
protected:
  std::vector<std::string> ParseOpNames( int NumOps );
};

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext);

// If present, remove Token from a string. Return true if removed
bool ExtractToken( std::string &Prefix, const std::string &Token );

// If present, remove integer preceded by Token from a string
void ExtractInteger( std::string &Prefix, bool &bHasValue, int &Value, const std::string Token );

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &s, bool &bHasTimeslice, int & Timeslice );

// Strip out DeltaT from a string if present
void ExtractDeltaT( std::string &Prefix, bool &bHasDeltaT, int &DeltaT );

// Append DeltaT to string
void AppendDeltaT( std::string &s, int DeltaT );

// Remove any gammas from Prefix
std::vector<Gamma::Algebra> ExtractGamma( std::string &Prefix );

// Append Gamma to string
void AppendGamma( std::string &s, Gamma::Algebra g, const std::string &Sep = Underscore );

// Append Gamma to string
inline void AppendGammaDeltaT( std::string &s, Gamma::Algebra g, int DeltaT )
{
  AppendGamma( s, g );
  AppendDeltaT( s, DeltaT );
}

// Remove any gammas from Prefix
void ReplaceGamma( std::string &Prefix, Gamma::Algebra gFrom, Gamma::Algebra gTo );

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

// Parse a string containing a selection from sOption, returning count of how many times specified
inline std::vector<int> ParseOptions( const char * pStringStart, std::size_t const Len, const std::string &sOption,
                                      bool bCaseInsensitive = true, bool bNoRepeats = true, bool bIgnoreSpace = true )
{
  std::vector<int> Count( sOption.length(), 0 );
  for( std::size_t i = 0; i < Len; ++i )
  {
    char c = pStringStart[i];
    if( !bIgnoreSpace || !std::isspace( c ) )
    {
      int idx = 0;
      while( idx < sOption.length() && Compare( c, sOption[idx] ) )
        ++idx;
      if( idx >= sOption.length() || ( Count[idx]++ && bNoRepeats ) )
      {
        const std::string s{ pStringStart, Len };
        throw std::runtime_error( "Option string '" + s + "' doesn't match \"" + sOption + "\"" );
      }
    }
  }
  return Count;
}

struct FoldProp
{
  RPS rps;
  bool t0Abs = true;
  bool Conjugate = false;
  bool Parse( const char * const pc, std::size_t const Len )
  {
    static const std::string sOption{ "RIEOPN0C" };
    std::vector<int> iCount{ ParseOptions( pc, Len, sOption ) };
    // Apply rules, i.e. first three pairs of options conflict
    bool bOK = true;
    switch( iCount[0] + iCount[1] )
    {
      case 0:
        rps.reality = Common::Reality::Unknown;
        break;
      case 1:
        rps.reality = iCount[0] ? Common::Reality::Real : Common::Reality::Imag;
        break;
      default:
        bOK = false;
    }
    switch( iCount[2] + iCount[3] )
    {
      case 0:
        rps.parity = Common::Parity::Unknown;
        break;
      case 1:
        rps.parity = iCount[2] ? Common::Parity::Even : Common::Parity::Odd;
        break;
      default:
        bOK = false;
    }
    switch( iCount[4] + iCount[5] )
    {
      case 0:
        rps.sign = Common::Sign::Unknown;
        break;
      case 1:
        rps.sign = iCount[4] ? Common::Sign::Positive : Common::Sign::Negative;
        break;
      default:
        bOK = false;
    }
    t0Abs = !iCount[6];
    Conjugate = iCount[7];
    // Chuck a wobbly (sorry, lapsed into Australian there for a minute!)
    if( !bOK )
    {
      const std::string s{ pc, Len };
      throw std::runtime_error( "Conflicting options '" + s + "'" );
    }
    bool bRealImagOnly = true;
    for( std::size_t i = 2; bRealImagOnly && i < sOption.length(); ++i )
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

// Traits for floating point types: float; double
template<typename ST> struct SampleTraits<ST, typename std::enable_if<std::is_floating_point<ST>::value>::type> : public std::true_type
{
  using scalar_type = ST;
  using value_type = ST;
  static constexpr bool is_complex = false;
  static constexpr int scalar_count = 1;
  static inline scalar_type * ScalarPtr( ST * p ) { return p; };
  static inline const scalar_type * ScalarPtr( const ST * p ) { return p; };
};

// Traits for std::complex<ST> types
template<typename ST> struct SampleTraits<std::complex<ST>, typename std::enable_if<SampleTraits<ST>::value>::type> : public std::true_type
{
  using scalar_type = typename SampleTraits<ST>::scalar_type;
  using value_type = std::complex<ST>;
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
      s << sep;
      Common::ValWithEr<T>::Header( n, s, sep );
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
  using value_type = typename Traits::value_type;
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
  void swap( CorrelatorFile &o )
  {
    using std::swap;
    swap( NumSnk_, o.NumSnk_ );
    swap( NumSrc_, o.NumSrc_ );
    swap( Nt_, o.Nt_ );
    swap( m_pData, o.m_pData );
    swap( AlgSnk_, o.AlgSnk_ );
    swap( AlgSrc_, o.AlgSrc_ );
    Name_.swap( o.Name_ );
    swap( bHasTimeslice, o.bHasTimeslice );
    swap( Timeslice_, o.Timeslice_ );
  }
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
             std::string *pGroupName = nullptr, const std::vector<bool> *pAlgSnkNeg = nullptr,
             const char * pDSName = nullptr );
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

template <typename T>
inline void swap( CorrelatorFile<T> &l, CorrelatorFile<T> &r )
{
  l.swap( r );
}

using CorrelatorFileC = CorrelatorFile<std::complex<double>>;
using CorrelatorFileD = CorrelatorFile<double>;

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void CorrelatorFile<T>::Read(const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk,
                             std::vector<Gamma::Algebra> &AlgSrc, const int * pTimeslice,
                             const char * PrintPrefix, std::string *pGroupName,
                             const std::vector<bool> *pAlgSnkNeg, const char * pDSName )
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
      ::H5::DataSet ds = g.openDataSet( ( pDSName && *pDSName ) ? pDSName : "correlator" );
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

// Generate random numbers
void GenerateRandom( std::vector<fint> &Random, SeedType Seed, std::size_t NumBoot, std::size_t NumSamples );

enum class SampleSource { Binned, Raw, Bootstrap };

std::istream& operator>>( std::istream& is, SampleSource &sampleSource );
std::ostream& operator<<( std::ostream& os, SampleSource sampleSource );

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
  static constexpr bool is_complex { Traits::is_complex };
  static constexpr int scalar_count { Traits::scalar_count };
  static constexpr int idxCentral{ BootRep<T>::idxCentral };
protected:
  int NumSamples_;
  int Nt_;
  std::unique_ptr<T[]> m_pData;
  std::unique_ptr<fint[]> m_pRandNum; // Random numbers used to generate the bootstrap sample
  std::vector<std::string> AuxNames;
  int NumSamplesRaw_;
  int NumSamplesBinned_;
  std::unique_ptr<T[]> m_pDataRaw;
  std::unique_ptr<T[]> m_pDataBinned;
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
  void Bootstrap();
public:
  FileNameAtt Name_;
  SeedType Seed_ = 0; // This is only ever set by the bootstrap code
  std::string SeedMachine_; // name of the machine that ran the bootstrap
  int binSize = 1;
  int SampleSize = 0; // Number of samples (after binning) used to create bootstrap (set during bootstrap)
  std::vector<Common::ConfigCount> ConfigCount; // Info on every config in the bootstrap in order
  std::vector<std::string> FileList; // Info on every config in the bootstrap in order
  bool bRawBootstrap = false; // Indicates the raw samples are bootstraps
  inline int NumSamples() const { return NumSamples_; }
  inline int GetNumExtraSamples() const { return NumExtraSamples; }
  inline int NumSamplesRaw() const { return NumSamplesRaw_; }
  inline int NumSamplesBinned() const { return NumSamplesBinned_; }
  inline int NumSamples( SampleSource ss ) const
  {
    return ss == SampleSource::Binned ? NumSamplesBinned()
         : ss == SampleSource::Bootstrap ? NumSamples()
         : NumSamplesRaw();
  }
  inline int Nt() const { return Nt_; }
  inline const fint * RandNum() const { return m_pRandNum.get(); };
  // Return the summary data for the nth type (0=central, 1=bias, etc)
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
  // Return the central replica summary data for the nth column
  inline ValWithEr<scalar_type> &getSummaryData( const std::string &ColumnName )
  {
    int idxField{ IndexIgnoreCase( GetColumnNames(), ColumnName ) };
    if( idxField >= GetColumnNames().size() )
      throw std::runtime_error( "Column " + ColumnName + " not available" );
    return getSummaryData()[idxField];
  }
  inline const ValWithEr<scalar_type> &getSummaryData( const std::string &ColumnName ) const
  {
    int idxField{ IndexIgnoreCase( GetColumnNames(), ColumnName ) };
    if( idxField >= GetColumnNames().size() )
      throw std::runtime_error( "Column " + ColumnName + " not available" );
    return getSummaryData()[idxField];
  }
  void WriteSummaryData( std::ostream &s, int idx = 0 ) const
  {
    s << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
    const ValWithEr<scalar_type> * p{ getSummaryData( idx ) };
    for( int t = 0; t < Nt_; t++ )
      s << ( t == 0 ? "" : " " ) << *p++;
  }
  const std::vector<std::string> &GetSummaryNames() const { return SummaryNames; };
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
      if( t )
        s << Sep;
      Common::ValWithEr<T>::Header( *ps, s, Sep );
    }
  }
  const std::vector<std::string> & GetColumnNames() const { return ColumnNames; }
  inline int GetColumnIndexNoThrow( const std::string & ColumnName ) const
  {
    assert( ColumnNames.size() <= std::numeric_limits<int>::max() && "Wow that's a lot of columns!" );
    for( int i = 0; i < ColumnNames.size(); i++ )
      if( EqualIgnoreCase( ColumnNames[i], ColumnName ) )
        return i;
    return -1;
  }
  inline int GetColumnIndexNoThrow( const std::string & ColumnName, int idx ) const
  {
    return GetColumnIndexNoThrow( ColumnName + std::to_string( idx ) );
  }
  inline int GetColumnIndex( const std::string & ColumnName ) const
  {
    int i = GetColumnIndexNoThrow( ColumnName );
    if( i < 0 )
      throw std::runtime_error( "Column " + ColumnName + " not found" );
    return i;
  }
  inline int GetColumnIndex( const std::string & ColumnName, int idx ) const
  {
    return GetColumnIndex( ColumnName + std::to_string( idx ) );
  }
  void resize( int NumSamples, int Nt, std::vector<std::string> * pAuxNames = nullptr )
  {
    const int NewNumExtraSamples{ static_cast<int>( pAuxNames && pAuxNames->size() ? pAuxNames->size() : 1 ) };
    if( Nt_ != Nt )
    {
      AllocSummaryBuffer( SummaryNames.size() * Nt );
      ColumnNames.clear();
      m_pDataRaw.reset( nullptr );
      NumSamplesRaw_ = 0;
      m_pDataBinned.reset( nullptr );
      NumSamplesBinned_ = 0;
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
  inline T * getRaw() { return m_pDataRaw.get(); }
  const inline T * getRaw() const { return m_pDataRaw.get(); }
  T * resizeRaw( int NumSamplesRaw )
  {
    if( !NumSamplesRaw )
      m_pDataRaw.reset( nullptr );
    else if( !Nt_ )
      throw std::runtime_error( "Cannot resizeRaw() when Nt==0" );
    else if( NumSamplesRaw_ != NumSamplesRaw )
    {
      m_pDataRaw.reset( new T[ static_cast<std::size_t>( NumSamplesRaw ) * Nt_ ] );
      NumSamplesRaw_ = NumSamplesRaw;
    }
    bRawBootstrap = false;
    return m_pDataRaw.get();
  }
  inline T * getBinned() { return m_pDataBinned.get(); }
  const inline T * getBinned() const { return m_pDataBinned.get(); }
  inline T * get( SampleSource ss, int idx = 0 ) { return ss == SampleSource::Bootstrap
    ? (*this)[idx] : ( ss == SampleSource::Raw ? getRaw() : getBinned() ) + Nt() * idx; }
  const inline T * get( SampleSource ss, int idx = 0 ) const { return ss == SampleSource::Bootstrap
    ? (*this)[idx] : ( ss == SampleSource::Raw ? getRaw() : getBinned() ) + Nt() * idx; }
  T * resizeBinned( int NumSamplesBinned )
  {
    if( !NumSamplesBinned )
      m_pDataBinned.reset( nullptr );
    else if( !Nt_ )
      throw std::runtime_error( "Cannot resizeBinned() when Nt==0" );
    else if( NumSamplesBinned_ != NumSamplesBinned )
      m_pDataBinned.reset( new T[ static_cast<std::size_t>( NumSamplesBinned ) * Nt_ ] );
    NumSamplesBinned_ = NumSamplesBinned;
    return m_pDataBinned.get();
  }
  bool IsFinite() { return Common::IsFinite( reinterpret_cast<scalar_type *>( m_pData.get() ),
      static_cast<size_t>( ( NumSamples_ + NumExtraSamples ) * ( SampleTraits<T>::is_complex ? 2 : 1 ) ) * Nt_ ); }
  template <typename U>
  void IsCompatible( const Sample<U> &o, int * pNumSamples = nullptr, unsigned int CompareFlags = COMPAT_DEFAULT ) const;
  template <typename U> void CopyAttributes( const Sample<U> &o );
  void Bin( SampleSource ssTo = SampleSource::Binned ); // Auto bin based on config count
  void Bin( int binSize_, SampleSource ssTo = SampleSource::Binned );
  void BootstrapBinned( int NumBoot, SampleSource ssFrom );
  void Bootstrap( SeedType Seed, const std::string * pMachineName = nullptr );
  template <typename U> void Bootstrap( const Sample<U> &Random );
  virtual void SetName( const std::string &FileName, std::vector<std::string> * pOpNames = nullptr )
  {
    Name_.Parse( FileName, pOpNames );
  }
  void SetName( FileNameAtt &&FileName ) { Name_ = std::move( FileName ); }
  void Read( const char *PrintPrefix = nullptr, std::string * pGroupName = nullptr );
  void Read( const std::string &FileName, const char *PrintPrefix = nullptr, std::vector<std::string> * pOpNames = nullptr,
             std::string * pGroupName = nullptr )
  {
    SetName( FileName, pOpNames );
    Read( PrintPrefix, pGroupName );
  }
  void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void MakeCorrSummary( const char * pAvgName );
  void WriteSummary( const std::string &sOutFileName, bool bVerboseSummary = false );
  explicit Sample(int NumSamples = 0, int Nt = 0, std::vector<std::string> * pAuxNames = nullptr,
                  std::vector<std::string> * pSummaryNames = nullptr)
    : NumSamples_{0}, Nt_{0}, NumExtraSamples{0}, NumSamplesRaw_{0}, NumSamplesBinned_{0}
  {
    resize( NumSamples, Nt, pAuxNames );
    if( pSummaryNames )
      SetSummaryNames( *pSummaryNames );
  }
  Sample(const std::string &FileName, const char *PrintPrefix = nullptr, std::vector<std::string> * pOpNames = nullptr,
         std::string * pGroupName = nullptr)
  : NumSamples_{0}, Nt_{0}, NumExtraSamples{0}, NumSamplesRaw_{0}, NumSamplesBinned_{0}
  {
    Read( FileName, PrintPrefix, pOpNames, pGroupName );
  }
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
  void MakeMean( int idxSlot = idxCentral, bool bBinned = false )
  {
    ZeroSlot( idxSlot );
    const int NumRows{ bBinned ? NumSamplesBinned_ : NumSamples_ };
    if( NumRows && Nt_ )
    {
      T * const dst{ (*this)[idxSlot] };
      const T * src{ bBinned ? getBinned() : (*this)[0] };
      for( int i = 0; i < NumRows; i++ )
        for( int t = 0; t < Nt_; t++ )
          dst[t] += *src++;
      for( int t = 0; t < Nt_; t++ )
        dst[t] /= NumRows;
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
        throw std::runtime_error( "Complex sample has imaginary component" );
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
        s << CommaSpace << cc;
      s << NewLine;
    }
    s << "# NumColumns: " << Nt();
    if( ColumnNames.empty() )
      s << " (rows t0 ... t" << ( Nt() - 1 ) << ")";
    else
    {
      s << NewLine << "# " << sColumnNames << ": ";
      for( std::size_t i = 0; i < ColumnNames.size(); ++i )
      {
        if( i )
          s << CommaSpace;
        s << ColumnNames[i];
      }
    }
    s << NewLine;
    if( bVerboseSummary )
    {
      s << "# FileCount: " << FileList.size() << NewLine;
      std::size_t i = 0;
      for( const std::string &f : FileList )
        s << "# File " << i++ << ": " << f << NewLine;
    }
  }
  virtual void SummaryColumnNames( std::ostream &os ) const
  {
    if( ColumnNames.empty() )
    {
      // First column is timeslice, then each of the summary columns per timeslice
      os << "t";
      for( const std::string &Name : SummaryNames )
      {
        os << Space;
        ValWithEr<T>::Header( Name, os, Space );
      }
    }
    else
    {
      // Show the bootstrap central replica for every column on one row
      bool bFirst{ true };
      for( const std::string &Name : ColumnNames )
      {
        if( bFirst )
          bFirst = false;
        else
          os << Space;
        ValWithEr<T>::Header( Name, os, Space );
      }
    }
  }
  virtual void SummaryContents( std::ostream &os ) const
  {
    os << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
    const ValWithEr<scalar_type> * p{ m_pSummaryData.get() };
    if( ColumnNames.empty() )
    {
      // One row per timeslice, showing each of the summary values
      for( int t = 0; t < Nt_; t++, p++ )
      {
        os << t;
        for( std::size_t f = 0; f < SummaryNames.size(); f++ )
          os << Space << p[f * Nt_];
        os << NewLine;
      }
    }
    else
    {
      // A single row, showing the central value for each column
      for( int t = 0; t < Nt_; t++ )
      {
        if( t )
          os << Space;
        os << *p++;
      }
      os << NewLine;
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
void Sample<T>::IsCompatible( const Sample<U> &o, int * pNumSamples, unsigned int CompareFlags ) const
{
  static const std::string sPrefix{ "Incompatible " + sBootstrap + " samples - " };
  std::string sSuffix{ ":\n  " + Name_.Filename + "\nvs " + o.Name_.Filename };
  if( !( CompareFlags & COMPAT_DISABLE_BASE ) && !Common::EqualIgnoreCase( o.Name_.Base, Name_.Base ) )
    throw std::runtime_error( sPrefix + "base " + o.Name_.Base + sNE + Name_.Base + sSuffix );
  if( !( CompareFlags & COMPAT_DISABLE_TYPE ) && !Common::EqualIgnoreCase( o.Name_.Type, Name_.Type ) )
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
  if( !( CompareFlags & COMPAT_DISABLE_NT ) && o.Nt() != Nt_ )
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
      /*if( !( CompareFlags & COMPAT_DISABLE_CONFIG_COUNT ) && r.Count != l.Count )
        throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                 ", NumTimeslices " + std::to_string(r.Count) + sNE +
                                 std::to_string(l.Count) + sSuffix );*/
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

// Auto bin: all measurements on each config binned into a single measurement
template <typename T> void Sample<T>::Bin( SampleSource ssTo )
{
  if( !NumSamplesRaw_ )
    throw std::runtime_error( "No raw data to bin" );
  if( bRawBootstrap )
    throw std::runtime_error( "Can't bin bootstrapped data" );
  int NumSamplesRaw{ 0 };
  int NewBinSize{ ConfigCount[0].Count };
  for( const Common::ConfigCount &cc : ConfigCount )
  {
    NumSamplesRaw += cc.Count;
    if( NewBinSize && NewBinSize != cc.Count )
      NewBinSize = 0; // Indicates the bin size varies per ConfigCount
  }
  if( NumSamplesRaw != NumSamplesRaw_ )
    throw std::runtime_error( "ConfigCount doesn't match raw data" );
  const T * pSrc = getRaw();
  std::unique_ptr<T[]> DeleteWhenOutOfScope;
  T * pDst;
  if( ssTo == SampleSource::Binned )
  {
    pDst = resizeBinned( static_cast<int>( ConfigCount.size() ) );
    binSize = NewBinSize;
  }
  else if( ssTo == SampleSource::Raw )
  {
    DeleteWhenOutOfScope = std::move( m_pDataRaw );
    NumSamplesRaw_ = 0;
    pDst = resizeRaw( static_cast<int>( ConfigCount.size() ) );
  }
  else
    throw std::runtime_error( "Can only store binned data in Binned or Raw samples" );
  for( const Common::ConfigCount &cc : ConfigCount )
  {
    for( int t = 0; t < Nt_; ++t )
      pDst[t] = *pSrc++;
    for( int i = 1; i < cc.Count; ++i )
      for( int t = 0; t < Nt_; ++t )
        pDst[t] += *pSrc++;
    for( int t = 0; t < Nt_; ++t )
      pDst[t] /= cc.Count;
    pDst += Nt_;
  }
}

template <typename T> void Sample<T>::Bin( int binSize_, SampleSource ssTo )
{
  if( binSize_ < 1 )
    throw std::runtime_error( "BinSize " + std::to_string( binSize_ ) + " invalid" );
  if( !NumSamplesRaw_ )
    throw std::runtime_error( "No raw data to bin" );
  if( bRawBootstrap )
    throw std::runtime_error( "Can't bin bootstrapped data" );
  int NumSamplesRaw{ 0 };
  for( const Common::ConfigCount &cc : ConfigCount )
  {
    if( cc.Count != ConfigCount[0].Count )
      throw std::runtime_error( "Can't manually bin - raw sample counts differ per config" );
    NumSamplesRaw += cc.Count;
  }
  if( NumSamplesRaw != NumSamplesRaw_ )
    throw std::runtime_error( "ConfigCount doesn't match raw data" );
  const bool PartialLastBin = NumSamplesRaw % binSize_;
  const int NewNumSamples{ static_cast<int>( NumSamplesRaw / binSize_ + ( PartialLastBin ? 1 : 0 ) ) };
  const T * pSrc = getRaw();
  std::unique_ptr<T[]> DeleteWhenOutOfScope;
  T * pDst;
  if( ssTo == SampleSource::Binned )
  {
    pDst = resizeBinned( NewNumSamples );
    binSize = binSize_;
  }
  else if( ssTo == SampleSource::Raw )
  {
    DeleteWhenOutOfScope = std::move( m_pDataRaw );
    NumSamplesRaw_ = 0;
    pDst = resizeRaw( NewNumSamples );
  }
  else
    throw std::runtime_error( "Can only store binned data in Binned or Raw samples" );
  for( int Bin = 0; Bin < NewNumSamples; ++Bin )
  {
    for( int t = 0; t < Nt_; ++t )
      pDst[t] = *pSrc++;
    const int ThisBinSize{ PartialLastBin && Bin==(NewNumSamples-1) ? NumSamplesRaw % binSize_ : binSize_ };
    for( int i = 1; i < ThisBinSize; ++i )
      for( int t = 0; t < Nt_; ++t )
        pDst[t] += *pSrc++;
    for( int t = 0; t < Nt_; ++t )
      pDst[t] /= ThisBinSize;
    pDst += Nt_;
  }
}

/**
 Perform bootstrap
 Compute the mean on this sample, saving this as the central value
 Returned sample has specified extra fields ... but only the central value is copied from original
 */
template <typename T> void Sample<T>::Bootstrap()
{
  if( bRawBootstrap )
    throw std::runtime_error( "Raw data already bootstrapped" );
  SampleSize = NumSamplesBinned_;
  MakeMean( idxCentral, true ); // Compute the mean from the binned data
  // Master thread can choose random samples while other threads bootstrap
  const fint * pRandom{ m_pRandNum.get() };
  T * dst{ (*this)[0] };
  for( int i = 0; i < NumSamples_; ++i, dst += Nt_ )
  {
    for( int s = 0; s < NumSamplesBinned_; ++s )
    {
      const T * src{ getBinned() + static_cast<std::ptrdiff_t>( *pRandom++ ) * Nt_ };
      for( int t = 0; t < Nt_; t++ )
      {
        if( !s )
          dst[t]  = *src++;
        else
          dst[t] += *src++;
      }
    }
    // Turn the sum into an average
    if( NumSamplesBinned_ > 1 )
      for( int t = 0; t < Nt_; ++t)
        dst[t] /= NumSamplesBinned_;
  }
}

// Bootstrap either the raw or binned data. Always store result in raw, marked as bootstrapped
template <typename T> void Sample<T>::BootstrapBinned( int NumBoot, SampleSource ssFrom )
{
  if( NumBoot < 2 )
    throw std::invalid_argument( "Bootstrap size " + std::to_string( NumBoot ) + " invalid" );
  if( ssFrom != SampleSource::Raw && ssFrom != SampleSource::Binned )
    throw std::runtime_error( "Can only bootstrap from Binned or Raw samples" );
  if( ssFrom == SampleSource::Raw && bRawBootstrap )
    throw std::runtime_error( "Raw data already bootstrapped" );
  const int NumSamplesFrom{ NumSamples( ssFrom ) };
  const T * src{ get( ssFrom ) };
  if( !Nt_ || !NumSamplesFrom || !src )
    throw std::runtime_error( "No data to bootstrap" );
  std::unique_ptr<T[]> DeleteWhenOutOfScope;
  if( ssFrom == SampleSource::Raw )
  {
    DeleteWhenOutOfScope = std::move( m_pDataRaw );
    NumSamplesRaw_ = 0;
  }
  T * dst{ resizeRaw( static_cast<std::size_t>( NumBoot ) * Nt_ ) };
  bRawBootstrap = true;
  std::mt19937                        engine( Seed_ );
  std::uniform_int_distribution<fint> random( 0, NumSamplesFrom - 1 );
  const T InvNumBinned{ static_cast<T>( 1 ) / static_cast<T>( NumSamplesFrom ) };
  for( int b = 0; b < NumBoot; ++b, dst += Nt_ )
  {
    for( int s = 0; s < NumSamplesFrom; ++s )
    {
      const T * ThisRow{ src + static_cast<std::ptrdiff_t>( random( engine ) ) * Nt_ };
      for( int t = 0; t < Nt_; t++ )
      {
        if( !s )
          dst[t]  = *ThisRow++;
        else
          dst[t] += *ThisRow++;
      }
    }
    // Turn the sum into an average
    if( NumSamplesFrom > 1 )
      for( int t = 0; t < Nt_; ++t)
        dst[t] *= InvNumBinned;
  }
}

// Generate a new set of random numbers to use for the seed
template <typename T> void Sample<T>::Bootstrap( SeedType Seed, const std::string * pMachineName )
{
  const std::size_t RandomLen{ static_cast<std::size_t>( NumSamples_ ) * static_cast<std::size_t>( NumSamplesBinned_ ) };
  assert( RandomLen && "Can't bootstrap nothing" );
  Seed_ = Seed;
  if( pMachineName && ! pMachineName->empty() )
    SeedMachine_ = *pMachineName;
  else
    SeedMachine_ = GetHostName();
  // Generate the random numbers
  fint * pRandom{ new fint[ RandomLen ] };
  m_pRandNum.reset( pRandom );
  std::mt19937                        engine( Seed );
  std::uniform_int_distribution<fint> random( 0, NumSamplesBinned_ - 1 );
  for( std::size_t i = 0; i < RandomLen; i++ )
    pRandom[i] = random( engine );
  Bootstrap();
}

template <typename T> template <typename U> void Sample<T>::Bootstrap( const Sample<U> &Random )
{
  const std::size_t RandomLen{ static_cast<std::size_t>( NumSamples_ ) * static_cast<std::size_t>( NumSamplesBinned_ ) };
  assert( RandomLen && "Can't bootstrap nothing" );
  const fint * pSrc{ Random.m_pRandNum.get() };
  if( !pSrc )
    throw std::runtime_error( "No random numbers in " + Random.Name_.Filename );
  if( Random.SampleSize != NumSamplesBinned_ )
    throw std::runtime_error( "Random numbers are for SampleSize " + std::to_string(Random.SampleSize)
                             + " not " + std::to_string( NumSamplesBinned_ ) );
  if( Random.NumSamples_ < NumSamples_ )
    throw std::runtime_error( "Not enough Random numbers for " + std::to_string( NumSamples_ )
                             + " samples. Only " + std::to_string( Random.NumSamples_ ) + " available" );
  Seed_ = Random.Seed_;
  SeedMachine_ = Random.SeedMachine_;
  // Copy the random numbers
  fint * pRandom{ new fint[ RandomLen ] };
  m_pRandNum.reset( pRandom );
  for( std::size_t i = 0; i < RandomLen; i++ )
    pRandom[i] = pSrc[i];
  Bootstrap();
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
  SummaryColumnNames( s );
  s << NewLine;
  SummaryContents( s );
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
void Sample<T>::Read( const char *PrintPrefix, std::string *pGroupName )
{
  if( !Name_.bSeedNum )
    throw std::runtime_error( "Seed missing from " + Name_.Filename );
  if( Name_.Type.empty() )
    throw std::runtime_error( "Type missing from " + Name_.Filename );
  ::H5::H5File f;
  ::H5::Group  g;
  H5::OpenFileGroup( f, g, Name_.Filename, PrintPrefix, pGroupName );
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
      std::int8_t i8RawBootstrap{ 0 };
      try
      {
        a = g.openAttribute( sRawBootstrap );
        a.read( ::H5::PredType::NATIVE_INT8, &i8RawBootstrap );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      std::vector<std::string> myFileList;
      // Try to load the FileList from dataSet first, falling back to attribute
      bool bGotFileList{ false };
      try
      {
        ::H5::DataSet ds = g.openDataSet( sFileList );
        myFileList = H5::ReadStrings( ds );
        ds.close();
        bGotFileList = true;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      if( !bGotFileList )
      {
        try
        {
          a = g.openAttribute(sFileList);
          myFileList = H5::ReadStrings( a );
          a.close();
          bGotFileList = true;
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
        }
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
                  try // to read the random numbers
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
                  try // to read the raw data
                  {
                    dsp.close();
                    ds.close();
                    ds = g.openDataSet( "samples_Raw" );
                    dsp = ds.getSpace();
                    if( dsp.getSimpleExtentNdims() == 2 )
                    {
                      dsp.getSimpleExtentDims( Dim );
                      if( Dim[0] < std::numeric_limits<int>::max() && Dim[1] == static_cast<unsigned int>( Nt_ ) )
                      {
                        ds.read( resizeRaw( static_cast<int>( Dim[0] ) ), H5::Equiv<T>::Type );
                      }
                    }
                  }
                  catch(const ::H5::Exception &)
                  {
                    ::H5::Exception::clearErrorStack();
                    resizeRaw( 0 );
                  }
                  // resizeRaw() resets this, so had to delay setting bRawBootstrap 'til now
                  bRawBootstrap = ( i8RawBootstrap != 0 );
                  try // to read the binned data
                  {
                    dsp.close();
                    ds.close();
                    ds = g.openDataSet( "samples_Binned" );
                    dsp = ds.getSpace();
                    if( dsp.getSimpleExtentNdims() == 2 )
                    {
                      dsp.getSimpleExtentDims( Dim );
                      if( Dim[0] < std::numeric_limits<int>::max() && Dim[1] == static_cast<unsigned int>( Nt_ ) )
                      {
                        ds.read( resizeBinned( static_cast<int>( Dim[0] ) ), H5::Equiv<T>::Type );
                      }
                    }
                  }
                  catch(const ::H5::Exception &)
                  {
                    ::H5::Exception::clearErrorStack();
                    resizeBinned( 0 );
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
                        const ::H5::DataType ErrType{ ds.getDataType() };
                        if( ErrType == H5::Equiv<ValWithErOldV1<scalar_type>>::Type )
                        {
                          const std::size_t Count{ SummaryNames.size() * Nt_ };
                          std::vector<ValWithErOldV1<scalar_type>> tmp( Count );
                          ds.read( &tmp[0], ErrType );
                          ValWithEr<scalar_type> * dst{ m_pSummaryData.get() };
                          for( std::size_t i = 0; i < Count; ++i )
                          {
                            dst[i].Central = tmp[i].Central;
                            dst[i].Low     = tmp[i].Low;
                            dst[i].High    = tmp[i].High;
                            dst[i].Check   = tmp[i].Check;
                            dst[i].Min     = 0;
                            dst[i].Max     = 0;
                          }
                          bOK = true;
                        }
                        else
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
  if( bFiniteError || !bOK )
  {
    resize( 0, 0 ); // We can use NumSamples() == 0 or Nt() == 0 to tell whether the file is open
    if( bFiniteError )
      throw std::runtime_error( "NANs in " + Name_.Filename );
    //if( !bOK )
      throw std::runtime_error( "Unable to read sample from " + Name_.Filename );
  }
}

template <typename T>
void Sample<T>::Write( const std::string &FileName, const char * pszGroupName )
{
  SetName( FileName );
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
    const std::int8_t i8{ 1 };
    if( bRawBootstrap )
    {
      a = g.createAttribute( sRawBootstrap, ::H5::PredType::STD_U8LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT8, &i8 );
      a.close();
      NumAttributes++;
    }
    NumAttributes += WriteAttributes( g );
    if( FileList.size() )
      H5::WriteStringData( g, sFileList, FileList );
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
    if( getRaw() )
    {
      Dims[0] = NumSamplesRaw_;
      Dims[1] = Nt_;
      dsp = ::H5::DataSpace( 2, Dims );
      ds = g.createDataSet( "samples_Raw", H5::Equiv<T>::Type, dsp );
      ds.write( getRaw(), H5::Equiv<T>::Type );
      ds.close();
      dsp.close();
    }
    if( getBinned() )
    {
      Dims[0] = NumSamplesBinned_;
      Dims[1] = Nt_;
      dsp = ::H5::DataSpace( 2, Dims );
      ds = g.createDataSet( "samples_Binned", H5::Equiv<T>::Type, dsp );
      ds.write( getBinned(), H5::Equiv<T>::Type );
      ds.close();
      dsp.close();
    }
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
      // The values in this table go from 0 ... NumSamples_ + 1, so choose a space-minimising size in the file
      ds = g.createDataSet( sRandom, SampleSize<=static_cast<int>(std::numeric_limits<std::uint8_t>::max())
        ? ::H5::PredType::STD_U8LE : SampleSize<=static_cast<int>(std::numeric_limits<std::uint16_t>::max())
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
  virtual void SummaryComments( std::ostream & s, bool bVerboseSummary = false ) const
  {
    Base::SummaryComments( s, bVerboseSummary );
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

template <typename T> std::size_t GetExtent( const std::vector<std::vector<T>> &v )
{
  std::size_t Extent{ 0 };
  for( const std::vector<T> & i : v )
    Extent += i.size();
  return Extent;
}

template <typename T>
struct Model : public Sample<T>
{
  using Base = Sample<T>;
  static const std::string EnergyPrefix;
  int dof = 0;
  bool CovarFrozen = false;
  using SS = SampleSource;
  SS CovarSource = SS::Binned;  // Where am I building the covariance from
  std::vector<int> CovarRebin;  // Did I rebin the data before computing covariances
  // 0 = build correlation from CovarSource, then scale by the variance of the data on each bootstrap replica
  // else this is the number of bootstrap replicas to use when estimating covariance
  int CovarNumBoot = 0;
  int CovarSampleSize = 0;      // What is the appropriate parameter for m to use in the T^2 distribution
  Vector<T> Guess;
  Matrix<T> CovarIn;      // Optional. From correlation source
  Matrix<T> Covar;        // As used in the fit
  Matrix<T> Correl;       // As used in the fit
  Matrix<T> CorrelCholesky;
  Matrix<T> CovarInv;
  Matrix<T> CorrelInv;
  Matrix<T> CorrelInvCholesky;
  Matrix<T> CovarInvCholesky;
  BootRep<T> StdErrorMean; // From all samples of binned data
  BootRep<T> FitInput;
  BootRep<T> ModelPrediction;
  BootRep<T> ErrorScaled;  // (Theory - Data) * CholeskyScale

  // These variables were added in current format
  Params params;
  std::vector<std::vector<int>> FitTimes;

  // Helper functions
  int GetExtent() { return static_cast<int>( ::Common::GetExtent( FitTimes ) ); };
  int NumFitParams() { return static_cast<int>( params.NumScalars( Param::Type::Variable ) ); };
  int NumParams() { return static_cast<int>( params.NumScalars( Param::Type::All ) ); };

  Model() : Base::Sample{} {}
  Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns );
  Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns,
         int CovarSampleSize_, bool CovarFrozen_, SampleSource CovarSource_,
         std::vector<int> CovarRebin_, int CovarNumBoot_ );
  const std::string & DefaultGroupName() override { return sModel; }
  inline bool NewParamsMorePrecise( bool covarFrozen_, int NumSamples ) const
  {
    // Unfreezing the covariance matrix makes a big difference!
    if( ( CovarFrozen && !covarFrozen_ ) || ( !CovarFrozen && covarFrozen_ ) )
      return CovarFrozen;
    return NumSamples > Base::NumSamples_;
  }
  /*void Read( const char *PrintPrefix = nullptr, std::string *pGroupName =nullptr )
  {
    Base::Read( PrintPrefix, pGroupName );
  }
  void Read( const std::string &FileName, const char *PrintPrefix = nullptr, std::string *pGroupName =nullptr )
  {
    SetName( FileName );
    Base::Read( PrintPrefix, pGroupName );
  }*/
  void ReadAttributes( ::H5::Group &g ) override;
  void ValidateAttributes() override;
  int WriteAttributes( ::H5::Group &g ) override;
  void SummaryComments( std::ostream & s, bool bVerboseSummary = false ) const override;
  void SummaryColumnNames( std::ostream &os ) const override;
  void SummaryContents( std::ostream &os ) const override;
  void WriteSummaryTD( const std::string &sOutFileName, bool bVerboseSummary = false );
protected:
  int OldFormatNumExponents;
  std::vector<std::string> OldFormatOpNames;
  void ReorderOldFormat( int NumOps, int NumExponents, std::unique_ptr<T[]> &pData, int Num );
  void CommonConstruct( const std::vector<std::string> &ExtraColumns );
};

struct ConstantSource
{
  std::size_t File; // Index of the constant file in the DataSet this comes from
  const Param &param;
  ConstantSource( std::size_t File_, const Param &param_ ) : File{File_}, param{param_} {}
};

template <typename T>
struct DataSet
{
  using SS = SampleSource;
  using ConstMap = std::map<Param::Key, ConstantSource, Param::Key::Less>;
  struct FixedParam
  {
    ConstantSource  src;   // Where to get values from
    std::size_t     Count; // How many to get
    int             idx;   // Where to put them (destination)
    FixedParam(const ConstantSource &src_, std::size_t Count_, int idx_):src{src_}, Count{Count_}, idx{idx_}{}
  };
  int NSamples;   // Number of samples we are using. These are guaranteed to exist
  int MaxSamples; // Maximum number of samples available. Guaranteed to exist. >= NSamples.
  int Extent = 0; // Number of data points in our fit (i.e. total number of elements in FitTimes)
  std::vector<Fold<T>>          corr;     // Correlator files
  std::vector<std::vector<int>> FitTimes; // The actual timeslices we are fitting to in each correlator
  ConstMap constMap;
  std::vector<int> RebinSize; // Bin sizes if I've rebinned the raw data
  int OriginalBinSize;        // Bin size before any rebinning (to deal with covariance matrices)
  // Cached data for selected FitTimes
  Vector<T> vCentral;
  Matrix<T> mRaw;
  Matrix<T> mBinned;
  Matrix<T> mBoot;
  const Matrix<T> &Cache( SampleSource ss ) const
  {
    const Matrix<T> &m{ ss == SampleSource::Bootstrap ? mBoot : ss == SampleSource::Binned ? mBinned : mRaw };
    if( m.size1 == 0 || m.size2 == 0 )
    {
      std::stringstream s;
      s << "Cached " << ss << " data unavailable";
      throw std::runtime_error( s.str().c_str() );
    }
    return m;
  }
protected:
  std::array<std::vector<fint>, 2> RandomCentralBuf; // No replacement. Saves having to rebuild continually
  std::array< MatrixView<fint>, 2> RandomViews;
  std::array<std::vector<fint>, 2> RandomBuffer;
  std::vector<Model<T>> constFile;// Each of the constant files (i.e. results from previous fits) I've loaded
  void AddConstant( const Param::Key &Key, std::size_t File, const Param &param );
  void SetValidatedFitTimes( std::vector<std::vector<int>> &&FitTimes );
  void MakeCentralReplicaNonRandom();
public:
  explicit DataSet( int nSamples = 0 ) : NSamples{ nSamples } {}
  inline bool empty() const { return corr.empty() && constFile.empty(); }
  void clear();
  int  LoadCorrelator( Common::FileNameAtt &&FileAtt, unsigned int CompareFlags = COMPAT_DEFAULT,
                       const char * PrintPrefix = "  " );
  void LoadModel     ( Common::FileNameAtt &&FileAtt, const std::string &Args );
  std::vector<std::string> GetModelFilenames() const;
  void SortOpNames( std::vector<std::string> &OpNames );
  void Rebin( const std::vector<int> &RebinSize );
  int NumSamplesBinned() const; // Minimum number of binned samples supported by all correlators
  void SetFitTimes( const std::vector<std::vector<int>> &FitTimes ); // A list of all the timeslices to include
  void SetFitTimes( int tMin, int tMax ); // All fit ranges are the same
  void GetData( int idx, VectorView<T> &vResult, SS ss = SS::Bootstrap ) const
  {
    if( idx == Fold<T>::idxCentral )
      vResult = vCentral;
    else
      vResult.MapRow( Cache( ss ), idx );
  }
  void GetFixed( int idx, Vector<T> &vResult, const std::vector<FixedParam> &Params ) const;
  bool InitRandomNumbers( SS ss ); // Make sure random numbers needed to build covariances are present
  void InitRandomNumbers() { InitRandomNumbers( SS::Raw ); InitRandomNumbers( SS::Binned ); }
  const std::vector<fint> &RandomCentral( SS ss ) const
  {
    const std::vector<fint> &r{ RandomCentralBuf[ ss == SS::Raw ? 0 : 1 ] };
    if( !( ss == SS::Binned || ss == SS::Raw ) || r.empty() )
    {
      std::stringstream os;
      os << "Central replica " << ss << " (non-) random numbers unavailable";
      throw std::runtime_error( os.str().c_str() );
    }
    return r;
  }
  const MatrixView<fint> &RandomNumbers( SS ss ) const
  {
    const MatrixView<fint> &r{ RandomViews[ ss == SS::Raw ? 0 : 1 ] };
    if( !( ss == SS::Binned || ss == SS::Raw ) || !r.data() )
    {
      std::stringstream os;
      os << ss << " random numbers unavailable";
      throw std::runtime_error( os.str().c_str() );
    }
    return r;
  }
  // All of these functions make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
  // From what is already a bootstrap replica, and therefore only the central replica
  void MakeCovarFromBootstrap( SS ss, Matrix<T> &Covar ) const;
  void MakeVarFromBootstrap( SS ss, Vector<T> &Var ) const;
  // From the underlying raw/(re)binned data (without bootstrapping)
  void MakeCovarFromNonBootstrap( int idx, SS ss, Matrix<T> &Covar ) const;
  void MakeVarFromNonBootstrap( int idx, SS ss, Vector<T> &Var ) const;
  // Calls the appropriate one of the previous two functions, depending on what the data are
  void MakeCovariance( int idx, SS ss, Matrix<T> &Covar ) const;
  void MakeVariance( int idx, SS ss, Vector<T> &Covar ) const;
  // Make covariance matrix using a secondary bootstrap - Random numbers must be provided by caller
  // I.e. unfrozen covariance matrix
  void MakeCovariance( int idx, SS ss, Matrix<T> &Covar, const MatrixView<fint> &Random)const;
  void SaveMatrixFile( const Matrix<T> &m, const std::string &Type, const std::string &FileName,
                       std::vector<std::string> &Abbreviations,
                       const std::vector<std::string> *FileComments = nullptr,
                       const char *pGnuplotExtra = nullptr ) const;
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
  static bool IsValuePresent( const char * & p );
  
public:
  void Parse( int argc, const char *argv[], const std::vector<SwitchDef> &defs );

  CommandLine() = default;
  CommandLine( int argc, const char *argv[], const std::vector<SwitchDef> &defs )
  { Parse( argc, argv, defs ); }

  inline bool GotSwitch( const std::string &SwitchName ) const {
    return Switches.find( SwitchName ) != Switches.end(); }

  inline int NumValues( const std::string &Switch ) const {
    int iNumValues{ 0 };
    SwitchMap::const_iterator it = Switches.find( Switch );
    if( it != Switches.end() ) {
      iNumValues = static_cast<int>( it->second.size() );
    }
    return iNumValues;
  }
  
  inline const std::vector<std::string> & SwitchStrings( const std::string &Switch ) const
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

  template<typename T> inline T SwitchValue( const std::string &Switch, int Subscript = 0 ) const
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
