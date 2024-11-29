/**

 Mike's lattice QCD utilities
 
 Source file: Utility.hpp
 
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

// MLU utilities (no dependencies other than c++ stdlib)

#ifndef MLU_Utility_hpp
#define MLU_Utility_hpp

#include <MLU/MLUFirst.hpp>
#include <MLU/FitRange.hpp>
#include <MLU/GSLVecMat.hpp>
#include <MLU/HDF5.hpp>
#include <MLU/Param.hpp>
#include <MLU/Posix.hpp>

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

// Default output file extension for text based summaries of binary data
#ifndef TEXT_EXT
#define TEXT_EXT "txt"
#endif
#ifndef DAT_EXT
#define DAT_EXT "dat"
#endif

extern "C" const char * MLUVersionInfoHuman();

BEGIN_MLU_NAMESPACE

// Compatibility flags for filename comparisons

static constexpr unsigned int COMPAT_DEFAULT{ 0 };
static constexpr unsigned int COMPAT_DISABLE_BASE{ 1 };
static constexpr unsigned int COMPAT_DISABLE_NT{ 2 };
static constexpr unsigned int COMPAT_DISABLE_TYPE{ 4 };
static constexpr unsigned int COMPAT_DISABLE_ENSEMBLE{ 8 };

template <typename T> int sgn( T x )
{
  return (T(0) < x) - (x < T(0));
}

template <typename T> inline void ApplyNegateStar( T &Value, NegateStar ns )
{
  if( ns == NegateStar::Negate || ns == NegateStar::NegateStar )
    Value = -Value;
  if( ns == NegateStar::Star || ns == NegateStar::NegateStar )
    Value = Conjugate( Value ); // Does nothing for real numbers
}

// Text required for summaries of correlators
namespace CorrSumm {
  extern const char sep[];
  extern const char Comment[];
  static constexpr int NumFields{ 3 };
  extern const char * FieldNames[NumFields];
};

extern const char pszCorrelGnuplot[];

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
extern const std::string sStdErrorMean;
extern const std::string sErrorScaled;
extern const std::string sCovariance;
extern const std::string sCovarianceIn;
extern const std::string sCovarianceInv;
extern const std::string sCovarianceInvCholesky;
extern const std::string sCorrelation;
extern const std::string sCorrelationInv;
extern const std::string sCorrelationCholesky;
extern const std::string sCorrelationInvCholesky;
extern const std::string sCorrelationParam;
extern const std::string sCorrelationParamNames;
extern const std::string sFitInput;
extern const std::string sModelPrediction;
extern const std::string sOperators;
extern const std::string sSummaryDSName;
extern const std::string sSummaryNames;
extern const std::string sColumnNames;
extern const std::string sSampleSize;
extern const std::string sEnsemble;
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
extern const std::string sModelTypeList;
extern const std::string sModelArgsList;
extern const std::string sNE;
extern const std::vector<std::string> sCorrSummaryNames;
extern const std::string sChiSqPerDof;
extern const std::string sPValue;
extern const std::string sPValueH;

// String containing switch name to enable alternate overlap normalisation
// (see definition of ModelOverlap in Fit/ModelCommon.hpp)
extern const std::string sOverlapAltNorm;

extern const std::vector<std::string> DefaultModelStats;

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
      MLU::Trim( v );
    return bOK;
  }
  template<class V = T> static typename std::enable_if<!std::is_same<V, std::string>::value, bool>::type
  GetValue( std::istringstream &iss, V &v ) { return iss >> v && MLU::StreamEmpty( iss ); }
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
    if( !MLU::StreamEmpty( iss ) )
    {
      Key k;
      if( !( iss >> k ) )
        throw std::runtime_error( "KeyValReader::ReadLine() bad key " + s );
      T t{};
      if( !bOptional || !MLU::StreamEmpty( iss ) )
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
    if( MLU::FileExists( Filename ) )
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

struct MomentumMap : std::map<std::string, Momentum, LessCaseInsensitive>
{
  void Parse( std::string &ExtractFrom );
  MomentumMap() = default;
  MomentumMap( std::string &ExtractFrom ) { Parse( ExtractFrom ); }
};
using MomentumPair = MomentumMap::value_type;

struct NamedMomentum
{
  std::string Name;
  Momentum p;
  NamedMomentum( const std::string &Name_ = Momentum::DefaultPrefix ) : Name(Name_) {}
  NamedMomentum( const std::string &Name_, const Momentum &p_ ) : Name{Name_}, p(p_) {}
  NamedMomentum( const MomentumPair &p ) : NamedMomentum( p.first, p.second ) {}
};

/*
 Attributes for filenames in form
    [Dir/]base1[_base2[...]][[...].extra1].extra0.[[...]op1_]op0.type[.seed].ext
 or
    [Dir/]base1[_base2[...]]_op1_op0.type[.seed].ext
 */
struct FileNameAtt
{
  std::string Filename; // Full (and unadulterated) original filename
  std::string NameNoExt;// Original filename, with path but without extension
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
  MomentumMap p;
  std::vector<int> DeltaT;
  inline bool GotDeltaT() const { return !DeltaT.empty(); }
  std::vector<Gamma::Algebra> Gamma; // Gamma algebras extracted from name
  // SpecDir and Spectator are either both set or neither is
  std::string SpecDir;
  std::string Spectator;
  bool bSpectatorGotSuffix = false;
  std::vector<std::string> Quark; // 0 is source, 1 is sink
  std::vector<std::string> Meson; // 0 is source, 1 is sink
  std::vector<std::string> MesonMom; // With momenta - if available
  std::vector<NamedMomentum> MesonP; // The momenta - if available
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
    DeltaT.clear();
    Gamma.clear();
    SpecDir.clear();
    Spectator.clear();
    bSpectatorGotSuffix = false;
    Quark.clear();
    Meson.clear();
    MesonMom.clear();
    MesonP.clear();
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
    std::swap( DeltaT, o.DeltaT );
    std::swap( Gamma, o.Gamma );
    std::swap( SpecDir, o.SpecDir );
    std::swap( Spectator, o.Spectator );
    std::swap( bSpectatorGotSuffix, o.bSpectatorGotSuffix );
    std::swap( Quark, o.Quark );
    std::swap( Meson, o.Meson );
    std::swap( MesonMom, o.MesonMom );
  }
  friend void swap( FileNameAtt &l, FileNameAtt &r )
  {
    l.swap( r );
  }
  void Parse( const std::string &DirBase, const std::string &Type, bool bHasSeed, SeedType Seed,              const std::string &Ext, std::vector<std::string> * pOpNames = nullptr,
              const std::vector<std::string> * pIgnoreMomenta = nullptr,
              const std::vector<std::string> * pIgnoreRegEx = nullptr,
              bool bPreBootstrap = false );
  void Parse( const std::string &DirBase, const std::string &Type, SeedType Seed,
              const std::string &Ext, std::vector<std::string> * pOpNames = nullptr,
              const std::vector<std::string> * pIgnoreMomenta = nullptr,
              const std::vector<std::string> * pIgnoreRegEx = nullptr,
              bool bPreBootstrap = false )
  {
    Parse( DirBase, Type, true, Seed, Ext, pOpNames, pIgnoreMomenta, pIgnoreRegEx, bPreBootstrap );
  }
  void Parse( const std::string &DirBase, const std::string &Type,
              const std::string &Ext, std::vector<std::string> * pOpNames = nullptr,
              const std::vector<std::string> * pIgnoreMomenta = nullptr,
              const std::vector<std::string> * pIgnoreRegEx = nullptr,
              bool bPreBootstrap = false )
  {
    Parse( DirBase, Type, false, 0, Ext, pOpNames, pIgnoreMomenta, pIgnoreRegEx, bPreBootstrap );
  }
  void Parse( const std::string &Filename, std::vector<std::string> * pOpNames = nullptr,
              const std::vector<std::string> * pIgnoreMomenta = nullptr,
              const std::vector<std::string> * pIgnoreRegEx = nullptr,
              bool bPreBootstrap = false );
protected:
  void Parse( std::vector<std::string> * pOpNames,
              const std::vector<std::string> * pIgnoreMomenta,
              const std::vector<std::string> * pIgnoreRegEx,
              bool bPreBootstrap );
public:
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
  /// Get the specified parts of the base, i.e. from first to last. Negative counts are relative to the end (-1=last, etc)
  std::string GetBaseShort( int First = 0, int Last = UINT_MAX ) const;
  /// As per GetBaseShort() but append all of the extra segments as well
  std::string GetBaseShortExtra( int First = 0, int Last = UINT_MAX ) const;
  /// Get all of the base, then append specified parts of Extra. Negative counts are relative to the end (-1=last, etc)
  std::string GetBaseExtra( int Last = 0, int First = INT_MAX ) const;
  /** Get complete base name of the file (i.e. everything before the type)
   without path, but with optional prefix (eg some other path) **/
  std::string GetBase( const std::string &Prefix = "" ) const;
  /// Get complete base name of the file (i.e. everything before the type) with path
  std::string GetBasePath() const { return GetBase( Dir ); }
  /// Get an alternate version of this name without path, but with optional prefix (eg some other path)
  std::string GetAlt( const std::string &Type, const std::string &Ext,
                      const std::string &Prefix = "" ) const;
  std::string GetAltPath( const std::string &Type, const std::string &Ext ) const
  { return GetAlt( Type, Ext, Dir ); }
  // Make a new name based on this one, overriding specified elements
  std::string DerivedName( const std::string &Suffix, const std::string &Snk, const std::string &Src,
                           const std::string &Ext ) const;
  const MomentumPair &GetMomentum( const std::string &Name = Momentum::DefaultPrefix ) const;
  bool HasNonZeroMomentum() const;
  const MomentumPair &GetFirstNonZeroMomentum() const;
  void AppendMomentum( std::string &s, const Momentum &p, const std::string &Name ) const
  { s.append( p.FileString( Name ) ); }
  void AppendMomentum( std::string &s, const MomentumPair &np ) const
  { AppendMomentum( s, np.second, np.first );}
  std::string MakeMesonName( const std::string &Quark ) const
  { return ::MLU::ConcatHeavyFirst( Quark, Spectator ); }
protected:
  std::vector<std::string> ParseOpNames( int NumOps );
};

/// Make a filename "Base(.* optional * Type)(.* optional * seed)(.* optional * Ext)"
std::string MakeFilename( const std::string &Base, const std::string &Type,
                          bool bHasSeed, SeedType Seed, const std::string &Ext );
/// Make a filename "Base(.* optional * Type).seed(.* optional * Ext)"
inline std::string MakeFilename( const std::string &Base, const std::string &Type,
                                 SeedType Seed, const std::string &Ext )
{
  return MakeFilename( Base, Type, true, Seed, Ext );
}
/// Make a filename "Base(.* optional * Type)(.* optional * Ext)" (no seed)
inline std::string MakeFilename( const std::string &Base, const std::string &Type,
                                 const std::string &Ext )
{
  return MakeFilename( Base, Type, false, 0, Ext );
}
/// Find existing file with specified Seed, but fall back to another seed if missing
std::string ExistingFilename( const std::string &Base, const std::string &Type, SeedType Seed,
                              const std::string &Ext );
/// Find existing file with specified Seed, but fall back to another seed if missing
std::string ExistingAnySeed( const std::string &sFilename );
/// Override Seed in filename to specified seed (if it exists)
std::string PreferSeed( const std::string &sFilename, SeedType NewSeed );
std::string PreferSeed( const std::string &sFilename ); // Uses default seed

// If present, remove Token from a string. Return true if removed
bool ExtractToken( std::string &Prefix, const std::string &Token );

// If present, remove integer preceded by Token from a string
void ExtractInteger( std::string &Prefix, bool &bHasValue, int &Value, const std::string &Token );

// If present, remove integer preceded by Token from a string
std::vector<int> ExtractIntegers( std::string &Prefix, const std::string &Token );

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &s, bool &bHasTimeslice, int & Timeslice );

// Strip out DeltaT from a string if present
std::vector<int> ExtractDeltaT( std::string &Prefix );

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
  if( l.reality == MLU::Reality::Imag && r.reality == MLU::Reality::Imag )
    rps.sign = rps.sign * MLU::Sign::Negative;
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
  /**
   Parse options describing how to fold
   - Returns:
   true if only real or imaginary specified, i.e. no instructions on how to fold.
   false if any of the other options specified ... which describe how to fold the forward and backward waves.
   */
  bool Parse( const char * const pc, std::size_t const Len )
  {
    static const std::string sOption{ "RIEOPN0C" };
    std::vector<int> iCount{ ParseOptions( pc, Len, sOption ) };
    // Apply rules, i.e. first three pairs of options conflict
    bool bOK = true;
    switch( iCount[0] + iCount[1] )
    {
      case 0:
        rps.reality = MLU::Reality::Unknown;
        break;
      case 1:
        rps.reality = iCount[0] ? MLU::Reality::Real : MLU::Reality::Imag;
        break;
      default:
        bOK = false;
    }
    switch( iCount[2] + iCount[3] )
    {
      case 0:
        rps.parity = MLU::Parity::Unknown;
        break;
      case 1:
        rps.parity = iCount[2] ? MLU::Parity::Even : MLU::Parity::Odd;
        break;
      default:
        bOK = false;
    }
    switch( iCount[4] + iCount[5] )
    {
      case 0:
        rps.sign = MLU::Sign::Unknown;
        break;
      case 1:
        rps.sign = iCount[4] ? MLU::Sign::Positive : MLU::Sign::Negative;
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
template<typename T> struct SampleTraits : public std::false_type{};

// Traits for float
template<> struct SampleTraits<float> : public std::true_type
{
  using scalar_type = float;
  using value_type = float;
  static constexpr bool is_complex = false;
  static constexpr int scalar_count = 1;
  static inline       scalar_type * ScalarPtr(       float * p ) { return p; };
  static inline const scalar_type * ScalarPtr( const float * p ) { return p; };
  static inline scalar_type Real( const float &v ) { return v; }
  static inline scalar_type Imag( const float &v ) {throw std::runtime_error("SampleTraits::Imag");}
  static inline scalar_type RealImag( const float &v, int i )
  {
    if( i )
      throw std::runtime_error( "SampleTraits<float>::Imag" );
    return v;
  }
};

// Traits for double
template<> struct SampleTraits<double> : public std::true_type
{
  using scalar_type = double;
  using value_type = double;
  static constexpr bool is_complex = false;
  static constexpr int scalar_count = 1;
  static inline       scalar_type * ScalarPtr(       double * p ) { return p; };
  static inline const scalar_type * ScalarPtr( const double * p ) { return p; };
  static inline scalar_type Real( const double &v ) { return v; }
  static inline scalar_type Imag( const double &v ) {throw std::runtime_error("SampleTraits::Imag");}
  static inline scalar_type RealImag( const double &v, int i )
  {
    if( i )
      throw std::runtime_error( "SampleTraits<double>::Imag" );
    return v;
  }
};

// Traits for std::complex<T> types
template<typename ST> struct SampleTraits<std::complex<ST>> : public std::true_type
{
  using scalar_type = ST;
  using value_type = std::complex<ST>;
  static constexpr bool is_complex = true;
  static constexpr int scalar_count = 2;
  static inline       scalar_type * ScalarPtr(       value_type * p )
      { return reinterpret_cast<      scalar_type *>( p ); };
  static inline const scalar_type * ScalarPtr( const value_type * p )
      { return reinterpret_cast<const scalar_type *>( p ); };
  static inline scalar_type Real( const value_type &v ) { return v.real(); }
  static inline scalar_type Imag( const value_type &v ) { return v.imag(); }
  static inline scalar_type RealImag( const value_type &v, int i )
  {
    if( i == 0 )
      return v.real();
    if( i == 1 )
      return v.imag();
    throw std::runtime_error( "SampleTraits<std::complex<T>>::RealImag " + std::to_string( i ) );
  }
};

template <typename T>
void SummaryHeader( std::ostream &s, const std::string & sOutFileName,
                    const char * pszHeaderComments = nullptr );

/// Make a summary of the data
/// Assume there is a central replica followed by nSample copies (i.e. nSample + 1 total)
template <typename T>
void SummaryHelper( const std::string & sOutFileName, const T * pData, const int nt,
                    const int nSample = 0, const char * pszHeaderComments = nullptr,
                    bool bFolded = false );

enum class SampleSource { Binned, Raw, Bootstrap };

std::istream& operator>>( std::istream& is, SampleSource &sampleSource );
std::ostream& operator<<( std::ostream& os, SampleSource sampleSource );

template <typename T> std::size_t GetExtent( const std::vector<std::vector<T>> &v )
{
  std::size_t Extent{ 0 };
  for( const std::vector<T> & i : v )
    Extent += i.size();
  return Extent;
}

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

END_MLU_NAMESPACE
#endif // MLU_Utility_hpp
