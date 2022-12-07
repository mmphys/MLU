/*************************************************************************************
 
 Create Ratios, e.g. R1, R2 and associated values such as Z_V
 Source file: Ratio.hpp
 Copyright (C) 2020-2022
 Author: Michael Marshall<Michael.Marshall@ed.ac.uk>
 
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

#pragma once
#include <MLU/Common.hpp>

//static const std::vector<std::string> vHeavy{ "h0", "h1", "h2", "h3" };
//static const int NumHeavy{ static_cast<int>( vHeavy.size() ) };
//static const std::vector<int> vDeltaT{ 12, 14, 16, 20 };
//static const int NumDeltaT{ static_cast<int>( vDeltaT.size() ) };
using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Model = Common::Model<Scalar>;
using Vector = Common::Vector<Scalar>;
using VectorView = Common::VectorView<Scalar>;

extern Model DummyModel;
extern const std::string DefaultColumnName;
extern const char LoadFilePrefix[];

enum class Freeze
{
  None,
  Central,
  Constant
};

inline void Append( std::string &s, const std::string &Add )
{
  s.append( Common::Underscore );
  s.append( Add );
}

// A quark with a sepecified momentum
struct QP
{
  std::string q;
  Common::FileNameMomentum p;
  QP( const std::string &Q, const Common::FileNameMomentum &P ) : q{Q}, p{P} {}
};

// Case insensitive compare of QP
struct LessQP
{
  bool operator()( const QP &lhs, const QP &rhs ) const
  {
    int i = Common::CompareIgnoreCase( lhs.q, rhs.q );
    if( i )
      return i < 0;
    return lhs.p < rhs.p;
  }
};

// Map from QP -> ModelInfo
/*struct QPMapModelInfo : public std::map<QP, ModelInfo, LessQP>
{
  int MaxModelSamples;
  Common::SeedType Seed;
  QPMapModelInfo( const std::string &Filename,
                  const std::string &modelBase, const std::string &C2Base );
};*/

// A quark with a sepecified DeltaT
struct QDT
{
  std::string q;
  int deltaT = 0;
  std::string opSnk;
  std::string opSrc;
  QDT() {}
  QDT( std::string q_, int deltaT_, const std::string &opSnk_, const std::string &opSrc_ )
  : q{q_}, deltaT{deltaT_}, opSnk{opSnk_}, opSrc{opSrc_} {}
};

// Case insensitive compare of QP
struct LessQDT
{
  bool operator()( const QDT &lhs, const QDT &rhs ) const
  {
    int i = Common::CompareIgnoreCase( lhs.q, rhs.q );
    if( i )
      return i < 0;
    if( lhs.deltaT != rhs.deltaT )
      return lhs.deltaT < rhs.deltaT;
    i = Common::CompareIgnoreCase( lhs.opSnk, rhs.opSnk );
    if( i )
      return i < 0;
    i = Common::CompareIgnoreCase( lhs.opSrc, rhs.opSrc );
    return i < 0;
  }
};

// Map from a DeltaT -> ModelInfo
/*struct QDTMapModelInfo : public std::map<QDT, ModelInfo, LessQDT>
{
  int MaxModelSamples;
  Common::SeedType Seed;
  QDTMapModelInfo( const std::string &FitListName, const std::string &modelBase, const std::string &C2Base );
};*/

// File cache
// BEWARE: Loading files, i.e. GetIndex(), can invalidate FileT references (by moving Files)
//         DO NOT keep FileT references around across GetIndex() calls
template<typename FileT>
struct FileCache
{
  static constexpr int BadIndex{ -1 };
  const char * PrintPrefix; // non-null to print filenames as they are loaded
protected:
  std::string Base;
  using TNameMap = std::map<std::string, int>;
  TNameMap NameMap;
  std::vector<FileT> Files;
public:
  void SetBase( const std::string &Base_ ) { Base = Base_; }
  std::vector<std::string> opNames; // These are operator names referred to by Files
  FileCache( const char * Prefix = nullptr ) : PrintPrefix{Prefix} {}
  void clear();
  // Add file to cache if not present - delay loading until accessed
  int GetIndex( const std::string &Filename, const char * printPrefix = nullptr );
  FileT &operator()( int iIndex, const char * printPrefix ); // This is what loads the file
  inline FileT &operator[]( int iIndex ) { return (*this)( iIndex, PrintPrefix ); };
  inline FileT &operator[]( const std::string &Filename ) { return (*this)[GetIndex(Filename)]; }
  inline const FileT &operator[]( const std::string &Filename ) const { return (*this)[GetIndex(Filename)]; }
};

// Allows a quark name to be read in as a string
// Extracts the momentum from the file name
struct QuarkReader : public std::string
{
  QP operator()( const std::string &Filename ) const;
};

// Static list of mappings from key to file. Can't be added to once read.
// Two-step construction (i.e. constructor then Read) so virtual functions in place
template<typename Key, typename LessKey, typename KeyRead = Key, typename LessKeyRead = LessKey,
         typename M = Model>
struct KeyFileCache
{
//  const std::string C2Base;
  Freeze freeze;
protected:
  using KeyMapT = std::map<Key, int, LessKey>;
  using FileCacheT = FileCache<M>;
  KeyMapT KeyMap;
  FileCacheT model;
  using FieldMapT = std::map<std::string, std::string>;
  using ConstantMapT = std::map<std::string, Scalar>;
  FieldMapT FieldMap;
  ConstantMapT ConstantMap;
protected:
  using KeyReadMapT = std::map<Key, std::string, LessKey>;
  template<typename K = Key, typename R = KeyRead>
  typename std::enable_if<std::is_same<K, R>::value, KeyReadMapT>::type
  ReadNameMap( const std::string &Filename, const char *pErrorMsg );
  template<typename K = Key, typename R = KeyRead>
  typename std::enable_if<!std::is_same<K, R>::value, KeyReadMapT>::type
  ReadNameMap( const std::string &Filename, const char *pErrorMsg );
  int GetIndex( const Key &k ) const; // Only call this of we know we're not frozen to a constant
  virtual void FrozenOptions( std::string &sOptions ) {}
public:
  static constexpr int BadIndex{ FileCacheT::BadIndex };
  virtual void clear();
  KeyFileCache() { clear(); }
  virtual ~KeyFileCache() {}
  // Use this to load a list of fit files
  void Read( const std::string &Filename, const char * PrintPrefix );
  M &operator()( const Key &key, const char * PrintPrefix ); // Make sure model loaded and print message
  inline M &operator[]( const Key &key ) { return operator()( key, nullptr ); }
  Vector GetVector( const Key &key, const std::string &ColumnName = DefaultColumnName );
};

using QDTModelMap = KeyFileCache<QDT, LessQDT>;

struct QPModelMap : public KeyFileCache<QP, LessQP, QuarkReader, Common::LessCaseInsensitive>
{
  using Base = KeyFileCache<QP, LessQP, QuarkReader, Common::LessCaseInsensitive>;
  std::string Spectator;
  std::string Get2ptName( const QP &key );
  void clear() override;
protected:
  void FrozenOptions( std::string &sOptions ) override;
};

// Base class for ratio construction
class Maker
{
public:
  const int MaxSamples;
  //const std::string inBase;
  //const std::string C2Base;
  const std::string modelBase;
  const std::string outBase;
  //const Common::Momentum p;
  const bool bSymmetrise;
protected:
  const bool bOverlapAltNorm;
  QPModelMap EFit;
  struct CorrT
  {
    std::string Name;
    int Handle = -1;
    Fold * Corr = nullptr;
  };
public:
  using CorrCache = FileCache<Fold>;
  CorrCache Cache2;
  CorrCache Cache3;
protected:
  std::string HeavyKey( const std::string &Heavy ) const;
public:
  Maker( const Common::CommandLine &cl );
  virtual ~Maker() {}
  virtual void Make( std::string &sFileName ) = 0;
};

class ZVRCommon : public Maker
{
protected:
  struct SSInfo
  {
    const std::string &op;
    const QP qp;
    Model &EModel;
    const Vector E;
    SSInfo( QPModelMap &efit, const std::string &op_, QP qp_ )
    : op{op_}, qp{qp_}, EModel{efit(qp,LoadFilePrefix)}, E{efit.GetVector(qp)} {}
  };
  inline void AppendOp( std::string &s, const std::string &Op ) { Append( s, Op ); }
  inline void AppendOps( std::string &s, const std::string &Snk, const std::string &Src)
  {
    AppendOp( s, Snk );
    AppendOp( s, Src );
  }
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const SSInfo &Snk, const SSInfo &Src ) = 0;
public:
  void Make( std::string &sFileName ) override;
  ZVRCommon( const Common::CommandLine &cl ) : Maker( cl ) {}
};

// Make Z_V: eq 2.7 pg 3 https://arxiv.org/pdf/1305.7217.pdf
class ZVMaker : public ZVRCommon
{
protected:
  void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
             const SSInfo &Snk, const SSInfo &Src ) override;
public:
  ZVMaker( const std::string &TypeParams, const Common::CommandLine &cl );
};

// Make R1 and R2: eq 2.8 pg 4 https://arxiv.org/pdf/1305.7217.pdf
class RMaker : public ZVRCommon
{
protected:
  const bool bAltR3;
  QDTModelMap ZVmi;
  void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
             const SSInfo &Snk, const SSInfo &Src ) override;
public:
  RMaker( const std::string &TypeParams, const Common::CommandLine &cl );
};

// Make F parralel, F perpendicular, F+ and F0
class FFitConstMaker : public Maker
{
  const std::string i3Base;
  const unsigned int N;
  const Scalar ap;
protected:
  bool bAdjustGammaSpatial;
  int RatioNum;
  std::string qSnk;
  std::string qSrc;
  Common::FileNameMomentum p;
  std::string Prefix;
  std::string Suffix;
  std::string FitType;
  std::vector<int> FitParts;
public:
  static int Weight( const std::string &Quark );
  FFitConstMaker( std::string TypeParams, const Common::CommandLine &cl );
  void Make( std::string &sFileName ) override;
};
