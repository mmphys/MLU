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

extern const std::string DefaultColumnName;

enum class Freeze
{
  None,
  Central,
  One
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
  const std::string Base;
  const char * PrintPrefix; // non-null to print filenames as they are loaded
protected:
  using TNameMap = std::map<std::string, int>;
  TNameMap NameMap;
  std::vector<FileT> Files;
public:
  std::vector<std::string> opNames; // These are operator names referred to by Files
  FileCache( const std::string &base, const char * Prefix = nullptr ) : Base{base}, PrintPrefix{Prefix} {}
  void clear();
  // Add file to cache if not present - delay loading until accessed
  int GetIndex( const std::string &Filename );
  FileT &operator[]( int iIndex );
  const FileT &operator[]( int iIndex ) const;
  inline FileT &operator[]( const std::string &Filename ) { return (*this)[GetIndex(Filename)]; }
  inline const FileT &operator[]( const std::string &Filename ) const { return (*this)[GetIndex(Filename)]; }
};

// Allows a quark name to be read in as a string
// Extracts the momentum from the file name
struct QuarkReader : public std::string
{
  QP Convert( const std::string &Filename ) const;
};

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
protected:
  using KeyReadMapT = std::map<Key, std::string, LessKey>;
  template<typename K = Key, typename R = KeyRead>
  typename std::enable_if<std::is_same<K, R>::value, KeyReadMapT>::type
  ReadNameMap( std::ifstream &s );
  template<typename K = Key, typename R = KeyRead>
  typename std::enable_if<!std::is_same<K, R>::value, KeyReadMapT>::type
  ReadNameMap( std::ifstream &s );
  virtual void FrozenOptions( std::string &sOptions ) {}
public:
  KeyFileCache( const std::string &modelBase );
  virtual ~KeyFileCache() {}
  inline       M &operator[]( const Key &key )       { return model[KeyMap[key]]; }
  inline const M &operator[]( const Key &key ) const { return model[KeyMap[key]]; }
  virtual void clear();
  void Read( const std::string &Filename, const char *pszPrintPrefix );
  Vector GetVector( const Key &key, const std::string &ColumnName = DefaultColumnName );
};

using QDTModelMap = KeyFileCache<QDT, LessQDT>;

struct QPModelMap : public KeyFileCache<QP, LessQP, QuarkReader, Common::LessCaseInsensitive>
{
  using Base = KeyFileCache<QP, LessQP, QuarkReader, Common::LessCaseInsensitive>;
  std::string Spectator;
  QPModelMap( const std::string &modelBase ) : Base( modelBase ) {}
  std::string Get2ptName( const QP &key );
  void clear() override;
protected:
  void FrozenOptions( std::string &sOptions ) override;
};

// Base class for ratio construction
class Maker
{
public:
  //const std::string inBase;
  //const std::string C2Base;
  const std::string modelBase;
  const std::string outBase;
  //const Common::Momentum p;
  const bool bSymmetrise;
protected:
  const bool RegExSwap;
  std::regex RegExExt;
  QPModelMap model;
public:
  using CorrCache = FileCache<Fold>;
  CorrCache Cache2;
  CorrCache Cache3;
protected:
  inline void AppendOp( std::string &s, const std::string &Op ) { Append( s, Op ); }
  inline void AppendOps( std::string &s, const std::string &Snk, const std::string &Src)
  {
    AppendOp( s, Snk );
    AppendOp( s, Src );
  }
  std::string HeavyKey( const std::string &Heavy ) const;
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const QP &QPSnk, const QP &QPSrc,
                     const Vector &ESnk, const Vector &EmiSrc,
                     const std::string &opSnk, const std::string &opSrc ) = 0;
public:
  Maker( std::string &TypeParams, const Common::CommandLine &cl );
  virtual ~Maker() {}
  void Make( std::string &sFileName );
};

// Make Z_V: eq 2.7 pg 3 https://arxiv.org/pdf/1305.7217.pdf
class ZVMaker : public Maker
{
protected:
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const QP &QPSnk, const QP &QPSrc,
                     const Vector &ESnk, const Vector &ESrc, const std::string &opSnk, const std::string &opSrc );
public:
  ZVMaker( std::string &TypeParams, const Common::CommandLine &cl ) : Maker( TypeParams, cl ) {}
};

// Make R1 and R2: eq 2.8 pg 4 https://arxiv.org/pdf/1305.7217.pdf
class R1R2Maker : public Maker
{
protected:
  QDTModelMap ZVmi;
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const QP &QPSnk, const QP &QPSrc,
                     const Vector &ESnk, const Vector &ESrc, const std::string &opSnk, const std::string &opSrc );
public:
  R1R2Maker( std::string &TypeParams, const Common::CommandLine &cl );
};
