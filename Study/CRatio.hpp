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

using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Model = Common::Model<Scalar>;
//using Vector = Common::Vector<Scalar>;
//using VectorView = Common::VectorView<Scalar>;
using Column = Common::JackBootColumn<Scalar>;

extern Model DummyModel;
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

// A meson with a named momentum
struct MP
{
  std::string Meson;
  Common::Momentum p;
  std::string pName;
  MP() {};
  MP( const std::string &Meson_, const Common::Momentum &P,
      const std::string &pName_ =  Common::Momentum::DefaultPrefix )
  : Meson{Meson_}, p{P}, pName{pName_} {}
  MP( const std::string &Meson_, const Common::NamedMomentum &NP ) : MP(Meson_,NP.p,NP.Name) {}
  MP( const std::string &Meson_, const Common::MomentumPair &NP ) : MP(Meson_,NP.second,NP.first) {}
  /// Get name of two-point correlator
  std::string C2Name() const { return Meson + p.FileString( Common::Momentum::DefaultPrefix ); }
  Common::Param::Key PK( const std::string &FieldName = Common::ModelBase::EnergyPrefix ) const
  { return Common::Param::Key( C2Name(), FieldName ); }
};

std::ostream &operator<<( std::ostream &os, const MP &mp );

// Case insensitive compare of MesonMom
struct LessMP
{
  bool operator()( const MP &lhs, const MP &rhs ) const
  {
    int i = Common::CompareIgnoreCase( lhs.Meson, rhs.Meson );
    if( i )
      return i < 0;
    return lhs.p < rhs.p;
  }
};

// Allows a meson with momentum to be read in as a string
struct MPReader : public std::string
{
  MP operator()( const std::string &Filename ) const;
};

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
  bool empty() const { return NameMap.empty(); }
  // Add file to cache if not present - delay loading until accessed
  int GetIndex( const std::string &Filename, const char * printPrefix = nullptr );
  FileT &operator()( int iIndex, const char * printPrefix ); // This is what loads the file
  inline FileT &operator[]( int iIndex ) { return (*this)( iIndex, PrintPrefix ); };
  inline FileT &operator[]( const std::string &Filename ) { return (*this)[GetIndex(Filename)]; }
  inline const FileT &operator[]( const std::string &Filename ) const { return (*this)[GetIndex(Filename)]; }
};

// Static list of mappings from key to file. Can't be added to once read.
// Two-step construction (i.e. constructor then Read) so virtual functions in place
template<typename Key, typename LessKey = std::less<Key>, typename KeyRead = Key,
         typename LessKeyRead = LessKey, typename Scalar = ::Scalar >
struct KeyFileCache
{
  using M = Common::Model<Scalar>;
  using ColumnT = Common::JackBootColumn<Scalar>;
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
public:
  static constexpr int BadIndex{ FileCacheT::BadIndex };
  virtual void clear();
  bool empty() const { return model.empty(); }
  KeyFileCache() { clear(); }
  virtual ~KeyFileCache() {}
  // Use this to load a list of fit files
  void Read( const std::string &Filename, const char * PrintPrefix );
  M &operator()( const Key &key, const char * PrintPrefix ); // Make sure model loaded and print message
  inline M &operator[]( const Key &key ) { return operator()( key, nullptr ); }
  ColumnT GetColumn( const Key &key, const Common::Param::Key &pKey, std::size_t Index = 0 );
};

using ZVModelMap = KeyFileCache<std::string, Common::LessCaseInsensitive>;

using MPModelMap = KeyFileCache<MP, LessMP, MPReader, Common::LessCaseInsensitive>;

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
  MPModelMap EFit;
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
  virtual void Run( std::size_t &NumOK, std::size_t &Total, std::vector<std::string> &FileList );
};

class ZVRCommon : public Maker
{
protected:
  struct SSInfo
  {
    const std::string &op;
    const MP mp;
    const std::string q;
    Model &EModel;
    const Column E;
    SSInfo( MPModelMap &efit, const std::string &op_, MP mp_, const std::string &q_ )
    : op{op_}, mp{mp_}, q{q_}, EModel{efit(mp,LoadFilePrefix)}, E{efit.GetColumn(mp, mp.PK())} {}
  };
  inline void AppendOp( std::string &s, const std::string &Op ) { Append( s, Op ); }
  inline void AppendOps( std::string &s, const std::string &Snk, const std::string &Src)
  {
    AppendOp( s, Snk );
    AppendOp( s, Src );
  }
  virtual void ZVRMake( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                        const SSInfo &Snk, const SSInfo &Src ) = 0;
public:
  void Make( std::string &sFileName ) override;
  ZVRCommon( const Common::CommandLine &cl ) : Maker( cl ) {}
};

// Make Z_V: eq 2.7 pg 3 https://arxiv.org/pdf/1305.7217.pdf
class ZVMaker : public ZVRCommon
{
protected:
  void ZVRMake( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                const SSInfo &Snk, const SSInfo &Src ) override;
public:
  ZVMaker( const std::string &TypeParams, const Common::CommandLine &cl );
};

// Make R1 and R2: eq 2.8 pg 4 https://arxiv.org/pdf/1305.7217.pdf
class RMaker : public ZVRCommon
{
protected:
  const bool bAltR3;
  ZVModelMap ZVmi;
  void ZVRMake( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
             const SSInfo &Snk, const SSInfo &Src ) override;
public:
  RMaker( const std::string &TypeParams, const Common::CommandLine &cl );
};

struct FormFactor
{
  static const std::vector<std::string> ParamNames;
  static constexpr int EL{ 0 };
  static constexpr int mL{ 1 };
  static constexpr int mH{ 2 };
  static constexpr int qSq{ 3 };
  static constexpr int kMu{ 4 };
  static constexpr int melV0{ 5 };
  static constexpr int melVi{ 6 };
  static constexpr int fPar{ 7 };
  static constexpr int fPerp{ 8 };
  static constexpr int fPlus{ 9 };
  static constexpr int f0{ 10 };
  static constexpr int ELLat{ 11 };  // EL,  but E_i derived from E_0 and lattice dispersion relation
  static constexpr int qSqLat{ 12 }; // qSq, but E_i derived from E_0 and lattice dispersion relation
  static constexpr int tPlus{ 13 };
  static constexpr int tMinus{ 14 }; // aka q^2_{max}
  static constexpr int t0{ 15 };
  static constexpr int z_re{ 16 };
  static constexpr int z_im{ 17 };
  static constexpr int idxZV{ 18 };
  static constexpr int melV0Raw{ 19 };
  static constexpr int melViRaw{ 20 };
  const unsigned int N;
  const Scalar ap;
  const Scalar apInv;
  const bool bAdjustGammaSpatial;
  FormFactor( std::string TypeParams );
  void Write( std::string &OutFileName, const Model &CopyAttributesFrom,
              std::vector<std::string> &&SourceFileNames, int NumSamples,
              const Column &MHeavy, const Column &ELight,
              const Column &vT, const Common::Momentum &p,
              const Column &ZVPrevSrc, const Column &ZVPrevSnk, const Column &ZVMixed,
             // These only required for non-zero momentum
              const Column *pMLight, const Column *pvXYZ );
protected:
  ZVModelMap ZVmi;
  std::array<Column, 2> ZV;
  Model mEnsembleInfo;
};

// Make F parralel, F perpendicular, F+ and F0
class FMaker : public Maker, public FormFactor
{
  struct RatioFile
  {
    Common::Momentum p;
    //std::string Prefix;
    std::string Source;
    std::string Sink;
    struct Less { bool operator()( const RatioFile &lhs, const RatioFile &rhs ) const; };
    RatioFile( const Common::Momentum &p_, const std::string &Source_, const std::string &Sink_ )
    : p{p_}, Source{Source_}, Sink{Sink_} {}
  };
  using FileMap = std::map<RatioFile, std::vector<std::string>, RatioFile::Less>;

  const std::string i3Base;
protected:
  std::array<FileMap,2> map;
public:
  static int Weight( const std::string &Quark );
  FMaker( std::string TypeParams, const Common::CommandLine &cl )
  : Maker( cl ), FormFactor{TypeParams}, i3Base{ cl.SwitchValue<std::string>("i3") } {}
  void Make( std::string &sFileName ) override;
  void Run( std::size_t &NumOK, std::size_t &Total, std::vector<std::string> &FileList ) override;
};

// Make F parralel, F perpendicular, F+ and F0
class FFitConstMaker : public Maker, public FormFactor
{
  const std::string i3Base;
protected:
  int RatioNum;
  std::string qSnk;
  std::string qSrc;
  Common::Momentum p; std::string pName;
  std::string Prefix;
  std::string Suffix;
  std::string FitType;
  std::vector<int> FitParts;
public:
  static int Weight( const std::string &Quark );
  FFitConstMaker( std::string TypeParams, const Common::CommandLine &cl )
  : Maker( cl ), FormFactor{TypeParams}, i3Base{ cl.SwitchValue<std::string>("i3") } {}
  void Make( std::string &sFileName ) override;
};
