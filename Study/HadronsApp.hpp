/*************************************************************************************
 
 Create a Hadrons Application, top-down building all required dependencies
 Source file: HadronsApp.hpp
 Copyright (C) 2021
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

#include <MLU/Common.hpp>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

extern const std::string Sep;    // used inside filenames
extern const std::string Space;  // whitespace as field separator / human readable info
extern const std::string defaultMomName;
extern const std::string SeqMomentumName;
extern const Common::Momentum p0;
extern const std::string GaugeFieldName;

enum class Precision { Double, Single };

inline void Append( std::string &sDest, const std::string &s )
{
  if( !s.empty() )
  {
    sDest.append( Sep );
    sDest.append( s );
  }
};

inline void Appendsz( std::string &sDest, const char * psz, std::size_t Len )
{
  if( psz && Len )
  {
    sDest.append( Sep );
    sDest.append( psz, Len );
  }
};

// Family

enum class Family : int { Z2, ZF, GF, GP, GR };
extern const std::array<std::string, 5> FamilyNames;
typename std::underlying_type<Family>::type FamilyIndex( Family family );

inline const std::string &FamilyName( Family family )
{
  return FamilyNames[FamilyIndex( family )];
}

inline std::ostream& operator<<(std::ostream& os, const Family family)
{
  return os << FamilyName( family );
}

std::istream& operator>>(std::istream& is, Family &family);

// Species

enum class Species : int { Point, Wall, Region };
extern const std::array<std::string, 3> SpeciesNames;
typename std::underlying_type<Species>::type SpeciesIndex( Species species );

inline const std::string &SpeciesName( Species species )
{
  return SpeciesNames[SpeciesIndex( species )];
}

inline char SpeciesNameShort( Species species )
{
  return SpeciesName( species )[0];
}

inline std::ostream& operator<<(std::ostream& os, const Species species)
{
  return os << SpeciesName( species );
}

std::istream& operator>>(std::istream& is, Species &species);

// Taxonomy

class Taxonomy
{
  friend class AppParams;
protected:
  static const std::string HissyFitDefault;
public:
  Family family;
  Species species;

  Taxonomy( Family family_, Species species_ ) : family{family_}, species{species_} {}
  Taxonomy() : Taxonomy( Family::GF, Species::Point ) {}

  [[noreturn]] void HissyFit( const std::string &Message = HissyFitDefault ) const
  {
    std::stringstream ss;
    ss << Message << ": " << family << Common::CommaSpace << species;
    throw std::runtime_error( ss.str() );
  }

  // Allow extraction of the relevant parts without needing to specify
  inline operator Family() const { return family; }
  inline operator Species() const { return species; }
  inline const std::string &FamilyName() const { return ::FamilyName( family ); }
  inline void SinkSourceType( std::string &s, bool bReverse = false ) const
  {
    char c1 = SpeciesNameShort( species );
    char c2;
    switch( family )
    {
      case Family::GF:
        c2 = 'W';
        break;
      default:
        c2 = 'P';
        break;
    }
    s.append( 1, bReverse ? c2 : c1 );
    s.append( 1, bReverse ? c1 : c2 );
  }

  inline std::string FilePrefix( bool bReverse = false ) const
  {
    std::string s{ FamilyName() };
    SinkSourceType( s, bReverse );
    return s;
  }
  inline bool GaugeFixed() const { return family != Family::Z2; }
  inline void AppendFixed( std::string &s, bool bSmeared ) const
  {
    if( GaugeFixed() )
      Append( s, "fixed" );
    if( bSmeared )
      Append( s, "stout" );
  }
  void Validate() const
  {
    // Validate family and species separately
    FamilyIndex( family );
    SpeciesIndex( species );
    // Weed out invalid combinations
    if( family == Family::Z2 && species == Species::Wall    // Shouldn't be asking for wall sinks without gauge-fixing
     || family == Family::GR && species != Species::Wall )  // These look like point source and wall sink
      HissyFit();
  }
};

inline std::ostream& operator<<(std::ostream& os, const Taxonomy &taxonomy)
{
  return os << taxonomy.family << Common::Space << taxonomy.species;
}

std::istream& operator>>(std::istream& is, Taxonomy &taxonomy);
bool operator<( const Taxonomy &lhs, const Taxonomy &rhs );
bool operator==( const Taxonomy &lhs, const Taxonomy &rhs );

// Quark

struct Quark: Serializable {
GRID_SERIALIZABLE_CLASS_MEMBERS(Quark,
                                std::string,  flavour,
                                double,       mass,
                                unsigned int, Ls,
                                double,       M5,
                                std::string,  boundary,
                                std::string,  twist,
                                double,       scale,
                                // Solver parameters
                                unsigned int, maxIteration,
                                unsigned int, maxOuterIteration,
                                double,       residual,
                                bool,         GaugeSmear,
                                std::string,  eigenPack,
                                unsigned int, size,
                                bool,         multiFile,
                                bool,         redBlack,
                                bool,         eigenSinglePrecision) // EigenPack single precision
public:
  template<typename T> void CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const;
  inline bool MixedPrecision() const { return maxOuterIteration != 0; }
};

// Append

inline void Append( std::string &sDest, const std::string &s1, const std::string &s2 )
{
  Append( sDest, s1 );
  Append( sDest, s2 );
};

inline void Append( std::string &sDest, const std::string &s, int i )
{
  Append( sDest, s, std::to_string( i ) );
};

inline void Append( std::string &sDest, char c, const std::string &s )
{
  Append( sDest, std::string( 1, c ), s );
};

inline void Append( std::string &sDest, char c, int i )
{
  Append( sDest, std::string( 1, c ), std::to_string( i ) );
};

// Use Grid's definitions of gamma names ... but abbreviate Gamma to 'g'
inline void Append( std::string &sDest, Gamma::Algebra gamma )
{
  if( gamma < 0 || gamma > Gamma::nGamma )
    throw std::runtime_error( "Undefined Gamma::Algebra" );
  const char * psz = Gamma::name[gamma];
  std::size_t Len = std::strlen( psz );
  while( Len && std::isspace( psz[Len - 1] ) )
    Len--;
  if( Len )
  {
    sDest.append( Sep );
    while( Len )
    {
      static constexpr int GammaLen = 5;
      if( Len >= GammaLen && !std::strncmp( psz, Gamma::name[Gamma::Algebra::Gamma5], GammaLen ) )
      {
        sDest.append( 1, 'g' );
        Len -= GammaLen;
        psz += GammaLen;
      }
      else
      {
        sDest.append( 1, *psz++ );
        Len--;
      }
    }
  }
};

inline void Append( std::string &sDest, const Taxonomy &taxonomy )
{
  Append( sDest, taxonomy.FilePrefix() );
};

/*inline void Append( std::string &sDest, Family family )
{
  std::stringstream ss;
  ss << family;
  Append( sDest, ss.str() );
}*/

inline void Append( std::string &sDest, Species species )
{
  std::stringstream ss;
  ss << species;
  Append( sDest, ss.str() );
}

inline void AppendP( std::string &sDest, const Common::Momentum &p, const std::string & sMomentumName = defaultMomName )
{
  Append( sDest, sMomentumName, p.to_string( Sep ) );
}

inline void AppendPSeq( std::string &sDest, const Common::Momentum &p )
{
  Append( sDest, "ps", p.to_string( Sep ) );
}

inline void AppendT( std::string &sDest, int t )
{
  Append( sDest, 't', t );
}

inline void AppendPT( std::string &sDest, int t, const Common::Momentum &p )
{
  AppendP( sDest, p );
  AppendT( sDest, t );
};

inline void AppendDeltaT( std::string &sDest, int t )
{
  Append( sDest, "dt", t );
}

/**************************
 Parameters that apply to the entire application being built
**************************/

struct AppParams
{
  struct databaseOptions: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(databaseOptions,
                                  bool,         enable,
                                  std::string,  resultDb,
                                  bool,         makeStatDb,
                                  std::string,  applicationDbPrefix )
    };

  struct RunPar: Serializable {
      GRID_SERIALIZABLE_CLASS_MEMBERS(RunPar,
                                      std::string,  runId,
        Grid::Hadrons::Application::TrajRange,      trajCounter,
                                  databaseOptions,  dbOptions,
        Grid::Hadrons::VirtualMachine::GeneticPar,  genetic,
                                      int,          Nt,
                                      std::string,  Gauge,
                                      std::string,  GaugeFixed,
                                      std::string,  GaugeFixedXform,
             Grid::Hadrons::MGauge::GaugeFix::Par,  GaugeFix,
                       MGauge::StoutSmearing::Par,  StoutSmear,
                                      std::string,  SpatialPos,
                                      std::string,  RegionSize,
                                      std::string,  OutputBase)
  };

  RunPar Run;
  inline int TimeBound( int t ) const
  { return t < 0 ? Run.Nt - ((-t) % Run.Nt) : t % Run.Nt; }
  AppParams( XmlReader &r );
};

/**************************
 Base class for my wrappers of Hadrons::Modules
 These objects know how to create their coprresponding Hadrons Module and all dependencies
 **************************/

class HModList;

class HMod
{
protected:
  std::string name; // The name is what makes each module unique
public:
  static const std::string sSinglePrec;
  const Taxonomy tax;
public:
  inline const std::string &Name() const { return name; };
  HMod( HModList &ModList, const Taxonomy &taxonomy, int NameLen = 80 ) : tax{ taxonomy }
  { name.reserve( NameLen ); }
  virtual ~HMod() = default;
  virtual void AddDependencies( HModList &ModList ) const = 0;
};

/**************************
 List of Hadrons Modules, i.e. a wrapper for Hadrons::Application
 **************************/

class HModList
{
protected:
  std::map<std::string,std::unique_ptr<HMod>> list;
public:
  // These are used by modules when adding dependencies
  Application &application;
  const AppParams &params;
public:
  HModList( Application &application_, const AppParams &params_ )
  : application{application_}, params{params_} {}
  const std::string TakeOwnership( HMod *pHMod );
};

/**************************
 Point sink
**************************/

class ModSink : public HMod
{
public:
  static const std::string Prefix;
  const Common::Momentum p;
  ModSink(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p);
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Wall sink
**************************/

class ModSinkSmear : public HMod
{
public:
  static const std::string Prefix;
  const Common::Momentum p;
  ModSinkSmear(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p);
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Source
**************************/

class ModSource : public HMod
{
public:
  static const std::string Prefix;
  const Common::Momentum p;
  const int t;
  ModSource(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p, int t);
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Stout smeared gauge
**************************/

class ModGauge : public HMod
{
public:
  const bool bSmeared;
  const Precision precision;
  ModGauge( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared, Precision precision );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Gauge transform
**************************/

class ModGaugeXform : public HMod
{
public:
  static const std::string Suffix;
  const bool bSmeared;
  const Precision precision;
  ModGaugeXform( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared, Precision precision );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Action
**************************/

class ModAction : public HMod
{
public:
  static const std::string Prefix;
  const Quark &q;
  const bool bSmeared;
  const Precision precision;
  ModAction( HModList &ModList, const Taxonomy &taxonomy, const Quark &q, Precision precision );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Solver
**************************/

class ModSolver : public HMod
{
public:
  static const std::string Prefix;
  const Quark &q;
  const bool bSmeared;
  ModSolver( HModList &ModList, const Taxonomy &taxonomy, const Quark &q );
  virtual void AddDependencies( HModList &ModList ) const;
protected:
  template<typename T> std::string LoadEigenPack( HModList &ModList, Precision epPres ) const;
};

/**************************
 Propagator
**************************/

class ModProp : public HMod
{
protected:
  std::string Suffix;
public:
  static const std::string Prefix;
  static const std::string PrefixConserved;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q, const Common::Momentum &p, int t );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Sliced Propagator / wall sink. For 2pt contractions only
**************************/

class ModSlicedProp : public HMod
{
public:
  static const std::string Prefix;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModSlicedProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q, const Common::Momentum &p, int t );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Sequential Source
**************************/

class ModSourceSeq : public HMod
{
public:
  static const std::string Prefix;
  const Gamma::Algebra Current;
  const int deltaT;
  const Common::Momentum pSeq;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModSourceSeq( HModList &ModList, const Taxonomy &taxonomy, Gamma::Algebra Current, int deltaT,
                const Common::Momentum &pSeq, const Quark &q, const Common::Momentum &p, int t );
  virtual void AddDependencies( HModList &ModList ) const;
protected:
  template<typename T> void AddDependenciesT( HModList &ModList, typename T::Par &seqPar ) const;
};

/**************************
 Sequential propagator
**************************/

class ModPropSeq : public HMod
{
public:
  static const std::string Prefix;
  const Quark &qSeq;
  const Gamma::Algebra Current;
  const int deltaT;
  const Common::Momentum pSeq;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModPropSeq( HModList &ModList, const Taxonomy &taxonomy, const Quark &qSeq, Gamma::Algebra Current, int deltaT,
              const Common::Momentum &pSeq, const Quark &q, const Common::Momentum &p, int t );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 2pt contraction
 By default, momentum +p1 at source and -p1 at sink
 However, can also specify p2 for quark 2
    i.e. source momentum = (p1 - p2), sink momentum = (p2 - p1)
**************************/

extern const std::string ContractionPrefix;

class ModContract2pt : public HMod
{
protected:
  std::string FileName;
public:
  static const std::string Prefix;
  const Quark &q1;
  const Quark &q2;
  const Common::Momentum p1;
  const Common::Momentum p2;
  const int t;
  const Common::Momentum pSource;
  ModContract2pt( HModList &ModList, const Taxonomy &taxonomy,
                  const Quark &q1, const Quark &q2, const Common::Momentum &p1, int t,
                  const Common::Momentum &p2 = p0 );
  virtual void AddDependencies( HModList &ModList ) const;
};

/**************************
 Three-point contraction
 Quark qSrc at time t, decaying to quark qDst via current insertion at time (t + deltaT), with spectator anti-quark
 NB: if bHeavyAnti is true, then qSrc and qDst are anti-quarks and spectator is quark
 Momentum p is for the source, -p at the current and 0 at the sink (this could change later)
**************************/

class ModContract3pt : public HMod
{
protected:
  std::string FileName;
public:
  const bool bReverse;
  const Quark &qSnk;
  const Quark &qSrc;
  const Quark &qSpec;
  const Common::Momentum pSeq;
  const Common::Momentum p;
  const Gamma::Algebra Current;
  const int deltaT;
  const int t;
  const bool bHeavyAnti;
  ModContract3pt( HModList &ModList, const Taxonomy &taxonomy, bool bReverse, const Quark &qSnk, const Quark &qSrc,
                  const Quark &qSpec, const Common::Momentum &pSeq_, const Common::Momentum &p,
                  Gamma::Algebra Current, int deltaT, int t, bool bHeavyAnti );
  virtual void AddDependencies( HModList &ModList ) const;
};

