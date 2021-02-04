/*************************************************************************************
 
 Create XML for a 3-pt current insertion study
 Initially use Z2 wall sources.
 Source file: xml3pt.cpp
 Copyright (C) 2020
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

#include "../Analyze/Common.hpp"
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

static const std::string Sep{ "_" };    // used inside filenames
static const std::string Space{ " " };  // whitespace as field separator / human readable info
static const Common::Momentum p0(0,0,0);

static const std::string GaugeFieldName{"gauge"};

inline void Append( std::string &sDest, const std::string &s )
{
  if( !s.empty() )
    sDest.append( Sep );
  sDest.append( s );
};

enum class Family : int { Z2, ZF, GF, GP };
static const std::array<std::string, 4> FamilyNames{ "Z2", "ZF", "GF", "GP" };

typename std::underlying_type<Family>::type FamilyIndex( Family family )
{
  const auto iFamily = static_cast<typename std::underlying_type<Family>::type>( family );
  if( iFamily < 0 || iFamily >= FamilyNames.size() )
    throw std::runtime_error( "Unknown family " + std::to_string( iFamily ) );
  return iFamily;
}

inline const std::string &FamilyName( Family family )
{
  return FamilyNames[FamilyIndex( family )];
}

inline std::ostream& operator<<(std::ostream& os, const Family family)
{
  return os << FamilyName( family );
}

std::istream& operator>>(std::istream& is, Family &family)
{
  std::string s;
  if( is >> s )
  {
    int i;
    for(i = 0; i < FamilyNames.size() && !Common::EqualIgnoreCase( s, FamilyNames[i] ); ++i)
      ;
    if( i < FamilyNames.size() )
      family = static_cast<Family>( i );
    else
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

enum class Species : int { Point, Wall, Region };
static const std::array<std::string, 3> SpeciesNames{ "Point", "Wall", "Region" };

typename std::underlying_type<Species>::type SpeciesIndex( Species species )
{
  const auto iSpecies = static_cast<typename std::underlying_type<Species>::type>( species );
  if( iSpecies < 0 || iSpecies >= ::SpeciesNames.size() )
    throw std::runtime_error( "Unknown species " + std::to_string( iSpecies ) );
  return iSpecies;
}

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

std::istream& operator>>(std::istream& is, Species &species)
{
  std::string s;
  if( is >> s )
  {
    int i;
    for(i = 0; i < SpeciesNames.size() && !Common::EqualIgnoreCase( s, SpeciesNames[i] ); ++i)
      ;
    if( i < SpeciesNames.size() )
      species = static_cast<Species>( i );
    else
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

//enum class Sink : int { Point, Wall };

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

  std::string FilePrefix() const
  {
    std::string s{ FamilyName() };
    s.append( 1, SpeciesNameShort( species ) );
    switch( family )
    {
      case Family::GF:
        s.append( 1, 'W' );
        break;
      default:
        s.append( 1, 'P' );
        break;
    }
    return s;
  }
  inline bool GaugeFixed() const { return family != Family::Z2; }
  inline void AppendFixed( std::string &s, bool bSmeared = false ) const
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
    if( family == Family::Z2 && species == Species::Wall ) // Shouldn't be asking for wall sinks without gauge-fixing
      HissyFit();
  }
};

const std::string Taxonomy::HissyFitDefault{ "Invalid Taxonomy" };

inline std::ostream& operator<<(std::ostream& os, const Taxonomy &taxonomy)
{
  return os << taxonomy.family << Common::Space << taxonomy.species;
}

std::istream& operator>>(std::istream& is, Taxonomy &taxonomy)
{
  Taxonomy local;
  std::string s;
  if( is >> local.family && is >> local.species )
  {
    local.Validate();
    taxonomy = local;
  }
  return is;
}

// These make taxa sortable
bool operator<( const Taxonomy &lhs, const Taxonomy &rhs )
{
  return       FamilyIndex ( lhs.family ) <  FamilyIndex ( rhs.family )
         || (  FamilyIndex ( lhs.family ) == FamilyIndex ( rhs.family )
            && SpeciesIndex( lhs.species) <  SpeciesIndex( rhs.species) );
}

bool operator==( const Taxonomy &lhs, const Taxonomy &rhs )
{
  return FamilyIndex( lhs.family ) == FamilyIndex ( rhs.family )
     && SpeciesIndex( lhs.species) == SpeciesIndex( rhs.species);
}

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
                                int,          maxIteration,
                                double,       residual,
                                bool,         GaugeSmear,
                                std::string,  EigenPackFilename )
public:
  void CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const;
};

void Quark::CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const
{
  MAction::ScaledDWF::Par Par;
  Par.gauge    = std::move( Gauge );
  Par.Ls       = Ls;
  Par.mass     = mass;
  Par.M5       = M5;
  Par.scale    = scale; // 1 = Shamir, 2 = Möbius
  Par.boundary = boundary;
  Par.twist    = twist;
  app.createModule<MAction::ScaledDWF>( name, Par );
}

static const Gamma::Algebra algInsert[] = {
  Gamma::Algebra::Gamma5,
  Gamma::Algebra::GammaTGamma5,
  //Gamma::Algebra::GammaT,
  //Gamma::Algebra::GammaX,
  //Gamma::Algebra::GammaY,
  //Gamma::Algebra::GammaZ,
};
static constexpr int NumInsert{ sizeof( algInsert ) / sizeof( algInsert[0] ) };
static const std::array<const std::string *, NumInsert> algInsertName {
  &Common::Gamma::nameShort[ static_cast<int>( Common::Gamma::Algebra::Gamma5 ) ],
  &Common::Gamma::nameShort[ static_cast<int>( Common::Gamma::Algebra::GammaTGamma5 ) ],
  //&Common::Gamma::nameShort[ static_cast<int>( Common::Gamma::Algebra::GammaT ) ],
  //&Common::Gamma::nameShort[ static_cast<int>( Common::Gamma::Algebra::GammaX ) ],
  //&Common::Gamma::nameShort[ static_cast<int>( Common::Gamma::Algebra::GammaY ) ],
  //&Common::Gamma::nameShort[ static_cast<int>( Common::Gamma::Algebra::GammaZ ) ],
};

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

inline void AppendP( std::string &sDest, const Common::Momentum &p )
{
  Append( sDest, 'p', p.to_string( Sep ) );
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
                                  std::string,  resultDb,
                                  bool,         makeStatDb,
                                  std::string,  applicationDbPrefix )
    };

  struct RunPar: Serializable {
      GRID_SERIALIZABLE_CLASS_MEMBERS(RunPar,
        Grid::Hadrons::Application::TrajRange,      trajCounter,
                                  databaseOptions,  dbOptions,
        Grid::Hadrons::VirtualMachine::GeneticPar,  genetic,
                                      int,          Nt,
        Grid::Hadrons::Application::TrajRange,      Timeslices,
                                      std::string,  runId,
                                      std::string,  Taxa,
                                      std::string,  Momenta,
                                      std::string,  deltaT,
                                      std::string,  SpatialPos,
                                      std::string,  RegionSize,
                                      std::string,  OutputBase,
                                      std::string,  Gauge,
             Grid::Hadrons::MGauge::GaugeFix::Par,  GaugeFix,
                       MGauge::StoutSmearing::Par,  StoutSmear,
                                      bool,         TwoPoint,
                                      bool,         HeavyQuark,
                                      bool,         HeavyAnti,
                                      bool,         DoNegativeMomenta,
                                      bool,         Run,
                                      std::string,  JobXmlPrefix)
  };

  RunPar Run;
  std::vector<Taxonomy> Taxa;
  std::vector<Quark> HeavyQuarks;
  std::vector<Quark> SpectatorQuarks;
  std::vector<Quark> SpectatorQuarks2pt;
  std::vector<Common::Momentum> Momenta;
  bool ThreePoint;
  std::vector<int> deltaT;
  inline int TimeBound( int t ) const
  { return t < 0 ? Run.Nt - ((-t) % Run.Nt) : t % Run.Nt; }
  AppParams( const std::string sXmlFilename, const std::string sRunPrefix );
  std::string RunID() const;
protected:
  static std::vector<Quark> ReadQuarks( XmlReader &r, const std::string &qType, std::size_t Min = 1 );
};

template <typename T> void NoDuplicates( const std::vector<T> &v, const std::string sErrorPrefix, bool bEmptyOK = false )
{
  if( v.empty() )
  {
    if( !bEmptyOK )
      throw std::runtime_error( sErrorPrefix + " is empty" );
  }
  else if( v.size() > 1 )
  {
    std::vector<T> vc{ v };
    std::sort( vc.begin(), vc.end() );
    const auto dup = std::adjacent_find( vc.begin(), vc.end() );
    if( dup != vc.end() )
    {
      std::stringstream ss;
      ss << sErrorPrefix << " contains duplicates, e.g. " << *dup;
      throw std::runtime_error( ss.str() );
    }
  }
}

AppParams::AppParams( const std::string sXmlFilename, const std::string sRunPrefix )
{
  static const std::string sXmlTopLevel{ "Study2" };
  static const std::string sXmlTagName{ "RunPar" };
  XmlReader r( sXmlFilename, false, sXmlTopLevel );
  read( r, sXmlTagName, Run );
  if( !sRunPrefix.empty() )
    Run.JobXmlPrefix = sRunPrefix + Run.JobXmlPrefix;
  HeavyQuarks        = ReadQuarks( r, "Heavy" );
  SpectatorQuarks    = ReadQuarks( r, "Spectator" );
  SpectatorQuarks2pt = ReadQuarks( r, "SpectatorExtra2pt", 0 );
  Momenta = Common::ArrayFromString<Common::Momentum>( Run.Momenta );
  // Check the taxonomy
  Taxa = Common::ArrayFromString<Taxonomy>( Run.Taxa );
  NoDuplicates( Taxa, "Taxa" );
  // Check parameters make sense
  //if( !( Run.TwoPoint || Run.HeavyQuark ||Run.HeavyAnti ) )
    //throw std::runtime_error( "At least one must be true of: TwoPoint; HeavyQuark; or HeavyAnti" );
  if( Momenta.empty() )
    throw std::runtime_error( "There should be at least one momentum" );
  ThreePoint = Run.HeavyQuark || Run.HeavyAnti;
  if( ThreePoint )
  {
    deltaT = Common::ArrayFromString<int>( Run.deltaT );
    if( deltaT.empty() )
      throw std::runtime_error( "There should be at least one deltaT" );
  }
}

std::vector<Quark> AppParams::ReadQuarks( XmlReader &r, const std::string &qType, std::size_t Min )
{
  std::vector<Quark> vq;
  r.readDefault( qType, vq );
  if( vq.size() < Min )
    throw std::runtime_error( "Only " + std::to_string( vq.size() ) + Common::Space + qType + " quarks" );
  return vq;
}

// Make a unique RunID string that completely describes the run
std::string AppParams::RunID() const
{
  std::ostringstream s;
  {
    // Start with all the Families and sink types
    std::vector<Taxonomy> taxSorted{ Taxa };
    std::sort( taxSorted.begin(), taxSorted.end() );
    for( std::size_t i = 0; i < taxSorted.size(); ++i )
    {
      if( i == 0 || taxSorted[i].family != taxSorted[i - 1].family )
      {
        if( i )
          s << Sep;
        s << taxSorted[i].family;
      }
      s << SpeciesNameShort( taxSorted[i].species );
    }
  }
  for( const Quark &q : SpectatorQuarks )
    s << Sep << q.flavour;
  if( Run.TwoPoint )
    s << Sep << "2pt";
  if( Run.HeavyQuark )
    s << Sep << "quark";
  if( Run.HeavyAnti )
    s << Sep << "anti";
  for( const Quark &q : HeavyQuarks )
    s << Sep << q.flavour;
  s << Sep << "t" << Sep << Run.Timeslices.start << Sep << Run.Timeslices.end << Sep << Run.Timeslices.step;
  if( !Run.DoNegativeMomenta )
    s << Sep << "pos";
  s << Sep << "p";
  for( const Common::Momentum &p : Momenta )
    s << Sep << p.to_string( Sep );
  if( ThreePoint )
  {
    s << Sep << "dT";
    for( int dT : deltaT )
      s << Sep << std::to_string( dT );
  }
  return s.str();
}

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

// Delete the module if it already exists, otherwise add it to the application and create all dependencies

const std::string HModList::TakeOwnership( HMod *pHMod )
{
  assert( pHMod && "HModList::TakeOwnership() : Null pointer" );
  std::string ModName{ pHMod->Name() };
  auto it{ list.find( ModName ) };
  if( it != list.end() )
    delete pHMod; // Already exists
  else
  {
    // I've not seen this yet - Make all its dependencies and add it to my list
    pHMod->AddDependencies( *this );
    list.emplace( std::make_pair( ModName, std::unique_ptr<HMod>( pHMod ) ) );
  }
  return ModName;
}

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

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p_)
: HMod(ModList, taxonomy), p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSink::AddDependencies( HModList &ModList ) const
{
  MSink::ScalarPoint::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
}

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

const std::string ModSinkSmear::Prefix{ "Sink_Smear" };

ModSinkSmear::ModSinkSmear(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p_)
: HMod(ModList, taxonomy), p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSinkSmear::AddDependencies( HModList &ModList ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::Point>(Name(), sinkPar);
}

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

const std::string ModSource::Prefix{ "Source" };

ModSource::ModSource(HModList &ModList, const Taxonomy &tx, const Common::Momentum &p_, int t_)
// Silently translate ZF sources into Z2 sources
: HMod(ModList, tx.family == Family::ZF ? Taxonomy(Family::Z2,tx.species) : tx), p{p_}, t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  AppendPT( name, t, p );
}

void ModSource::AddDependencies( HModList &ModList ) const
{
  switch( tax.family )
  {
    case Family::Z2:
      if( p )
      {
        MSource::MomentumPhase::Par par;
        par.src = ModList.TakeOwnership( new ModSource( ModList, tax, Common::Momentum(0,0,0), t ) );
        par.mom = p.to_string4d( Common::Space );
        ModList.application.createModule<MSource::MomentumPhase>(name, par);
      }
      else
      {
        MSource::Z2::Par par;
        par.tA = t;
        par.tB = t;
        ModList.application.createModule<MSource::Z2>(name, par);
      }
      break;
    case Family::GP:
      {
        MSource::Point::Par par;
        par.position = ModList.params.Run.SpatialPos + Common::Space + std::to_string( t );
        ModList.application.createModule<MSource::Point>(name, par);
      }
      break;
    case Family::GF:
      {
        MSource::Wall::Par par;
        par.tW = t;
        par.mom = p.to_string4d( Common::Space );
        ModList.application.createModule<MSource::Wall>(name, par);
      }
      break;
    default:
      {
        std::stringstream ss;
        ss << "Unrecognised family \"" << tax.family << "\"";
        throw std::runtime_error( ss.str() );
      }
  }
}

/**************************
 Stout smeared gauge
**************************/

class ModGauge : public HMod
{
public:
  const bool bSmeared;
  ModGauge( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared );
  virtual void AddDependencies( HModList &ModList ) const;
};

ModGauge::ModGauge( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared_ )
: HMod(ModList, taxonomy), bSmeared{ bSmeared_ }
{
  name = GaugeFieldName;
  tax.AppendFixed( name, bSmeared );
}

void ModGauge::AddDependencies( HModList &ModList ) const
{
  if( bSmeared )
  {
    // Stout smeared field
    MGauge::StoutSmearing::Par stoutPar( ModList.params.Run.StoutSmear );
    stoutPar.gauge = tax.GaugeFixed() ? ModList.TakeOwnership( new ModGauge( ModList, tax, false ) ) : GaugeFieldName;
    ModList.application.createModule<MGauge::StoutSmearing>(name, stoutPar);
  }
  else if( tax.GaugeFixed() )
  {
    // Gauge-fixed
    MGauge::GaugeFix::Par Par( ModList.params.Run.GaugeFix );
    Par.gauge = GaugeFieldName;
    ModList.application.createModule<MGauge::GaugeFix>(name, Par);
  }
}

/**************************
 Action
**************************/

class ModAction : public HMod
{
public:
  static const std::string Prefix;
  const Quark &q;
  const bool bSmeared;
  ModAction( HModList &ModList, const Taxonomy &taxonomy, const Quark &q );
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModAction::Prefix{ "DWF" };

ModAction::ModAction( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_ )
: HMod( ModList, taxonomy ), q{q_}, bSmeared{q_.GaugeSmear && !ModList.params.Run.Gauge.empty()}
{
  name = Prefix;
  tax.AppendFixed( name, bSmeared );
  Append( name, q.flavour );
}

void ModAction::AddDependencies( HModList &ModList ) const
{
  std::string Gauge{ModList.TakeOwnership(new ModGauge( ModList, tax, bSmeared))};
  q.CreateAction( ModList.application, name, std::move( Gauge ) );
}

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
};

const std::string ModSolver::Prefix{ "CG" };

ModSolver::ModSolver( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_ )
: HMod( ModList, taxonomy ), q( q_ ), bSmeared{q_.GaugeSmear && !ModList.params.Run.Gauge.empty()}
{
  name = Prefix;
  tax.AppendFixed( name, bSmeared );
  Append( name, q.flavour );
}

void ModSolver::AddDependencies( HModList &ModList ) const
{
  // solvers
  MSolver::RBPrecCG::Par solverPar;
  if( ModList.params.Run.Gauge.length() && q.EigenPackFilename.length() )
  {
    // eigenpacks for deflation
    MIO::LoadFermionEigenPack::Par epPar;
    epPar.filestem = q.EigenPackFilename;
    epPar.multiFile = false;
    epPar.size = 600;
    epPar.Ls = q.Ls;
    if( tax.GaugeFixed() )
      epPar.gaugeXform = ModList.TakeOwnership( new ModGauge( ModList, tax, false ) ) + "_xform";
    solverPar.eigenPack = "epack_" + q.flavour;
    ModList.application.createModule<MIO::LoadFermionEigenPack>(solverPar.eigenPack, epPar);
  }
  solverPar.action       = ModList.TakeOwnership( new ModAction( ModList, tax, q ) );
  solverPar.residual     = q.residual;
  solverPar.maxIteration = q.maxIteration;
  ModList.application.createModule<MSolver::RBPrecCG>(name, solverPar);
}

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

const std::string ModProp::Prefix{ "Prop" };
const std::string ModProp::PrefixConserved{ "Ward" };

ModProp::ModProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod( ModList, taxonomy ), q{q_}, p{p_}, t{t_}
{
  Suffix = tax.FamilyName();
  Append( Suffix, q.flavour );
  AppendPT( Suffix, t, p );
  // Name of the propagator
  name = Prefix;
  Append( name, Suffix );
}

void ModProp::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par par;
  par.source = ModList.TakeOwnership( new ModSource( ModList, tax, p, t ) );
  par.solver = ModList.TakeOwnership( new ModSolver( ModList, tax, q ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, par);
  // Check residual mass for zero-momentum propagators
  if( !p )
  {
    std::string WardName{ PrefixConserved };
    Append( WardName, Suffix );
    // Check residual mass
    MContraction::WardIdentity::Par WIP;
    WIP.prop = name + "_5d";
    WIP.action = ModList.TakeOwnership( new ModAction( ModList, tax, q ) );
    WIP.source = par.source;
    WIP.mass = q.mass;
    WIP.output = ModList.params.Run.OutputBase + PrefixConserved + "/" + Suffix;
    ModList.application.createModule<MContraction::WardIdentity>(WardName, WIP);
  }
}

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

const std::string ModSlicedProp::Prefix{ "SlicedProp" };

ModSlicedProp::ModSlicedProp( HModList &ML, const Taxonomy &tax_, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod( ML, tax_ ), q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

void ModSlicedProp::AddDependencies( HModList &ModList ) const
{
  MSink::Smear::Par smearPar;
  smearPar.q = ModList.TakeOwnership( new ModProp( ModList, tax, q, p, t ) );
  smearPar.sink = ModList.TakeOwnership( new ModSinkSmear( ModList, tax, p0 ) ); // TODO: Have I put momentum in properly?
  ModList.application.createModule<MSink::Smear>(name, smearPar);
}

/**************************
 Sequential Source
**************************/

class ModSourceSeq : public HMod
{
public:
  static const std::string Prefix;
  const int Current;
  const int deltaT;
  const Common::Momentum pSeq;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModSourceSeq( HModList &ModList, const Taxonomy &taxonomy, int Current, int deltaT, const Common::Momentum &pSeq,
                const Quark &q, const Common::Momentum &p, int t );
  virtual void AddDependencies( HModList &ModList ) const;
protected:
  template<typename T> void AddDependenciesT( HModList &ModList, typename T::Par &seqPar ) const;
};

const std::string ModSourceSeq::Prefix{ "Seq" + ModSource::Prefix };

ModSourceSeq::ModSourceSeq( HModList &ML, const Taxonomy &tax_, int current_, int deltaT_, const Common::Momentum &pSeq_,
                            const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod(ML,tax_), Current{current_}, deltaT{deltaT_}, pSeq{pSeq_}, q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, tax.species );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

template<typename T> void ModSourceSeq::AddDependenciesT( HModList &ModList, typename T::Par &seqPar ) const
{
  seqPar.q = ModList.TakeOwnership( new ModProp( ModList, tax, q, p, t ) );
  seqPar.mom = pSeq.to_string4d( Space );
  seqPar.gamma = algInsert[Current];
  ModList.application.createModule<T>(name, seqPar);
}

void ModSourceSeq::AddDependencies( HModList &ModList ) const
{
  const int tSink{ModList.params.TimeBound( t + deltaT )};
  switch( tax.species )
  {
    case Species::Point:
    {
      using M = MSource::SeqGamma;
      typename M::Par seqPar;
      seqPar.tA  = tSink;
      seqPar.tB  = tSink;
      AddDependenciesT<M>( ModList, seqPar );
    }
      break;
    case Species::Wall:
    {
      using M = MSource::SeqGammaWall;
      typename M::Par seqPar;
      seqPar.tA  = tSink;
      seqPar.tB  = tSink;
      AddDependenciesT<M>( ModList, seqPar );
    }
      break;
    case Species::Region:
    {
      using M = MSource::SeqGammaRegion;
      typename M::Par seqPar;
      seqPar.LowerLeft  = ModList.params.Run.SpatialPos + Common::Space + std::to_string( tSink );
      seqPar.RegionSize = "1 1 1 1";
      AddDependenciesT<M>( ModList, seqPar );
    }
      break;
  }
}

/**************************
 Sequential propagator
**************************/

class ModPropSeq : public HMod
{
public:
  static const std::string Prefix;
  const Quark &qSeq;
  const int Current;
  const int deltaT;
  const Common::Momentum pSeq;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModPropSeq( HModList &ModList, const Taxonomy &taxonomy, const Quark &qSeq, int current, int deltaT,
              const Common::Momentum &pSeq, const Quark &q, const Common::Momentum &p, int t );
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModPropSeq::Prefix{ "Seq" + ModProp::Prefix };

ModPropSeq::ModPropSeq( HModList &ML, const Taxonomy &tax_, const Quark &qSeq_, int current_, int deltaT_,
                        const Common::Momentum &pSeq_, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod(ML,tax_),qSeq{qSeq_},Current{current_},deltaT{deltaT_},pSeq{pSeq_},q{q_},p{p_},t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, tax.species );
  Append( name, qSeq.flavour );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

void ModPropSeq::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par quarkPar;
  quarkPar.source = ModList.TakeOwnership( new ModSourceSeq( ModList, tax, Current, deltaT, pSeq, q, p, t ) );
  quarkPar.solver = ModList.TakeOwnership( new ModSolver( ModList, tax, qSeq ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, quarkPar);
}

/**************************
 2pt contraction
**************************/

const std::string ContractionPrefix{ "meson" };

class ModContract2pt : public HMod
{
protected:
  std::string FileName;
public:
  static const std::string Prefix;
  const Quark &q1;
  const Quark &q2;
  const Common::Momentum p;
  const int t;
  ModContract2pt( HModList &ModList, const Taxonomy &taxonomy,
                  const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_);
  virtual void AddDependencies( HModList &ModList ) const;
};

ModContract2pt::ModContract2pt( HModList &ModList, const Taxonomy &taxonomy,
                                const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_)
// Two-point Region-sinks don't really exist - silently translate to point-sinks
: HMod(ModList, taxonomy == Species::Region ? Taxonomy(taxonomy.family, Species::Point) : taxonomy),
  q1{q1_}, q2{q2_}, p{p_}, t{t_}
{
  std::string Prefix{ tax.FilePrefix() };
  std::string s{ q2.flavour };
  Append( s, q1.flavour );
  AppendPT( s, t, p );
  static const std::string Prefixfamily{ "2pt" };
  FileName = ModList.params.Run.OutputBase;
  FileName.append( Prefixfamily );
  FileName.append( 1, '/' );
  FileName.append( Prefix );
  FileName.append( 1, '/' );
  FileName.append( Prefix ); // Repetition not an error. Part of directory hierarchy and filename
  Append( FileName, s );
  name = ContractionPrefix;
  Append( name, Prefixfamily );
  Append( name, Prefix );
  Append( name, s );
}

void ModContract2pt::AddDependencies( HModList &ModList ) const
{
  if( name.empty() )
    return; // I can't do this for two-point. Just silently ignore
  MContraction::Meson::Par mesPar;
  mesPar.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  mesPar.gammas = "all";
  if( tax == Species::Wall )
  {
    mesPar.q1 = ModList.TakeOwnership( new ModSlicedProp( ModList, tax, q1, p, t ) );
    mesPar.q2 = ModList.TakeOwnership( new ModSlicedProp( ModList, tax, q2, p, t ) );
  }
  else
  {
    mesPar.q1 = ModList.TakeOwnership( new ModProp( ModList, tax, q1, p, t ) );
    mesPar.q2 = ModList.TakeOwnership( new ModProp( ModList, tax, q2, p, t ) );
  }
  mesPar.sink = ModList.TakeOwnership( new ModSink( ModList, tax, p0 ) );
  ModList.application.createModule<MContraction::Meson>(name, mesPar);
}

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
  const Quark &qSnk;
  const Quark &qSrc;
  const Quark &qSpec;
  const Common::Momentum p;
  const int t;
  const bool bHeavyAnti;
  const int Current;
  const int deltaT;
  const bool bSinkMom;
  ModContract3pt( HModList &ModList, const Taxonomy &taxonomy, const Quark &qSnk, const Quark &qSrc, const Quark &qSpec,
                  const Common::Momentum &p, int t, bool bHeavyAnti, int Current, int deltaT, bool bSinkMom );
  virtual void AddDependencies( HModList &ModList ) const;
};

ModContract3pt::ModContract3pt(HModList &ML,const Taxonomy &tax_,const Quark &Snk_,const Quark &Src_,const Quark &Spec_,
                      const Common::Momentum &p_, int t_, bool bHA_, int Cur_, int dT_, bool bSM_ )
: HMod(ML, tax_),qSnk{Snk_},qSrc{Src_},qSpec{Spec_},p{p_},t{t_},bHeavyAnti{bHA_},Current{Cur_},deltaT{dT_},bSinkMom{bSM_}
{
  std::string Prefix1{ "3pt" };
  Append( Prefix1, qSpec.flavour );
  std::string Prefix2{ tax.FilePrefix() };
  Append( Prefix2, bHeavyAnti ? "anti" : "quark" );
  std::string s{ qSnk.flavour };
  Append( s, qSrc.flavour );
  Append( s, *algInsertName[Current] );
  AppendDeltaT( s, deltaT );
  AppendPT( s, t, p );
  // File name
  FileName = ML.params.Run.OutputBase + Prefix1;
  if( bSinkMom )
    Append( FileName, "ms" );
  FileName.append( 1, '/' );
  FileName.append( Prefix2 );
  FileName.append( 1, '/' );
  FileName.append( Prefix2 ); // Repetition not a mistake. This is part of directory hierarchy and filename
  Append( FileName, s );
  // Object name
  name = ContractionPrefix;
  Append( name, Prefix1 );
  Append( name, Prefix2 );
  Append( name, s );
}

void ModContract3pt::AddDependencies( HModList &ML ) const
{
  MContraction::Meson::Par par;
  par.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  par.gammas = "all";
  const bool bInvertSeq{ !bHeavyAnti };
  if( bInvertSeq )
  {
    par.q1 = ML.TakeOwnership( new ModProp( ML, tax, qSrc, p, t ) );
    par.q2 = ML.TakeOwnership(new ModPropSeq(ML,tax, qSnk, Current, deltaT, bSinkMom ? p : p0, qSpec, p0, t ) );
  }
  else
  {
    par.q1 = ML.TakeOwnership(new ModPropSeq(ML,tax, qSnk, Current, deltaT, bSinkMom ? -p : p0, qSpec, p , t ) );
    par.q2 = ML.TakeOwnership( new ModProp( ML, tax, qSrc, p0, t ) );
  }
  par.sink = ML.TakeOwnership( new ModSink( ML, tax, bSinkMom ? p0 : -p ) );
  ML.application.createModule<MContraction::Meson>(name, par);
}

/**************************
 Make application
**************************/

class AppMaker
{
public:
  Application application;
  HModList l;
protected:
  Application Setup( const AppParams &params );
public:
  explicit AppMaker( const AppParams &params )
  : application{ Setup( params ) }, l( application, params ) {}
  void Make();
};

// One-time initialisation
Application AppMaker::Setup( const AppParams &params )
{
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter  = params.Run.trajCounter;
  globalPar.genetic      = params.Run.genetic;
  globalPar.runId        = params.Run.runId;
  // database options
  globalPar.database.resultDb   = params.Run.dbOptions.resultDb;
  globalPar.database.makeStatDb = params.Run.dbOptions.makeStatDb;
  static const std::string RunID{ params.RunID() };
  globalPar.database.applicationDb = params.Run.dbOptions.applicationDbPrefix + RunID + ".db";
  Grid::Hadrons::makeFileDir( globalPar.database.applicationDb );
  globalPar.database.restoreModules = false;
  globalPar.database.restoreMemoryProfile = false;
  globalPar.database.restoreSchedule = false;
  globalPar.scheduleFile = params.Run.JobXmlPrefix + RunID + ".sched";
  globalPar.saveSchedule = true;
  Application application( globalPar );
  // gauge field
  if( params.Run.Gauge.empty() )
  {
    application.createModule<MGauge::Unit>( GaugeFieldName );
  }
  else
  {
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = params.Run.Gauge;
    application.createModule<MIO::LoadNersc>( GaugeFieldName, gaugePar );
  }
  return application;
}

void AppMaker::Make()
{
  // Unitary (sea mass) spectators
  bool bFirstSpec{ true };
  for( const Quark &qSpectator : l.params.SpectatorQuarks )
  {
    for(unsigned int t = l.params.Run.Timeslices.start; t < l.params.Run.Timeslices.end; t += l.params.Run.Timeslices.step)
    {
      for( const Taxonomy &tax : l.params.Taxa )
      {
        for( Common::Momentum p : l.params.Momenta )
        {
          for( int pDoNeg = 0; pDoNeg < ( ( l.params.Run.DoNegativeMomenta && p ) ? 2 : 1 ); ++pDoNeg )
          {
            for( const Quark &qH1 : l.params.HeavyQuarks )
            {
              if( l.params.Run.TwoPoint )
              {
                l.TakeOwnership( new ModContract2pt( l, tax, qSpectator, qH1, p, t ) );
                l.TakeOwnership( new ModContract2pt( l, tax, qH1, qSpectator, p, t ) );
                if( bFirstSpec )
                {
                  // Additional spectators for 2pt functions only
                  for( const Quark &qSpec2 : l.params.SpectatorQuarks2pt )
                  {
                    l.TakeOwnership( new ModContract2pt( l, tax, qSpec2, qH1, p, t ) );
                    l.TakeOwnership( new ModContract2pt( l, tax, qH1, qSpec2, p, t ) );
                  }
                }
              }
              else if( !l.params.ThreePoint )
              {
                // We are only performing residual-mass checks on each propagator
                l.TakeOwnership( new ModProp( l, tax, qH1, p, t ) );
                //l.TakeOwnership( new ModProp( l, tax, qSpectator, p,t) );
              }
              if( l.params.ThreePoint )
              {
                for( const Quark &qH2 : l.params.HeavyQuarks )
                {
                  for( int iHeavy  = l.params.Run.HeavyQuark ? 0 : 1;
                           iHeavy <= l.params.Run.HeavyAnti  ? 1 : 0; ++iHeavy )
                  {
                    const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
                    if( !p || qH1.mass >= qH2.mass )
                    {
                      for( int deltaT : l.params.deltaT )
                      {
                        for( int j = 0; j < NumInsert; j++ )
                        {
                          l.TakeOwnership( new ModContract3pt( l, tax, qH1, qH2, qSpectator, p,
                                                               t, bHeavyAnti, j, deltaT, false ) );
                          if( p && qH1.mass == qH2.mass )
                          {
                            l.TakeOwnership( new ModContract3pt( l, tax, qH1, qH2, qSpectator, p,
                                                                 t, bHeavyAnti, j, deltaT, true ) );
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            p = -p;
          }
        }
      }
    }
    bFirstSpec = false;
  }
}

int main(int argc, char *argv[])
{
  int iReturn = EXIT_SUCCESS;
  #ifdef TEST_DEBUG_NEW
  std::cout << "main() :: before Debug()" << std::endl;
  if( Debug() ) return iReturn;
  #endif

  // See whether parameter file exists
  if( argc < 2 || !Common::FileExists( argv[1] ) )
  {
    std::cout << "Usage: xml3pt in_file [out_prefix]\n"
              << "where:\n"
              << "  in_file is the name of the Parameter.xml driving this job\n"
              << "  out_prefix+JobXmlPrefix (from .xml) is location for job and schedule files"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string sXmlFilename{ argv[1] };
  // 2nd parameter is optional. If present (and not a switch) it's a prefix for current run
  const std::string sRunPrefix( argc >= 3 && argv[2][0] != '-' ? argv[2] : "" );

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  try
  {
    const AppParams params( sXmlFilename, sRunPrefix );
    AppMaker x( params );
    x.Make();
    // Run or save the job
    x.application.saveParameterFile( params.Run.JobXmlPrefix + params.RunID() + ".xml" );
    if( params.Run.Run )
      x.application.run();
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }

  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return iReturn;
}
