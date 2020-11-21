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

enum SourceT : int { Z2, GF };
static const std::array<std::string, 2> SourceTName{ "Z2", "GF" };

std::ostream& operator<<(std::ostream& os, const SourceT Type)
{
  if( Type >= 0 && Type < SourceTName.size() )
    os << SourceTName[Type];
  else
    os << "SourceTUnknown" << std::to_string( Type );
  return os;
}

std::istream& operator>>(std::istream& is, SourceT &Type)
{
  std::string s;
  if( is >> s )
  {
    int i;
    for(i = 0; i < SourceTName.size() && !Common::EqualIgnoreCase( s, SourceTName[i] ); i++)
      ;
    if( i < SourceTName.size() )
      Type = static_cast<SourceT>( i );
    else
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

enum SinkT : int { Point, Wall };
static const std::array<std::string, 2> SinkTName{ "point", "wall" };
//const std::string &sSinkPoint{ SinkTName[SinkT::Point] };
//const std::string &sSinkWall { SinkTName[SinkT::Wall ] };

std::ostream& operator<<(std::ostream& os, const SinkT Type)
{
  if( Type >= 0 && Type < SinkTName.size() )
    os << SinkTName[Type];
  else
    os << "SinkTUnknown" << std::to_string( Type );
  return os;
}

std::istream& operator>>(std::istream& is, SinkT &Type)
{
  std::string s;
  if( is >> s )
  {
    int i;
    for(i = 0; i < SinkTName.size() && !Common::EqualIgnoreCase( s, SinkTName[i] ); i++)
      ;
    if( i < SinkTName.size() )
      Type = static_cast<SinkT>( i );
    else
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

inline const std::string &FilePrefix( const SourceT Type, const SinkT Sink )
{
  if( Type == SourceT::Z2 )
  {
    if( Sink == SinkT::Point )
    {
      static const std::string s{ "Z2PP" };
      return s;
    }
    else if( Sink == SinkT::Wall )
    {
      static const std::string s{ "Z2WP" };
      return s;
    }
  }
  else if( Type == SourceT::GF )
  {
    if( Sink == SinkT::Point )
    {
      static const std::string s{ "GFPW" };
      return s;
    }
    else if( Sink == SinkT::Wall )
    {
      static const std::string s{ "GFWW" };
      return s;
    }
  }
  assert( 0 && "Unknown sink/source" );
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
protected:
  template<typename Action> void ActionCommon( Application &app, const std::string &name, typename Action::Par &Par ) const;
public:
  void CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const;
};

template<typename Action> void Quark::ActionCommon( Application &app, const std::string &name, typename Action::Par &Par ) const
{
  Par.Ls       = Ls;
  Par.M5       = M5;
  Par.mass     = mass;
  Par.boundary = boundary;
  Par.twist    = twist;
  app.createModule<Action>( name, Par );
}

void Quark::CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const
{
  if( scale <= 0 )
  {
    MAction::DWF::Par Par;
    Par.gauge = std::move( Gauge );
    ActionCommon<MAction::DWF>( app, name, Par );
  }
  else
  {
    MAction::ScaledDWF::Par Par;
    Par.gauge = std::move( Gauge );
    Par.scale = scale;
    ActionCommon<MAction::ScaledDWF>( app, name, Par );
  }
}

class QuarkType
{
public:
  const SourceT Type;
  const SinkT Sink;
  const Quark &q;
  static bool FixEverything;
  QuarkType( const SourceT type_, const SinkT sink_, const Quark &Q ) : Type{type_}, Sink{sink_}, q{Q}{}
  QuarkType( const QuarkType &qt ) : QuarkType( qt.Type, qt.Sink, qt.q ) {}
  inline bool GaugeFixed() const { return FixEverything || Type == SourceT::GF || Sink == SinkT::Wall; }
  inline std::string GaugeFixedName( const std::string &s ) const
  {
    std::string sAdjusted{ s };
    if( GaugeFixed() )
      sAdjusted.append( 1, 'f' );
    return sAdjusted;
  }
  inline std::string Flavour() const { return GaugeFixedName( q.flavour ); };
};

bool QuarkType::FixEverything{ false };

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

inline void Append( std::string &sDest, const std::string &s )
{
  if( !s.empty() )
  {
    sDest.append( Sep );
    sDest.append( s );
  }
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

inline void Append( std::string &sDest, SourceT Type )
{
  Append( sDest, SourceTName[Type] );
}

inline void Append( std::string &sDest, SinkT Type )
{
  Append( sDest, SinkTName[Type] );
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
        Grid::Hadrons::VirtualMachine::GeneticPar,  genetic,
                                   databaseOptions, dbOptions,
                                      std::string,  ScheduleFile,
                                      int,          Nt,
        Grid::Hadrons::Application::TrajRange,      Timeslices,
                                      std::string,  Type,
                                      std::string,  Momenta,
                                      std::string,  deltaT,
                                      std::string,  OutputBase,
                                      std::string,  Gauge,
             Grid::Hadrons::MGauge::GaugeFix::Par,  GaugeFix,
                       MGauge::StoutSmearing::Par,  StoutSmear,
                                      bool,         TwoPoint,
                                      bool,         HeavyQuark,
                                      bool,         HeavyAnti,
                                      bool,         SinkPoint,
                                      bool,         SinkWall,
                                      bool,         DoNegativeMomenta,
                                      bool,         Run)
  };

  RunPar Run;
  SourceT Type;
  std::vector<Quark> HeavyQuarks;
  std::vector<Quark> SpectatorQuarks;
  std::vector<Quark> SpectatorQuarks2pt;
  std::vector<Common::Momentum> Momenta;
  std::vector<int> deltaT;
  inline int TimeBound( int t ) const
  { return t < 0 ? Run.Nt - ((-t) % Run.Nt) : t % Run.Nt; }
  AppParams( const std::string sXmlFilename );
  std::string ShortID() const;
  std::string RunID() const;
protected:
  static std::vector<Quark> ReadQuarks( XmlReader &r, const std::string &qType, int Min = 1 );
};

AppParams::AppParams( const std::string sXmlFilename )
{
  static const std::string sXmlTopLevel{ "Study2" };
  static const std::string sXmlTagName{ "RunPar" };
  XmlReader r( sXmlFilename, false, sXmlTopLevel );
  read( r, sXmlTagName, Run );
  HeavyQuarks     = ReadQuarks( r, "Heavy" );
  SpectatorQuarks = ReadQuarks( r, "Spectator" );
  SpectatorQuarks2pt = ReadQuarks( r, "Spectator2pt", 0 );
  Momenta = Common::ArrayFromString<Common::Momentum>( Run.Momenta );
  // Check the type
  std::istringstream ss( Run.Type );
  if( ! ( ss >> Type && Common::StreamEmpty( ss ) ) )
    throw std::runtime_error( "Unrecognised type \"" + Run.Type + "\"" );
  // Check parameters make sense
  if( !( Run.TwoPoint || Run.HeavyQuark ||Run.HeavyAnti ) )
    throw std::runtime_error( "At least one must be true of: TwoPoint; HeavyQuark; or HeavyAnti" );
  if( !Run.SinkPoint && !Run.SinkWall )
    throw std::runtime_error( "At least one of SinkPoint and SinkWall must be true" );
  if( Momenta.empty() )
    throw std::runtime_error( "There should be at least one momentum" );
  if( Run.HeavyQuark || Run.HeavyAnti )
  {
    deltaT = Common::ArrayFromString<int>( Run.deltaT );
    if( deltaT.empty() )
      throw std::runtime_error( "There should be at least one deltaT" );
  }
}

std::vector<Quark> AppParams::ReadQuarks( XmlReader &r, const std::string &qType, int Min )
{
  static const std::string sQuark{ "Quark" };
  std::string TagName{ "Num" };
  TagName.append( qType );
  TagName.append( sQuark );
  int NumQuarks;
  read( r, TagName, NumQuarks );
  if( NumQuarks < Min )
    throw std::runtime_error( std::to_string( NumQuarks ) + Common::Space + qType + " quarks" );
  TagName = qType;
  TagName.append( sQuark );
  const std::size_t PrefixLen{ TagName.length() };
  Quark q;
  std::vector<Quark> vq;
  vq.reserve( NumQuarks );
  for( int i = 0; i < NumQuarks; ++i )
  {
    TagName.resize( PrefixLen );
    TagName.append( std::to_string( i ) );
    read( r, TagName, q );
    vq.push_back( q );
  }
  return vq;
}

// Make a brief name for the job, not necessarily unique, but likely to identify the job
// ... for use as a base filename
std::string AppParams::ShortID() const
{
  std::ostringstream s;
  s << "S2" << Type;
  // I originally defined the GF type as GFPW
  // Make sure the runID doesn't change, because it's the seed for random number generation
  if( Type == SourceT::GF )
    s << "PW";
  for( const Quark &q : SpectatorQuarks )
    s << q.flavour;
  return s.str();
}

// Make a unique RunID string that completely describes the run
std::string AppParams::RunID() const
{
  std::ostringstream s;
  s << Type;
  if( Run.SinkPoint )
    s << "P";
  if( Run.SinkWall )
    s << "W";
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
  if( Run.HeavyQuark || Run.HeavyAnti )
  {
    s << Sep << "dT";
    for( int dT : deltaT )
      s << Sep << std::to_string( dT );
  }
  return s.str();
}

/**************************
 Generic classes for a list of dependencies
 **************************/

class HModList;

class HMod
{
protected:
  std::string name;
public:
  inline const std::string &Name() const { return name; };
  HMod() { name.reserve( 80 ); }
  virtual ~HMod() = default;
  virtual bool AddDependencies( HModList &ModList ) const { return true; };
};

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

const std::string HModList::TakeOwnership( HMod *pHMod )
{
  assert( pHMod && "HModList::TakeOwnership() : Null pointer" );
  std::string ModName{ pHMod->Name() };
  auto it{ list.find( ModName ) };
  if( it == list.end() && pHMod->AddDependencies( *this ) )
    list.emplace( std::make_pair( ModName, std::unique_ptr<HMod>( pHMod ) ) );
  else
    delete pHMod;
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
  ModSink(const Common::Momentum &p);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(const Common::Momentum &p_) : p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

bool ModSink::AddDependencies( HModList &ModList ) const
{
  MSink::ScalarPoint::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
  return true;
}

/**************************
 Point sink
**************************/

class ModSinkSmear : public HMod
{
public:
  static const std::string Prefix;
  const Common::Momentum p;
  ModSinkSmear(const Common::Momentum &p);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModSinkSmear::Prefix{ "Sink_Smear" };

ModSinkSmear::ModSinkSmear(const Common::Momentum &p_) : p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

bool ModSinkSmear::AddDependencies( HModList &ModList ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::Point>(Name(), sinkPar);
  return true;
}

/**************************
 Source
**************************/

class ModSource : public HMod
{
public:
  static const std::string Prefix;
  const SourceT Type;
  const Common::Momentum p;
  const int t;
  ModSource(const SourceT Type, const Common::Momentum &p, int t);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModSource::Prefix{ "Source" };

ModSource::ModSource(const SourceT type_, const Common::Momentum &p_, int t_) : Type{type_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, SourceTName[Type] );
  AppendPT( name, t, p );
}

bool ModSource::AddDependencies( HModList &ModList ) const
{
  switch( Type )
  {
    case Z2:
      if( p )
      {
        MSource::MomentumPhase::Par par;
        par.src = ModList.TakeOwnership( new ModSource( Type, Common::Momentum(0,0,0), t ) );
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
    case GF:
      {
        MSource::Wall::Par par;
        par.tW = t;
        par.mom = p.to_string4d( Common::Space );
        ModList.application.createModule<MSource::Wall>(name, par);
      }
      break;
    default:
      return false;
  }
  return true;
}

/**************************
 Stout smeared gauge
**************************/

class ModGauge : public HMod
{
public:
  const bool bFixed;
  const bool bSmeared;
  ModGauge( bool bFixed, bool bSmeared );
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModGauge::ModGauge( bool bFixed_, bool bSmeared_ ) : bFixed{bFixed_}, bSmeared{ bSmeared_ }
{
  name = GaugeFieldName;
  if( bFixed )
    Append( name, "fixed" );
  if( bSmeared )
    Append( name, "stout" );
}

bool ModGauge::AddDependencies( HModList &ModList ) const
{
  if( bSmeared )
  {
    // Stout smeared field
    MGauge::StoutSmearing::Par stoutPar( ModList.params.Run.StoutSmear );
    stoutPar.gauge = bFixed ? ModList.TakeOwnership(new ModGauge( bFixed, false )) : GaugeFieldName;
    ModList.application.createModule<MGauge::StoutSmearing>(name, stoutPar);
  }
  else if( bFixed )
  {
    // Gauge-fixed
    MGauge::GaugeFix::Par Par( ModList.params.Run.GaugeFix );
    Par.gauge = GaugeFieldName;
    ModList.application.createModule<MGauge::GaugeFix>(name, Par);
  }
  return true;
}

/**************************
 Action
**************************/

class ModAction : public HMod
{
public:
  static const std::string Prefix;
  const QuarkType qt;
  ModAction(const QuarkType &qt);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModAction::Prefix{ "DWF" };

ModAction::ModAction(const QuarkType &qt_) : qt{qt_}
{
  name = Prefix;
  Append( name, qt.Flavour() );
}

bool ModAction::AddDependencies( HModList &ModList ) const
{
  std::string Gauge = ModList.TakeOwnership( new ModGauge( qt.GaugeFixed(), qt.q.GaugeSmear &&
                                                           !ModList.params.Run.Gauge.empty() ) );
  qt.q.CreateAction( ModList.application, name, std::move( Gauge ) );
  return true;
}

/**************************
 Solver
**************************/

class ModSolver : public HMod
{
public:
  static const std::string Prefix;
  const QuarkType qt;
  ModSolver( const QuarkType &qt );
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModSolver::Prefix{ "CG" };

ModSolver::ModSolver( const QuarkType &QT ) : qt( QT )
{
  name = Prefix;
  Append( name, qt.Flavour() );
}

bool ModSolver::AddDependencies( HModList &ModList ) const
{
  // solvers
  MSolver::RBPrecCG::Par solverPar;
  if( ModList.params.Run.Gauge.length() && qt.q.EigenPackFilename.length() )
  {
    // eigenpacks for deflation
    MIO::LoadFermionEigenPack::Par epPar;
    epPar.filestem = qt.q.EigenPackFilename;
    epPar.multiFile = false;
    epPar.size = 600;
    epPar.Ls = qt.q.Ls;
    if( qt.GaugeFixed() )
      epPar.gaugeXform = ModList.TakeOwnership( new ModGauge( true, false ) ) + "_xform";
    solverPar.eigenPack = "epack_" + qt.Flavour();
    ModList.application.createModule<MIO::LoadFermionEigenPack>(solverPar.eigenPack, epPar);
  }
  solverPar.action       = ModList.TakeOwnership( new ModAction( qt ) );
  solverPar.residual     = qt.q.residual;
  solverPar.maxIteration = qt.q.maxIteration;
  ModList.application.createModule<MSolver::RBPrecCG>(name, solverPar);
  return true;
}

/**************************
 Propagator
**************************/

class ModProp : public HMod
{
protected:
  std::string FileName;
public:
  static const std::string Prefix;
  static const std::string PrefixConserved;
  static const std::string PrefixConservedSep;
  const QuarkType &qt;
  const Common::Momentum p;
  const int t;
  const SinkT SinkSmearing;
  // bForcePointSink=true to make the underlying propagator for wall-sinks
  ModProp(const std::string &OutputBase, const QuarkType &qt, const Common::Momentum &p, int t, bool bForcePointSink=false);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModProp::Prefix{ "Prop" };
const std::string ModProp::PrefixConserved{ "Ward" };
const std::string ModProp::PrefixConservedSep{ PrefixConserved + Sep };

ModProp::ModProp(const std::string &OutputBase, const QuarkType &qt_, const Common::Momentum &p_, int t_,bool bForcePoint)
: qt{qt_}, p{p_}, t{t_}, SinkSmearing{ bForcePoint ? SinkT::Point : qt_.Sink }
{
  name = Prefix;
  Append( name, qt.Type );
  if( SinkSmearing != SinkT::Point )
    Append( name, SinkSmearing );
  Append( name, qt.Flavour() );
  AppendPT( name, t, p );
  if( SinkSmearing == SinkT::Point && !p )
  {
    FileName = OutputBase + PrefixConserved + "/";
    FileName.append( name );
  }
}

bool ModProp::AddDependencies( HModList &ModList ) const
{
  switch( SinkSmearing )
  {
    case SinkT::Point:
    {
      MFermion::GaugeProp::Par par;
      par.source = ModList.TakeOwnership( new ModSource( qt.Type, p, t ) );
      par.solver = ModList.TakeOwnership( new ModSolver( qt ) );
      ModList.application.createModule<MFermion::GaugeProp>(name, par);
      if( !FileName.empty() )
      {
        // Check residual mass for zero-momentum propagators
        MContraction::WardIdentity::Par WIP;
        WIP.q = name;
        WIP.action = ModList.TakeOwnership( new ModAction( qt ) );
        WIP.mass = qt.q.mass;
        WIP.test_axial = true;
        WIP.output = FileName;
        ModList.application.createModule<MContraction::WardIdentity>(PrefixConservedSep + name, WIP);
      }
    }
      break;
    case SinkT::Wall:
    {
      MSink::Smear::Par smearPar;
      smearPar.q = ModList.TakeOwnership( new ModProp( ModList.params.Run.OutputBase, qt, p, t, true ) );
      smearPar.sink = ModList.TakeOwnership( new ModSinkSmear( p0 ) );
      ModList.application.createModule<MSink::Smear>(name, smearPar);
    }
      break;
    default:
      assert( 0 && "Unknown Sink type" );
  }
  return true;
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
  const QuarkType &qt;
  const Common::Momentum p;
  const int t;
  ModSourceSeq( int Current, int deltaT, const Common::Momentum &pSeq,
                const QuarkType &qt, const Common::Momentum &p, int t );
  virtual bool AddDependencies( HModList &ModList ) const;
protected:
  template<typename T> void AddDependenciesT( HModList &ModList ) const;
};

const std::string ModSourceSeq::Prefix{ ModSource::Prefix + "Seq" };

ModSourceSeq::ModSourceSeq( int current_, int deltaT_, const Common::Momentum &pSeq_,
                            const QuarkType &qt_, const Common::Momentum &p_, int t_ )
: Current{current_}, deltaT{deltaT_}, pSeq{pSeq_}, qt{qt_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, qt.Type );
  Append( name, qt.Sink );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, qt.Flavour() );
  AppendPT( name, t, p );
}

template<typename T> void ModSourceSeq::AddDependenciesT( HModList &ModList ) const
{
  typename T::Par seqPar;
  seqPar.q = ModList.TakeOwnership( new ModProp( ModList.params.Run.OutputBase, qt, p, t, true ) );
  seqPar.tA  = ModList.params.TimeBound( t + deltaT );
  seqPar.tB  = seqPar.tA;
  seqPar.mom = pSeq.to_string4d( Space );
  seqPar.gamma = algInsert[Current];
  ModList.application.createModule<T>(name, seqPar);
}

bool ModSourceSeq::AddDependencies( HModList &ModList ) const
{
  switch( qt.Sink )
  {
    case SinkT::Point:
      AddDependenciesT<MSource::SeqGamma>( ModList );
      break;
    case SinkT::Wall:
      AddDependenciesT<MSource::SeqGammaWall>( ModList );
      break;
    default:
      assert( 0 && "Unknown sink type" );
  }
  return true;
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
  const QuarkType &qt;
  const Common::Momentum p;
  const int t;
  ModPropSeq( const Quark &qSeq, int current, int deltaT, const Common::Momentum &pSeq,
              const QuarkType &qt, const Common::Momentum &p, int t );
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModPropSeq::Prefix{ ModProp::Prefix + "Seq" };

ModPropSeq::ModPropSeq( const Quark &qSeq_, int current_, int deltaT_, const Common::Momentum &pSeq_,
                        const QuarkType &qt_, const Common::Momentum &p_, int t_ )
: qSeq{qSeq_},Current{current_},deltaT{deltaT_},pSeq{pSeq_},qt{qt_},p{p_},t{t_}
{
  name = Prefix;
  Append( name, qt.Type );
  Append( name, qt.Sink );
  Append( name, qt.GaugeFixedName( qSeq.flavour ) );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, qt.Flavour() );
  AppendPT( name, t, p );
}

bool ModPropSeq::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par quarkPar;
  quarkPar.source = ModList.TakeOwnership( new ModSourceSeq( Current, deltaT, pSeq, qt, p, t ) );
  quarkPar.solver = ModList.TakeOwnership( new ModSolver( QuarkType( qt.Type, qt.Sink, qSeq ) ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, quarkPar);
  return true;
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
  const SourceT Type;
  const SinkT Sink;
  const Quark &q1;
  const Quark &q2;
  const Common::Momentum p;
  const int t;
  ModContract2pt(const std::string &OutputBase, const SourceT Type, const SinkT Sink,
                 const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_);
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModContract2pt::ModContract2pt(const std::string &OutputBase, const SourceT type_, const SinkT sink_,
                               const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_)
: Type{type_}, Sink{sink_}, q1{q1_}, q2{q2_}, p{p_}, t{t_}
{
  std::string s{ FilePrefix( Type, Sink ) };
  Append( s, q1.flavour );
  Append( s, q2.flavour );
  AppendPT( s, t, p );
  const std::string PrefixType{ "2pt" };
  FileName = OutputBase + PrefixType + "/" + s;
  name = ContractionPrefix;
  Append( name, PrefixType );
  Append( name, s );
}

bool ModContract2pt::AddDependencies( HModList &ModList ) const
{
  MContraction::Meson::Par mesPar;
  mesPar.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  mesPar.gammas = "all";
  const Common::Momentum p_q2  { Sink == SinkT::Wall ? p  : p0 };
  const Common::Momentum p_sink{ Sink == SinkT::Wall ? p0 : -p };
  mesPar.q1 = ModList.TakeOwnership(new ModProp( ModList.params.Run.OutputBase, QuarkType( Type, Sink, q1 ), p, t ));
  mesPar.q2 = ModList.TakeOwnership(new ModProp( ModList.params.Run.OutputBase, QuarkType( Type, Sink, q2 ), p_q2, t ));
  mesPar.sink = ModList.TakeOwnership(new ModSink( p_sink ));
  ModList.application.createModule<MContraction::Meson>(name, mesPar);
  return true;
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
  const SourceT Type;
  const SinkT Sink;
  const Quark &qSnk;
  const Quark &qSrc;
  const Quark &qSpectator;
  const Common::Momentum p;
  const int t;
  const bool bHeavyAnti;
  const int Current;
  const int deltaT;
  const bool bSinkMom;
  ModContract3pt( const std::string &OutputBase, const SourceT Type, const SinkT Sink,
                  const Quark &qSnk, const Quark &qSrc, const Quark &qSpectator,
                  const Common::Momentum &p, int t, bool bHeavyAnti, int Current, int deltaT, bool bSinkMom );
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModContract3pt::ModContract3pt(const std::string &OutputBase, const SourceT type_, const SinkT sink_,
                               const Quark &qSnk_, const Quark &qSrc_, const Quark &qSpectator_,
                               const Common::Momentum &p_, int t_, bool bHeavyAnti_, int Current_, int deltaT_,
                               bool bSinkMom_ )
: Type{type_}, Sink{sink_}, qSnk{qSnk_}, qSrc{qSrc_}, qSpectator{qSpectator_},
  p{p_}, t{t_}, bHeavyAnti{bHeavyAnti_}, Current{Current_}, deltaT(deltaT_), bSinkMom{ bSinkMom_ }
{
  std::string s{ FilePrefix( Type, Sink ) };
  Append( s, bHeavyAnti ? "anti" : "quark" );
  Append( s, qSnk.flavour );
  Append( s, qSrc.flavour );
  Append( s, *algInsertName[Current_] );
  AppendDeltaT( s, deltaT );
  AppendPT( s, t, p );
  if( bSinkMom )
    Append( s, "snkmom" );
  const std::string PrefixType{ "3pt" + Sep + qSpectator.flavour };
  FileName = OutputBase + PrefixType + "/" + s;
  name = ContractionPrefix;
  Append( name, PrefixType );
  Append( name, s );
}

bool ModContract3pt::AddDependencies( HModList &ModList ) const
{
  MContraction::Meson::Par par;
  par.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  par.gammas = "all";
  const bool bInvertSeq{ !bHeavyAnti };
  const QuarkType qtSrc( Type, Sink, qSrc );
  const QuarkType qtSpec( Type, Sink, qSpectator );
  if( bInvertSeq )
  {
    par.q1 = ModList.TakeOwnership( new ModProp( ModList.params.Run.OutputBase, qtSrc, p, t, true ) );
    par.q2 = ModList.TakeOwnership(new ModPropSeq( qSnk, Current, deltaT, bSinkMom ?  p : p0, qtSpec, p0, t ) );
  }
  else
  {
    par.q1 = ModList.TakeOwnership(new ModPropSeq( qSnk, Current, deltaT, bSinkMom ? -p : p0, qtSpec, p , t ) );
    par.q2 = ModList.TakeOwnership( new ModProp( ModList.params.Run.OutputBase, qtSrc, p0, t, true ) );
  }
  par.sink = ModList.TakeOwnership( new ModSink( bSinkMom ? p0 : -p ) );
  ModList.application.createModule<MContraction::Meson>(name, par);
  return true;
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
  globalPar.runId        = params.ShortID();
  // database options
  globalPar.database.resultDb   = params.Run.dbOptions.resultDb;
  globalPar.database.makeStatDb = params.Run.dbOptions.makeStatDb;
  globalPar.database.applicationDb = params.Run.dbOptions.applicationDbPrefix + params.RunID() + ".db";
  globalPar.database.restoreModules = false;
  globalPar.database.restoreMemoryProfile = false;
  globalPar.database.restoreSchedule = false;
  globalPar.scheduleFile = globalPar.runId + ".sched";
  globalPar.saveSchedule = true;
  Application application( globalPar );
  // gauge field
  QuarkType::FixEverything = params.Run.SinkWall; // If I'm doing wall sink, gauge-fixing point sink saves inversions
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
    for( unsigned int t = l.params.Run.Timeslices.start; t < l.params.Run.Timeslices.end; t += l.params.Run.Timeslices.step )
    {
      for(int iSink = l.params.Run.SinkPoint ? 0 : 1; iSink <= l.params.Run.SinkWall ? 1  : 0; iSink++)
      {
        const SinkT Sink{ iSink ? SinkT::Wall : SinkT::Point };
        for( Common::Momentum p : l.params.Momenta )
        {
          for( int pDoNeg = 0; pDoNeg < ( ( l.params.Run.DoNegativeMomenta && p ) ? 2 : 1 ); ++pDoNeg )
          {
            for( const Quark &qH1 : l.params.HeavyQuarks )
            {
              if( bFirstSpec )
              {
                // Physical point spectators for 2pt functions only
                for( const Quark &qSpec2 : l.params.SpectatorQuarks2pt )
                {
                  l.TakeOwnership( new ModContract2pt( l.params.Run.OutputBase, l.params.Type, Sink, qSpec2, qH1, p, t ) );
                  l.TakeOwnership( new ModContract2pt( l.params.Run.OutputBase, l.params.Type, Sink, qH1, qSpec2, p, t ) );
                }
              }
              if( l.params.Run.TwoPoint )
              {
                l.TakeOwnership( new ModContract2pt( l.params.Run.OutputBase, l.params.Type, Sink,
                                                       qSpectator, qH1, p, t ) );
                l.TakeOwnership( new ModContract2pt( l.params.Run.OutputBase, l.params.Type, Sink,
                                                       qH1, qSpectator, p, t ) );
              }
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
                        l.TakeOwnership( new ModContract3pt( l.params.Run.OutputBase, l.params.Type, Sink,
                                                             qH1, qH2, qSpectator, p,
                                                             t, bHeavyAnti, j, deltaT, false ) );
                      }
                    }
                  }
                  if( p && qH1.mass == qH2.mass )
                  {
                    for( int deltaT : l.params.deltaT )
                    {
                      for( int j = 0; j < NumInsert; j++ )
                      {
                        l.TakeOwnership( new ModContract3pt( l.params.Run.OutputBase, l.params.Type, Sink,
                                                             qH1, qH2, qSpectator, p,
                                                             t, bHeavyAnti, j, deltaT, true ) );
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
  if( argc < 2 )
  {
    std::cout << "1st argument should be Parameter.xml filename" << std::endl;
    return EXIT_FAILURE;
  }
  const std::string sXmlFilename{ argv[1] };

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  try
  {
    const AppParams params( sXmlFilename );
    AppMaker x( params );
    x.Make();
    // Run or save the job
    if( params.Run.SinkWall && !params.Run.SinkPoint )
    {
      LOG(Warning) << "Creation of wall sinks necessitates creation of point sinks." << std::endl;
      LOG(Warning) << "Are you sure you don't want to save point sinks?" << std::endl;
    }
    if( !params.Run.Run )
      x.application.saveParameterFile( params.RunID() + ".xml" );
    else
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
