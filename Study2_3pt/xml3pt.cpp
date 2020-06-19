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

enum SourceT : int { Z2, GFPW };
static const std::array<std::string, 2> SourceTName{ "Z2", "GFPW" };

std::ostream& operator<<(std::ostream& os, const SourceT Type)
{
  if( Type >= 0 && Type < SourceTName.size() )
    os << SourceTName[Type];
  else
    os << "SourceT_Unknown_" << std::to_string( Type );
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

struct Quark: Serializable {
GRID_SERIALIZABLE_CLASS_MEMBERS(Quark,
                                std::string,  flavour,
                                double,       mass,
                                unsigned int, Ls,
                                double,       M5,
                                // Solver parameters
                                int,          maxIteration,
                                double,       residual,
                                bool,         GaugeSmear,
                                std::string,  EigenPackFilename )
public:
  Quark() = default;
  Quark(std::string flavour_, double mass_, unsigned int Ls_, double M5_, int maxIteration_,
        double residual_, bool GaugeSmear_ = false, const char * EigenPackFilename_ = nullptr)
  : flavour{flavour_}, mass{mass_}, Ls{Ls_}, M5{M5_}, maxIteration{maxIteration_},
    residual{residual_}, GaugeSmear{GaugeSmear_},
    EigenPackFilename{EigenPackFilename_ ? EigenPackFilename_ : "" } {}
};
/*
static const std::string GaugePathName{ "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat" };
static const char EigenPackC1_l[] = "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine";
static const Quark Quark_l {"l" , 0.005, 16, 1.8, 5000, 1e-8, false, EigenPackC1_l};
static const Quark Quark_s {"s" , 0.04, 16, 1.8, 5000, 1e-12};
static const Quark Quark_h0{"h0", 0.50, 12, 1.0, 5000, 1e-12, true};
static const Quark Quark_h1{"h1", 0.58, 12, 1.0, 5000, 1e-12, true};
static const Quark Quark_h2{"h2", 0.64, 12, 1.0, 5000, 1e-12, true};
static const Quark Quark_h3{"h3", 0.69, 12, 1.0, 5000, 1e-12, true};
*/
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
                                      std::string,  Lattice,
                                      std::string,  Gauge,
                                      bool,         HeavyQuark,
                                      bool,         HeavyAnti)
  };

  RunPar Run;
  SourceT Type;
  std::vector<Quark> HeavyQuarks;
  std::vector<Quark> SpectatorQuarks;
  std::vector<Common::Momentum> Momenta;
  std::vector<int> deltaT;
  inline int TimeBound( int t ) const
  { return t < 0 ? Run.Nt - ((-t) % Run.Nt) : t % Run.Nt; }
  AppParams( const std::string sXmlFilename );
  std::string ShortID() const;
  std::string RunID() const;
protected:
  static std::vector<Quark> ReadQuarks( XmlReader &r, const std::string &qType );
};

AppParams::AppParams( const std::string sXmlFilename )
{
  static const std::string sXmlTopLevel{ "Study2" };
  static const std::string sXmlTagName{ "RunPar" };
  XmlReader r( sXmlFilename, false, sXmlTopLevel );
  read( r, sXmlTagName, Run );
  HeavyQuarks     = ReadQuarks( r, "Heavy" );
  SpectatorQuarks = ReadQuarks( r, "Spectator" );
  Momenta = Common::ArrayFromString<Common::Momentum>( Run.Momenta );
  deltaT = Common::ArrayFromString<int>( Run.deltaT );
  // Check the type
  std::istringstream ss( Run.Type );
  if( ! ( ss >> Type && Common::StreamEmpty( ss ) ) )
    throw std::runtime_error( "Unrecognised type \"" + Run.Type + "\"" );
  // Check parameters make sense
  if( !Run.HeavyQuark && !Run.HeavyAnti )
    throw std::runtime_error( "At least one of HeavyQuark and HeavyAnti must be true" );
  if( Momenta.empty() )
    throw std::runtime_error( "There should be at least one momentum" );
  if( deltaT.empty() )
    throw std::runtime_error( "There should be at least one deltaT" );
}

std::vector<Quark> AppParams::ReadQuarks( XmlReader &r, const std::string &qType )
{
  static const std::string sQuark{ "Quark" };
  std::string TagName{ "Num" };
  TagName.append( qType );
  TagName.append( sQuark );
  int NumQuarks;
  read( r, TagName, NumQuarks );
  if( NumQuarks < 1 )
    throw std::runtime_error( "At least one " + qType + " quark must be specified" );
  TagName = qType;
  TagName.append( sQuark );
  const std::size_t PrefixLen{ TagName.length() };
  Quark q;
  std::vector<Quark> vq;
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
  for( const Quark &q : SpectatorQuarks )
    s << q.flavour;
  return s.str();
}

// Make a unique RunID string that completely describes the run
std::string AppParams::RunID() const
{
  std::ostringstream s;
  s << Type;
  for( const Quark &q : SpectatorQuarks )
    s << Sep << q.flavour;
  if( Run.HeavyQuark )
    s << Sep << "quark";
  if( Run.HeavyAnti )
    s << Sep << "anti";
  for( const Quark &q : HeavyQuarks )
    s << Sep << q.flavour;
  s << Sep << "t" << Sep << Run.Timeslices.start << Sep << Run.Timeslices.end << Sep << Run.Timeslices.step;
  s << Sep << "p";
  for( const Common::Momentum &p : Momenta )
    s << Sep << p.to_string( Sep );
  s << Sep << "dT";
  for( int dT : deltaT )
    s << Sep << std::to_string( dT );
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
  const SourceT Type;
  const Common::Momentum p;
  ModSink(const SourceT Type, const Common::Momentum &p);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(const SourceT type_, const Common::Momentum &p_) : Type{type_}, p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

bool ModSink::AddDependencies( HModList &ModList ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
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
    case GFPW:
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

class ModStoutGauge : public HMod
{
public:
  ModStoutGauge();
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModStoutGauge::ModStoutGauge()
{
  name = GaugeFieldName;
  Append( name, "stout" );
}

bool ModStoutGauge::AddDependencies( HModList &ModList ) const
{
  // Stout smeared field
  MGauge::StoutSmearing::Par stoutPar;
  stoutPar.gauge = GaugeFieldName;
  stoutPar.steps = 3;
  stoutPar.rho = 0.1;
  ModList.application.createModule<MGauge::StoutSmearing>(name, stoutPar);
  return true;
}

/**************************
 Action
**************************/

class ModAction : public HMod
{
public:
  static const std::string Prefix;
  const Quark &q;
  ModAction(const Quark &q);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModAction::Prefix{ "DWF" };

ModAction::ModAction(const Quark &q_) : q{q_}
{
  name = Prefix;
  Append( name, q.flavour );
}

bool ModAction::AddDependencies( HModList &ModList ) const
{
  MAction::DWF::Par actionPar;
  actionPar.gauge = q.GaugeSmear ? ModList.TakeOwnership(new ModStoutGauge()) : GaugeFieldName;
  actionPar.Ls    = q.Ls;
  actionPar.M5    = q.M5;
  actionPar.mass  = q.mass;
  // set fermion boundary conditions to be periodic space, antiperiodic time.
  static const std::string boundary{"1 1 1 -1"};
  static const std::string twist{"0. 0. 0. 0."};
  actionPar.boundary = boundary;
  actionPar.twist = twist;
  ModList.application.createModule<MAction::DWF>(name, actionPar);
  return true;
}

/**************************
 Solver
**************************/

class ModSolver : public HMod
{
public:
  static const std::string Prefix;
  const Quark &q;
  ModSolver(const Quark &q);
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModSolver::Prefix{ "CG" };

ModSolver::ModSolver(const Quark &q_) : q{q_}
{
  name = Prefix;
  Append( name, q.flavour );
}

bool ModSolver::AddDependencies( HModList &ModList ) const
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
    solverPar.eigenPack = "epack_" + q.flavour;
    ModList.application.createModule<MIO::LoadFermionEigenPack>(solverPar.eigenPack, epPar);
  }
  solverPar.action       = ModList.TakeOwnership( new ModAction( q ) );
  solverPar.residual     = q.residual;
  solverPar.maxIteration = ModList.params.Run.Gauge.empty() ? q.maxIteration / 10 : q.maxIteration;
  ModList.application.createModule<MSolver::RBPrecCG>(name, solverPar);
  return true;
}

/**************************
 Propagator
**************************/

class ModProp : public HMod
{
public:
  static const std::string Prefix;
  const SourceT Type;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModProp( const SourceT Type, const Quark &q, const Common::Momentum &p, int t );
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModProp::Prefix{ "Prop" };

ModProp::ModProp( const SourceT type_, const Quark &q_, const Common::Momentum &p_, int t_ )
: Type{type_}, q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, Type );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

bool ModProp::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par par;
  par.source = ModList.TakeOwnership( new ModSource( Type, p, t ) );
  par.solver = ModList.TakeOwnership( new ModSolver( q ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, par);
  return true;
}

/**************************
 Sequential Source
**************************/

class ModSourceSeq : public HMod
{
public:
  static const std::string Prefix;
  const SourceT Type;
  const int Current;
  const int deltaT;
  const Common::Momentum pSeq;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModSourceSeq( const SourceT Type, int Current, int deltaT, const Common::Momentum &pSeq,
                const Quark &q, const Common::Momentum &p, int t );
  virtual bool AddDependencies( HModList &ModList ) const;
protected:
  template<typename T> void AddDependenciesT( HModList &ModList ) const;
};

const std::string ModSourceSeq::Prefix{ ModSource::Prefix + "Seq" };

ModSourceSeq::ModSourceSeq( const SourceT type_, int current_, int deltaT_, const Common::Momentum &pSeq_,
                            const Quark &q_, const Common::Momentum &p_, int t_ )
: Type{type_}, Current{current_}, deltaT{deltaT_}, pSeq{pSeq_}, q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, Type );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

template<typename T> void ModSourceSeq::AddDependenciesT( HModList &ModList ) const
{
  typename T::Par seqPar;
  seqPar.q = ModList.TakeOwnership( new ModProp( Type, q, p, t ) );
  seqPar.tA  = ModList.params.TimeBound( t + deltaT );
  seqPar.tB  = seqPar.tA;
  seqPar.mom = pSeq.to_string4d( Space );
  seqPar.gamma = algInsert[Current];
  ModList.application.createModule<T>(name, seqPar);
}

bool ModSourceSeq::AddDependencies( HModList &ModList ) const
{
  //if( Type == SourceT::Z2 )
    AddDependenciesT<MSource::SeqGamma>( ModList );
  //else
    //AddDependenciesT<MSource::SeqGammaWall>( ModList );
  return true;
}

/**************************
 Sequential propagator
**************************/

class ModPropSeq : public HMod
{
public:
  static const std::string Prefix;
  const SourceT Type;
  const Quark &qSeq;
  const int Current;
  const int deltaT;
  const Common::Momentum pSeq;
  const Quark &q;
  const Common::Momentum p;
  const int t;
  ModPropSeq( const SourceT type, const Quark &qSeq, int current, int deltaT, const Common::Momentum &pSeq,
              const Quark &q, const Common::Momentum &p, int t );
  virtual bool AddDependencies( HModList &ModList ) const;
};

const std::string ModPropSeq::Prefix{ ModProp::Prefix + "Seq" };

ModPropSeq::ModPropSeq( const SourceT type_, const Quark &qSeq_, int current_, int deltaT_, const Common::Momentum &pSeq_,
                        const Quark &q_, const Common::Momentum &p_, int t_ )
: Type{type_}, qSeq{qSeq_}, Current{current_}, deltaT{deltaT_}, pSeq{pSeq_}, q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, Type );
  Append( name, qSeq.flavour );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

bool ModPropSeq::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par quarkPar;
  quarkPar.source = ModList.TakeOwnership( new ModSourceSeq( Type, Current, deltaT, pSeq, q, p, t ) );
  quarkPar.solver = ModList.TakeOwnership( new ModSolver( qSeq ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, quarkPar);
  return true;
}

/**************************
 2pt contraction
**************************/

const std::string ContractionPrefix{ "meson" };
const std::string ContractionBaseOutput{ "C1/meson/" };

class ModContract2pt : public HMod
{
protected:
  std::string FileName;
public:
  const SourceT Type;
  const Quark &q1;
  const Quark &q2;
  const Common::Momentum p;
  const int t;
  ModContract2pt(const SourceT Type, const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_);
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModContract2pt::ModContract2pt(const SourceT type_, const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_)
: Type{type_}, q1{q1_}, q2{q2_}, p{p_}, t{t_}
{
  std::string s{SourceTName[type_]};
  Append( s, q1.flavour );
  Append( s, q2.flavour );
  AppendPT( s, t, p );
  const std::string PrefixType{ "2pt" };
  FileName = ContractionBaseOutput + PrefixType + "/" + s;
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
  mesPar.q1 = ModList.TakeOwnership(new ModProp( Type, q1, p, t ));
  mesPar.q2 = ModList.TakeOwnership(new ModProp( Type, q2, p0, t ));
  mesPar.sink = ModList.TakeOwnership(new ModSink( Type, -p ));
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
  const Quark &qSnk;
  const Quark &qSrc;
  const Quark &qSpectator;
  const Common::Momentum p;
  const int t;
  const bool bHeavyAnti;
  const int Current;
  const int deltaT;
  ModContract3pt( const SourceT Type, const Quark &qSnk, const Quark &qSrc, const Quark &qSpectator,
                  const Common::Momentum &p, int t, bool bHeavyAnti, int Current, int deltaT );
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModContract3pt::ModContract3pt( const SourceT type_, const Quark &qSnk_, const Quark &qSrc_, const Quark &qSpectator_,
                                const Common::Momentum &p_, int t_, bool bHeavyAnti_, int Current_, int deltaT_)
: Type{type_}, qSnk{qSnk_}, qSrc{qSrc_}, qSpectator{qSpectator_},
  p{p_}, t{t_}, bHeavyAnti{bHeavyAnti_}, Current{Current_}, deltaT(deltaT_)
{
  std::string s{SourceTName[type_]};
  Append( s, bHeavyAnti ? "anti" : "quark" );
  Append( s, qSnk.flavour );
  Append( s, qSrc.flavour );
  Append( s, *algInsertName[Current_] );
  AppendDeltaT( s, deltaT );
  AppendPT( s, t, p );
  const std::string PrefixType{ "3pt" + Sep + qSpectator.flavour };
  FileName = ContractionBaseOutput + PrefixType + "/" + s;
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
  if( bInvertSeq )
  {
    par.q1 = ModList.TakeOwnership( new ModProp( Type, qSrc, p, t ) );
    par.q2 = ModList.TakeOwnership(new ModPropSeq( Type, qSnk, Current, deltaT, p0, qSpectator,p0, t ));
  }
  else
  {
    par.q1 = ModList.TakeOwnership(new ModPropSeq( Type, qSnk, Current, deltaT, p0, qSpectator, p, t ));
    par.q2 = ModList.TakeOwnership( new ModProp( Type, qSrc, p0, t ) );
  }
  par.sink = ModList.TakeOwnership( new ModSink( Type, -p ) );
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
  globalPar.saveSchedule = true;
  Application application( globalPar );
  // gauge field
  if( params.Run.Gauge.empty() )
  {
    application.createModule<MGauge::Random>( GaugeFieldName );
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
  for( const Quark &qSpectator : l.params.SpectatorQuarks )
  {
    for( unsigned int t = l.params.Run.Timeslices.start; t < l.params.Run.Timeslices.end; t += l.params.Run.Timeslices.step )
    {
      for( Common::Momentum p : l.params.Momenta )
      {
        for( int pDoNeg = 0; pDoNeg < ( p ? 2 : 1 ); ++pDoNeg )
        {
          for( const Quark &qH1 : l.params.HeavyQuarks )
          {
            for( const Quark &qH2 : l.params.HeavyQuarks )
            {
              for( int iHeavy  = l.params.Run.HeavyQuark ? 0 : 1;
                       iHeavy <= l.params.Run.HeavyAnti  ? 1 : 0; ++iHeavy )
              {
                const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
                if( !p || qH1.mass >= qH2.mass )
                {
                  // 2pt functions
                  //if( bHeavyAnti )
                  //{
                    l.TakeOwnership( new ModContract2pt( l.params.Type, qSpectator, qH1, p, t ) );
                    l.TakeOwnership( new ModContract2pt( l.params.Type, qSpectator, qH2, p, t ) );
                  //}
                  //else
                  //{
                    l.TakeOwnership( new ModContract2pt( l.params.Type, qH1, qSpectator, p, t ) );
                    l.TakeOwnership( new ModContract2pt( l.params.Type, qH2, qSpectator, p, t ) );
                  //}
                  for( int deltaT : l.params.deltaT )
                  {
                    for( int j = 0; j < NumInsert; j++ )
                    {
                      l.TakeOwnership( new ModContract3pt( l.params.Type, qH1, qH2, qSpectator, p,
                                                           t, bHeavyAnti, j, deltaT ) );
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
  const Grid::Coordinate &lat{GridDefaultLatt()};
  try
  {
    const AppParams params( sXmlFilename );
    AppMaker x( params );
    x.Make();
    // Is the lattice the same size as specified in the job?
    std::vector<int> Lattice = Common::ArrayFromString<int>( params.Run.Lattice );
    bool bRun = Lattice.size() == lat.size();
    for( int i = 0; bRun && i < Lattice.size(); ++i )
      if( lat[i] != Lattice[i] )
        bRun = false;
    // Run or save the job
    if( !bRun )
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