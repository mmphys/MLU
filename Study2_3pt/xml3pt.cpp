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

struct Quark
{
  std::string flavour;
  double mass;
  unsigned int Ls;
  double M5;
  // Solver parameters
  int maxIteration;
  double residual;
  const bool bGaugeSmear = false;
  const char * EigenPackFilename = nullptr;
};

static const std::string GaugePathName{ "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat" };
static const char EigenPackC1_l[] = "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine";
static const Quark Quark_l {"l" , 0.005, 16, 1.8, 5000, 1e-8, false, EigenPackC1_l};
static const Quark Quark_s {"s" , 0.04, 16, 1.8, 5000, 1e-12};
static const Quark Quark_h0{"h0", 0.50, 12, 1.0, 5000, 1e-12, true};
static const Quark Quark_h1{"h1", 0.58, 12, 1.0, 5000, 1e-12, true};
static const Quark Quark_h2{"h2", 0.64, 12, 1.0, 5000, 1e-12, true};
static const Quark Quark_h3{"h3", 0.69, 12, 1.0, 5000, 1e-12, true};

static const std::array<const Quark *, 1> SpectatorQuarks{ &Quark_s };
static const std::vector<const Quark *> HeavyQuarks{ &Quark_h0, &Quark_h1, &Quark_h2, &Quark_h3 };

static const std::vector<Common::Momentum> Momenta{
  { 0, 0, 0 },
  //{ 1, 0, 0 },
  //{ 1, 1, 0 },
  //{ 1, 1, 1 },
  //{ 2, 0, 0 },
};

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
  unsigned int nt;
  bool bEigenpackEnable;
  bool bRandom;
  AppParams( unsigned int nt_, bool bEigenpackEnable_, bool bRandom_ )
  : nt{nt_}, bEigenpackEnable{bEigenpackEnable_}, bRandom{bRandom_} {}
};

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
  actionPar.gauge = q.bGaugeSmear ? ModList.TakeOwnership(new ModStoutGauge()) : GaugeFieldName;
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
  if( ModList.params.bEigenpackEnable && q.EigenPackFilename )
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
  solverPar.maxIteration = ModList.params.bRandom ? q.maxIteration / 10 : q.maxIteration;
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

bool ModSourceSeq::AddDependencies( HModList &ModList ) const
{
  MSource::SeqGamma::Par seqPar;
  seqPar.q = ModList.TakeOwnership( new ModProp( Type, q, p, t ) );
  seqPar.tA  = ( t + deltaT ) % ModList.params.nt;
  seqPar.tB  = seqPar.tA;
  seqPar.mom = pSeq.to_string4d( Space );
  seqPar.gamma = algInsert[Current];
  ModList.application.createModule<MSource::SeqGamma>(name, seqPar);
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
const std::string ContractionBaseOutput{ "mesons/C1/" };

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
  std::string s{q1.flavour};
  Append( s, q2.flavour );
  AppendPT( s, t, p );
  const std::string PrefixType{ SourceTName[type_] + Sep + "2pt" };
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
  mesPar.sink = ModList.TakeOwnership(new ModSink( Type, -p ));
  mesPar.q1 = ModList.TakeOwnership(new ModProp( Type, q1, p, t ));
  mesPar.q2 = ModList.TakeOwnership(new ModProp( Type, q2, Common::Momentum(0,0,0), t ));
  ModList.application.createModule<MContraction::Meson>(name, mesPar);
  return true;
}

/**************************
 Three-point contraction
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
  const bool bCurrentAnti;
  const int Current;
  const int deltaT;
  ModContract3pt( const SourceT Type, const Quark &qSnk, const Quark &qSrc, const Quark &qSpectator,
                  const Common::Momentum &p, int t, bool bCurrentAnti, int Current, int deltaT );
  virtual bool AddDependencies( HModList &ModList ) const;
};

ModContract3pt::ModContract3pt( const SourceT type_, const Quark &qSnk_, const Quark &qSrc_, const Quark &qSpectator_,
                                const Common::Momentum &p_, int t_, bool bCurrentAnti_, int Current_, int deltaT_)
: Type{type_}, qSnk{qSnk_}, qSrc{qSrc_}, qSpectator{qSpectator_},
  p{p_}, t{t_}, bCurrentAnti{bCurrentAnti_}, Current{Current_}, deltaT(deltaT_)
{
  std::string s{bCurrentAnti ? "anti" : "quark" };
  Append( s, qSnk.flavour );
  Append( s, qSrc.flavour );
  Append( s, *algInsertName[Current_] );
  AppendDeltaT( s, deltaT );
  AppendPT( s, t, p );
  const std::string PrefixType{ SourceTName[type_] + Sep + "3pt" + Sep + qSpectator.flavour };
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
  par.sink = ModList.TakeOwnership(new ModSink( Type, p ));
  std::string qSeq{ ModList.TakeOwnership( new ModPropSeq( Type, qSnk, Current, deltaT, -p,
                                                           qSpectator, Common::Momentum(0,0,0), t ) ) };
  std::string q{ ModList.TakeOwnership( new ModProp( Type, qSrc, Common::Momentum(0,0,0), t ) ) };
  par.q1 = bCurrentAnti ? qSeq : q;
  par.q2 = bCurrentAnti ? q    : qSeq;
  ModList.application.createModule<MContraction::Meson>(name, par);
  return true;
}

/**************************
 Make application
**************************/

class AppMaker
{
public:
  const AppParams params;
  Application application;
  HModList l;
protected:
  Application Setup( const std::string &RunID );
public:
  explicit AppMaker( unsigned int nt_, bool bEigenpackEnable_, bool bRandom_, const std::string &RunID )
  : params(nt_, bEigenpackEnable_, bRandom_), application{Setup( RunID )}, l(application, params) {}
  void Make( SourceT Type, const Quark &qSpectator, const std::vector<Common::Momentum> &mom );
};

// One-time initialisation
Application AppMaker::Setup( const std::string &RunID )
{
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = 3000;
  globalPar.trajCounter.end      = 3001;
  globalPar.trajCounter.step     = 40;
  globalPar.runId                = RunID;
  globalPar.genetic.maxGen       = params.bRandom ? 10 : 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  globalPar.saveSchedule = false;
  Application application(globalPar);
  // gauge field
  if( params.bRandom )
  {
    application.createModule<MGauge::Random>(GaugeFieldName);
  }
  else
  {
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = GaugePathName;
    application.createModule<MIO::LoadNersc>(GaugeFieldName, gaugePar);
  }
  return application;
}

void AppMaker::Make( SourceT Type, const Quark &qSpectator, const std::vector<Common::Momentum> &mom )
{
  for( unsigned int t = 0; t < params.nt; t+=4 )
  //for( unsigned int t = 44; t < nt; t+=4 )
  //unsigned int t = 8;
  {
    for( const Quark *qSnk : HeavyQuarks )
    {
      for( const Quark *qSrc : HeavyQuarks )
      {
        for( const Common::Momentum &p : mom )
        {
          for( int qSpectatorAnti = 0; qSpectatorAnti < 2; ++qSpectatorAnti )
          {
            // 2pt functions
            if( qSpectatorAnti )
            {
              l.TakeOwnership( new ModContract2pt( Type, *qSrc, qSpectator, p, t ) );
              l.TakeOwnership( new ModContract2pt( Type, *qSnk, qSpectator, p, t ) );
            }
            else
            {
              l.TakeOwnership( new ModContract2pt( Type, qSpectator, *qSrc, p, t ) );
              l.TakeOwnership( new ModContract2pt( Type, qSpectator, *qSnk, p, t ) );
            }
            static const std::vector<int> deltaTList{ 12, 14, 16, 20 };
            for( int deltaT : deltaTList )
            {
              for( int j = 0; j < NumInsert; j++ )
              {
                l.TakeOwnership( new ModContract3pt( Type, *qSnk, *qSrc, qSpectator, p, t, qSpectatorAnti, j, deltaT ) );
              }
            }
          }
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

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  const Grid::Coordinate &lat{GridDefaultLatt()};
  const bool bRandom{ GridCmdOptionExists( argv, argv + argc, "-r" ) };
  const bool bEigenpackEnable{ !bRandom && !GridCmdOptionExists( argv, argv + argc, "-e" ) };
  const unsigned int nt{ lat.size() < 4 ? 64 : static_cast<unsigned int>( lat[3] ) };
  try
  {
    SourceT Type;
    {
      std::string s{ GridCmdOptionPayload( argv, argv + argc, "-t" ) };
      int i{ Z2 };
      GridCmdOptionInt( s, i );
      Type = static_cast<SourceT>( i );
    }
    if( Type < 0 || Type >= SourceTName.size() )
      throw std::runtime_error( "Unrecognised type " + std::to_string( Type ) );
    LOG(Message)
      << std::boolalpha
      << "nt=" << std::to_string( nt )
      << ", Random gauge " << bRandom
      << ", Eigen packs " << bEigenpackEnable
      << ", Type " << Type
      << std::endl;
    for( const Quark *qSpectator : SpectatorQuarks )
    {
      //for( int Type = 0; Type < SourceTName.size(); ++Type )
      {
        const std::string RunID{ "S2" + SourceTName[Type] + Sep + qSpectator->flavour };
        AppMaker x( nt, bEigenpackEnable, bRandom, RunID );
        x.Make( Type, *qSpectator, Momenta );
        // save
        static const std::string XmlFileName{ RunID + ".xml" };
        x.application.saveParameterFile( XmlFileName );
        // execute
        if( bRandom || ( lat.size() == 4 && lat[0] == 24 && lat[1] == 24 && lat[2] == 24 && lat[3] == 64 ) )
          x.application.run();
        else
          LOG(Warning) << "The parameters in " << XmlFileName << " are designed for --grid 24.24.24.64" << std::endl;
      }
    }
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
