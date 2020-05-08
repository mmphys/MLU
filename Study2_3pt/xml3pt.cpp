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

enum SourceT { Z2, GFPW };
static const std::array<std::string, 2> SourceTName{ "Z2", "GFPW" };

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


static const Quark Quarks[] = {
  // {"l" , 0.005, 16, 1.8, 5000, 1e-8, false, "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine"},
  {"s" , 0.04, 16, 1.8, 5000, 1e-12},
  {"h1", 0.58, 12, 1.0, 5000, 1e-12, true}, // charm
  // {"h2", 0.64, 12, 1.0, 5000, 1e-12, true},  // charm
};
static constexpr int NumQuarks{ sizeof( Quarks ) / sizeof( Quarks[0] ) };

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
  virtual void AddDependencies( HModList &ModList ) const {};
};

class HModList
{
protected:
  std::map<std::string,std::unique_ptr<HMod>> list;
public:
  // These are used by modules when adding dependencies
  Application &application;
  const int nt;
  bool bEigenpackEnable;
public:
  HModList( Application &application, int nt, bool bEigenpackEnable );
  const std::string TakeOwnership( HMod *pHMod );
};

HModList::HModList( Application &application_, int nt_, bool bEigenpackEnable_ )
: application{application_}, nt{nt_}, bEigenpackEnable{bEigenpackEnable_} {}

const std::string HModList::TakeOwnership( HMod *pHMod )
{
  assert( pHMod && "HModList::TakeOwnership() : Null pointer" );
  std::string ModName{ pHMod->Name() };
  auto it{ list.find( ModName ) };
  if( it == list.end() )
  {
    pHMod->AddDependencies( *this );
    list.emplace( std::make_pair( ModName, std::unique_ptr<HMod>( pHMod ) ) );
  }
  else
    delete pHMod;
  return ModName;
}

/**************************
 Point sink
**************************/

class ModSink : public HMod
{
protected:
  static const std::string Prefix;
  SourceT Type;
  Common::Momentum p;
public:
  ModSink(const SourceT Type, const Common::Momentum &p);
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(const SourceT type_, const Common::Momentum &p_) : Type{type_}, p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSink::AddDependencies( HModList &ModList ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
}

/**************************
 Source
**************************/

class ModSource : public HMod
{
protected:
  static const std::string Prefix;
  SourceT Type;
  int t;
public:
  ModSource(const SourceT Type, int t);
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModSource::Prefix{ "Source" };

ModSource::ModSource(const SourceT type_, int t_) : Type{type_}, t{t_}
{
  name = Prefix;
  Append( name, SourceTName[Type] );
  AppendT( name, t );
}

void ModSource::AddDependencies( HModList &ModList ) const
{
  MSource::Z2::Par z2Par;
  z2Par.tA = t;
  z2Par.tB = t;
  ModList.application.createModule<MSource::Z2>(name, z2Par);
}

/**************************
 Stout smeared gauge
**************************/

class ModStoutGauge : public HMod
{
public:
  ModStoutGauge();
  virtual void AddDependencies( HModList &ModList ) const;
};

ModStoutGauge::ModStoutGauge()
{
  name = GaugeFieldName;
  Append( name, "stout" );
}

void ModStoutGauge::AddDependencies( HModList &ModList ) const
{
  // Stout smeared field
  MGauge::StoutSmearing::Par stoutPar;
  stoutPar.gauge = GaugeFieldName;
  stoutPar.steps = 3;
  stoutPar.rho = 0.1;
  ModList.application.createModule<MGauge::StoutSmearing>(name, stoutPar);
}

/**************************
 Action
**************************/

class ModAction : public HMod
{
protected:
  static const std::string Prefix;
  const Quark &q;
public:
  ModAction(const Quark &q);
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModAction::Prefix{ "DWF" };

ModAction::ModAction(const Quark &q_) : q{q_}
{
  name = Prefix;
  Append( name, q.flavour );
}

void ModAction::AddDependencies( HModList &ModList ) const
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
}

/**************************
 Solver
**************************/

class ModSolver : public HMod
{
protected:
  static const std::string Prefix;
  const Quark &q;
public:
  ModSolver(const Quark &q);
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModSolver::Prefix{ "CG" };

ModSolver::ModSolver(const Quark &q_) : q{q_}
{
  name = Prefix;
  Append( name, q.flavour );
}

void ModSolver::AddDependencies( HModList &ModList ) const
{
  // solvers
  MSolver::RBPrecCG::Par solverPar;
  if( ModList.bEigenpackEnable && q.EigenPackFilename )
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
  solverPar.maxIteration = q.maxIteration;
  ModList.application.createModule<MSolver::RBPrecCG>(name, solverPar);
}

/**************************
 Sequential Source
**************************/

class ModSeqSource : public HMod
{
protected:
  static const std::string Prefix;
  SourceT Type;
  int Current;
  int deltaT;
  const Quark &q;
  Common::Momentum p;
  int t;
public:
  ModSeqSource(const SourceT Type, int Current, int deltaT, const Quark &q, const Common::Momentum &p, int t);
  virtual void AddDependencies( HModList &ModList ) const;
};

class ModProp : public HMod
{
protected:
  static const std::string Prefix;
  SourceT Type;
  const Quark &q;
  Common::Momentum p;
  int t;
  int Current;
  int deltaT;
public:
  ModProp(const SourceT Type, const Quark &q, const Common::Momentum &p, int t, int Current, int deltaT);
  ModProp(const SourceT Type, const Quark &q, const Common::Momentum &p, int t)
  : ModProp(Type, q, p, t, std::numeric_limits<int>::max(), 0) {}
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModSeqSource::Prefix{ "SeqSource" };

ModSeqSource::ModSeqSource(const SourceT type_, int current_, int deltaT_,
                           const Quark &q_, const Common::Momentum &p_, int t_)
: Type{type_}, Current{current_}, deltaT{deltaT_}, q{q_}, p{p_}, t{t_}
{
  name.reserve( 80 );
  name = Prefix;
  Append( name, Type );
  Append( name, *algInsertName[Current] );
  AppendDeltaT( name, deltaT );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

void ModSeqSource::AddDependencies( HModList &ModList ) const
{
  MSource::SeqGamma::Par seqPar;
  seqPar.q = ModList.TakeOwnership( new ModProp( Type, q, p, t ) );
  seqPar.tA  = ( t + deltaT ) % ModList.nt;
  seqPar.tB  = seqPar.tA;
  seqPar.mom = p.to_string4d( Space );
  seqPar.gamma = algInsert[Current];
  ModList.application.createModule<MSource::SeqGamma>(name, seqPar);
}

/**************************
 Propagator
**************************/

const std::string ModProp::Prefix{ "Prop" };

ModProp::ModProp(const SourceT type_, const Quark &q_, const Common::Momentum &p_, int t_, int current_, int deltaT_)
: Type{type_}, q{q_}, t{t_}, p{p_}, Current{current_}, deltaT(deltaT_)
{
  name.reserve( 64 );
  name = Prefix;
  Append( name, Type );
  Append( name, q.flavour );
  if( Current != std::numeric_limits<int>::max() )
  {
    Append( name, *algInsertName[Current] );
    AppendDeltaT( name, deltaT );
  }
  AppendPT( name, t, p );
}

void ModProp::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par quarkPar;
  quarkPar.source = ModList.TakeOwnership( ( Current != std::numeric_limits<int>::max() )
                    ? static_cast<HMod *>( new ModSeqSource( Type, Current, deltaT, q, p, t ) )
                    : static_cast<HMod *>( new ModSource( Type, t ) ) );
  quarkPar.solver = ModList.TakeOwnership( new ModSolver( q ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, quarkPar);
}

/**************************
 Generic contraction
**************************/

class ModContract : public HMod
{
protected:
  static const std::string Prefix;
  static const std::string PrefixType2pt;
  static const std::string PrefixType3pt;
  static const std::string baseOutput;
  const bool b3pt;
  const std::string PrefixType;
  const SourceT Type;
  const Quark &q1;
  const Quark &q2;
  const Common::Momentum p;
  const int t;
  const int Current;
  const int deltaT;
  std::string FileName;
public:
  ModContract(const SourceT Type, const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_, int Current_, int deltaT);
  ModContract(const SourceT Type, const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_)
  : ModContract(Type, q1_, q2_, p_, t_, std::numeric_limits<int>::max(), 0) {}
  virtual void AddDependencies( HModList &ModList ) const;
};

const std::string ModContract::Prefix{ "meson" };
const std::string ModContract::PrefixType2pt{ "2pt" };
const std::string ModContract::PrefixType3pt{ "3pt" };
const std::string ModContract::baseOutput{ "mesons/C1/" };

ModContract::ModContract(const SourceT type_, const Quark &q1_, const Quark &q2_,
                         const Common::Momentum &p_, int t_, int Current_, int deltaT_)
: b3pt{Current_ != std::numeric_limits<int>::max()},
  PrefixType{SourceTName[type_] + Sep + ( b3pt ? PrefixType3pt : PrefixType2pt ) },
  Type{type_}, q1{q1_}, q2{q2_}, p{p_}, t{t_}, Current{Current_}, deltaT(deltaT_)
{
  std::string s{q1.flavour};
  Append( s, q2.flavour );
  if( b3pt )
  {
    Append( s, *algInsertName[Current_] );
    AppendDeltaT( s, deltaT );
  }
  AppendPT( s, t, p );

  FileName = baseOutput + PrefixType + "/" + s;
  name = Prefix;
  Append( name, PrefixType );
  Append( name, s );
}

void ModContract::AddDependencies( HModList &ModList ) const
{
  MContraction::Meson::Par mesPar;
  mesPar.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  mesPar.gammas = "all";
  mesPar.sink = ModList.TakeOwnership(new ModSink( Type, -p ));
  mesPar.q1 = ModList.TakeOwnership(new ModProp( Type, q1, p, t ));
  mesPar.q2 = ModList.TakeOwnership(new ModProp( Type, q2, Common::Momentum(0,0,0), t ));
  ModList.application.createModule<MContraction::Meson>(name, mesPar);
}

/**************************
 Three-point contraction
**************************/

class ModContractCurrent : public ModContract
{
public:
  ModContractCurrent(const SourceT Type, const Quark &q1, const Quark &q2, const Common::Momentum &p, int t, int Current, int deltaT)
  : ModContract(Type, q1, q2, p, t, Current, deltaT) {}
  virtual void AddDependencies( HModList &ModList ) const;
};

void ModContractCurrent::AddDependencies( HModList &ModList ) const
{
  if( q1.mass == q2.mass )
  {
    LOG(Error) << "Mass of " << q1.flavour << " and " << q2.flavour << " are both "
      << q1.mass << ". Skipping 3pt contraction"
      << std::endl;
  }
  else
  {
    MContraction::Meson::Par mesPar;
    mesPar.output = FileName;
    //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
    mesPar.gammas = "all";
    mesPar.sink = ModList.TakeOwnership(new ModSink( Type, -p ));
    bool q1Heavy{ q1.mass >= q2.mass };
    mesPar.q1 = ModList.TakeOwnership(q1Heavy ? new ModProp( Type, q1, Common::Momentum(0,0,0), t )
                                      : new ModProp( Type, q1, p, t, Current, deltaT ));
    mesPar.q2 = ModList.TakeOwnership(q1Heavy ? new ModProp( Type, q2, p, t, Current, deltaT )
                                      : new ModProp( Type, q2, Common::Momentum(0,0,0), t ));
    ModList.application.createModule<MContraction::Meson>(name, mesPar);
  }
}

/**************************
 Make application
**************************/

class xml3pt
{
protected:
  static constexpr int iL{ 0 }; // Offset within internal structures indexed by quark
  static constexpr int iH{ 1 };
  const Quark &qh;
  const Quark &ql;
  const int p0; // This is the offset into Momenta for momentum 0
  std::vector<Common::Momentum> Momenta;
  int NumMomenta;
public:
  Application application;
protected:
  Application Setup();
  std::array<std::string, NumQuarks> solverName;
public:
  unsigned int nt;
  bool bRandom;
  bool bEigenpackEnable;
  xml3pt( const Quark &qH, const Quark &qL, const std::initializer_list<Common::Momentum> &mom )
  : qh{qH}, ql{qL}, p0{0}, Momenta{mom}, application{Setup()} {}
  void Make();
  void NewMake();
};

Application xml3pt::Setup()
{
  // Make sure momentum list starts with 0
  NumMomenta = static_cast<int>( Momenta.size() );
  if( NumMomenta == 0 )
  {
    Momenta.emplace_back( 0, 0, 0 );
    ++NumMomenta;
  }
  else if( Momenta[p0].p2() )
    throw std::runtime_error( "The first momentum must be 0" );

  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = 3000;
  globalPar.trajCounter.end      = 3001;
  globalPar.trajCounter.step     = 40;
  globalPar.runId                = qh.flavour + Sep + ql.flavour + Sep + "3pt";
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  globalPar.saveSchedule = false;
  Application application(globalPar);
  
  // gauge field
  if( bRandom )
  {
    application.createModule<MGauge::Unit>(GaugeFieldName);
  }
  else
  {
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat";
    application.createModule<MIO::LoadNersc>(GaugeFieldName, gaugePar);
  }

  return application;
}

void xml3pt::NewMake()
{
  HModList l( application, nt, bEigenpackEnable );
  for( unsigned int t = 0; t < 4/*nt*/; t+=4 )
  {
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      static const SourceT Type{ Z2 };
      l.TakeOwnership( new ModContract( Type, qh, ql, Momenta[p], t ) );
      for( int deltaT = 12; deltaT <= 20; deltaT += 4 )
      {
        for( int j = 0; j < NumInsert; j++ )
        {
          l.TakeOwnership( new ModContractCurrent( Type, qh, ql, Momenta[p], t, j, deltaT ) );
        }
      }
    }
  }
}

void xml3pt::Make()
{
  // run setup ///////////////////////////////////////////////////////////////
  //static constexpr int NumQuarks{ 4 };
  //const std::array<std::string, NumQuarks> flavour{"l", "s", "c1", "c2"};
  //const std::array<double, NumQuarks>      mass   {.005, .04, .58  , .64};
  //const unsigned int nt{static_cast<unsigned int>(GridDefaultLatt()[Tp])};
  
  // sink for each momentum (negative)
  std::vector<std::string> SinkPos( NumMomenta );
  std::vector<std::string> SinkNeg( NumMomenta );
  for( unsigned int p = 0; p < NumMomenta; p++ )
  {
    const std::string SinkPrefix{ "sink_p_" };
    MSink::Point::Par sinkPar;
    sinkPar.mom = Momenta[p].to_string( Space );
    SinkPos[p] = SinkPrefix + Momenta[p].to_string( Sep );
    application.createModule<MSink::ScalarPoint>(SinkPos[p], sinkPar);
    if( Momenta[p].p2() )
    {
      sinkPar.mom = Momenta[p].to_string( Space, true );
      SinkNeg[p] = SinkPrefix + Momenta[p].to_string( Sep, true );
      application.createModule<MSink::ScalarPoint>(SinkNeg[p], sinkPar);
    }
    else
      SinkNeg[p] = SinkPos[p];
  }
  
  // Loop through all timeslices
  static const std::string sPropPrefix{ "prop" + Sep };
  const unsigned int deltaT{ nt / 4 };
  for (unsigned int t = 0; t < 1/*nt*/; t+=4)
  {
    const std::string timeSuffix{ Sep + "t" + Sep + std::to_string( t ) };
    const std::string sZ2{ "Z2" };
    // I will always need the momentum 0 Z_2 wall source for this timeslice
    const std::string Suffix0{ timeSuffix + Sep + "p" + Sep + Common::Momentum( 0, 0, 0 ).to_string( Sep ) };
    const std::string srcName0{ sZ2 + Suffix0 };
    {
      MSource::Z2::Par z2Par;
      z2Par.tA = t;
      z2Par.tB = t;
      application.createModule<MSource::Z2>(srcName0, z2Par);
    }
    // I will always need the momentum 0 propagator for this timeslice
    std::array<std::string, NumQuarks> propName0;
    std::array<std::array<std::string, NumInsert>, NumQuarks> SeqSrcName0;
    std::array<std::array<std::string, NumInsert>, NumQuarks> SeqPropName0;
    for (unsigned int i = 0; i < NumQuarks; ++i)
    {
      
      propName0[i] = sPropPrefix + Quarks[i].flavour + Suffix0;
      MFermion::GaugeProp::Par quarkPar;
      quarkPar.solver = solverName[i];
      quarkPar.source = srcName0;
      application.createModule<MFermion::GaugeProp>(propName0[i], quarkPar);
      // I also need a sequential source and propagator for each quark
      MSource::SeqGamma::Par seqPar;
      seqPar.q   = propName0[i];
      seqPar.tA  = ( t + deltaT ) % nt;
      seqPar.tB  = seqPar.tA;
      seqPar.mom = Common::Momentum( 0, 0, 0 ).to_string4d( Space );
      for( int j = 0; j < NumInsert; j++ )
      {
        /*switch( j )
        {
          case 0:
            seqPar.mom = "0 0 0 1";
            break;
          case 1:
            seqPar.mom = Common::Momentum( 1, 0, 0 ).to_string4d( Space );
            break;
          case 2:
            seqPar.mom = Common::Momentum( 0, 1, 0 ).to_string4d( Space );
            break;
          case 3:
            seqPar.mom = Common::Momentum( 0, 0, 1 ).to_string4d( Space );
            break;
        }*/
        seqPar.gamma = algInsert[j];
        SeqSrcName0[i][j] = sZ2 + Sep + *(algInsertName[j]) + Sep + Quarks[i].flavour + Suffix0;
        application.createModule<MSource::SeqGamma>(SeqSrcName0[i][j], seqPar);
        quarkPar.source = SeqSrcName0[i][j];
        SeqPropName0[i][j] = sPropPrefix + SeqSrcName0[i][j];
        application.createModule<MFermion::GaugeProp>(SeqPropName0[i][j], quarkPar);
      }
    }
    // Loop through all momenta
    for( unsigned int p = 0; p < NumMomenta; ++p ) {
      const std::string momentumSuffix{ Sep + "p" + Sep + Momenta[p].to_string( Sep ) };
      const std::string Suffix{ momentumSuffix + timeSuffix };

      // Z2 source with this momenta (doesn't need creating for momentum 0)
      const std::string srcName{ Momenta[p] ? sZ2 + Suffix : srcName0 };
      if( Momenta[p] ) {
        MSource::MomentumPhase::Par z2Par;
        z2Par.src = srcName0;
        z2Par.mom = Momenta[p].to_string( Space ) + Space + "0";
        application.createModule<MSource::MomentumPhase>(srcName, z2Par);
      }
      // Make propagators
      std::array<std::string, NumQuarks> propName;
      for (unsigned int i = 0; i < NumQuarks; ++i)
      {
        if( !Momenta[p] )
          propName[i] = propName0[i];
        else {
          propName[i] = sPropPrefix + Quarks[i].flavour + Suffix;
          MFermion::GaugeProp::Par quarkPar;
          quarkPar.solver = solverName[i];
          quarkPar.source = srcName;
          application.createModule<MFermion::GaugeProp>(propName[i], quarkPar);
        }
      }

      // contractions
      MContraction::Meson::Par mesPar;
      std::string MesonSuffix;
      for (unsigned int i = 0; i < NumQuarks; ++i)
        for (unsigned int j = 0; j < NumQuarks; ++j)
        {
          static const std::string sMeson{ "meson" };
          static const std::string baseOutput{ "mesons/C1/Z2/" };
          //const std::string MesonSuffix{ propName[i] + Sep + propName0[j] + Sep + SinkNeg[p] };
          const std::string MesonSuffix{ Quarks[i].flavour + Sep + Quarks[j].flavour + Suffix };
          mesPar.output = baseOutput + "2pt/" + MesonSuffix;
          mesPar.q1     = propName[i];
          mesPar.q2     = propName0[j];
          mesPar.gammas = "all";
          //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
          mesPar.sink   = SinkNeg[p];
          application.createModule<MContraction::Meson>(sMeson + Sep + MesonSuffix, mesPar);
          for (unsigned int k = 0; k < NumInsert; ++k)
          {
            const std::string MesonSuffix3pt{Quarks[i].flavour+Sep+Quarks[j].flavour+Sep+*(algInsertName[k])+Suffix};
            mesPar.output = baseOutput + "3pt/" + MesonSuffix3pt;
            mesPar.q2     = SeqPropName0[j][k];
            application.createModule<MContraction::Meson>(sMeson + Sep + MesonSuffix3pt, mesPar);
          }
        }
    }
  }
}

#ifdef TEST_DEBUG_NEW
static thread_local bool bUseDevice{ false };

struct DebugNewHeader
{
  std::size_t szBufSize;
  bool bOnDevice;
};
static constexpr std::size_t DebugNewHeaderSize = ( sizeof( DebugNewHeader ) + 15 ) & (~15);

std::ostream& operator<<(std::ostream& os, const DebugNewHeader &h)
{
  return os << &h << ": " << h.szBufSize << " bytes of " << ( h.bOnDevice ? "device" : "host" ) << " memory";
}

void * DebugNew( std::size_t size )
{
  DebugNewHeader * pHeader = static_cast<DebugNewHeader *>( std::malloc( size + DebugNewHeaderSize ) );
  pHeader->szBufSize = size + DebugNewHeaderSize;
  pHeader->bOnDevice = bUseDevice;
  std::cout << "DebugNew( " << size << " ) = " << *pHeader << std::endl;
  return reinterpret_cast<char *>( pHeader ) + DebugNewHeaderSize;
}

void DebugDelete( void * mem )
{
  DebugNewHeader * pHeader = reinterpret_cast<DebugNewHeader *>( static_cast<char *>( mem ) - DebugNewHeaderSize );
  std::cout << "DebugDelete( " << mem << " ) = " << *pHeader << std::endl;
  std::free( pHeader );
}

void * operator new( std::size_t size ) { return DebugNew( size ); }
void operator delete( void * mem ) { DebugDelete( mem ); }

class CTest
{
public:
  using Scalar = Grid::ComplexD;
  using ET = Eigen::Tensor<Scalar, 6, Eigen::RowMajor>;
  std::vector<Scalar> v;
  //std::vector<Scalar,Grid::alignedAllocator<Scalar>> vGrid;
  std::vector<Scalar> vGrid;
  ET tensor;
  Eigen::TensorMap<ET> tMap;
  CTest();
  static void * operator new( std::size_t size ) { return DebugNew( size ); }
  static void operator delete( void * mem ) { DebugDelete( mem ); }
};

CTest::CTest() : v(2), vGrid(8*1*2*3*4*5), tensor(8,1,2,3,4,5), tMap(&vGrid[0],8,1,2,3,4,5)
{
  std::cout << "CTest::CTest" << std::endl;
}

bool Debug()
{
  CTest t1;
  bUseDevice = true;
  std::unique_ptr<CTest> pBingo{ new CTest };
  t1.tensor.resize(1,2,1,2,1,2);
  bUseDevice = false;
  t1.tensor.resize(2,2,2,2,2,2);
  return true;
}
#endif

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
  static const std::initializer_list<Common::Momentum> mom{
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 2, 0, 0 },
    //{ 1, 1, 0 },
    //{ 1, 1, 1 },
  };
  xml3pt x( Quarks[1], Quarks[0], mom );
  x.bRandom = GridCmdOptionExists( argv, argv + argc, "-r" );
  x.bEigenpackEnable = !x.bRandom && !GridCmdOptionExists( argv, argv + argc, "-e" );
  x.nt = lat.size() < 4 ? 64 : static_cast<unsigned int>( lat[3] );
  LOG(Message)
    << std::boolalpha
    << "nt=" << std::to_string( x.nt )
    << ", Random gauge " << x.bRandom
    << ", Eigen packs " << x.bEigenpackEnable
    << std::endl;

  try
  {
    x.NewMake();

    // execution
    static const std::string XmlFileName{ x.application.getPar().runId + ".xml" };
    x.application.saveParameterFile( XmlFileName );
    if( x.bRandom || ( lat.size() == 4 && lat[0] == 24 && lat[1] == 24 && lat[2] == 24 && lat[3] == 64 ) )
      x.application.run();
    else
      LOG(Warning) << "The parameters in " << XmlFileName << " are designed for --grid 24.24.24.64" << std::endl;
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
