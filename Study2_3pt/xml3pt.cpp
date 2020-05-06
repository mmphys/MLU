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
static const std::string GaugeFieldNameSmeared{"gauge_stout"};

struct Quark
{
  std::string flavour;
  double mass;
  unsigned int Ls;
  double M5;
  // Solver parameters
  int maxIteration;
  double residual;
  const std::string & GaugeField;
  const char * EigenPackFilename;
};


static const Quark Quarks[] = {
  // {"l" , 0.005, 16, 1.8, 5000, 1e-8, GaugeFieldName, "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine"},
  {"s" , 0.04, 16, 1.8, 5000, 1e-12, GaugeFieldName, nullptr},
  {"h1", 0.58, 12, 1.0, 5000, 1e-12, GaugeFieldNameSmeared, nullptr}, // charm
  // {"h2", 0.64, 12, 1.0, 5000, 1e-12, GaugeFieldNameSmeared, nullptr},  // charm
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
  sDest.append( Sep );
  sDest.append( s );
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

inline void AppendP( std::string &sDest, const Common::Momentum &p, bool bNegative = false )
{
  Append( sDest, 'p', p.to_string( Sep, bNegative ) );
}

inline void AppendT( std::string &sDest, int t )
{
  Append( sDest, 't', t );
}

inline void AppendPT( std::string &sDest, int t, const Common::Momentum &p, bool bNegative = false )
{
  AppendP( sDest, p, bNegative );
  AppendT( sDest, t );
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
  inline const std::string & Name() const { return name; };
  HMod() = default;
  HMod( const std::string &Name ) : name{Name} {}
  virtual ~HMod() = default;
  virtual void AddDependencies( HModList &ModList, Application &application ) const {};
};

class HModList
{
protected:
  std::map<std::string,std::unique_ptr<HMod>> list;
public:
  void TakeOwnership( HMod *pHMod, Application &application );
};

void HModList::TakeOwnership( HMod *pHMod, Application &application )
{
  std::string Name{ pHMod->Name() };
  auto it{ list.find( Name ) };
  if( it == list.end() )
  {
    pHMod->AddDependencies( *this, application );
    list.emplace( std::make_pair( std::move( Name ), std::unique_ptr<HMod>( pHMod ) ) );
  }
  else
    delete pHMod;
}

/**************************
 Point sink
**************************/

class ModSink : public HMod
{
protected:
  static const std::string Prefix;
  Common::Momentum p;
  bool bNegative;
public:
  ModSink(const Common::Momentum &p, bool bNegative);
  virtual void AddDependencies( HModList &ModList, Application &application ) const;
};

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(const Common::Momentum &p_, bool bNegative_) : HMod( Prefix ), p{p_}, bNegative{bNegative_}
{
  AppendP( name, p, bNegative );
}

void ModSink::AddDependencies( HModList &ModList, Application &application ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string( Space, bNegative );
  application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
}

/**************************
 Sequential propagator
**************************/

class SeqProp : public HMod
{
protected:
  static const std::string Prefix;
  const HMod &q;
  int tA, tB;
  int Current;
  Common::Momentum p;
public:
  SeqProp(const HMod &q, int tA, int tB, int Current, Common::Momentum &p);
  virtual void AddDependencies( HModList &ModList, Application &application ) const;
};

const std::string SeqProp::Prefix{ "SeqProp" };

SeqProp::SeqProp(const HMod &q_, int ta_, int tb_, int current_, Common::Momentum &p_)
: q{q_}, tA{ta_}, tB{tb_}, Current{current_}, p{p_}
{
  name.reserve( 64 );
  name = Prefix;
  Append( name, 'a', tA );
  Append( name, 'b', tB );
  Append( name, q.Name() );
}

void SeqProp::AddDependencies( HModList &ModList, Application &application ) const
{
  MSource::SeqGamma::Par seqPar;
  seqPar.q   = q.Name();
  seqPar.tA  = tA;
  seqPar.tB  = tB;
  seqPar.mom = p.to_string4d( Space );
  seqPar.gamma = algInsert[Current];
  application.createModule<MSource::SeqGamma>(name, seqPar);
}

/**************************
 Generic contraction
**************************/

class ModContract : public HMod
{
protected:
  static const std::string Prefix;
  std::string qH;
  std::string qL;
  int Current;
  Common::Momentum p;
  int t;
  std::string FileName;
public:
  ModContract(const Quark &qH_, const Quark &qL_, Common::Momentum &p_, int t_,
              int Current_ = std::numeric_limits<int>::max());
  virtual void AddDependencies( HModList &ModList, Application &application ) const;
protected:
  void AddDependencies( HModList &ModList, Application &application, const std::string &output ) const;
};

const std::string ModContract::Prefix{ "meson" };

ModContract::ModContract(const Quark &qH_, const Quark &qL_, Common::Momentum &p_, int t_, int Current_)
: qH{qH_.flavour}, qL{qL_.flavour}, Current{Current_}, p{p_}, t{t_}
{
  FileName.reserve( 64 );
  FileName.append( qH );
  FileName.append( Sep );
  FileName.append( qL );
  if( Current_ != std::numeric_limits<int>::max() )
  {
    FileName.append( Sep );
    FileName.append( *algInsertName[Current_] );
  }
  AppendPT( FileName, t, p );
  name = Prefix;
  Append( name, FileName );
}

void ModContract::AddDependencies( HModList &ModList, Application &application, const std::string &output ) const
{
  // Make all the objects I depend on
  HMod * pSink = new ModSink( p, true );
  std::string Sink = pSink->Name();
  ModList.TakeOwnership(pSink, application);
  // Now make my object
  static const std::string baseOutput{ "mesons/C1/Z2/" };
  MContraction::Meson::Par mesPar;
  mesPar.output = baseOutput + output + FileName;
  mesPar.q1     = qH;
  mesPar.q2     = qL;
  mesPar.gammas = "all";
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  mesPar.sink   = Sink;
  application.createModule<MContraction::Meson>(name, mesPar);
}

void ModContract::AddDependencies( HModList &ModList, Application &application ) const
{
  AddDependencies( ModList, application, "2pt/" );
}

/**************************
 Three-point contraction
**************************/

class ModContractCurrent : public ModContract
{
public:
  ModContractCurrent(const Quark &qH_, const Quark &qL_, Common::Momentum &p_, int t_, int Current_)
  : ModContract(qH_, qL_, p_, t_, Current_) {}
  virtual void AddDependencies( HModList &ModList, Application &application ) const;
};

void ModContractCurrent::AddDependencies( HModList &ModList, Application &application ) const
{
  ModContract::AddDependencies( ModList, application, "3pt/" );
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
  
  LOG(Message)
    << std::boolalpha
    << "runId " << globalPar.runId
    << ", nt=" << std::to_string( nt )
    << ", Random gauge " << bRandom
    << ", Eigen packs " << bEigenpackEnable
    << std::endl;

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
  HModList l;
  for( unsigned int t = 0; t < nt; t+=4 )
  {
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      static const std::string sink{"sink_kitchen"};
      HMod * ptr = new ModContract( qh, ql, Momenta[p], t );
      l.TakeOwnership(ptr, application);
      for( int j = 0; j < NumInsert; j++ )
      {
        HMod * ptr = new ModContractCurrent( qh, ql, Momenta[p], t, j );
        l.TakeOwnership(ptr, application);
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
  
  // Stout smeared field
  MGauge::StoutSmearing::Par stoutPar;
  stoutPar.gauge = GaugeFieldName;
  stoutPar.steps = 3;
  stoutPar.rho = 0.1;
  application.createModule<MGauge::StoutSmearing>(GaugeFieldNameSmeared, stoutPar);
  
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
  
  // Action and a solver for each quark
  std::array<std::string, NumQuarks> solverName;
  for( int i = 0; i < NumQuarks; ++i )
  {
    const Quark & q{ Quarks[i] };
    // actions
    MAction::DWF::Par actionPar;
    actionPar.gauge = q.GaugeField;
    actionPar.Ls    = q.Ls;
    actionPar.M5    = q.M5;
    actionPar.mass  = q.mass;
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    static const std::string boundary{"1 1 1 -1"};
    static const std::string twist{"0. 0. 0. 0."};
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    const std::string ActionName{ "DWF" + Sep + q.flavour };
    application.createModule<MAction::DWF>(ActionName, actionPar);
    
    // solvers
    MSolver::RBPrecCG::Par solverPar;
    if( !bRandom && bEigenpackEnable && q.EigenPackFilename ) {
      // eigenpacks for deflation
      MIO::LoadFermionEigenPack::Par epPar;
      epPar.filestem = q.EigenPackFilename;
      epPar.multiFile = false;
      epPar.size = 600;
      epPar.Ls = q.Ls;
      const std::string epackObjname{ "epack_" + q.flavour };
      application.createModule<MIO::LoadFermionEigenPack>(epackObjname, epPar);
      solverPar.eigenPack = epackObjname;
    }
    solverPar.action       = ActionName;
    solverPar.residual     = q.residual;
    solverPar.maxIteration = q.maxIteration;
    solverName[i] = "CG" + Sep + q.flavour;
    application.createModule<MSolver::RBPrecCG>(solverName[i], solverPar);
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
  xml3pt x( Quarks[1], Quarks[0], {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 2, 0, 0 },
    //{ 1, 1, 0 },
    //{ 1, 1, 1 },
  } );
  x.bRandom = GridCmdOptionExists( argv, argv + argc, "-r" );
  x.bEigenpackEnable = !GridCmdOptionExists( argv, argv + argc, "-e" );
  x.nt = lat.size() < 4 ? 64 : static_cast<unsigned int>( lat[3] );

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
