/**
 
 Create XML for a 2-pt mass plateau study using Z2 wall sources.

 Source file: xmlZ2.cpp
 
 Based on an original Grid test program:
 Grid physics library, www.github.com/paboyle/Grid
 Source file: Tests/Hadrons/Test_hadrons_meson_3pt.cc
 Copyright (C) 2015-2018
 Author: Antonin Portelli <antonin.portelli@me.com>
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Semileptonic Data Generation (SemiLep).
 
 SemiLep is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 SemiLep is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <SemiLepConfig.h>
#include <MLU/MLU.hpp>
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
  {"l", 0.005, 16, 1.8, 5000, 1e-8, GaugeFieldName,
    "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine"},
  // {"s", 0.04, 16, 1.8, 5000, 1e-12, GaugeFieldName, nullptr},
  {"h1", 0.58, 12, 1.0, 5000, 1e-12, GaugeFieldNameSmeared, nullptr}, // charm
  // {"h2", 0.64, 12, 1.0, 5000, 1e-12, GaugeFieldNameSmeared, nullptr},  // charm
};
static constexpr int NumQuarks{ sizeof( Quarks ) / sizeof( Quarks[0] ) };

static const MLU::Momentum Momenta[] = {
  { 0, 0, 0 },
  { 1, 0, 0 },
  { 1, 1, 0 },
  { 1, 1, 1 },
  { 2, 0, 0 },
};
static constexpr int NumMomenta{ sizeof( Momenta ) / sizeof( Momenta[0] ) };

void Make( Application &application, unsigned int nt, bool bRandom, bool bEigenpackEnable )
{
  // run setup ///////////////////////////////////////////////////////////////
  //static constexpr int NumQuarks{ 4 };
  //const std::array<std::string, NumQuarks> flavour{"l", "s", "c1", "c2"};
  //const std::array<double, NumQuarks>      mass   {.005, .04, .58  , .64};
  //const unsigned int nt{static_cast<unsigned int>(GridDefaultLatt()[Tp])};
  
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = 3000;
  globalPar.trajCounter.end      = 3001;
  globalPar.trajCounter.step     = 40;
  globalPar.runId                = "Z2Plateau";
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  globalPar.saveSchedule = false;
  application.setPar(globalPar);
  
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

  // Stout smeared field
  MGauge::StoutSmearing::Par stoutPar;
  stoutPar.gauge = GaugeFieldName;
  stoutPar.steps = 3;
  stoutPar.rho = 0.1;
  application.createModule<MGauge::StoutSmearing>(GaugeFieldNameSmeared, stoutPar);
  
  // sink for each momentum (negative)
  std::array<std::string, NumMomenta> SinkName;
  for( unsigned int p = 0; p < NumMomenta; p++ ) {
    MSink::Point::Par sinkPar;
    sinkPar.mom = Momenta[p].to_string3d( Space, true );
    SinkName[p] = "sink_p_" + Momenta[p].to_string3d( Sep, true );
    application.createModule<MSink::ScalarPoint>(SinkName[p], sinkPar);
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
#ifdef MLU_HADRONS_HAS_GUESSERS
      epPar.redBlack = true;
#endif
      const std::string epackObjname{ "epack_" + q.flavour };
      application.createModule<MIO::LoadFermionEigenPack>(epackObjname, epPar);
#ifdef MLU_HADRONS_HAS_GUESSERS
      const std::string GuesserName{ "guesser_" + q.flavour };
      MGuesser::ExactDeflationPar gPar;
      gPar.eigenPack = epackObjname;
      gPar.size = epPar.size;
      application.createModule<MGuesser::ExactDeflation>(GuesserName, gPar);
      solverPar.guesser = GuesserName;
#else
      solverPar.eigenPack = epackObjname;
#endif
    }
    solverPar.action       = ActionName;
    solverPar.residual     = q.residual;
    solverPar.maxIteration = q.maxIteration;
    solverName[i] = "CG" + Sep + q.flavour;
    application.createModule<MSolver::RBPrecCG>(solverName[i], solverPar);
  }

  // Loop through all timeslices
  static const std::string sPropPrefix{ "prop" + Sep };
  for (unsigned int t = 0; t < nt; t+=4) {
    const std::string timeSuffix{ Sep + "t" + Sep + std::to_string( t ) };
    const std::string sZ2{ "Z2" };
    // I will always need the momentum 0 Z_2 wall source for this timeslice
    const std::string Suffix0{ timeSuffix + Sep + "p" + Sep + MLU::p0.to_string3d( Sep ) };
    const std::string srcName0{ sZ2 + Suffix0 };
    {
      MSource::Z2::Par z2Par;
      z2Par.tA = t;
      z2Par.tB = t;
      application.createModule<MSource::Z2>(srcName0, z2Par);
    }
    // I will always need the momentum 0 propagator for this timeslice
    std::array<std::string, NumQuarks> propName0;
    for (unsigned int i = 0; i < NumQuarks; ++i)
    {
      propName0[i] = sPropPrefix + Quarks[i].flavour + Suffix0;
      MFermion::GaugeProp::Par quarkPar;
      quarkPar.solver = solverName[i];
      quarkPar.source = srcName0;
      application.createModule<MFermion::GaugeProp>(propName0[i], quarkPar);
    }
    // Loop through all momenta
    for( unsigned int p = 0; p < NumMomenta; ++p ) {
      const std::string momentumSuffix{ Sep + "p" + Sep + Momenta[p].to_string3d( Sep ) };
      const std::string Suffix{ momentumSuffix + timeSuffix };

      // Z2 source with this momenta (doesn't need creating for momentum 0)
      const std::string srcName{ Momenta[p] ? sZ2 + Suffix : srcName0 };
      if( Momenta[p] ) {
        MSource::MomentumPhase::Par z2Par;
        z2Par.src = srcName0;
        z2Par.mom = Momenta[p].to_string3d( Space ) + Space + "0";
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
          //const std::string MesonSuffix{ propName[i] + Sep + propName0[j] + Sep + SinkName[p] };
          const std::string MesonSuffix{ Quarks[i].flavour + Sep + Quarks[j].flavour + Suffix };
          mesPar.output = "mesons/C1/Z2/" + MesonSuffix;
          mesPar.q1     = propName[i];
          mesPar.q2     = propName0[j];
          mesPar.gammas = "all";
          //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
          mesPar.sink   = SinkName[p];
          application.createModule<MContraction::Meson>("meson" + Sep + MesonSuffix, mesPar);
        }
    }
  }
}

int main(int argc, char *argv[])
{
  int iReturn = EXIT_SUCCESS;
  std::cout << "main() :: before Debug()" << std::endl;

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());

  const Grid::Coordinate &lat{GridDefaultLatt()};
  const unsigned int nt{ lat.size() < 4 ? 64 : static_cast<unsigned int>( lat[3] ) };
  const bool bRandom{ GridCmdOptionExists( argv, argv + argc, "-r" ) };
  const bool bEigenpackEnable{ !GridCmdOptionExists( argv, argv + argc, "-e" ) };
  LOG(Message)
    << std::boolalpha
    << "nt=" << std::to_string( nt )
    << ", Random gauge " << bRandom
    << ", Eigen packs " << bEigenpackEnable
    << std::endl;

  try
  {
    Application application;
    Make( application, nt, bRandom, bEigenpackEnable );

    // execution
    static const std::string XmlFileName{ "Z2.template.xml" };
    application.saveParameterFile( XmlFileName );
    if( bRandom || ( lat.size() == 4 && lat[0] == 24 && lat[1] == 24 && lat[2] == 24 && lat[3] == 64 ) )
      application.run();
    else
      LOG(Warning) << "The parameters in " << XmlFileName << " are designed for --grid 24.24.24.64\nOn 16 nodes each config takes about 40 hours on Tesseract" << std::endl;
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
