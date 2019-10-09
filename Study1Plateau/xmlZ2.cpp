/*************************************************************************************
 
 Create XML for a 2-pt mass plateau study using Z2 wall sources.
 Copyright (C) 2019
 Author: Michael Marshall<Michael.Marshall@ed.ac.uk>
 
 Based on an original Grid test program:
 Grid physics library, www.github.com/paboyle/Grid
 Source file: Tests/Hadrons/Test_hadrons_meson_3pt.cc
 Copyright (C) 2015-2018
 Author: Antonin Portelli <antonin.portelli@me.com>
 
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

struct Momentum
{
  const int x;
  const int y;
  const int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
  inline bool isZero() const { return x==0 && y==0 && z==0; }
  std::string to_string( const std::string &separator = Sep, bool bNegative = false ) const {
    return std::to_string( bNegative ? -x : x ) + separator
    + std::to_string( bNegative ? -y : y ) + separator
    + std::to_string( bNegative ? -z : z );
  }
};

static const Momentum Momenta[] = {
  { 0, 0, 0 },
  { 1, 0, 0 },
  { 1, 1, 0 },
  { 1, 1, 1 },
  { 2, 0, 0 },
};
static constexpr int NumMomenta{ sizeof( Momenta ) / sizeof( Momenta[0] ) };

int main(int argc, char *argv[])
{
  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  LOG(Message) << "Grid initialized" << std::endl;
  
  // run setup ///////////////////////////////////////////////////////////////
  Application              application;
  const unsigned int nt{ 64 };
  //static constexpr int NumQuarks{ 4 };
  //const std::array<std::string, NumQuarks> flavour{"l", "s", "c1", "c2"};
  //const std::array<double, NumQuarks>      mass   {.005, .04, .58  , .64};
  //const unsigned int nt{static_cast<unsigned int>(GridDefaultLatt()[Tp])};
  
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = 3000;
  globalPar.trajCounter.end      = 4000;
  globalPar.trajCounter.step     = 40;
  globalPar.runId                = "Z2Plateau";
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  application.setPar(globalPar);
  
  // gauge field
  MIO::LoadNersc::Par gaugePar;
  gaugePar.file = "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat";
  application.createModule<MIO::LoadNersc>(GaugeFieldName, gaugePar);
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
    sinkPar.mom = Momenta[p].to_string( Space, true );
    SinkName[p] = "sink_p_" + Momenta[p].to_string( Sep, true );
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
    if( q.EigenPackFilename ) {
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
  for (unsigned int t = 0; t < nt; ++t) {
    const std::string timeSuffix{ Sep + "t" + std::to_string( t ) };
    const std::string sZ2{ "Z2" };
    // I will always need the momentum 0 Z_2 wall source for this timeslice
    const std::string Suffix0{ timeSuffix + Sep + "p" + Sep + Momentum( 0, 0, 0 ).to_string() };
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
      propName0[i] = "prop" + Sep + Quarks[i].flavour + Suffix0;
      MFermion::GaugeProp::Par quarkPar;
      quarkPar.solver = solverName[i];
      quarkPar.source = srcName0;
      application.createModule<MFermion::GaugeProp>(propName0[i], quarkPar);
    }
    // Loop through all momenta
    for( unsigned int p = 0; p < NumMomenta; ++p ) {
      const std::string momentumSuffix{ Sep + "p" + Sep + Momenta[p].to_string() };
      const std::string Suffix{ timeSuffix + momentumSuffix };

      // Z2 source with this momenta (doesn't need creating for momentum 0)
      const std::string srcName{ Momenta[p].isZero() ? srcName0 : sZ2 + Suffix };
      if( !Momenta[p].isZero() ) {
        MSource::MomentumPhase::Par z2Par;
        z2Par.src = srcName0;
        z2Par.mom = Momenta[p].to_string( Space ) + Space + "0";
        application.createModule<MSource::MomentumPhase>(srcName, z2Par);
      }
      // Make propagators
      std::array<std::string, NumQuarks> propName;
      for (unsigned int i = 0; i < NumQuarks; ++i)
      {
        if( Momenta[p].isZero() )
          propName[i] = propName0[i];
        else {
          propName[i] = "prop" + Sep + Quarks[i].flavour + Suffix;
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
          const std::string MesonSuffix{ propName[i] + Sep + propName0[j] + Sep + SinkName[p] };
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

  // execution
  application.saveParameterFile("Z2Plateau.xml");
#ifndef DEBUG
  application.run(); // I only want to run this on Tesseract
#endif
  
  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}