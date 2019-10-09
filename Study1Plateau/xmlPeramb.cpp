/*************************************************************************************
 
 Create XML for a 2-pt mass plateau study using distillation.
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
static const std::string GaugeFieldName{ "gauge" };
static const std::string GaugeFieldNameSmeared{ "gauge" + Sep + "stout" };
static const std::string UnsmearedSinkSuffix{ Sep + "unsmeared" + Sep + "sink" };
static const std::array<std::string, 2> RhoPhi{ "rho", "phi" };
static const std::string PerambPrefix{ "Peramb" + Sep };

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

  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = 3000;
  globalPar.trajCounter.step     = 160;
  globalPar.trajCounter.end      = globalPar.trajCounter.start + globalPar.trajCounter.step + 1;
  globalPar.runId                = "HLPeramb";
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

  // Eigenvectors of the Laplacian
  static const int Nvec{ 60 };
  static const std::string LapEvecName{ "LapEvec" };
  MDistil::LapEvec::Par lePar;
  lePar.gauge = GaugeFieldName;
  lePar.Stout.steps = 3;
  lePar.Stout.rho = 0.2;
  lePar.Cheby.PolyOrder = 11;
  lePar.Cheby.alpha = 0.5;
  lePar.Cheby.beta = 12.5;
  lePar.Lanczos.Nvec = Nvec;
  lePar.Lanczos.Nk = Nvec * 4 / 3;
  lePar.Lanczos.Np = lePar.Lanczos.Nk - Nvec;
  lePar.Lanczos.MaxIt = 500;
  lePar.Lanczos.resid = 1e-8;
  lePar.Lanczos.IRLLog = 0;
  application.createModule<MDistil::LapEvec>( LapEvecName, lePar );

  static const std::string NoiseName{ "Noise" };
  MDistil::Noises::Par noisePar;
  noisePar.nnoise = 1;
  noisePar.nvec = Nvec;
  application.createModule<MDistil::Noises>( NoiseName, noisePar );

  // Loop through all quarks
  for (unsigned int i = 0; i < NumQuarks; ++i) {
    const std::string Prefix{ PerambPrefix + Quarks[i].flavour + Sep };
    // Loop through all timeslices
    for (unsigned int t = 0; t < nt; t += 4) {
      const std::string PerambName{ Prefix + std::to_string( t ) };
      MDistil::Perambulator::Par pPar;
      pPar.lapevec = LapEvecName;
      pPar.solver = solverName[i];
      pPar.noise = NoiseName;
      pPar.nvec = Nvec;
      pPar.Distil.nnoise = noisePar.nnoise;
      pPar.Distil.tsrc = t;
      pPar.UnsmearedSinkFileName = PerambName + UnsmearedSinkSuffix;
      application.createModule<MDistil::Perambulator>(PerambName, pPar);

      const std::string DVName{ "DV" + Sep + Quarks[i].flavour + Sep + std::to_string( t ) };
      const std::string SourceName{ RhoPhi[0] + Sep + Quarks[i].flavour + Sep + std::to_string( t ) };
      const std::string SinkName{ RhoPhi[1] + Sep + Quarks[i].flavour + Sep + std::to_string( t ) };
      MDistil::DistilVectors::Par dvPar;
      dvPar.noise = NoiseName;
      dvPar.perambulator = PerambName;
      dvPar.lapevec = LapEvecName;
      dvPar.source = SourceName;
      dvPar.sink = SinkName;
      dvPar.tsrc = t;
      application.createModule<MDistil::DistilVectors>(DVName, dvPar);
    }
  }

  // Make momenta
  std::vector<std::string> sMomenta;
  {
    static constexpr int p{ 2 };
    for( int px = -p; px <= p; ++px )
      for( int py = -p; py <= p; ++py )
        for( int pz = -p; pz <= p; ++pz )
          if( px * px + py * py + pz * pz <= p * p)
            sMomenta.push_back( std::to_string(px) + Space + std::to_string(py) + Space + std::to_string(pz));
  }

  // For each vector type
  for( const std::string &vec : RhoPhi ) {
    // Loop through all quarks on right / source
    for (unsigned int qR = 0; qR < NumQuarks; ++qR) {
      // Loop through all quarks on left / sink
      for (unsigned int qL = 0; qL < NumQuarks; ++qL) {
        // Loop through all timeslices
        for (unsigned int t = 0; t < nt; t += 4) {
          const std::string st{ Sep + std::to_string(t) };
          //const std::string PerambLeft { PerambPrefix + Quarks[qL].flavour + st + Sep + vec};
          //const std::string PerambRight{ PerambPrefix + Quarks[qR].flavour + st + Sep + vec};
          const std::string Left { vec + Sep + Quarks[qL].flavour + Sep + std::to_string(t) };
          const std::string Right{ vec + Sep + Quarks[qR].flavour + Sep + std::to_string(t) };
          static const std::string MesonDir{ "mesons/C1/Distil/" };
          const std::string MesonFieldName{ Left + Sep + Right };
          MContraction::A2AMesonField::Par mfPar;
          mfPar.cacheBlock = 4;
          mfPar.block = 100;
          mfPar.left =  /*Peramb*/Left;
          mfPar.right = /*Peramb*/Right;
          mfPar.output = MesonDir + MesonFieldName;
          mfPar.gammas = "all";
          mfPar.mom = sMomenta;
          application.createModule<MContraction::A2AMesonField>(MesonFieldName, mfPar);
        }
      }
    }
  }

  // execution
  application.saveParameterFile("Peramb.xml");
#ifndef DEBUG
  application.run(); // I only want to run this on Tesseract
#endif
  
  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
