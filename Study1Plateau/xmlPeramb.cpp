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
static const std::array<std::string, 2> RhoPhi{ "Rho", "Phi" };
static const std::string PerambPrefix{ "Peramb" + Sep };

static constexpr bool bMakePeramb{ false }; // Alternative is to load the perambulators
static constexpr bool bRhoPhi{ true }; // Should we make the rho and phi vectors

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

struct Contraction
{
  std::vector<std::string> sMomenta;
  
  explicit Contraction( int p )
  {
    // Make momenta
    for( int px = -p; px <= p; ++px )
      for( int py = -p; py <= p; ++py )
        for( int pz = -p; pz <= p; ++pz )
          if( px * px + py * py + pz * pz <= p * p)
            sMomenta.push_back( std::to_string(px) + Space + std::to_string(py) + Space + std::to_string(pz));
  }

  void Make( Application &application, const std::string &Left, const std::string &Right)
  {
    static const std::string MesonDir{ "mesons/C1/Distil/" };
    const std::string MesonFieldName{ Left + Sep + Right };
    MContraction::A2AMesonField::Par mfPar;
    mfPar.cacheBlock = 4;
    mfPar.block = 100;
    mfPar.left =  Left;
    mfPar.right = Right;
    mfPar.output = MesonDir + MesonFieldName;
    mfPar.gammas = "all";
    mfPar.mom = sMomenta;
    application.createModule<MContraction::A2AMesonField>(MesonFieldName, mfPar);
  }
};

int main(int argc, char *argv[])
{
  // We're not using C-style I/O
  std::ios_base::sync_with_stdio(false);
  if( !bMakePeramb && !bRhoPhi ) {
    std::cout << "Error: This would create a file that loads perambulators, but does nothing with them\n";
    return 1;
  }
  std::string RunID{ bMakePeramb ? "Make" : "Load" };
  if( bRhoPhi )
    RunID.append( "RhoPhi");

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
  globalPar.trajCounter.end      = globalPar.trajCounter.start + 1;
  globalPar.runId                = RunID;
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  application.setPar(globalPar);
  
  // gauge field
  MIO::LoadNersc::Par gaugePar;
  gaugePar.file = "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat";
  application.createModule<MIO::LoadNersc>(GaugeFieldName, gaugePar);
  // Stout smeared gauge field
  if( bMakePeramb )
  {
    MGauge::StoutSmearing::Par stoutPar;
    stoutPar.gauge = GaugeFieldName;
    stoutPar.steps = 3;
    stoutPar.rho = 0.1;
    application.createModule<MGauge::StoutSmearing>(GaugeFieldNameSmeared, stoutPar);
  }

  // Action and a solver for each quark
  std::array<std::string, NumQuarks> solverName;
  for( int i = 0; bMakePeramb && i < NumQuarks; ++i )
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
    const std::string Prefix{ PerambPrefix + Quarks[i].flavour };
    // Loop through all timeslices
    for (unsigned int t = 0; t < nt; t += 4) {
      const std::string st{ Sep + std::to_string( t ) };
      const std::string PerambName{ Prefix + st };
      if( bMakePeramb ) {
        MDistil::Perambulator::Par pPar;
        pPar.lapevec = LapEvecName;
        pPar.solver = solverName[i];
        pPar.noise = NoiseName;
        pPar.nvec = Nvec;
        pPar.Distil.nnoise = noisePar.nnoise;
        pPar.Distil.tsrc = t;
        pPar.UnsmearedSinkFileName = PerambName + UnsmearedSinkSuffix;
        application.createModule<MDistil::Perambulator>(PerambName, pPar);
      } else {
        static const std::string PerambDir{ "Peramb/C1/" };
        MIO::LoadPerambulator::Par pPar;
        pPar.nvec = Nvec;
        pPar.Distil.nnoise = noisePar.nnoise;
        pPar.Distil.tsrc = t;
        pPar.PerambFileName = PerambDir + PerambName;
        application.createModule<MIO::LoadPerambulator>(PerambName, pPar);
      }
      if( bRhoPhi ) {
        const std::string DVName{ "DV" + Sep + Quarks[i].flavour + st };
        MDistil::DistilVectors::Par dvPar;
        dvPar.noise = NoiseName;
        dvPar.perambulator = PerambName;
        dvPar.lapevec = LapEvecName;
        if( Quarks[i].flavour != "l" )
          dvPar.source = RhoPhi[0] + st;
        dvPar.sink = RhoPhi[1] + Sep + Quarks[i].flavour + st;
        dvPar.tsrc = t;
        application.createModule<MDistil::DistilVectors>(DVName, dvPar);
      }
    }
  }
  if( bRhoPhi ) {
    Contraction contraction{ 2 };

    // 16 x Rho_t_Rho_t
    for (unsigned int t = 0; t < nt; t += 4) {
      const std::string Left{ RhoPhi[0] + Sep + std::to_string(t) };
      contraction.Make( application, Left, Left );
    }

    // 64 x Phi_qL_t_Phi_qR_t
    for (unsigned int qL = 0; qL < NumQuarks; ++qL) {
      const std::string Left{ RhoPhi[1] + Sep + Quarks[qL].flavour + Sep };
      for (unsigned int qR = 0; qR < NumQuarks; ++qR) {
        const std::string Right{ RhoPhi[1] + Sep + Quarks[qR].flavour + Sep };
        for (unsigned int t = 0; t < nt; t += 4) {
          const std::string st{ std::to_string(t) };
          contraction.Make( application, Left + st, Right + st );
        }
      }
    }

    // 128 x Phi_q_dt_Rho_t
    for (unsigned int q = 0; q < NumQuarks; ++q) {
      const std::string LeftBase{ RhoPhi[1] + Sep + Quarks[q].flavour + Sep };
      for (unsigned int t = 0; t < nt; t += 4) {
        const std::string Right{ RhoPhi[0] + Sep + std::to_string(t) };
        for (unsigned int dt = 12; dt <= 24; dt += 4)
          contraction.Make( application, LeftBase + std::to_string((t+dt)%nt), Right );
      }
    }

    // 128 x Rho_dt_Phi_q_t
    for (unsigned int q = 0; q < NumQuarks; ++q) {
      const std::string RightBase{ RhoPhi[1] + Sep + Quarks[q].flavour + Sep };
      for (unsigned int t = 0; t < nt; t += 4) {
        const std::string Right{ RightBase + std::to_string(t) };
        for (unsigned int dt = 12; dt <= 24; dt += 4)
          contraction.Make( application, RhoPhi[0] + Sep + std::to_string((t+dt)%nt), Right );
      }
    }
  }

  // Save the parameter file
  application.saveParameterFile( RunID + ".xml" );

  // I only want to run this if lattice size has been specified and larger than default
  const Grid::Coordinate &lat{GridDefaultLatt()};
  if( lat.size() == 4 && lat[0] > 8 && lat[1] > 8 && lat[2] > 8 && lat[3] > 8 )
    application.run();
  
  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
