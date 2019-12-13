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

#include "../Analyze/Common.hpp"
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

static const std::string RunName{ "GFWall" };

using namespace Grid;
using namespace Hadrons;

static const std::string Sep{ "_" };    // used inside filenames
static const std::string SepBig{ "." };    // used inside filenames
static const std::string Space{ " " };  // whitespace as field separator / human readable info

// Should I stout smear the heavy quark gauge field?
#define SMEAR_HEAVY

static const std::string GaugeFieldNameUnfixed{"gauge"};
static const std::string GaugeFieldName{"gauge_fixed"};
#ifdef SMEAR_HEAVY
static const std::string GaugeFieldNameUnfixedSmeared{"gauge_stout"};
static const std::string GaugeFieldNameSmeared{"gauge_stout_fixed"};
#define HEAVY_GAUGE_NAME GaugeFieldNameSmeared
#else
#define HEAVY_GAUGE_NAME GaugeFieldName
#endif
//static const std::string GaugeFieldNameSingle{"gauge_fixed_float"};
//static const std::string GaugeFieldNameSingleSmeared{"gauge_fixed_stout_float"};

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
  {"l", 0.005, 16, 1.8, 5000, 1e-8, GaugeFieldName, \
    "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine"},
  // {"s", 0.04, 16, 1.8, 5000, 1e-12, GaugeFieldName, nullptr},
  {"h1", 0.58, 12, 1.0, 5000, 1e-12, HEAVY_GAUGE_NAME, nullptr}, // charm
  // {"h2", 0.64, 12, 1.0, 5000, 1e-12, HEAVY_GAUGE_NAME, nullptr},  // charm
};
static constexpr int NumQuarks{ sizeof( Quarks ) / sizeof( Quarks[0] ) };

static const Common::Momentum Momentum0{ 0, 0, 0 };

static const Common::Momentum Momenta[] = {
  { 0, 0, 0 },
  { 1, 0, 0 },
  { 0, 1, 0 },
  { 0, 0, 1 },
  { 1, 1, 0 },
  { 1, 0, 1 },
  { 0, 1, 1 },
  { 1, 1, 1 },
  { 2, 0, 0 },
  { 0, 2, 0 },
  { 0, 0, 2 },
};
static constexpr int NumMomenta{ sizeof( Momenta ) / sizeof( Momenta[0] ) };

void CreateApp( Application &application )
{
  const unsigned int nt{ 64 };
  //static constexpr int NumQuarks{ 4 };
  //const std::array<std::string, NumQuarks> flavour{"l", "s", "c1", "c2"};
  //const std::array<double, NumQuarks>      mass   {.005, .04, .58  , .64};
  //const unsigned int nt{static_cast<unsigned int>(GridDefaultLatt()[Tp])};
  
  // gauge field
  MIO::LoadNersc::Par gaugePar;
  gaugePar.file = "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat";
  application.createModule<MIO::LoadNersc>(GaugeFieldNameUnfixed, gaugePar);

  // Should I stout smear the heavy quark gauge field?
#ifdef SMEAR_HEAVY
  MGauge::StoutSmearing::Par stoutPar;
  stoutPar.gauge = GaugeFieldNameUnfixed;
  stoutPar.steps = 3;
  stoutPar.rho = 0.1;
  application.createModule<MGauge::StoutSmearing>(GaugeFieldNameUnfixedSmeared, stoutPar);
#endif

  // Gauge fixed versions
  MGauge::GaugeFix::Par gfPar;
  gfPar.alpha = 0.05;
  gfPar.maxiter = 1e6;
  gfPar.Omega_tol = 1e-8;
  gfPar.Phi_tol = 1e-8;
  gfPar.gaugeFix = MGauge::Fix::coulomb;
  gfPar.Fourier = true;
  gfPar.gauge = GaugeFieldNameUnfixed;
  application.createModule<MGauge::GaugeFix>(GaugeFieldName, gfPar);
#ifdef SMEAR_HEAVY
  gfPar.gauge = GaugeFieldNameUnfixedSmeared;
  application.createModule<MGauge::GaugeFix>(GaugeFieldNameSmeared, gfPar);
#endif

  // Single-precision versions
  //MUtilities::GaugeSinglePrecisionCast::Par g1Par;
  //g1Par.field = GaugeFieldName;
  //application.createModule<MUtilities::GaugeSinglePrecisionCast>(GaugeFieldNameSingle, g1Par);
  //g1Par.field = GaugeFieldNameSmeared;
  //application.createModule<MUtilities::GaugeSinglePrecisionCast>(GaugeFieldNameSingleSmeared, g1Par);

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

  // The contraction sink will always have zero momentum
  static const std::string ContractionSinkName0{ "contraction_sink_p_" + Momentum0.to_string( Sep ) };
  {
    MSink::ScalarPoint::Par sinkPar;
    sinkPar.mom = Momentum0.to_string( Space );
    application.createModule<MSink::ScalarPoint>(ContractionSinkName0, sinkPar);
  }

  // Make a propagator sink for each of the non-zero momentum
  static std::array<std::string,NumMomenta> PropSinkName;
  {
    MSink::Point::Par sinkPar;
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      if( Momenta[p] )
      {
        PropSinkName[p] = "prop_sink_p_" + Momenta[p].to_string( Sep );
        sinkPar.mom = Momenta[p].to_string( Space );
        application.createModule<MSink::Point>(PropSinkName[p], sinkPar);
      }
    }
  }
  
  // Loop through all timeslices
  for (unsigned int t = 0; t < nt; t += 4 )
  {
    const std::string TimeSuffix{ Sep + "t" + Sep + std::to_string( t ) };

    // Zero-momentum wall-source with this momenta
    const std::string srcName{ "wallsrc" + TimeSuffix };
    {
      MSource::Wall::Par srcPar;
      srcPar.tW = t;
      srcPar.mom = Momentum0.to_string( Space );
      application.createModule<MSource::Wall>(srcName, srcPar);
    }

    // Zero-momentum propagators
    std::array<std::string, NumQuarks> propNameShort;
    std::array<std::string, NumQuarks> propName;
    for( unsigned int i = 0; i < NumQuarks; ++i )
    {
      propNameShort[i] = "prop" + Sep + Quarks[i].flavour;
      propName[i] = propNameShort[i] + TimeSuffix;
      MFermion::GaugeProp::Par quarkPar;
      quarkPar.solver = solverName[i];
      quarkPar.source = srcName;
      application.createModule<MFermion::GaugeProp>(propName[i], quarkPar);
    }

    // Loop through all momenta
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      const std::string Suffix{ Sep + "p" + Sep + Momenta[p].to_string( Sep ) + TimeSuffix};

      // Make smeared propagators
      std::array<std::string, NumQuarks> propNameSmearedShort;
      std::array<std::string, NumQuarks> propNameSmeared;
      for (unsigned int i = 0; i < NumQuarks; ++i)
      {
        if( !Momenta[p] ) {
          propNameSmearedShort[i] = propNameShort[i];
          propNameSmeared[i] = propName[i];
        } else {
          propNameSmearedShort[i] = "prop" + Sep + Quarks[i].flavour;
          propNameSmeared[i] = propNameSmearedShort[i] + Suffix;
          MSink::Smear::Par smPar;
          smPar.q = propName[i];
          smPar.sink = PropSinkName[p];
          application.createModule<MSink::Smear>(propNameSmeared[i], smPar);
        }
      }
      // contractions
      MContraction::Meson::Par mesPar;
      for (unsigned int i = 0; i < NumQuarks; ++i)
        for (unsigned int j = 0; j < NumQuarks; ++j)
        {
          const std::string MesonSuffix{ Quarks[i].flavour + Sep + Quarks[j].flavour + Suffix };
          mesPar.output = "mesons/C1/" + RunName + '/' + MesonSuffix;
          mesPar.q1     = propNameSmeared[i];
          mesPar.q2     = propNameSmeared[j];
          //mesPar.gammas = "all";
          mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
          mesPar.sink   = ContractionSinkName0;
          application.createModule<MContraction::Meson>(MesonSuffix, mesPar);
        }
    }
  }
}

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
  
  const Grid::Coordinate &lat{GridDefaultLatt()};
  const bool bRun{ lat.size() == 4 && lat[0] == 24 && lat[1] == 24 && lat[2] == 24 && lat[3] == 64 };
  if( !bRun )
  {
    LOG(Warning) << "The parameters in " << RunName << " are designed for --grid 24.24.24.64" << std::endl;
  }

  // global parameters
  static const int Config{ 3000 };
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = Config;
  globalPar.trajCounter.end      = Config + 1;
  globalPar.trajCounter.step     = 40;
  globalPar.runId                = RunName + std::to_string( globalPar.trajCounter.start );
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  globalPar.saveSchedule         = false;

  Application application( globalPar );
  CreateApp( application );
  application.saveParameterFile( globalPar.runId + ".xml" );
  if( bRun )
  {
    application.run();
  }
  else
  {
    for( int i = 1; i < 10; i++ )
    {
      globalPar.trajCounter.start += globalPar.trajCounter.step;
      globalPar.trajCounter.end   += globalPar.trajCounter.step;
      globalPar.runId              = RunName + std::to_string( globalPar.trajCounter.start );
      application.setPar(globalPar);
      application.saveParameterFile( globalPar.runId + ".xml" );
    }
  }

  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
