/*************************************************************************************
 
 Create XML for a 2-pt mass plateau study using Z2 wall sources.
 Source file: xmlGFWall.cpp
 Copyright (C) 2019-2020
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

static const std::string RunName{ "GFWall" };

using namespace Grid;
using namespace Hadrons;

static const std::string Sep{ "_" };    // used inside filenames
static const std::string SepBig{ "." };    // used inside filenames
static const std::string Space{ " " };  // whitespace as field separator / human readable info

static const std::string GaugeFieldNameUnfixed{"gauge"};
static const std::string GaugeFieldName{"gauge_fixed"};
static const std::string GaugeFieldNameUnfixedSmeared{"gauge_stout"};
static const std::string GaugeFieldNameSmeared{"gauge_stout_fixed"};

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

// Values used by Fionn on 48^3
//#define LS_LIGHT    10
//#define LS_STRANGE  24
//#define LS_HEAVY    LS_LIGHT
//#define RESID_STRANGE  1e-8
//#define LIGHT_ITERATION 30000
//#define LIGHT_EIGENPACK nullptr

// Values for 24^3
#define LS_LIGHT    16
#define LS_STRANGE  LS_LIGHT
#define LS_HEAVY    12
#define RESID_STRANGE  1e-12
#define LIGHT_ITERATION 5000
#define LIGHT_EIGENPACK C1EigenPack

#define RESID_LIGHT   1e-8
#define RESID_HEAVY   1e-12
#define MAX_ITERATION 5000

// This is the correct configuration - C1
static const char * GaugePrefix{ "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat" };
static const char * C1EigenPack{ "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine" };
static const std::string OutputPrefix{ "" };

static const std::string Strange{ "s" };
static const Quark Quarks[] = {
  {"l", 0.005, LS_LIGHT, 1.8, LIGHT_ITERATION, RESID_LIGHT, GaugeFieldName, LIGHT_EIGENPACK }, // light
  {Strange, 0.04, LS_STRANGE, 1.8, MAX_ITERATION, RESID_STRANGE, GaugeFieldName, nullptr}, // strange
  {"h1", 0.58, LS_HEAVY, 1.0, MAX_ITERATION, RESID_HEAVY, GaugeFieldNameSmeared, nullptr}, // charm
  //{"h2", 0.64, LS_HEAVY, 1.0, MAX_ITERATION, RESID_HEAVY, GaugeFieldNameSmeared, nullptr}, // charm
};
static constexpr int NumQuarks{ sizeof( Quarks ) / sizeof( Quarks[0] ) };

static const Common::Momentum Momentum0{ 0, 0, 0 };

static const Common::Momentum Momenta[] = {
  { 0, 0, 0 },
  { 1, 0, 0 },
  { 1, 1, 0 },
  { 1, 1, 1 },
  { 2, 0, 0 },/*
  { 0, 1, 0 },
  { 0, 0, 1 },
  { 1, 0, 1 },
  { 0, 1, 1 },
  { 0, 2, 0 },
  { 0, 0, 2 },*/
};
static constexpr int NumMomenta{ sizeof( Momenta ) / sizeof( Momenta[0] ) };

void CreateApp( Application &application )
{
  const unsigned int NtStart{ 0 };
  const unsigned int Nt{ 32 };//64 };
  const unsigned int NtIncrement{ 4 };

  // gauge field
  MIO::LoadNersc::Par gaugePar;
  gaugePar.file = GaugePrefix;
  application.createModule<MIO::LoadNersc>(GaugeFieldNameUnfixed, gaugePar);

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

  // Do I need the stout smeared heavy quark gauge field?
  for( int q = 0; q < NumQuarks; q++ )
  {
    if( Common::EqualIgnoreCase( Quarks[q].GaugeField, GaugeFieldNameSmeared ) )
    {
      MGauge::StoutSmearing::Par stoutPar;
      stoutPar.steps = 3;
      stoutPar.rho = 0.1;
      stoutPar.gauge = GaugeFieldName;
      application.createModule<MGauge::StoutSmearing>(GaugeFieldNameSmeared, stoutPar);
      break;
    }
  }

  // Action and a solver for each quark
  std::array<std::string, NumQuarks> solverName;
  std::array<std::string, NumQuarks> propName;
  for( int i = 0; i < NumQuarks; ++i )
  {
    const Quark & q{ Quarks[i] };
    propName[i] = "prop" + Sep + q.flavour;
    // actions
    const std::string ActionName{ "Action" + Sep + q.flavour };
    //const std::string ActionNameSingle{ ActionName + suffixSingle };
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    static const std::string boundary{"1 1 1 -1"};
    static const std::string twist{"0. 0. 0. 0."};
    MAction::ScaledDWF::Par actionPar;
    actionPar.gauge = q.GaugeField;
    actionPar.Ls    = q.Ls;
    actionPar.M5    = q.M5;
    actionPar.mass  = q.mass;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    actionPar.scale = 1.0;
    application.createModule<MAction::ScaledDWF>(ActionName, actionPar);

    // eigenpacks for deflation
    std::string epackObjname;
    if( q.EigenPackFilename )
    {
      epackObjname = "epack_" + q.flavour;
      MIO::LoadFermionEigenPack::Par epPar;
      epPar.filestem = q.EigenPackFilename;
      epPar.multiFile = false;
      epPar.size = 600;
      epPar.Ls = q.Ls;
      application.createModule<MIO::LoadFermionEigenPack>(epackObjname, epPar);
    }
    // solvers
    solverName[i] = "CG" + Sep + q.flavour;
    MSolver::RBPrecCG::Par solverPar;
    solverPar.eigenPack    = epackObjname;
    solverPar.action       = ActionName;
    solverPar.residual     = q.residual;
    solverPar.maxIteration = q.maxIteration;
    application.createModule<MSolver::RBPrecCG>(solverName[i], solverPar);
  }

  // Contraction sink for all momentum
  static const std::string ContractionSinkPrefix{ "sink_contract_p_" };
  static const std::string ContractionSinkName0{ ContractionSinkPrefix + Momentum0.to_string( Sep ) };
  std::array<std::string,NumMomenta> ContractionSinkName;
  {
    MSink::ScalarPoint::Par sinkPar;
    sinkPar.mom = Momentum0.to_string( Space );
    application.createModule<MSink::ScalarPoint>(ContractionSinkName0, sinkPar);
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      if( !Momenta[p] )
        ContractionSinkName[p] = ContractionSinkName0;
      else
      {
        ContractionSinkName[p] = ContractionSinkPrefix + Momenta[p].to_string( Sep );
        sinkPar.mom = Momenta[p].to_string( Space );
        application.createModule<MSink::ScalarPoint>(ContractionSinkName[p], sinkPar);
      }
    }
  }
  
  // Smearing sink for all momentum
  static const std::string SmearSinkPrefix{ "sink_smear_p_" };
  static const std::string SmearSinkName0{ SmearSinkPrefix + Momentum0.to_string( Sep ) };
  std::array<std::string,NumMomenta> SmearSinkName;
  {
    MSink::Point::Par sinkPar;
    sinkPar.mom = Momentum0.to_string( Space );
    application.createModule<MSink::Point>(SmearSinkName0, sinkPar);
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      if( !Momenta[p] )
        SmearSinkName[p] = SmearSinkName0;
      else
      {
        SmearSinkName[p] = SmearSinkPrefix + Momenta[p].to_string( Sep );
        sinkPar.mom = Momenta[p].to_string( Space );
        application.createModule<MSink::Point>(SmearSinkName[p], sinkPar);
      }
    }
  }
  
  // Loop through all timeslices
  for (unsigned int t = NtStart; t < Nt; t += NtIncrement )
  {
    const std::string TimeSuffix{ Sep + "t" + Sep + std::to_string( t ) };
    const std::string Suffix0{ Sep + "p" + Sep + Momentum0.to_string( Sep ) + TimeSuffix};

    // Make zero-momentum wall-source
    const std::string srcName0{ "wallsrc" + Suffix0 };
    {
      MSource::Wall::Par srcPar;
      srcPar.tW = t;
      srcPar.mom = Momentum0.to_string4d( Space );
      application.createModule<MSource::Wall>(srcName0, srcPar);
    }

    // Make zero-momentum propagators
    std::array<std::string, NumQuarks> propNameUnsmeared0;
    std::array<std::string, NumQuarks> propNameSmeared0;
    for (unsigned int i = 0; i < NumQuarks; ++i)
    {
      // This version used for wall-point
      propNameUnsmeared0[i] = propName[i] + Sep + srcName0;
      MFermion::GaugeProp::Par quarkPar;
      quarkPar.solver = solverName[i];
      quarkPar.source = srcName0;
      application.createModule<MFermion::GaugeProp>(propNameUnsmeared0[i], quarkPar);
      // And smeared version for wall-wall
      propNameSmeared0[i] = propName[i] + Sep + "sm" + Sep + srcName0;
      MSink::Smear::Par smearPar;
      smearPar.q = propNameUnsmeared0[i];
      smearPar.sink = SmearSinkName0;
      application.createModule<MSink::Smear>(propNameSmeared0[i], smearPar);
    }

    // Loop through all momenta
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      const std::string Suffix{ Sep + "p" + Sep + Momenta[p].to_string( Sep ) + TimeSuffix};

      // Wall-source(s) with this momentum
      const std::string srcName{ Momenta[p] ? "wallsrc" + Suffix : srcName0 };
      if( Momenta[p] )
      {
        MSource::Wall::Par srcPar;
        srcPar.tW = t;
        srcPar.mom = Momenta[p].to_string4d( Space );
        application.createModule<MSource::Wall>(srcName, srcPar);
      }

      // Make propagators
      std::array<std::string, NumQuarks> propNameUnsmeared;
      std::array<std::string, NumQuarks> propNameSmeared;
      for (unsigned int i = 0; i < NumQuarks; ++i)
      {
        if( !Momenta[p] )
        {
          propNameUnsmeared[i] = propNameUnsmeared0[i];
          propNameSmeared[i]   = propNameSmeared0[i];
        }
        else
        {
          // This version used for wall-point
          propNameUnsmeared[i] = propName[i] + Sep + srcName;
          MFermion::GaugeProp::Par quarkPar;
          quarkPar.solver = solverName[i];
          quarkPar.source = srcName;
          application.createModule<MFermion::GaugeProp>(propNameUnsmeared[i], quarkPar);
          // And smeared version for wall-wall
          propNameSmeared[i] = propName[i] + Sep + "sm" + Sep + srcName;
          MSink::Smear::Par smearPar;
          smearPar.q = propNameUnsmeared[i];
          smearPar.sink = SmearSinkName0; // we don't need momentum here
          application.createModule<MSink::Smear>(propNameSmeared[i], smearPar);
        }
      }

      // contractions
      MContraction::Meson::Par mesPar;
      mesPar.gammas = "all";
      for (unsigned int i = 0; i < NumQuarks; ++i)
        for (unsigned int j = 0; j < NumQuarks; ++j)
        {
          static const std::string MesonDir{ OutputPrefix + "mesons/C1/" }; //+ RunName + "/" };
          // Wall-point
          std::string MesonSuffix{ Quarks[i].flavour + Sep + Quarks[j].flavour + Suffix };
          std::string TypeMesonSuffix{ "GFWP/" + MesonSuffix };
          mesPar.output = MesonDir + TypeMesonSuffix;
          mesPar.q1     = propNameUnsmeared0[i];
          mesPar.q2     = propNameUnsmeared[j];
          mesPar.sink   = ContractionSinkName[p];
          application.createModule<MContraction::Meson>(TypeMesonSuffix, mesPar);
          // Wall-wall
          TypeMesonSuffix = "GFWW/" + MesonSuffix;
          mesPar.output = MesonDir + TypeMesonSuffix;
          mesPar.q1     = propNameSmeared[i];
          mesPar.q2     = propNameSmeared[j];
          mesPar.sink   = ContractionSinkName0; // NB: sink not actually used (with sliced propagators)
          application.createModule<MContraction::Meson>(TypeMesonSuffix, mesPar);
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
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = 3000;
  globalPar.trajCounter.step     = 40;
  globalPar.trajCounter.end      = 3001;
  globalPar.runId                = RunName;
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  globalPar.saveSchedule         = false;

  Application application( globalPar );
  CreateApp( application );
  application.saveParameterFile( RunName + ".template.xml" );
  if( bRun )
    application.run();

  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
