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

static const std::string GaugeFieldNameUnfixed{"gauge"};
static const std::string GaugeFieldName{"gauge_fixed"};
static const std::string GaugeFieldNameUnfixedSmeared{"gauge_stout"};
static const std::string GaugeFieldNameSmeared{"gauge_stout_fixed"};
// Should I stout smear the heavy quark gauge field?
//#define HEAVY_GAUGE_NAME GaugeFieldName
#define HEAVY_GAUGE_NAME GaugeFieldNameSmeared
//#define SMEAR_THEN_FIX  // Otherwise, fix the gauge, then smear

//static const std::string GaugeFieldNameSingle{"gauge_fixed_float"};
//static const std::string GaugeFieldNameSingleSmeared{"gauge_fixed_stout_float"};

// Comment this out to use the actions inherited from Z2 noise and distillation
#define USE_ZMOBIUS
// NB: This is only partially implemented - would need to be finished if required
//#define USE_MIXED_PRECISION

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

#define RESID_LIGHT   1e-8
#define RESID_HEAVY   1e-12
#define MAX_ITERATION 5000

#ifdef USE_ZMOBIUS
#define LS_LIGHT    10
#define LS_STRANGE  24
#define LS_HEAVY    LS_LIGHT
#define RESID_STRANGE  1e-8
#define LIGHT_ITERATION 30000
#else
#define LS_LIGHT    16
#define LS_STRANGE  LS_LIGHT
#define LS_HEAVY    12
#define RESID_STRANGE  1e-12
#define LIGHT_ITERATION 5000
#endif

static const std::string Strange{ "s" };
static const Quark Quarks[] = {
  {"l", 0.005, LS_LIGHT, 1.8, LIGHT_ITERATION, RESID_LIGHT, GaugeFieldName,
    "/tessfs1/work/dp008/dp008/shared/data/eigenpack/C1/vec_fine"}, // light
  {Strange, 0.04, LS_STRANGE, 1.8, MAX_ITERATION, RESID_STRANGE, GaugeFieldName, nullptr}, // strange
  {"h1", 0.58, LS_HEAVY, 1.0, MAX_ITERATION, RESID_HEAVY, HEAVY_GAUGE_NAME, nullptr}, // charm
  // {"h2", 0.64, LS_HEAVY, 1.0, MAX_ITERATION, RESID_HEAVY, HEAVY_GAUGE_NAME, nullptr}, // charm
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
  const unsigned int Nt{ 2 };//64 };
  const unsigned int NtIncrement{ 1 };

  // gauge field
  MIO::LoadNersc::Par gaugePar;
  gaugePar.file = "/tessfs1/work/dp008/dp008/shared/dwf_2+1f/C1/ckpoint_lat";
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
  bool bNeedHeavyGaugeField = false;
  for( int q = 0; q < NumQuarks; q++ )
  {
    if( Common::EqualIgnoreCase( Quarks[q].GaugeField, GaugeFieldNameSmeared ) )
    {
      MGauge::StoutSmearing::Par stoutPar;
      stoutPar.steps = 3;
      stoutPar.rho = 0.1;
#ifdef SMEAR_THEN_FIX
      stoutPar.gauge = GaugeFieldNameUnfixed;
      application.createModule<MGauge::StoutSmearing>(GaugeFieldNameUnfixedSmeared, stoutPar);
      gfPar.gauge = GaugeFieldNameUnfixedSmeared;
      gfPar.Omega_tol = 1e-16;
      gfPar.Phi_tol = 1e-16;
      application.createModule<MGauge::GaugeFix>(GaugeFieldNameSmeared, gfPar);
#else
      stoutPar.gauge = GaugeFieldName;
      application.createModule<MGauge::StoutSmearing>(GaugeFieldNameSmeared, stoutPar);
#endif
      bNeedHeavyGaugeField = true;
      break;
    }
  }

  // Single-precision versions
#ifdef USE_MIXED_PRECISION
  static const std::string suffixSingle{ "_single" };
  static const std::string GaugeFieldNameSingle{ GaugeFieldName + suffixSingle };
  static const std::string GaugeFieldNameSmearedSingle{ GaugeFieldNameSmeared + suffixSingle };
  MUtilities::GaugeSinglePrecisionCast::Par g1Par;
  g1Par.field = GaugeFieldName;
  application.createModule<MUtilities::GaugeSinglePrecisionCast>( GaugeFieldNameSingle, g1Par );
  if( bNeedHeavyGaugeField )
  {
    g1Par.field = GaugeFieldNameSmeared;
    application.createModule<MUtilities::GaugeSinglePrecisionCast>(GaugeFieldNameSmearedSingle, g1Par);
  }
#endif

  // Action and a solver for each quark
  std::array<std::string, NumQuarks> solverName;
  for( int i = 0; i < NumQuarks; ++i )
  {
    const Quark & q{ Quarks[i] };
    // actions
    const std::string ActionName{ "Action" + Sep + q.flavour };
    //const std::string ActionNameSingle{ ActionName + suffixSingle };
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    static const std::string boundary{"1 1 1 -1"};
    static const std::string twist{"0. 0. 0. 0."};
#ifdef USE_ZMOBIUS
    if( Common::EqualIgnoreCase( q.flavour, Strange ) )
    {
      MAction::ScaledDWF::Par actionPar;
      actionPar.gauge = q.GaugeField;
      actionPar.Ls    = q.Ls;
      actionPar.M5    = q.M5;
      actionPar.mass  = q.mass;
      actionPar.boundary = boundary;
      actionPar.twist = twist;
      actionPar.scale = 2.0;
      application.createModule<MAction::ScaledDWF>(ActionName, actionPar);
    }
    else
    {
      static const std::vector<std::complex<double>> omega {
        {1.458064389850479e+00,-0.000000000000000e+00},
        {1.182313183893475e+00,-0.000000000000000e+00},
        {8.309511666859551e-01,-0.000000000000000e+00},
        {5.423524091567911e-01,-0.000000000000000e+00},
        {3.419850204537295e-01,-0.000000000000000e+00},
        {2.113790261902896e-01,-0.000000000000000e+00},
        {1.260742995029118e-01,-0.000000000000000e+00},
        {9.901366519626265e-02,-0.000000000000000e+00},
        {6.863249884465925e-02,5.506585308274019e-02},
        {6.863249884465925e-02,-5.506585308274019e-02},
      };
      MAction::ZMobiusDWF::Par actionPar;
      actionPar.gauge = q.GaugeField;
      actionPar.Ls    = q.Ls;
      actionPar.M5    = q.M5;
      actionPar.mass  = q.mass;
      actionPar.boundary = boundary;
      actionPar.b = 1.;
      actionPar.c = 0.;
      actionPar.omega = omega;
      actionPar.twist = twist;
      application.createModule<MAction::ZMobiusDWF>(ActionName, actionPar);
    }
#else
    MAction::DWF::Par actionPar;
    actionPar.gauge = q.GaugeField;
    actionPar.Ls    = q.Ls;
    actionPar.M5    = q.M5;
    actionPar.mass  = q.mass;
    actionPar.boundary = boundary;
    actionPar.twist = twist;
    application.createModule<MAction::DWF>(ActionName, actionPar);
#endif

    // eigenpacks for deflation
    std::string epackObjname;
    if( q.EigenPackFilename && q.Ls == 16 ) // The eigenpacks expect Ls == 16
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
#ifdef USE_ZMOBIUS
    if( !Common::EqualIgnoreCase( q.flavour, Strange ) )
    {
      MSolver::ZRBPrecCG::Par solverPar;
      solverPar.eigenPack    = epackObjname;
      solverPar.action       = ActionName;
      solverPar.residual     = q.residual;
      solverPar.maxIteration = q.maxIteration;
      application.createModule<MSolver::ZRBPrecCG>(solverName[i], solverPar);
    }
    else
#endif
    {
      MSolver::RBPrecCG::Par solverPar;
      solverPar.eigenPack    = epackObjname;
      solverPar.action       = ActionName;
      solverPar.residual     = q.residual;
      solverPar.maxIteration = q.maxIteration;
      application.createModule<MSolver::RBPrecCG>(solverName[i], solverPar);
    }
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
  static std::array<std::string,NumMomenta> PropSinkNameNeg;
  {
    MSink::Point::Par sinkPar;
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      if( Momenta[p] )
      {
        PropSinkName[p] = "prop_sink_p_" + Momenta[p].to_string( Sep );
        sinkPar.mom = Momenta[p].to_string( Space );
        application.createModule<MSink::Point>(PropSinkName[p], sinkPar);
        PropSinkNameNeg[p] = "prop_sink_p_" + Momenta[p].to_string( Sep, true );
        sinkPar.mom = Momenta[p].to_string( Space, true );
        application.createModule<MSink::Point>(PropSinkNameNeg[p], sinkPar);
      }
    }
  }
  
  // Loop through all timeslices
  for (unsigned int t = NtIncrement; t < Nt; t += NtIncrement )
  {
    const std::string TimeSuffix{ Sep + "t" + Sep + std::to_string( t ) };

    // Zero-momentum wall-source with this momenta
    const std::string srcName{ "wallsrc" + TimeSuffix };
    {
      MSource::Wall::Par srcPar;
      srcPar.tW = t;
      srcPar.mom = Momentum0.to_string4d( Space );
      application.createModule<MSource::Wall>(srcName, srcPar);
    }

    // Zero-momentum propagators
    std::array<std::string, NumQuarks> propNameShort;
    std::array<std::string, NumQuarks> propName;
    for( unsigned int i = 0; i < NumQuarks; ++i )
    {
      propNameShort[i] = "prop" + Sep + Quarks[i].flavour;
      propName[i] = propNameShort[i] + TimeSuffix;
#ifdef USE_ZMOBIUS
      if( !Common::EqualIgnoreCase( Quarks[i].flavour, Strange ) )
      {
        MFermion::ZGaugeProp::Par quarkPar;
        quarkPar.solver = solverName[i];
        quarkPar.source = srcName;
        application.createModule<MFermion::ZGaugeProp>(propName[i], quarkPar);
      }
      else
#endif
      {
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = solverName[i];
        quarkPar.source = srcName;
        application.createModule<MFermion::GaugeProp>(propName[i], quarkPar);
      }
    }

    // Loop through all momenta
    for( unsigned int p = 0; p < NumMomenta; ++p )
    {
      const std::string Suffix   { Sep + "p" + Sep + Momenta[p].to_string( Sep ) + TimeSuffix};
      const std::string SuffixNeg{ Sep + "p" + Sep + Momenta[p].to_string( Sep, true ) + TimeSuffix};

      // Make smeared propagators
      std::array<std::string, NumQuarks> propNameSmearedShort;
      std::array<std::string, NumQuarks> propNameSmeared;
      std::array<std::string, NumQuarks> propNameSmearedNeg;
      for (unsigned int i = 0; i < NumQuarks; ++i)
      {
        if( !Momenta[p] ) {
          propNameSmearedShort[i] = propNameShort[i];
          propNameSmeared[i] = propName[i];
          propNameSmearedNeg[i] = propName[i];
        } else {
          propNameSmearedShort[i] = "prop" + Sep + Quarks[i].flavour;
          propNameSmeared[i] = propNameSmearedShort[i] + Suffix;
          MSink::Smear::Par smPar;
          smPar.q = propName[i];
          smPar.sink = PropSinkName[p];
          application.createModule<MSink::Smear>(propNameSmeared[i], smPar);
          propNameSmearedNeg[i] = propNameSmearedShort[i] + SuffixNeg;
          smPar.sink = PropSinkNameNeg[p];
          application.createModule<MSink::Smear>(propNameSmearedNeg[i], smPar);
        }
      }
      // contractions
      MContraction::Meson::Par mesPar;
      for (unsigned int i = 0; i < NumQuarks; ++i)
        for (unsigned int j = 0; j < NumQuarks; ++j)
        {
          static const std::string MyGammas{ "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)" }; // used to be "all"
          static const std::string MesonDir{ "mesons/C1/" + RunName
            //+ "/Unsm"
          //#ifdef USE_ZMOBIUS
                                           //+ "/ZMob"
          //#endif
                                           + "/" };
          std::string MesonSuffix{ Quarks[i].flavour + Sep + Quarks[j].flavour + Suffix };
          mesPar.output = MesonDir + MesonSuffix;
          mesPar.q1     = propNameSmeared[i];
          mesPar.q2     = propNameSmearedNeg[j];
          mesPar.gammas = MyGammas;
          mesPar.sink   = ContractionSinkName0;
          application.createModule<MContraction::Meson>(MesonSuffix, mesPar);
          if( Momenta[p] )
          {
            MesonSuffix = Quarks[i].flavour + Sep + Quarks[j].flavour + SuffixNeg;
            mesPar.output = MesonDir + MesonSuffix;
            mesPar.q1     = propNameSmearedNeg[i];
            mesPar.q2     = propNameSmeared[j];
            application.createModule<MContraction::Meson>(MesonSuffix, mesPar);
          }
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
  static const int NumConfigs{ 1 };
  assert( NumConfigs > 0 );
  Application::GlobalPar globalPar;
  globalPar.trajCounter.start    = Config;
  globalPar.trajCounter.step     = 40;
  globalPar.trajCounter.end      = Config + (NumConfigs - 1) * globalPar.trajCounter.step + 1;
  globalPar.runId                = RunName;
  globalPar.genetic.maxGen       = 1000;
  globalPar.genetic.maxCstGen    = 200;
  globalPar.genetic.popSize      = 20;
  globalPar.genetic.mutationRate = .1;
  globalPar.saveSchedule         = false;

  Application application( globalPar );
  CreateApp( application );
  application.saveParameterFile( RunName + std::to_string( globalPar.trajCounter.start ) + ".xml" );
  if( bRun )
  {
    application.run();
  }
  else if(( 0 ))
  {
    for( int i = 1; i < 10; i++ )
    {
      globalPar.trajCounter.start += globalPar.trajCounter.step;
      globalPar.trajCounter.end   += globalPar.trajCounter.step;
      application.setPar( globalPar );
      application.saveParameterFile( RunName + std::to_string( globalPar.trajCounter.start ) + ".xml" );
    }
  }

  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return EXIT_SUCCESS;
}
