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

#include "HadronsApp.hpp"

/**************************
 Make application
**************************/

class AppMaker
{
public:
  Application application;
  HModList l;
protected:
  Application Setup( const AppParams &params );
public:
  explicit AppMaker( const AppParams &params )
  : application{ Setup( params ) }, l( application, params ) {}
  void Make();
  void MakeStudy3(const Quark &qf, const Quark &qs);
};

// One-time initialisation
Application AppMaker::Setup( const AppParams &params )
{
  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter  = params.Run.trajCounter;
  globalPar.genetic      = params.Run.genetic;
  globalPar.runId        = params.Run.runId;
  // database options
  globalPar.database.resultDb   = params.Run.dbOptions.resultDb;
  globalPar.database.makeStatDb = params.Run.dbOptions.makeStatDb;
  static const std::string RunID{ params.RunID() };
  globalPar.database.applicationDb = params.Run.dbOptions.applicationDbPrefix + RunID + params.sRunSuffix + ".db";
  Grid::Hadrons::makeFileDir( globalPar.database.applicationDb );
  globalPar.database.restoreModules = false;
  globalPar.database.restoreMemoryProfile = false;
  globalPar.database.restoreSchedule = false;
  globalPar.scheduleFile = params.Run.JobXmlPrefix + RunID + params.sRunSuffix + ".sched";
  globalPar.saveSchedule = true;
  Application application( globalPar );
  // gauge field
  if( params.Run.Gauge.empty() )
  {
    application.createModule<MGauge::Unit>( GaugeFieldName );
  }
  else
  {
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = params.Run.Gauge;
    application.createModule<MIO::LoadNersc>( GaugeFieldName, gaugePar );
  }
  return application;
}

void AppMaker::Make()
{
  // Unitary (sea mass) spectators
  bool bFirstSpec{ true };
  for( const Quark &qSpectator : l.params.SpectatorQuarks )
  {
    for(unsigned int t = l.params.Run.Timeslices.start; t < l.params.Run.Timeslices.end; t += l.params.Run.Timeslices.step)
    {
      for( const Taxonomy &tax : l.params.Taxa )
      {
        for( Common::Momentum p : l.params.Momenta )
        {
          for( int pDoNeg = 0; pDoNeg < ( ( l.params.Run.DoNegativeMomenta && p ) ? 2 : 1 ); ++pDoNeg )
          {
            for( const Quark &qH1 : l.params.HeavyQuarks )
            {
              bool bDidSomething{ false };
              if( l.params.Run.TwoPoint )
              {
                bDidSomething = true;
                l.TakeOwnership( new ModContract2pt( l, tax, qSpectator, qH1, p, t ) );
                l.TakeOwnership( new ModContract2pt( l, tax, qH1, qSpectator, p, t ) );
                if( bFirstSpec )
                {
                  // Additional spectators for 2pt functions only
                  for( const Quark &qSpec2 : l.params.SpectatorQuarks2pt )
                  {
                    l.TakeOwnership( new ModContract2pt( l, tax, qSpec2, qH1, p, t ) );
                    l.TakeOwnership( new ModContract2pt( l, tax, qH1, qSpec2, p, t ) );
                  }
                }
              }
              if( l.params.ThreePoint )
              {
                bDidSomething = true;
                for( const Quark &qH2 : l.params.HeavyQuarks )
                {
                  for( int iHeavy  = l.params.Run.HeavyQuark ? 0 : 1;
                           iHeavy <= l.params.Run.HeavyAnti  ? 1 : 0; ++iHeavy )
                  {
                    const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
                    if( !p || qH1.mass >= qH2.mass )
                    {
                      for( int deltaT : l.params.deltaT )
                      {
                        for( int j = 0; j < NumInsert; j++ )
                        {
                          l.TakeOwnership( new ModContract3pt( l, tax, qH2, qH1, qSpectator, p0, p,
                                                               t, bHeavyAnti, j, deltaT ) );
                          if( p && qH1.mass == qH2.mass )
                          {
                            l.TakeOwnership( new ModContract3pt( l, tax, qH2, qH1, qSpectator, -p, p,
                                                                 t, bHeavyAnti, j, deltaT ) );
                          }
                        }
                      }
                    }
                  }
                }
              }
              if( !bDidSomething )
              {
                // We are only performing residual-mass checks on each propagator
                l.TakeOwnership( new ModProp( l, tax, qH1, p, t ) );
                //l.TakeOwnership( new ModProp( l, tax, qSpectator, p,t) );
              }
            }
            p = -p;
          }
        }
      }
    }
    bFirstSpec = false;
  }
}

void AppMaker::MakeStudy3(const Quark &qf, const Quark &qSpectator)
{
  for(unsigned int t = l.params.Run.Timeslices.start; t < l.params.Run.Timeslices.end; t += l.params.Run.Timeslices.step)
  {
    for( const Taxonomy &tax : l.params.Taxa )
    {
      for( const Quark &qi : l.params.HeavyQuarks )
      {
        for( int iPx = 0; iPx < ( qi.flavour[1] - '0' + 5 ) ; ++iPx )
        {
          const Common::Momentum p( iPx, 0, 0 );
          bool bDidSomething{ false };
          if( l.params.Run.TwoPoint )
          {
            bDidSomething = true;
            l.TakeOwnership( new ModContract2pt( l, tax, qSpectator, qf, p, t ) );
            l.TakeOwnership( new ModContract2pt( l, tax, qf, qSpectator, p, t ) );
            l.TakeOwnership( new ModContract2pt( l, tax, qSpectator, qi, p, t ) );
            l.TakeOwnership( new ModContract2pt( l, tax, qi, qSpectator, p, t ) );
          }
          if( l.params.ThreePoint )
          {
            bDidSomething = true;
            for( int iHeavy  = l.params.Run.HeavyQuark ? 0 : 1;
                iHeavy <= l.params.Run.HeavyAnti  ? 1 : 0; ++iHeavy )
            {
              const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
              for( int deltaT : l.params.deltaT )
              {
                for( int j = 0; j < NumInsert; j++ )
                {
                  l.TakeOwnership( new ModContract3pt( l, tax, qf, qi, qSpectator, -p, p0, t, bHeavyAnti, j, deltaT ) );
                  l.TakeOwnership( new ModContract3pt( l, tax, qi, qf, qSpectator, p0, p , t, bHeavyAnti, j, deltaT ) );
                  l.TakeOwnership( new ModContract3pt( l, tax, qf, qf, qSpectator, -p, p , t, bHeavyAnti, j, deltaT ) );
                  if( iPx == 0 ) // Don't repeat this for non-zero momenta (if job being performed separately)
                    l.TakeOwnership( new ModContract3pt( l, tax, qi, qi, qSpectator, p0, p0, t, bHeavyAnti, j, deltaT ) );
                }
              }
            }
          }
          if( !bDidSomething )
          {
            // We are only performing residual-mass checks on each propagator
            l.TakeOwnership( new ModProp( l, tax, qi, p, t ) );
            //l.TakeOwnership( new ModProp( l, tax, qSpectator, p,t) );
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

  // See whether parameter file exists
  if( argc < 2 || !Common::FileExists( argv[1] ) )
  {
    std::cout << "Usage: xml3pt in_file [job_suffix]\n"
              << "where:\n"
              << "  in_file is the name of the Parameter.xml driving this job\n"
              << "  job_suffix (if present) will be appended to .db, .xml and .sched files"
              << std::endl;
    return EXIT_FAILURE;
  }
  const std::string sXmlFilename{ argv[1] };
  // 2nd parameter is optional. If present (and not a switch) it's a prefix for current run
  const std::string sRunSuffix( argc >= 3 && argv[2][0] != '-' ? ( Common::Period + argv[2] ) : "" );

  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  try
  {
    const AppParams params( sXmlFilename, sRunSuffix );
    AppMaker x( params );
    //x.Make();
    static constexpr int iLight = 0;
    static constexpr int iStrange = 1;
    x.MakeStudy3( params.SpectatorQuarks[iStrange], params.SpectatorQuarks[iLight]   ); //  K <- D
    x.MakeStudy3( params.SpectatorQuarks[iLight],   params.SpectatorQuarks[iLight]   ); // pi <- D
    x.MakeStudy3( params.SpectatorQuarks[iLight],   params.SpectatorQuarks[iStrange] ); //  K <- D_s
    // Run or save the job
    x.application.saveParameterFile( params.Run.JobXmlPrefix + params.RunID() + sRunSuffix + ".xml" );
    if( params.Run.Run )
      x.application.run();
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
