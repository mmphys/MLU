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

#include "xml3pt.hpp"

/**************************
 Base class for my application maker
**************************/

void AppMaker::WriteTaxonomy( std::ostringstream &s, const std::vector<Taxonomy> &Taxa ) const
{
  // Start with all the Families and sink types
  std::vector<Taxonomy> taxSorted{ Taxa };
  std::sort( taxSorted.begin(), taxSorted.end() );
  for( std::size_t i = 0; i < taxSorted.size(); ++i )
  {
    if( i == 0 || taxSorted[i].family != taxSorted[i - 1].family )
    {
      if( i )
        s << Sep;
      s << taxSorted[i].family;
    }
    s << SpeciesNameShort( taxSorted[i].species );
  }
}

void AppMaker::ReadQuarks( XmlReader &r )
{
  std::string sQuarks;
  read( r, "Quarks", sQuarks );
  std::vector<std::string> Quarks = Common::ArrayFromString( sQuarks, Common::WhiteSpace );
  for( const std::string &sQuark : Quarks )
  {
    LOG(Message) << "Reading quarks from " << sQuark << std::endl;
    std::vector<Quark> vq;
    r.readDefault( sQuark, vq );
    for( Quark &q : vq )
    {
      LOG(Message) << " Loaded quark " << sQuark << " : " << q.flavour << std::endl;
      Q.emplace( q.flavour, q );
    }
  }
}

std::vector<Quark *> AppMaker::ValidateQuarkList( const std::string &sList, const std::string &sName, std::size_t Min )
{
  LOG(Message) << sName << " quarks: " << sList << std::endl;
  std::vector<std::string> vs = Common::ArrayFromString( sList, Common::WhiteSpace );
  Common::NoDuplicates( vs, sName, Min );
  std::vector<Quark *> v;
  v.reserve( vs.size() );
  for( const std::string &s : vs )
    v.emplace_back( &Q.at( s ) );
  return v;
}

// One-time initialisation
void AppMaker::Setup( XmlReader &r, const std::string &sRunSuffix )
{
  // Read my parameters from Xml
  read( r, "JobFilePrefix", JobFilePrefix );
  read( r, "Run", Run );
  ReadQuarks( r );

  // Allow the derived class to initialise
  Setup( r );

  std::string RunIDSuffix{ RunID() };
  RunIDSuffix.append( sRunSuffix );
  JobFilePrefix.append( RunIDSuffix );

  // global parameters
  Application::GlobalPar globalPar;
  globalPar.trajCounter  = appPar.Run.trajCounter;
  globalPar.genetic      = appPar.Run.genetic;
  globalPar.runId        = appPar.Run.runId;
  // database options
  globalPar.database.resultDb   = appPar.Run.dbOptions.resultDb;
  globalPar.database.makeStatDb = appPar.Run.dbOptions.makeStatDb;
  globalPar.database.applicationDb = appPar.Run.dbOptions.applicationDbPrefix + RunIDSuffix + ".db";
  Grid::Hadrons::makeFileDir( globalPar.database.applicationDb );
  globalPar.database.restoreModules = false;
  globalPar.database.restoreMemoryProfile = false;
  globalPar.database.restoreSchedule = false;
  globalPar.scheduleFile = JobFilePrefix + ".sched";
  globalPar.saveSchedule = true;
  Application application( globalPar );
  // gauge field
  if( appPar.Run.Gauge.empty() )
  {
    application.createModule<MGauge::Unit>( GaugeFieldName );
  }
  else
  {
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = appPar.Run.Gauge;
    application.createModule<MIO::LoadNersc>( GaugeFieldName, gaugePar );
  }
}

#define APP_MAKER_GLUE( ClassName ) \
{ \
  if( reader->push( ClassName::sXmlTagName ) ) \
  { \
    reader->pop(); \
    LOG(Message) << "Making " << ClassName::sXmlTagName << std::endl; \
    x.reset( new ClassName( appParams ) ); \
  } \
}

int AppMaker::MakeThreePoint( int argc, char *argv[], const std::string &sXmlFilename, const std::string &sRunSuffix )
{
  // initialization //////////////////////////////////////////////////////////
  Grid_init(&argc, &argv);
  HadronsLogError.Active(GridLogError.isActive());
  HadronsLogWarning.Active(GridLogWarning.isActive());
  HadronsLogMessage.Active(GridLogMessage.isActive());
  HadronsLogIterative.Active(GridLogIterative.isActive());
  HadronsLogDebug.Active(GridLogDebug.isActive());
  int iReturn{ EXIT_FAILURE };
  try
  {
    static const std::string sXmlTopLevel{ "xml3pt" };
    std::unique_ptr<XmlReader> reader( new XmlReader( sXmlFilename, false, sXmlTopLevel ) );
    const AppParams appParams( *reader );
    std::unique_ptr<AppMaker> x;
    for( int i = 0; i < 2 && !x.get(); ++i )
    {
      switch( i )
      {
        case 0:
          APP_MAKER_GLUE( Study2 )
          break;
        case 1:
          APP_MAKER_GLUE( Study3 )
          break;
      }
    }
    if( !x.get() )
      throw std::runtime_error( "xml format unrecognised" );
    x->Setup( *reader, sRunSuffix );
    reader.reset( nullptr );
    x->Make();
    // Run or save the job
    x->application.saveParameterFile( x->JobFilePrefix + ".xml" );
    if( x->Run )
      x->application.run();
    iReturn = EXIT_SUCCESS;
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
  }

  // epilogue
  LOG(Message) << "Grid is finalizing now" << std::endl;
  Grid_finalize();
  
  return iReturn;
}

/**************************
 Study 2: generic study for heavy-light meson decays
**************************/

const std::string Study2::sXmlTagName{ "Study2" };

// One-time initialisation
void Study2::Setup( XmlReader &r )
{
  // Read parameters from Xml
  read( r, sXmlTagName, makePar );
  // Check the taxonomy
  Taxa = Common::ArrayFromString<Taxonomy>( makePar.Taxa );
  Common::NoDuplicates( Taxa, "Taxa", 1 );
  for( const Taxonomy &t : Taxa )
    if( t == Family::GR )
      LOG(Warning) << t << " don't swap the gamma structures at sink and source (18-Feb-2021)" << std::endl;
  // Check parameters make sense
  //if( !( Run.TwoPoint || Run.HeavyQuark ||Run.HeavyAnti ) )
    //throw std::runtime_error( "At least one must be true of: TwoPoint; HeavyQuark; or HeavyAnti" );
  Momenta = Common::ArrayFromString<Common::Momentum>( makePar.Momenta );
  Common::NoDuplicates( Momenta, "Momenta", 1 );
  ThreePoint = makePar.HeavyQuark || makePar.HeavyAnti;
  if( ThreePoint )
  {
    deltaT = Common::ArrayFromString<int>( makePar.deltaT );
    Common::NoDuplicates( deltaT, "deltaT", 1 );
    gamma = Common::ArrayFromString<Gamma::Algebra>( makePar.gamma );
    Common::NoDuplicates( gamma, "gamma", 1 );
  }
  HeavyQuarks        = ValidateQuarkList( makePar.Heavy, "Heavy", 1 );
  SpectatorQuarks    = ValidateQuarkList( makePar.Spectator, "Spectator", 1 );
  SpectatorQuarks2pt = ValidateQuarkList( makePar.SpectatorExtra2pt, "SpectatorExtra2pt", 0 );
}

// Make a unique RunID string that completely describes the run
std::string Study2::RunID() const
{
  std::ostringstream s;
  WriteTaxonomy( s, Taxa );
  for( const Quark *q : SpectatorQuarks )
    s << Sep << q->flavour;
  if( makePar.TwoPoint )
    s << Sep << "2pt";
  if( makePar.HeavyQuark )
    s << Sep << "quark";
  if( makePar.HeavyAnti )
    s << Sep << "anti";
  for( const Quark *q : HeavyQuarks )
    s << Sep << q->flavour;
  s << Sep << "t" << Sep << makePar.Timeslices.start << Sep << makePar.Timeslices.end << Sep << makePar.Timeslices.step;
  if( !makePar.DoNegativeMomenta )
    s << Sep << "pos";
  s << Sep << "p";
  for( const Common::Momentum &p : Momenta )
    s << Sep << p.to_string( Sep );
  if( ThreePoint )
  {
    s << Sep << "dT";
    for( int dT : deltaT )
      s << Sep << std::to_string( dT );
  }
  return s.str();
}

void Study2::Make()
{
  // Unitary (sea mass) spectators
  bool bFirstSpec{ true };
  for( const Quark *qSpectator : SpectatorQuarks )
  {
    for(unsigned int t = makePar.Timeslices.start; t < makePar.Timeslices.end; t += makePar.Timeslices.step)
    {
      for( const Taxonomy &tax : Taxa )
      {
        for( Common::Momentum p : Momenta )
        {
          for( int pDoNeg = 0; pDoNeg < ( ( makePar.DoNegativeMomenta && p ) ? 2 : 1 ); ++pDoNeg )
          {
            for( const Quark *qH1 : HeavyQuarks )
            {
              bool bDidSomething{ false };
              if( makePar.TwoPoint )
              {
                bDidSomething = true;
                l.TakeOwnership( new ModContract2pt( l, tax, *qSpectator, *qH1, p, t ) );
                l.TakeOwnership( new ModContract2pt( l, tax, *qH1, *qSpectator, p, t ) );
                if( bFirstSpec )
                {
                  // Additional spectators for 2pt functions only
                  for( const Quark *qSpec2 : SpectatorQuarks2pt )
                  {
                    l.TakeOwnership( new ModContract2pt( l, tax, *qSpec2, *qH1, p, t ) );
                    l.TakeOwnership( new ModContract2pt( l, tax, *qH1, *qSpec2, p, t ) );
                  }
                }
              }
              if( ThreePoint )
              {
                bDidSomething = true;
                for( const Quark *qH2 : HeavyQuarks )
                {
                  for( int iHeavy  = makePar.HeavyQuark ? 0 : 1;
                           iHeavy <= makePar.HeavyAnti  ? 1 : 0; ++iHeavy )
                  {
                    const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
                    if( !p || qH1->mass >= qH2->mass )
                    {
                      for( int deltaT : deltaT )
                      {
                        for( Gamma::Algebra j : gamma )
                        {
                          l.TakeOwnership( new ModContract3pt( l, tax, false, *qH2, *qH1, *qSpectator, p0, p,
                                                               j, deltaT, t, bHeavyAnti ) );
                          if( p && qH1->mass == qH2->mass )
                          {
                            l.TakeOwnership( new ModContract3pt( l, tax, false, *qH2, *qH1, *qSpectator, -p, p,
                                                                 j, deltaT, t, bHeavyAnti ) );
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
                l.TakeOwnership( new ModProp( l, tax, *qH1, p, t ) );
                //l.TakeOwnership( new ModProp( l, tax, *qSpectator, p,t) );
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

/**************************
 Study 3: compute R2 and R1 ratios for meson decays
**************************/

const std::string Study3::sXmlTagName{ "Study3" };

// One-time initialisation
void Study3::Setup( XmlReader &r )
{
  // Read parameters from Xml
  read( r, sXmlTagName, makePar );
  // Check the taxonomy
  Taxa = Common::ArrayFromString<Taxonomy>( makePar.Taxa );
  Common::NoDuplicates( Taxa, "Taxa", 1 );
  for( const Taxonomy &t : Taxa )
    if( t == Family::GR )
      LOG(Warning) << t << " don't swap the gamma structures at sink and source (18-Feb-2021)" << std::endl;
  // Check parameters make sense
  //if( !( Run.TwoPoint || Run.HeavyQuark ||Run.HeavyAnti ) )
    //throw std::runtime_error( "At least one must be true of: TwoPoint; HeavyQuark; or HeavyAnti" );
  ThreePoint = makePar.HeavyQuark || makePar.HeavyAnti;
  if( ThreePoint )
  {
    deltaT = Common::ArrayFromString<int>( makePar.deltaT );
    Common::NoDuplicates( deltaT, "deltaT", 1 );
    gamma = Common::ArrayFromString<Gamma::Algebra>( makePar.gamma );
    Common::NoDuplicates( gamma, "gamma", 1 );
  }
  HeavyQuarks = ValidateQuarkList( makePar.Heavy, "Heavy", 1 );
  CountHeavyMomenta = 0;
  for( const Decay &d : makePar.Decays )
  {
    LOG(Message) << "Decay " << d.name << ": light=" << d.qLight << ", spectator=" << d.qSpectator
                 << ", " << d.HeavyMom.size() << " heavy quarks" << std::endl;
    Q.at( d.qLight );
    Q.at( d.qSpectator );
    if( d.HeavyMom.empty() )
      throw std::runtime_error( "There should be at least 1 heavy quark" );
    for( const HeavyMomenta &hp : d.HeavyMom )
    {
      LOG(Message) << " heavy " << hp.qHeavy << ": momenta " << hp.Momenta << std::endl;
      Q.at( hp.qHeavy );
      std::vector<Common::Momentum> Momenta = Common::ArrayFromString<Common::Momentum>( hp.Momenta );
      Common::NoDuplicates( Momenta, "HeavyMomenta", 1 );
      CountHeavyMomenta += Momenta.size();
    }
  }
}

// Make a unique RunID string that completely describes the run
std::string Study3::RunID() const
{
  std::ostringstream s;
  WriteTaxonomy( s, Taxa );
  for( const Decay &d : makePar.Decays )
    s << Sep << d.name;
  if( makePar.TwoPoint )
    s << Sep << "2pt";
  if( makePar.HeavyQuark )
    s << Sep << "quark";
  if( makePar.HeavyAnti )
    s << Sep << "anti";
  for( const Quark *q : HeavyQuarks )
    s << Sep << q->flavour;
  s << Sep << "t" << Sep << makePar.Timeslices.start << Sep << makePar.Timeslices.end << Sep << makePar.Timeslices.step;
  if( makePar.DoNegativeMomenta )
    s << Sep << "neg";
  s << Sep << "hp" << Sep << CountHeavyMomenta;
  if( ThreePoint )
  {
    s << Sep << "dT";
    for( int dT : deltaT )
      s << Sep << std::to_string( dT );
  }
  return s.str();
}

void Study3::MakeStudy3( const Decay &d )
{
  const Quark &ql{ Q.at( d.qLight ) };
  const Quark &qSpectator{ Q.at( d.qSpectator ) };
  for(unsigned int t = makePar.Timeslices.start; t < makePar.Timeslices.end; t += makePar.Timeslices.step)
  {
    for( const Taxonomy &tax : Taxa )
    {
      for( const HeavyMomenta &hp : d.HeavyMom )
      {
        const Quark &qh{ Q.at( hp.qHeavy ) };
        std::vector<Common::Momentum> Momenta = Common::ArrayFromString<Common::Momentum>( hp.Momenta );
        for( Common::Momentum p : Momenta )
        {
          for( int pDoNeg = 0; pDoNeg < ( ( makePar.DoNegativeMomenta && p ) ? 2 : 1 ); ++pDoNeg )
          {
            bool bDidSomething{ false };
            if( makePar.TwoPoint )
            {
              // Create the two-point functions I need for R1. Second pair not really needed, but doesn't take much time
              bDidSomething = true;
              l.TakeOwnership( new ModContract2pt( l, tax, qSpectator, ql, p, t ) );
              l.TakeOwnership( new ModContract2pt( l, tax, ql, qSpectator, p, t ) );
              // These two only needed at zero momentum ... but doesn't take much time
              l.TakeOwnership( new ModContract2pt( l, tax, qSpectator, qh, p, t ) );
              l.TakeOwnership( new ModContract2pt( l, tax, qh, qSpectator, p, t ) );
            }
            if( ThreePoint )
            {
              // Create the three-point functions I need for R2 (which includes what's required for R1)
              bDidSomething = true;
              for( int iHeavy  = makePar.HeavyQuark ? 0 : 1;
                       iHeavy <= makePar.HeavyAnti  ? 1 : 0; ++iHeavy )
              {
                const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
                for( int deltaT : deltaT )
                {
                  for( Gamma::Algebra j : gamma )
                  {
                    const bool bHLB{ makePar.HeavyLightBackwards };
                    const unsigned int tHLB{ bHLB ? l.params.TimeBound( t - deltaT ) : t };
                    const Common::Momentum pHLB{ bHLB ? p : -p };
                    l.TakeOwnership(new ModContract3pt(l,tax, bHLB,ql,qh, qSpectator,pHLB,p0,j, deltaT, tHLB, bHeavyAnti));
                    l.TakeOwnership(new ModContract3pt(l,tax,false,qh,ql, qSpectator, p0, p, j, deltaT, t, bHeavyAnti ) );
                    if( makePar.R2Ratio ) // In practice, M2Ratios should equal !HeavyLightBackwards
                    {
                      l.TakeOwnership(new ModContract3pt(l,tax,false,ql,ql, qSpectator,-p, p, j, deltaT, t, bHeavyAnti ) );
                      if( !p ) // Don't repeat this for non-zero momenta (if job being performed separately)
                        l.TakeOwnership(new ModContract3pt(l,tax,false,qh,qh,qSpectator,p0,p0,j,deltaT,t, bHeavyAnti ) );
                    }
                  }
                }
              }
            }
            if( !bDidSomething )
            {
              // We are only performing residual-mass checks on each propagator
              l.TakeOwnership( new ModProp( l, tax, qh, p, t ) );
              //l.TakeOwnership( new ModProp( l, tax, qSpectator, p,t) );
            }
            p = -p;
          }
        }
      }
    }
  }
}

void Study3::Make()
{
  for( const Decay &d : makePar.Decays )
    MakeStudy3( d );
}

int main(int argc, char *argv[])
{
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
  return AppMaker::MakeThreePoint( argc, argv, sXmlFilename, sRunSuffix );
}
