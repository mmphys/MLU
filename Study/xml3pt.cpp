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
      //q.Validate();
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
void AppMaker::SetupBase( XmlReader &r, const std::string &sRunSuffix, Grid::GridCartesian * grid )
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
  // result database
  globalPar.database.resultDb = appPar.Run.dbOptions.resultDb;
  Grid::Hadrons::makeFileDir( globalPar.database.resultDb );
  // Specify non-zero sampling period to make statistics dB
  globalPar.database.makeStatDb     = appPar.Run.dbOptions.statDbPeriodMs;
  globalPar.database.statDbPeriodMs = appPar.Run.dbOptions.statDbPeriodMs;
  // Make the application database if requested
  if( appPar.Run.dbOptions.enable )
  {
    globalPar.database.applicationDb = appPar.Run.dbOptions.applicationDbPrefix + RunIDSuffix + ".db";
    Grid::Hadrons::makeFileDir( globalPar.database.applicationDb, grid );
  }
  globalPar.database.restoreModules = false;
  globalPar.database.restoreMemoryProfile = false;
  globalPar.database.restoreSchedule = false;
  globalPar.scheduleFile = JobFilePrefix + ".sched";
  globalPar.graphFile = JobFilePrefix + ".graph";
  globalPar.saveSchedule = true;
  application.setPar ( globalPar );
}

#define APP_MAKER_GLUE( ClassName, StudyName ) \
{ \
  if( reader->push( StudyName ) ) \
  { \
    reader->pop(); \
    LOG(Message) << "Making " << StudyName << std::endl; \
    x.reset( new ClassName( appParams, StudyName ) ); \
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
  LOG(Message) << MLUVersionInfoHuman() << std::endl;
  // Make a default spacetime grid ... just in case we need to communicate
  const Grid::Coordinate &gridDefaultLatt{ GridDefaultLatt() };
  Grid::GridCartesian * grid = Grid::SpaceTimeGrid::makeFourDimGrid( gridDefaultLatt,
                                                                     GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                     GridDefaultMpi() );
  int iReturn{ EXIT_FAILURE };
  try
  {
    static const std::string sXmlTopLevel{ "xml3pt" };
    std::unique_ptr<XmlReader> reader( new XmlReader( sXmlFilename, false, sXmlTopLevel ) );
    const AppParams appParams( *reader );
    for( const std::string &warning : appParams.GetWarnings() )
      LOG(Warning) << warning << std::endl;
    int StudyType{ -1 };
    std::string StudyName;
    read( *reader, "StudyType", StudyType );
    read( *reader, "StudyName", StudyName );
    std::string StudyDescription;
    {
      std::ostringstream os;
      os << "Study " << StudyName << " [Type " << StudyType << "]";
      StudyDescription = os.str();
    }
    std::unique_ptr<AppMaker> x;
    switch( StudyType )
    {
      case 2:
        APP_MAKER_GLUE( Study2, StudyName )
        break;
      case 3:
        APP_MAKER_GLUE( Study3, StudyName )
        break;
    }
    if( !x.get() )
      throw std::runtime_error( StudyDescription + " not recognised" );
    x->SetupBase( *reader, sRunSuffix, grid );
    reader.reset( nullptr );
    x->Make();
    // Run or save the job
    x->application.saveParameterFile( x->JobFilePrefix + ".xml" );
    if( x->Run )
    {
      const int GridNt = gridDefaultLatt[Tdir];
      if( x->appPar.Run.Nt != GridNt )
      {
        std::ostringstream os;
        os << "Cannot run: Job Nt " << x->appPar.Run.Nt << " != Grid Nt " << GridNt;
        throw std::runtime_error( os.str() );
      }
      x->application.run();
    }
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
                l.TakeOwnership( new ModContract2pt( l, tax, *qSpectator, *qH1, p, t, p0, 1 ) );
                l.TakeOwnership( new ModContract2pt( l, tax, *qH1, *qSpectator, p, t, p0, 1 ) );
                if( bFirstSpec )
                {
                  // Additional spectators for 2pt functions only
                  for( const Quark *qSpec2 : SpectatorQuarks2pt )
                  {
                    l.TakeOwnership( new ModContract2pt( l, tax, *qSpec2, *qH1, p, t, p0, 1 ) );
                    l.TakeOwnership( new ModContract2pt( l, tax, *qH1, *qSpec2, p, t, p0, 1 ) );
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
  // Number of hits is optional - because I added this after performing production runs
  Taxonomy::NumHits( Common::FromString<int>( makePar.NumHits, 1 ) );
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
  CountHeavyMomenta = 0;
  for( Decay &d : makePar.Decays )
  {
    Common::Trim( d.qLight );
    Common::Trim( d.qSpectator );
    LOG(Message) << "Decay " << d.name << ": light=" << d.qLight << ", spectator=" << d.qSpectator
                 << ", " << d.HeavyMom.size() << " heavy quarks" << std::endl;
    Q.at( d.qLight );
    Q.at( d.qSpectator );
    if( d.HeavyMom.empty() )
      throw std::runtime_error( "There should be at least 1 heavy quark" );
    for( HeavyMomenta &hp : d.HeavyMom )
    {
      Common::Trim( hp.qHeavy );
      Common::Trim( hp.Momenta );
      LOG(Message) << " heavy " << hp.qHeavy << ": momenta " << hp.Momenta << std::endl;
      const std::string & qh{ Q.at( hp.qHeavy ).flavour };
      if( HeavyQuarks.find( qh ) == HeavyQuarks.end() )
        HeavyQuarks.emplace( qh );
      std::vector<Common::Momentum> Momenta = Common::ArrayFromString<Common::Momentum>( hp.Momenta );
      Check0Negate( Momenta, makePar.DoNegativeMomenta );
      Common::NoDuplicates( Momenta, "HeavyMomenta", 1 );
      CountHeavyMomenta += Momenta.size();
      for( const Common::Momentum &p : Momenta )
      {
        if( UniqueP2.find( p.p2() ) == UniqueP2.end() )
          UniqueP2.emplace( p.p2() );
      }
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
  for( const auto &q : HeavyQuarks )
    s << Sep << q;
  s << Sep << "t" << Sep << makePar.Timeslices.start << Sep << makePar.Timeslices.end << Sep << makePar.Timeslices.step;
  if( makePar.DoNegativeMomenta )
    s << Sep << "neg";
  s << Sep << "hp" << Sep << CountHeavyMomenta << Sep << "p2";
  for( int p2 : UniqueP2 )
    s << Sep << p2;
  if( ThreePoint )
  {
    s << Sep << "dT";
    for( int dT : deltaT )
      s << Sep << std::to_string( dT );
  }
  return s.str();
}

// Get the list of momenta to use for Q2 when we perform two-point contractions
// Copy the momenta, add zero-momentum and (optionally) add the negative momenta. De-duplicate the list
bool Study3::Check0Negate( std::vector<Common::Momentum> &Momenta, bool bNegate )
{
  bool bGot0{ false };
  const std::size_t Num{ Momenta.size() };
  if( bNegate )
    Momenta.reserve( 2 * Num );
  for( std::size_t i = 0; i < Num; ++i )
  {
    if( !Momenta[i] )
      bGot0 = true;
    else if( bNegate )
      Momenta.emplace_back( -Momenta[i] );
  }
  return bGot0;
}

void Study3::Contract2pt( const Taxonomy &tax, const Quark &q1, const Quark &q2, int idxMomentum, int t,
                          const std::vector<Common::Momentum> &Momenta, bool bGotp0 )
{
  // Only the Z2 sources are protected by the gauge-average delta
  const bool bMultiP2{ tax == Family::Z2 || tax == Family::ZF };
  const Common::Momentum &p{ Momenta[idxMomentum] };
  const bool bFlavourDiagonal{ !Common::CompareIgnoreCase( q1.flavour, q2.flavour ) };
  const int idxMomentumStart{ !bMultiP2 || !bGotp0 ? -1 : 0 };
  const int idxMomentumLimit{ !bMultiP2 ?  0
                              : bFlavourDiagonal ? idxMomentum + 1 : static_cast<int>( Momenta.size() ) };
  for( int idxMomentum2 = idxMomentumStart; idxMomentum2 < idxMomentumLimit; ++idxMomentum2 )
  {
    for( int hit = 1; hit <= tax.NumHits(); ++hit )
    {
      const Common::Momentum &p2{ idxMomentum2 < 0 ? p0 : Momenta[idxMomentum2] };
      l.TakeOwnership( new ModContract2pt( l, tax, q1, q2, p, t, p2, hit ) );
      if( !bFlavourDiagonal )
        l.TakeOwnership( new ModContract2pt( l, tax, q2, q1, p, t, p2, hit ) );
    }
  }
}

void Study3::MakeStudy3( unsigned int t, const Decay &d )
{
  const Quark &ql{ Q.at( d.qLight ) };
  const Quark &qSpectator{ Q.at( d.qSpectator ) };
  {
    for( const Taxonomy &tax : Taxa )
    {
      for( std::size_t idxHP = 0; idxHP < d.HeavyMom.size(); ++idxHP )
      {
        const Quark &qh{ Q.at( d.HeavyMom[idxHP].qHeavy ) };
        std::vector<Common::Momentum> Momenta = Common::ArrayFromString<Common::Momentum>( d.HeavyMom[idxHP].Momenta );
        const bool bGotp0{ Check0Negate( Momenta, makePar.DoNegativeMomenta ) };
        for( int idxMomentum = 0; idxMomentum < Momenta.size(); ++idxMomentum )
        {
            const Common::Momentum &p{ Momenta[idxMomentum] };
            if( makePar.TwoPoint )
            {
              // Create the two-point functions I need for R1. Second pair not really needed, but doesn't take much time
              Contract2pt( tax, ql, qSpectator, idxMomentum, t, Momenta, bGotp0 );
              // heavy-spec only needed at 0 momentum. Causes p!=0 propagators to be built, but doesn't take much time
              Contract2pt( tax, qh, qSpectator, idxMomentum, t, Momenta, bGotp0 );
              // Discussion Tobi 9 Jun 2021 - just contract every combination
              Contract2pt( tax, qSpectator, qSpectator, idxMomentum, t, Momenta, bGotp0 );
              Contract2pt( tax, ql, ql, idxMomentum, t, Momenta, bGotp0 );
              for( std::size_t idxHP2 = 0; idxHP2 <= idxHP; ++idxHP2 )
                Contract2pt( tax, qh, Q.at( d.HeavyMom[idxHP2].qHeavy ), idxMomentum, t, Momenta, bGotp0 );
            }
            if( ThreePoint )
            {
              // Create the three-point functions I need for R2 (which includes what's required for R1)
              for( int iHeavy  = makePar.HeavyQuark ? 0 : 1;
                       iHeavy <= makePar.HeavyAnti  ? 1 : 0; ++iHeavy )
              {
                const bool bHeavyAnti{ static_cast<bool>( iHeavy ) };
                for( int deltaT : deltaT )
                {
                  for( Gamma::Algebra j : gamma )
                  {
                    if( makePar.R1Term1 )
                    {
                      const bool bHLB{ makePar.R1Term1Backwards };
                      const unsigned int tHLB{ bHLB ? l.params.TimeBound( t - deltaT ) : t };
                      const Common::Momentum pHLB{ makePar.R2Terms ? (bHLB ? -p : -p) : ( bHLB ? -p : -p ) };
                      l.TakeOwnership(new ModContract3pt(l,tax,bHLB,ql,qh,qSpectator,pHLB,p0,j,deltaT,tHLB,bHeavyAnti));
                    }
                    if( makePar.R1Term2 ) // In practice, R2Terms should equal !R1Term1Backwards
                      l.TakeOwnership(new ModContract3pt(l,tax,false,qh,ql, qSpectator, p0, p, j, deltaT, t, bHeavyAnti ));
                    // I'll always need the zero-momentum versions of l-l and h-h for Z_V ... in addition to R2
                    if( !p || makePar.R2Terms ) // In practice, R2Terms should equal !R1Term1Backwards
                      l.TakeOwnership(new ModContract3pt(l,tax,false,ql,ql, qSpectator,-p, p, j, deltaT, t, bHeavyAnti ) );
                    if( !p ) // Don't repeat this for non-zero momenta if job being performed separately
                    {
                      l.TakeOwnership(new ModContract3pt(l,tax,false,qh,qh, qSpectator,p0, p0,j, deltaT, t, bHeavyAnti ) );
                      l.TakeOwnership(new ModContract3pt(l,tax,false,qh,qh, qh        ,p0, p0,j, deltaT, t, bHeavyAnti));
                      l.TakeOwnership(new ModContract3pt(l,tax,false,ql,ql, qh        ,p0, p0,j, deltaT, t, bHeavyAnti));
                      l.TakeOwnership(new ModContract3pt(l,tax,false,qSpectator,qSpectator,qh,p0,p0,j,deltaT,t,bHeavyAnti));
                    }
                  }
                }
              }
            }
            if( !makePar.TwoPoint && !ThreePoint && !p )
            {
              // We are only performing residual-mass checks on each propagator
              l.TakeOwnership( new ModProp( l, tax, qh, p, t ) );
              l.TakeOwnership( new ModProp( l, tax, ql, p, t ) );
              l.TakeOwnership( new ModProp( l, tax, qSpectator, p,t) );
            }
        }
      }
    }
  }
}

void Study3::Make()
{
  for(unsigned int t = makePar.Timeslices.start; t < makePar.Timeslices.end; t += makePar.Timeslices.step)
  {
    for( const Decay &d : makePar.Decays )
      MakeStudy3( t, d );
    l.Write();
  }
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
