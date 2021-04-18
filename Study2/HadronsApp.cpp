/*************************************************************************************
 
 Create a Hadrons Application, top-down building all required dependencies
 Source file: HadronsApp.hpp
 Copyright (C) 2021
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

const std::string Sep{ "_" };    // used inside filenames
const std::string Space{ " " };  // whitespace as field separator / human readable info
const std::string defaultMomName{ "p" };
const std::string SeqMomentumName{ "ps" };
const Common::Momentum p0(0,0,0);
const std::string GaugeFieldName{"gauge"};

const std::array<std::string, 5> FamilyNames{ "Z2", "ZF", "GF", "GP", "GR" };

typename std::underlying_type<Family>::type FamilyIndex( Family family )
{
  const auto iFamily = static_cast<typename std::underlying_type<Family>::type>( family );
  if( iFamily < 0 || iFamily >= FamilyNames.size() )
    throw std::runtime_error( "Unknown family " + std::to_string( iFamily ) );
  return iFamily;
}

std::istream& operator>>(std::istream& is, Family &family)
{
  std::string s;
  if( is >> s )
  {
    int i;
    for(i = 0; i < FamilyNames.size() && !Common::EqualIgnoreCase( s, FamilyNames[i] ); ++i)
      ;
    if( i < FamilyNames.size() )
      family = static_cast<Family>( i );
    else
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

const std::array<std::string, 3> SpeciesNames{ "Point", "Wall", "Region" };

typename std::underlying_type<Species>::type SpeciesIndex( Species species )
{
  const auto iSpecies = static_cast<typename std::underlying_type<Species>::type>( species );
  if( iSpecies < 0 || iSpecies >= ::SpeciesNames.size() )
    throw std::runtime_error( "Unknown species " + std::to_string( iSpecies ) );
  return iSpecies;
}

std::istream& operator>>(std::istream& is, Species &species)
{
  std::string s;
  if( is >> s )
  {
    int i;
    for(i = 0; i < SpeciesNames.size() && !Common::EqualIgnoreCase( s, SpeciesNames[i] ); ++i)
      ;
    if( i < SpeciesNames.size() )
      species = static_cast<Species>( i );
    else
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

// Taxonomy

const std::string Taxonomy::HissyFitDefault{ "Invalid Taxonomy" };

std::istream& operator>>(std::istream& is, Taxonomy &taxonomy)
{
  Taxonomy local;
  std::string s;
  if( is >> local.family && is >> local.species )
  {
    local.Validate();
    taxonomy = local;
  }
  return is;
}

// These make taxa sortable
bool operator<( const Taxonomy &lhs, const Taxonomy &rhs )
{
  return       FamilyIndex ( lhs.family ) <  FamilyIndex ( rhs.family )
         || (  FamilyIndex ( lhs.family ) == FamilyIndex ( rhs.family )
            && SpeciesIndex( lhs.species) <  SpeciesIndex( rhs.species) );
}

bool operator==( const Taxonomy &lhs, const Taxonomy &rhs )
{
  return FamilyIndex( lhs.family ) == FamilyIndex ( rhs.family )
     && SpeciesIndex( lhs.species) == SpeciesIndex( rhs.species);
}

// Quark

void Quark::CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const
{
  MAction::ScaledDWF::Par Par;
  Par.gauge    = std::move( Gauge );
  Par.Ls       = Ls;
  Par.mass     = mass;
  Par.M5       = M5;
  Par.scale    = scale; // 1 = Shamir, 2 = Möbius
  Par.boundary = boundary;
  Par.twist    = twist;
  app.createModule<MAction::ScaledDWF>( name, Par );
}

/**************************
 Parameters that apply to the entire application being built
**************************/

AppParams::AppParams( const std::string sXmlFilename, const std::string sRunSuffix_ )
: sRunSuffix{sRunSuffix_}
{
  static const std::string sXmlTopLevel{ "Study2" };
  static const std::string sXmlTagName{ "RunPar" };
  XmlReader r( sXmlFilename, false, sXmlTopLevel );
  read( r, sXmlTagName, Run );
  HeavyQuarks        = ReadQuarks( r, "Heavy" );
  SpectatorQuarks    = ReadQuarks( r, "Spectator" );
  SpectatorQuarks2pt = ReadQuarks( r, "SpectatorExtra2pt", 0 );
  Momenta = Common::ArrayFromString<Common::Momentum>( Run.Momenta );
  // Check the taxonomy
  Taxa = Common::ArrayFromString<Taxonomy>( Run.Taxa );
  Common::NoDuplicates( Taxa, "Taxa", 1 );
  for( const Taxonomy &t : Taxa )
    if( t == Family::GR )
      LOG(Warning) << t << " don't swap the gamma structures at sink and source (18-Feb-2021)" << std::endl;
  // Check parameters make sense
  //if( !( Run.TwoPoint || Run.HeavyQuark ||Run.HeavyAnti ) )
    //throw std::runtime_error( "At least one must be true of: TwoPoint; HeavyQuark; or HeavyAnti" );
  if( Momenta.empty() )
    throw std::runtime_error( "There should be at least one momentum" );
  ThreePoint = Run.HeavyQuark || Run.HeavyAnti;
  if( ThreePoint )
  {
    deltaT = Common::ArrayFromString<int>( Run.deltaT );
    Common::NoDuplicates( deltaT, "deltaT", 1 );
    gamma = Common::ArrayFromString<Gamma::Algebra>( Run.gamma );
    Common::NoDuplicates( gamma, "gamma", 1 );
  }
}

std::vector<Quark> AppParams::ReadQuarks( XmlReader &r, const std::string &qType, std::size_t Min )
{
  std::vector<Quark> vq;
  r.readDefault( qType, vq );
  if( vq.size() < Min )
    throw std::runtime_error( "Only " + std::to_string( vq.size() ) + Common::Space + qType + " quarks" );
  return vq;
}

// Make a unique RunID string that completely describes the run
std::string AppParams::RunID() const
{
  std::ostringstream s;
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
  for( const Quark &q : SpectatorQuarks )
    s << Sep << q.flavour;
  if( Run.TwoPoint )
    s << Sep << "2pt";
  if( Run.HeavyQuark )
    s << Sep << "quark";
  if( Run.HeavyAnti )
    s << Sep << "anti";
  for( const Quark &q : HeavyQuarks )
    s << Sep << q.flavour;
  s << Sep << "t" << Sep << Run.Timeslices.start << Sep << Run.Timeslices.end << Sep << Run.Timeslices.step;
  if( !Run.DoNegativeMomenta )
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

/**************************
 List of Hadrons Modules, i.e. a wrapper for Hadrons::Application
 **************************/

// Delete the module if it already exists, otherwise add it to the application and create all dependencies

const std::string HModList::TakeOwnership( HMod *pHMod )
{
  assert( pHMod && "HModList::TakeOwnership() : Null pointer" );
  std::string ModName{ pHMod->Name() };
  auto it{ list.find( ModName ) };
  if( it != list.end() )
    delete pHMod; // Already exists
  else
  {
    // I've not seen this yet - Make all its dependencies and add it to my list
    pHMod->AddDependencies( *this );
    list.emplace( std::make_pair( ModName, std::unique_ptr<HMod>( pHMod ) ) );
  }
  return ModName;
}

/**************************
 Point sink
**************************/

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p_)
: HMod(ModList, taxonomy), p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSink::AddDependencies( HModList &ModList ) const
{
  MSink::ScalarPoint::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
}

/**************************
 Wall sink
**************************/

const std::string ModSinkSmear::Prefix{ "Sink_Smear" };

ModSinkSmear::ModSinkSmear(HModList &ModList, const Taxonomy &taxonomy, const Common::Momentum &p_)
: HMod(ModList, taxonomy), p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSinkSmear::AddDependencies( HModList &ModList ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string( Space );
  ModList.application.createModule<MSink::Point>(Name(), sinkPar);
}

/**************************
 Source
**************************/

const std::string ModSource::Prefix{ "Source" };

ModSource::ModSource(HModList &ModList, const Taxonomy &tx, const Common::Momentum &p_, int t_)
// Silently translate ZF sources into Z2 sources
: HMod(ModList, tx.family == Family::ZF ? Taxonomy(Family::Z2,tx.species) : tx), p{p_}, t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  AppendPT( name, t, p );
}

void ModSource::AddDependencies( HModList &ModList ) const
{
  switch( tax.family )
  {
    case Family::Z2:
      if( p )
      {
        MSource::MomentumPhase::Par par;
        par.src = ModList.TakeOwnership( new ModSource( ModList, tax, Common::Momentum(0,0,0), t ) );
        par.mom = p.to_string4d( Common::Space );
        ModList.application.createModule<MSource::MomentumPhase>(name, par);
      }
      else
      {
        MSource::Z2::Par par;
        par.tA = t;
        par.tB = t;
        ModList.application.createModule<MSource::Z2>(name, par);
      }
      break;
    case Family::GP:
      {
        MSource::Point::Par par;
        par.position = ModList.params.Run.SpatialPos + Common::Space + std::to_string( t );
        ModList.application.createModule<MSource::Point>(name, par);
      }
      break;
    case Family::GF:
      {
        MSource::Wall::Par par;
        par.tW = t;
        par.mom = p.to_string4d( Common::Space );
        ModList.application.createModule<MSource::Wall>(name, par);
      }
      break;
    default:
      {
        std::stringstream ss;
        ss << "Unrecognised family \"" << tax.family << "\"";
        throw std::runtime_error( ss.str() );
      }
  }
}

/**************************
 Stout smeared gauge
**************************/

ModGauge::ModGauge( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared_ )
: HMod(ModList, taxonomy), bSmeared{ bSmeared_ }
{
  name = GaugeFieldName;
  tax.AppendFixed( name, bSmeared );
}

void ModGauge::AddDependencies( HModList &ModList ) const
{
  if( bSmeared )
  {
    // Stout smeared field
    MGauge::StoutSmearing::Par stoutPar( ModList.params.Run.StoutSmear );
    stoutPar.gauge = tax.GaugeFixed() ? ModList.TakeOwnership( new ModGauge( ModList, tax, false ) ) : GaugeFieldName;
    ModList.application.createModule<MGauge::StoutSmearing>(name, stoutPar);
  }
  else if( tax.GaugeFixed() )
  {
    // Gauge-fixed
    MGauge::GaugeFix::Par Par( ModList.params.Run.GaugeFix );
    Par.gauge = GaugeFieldName;
    ModList.application.createModule<MGauge::GaugeFix>(name, Par);
  }
}

/**************************
 Action
**************************/

const std::string ModAction::Prefix{ "DWF" };

ModAction::ModAction( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_ )
: HMod( ModList, taxonomy ), q{q_}, bSmeared{q_.GaugeSmear && !ModList.params.Run.Gauge.empty()}
{
  name = Prefix;
  tax.AppendFixed( name, bSmeared );
  Append( name, q.flavour );
}

void ModAction::AddDependencies( HModList &ModList ) const
{
  std::string Gauge{ModList.TakeOwnership(new ModGauge( ModList, tax, bSmeared))};
  q.CreateAction( ModList.application, name, std::move( Gauge ) );
}

/**************************
 Solver
**************************/

const std::string ModSolver::Prefix{ "CG" };

ModSolver::ModSolver( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_ )
: HMod( ModList, taxonomy ), q( q_ ), bSmeared{q_.GaugeSmear && !ModList.params.Run.Gauge.empty()}
{
  name = Prefix;
  tax.AppendFixed( name, bSmeared );
  Append( name, q.flavour );
}

void ModSolver::AddDependencies( HModList &ModList ) const
{
  // solvers
  MSolver::RBPrecCG::Par solverPar;
  if( ModList.params.Run.Gauge.length() && q.EigenPackFilename.length() )
  {
    // eigenpacks for deflation
    MIO::LoadFermionEigenPack::Par epPar;
    epPar.filestem = q.EigenPackFilename;
    epPar.multiFile = false;
    epPar.size = 600;
    epPar.Ls = q.Ls;
    if( tax.GaugeFixed() )
      epPar.gaugeXform = ModList.TakeOwnership( new ModGauge( ModList, tax, false ) ) + "_xform";
    solverPar.eigenPack = "epack_" + q.flavour;
    ModList.application.createModule<MIO::LoadFermionEigenPack>(solverPar.eigenPack, epPar);
  }
  solverPar.action       = ModList.TakeOwnership( new ModAction( ModList, tax, q ) );
  solverPar.residual     = q.residual;
  solverPar.maxIteration = q.maxIteration;
  ModList.application.createModule<MSolver::RBPrecCG>(name, solverPar);
}

/**************************
 Propagator
**************************/

const std::string ModProp::Prefix{ "Prop" };
const std::string ModProp::PrefixConserved{ "Ward" };

ModProp::ModProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod( ModList, taxonomy ), q{q_}, p{p_}, t{t_}
{
  Suffix = tax.FamilyName();
  Append( Suffix, q.flavour );
  AppendPT( Suffix, t, p );
  // Name of the propagator
  name = Prefix;
  Append( name, Suffix );
}

void ModProp::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par par;
  par.source = ModList.TakeOwnership( new ModSource( ModList, tax, p, t ) );
  par.solver = ModList.TakeOwnership( new ModSolver( ModList, tax, q ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, par);
  // Check residual mass for zero-momentum propagators
  if( !p )
  {
    std::string WardName{ PrefixConserved };
    Append( WardName, Suffix );
    // Check residual mass
    MContraction::WardIdentity::Par WIP;
    WIP.prop = name + "_5d";
    WIP.action = ModList.TakeOwnership( new ModAction( ModList, tax, q ) );
    WIP.source = par.source;
    WIP.mass = q.mass;
    WIP.output = ModList.params.Run.OutputBase + PrefixConserved + "/" + Suffix;
    ModList.application.createModule<MContraction::WardIdentity>(WardName, WIP);
  }
}

/**************************
 Sliced Propagator / wall sink. For 2pt contractions only
**************************/

const std::string ModSlicedProp::Prefix{ "SlicedProp" };

ModSlicedProp::ModSlicedProp( HModList &ML, const Taxonomy &tax_, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod( ML, tax_ ), q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

void ModSlicedProp::AddDependencies( HModList &ModList ) const
{
  MSink::Smear::Par smearPar;
  smearPar.q = ModList.TakeOwnership( new ModProp( ModList, tax, q, p, t ) );
  smearPar.sink = ModList.TakeOwnership( new ModSinkSmear( ModList, tax, p0 ) ); // TODO: Have I put momentum in properly?
  ModList.application.createModule<MSink::Smear>(name, smearPar);
}

/**************************
 Sequential Source
**************************/

const std::string ModSourceSeq::Prefix{ "Seq" + ModSource::Prefix };

ModSourceSeq::ModSourceSeq( HModList &ML, const Taxonomy &tax_, Gamma::Algebra current_, int deltaT_,
                            const Common::Momentum &pSeq_, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod(ML,tax_), Current{current_}, deltaT{deltaT_}, pSeq{pSeq_}, q{q_}, p{p_}, t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, tax.species );
  Append( name, Current );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

template<typename T> void ModSourceSeq::AddDependenciesT( HModList &ModList, typename T::Par &seqPar ) const
{
  seqPar.q = ModList.TakeOwnership( new ModProp( ModList, tax, q, p, t ) );
  seqPar.mom = pSeq.to_string4d( Space );
  seqPar.gamma = Current;
  ModList.application.createModule<T>(name, seqPar);
}

void ModSourceSeq::AddDependencies( HModList &ModList ) const
{
  const int tSink{ModList.params.TimeBound( t + deltaT )};
  switch( tax.species )
  {
    case Species::Point:
    {
      using M = MSource::SeqGamma;
      typename M::Par seqPar;
      seqPar.tA  = tSink;
      seqPar.tB  = tSink;
      AddDependenciesT<M>( ModList, seqPar );
    }
      break;
    case Species::Wall:
    {
      using M = MSource::SeqGammaWall;
      typename M::Par seqPar;
      seqPar.tA  = tSink;
      seqPar.tB  = tSink;
      AddDependenciesT<M>( ModList, seqPar );
    }
      break;
    case Species::Region:
    {
      using M = MSource::SeqGammaRegion;
      typename M::Par seqPar;
      seqPar.LowerLeft  = ModList.params.Run.SpatialPos + Common::Space + std::to_string( tSink );
      seqPar.RegionSize = ModList.params.Run.RegionSize;
      AddDependenciesT<M>( ModList, seqPar );
    }
      break;
  }
}

/**************************
 Sequential propagator
**************************/

const std::string ModPropSeq::Prefix{ "Seq" + ModProp::Prefix };

ModPropSeq::ModPropSeq( HModList &ML, const Taxonomy &tax_, const Quark &qSeq_, Gamma::Algebra current_, int deltaT_,
                        const Common::Momentum &pSeq_, const Quark &q_, const Common::Momentum &p_, int t_ )
: HMod(ML,tax_),qSeq{qSeq_},Current{current_},deltaT{deltaT_},pSeq{pSeq_},q{q_},p{p_},t{t_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, tax.species );
  Append( name, qSeq.flavour );
  Append( name, Current );
  AppendDeltaT( name, deltaT );
  AppendPSeq( name, pSeq );
  Append( name, q.flavour );
  AppendPT( name, t, p );
}

void ModPropSeq::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par quarkPar;
  quarkPar.source = ModList.TakeOwnership( new ModSourceSeq( ModList, tax, Current, deltaT, pSeq, q, p, t ) );
  quarkPar.solver = ModList.TakeOwnership( new ModSolver( ModList, tax, qSeq ) );
  ModList.application.createModule<MFermion::GaugeProp>(name, quarkPar);
}

/**************************
 2pt contraction
**************************/

const std::string ContractionPrefix{ "meson" };

ModContract2pt::ModContract2pt( HModList &ModList, const Taxonomy &taxonomy,
                                const Quark &q1_, const Quark &q2_, const Common::Momentum &p_, int t_)
// Two-point Region-sinks don't really exist - silently translate to point-sinks
: HMod(ModList, taxonomy == Species::Region ? Taxonomy(taxonomy.family, Species::Point) : taxonomy),
  q1{q1_}, q2{q2_}, p{p_}, t{t_}
{
  std::string Prefix{ tax.FilePrefix() };
  std::string s{ q2.flavour };
  Append( s, q1.flavour );
  AppendPT( s, t, p );
  static const std::string Prefixfamily{ "2pt" };
  FileName = ModList.params.Run.OutputBase;
  FileName.append( Prefixfamily );
  FileName.append( 1, '/' );
  FileName.append( Prefix );
  FileName.append( 1, '/' );
  FileName.append( Prefix ); // Repetition not an error. Part of directory hierarchy and filename
  Append( FileName, s );
  name = ContractionPrefix;
  Append( name, Prefixfamily );
  Append( name, Prefix );
  Append( name, s );
}

void ModContract2pt::AddDependencies( HModList &ModList ) const
{
  // Two-point Gauge-fixed reversed don't really exist ... just silently ignore
  if( tax == Family::GR )
    return; // I can't do this for two-point. Just silently ignore
  MContraction::Meson::Par mesPar;
  mesPar.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  mesPar.gammas = "all";
  if( tax == Species::Wall )
  {
    mesPar.q1 = ModList.TakeOwnership( new ModSlicedProp( ModList, tax, q1, p, t ) );
    mesPar.q2 = ModList.TakeOwnership( new ModSlicedProp( ModList, tax, q2, p, t ) );
  }
  else
  {
    mesPar.q1 = ModList.TakeOwnership( new ModProp( ModList, tax, q1, p, t ) );
    mesPar.q2 = ModList.TakeOwnership( new ModProp( ModList, tax, q2, p, t ) );
  }
  mesPar.sink = ModList.TakeOwnership( new ModSink( ModList, tax, p0 ) );
  ModList.application.createModule<MContraction::Meson>(name, mesPar);
}

/**************************
 Three-point contraction
 Quark qSrc at time t, decaying to quark qDst via current insertion at time (t + deltaT), with spectator anti-quark
 NB: if bHeavyAnti is true, then qSrc and qDst are anti-quarks and spectator is quark
 Momentum p is for the source, -p at the current and 0 at the sink (this could change later)
**************************/

ModContract3pt::ModContract3pt(HModList &ML,const Taxonomy &tax_,const Quark &Snk_,const Quark &Src_,const Quark &Spec_,
              const Common::Momentum &pSeq_, const Common::Momentum &p_, int t_, bool bHA_, Gamma::Algebra Cur_, int dT_ )
: HMod(ML, tax_),qSnk{Snk_},qSrc{Src_},qSpec{Spec_},pSeq{pSeq_},p{p_},t{t_},bHeavyAnti{bHA_},Current{Cur_},deltaT{dT_}
{
  std::string Prefix1{ "3pt" };
  Append( Prefix1, qSpec.flavour );
  std::string Prefix2{ tax.FilePrefix() };
  Append( Prefix2, bHeavyAnti ? "anti" : "quark" );
  std::string s{ qSnk.flavour };
  Append( s, qSrc.flavour );
  Append( s, Current );
  AppendDeltaT( s, deltaT );
  AppendP( s, pSeq, SeqMomentumName );
  AppendPT( s, t, p );
  // File name
  FileName = ML.params.Run.OutputBase + Prefix1;
  FileName.append( 1, '/' );
  FileName.append( Prefix2 );
  FileName.append( 1, '/' );
  FileName.append( Prefix2 ); // Repetition not a mistake. This is part of directory hierarchy and filename
  Append( FileName, s );
  // Object name
  name = ContractionPrefix;
  Append( name, Prefix1 );
  Append( name, Prefix2 );
  Append( name, s );
}

void ModContract3pt::AddDependencies( HModList &ML ) const
{
  MContraction::Meson::Par par;
  par.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  par.gammas = "all";
  const bool bInvertSeq{ !bHeavyAnti };
  const bool bRev{ tax == Family::GR };
  // The "s" prefix for the variable names is short for "sub", i.e. the values to use for dependencies
  Taxonomy stax( bRev ? Family::GF : tax.family, bRev ? Species::Region : tax.species );
  const Quark &sqSrc{ bRev ? qSnk : qSrc };
  const Quark &sqSnk{ bRev ? qSrc : qSnk };
  const int st     { bRev ? ML.params.TimeBound( t + deltaT ) : t      };
  const int sdeltaT{ bRev ? ML.params.TimeBound(   - deltaT ) : deltaT };
  const Common::Momentum &pSource( bRev ? pSeq : p );
  const Common::Momentum &pSink( bRev ? p : pSeq );
  Common::Momentum pCurrent( - ( pSource + pSink ) );
  if( bInvertSeq )
  {
    par.q1 = ML.TakeOwnership( new ModProp( ML, stax, sqSrc, pSource, st ) );
    par.q2 = ML.TakeOwnership(new ModPropSeq(ML,stax, sqSnk, Current, sdeltaT, -pSink, qSpec, p0,       st ) );
  }
  else
  {
    par.q1 = ML.TakeOwnership(new ModPropSeq(ML,stax, sqSnk, Current, sdeltaT,  pSink, qSpec, pSource , st ) );
    par.q2 = ML.TakeOwnership( new ModProp( ML, stax, sqSrc, p0, st ) );
  }
  par.sink = ML.TakeOwnership( new ModSink( ML, stax, pCurrent ) );
  ML.application.createModule<MContraction::Meson>(name, par);
}
