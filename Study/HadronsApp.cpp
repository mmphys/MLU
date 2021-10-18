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
const std::string HMod::sSinglePrec{ "float" };

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
int Taxonomy::numHits{ 1 };

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

template<typename T> void Quark::CreateAction( Application &app, const std::string &name, std::string &&Gauge ) const
{
  typename T::Par Par;
  Par.gauge    = std::move( Gauge );
  Par.Ls       = Ls;
  Par.mass     = mass;
  Par.M5       = M5;
  Par.scale    = scale; // 1 = Shamir, 2 = Möbius
  Par.boundary = boundary;
  Par.twist    = twist;
  app.createModule<T>( name, Par );
}

/**************************
 Parameters that apply to the entire application being built
**************************/

AppParams::AppParams( XmlReader &r )
{
  static const std::string sXmlTagName{ "RunPar" };
  read( r, sXmlTagName, Run );
  if( Run.GaugeFixed.empty() && !Run.GaugeFixedXform.empty() )
    throw std::runtime_error( "GaugeFixedXform specified without GaugeFixed" );
  if( !Run.GaugeFixed.empty() && Run.GaugeFixedXform.empty() )
    throw std::runtime_error( "GaugeFixed specified without GaugeFixedXform" );
}

std::vector<std::string> AppParams::GetWarnings() const
{
  std::vector<std::string> w;
  if(  Run.Gauge.empty() &&  Run.GaugeFixed.empty() )
    w.push_back( "Neither Gauge nor GaugeFixed specified. Using Unit gauge" );
  if( !Run.Gauge.empty() && !Run.GaugeFixed.empty() )
    w.push_back( "Both Gauge and GaugeFixed specified. Using GaugeFixed" );
  return w;
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

ModSource::ModSource(HModList &ModList, const Taxonomy &tx, const Common::Momentum &p_, int t_, int hit_)
// Silently translate ZF sources into Z2 sources
: HMod(ModList, tx.family == Family::ZF ? Taxonomy(Family::Z2,tx.species) : tx), p{p_}, t{t_}, hit{hit_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  AppendPT( name, t, p );
  tax.AppendHit( name, hit );
}

void ModSource::AddDependencies( HModList &ModList ) const
{
  switch( tax.family )
  {
    case Family::Z2:
      if( p )
      {
        MSource::MomentumPhase::Par par;
        par.src = ModList.TakeOwnership( new ModSource( ModList, tax, p0, t, hit ) );
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

ModGauge::ModGauge( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared_, Precision precision_ )
: HMod(ModList, taxonomy), bSmeared{ bSmeared_ }, precision{ precision_ }
{
  name = GaugeFieldName;
  tax.AppendFixed( name, bSmeared );
  if( precision == Precision::Single )
    Append( name, HMod::sSinglePrec );
}

void ModGauge::AddDependencies( HModList &ModList ) const
{
  if( precision == Precision::Single )
  {
    MUtilities::GaugeSinglePrecisionCast::Par castPar;
    castPar.field = ModList.TakeOwnership( new ModGauge( ModList, tax, bSmeared, Precision::Double ) );
    ModList.application.createModule<MUtilities::GaugeSinglePrecisionCast>(name, castPar);
  }
  else if( bSmeared )
  {
    // Stout smeared field
    MGauge::StoutSmearing::Par stoutPar( ModList.params.Run.StoutSmear );
    stoutPar.gauge = ModList.TakeOwnership( new ModGauge( ModList, tax, false, precision ) );
    ModList.application.createModule<MGauge::StoutSmearing>(name, stoutPar);
  }
  else if( tax.GaugeFixed() )
  {
    // Gauge-fixed
    if( !ModList.params.Run.GaugeFixed.empty() )
    {
      MIO::LoadNersc::Par gaugePar;
      gaugePar.file = ModList.params.Run.GaugeFixed;
      ModList.application.createModule<MIO::LoadNersc>(name, gaugePar);
    }
    else
    {
      MGauge::GaugeFix::Par Par( ModList.params.Run.GaugeFix );
      Par.gauge = ModList.TakeOwnership( new ModGauge( ModList, Taxonomy( Family::Z2, Species::Point ),
                                                       bSmeared, precision ) );
      ModList.application.createModule<MGauge::GaugeFix>(name, Par);
    }
  }
  else if( !ModList.params.Run.Gauge.empty() )
  {
    MIO::LoadNersc::Par gaugePar;
    gaugePar.file = ModList.params.Run.Gauge;
    ModList.application.createModule<MIO::LoadNersc>( name, gaugePar );
  }
  else
  {
    ModList.application.createModule<MGauge::Unit>( name );
  }
}

/**************************
 Gauge transform
**************************/

const std::string ModGaugeXform::Suffix{ "xform" };

ModGaugeXform::ModGaugeXform( HModList &ModList, const Taxonomy &taxonomy, bool bSmeared_, Precision precision_ )
: HMod(ModList, taxonomy), bSmeared{ bSmeared_ }, precision{ precision_ }
{
  name = GaugeFieldName;
  tax.AppendFixed( name, bSmeared );
  if( precision == Precision::Single )
    Append( name, HMod::sSinglePrec );
  Append( name, Suffix );
}

void ModGaugeXform::AddDependencies( HModList &ModList ) const
{
  if( precision == Precision::Single )
  {
    MUtilities::ColourMatrixSinglePrecisionCast::Par xformPar;
    xformPar.field = ModList.TakeOwnership( new ModGaugeXform( ModList, tax, bSmeared, Precision::Double ) );
    ModList.application.createModule<MUtilities::ColourMatrixSinglePrecisionCast>(name, xformPar);
  }
  else if( bSmeared )
  {
    throw std::runtime_error( "Don't know how to make gauge transform for smeared field" );
  }
  else if( tax.GaugeFixed() )
  {
    const std::string GaugeName{ ModList.TakeOwnership( new ModGauge( ModList, tax, bSmeared, precision ) ) };
    if( !ModList.params.Run.GaugeFixed.empty() )
    {
      // Since we loaded the gauge field, we'll need to load the transform as well
      MIO::LoadColourMatrixField::Par xformPar;
      xformPar.name = name;
      xformPar.Ls = 1;
      xformPar.fileStem = ModList.params.Run.GaugeFixedXform;
      ModList.application.createModule<MIO::LoadColourMatrixField>(name + "_load", xformPar);
    }
  }
  else
  {
    throw std::runtime_error( "No such thing as a gauge transform for un-fixed field" );
  }
}

/**************************
 Action
**************************/

const std::string ModAction::Prefix{ "DWF" };

ModAction::ModAction( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_, Precision precision_ )
: HMod( ModList, taxonomy ), q{q_}, bSmeared{q_.GaugeSmear}, precision{precision_}
{
  name = Prefix;
  tax.AppendFixed( name, bSmeared );
  Append( name, q.flavour );
  if( precision == Precision::Single )
    Append( name, HMod::sSinglePrec );
}

void ModAction::AddDependencies( HModList &ModList ) const
{
  std::string Gauge{ ModList.TakeOwnership( new ModGauge( ModList, tax, bSmeared, precision ) ) };
  if( precision == Precision::Single )
    q.CreateAction<MAction::ScaledDWFF>( ModList.application, name, std::move( Gauge ) );
  else
    q.CreateAction<MAction::ScaledDWF >( ModList.application, name, std::move( Gauge ) );
}

/**************************
 Solver
**************************/

const std::string ModSolver::Prefix{ "CG" };

ModSolver::ModSolver( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_ )
: HMod( ModList, taxonomy ), q( q_ ), bSmeared{q_.GaugeSmear}
{
  name = Prefix;
  tax.AppendFixed( name, bSmeared );
  Append( name, q.flavour );
}

template<typename TEPLoad, typename TGuesser>
std::string ModSolver::LoadEigenPack( HModList &ModList, Precision XformPres ) const
{
  std::string EigenPackName{ "epack_" + name };
  typename TEPLoad::Par epPar;
  epPar.filestem = q.eigenPack;
  epPar.multiFile = q.multiFile;
  epPar.redBlack = q.redBlack;
  epPar.size = q.size;
  epPar.Ls = q.Ls;
  if( tax.GaugeFixed() )
  {
    epPar.gaugeXform = ModList.TakeOwnership( new ModGaugeXform( ModList, tax, q.GaugeSmear, XformPres ) );
  }
#ifdef MLU_HADRONS_HAS_GUESSERS
  epPar.redBlack = true;
#endif
  ModList.application.createModule<TEPLoad>( EigenPackName, epPar );
#ifdef MLU_HADRONS_HAS_GUESSERS
  std::string GuesserName{ "guesser_" + name };
  typename TGuesser::Par guessPar;
  guessPar.eigenPack = EigenPackName;
  guessPar.size = epPar.size;
  ModList.application.createModule<TGuesser>( GuesserName, guessPar );
  return GuesserName;
#else
  return EigenPackName;
#endif
}

void ModSolver::AddDependencies( HModList &ModList ) const
{
  std::string EigenPackName;
  if( ModList.params.Run.Gauge.length() && q.eigenPack.length() )
  {
    using TDoubleEP             = MIO::LoadFermionEigenPack; // Double-precision Eigen packs can be used anywhere
    using TSingleEP             = MIO::LoadFermionEigenPackF; // Single-precision Eigen packs for MP solver
    using TLoadSingleAsDoubleEP = MIO::LoadFermionEigenPackIo32; // Load single-precision as double-
#ifdef MLU_HADRONS_HAS_GUESSERS
    using TGuesserD = MGuesser::ExactDeflation;
    using TGuesserF = MGuesser::ExactDeflationF;
#else
    using TGuesserD = void;
    using TGuesserF = void;
#endif
    if( !q.eigenSinglePrecision )
      EigenPackName = LoadEigenPack<TDoubleEP,             TGuesserD>( ModList, Precision::Double );
    else if( q.MixedPrecision() )
      EigenPackName = LoadEigenPack<TSingleEP,             TGuesserF>( ModList, Precision::Single );
    else
      EigenPackName = LoadEigenPack<TLoadSingleAsDoubleEP, TGuesserD>( ModList, Precision::Double );
  }
  if( q.MixedPrecision() )
  {
    using T = MSolver::MixedPrecisionRBPrecCG;
    typename T::Par solverPar;
#ifdef MLU_HADRONS_HAS_GUESSERS
    solverPar.outerGuesser      = EigenPackName;
#else
    solverPar.eigenPack         = EigenPackName;
#endif
    solverPar.residual          = q.residual;
    solverPar.maxInnerIteration = q.maxIteration;
    solverPar.maxOuterIteration = q.maxOuterIteration;
    solverPar.innerAction       = ModList.TakeOwnership( new ModAction( ModList, tax, q, Precision::Single ) );
    solverPar.outerAction       = ModList.TakeOwnership( new ModAction( ModList, tax, q, Precision::Double ) );
    ModList.application.createModule<T>(name, solverPar);
  }
  else
  {
    using T = MSolver::RBPrecCG;
    typename T::Par solverPar;
#ifdef MLU_HADRONS_HAS_GUESSERS
    solverPar.guesser      = EigenPackName;
#else
    solverPar.eigenPack    = EigenPackName;
#endif
    solverPar.residual     = q.residual;
    solverPar.maxIteration = q.maxIteration;
    solverPar.action       = ModList.TakeOwnership( new ModAction( ModList, tax, q, Precision::Double ) );
    ModList.application.createModule<T>(name, solverPar);
  }
}

/**************************
 Propagator
**************************/

const std::string ModProp::Prefix{ "Prop" };
const std::string ModProp::PrefixConserved{ "Ward" };

ModProp::ModProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_, const Common::Momentum &p_,
                  int t_, int hit_ )
: HMod( ModList, taxonomy ), q{q_}, p{p_}, t{t_}, hit{hit_}
{
  Suffix = tax.FamilyName();
  Append( Suffix, q.flavour );
  AppendPT( Suffix, t, p );
  tax.AppendHit( Suffix, hit );
  // Name of the propagator
  name = Prefix;
  Append( name, Suffix );
}

void ModProp::AddDependencies( HModList &ModList ) const
{
  MFermion::GaugeProp::Par par;
  par.source = ModList.TakeOwnership( new ModSource( ModList, tax, p, t, hit ) );
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
    WIP.action = ModList.TakeOwnership( new ModAction( ModList, tax, q, Precision::Double ) );
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

ModSlicedProp::ModSlicedProp( HModList &ML, const Taxonomy &tax_, const Quark &q_, const Common::Momentum &p_,
                              int t_, int hit_ )
: HMod( ML, tax_ ), q{q_}, p{p_}, t{t_}, hit{hit_}
{
  name = Prefix;
  Append( name, tax.FamilyName() );
  Append( name, q.flavour );
  AppendPT( name, t, p );
  tax.AppendHit( name, hit );
}

void ModSlicedProp::AddDependencies( HModList &ModList ) const
{
  MSink::Smear::Par smearPar;
  smearPar.q = ModList.TakeOwnership( new ModProp( ModList, tax, q, p, t, hit ) );
  smearPar.sink = ModList.TakeOwnership( new ModSinkSmear( ModList, tax, -p ) );
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
 By default, momentum +p1 at source and -p1 at sink
 However, can also specify p2 for quark 2
    i.e. source momentum = (p1 - p2), sink momentum = (p2 - p1)
**************************/

const std::string ContractionPrefix{ "meson" };

ModContract2pt::ModContract2pt( HModList &ModList, const Taxonomy &taxonomy,
                                const Quark &q1_, const Quark &q2_, const Common::Momentum &p1_, int t_,
                                const Common::Momentum &p2_, int hit_ )
// Two-point Region-sinks don't really exist - silently translate to point-sinks
: HMod(ModList, taxonomy == Species::Region ? Taxonomy(taxonomy.family, Species::Point) : taxonomy),
q1{q1_}, q2{q2_}, p1{p1_}, t{t_}, p2{p2_}, pSource{p1_ - p2_}, hit{hit_}
{
  std::string Prefix2a( tax.FamilyName() );
  std::string Prefix2b;
  tax.SinkSourceType( Prefix2b );
  std::string s{ q2.flavour };
  Append( s, q1.flavour );
  AppendPT( s, t, pSource );
  if( p2 )
    AppendP( s, p2, "pq2" );
  tax.AppendHit( s, hit );
  static const std::string Prefixfamily{ "2pt" };
  FileName = ModList.params.Run.OutputBase;
  FileName.append( Prefixfamily );
  FileName.append( 1, '/' );
  FileName.append( Prefix2a );
  FileName.append( Prefix2b );
  FileName.append( 1, '/' );
  FileName.append( Prefix2b ); // Repetition not an error. Part of directory hierarchy and filename
  Append( FileName, s );
  name = ContractionPrefix;
  Append( name, Prefixfamily );
  Append( name, Prefix2a );
  name.append( Prefix2b );
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
  const bool bWallSink{ tax == Species::Wall };
  if( bWallSink )
  {
    mesPar.q1 = ModList.TakeOwnership( new ModSlicedProp( ModList, tax, q1, p1, t, hit ) );
    mesPar.q2 = ModList.TakeOwnership( new ModSlicedProp( ModList, tax, q2, p2, t, hit ) );
  }
  else
  {
    mesPar.q1 = ModList.TakeOwnership( new ModProp( ModList, tax, q1, p1, t, hit ) );
    mesPar.q2 = ModList.TakeOwnership( new ModProp( ModList, tax, q2, p2, t, hit ) );
  }
  mesPar.sink = ModList.TakeOwnership( new ModSink( ModList, tax, bWallSink ? p0 : -pSource ) );
  ModList.application.createModule<MContraction::Meson>(name, mesPar);
}

/**************************
 Three-point contraction
 Quark qSrc at time t, decaying to quark qDst via current insertion at time (t + deltaT), with spectator anti-quark
 NB: if bHeavyAnti is true, then qSrc and qDst are anti-quarks and spectator is quark
 Momentum p is for the source, pSeq at the sink, and p_current set so that total is zero
 bRev true means computed starting at sink, not source, so gamma insertion is for source, not sink
  ... otherwise everything else about the filename and parameters has the same meaning
 Design requirement: bRev=true should compute exactly the same function as bRev=false
**************************/

ModContract3pt::ModContract3pt( HModList &ML,const Taxonomy &tax_, bool bRev_,const Quark &Snk_,const Quark &Src_,
                                const Quark &Spec_, const Common::Momentum &pSeq_, const Common::Momentum &p_,
                                Gamma::Algebra Cur_, int dT_, int t_, bool bHA_ )
: HMod(ML, tax_), bReverse{tax == Family::GR ? !bRev_ : bRev_}, qSnk{Snk_}, qSrc{Src_}, qSpec{Spec_},
  pSeq{pSeq_}, p{p_}, Current{Cur_}, deltaT{dT_}, t{t_}, bHeavyAnti{bHA_}
{
  std::string Prefix1{ "3pt" };
  Append( Prefix1, qSpec.flavour );
  std::string Prefix2a( tax.FamilyName() );
  std::string Prefix2b;
  tax.SinkSourceType( Prefix2b, bReverse );
  Append( Prefix2b, bHeavyAnti ? "anti" : "quark" );
  if( bReverse )
    Append( Prefix2b, "rev" );
  std::string s{ qSnk.flavour };
  Append( s, qSrc.flavour );
  Append( s, Current );
  AppendDeltaT( s, deltaT );
  AppendP( s, pSeq, SeqMomentumName );
  AppendPT( s, t, p );
  // File name
  FileName = ML.params.Run.OutputBase + Prefix1;
  FileName.append( 1, '/' );
  FileName.append( Prefix2a );
  FileName.append( Prefix2b );
  FileName.append( 1, '/' );
  FileName.append( Prefix2b ); // Repetition not a mistake. This is part of directory hierarchy and filename
  Append( FileName, s );
  // Object name
  name = ContractionPrefix;
  Append( name, Prefix1 );
  Append( name, Prefix2a );
  name.append(  Prefix2b );
  Append( name, s );
}

void ModContract3pt::AddDependencies( HModList &ML ) const
{
  MContraction::Meson::Par par;
  par.output = FileName;
  //mesPar.gammas = "(Gamma5 Gamma5)(Gamma5 GammaTGamma5)(GammaTGamma5 Gamma5)(GammaTGamma5 GammaTGamma5)";
  par.gammas = "all";
  const bool bInvertSeq{ bReverse ? bHeavyAnti : !bHeavyAnti };
  // The "s" prefix for the variable names is short for "sub", i.e. the values to use for dependencies
  Taxonomy stax( tax == Family::GR ? Family::GF : tax.family, tax == Family::GR ? Species::Region : tax.species );
  const Quark &sqSrc{ bReverse ? qSnk : qSrc };
  const Quark &sqSnk{ bReverse ? qSrc : qSnk };
  const int st      { bReverse ? ML.params.TimeBound( t + deltaT ) : t      };
  const int sdeltaT { bReverse ? ML.params.TimeBound(   - deltaT ) : deltaT };
  const Common::Momentum &pSource( bReverse ? pSeq : p );
  const Common::Momentum &pSink  ( bReverse ? p : pSeq );
  Common::Momentum pCurrent( - ( pSource + pSink ) );
  if( bInvertSeq )
  {
    par.q1 = ML.TakeOwnership( new ModProp( ML, stax, sqSrc, pSource, st ) );
    par.q2 = ML.TakeOwnership(new ModPropSeq(ML,stax, sqSnk, Current, sdeltaT, -pSink, qSpec, p0,  st ) );
  }
  else
  {
    par.q1 = ML.TakeOwnership(new ModPropSeq(ML,stax, sqSnk, Current, sdeltaT,  pSink, qSpec, p0 , st ) );
    par.q2 = ML.TakeOwnership( new ModProp( ML, stax, sqSrc, -pSource, st ) );
  }
  par.sink = ML.TakeOwnership( new ModSink( ML, stax, pCurrent ) );
  ML.application.createModule<MContraction::Meson>(name, par);
}
