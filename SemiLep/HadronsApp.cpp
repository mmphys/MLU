/**
 
 Create a Hadrons Application, top-down building all required dependencies

 Source file: HadronsApp.cpp

 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Semileptonic Data Generation (SemiLep).
 
 SemiLep is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 SemiLep is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <SemiLepConfig.h>
#include "HadronsApp.hpp"

const std::string Sep{ "_" };    // used inside filenames
const std::string Space{ " " };  // whitespace as field separator / human readable info
const std::string defaultMomName{ "p" };
const std::string SeqMomentumName{ "ps" };
const MLU::Momentum p0(0,0,0);
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
    for(i = 0; i < FamilyNames.size() && !MLU::EqualIgnoreCase( s, FamilyNames[i] ); ++i)
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
    for(i = 0; i < SpeciesNames.size() && !MLU::EqualIgnoreCase( s, SpeciesNames[i] ); ++i)
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
#ifdef MLU_HADRONS_HAS_GUESSERS
  if( Run.BatchSize )
  {
    std::string sError{ "BatchSize = " + std::to_string( Run.BatchSize ) + " but " };
    if( !Run.evBatchSize )
      throw std::runtime_error( sError + "evBatchSize = 0" );
    if( !Run.sourceBatchSize )
      throw std::runtime_error( sError + "sourceBatchSize = 0" );
  }
  else
  {
    std::string sError{ "BatchSize = 0 but " };
    if( Run.evBatchSize )
      throw std::runtime_error( sError + "evBatchSize = " + std::to_string( Run.evBatchSize ) );
    if( Run.sourceBatchSize )
      throw std::runtime_error( sError + "sourceBatchSize = " + std::to_string( Run.sourceBatchSize ) );
  }
#endif
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
 Multi-propagator, i.e. a list of propagators to solve at the same time
 **************************/

#ifdef MLU_HADRONS_HAS_GUESSERS
const std::string MultiProp::sPack{ "Pack" };
const std::string MultiProp::sProp{ "MultiProp" };
const std::string MultiProp::sUnpack{ "Unpack" };

MultiProp::MultiProp(const std::string &Name, const Taxonomy &taxonomy, const Quark &quark, std::string solver)
: HMod(taxonomy), q{quark}, Solver{solver}
{
  name = Name;
  q.AppendName( name, taxonomy );
}

// Make a propagator as part of a multi-propagator

std::string MultiProp::Add( HModList &ModList, const std::string &Source, const std::string &PropName )
{
  Map::iterator it{ m.find( Source ) };
  if( it == m.end() )
  {
    // Add this to the list
    if( !bDirty )
    {
      SourceList.resize( SourceList.size() + 1 );
      bDirty = true;
    }
    ID NewID;
    NewID.Batch = static_cast<unsigned int>( SourceList.size() - 1 );
    NewID.Item  = static_cast<unsigned int>( SourceList[NewID.Batch].size() );
    SourceList[NewID.Batch].push_back( SourceProp( Source, PropName ) );
    it = m.emplace( std::make_pair( Source, NewID ) ).first;
    if( NewID.Item >= ModList.params.Run.BatchSize - 1 )
      Write( ModList ); // Write out this latest batch
  }
  return SourceList[it->second.Batch][it->second.Item].Prop;
}

// Write out the last Num entries as a packed solver
void MultiProp::Write( HModList &ModList )
{
  if( bDirty )
  {
    bDirty = false;
    int Batch = static_cast<unsigned int>( SourceList.size() - 1 );
    std::string Suffix{ name };
    Append( Suffix, Batch );
    // Write the vector packer
    const std::string namePack{ sPack + Suffix };
    MUtilities::PropagatorVectorPackRef::Par parPack;
    for( const SourceProp &sp : SourceList.back() )
      parPack.fields.push_back( sp.Source );
    ModList.application.createModule<MUtilities::PropagatorVectorPackRef>( namePack, parPack );
    // Write the multi RHS propagator
    const std::string nameProp{ sProp + Suffix };
    MFermion::GaugeProp::Par parProp;
    parProp.source = namePack;
    parProp.solver = Solver;
    ModList.application.createModule<MFermion::GaugeProp>( nameProp, parProp );
    // Write the unpacker
    const std::string nameUnpack{ sUnpack + Suffix };
    MUtilities::PropagatorVectorUnpack::Par parUnpack;
    parUnpack.input = nameProp;
    parUnpack.size = static_cast<unsigned int>( SourceList.back().size() );
    for( const SourceProp &sp : SourceList.back() )
      parUnpack.fields.push_back( sp.Prop );
    ModList.application.createModule<MUtilities::PropagatorVectorUnpack>( nameUnpack, parUnpack );
  }
}

/**************************
 Multi-propagator Map, One Multi-propagator per solver
 **************************/

std::string MultiPropMap::Add( HModList &ML, const Taxonomy &tax, const Quark &q, const std::string &Source,
                               const std::string &PropName )
{
  std::string name;
  q.AppendName( name, tax );
  Map::iterator it{ find( name ) };
  if( it == end() )
  {
    std::string sSolver{ ML.TakeOwnership( new ModSolver( ML, tax, q ) ) };
    it = emplace( std::make_pair( name, MultiProp( Name, tax, q, sSolver ) ) ).first;
  }
  return it->second.Add( ML, Source, PropName );
}

void MultiPropMap::Write( HModList &ModList )
{
  for( Map::iterator it = begin(); it != end(); ++it )
    it->second.Write( ModList );
}
#endif // MLU_HADRONS_HAS_GUESSERS

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

std::string HModList::MakeProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q,
                                const MLU::Momentum &p, int t, int hit )
{
  std::string PropName;
#ifdef MLU_HADRONS_HAS_GUESSERS
  if( !q.eigenPack.empty() && ModList.params.Run.BatchSize )
  {
    std::string Source{ ModList.TakeOwnership( new ModSource( ModList, taxonomy, p, t, hit ) ) };
    PropName = Prop.Add( ModList, taxonomy, q, Source, ModProp( ModList, taxonomy, q, p, t, hit ).Name() );
  }
  else
#endif
  {
    PropName = ModList.TakeOwnership( new ModProp( ModList, taxonomy, q, p, t, hit ) );
  }
  return PropName;
}

std::string HModList::MakePropSeq( HModList &ML, const Taxonomy &tax, const Quark &qSeq,
                                   Gamma::Algebra current, int deltaT, const MLU::Momentum &pSeq,
                                   const Quark &q, const MLU::Momentum &p, int t )
{
  std::string PropName;
#ifdef MLU_HADRONS_HAS_GUESSERS
  if( !qSeq.eigenPack.empty() && ML.params.Run.BatchSize )
  {
    std::string Source{ ML.TakeOwnership( new ModSourceSeq( ML, tax, current, deltaT, pSeq, q, p, t ) ) };
    PropName = PropSeq.Add( ML, tax, qSeq, Source, ModPropSeq(ML,tax,qSeq,current,deltaT,pSeq,q,p,t).Name() );
  }
  else
#endif
  {
    PropName = ML.TakeOwnership( new ModPropSeq( ML, tax, qSeq, current, deltaT, pSeq, q, p, t ) );
  }
  return PropName;
}

/**************************
 Point sink
**************************/

const std::string ModSink::Prefix{ "Sink" };

ModSink::ModSink(HModList &ModList, const Taxonomy &taxonomy, const MLU::Momentum &p_)
: HMod(taxonomy), p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSink::AddDependencies( HModList &ModList ) const
{
  MSink::ScalarPoint::Par sinkPar;
  sinkPar.mom = p.to_string3d( Space );
  ModList.application.createModule<MSink::ScalarPoint>(Name(), sinkPar);
}

/**************************
 Wall sink
**************************/

const std::string ModSinkSmear::Prefix{ "Sink_Smear" };

ModSinkSmear::ModSinkSmear(HModList &ModList, const Taxonomy &taxonomy, const MLU::Momentum &p_)
: HMod(taxonomy), p{p_}
{
  name = Prefix;
  AppendP( name, p );
}

void ModSinkSmear::AddDependencies( HModList &ModList ) const
{
  MSink::Point::Par sinkPar;
  sinkPar.mom = p.to_string3d( Space );
  ModList.application.createModule<MSink::Point>(Name(), sinkPar);
}

/**************************
 Source
**************************/

const std::string ModSource::Prefix{ "Source" };

ModSource::ModSource(HModList &ModList, const Taxonomy &tx, const MLU::Momentum &p_, int t_, int hit_)
// Silently translate ZF sources into Z2 sources
: HMod(tx.family == Family::ZF ? Taxonomy(Family::Z2,tx.species) : tx), p{p_}, t{t_}, hit{hit_}
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
        par.mom = p.to_string4d( MLU::Space );
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
        par.position = ModList.params.Run.SpatialPos + MLU::Space + std::to_string( t );
        ModList.application.createModule<MSource::Point>(name, par);
      }
      break;
    case Family::GF:
      {
        MSource::Wall::Par par;
        par.tW = t;
        par.mom = p.to_string4d( MLU::Space );
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
: HMod(taxonomy), bSmeared{ bSmeared_ }, precision{ precision_ }
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
: HMod(taxonomy), bSmeared{ bSmeared_ }, precision{ precision_ }
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
: HMod(taxonomy), q{q_}, precision{precision_}
{
  name = Prefix;
  q.AppendName( name, taxonomy );
  if( precision == Precision::Single )
    Append( name, HMod::sSinglePrec );
}

void ModAction::AddDependencies( HModList &ModList ) const
{
  std::string Gauge{ ModList.TakeOwnership( new ModGauge( ModList, tax, q.GaugeSmear, precision ) ) };
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
: HMod(taxonomy), q(q_)
{
  name = Prefix;
  q.AppendName( name, taxonomy );
}

template<typename TPar>
void ModSolver::LoadEigenPar( TPar &epPar, HModList &ModList, Precision XformPres ) const
{
  epPar.filestem = q.eigenPack;
  epPar.multiFile = q.multiFile;
  epPar.redBlack = q.redBlack;
  epPar.size = q.size;
  epPar.Ls = q.Ls;
  if( tax.GaugeFixed() )
  {
    epPar.gaugeXform = ModList.TakeOwnership( new ModGaugeXform( ModList, tax, q.GaugeSmear, XformPres ) );
  }
}

template<typename TEPLoad>
std::string ModSolver::LoadEigenPack( HModList &ModList, Precision XformPres ) const
{
  std::string EigenPackName{ "epack_" + name };
  typename TEPLoad::Par epPar;
  LoadEigenPar( epPar, ModList, XformPres );
  ModList.application.createModule<TEPLoad>( EigenPackName, epPar );
  return EigenPackName;
}

#ifdef MLU_HADRONS_HAS_GUESSERS
template<typename TGuesser>
void ModSolver::BatchGuessLoad( HModList &ModList, const std::string &GuesserName ) const
{
  typename TGuesser::Par guessPar;
  LoadEigenPar( guessPar.eigenPack, ModList, Precision::Double );
  guessPar.evBatchSize = ModList.params.Run.evBatchSize;
  guessPar.sourceBatchSize = ModList.params.Run.sourceBatchSize;
  ModList.application.createModule<TGuesser>( GuesserName, guessPar );
}

template<typename TGuesser>
void ModSolver::BatchGuessPreload( HModList &ModList, const std::string &GuesserName,
                                  const std::string &epName ) const
{
  typename TGuesser::Par guessPar;
  guessPar.eigenPack = epName;
  guessPar.epSize = q.size;
  guessPar.evBatchSize = ModList.params.Run.evBatchSize;
  guessPar.sourceBatchSize = ModList.params.Run.sourceBatchSize;
  ModList.application.createModule<TGuesser>( GuesserName, guessPar );
}
#endif

void ModSolver::AddDependencies( HModList &ModList ) const
{
  std::string EigenGuessName;
  if( q.eigenPack.length() )
  {
#ifdef MLU_HADRONS_HAS_GUESSERS
    EigenGuessName = "guesser_" + name;
    if( ModList.params.Run.BatchSize )
    {
      if( ModList.params.Run.PreLoadEigen )
      {
        // Double-precision batch guesser together with pre-loaded eigen pack
        std::string epName;
        if( q.eigenSinglePrecision )
        {
          epName = LoadEigenPack<MIO::LoadFermionEigenPackF>( ModList, Precision::Single );
          BatchGuessPreload<MGuesser::BatchExactDeflationEPackF>( ModList, EigenGuessName, epName );
        }
        else
        {
          epName = LoadEigenPack<MIO::LoadFermionEigenPack>( ModList, Precision::Double );
          BatchGuessPreload<MGuesser::BatchExactDeflation>( ModList, EigenGuessName, epName );
        }
      }
      else
      {
        // Double-precision batch guesser, loading from single or double-precision
        if( q.eigenSinglePrecision )
          BatchGuessLoad<MGuesser::BatchExactDeflationLoadDIoF>( ModList, EigenGuessName );
        else
          BatchGuessLoad<MGuesser::BatchExactDeflationLoad>( ModList, EigenGuessName );
      }
    }
    else
    {
      // Batch size == 0 -> Exact deflation. Always double-precision
      typename MGuesser::ExactDeflation::Par guessPar;
      if( q.eigenSinglePrecision )
        guessPar.eigenPack = LoadEigenPack<MIO::LoadFermionEigenPackIo32>( ModList, Precision::Double );
      else
        guessPar.eigenPack = LoadEigenPack<MIO::LoadFermionEigenPack>    ( ModList, Precision::Double );
      guessPar.size = q.size;
      ModList.application.createModule<MGuesser::ExactDeflation>( EigenGuessName, guessPar );
    }
#else
    if( !q.eigenSinglePrecision ) // Double-precision Eigen packs can be used anywhere
      EigenGuessName = LoadEigenPack<MIO::LoadFermionEigenPack>( ModList, Precision::Double );
    else if( q.MixedPrecision() ) // Single-precision Eigen packs for MP solver
      EigenGuessName = LoadEigenPack<MIO::LoadFermionEigenPackF>( ModList, Precision::Single );
    else // Load single-precision as double-
      EigenGuessName = LoadEigenPack<MIO::LoadFermionEigenPackIo32>( ModList, Precision::Double );
#endif
  }
  if( q.MixedPrecision() )
  {
    using T = MSolver::MixedPrecisionRBPrecCG;
    typename T::Par solverPar;
#ifdef MLU_HADRONS_HAS_GUESSERS
    solverPar.outerGuesser      = EigenGuessName;
#else
    solverPar.eigenPack         = EigenGuessName;
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
    solverPar.guesser      = EigenGuessName;
#else
    solverPar.eigenPack    = EigenGuessName;
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

ModProp::ModProp( HModList &ModList, const Taxonomy &taxonomy, const Quark &q_, const MLU::Momentum &p_,
                  int t_, int hit_ )
: HMod(taxonomy), q{q_}, p{p_}, t{t_}, hit{hit_}
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

ModSlicedProp::ModSlicedProp( HModList &ML, const Taxonomy &tax_, const Quark &q_, const MLU::Momentum &p_,
                              int t_, int hit_ )
: HMod(tax_), q{q_}, p{p_}, t{t_}, hit{hit_}
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
  smearPar.q = ModList.MakeProp( ModList, tax, q, p, t, hit );
  smearPar.sink = ModList.TakeOwnership( new ModSinkSmear( ModList, tax, -p ) );
  ModList.application.createModule<MSink::Smear>(name, smearPar);
}

/**************************
 Sequential Source
**************************/

const std::string ModSourceSeq::Prefix{ "Seq" + ModSource::Prefix };

ModSourceSeq::ModSourceSeq( HModList &ML, const Taxonomy &tax_, Gamma::Algebra current_, int deltaT_,
                            const MLU::Momentum &pSeq_, const Quark &q_, const MLU::Momentum &p_, int t_ )
: HMod(tax_), Current{current_}, deltaT{deltaT_}, pSeq{pSeq_}, q{q_}, p{p_}, t{t_}
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
  seqPar.q = ModList.MakeProp( ModList, tax, q, p, t );
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
      seqPar.LowerLeft  = ModList.params.Run.SpatialPos + MLU::Space + std::to_string( tSink );
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
                        const MLU::Momentum &pSeq_, const Quark &q_, const MLU::Momentum &p_, int t_ )
: HMod(tax_),qSeq{qSeq_},Current{current_},deltaT{deltaT_},pSeq{pSeq_},q{q_},p{p_},t{t_}
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
                                const Quark &q1_, const Quark &q2_, const MLU::Momentum &p1_, int t_,
                                const MLU::Momentum &p2_, int hit_ )
// Two-point Region-sinks don't really exist - silently translate to point-sinks
: HMod(taxonomy == Species::Region ? Taxonomy(taxonomy.family, Species::Point) : taxonomy),
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
    mesPar.q1 = ModList.MakeProp( ModList, tax, q1, p1, t, hit );
    mesPar.q2 = ModList.MakeProp( ModList, tax, q2, p2, t, hit );
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
                                const Quark &Spec_, const MLU::Momentum &pSeq_, const MLU::Momentum &p_,
                                Gamma::Algebra Cur_, int dT_, int t_, bool bHA_ )
: HMod(tax_), bReverse{tax == Family::GR ? !bRev_ : bRev_}, qSnk{Snk_}, qSrc{Src_}, qSpec{Spec_},
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
  const MLU::Momentum &pSource( bReverse ? pSeq : p );
  const MLU::Momentum &pSink  ( bReverse ? p : pSeq );
  MLU::Momentum pCurrent( - ( pSource + pSink ) );
  if( bInvertSeq )
  {
    par.q1 = ML.MakeProp   ( ML, stax, sqSrc, pSource, st );
    par.q2 = ML.MakePropSeq( ML, stax, sqSnk, Current, sdeltaT, -pSink, qSpec, p0, st );
  }
  else
  {
    par.q1 = ML.MakePropSeq( ML, stax, sqSnk, Current, sdeltaT,  pSink, qSpec, p0, st );
    par.q2 = ML.MakeProp   ( ML, stax, sqSrc,-pSource, st );
  }
  par.sink = ML.TakeOwnership( new ModSink( ML, stax, pCurrent ) );
  ML.application.createModule<MContraction::Meson>(name, par);
}
