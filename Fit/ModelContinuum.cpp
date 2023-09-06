/*************************************************************************************
 
 Chiral continuum fit

 Source file: ModelContinuum.cpp
 
 Copyright (C) 2023
 
 Author: Michael Marshall <Mike@lqcd.me>
 
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

// Keep my models separate from my fitting code

#include "ModelContinuum.hpp"

const std::string sFF{ "ff" };

// This is the glue between abstract models and specific implementations
// Other than here, these must be kept strictly separate

static const std::string sModelTypeUnknown{ "Unknown" };
static const std::string sModelTypeContinuum{ "Continuum" };

std::ostream & operator<<( std::ostream &os, const eModelType m )
{
  switch( m )
  {
    case eModelType::Continuum:
      os << sModelTypeContinuum;
      break;
    default:
      os << sModelTypeUnknown;
      break;
  }
  return os;
}

std::ostream & operator<<( std::ostream &os, const ModelType m )
{
  return os << static_cast<eModelType>( m.t );
}

std::istream & operator>>( std::istream &is, eModelType &m )
{
  std::string s;
  if( is >> s )
  {
    if( Common::EqualIgnoreCase( s, sModelTypeContinuum ) )
      m = eModelType::Continuum;
    else
    {
      m = eModelType::Unknown;
      is.setstate( std::ios_base::failbit );
    }
  }
  return is;
}

std::istream & operator>>( std::istream &is, ModelType &m )
{
  eModelType t;
  is >> t;
  m.t = static_cast<int>( t );
  return is;
}

// Create a model of the appropriate type - this is the only place with knowledge of this mapping
ModelPtr Model::MakeModel( const Model::CreateParams &cp, Model::Args &Args )
{
  // Now work out what type of model we are creating
  const eModelType modelType{ eModelType::Continuum };
  const Common::FormFactor ff{ Common::FromString<Common::FormFactor>( Args.Remove( sFF ) ) };
  // Make the model
  std::cout << "  " << modelType;
  for( typename Model::Args::value_type v : Args )
  {
    std::cout << Common::CommaSpace << v.first;
    if( !v.second.empty() )
      std::cout << '=' << v.second;
  }
  ModelPtr model;
  try
  {
    model.reset( new ModelContinuum( cp, Args, ff ) );
  }
  catch(...)
  {
    std::cout << Common::Space;
    throw;
  }
  std::cout << " => " << model->Description() << Common::NewLine;
  // 2 part construction - now that virtual functions in place
  return model;
}

const EnsembleInfo &ModelContinuum::GetEnsembleInfo() const
{
  if( Ensemble.empty() )
    throw std::runtime_error( "ModelContinuum(): Ensemble unavailable - specify manually" );
  typename EnsembleMapT::const_iterator cit{ Parent.EnsembleMap.find( Ensemble ) };
  if( cit == Parent.EnsembleMap.cend() )
    throw std::runtime_error( "ModelContinuum() Ensemble not supported " + Ensemble );
  return cit->second;
}

const std::string ModelContinuum::FieldQSq{ "qSq" };
const std::string ModelContinuum::FieldEL{ "EL" };

ModelContinuum::ModelContinuum( const Model::CreateParams &mcp, Model::Args &Args,
                                Common::FormFactor ff_ )
: Model( mcp, 1 ),
  ff{ff_},
  Ensemble{ Args.Remove( "Ensemble", mcp.pCorr->Ensemble ) },
  ei{ GetEnsembleInfo() },
  Parent{ dynamic_cast<const ::CreateParams &>( mcp ).Parent },
  mf{ * dynamic_cast<const ModelFile *>( mcp.pCorr ) },
  idxFitColumn{ mf.params.Find( Common::Param::Key( Common::GetFormFactorString( ff ) ),
                                "ModelContinuum::GetFitColumn()" )->second() },
  idxFF{ Parent.ffIndex( ff ) },
  p{ mcp.pCorr->Name_.MesonP[1].p }
{
  // Check source and sink names match the parent
  const Common::FileNameAtt &fna{ mcp.pCorr->Name_ };
  if( fna.BaseShortParts.size() < 3 || fna.Spectator.empty() )
    throw std::runtime_error( "Meson names not included in model filename " + fna.Filename );
  std::string MesonSnk{ Common::MesonName( fna.BaseShortParts[1], fna.Spectator ) };
  std::string MesonSrc{ Common::MesonName( fna.BaseShortParts[2], fna.Spectator ) };;
  if( !Common::EqualIgnoreCase( MesonSnk, Parent.Meson[idxSnk] ) )
    throw std::runtime_error( "ModelContinuum(): Snk " + MesonSnk
                             + " doesn't match parent " + Parent.Meson[idxSnk] );
  if( !Common::EqualIgnoreCase( MesonSrc, Parent.Meson[idxSrc] ) )
    throw std::runtime_error( "ModelContinuum(): Src " + MesonSrc
                             + " doesn't match parent " + Parent.Meson[idxSrc] );
  const ::CreateParams &cp{ dynamic_cast<const ::CreateParams &>( mcp ) };
  if( cp.pCorr->Name_.MesonP[0].p )
  {
    std::ostringstream os;
    os << "ModelContinuum(): Src momentum " << cp.pCorr->Name_.MesonP[0].p << " != 0";
    throw std::runtime_error( os.str().c_str() );
  }
  // General parameters
  // X Vector
  EL.Key.Object = { Ensemble, std::to_string( p.p2() ) };
  EL.Key.Name = FieldEL;
  mH.Key.Object = { Ensemble };
  mH.Key.Name = "mH";
  mL.Key.Object = mH.Key.Object;
  mL.Key.Name = "mL";
  qSq.Key.Object = EL.Key.Object;
  qSq.Key.Name = FieldQSq;
}

void ModelContinuum::DefineXVector( DataSet &ds, int i )
{
  aEL.Key.Object = { Ensemble, std::to_string( p.p2() ) };
  aEL.Key.Name = "a" + FieldEL;
  ds.AddConstant( aEL.Key, i, Param::Key{EL.Key.Name} );
  amH.Key.Object = { Ensemble };
  amH.Key.Name = "amH";
  ds.AddConstant( amH.Key, i, Param::Key{mH.Key.Name} );
  amL.Key.Object = mH.Key.Object;
  amL.Key.Name = "amL";
  ds.AddConstant( amL.Key, i, Param::Key{mL.Key.Name} );
  aqSq.Key.Object = aEL.Key.Object;
  aqSq.Key.Name = "a" + FieldQSq;
  ds.AddConstant( aqSq.Key, i, Param::Key{qSq.Key.Name} );
}

const std::string &ModelContinuum::XVectorKeyName() const
{
  static const std::string MyKey{ "nSq" };
  return MyKey;
}

std::string ModelContinuum::XVectorKey() const
{
  return std::to_string( p.p2() );
}

std::vector<Param::Key> ModelContinuum::XVectorKeyNames() const
{
  std::vector<Param::Key> vk;
  vk.reserve( 20 );
  vk.push_back( aEL.Key );
  vk.push_back( amH.Key );
  vk.push_back( amL.Key );
  vk.push_back( aqSq.Key );
  vk.push_back( EL.Key );
  vk.push_back( mH.Key );
  vk.push_back( mL.Key );
  vk.push_back( qSq.Key );
  // Do I need these? practical effect is to write these fields to the summary
  vk.push_back( Parent.kaInv[ei.idx] );
  vk.push_back( Parent.kPDGH );
  vk.push_back( Parent.kPDGL );
  vk.push_back( Parent.kfPi );
  vk.push_back( Parent.kmPDGPi );
  vk.push_back( Parent.kmPi[ei.idx] );
  vk.push_back( Parent.kDelta[idxFF] );
  if( Parent.FVEnabled[idxFF] )
  {
    vk.push_back( Parent.kFVSim[ei.idx] );
    vk.push_back( Parent.kFVPhys[ei.idx] );
  }
  if( Parent.ChiEnabled[idxFF] )
  {
    vk.push_back( Parent.kChiSim[ei.idx] );
    vk.push_back( Parent.kChiPhys );
  }
  if( Parent.FVEnabled[idxFF] || Parent.ChiEnabled[idxFF] )
    vk.push_back( Parent.kChiFV[ei.idx] );
  return vk;
}

void ModelContinuum::AddParameters( Params &mp )
{
  AddParam( mp, aEL );
  AddParam( mp, amH );
  AddParam( mp, amL );
  AddParam( mp, aqSq );
  AddParam( mp, EL, 1, false, Common::Param::Type::Derived );
  AddParam( mp, mH, 1, false, Common::Param::Type::Derived );
  AddParam( mp, mL, 1, false, Common::Param::Type::Derived );
  AddParam( mp, qSq, 1, false, Common::Param::Type::Derived );
}

std::size_t ModelContinuum::GetFitColumn() const
{
  return idxFitColumn;
}

void ModelContinuum::SaveParameters( const Params &mp )
{
  aEL.idx = mp.at( aEL.Key )();
  amH.idx = mp.at( amH.Key )();
  amL.idx = mp.at( amL.Key )();
  aqSq.idx = mp.at( aqSq.Key )();
  EL.idx = mp.at( EL.Key )();
  mH.idx = mp.at( mH.Key )();
  mL.idx = mp.at( mL.Key )();
  qSq.idx = mp.at( qSq.Key )();
}

// Get a descriptive string for the model
std::string ModelContinuum::Description() const
{
  std::ostringstream ss;
  ss << sModelTypeContinuum << '(' << ff << ',' << Ensemble << ',' << p << ',' << ei.aInv_L << ')';
  return ss.str().c_str();
}

void ModelContinuum::Guessable( ParamsPairs &PP ) const
{
  if( Parent.c0Enabled[idxFF] )
    PP.SetState( ParamsPairs::State::Known, Parent.kC0[idxFF], 1 );
  if( Parent.c1Enabled[idxFF] )
    PP.SetState( ParamsPairs::State::Known, Parent.kC1[idxFF], 1 );
  for( unsigned int i = 0; i < Parent.dEnabled[idxFF].size(); ++i )
    if( Parent.dEnabled[idxFF][i] )
      PP.SetState( ParamsPairs::State::Known, Parent.kD[idxFF][i], 1 );
  for( unsigned int i = 0; i < Parent.eEnabled[idxFF].size(); ++i )
    if( Parent.eEnabled[idxFF][i] )
      PP.SetState( ParamsPairs::State::Known, Parent.kE[idxFF][i], 1 );
}

std::size_t ModelContinuum::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                   const VectorView &FitData, std::vector<int> FitTimes,
                   bool bLastChance ) const
{
  if( Parent.c0Enabled[idxFF] && !bKnown[Parent.idxC0[idxFF]] )
  {
    Guess[Parent.idxC0[idxFF]] = 1;
    bKnown[Parent.idxC0[idxFF]] = true;
  }
  if( Parent.c1Enabled[idxFF] && !bKnown[Parent.idxC1[idxFF]] )
  {
    Guess[Parent.idxC1[idxFF]] = 1;
    bKnown[Parent.idxC1[idxFF]] = true;
  }
  for( unsigned int i = 0; i < Parent.dEnabled[idxFF].size(); ++i )
    if( Parent.dEnabled[idxFF][i] && !bKnown[Parent.idxD[idxFF][i]] )
    {
      Guess[Parent.idxD[idxFF][i]] = 1;
      bKnown[Parent.idxD[idxFF][i]] = true;
    }
  for( unsigned int i = 0; i < Parent.eEnabled[idxFF].size(); ++i )
    if( Parent.eEnabled[idxFF][i] && !bKnown[Parent.idxE[idxFF][i]] )
    {
      Guess[Parent.idxE[idxFF][i]] = 1;
      bKnown[Parent.idxE[idxFF][i]] = true;
    }
  return 0;
}

ModelType ModelContinuum::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::Continuum );
  return m;
}

scalar ModelContinuum::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  // Compute the chiral term
  // Chiral log term
  const scalar sChiFV{ Parent.ChiEnabled[idxFF] || Parent.FVEnabled[idxFF]
                       ? ModelParams[Parent.idxChiFV[ei.idx]] : 0 };
  const scalar TermChiral = Parent.c0Enabled[idxFF]
                            ? ModelParams[Parent.idxC0[idxFF]] * ( 1 + sChiFV ) : 0;
  // Compute the MPi term
  const scalar TermMPi{ Parent.c1Enabled[idxFF]
    ? ModelParams[Parent.idxC1[idxFF]] * ModelParams[Parent.idxDeltaMPi[ei.idx]] : 0 };
  // Compute E_L / Lambda
  const scalar ELOnLambda{ ModelParams[EL.idx] * LambdaInv };
  scalar Factor = ELOnLambda;
  scalar TermEL = 0;
  for( unsigned int i = 0; i < Parent.eEnabled[idxFF].size(); ++i )
  {
    if( Parent.eEnabled[idxFF][i] )
      TermEL += ModelParams[Parent.idxE[idxFF][i]] * Factor;
    Factor *= ELOnLambda;
  }
  // Compute the Discretisation term
  const scalar aInv{ ModelParams[Parent.idxaInv[ei.idx]] };
  const scalar aLambda{ Lambda / aInv };
  const scalar aLambdaSq{ aLambda * aLambda };
  Factor = aLambdaSq;
  scalar TermDisc = 0;
  for( unsigned int i = 0; i < Parent.dEnabled[idxFF].size(); ++i )
  {
    if( Parent.dEnabled[idxFF][i] )
      TermDisc += ModelParams[Parent.idxD[idxFF][i]] * Factor;
    Factor *= aLambdaSq;
  }
  // Return model result
  const scalar PoleTerm{ Lambda / ( ModelParams[EL.idx] + ModelParams[Parent.idxDelta[idxFF]] ) };
  const scalar Result = PoleTerm * ( TermChiral + TermMPi + TermEL + TermDisc );
  return Result;
}

// Cache values based solely on the model parameters (to speed up computation)
void ModelContinuum::SetReplica( Vector &ScratchPad, Vector &ModelParams ) const
{
  // At zero momentum, check that the energy is the mass
  if( !p.p2() && ModelParams[amL.idx] != ModelParams[aEL.idx] )
    throw std::runtime_error( "ModelContinuum amL != aEL(0) on Ensemble " + Ensemble );
  // Scale the data in lattice units by a
  const scalar aInv{ ModelParams[Parent.idxaInv[ei.idx]] };
  ModelParams[EL.idx] = ModelParams[aEL.idx] * aInv;
  ModelParams[mH.idx] = ModelParams[amH.idx] * aInv;
  ModelParams[mL.idx] = ModelParams[amL.idx] * aInv;
  ModelParams[qSq.idx] = ModelParams[aqSq.idx] * aInv * aInv;
}
