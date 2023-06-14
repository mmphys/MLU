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
#include "Continuum.hpp"

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

ModelContinuum::ModelContinuum( const Model::CreateParams &mcp, Model::Args &Args,
                                Common::FormFactor ff_ )
: Model( mcp, 1 ),
  mf{ * dynamic_cast<const ModelFile *>( mcp.pCorr ) },
  ff{ff_},
  p{ mcp.pCorr->Name_.MesonP[1].p },
  Ensemble{ Args.Remove( "Ensemble", mcp.pCorr->Ensemble ) }
{
  const ::CreateParams &cp{ dynamic_cast<const ::CreateParams &>( mcp ) };
  if( Ensemble.empty() )
    throw std::runtime_error( "ModelContinuum(): Ensemble unavailable - specify manually" );
  typename EnsembleMapT::const_iterator cit{ cp.EnsembleMap.find( Ensemble ) };
  if( cit == cp.EnsembleMap.cend() )
    throw std::runtime_error( "ModelContinuum() Ensemble not supported " + Ensemble );
  const EnsembleInfo &ei{ cit->second };
  aInv_L = ei.aInv_L;
  // Now remove source and sink names
  bool bManualObjectID[2];
  Meson[0] = Args.Remove( "Src", &bManualObjectID[0] );
  Meson[1] = Args.Remove( "Snk", &bManualObjectID[1] );
  if( !bManualObjectID[0] )
  {
    if( cp.pCorr->Name_.Quark.empty() )
      throw std::runtime_error( "ModelContinuum(): Src unavailable - specify manually" );
    Meson[0] = Common::MesonName( cp.pCorr->Name_.Quark[0], cp.pCorr->Name_.Spectator );
  }
  if( !bManualObjectID[1] )
  {
    if( cp.pCorr->Name_.Quark.size() < 2 )
      throw std::runtime_error( "ModelContinuum(): Snk unavailable - specify manually" );
    Meson[1] = Common::MesonName( cp.pCorr->Name_.Quark[1], cp.pCorr->Name_.Spectator );
  }
  if( Common::EqualIgnoreCase( Meson[0], Meson[1] ) )
    throw std::runtime_error( "ModelContinuum(): Snk and Src are both " + Meson[0] );
  if( cp.pCorr->Name_.MesonP[0].p )
  {
    std::ostringstream os;
    os << "ModelContinuum(): Src momentum " << cp.pCorr->Name_.MesonP[0].p << " != 0";
    throw std::runtime_error( os.str().c_str() );
  }
  // General parameters
  aInv.Key.Object = { Ensemble };
  aInv.Key.Name = "aInv";
  fPi.Key.Name = "fPi";
  mPi.Key.Object = aInv.Key.Object;
  mPi.Key.Name = "mPi"; // Simulated pion mass on each ensemble
  mPDGPi.Key.Name = "PDGPi";
  // The excited meson in the pole depends on the charm decay final quark
  const bool bVector{ ff == Common::FormFactor::fplus || ff == Common::FormFactor::fperp };
  if( std::toupper( cp.pCorr->Name_.Quark[1][0] ) == 'L' )
    mPDGDStar.Key.Name = bVector ? "PDGDStar" : "PDGD0Star";
  else
    mPDGDStar.Key.Name = bVector ? "PDGDsStar" : "PDGDs0Star";
  mPDGH.Key.Name = std::toupper( cp.pCorr->Name_.Spectator[0] ) == 'S' ? "PDGDs" : "PDGD";
  mPDGL.Key.Name = std::toupper( Meson[1][0] ) == 'K' ? "PDGK" : "PDGP";
  Delta.Key.Name = "Delta";
  // X Vector
  EL.Key.Object = { Ensemble, std::to_string( p.p2() ) };
  EL.Key.Name = "EL";
  kMu.Key.Object = EL.Key.Object;
  kMu.Key.Name = "kMu";
  mH.Key.Object = { Ensemble };
  mH.Key.Name = "mH";
  mL.Key.Object = mH.Key.Object;
  mL.Key.Name = "mL";
  qSq.Key.Object = kMu.Key.Object;
  qSq.Key.Name = "qSq";
  // This is my model
  c0.Key.Object = { GetFormFactorString( ff ) };
  c0.Key.Name = "c0";
  c1.Key.Object = c0.Key.Object;
  c1.Key.Name = "c1";
  c2.Key.Object = c0.Key.Object;
  c2.Key.Name = "c2";
  c3.Key.Object = c0.Key.Object;
  c3.Key.Name = "c3";
  c4.Key.Object = c0.Key.Object;
  c4.Key.Name = "c4";
}

void ModelContinuum::DefineXVector( DataSet &ds, int i )
{
  aEL.Key.Object = { Ensemble, std::to_string( p.p2() ) };
  aEL.Key.Name = "aEL";
  ds.AddConstant( aEL.Key, i, Param::Key{EL.Key.Name} );
  akMu.Key.Object = aEL.Key.Object;
  akMu.Key.Name = "akMu";
  ds.AddConstant( akMu.Key, i, Param::Key{kMu.Key.Name} );
  amH.Key.Object = { Ensemble };
  amH.Key.Name = "amH";
  ds.AddConstant( amH.Key, i, Param::Key{mH.Key.Name} );
  amL.Key.Object = mH.Key.Object;
  amL.Key.Name = "amL";
  ds.AddConstant( amL.Key, i, Param::Key{mL.Key.Name} );
  aqSq.Key.Object = kMu.Key.Object;
  aqSq.Key.Name = "aqSq";
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
  vk.reserve( 18 );
  vk.push_back( aInv.Key );
  vk.push_back( fPi.Key );
  vk.push_back( mPi.Key );
  vk.push_back( mPDGPi.Key );
  vk.push_back( aEL.Key );
  vk.push_back( akMu.Key );
  vk.push_back( amH.Key );
  vk.push_back( amL.Key );
  vk.push_back( aqSq.Key );
  vk.push_back( EL.Key );
  vk.push_back( kMu.Key );
  vk.push_back( mH.Key );
  vk.push_back( mL.Key );
  vk.push_back( qSq.Key );
  vk.push_back( mPDGDStar.Key );
  vk.push_back( mPDGH.Key );
  vk.push_back( mPDGL.Key );
  vk.push_back( Delta.Key );
  return vk;
}

void ModelContinuum::AddParameters( Params &mp )
{
  AddParam( mp, aInv );
  AddParam( mp, fPi );
  AddParam( mp, mPi );
  AddParam( mp, mPDGPi );
  AddParam( mp, mPDGDStar );
  AddParam( mp, mPDGH );
  AddParam( mp, mPDGL );
  AddParam( mp, Delta, 1, false, Common::Param::Type::Derived );
  AddParam( mp, aEL );
  AddParam( mp, akMu );
  AddParam( mp, amH );
  AddParam( mp, amL );
  AddParam( mp, aqSq );
  AddParam( mp, EL, 1, false, Common::Param::Type::Derived );
  AddParam( mp, kMu, 1, false, Common::Param::Type::Derived );
  AddParam( mp, mH, 1, false, Common::Param::Type::Derived );
  AddParam( mp, mL, 1, false, Common::Param::Type::Derived );
  AddParam( mp, qSq, 1, false, Common::Param::Type::Derived );
  AddParam( mp, c0 );
  AddParam( mp, c1 );
  AddParam( mp, c2 );
  AddParam( mp, c3 );
  AddParam( mp, c4 );
}

std::size_t ModelContinuum::GetFitColumn() const
{
  Common::Param::Key k;
  k.Name = Common::GetFormFactorString( ff );
  Params::const_iterator it{ mf.params.Find( k, "ModelContinuum::GetFitColumn()" ) };
  return it->second();
}

void ModelContinuum::SaveParameters( const Params &mp )
{
  aInv.idx = mp.at( aInv.Key )();
  fPi.idx = mp.at( fPi.Key )();
  mPi.idx = mp.at( mPi.Key )();
  mPDGPi.idx = mp.at( mPDGPi.Key )();
  mPDGDStar.idx = mp.at( mPDGDStar.Key )();
  mPDGH.idx = mp.at( mPDGH.Key )();
  mPDGL.idx = mp.at( mPDGL.Key )();
  Delta.idx = mp.at( Delta.Key )();
  aEL.idx = mp.at( aEL.Key )();
  akMu.idx = mp.at( akMu.Key )();
  amH.idx = mp.at( amH.Key )();
  amL.idx = mp.at( amL.Key )();
  aqSq.idx = mp.at( aqSq.Key )();
  EL.idx = mp.at( EL.Key )();
  kMu.idx = mp.at( kMu.Key )();
  mH.idx = mp.at( mH.Key )();
  mL.idx = mp.at( mL.Key )();
  qSq.idx = mp.at( qSq.Key )();
  c0.idx = mp.at( c0.Key )();
  c1.idx = mp.at( c1.Key )();
  c2.idx = mp.at( c2.Key )();
  c3.idx = mp.at( c3.Key )();
  c4.idx = mp.at( c4.Key )();
}

// Get a descriptive string for the model
std::string ModelContinuum::Description() const
{
  std::ostringstream ss;
  ss << sModelTypeContinuum << '(' << ff << ',' << Ensemble << ',' << p << ',' << aInv_L << ')';
  return ss.str().c_str();
}

void ModelContinuum::Guessable( ParamsPairs &PP ) const
{
  PP.SetState( ParamsPairs::State::Known, c0.Key, c0.param->size );
  PP.SetState( ParamsPairs::State::Known, c1.Key, c1.param->size );
  PP.SetState( ParamsPairs::State::Known, c2.Key, c2.param->size );
  PP.SetState( ParamsPairs::State::Known, c3.Key, c3.param->size );
  PP.SetState( ParamsPairs::State::Known, c4.Key, c4.param->size );
}

std::size_t ModelContinuum::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                   const VectorView &FitData, std::vector<int> FitTimes,
                   bool bLastChance ) const
{
  if( !bKnown[c0.idx] )
  {
    bKnown[c0.idx] = true;
    Guess[c0.idx] = 1;
  }
  if( !bKnown[c1.idx] )
  {
    bKnown[c1.idx] = true;
    Guess[c1.idx] = 1;
  }
  if( !bKnown[c2.idx] )
  {
    bKnown[c2.idx] = true;
    Guess[c2.idx] = 1;
  }
  if( !bKnown[c3.idx] )
  {
    bKnown[c3.idx] = true;
    Guess[c3.idx] = 1;
  }
  if( !bKnown[c4.idx] )
  {
    bKnown[c4.idx] = true;
    Guess[c4.idx] = 1;
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
  static constexpr scalar FourPi{ 4. * M_PI };
  // At zero momentum, check that the energy is the mass
  if( !p.p2() && ModelParams[amL.idx] != ModelParams[aEL.idx] )
    throw std::runtime_error( "ModelContinuum amL != aEL(0) on Ensemble " + Ensemble );
  // Scale the data in lattice units by a
  const scalar &_aInv{ ModelParams[aInv.idx] };
  ModelParams[EL.idx] = ModelParams[aEL.idx] * _aInv;
  ModelParams[kMu.idx] = ModelParams[akMu.idx] * _aInv;
  ModelParams[mH.idx] = ModelParams[amH.idx] * _aInv;
  ModelParams[mL.idx] = ModelParams[amL.idx] * _aInv;
  ModelParams[qSq.idx] = ModelParams[aqSq.idx] * _aInv * _aInv;
  // Compute the pole location and pole term
  ModelParams[Delta.idx] = 0.5 * ( ( ModelParams[mPDGDStar.idx] * ModelParams[mPDGDStar.idx]
                                    - ModelParams[mPDGL.idx] * ModelParams[mPDGL.idx] )
                                  / ModelParams[mPDGH.idx] - ModelParams[mPDGH.idx] );
  const scalar PoleTerm{ Lambda / ( ModelParams[EL.idx] + ModelParams[Delta.idx] ) };
  // Compute the c0 term
  scalar c0Num = DeltaF( ModelParams[mPi.idx] ) - DeltaF( ModelParams[mPDGPi.idx] );
  scalar c0Denom = FourPi * ModelParams[fPi.idx];
  c0Denom *= c0Denom;
  scalar c0Term = ModelParams[c0.idx] * ( 1. + c0Num / c0Denom );
  // Compute the c1 term
  scalar DeltaMPiSq = ModelParams[mPi.idx] * ModelParams[mPi.idx];
  DeltaMPiSq -= ModelParams[mPDGPi.idx] * ModelParams[mPDGPi.idx];
  scalar c1Term = ModelParams[c1.idx] * DeltaMPiSq * LambdaInv * LambdaInv;
  // Compute E_L / Lambda
  scalar ELOnLambda = ModelParams[EL.idx] * LambdaInv;
  // Compute the c2 term
  scalar c2Term = ModelParams[c2.idx] * ELOnLambda;
  // Compute the c3 term
  scalar c3Term = ModelParams[c3.idx] * ELOnLambda * ELOnLambda;
  // Compute the c4 term
  scalar aLambda = Lambda / _aInv;
  scalar c4Term = ModelParams[c4.idx] * aLambda * aLambda;
  // Return model result
  scalar Result = PoleTerm * ( c0Term + c1Term + c2Term + c3Term + c4Term );
  return Result;
}

scalar ModelContinuum::DeltaF( scalar M ) const
{
  const scalar MOnLambda{ M * LambdaInv };
  return -0.75 * M * M * std::log( MOnLambda * MOnLambda );
}
