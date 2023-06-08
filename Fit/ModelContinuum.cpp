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

ModelContinuum::ModelContinuum( const Model::CreateParams &cp, Model::Args &Args,
                                Common::FormFactor ff_ )
: Model( cp, 1 ),
  mf{ * dynamic_cast<const ModelFile *>( cp.pCorr ) },
  ff{ff_},
  p{ cp.pCorr->Name_.MesonP[1].p },
  Ensemble{ Args.Remove( "Ensemble", cp.pCorr->Ensemble ) }
{
  if( Ensemble.empty() )
    throw std::runtime_error( "ModelContinuum(): Ensemble unavailable - specify manually" );
  // Now remove source and sink names
  bool bManualObjectID[2];
  Meson[0] = Args.Remove( "Src", &bManualObjectID[0] );
  Meson[1] = Args.Remove( "Snk", &bManualObjectID[1] );
  if( !bManualObjectID[0] )
  {
    if( cp.pCorr->Name_.Meson.empty() )
      throw std::runtime_error( "ModelContinuum(): Src unavailable - specify manually" );
    Meson[0] = Common::MesonName( cp.pCorr->Name_.BaseShortParts[2], cp.pCorr->Name_.Spectator );
  }
  if( !bManualObjectID[1] )
  {
    if( cp.pCorr->Name_.Meson.size() < 2 )
      throw std::runtime_error( "ModelContinuum(): Snk unavailable - specify manually" );
    Meson[1] = Common::MesonName( cp.pCorr->Name_.BaseShortParts[1], cp.pCorr->Name_.Spectator );
  }
  if( Common::EqualIgnoreCase( Meson[0], Meson[1] ) )
    throw std::runtime_error( "ModelContinuum(): Snk and Src are both " + Meson[0] );
  if( cp.pCorr->Name_.MesonP[0].p )
  {
    std::ostringstream os;
    os << "ModelContinuum(): Src momentum " << cp.pCorr->Name_.MesonP[0].p << " != 0";
    throw std::runtime_error( os.str().c_str() );
  }
}

void ModelContinuum::DefineXVector( DataSet &ds, int i )
{
  EL.Key.Object = { Ensemble, std::to_string( p.p2() ) };
  EL.Key.Name = "EL";
  ds.AddConstant( EL.Key, i, Param::Key{EL.Key.Name} );
  kMu.Key.Object = EL.Key.Object;
  kMu.Key.Name = "kMu";
  ds.AddConstant( kMu.Key, i, Param::Key{kMu.Key.Name} );
  mH.Key.Object = { Ensemble };
  mH.Key.Name = "mH";
  ds.AddConstant( mH.Key, i, Param::Key{mH.Key.Name} );
  mL.Key.Object = mH.Key.Object;
  mL.Key.Name = "mL";
  ds.AddConstant( mL.Key, i, Param::Key{"mL"} );
  qSq.Key.Object = kMu.Key.Object;
  qSq.Key.Name = "qSq";
  ds.AddConstant( qSq.Key, i, Param::Key{qSq.Key.Name} );
  // For now I use a per-ensemble linear model
  Slope.Key.Object = { Ensemble };
  Slope.Key.Name = "Slope";
  YIntercept.Key.Object = Slope.Key.Object;
  YIntercept.Key.Name = "YIntercept";
}

void ModelContinuum::AddParameters( Params &mp )
{
  AddParam( mp, EL );
  AddParam( mp, kMu );
  AddParam( mp, mH );
  AddParam( mp, mL );
  AddParam( mp, qSq );
  AddParam( mp, Slope );
  AddParam( mp, YIntercept );
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
  EL.idx = mp.at( EL.Key )();
  kMu.idx = mp.at( kMu.Key )();
  mH.idx = mp.at( mH.Key )();
  mL.idx = mp.at( mL.Key )();
  qSq.idx = mp.at( qSq.Key )();
  Slope.idx = mp.at( Slope.Key )();
  YIntercept.idx = mp.at( YIntercept.Key )();
}

// Get a descriptive string for the model
std::string ModelContinuum::Description() const
{
  std::ostringstream ss;
  ss << sModelTypeContinuum << '(' << ff << ',' << Ensemble << ',' << p << ')';
  return ss.str().c_str();
}

void ModelContinuum::Guessable( ParamsPairs &PP ) const
{
  PP.SetState( ParamsPairs::State::Known, Slope.Key, Slope.param->size );
  PP.SetState( ParamsPairs::State::Known, YIntercept.Key, YIntercept.param->size );
}

std::size_t ModelContinuum::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                   const VectorView &FitData, std::vector<int> FitTimes,
                   bool bLastChance ) const
{
  if( !bKnown[Slope.idx] )
  {
    bKnown[Slope.idx] = true;
    Guess[Slope.idx] = 7;
  }
  if( !bKnown[YIntercept.idx] )
  {
    bKnown[YIntercept.idx] = true;
    Guess[YIntercept.idx] = 0.7;
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
  if( !p.p2() && ModelParams[mL.idx] != ModelParams[EL.idx] )
    throw std::runtime_error( "ModelContinuum mL != EL(0) on Ensemble " + Ensemble );
  return ModelParams[Slope.idx] * ModelParams[qSq.idx] + ModelParams[YIntercept.idx];
}
