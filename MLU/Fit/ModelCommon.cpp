/*************************************************************************************
 
 Common components for model implementations

 Source file: ModelCommon.cpp
 
 Copyright (C) 2019-2022
 
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

#include <MLUconfig.h>
#include "ModelCommon.hpp"

// This is the glue between abstract models and specific implementations
// Other than here, these must be kept strictly separate

#include "Model2pt.hpp"
#include "Model3pt.hpp"
#include "ModelConstant.hpp"
#include "ModelRatio.hpp"
#include "ModelZV.hpp"
#include "ModelMisc.hpp"

static const std::string sModelTypeUnknown{ "Unknown" };
static const std::string sModelTypeExp{ "Exp" };
static const std::string sModelTypeCosh{ "Cosh" };
static const std::string sModelTypeSinh{ "Sinh" };
static const std::string sModelTypeThreePoint{ "3pt" };
static const std::string sModelTypeConstant{ "Const" };
static const std::string sModelTypeRatio{ "R3" };
static const std::string sModelTypeZV{ "ZV" };
static const std::string sModelTypeMisc{ "Misc" };

std::ostream & operator<<( std::ostream &os, const eModelType m )
{
  switch( m )
  {
    case eModelType::Exp:
      os << sModelTypeExp;
      break;
    case eModelType::Cosh:
      os << sModelTypeCosh;
      break;
    case eModelType::Sinh:
      os << sModelTypeSinh;
      break;
    case eModelType::ThreePoint:
      os << sModelTypeThreePoint;
      break;
    case eModelType::Constant:
      os << sModelTypeConstant;
      break;
    case eModelType::R3:
      os << sModelTypeRatio;
      break;
    case eModelType::ZV:
      os << sModelTypeZV;
      break;
    case eModelType::Misc:
      os << sModelTypeMisc;
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
    if( MLU::EqualIgnoreCase( s, sModelTypeExp ) )
      m = eModelType::Exp;
    else if( MLU::EqualIgnoreCase( s, sModelTypeCosh ) )
      m = eModelType::Cosh;
    else if( MLU::EqualIgnoreCase( s, sModelTypeSinh ) )
      m = eModelType::Sinh;
    else if( MLU::EqualIgnoreCase( s, sModelTypeThreePoint ) )
      m = eModelType::ThreePoint;
    else if( MLU::EqualIgnoreCase( s, sModelTypeConstant ) )
      m = eModelType::Constant;
    else if( MLU::EqualIgnoreCase( s, sModelTypeRatio ) )
      m = eModelType::R3;
    else if( MLU::EqualIgnoreCase( s, sModelTypeZV ) )
      m = eModelType::ZV;
    else if( MLU::EqualIgnoreCase( s, sModelTypeMisc ) )
      m = eModelType::Misc;
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

MultiFitCreateParams::MultiFitCreateParams( const std::vector<std::string> &OpNames_,
                                            const MLU::CommandLine &cl_ )
: CreateParams( OpNames_, cl_),
  NumExponents{ cl.SwitchValue<int>( "e" ) },
  N{ cl.SwitchValue<int>( "N" ) },
  dispType{ cl.GotSwitch( "dispersion" ) ? cl.SwitchValue<MLU::DispersionType>( "dispersion" )
                                         : MLU::DispersionType::LatFreeScalar }
{
  if( N < 0 )
    throw std::runtime_error( "L/a " + std::to_string( N ) + " < 0");
}

std::string MultiFitCreateParams::Description() const
{
  std::string s{ Model::CreateParams::Description() };
  if( N )
  {
    if( !s.empty() )
      s.append( MLU::CommaSpace );
    s.append( "using dispersion relation N=" );
    s.append( std::to_string( N ) );
    s.append( " with " );
    s.append( MLU::GetDispersionString( dispType ) );
  }
  return s;
}

// Create a model of the appropriate type - this is the only place with knowledge of this mapping
ModelPtr Model::MakeModel( int ModelNum, const Model::CreateParams &mcp, Model::Args &Args )
{
  const MultiFitCreateParams &cp{ dynamic_cast<const MultiFitCreateParams &>( mcp ) };
  // Now work out what type of model we are creating
  eModelType modelType{ eModelType::Unknown };
  bool bGotModelType{ false };
  {
    std::string s{ Args.Remove( "Model", &bGotModelType ) };
    if( bGotModelType )
    {
      std::istringstream ss( s );
      if( !( ss >> modelType ) || !MLU::StreamEmpty( ss ) )
        throw std::runtime_error( "Unknown ModelType \"" + s + "\"" );
    }
  }
  if( !bGotModelType )
  {
    // We haven't been told which model to use. Choose a suitable default
    const std::size_t NumParts{ cp.pCorr->Name_.BaseShortParts.size() };
    const bool bGotDeltaT{ cp.pCorr->Name_.GotDeltaT() };
    const bool b3pt{ bGotDeltaT && cp.pCorr->Name_.Gamma.size() == 1 };
    if( b3pt )
    {
      if( NumParts && MLU::EqualIgnoreCase( cp.pCorr->Name_.BaseShortParts[0], "R3" ) )
        modelType = eModelType::R3;
      else
        modelType = eModelType::ThreePoint;
    }
    else if( bGotDeltaT && MLU::EqualIgnoreCase( cp.pCorr->Name_.BaseShortParts[0], "ZV" )
            && NumParts >= 2 && !cp.pCorr->Name_.HasNonZeroMomentum() )
    {
      modelType = eModelType::ZV;
    }
    else switch( dynamic_cast<const Fold *>( cp.pCorr )->parity )
    {
      case MLU::Parity::Even:
        modelType = eModelType::Cosh;
        break;
      case MLU::Parity::Odd:
        modelType = eModelType::Sinh;
        break;
      default:
        modelType = eModelType::Exp;
        break;
    }
  }
  // Say which model we're making and what the parameters are
  std::cout << std::setw( 5 ) << ModelNum << MLU::Space << modelType;
  for( typename Model::Args::value_type v : Args )
  {
    std::cout << MLU::CommaSpace << v.first;
    if( !v.second.empty() )
      std::cout << '=' << v.second;
  }
  // Make the model
  const int nExp{ Args.Remove( "e", cp.NumExponents ) };
  ModelPtr model;
  try
  {
    switch( modelType )
    {
      case eModelType::Exp:
        model.reset( new ModelExp( cp, Args, nExp ) );
        break;
      case eModelType::Cosh:
        model.reset( new ModelCosh( cp, Args, nExp ) );
        break;
      case eModelType::Sinh:
        model.reset( new ModelSinh( cp, Args, nExp ) );
        break;
      case eModelType::ThreePoint:
        model.reset( new Model3pt( cp, Args, nExp ) );
        break;
      case eModelType::Constant:
        model.reset( new ModelConstant( cp, Args, nExp ) );
        break;
      case eModelType::R3:
        model.reset( new ModelRatio( cp, Args, nExp ) );
        break;
      case eModelType::ZV:
        model.reset( new ModelZV( cp, Args, nExp ) );
        break;
      case eModelType::Misc:
        model.reset( new ModelMisc( cp, Args, nExp ) );
        break;
      default:
        throw std::runtime_error( "Unrecognised ModelType for " + cp.pCorr->Name_.Filename );
    }
  }
  catch(...)
  {
    std::cout << MLU::Space;
    throw;
  }
  std::cout << " => " << model->Description() << MLU::NewLine;
  // 2 part construction - now that virtual functions in place
  return model;
}

std::vector<std::string> ModelOverlap::GetOpNames( const Model::CreateParams &cp, Model::Args &Args )
{
  bool bForceSrcSnkDifferent;
  Args.Remove( "srcsnk", &bForceSrcSnkDifferent );
  std::vector<std::string> v( 2 );
  for( int i = 0; i < v.size(); ++i )
  {
    v[i] = Args.Remove( pSrcSnk[i], cp.GetOpName( i ) );
    if( bForceSrcSnkDifferent )
      v[i].append( pSrcSnk[i] );
  }
  return v;
}

// 2nd phase of construction (now that virtual functions in place)
ModelOverlap::ModelOverlap( const Model::CreateParams &cp, Model::Args &Args,
                            int NumExponents_, std::size_t NumOverlapExp_,
                            std::vector<std::string> &&OID, std::vector<std::string> &&opNames )
: Model( cp, NumExponents_ ),
  Object( std::move( OID ) ),
  bOverlapAltNorm{ cp.bOverlapAltNorm },
  NumOverlapExp{ NumOverlapExp_ },
  NumOverlapExpDual{ NumOverlapExp_ },
  vOverlap( 2 )
{
  for( int i = idxSrc; i <= idxSnk; ++i )
  {
    vOverlap[i].Key.Object.push_back( ObjectID( i ) );
    vOverlap[i].Key.Name = opNames[i];
    if( cp.OpMomIndependent( i ) )
    {
      MLU::Momentum pIgnore;
      pIgnore.Extract( vOverlap[i].Key.Object.back(), MLU::Momentum::DefaultPrefix );
    }
  }
  if( vOverlap[0].Key == vOverlap[1].Key )
  {
    vOverlap.resize( 1 );
    NumOverlapExpDual = 0;
  }
}

void ModelOverlap::AddParameters( Params &mp )
{
  if( NumOverlapExpDual )
  {
    AddParam( mp, vOverlap[idxSrc], NumOverlapExpDual );
    AddParam( mp, vOverlap[idxSnk], NumOverlapExpDual );
    if( NumOverlapExp > NumOverlapExpDual )
      AddParam( mp, vOverlap[2], NumOverlapExp - NumOverlapExpDual );
  }
  else
    AddParam( mp, vOverlap[0], NumOverlapExp );
}

void ModelOverlap::SaveParameters( const Params &mp )
{
  for( ModelParam &p : vOverlap )
    p.idx = mp.at( p.Key )();
}

// Get a descriptive string for the model
std::string ModelOverlap::Description() const
{
  std::string s;
  for( std::size_t i = vOverlap.size(); i-- > 0; )
  {
    s.append( 1, ',' );
    s.append( vOverlap[i].Key.Name );
  }
  return s;
}

void ModelOverlap::Guessable( ParamsPairs &PP ) const
{
  // Dual overlap factors are products
  if( NumOverlapExpDual )
    PP.KnowProduct( vOverlap[0].Key, vOverlap[1].Key, NumOverlapExpDual );
  // If I only have one overlap factor, I can guess it from the data (but not its sign)
  if( NumOverlapExp - NumOverlapExpDual )
    PP.SetState( ParamsPairs::State::AmbiguousSign, vOverlap[NumOverlapExpDual ? 2 : 0].Key,
                 NumOverlapExp - NumOverlapExpDual );
}

// If I can't work out what the overlap coefficients are, I might be able to determine their product
void ModelOverlap::ReduceUnknown( const ParamsPairs &PP )
{
  if( NumOverlapExpDual )
  {
    // I have at least one differing overlap coefficients
    std::size_t NewDual{ PP.NumKnownProducts( vOverlap[0].Key, vOverlap[1].Key, NumOverlapExpDual ) };
    if( NumOverlapExpDual != NewDual )
    {
      // I need to reduce my number of overlap coefficients
      if( NumOverlapExpDual == NumOverlapExp )
      {
        // I need to make a new name
        if( vOverlap[0].Key.Object.size() != 1 || vOverlap[1].Key.Object.size() != 1 )
          throw std::runtime_error( "ModelOverlap::ReduceUnknown() : Class invariant breached" );
        vOverlap.resize( 3 );
        vOverlap[2].Key = vOverlap[1].Key;
        if( !MLU::EqualIgnoreCase( vOverlap[0].Key.Object[0], vOverlap[1].Key.Object[0] ) )
          vOverlap[2].Key.Object.push_back( vOverlap[0].Key.Object[0] );
        if( !MLU::EqualIgnoreCase( vOverlap[0].Key.Name, vOverlap[1].Key.Name ) )
          vOverlap[2].Key.Name.append( vOverlap[0].Key.Name );
      }
      if( NewDual == 0 )
      {
        // I know the name of the combined operator and no longer have individual operators
        vOverlap[0] = std::move( vOverlap[2] );
        vOverlap.resize( 1 );
      }
      NumOverlapExpDual = NewDual;
    }
  }
}

/*std::size_t ModelOverlap::NumUnknown( std::vector<bool> &bKnown ) const
{
  std::size_t UnKnown{ 0 };
  for( const ModelParam &p : Overlap )
  {
    for( std::size_t i = 0; i < NumOverlapExp; ++i )
    {
      if( !bKnown[p.idx + i] )
        ++UnKnown;
    }
  }
  return UnKnown;
}*/
/***********************************
 
 Guess

{
  for( int i = 0; i < NumModelParams; ++i )
    modelParams[i] = 0;
  // Get all the constants and start keeping track of which parameters we have managed to take a guess for
  std::vector<bool> bKnown( NumModelParams, false );
  ds.GetFixed( Fold::idxCentral , modelParams, ParamFixed );
  for( const DataSet::FixedParam &p : ParamFixed )
    bKnown[p.idx] = true;
  const std::vector<bool> bFixed( bKnown );
  // Take a guess for E0 if we didn't manage to load it
  if( !bKnown[0] )
  {
    int E0Count{ 0 };
    Vector vCorr;
    for( int f = 0; f < NumFiles; ++f )
    {
      scalar z;
      vCorr.MapView( const_cast<scalar *>( ds.corr[f][Fold::idxCentral] ), ds.corr[f].Nt() );
      if( model[f]->GuessE0( z, vCorr ) )
      {
        modelParams[0] += z;
        E0Count++;
      }
    }
    assert( E0Count && "Nothing guessed E0" );
    modelParams[0] /= E0Count;
    bKnown[0] = true;
  }
  // Now guess excited state energies (increment of half the previous difference each time)
  {
    scalar EPrior{ 0 };
    for( int e = 1; e < NumExponents; ++e )
    {
      const int idxE{ e * NumPerExp };
      const int idxELast{ idxE - NumPerExp };
      if( !bKnown[idxE] )
      {
        scalar PriorDiff{ modelParams[idxELast] - EPrior };
        modelParams[idxE] = modelParams[idxELast] + PriorDiff / 2;
        EPrior = modelParams[idxELast];
        bKnown[idxE] = true;
      }
    }
  }
  // Now guess everything else
  {
    Vector vCorr;
    int pass{ 0 };
    for( bool bNeedAnotherPass = true; bNeedAnotherPass; ++pass )
    {
      bNeedAnotherPass = false;
      for( int f = 0; f < NumFiles; ++f )
      {
        vCorr.MapView( const_cast<scalar *>( ds.corr[f][Fold::idxCentral] ), ds.corr[f].Nt() );
        if( model[f]->Guess( modelParams, bKnown, pass, vCorr ) )
          bNeedAnotherPass = true;
      }
    }
  }
  // Make sure we've taken a guess for everything
  int iVar = 0;
  for( int i = 0; i < NumModelParams; ++i )
  {
    assert( bKnown[i] && "Bad model: unable to guess all parameters" );
    if( !bFixed[i] )
      parGuess.Add( ParamNames[i], vGuess.empty() ? modelParams[i] : vGuess[iVar++], 0 );
  }
}

 ***********************************/
