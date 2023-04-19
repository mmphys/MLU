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

#include "ModelCommon.hpp"

std::vector<std::string> Object::GetObjectNameSingle( const Model::CreateParams &cp, Model::Args &Args )
{
  bool bManualObjectID;
  std::string ObjectID{ Args.Remove( "ObjectID", &bManualObjectID ) };
  if( !bManualObjectID )
  {
    ObjectID = cp.pCorr->Name_.Base;
    std::size_t Len{ ObjectID.find( '.' ) };
    if( Len != std::string::npos )
      ObjectID.resize( Len );
  }
  return { ObjectID };
}

std::vector<std::string> Object::GetObjectNameSnkSrc( const Model::CreateParams &cp, Model::Args &Args )
{
  // Now remove object names
  bool bManualObjectID[2];
  std::vector<std::string> ObjectID( 2 );
  ObjectID[0] = Args.Remove( "Src", &bManualObjectID[0] );
  ObjectID[1] = Args.Remove( "Snk", &bManualObjectID[1] );
  if( !bManualObjectID[0] )
  {
    if( cp.pCorr->Name_.MesonMom.empty() )
      throw std::runtime_error( "GetObjectNameSnkSrc(): Src unavailable - specify manually" );
    ObjectID[0] = cp.pCorr->Name_.MesonMom[0];
  }
  if( !bManualObjectID[1] )
  {
    if( cp.pCorr->Name_.Meson.size() < 2 )
      throw std::runtime_error( "GetObjectNameSnkSrc(): Snk unavailable - specify manually" );
    ObjectID[1] = cp.pCorr->Name_.MesonMom[1];
  }
  if( Common::EqualIgnoreCase( ObjectID[0], ObjectID[1] ) )
    ObjectID.resize( 1 );
  return { ObjectID };
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
  {
    ParamsPairs::Key key0{ vOverlap[0].Key, 0 };
    ParamsPairs::Key key1{ vOverlap[1].Key, 0 };
    for( key0.Index = 0; key0.Index < NumOverlapExpDual; ++key0.Index )
    {
      key1.Index = key0.Index;
      PP.KnowProduct( key0, key1 );
    }
  }
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
        if( !Common::EqualIgnoreCase( vOverlap[0].Key.Object[0], vOverlap[1].Key.Object[0] ) )
          vOverlap[2].Key.Object.push_back( vOverlap[0].Key.Object[0] );
        if( !Common::EqualIgnoreCase( vOverlap[0].Key.Name, vOverlap[1].Key.Name ) )
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
