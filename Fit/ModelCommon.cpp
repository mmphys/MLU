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
  if( !bManualObjectID[0] || !bManualObjectID[1] )
    throw std::runtime_error( "Model3pt constructor: Snk and Src parameters are required (TODO: fix this)" );
  if( Common::EqualIgnoreCase( ObjectID[0], ObjectID[1] ) )
    ObjectID.resize( 1 );
  return { ObjectID };
}

// 2nd phase of construction (now that virtual functions in place)
ModelOverlap::ModelOverlap( const Model::CreateParams &cp, Model::Args &Args,
                            std::vector<std::string> &&OID, std::size_t NumOverlapExp_ )
: Model( cp, Args ),
  Object( std::move( OID ) ),
  bNormalisationByEnergy{ Args.Remove( "ENorm", true ) },
  NumOverlapExp{ NumOverlapExp_ },
  Overlap( 2 )
{
  // Overlap coefficient names come from the correlator file name
  bool bForceSrcSnkDifferent;
  Args.Remove( "srcsnk", &bForceSrcSnkDifferent );
  for( int i = idxSrc; i <= idxSnk; ++i )
  {
    Overlap[i].Key.Object.push_back( ObjectID( i ) );
    Overlap[i].Key.Name = Args.Remove( pSrcSnk[i], cp.OpNames[cp.pCorr->Name_.op[i]] );
    if( bForceSrcSnkDifferent )
      Overlap[i].Key.Name.append( pSrcSnk[i] );
  }
  if( Overlap[0].Key == Overlap[1].Key )
    Overlap.resize( 1 );
}

void ModelOverlap::AddParameters( Params &mp )
{
  for( ModelParam &p : Overlap )
    AddParam( mp, p, NumOverlapExp );
}

void ModelOverlap::SaveParameters( const Params &mp )
{
  for( ModelParam &p : Overlap )
    p.idx = mp.at( p.Key )();
}

// Get a descriptive string for the model
std::string ModelOverlap::Description() const
{
  std::string s;
  for( std::size_t i = Overlap.size(); i-- > 0; )
  {
    s.append( 1, ',' );
    s.append( Overlap[i].Key.Name );
  }
  return s;
}

std::size_t ModelOverlap::Guessable( std::vector<bool> &bKnown, bool bLastChance ) const
{
  // Assume the energies are known (or we wouldn't be called)
  std::size_t NumUnknown{ 0 };
  if( Overlap.size() == 1 )
    for( std::size_t i = 0; i < Overlap[0].param->size; ++i )
      bKnown[Overlap[0].idx + i] = true;
  else if( bLastChance || bKnown[Overlap[0].idx] || bKnown[Overlap[1].idx] )
  {
    // I know at least one of them, so I in fact know both
    for( std::size_t i = 0; i < Overlap[0].param->size; ++i )
    {
      bKnown[Overlap[0].idx + i] = true;
      bKnown[Overlap[1].idx + i] = true;
    }
  }
  else
    NumUnknown = 2;
  return NumUnknown;
}

// If I can't work out what the overlap coefficients are, I might be able to determine their product
void ModelOverlap::ReduceUnknown()
{
  if( Overlap.size() > 1 )
  {
    if( Overlap[0].Key.Object.size() != 1 || Overlap[1].Key.Object.size() != 1 )
      throw std::runtime_error( "ModelOverlap::ReduceUnknown() : Class invariant breached" );
    if( !Common::EqualIgnoreCase( Overlap[0].Key.Object[0], Overlap[1].Key.Object[0] ) )
      Overlap[1].Key.Object.push_back( std::move( Overlap[0].Key.Object[0] ) );
    if( !Common::EqualIgnoreCase( Overlap[0].Key.Name, Overlap[1].Key.Name ) )
      Overlap[1].Key.Name.append( std::move( Overlap[0].Key.Name ) );
    Overlap[0].Key = std::move( Overlap[1].Key );
    Overlap.resize( 1 );
  }
}

std::size_t ModelOverlap::NumUnknown( std::vector<bool> &bKnown ) const
{
  std::size_t UnKnown{ 0 };
  for( const ModelParam &p : Overlap )
  {
    for( std::size_t i = 0; i < p.param->size; ++i )
    {
      if( !bKnown[p.idx + i] )
        ++UnKnown;
    }
  }
  return UnKnown;
}
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
