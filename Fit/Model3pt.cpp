/*************************************************************************************
 
 Models for 3pt functions

 Source file: Model3pt.cpp
 
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

#include "Model3pt.hpp"

Model3pt::Model3pt( const Model::CreateParams &cp, Model::Args &Args )
: ModelOverlap( cp, Args, GetObjectNameSnkSrc( cp, Args ),
                Args.Remove( "e", cp.NumExponents, true ) > 1 ? 2 : 1 )
{
  if( !cp.pCorr->Name_.bGotDeltaT )
    throw std::runtime_error( "DeltaT not available in " + cp.pCorr->Name_.Filename );
  DeltaT = cp.pCorr->Name_.DeltaT;
  if( NumExponents > 3 )
    throw std::runtime_error( "Model3pt supports a maximum of 3 exponentials: "
                              "Gnd-Gnd, Gnd-Ex (and ^\\dag), Ex-Ex" );
  std::string ESrc = Args.Remove( "ESrc", ::E );
  std::string ESnk = Args.Remove( "ESnk", ::E );
  const bool bESame{ objectID.size() == 1 && Common::EqualIgnoreCase( ESrc, ESnk ) };
  E.resize( bESame ? 1 : 2 );
  E[idxSrc].Key.Object = { ObjectID( idxSrc ) };
  E[idxSrc].Key.Name = std::move( ESrc );
  if( !bESame )
  {
    E[idxSnk].Key.Object = { ObjectID( idxSnk ) };
    E[idxSnk].Key.Name = std::move( ESnk );
  }
  // Now create the matrix element
  MEL.Key.Object = objectID;
  static const std::string sMEL{ "MEL" };
  MEL.Key.Name = Args.Remove( sMEL, sMEL );
}

void Model3pt::AddParameters( Params &mp )
{
  for( ModelParam &p : E )
    AddParam( mp, p, NumOverlapExp, true );
  AddParam( mp, MEL, NumExponents, false );
  ModelOverlap::AddParameters( mp );
}

void Model3pt::SaveParameters( const Params &mp )
{
  for( ModelParam &p : E )
    p.idx = mp.at( p.Key )();
  MEL.idx = mp.at( MEL.Key )( 0, 0 );
  ModelOverlap::SaveParameters( mp );
}

// Get a descriptive string for the model
std::string Model3pt::Description() const
{
  std::string s{ "C3(" };
  s.append( ObjectID( idxSrc ) );
  if( objectID.size() > 1 )
  {
    s.append( 1, ',' );
    s.append( ObjectID( idxSnk ) );
  }
  s.append( ModelOverlap::Description() );
  s.append( 1, ')' );
  return s;
}

// For now, all we can guess are matrix elements - and we assume they are 1
std::size_t Model3pt::Guessable( std::vector<bool> &bKnown, bool bLastChance ) const
{
  // I can guess the matrix element
  for( std::size_t i = 0; i < MEL.param->size; ++i )
    bKnown[MEL.idx + i] = true;
  return 0; // TODO: Implement guessing
}

std::size_t Model3pt::Guess( Vector &Guess, std::vector<bool> &bKnown,
                       const VectorView &FitData, std::vector<int> FitTimes, bool bLastChance ) const
{
  // I can guess the matrix element
  for( std::size_t i = 0; i < MEL.param->size; ++i )
  {
    if( !bKnown[MEL.idx + i] )
    {
      bKnown[MEL.idx + i] = true;
      Guess[MEL.idx + i] = 1;
    }
  }
  return 0; // TODO: Implement guessing
}

// TODO: Make sure we add partial derivatives for all models
double Model3pt::Derivative( int t, int p ) const
{
  return 0;
}

scalar Model3pt::operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const
{
  scalar z = 0;
  for( int eSnk = 0; eSnk < NumOverlapExp; ++eSnk )
  {
    const scalar &ESnk{ ModelParams[E.back().idx + eSnk] };
    const scalar ASnk = ModelParams[Overlap.back().idx + eSnk];
    const scalar expSnk{ - ESnk * ( DeltaT - t ) };
    for( int eSrc = 0; eSrc < NumOverlapExp; ++eSrc )
    {
      if( !eSnk || !eSrc || NumExponents >= 3 )
      {
        const scalar &ESrc{ ModelParams[E[idxSrc].idx + eSrc] };
        const scalar ASrc = ModelParams[Overlap[idxSrc].idx + eSrc];
        const scalar expSrc{ - ESrc * t };
        const std::size_t idxMEL = ( *MEL.param )( eSnk, eSrc );
        scalar d{ ASnk * ASrc * ModelParams[idxMEL] * std::exp( expSnk + expSrc ) };
        if( bNormalisationByEnergy )
          d /= 4 * ESnk * ESrc;
        z += d;
      }
    }
  }
  return z;
}

std::size_t Model3pt::ParamIndex( std::size_t idxSnk, std::size_t idxSrc ) const
{
  if( idxSnk > 1 || idxSrc > 1 )
    throw std::runtime_error( "Model3pt::ParamIndex index out of bounds" );
  std::size_t idx{ idxSnk + idxSrc };
  if( idxSnk && E.size() > 1 )
    idx++;
  return idx;
}

// TODO: Not sure I use this. Inconsistent definitions of how many parameters
std::size_t Model3pt::NumUnknown( std::vector<bool> &bKnown ) const
{
  std::size_t UnKnown{ ModelOverlap::NumUnknown( bKnown ) };
  for( const ModelParam &p : E )
  {
    for( std::size_t i = 0; i < p.param->size; ++i )
    {
      if( !bKnown[p.idx + i] )
        ++UnKnown;
    }
  }
  return UnKnown;
}
