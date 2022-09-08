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
: ModelOverlap( cp, Args, GetObjectNameSnkSrc( cp, Args ), NumExponents > 1 ? 2 : 1 )
{
  if( NumExponents > 3 )
    throw std::runtime_error( "Model3pt supports a maximum of 3 exponentials: "
                              "Gnd-Gnd, Gnd-Ex (and ^\\dag), Ex-Ex" );
  std::string ESrc = Args.Remove( "ESrc", ::E );
  std::string ESnk = Args.Remove( "ESnk", ::E );
  const bool bESame{ objectID.size() == 1 && Common::EqualIgnoreCase( ESrc, ESnk ) };
  E.resize( bESame ? 1 : 2 );
  E[idxSrc].Key.Object = ObjectID( idxSrc );
  E[idxSrc].Key.Name = std::move( ESrc );
  if( !bESame )
  {
    E[idxSnk].Key.Object = ObjectID( idxSnk );
    E[idxSnk].Key.Name = std::move( ESnk );
  }
  // Now create the matrix element
  for( std::size_t i = 0; i < objectID.size(); ++i )
  {
    if( i )
      MEL.Key.Object.append( 1, '-' );
    MEL.Key.Object.append( objectID[i] );
  }
  static const std::string sMEL{ "MEL" };
  MEL.Key.Name = Args.Remove( sMEL, sMEL );
}

void Model3pt::AddParameters( struct Params &mp )
{
  for( ModelParam &p : E )
    p.it = mp.Add( p.Key, NumOverlapExp, true );
  MEL.it = mp.Add( MEL.Key, NumExponents, true );
  ModelOverlap::AddParameters( mp );
}

void Model3pt::SaveParameters( const struct Params &mp )
{
  for( ModelParam &p : E )
    p.idx = p.it->second();
  MEL.idx = MEL.it->second();
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
  for( std::size_t i = 0; i < MEL.it->second.size; ++i )
  {
    if( !bKnown[MEL.idx + i] )
      bKnown[MEL.idx + i] = true;
  }
  return 0; // TODO: Implement guessing
}

std::size_t Model3pt::Guess( Vector &Guess, std::vector<bool> &bKnown,
                       const Vector &FitData, std::vector<int> FitTimes, bool bLastChance ) const
{
  // I can guess the matrix element
  for( std::size_t i = 0; i < MEL.it->second.size; ++i )
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
  return 0;
}

// TODO: Not sure I use this. Inconsistent definitions of how many parameters
std::size_t Model3pt::NumUnknown( std::vector<bool> &bKnown ) const
{
  std::size_t UnKnown{ ModelOverlap::NumUnknown( bKnown ) };
  for( const ModelParam &p : E )
  {
    for( std::size_t i = 0; i < p.it->second.size; ++i )
    {
      if( !bKnown[p.idx + i] )
        ++UnKnown;
    }
  }
  return UnKnown;
}
