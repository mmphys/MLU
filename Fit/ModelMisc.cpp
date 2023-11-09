/*************************************************************************************
 
 A place to implement miscellaneous models
 (This version for regression testing with Taku Izubuchi Nov 2023)

 Source file: ModelMisc.cpp
 
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

#include "ModelMisc.hpp"

ModelMisc::ModelMisc( const Model::CreateParams &cp, Model::Args &Args, int NumExponents_ )
: Model( cp, NumExponents_ ),
  Object( GetObjectNameSingle( cp, Args ) )
{
  for( std::size_t i = 0; i < ap.size(); ++i )
  {
    ap[i].Key.Object.push_back( objectID[0] );
    ap[i].Key.Name = std::string( 1, 'A' + i );
  }
}

void ModelMisc::AddParameters( Params &mp )
{
  for( ModelParam &p : ap )
    AddParam( mp, p, 1 );
}

void ModelMisc::SaveParameters( const Params &mp )
{
  for( ModelParam &p : ap )
    p.idx = mp.at( p.Key )();
}

// Get a descriptive string for the model
std::string ModelMisc::Description() const
{
  return objectID[0] + ",A,B,C";
}

void ModelMisc::Guessable( ParamsPairs &PP ) const
{
}

std::size_t ModelMisc::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                              const VectorView &FitData, std::vector<int> FitTimes,
                              bool bLastChance ) const
{
  std::size_t NumUnknown{};
  for( const ModelParam &p : ap )
  {
    if( !bKnown[p.idx] )
      NumUnknown++;
  }
  return NumUnknown;
}

ModelType ModelMisc::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::Exp );
  return m;
}

scalar ModelMisc::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  const scalar A{ ModelParams[ap[0].idx] };
  const scalar B{ ModelParams[ap[1].idx] };
  const scalar C{ ModelParams[ap[2].idx] };
  return A * ( 1 + B * std::exp( - C * t ) );
}
