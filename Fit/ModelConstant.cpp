/*************************************************************************************
 
 Fit to a constant
 
 Source file: ModelConstant.cpp

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

#include "ModelConstant.hpp"

ModelConstant::ModelConstant( const Model::CreateParams &cp, Model::Args &Args )
: Model( cp, Args ), Object( GetObjectNameSingle( cp, Args ) )
{
  Constant.Key.Object = { ObjectID( idxSrc ) };
  Constant.Key.Name = Args.Remove( "const", Common::ModelBase::ConstantPrefix );
}

void ModelConstant::AddParameters( Params &mp )
{
  AddParam( mp, Constant, 1 );
}

void ModelConstant::SaveParameters( const Params &mp )
{
  Constant.idx = mp.at( Constant.Key )();
}

std::string ModelConstant::Description() const
{
  std::string s{ "K(" };
  s.append( ObjectID( idxSrc ) );
  s.append( 1, ',' );
  s.append( Constant.Key.Name );
  s.append( 1, ')' );
  return s;
}

void ModelConstant::Guessable( ParamsPairs &PP ) const
{
  PP.SetState( ParamsPairs::State::Known, Constant.Key, Constant.param->size );
}

std::size_t ModelConstant::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                                  const VectorView &FitData, std::vector<int> FitTimes,
                                  bool bLastChance ) const
{
  if( !bKnown[Constant.idx] )
  {
    scalar Sum{ 0 };
    for( scalar z : FitTimes )
      Sum += z;
    Guess[Constant.idx] = Sum / FitTimes.size();
    bKnown[Constant.idx] = true;
  }
  return 0;
}

scalar ModelConstant::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  return ModelParams[Constant.idx];
}

double ModelConstant::Derivative( int t, int p ) const
{
  return p == Constant.idx ? 1 : 0;
}
