/**
 
 Fit to a constant
 
 Source file: ModelConstant.cpp

 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "ModelConstant.hpp"

ModelConstant::ModelConstant( const Model::CreateParams &cp, Model::Args &Args, int NumExponents )
: Model( cp, 1 ), Object( GetObjectNameSingle( cp, Args ) ) // Silently ignore NumExponents
{
  MLU::Param::Key kDefault{ ObjectID( idxSrc ), MLU::ModelBase::ConstantPrefix };
  Constant.Key = Args.Remove<MLU::Param::Key>( "const", kDefault );
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
    for( scalar z : FitData )
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

ModelType ModelConstant::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::Constant );
  return m;
}

double ModelConstant::Derivative( int t, int p ) const
{
  return p == Constant.idx ? 1 : 0;
}
