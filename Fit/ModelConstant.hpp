/*************************************************************************************
 
 Fit to a constant
 
 Source file: ModelConstant.hpp
 
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

#ifndef ModelConstant_hpp
#define ModelConstant_hpp

#include "ModelCommon.hpp"

struct ModelConstant : Model, Object
{
  ModelConstant( const Model::CreateParams &cp, Model::Args &Args );
  void AddParameters( struct Params &mp ) override;
  void SaveParameters( const struct Params &mp ) override;
  std::string Description() const override;
  std::size_t Guessable( std::vector<bool> &bKnown, bool bLastChance ) const override;
  std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown,
               const Vector &FitData, std::vector<int> FitTimes, bool bLastChance ) const override;
  ModelType Type() const override { return ModelType::Constant; }
  scalar operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const override;
  double Derivative( int t, int p ) const override;
protected:
  ModelParam Constant;
};

#endif // ModelConstant_hpp
