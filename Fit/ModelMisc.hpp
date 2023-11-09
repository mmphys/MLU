/*************************************************************************************
 
 A place to implement miscellaneous models
 (This version for regression testing with Taku Izubuchi Nov 2023)

 Source file: ModelMisc.hpp
 
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

#ifndef ModelMisc_hpp
#define ModelMisc_hpp

#include "ModelCommon.hpp"

/// A place to build miscellaneous models while testing
struct ModelMisc : Model, Object
{
  ModelMisc( const Model::CreateParams &cp, Model::Args &Args, int NumExponents );
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                     const VectorView &FitData, std::vector<int> FitTimes,
                     bool bLastChance ) const override;
  ModelType Type() const override;
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;
private:
  std::array<ModelParam, 3> ap;
};

#endif // ModelMisc_hpp
