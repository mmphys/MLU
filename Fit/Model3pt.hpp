/*************************************************************************************
 
 Models for 3pt functions

 Source file: Model3pt.hpp
 
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

#ifndef Model3pt_hpp
#define Model3pt_hpp

#include "ModelCommon.hpp"

struct Model3pt : public ModelOverlap
{
  Model3pt( const Model::CreateParams &cp, Model::Args &Args,
            std::vector<std::string> &&objectID, std::vector<std::string> &&opNames,
            std::size_t NumOverlapExp );
  Model3pt( const Model::CreateParams &cp, Model::Args &Args )
  : Model3pt( cp, Args, GetObjectNameSnkSrc( cp, Args ), ModelOverlap::GetOpNames( cp, Args ),
              Args.Remove( "e", cp.NumExponents, true ) ) {}
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown,
               const VectorView &FitData, std::vector<int> FitTimes, bool bLastChance ) const override;
  double Derivative( int t, int p ) const override;
  ModelType Type() const override { return ModelType::ThreePoint; }
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;

protected:
  int DeltaT;
  //std::size_t ParamIndex( std::size_t idxSnk, std::size_t idxSrc ) const;
  //std::size_t NumUnknown( std::vector<bool> &bKnown ) const;
  std::vector<ModelParam> E;
  ModelParam MEL;
  // Derived parameters
  ModelParam EDiff;
};

#endif // Model3pt_hpp
