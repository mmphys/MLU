/**
 
 Models for 3pt functions

 Source file: Model3pt.hpp
 
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

#ifndef Model3pt_hpp
#define Model3pt_hpp

#include "ModelCommon.hpp"

struct Model3pt : public ModelOverlap
{
  Model3pt( const Model::CreateParams &cp, Model::Args &Args, int NumExponents,
            std::vector<std::string> &&objectID, std::vector<std::string> &&opNames );
  Model3pt( const Model::CreateParams &cp, Model::Args &Args, int NumExponents )
  : Model3pt( cp, Args, NumExponents, GetObjectNameSnkSrc( cp, Args ),
              ModelOverlap::GetOpNames( cp, Args ) ) {}
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                     const VectorView &FitData, std::vector<int> FitTimes,
                     bool bLastChance ) const override;
  double Derivative( int t, int p ) const override;
  ModelType Type() const override;
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;

protected:
  const int N;
  const MLU::DispersionType dispType;
  int DeltaT;
  //std::size_t ParamIndex( std::size_t idxSnk, std::size_t idxSrc ) const;
  //std::size_t NumUnknown( std::vector<bool> &bKnown ) const;
  std::vector<ModelParam> E;
  ModelParam MEL;
  // Derived parameters
  ModelParam EDiff;
};

#endif // Model3pt_hpp
