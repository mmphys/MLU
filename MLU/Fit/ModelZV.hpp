/**
 
 Models for ZV extraction

 Source file: ModelZV.hpp
 
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

#ifndef ModelZV_hpp
#define ModelZV_hpp

#include "Model2pt.hpp"
#include "Model3pt.hpp"

struct ModelZV : public Model3pt
{
  // static const std::string sRaw;
  static const std::string sZV;
  static const std::string sChildModel;
  ModelZV( const Model::CreateParams &cp, Model::Args &Args, int NumExponents,
              std::vector<std::string> &&objectID, std::vector<std::string> &&opNames );
  ModelZV( const Model::CreateParams &cp, Model::Args &Args, int NumExp );
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  double Derivative( int t, int p ) const override;
  ModelType Type() const override;
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;

protected:
  const eModelType ChildType;
  std::unique_ptr<Model2pt> C2;
  // Derived parameters
  ModelParam ZV;
  static std::vector<std::string> GetObjectNameZV( const Model::CreateParams &cp, Model::Args &Args );
};

#endif // ModelZV_hpp
