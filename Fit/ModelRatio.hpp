/*************************************************************************************
 
 Models for Ratios

 Source file: ModelRatio.hpp
 
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

#ifndef ModelRatio_hpp
#define ModelRatio_hpp

#include "Model2pt.hpp"
#include "Model3pt.hpp"

struct ModelRatio : public Model3pt
{
  static const std::string sRaw;
  static const std::string sR3;
  static const std::string sChildModel;
  ModelRatio( const Model::CreateParams &cp, Model::Args &Args, int NumExponents,
              std::vector<std::string> &&objectID, std::vector<std::string> &&opNames );
  ModelRatio( const Model::CreateParams &cp, Model::Args &Args, int NumExp )
  : ModelRatio(cp, Args, NumExp, GetObjectNameSnkSrc(cp, Args), ModelOverlap::GetOpNames(cp, Args)){}
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  double Derivative( int t, int p ) const override;
  ModelType Type() const override;
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;

protected:
  const eModelType ChildType;
  bool bRaw;
  std::vector<std::unique_ptr<Model2pt>> C2;
  // Derived parameters
  ModelParam R3Raw;
};

#endif // ModelRatio_hpp
