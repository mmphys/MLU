/*************************************************************************************
 
 Models for 2pt functions
 
 Source file: Model2pt.hpp
 
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

#ifndef Model2pt_hpp
#define Model2pt_hpp

#include "ModelCommon.hpp"

struct Model2pt : public ModelOverlap
{
  Model2pt( const Model::CreateParams &cp, Model::Args &Args )
  : Model2pt( cp, Args, GetObjectNameSingle( cp, Args ), Args.Remove( "e", cp.NumExponents, true ) ) {}
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown,
               const VectorView &FitData, std::vector<int> FitTimes, bool bLastChance ) const override;
  double Derivative( int t, int p ) const override;

protected:
  Model2pt( const Model::CreateParams &cp, Model::Args &Args,
            std::vector<std::string> &&ObjectID, std::size_t NumOverlapExp );
  scalar Estimate( Vector &Guess, const VectorView &FitData, std::vector<int> FitTimes,
                   std::size_t NumExp, std::size_t Timeslice ) const;
  ModelParam E;
};

struct ModelExp : public Model2pt
{
  ModelExp( const Model::CreateParams &cp, Model::Args &Args,
            std::vector<std::string> &&objectID, std::size_t NumOverlapExp )
  : Model2pt( cp, Args, std::move( objectID ), NumOverlapExp ) {}
  ModelExp( const Model::CreateParams &cp, Model::Args &Args )
  : ModelExp(cp, Args, GetObjectNameSingle( cp, Args ), Args.Remove( "e", cp.NumExponents, true )) {}
  ModelType Type() const override { return ModelType::Exp; }
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;
};

struct ModelCosh : public Model2pt
{
  ModelCosh( const Model::CreateParams &cp, Model::Args &Args,
             std::vector<std::string> &&objectID, std::size_t NumOverlapExp )
  : Model2pt( cp, Args, std::move( objectID ), NumOverlapExp ) {}
  ModelCosh( const Model::CreateParams &cp, Model::Args &Args )
  : ModelCosh(cp, Args, GetObjectNameSingle( cp, Args ), Args.Remove( "e", cp.NumExponents, true )) {}
  ModelType Type() const override { return ModelType::Cosh; }
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;
};

struct ModelSinh : public Model2pt
{
  ModelSinh( const Model::CreateParams &cp, Model::Args &Args,
             std::vector<std::string> &&objectID, std::size_t NumOverlapExp )
  : Model2pt( cp, Args, std::move( objectID ), NumOverlapExp ) {}
  ModelSinh( const Model::CreateParams &cp, Model::Args &Args )
  : ModelSinh(cp, Args, GetObjectNameSingle( cp, Args ), Args.Remove( "e", cp.NumExponents, true )) {}
  ModelType Type() const override { return ModelType::Sinh; }
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;
};

#endif // Model2pt_hpp
