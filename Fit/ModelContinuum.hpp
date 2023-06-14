/*************************************************************************************
 
 Chiral continuum fit

 Source file: ModelContinuum.hpp
 
 Copyright (C) 2023

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

#ifndef ModelContinuum_hpp
#define ModelContinuum_hpp

#include "Model.hpp"

extern const std::string sFF;

enum class eModelType{ Unknown, Continuum };
std::ostream & operator<<( std::ostream &os, const eModelType m );
std::istream & operator>>( std::istream &is,       eModelType &m );

/// Chiral continuum fit model
struct ModelContinuum : Model
{
  ModelContinuum( const Model::CreateParams &cp, Model::Args &Args, Common::FormFactor ff );
  void DefineXVector( DataSet &ds, int i ) override;
  const std::string &XVectorKeyName() const override;
  std::string XVectorKey() const override;
  std::vector<Param::Key> XVectorKeyNames() const override;
  void AddParameters( Params &mp ) override;
  std::size_t GetFitColumn() const override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                     const VectorView &FitData, std::vector<int> FitTimes,
                     bool bLastChance ) const override;
  ModelType Type() const override;
  scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const override;
protected:
  static constexpr scalar Lambda{ 1e9 }; // Units GeV (for now this is hardwired)
  static constexpr scalar LambdaInv{ 1. / Lambda }; // For now this is hardwired
  const ModelFile &mf;
  const Common::FormFactor ff;
  const Common::Momentum p;
  const std::string Ensemble;
  unsigned int aInv_L;
  std::array<std::string, 2> Meson;
  ModelParam aInv, fPi, mPi, mPDGPi, mPDGDStar, mPDGH, mPDGL, Delta;
  ModelParam EL, kMu, mH, mL, qSq;
  ModelParam aEL, akMu, amH, amL, aqSq;
  ModelParam c0, c1, c2, c3, c4;
  scalar DeltaF( scalar M ) const;
};

#endif // ModelContinuum_hpp
