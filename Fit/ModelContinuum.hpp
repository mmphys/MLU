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

struct ContinuumFit;

/// Chiral continuum fit model
struct ModelContinuum : Model
{
  static constexpr int NumConst{ 6 };
  static constexpr int CChiral{ 0 };
  static constexpr int CMPi{ 1 };
  static constexpr int CDiscret{ 4 };
  static constexpr int CEOnL{ 2 };
  static constexpr int CEOnL2{ 3 };
  static constexpr int CEOnL3{ 5 };
  static constexpr scalar FourPi{ 4. * M_PI };
  friend class ContinuumFit;
  static const std::string FieldQSq;
  static const std::string FieldEL;
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
  // Cache values based solely on the model parameters (to speed up computation)
  int GetScratchPadSize() const override { return 1; }
  void ModelParamsChanged( Vector &ScratchPad, const Vector &ModelParams ) const override;
protected:
  static constexpr scalar Lambda{ 1e9 }; // Units GeV (for now this is hardwired)
  static constexpr scalar LambdaInv{ 1. / Lambda }; // For now this is hardwired
  const ContinuumFit &Parent;
  const ModelFile &mf;
  const Common::FormFactor ff;
  const int idxFF;
  const Common::Momentum p;
  const std::string Ensemble;
  unsigned int aInv_L;
  std::array<std::string, 2> Meson;
  ModelParam aInv, fPi, mPi, mPDGPi, FVSim, FVPhys, mPDGH, mPDGL, Delta;
  ModelParam EL, kMu, mH, mL, qSq;
  ModelParam aEL, akMu, amH, amL, aqSq;
  std::array<ModelParam, NumConst> c;
  inline scalar DeltaF( scalar M, scalar FV ) const
  {
    const scalar ChiralLog{ 2. * M * M * std::log( std::abs( M * LambdaInv ) ) };
    return -0.75 * ( ChiralLog + FV );
  }
};

#endif // ModelContinuum_hpp
