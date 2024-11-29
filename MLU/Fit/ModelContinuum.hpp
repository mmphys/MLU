/**
 
 Chiral continuum fit

 Source file: ModelContinuum.hpp
 
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

#ifndef ModelContinuum_hpp
#define ModelContinuum_hpp

#include "Continuum.hpp"

extern const std::string sFF;

enum class eModelType{ Unknown, Continuum };
std::ostream & operator<<( std::ostream &os, const eModelType m );
std::istream & operator>>( std::istream &is,       eModelType &m );

struct ContinuumFit;

/// Chiral continuum fit model
struct ModelContinuum : Model
{
  //friend class ContinuumFit;
  static constexpr scalar Lambda{ ContinuumFit::Lambda }; // Units GeV (for now this is hardwired)
  static constexpr scalar LambdaInv{ ContinuumFit::LambdaInv };
  static constexpr scalar LambdaInvSq{ ContinuumFit::LambdaInvSq };
  static const std::string FieldQSq;
  static const std::string FieldEL;
  const MLU::FormFactor ff;
  const ContinuumFit &Parent;
  const ModelFile &mf;
  const std::size_t idxFitColumn;
  const int idxFF;
  const MLU::Momentum p;
  const std::string Ensemble;
  const EnsembleInfo &ei;
  ModelParam EL, mH, mL, qSq;
protected:
  ModelParam aEL, amH, amL, aqSq;
  const EnsembleInfo &GetEnsembleInfo() const;
public:
  ModelContinuum( const Model::CreateParams &cp, Model::Args &Args, MLU::FormFactor ff );
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
  /// This is called once at the start of each replica
  void SetReplica( Vector &ScratchPad, Vector &ModelParams ) const override;
};

#endif // ModelContinuum_hpp
