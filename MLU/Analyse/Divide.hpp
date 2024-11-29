/**
 
 For now this divides arg 1 by arg 2. Might become generic mathematical manipulation

 Source file: Divide.hpp

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

#pragma once
#include <MLU/MLU.hpp>

using Scalar = double;
using Fold = MLU::Fold<Scalar>;
using Boot = MLU::Sample<std::complex<Scalar>>;

class Params
{
public:
  const bool bForceOverwrite;
  const bool bSaveHdf5;
  const std::string InBase;
  const std::string OutBase;
  using VFMap = std::map<std::string, double>;

protected:
  VFMap VolFactor;
  bool bFold = false;
  Fold fDenominator, fNumerator;
  Boot bDenominator, bNumerator;

  void DivideCommon( const Boot &Denominator, Boot &Numerator );
  void ReadNumerator( const std::string &FileName );
  template <typename T>
  void Normalise( MLU::JackBoot<T> &jb, const std::vector<int> &opNum,
             const std::vector<std::string> &opNames );
  template <typename FoldBoot>
  void DivideBoot( FoldBoot &Numerator, const FoldBoot &Denominator );
  void DivideFold( Fold &Numerator, const Fold &Denominator );

public:
  Params( const MLU::CommandLine &cl );
  void Run( const MLU::CommandLine &cl, bool &bShowUsage );
};
