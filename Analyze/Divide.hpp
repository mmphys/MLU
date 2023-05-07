/*************************************************************************************
 
 For now this divides arg 1 by arg 2. Might become generic mathematical manipulation
 Source file: Divide.cpp
 Copyright (C) 2021
 Author: Michael Marshall<Michael.Marshall@ed.ac.uk>
 
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

#pragma once
#include <MLU/Common.hpp>

using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Boot = Common::Sample<std::complex<Scalar>>;

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
  void Normalise( Common::JackBoot<T> &jb, const std::vector<int> &opNum,
             const std::vector<std::string> &opNames );
  template <typename FoldBoot>
  void DivideBoot( FoldBoot &Numerator, const FoldBoot &Denominator );
  void DivideFold( Fold &Numerator, const Fold &Denominator );

public:
  Params( const Common::CommandLine &cl );
  void Run( const Common::CommandLine &cl, bool &bShowUsage );
};
