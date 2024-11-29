/**
 
 Mike's lattice QCD utilities: Import correlators
 
 Source file: ImportCorr.hpp
 
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

#ifndef ImportCorr_hpp
#define ImportCorr_hpp

#include <cmath>
#include <iomanip>
#include <MLU/MLU.hpp>

using Scalar = double;
using Fold = MLU::Fold<Scalar>;
using Algebra = MLU::Gamma::Algebra;
using Vector = MLU::Vector<Scalar>;
using Matrix = MLU::Matrix<Scalar>;

struct CorrInfo
{
  std::string id;
  std::string opSrc;
  std::string opSnk;
  MLU::Parity parity;
};

/*****************************************************************

Importer

*****************************************************************/

struct Importer
{
  const MLU::SeedType Seed;
  const std::string outStem;
  const bool bDebug;
protected:
  int Nt;
  int NtOut;
  std::string GroupName;
  static const std::vector<CorrInfo> corrInfo; // The correlators I want to import - in order
  std::vector<std::vector<int>> vFitRanges;
  std::vector<std::vector<Matrix>> vRawData;
  Vector vBinnedCentral;
  Matrix mBinnedData;
  Vector vUnbinnedCentral;
  Matrix mUnbinnedData;
  template <typename T>
  void ReadCorrVector( ::H5::Group &g, const std::string &Name, std::vector<T> &v );
  void DebugDump( const std::vector<Matrix> &v ) const;
  void DebugDump( const std::vector<std::vector<Matrix>> &v ) const;
  void ReadInput( const std::string &Filename );
  void SpreadData( std::vector<Fold> &out, const Scalar *pSource, int idx, int NumSamples );
  void SpreadDataRawBoot( std::vector<Fold> &out, const Matrix &RawData );
  void SaveRawData( std::vector<Fold> &out, bool bPreserveSign );
  void Write( const std::string &Base, bool bPreserveSign );
public:
  Importer( const MLU::CommandLine &cl );
  void Import( const std::string &Filename, bool bPreserveSign ); // Tobi format
  void Import( const std::string &Filename, const std::string &Group, const std::string &DS );
};

#endif // ImportCorr_hpp
