/**
 
 Mike's lattice QCD utilities: Import correlators
 
 Source file: ImportCorr.hpp
 
 Copyright (C) 2022

 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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
**/

#ifndef ImportCorr_hpp
#define ImportCorr_hpp

#include <cmath>
#include <iomanip>
#include <MLU/Common.hpp>

using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Algebra = Common::Gamma::Algebra;
using Vector = Common::Vector<Scalar>;
using Matrix = Common::Matrix<Scalar>;

struct CorrInfo
{
  std::string id;
  std::string opSrc;
  std::string opSnk;
  Common::Parity parity;
};

/*****************************************************************

Importer

*****************************************************************/

struct Importer
{
  const Common::SeedType Seed;
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
  Importer( const Common::CommandLine &cl );
  void Import( const std::string &Filename, bool bPreserveSign );
};

#endif // ImportCorr_hpp
