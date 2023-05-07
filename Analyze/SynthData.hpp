/**
 
 Mike's lattice QCD utilities: Synthetic data generator (e.g. for fit test)
 
 Source file: SynthData.hpp
 
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

#ifndef SynthData_hpp
#define SynthData_hpp

#include <MLU/Common.hpp>
//using MomentumMap = Common::FileNameAtt::MomentumMap;
//using MomentumMapValue = typename MomentumMap::value_type;
using Algebra = Common::Gamma::Algebra;
//using namespace Common;

using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Matrix = Common::Matrix<Scalar>;

#include <cmath>
#include <iomanip>
//#include <mutex> // Apparently empty under __INTEL_COMPILER
//#include <set>
//#include <filesystem> // C++17
//#include <LatAnalyze/Core/OptParser.hpp>
//#include <LatAnalyze/Statistics/Dataset.hpp>
//#include <LatAnalyze/Io/Io.hpp>
//#include <LatAnalyze/Io/Hdf5File.hpp>

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif

struct MeanSigma
{
  Scalar Mean;
  Scalar Sigma;
  MeanSigma( Scalar mean, Scalar sigma ) : Mean{mean}, Sigma{sigma} {}
  MeanSigma() : MeanSigma( 1, 0 ) {}
};

std::ostream & operator<<( std::ostream &os, const MeanSigma &mv );
std::istream & operator>>( std::istream &is, MeanSigma &mv );

struct Synth
{
  static constexpr int NumOps{2};
  static const std::array<std::string, NumOps> opNames;
  static constexpr int NumParams{ NumOps + 1 };
  static constexpr int ParamE{ 0 };
  static constexpr int ParamA{ ParamE + 1 };
  static constexpr int MaxExp{2};

  const int NumExp;
  const Scalar CorrelatorJitter;
  const std::string MachineName;
  const int nSample;
  const int nBootSample;
  const std::string outStem;
  const Common::SeedType seed;
  const bool bVerbose;
  const bool bWarnIfExists;
  Synth( const Common::CommandLine &cl, const std::string MachineNameActual );
  void Make( std::string Basename ) const;
  std::ostream & DumpParams( std::ostream &os ) const;
protected:
  std::array<std::array<MeanSigma, NumParams>, MaxExp> ParamMS;
  Fold RandomSample;
  Common::SeedType GetSeedType( const Common::CommandLine &cl );
  std::string GetParamName( int e, int p ) const;
};

#endif // SynthData_hpp
