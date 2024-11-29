/**
 
 Mike's lattice QCD utilities: Synthetic data generator (e.g. for fit test)
 
 Source file: SynthData.hpp
 
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

#ifndef SynthData_hpp
#define SynthData_hpp

#include <MLU/MLU.hpp>
//using MomentumMap = MLU::FileNameAtt::MomentumMap;
//using MomentumMapValue = typename MomentumMap::value_type;
using Algebra = MLU::Gamma::Algebra;

using Scalar = double;
using Fold = MLU::Fold<Scalar>;
using Matrix = MLU::Matrix<Scalar>;

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
  const MLU::SeedType seed;
  const bool bVerbose;
  const bool bWarnIfExists;
  Synth( const MLU::CommandLine &cl, const std::string MachineNameActual );
  void Make( std::string Basename ) const;
  std::ostream & DumpParams( std::ostream &os ) const;
protected:
  std::array<std::array<MeanSigma, NumParams>, MaxExp> ParamMS;
  Fold RandomSample;
  MLU::SeedType GetSeedType( const MLU::CommandLine &cl );
  std::string GetParamName( int e, int p ) const;
};

#endif // SynthData_hpp
