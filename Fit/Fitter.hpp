/*************************************************************************************
 
 Base class for fitter
 
 Source file: Fitter.cpp

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

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Covar.hpp"
#include "ModelBase.hpp"

// Forward declarations to minimise coupling
class FitterThread;       // OpenMP thread that will perform fit on each replica
class ModelDefaultParams; // Default parameters passed to models when creating

class Fitter
{
public:
  // Simple command-line options
  const bool bAnalyticDerivatives;
  const double HotellingCutoff;
  const double RelEnergySep;
  const int MinDof;
  const int Retry;
  const int MaxIt;
  const double Tolerance;
  const bool bSaveCorr;
  const bool bSaveCMat;
  const int Verbosity;
  const bool bForceSrcSnkDifferent;
  const std::vector<scalar> vGuess;
  // More complex command-line options
  const DataSet &ds;
  const int NumFiles; // Number of correlator files in the dataset
  const std::vector<std::string> &OpNames;
  const int NumOps;
  std::vector<ModelPtr> model;      // Model for each correlator
  const int NumExponents;
  const std::vector<std::string> PerExpNames; // Names of all the per-energy-level parameters ("E" is first, rest sorted)
  const std::vector<std::string> ParamNames; // Names of all parameters, with a numeric suffix for per-energy-level params
  const std::vector<DataSet::FixedParam> ParamFixed; // Map from constants to parameters
  const std::vector<int> ParamVariable;     // Map from fitting engine parameters to parameters
  const int NumModelParams; // Total number of parameters = NumFixed + NumVariable = NumExponents * NumPerExp +
  const int NumPerExp;      // Number of parameters per exponent (Needed when sorting)
  const int NumOneOff;      // Number of one-off (i.e. not per-exponent) parameters
  const int NumFixed;       // Number of fixed parameters
  const int NumVariable;    // Number of variable parameters (used by the fitting engines)
  CovarParams cp;

  // These variables set once at the start of each fit
  int dof;
  bool bCorrelated;

  explicit Fitter( const Common::CommandLine &cl, const DataSet &ds_,
                   const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames_,
                   CovarParams &&cp );
  virtual ~Fitter() {}

protected:
  // Used during construction (so that we can make the results const)
  std::vector<ModelPtr> CreateModels( const std::vector<std::string> &ModelArgs,
                                      const ModelDefaultParams &modelDefault );
  int GetNumExponents();
  std::size_t EnsureModelsSolubleHelper( UniqueNames &Names, std::size_t &NumWithUnknowns );
  std::vector<std::string> MakePerExpNames();
  std::vector<std::string> MakeParamNames();
  std::vector<DataSet::FixedParam> MakeParamFixed();
  std::vector<int> MakeParamVariable();
  virtual FitterThread * MakeThread( bool bCorrelated, ModelFile &OutputModel, vCorrelator &CorrSynthetic ) = 0;

public:
  static void SayNumThreads( std::ostream &os );
  bool Dump( int idx ) const;
  void Dump( int idx, const std::string &Name, const Matrix &m ) const;
  void Dump( int idx, const std::string &Name, const Vector &v ) const;
  void SaveMatrixFile( const Matrix &m, const std::string &Type, const std::string &Filename,
                         const char *pGnuplotExtra = nullptr ) const;
  std::vector<Common::ValWithEr<scalar>>
  PerformFit( bool bCorrelated, double &ChiSq, int &dof, const std::string &OutBaseName,
              const std::string &ModelSuffix, Common::SeedType Seed );
  virtual const std::string &Type() const = 0;
};

#endif // Fitter_hpp
