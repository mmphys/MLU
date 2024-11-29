/**
 
 Base class for fitter
 
 Source file: Fitter.hpp

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

#ifndef Fitter_hpp
#define Fitter_hpp

#include "Covar.hpp"
#include "Model.hpp"

// Forward declarations to minimise coupling
struct FitterThread;       // OpenMP thread that will perform fit on each replica
struct Fitter;

// Callbacks to the Fit Controller
struct FitController
{
  virtual ~FitController() {}
  /// Models have agreed parameters - controller can adjust the parameters before offsets assigned
  virtual void ParamsAdjust( MLU::Params &mp, const Fitter &f ) {}
  /// Offsets have been assigned - controller can save indices
  virtual void SaveParameters( MLU::Params &mp, const Fitter &f ) {}
  /// FitController gets to participate in guesses
  virtual void Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                      const VectorView &FitData, bool bLastChance ) const {};
  /// Called once at the start of each replica, immediately after constants loaded. Compute per-replica derived values
  virtual void SetReplica( Vector &ModelParams ) const {}
  /// Called every time the fitter makes a new guess on each replica - allow controller to compute derived properties
  virtual void ComputeDerived( Vector &ModelParams ) const {}
  /// Allow the fit controller to adjust the list of parameters in the parameter covariance matrix
  virtual void ParamCovarList( MLU::Params &paramsCovar ) const {}
  static FitController None;
};

struct Fitter
{
  // Simple command-line options
  const bool bAnalyticDerivatives;
  const bool bTestRun;
  const bool bCentralGuess; // Use fit on central replica as guess for other replicas
  const bool bOverwrite;
  const double HotellingCutoff;
  const double ChiSqDofCutoff;
  const double RelEnergySep;
  const double ShrinkFactor;
  const int MinDof;
  const int Retry;
  const int MaxIt;
  const double Tolerance;
  const int SummaryLevel;
  const bool bSaveCMat;
  const int Verbosity;
  const bool UserGuess;
  const int Strictness;
  const scalar MonotonicUpperLimit;
  const int ErrorDigits;
  // More complex command-line options
  DataSet &ds;
  const bool bFitCorr; // true if we are fitting correlators, false = fitting models
  FitController &fitController;
  const int NumFiles; // Number of correlator files in the dataset
  //const std::vector<std::string> &OpNames;
  const std::vector<Model::Args> ModelArgs;
  //const int NumOps;
  std::vector<ModelPtr> model;      // Model for each correlator
  const int NumExponents;
  std::vector<DataSet::FixedParam> ParamFixed; // Map from constants in DataSet to parameters
  const Params mp; // Model Parameters
  const bool bAllParamsKnown; // Do I already know everything? No fit - just produce theory - data
  const CovarParams cp;
  Vector Guess; // This is the guess all replicas should use
  ModelFile OutputModel;

  // These variables set once at the start of each fit
  int dof;
  bool bCorrelated;

  explicit Fitter( Model::CreateParams &mcp, DataSet &ds_, std::vector<Model::Args> &&ModelArgs,
                   CovarParams &&cp, bool bFitCorr, FitController &fitController );
  virtual ~Fitter() {}

protected:
  // Used during construction (so that we can make the results const)
  std::vector<ModelPtr> CreateModels( Model::CreateParams &mcp, std::vector<Model::Args> ModelArgs );
  int GetNumExponents();
  Params MakeModelParams();
  void MakeGuess();
  /// Save the parameter correlation matrix (at the end of the fit)
  void SaveParamCorrel( const std::string &sFileName );
  virtual FitterThread * MakeThread( bool bCorrelated, ModelFile &OutputModel ) = 0;

public:
  static Fitter * Make( Model::CreateParams &&mcp, DataSet &ds,
                        std::vector<Model::Args> &&ModelArgs, CovarParams &&cp, bool bFitCorr,
                        FitController &fitController = FitController::None );
  static void SayNumThreads( std::ostream &os );
  bool Dump( int idx ) const;
  void Dump( int idx, const std::string &Name, const Matrix &m ) const;
  void Dump( int idx, const std::string &Name, const Vector &v ) const;
  void SaveMatrixFile( const Matrix &m, const std::string &Type, const std::string &Filename,
                         const char *pGnuplotExtra = nullptr ) const;
  bool PerformFit( bool bCorrelated, double &ChiSq, int &dof, const std::string &OutBaseName,
              const std::string &ModelSuffix );
  void Show( Param::Type type ) const;
  virtual const std::string &Type() const = 0;
  std::vector<std::string> GetModelTypes() const;
  std::vector<std::string> GetModelArgs() const;
  void WriteSummaryTD( const std::string &sOutFileName, bool bVerboseSummary );
};

#endif // Fitter_hpp
