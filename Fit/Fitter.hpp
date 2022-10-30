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
#include "Model.hpp"

// Forward declarations to minimise coupling
class FitterThread;       // OpenMP thread that will perform fit on each replica

class Fitter
{
public:
  // Simple command-line options
  const bool bAnalyticDerivatives;
  const bool bTestRun;
  const bool bOverwrite;
  const double HotellingCutoff;
  const double ChiSqDofCutoff;
  const double RelEnergySep;
  const int MinDof;
  const int Retry;
  const int MaxIt;
  const double Tolerance;
  const int SummaryLevel;
  const bool bSaveCMat;
  const int Verbosity;
  const bool UserGuess;
  // More complex command-line options
  const DataSet &ds;
  const int NumFiles; // Number of correlator files in the dataset
  const std::vector<std::string> &OpNames;
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

  explicit Fitter( const Common::CommandLine &cl, const DataSet &ds_,
                   std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames_,
                   CovarParams &&cp );
  virtual ~Fitter() {}

protected:
  // Used during construction (so that we can make the results const)
  std::vector<ModelPtr> CreateModels(const Common::CommandLine &cl, std::vector<Model::Args> &ModelArgs);
  int GetNumExponents();
  Params MakeModelParams();
  void MakeGuess();
  virtual FitterThread * MakeThread( bool bCorrelated, ModelFile &OutputModel ) = 0;

public:
  static void SayNumThreads( std::ostream &os );
  bool Dump( int idx ) const;
  void Dump( int idx, const std::string &Name, const Matrix &m ) const;
  void Dump( int idx, const std::string &Name, const Vector &v ) const;
  void SaveMatrixFile( const Matrix &m, const std::string &Type, const std::string &Filename,
                         const char *pGnuplotExtra = nullptr ) const;
  void PerformFit( bool bCorrelated, double &ChiSq, int &dof, const std::string &OutBaseName,
              const std::string &ModelSuffix, Common::SeedType Seed );
  virtual const std::string &Type() const = 0;
};

#endif // Fitter_hpp
