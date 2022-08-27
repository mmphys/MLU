/*************************************************************************************
 
 OpenMP thread that will perform fit on each replica
 
 Source file: FitterThread.cpp

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

#ifndef FitterThread_hpp
#define FitterThread_hpp

#include "MultiFit.hpp"

// Forward declaration of fitter for multi-exponential fit
class Fitter;
class Parameters;
class ParamState;

// Several of these will be running at the same time on different threads during a fit
class FitterThread
{
public:
  const Fitter &parent;
  const int Extent;
  static constexpr int FirstIndex{ std::numeric_limits<int>::lowest() };
  // These variables change on each iteration
  int idx;
  Matrix Covar; // Covariance matrix used in the fit
  Matrix Correl; // Correlation matrix
  Vector StdErrorMean; // sqrt of diagonals of Covar
  // GSL:     Cholesky decomposition of InverseCorrelationMatrix
  // Minuit2: Cholesky decomposition of CorrelationMatrix
  Matrix Cholesky;
  Vector CholeskyDiag; // 1 / StdErrorMean
  VectorView Data;        // Data points we are fitting
  Vector Theory;  // Theory prediction of our model on this replica
  mutable Vector Error;   // (Theory - Data) * CholeskyScale (mutable because of Minuit2)
  // These next are for model parameters
  Vector ModelParams;
  std::vector<std::pair<double,int>> SortingHat;
protected:
  bool bCorrelated;
  ModelFile &OutputModel; // Fill this with parameters for each replica
  vCorrelator &CorrSynthetic; // Fill this with data for model correlator
  std::vector<Vector> ModelBuffer;
  // Helper functions
  virtual ParamState * MakeParamState( const Parameters &Params ) = 0;
  void SaveStdError();
public:
  FitterThread( const Fitter &fitter, bool bCorrelated, ModelFile &OutputModel, vCorrelator &CorrSynthetic );
  virtual ~FitterThread() {}
  // Switch to this index
  void SetReplica( int idx, bool bShowOutput = false, bool bSaveMatrices = false,
                   const std::string *pBaseName = nullptr );
  inline std::string ReplicaString( int iFitNum ) const
  {
    std::stringstream ss;
    ss << ( bCorrelated ? "C" : "Unc" ) << "orrelated fit " << iFitNum << " on ";
    if( idx == Fold::idxCentral )
      ss << "central replica";
    else
      ss << "replica " << idx;
    return ss.str();
  }
  void ReplicaMessage( const ParamState &state, int iFitNum ) const;
  bool SaveError( Vector &Error, const scalar * FitterParams, std::size_t Size, std::size_t Stride = 1 );
  bool AnalyticJacobian( Matrix &Jacobian ) const;
  scalar RepeatFit( ParamState &Guess, int MaxGuesses );
  void UpdateGuess( Parameters &parGuess );
  scalar FitOne( const Parameters &parGuess );
  // Implement this to support a new type of fitter
  virtual void Minimise( ParamState &Guess, int iNumGuesses ) = 0;
  virtual bool CholeskyAdjust() { return false; } // true to indicate Cholesky matrix needs inversion
  virtual int NumRetriesGuess() const = 0;
  virtual int NumRetriesFit() const = 0;
  virtual std::string Description() const { return std::string(); }
};

#endif // FitterThread_hpp
