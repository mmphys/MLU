/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: MultiFit.hpp
 
 Copyright (C) 2019
 
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
 *************************************************************************************/
/*  END LEGAL */

#ifndef MultiFit_hpp
#define MultiFit_hpp

//#include "CommonLatAn.hpp"
#include "Common.hpp"

// Uncomment the next line if your cmath doesn't define M_PI etc by default
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <ios>
//#include <sys/stat.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/VariableMetricMinimizer.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MnPrint.h>

#include <gsl/gsl_multifit_nlinear.h>

using scalar = double;
using Fold = Common::Fold<scalar>;
using Correlators = std::vector<Fold>;
using ModelFile = Common::Model<scalar>;
//using Matrix = Eigen::MatrixXd; // dynamic sized matrix of complex double
using Matrix = Common::Matrix<scalar>;
using Vector = Common::Vector<scalar>;

enum class FitterType{ Minuit2, GSL };
inline std::ostream & operator<<( std::ostream &os, const FitterType f );

class Parameters
{
public:
  struct Parameter
  {
    std::string Name;
    scalar      Value;
    scalar      Error;
    Parameter() {} // Random Value and Error is ok - would like this to show up
    Parameter( const std::string &name_, scalar value_, scalar error_ ) : Name{name_}, Value{value_}, Error{error_} {}
  };
  std::size_t MaxLen;
  std::vector<Parameter>  Params;
protected:
  inline void CheckLength( const std::string &name_ )
  {
    std::size_t Len = name_.length();
    if( MaxLen < Len )
      MaxLen = Len;
  }
public:
  // constructors
  explicit Parameters( std::size_t Num = 0 ) : MaxLen{ 0 }, Params( Num ){};
  Parameters( const Parameters &p ) : Params( p.Params ) { MaxLen = p.MaxLen; }
  Parameters( const ROOT::Minuit2::MnUserParameters &Minuit2Params ) : MaxLen{ 0 } { *this = Minuit2Params; }
  Parameters( const ROOT::Minuit2::MnUserParameterState &Minuit2State ) : Parameters( Minuit2State.Parameters() ) {}
  inline Parameters& operator=( const ROOT::Minuit2::MnUserParameters &Minuit2Params )
  {
    Params.clear();
    MaxLen = 0;
    const std::vector<ROOT::Minuit2::MinuitParameter> & params{ Minuit2Params.Parameters() };
    Params.reserve( params.size() );
    for( const ROOT::Minuit2::MinuitParameter &p : params )
      Add( p.GetName(), p.Value(), p.Error() );
    return *this;
  }
  inline Parameters& operator=( const ROOT::Minuit2::MnUserParameterState &State ) { return *this=State.Parameters(); };
  inline Parameter & operator[]( std::size_t i ) { return Params[i]; }
  inline void Add( const std::string &name_, scalar value_, scalar error_ )
  {
    Params.emplace_back( name_, value_, error_ );
    CheckLength( name_ );
  }
  inline scalar Value( std::size_t i ) const { return Params[i].Value; }
};

std::ostream & operator<<( std::ostream &os, const Parameters &Params );

struct GSLState
{
  int ConvergeReason;
};

struct ParamState
{
  bool bValid;
  Parameters parameters;
  scalar TestStat;
  unsigned int NumCalls;
  scalar edm;
  bool bGotMinuit2State;
  ROOT::Minuit2::MnUserParameterState Minuit2State;
  bool bGotGSLState;
  GSLState gslState;
  ParamState( Parameters parameters_ ) : bValid{false}, bGotMinuit2State{false}, bGotGSLState{false}, parameters( parameters_ ) {}
  ParamState( const ROOT::Minuit2::MnUserParameterState &State )
  : bValid{State.IsValid()}, parameters( State ), TestStat{State.Fval()}, NumCalls{State.NFcn()}, edm{State.Edm()},
    bGotMinuit2State{true}, Minuit2State{State} {}
  bool IsValid() const { return bValid; }
  inline ParamState& operator=( const ROOT::Minuit2::MnUserParameterState &State )
  {
    bValid = State.IsValid();
    if( bValid )
    {
      parameters = State.Parameters();
      TestStat = State.Fval();
      NumCalls = State.NFcn();
      edm = State.Edm();
    }
    bGotMinuit2State = true;
    Minuit2State = State;
    return *this;
  };
  const Parameters &Parameters() const { return parameters; }
  scalar Fval() const { return TestStat; }
  unsigned int NFcn() const { return NumCalls; }
  scalar Edm() const { return edm; }
  void Add( const std::string &name_, scalar value_, scalar error_ ) { parameters.Add( name_, value_, error_ ); }
};

std::ostream & operator<<( std::ostream &os, const ParamState &State );

// Forward declaration of fitter for multi-exponential fit
class Fitter;

// This represents the model I'm fitting to
// When I introduce multiple models, this needs to split into base and derived classes
// NB: One of these will be instantiated per ThreadFitter (for thread safety)
class Model
{
public:
  const int Nt;
  const int HalfNt;
  const int NumFiles;
  const int NumExponents;
  const int NumOps;
protected:
  std::vector<std::vector<int>> snk;
  std::vector<std::vector<int>> src;
  std::vector<Common::Parity> parity;
  std::vector<double> SinhCoshAdjust;
  Vector Energy;
  Vector Coeff;
public:
  Model( const Fitter &parent );
  void Init( const scalar * ModelParameters, std::size_t NumParams, std::size_t Stride = 1 );
  inline void Init( const Vector &x ) { Init( x.data, x.size, x.stride ); }
  inline double operator()( int f, int t ) const;
  inline double Derivative( int f, int t, int p ) const;
};

// Several of these will be running at the same time on different threads during a fit
class FitterThread
{
public:
  const Fitter &parent;
  // These variables change on each iteration
  int idx;
  Matrix Covar;
  Matrix Cholesky;
  Vector CholeskyDiag;
  mutable Vector Error;
  std::vector<std::vector<double>> SortingHat;
protected:
  bool bCorrelated;
  ModelFile &ModelParams; // Fill this with parameters for each replica
  Correlators &CorrSynthetic; // Fill this with data for model correlator
  mutable Model model; // mutable only because Minuit2 operator() is const
  // Helper functions
public:
  FitterThread( const Fitter &fitter, bool bCorrelated, ModelFile &ModelParams, Correlators &CorrSynthetic );
  virtual ~FitterThread() {}
  void MakeCovar( int idx, bool bShowOutput = false ); // Switch to this index
  inline std::string ReplicaString( int iFitNum ) const
  {
    std::string sError{ bCorrelated ? "C" : "Unc" };
    sError.append( "orrelated fit " + std::to_string( iFitNum ) + " on " );
    if( idx == Fold::idxCentral )
      sError.append( "central replica" );
    else
      sError.append( "replica " + std::to_string( idx ) );
    return sError;
  }
  void ReplicaMessage( const ParamState &state, int iFitNum ) const;
  bool SaveError( Vector &ModelError ) const;
  bool SaveJacobian( Matrix &Jacobian ) const;
  scalar RepeatFit( ParamState &Guess, int MaxGuesses );
  void UpdateGuess( Parameters &parGuess );
  scalar FitOne( const Parameters &parGuess, const std::string &SaveCorMatFileName );
  // Implement this to support a new type of fitter
  virtual void Minimise( ParamState &Guess, int iNumGuesses ) = 0;
  virtual void MakeCovarCorrelated() = 0;
  virtual int NumRetriesGuess() const = 0;
  virtual int NumRetriesFit() const = 0;
};

class Fitter
{
public:
  const FitterType fitType;
  //static constexpr double pi{ M_PI };
  const int NumOps;
  const std::vector<std::string> &OpNames;
  const bool bFactor;
  const std::string &sOpNameConcat;
  const int Verbosity;
  const std::string OutputBaseName;
  const Common::SeedType Seed;
  const bool bFreezeCovar;
  const bool bSaveCorr;
  const bool bSaveCMat;
  const int Retry;
  const int MaxIt;
  const double Tolerance;
  const double RelEnergySep;
  const int NumFiles;
  const int NSamples;
  const int Nt;
  const int NumExponents;
  const int NumParams;
  const Correlators Corr;
  const std::vector<std::string> ParamNames;

  // These variables set once at the start of each fit
  int tMin = -1;
  int tMax = -1;
  int NtCorr;
  int Extent;
  int dof;
  bool bCorrelated;

  // Helper functions
  //inline int AlphaIndex( int snk, int src) const { return snk * NumOps + src; };
  std::vector<std::string> MakeParamNames() const;
  inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; };

public:
  explicit Fitter( FitterType fitType, Correlators &&Corr, const std::vector<std::string> &OpNames, bool bFactor,
                   const std::string &sOpNameConcat, int NumExponents, int Verbosity, const std::string &OutputBaseName,
                   Common::SeedType Seed, bool bFreezeCovar, bool bSaveCorr, bool bSaveCMat,
                   int Retry, int MaxIt, double Tolerance, double RelEnergySep, int NSamples );
  virtual ~Fitter() {}
  std::vector<Common::ValWithEr<scalar>>
  PerformFit(bool bCorrelated, int tMin, int tMax, double &ChiSq, int &dof );
};

// Several of these will be running at the same time on different threads during a fit
class FitterThreadMinuit2 : public FitterThread, ROOT::Minuit2::FCNBase
{
protected:
  static constexpr unsigned int StrategyLevel{ 1 }; // for parameter ERRORS (see MnStrategy) 0=low, 1=medium, 2=high
  static const ROOT::Minuit2::MnStrategy Strategy;
  ROOT::Minuit2::VariableMetricMinimizer Minimiser;
  // Helper functions
public:
  FitterThreadMinuit2( const Fitter &fitter_, bool bCorrelated_, ModelFile &modelParams_, Correlators &CorrSynthetic_ )
  : FitterThread( fitter_, bCorrelated_, modelParams_, CorrSynthetic_ ) {}
  virtual ~FitterThreadMinuit2() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const std::vector<double> &ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  virtual void Minimise( ParamState &Guess, int iNumGuesses );
  virtual void MakeCovarCorrelated();
  virtual int NumRetriesGuess() const { return parent.Retry ? parent.Retry + 10 : 15; };
  virtual int NumRetriesFit() const { return parent.Retry ? parent.Retry : 5; };
};

// Several of these will be running at the same time on different threads during a fit
class FitterThreadGSL : public FitterThread
{
protected:
  Vector vGuess;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_workspace * ws;
  int  f( const Vector &x, Vector &f_ );
  int df( const Vector &x, Matrix &J );
  using Me = FitterThreadGSL;
  static int  sf( const gsl_vector * x, void *data, gsl_vector * f_ )
  { return reinterpret_cast<Me*>(data)->f( *reinterpret_cast<const Vector*>(x), *reinterpret_cast<Vector*>(f_) ); }
  static int sdf( const gsl_vector * x, void *data, gsl_matrix * J  )
  { return reinterpret_cast<Me*>(data)->df( *reinterpret_cast<const Vector*>(x), *reinterpret_cast<Matrix*>(J) ); }
public:
  FitterThreadGSL( const Fitter &Fitter, bool bCorrelated, ModelFile &ModelParams, Correlators &CorrSynthetic );
  virtual ~FitterThreadGSL();
  virtual void Minimise( ParamState &Guess, int iNumGuesses );
  virtual void MakeCovarCorrelated();
  virtual int NumRetriesGuess() const { return parent.Retry; };
  virtual int NumRetriesFit() const { return parent.Retry; };
};

#endif // MultiFit_hpp
