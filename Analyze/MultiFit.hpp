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
#include <set>
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
using vCorrelator = std::vector<Fold>;
using ModelFile = Common::Model<scalar>;
//using Matrix = Eigen::MatrixXd; // dynamic sized matrix of complex double
using Matrix = Common::Matrix<scalar>;
using Vector = Common::Vector<scalar>;
using vString = std::vector<std::string>;
using vInt = std::vector<int>;
// A list of case insensitive, but unique names, each mapped to an int
using UniqueNames = std::map<std::string, int, Common::LessCaseInsensitive>;

constexpr int idxSrc{ 0 };
constexpr int idxSnk{ 1 };
extern const char * pSrcSnk[];

// This is the name of an energy level
extern const std::string E;

enum class FitterType{ Minuit2, GSL };
std::ostream & operator<<( std::ostream &os, const FitterType f );

enum class ModelType{ Unknown, Exp, Cosh, Sinh, ThreePoint, Constant };
std::ostream & operator<<( std::ostream &os, const ModelType m );
std::istream & operator>>( std::istream &is, ModelType &m );

// This is the per-thread instance of each model
/*struct ThreadModel
{
public:
  const int Nt;
  const int HalfNt;
  const int NumExponents;
  //protected:
  Vector Params;
public:
  ThreadModel( int Nt_, int HalfNt_, int NumExponents_ ) : Nt{Nt_}, HalfNt{HalfNt_}, NumExponents{NumExponents_} {}
  virtual void Init( const scalar * ModelParameters, std::size_t NumParams, std::size_t Stride = 1 );
  inline void Init( const Vector &x ) { Init( x.data, x.size, x.stride ); }
  virtual double operator()( int t ) const = 0;
  virtual double Derivative( int t, int p ) const { return 0; };
  virtual void ParamsChanged() = 0;
};*/

// Default parameters for model creation
struct ModelDefaultParams
{
  int NumExponents;
  int NumOps;
  bool bFactor;
};

// This represents the model I'm fitting to
class Model;
using ModelPtr = std::unique_ptr<Model>;

class ThreadModel
{
  const Model &model;
  friend class Model;
protected:
  explicit ThreadModel( const Model &Model_ ) : model{ Model_ } {}
};

class Model
{
  friend class ModelSet;
public:
  int Nt;
  int HalfNt;
  int NumExponents;
  vString ParamNames;
  vInt    ParamIdx;
  vString ParamNamesPerExp;
  vInt    ParamIdxPerExp;
  //protected:
  Vector Params;
protected:
  inline void SetParams( int Nt_, int HalfNt_, int NumExponents_ ) { Nt=Nt_; HalfNt=HalfNt_; NumExponents=NumExponents_; }
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName )
  { ParamNamesPerExp.push_back( E ); } // Make sure to call base if you override (because everything uses Energy)
public:
  virtual ~Model() {}
  // Save a per-thread buffer for the model to use during a fit
  virtual ThreadModel * MakeThreadModel() const { throw std::runtime_error( "Implement MakeThreadModel()" ); };
  // This is how to make a model
  static ModelPtr MakeModel( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  // Called while validating models. This should really be cleaned up
  virtual std::size_t UnknownParameterCount( UniqueNames &Names ) const { return 0; }
  // There are more unknowns than correlators - combine overlap coefficients into single operator
  virtual void ReduceUnknown( const UniqueNames &Names ) {}
};

struct DataSet
{
  struct ConstantSource
  {
    std::size_t File;
    std::size_t idx;
    ConstantSource( std::size_t File_, std::size_t idx_ ) : File{File_}, idx{idx_} {}
  };
  int NSamples; // Number of samples we are using. These are guaranteed to exist
  int Extent = 0;   // Number of data points in our fit (i.e. total number of elements in FitTimes)
  int MinExponents = 0;
  int MaxExponents = 0;
  std::vector<Fold>             corr;     // Correlator files
  std::vector<std::vector<int>> FitTimes; // The actual timeslices we are fitting to in each correlator
  UniqueNames ConstantNames;
  UniqueNames ConstantNamesPerExp;
  using ConstMap = std::map<std::string, ConstantSource, Common::LessCaseInsensitive>;
  ConstMap constMap;
protected:
  std::vector<ModelFile>        constFile;// Each of the constant files (i.e. results from previous fits) I've loaded
  void AddConstant( const std::string &Name, std::size_t File, std::size_t idx );
  void AddConstant( const std::string &Name, std::size_t File, std::size_t idx, int e );
public:
  explicit DataSet( int nSamples = 0 ) : NSamples{ nSamples } {}
  inline bool empty() const { return corr.empty() && constFile.empty(); }
  void clear();
  void LoadFile( const std::string &sFileName, std::vector<std::string> &OpNames, std::vector<std::string> &ModelArgs,
                 const std::string &Args );
  void SetFitTimes( const std::vector<std::vector<int>> &FitTimes );
  void SetFitTimes( int tMin, int tMax );
  void GetData( int idx, Vector &vResult ) const;
  void MakeInvErr( int idx, Vector &Var ) const;
  void MakeCovariance( int idx, Matrix &Covar ) const;
  //void MakePredictions( Vector &Prediction ) const;
  void SaveCovariance( const std::string FileName, const Matrix &Covar, const std::vector<ModelPtr> *pModels=nullptr ) const;
};

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
  std::size_t nevalf;
  std::size_t nevaldf;
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
  vCorrelator &CorrSynthetic; // Fill this with data for model correlator
  mutable std::vector<std::unique_ptr<ThreadModel>> model; // mutable only because Minuit2 operator() is const
  // Helper functions
public:
  FitterThread( const Fitter &fitter, bool bCorrelated, ModelFile &ModelParams, vCorrelator &CorrSynthetic );
  virtual ~FitterThread() {}
  void SetReplica( int idx, bool bShowOutput ); // Switch to this index
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
  bool SaveError( Vector &ModelError ) const;
  bool AnalyticJacobian( Matrix &Jacobian ) const;
  scalar RepeatFit( ParamState &Guess, int MaxGuesses );
  void UpdateGuess( Parameters &parGuess );
  scalar FitOne( const Parameters &parGuess, const std::string &SaveCorMatFileName );
  // Implement this to support a new type of fitter
  virtual void Minimise( ParamState &Guess, int iNumGuesses ) = 0;
  virtual void MakeCovarCorrelated() {} // Override to adjust covariance matrix
  virtual int NumRetriesGuess() const = 0;
  virtual int NumRetriesFit() const = 0;
  virtual std::string Description() const { return std::string(); }
};

class Fitter
{
protected:
  struct FixedParam
  {
    int                     idx;
    DataSet::ConstantSource src;
    FixedParam( int idx_, const DataSet::ConstantSource &src_ ) : idx{idx_}, src{src_} {}
  };
public:
  // These variables set in the constructor
  const FitterType fitType;
  DataSet &ds;  // Non-const because I need to set fit ranges
  const bool bAnalyticDerivatives;
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
  const int NumExponents;
  const std::vector<std::string> ParamNames;
  const int NumModelParams; // Total number of parameters across all models (fixed and variable)
  const int NumFixed;       // Number of fixed parameters
  const int NumVariable;    // Number of variable parameters

  // These variables set once at the start of each fit
  int dof;
  bool bCorrelated;

protected:
  std::vector<ModelPtr> model;      // Model for each correlator
  std::vector<FixedParam> ParamFixed; // Map from constants to parameters
  std::vector<int> ParamVariable;     // Map from fitting engine parameters to parameters

protected:
  // Used during construction (so that we can make the results const)
  int CreateModels( const std::vector<std::string> &ModelArgs, const ModelDefaultParams &modelDefault );
  std::size_t EnsureModelsSolubleHelper( UniqueNames &Names, std::size_t &NumWithUnknowns );
  std::vector<std::string> MakeParamNames();
  int FixedVsVariableParams();
  // Helper functions
  //inline int AlphaIndex( int snk, int src) const { return snk * NumOps + src; };
public:
  inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; }; // TODO: Delete

public:
  explicit Fitter( FitterType fitType, const DataSet &ds_, const std::vector<std::string> &ModelArgs,
                   const ModelDefaultParams &modelDefault, const std::vector<std::string> &opNames_,
                   const std::string &sOpNameConcat, int Verbosity, const std::string &OutputBaseName,
                   Common::SeedType Seed, bool bFreezeCovar, bool bSaveCorr, bool bSaveCMat,
                   int Retry, int MaxIt, double Tolerance, double RelEnergySep, bool bNumericDerivatives );
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
  FitterThreadMinuit2( const Fitter &fitter_, bool bCorrelated_, ModelFile &modelParams_, vCorrelator &CorrSynthetic_ )
  : FitterThread( fitter_, bCorrelated_, modelParams_, CorrSynthetic_ ) {}
  virtual ~FitterThreadMinuit2() {}
  // These are part of the FCNBase interface
  virtual double Up() const { return 1.; }
  virtual double operator()( const std::vector<double> &ModelParameters ) const;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  virtual void Minimise( ParamState &Guess, int iNumGuesses );
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
  FitterThreadGSL( const Fitter &Fitter, bool bCorrelated, ModelFile &ModelParams, vCorrelator &CorrSynthetic );
  virtual ~FitterThreadGSL();
  virtual void Minimise( ParamState &Guess, int iNumGuesses );
  virtual void MakeCovarCorrelated();
  virtual int NumRetriesGuess() const { return parent.Retry; };
  virtual int NumRetriesFit() const { return parent.Retry; };
  virtual std::string Description() const;
};

#endif // MultiFit_hpp
