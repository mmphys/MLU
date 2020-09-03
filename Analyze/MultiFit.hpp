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

class Model
{
  friend class ModelSet;
public:
  int Nt;
  int HalfNt;
  int NumExponents;
  vString ParamNames;
  vString ParamNamesPerExp;                     // Be careful: these names EXCLUDE the energy
  vInt    ParamIdx;
  std::vector<std::vector<int>> ParamIdxPerExp; // Be careful: Zero'th element of this IS energy
  //protected:
  Vector Params;
protected:
  inline void SetParams( int Nt_, int HalfNt_, int NumExponents_ ) { Nt=Nt_; HalfNt=HalfNt_; NumExponents=NumExponents_; }
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName )
  { ParamNamesPerExp.push_back( E ); } // Make sure to call base if you override (because everything uses Energy)
public:
  virtual ~Model() {}
  // This is how to make a model
  static ModelPtr MakeModel( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  // NOT COUNTING ENERGY, put names of parameters that can be determined in list, return number remaining to determine
  // For now, synonymous with overlap coefficients
  virtual std::size_t UnknownParameterCount( UniqueNames &Names ) const { return 0; }
  // There are more unknowns than correlators - combine overlap coefficients into single operator
  virtual void ReduceUnknown( const UniqueNames &Names ) {}
  // Take a guess as to E0 (guesses will be averaged, and guesses for higher energies based on E0)
  virtual bool GuessE0( scalar &E0, const Vector &Corr ) const { return false; }
  // Take a guess for all parameters other than energies (guaranteed to be valid on entry). Return true for another pass
  virtual bool Guess( Vector &Parameters, std::vector<bool> &bKnown, int pass, const Vector &Corr ) const { return false; }
  // How big a scratchpad is required for each model?
  virtual int GetScratchPadSize() const { return 0; }
  // Cache values based solely on the model parameters (to speed up computation)
  virtual void UpdateScratchPad( Vector &ScratchPad, const Vector &ModelParams ) const {};
  // This is where the actual computation is performed
  virtual scalar operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const = 0;
};

struct DataSet
{
  struct ConstantSource
  {
    std::size_t File; // Index of the constant file in the DataSet this comes from
    std::size_t idx;  // Index of the parameter in this file
    ConstantSource( std::size_t File_, std::size_t idx_ ) : File{File_}, idx{idx_} {}
  };
  struct FixedParam
  {
    int             idx;  // Index of the parameter (relative to Fitter::ParamNames)
    ConstantSource  src;  // Where to get the value from
    FixedParam( int idx_, const ConstantSource &src_ ) : idx{idx_}, src{src_} {}
  };
  int NSamples;   // Number of samples we are using. These are guaranteed to exist
  int MaxSamples; // Maximum number of samples available - used for covariance. Guaranteed to exist. >= NSamples.
  int Extent = 0; // Number of data points in our fit (i.e. total number of elements in FitTimes)
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
  void SortOpNames( std::vector<std::string> &OpNames );
  void SetFitTimes( const std::vector<std::vector<int>> &FitTimes );
  void SetFitTimes( int tMin, int tMax );
  void GetData( int idx, Vector &vResult ) const;
  void GetFixed( int idx, Vector &vResult, const std::vector<FixedParam> &Params ) const;
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
  const int Extent;
  static constexpr int FirstIndex{ std::numeric_limits<int>::lowest() };
  // These variables change on each iteration
  int idx;
  Matrix Covar;
  Matrix Cholesky;
  Vector CholeskyDiag;
  Vector Data;            // Data points we are fitting
  Vector Error;           // Difference between model predictions (theory) and Data
  // These next are for model parameters
  Vector ModelParams;
  Vector SortingHat;
protected:
  bool bCorrelated;
  ModelFile &OutputModel; // Fill this with parameters for each replica
  vCorrelator &CorrSynthetic; // Fill this with data for model correlator
  std::vector<Vector> ModelBuffer;
  // Helper functions
public:
  FitterThread( const Fitter &fitter, bool bCorrelated, ModelFile &OutputModel, vCorrelator &CorrSynthetic );
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
  bool SaveError( Vector &Error, const scalar * FitterParams, std::size_t Size, std::size_t Stride = 1 );
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
public:
  // These variables set in the constructor
  const FitterType fitType;
  const DataSet &ds;  // Non-const because I need to set fit ranges
  const bool bAnalyticDerivatives;
  const int NumOps;
  const std::vector<std::string> &OpNames;
  const bool bFactor;
  //const std::string &sOpNameConcat;
  const int Verbosity;
  //const std::string OutputBaseName;
  const bool bFreezeCovar;
  const bool bSaveCorr;
  const bool bSaveCMat;
  const int Retry;
  const int MaxIt;
  const double Tolerance;
  const double RelEnergySep;
  const int NumFiles;
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

  // These variables set once at the start of each fit
  int dof;
  bool bCorrelated;

protected:
  // Used during construction (so that we can make the results const)
  std::vector<ModelPtr> CreateModels( const std::vector<std::string> &ModelArgs, const ModelDefaultParams &modelDefault );
  int GetNumExponents();
  std::size_t EnsureModelsSolubleHelper( UniqueNames &Names, std::size_t &NumWithUnknowns );
  std::vector<std::string> MakePerExpNames();
  std::vector<std::string> MakeParamNames();
  std::vector<DataSet::FixedParam> MakeParamFixed();
  std::vector<int> MakeParamVariable();
  // Helper functions
  //inline int AlphaIndex( int snk, int src) const { return snk * NumOps + src; };
public:
  //inline int MELIndex( int op, int EnergyLevel) const { return EnergyLevel * NumOps + op; }; // TODO: Delete

public:
  explicit Fitter( FitterType fitType, const DataSet &ds_, const std::vector<std::string> &ModelArgs,
                   const ModelDefaultParams &modelDefault, const std::vector<std::string> &opNames_,
                   int Verbosity, bool bFreezeCovar, bool bSaveCorr, bool bSaveCMat,
                   int Retry, int MaxIt, double Tolerance, double RelEnergySep, bool bNumericDerivatives );
  virtual ~Fitter() {}
  std::vector<Common::ValWithEr<scalar>>
  PerformFit( bool bCorrelated, double &ChiSq, int &dof, const std::string &OutBaseName, const std::string &ModelSuffix,
              Common::SeedType Seed );
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
  FitterThreadMinuit2( const Fitter &fitter_, bool bCorrelated_, ModelFile &OutputModel, vCorrelator &CorrSynthetic_ )
  : FitterThread( fitter_, bCorrelated_, OutputModel, CorrSynthetic_ ) {}
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
  int  f( const Vector &FitParams, Vector &Errors );
  int df( const Vector &x, Matrix &J );
  using Me = FitterThreadGSL;
  static int  sf( const gsl_vector * x, void *data, gsl_vector * f_ )
  { return reinterpret_cast<Me*>(data)->f( *reinterpret_cast<const Vector*>(x), *reinterpret_cast<Vector*>(f_) ); }
  static int sdf( const gsl_vector * x, void *data, gsl_matrix * J  )
  { return reinterpret_cast<Me*>(data)->df( *reinterpret_cast<const Vector*>(x), *reinterpret_cast<Matrix*>(J) ); }
public:
  FitterThreadGSL( const Fitter &Fitter, bool bCorrelated, ModelFile &OutputModel, vCorrelator &CorrSynthetic );
  virtual ~FitterThreadGSL();
  virtual void Minimise( ParamState &Guess, int iNumGuesses );
  virtual void MakeCovarCorrelated();
  virtual int NumRetriesGuess() const { return parent.Retry; };
  virtual int NumRetriesFit() const { return parent.Retry; };
  virtual std::string Description() const;
};

#endif // MultiFit_hpp
