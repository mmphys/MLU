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

#include <MLU/Common.hpp>

// Uncomment the next line if your cmath doesn't define M_PI etc by default
//#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <ios>
#include <set>
//#include <sys/stat.h>

using scalar = double;
using Matrix = Common::Matrix<scalar>;
using Vector = Common::Vector<scalar>;
using MatrixView = Common::MatrixView<scalar>;
using VectorView = Common::VectorView<scalar>;

using Fold = Common::Fold<scalar>;
using vCorrelator = std::vector<Fold>;
using ModelFile = Common::Model<scalar>;
using DataSet = Common::DataSet<scalar>;
using vString = std::vector<std::string>;
using vInt = std::vector<int>;
using UniqueNames = Common::UniqueNames;

constexpr int idxSrc{ 0 };
constexpr int idxSnk{ 1 };
extern const char * pSrcSnk[];

// This is the name of an energy level
extern const std::string E;

enum class ModelType{ Unknown, Exp, Cosh, Sinh, ThreePoint, Constant };
std::ostream & operator<<( std::ostream &os, const ModelType m );
std::istream & operator>>( std::istream &is, ModelType &m );

// Default parameters for model creation
struct ModelDefaultParams
{
  const bool bForceSrcSnkDifferent;
  const int NumExponents;
  //int NumOps;
  ModelDefaultParams( bool bForceSrcSnkDifferent_, const Common::CommandLine &cl )
  : bForceSrcSnkDifferent{ bForceSrcSnkDifferent_ }, NumExponents{ cl.SwitchValue<int>( "e" ) } {}
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
  virtual void ModelParamsChanged( Vector &ScratchPad, const Vector &ModelParams ) const {};
  // This is where the actual computation is performed
  virtual scalar operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const = 0;
  // Partial derivative of the model function on a timeslice with respect to each parameter
  virtual double Derivative( int t, int p ) const { return 0; }
};

struct Parameters
{
  struct Parameter
  {
    std::string Name;
    scalar      Value;
    scalar      Error;
    Parameter() {} // Random Value and Error is ok - would like this to show up
    Parameter( const std::string &name_, scalar value_, scalar error_ ) : Name{name_}, Value{value_}, Error{error_} {}
  };
  using iterator = std::vector<Parameter>::iterator;
  using const_iterator = std::vector<Parameter>::const_iterator;
protected:
  std::size_t maxLen = 0;
  std::vector<Parameter> Params;
public:
  // Helpers
  inline void clear() { maxLen = 0; Params.clear(); }
  inline std::size_t size() const { return Params.size(); }
  inline std::size_t MaxLen() const { return maxLen; }
  //inline       Parameter & operator[]( std::size_t i )       { return Params[i]; }
  inline const Parameter & operator[]( std::size_t i ) const { return Params[i]; }
  inline void Add( const std::string &Name, scalar Value, scalar Error )
  {
    Params.emplace_back( Name, Value, Error );
    std::size_t Len = Name.length();
    if( maxLen < Len )
      maxLen = Len;
  }
  inline iterator begin() { return Params.begin(); }
  inline iterator end()   { return Params.end(); }
  inline const_iterator begin() const { return Params.begin(); }
  inline const_iterator end()   const { return Params.end(); }
};

std::ostream & operator<<( std::ostream &os, const Parameters &Params );

struct ParamState
{
  Parameters parameters;
  bool bValid;
  scalar TestStat;
  unsigned int NumCalls;
        scalar getTestStat() const { return bValid ? TestStat : 0; }
  unsigned int getNumCalls() const { return bValid ? NumCalls : 0; }
  ParamState( const Parameters &parameters_, bool bValid_=false, scalar TestStat_=0, unsigned int NumCalls_=0 )
  : parameters{parameters_}, bValid{bValid_}, TestStat{TestStat_}, NumCalls{NumCalls_} {}
  ParamState( Parameters &&parameters_, bool bValid_=false, scalar TestStat_=0, unsigned int NumCalls_=0 )
  : parameters{std::move( parameters_ )}, bValid{bValid_}, TestStat{TestStat_}, NumCalls{NumCalls_} {}
  virtual ~ParamState(){}
  virtual void StandardOut( std::ostream &os ) const = 0; //TODO: Deuglify
  virtual void ReplicaMessage( std::ostream &os ) const = 0; //TODO: Deuglify
};

std::ostream & operator<<( std::ostream &os, const ParamState &State );

// Fitter parameters
struct CovarParams
{
  const DataSet &ds;
  // Whether to rebin the data before computing covariances
  std::vector<int> RebinSize;
  // Where am I building the covariance from
  using SS = Common::SampleSource;
  SS Source;
  // 0 = build correlation from CovarSource, then scale by the variance of the data on each bootstrap replica
  // else this is the number of bootstrap replicas to use when estimating covariance
  int CovarNumBoot;
  // Random numbers to get the covariance from a bootstrap
  using fint = Common::fint;
  Common::MatrixView<fint> CovarRandom;
  std::vector<fint> vCovarRandom; // Contains the random numbers iff I generated them
  // Freeze covariance to the central replica
  bool bFreeze;
  // Manually loaded correlation matrix
  Matrix Covar;
  // Constructor
  CovarParams( const Common::CommandLine &cl, DataSet &ds );
  // Can we compute (co)variance on non-central replicas - i.e. do we have rebinned data
  inline bool SupportsUnfrozen() const { return ds.corr[0].NumSamplesBinned(); }
  // Source is a bootstrap - only valid to compute variance on central replica
  inline bool SourceIsBootstrap() const
  { return Source == SS::Bootstrap || ( Source == SS::Raw && ds.corr[0].bRawBootstrap ); }
  // This is the appropriate parameter for m to use in the T^2 distribution
  inline int CovarSampleSize() const
  {
    return SourceIsBootstrap() ? ds.corr[0].SampleSize
                               : ds.corr[0].NumSamples( Source==SS::Raw ? SS::Raw : SS::Binned );
  }
  // Make correlation matrix, then scale it by variance of data on each replica
  inline bool OneStep( int Replica ) const
  {
    // TODO: support covar from replica data
    return CovarNumBoot || (Replica == Fold::idxCentral && (Source == SS::Binned || Source == SS::Bootstrap) );
  }
  // How many covariance samples are there
  inline int CovarCount() const { return ds.corr[0].NumSamples( Source ); }
};

std::ostream & operator<<( std::ostream &os, const CovarParams &cp );

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
  virtual const std::string &Type() const = 0;
  explicit Fitter( const Common::CommandLine &cl, const DataSet &ds_,
                   const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames_,
                   CovarParams &&cp );
  virtual ~Fitter() {}
  inline bool Dump( int idx ) const { return Verbosity > 2 || ( Verbosity >= 1 && idx == Fold::idxCentral ); }
  inline void Dump( int idx, const std::string &Name, const Matrix &m ) const
  {
    if( Verbosity > 2 || ( Verbosity >= 1 && idx == Fold::idxCentral ) )
      std::cout << Name << Common::Space << m << Common::NewLine;
  }
  void Dump( int idx, const std::string &Name, const Vector &v ) const
  {
    if( Verbosity > 2 || ( Verbosity >= 1 && idx == Fold::idxCentral ) )
      std::cout << Name << Common::Space << v << Common::NewLine;
  }
  std::vector<Common::ValWithEr<scalar>>
  PerformFit( bool bCorrelated, double &ChiSq, int &dof, const std::string &OutBaseName,
              const std::string &ModelSuffix, Common::SeedType Seed );
  void SaveMatrixFile( const Matrix &m, const std::string &Type, const std::string &Filename,
                         const char *pGnuplotExtra = nullptr ) const;
};

#endif // MultiFit_hpp
