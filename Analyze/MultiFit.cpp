/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: FitModel.cpp
 
 Copyright (C) 2020
 
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

#include "MultiFit.hpp"

//#define DEBUG_DISABLE_OMP

// Indices for operators in correlator names
const char * pSrcSnk[] = { "src", "snk" };
const std::string E{ "E" };

static const std::string sFitterTypeMinuit2{ "Minuit2" };
static const std::string sFitterTypeGSL{ "GSL" };

inline std::ostream & operator<<( std::ostream &os, const FitterType f )
{
  switch( f )
  {
    case FitterType::Minuit2:
      os << sFitterTypeMinuit2;
      break;
    case FitterType::GSL:
      os << sFitterTypeGSL;
      break;
    default:
      os.setstate( std::ios_base::failbit );
      break;
  }
  return os;
}

std::ostream & operator<<( std::ostream &os, const Parameters &Params )
{
  using Parameter = Parameters::Parameter;
  for( const Parameter &p : Params.Params )
    os << std::string( Params.MaxLen - p.Name.length() + 2, ' ' ) << p.Name
       << Common::Space << p.Value << "\t+/- " << p.Error << Common::NewLine;
  return os;
}

std::ostream & operator<<( std::ostream &os, const ParamState &State )
{
  if( State.bGotMinuit2State )
    os << State.Minuit2State;
  else
  {
    os << State.Parameters();
    if( State.bGotGSLState )
    {
      os << "";
    }
  }
  return os;
}

std::ostream & operator<<( std::ostream &os, const FitRange &fr )
{
  char colon{ ':' };
  os << fr.ti << colon << fr.tf;
  if( fr.dti != fr.dtf )
    os << colon << fr.dti << colon << fr.dtf;
  else if( fr.dti != 1 )
    os << colon << fr.dti;
  return os;
}

std::istream & operator>>( std::istream &is, FitRange &fr_ )
{
  int Numbers[4];
  int i{ 0 };
  for( bool bMore = true; bMore && i < 4 && is >> Numbers[i]; )
  {
    if( ++i < 4 && !is.ios_base::eof() && is.peek() == ':' )
      is.get();
    else
      bMore = false;
  }
  if( i < 2 )
    is.setstate( std::ios_base::failbit );
  else
  {
    FitRange fr( Numbers[0], Numbers[1] );
    if( i >= 3 )
    {
      fr.dti = Numbers[2];
      fr.dtf = ( i == 3 ? Numbers[2] : Numbers[3] );
    }
    fr_ = fr;
  }
  return is;
}

FitRangesIterator FitRanges::begin()
{
  FitRangesIterator it( *this );
  std::vector<FitRange> &fitRange( *this );
  for( std::size_t i = 0; i < fitRange.size(); ++i )
  {
    it[i].ti = fitRange[i].ti;
    it[i].tf = fitRange[i].tf;
  }
  return it;
}

FitRangesIterator FitRanges::end()
{
  FitRangesIterator it{ begin() };
  std::vector<FitRange> &fitRange( *this );
  it.back().tf = fitRange.back().tf + fitRange.back().dtf;
  return it;
}

FitRangesIterator &FitRangesIterator::operator++()
{
  // Don't increment if we are already at the end
  if( !AtEnd() )
  {
    std::vector<FitTime> &fitTime( *this );
    bool bOverflow{ true };
    for( std::size_t i = 0; bOverflow && i < fitTime.size(); ++i )
    {
      bOverflow = false;
      if( ++fitTime[i].ti >= Ranges[i].ti + Ranges[i].dti )
      {
        fitTime[i].ti = Ranges[i].ti;
        if( ++fitTime[i].tf >= Ranges[i].tf + Ranges[i].dtf )
        {
          fitTime[i].tf = Ranges[i].tf;
          bOverflow = true;
        }
      }
    }
    if( bOverflow )
      back().tf = Ranges.back().tf + Ranges.back().dtf;
  }
  return *this;
}

std::string FitRangesIterator::to_string( const std::string &Sep1, const std::string &Sep2 )
{
  std::string s;
  std::vector<FitTime> &fitTime( *this );
  for( std::size_t i = 0; i < fitTime.size(); ++i )
  {
    if( i )
      s.append( Sep2 );
    s.append( std::to_string( fitTime[i].ti ) );
    s.append( Sep1 );
    s.append( std::to_string( fitTime[i].tf ) );
  }
  return s;
}

FitterThread::FitterThread( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_, vCorrelator &CorrSynthetic_ )
: parent{ fitter_ },
  Extent{ fitter_.ds.Extent },
  idx{ FirstIndex }, // So we'll notice a change on the first call to SetReplica
  CholeskyDiag( Extent ),
  Data( Extent ),
  ModelParams( fitter_.NumModelParams ),
  SortingHat( parent.NumExponents ),
  bCorrelated{ bCorrelated_ },
  OutputModel{ outputModel_ },
  CorrSynthetic( CorrSynthetic_ ),
  ModelBuffer( fitter_.NumFiles )
{
  if( bCorrelated )
  {
    Covar.resize( Extent, Extent );
    Cholesky.resize( Extent, Extent );
  }
  // Make a buffer for each model to use as a scratchpad
  for( int f = 0; f < parent.NumFiles; f++ )
    ModelBuffer[f].resize( parent.model[f]->GetScratchPadSize() );
}

// Make the covariance matrix, for only the timeslices we are interested in
void FitterThread::SetReplica( int idx_, bool bShowOutput )
{
  // Don't make the same covariance matrix twice in a row
  if( idx_ != idx )
  {
    // Switch to the requested replica
    const bool bFirstTime{ idx == FirstIndex };
    idx = idx_;
    // Don't build covariance matrix if frozen (unless it's the first time)
    if( !parent.bFreezeCovar || bFirstTime )
    {
      // Covariance matrix always comes from the central replica if frozen
      if( parent.bFreezeCovar )
        idx_ = Fold::idxCentral;
      // Make covariance or inverse error
      if( !bCorrelated )
        parent.ds.MakeInvErr( idx_, CholeskyDiag );
      else
      {
        parent.ds.MakeCovariance( idx_, Covar );
        // Cholesky decompose the covariance matrix, extracting diagonals for condition number
        Cholesky = Covar.Cholesky( CholeskyDiag );
        if( !Cholesky.IsFinite() ) // I don't think this is needed because GSL checks
          throw std::runtime_error( "Cholesky decomposition of covariance matrix isn't finite" );
        if( idx_ == Fold::idxCentral && bShowOutput )
        {
          const double CondNumber{ Cholesky.CholeskyRCond() };
          const int CondDigits{ static_cast<int>( 0.5 - std::log10( CondNumber ) ) };
          //std::cout << Common::NewLine;
          if( CondDigits >= 12 )
            std::cout << "WARNING see https://www.gnu.org/software/gsl/doc/html/linalg.html#cholesky-decomposition\n";
          std::cout << "Covariance reciprocal condition number " << CondNumber
                    << ", ~" << CondDigits << " digits\n";
        }
        // Allow GSL or Minuit2 to tailor the covariance matrix to their requirements
        MakeCovarCorrelated();
      }
    }
    // Now get the data for this replica
    parent.ds.GetData( idx, Data );
    // And get associated constants
    parent.ds.GetFixed( idx, ModelParams, parent.ParamFixed );
  }
}

void FitterThread::ReplicaMessage( const ParamState &state, int iFitNum ) const
{
  double ChiSq{ state.Fval() };
  std::cout << ReplicaString( iFitNum ) << ", calls " << state.NFcn() << ", chi^2 " << ChiSq << Common::NewLine;
  if( state.bGotMinuit2State )
    std::cout << "edm " << state.Edm() << ", ";
  if( state.bGotGSLState )
  {
    std::cout << "Stop: " << ( state.gslState.ConvergeReason == 1 ? "step size" : "gradient" )
              << ", f()=" << state.gslState.nevalf << ", df()=" << state.gslState.nevaldf << ", ";
  }
  std::cout << "dof " << parent.dof << ", chi^2/dof " << ( ChiSq / parent.dof ) << Common::NewLine;
  if( parent.Verbosity > 1 )
    std::cout << state;
  else
    std::cout << state.Parameters();
}

// Compute Cholesky scaled (theory - data) based on parameters from the fitting engine
bool FitterThread::SaveError( Vector &Error, const scalar * FitterParams, std::size_t Size, std::size_t Stride )
{
  // Put the fitter parameters into our model parameters
  assert( Size == parent.NumVariable && "Fitter parameters don't match our variable parameters" );
  for( int i = 0; i < parent.NumVariable; ++i )
    ModelParams[parent.ParamVariable[i]] = FitterParams[i * Stride];
  // Compute theory - error for each timeslice
  int i{0};
  for( int f = 0; f < parent.NumFiles; ++f )
  {
    // Allow each model to cache any common computations based solely on the model parameters
    const Model &m{ *parent.model[f] };
    if( ModelBuffer[f].size )
      m.ModelParamsChanged( ModelBuffer[f], ModelParams );
    const vInt &FitTimes{ parent.ds.FitTimes[f] };
    for( int t : FitTimes )
    {
      Error[i] = ( m( t, ModelBuffer[f], ModelParams ) - Data[i] ) * CholeskyDiag[i];
      if( !std::isfinite( Error[i++] ) )
        return false;
    }
  }
  return true;
}

bool FitterThread::AnalyticJacobian( Matrix &Jacobian ) const
{
  return false;
  /*for( int i = 0; i < parent.Extent; ++i )
  {
    const int f{ i / parent.NtCorr };
    const int t{ i % parent.NtCorr + parent.tMin };
    for( int p = 0; p < parent.NumParams; ++p )
    {
      double z = 1;// TODO: (*model[f]).Derivative( t, p ) * CholeskyDiag[i];
      if( !std::isfinite( z ) )
        return false;
      Jacobian( i, p ) = z;
    }
  }
  if( bCorrelated )
  {
    Jacobian.blas_trmm( CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1, Cholesky );
    //std::cout << Covar << Common::NewLine << Common::NewLine << Common::NewLine << CovarInv << Common::NewLine;
  }
  return true;*/
}

scalar FitterThread::RepeatFit( ParamState &Guess, int MaxGuesses )
{
  // Call the minimiser until it provides the same answer twice
  double dTestStat = -747;
  int iNumGuesses{ 0 };
  bool bFinished{ false };
  while( !bFinished && iNumGuesses++ <= MaxGuesses )
  {
    Minimise( Guess, iNumGuesses );
    double dNewTestStat{ Guess.Fval() };
    if( MaxGuesses == 0 )
      bFinished = true;
    else if( iNumGuesses != 1 )
      bFinished = ( dTestStat == dNewTestStat );
    dTestStat = dNewTestStat;
    if( idx == Fold::idxCentral && ( bFinished || parent.Verbosity ) )
      ReplicaMessage( Guess, iNumGuesses );
  }
  if( !bFinished )
    throw std::runtime_error( "Fit did not converge within " + std::to_string( MaxGuesses ) + " retries" );
  return dTestStat;
}

// Perform an uncorrelated fit on the central timeslice to update parameters guessed
void FitterThread::UpdateGuess( Parameters &Guess )
{
  bool bSaveCorrelated{ bCorrelated };
  bCorrelated = false;
  ParamState ThisGuess( Guess );
  try
  {
    RepeatFit( ThisGuess, NumRetriesGuess() );
  }
  catch(...)
  {
    bCorrelated = bSaveCorrelated;
    throw;
  }
  bCorrelated = bSaveCorrelated;
  Guess = ThisGuess.Parameters();
}

// Perform a fit. NB: Only the fit on central replica updates ChiSq
scalar FitterThread::FitOne( const Parameters &parGuess, const std::string &SaveCorMatFileName )
{
  // Perform fit
  ParamState Result{ parGuess };
  scalar dTestStat = RepeatFit( Result, NumRetriesFit() );
  // Put the variable fit parameters into the full model parameter set (constants are already there)
  for( int i = 0; i < parent.NumVariable; ++i )
    ModelParams[parent.ParamVariable[i]] = Result.parameters[i].Value;
  // Save the fit parameters for this replica, sorted by E_0
  for( int e = 0; e < parent.NumExponents; ++e )
  {
    SortingHat[e].first = ModelParams[e * parent.NumPerExp];
    SortingHat[e].second = e;
  }
  std::sort( SortingHat.begin(), SortingHat.end() );

  // Copy results into OutputModel
  scalar * const OutputData{ OutputModel[idx] };
  for( int e = 0; e < parent.NumExponents; ++e )
  {
    const int src{ SortingHat[e].second * parent.NumPerExp };
    const int dst{            e         * parent.NumPerExp };
    for( int i = 0; i < parent.NumPerExp; ++i )
      OutputData[dst + i] = ModelParams[src + i];
  }
  // Now copy the one-off parameters and append Chi^2 per dof
  int OutputDataSize{ parent.NumExponents * parent.NumPerExp };
  for( int i = 0; i < parent.NumOneOff; ++i, ++OutputDataSize )
    OutputData[OutputDataSize] = ModelParams[OutputDataSize];
  OutputData[OutputDataSize++] = dTestStat / parent.dof;

  // Check whether energy levels are separated by minimum separation
  if( idx == Fold::idxCentral )
  {
    for( int e = 1; e < parent.NumExponents; e++ )
    {
      double RelSep = OutputData[e * parent.NumPerExp] / OutputData [ ( e - 1 ) * parent.NumPerExp];
      if( ( RelSep - 1 ) < parent.RelEnergySep || RelSep * parent.RelEnergySep > 1 )
        throw std::runtime_error( "Fit failed energy separation criteria: "
                                 + std::to_string( 1 + parent.RelEnergySep ) + " < E_n / E_{n-1} < "
                                 + std::to_string( 1 / parent.RelEnergySep ) );
    }
  }

  // Save correlation matrix for central replica
  // NB: file name will only be empty on the central replica, so only one thread will do this
  if( bCorrelated && !SaveCorMatFileName.empty() )
  {
    if( parent.model.size() != parent.ds.corr.size() )
      throw std::runtime_error( "ModelSet doesn't match DataSet" );
    std::vector<std::vector<std::string>> ExtraStrings( parent.ds.corr.size() );
    for( std::size_t f = 0; f < parent.model.size(); ++f )
    {
      ExtraStrings[f].resize( 3 );
      // Name of each operator THIS ONLY WORKS FOR MODELS WITH OVERLAP CONSTANTS - NEEDS TO BE FIXED
      for( std::size_t op = 1; op < parent.model[f]->ParamNamesPerExp.size(); ++op )
        ExtraStrings[f][0].append( parent.model[f]->ParamNamesPerExp[op] );
      if( parent.model[f]->ParamNamesPerExp.size() == 2 )
        ExtraStrings[f][0].append( parent.model[f]->ParamNamesPerExp[1] );
      // Operators - one-off
      for( const std::string &op : parent.model[f]->ParamNames )
      {
        if( !ExtraStrings[f][1].empty() )
          ExtraStrings[f][1].append( Common::CommaSpace );
        ExtraStrings[f][1].append( op );
      }
      // Operators per energy
      for( const std::string &op : parent.model[f]->ParamNamesPerExp )
      {
        if( !ExtraStrings[f][2].empty() )
          ExtraStrings[f][2].append( Common::CommaSpace );
        ExtraStrings[f][2].append( op );
      }
    }
    // Before we plot the matrix, we need to: a) remove the scaling b) make sure the upper triangle is valid
    Matrix CovarScaled{ Covar };
    CovarScaled.CholeskyScaleApply( CholeskyDiag );
    for( int i = 1; i < CovarScaled.size1; ++i )
      for( int j = 0; j < i; ++j )
        CovarScaled( j, i ) = CovarScaled( i, j );
    parent.ds.SaveCovariance( SaveCorMatFileName, CovarScaled, &ExtraStrings );
  }

  // Save the reconstructed correlator values for this replica
  // Put the fitter parameters into our model parameters
  for( int i = 0; i < parent.NumModelParams; ++i )
    ModelParams[i] = OutputData[i];
  // Compute theory values for each timeslice of synthetic correlators
  for( int f = 0; f < parent.NumFiles; ++f )
  {
    // Allow each model to cache any common computations based solely on the model parameters
    const Model &m{ *parent.model[f] };
    if( ModelBuffer[f].size )
      m.ModelParamsChanged( ModelBuffer[f], ModelParams );
    scalar * SyntheticData{ CorrSynthetic[f][idx] };
    for( int t = 0; t < m.Nt; ++t )
      *SyntheticData++ = m( t, ModelBuffer[f], ModelParams );
  }
  return dTestStat;
}

const ROOT::Minuit2::MnStrategy FitterThreadMinuit2::Strategy( FitterThreadMinuit2::StrategyLevel );

// Painful, but Minuit2 calls operator() const
double FitterThreadMinuit2::operator()( const std::vector<double> & par ) const
{
  scalar localBuffer[Extent];
  Vector Error;
  Error.MapView( localBuffer, Extent );
  FitterThreadMinuit2 * pMe{ const_cast<FitterThreadMinuit2 *>( this ) };
  if( !pMe->SaveError( Error, &par[0], par.size() ) )
    return std::numeric_limits<double>::max();
  double chi2;
  if( bCorrelated )
    chi2 = Error.Dot( Cholesky.CholeskySolve( Error ) );
  else
  {
    chi2 = 0;
    for( int i = 0; i < Extent; ++i )
      chi2 += Error[i] * Error[i];
  }
  if( chi2 < 0 )
    throw std::runtime_error( "Chi^2 < 0 on replica " + std::to_string( idx ) );
  return chi2;
}

void FitterThreadMinuit2::Minimise( ParamState &Guess, int iNumGuesses )
{
  ROOT::Minuit2::MnUserParameters Minuit2Par;
  if( !Guess.bGotMinuit2State )
    for( const Parameters::Parameter & p : Guess.parameters.Params )
      Minuit2Par.Add( p.Name, p.Value, p.Error );
  ROOT::Minuit2::FunctionMinimum min = Minimiser.Minimize( *this, Guess.bGotMinuit2State ? Guess.Minuit2State : Minuit2Par,
                                                           Strategy, parent.MaxIt, parent.Tolerance * 1000 );
  const ROOT::Minuit2::MnUserParameterState &state{ min.UserState() };
  if( !state.IsValid() )
    throw std::runtime_error( ReplicaString( iNumGuesses ) + " did not converge" );
  Guess = state;
}

FitterThreadGSL::FitterThreadGSL( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_,
                                  vCorrelator &CorrSynthetic_ )
: FitterThread( fitter_, bCorrelated_, outputModel_, CorrSynthetic_ ), vGuess( parent.NumVariable )
{
  // Define my finite difference function
  std::memset( &fdf, 0, sizeof( fdf ) );
  /* define the function to be minimized */
  fdf.f = &sf;
  if( parent.bAnalyticDerivatives )
    fdf.df = &sdf; // Analytic derivatives
  fdf.n = Extent;
  fdf.p = parent.NumVariable;
  fdf.params = this;

  /* allocate workspace with default parameters */
  gsl_multifit_nlinear_parameters fdf_params{ gsl_multifit_nlinear_default_parameters() };
  ws = gsl_multifit_nlinear_alloc( gsl_multifit_nlinear_trust, &fdf_params, fdf.n, fdf.p );
}

FitterThreadGSL::~FitterThreadGSL()
{
  if( ws )
    gsl_multifit_nlinear_free( ws );
}

void FitterThreadGSL::MakeCovarCorrelated()
{
  // Cholesky gives LL^T, so inverting lower triangle (L) gives L^{-1} in lower triangle
  Cholesky.CholeskyInvert();
  Cholesky.Cholesky();
}

int FitterThreadGSL::f( const Vector &FitParams, Vector &Error )
{
  assert( FitParams.size == parent.NumVariable && "Parameter vector is not the right size" );
  assert( Error.size == Extent && "Result vector is not the right size" );
  if( !SaveError( Error, FitParams.data, FitParams.size, FitParams.stride ) )
    throw std::runtime_error( "Error computing residuals" );
  if( bCorrelated )
  {
    Error.blas_trmv( CblasLower, CblasTrans, CblasNonUnit, Cholesky );
    //std::cout << Covar << Common::NewLine << Common::NewLine << Common::NewLine << CovarInv << Common::NewLine;
  }
  return 0;
}

int FitterThreadGSL::df( const Vector &x, Matrix &J )
{
  assert( x.size == parent.NumVariable && "Parameter vector is not the right size" );
  assert( J.size1 == Extent && "Jacobian rows != data points" );
  assert( J.size2 == parent.NumVariable && "Parameter columns != parameters" );
  for( int f = 0; f < parent.NumFiles; f++ )
    ;// TODO: (*model[f]).Init( x );
  if( !AnalyticJacobian( J ) )
    throw std::runtime_error( "Error computing Jacobian" );
  return 0;
}

void FitterThreadGSL::Minimise( ParamState &Guess, int iNumGuesses )
{
  std::size_t i = 0;
  for( const Parameters::Parameter & p : Guess.parameters.Params )
    vGuess[i++] = p.Value;

  /* initialize solver with starting point and weights */
  gsl_multifit_nlinear_init( &vGuess, &fdf, ws );

  /* compute initial cost function */
  if( idx == Fold::idxCentral && parent.Verbosity > 1 )
  {
    const Vector &vResidual{ * reinterpret_cast<Vector *>( gsl_multifit_nlinear_residual( ws ) ) };
    gsl_blas_ddot( &vResidual, &vResidual, &Guess.TestStat );
    std::cout << "Guess chi^2=" << Guess.TestStat << Common::NewLine;
  }

  /* solve the system with a maximum of parent.MaxIt iterations */
  const auto tol = parent.Tolerance;
  Guess.bValid = false;
  gsl_multifit_nlinear_driver( parent.MaxIt ? parent.MaxIt : std::numeric_limits<std::size_t>::max(), // Infinite
                               tol, tol, tol, nullptr, nullptr, &Guess.gslState.ConvergeReason, ws );
  Guess.bValid = true;
  Guess.bGotGSLState = true;
  {
    const std::size_t nIter{ gsl_multifit_nlinear_niter( ws ) };
    Guess.NumCalls = nIter < std::numeric_limits<unsigned int>::max() ? static_cast<unsigned int>( nIter ) : std::numeric_limits<unsigned int>::max();
    Guess.gslState.nevalf = fdf.nevalf;
    Guess.gslState.nevaldf = fdf.nevaldf;
  }
  const Vector &vResidual{ * reinterpret_cast<Vector *>( gsl_multifit_nlinear_residual( ws ) ) };
  gsl_blas_ddot( &vResidual, &vResidual, &Guess.TestStat );
  i = 0;
  const Vector &vResult{ * reinterpret_cast<Vector *>( gsl_multifit_nlinear_position( ws ) ) };
  Matrix &mJacobian{ * reinterpret_cast<Matrix *>( gsl_multifit_nlinear_jac( ws ) ) };
  Matrix mErrors( parent.NumVariable, parent.NumVariable );
  gsl_multifit_nlinear_covar( &mJacobian, 0, &mErrors );
  for( Parameters::Parameter & p : Guess.parameters.Params )
  {
    p.Error = std::sqrt( mErrors( i, i ) );
    p.Value = vResult[i++];
  }
  if( parent.Verbosity > 1 && !fdf.df && idx == Fold::idxCentral )
  {
    // Compare numeric derivatives to the analytic ones I would have computed
    Matrix MyJacobian( parent.ds.Extent, parent.NumVariable );
    AnalyticJacobian( MyJacobian );
    std::cout << Common::NewLine << "GSL Jacobian:\n" << mJacobian << Common::NewLine
              << "My Jacobian:\n" << MyJacobian << Common::NewLine;
    MyJacobian.cols(); // Debug breakpoint here
  }
}

std::string FitterThreadGSL::Description() const
{
  std::stringstream ss;
  ss << gsl_multifit_nlinear_name( ws ) << "/" << gsl_multifit_nlinear_trs_name( ws ) << " with ";
  if( fdf.df )
    ss << "analytic";
  else
    ss << "numeric";
  ss << " derivatives";
  return ss.str();
}

Fitter::Fitter( FitterType fitType_, const DataSet &ds_, const std::vector<std::string> &ModelArgs,
                const ModelDefaultParams &modelDefault, const std::vector<std::string> &opNames_,
                int verbosity_, bool bFreezeCovar_, bool bSaveCorr_, bool bSaveCMat_,
                int Retry_, int MaxIt_, double Tolerance_, double RelEnergySep_, bool bAnalyticDerivatives_ )
  : fitType{fitType_},
    ds{ std::move( ds_ ) },
    bAnalyticDerivatives{bAnalyticDerivatives_},
    NumOps{ static_cast<int>( opNames_.size() ) },
    OpNames{ opNames_ },
    bForceSrcSnkDifferent{ modelDefault.bForceSrcSnkDifferent },
    Verbosity{ verbosity_ },
    bFreezeCovar{bFreezeCovar_},
    bSaveCorr{bSaveCorr_},
    bSaveCMat{bSaveCMat_},
    Retry{Retry_},
    MaxIt{MaxIt_},
    Tolerance{Tolerance_},
    RelEnergySep{RelEnergySep_},
    NumFiles{ static_cast<int>( ds.corr.size() ) },
    model{ CreateModels( ModelArgs, modelDefault ) },
    NumExponents{ GetNumExponents() },
    PerExpNames{ MakePerExpNames() },
    ParamNames( MakeParamNames() ),
    ParamFixed{ MakeParamFixed() },
    ParamVariable{ MakeParamVariable() },
    NumModelParams{ static_cast<int>( ParamNames.size() ) },
    NumPerExp{ static_cast<int>( PerExpNames.size() ) },
    NumOneOff{ NumModelParams - NumExponents * NumPerExp },
    NumFixed{ static_cast<int>( ParamFixed.size() ) },
    NumVariable{ static_cast<int>( ParamVariable.size() ) }
{
  assert( ds.corr.size() == NumFiles && "Number of files extremely large" );
  assert( NumModelParams == NumFixed + NumVariable && "NumModelParams doesn't match fixed and variable" );
  assert( NumModelParams == NumExponents * NumPerExp + NumOneOff && "NumModelParams doesn't match NumPerExp and NumOneOff" );
  switch( fitType )
  {
    case FitterType::Minuit2:
    case FitterType::GSL:
      break;
    default:
      throw std::runtime_error( "Unknown FitterType " + std::to_string( static_cast<int>( fitType ) ) );
  }
}

// Create all the models and return the number of exponents in the fit
std::vector<ModelPtr> Fitter::CreateModels(const std::vector<std::string> &ModelArgs, const ModelDefaultParams &modelDefault)
{
  const int NumModels{ static_cast<int>( ds.corr.size() ) };
  if( NumModels == 0 )
    throw std::runtime_error( "Can't construct a ModelSet for an empty data set" );
  if( NumModels != ModelArgs.size() )
    throw std::runtime_error( "ModelSet for " + std::to_string( NumModels ) + " correrlators, but "
                     + std::to_string( ModelArgs.size() ) + " arguments" );
  // Create this model
  std::vector<ModelPtr> model;
  model.reserve( ModelArgs.size() );
  for( int i = 0; i < ModelArgs.size(); ++i )
  {
    std::vector<std::string> vThisArg{ Common::ArrayFromString( ModelArgs[i] ) };
    model.emplace_back( Model::MakeModel( vThisArg, modelDefault, ds.corr[i], OpNames ) );
    if( !vThisArg.empty() )
      throw std::runtime_error( "Model " + std::to_string( i ) + " has " + std::to_string( vThisArg.size() )
                               + " leftover parameters \"" + ModelArgs[i] + "\"" );
  }
  return model;
}

// Create all the models and return the number of exponents in the fit
int Fitter::GetNumExponents()
{
  int MaxExponents{0};
  for( int i = 0; i < model.size(); ++i )
  {
    if( i == 0 )
      MaxExponents = model[i]->NumExponents;
    else if( MaxExponents < model[i]->NumExponents )
      MaxExponents = model[i]->NumExponents;
  }
  return MaxExponents;
}

// The assumption is that only per energy level constants might or might not be soluble
std::size_t Fitter::EnsureModelsSolubleHelper( UniqueNames &Names, std::size_t &NumWithUnknowns )
{
  std::size_t NumUnknown{ 0 };
  std::size_t LastUnknown;
  do
  {
    LastUnknown = NumUnknown;
    NumUnknown = 0;
    NumWithUnknowns = 0;
    for( ModelPtr &m : model )
    {
      std::size_t z{ m->UnknownParameterCount( Names ) };
      if( z )
      {
        NumUnknown += z;
        NumWithUnknowns++;
      }
    }
  }
  while( NumUnknown && NumUnknown != LastUnknown );
  return NumUnknown;
}

// Finalise the parameter lists the models will fit, then build a list of per-exponential parameter names
// Tell each model which per-exponential parameters they are interested in
std::vector<std::string> Fitter::MakePerExpNames()
{
  // Loop through all the models making sure we can solve for the parameters we've been given
  // Throws an error if problematic
  UniqueNames Names{ ds.ConstantNamesPerExp };
  std::size_t NumWithUnknowns;
  std::size_t NumUnknown{ EnsureModelsSolubleHelper( Names, NumWithUnknowns ) };
  // See whether we can live with any remaining unknowns
  if( NumUnknown > NumWithUnknowns )
  {
    for( ModelPtr &m : model )
      m->ReduceUnknown( Names ); // Which of course invalidates our name list
    Names = ds.ConstantNamesPerExp;
    NumUnknown = EnsureModelsSolubleHelper( Names, NumWithUnknowns );
    if( NumUnknown > NumWithUnknowns )
      throw std::runtime_error( "Insoluble: " + std::to_string( NumUnknown ) + " unknowns > "
                               + std::to_string( NumUnknown ) + " correlators" );
  }
  // Make a sorted list of all the per exponential parameters but WITHOUT energy, as this must be first
  Names.clear();
  for( ModelPtr &m : model )
    for( const std::string &s : m->ParamNamesPerExp )
      if( !Common::EqualIgnoreCase( E, s ) )
        Names[s];
  // Sort per-exponential parameters by name, adding energy as parameter 0
  {
    int PerExpSort{ 1 };
    for( UniqueNames::iterator it = Names.begin(); it != Names.end(); ++it )
      it->second = PerExpSort++;
    Names.insert( { E, 0 } );
  }
  // Now tell every model which parameter to use for per-exponential parameters
  const int NumParamsPerExp{ static_cast<int>( Names.size() ) };
  for( ModelPtr &m : model )
  {
    m->ParamIdxPerExp.resize( m->NumExponents );
    const int ModelPerExp{ static_cast<int>( m->ParamNamesPerExp.size() ) };
    for( int i = 0; i < ModelPerExp; ++i )
    {
      UniqueNames::iterator it = Names.find( m->ParamNamesPerExp[i] );
      assert( it != Names.end() && "ModelSet constructor buggy" );
      int idx{ it->second };
      for( int e = 0; e < m->NumExponents; ++e )
      {
        if( i == 0 )
          m->ParamIdxPerExp[e].resize( ModelPerExp );
        m->ParamIdxPerExp[e][i] = idx;
        idx += NumParamsPerExp;
      }
    }
  }
  // Now Make our list of parameter names per exponent
  std::vector<std::string> PerExpNames( NumParamsPerExp );
  for( UniqueNames::iterator it = Names.begin(); it != Names.end(); ++it )
    PerExpNames[it->second] = it->first;
  return PerExpNames;
}

// Make the full list of all parameters. Tell models about one-off parameters they are interested in
std::vector<std::string> Fitter::MakeParamNames()
{
  // Make a sorted list of all the per exponential parameters but WITHOUT energy, as this must be first
  int idxParam{ 0 };
  UniqueNames Names;
  std::vector<std::string> ParamNames( PerExpNames.size() * NumExponents );
  for( int e = 0; e < NumExponents; ++e )
  {
    const std::string eString{ std::to_string( e ) };
    for( std::string s : PerExpNames )
    {
      s.append( eString );
      ParamNames[idxParam] = s;
      Names.insert( { s, idxParam++ } );
    }
  }
  // Add one-off parameters from every model to the parameter list - if not there already
  // While we're at it, tell each model the indices of the parameters they are interested in
  for( ModelPtr &m : model )
  {
    const std::size_t peSize{ m->ParamNames.size() };
    m->ParamIdx.resize( peSize );
    for( int i = 0; i < peSize; ++i )
    {
      const std::string &s{ m->ParamNames[i] };
      UniqueNames::iterator it = Names.find( s );
      if( it == Names.end() )
      {
        m->ParamIdx[i] = idxParam;
        Names.insert( { s, idxParam++ } );
        ParamNames.push_back( s );
      }
      else
        m->ParamIdx[i] = it->second;
    }
  }
  return ParamNames;
}

// Now work out which parameters are fixed and which variable. Save info for back-references
std::vector<DataSet::FixedParam> Fitter::MakeParamFixed()
{
  std::vector<DataSet::FixedParam> vFixed;
  const int NumParams_{ static_cast<int>( ParamNames.size() ) };
  assert( NumParams_ <= std::numeric_limits<int>::max() );
  for( int i = 0; i < NumParams_; ++i )
  {
    DataSet::ConstMap::const_iterator it{ ds.constMap.find( ParamNames[i] ) };
    if( it != ds.constMap.end() )
      vFixed.emplace_back( i, it->second ); // fixed parameter
  }
  return vFixed;
}

// Now work out which parameters are fixed and which variable. Save info for back-references
std::vector<int> Fitter::MakeParamVariable()
{
  const int NumParams_{ static_cast<int>( ParamNames.size() ) };
  assert( NumParams_ <= std::numeric_limits<int>::max() );
  const int NumFixed_{ static_cast<int>( ParamFixed.size() ) };
  std::vector<int> vFree;
  vFree.reserve( NumParams_ - NumFixed_ );
  for( int i = 0; i < NumParams_; ++i )
  {
    DataSet::ConstMap::const_iterator it{ ds.constMap.find( ParamNames[i] ) };
    if( it == ds.constMap.end() )
      vFree.push_back( i );
  }
  if( vFree.empty() )
    throw std::runtime_error( "Model has no parameters to fit" );
  return vFree;
}

// Perform a fit - assume fit ranges have been set on the DataSet prior to the call
std::vector<Common::ValWithEr<scalar>>
Fitter::PerformFit( bool Bcorrelated, double &ChiSq, int &dof_, const std::string &OutBaseName,
                    const std::string &ModelSuffix, Common::SeedType Seed )
{
  bCorrelated = Bcorrelated;
  dof = ds.Extent - NumVariable;
  if( dof <= 0 )
    throw std::runtime_error( "Fit has " + std::to_string( dof ) + " degrees of freedom" );
  dof_ = dof;

  // Make somewhere to store the results of the fit for each bootstrap sample
  const int tMin{ ds.FitTimes[0][0] };
  const int tMax{ ds.FitTimes[0].back() };
  ModelFile OutputModel( OpNames, NumExponents, NumFiles, tMin, tMax, dof, !bForceSrcSnkDifferent, bFreezeCovar, ds.NSamples,
                         NumModelParams + 1 );
  for( const Fold &f : ds.corr )
    OutputModel.FileList.emplace_back( f.Name_.Filename );
  OutputModel.CopyAttributes( ds.corr[0] );
  {
    std::vector<std::string> ColNames{ ParamNames };
    ColNames.push_back( "ChiSqPerDof" );
    OutputModel.SetColumnNames( ColNames );
  }

  // See whether this fit already exists
  bool bPerformFit{ true };
  const std::string sModelBase{ OutBaseName + ModelSuffix };
  const std::string ModelFileName{ Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ) };
  if( Common::FileExists( ModelFileName ) )
  {
    ModelFile PreBuilt;
    PreBuilt.Read( ModelFileName, "Pre-built: " );
    if( std::tie( NumExponents, /*dof,*/ tMin, tMax, OutputModel.GetColumnNames() )
        != std::tie( PreBuilt.NumExponents, /*PreBuilt.dof,*/ PreBuilt.ti, PreBuilt.tf, PreBuilt.GetColumnNames() ) )
    {
      throw std::runtime_error( "Pre-built fit not compatible with parameters from this run" );
    }
    if( dof != PreBuilt.dof )
      throw std::runtime_error( "Pre-built fit had different constraints" );
    bPerformFit = PreBuilt.NewParamsMorePrecise( bFreezeCovar, ds.NSamples );
    if( !bPerformFit )
    {
      ChiSq = dof * PreBuilt.getSummaryData()[NumModelParams].Central; // Last summary has chi squared per dof
      OutputModel = std::move( PreBuilt );
    }
    else
      std::cout << "Overwriting\n";
  }

  if( bPerformFit )
  {
    // If we're saving the correlation matrix, it should have a similar name to the model
    const std::string SaveCorrMatrixFileName{ bCorrelated && bSaveCMat
                            ? Common::MakeFilename( sModelBase, Common::sCormat, Seed, TEXT_EXT ) : "" };
    // Make somewhere to hold the correlators corresponding to the fitted model
    vCorrelator CorrSynthetic( NumFiles ); // correlators resulting from the fit params
    for( int f = 0; f < NumFiles; f++ )
    {
      CorrSynthetic[f].resize( ds.NSamples, ds.corr[f].Nt() );
      CorrSynthetic[f].FileList.push_back( ModelFileName );
      CorrSynthetic[f].CopyAttributes( ds.corr[f] );
    }

    // Build Parameter lists and take initial guesses
    Parameters parGuess; // Only the variable parameters - used by fitting engines
    Vector modelParams( NumModelParams ); // All parameters - used by models
    {
      for( int i = 0; i < NumModelParams; ++i )
        modelParams[i] = 0;
      // Get all the constants and start keeping track of which parameters we have managed to take a guess for
      std::vector<bool> bKnown( NumModelParams, false );
      ds.GetFixed( Fold::idxCentral , modelParams, ParamFixed );
      for( const DataSet::FixedParam &p : ParamFixed )
        bKnown[p.idx] = true;
      const std::vector<bool> bFixed( bKnown );
      // Take a guess for E0 if we didn't manage to load it
      if( !bKnown[0] )
      {
        int E0Count{ 0 };
        Vector vCorr;
        for( int f = 0; f < NumFiles; ++f )
        {
          scalar z;
          vCorr.MapView( const_cast<scalar *>( ds.corr[f][Fold::idxCentral] ), ds.corr[f].Nt() );
          if( model[f]->GuessE0( z, vCorr ) )
          {
            modelParams[0] += z;
            E0Count++;
          }
        }
        assert( E0Count && "Nothing guessed E0" );
        modelParams[0] /= E0Count;
        bKnown[0] = true;
      }
      // Now guess excited state energies (increment of half the previous difference each time)
      {
        scalar EPrior{ 0 };
        for( int e = 1; e < NumExponents; ++e )
        {
          const int idxE{ e * NumPerExp };
          const int idxELast{ idxE - NumPerExp };
          if( !bKnown[idxE] )
          {
            scalar PriorDiff{ modelParams[idxELast] - EPrior };
            modelParams[idxE] = modelParams[idxELast] + PriorDiff / 2;
            EPrior = modelParams[idxELast];
            bKnown[idxE] = true;
          }
        }
      }
      // Now guess everything else
      {
        Vector vCorr;
        int pass{ 0 };
        for( bool bNeedAnotherPass = true; bNeedAnotherPass; ++pass )
        {
          bNeedAnotherPass = false;
          for( int f = 0; f < NumFiles; ++f )
          {
            vCorr.MapView( const_cast<scalar *>( ds.corr[f][Fold::idxCentral] ), ds.corr[f].Nt() );
            if( model[f]->Guess( modelParams, bKnown, pass, vCorr ) )
              bNeedAnotherPass = true;
          }
        }
      }
      // Make sure we've taken a guess for everything
      for( int i = 0; i < NumModelParams; ++i )
      {
        assert( bKnown[i] && "Bad model: unable to guess all parameters" );
        if( !bFixed[i] )
          parGuess.Add( ParamNames[i], modelParams[i], modelParams[i] * 0.1 );
      }
    }
    // Protected by CancelCritSec
    {
      volatile bool Abort{ false };
      std::string sError;
#ifndef DEBUG_DISABLE_OMP
      #pragma omp parallel default(shared)
#endif
      {
        try
        {
          // fitThread only throws an exception if the correlation matrix is frozen and bad
          // ... which will either happen for all threads or none
          std::unique_ptr<FitterThread> ft;
          switch( fitType )
          {
            case FitterType::Minuit2:
              ft.reset( new FitterThreadMinuit2( *this, bCorrelated, OutputModel, CorrSynthetic ) );
              break;
            case FitterType::GSL:
              ft.reset( new FitterThreadGSL( *this, bCorrelated, OutputModel, CorrSynthetic ) );
              break;
          }
          FitterThread &fitThread( *ft.get() );
#ifndef DEBUG_DISABLE_OMP
          #pragma omp single
#endif
          {
            const std::string sDescription{ fitThread.Description() };
            if( !sDescription.empty() )
              std::cout << sDescription << Common::NewLine;
            std::cout << "Tolerance " << Tolerance << ". Using uncorrelated fit as guess for each replica.\n";
            if( NumFixed )
            {
              std::size_t MaxLen{ 0 };
              for( const auto & p : ParamFixed )
              {
                std::size_t ThisLen{ ParamNames[p.idx].size() };
                if( MaxLen < ThisLen )
                  MaxLen = ThisLen;
              }
              std::cout << " Fixed parameters:\n";
              for( const auto & p : ParamFixed )
                std::cout << std::string( 2 + MaxLen - ParamNames[p.idx].size(), ' ' ) << ParamNames[p.idx]
                          << Common::Space << modelParams[p.idx] << Common::NewLine;
            }
            std::cout << " Initial guess:\n" << parGuess;
            if( bCorrelated )
            {
              // Perform an uncorrelated fit on the central replica, and use that as the guess for every replica
              try
              {
                fitThread.SetReplica( Fold::idxCentral, true );
                fitThread.UpdateGuess( parGuess );
              }
              catch( const std::exception &e )
              {
                bool WasAbort;
#ifndef DEBUG_DISABLE_OMP
                #pragma omp atomic capture
#endif
                {
                  WasAbort = Abort;
                  Abort = true;
                }
                if( !WasAbort )
                  sError = e.what();
              }
            }
          }
#ifndef DEBUG_DISABLE_OMP
          #pragma omp for schedule(dynamic)
#endif
          for( int idx = Fold::idxCentral; idx < ds.NSamples; ++idx )
          {
            try
            {
              // Use uncorrelated fit as guess for correlated fit
              bool WasAbort;
#ifndef DEBUG_DISABLE_OMP
              #pragma omp atomic read
#endif
                WasAbort = Abort;
              if( !WasAbort )
              {
                fitThread.SetReplica( idx, false );
                scalar z{ fitThread.FitOne( parGuess, idx == Fold::idxCentral ? SaveCorrMatrixFileName : "" ) };
                if( idx == Fold::idxCentral )
                  ChiSq = z;
              }
            }
            catch( const std::exception &e )
            {
              bool WasAbort;
#ifndef DEBUG_DISABLE_OMP
              #pragma omp atomic capture
#endif
              {
                WasAbort = Abort;
                Abort = true;
              }
              if( !WasAbort )
                sError = e.what();
            }
          }
        }
        catch( const std::exception &e )
        {
          bool WasAbort;
#ifndef DEBUG_DISABLE_OMP
          #pragma omp atomic capture
#endif
          {
            WasAbort = Abort;
            Abort = true;
          }
          if( !WasAbort )
            sError = e.what();
        }
      }
      if( Abort )
        throw std::runtime_error( sError );
    }
    OutputModel.MakeCorrSummary( "Params" );
    OutputModel.Write( ModelFileName );
    //OutputModel.WriteSummary( Common::MakeFilename( sModelBase, Common::sModel, Seed, TEXT_EXT ) );
    for( int f = 0; f < NumFiles; f++ )
    {
      const int snk{ ds.corr[f].Name_.op[idxSnk] };
      const int src{ ds.corr[f].Name_.op[idxSrc] };
      std::string sSink{ OpNames[snk] };
      std::size_t pos = sSink.find_last_of( '_' );
      if( pos != std::string::npos )
        sSink.resize( pos );
      std::string sSrc{ OpNames[src] };
      pos = sSrc.find_last_of( '_' );
      if( pos != std::string::npos )
        sSrc.resize( pos );
      const std::string SummaryBase{ OutBaseName + sSink + '_' + sSrc };
      if( bSaveCorr )
        CorrSynthetic[f].Write( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, DEF_FMT ) );
      CorrSynthetic[f].MakeCorrSummary( nullptr );
      CorrSynthetic[f].WriteSummary( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, TEXT_EXT ));
    }
  }
  // Return the statistics on the fit results
  const int NumSummaries{ OutputModel.NumSamples() }; // because we might have read back an old fit
  std::vector<Common::ValWithEr<scalar>> Results( NumModelParams );
  std::vector<double> data( NumSummaries );
  for( int p = 0; p < NumModelParams; p++ ) {
    double Central = OutputModel[Fold::idxCentral][p];
    std::size_t Count{ 0 };
    for( int j = 0; j < NumSummaries; ++j ) {
      double d = OutputModel[j][p];
      if( std::isfinite( d ) )
        data[Count++] = d;
    }
    Results[p].Get( Central, data, Count );
  }
  return Results;
}

int main(int argc, const char *argv[])
{
  // Can't do this because of Minuit2     std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"sep", CL::SwitchType::Single, "0.2"},
      {"delta", CL::SwitchType::Single, "3"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"t", CL::SwitchType::Single, nullptr},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"e", CL::SwitchType::Single, "1"},
      {"n", CL::SwitchType::Single, "0"},
      {"v", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"savecmat", CL::SwitchType::Flag, nullptr},
      {"minuit2", CL::SwitchType::Flag, nullptr},
      {"analytic", CL::SwitchType::Flag, nullptr},
      {"srcsnk", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      const double RelEnergySep{ cl.SwitchValue<double>("sep") };
      const int delta{ cl.SwitchValue<int>("delta") };
      const int Retry{ cl.SwitchValue<int>("retry") };
      const int MaxIterations{ cl.SwitchValue<int>("iter") }; // Max iteration count, 0=unlimited
      const double Tolerance{ cl.SwitchValue<double>("tol") };
      const int Verbosity{ cl.SwitchValue<int>("v") };
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      const int NSamples{ cl.SwitchValue<int>("n") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bFreezeCovar{ cl.GotSwitch( "freeze" ) };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };
      const bool bSaveCMat{ cl.GotSwitch("savecmat") };
      const FitterType fitType{ cl.GotSwitch("minuit2") ? FitterType::Minuit2 : FitterType::GSL };
      const bool bAnalyticDerivatives{ cl.GotSwitch("analytic") };
      ModelDefaultParams modelDefault;
      modelDefault.NumExponents = cl.SwitchValue<int>( "e" );
      modelDefault.bForceSrcSnkDifferent = cl.GotSwitch( "srcsnk" );

      if( Retry < 0 )
        throw std::invalid_argument( "Retry must be >= 0" );
      if( MaxIterations < 0 )
        throw std::invalid_argument( "MaxIterations must be >= 0" );

      FitRanges fitRanges;
      if( cl.GotSwitch( "t") )
        fitRanges = Common::ArrayFromString<FitRange>( cl.SwitchValue<std::string>( "t" ) );
      if( !fitRanges.size() )
        throw std::runtime_error( "No fit ranges specified" );

      // Walk the list of parameters on the command-line, loading correlators and making models
      bShowUsage = false;
      std::vector<std::string> OpName;
      std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ ) << "Loading folded correlators\n";
      // Split each argument at the first comma (so the first part can be treated as a filename to glob
      const std::size_t NumArgs{ cl.Args.size() };
      DataSet ds( NSamples );
      std::vector<std::string> ModelArgs;
      std::vector<int> ModelFitRange;
      std::vector<int> ModelFitRangeCount( fitRanges.size() );
      for( std::size_t ArgNum = 0; ArgNum < NumArgs; ++ArgNum )
      {
        // First parameter is the filename we're looking for
        std::string FileToGlob{ Common::ExtractToSeparator( cl.Args[ArgNum] ) };
        for( const std::string &sFileName : Common::glob( &FileToGlob, &FileToGlob + 1, inBase.c_str() ) )
        {
          Common::FileNameAtt Att( sFileName );
          const bool bIsCorr{ Common::EqualIgnoreCase( Att.Type, Common::sFold ) };
          if( bIsCorr )
          {
            // This is a correlator - load it
            Att.ParseOpNames( OpName );
            ds.LoadCorrelator( std::move( Att ) );
            // Now see whether the first parameter is a fit range (i.e. integer index of a defined fit range)
            std::string sModelArgs{ cl.Args[ArgNum] };
            int ThisFitRange{ 0 };
            {
              int tmp;
              std::string sPossibleFitRange{ Common::ExtractToSeparator( sModelArgs ) };
              std::istringstream ss{ sPossibleFitRange };
              if( ss >> tmp && Common::StreamEmpty( ss ) )
              {
                if( tmp < 0 || tmp >= fitRanges.size() )
                  throw std::runtime_error( "Fit range " + sPossibleFitRange + " not defined" );
                if( !fitRanges[tmp].Validate( ds.corr.back().Nt() ) )
                {
                  std::stringstream oss;
                  oss << "Fit range " << fitRanges[tmp] << " not valid";
                  throw std::runtime_error( oss.str() );
                }
                ThisFitRange = tmp;
              }
              else
                sModelArgs = cl.Args[ArgNum]; // Unadulterated parameters
            }
            ModelArgs.push_back( sModelArgs );
            ModelFitRange.push_back( ThisFitRange );
            ++ModelFitRangeCount[ThisFitRange];
          }
          else
          {
            ds.LoadModel( std::move( Att ), cl.Args[ArgNum] );
          }
        }
      }
      if( ds.corr.empty() )
        throw std::runtime_error( "At least one correlator must be loaded to perform a fit" );
      for( int i = 0; i < ModelFitRangeCount.size(); ++i )
        if( ModelFitRangeCount[i] == 0 )
          throw std::runtime_error( "Models don't refer to all fit ranges" );
      ds.SortOpNames( OpName );
      // Describe the number of replicas
      std::cout << "Using ";
      if( ds.NSamples == ds.MaxSamples )
        std::cout << "all ";
      else
        std::cout << "first " << ds.NSamples << " of ";
      std::cout << ds.MaxSamples << " bootstrap replicas";
      if( ds.NSamples != ds.MaxSamples )
        std::cout << " (all " << ds.MaxSamples << " for var/covar)";
      std::cout << Common::NewLine;

      // I'll need a sorted, concatenated list of the operators in the fit for filenames
      std::string sOpNameConcat;
      {
        sOpNameConcat = OpName[0];
        for( std::size_t i = 1; i < OpName.size(); i++ )
        {
          sOpNameConcat.append( 1, '_' );
          sOpNameConcat.append( OpName[i] );
        }
      }

      // Now make base filenames for output
      outBaseFileName.append( ds.corr[0].Name_.Base );
      outBaseFileName.append( 1, '.' );
      outBaseFileName.append( doCorr ? "corr" : "uncorr" );
      const std::string sSummaryBase{ outBaseFileName + Common::Period + sOpNameConcat };
      outBaseFileName.append( 1, '_' );
      const std::size_t outBaseFileNameLen{ outBaseFileName.size() };
      const Common::SeedType Seed{ ds.corr[0].Name_.Seed };

      // All the models are loaded
      Fitter m( fitType, ds, ModelArgs, modelDefault, OpName, Verbosity, bFreezeCovar, bSaveCorr, bSaveCMat, Retry,
                MaxIterations, Tolerance, RelEnergySep, bAnalyticDerivatives );
      const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
      std::ofstream s;
      for( FitRangesIterator it = fitRanges.begin(); !it.AtEnd(); ++it )
      {
        // Set fit ranges
        {
          std::vector<std::vector<int>> fitTimes( ds.corr.size() );
          for( int i = 0; i < ds.corr.size(); ++i )
          {
            const FitTime &ft{ it[ModelFitRange[i]] };
            const int Extent{ ft.tf - ft.ti + 1 };
            fitTimes[i].resize( Extent );
            for( int j = 0; j < Extent; ++j )
              fitTimes[i][j] = ft.ti + j;
          }
          ds.SetFitTimes( fitTimes );
        }
        if( ds.Extent >= delta )
        {
          // Log what file we're processing and when we started
          std::time_t then;
          std::time( &then );
          try
          {
            {
              std::stringstream ss;
              ss << ( doCorr ? "C" : "unc" ) << "orrelated " << fitType << " fit on timeslices " << it.to_string( "-", ", " )
          #ifdef DEBUG_DISABLE_OMP
                 << " with Open MP disabled";
          #else
                 << " using " << omp_get_max_threads() << " Open MP threads";
          #endif
              const std::string &sMsg{ ss.str() };
              std::cout << std::string( sMsg.length(), '=' ) << Common::NewLine << sMsg << Common::NewLine;
            }
            double ChiSq;
            int dof;
            outBaseFileName.resize( outBaseFileNameLen );
            outBaseFileName.append( it.to_string( Common::Underscore ) );
            outBaseFileName.append( 1, '.' );
            auto params = m.PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat, Seed );
          }
          catch(const std::exception &e)
          {
            std::cout << "Error: " << e.what() << "\n";
          }
          // Mention that we're finished, what the time is and how long it took
          std::time_t now;
          std::time( &now );
          double dNumSecs = std::difftime( now, then );
          std::string sNow{ std::ctime( &now ) };
          while( sNow.length() && sNow[sNow.length() - 1] == '\n' )
            sNow.resize( sNow.length() - 1 );
          std::stringstream ss;
          ss << sNow << ". Total duration " << std::fixed << std::setprecision(1)
                     << dNumSecs << " seconds.\n";
          std::cout << ss.str();
        }
      }
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> Bootstrap1[,model1[,params1]] [Bootstrap2[,model2[,params2]] ...]\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, where:\n"
    " modeln  is one of {Exp, Cosh, Sinh, ThreePoint, Constant}\n"
    " params1 t<n> chooses the n'th fit range (default=0, see -t)\n"
    " paramsn Depends on the model\n"
    "and <options> are:\n"
    "--sep  Minimum relative separation between energy levels (default 0.2)\n"
    "--delta Minimum number of timeslices in fit range (default 3)\n"
    "--retry Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter Max iteration count, 0 (default) = unlimited\n"
    "--tol  Tolerance of required fits (default 1e-7)\n"
    "-t     Fit range1[,range2[,...]] (start:stop[:numstart=1[:numstop=numstart]])\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-v     Verbosity, 0 (default)=central fit results, 1=all fits, 2=detailed\n"
    "Flags:\n"
    "--minuit2  Use Minuit2 fitter (default GSL Levenberg-Marquardt)\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--savecmat Save correlation matrix\n"
    "--analytic Analytic derivatives for GSL (default: numeric)\n"
    "--srcsnk   Append _src and _snk to overlap coefficients (ie force different)\n"
    "--help     This message\n";
  }
  return iReturn;
}
