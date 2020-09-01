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

FitterThread::FitterThread( const Fitter &fitter_, bool bCorrelated_, ModelFile &modelParams_, vCorrelator &CorrSynthetic_ )
: parent{ fitter_ },
  idx{ Fold::idxCentral },
  Error( fitter_.Extent ),
  SortingHat( fitter_.NumExponents ),
  bCorrelated{bCorrelated_},
  ModelParams{modelParams_},
  CorrSynthetic( CorrSynthetic_ )
{
  // Make each model
  // TODO: model.reserve( fitter_.ds.model.size() );
  // TODO: for( const auto &m : fitter_.ds.model )
    // TODO: model.emplace_back( m->MakeThreadModel() );
  // Make somewhere to sort the results of each fit by energy level
  for( int e = 0; e < parent.NumExponents; ++e )
    SortingHat[e].resize( parent.NumOps + 1 );
}

void FitterThread::ReplicaMessage( const ParamState &state, int iFitNum ) const
{
  double ChiSq{ state.Fval() };
  std::cout << ReplicaString( iFitNum ) << ", calls " << state.NFcn() << ", chi^2 " << ChiSq << Common::NewLine;
  if( state.bGotMinuit2State )
    std::cout << "edm " << state.Edm() << ", ";
  std::cout << "dof " << parent.dof << ", chi^2/dof " << ( ChiSq / parent.dof );
  if( state.bGotGSLState )
  {
    std::cout << ", Stop: " << ( state.gslState.ConvergeReason == 1 ? "step size" : "gradient" )
              << ", f()=" << state.gslState.nevalf << ", df()=" << state.gslState.nevaldf;
  }
  std::cout << Common::NewLine;
  if( parent.Verbosity > 1 )
    std::cout << state;
  else
    std::cout << state.Parameters();
}

bool FitterThread::SaveError( Vector &ModelError ) const
{
  const scalar * CorrData{ nullptr };
  for( int i = 0; i < parent.Extent; ++i )
  {
    const int f{ i / parent.NtCorr };
    const int t{ i % parent.NtCorr + parent.tMin };
    if( i % parent.NtCorr == 0 )
      ;// TODO: CorrData = parent.Corr[f].Corr[idx];
    double z = 1;// TODO: ( (*model[f])( t ) - CorrData[t] ) * CholeskyDiag[i];
    if( !std::isfinite( z ) )
      return false;
    ModelError[i] = z;
  }
  return true;
}

bool FitterThread::AnalyticJacobian( Matrix &Jacobian ) const
{
  for( int i = 0; i < parent.Extent; ++i )
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
  return true;
}

// Make the covariance matrix, for only the timeslices we are interested in
void FitterThread::SetReplica( int idx_, bool bShowOutput )
{
  const bool bFirstTime{ CholeskyDiag.size == 0 };
  if( bFirstTime )
  {
    CholeskyDiag.resize( parent.Extent );
    if( bCorrelated )
    {
      Covar.resize( parent.Extent, parent.Extent );
      Cholesky.resize( parent.Extent, parent.Extent );
    }
  }
  else if( idx_ == idx )
    return; // Don't make the same covariance matrix twice in a row
  // Switch to the requested covariance matrix ... but Freeze on the central replica
  idx = idx_;
  if( bFirstTime && parent.bFreezeCovar )
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
  // Save the fit parameters for this replica, sorted by E_0
  for( int e = 0; e < parent.NumExponents; ++e )
  {
    SortingHat[e][0] = Result.parameters[e].Value;
    for( int o = 0; o < parent.NumOps; ++o )
      SortingHat[e][o + 1] = Result.parameters[parent.MELIndex(o, e) + parent.NumExponents].Value;
  }
  std::sort( SortingHat.begin(), SortingHat.end(),
            []( const std::vector<double> &l, const std::vector<double> &r )
            { return l[0] < r[0]; } );
  double * const FitParams{ ModelParams[idx] };
  for( int e = 0; e < parent.NumExponents; ++e )
  {
    FitParams[e] = SortingHat[e][0];
    for( int o = 0; o < parent.NumOps; ++o )
      FitParams[parent.MELIndex(o, e) + parent.NumExponents] = SortingHat[e][o + 1];
  }
  FitParams[parent.MELIndex(0, parent.NumExponents) + parent.NumExponents] = dTestStat / parent.dof;
  // Check whether energy levels are separated by minimum separation
  if( idx == Fold::idxCentral )
  {
    for( int e = 1; e < parent.NumExponents; e++ )
    {
      double RelSep = FitParams[e] / FitParams[e - 1];
      if( ( RelSep - 1 ) < parent.RelEnergySep || RelSep * parent.RelEnergySep > 1 )
        throw std::runtime_error( "Fit failed energy separation criteria: "
                                 + std::to_string( 1 + parent.RelEnergySep ) + " < E_n / E_{n-1} < "
                                 + std::to_string( 1 / parent.RelEnergySep ) );
    }
  }
  
  // Save correlation matrix for central replica
  // NB: file name will only be empty on the central replica, so only one thread will do this
  if( bCorrelated && !SaveCorMatFileName.empty() )
    parent.ds.SaveCovariance( SaveCorMatFileName, Covar );

  // Save the reconstructed correlator values for this replica
  for( int f = 0; f < parent.NumFiles; f++ )
  {
    // TODO:
    /*(*model[f]).Init( FitParams, parent.NumParams );
    scalar * mc{ CorrSynthetic[f][idx] };
    for( int t = 0; t < parent.Nt; t++ )
    {
      *mc++ = (*model[f])( t );
    }*/
  }
  return dTestStat;
}

const ROOT::Minuit2::MnStrategy FitterThreadMinuit2::Strategy( FitterThreadMinuit2::StrategyLevel );

// Compute chi-squared given parameters for multi-exponential fit
  // Energy Levels        Parameters 0 ... NumExponents - 1
  //  NB: second and subsequent energy levels are deltas
  // Overlap Coefficients Parameters NumExponents ... NumExponents ( NumOps + 1 ) - 1

double FitterThreadMinuit2::operator()( const std::vector<double> & par ) const
{
  for( int f = 0; f < parent.NumFiles; f++ )
    ;// TODO: (*model[f]).Init( par.data(), par.size() );
  if( !SaveError( Error ) )
    return std::numeric_limits<double>::max();
  double chi2;
  if( bCorrelated )
    chi2 = Error.Dot( Cholesky.CholeskySolve( Error ) );
  else
  {
    chi2 = 0;
    for( int i = 0; i < parent.Extent; ++i )
    {
      double z = Error[i];
      chi2 += z * z;
    }
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

FitterThreadGSL::FitterThreadGSL( const Fitter &fitter_, bool bCorrelated_, ModelFile &modelParams_,
                                  vCorrelator &CorrSynthetic_ )
: FitterThread( fitter_, bCorrelated_, modelParams_, CorrSynthetic_ ), vGuess( parent.NumParams )
{
  // Define my finite difference function
  std::memset( &fdf, 0, sizeof( fdf ) );
  /* define the function to be minimized */
  fdf.f = &sf;
  if( parent.bAnalyticDerivatives )
    fdf.df = &sdf; // Analytic derivatives
  fdf.n = parent.Extent;
  fdf.p = parent.NumParams;
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

int FitterThreadGSL::f( const Vector &x, Vector &f )
{
  assert( x.size == parent.NumParams && "Parameter vector is not the right size" );
  assert( f.size == parent.Extent && "Result vector is not the right size" );
  for( int file = 0; file < parent.NumFiles; file++ )
    ;// TODO: (*model[file]).Init( x );
  if( !SaveError( f ) )
    throw std::runtime_error( "Error computing residuals" );
  if( bCorrelated )
  {
    f.blas_trmv( CblasLower, CblasTrans, CblasNonUnit, Cholesky );
    //std::cout << Covar << Common::NewLine << Common::NewLine << Common::NewLine << CovarInv << Common::NewLine;
  }
  return 0;
}

int FitterThreadGSL::df( const Vector &x, Matrix &J )
{
  assert( x.size == parent.NumParams && "Parameter vector is not the right size" );
  assert( J.size1 == parent.Extent && "Jacobian rows != data points" );
  assert( J.size2 == parent.NumParams && "Parameter columns != parameters" );
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

  /* solve the system with a maximum of 100 iterations */
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
  Matrix mErrors( parent.NumParams, parent.NumParams );
  gsl_multifit_nlinear_covar( &mJacobian, 0, &mErrors );
  for( Parameters::Parameter & p : Guess.parameters.Params )
  {
    p.Error = std::sqrt( mErrors( i, i ) );
    p.Value = vResult[i++];
  }
  if( parent.Verbosity > 1 && !fdf.df && idx == Fold::idxCentral )
  {
    // Compare numeric derivatives to the analytic ones I would have computed
    Matrix MyJacobian( parent.Extent, parent.NumParams );
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
                const std::string &sopNameConcat_, int verbosity_, const std::string &outputBaseName_,
                Common::SeedType seed_, bool bFreezeCovar_, bool bSaveCorr_, bool bSaveCMat_,
                int Retry_, int MaxIt_, double Tolerance_, double RelEnergySep_, bool bAnalyticDerivatives_ )
  : fitType{fitType_},
    ds{ std::move( ds_ ) },
    bAnalyticDerivatives{bAnalyticDerivatives_},
    NumOps{ static_cast<int>( opNames_.size() ) },
    OpNames{ opNames_ },
    bFactor{ modelDefault.bFactor },
    sOpNameConcat{ sopNameConcat_ },
    Verbosity{ verbosity_ },
    OutputBaseName{ outputBaseName_ },
    Seed{ seed_ },
    bFreezeCovar{bFreezeCovar_},
    bSaveCorr{bSaveCorr_},
    bSaveCMat{bSaveCMat_},
    Retry{Retry_},
    MaxIt{MaxIt_},
    Tolerance{Tolerance_},
    RelEnergySep{RelEnergySep_},
    NumFiles{ static_cast<int>( ds.corr.size() ) },
    NumExponents{ CreateModels( ModelArgs, modelDefault ) },
    ParamNames( MakeParamNames() ),
    NumModelParams{ FixedVsVariableParams() },
    NumFixed{ static_cast<int>( ParamFixed.size() ) },
    NumVariable{ static_cast<int>( ParamVariable.size() ) }
{
  assert( ds.corr.size() == NumFiles && "Number of files extremely large" );
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
int Fitter::CreateModels( const std::vector<std::string> &ModelArgs, const ModelDefaultParams &modelDefault )
{
  const int NumModels{ static_cast<int>( ds.corr.size() ) };
  if( NumModels == 0 )
    throw std::runtime_error( "Can't construct a ModelSet for an empty data set" );
  if( NumModels != ModelArgs.size() )
    throw std::runtime_error( "ModelSet for " + std::to_string( NumModels ) + " correrlators, but "
                     + std::to_string( ModelArgs.size() ) + " arguments" );
  // Create this model
  model.clear();
  model.reserve( ModelArgs.size() );
  int MinExponents{0};
  int MaxExponents{0};
  for( int i = 0; i < ModelArgs.size(); ++i )
  {
    std::vector<std::string> vThisArg{ Common::ArrayFromString( ModelArgs[i] ) };
    model.emplace_back( Model::MakeModel( vThisArg, modelDefault, ds.corr[i], OpNames ) );
    if( !vThisArg.empty() )
      throw std::runtime_error( "Model " + std::to_string( i ) + " has " + std::to_string( vThisArg.size() )
                               + " leftover parameters \"" + ModelArgs[i] + "\"" );
    if( i == 0 )
    {
      MinExponents = model[i]->NumExponents;
      MaxExponents = model[i]->NumExponents;
    }
    else
    {
      if( MinExponents > model[i]->NumExponents )
        MinExponents = model[i]->NumExponents;
      if( MaxExponents < model[i]->NumExponents )
        MaxExponents = model[i]->NumExponents;
    }
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

std::vector<std::string> Fitter::MakeParamNames()
{
  std::vector<std::string> ParamNames;
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
  UniqueNames NamesPerExp;
  for( ModelPtr &m : model )
    for( const std::string &s : m->ParamNamesPerExp )
      if( !Common::EqualIgnoreCase( E, s ) )
        NamesPerExp[s];
  int NumParamsPerExp = static_cast<int>( NamesPerExp.size() + 1 );
  // Build our actual parameters per-exponential
  int idxParam{ 0 };
  Names.clear();
  ParamNames.reserve( NumParamsPerExp * NumExponents );
  for( int e = 0; e < NumExponents; ++e )
  {
    std::string eString( std::to_string( e ) );
    std::string s{ E };
    s.append( eString );
    Names.insert( { s, idxParam++ } );
    ParamNames.push_back( s );
    for( UniqueNames::iterator it = NamesPerExp.begin(); it != NamesPerExp.end(); ++it )
    {
      s = it->first;
      s.append( eString );
      Names.insert( { s, idxParam++ } );
      ParamNames.push_back( s );
    }
  }
  // For all the one-off parameters, add it to the parameter list if not there already
  // While we're at it, tell each model the indices of the parameters they are interested in
  int NumParamsSingle = 0;
  for( ModelPtr &m : model )
  {
    const std::size_t pnSize{ m->ParamNames.size() };
    m->ParamIdx.resize( pnSize );
    for( int i = 0; i < pnSize; ++i )
    {
      const std::string &s{ m->ParamNames[i] };
      UniqueNames::iterator it = Names.find( s );
      if( it == Names.end() )
      {
        m->ParamIdx[i] = idxParam;
        Names.insert( { s, idxParam++ } );
        ParamNames.push_back( s );
        NumParamsSingle++;
      }
      else
        m->ParamIdx[i] = it->second;
    }
  }
  // Now tell every model which parameter to use for per-exponential parameters
  int PerExpSort{ 1 };
  for( UniqueNames::iterator it = NamesPerExp.begin(); it != NamesPerExp.end(); ++it )
    it->second = PerExpSort++;
  NamesPerExp.insert( { E, 0 } );
  for( ModelPtr &m : model )
  {
    const std::size_t pnSize{ m->ParamNamesPerExp.size() };
    m->ParamIdxPerExp.resize( pnSize );
    for( int i = 0; i < pnSize; ++i )
    {
      const std::string &s{ m->ParamNamesPerExp[i] };
      UniqueNames::iterator it = NamesPerExp.find( s );
      assert( it != Names.end() && "ModelSet constructor buggy" );
      m->ParamIdxPerExp[i] = it->second;
    }
  }
  return ParamNames;
}

// Now work out which parameters are fixed and which variable. Save info for back-references
int Fitter::FixedVsVariableParams()
{
  const int NumParams_{ static_cast<int>( ParamNames.size() ) };
  assert( NumParams_ <= std::numeric_limits<int>::max() );
  for( int i = 0; i < NumParams_; ++i )
  {
    DataSet::ConstMap::const_iterator it{ ds.constMap.find( ParamNames[i] ) };
    if( it == ds.constMap.end() )
      ParamVariable.push_back( i ); // free parameter
    else
      ParamFixed.emplace_back( i, it->second ); // fixed parameter
  }
  if( ParamVariable.empty() )
    throw std::runtime_error( "Model has no parameters to fit" );
  return NumParams_;
}

// Perform a fit
std::vector<Common::ValWithEr<scalar>> Fitter::PerformFit( bool Bcorrelated_, int tMin, int tMax, double &ChiSq, int &dof_ )
{
  {
    std::stringstream ss;
    ss << ( Bcorrelated_ ? "C" : "unc" ) << "orrelated " << fitType << " fit on timeslices " << TMin_ << " to " << TMax_
#ifdef DEBUG_DISABLE_OMP
       << " with Open MP disabled";
#else
       << " using " << omp_get_max_threads() << " Open MP threads";
#endif
    const std::string &sMsg{ ss.str() };
    std::cout << std::string( sMsg.length(), '=' ) << Common::NewLine << sMsg << Common::NewLine;
  }
  bCorrelated = Bcorrelated_;
  dof = ds.Extent - NumVariable;
  if( dof <= 0 )
    throw std::runtime_error( "Fit has " + std::to_string( dof ) + " degrees of freedom" );
  dof_ = dof;

  // Make somewhere to store the results of the fit for each bootstrap sample
  ds.SetFitTimes( tMin, tMax );
  ModelFile ModelParams(OpNames, NumExponents, NumFiles, tMin, tMax, dof, bFactor, bFreezeCovar, NSamples, NumModelParams+1);
  for( const Fold &f : ds.corr )
    ModelParams.FileList.emplace_back( f.Name_.Filename );
  ModelParams.CopyAttributes( ds.corr[0] );
  {
    std::vector<std::string> ColNames{ ParamNames };
    ColNames.push_back( "ChiSqPerDof" );
    ModelParams.SetColumnNames( ColNames );
  }

  // See whether this fit already exists
  bool bPerformFit{ true };
  std::string OutputRunBase{ OutputBaseName };
  OutputRunBase.append( 1, '.' );
  OutputRunBase.append( bCorrelated ? "corr" : "uncorr" );
  OutputRunBase.append( 1, '_' );
  OutputRunBase.append( std::to_string( tMin ) );
  OutputRunBase.append( 1, '_' );
  OutputRunBase.append( std::to_string( tMax ) );
  std::string sModelBase{ OutputRunBase };
  sModelBase.append( 1, '.' );
  sModelBase.append( sOpNameConcat );
  const std::string ModelFileName{ Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ) };
  if( Common::FileExists( ModelFileName ) )
  {
    ModelFile PreBuilt;
    PreBuilt.Read( ModelFileName, "\nPre-built: " );
    if( !PreBuilt.Compatible( NumExponents, dof, tMin, tMax, OpNames ) )
      throw std::runtime_error( "Pre-existing fits not compatible with parameters from this run" );
    bPerformFit = PreBuilt.NewParamsMorePrecise( bFreezeCovar, NSamples );
    if( !bPerformFit )
    {
      ChiSq = dof * PreBuilt.getSummaryData()[NumParams].Central; // Last summary has chi squared per dof
      ModelParams = std::move( PreBuilt );
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
      CorrSynthetic[f].resize( NSamples, Nt );
      CorrSynthetic[f].FileList.push_back( ModelFileName );
      CorrSynthetic[f].CopyAttributes( Corr[f].Corr );
    }
    // Make initial guesses for the parameters
    // For each Exponent, I need the delta_E + a constant for each operator
    Parameters parGuess;
    {
      static constexpr double ErrorFactor = 0.1; // I actually have no idea what the parameter errors are
      std::vector<double> GuessValue( NumParams );
      std::vector<double> GuessError( NumParams );
      // Take a starting guess for the parameters - same as LatAnalyze
      double * const Energy{ GuessValue.data() };
      double * const Coeff{ Energy + NumExponents };
      double * const EnergyErr{ GuessError.data() };
      double * const CoeffErr{ EnergyErr + NumExponents };
      const int tGuess = 11;// TODO: { Nt / 4 };
      Energy[0] = 0;
      bool bOpDone[NumOps];
      for( int i = 0; i < NumOps; ++i )
        bOpDone[i] = false;
      for( int f = 0; f < NumFiles; f++ )
      {
        // Guess for energy is the log of the ratios of the exponentials at tGuess ... averaged over all correlators
        const double * const c{ nullptr }; // TODO: Corr[f].Corr[Fold::idxCentral] };
        double E0 = std::log( c[tGuess] / c[tGuess + 1] );
        Energy[0] += E0;
        // Guess for matrix elements in two passes, first where source and sink same, then different
        for( int Pass = 0; Pass < 2; ++Pass )
        {
          const int iOpSrc{ 0};// TODO: MELIndex( Corr[f].Corr.Name_.op[idxSrc], 0 ) };
          const int iOpSnk{ 0};// TODO: MELIndex( Corr[f].Corr.Name_.op[idxSnk], 0 ) };
          if( Pass == 0 && iOpSrc == iOpSnk && !bOpDone[iOpSrc] )
          {
            double MELGuess = std::sqrt( std::abs( c[tGuess] ) ) * std::exp( E0 * tGuess / 2 );
            Coeff[iOpSrc] = MELGuess;
            CoeffErr[iOpSrc] = MELGuess * ErrorFactor;
            bOpDone[iOpSrc] = true;
          }
          else if( Pass == 1 && iOpSrc != iOpSnk )
          {
            double MELGuess;
            if( !bOpDone[iOpSrc] && !bOpDone[iOpSnk] )
            {
              // Different ops, but this is only correlator we see them in
              // TODO: std::cout << "Warning: can't distinguish " << OpNames[Corr[f].Corr.Name_.op[idxSnk]]
                        // TODO: << " from " << OpNames[Corr[f].Corr.Name_.op[idxSrc]] << Common::NewLine;
              MELGuess = std::sqrt( std::abs( c[tGuess] ) ) * std::exp( E0 * tGuess / 2 );
            }
            else
            {
              MELGuess = c[tGuess] * std::exp( E0 * tGuess ) / Coeff[bOpDone[iOpSrc] ? iOpSrc : iOpSnk];
            }
            if( !bOpDone[iOpSrc] )
            {
              Coeff[iOpSrc] = MELGuess;
              CoeffErr[iOpSrc] = MELGuess * ErrorFactor;
              bOpDone[iOpSrc] = true;
            }
            if( !bOpDone[iOpSnk] )
            {
              Coeff[iOpSnk] = MELGuess;
              CoeffErr[iOpSnk] = MELGuess * ErrorFactor;
              bOpDone[iOpSnk] = true;
            }
          }
        }
      }
      Energy[0] /= NumFiles;
      EnergyErr[0] = Energy[0] * ErrorFactor;
      // Now guess Higher exponents - same as LatAnalyze
      //static const double MELFactor{ std::sqrt( 0.5 ) };
      static const double MELFactor{ std::sqrt( 2 ) };
      for( int e = 1; e < NumExponents; ++e )
      {
        Energy[e] = Energy[e - 1] + ( Energy[e - 1] - ( e > 1 ? Energy[e - 2] : 0 ) ) * 0.5;
        EnergyErr[e] = Energy[e] * ErrorFactor;
        for( int o = 0; o < NumOps; ++o )
        {
          int i = MELIndex( o, e );
          Coeff[i] = Coeff[ MELIndex( o, e - 1 ) ] * MELFactor;
          CoeffErr[i] = Coeff[i] * ErrorFactor;
        }
      }
      for( int i = 0; i < NumParams; ++i )
      {
        parGuess.Add( ParamNames[i], GuessValue[i], GuessError[i] );
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
              // TODO: ft.reset( new FitterThreadMinuit2( *this, bCorrelated, ModelParams, CorrSynthetic ) );
              break;
            case FitterType::GSL:
              // TODO: ft.reset( new FitterThreadGSL( *this, bCorrelated, ModelParams, CorrSynthetic ) );
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
            std::cout << "Tolerance " << Tolerance << ". Using uncorrelated fit as guess for each replica. Initial guess:\n"
                      << parGuess;
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
    // TODO: ModelParams.MakeCorrSummary( "Params" );
    // TODO: ModelParams.Write( ModelFileName );
    //ModelParams.WriteSummary( Common::MakeFilename( sModelBase, Common::sModel, Seed, TEXT_EXT ) );
    for( int f = 0; f < NumFiles; f++ )
    {
      const int snk{ 0};// TODO: Corr[f].Corr.Name_.op[idxSnk] };
      const int src{ 0};// TODO: Corr[f].Corr.Name_.op[idxSrc] };
      std::string sSink{ OpNames[snk] };
      std::size_t pos = sSink.find_last_of( '_' );
      if( pos != std::string::npos )
        sSink.resize( pos );
      std::string sSrc{ OpNames[src] };
      pos = sSrc.find_last_of( '_' );
      if( pos != std::string::npos )
        sSrc.resize( pos );
      const std::string SummaryBase{ OutputRunBase + '.' + sSink + '_' + sSrc };
      if( bSaveCorr )
        CorrSynthetic[f].Write( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, DEF_FMT ) );
      CorrSynthetic[f].MakeCorrSummary( nullptr );
      CorrSynthetic[f].WriteSummary( Common::MakeFilename( SummaryBase, Common::sBootstrap, Seed, TEXT_EXT ));
    }
  }
  // Return the statistics on the fit results
  // TODO:
  /*
  const int NumSummaries{ ModelParams.NumSamples() }; // because we might have read back an old fit
  std::vector<Common::ValWithEr<scalar>> Results( NumParams );
  std::vector<double> data( NumSummaries );
  for( int p = 0; p < NumParams; p++ ) {
    double Central = ModelParams[Fold::idxCentral][p];
    std::size_t Count{ 0 };
    for( int j = 0; j < NumSummaries; ++j ) {
      double d = ModelParams[j][p];
      if( std::isfinite( d ) )
        data[Count++] = d;
    }
    Results[p].Get( Central, data, Count );
  }*/
  return {};// TODO: Results;
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
      {"ti", CL::SwitchType::Single, nullptr},
      {"tf", CL::SwitchType::Single, nullptr},
      {"dti", CL::SwitchType::Single, "1"},
      {"dtf", CL::SwitchType::Single, "1"},
      {"sep", CL::SwitchType::Single, "0.2"},
      {"delta", CL::SwitchType::Single, "3"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"v", CL::SwitchType::Single, "0"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"e", CL::SwitchType::Single, "1"},
      {"n", CL::SwitchType::Single, "0"},
      {"f", CL::SwitchType::Flag, nullptr},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"savecmat", CL::SwitchType::Flag, nullptr},
      {"minuit2", CL::SwitchType::Flag, nullptr},
      {"analytic", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      const int ti_start{ cl.SwitchValue<int>("ti") };
      const int tf_start{ cl.SwitchValue<int>("tf") };
      const int dti_max{ cl.SwitchValue<int>("dti") };
      const int dtf_max{ cl.SwitchValue<int>("dtf") };
      const double RelEnergySep{ cl.SwitchValue<double>("sep") };
      const int delta{ cl.SwitchValue<int>("delta") };
      const int Retry{ cl.SwitchValue<int>("retry") };
      const int MaxIterations{ cl.SwitchValue<int>("iter") }; // Max iteration count, 0=unlimited
      const double Tolerance{ cl.SwitchValue<double>("tol") };
      const int Verbosity{ cl.SwitchValue<int>("v") };
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      const int NSamples{ cl.SwitchValue<int>("n") };
      //const std::string model{ opt.optionValue("m") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bFreezeCovar{ cl.GotSwitch( "freeze" ) };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };
      const bool bSaveCMat{ cl.GotSwitch("savecmat") };
      const FitterType fitType{ cl.GotSwitch("minuit2") ? FitterType::Minuit2 : FitterType::GSL };
      const bool bAnalyticDerivatives{ cl.GotSwitch("analytic") };
      ModelDefaultParams modelDefault;
      modelDefault.NumExponents = cl.SwitchValue<int>( "e" );
      modelDefault.bFactor = cl.GotSwitch( "f" );

      if( Retry < 0 )
        throw std::invalid_argument( "Retry must be >= 0" );
      if( MaxIterations < 0 )
        throw std::invalid_argument( "MaxIterations must be >= 0" );

      // Walk the list of parameters on the command-line, loading correlators and making models
      bShowUsage = false;
      Common::SeedType Seed = 0;
      std::vector<std::string> OpNameFile;
      std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ ) << "Loading folded correlators\n";
      // Split each argument at the first comma (so the first part can be treated as a filename to glob
      const std::size_t NumArgs{ cl.Args.size() };
      DataSet ds( NSamples );
      std::vector<std::string> ModelArgs;
      for( std::size_t ArgNum = 0; ArgNum < NumArgs; ++ArgNum )
      {
        std::string FileToGlob{ Common::ExtractToSeparator( cl.Args[ArgNum] ) };
        for( const std::string &sFileName : Common::glob( &FileToGlob, &FileToGlob + 1, inBase.c_str() ) )
          ds.LoadFile( sFileName, OpNameFile, ModelArgs, cl.Args[ArgNum] );
      }

      // At least one correlator must be loaded
      if( ds.corr.empty() )
        throw std::runtime_error( "At least one correlator must be loaded to perform a fit" );
      outBaseFileName.append( ds.corr[0].Name_.Base );
      Seed = ds.corr[0].Name_.Seed;

      // All the models are loaded
      std::sort( OpNameFile.begin(), OpNameFile.end() );
      std::string sOpNameConcat{ OpNameFile[0] };
      for( std::size_t i = 1; i < OpNameFile.size(); i++ )
      {
        sOpNameConcat.append( 1, '_' );
        sOpNameConcat.append( OpNameFile[i] );
      }
      Fitter m( fitType, ds, ModelArgs, modelDefault, OpNameFile, sOpNameConcat, Verbosity,
                outBaseFileName, Seed, bFreezeCovar, bSaveCorr, bSaveCMat, Retry, MaxIterations, Tolerance, RelEnergySep,
                bAnalyticDerivatives );
      std::string sSummaryBase{ outBaseFileName };
      sSummaryBase.append( 1, '.' );
      if( doCorr )
        sSummaryBase.append( "corr" );
      else
        sSummaryBase.append( "uncorr" );
      sSummaryBase.append( 1, '.' );
      sSummaryBase.append( sOpNameConcat );
      static const char Sep[] = " ";
      const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
      std::ofstream s;
      for( int tf = tf_start; tf - tf_start < dtf_max; tf++ )
      {
        bool bNeedHeader = true;
        for( int ti = ti_start; ti - ti_start < dti_max; ti++ )
        {
          if( tf - ti + 1 >= delta )
          {
            // Log what file we're processing and when we started
	    std::time_t then;
	    std::time( &then );
            try
            {
              double ChiSq;
              int dof;
              auto params = m.PerformFit( doCorr, ti, tf, ChiSq, dof );
              if( bNeedHeader )
              {
                bNeedHeader = false;
                if( !s.is_open() )
                {
                  s.open( sFitFilename );
                  Common::SummaryHeader<scalar>( s, sSummaryBase );
                  s << "# Seed " << Seed << std::endl;
                }
                else
                {
                  // two blank lines at start of new data block
                  s << "\n" << std::endl;
                }
                // Name the data series
                s << "# [tf=" << tf << "]" << std::endl;
                // Column names, with the series value embedded in the column header (best I can do atm)
                s << "tf=" << tf << Sep << "ti" << Sep << "dof";
                for( int p = 0; p < m.NumParams; p++ )
                  s << Sep << m.ParamNames[p] << Sep << m.ParamNames[p] << "_low" << Sep << m.ParamNames[p] << "_high" << Sep << m.ParamNames[p] << "_check";
                s << " ChiSqPerDof" << std::endl;
              }
              s << tf << Sep << ti << Sep << dof;
              for( int p = 0; p < m.NumParams; p++ )
                s << Sep << params[p];
              s << Sep << ( ChiSq / dof ) << std::endl;
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
    "model<n> is one of {2pt, 3pt, const}, followed by optional params\n"
    "and <options> are:\n"
    "--ti   Initial fit time\n"
    "--tf   Final   fit time\n"
    "--dti  Number of initial fit times (default 1)\n"
    "--dtf  Number of final   fit times (default 1)\n"
    "--sep  Minimum relative separation between energy levels (default 0.2)\n"
    "--delta Minimum number of timeslices in fit range (default 3)\n"
    "--retry Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter Max iteration count, 0 (default) = unlimited\n"
    "--tol  Tolerance of required fits (default 1e-7)\n"
    "-v     Verbosity, 0 (default)=central fit results, 1=all fits, 2=detailed\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "Flags:\n"
    "-f         Factorising operators (default non-factorising)\n"
    "--minuit2  Use Minuit2 fitter (default GSL Levenberg-Marquardt)\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--savecmat Save correlation matrix\n"
    "--analytic Analytic derivatives for GSL (default: numeric)\n"
    "--help     This message\n";
  }
  return iReturn;
}
