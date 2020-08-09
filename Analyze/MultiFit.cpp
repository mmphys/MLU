/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: MultiFit.cpp
 
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

#include "MultiFit.hpp"

//#define DEBUG_DISABLE_OMP

// Indices for operators in correlator names
static constexpr int idxSrc{ 0 };
static constexpr int idxSnk{ 1 };
static const char * pSrcSnk[] = { "src", "snk" };

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

Model::Model( const Fitter &parent )
: Nt{parent.Nt},
  HalfNt{ Nt - 1 }, // NB: This is half of the ORIGINAL correlator time dimension
  NumFiles{parent.NumFiles},
  NumExponents{parent.NumExponents},
  NumOps{parent.NumOps},
  snk(NumFiles),
  src(NumFiles),
  parity(NumFiles),
  SinhCoshAdjust(NumExponents)
{
  for( int f = 0; f < NumFiles; f++ )
  {
    auto &Corr{ parent.Corr[f] };
    parity[f] = Corr.parity;
    auto &Op{ Corr.Name_.op };
    snk[f].resize( NumExponents );
    src[f].resize( NumExponents );
    for( int e = 0; e < NumExponents; ++e )
    {
      snk[f][e] = parent.MELIndex( Op[idxSnk], e );
      src[f][e] = parent.MELIndex( Op[idxSrc], e );
    }
  }
}

void Model::Init( const scalar * ModelParameters, std::size_t NumParams, std::size_t Stride )
{
  Vector NewParams;
  NewParams.MapView( const_cast<scalar *>( ModelParameters ), NumParams, Stride );
  if( NewParams != Params )
  {
    Params = NewParams;
    Energy.MapView( &Params[0], NumExponents, Params.stride );
    Coeff.MapView( &Params[NumExponents], NumExponents * NumOps, Params.stride );
    for( int e = 0; e < NumExponents; ++e )
      SinhCoshAdjust[e] = std::exp( - Energy[e] * HalfNt );
    /*std::cout << "E[0]=" << Energy[0] << ", E[1]=" << Energy[1]
              << " A[0]=" << Coeff[0] << ", A[1]=" << Coeff[1]
              << " A[2]=" << Coeff[2] << ", A[3]=" << Coeff[3] << "\n";*/
  }
}

double Model::operator()( int f, int t ) const
{
  double z = 0;
  for( int e = 0; e < NumExponents; ++e )
  {
    double d;
    switch( parity[f] )
    {
      case Common::Parity::Even:
        d = SinhCoshAdjust[e] * std::cosh( - Energy[e] * ( t - HalfNt ) );
        break;
      case Common::Parity::Odd:
        d = SinhCoshAdjust[e] * std::sinh( - Energy[e] * ( t - HalfNt ) );
        break;
      default:
        d = std::exp( - Energy[e] * ( t ) );
        break;
    }
    z += d * Coeff[src[f][e]] * Coeff[snk[f][e]];
  }
  /*if( !std::isfinite( z ) )
  {
    std::cout << "Energy[0]=" << Energy[0] << ", Energy[1]=" << Energy[1] << "\n";
    std::cout << "Coeff[0]=" << Coeff[0] << ", Coeff[1]=" << Coeff[1]
              << "Coeff[2]=" << Coeff[2] << ", Coeff[3]=" << Coeff[3] << "\n";
    std::cout << "Debug here\n";
  }*/
  return z;
}

double Model::Derivative( int f, int t, int p ) const
{
  double d = 0;
  if( p < NumExponents )
  {
    // Derivative wrt Energy
    const int e{p};
    d = t * std::exp( - Energy[e] * t );
    if( parity[f] == Common::Parity::Even || parity[f] == Common::Parity::Odd )
    {
      int NtMinusT = parity[f] == Common::Parity::Even ? Nt - t : t - Nt;
      d += NtMinusT * std::exp( - Energy[e] * ( Nt - t ) );
      d *= 0.5;
    }
    d *= - Coeff[src[f][e]] * Coeff[snk[f][e]];
  }
  else
  {
    // Derivative wrt overlap coefficient
    p -= NumExponents;
    const int e{ p / NumOps };
    // Work out whether this operator is source or sink or both
    double OtherOverlap = 0;
    int Factor{ 0 };
    if( snk[f][e] == p )
    {
      OtherOverlap = Coeff[src[f][e]];
      ++Factor;
    }
    if( src[f][e] == p )
    {
      OtherOverlap = Coeff[snk[f][e]]; // Doesn't matter if we do this twice ... it's the same operator!
      ++Factor;
    }
    if( Factor )
    {
      switch( parity[f] )
      {
        case Common::Parity::Even:
          d = SinhCoshAdjust[e] * std::cosh( - Energy[e] * ( t - HalfNt ) );
          break;
        case Common::Parity::Odd:
          d = SinhCoshAdjust[e] * std::sinh( - Energy[e] * ( t - HalfNt ) );
          break;
        default:
          d = std::exp( - Energy[e] * t );
          break;
      }
      d *= Factor * OtherOverlap;
    }
  }
  return d;
}

FitterThread::FitterThread( const Fitter &fitter_, bool bCorrelated_, ModelFile &modelParams_, Correlators &CorrSynthetic_ )
: parent{ fitter_ },
  idx{ Fold::idxCentral },
  Error( fitter_.Extent ),
  SortingHat( fitter_.NumExponents ),
  bCorrelated{bCorrelated_},
  ModelParams{modelParams_},
  CorrSynthetic( CorrSynthetic_ ),
  model( fitter_ )
{
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
      CorrData = parent.Corr[f][idx];
    double z = ( model( f, t ) - CorrData[t] ) * CholeskyDiag[i];
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
      double z = model.Derivative( f, t, p ) * CholeskyDiag[i];
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
void FitterThread::MakeCovar( int idx_, bool bShowOutput )
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
  // Now make covariance matrix for replica idx_
  const int NumBoot{ parent.Corr[0].NumSamples() };
  // Make variance / covariance
  const scalar * CentralX = nullptr;
  const scalar * CentralY = nullptr;
  const scalar * ReplicaX = nullptr;
  const scalar * ReplicaY = nullptr;
  for( int x = 0; x < parent.Extent; x++ )
  {
    const int t1{ x % parent.NtCorr };
    if( !t1 )
    {
      CentralX = parent.Corr[x / parent.NtCorr][idx_] + parent.tMin;
      ReplicaX = parent.Corr[x / parent.NtCorr][ 0 ] + parent.tMin;
    }
    for( int y = bCorrelated ? 0 : x; y <= x; y++ )
    {
      const int t2{ y % parent.NtCorr };
      if( !t2 || !bCorrelated )
      {
        CentralY = parent.Corr[y / parent.NtCorr][idx_] + parent.tMin;
        ReplicaY = parent.Corr[y / parent.NtCorr][ 0 ] + parent.tMin;
      }
      double z = 0;
      const scalar * DataX{ ReplicaX };
      const scalar * DataY{ ReplicaY };
      for( int i = 0; i < NumBoot; i++ )
      {
        z += ( DataX[t1] - CentralX[t1] ) * ( DataY[t2] - CentralY[t2] );
        DataX += parent.Nt;
        DataY += parent.Nt;
      }
      z /= NumBoot;
      if( bCorrelated )
      {
        Covar( x, y ) = z;
        if( x != y )
          Covar( y, x ) = z;
      }
      else//if( x == y )
      {
        z = 1. / sqrt( z );
        CholeskyDiag[x] = z;
      }
      if( !std::isfinite( z ) )
        throw std::runtime_error( "Variance / covariance isn't finite" );
    }
  }
  if( bCorrelated )
  {
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
  {
    const std::string Sep{ " " };
    const std::string NewLine{ "\n" };
    std::ofstream s{ SaveCorMatFileName };
    s << "# Correlation matrix\n# Files: " << parent.NumFiles << "\n# NtCorr: " << parent.NtCorr
      << "\n# gnuplot: plot '" << SaveCorMatFileName
      << "' matrix columnheaders rowheaders with image pixels" << NewLine << parent.Extent;
    for( int f = 0; f < parent.NumFiles; f++ )
      for( int t = 0; t < parent.NtCorr; t++ )
        s << Sep << parent.OpNames[f] << t + parent.tMin;
    s << NewLine;
    for( int f = 0; f < parent.NumFiles; f++ )
      for( int t = 0; t < parent.NtCorr; t++ )
      {
        int i = f * parent.NtCorr + t;
        s << parent.OpNames[f] << t + parent.tMin;
        for( int j = 0; j < parent.Extent; j++ )
          s << Sep << Covar( i, j );
        s << NewLine;
      }
  }

  // Save the reconstructed correlator values for this replica
  model.Init( FitParams, parent.NumParams );
  for( int f = 0; f < parent.NumFiles; f++ )
  {
    scalar * mc{ CorrSynthetic[f][idx] };
    for( int t = 0; t < parent.Nt; t++ )
    {
      *mc++ = model( f, t );
    }
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
  model.Init( par.data(), par.size() );
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
                                  Correlators &CorrSynthetic_ )
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
  model.Init( x );
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
  model.Init( x );
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

Fitter::Fitter( FitterType fitType_, Correlators &&Corr_, const std::vector<std::string> &opNames_, bool bfactor_,
                const std::string &sopNameConcat_, int numExponents_, int verbosity_, const std::string &outputBaseName_,
                Common::SeedType seed_, bool bFreezeCovar_, bool bSaveCorr_, bool bSaveCMat_,
                int Retry_, int MaxIt_, double Tolerance_, double RelEnergySep_, int nSamples_, bool bAnalyticDerivatives_ )
  : fitType{fitType_},
    bAnalyticDerivatives{bAnalyticDerivatives_},
    NumOps{ static_cast<int>( opNames_.size() ) },
    OpNames{ opNames_ },
    bFactor{ bfactor_ },
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
    NumFiles{ static_cast<int>( Corr_.size() ) },
    NSamples{ ( nSamples_ <= 0 || nSamples_ > Corr_[0].NumSamples() ) ? Corr_[0].NumSamples() : nSamples_ },
    Nt      { Corr_[0].Nt() },
    NumExponents{ numExponents_ },
    NumParams{ NumExponents * ( 1 + NumOps ) },
    Corr{ std::move( Corr_ ) },
    ParamNames( MakeParamNames() )
{
  switch( fitType )
  {
    case FitterType::Minuit2:
    case FitterType::GSL:
      break;
    default:
      throw std::runtime_error( "Unknown FitterType " + std::to_string( static_cast<int>( fitType ) ) );
  }
}

std::vector<std::string> Fitter::MakeParamNames() const
{
  assert( NumExponents > 0 && "Number of exponents in the fit should be positive" );
  std::vector<std::string> Names( NumParams );
  for( int e = 0; e < NumExponents; e++ ) {
    Names[e] = "E" + std::to_string( e );
    for( int op = 0; op < NumOps; op++ )
      Names[NumExponents + MELIndex( op, e )] = OpNames[op] + std::to_string( e );
  }
  return Names;
}

// Perform a fit
std::vector<Common::ValWithEr<scalar>> Fitter::PerformFit( bool Bcorrelated_, int TMin_, int TMax_, double &ChiSq, int &dof_ )
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
  tMin = TMin_;
  tMax = TMax_;
  NtCorr = TMax_ - TMin_ + 1;
  Extent = NumFiles * NtCorr;
  dof = Extent - NumParams;
  if( dof <= 0 )
    throw std::runtime_error( "Fit from " + std::to_string( tMin ) + " to " + std::to_string( tMax ) + " has "
                              + std::to_string( dof ) + " degrees of freedom" );
  dof_ = dof;

  // Make somewhere to store the results of the fit for each bootstrap sample
  ModelFile ModelParams( OpNames, NumExponents, NumFiles, tMin, tMax, dof, bFactor, bFreezeCovar, NSamples, NumParams+1 );
  for( const Fold &f : Corr )
    ModelParams.FileList.emplace_back( f.Name_.Filename );
  ModelParams.CopyAttributes( Corr[0] );
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
    Correlators CorrSynthetic( NumFiles ); // correlators resulting from the fit params
    for( int f = 0; f < NumFiles; f++ )
    {
      CorrSynthetic[f].resize( NSamples, Nt );
      CorrSynthetic[f].FileList.push_back( ModelFileName );
      CorrSynthetic[f].CopyAttributes( Corr[f] );
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
      const int tGuess{ Nt / 4 };
      Energy[0] = 0;
      bool bOpDone[NumOps];
      for( int i = 0; i < NumOps; ++i )
        bOpDone[i] = false;
      for( int f = 0; f < NumFiles; f++ )
      {
        // Guess for energy is the log of the ratios of the exponentials at tGuess ... averaged over all correlators
        const double * const c{ Corr[f][Fold::idxCentral] };
        double E0 = std::log( c[tGuess] / c[tGuess + 1] );
        Energy[0] += E0;
        // Guess for matrix elements in two passes, first where source and sink same, then different
        for( int Pass = 0; Pass < 2; ++Pass )
        {
          const int iOpSrc{ MELIndex( Corr[f].Name_.op[idxSrc], 0 ) };
          const int iOpSnk{ MELIndex( Corr[f].Name_.op[idxSnk], 0 ) };
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
              std::cout << "Warning: can't distinguish " << OpNames[Corr[f].Name_.op[idxSnk]]
                        << " from " << OpNames[Corr[f].Name_.op[idxSrc]] << Common::NewLine;
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
              ft.reset( new FitterThreadMinuit2( *this, bCorrelated, ModelParams, CorrSynthetic ) );
              break;
            case FitterType::GSL:
              ft.reset( new FitterThreadGSL( *this, bCorrelated, ModelParams, CorrSynthetic ) );
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
                fitThread.MakeCovar( Fold::idxCentral, true );
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
          for( int idx = Fold::idxCentral; idx < NSamples; ++idx )
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
                fitThread.MakeCovar( idx );
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
    ModelParams.MakeCorrSummary( "Params" );
    ModelParams.Write( ModelFileName );
    //ModelParams.WriteSummary( Common::MakeFilename( sModelBase, Common::sModel, Seed, TEXT_EXT ) );
    for( int f = 0; f < NumFiles; f++ )
    {
      const int snk{ Corr[f].Name_.op[idxSnk] };
      const int src{ Corr[f].Name_.op[idxSrc] };
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
      const int NumExponents{ cl.SwitchValue<int>( "e" ) };
      const int NSamples{ cl.SwitchValue<int>("n") };
      //const std::string model{ opt.optionValue("m") };
      const bool bFactor{ cl.GotSwitch( "f" ) };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bFreezeCovar{ cl.GotSwitch( "freeze" ) };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };
      const bool bSaveCMat{ cl.GotSwitch("savecmat") };
      const FitterType fitType{ cl.GotSwitch("minuit2") ? FitterType::Minuit2 : FitterType::GSL };
      const bool bAnalyticDerivatives{ cl.GotSwitch("analytic") };

      if( Retry < 0 )
        throw std::invalid_argument( "Retry must be >= 0" );
      if( MaxIterations < 0 )
        throw std::invalid_argument( "MaxIterations must be >= 0" );

      bShowUsage = false;
      std::size_t i = 0;
      Common::SeedType Seed = 0;
      std::vector<std::string> OpNameFile;
      Correlators Corr;
      std::cout << std::setprecision( 13 /*std::numeric_limits<double>::max_digits10*/ ) << "Loading folded correlators\n";
      for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
      {
        Corr.emplace_back();
        Corr[i].Read( sFileName, "  ", &OpNameFile );
        if( i == 0 )
        {
          outBaseFileName.append( Corr[0].Name_.Base );
          Seed = Corr[0].Name_.Seed;
        }
        else
        {
          Corr[0].IsCompatible( Corr[i] );
        }
        i++;
      }
      std::vector<std::string> OpNames; // List of all the operator names referred to in the file
      if( bFactor )
        OpNames = OpNameFile;
      else
      {
        // OpNames become a list of all the Op_src or Op_snk combinations actually used
        OpNames.clear();
        for( Fold & f : Corr )
        {
          using It = const typename std::vector<std::string>::iterator;
          for( int i = 0; i < 2; i++ )
          {
            const std::string Op{ OpNameFile[f.Name_.op[i]] + "_" + pSrcSnk[i] };
            It it{ std::find( OpNames.begin(), OpNames.end(), Op ) };
            if( it == OpNames.end() )
            {
              f.Name_.op[i] = static_cast<int>( OpNames.size() );
              OpNames.emplace_back( std::move( Op ) );
            }
            else
              f.Name_.op[i] = static_cast<int>( it - OpNames.begin() );
          }
        }
      }
      std::sort( OpNameFile.begin(), OpNameFile.end() );
      std::string sOpNameConcat{ OpNameFile[0] };
      for( std::size_t i = 1; i < OpNameFile.size(); i++ )
      {
        sOpNameConcat.append( 1, '_' );
        sOpNameConcat.append( OpNameFile[i] );
      }
      Fitter m( fitType, std::move( Corr ), OpNames, bFactor, sOpNameConcat, NumExponents, Verbosity, outBaseFileName,
                Seed, bFreezeCovar, bSaveCorr, bSaveCMat, Retry, MaxIterations, Tolerance, RelEnergySep, NSamples,
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
    " <options> Bootstrap1 [Bootstrap2 ...]\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, where <options> are:\n"
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
