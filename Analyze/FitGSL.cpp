/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: FitGSL.cpp
 
 Copyright (C) 2021
 
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

#include "FitGSL.hpp"

void ParamStateGSL::StandardOut( std::ostream &os ) const
{
  os << parameters;
  //if( bValid )
  //{
    //os << "";
  //}
}

void ParamStateGSL::ReplicaMessage( std::ostream &os ) const
{
  if( bValid )
  {
    os << "Stop: "  <<(ConvergeReason == 1 ? "step size" : "gradient" )
       << ", f()="  << nevalf
       << ", df()=" << nevaldf << ", ";
  }
}

ParamState * FitterThreadGSL::MakeParamState( const Parameters &Params )
{
  return new ParamStateGSL( Params );
}

FitterThreadGSL::FitterThreadGSL( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_,
                                  vCorrelator &CorrSynthetic_ )
: FitterThread( fitter_, bCorrelated_, outputModel_, CorrSynthetic_ ), vGuess( parent.NumVariable )
{
  const FitterGSL &parentGSL{ *dynamic_cast<const FitterGSL*>( &parent ) };
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
  switch( parentGSL.trs )
  {
    case FitterGSL::TRS::lm:
      fdf_params.trs = gsl_multifit_nlinear_trs_lm;
      break;
    case FitterGSL::TRS::lmaccel:
      fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
      break;
    case FitterGSL::TRS::dogleg:
      fdf_params.trs = gsl_multifit_nlinear_trs_dogleg;
      break;
    case FitterGSL::TRS::ddogleg:
      fdf_params.trs = gsl_multifit_nlinear_trs_ddogleg;
      break;
    case FitterGSL::TRS::subspace2D:
      fdf_params.trs = gsl_multifit_nlinear_trs_subspace2D;
      break;
  }
  ws = gsl_multifit_nlinear_alloc( gsl_multifit_nlinear_trust, &fdf_params, fdf.n, fdf.p );
}

FitterThreadGSL::~FitterThreadGSL()
{
  if( ws )
    gsl_multifit_nlinear_free( ws );
}

void FitterThreadGSL::MakeCovarCorrelated()
{
  // Cholesky is Cholesky decomposition of (i.e. LL^T=) CORRELATION matrix
  // CholeskyInvert() finds the inverse correlation matrix
  Cholesky.CholeskyInvert();
  Dump( "Inverse Covar", Cholesky );
  // Cholesky() finds the Cholesky decomposition of the inverse correlation matrix
  Cholesky.Cholesky();
  Dump( "Inverse Covar Cholesky", Cholesky );
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
  if( !AnalyticJacobian( J ) )
    throw std::runtime_error( "Error computing Jacobian" );
  return 0;
}

void FitterThreadGSL::Minimise( ParamState &Guess_, int iNumGuesses )
{
  ParamStateGSL &Guess{ *dynamic_cast<ParamStateGSL*>( &Guess_ ) };
  std::size_t i = 0;
  for( const Parameters::Parameter & p : Guess.parameters )
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
  // Use the requested number of iterations ... or infinite if zero
  const std::size_t MaxIt{ parent.MaxIt ? parent.MaxIt : std::numeric_limits<std::size_t>::max() };
  gsl_multifit_nlinear_driver( MaxIt, tol, tol, tol, nullptr, nullptr, &Guess.ConvergeReason, ws );
  Guess.bValid = true;
  const std::size_t nIter{ gsl_multifit_nlinear_niter( ws ) };
  Guess.NumCalls = nIter < std::numeric_limits<unsigned int>::max() ? static_cast<unsigned int>( nIter ) : std::numeric_limits<unsigned int>::max();
  Guess.nevalf = fdf.nevalf;
  Guess.nevaldf = fdf.nevaldf;
  const Vector &vResidual{ * reinterpret_cast<Vector *>( gsl_multifit_nlinear_residual( ws ) ) };
  gsl_blas_ddot( &vResidual, &vResidual, &Guess.TestStat );
  i = 0;
  const Vector &vResult{ * reinterpret_cast<Vector *>( gsl_multifit_nlinear_position( ws ) ) };
  Matrix &mJacobian{ * reinterpret_cast<Matrix *>( gsl_multifit_nlinear_jac( ws ) ) };
  Matrix mErrors( parent.NumVariable, parent.NumVariable );
  gsl_multifit_nlinear_covar( &mJacobian, 0, &mErrors );
  for( Parameters::Parameter & p : Guess.parameters )
  {
    p.Error = std::sqrt( mErrors( i, i ) );
    p.Value = vResult[i++];
  }
  // TODO: Should I reinstate this? Or is it debugging?
  /* if( parent.Verbosity > 1 && !fdf.df && idx == Fold::idxCentral )
  {
    // Compare numeric derivatives to the analytic ones I would have computed
    Matrix MyJacobian( parent.ds.Extent, parent.NumVariable );
    AnalyticJacobian( MyJacobian );
    std::cout << Common::NewLine << "GSL Jacobian:\n" << mJacobian << Common::NewLine
              << "My Jacobian:\n" << MyJacobian << Common::NewLine;
    MyJacobian.cols(); // Debug breakpoint here
  }*/
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

const std::string & FitterGSL::Type() const
{
  static const std::string MyType{ "GSL" };
  return MyType;
}

FitterGSL * FitterGSL::Make(const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                            const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames)
{
  TRS trs{ TRS::lm };
  static const std::vector<std::string> Algorithms{ "lm", "lmaccel", "dogleg", "ddogleg", "subspace2d" };
  if( !FitterArgs.empty() )
  {
    int idx{ Common::IndexIgnoreCase( Algorithms, FitterArgs ) };
    if( idx >= Algorithms.size() )
      return nullptr;
    trs = static_cast<TRS>( idx );
  }
  return new FitterGSL( trs, cl, ds, ModelArgs, opNames );
}

Fitter * MakeFitterGSL( const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                        const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames )
{
  return FitterGSL::Make( FitterArgs, cl, ds, ModelArgs, opNames );
}
