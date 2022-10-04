/*************************************************************************************
 
 Use GSL as fitting exgine

 Source file: FitGSL.cpp
 
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

#include "FitGSL.hpp"

void FitterThreadGSL::DumpParamsFitter( std::ostream &os ) const
{
  parent.mp.Dump( os, ModelParams, Param::Type::Variable, state.bValid ? &state.ModelErrors : nullptr ); // TODO: Dump more of the GSL state
}

void FitterThreadGSL::ReplicaMessage( std::ostream &os ) const
{
  if( state.bValid )
  {
    os << "Stop: ";
    SayConvergeReason( os, ConvergeReason );
    os << ", f()="  << nevalf
       << ", df()=" << nevaldf << ", ";
  }
}

void FitterThreadGSL::InitialiseGSL()
{
  const FitterGSL &parentGSL{ *dynamic_cast<const FitterGSL*>( &parent ) };
  // Define my finite difference function
  std::memset( &fdf, 0, sizeof( fdf ) );
  /* define the function to be minimized */
  fdf.f = &sf;
  if( parent.bAnalyticDerivatives )
    fdf.df = &sdf; // Analytic derivatives
  fdf.n = Extent;
  fdf.p = parent.mp.NumScalars( Param::Type::Variable );
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
  if( !ws )
    throw std::runtime_error( "gsl_multifit_nlinear_alloc() failed" );
}

void FitterThreadGSL::SayConvergeReason( std::ostream &os, int ConvergeReason )
{
  switch( ConvergeReason )
  {
    case 1:
      os << "step size";
      break;
    case 2:
      os << "gradient";
      break;
    case 0:
      os << "none (iterations?)";
      break;
    case GSL_EMAXITER:
      os << "iterations";
      break;
    case GSL_ENOPROG:
      os << "no progress (no new step delta)";
      break;
    default:
      os << "unknown " << ConvergeReason;
  }
}

FitterThreadGSL::FitterThreadGSL( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_ )
: FitterThread( fitter_, bCorrelated_, outputModel_ )
{
  InitialiseGSL();
}

FitterThreadGSL::FitterThreadGSL( const FitterThreadGSL &ftGSL ) : FitterThread( ftGSL )
{
  InitialiseGSL();
}

FitterThreadGSL::~FitterThreadGSL()
{
  if( ws )
    gsl_multifit_nlinear_free( ws );
}

FitterThread * FitterThreadGSL::Clone() const
{
  return new FitterThreadGSL( *this );
}

int FitterThreadGSL::f( const Vector &FitParams, Vector &Error )
{
  assert( FitParams.size == parent.mp.NumScalars( Param::Type::Variable ) && "Parameter vector is not the right size" );
  assert( Error.size == Extent && "Result vector is not the right size" );
  if( !SaveError( Error, FitParams.data, FitParams.size, FitParams.stride ) )
    throw std::runtime_error( "Overflow computing residuals (usually caused by bad guess)" );
  return 0;
}

int FitterThreadGSL::df( const Vector &x, Matrix &J )
{
  assert( x.size == parent.mp.NumScalars( Param::Type::Variable ) && "Parameter vector is not the right size" );
  assert( J.size1 == Extent && "Jacobian rows != data points" );
  assert( J.size2 == parent.mp.NumScalars( Param::Type::Variable ) && "Parameter columns != parameters" );
  if( !AnalyticJacobian( J ) )
    throw std::runtime_error( "Error computing Jacobian" );
  return 0;
}

void FitterThreadGSL::Minimise( int )
{
  // First time around,  initialize solver with starting point and weights
  int gsl_e;
  if( !state.bValid )
  {
    gsl_e = gsl_multifit_nlinear_init( &FitterParams, &fdf, ws );
    if( gsl_e )
      GSLLG::Error( "FitterThreadGSL::Minimise() unable to start fit using gsl_multifit_nlinear_init()",
                   __FILE__, __LINE__, gsl_e );
  }

  // compute initial cost function
  if( idx == Fold::idxCentral )
  {
    double TestStat;
    const Vector &vResidual{ * reinterpret_cast<Vector *>( gsl_multifit_nlinear_residual( ws ) ) };
    gsl_e = gsl_blas_ddot( &vResidual, &vResidual, &TestStat );
    if( gsl_e )
      GSLLG::Error( "FitterThreadGSL::Minimise() unable to compute starting residual using gsl_blas_ddot()",
                   __FILE__, __LINE__, gsl_e );
    std::cout << (state.bValid ? "Intermediate" : "Guess") << " chi^2=" << TestStat << Common::NewLine;
  }

  // solve the system with a maximum of parent.MaxIt iterations
  const auto tol = parent.Tolerance;
  // Use the requested number of iterations ... or infinite if zero
  const std::size_t MaxIt{ parent.MaxIt ? parent.MaxIt : std::numeric_limits<std::size_t>::max() };
  gsl_e = gsl_multifit_nlinear_driver( MaxIt, tol, tol, tol, nullptr, nullptr, &ConvergeReason, ws );
  if( gsl_e )
  {
    std::ostringstream es;
    es << "Convergence reason ";
    SayConvergeReason( es, ConvergeReason );
    es << " during gsl_multifit_nlinear_driver()";
    GSLLG::Error( es.str().c_str(), __FILE__, __LINE__, gsl_e );
  }
  nevalf = fdf.nevalf;
  nevaldf = fdf.nevaldf;
  const std::size_t nIter{ gsl_multifit_nlinear_niter( ws ) };
  state.NumCalls = nIter < std::numeric_limits<unsigned int>::max()
                          ? static_cast<unsigned int>( nIter ) : std::numeric_limits<unsigned int>::max();
  FitterParams = * reinterpret_cast<Vector *>( gsl_multifit_nlinear_position( ws ) );
  f( FitterParams, Error );
  gsl_e = gsl_blas_ddot( &Error, &Error, &state.TestStat );
  GSLLG::Error( "FitterThreadGSL::Minimise() unable to compute residual using gsl_blas_ddot()",
               __FILE__, __LINE__, gsl_e );
  Matrix &mJacobian{ * reinterpret_cast<Matrix *>( gsl_multifit_nlinear_jac( ws ) ) };
  Matrix mErrors(parent.mp.NumScalars( Param::Type::Variable ), parent.mp.NumScalars(Param::Type::Variable));
  gsl_e = gsl_multifit_nlinear_covar( &mJacobian, 0, &mErrors );
  if( gsl_e )
    GSLLG::Error( "FitterThreadGSL::Minimise() gsl_multifit_nlinear_covar() failed",
                 __FILE__, __LINE__, gsl_e );
  for( std::size_t i = 0; i < parent.mp.NumScalars( Param::Type::Variable ); ++i )
    state.FitterErrors[i] = std::sqrt( mErrors( i, i ) );
  state.bValid = true;
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
                            std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                            CovarParams &&cp )
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
  return new FitterGSL( trs, cl, ds, ModelArgs, opNames, std::move( cp ) );
}

Fitter * MakeFitterGSL( const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                        std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                        CovarParams &&cp )
{
  return FitterGSL::Make( FitterArgs, cl, ds, ModelArgs, opNames, std::move( cp ) );
}
