/**
 
 OpenMP thread that will perform fit on each replica
 
 Source file: FitterThread.cpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "Fitter.hpp"
#include "FitterThread.hpp"

FitterThread::FitterThread( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_ )
: parent{ fitter_ },
  Extent{ fitter_.ds.Extent },
  idx{ FirstIndex }, // So we'll notice a change on the first call to SetReplica
  CholeskyDiag( Extent ),
  ModelParams( fitter_.mp.NumScalars( Param::Type::All ) ),
  FitterParams( fitter_.mp.NumScalars( Param::Type::Variable ) ),
  state( fitter_.mp.NumScalars( Param::Type::All ), fitter_.mp.NumScalars( Param::Type::Variable ) ),
  bCorrelated{ bCorrelated_ },
  OutputModel{ outputModel_ },
  ModelBuffer( fitter_.NumFiles )
{
  if( bCorrelated )
  {
    Covar.resize( Extent, Extent );
    Correl.resize( Extent, Extent );
    Cholesky.resize( Extent, Extent );
  }
  // Make a buffer for each model to use as a scratchpad
  for( int f = 0; f < parent.NumFiles; f++ )
    ModelBuffer[f].resize( parent.model[f]->GetScratchPadSize() );
}

void FitterThread::SaveStdError()
{
  for( int i = 0; i < Extent; ++i )
  {
    StdErrorMean[i] = std::sqrt( StdErrorMean[i] );
    CholeskyDiag[i] = 1. / StdErrorMean[i];
  }
  parent.Dump( idx, "StdErrorMean", StdErrorMean );
  parent.Dump( idx, "CholeskyDiag", CholeskyDiag );
}

/**
 Given: Covar
 1) Extract StdErrorMean = sqrt(diagonals) and CholeskyDiag = 1/StdErrorMean
 2) Create Correl = Covar with CholeskyDiag applied
 3) Cholesky decompose and invert Correl
   NB: If that fails, try pivoted Cholesky decomposition
 4) Cholesky decompose the inverse Correl
 Save/dump the matrices on the central replica
 */
void FitterThread::InitialiseCovar( const std::string *pBaseName )
{
  if( pBaseName && pBaseName->empty() )
    pBaseName = nullptr;
  const MLU::SeedType Seed{ OutputModel.Name_.Seed };
  // Extract StdErrorMean = sqrt(diagonals) and CholeskyDiag = 1/StdErrorMean
  for( std::size_t i = 0; i < Extent; ++i )
    StdErrorMean[i] = Covar( i, i );
  SaveStdError();
  // Create Correl = Covar with CholeskyDiag applied
  Correl = Covar;
  Correl.CholeskyScaleApply( CholeskyDiag, true );
  LedoitWolfShrink( Correl, parent.ShrinkFactor );
  parent.Dump( idx, "Covariance", Covar );
  parent.Dump( idx, "Correlation", Correl );
  if( idx == Fold::idxCentral )
  {
    OutputModel.Covar = Covar;
    OutputModel.Correl = Correl;
    if( pBaseName )
    {
      parent.SaveMatrixFile( Covar, MLU::sCovariance,
            MLU::MakeFilename( *pBaseName, MLU::sCovmat, Seed, TEXT_EXT ) );
      parent.SaveMatrixFile( Correl, MLU::sCorrelation,
            MLU::MakeFilename( *pBaseName, MLU::sCormat, Seed, TEXT_EXT ),
                                  MLU::pszCorrelGnuplot );
    }
  }
  // Cholesky decompose the correlation matrix
  // CholeskyDiag[i] = 1/sqrt( \sigma_{ii} ) i.e. inverse of square root of variances
  // Cholesky is generated by gsl_linalg_cholesky_decomp2
  //   So lower triangle and diagonal is cholesky decomposition of CORRELATION matrix
  //   Upper right triangle (excluding diagonal - assume 1) is CORRELATION matrix. NB: Not documentated
  Cholesky = Correl;
  Matrix CorrelInv( Cholesky.size1, Cholesky.size2 );
  double CondNumber = 1.;
  // Invert correlation matrix. Set:
  //  CondNumber  Condition number of correlation matrix
  //  Cholesky    Cholesky decomposition of correlation matrix
  //  CorelInv    Inverse correlation matrix
  if( Cholesky.Cholesky( true ) )
  {
    // Cholesky decomposition of correlation matrix worked
    if( idx == Fold::idxCentral )
      CondNumber = Cholesky.CholeskyRCond();
    // Cholesky is Cholesky decomposition of (i.e. LL^T=) CORRELATION matrix
    // CholeskyInvert() finds the inverse correlation matrix
    CorrelInv = Cholesky;
    CorrelInv.CholeskyInvert();
  }
  else
  {
    // Cholesky decomposition failed - try pivoted Cholesky decomposition
    // Dump eigenvalues
    Vector S( Correl.GetEigenValues() );
    std::cout << "Cholesky decomp of correlation matrix failed on ";
    if( idx == Fold::idxCentral )
      std::cout << "central replica";
    else
      std::cout << "replica " << idx;
    std::cout << ". Trying pivoted\n";
    if( parent.Dump( idx ) )
    {
      std::cout << "Correlation matrix Eigenvalues:\n";
      for( std::size_t i = 0; i < S.size; ++i )
        std::cout << '\t' << i << '\t' << S[i] << '\n';
    }
    // Try pivoted Cholesky decomposition
    Cholesky = Correl;
    gsl_permutation *Perm{ gsl_permutation_alloc( Cholesky.size1 ) };
    try
    {
      if( gsl_linalg_pcholesky_decomp2( &Cholesky, Perm, &S ) )
        throw std::runtime_error( "Unable to pivoted Cholesky decompose correlation matrix" );
      if( idx == Fold::idxCentral )
      {
        Vector CondBuffer( 3 * Cholesky.size1 );
        const int e2{ gsl_linalg_pcholesky_rcond( &Cholesky, Perm, &CondNumber, &CondBuffer ) };
        if( e2 )
          MLU::GSLLibraryGlobal::Error( "Unable to get condition number using pivoted Cholesky"
                                          " decomposition", __FILE__, __LINE__, e2 );
      }
      if( gsl_linalg_pcholesky_invert( &Cholesky, Perm, &CorrelInv ) )
        throw std::runtime_error( "Unable to pivoted Cholesky invert correlation matrix" );
    }
    catch(...)
    {
      gsl_permutation_free( Perm );
      throw;
    }
    gsl_permutation_free( Perm );
  }
  // Here we have:
  //  CondNumber  Condition number of correlation matrix
  //  Cholesky    Cholesky decomposition of correlation matrix
  //  CorelInv    Inverse correlation matrix
  if( idx == Fold::idxCentral )
  {
    // Show condition number
    const int CondDigits{ static_cast<int>( 0.5 - std::log10( CondNumber ) ) };
    if( CondDigits >= 12 )
      std::cout << "WARNING see https://www.gnu.org/software/gsl/doc/html/linalg.html#cholesky-decomposition\n";
    std::cout << "Covariance matrix condition number " << CondNumber
              << ", ~" << CondDigits << " digits\n";
  }
  if( !CorrelInv.IsFinite() )
    throw std::runtime_error( "Inversion of correlation matrix produced NaNs" );
  if( parent.Dump( idx ) )
  {
    parent.Dump( idx, "Correlation Cholesky", Cholesky );
    parent.Dump( idx, "Inverse Correl", CorrelInv );
    Matrix Copy{ Cholesky };
    Copy.blas_trmm( CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1, Cholesky );
    parent.Dump( idx, "L L^T = Correlation", Copy );
    Matrix One{ Copy.size1, Copy.size2 };
    One = 0.;
    One.blas_symm( CblasLeft, CblasLower, 1, Copy, CorrelInv, 1 );
    parent.Dump( idx, "L L^T C^{-1} = 1", One );
  }
  if( idx == Fold::idxCentral )
  {
    OutputModel.CorrelCholesky = Cholesky;
    OutputModel.CorrelInv = CorrelInv;
    OutputModel.CovarInv = CorrelInv;
    OutputModel.CovarInv.CholeskyScaleApply( CholeskyDiag );
    parent.Dump( idx, "Inverse Covar", OutputModel.CovarInv );
    if( pBaseName )
    {
      parent.SaveMatrixFile( Cholesky, MLU::sCorrelationCholesky,
         MLU::MakeFilename(*pBaseName, MLU::sCormatCholesky, Seed, TEXT_EXT));
      parent.SaveMatrixFile( OutputModel.CovarInv, MLU::sCovarianceInv,
         MLU::MakeFilename( *pBaseName, MLU::sCovmatInv, Seed,TEXT_EXT ) );
    }
  }
  // Set Cholesky = Cholesky decomposition of the inverse correlation matrix
  Cholesky = CorrelInv;
  Cholesky.Cholesky( true );
  parent.Dump( idx, "Inverse Correl Cholesky", Cholesky );
  if( idx == Fold::idxCentral )
  {
    OutputModel.CorrelInvCholesky = Cholesky;
    OutputModel.CovarInvCholesky = Cholesky;
    // Apply scale on the left, i.e. S L
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        OutputModel.CovarInvCholesky( i, j ) *= CholeskyDiag[i];
    if( pBaseName )
    {
      parent.SaveMatrixFile( Cholesky, MLU::sCorrelationInvCholesky,
                MLU::MakeFilename( *pBaseName, MLU::sCormatInvCholesky, Seed, TEXT_EXT ) );
      parent.SaveMatrixFile( OutputModel.CovarInvCholesky, MLU::sCovarianceInvCholesky,
                MLU::MakeFilename( *pBaseName, MLU::sCovmatInvCholesky, Seed, TEXT_EXT ) );
    }
  }
}

void FitterThread::SetReplicaVars( int idx_ )
{
  // Switch to the requested replica
  idx = idx_;
  state.bValid = false;
  // Now get the data for this replica
  parent.ds.mFitData.MapRow( Data, idx );
  OutputModel.ModelPrediction.MapRow( Theory, idx );
  OutputModel.ErrorScaled.MapRow( Error, idx );
  // StdErrorMean for every replica - unless frozen, in which case all use central replica
  if( idx == Fold::idxCentral || !parent.ds.mBinned.Empty() )
    OutputModel.StdErrorMean.MapRow( StdErrorMean, idx );
  // All replicas start from the same Guess
  ModelParams = parent.Guess;
  // Load the constants for this replica into the fixed portion of the guess
  parent.ds.GetFixed( idx, ModelParams, parent.ParamFixed );
  // Give the controller the chance to initialise on each replica
  parent.fitController.SetReplica( ModelParams );
  // Allow each model to initialise on each replica
  for( int f = 0; f < parent.NumFiles; ++f )
    if( ModelBuffer[f].size )
      parent.model[f]->SetReplica( ModelBuffer[f], ModelParams );
  // Export variable parameters from the guess for the fitter to start with
  parent.mp.AllToType( FitterParams, ModelParams, Param::Type::Variable );
}

void FitterThread::Initialise( const std::string *pBaseName )
{
  SetReplicaVars( Fold::idxCentral );
  if( parent.bAllParamsKnown )
  {
    // No point inverting covariance - we're not doing a fit
  }
  if( !parent.cp.Covar.Empty() )
  {
    // Use the inverse covariance matrix loaded from file
    if( parent.cp.Covar.size2 != Extent )
    {
      std::ostringstream os;
      os << "Loaded covariance matrix for " << parent.cp.Covar.size2
         << " data points but fit is for " << Extent << " data points";
      throw std::runtime_error( os.str().c_str() );
    }
    // Extract the diagonal (inverse of the variance^2) - ie turn into inverse correlation
    Cholesky = parent.cp.Covar;
    for( int i = 0; i < Extent; ++i )
    {
      CholeskyDiag[i] = std::sqrt( Cholesky( i, i ) );
      StdErrorMean[i] = 1. / CholeskyDiag[i];
    }
    Cholesky.CholeskyScaleApply( StdErrorMean, true );
    // Cholesky decompose
    Cholesky.Cholesky( true );
    parent.Dump( idx, "Inverse Correl Cholesky", Cholesky );
    // Save matrices
    OutputModel.CovarInv = parent.cp.Covar;
    OutputModel.CorrelInvCholesky = Cholesky;
    OutputModel.CovarInvCholesky = Cholesky;
    // Apply scale on the left, i.e. S L
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        OutputModel.CovarInvCholesky( i, j ) *= CholeskyDiag[i];
  }
  else if( parent.ds.mCovar.Empty() )
  {
    // Fully unfrozen. Make (co)variance on central replica
    Vector Mean;
    JackBoot::MakeMean( Mean, parent.ds.mBinned );
    if( bCorrelated )
    {
      JackBoot::MakeCovar( Covar, Mean, parent.ds.mBinned, JackBoot::Norm::RawBinned );
      InitialiseCovar( pBaseName );
    }
    else
    {
      JackBoot::MakeVar( StdErrorMean, Mean, parent.ds.mBinned, JackBoot::Norm::RawBinned );
      SaveStdError();
    }
  }
  else
  {
    // (semi-)frozen - Use the correlation matrix we've been given for the central replica
    if( bCorrelated )
    {
      Covar = parent.ds.mCovar;
      InitialiseCovar( pBaseName );
    }
    else
    {
      for( std::size_t i = 0; i < Extent; ++i )
        StdErrorMean[i] = parent.ds.mCovar( i, i );
      SaveStdError();
    }
  }
}

void FitterThread::SwitchReplica( int idx_ )
{
  if( idx_ == Fold::idxCentral )
    throw std::runtime_error( "FitterThread::SwitchReplica() called on central replica" );
  // Don't make the same covariance matrix twice in a row
  if( idx_ != idx )
  {
    SetReplicaVars( idx_ );
    // Speed optimisation - don't re-calculate covariance if we know all parameters
    if( !parent.ds.mBinned.Empty() && !parent.bAllParamsKnown )
    {
      // Unfrozen
      if( parent.ds.mCovar.Empty() )
      {
        // Fully unfrozen
        if( bCorrelated )
        {
          JackBoot::MakeCovar( Covar, parent.ds.mBinned, idx );
          InitialiseCovar();
        }
        else
        {
          JackBoot::MakeVar( StdErrorMean, parent.ds.mBinned, idx );
          SaveStdError();
        }
      }
      else
      {
        // Semi frozen - apply variance for this replica to inverse correl
        JackBoot::MakeVar( StdErrorMean, parent.ds.mBinned, idx );
        SaveStdError();
      }
    }
  }
}

std::string FitterThread::ReplicaString( int iFitNum ) const
{
  std::stringstream ss;
  ss << ( bCorrelated ? "C" : "Unc" ) << "orrelated fit " << iFitNum << " on ";
  if( idx == Fold::idxCentral )
    ss << "central replica";
  else
    ss << "replica " << idx;
  return ss.str();
}

void FitterThread::ShowReplicaMessage( int iFitNum ) const
{
  double ChiSq{ getTestStat() };
  std::cout << ReplicaString( iFitNum ) << ", calls " << getNumCalls() << ", chi^2 " << ChiSq << MLU::NewLine;
  ReplicaMessage( std::cout );
  std::cout << "dof " << parent.dof << ", chi^2/dof " << ( ChiSq / ( parent.dof ? parent.dof : 1 ) ) << MLU::NewLine;
  // TODO: I don't think this first case is correct?
  if( parent.Verbosity > 2 )
    DumpParamsFitter( std::cout );
  else
    parent.mp.Dump( std::cout, ModelParams, Param::Type::Variable, state.bValid ? &state.ModelErrors : nullptr );
}

// Compute Cholesky scaled (theory - data) based on parameters from the fitting engine
bool FitterThread::SaveError( Vector &Error, const scalar * FitterParams, std::size_t Size, std::size_t Stride )
{
  // Put the fitter parameters into our model parameters
  if( Size )
  {
    if( Size != parent.mp.NumScalars( Param::Type::Variable ) )
      throw std::runtime_error( "Fitter parameters don't match our variable parameters" );
    parent.mp.TypeToAll( ModelParams, VectorView( FitterParams, Size, Stride ),
                         Param::Type::Variable, true );
  }
  parent.mp.PropagateEnergy( ModelParams );
  // Give the controller the option to manipulate parameters
  parent.fitController.ComputeDerived( ModelParams );
  // Compute theory - error for each timeslice
  int i{0};
  bool bOK{ true };
  for( int f = 0; f < parent.NumFiles; ++f )
  {
    // Allow each model to cache any common computations based solely on the model parameters
    const Model &m{ *parent.model[f] };
    if( ModelBuffer[f].size )
      m.ModelParamsChanged( ModelBuffer[f], ModelParams );
    const vInt &FitTimes{ parent.ds.FitTimes[f] };
    for( int t : FitTimes )
    {
      Theory[i] = m( t, ModelBuffer[f], ModelParams );
      Error[i] = ( Theory[i] - Data[i] ) * CholeskyDiag[i];
      if( !std::isfinite( Error[i++] ) )
        bOK = false;
    }
  }
  if( bOK && bCorrelated )
  {
    Error.blas_trmv( CblasLower, CblasTrans, CblasNonUnit, Cholesky );
    bOK = Error.IsFinite();
  }
  return bOK;
}

bool FitterThread::AnalyticJacobian( Matrix &Jacobian ) const
{
  int i{0};
  for( int f = 0; f < parent.NumFiles; ++f )
  {
    const Model &m{ *parent.model[f] };
    const vInt &FitTimes{ parent.ds.FitTimes[f] };
    for( int t : FitTimes )
    {
      for( int p = 0; p < parent.mp.NumScalars( Param::Type::Variable ); ++p )
      {
        double z = m.Derivative( t, p ) * CholeskyDiag[i];
        if( !std::isfinite( z ) )
          return false;
        Jacobian( i, p ) = z;
      }
      ++i;
    }
  }
  if( bCorrelated )
  {
    parent.Dump( idx, "Jacobian", Jacobian );
    Jacobian.blas_trmm( CblasLeft, CblasLower, CblasTrans, CblasNonUnit, 1, Cholesky );
    parent.Dump( idx, "Jacobian Scaled", Jacobian );
  }
  return true;
}

scalar FitterThread::RepeatFit( int MaxGuesses )
{
  // Call the minimiser until it provides the same answer twice
  state.bValid = false;
  double dTestStat = -747;
  int iNumGuesses{ 0 };
  bool bFinished{ false };
  while( !bFinished && iNumGuesses++ <= MaxGuesses )
  {
    if( parent.bAllParamsKnown )
    {
      if( !SaveError( Error, nullptr, 0, 0 ) )
        throw std::runtime_error( "NaN values on replica " + std::to_string( idx ) );
      int gsl_e = gsl_blas_ddot( &Error, &Error, &dTestStat );
      if( gsl_e )
        GSLLG::Error( "FitterThread::RepeatFit() unable to compute test statistic"
                      " using gsl_blas_ddot()", __FILE__, __LINE__, gsl_e );
      bFinished = true;
      state.bValid = true;
      state.TestStat = dTestStat;
    }
    else
    {
      Minimise( iNumGuesses );
      double dNewTestStat{ getTestStat() };
      if( MaxGuesses == 0 )
        bFinished = true;
      else if( iNumGuesses != 1 )
        bFinished = ( dTestStat == dNewTestStat );
      dTestStat = dNewTestStat;
    }
    if( !bFinished && ( ( idx == Fold::idxCentral && parent.Verbosity ) || parent.Verbosity > 2 ) )
      ShowReplicaMessage( iNumGuesses );
  }
  if( !bFinished )
  {
    std::ostringstream es;
    es << "Fit on replica " << idx << " did not converge within " << MaxGuesses << " retries";
    throw std::runtime_error( es.str().c_str() );
  }
  if( !parent.bAllParamsKnown )
  {
    if( !SaveError( Error, FitterParams.data, FitterParams.size, FitterParams.stride ) )
      throw std::runtime_error( "NaN values on replica " + std::to_string( idx ) );
    parent.mp.TypeToAll<scalar>( state.ModelErrors, state.FitterErrors, Param::Type::Variable,
                                 true, &ModelParams );
  }
  if( idx == Fold::idxCentral || parent.Verbosity >= 2 )
    ShowReplicaMessage( iNumGuesses );
  return dTestStat;
}

// Perform an uncorrelated fit on the central timeslice to update parameters guessed
const Vector &FitterThread::UncorrelatedFit()
{
  bool bSaveCorrelated{ bCorrelated };
  bCorrelated = false;
  try
  {
    RepeatFit( NumRetriesGuess() );
  }
  catch(...)
  {
    bCorrelated = bSaveCorrelated;
    throw;
  }
  bCorrelated = bSaveCorrelated;
  return ModelParams;
}

// Perform a fit. NB: Only the fit on central replica updates ChiSq
scalar FitterThread::FitOne()
{
  // Perform fit
  scalar dTestStat = RepeatFit( NumRetriesFit() );
  // This is a convenience to make regression testing & comparison with errors generated by other codes easier
  //for( int i = 0; i < Extent; ++i )
    //Error[i] = std::abs( Error[i] );
  // Get test statistic
  scalar qValue, qValueH;
  scalar ChiSqPerDof{ dTestStat };
  const int FDist_p{ parent.dof };
  const int FDist_m{ parent.cp.CovarSampleSize() - 1 };
  if( parent.dof == 0 )
  {
    qValue = 1;
    qValueH = 1;
  }
  else
  {
    ChiSqPerDof /= parent.dof;
    qValue = MLU::qValueChiSq( dTestStat, parent.dof );
    if( MLU::HotellingDist::Usable( FDist_p, FDist_m ) )
      qValueH = MLU::HotellingDist::qValue( dTestStat, FDist_p, FDist_m );
    else
      qValueH = qValue;
  }

  // Make sure the test statistic is within acceptable bounds
  if( idx == Fold::idxCentral )
  {
    std::ostringstream ss;
    // Now add each component of the error vector
    ss << (bCorrelated ? "C" : "Unc") << "orrelated ( Theory - Data ) / sigma components:\n";
    std::size_t idxError = 0;
    for( std::size_t Corr = 0; Corr < parent.ds.FitTimes.size(); Corr++ )
    {
      ss << "C(" << Corr << ") " << parent.model[Corr]->Description() << ": ";
      for( std::size_t idxT = 0; idxT < parent.ds.FitTimes[Corr].size(); idxT++ )
      {
        if( idxT )
          ss << ", ";
        ss << '[' << parent.ds.FitTimes[Corr][idxT] << "]=" << Error[idxError++];
      }
      ss << "\n";
    }
    std::cout << ss.str();
    ss.str(std::string());
    ss << "chi^2/dof " << ChiSqPerDof << ", dof " << parent.dof << MLU::CommaSpace;
    if( parent.dof == 0 )
      ss << "Extrapolation => qValue 1";
    else
      ss << "Chi^2 qValue " << qValue;
    ss << ", p=" << FDist_p << ", m=" << FDist_m;
    if( parent.dof )
    {
      ss << ", Hotelling qValue ";
      if( MLU::HotellingDist::Usable( FDist_p, FDist_m ) )
        ss << qValueH;
      else
        ss << "invalid (m <= p)";
    }
    bool bOK{ true };
    if( parent.HotellingCutoff )
    {
      bOK = qValueH >= parent.HotellingCutoff;
      ss << ( bOK ? " >=" : " <" ) << " cutoff " << parent.HotellingCutoff;
    }
    if( parent.ChiSqDofCutoff )
    {
      const bool bOKC{ ChiSqPerDof <= parent.ChiSqDofCutoff };
      ss << ( bOKC ? " <=" : " >" ) << " cutoff " << parent.ChiSqDofCutoff;
      if( !bOKC )
        bOK = false;
    }
    if( !bOK )
      throw std::runtime_error( ss.str() );
    // See whether any of the energies is so large it's effectively indeterminate
    for( const Params::value_type &it : parent.mp )
    {
      const Param &p{ it.second };
      if( p.type == Param::Type::Variable && p.bMonotonic )
      {
        const std::size_t Offset{ p.GetOffset( 0, Param::Type::All ) };
        for( std::size_t i = 0; i < p.size; ++i )
        {
          if( ModelParams[Offset + i] > parent.MonotonicUpperLimit )
          {
            const Param::Key &pk{ it.first };
            ss << "\n" << pk.FullName( i, p.size ) << " > " << parent.MonotonicUpperLimit;
            throw std::runtime_error( ss.str() );
          }
        }
      }
    }
    std::cout << "OK: " << ss.str() << MLU::NewLine;
  }
  // Make pairs of parameters the right sign
  parent.mp.AdjustSigns( ModelParams );
  // Copy results into OutputModel
  for( typename Params::value_type it : parent.mp )
  {
    const Param::Key key{ it.first };
    const Param p{ it.second };
    const std::size_t Offset{ p.GetOffset( 0, Param::Type::All ) };
    for( std::size_t i = 0; i < p.size; ++i )
    {
      OutputModel(idx,Offset + i) = ModelParams[Offset + i];
      // Check whether energy levels are separated by minimum separation
      if( i && p.bMonotonic && idx == Fold::idxCentral && parent.RelEnergySep )
      {
        double RelSep = ModelParams[ Offset + i ] / ModelParams [ Offset + i - 1 ];
        if( ( RelSep - 1 ) < parent.RelEnergySep || RelSep * parent.RelEnergySep > 1 )
        {
          std::ostringstream es;
          es << "Fit failed energy separation criteria: " << ( 1 + parent.RelEnergySep )
             << " < " << key << i << " / " << key << i - 1 << " < " << ( 1 / parent.RelEnergySep );
          throw std::runtime_error( es.str().c_str() );
        }
      }
    }
  }
  std::size_t NumScalars{ parent.mp.NumScalars( Param::Type::All ) };
  OutputModel(idx,NumScalars++) = parent.dof > 1 ? ( dTestStat / parent.dof ) : dTestStat;
  OutputModel(idx,NumScalars++) = qValue;
  OutputModel(idx,NumScalars++) = qValueH;
  return dTestStat;
}
