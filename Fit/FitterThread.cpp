/*************************************************************************************
 
 OpenMP thread that will perform fit on each replica
 
 Source file: Fitter.cpp
 
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

#include "Param.hpp"
#include "Fitter.hpp"
#include "FitterThread.hpp"

FitterThread::FitterThread( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_, vCorrelator &CorrSynthetic_ )
: parent{ fitter_ },
  Extent{ fitter_.ds.Extent },
  idx{ FirstIndex }, // So we'll notice a change on the first call to SetReplica
  CholeskyDiag( Extent ),
  ModelParams( fitter_.mp.NumScalars( Param::Type::All ) ),
  FitterParams( fitter_.mp.NumScalars( Param::Type::Variable ) ),
  state( fitter_.mp.NumScalars( Param::Type::All ), fitter_.mp.NumScalars( Param::Type::Variable ) ),
  bCorrelated{ bCorrelated_ },
  OutputModel{ outputModel_ },
  CorrSynthetic( CorrSynthetic_ ),
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

// Make the covariance matrix, for only the timeslices we are interested in
void FitterThread::SetReplica( int idx_, bool bShowOutput, bool bSaveMatrices, const std::string *pBaseName )
{
  static const char pszCorrelGnuplot[] = "set cbrange[-1:1]";
  if( !bSaveMatrices && pBaseName )
    throw std::runtime_error( "Can only write matrices to file on central replica" );
  // Don't make the same covariance matrix twice in a row
  if( idx_ != idx )
  {
    // Switch to the requested replica
    const bool bFirstTime{ idx == FirstIndex };
    idx = idx_;
    // All replicas start from the same Guess
    ModelParams = parent.Guess;
    // Load the constants for this replica into the fixed portion of the guess
    parent.ds.GetFixed( idx, ModelParams, parent.ParamFixed );
    // Export variable parameters from the guess for the fitter to start with
    parent.mp.Export( FitterParams, ModelParams, Param::Type::Variable );
    state.bValid = false;
    // Now get the data for this replica
    parent.ds.GetData( idx, Data );
    OutputModel.ModelPrediction.MapView( Theory, idx );
    OutputModel.ErrorScaled.MapView( Error, idx );
    // Don't build covariance matrix if frozen (unless it's the first time)
    if( !parent.cp.bFreeze || bFirstTime )
    {
      OutputModel.StdErrorMean.MapView( StdErrorMean, idx ); // Each FitterThread needs it's own work area
      if( parent.cp.Covar.size1 == Extent )
      {
        // I've loaded the inverse covariance matrix already
        if( bFirstTime )
        {
          if( !CholeskyAdjust() )
            throw std::runtime_error( "Can't load inverse covariance matrix with this fitter" );
          Cholesky = parent.cp.Covar;
          for( int i = 0; i < Extent; ++i )
          {
            CholeskyDiag[i] = std::sqrt( Cholesky( i, i ) );
            StdErrorMean[i] = 1. / CholeskyDiag[i];
          }
          Cholesky.CholeskyScaleApply( StdErrorMean, true );
          Cholesky.Cholesky();
          Cholesky.ZeroUpperTriangle();
          parent.Dump( idx, "Inverse Correl Cholesky", Cholesky );
          if( bSaveMatrices )
          {
            OutputModel.CovarInv = parent.cp.Covar;
            OutputModel.CorrelInvCholesky = Cholesky;
            OutputModel.CovarInvCholesky = Cholesky;
            // Apply scale on the left, i.e. S L
            for( int i = 0; i < Extent; ++i )
              for( int j = 0; j <= i; ++j )
                OutputModel.CovarInvCholesky( i, j ) *= CholeskyDiag[i];
          }
        }
      }
      else
      {
        // (Co)variance always comes from the central replica if frozen
        if( parent.cp.bFreeze )
          idx_ = Fold::idxCentral;
        if( !bCorrelated || !parent.cp.OneStep( idx_ ) )
        {
          // Make the variance from the binned data on this replica
          parent.ds.MakeVariance( idx_, Common::SampleSource::Binned, StdErrorMean );
          SaveStdError();
        }
        if( bCorrelated )
        {
          if( parent.cp.OneStep( idx_ ) )
          {
            // I can construct the covariance matrix in one step
            if( parent.cp.Covar.size1 == Extent && parent.cp.bFreeze )
              Covar = parent.cp.Covar;
            else if( parent.cp.CovarNumBoot )
              parent.ds.MakeCovariance( idx_, parent.cp.Source, Covar, parent.cp.CovarRandom );
            else
              parent.ds.MakeCovariance( idx_, parent.cp.Source, Covar );
            // Correlated fit with correlation matrix from bootstrap. Extract inverse variance from diagonals
            for( int i = 0; i < Extent; ++i )
              StdErrorMean[i] = Covar( i, i );
            SaveStdError();
            Correl = Covar;
            Correl.CholeskyScaleApply( CholeskyDiag, true );
          }
          else
          {
            // Two step: correlation matrix from central replica + variance from binned data on each replica
            if( bFirstTime )
            {
              // Make the correlation matrix from the central replica
              if( parent.cp.Covar.size1 == Extent )
                Correl = parent.cp.Covar;
              else
                parent.ds.MakeCovariance( Fold::idxCentral, parent.cp.Source, Correl );
              parent.Dump( idx, "CovarianceIn", Correl );
              if( bSaveMatrices )
                OutputModel.CovarIn = Correl;
              if( pBaseName && !pBaseName->empty() )
                parent.SaveMatrixFile( Correl, Common::sCovarianceIn,
                    Common::MakeFilename( *pBaseName, Common::sCovmatIn, OutputModel.Name_.Seed, TEXT_EXT ) );
              // Make this a correlation matrix
              Correl.CholeskyExtract();
            }
            // Now apply the variance from the binned data on this replica
            Covar = Correl;
            Covar.CholeskyScaleApply( StdErrorMean );
          }
          parent.Dump( idx, "Covariance", Covar );
          parent.Dump( idx, "Correlation", Correl );
          if( bSaveMatrices )
          {
            OutputModel.Covar = Covar;
            OutputModel.Correl = Correl;
          }
          if( pBaseName && !pBaseName->empty() )
          {
            parent.SaveMatrixFile( Covar, Common::sCovariance,
                      Common::MakeFilename( *pBaseName, Common::sCovmat, OutputModel.Name_.Seed, TEXT_EXT ) );
            parent.SaveMatrixFile( Correl, Common::sCorrelation,
                      Common::MakeFilename( *pBaseName, Common::sCormat, OutputModel.Name_.Seed, TEXT_EXT ),
                                  pszCorrelGnuplot );
          }
          // Cholesky decompose the covariance matrix
          // CholeskyDiag[i] = 1/sqrt( \sigma_{ii} ) i.e. inverse of square root of variances
          // Cholesky is generated by gsl_linalg_cholesky_decomp2
          //   So lower triangle and diagonal is cholesky decomposition of CORRELATION matrix
          //   Upper right triangle (excluding diagonal - assume 1) is CORRELATION matrix. NB: Not documentated
          Cholesky = Correl;
          Cholesky.Cholesky();
          for( int i = 0; i < Cholesky.size1; ++i )
            for( int j = i + 1; j < Cholesky.size2; ++j )
              Cholesky( i, j ) = 0;
          parent.Dump( idx, "Cholesky", Cholesky );
          if( bSaveMatrices )
            OutputModel.CorrelCholesky = Cholesky;
          if( pBaseName && !pBaseName->empty() )
            parent.SaveMatrixFile( Cholesky, Common::sCorrelationCholesky,
                  Common::MakeFilename( *pBaseName, Common::sCormatCholesky, OutputModel.Name_.Seed, TEXT_EXT ) );
          if( parent.Dump( idx ) )
          {
            Matrix Copy{ Cholesky };
            Copy.blas_trmm( CblasRight, CblasLower, CblasTrans, CblasNonUnit, 1, Cholesky );
            parent.Dump( idx, "L L^T", Copy );
          }
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
          // Allow GSL or Minuit2 to use the cholesky decomposition of the inverse covariance matrix
          if( CholeskyAdjust() )
          {
            // Cholesky is Cholesky decomposition of (i.e. LL^T=) CORRELATION matrix
            // CholeskyInvert() finds the inverse correlation matrix
            Cholesky.CholeskyInvert();
            parent.Dump( idx, "Inverse Correl", Cholesky );
            if( bSaveMatrices )
            {
              OutputModel.CovarInv = Cholesky;
              OutputModel.CovarInv.CholeskyScaleApply( CholeskyDiag );
              parent.Dump( idx, "Inverse Covar", OutputModel.CovarInv );
              if( pBaseName && !pBaseName->empty() )
                parent.SaveMatrixFile( OutputModel.CovarInv, Common::sCovarianceInv,
                                      Common::MakeFilename( *pBaseName, Common::sCovmatInv, OutputModel.Name_.Seed,TEXT_EXT));
            }
            // Cholesky() finds the Cholesky decomposition of the inverse correlation matrix
            Cholesky.Cholesky();
            for( int i = 0; i < Extent; ++i )
              for( int j = i + 1; j < Extent; ++j )
                Cholesky( i, j ) = 0;
            parent.Dump( idx, "Inverse Correl Cholesky", Cholesky );
            if( bSaveMatrices )
            {
              OutputModel.CorrelInvCholesky = Cholesky;
              OutputModel.CovarInvCholesky = Cholesky;
              // Apply scale on the left, i.e. S L
              for( int i = 0; i < Extent; ++i )
                for( int j = 0; j <= i; ++j )
                  OutputModel.CovarInvCholesky( i, j ) *= CholeskyDiag[i];
              if( pBaseName && !pBaseName->empty() )
              {
                parent.SaveMatrixFile( Cholesky, Common::sCorrelationInvCholesky,
                Common::MakeFilename( *pBaseName, Common::sCormatInvCholesky, OutputModel.Name_.Seed, TEXT_EXT));
                parent.SaveMatrixFile( OutputModel.CovarInvCholesky, Common::sCovarianceInvCholesky,
                Common::MakeFilename( *pBaseName, Common::sCovmatInvCholesky, OutputModel.Name_.Seed, TEXT_EXT));
              }
            }
          }
        }
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
  std::cout << ReplicaString( iFitNum ) << ", calls " << getNumCalls() << ", chi^2 " << ChiSq << Common::NewLine;
  ReplicaMessage( std::cout );
  std::cout << "dof " << parent.dof << ", chi^2/dof " << ( ChiSq / ( parent.dof ? parent.dof : 1 ) ) << Common::NewLine;
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
  if( Size != parent.mp.NumScalars( Param::Type::Variable ) )
    throw std::runtime_error( "Fitter parameters don't match our variable parameters" );
  parent.mp.Import( ModelParams, VectorView( FitterParams, Size, Stride ) );
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
    Minimise( iNumGuesses );
    double dNewTestStat{ getTestStat() };
    if( MaxGuesses == 0 )
      bFinished = true;
    else if( iNumGuesses != 1 )
      bFinished = ( dTestStat == dNewTestStat );
    dTestStat = dNewTestStat;
    if( !bFinished && ( ( idx == Fold::idxCentral && parent.Verbosity ) || parent.Verbosity > 2 ) )
      ShowReplicaMessage( iNumGuesses );
  }
  if( !bFinished )
  {
    std::ostringstream es;
    es << "Fit on replica " << idx << " did not converge within " << MaxGuesses << " retries";
    throw std::runtime_error( es.str().c_str() );
  }
  if( !SaveError( Error, FitterParams.data, FitterParams.size, FitterParams.stride ) )
    throw std::runtime_error( "NaN values on replica " + std::to_string( idx ) );
  parent.mp.Import<scalar>( state.ModelErrors, state.FitterErrors, Param::Type::Variable, &ModelParams );
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
scalar FitterThread::FitOne( const std::vector<std::size_t> &Reorder )
{
  // Perform fit
  scalar dTestStat = RepeatFit( NumRetriesFit() );
  // This is a convenience to make regression testing & comparison with errors generated by other codes easier
  for( int i = 0; i < Extent; ++i )
    Error[i] = std::abs( Error[i] );
  // Make sure the test statistic is within acceptable bounds
  if( idx == Fold::idxCentral )
  {
    std::ostringstream ss;
    scalar qValue;
    const int FDist_p{ parent.dof };
    const int FDist_m{ parent.cp.CovarSampleSize() - 1 };
    ss << "p=" << FDist_p << ", m=" << FDist_m << ", ";
    if( parent.dof == 0 )
    {
      qValue = 1;
      ss << "Extrapolation";
    }
    else if( Common::HotellingDist::Usable( FDist_p, FDist_m ) )
    {
      qValue = Common::HotellingDist::qValue( dTestStat, FDist_p, FDist_m );
      ss << "Hotelling";
    }
    else
    {
      // Can't use Hotelling distribution. Use chi-squared distribution instead
      qValue = Common::qValueChiSq( dTestStat, parent.dof );
      ss << "m <= p, Chi^2";
    }
    ss << " qValue " << qValue;
    if( parent.HotellingCutoff )
    {
      const bool bOK{ qValue >= parent.HotellingCutoff };
      ss << ( bOK ? " >=" : " <" ) << " cutoff " << parent.HotellingCutoff;
      if( !bOK )
      {
        // Now add each component of the error vector
        ss << std::setprecision( std::numeric_limits<scalar>::max_digits10 ) << "\nTheory - Data ("
           << ( !bCorrelated || !CholeskyAdjust() ? "un" : "" ) << "correlated):";
        for( std::size_t i = 0; i < Error.size; ++i )
          ss << ( i ? Common::CommaSpace : Common::Space ) << Error[i];
        throw std::runtime_error( ss.str() );
      }
    }
    std::cout << "OK: " << ss.str() << Common::NewLine;
  }
  // Copy results into OutputModel
  scalar * const OutputData{ OutputModel[idx] };
  for( typename Params::value_type it : parent.mp )
  {
    const Param::Key key{ it.first };
    const Param p{ it.second };
    const std::size_t Offset{ p.GetOffset( 0, Param::Type::All ) };
    for( std::size_t i = 0; i < p.size; ++i )
    {
      OutputData[ Reorder[ Offset + i ] ] = ModelParams[ Offset + i ];
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
  OutputData[parent.mp.NumScalars( Param::Type::All )] = dTestStat / ( parent.dof ? parent.dof : 1 );

  // Save the reconstructed correlator values for this replica
  for( int f = 0; f < parent.NumFiles; ++f )
  {
    const Model &m{ *parent.model[f] };
    scalar * SyntheticData{ CorrSynthetic[f][idx] };
    for( int t = 0; t < CorrSynthetic[f].Nt(); ++t )
      *SyntheticData++ = m( t, ModelBuffer[f], ModelParams );
  }
  return dTestStat;
}
