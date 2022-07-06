/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: MultiFit.cpp
 
 Copyright (C) 2020-2021
 
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

// Uncomment this next line to debug without OpenMP
#define DEBUG_DISABLE_OMP

#ifdef HAVE_MINUIT2
// Comment out the next line to disable Minuit2 (for testing)
#define HAVE_MINUIT2_SO_BUILD_IT
#endif

// Indices for operators in correlator names
const char * pSrcSnk[] = { "src", "snk" };
const std::string E{ "E" };

std::ostream & operator<<( std::ostream &os, const Parameters &Params )
{
  for( const Parameters::Parameter &p : Params )
    os << std::string( Params.MaxLen() - p.Name.length() + 2, ' ' ) << p.Name
       << Common::Space << p.Value << "\t+/- " << p.Error << Common::NewLine;
  return os;
}

std::ostream & operator<<( std::ostream &os, const ParamState &State )
{
  State.StandardOut( os );
  return os;
}

// Work out where the covariance matrix comes from
CovarParams::CovarParams( const Common::CommandLine &cl, DataSet &dsrw ) : ds{ dsrw }
{
  // Are we using a frozen covariance matrix? Or does it vary for each sample?
  bFreeze = SupportsUnfrozen() ? cl.GotSwitch( "freeze" ) : true;

  // Decode the user's choice of covariance source
  std::string sCovOptions{ cl.SwitchValue<std::string>( "covsrc" ) };
  std::string sCovSrc{ Common::ExtractToSeparator( sCovOptions ) };
  if( Common::EqualIgnoreCase( sCovSrc, "Rebin" ) )
  {
    // Rebin the raw data. Overwrites the raw data in the correlators
    dsrw.Rebin( Common::ArrayFromString<int>( sCovOptions ) );
    RebinSize = ds.RebinSize;
    std::cout << "Rebinned raw data:";
    for( int i : RebinSize )
      std::cout << Common::Space << i;
    std::cout << Common::NewLine;
    Source = SS::Raw;
  }
  else if( Common::EqualIgnoreCase( sCovSrc, "h5" ) )
  {
    // Load inverse covariance matrix from hdf5 - for debugging and comparison with other fit code
    std::vector<std::string> Opts{ Common::ArrayFromString( sCovOptions ) };
    if( Opts.size() < 2 || Opts.size() > 3 )
      throw std::runtime_error( "Options should contain file[,group],dataset. Bad options: " + sCovOptions );
    std::string sRoot( "/" );
    const int iHaveGroupName{ Opts.size() == 2 ? 0 : 1 };
    std::string &GroupName{ iHaveGroupName ? Opts[1] : sRoot };
    std::string &DataSetName{ Opts[1 + iHaveGroupName] };
    ::H5::H5File f;
    ::H5::Group  g;
    Common::H5::OpenFileGroup( f, g, Opts[0], "Loading covariance matrix from ", &GroupName );
    try
    {
      Common::H5::ReadMatrix( g, DataSetName, Covar );
    }
    catch(const ::H5::Exception &)
    {
      ::H5::Exception::clearErrorStack();
      throw std::runtime_error( "Unable to load covariance matrix " + sCovOptions );
    }
    Source = SS::Bootstrap;
    bFreeze = true; // because we only have inv_cov. TODO: invert inv_cov if unfrozen? (Better to load cov)
  }
  else
  {
    // See what the user has asked for
    Source = Common::FromString<SS>( sCovSrc );
    if( !sCovOptions.empty() )
      throw std::runtime_error( "Covariance source " + sCovSrc + " unexpected parameters: " + sCovOptions );
  }

  // Check that all the input files match for raw/rebinned (checked for binned data on load)
  if( Source == SS::Raw || Source == SS::Binned )
  {
    // Using raw samples - check all files have the same number
    for( std::size_t i = 0; i < ds.corr.size(); ++i )
    {
      const int Count{ ds.corr[i].NumSamples( Source ) };
      if( Count == 0 || Count != ds.corr[0].NumSamples( Source ) )
      {
        std::ostringstream os;
        os << "Can't use " << Source << " samples for covariance. "
           << ds.corr[0].Name_.Filename << " has " << ds.corr[0].NumSamples( Source )
        << " samples but " << ds.corr[i].Name_.Filename << " has " << ds.corr[i].NumSamples( Source );
        throw std::runtime_error( os.str().c_str() );
      }
    }
  }

  // Has the user requested we get covariance using a bootstrap
  if( !cl.GotSwitch( "covboot" ) )
    CovarNumBoot = 0;
  else
  {
    CovarNumBoot = cl.SwitchValue<int>( "covboot" );
    if( CovarNumBoot < 0 || CovarNumBoot == 1 )
      throw std::runtime_error( "Rebootstrapping with " + std::to_string( CovarNumBoot ) + " replicas" );
    if( SourceIsBootstrap() )
    {
      std::ostringstream os;
      os << "Can't use ";
      if( Source == SS::Raw )
        os << "bootstrapped ";
      os << Source << " samples when bootstrapping covariance matrix";
      throw std::runtime_error( os.str().c_str() );
    }
    // 0 => Use the maximum number of samples available
    if( CovarNumBoot == 0 || CovarNumBoot > ds.MaxSamples )
      CovarNumBoot = ds.MaxSamples;
    // If I'm bootstrapping, I will need some random numbers
    if(   ds.corr[0].RandNum() // Do I have original bootstrap random numbers
       && ds.corr[0].SampleSize == CovarCount() // Are the random numbers over same range (0...SampleSize)
       && CovarNumBoot <= ds.corr[0].NumSamples() ) // Are there enough replicas available
    {
      // Re-use existing random numbers
      vCovarRandom.clear();
      CovarRandom.Map( ds.corr[0].RandNum(), CovarNumBoot, CovarCount() );
    }
    else
    {
      Common::GenerateRandom( vCovarRandom, ds.corr[0].Seed_, CovarNumBoot, CovarCount() );
      CovarRandom.Map( vCovarRandom.data(), CovarNumBoot, CovarCount() );
    }
  }

  // Make sure the random numbers I'll need are built once and only once
  for( int i = 0; i < 2; ++i )
  {
    const SS ss{ i == 0 ? SS::Raw : SS::Binned };
    if( Source == ss || ( i == 1 && !bFreeze ) )
    {
      if( dsrw.InitRandomNumbers( ss ) )
      {
        std::stringstream os;
        os << "Generated bootstrap random numbers for " << ss << " data using Mersenne Twister (C++ 11 std::mt19937) with seed " << ds.corr[0].Name_.Seed;
        std::cout << os.str() << Common::NewLine;
      }
    }
  }
}

std::ostream & operator<<( std::ostream &os, const CovarParams &cp )
{
  if( cp.bFreeze )
    os << "F";
  else
    os << "Unf";
  os << "rozen covariance matrix constructed from " << cp.CovarCount() << Common::Space;
  if( cp.CovarCount() != cp.CovarSampleSize() )
    os << "(" << cp.CovarSampleSize() << " independent) ";
  if( cp.SourceIsBootstrap() )
    os << "pre-bootstrapped ";
  os << cp.Source << " samples";
  return os;
}

FitterThread::FitterThread( const Fitter &fitter_, bool bCorrelated_, ModelFile &outputModel_, vCorrelator &CorrSynthetic_ )
: parent{ fitter_ },
  Extent{ fitter_.ds.Extent },
  idx{ FirstIndex }, // So we'll notice a change on the first call to SetReplica
  CholeskyDiag( Extent ),
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
  // TODO: Correlated & Uncorrelated happen on different threads in MPI
  // so same matrix built twice on central replica
  if( idx_ != idx )
  {
    // Switch to the requested replica
    const bool bFirstTime{ idx == FirstIndex };
    idx = idx_;
    // Now get the data for this replica
    if( idx == Fold::idxCentral )
      Data = parent.ds.vCentral;
    else
      Data.MapRow( parent.ds.mBoot, idx);
    OutputModel.ModelPrediction.MapView( Theory, idx );
    OutputModel.ErrorScaled.MapView( Error, idx );
    parent.ds.GetData( idx, Data );
    // And get associated constants
    parent.ds.GetFixed( idx, ModelParams, parent.ParamFixed );
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

void FitterThread::ReplicaMessage( const ParamState &state, int iFitNum ) const
{
  double ChiSq{ state.getTestStat() };
  std::cout << ReplicaString( iFitNum ) << ", calls " << state.getNumCalls() << ", chi^2 " << ChiSq << Common::NewLine;
  state.ReplicaMessage( std::cout );
  std::cout << "dof " << parent.dof << ", chi^2/dof " << ( ChiSq / ( parent.dof ? parent.dof : 1 ) ) << Common::NewLine;
  if( parent.Verbosity > 2 )
    std::cout << state;
  else
    std::cout << state.parameters;
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
      Theory[i] = m( t, ModelBuffer[f], ModelParams );
      Error[i] = ( Theory[i] - Data[i] ) * CholeskyDiag[i];
      if( !std::isfinite( Error[i++] ) )
        return false;
    }
  }
  return true;
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
      for( int p = 0; p < parent.NumVariable; ++p )
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

scalar FitterThread::RepeatFit( ParamState &Guess, int MaxGuesses )
{
  // Call the minimiser until it provides the same answer twice
  double dTestStat = -747;
  int iNumGuesses{ 0 };
  bool bFinished{ false };
  while( !bFinished && iNumGuesses++ <= MaxGuesses )
  {
    Minimise( Guess, iNumGuesses );
    double dNewTestStat{ Guess.getTestStat() };
    if( MaxGuesses == 0 )
      bFinished = true;
    else if( iNumGuesses != 1 )
      bFinished = ( dTestStat == dNewTestStat );
    dTestStat = dNewTestStat;
    if( ( idx == Fold::idxCentral && ( bFinished || parent.Verbosity ) )
       || ( parent.Verbosity >= 2 && bFinished ) || parent.Verbosity > 2 )
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
  std::unique_ptr<ParamState> ThisGuess( MakeParamState( Guess ) );
  try
  {
    RepeatFit( *ThisGuess, NumRetriesGuess() );
  }
  catch(...)
  {
    bCorrelated = bSaveCorrelated;
    throw;
  }
  bCorrelated = bSaveCorrelated;
  Guess = ThisGuess->parameters;
}

// Perform a fit. NB: Only the fit on central replica updates ChiSq
scalar FitterThread::FitOne( const Parameters &parGuess )
{
  // Perform fit
  std::unique_ptr<ParamState> Result( MakeParamState( parGuess ) );
  scalar dTestStat = RepeatFit( *Result, NumRetriesFit() );
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
    if( Common::HotellingDist::Usable( FDist_p, FDist_m ) )
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
  // Put the variable fit parameters into the full model parameter set (constants are already there)
  for( int i = 0; i < parent.NumVariable; ++i )
    ModelParams[parent.ParamVariable[i]] = Result->parameters[i].Value;
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
  OutputData[OutputDataSize++] = dTestStat / ( parent.dof ? parent.dof : 1 );

  // Check whether energy levels are separated by minimum separation
  if( idx == Fold::idxCentral && parent.RelEnergySep )
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

Fitter::Fitter( const Common::CommandLine &cl, const DataSet &ds_,
                const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames_,
                CovarParams &&cp_ )
  : bAnalyticDerivatives{ cl.GotSwitch("analytic") },
    HotellingCutoff{ cl.SwitchValue<double>( "Hotelling" ) },
    RelEnergySep{ cl.SwitchValue<double>("sep") },
    MinDof{ cl.SwitchValue<int>("mindof") },
    Retry{ cl.SwitchValue<int>("retry") },
    MaxIt{ cl.SwitchValue<int>("iter") },
    Tolerance{ cl.SwitchValue<double>("tol") },
    bSaveCorr{ cl.GotSwitch("savecorr") },
    bSaveCMat{ cl.GotSwitch("savecmat") },
    Verbosity{ cl.SwitchValue<int>("v") },
    bForceSrcSnkDifferent{ cl.GotSwitch( "srcsnk" ) },
    vGuess{ Common::ArrayFromString<scalar>( cl.SwitchValue<std::string>( "guess" ) ) },
    ds{ std::move( ds_ ) },
    NumFiles{ static_cast<int>( ds.corr.size() ) },
    OpNames{ opNames_ },
    NumOps{ static_cast<int>( OpNames.size() ) },
    model{ CreateModels( ModelArgs, ModelDefaultParams( bForceSrcSnkDifferent, cl ) ) },
    NumExponents{ GetNumExponents() },
    PerExpNames{ MakePerExpNames() },
    ParamNames( MakeParamNames() ),
    ParamFixed{ MakeParamFixed() },
    ParamVariable{ MakeParamVariable() },
    NumModelParams{ static_cast<int>( ParamNames.size() ) },
    NumPerExp{ static_cast<int>( PerExpNames.size() ) },
    NumOneOff{ NumModelParams - NumExponents * NumPerExp },
    NumFixed{ static_cast<int>( ParamFixed.size() ) },
    NumVariable{ static_cast<int>( ParamVariable.size() ) },
    cp{ std::move( cp_ ) }
{
  assert( ds.corr.size() == NumFiles && "Number of files extremely large" );
  assert( NumModelParams == NumFixed + NumVariable && "NumModelParams doesn't match fixed and variable" );
  assert( NumModelParams == NumExponents * NumPerExp + NumOneOff && "NumModelParams doesn't match NumPerExp and NumOneOff" );
  if( MinDof < 0 )
    throw std::invalid_argument( "Degrees of freedom (mindof) must be >= 0" );
  if( Retry < 0 )
    throw std::invalid_argument( "Number of fit retries (retry) must be >= 0" );
  if( MaxIt < 0 )
    throw std::invalid_argument( "Maximum fitter iterations (iter) must be >= 0" );
  if( vGuess.size() && vGuess.size() != NumVariable )
    throw std::invalid_argument( "Guess contains " + std::to_string( vGuess.size() ) + " parameters, but "
                                + std::to_string( NumVariable ) + " variable parameters" );
}

// Create all the models and return the number of exponents in the fit
std::vector<ModelPtr> Fitter::CreateModels( const std::vector<std::string> &ModelArgs,
                                            const ModelDefaultParams &modelDefault )
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
  std::cout << "Making models\n";
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

void Fitter::SaveMatrixFile( const Matrix &m, const std::string &Type, const std::string &Filename,
                               const char *pGnuplotExtra ) const
{
  if( model.size() != ds.corr.size() )
    throw std::runtime_error( "ModelSet doesn't match DataSet" );
  std::vector<std::string> Abbreviations( ds.corr.size() );
  std::vector<std::string> FileComments;
  FileComments.reserve( ds.corr.size() );
  for( std::size_t f = 0; f < model.size(); ++f )
  {
    // Name of each operator THIS ONLY WORKS FOR MODELS WITH OVERLAP CONSTANTS - NEEDS TO BE FIXED
    for( std::size_t op = 1; op < model[f]->ParamNamesPerExp.size(); ++op )
      Abbreviations[f].append( model[f]->ParamNamesPerExp[op] );
    if( model[f]->ParamNamesPerExp.size() == 2 )
      Abbreviations[f].append( model[f]->ParamNamesPerExp[1] );
    // Operators - one-off
    std::ostringstream s;
    if( !model[f]->ParamNames.empty() )
    {
      for( std::size_t i = 0; i < model[f]->ParamNames.size(); ++i )
      {
        if( i == 0 )
          s << "# OpsOneOff" << f << ": ";
        else
          s << Common::CommaSpace;
        s << model[f]->ParamNames[i];
      }
      s << Common::NewLine;
    }
    // Operators per energy
    if( !model[f]->ParamNamesPerExp.empty() )
    {
      for( std::size_t i = 0; i < model[f]->ParamNamesPerExp.size(); ++i )
      {
        if( i == 0 )
          s << "# OpsPerExp" << f << ": ";
        else
          s << Common::CommaSpace;
        s << model[f]->ParamNamesPerExp[i];
      }
      s << Common::NewLine;
    }
    FileComments.emplace_back( s.str() );
  }
  // If we're saving the correlation matrix, it should have a similar name to the model
  ds.SaveMatrixFile( m, Type, Filename, Abbreviations, &FileComments, pGnuplotExtra );
}

// Perform a fit - assume fit ranges have been set on the DataSet prior to the call
std::vector<Common::ValWithEr<scalar>>
Fitter::PerformFit( bool Bcorrelated, double &ChiSq, int &dof_, const std::string &OutBaseName,
                    const std::string &ModelSuffix, Common::SeedType Seed )
{
  bCorrelated = Bcorrelated;
  dof = ds.Extent - NumVariable;
  if( dof < MinDof )
  {
    std::string Message{};
    if( dof )
      Message = "Fit has " + std::to_string( dof ) + " degrees of freedom";
    else
      Message = "Fit is an extrapolation (0 degrees of freedom)";
    throw std::runtime_error( Message );
  }
  dof_ = dof;

  // Make somewhere to store the results of the fit for each bootstrap sample
  const int tMin{ ds.FitTimes[0][0] };
  const int tMax{ ds.FitTimes[0].back() };
  ModelFile OutputModel( OpNames, NumExponents, NumFiles, tMin, tMax, dof, !bForceSrcSnkDifferent,
                         cp.bFreeze, ds.NSamples, NumModelParams + 1 );
  for( const Fold &f : ds.corr )
    OutputModel.FileList.emplace_back( f.Name_.Filename );
  OutputModel.CopyAttributes( ds.corr[0] );
  {
    std::vector<std::string> ColNames{ ParamNames };
    //TODO: Work out how to copy parameter names into output file properly
    if( ColNames.size() == 1 )
      OutputModel.OpNames.clear();
    ColNames.push_back( "ChiSqPerDof" );
    OutputModel.SetColumnNames( ColNames );
  }
  OutputModel.Name_.Seed = Seed;
  OutputModel.binSize = ds.OriginalBinSize;
  OutputModel.CovarFrozen = cp.bFreeze;
  OutputModel.CovarSource = cp.Source;
  OutputModel.CovarRebin = cp.RebinSize;
  OutputModel.CovarNumBoot = cp.CovarNumBoot;
  OutputModel.CovarSampleSize = cp.CovarSampleSize();
  OutputModel.StdErrorMean.resize( ds.Extent, ds.NSamples ); // Each thread needs to use its own replica
  OutputModel.ModelPrediction.resize( ds.Extent, ds.NSamples );
  OutputModel.ErrorScaled.resize( ds.Extent, ds.NSamples );

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
    bPerformFit = PreBuilt.NewParamsMorePrecise( cp.bFreeze, ds.NSamples );
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
      int iVar = 0;
      for( int i = 0; i < NumModelParams; ++i )
      {
        assert( bKnown[i] && "Bad model: unable to guess all parameters" );
        if( !bFixed[i] )
          parGuess.Add( ParamNames[i], vGuess.empty() ? modelParams[i] : vGuess[iVar++], 0 );
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
          std::unique_ptr<FitterThread> ft{ MakeThread( bCorrelated, OutputModel, CorrSynthetic ) };
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
                fitThread.SetReplica( Fold::idxCentral, true, true, bSaveCMat ? &sModelBase : nullptr );
                if( vGuess.empty() )
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
                fitThread.SetReplica( idx );
                scalar z{ fitThread.FitOne( parGuess ) };
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
    if( cp.bFreeze )
      OutputModel.StdErrorMean.Replica.clear();
    OutputModel.FitInput.Central = ds.vCentral;
    OutputModel.FitInput.Replica = ds.Cache( cp.Source );
    OutputModel.MakeCorrSummary( "Params" );
    {
      // Show parameters - in neat columns
      std::cout << OutputModel.GetSummaryNames()[0] << " after " << OutputModel.NumSamples()
                << " bootstrap replicas" << Common::NewLine;
      const Common::ValWithEr<scalar> * const v{ OutputModel.getSummaryData() };
      const std::vector<std::string> &Cols{ OutputModel.GetColumnNames() };
      const int maxLen{ static_cast<int>( std::max_element( Cols.begin(), Cols.end(),
                                          [](const std::string &a, const std::string &b)
                                          { return a.length() < b.length(); } )->length() ) };
      const int NumWidth{ static_cast<int>( std::cout.precision() ) + 7 };
      for( int i = 0; i < OutputModel.Nt(); ++i )
      {
        std::cout << std::setw(maxLen) << Cols[i] << std::left
                  << Common::Space << std::setw(NumWidth) << v[i].Central
                  << " +"    << std::setw(NumWidth) << (v[i].High - v[i].Central)
                  << " -"    << std::setw(NumWidth) << (v[i].Central - v[i].Low)
                  << " Max " << std::setw(NumWidth) << v[i].Max
                  << " Min "                        << v[i].Min
                  << Common::NewLine << std::right; // Right-aligned is the default
      }
    }
    // Save the file
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

// This should be the only place which knows about different fitters

Fitter * MakeFitterGSL( const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                        const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames,
                        CovarParams &&cp );
#ifdef HAVE_MINUIT2_SO_BUILD_IT
Fitter * MakeFitterMinuit2(const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                           const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames,
                           CovarParams &&cp );
#endif

Fitter * MakeFitter( const Common::CommandLine &cl, const DataSet &ds,
                     const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames,
                     CovarParams &&cp )
{
  Fitter * f;
  std::string FitterArgs{ cl.SwitchValue<std::string>( "fitter" ) };
  std::string FitterType{ Common::ExtractToSeparator( FitterArgs ) };
  if( Common::EqualIgnoreCase( FitterType, "GSL" ) )
    f = MakeFitterGSL( FitterArgs, cl, ds, ModelArgs, opNames, std::move( cp ) );
#ifdef HAVE_MINUIT2_SO_BUILD_IT
  else if( Common::EqualIgnoreCase( FitterType, "Minuit2" ) )
    f = MakeFitterMinuit2( FitterArgs, cl, ds, ModelArgs, opNames, std::move( cp ) );
#endif
  else
    throw std::runtime_error( "Unrecognised fitter: " + FitterType );
  if( !f )
    throw std::runtime_error( "Unrecognised " + FitterType + " fitter options: " + FitterArgs );
  return f;
}

int main(int argc, const char *argv[])
{
#ifndef HAVE_MINUIT2_SO_BUILD_IT
  std::ios_base::sync_with_stdio( false );
#endif
  static const char DefaultEnergySep[] = "0"; // Default used to be 0.2 until 10 Jul 2021
  static const char DefaultHotelling[] = "0.05";
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"analytic", CL::SwitchType::Flag, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"mindof", CL::SwitchType::Single, "1"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"savecmat", CL::SwitchType::Flag, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      {"v", CL::SwitchType::Single, "0"},
      {"srcsnk", CL::SwitchType::Flag, nullptr},
      {"guess", CL::SwitchType::Single, ""},
      // ModelDefaultParams
      {"e", CL::SwitchType::Single, "1"},
      // Other params
      {"fitter", CL::SwitchType::Single, "GSL"},
      {"delta", CL::SwitchType::Single, "3"},
      {"t", CL::SwitchType::Single, nullptr},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"opnames", CL::SwitchType::Flag, nullptr},
      {"covsrc", CL::SwitchType::Single, "Bootstrap"},
      {"covboot", CL::SwitchType::Single, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      const int delta{ cl.SwitchValue<int>("delta") };
      const std::string inBase{ cl.SwitchValue<std::string>("i") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      Common::MakeAncestorDirs( outBaseFileName );
      const int NSamples{ cl.SwitchValue<int>("n") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };
      const bool bOpSort{ !cl.GotSwitch("opnames") };

      if( !cl.GotSwitch( "t") )
        throw std::runtime_error( "No fit ranges specified" );
      Common::FitRanges fitRanges( Common::ArrayFromString( cl.SwitchValue<std::string>( "t" ) ) );

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
            std::string PrintPrefix( 2, ' ' );
            if( !cl.Args[ArgNum].empty() )
            {
              PrintPrefix.append( cl.Args[ArgNum] );
              PrintPrefix.append( 1, ' ' );
            }
            ds.LoadCorrelator( std::move( Att ), Common::COMPAT_DEFAULT, PrintPrefix.c_str() );
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

      // Work out how covariance matrix should be constructed and tell user
      CovarParams cp{ cl, ds };
      std::cout << cp << Common::NewLine;

      // I'll need a sorted, concatenated list of the operators in the fit for filenames
      std::string sOpNameConcat;
      if( bOpSort )
      {
        sOpNameConcat = OpName[0];
        for( std::size_t i = 1; i < OpName.size(); i++ )
        {
          sOpNameConcat.append( 1, '_' );
          sOpNameConcat.append( OpName[i] );
        }
      }
      else
      {
        for( std::size_t i = 0; i < ds.corr.size(); ++i )
        {
          const std::size_t ThisSize{ ds.corr[i].Name_.op.size() };
          for( std::size_t j = ThisSize; j; )
          {
            if( i || j != ThisSize )
              sOpNameConcat.append( Common::Underscore );
            sOpNameConcat.append( OpName[ds.corr[i].Name_.op[--j]] );
          }
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
      std::unique_ptr<Fitter> m{ MakeFitter( cl, ds, ModelArgs, OpName, std::move( cp ) ) };
      const std::string sFitFilename{ Common::MakeFilename( sSummaryBase, "params", Seed, TEXT_EXT ) };
      std::ofstream s;
      for( Common::FitRangesIterator it = fitRanges.begin(); !it.PastEnd(); ++it )
      {
        // Set fit ranges
        int MinExtent = delta;
        std::vector<std::vector<int>> fitTimes( ds.corr.size() );
        for( int i = 0; MinExtent >= delta && i < ds.corr.size(); ++i )
        {
          const Common::FitTime &ft{ it[ModelFitRange[i]] };
          const int Extent{ ft.tf - ft.ti + 1 };
          if( i == 0 || MinExtent > Extent )
            MinExtent = Extent;
          fitTimes[i].resize( Extent );
          for( int j = 0; j < Extent; ++j )
            fitTimes[i][j] = ft.ti + j;
        }
        if( MinExtent >= delta )
        {
          ds.SetFitTimes( fitTimes );
          // Log what file we're processing and when we started
          std::time_t then;
          std::time( &then );
          try
          {
            {
              std::stringstream ss;
              ss << ( doCorr ? "C" : "unc" ) << "orrelated " << m->Type() << " fit on timeslices " << it.to_string( "-", ", " )
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
            auto params = m->PerformFit( doCorr, ChiSq, dof, outBaseFileName, sOpNameConcat, Seed );
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
    " <options> Bootstrap1[[,FitRange1],model1[,params1]] [Bootstrap2[[,FitRange1],model2[,params2]] ...]\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, where:\n"
    " FitRange chooses the n'th fit range (default=0)\n"
    " model    is one of {Exp, Cosh, Sinh, ThreePoint, Const}\n"
    " params   Depend on the model:\n"
    " Exp/Cosh/Sinh\n"
    "   Param1 Snk_Src, sink and source names\n"
    "   Param2 'n' to normalise by energy\n"
    "<options> are:\n"
    "--Hotelling Minimum Hotelling Q-value on central replica (default " << DefaultHotelling << ")\n"
    "--sep    Minimum relative separation between energy levels (default " << DefaultEnergySep << ")\n"
    "--delta  Minimum number of timeslices in fit range (default 3)\n"
    "--retry  Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter   Max iteration count, 0 (default) = unlimited\n"
    "--tol    Tolerance of required fits (default 1e-7)\n"
    "--mindof Minimum degrees of freedom (default 1)\n"
    "--fitter (GSL|Minuit2)[,options] fitter (default GSL,Levenberg-Marquardt)\n"
    "         GSL options: lm, lmaccel, dogleg, ddogleg, subspace2D\n"
    "--covsrc source[,options[,...]] build (co)variance from source, i.e. one of\n"
    "         Binned    Use the already binned data\n"
    "         Raw       Use the raw (unbinned) data\n"
    "         Rebin     Rebin the raw data using bin size(s) specified\n"
    "         Bootstrap Use bootstrap replicas (default)\n"
    "         H5,f[,g],d  Load INVERSE covariance from hdf5 file f, group g, dataset d"
    "--covboot How many bootstrap replicas in covariance (-1=no bootstrap)\n"
    "--guess  List of specific values to use for inital guess"
    "-t     Fit range1[,range2[,...]] (start:stop[:numstart=1[:numstop=numstart]])\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-e     number of Exponents (default 1)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--freeze   Freeze the covariance matrix/variance on the central replica\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--savecmat Save correlation matrix\n"
    "--analytic Analytic derivatives for GSL (default: numeric)\n"
    "--srcsnk   Append _src and _snk to overlap coefficients (ie force different)\n"
    "--opnames  Disable sorting and deduplicating operator name list\n"
    "--help     This message\n";
  }
  return iReturn;
}
