/*************************************************************************************
 
 Base class for fitter
 
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

#include "Fitter.hpp"
#include "FitterThread.hpp"
#include "Param.hpp"

// Uncomment this next line to debug without OpenMP
//#define DEBUG_DISABLE_OMP

Fitter::Fitter( const Common::CommandLine &cl, const DataSet &ds_,
                std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames_,
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
    UserGuess{ cl.GotSwitch( "guess" ) },
    ds{ std::move( ds_ ) },
    NumFiles{ static_cast<int>( ds.corr.size() ) },
    OpNames{ opNames_ },
    //NumOps{ static_cast<int>( OpNames.size() ) },
    model{ CreateModels( cl, ModelArgs ) },
    NumExponents{ GetNumExponents() },
    mp{ MakeModelParams() },
    /*ParamNames( MakeParamNames() ),
    ParamFixed{ MakeParamFixed() },
    ParamVariable{ MakeParamVariable() },
    NumModelParams{ static_cast<int>( ParamNames.size() ) },
    NumPerExp{ static_cast<int>( PerExpNames.size() ) },
    NumOneOff{ NumModelParams - NumExponents * NumPerExp },
    NumFixed{ static_cast<int>( ParamFixed.size() ) },
    NumVariable{ static_cast<int>( ParamVariable.size() ) },*/
    cp{ std::move( cp_ ) },
    Guess{ mp.NumScalars( Param::Type::All ) }
{
  assert( ds.corr.size() == NumFiles && "Number of files extremely large" );
  /*assert( NumModelParams == NumFixed + NumVariable && "NumModelParams doesn't match fixed and variable" );
  assert( NumModelParams == NumExponents * NumPerExp + NumOneOff && "NumModelParams doesn't match NumPerExp and NumOneOff" );*/
  if( MinDof < 0 )
    throw std::invalid_argument( "Degrees of freedom (mindof) must be >= 0" );
  if( Retry < 0 )
    throw std::invalid_argument( "Number of fit retries (retry) must be >= 0" );
  if( MaxIt < 0 )
    throw std::invalid_argument( "Maximum fitter iterations (iter) must be >= 0" );
  // If we've been given a guess, initialise variable parameters with the user-supplied guess
  if( UserGuess )
  {
    const std::vector<scalar> vGuess{Common::ArrayFromString<scalar>(cl.SwitchValue<std::string>("guess"))};
    if( vGuess.size() != mp.NumScalars( Param::Type::Variable ) )
      throw std::invalid_argument( "Guess contains " + std::to_string( vGuess.size() )
                                 + " parameters, but there are "
                                 + std::to_string( mp.NumScalars( Param::Type::Variable ) )
                                 + " variable parameters" );
    mp.Import( Guess, vGuess );
  }
  // Load the constants into the fixed portion of the guess
  ds.GetFixed( Fold::idxCentral, Guess, ParamFixed );
}

// Create all the models and return the number of exponents in the fit
std::vector<ModelPtr> Fitter::CreateModels( const Common::CommandLine &cl,
                                            std::vector<Model::Args> &ModelArgs )
{
  Model::CreateParams cp( OpNames, cl );
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
    cp.pCorr = &ds.corr[i];
    model.emplace_back( Model::MakeModel( cp, ModelArgs[i] ) );
    if( !ModelArgs[i].empty() )
    {
      std::ostringstream os;
      os << "Model " << i << " has " << ModelArgs[i].size() << " invalid parameters: ";
      bool bFirst{ true };
      for( typename Model::Args::value_type v : ModelArgs[i] )
      {
        if( bFirst )
          bFirst = false;
        else
          os << Common::CommaSpace;
        os << v.first;
        if( !v.second.empty() )
          os << Common::EqualSign << v.second;
      }
      throw std::runtime_error( os.str().c_str() );
    }
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

// Finalise the parameter lists the models will fit
// Build complete list of fixed and variable parameters
// Tell each model where to get the parameters they are interested in
Params Fitter::MakeModelParams()
{
  Params mp;
  bool bSoluble{ false };
  for( int pass = 0; pass < 2 && !bSoluble; ++pass )
  {
    if( pass )
      for( ModelPtr &m : model )
        m->ReduceUnknown();
    // Ask all the models what parameters they need
    mp.clear();
    for( ModelPtr &m : model )
      m->AddParameters( mp );
    mp.AssignOffsets();
    for( ModelPtr &m : model )
      m->SaveParameters( mp );
    // Loop through parameters - see if they are available as constants
    DataSet::ConstantSource::Key key;
    ParamFixed.clear();
    for( typename Params::value_type it : mp )
    {
      const Param::Key &pk{ it.first };
      const Param &p{ it.second };
      key.Object = pk.Object;
      key.Name = pk.Name;
      typename DataSet::ConstMap::const_iterator cit{ ds.constMap.find( key ) };
      if( cit != ds.constMap.end() )
      {
        // This is available as a constant
        const DataSet::ConstantSource &cs{ cit->second };
        if( cs.idx.size() < p.size )
        {
          std::ostringstream os;
          os << "AreModelsSoluble " << pk.Object << "/" << pk.Name << "[" << p.size << "] has only "
          << cs.idx.size() << " constants available";
          throw std::runtime_error( os.str().c_str() );
        }
        mp.Add( pk, p.size, false, Param::Type::Fixed );
        // Now remember where this comes from
        ParamFixed.emplace_back( cs, p.size, static_cast<int>( p() ) );
      }
    }
    // Create list of parameters we know - starting from constants
    const std::size_t NumParams{ mp.NumScalars( Param::Type::All ) };
    std::vector<bool> ParamKnown( NumParams, false );
    for( typename Params::value_type it : mp )
    {
      const Param &p{ it.second };
      if( p.type == Param::Type::Fixed )
      {
        for( std::size_t i = 0; i < p.size; ++i )
          ParamKnown[ p( i ) ] = true;
      }
    }
    // Now ask the models whether they can guess parameters
    std::size_t LastUnknown;
    std::size_t NumUnknown{ 0 };
    std::size_t NumWithUnknowns;
    do
    {
      LastUnknown = NumUnknown;
      NumUnknown = 0;
      NumWithUnknowns = 0;
      for( std::size_t i = 0; i < model.size(); ++i )
      {
        std::size_t z{ model[i]->Guessable( ParamKnown, false ) };
        if( z )
        {
          NumUnknown += z;
          NumWithUnknowns++;
        }
      }
    }
    while( NumUnknown && NumUnknown != LastUnknown );
    bSoluble = ( NumUnknown <= NumWithUnknowns );
    // If the model is soluble, but we've not guessed all parameters, force the models to do so
    if( bSoluble && NumUnknown )
    {
      for( std::size_t i = 0; i < model.size(); ++i )
        model[i]->Guessable( ParamKnown, true );
      for( bool bKnown : ParamKnown )
        if( !bKnown )
          throw std::runtime_error( "Model is not solvable" );
    }
  }
  if( mp.NumScalars( Param::Type::Variable ) == 0 )
    throw std::runtime_error( "Model has no parameters to fit" );
  return mp;
}

void Fitter::MakeGuess()
{
  // Create list of parameters we know - starting from constants
  const std::size_t NumParams{ mp.NumScalars( Param::Type::All ) };
  if( Guess.size != NumParams )
    throw std::runtime_error( "Fitter::MakeGuess() guess is wrong size" );
  std::vector<bool> ParamKnown( NumParams, false );
  for( typename Params::value_type it : mp )
  {
    const Param &p{ it.second };
    if( p.type == Param::Type::Fixed )
    {
      for( std::size_t i = 0; i < p.size; ++i )
        ParamKnown[ p( i ) ] = true;
    }
  }
  ds.GetFixed( Fold::idxCentral, Guess, ParamFixed );
  // Now ask the models to guess parameters
  std::size_t LastUnknown;
  std::size_t NumUnknown{ 0 };
  std::size_t NumWithUnknowns;
  do
  {
    LastUnknown = NumUnknown;
    NumUnknown = 0;
    NumWithUnknowns = 0;
    for( std::size_t i = 0; i < model.size(); ++i )
    {
      std::size_t z{ model[i]->Guess( Guess, ParamKnown, ds.vCentral, ds.FitTimes[i], false ) };
      if( z )
      {
        NumUnknown += z;
        NumWithUnknowns++;
      }
    }
  }
  while( NumUnknown && NumUnknown != LastUnknown );
  // Since the model is soluble, but we've not guessed all parameters, give the models one last chance
  if( NumUnknown )
  {
    for( std::size_t i = 0; i < model.size(); ++i )
      model[i]->Guess( Guess, ParamKnown, ds.vCentral, ds.FitTimes[i], true );
    for( bool bKnown : ParamKnown )
      if( !bKnown )
        throw std::runtime_error( "Bad model: unable to guess all parameters" );
  }
}

void Fitter::SayNumThreads( std::ostream &os )
{
#ifdef DEBUG_DISABLE_OMP
  os << " with Open MP disabled";
#else
  os << " using " << omp_get_max_threads() << " Open MP threads";
#endif
}

bool Fitter::Dump( int idx ) const
{
  return Verbosity > 2 || ( Verbosity >= 1 && idx == Fold::idxCentral );
}

void Fitter::Dump( int idx, const std::string &Name, const Matrix &m ) const
{
  if( Verbosity > 2 || ( Verbosity >= 1 && idx == Fold::idxCentral ) )
    std::cout << Name << Common::Space << m << Common::NewLine;
}

void Fitter::Dump( int idx, const std::string &Name, const Vector &v ) const
{
  if( Verbosity > 2 || ( Verbosity >= 1 && idx == Fold::idxCentral ) )
    std::cout << Name << Common::Space << v << Common::NewLine;
}

void Fitter::SaveMatrixFile( const Matrix &m, const std::string &Type, const std::string &Filename,
                               const char *pGnuplotExtra ) const
{
  if( model.size() != ds.corr.size() )
    throw std::runtime_error( "ModelSet doesn't match DataSet" );
  std::vector<std::string> Abbreviations;
  Abbreviations.reserve( ds.corr.size() );
  std::vector<std::string> FileComments;
  FileComments.reserve( ds.corr.size() );
  for( std::size_t f = 0; f < model.size(); ++f )
  {
    Abbreviations.emplace_back( model[f]->Description() );
    // Save type of each model and parameters
    std::ostringstream s;
    s << "# Model" << f << ": " << model[f]->Type();
    for( std::size_t i = 0; i < model[f]->param.size(); ++i )
      s << Common::CommaSpace << model[f]->param[i];
    s << Common::NewLine;
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
  dof = ds.Extent - mp.NumScalars( Param::Type::Variable );
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
  // TODO: This is temporary - makes it look like existing format
  // TODO: Update to save the parameter list + additional columns
  std::vector<std::string> OldFormatOpNames;
  for( typename Params::value_type it : mp )
  {
    if( it.second.size > 1 && !Common::EqualIgnoreCase( E, it.first.Name ) )
    {
      std::ostringstream os;
      os << it.first;
      OldFormatOpNames.push_back( os.str() );
    }
  }
  const std::size_t NumParams{ mp.NumScalars( Param::Type::All ) };
  ModelFile OutputModel( OldFormatOpNames, NumExponents, NumFiles, tMin, tMax, dof,
                         true, //TODO: !bForceSrcSnkDifferent
                         cp.bFreeze, ds.NSamples, NumParams + 1 );
  for( const Fold &f : ds.corr )
    OutputModel.FileList.emplace_back( f.Name_.Filename );
  OutputModel.CopyAttributes( ds.corr[0] );
  {
    std::vector<std::string> ColNames{ NumParams + 1 };
    //TODO: Work out how to copy parameter names into output file properly
    if( ColNames.size() == 1 )
      OutputModel.OpNames.clear();
    ColNames.push_back( Common::sChiSqPerDof );
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
      ChiSq = dof * PreBuilt.getSummaryData()[NumParams].Central; // Last summary has chi squared per dof
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

    // If there's no user-supplied guess, guess based on the data we're fitting
    if( !UserGuess )
      MakeGuess();

    // Fit the central replica, then use this thread as template for all others
    std::unique_ptr<FitterThread> ftC{ MakeThread( bCorrelated, OutputModel, CorrSynthetic ) };
    ftC->SetReplica( Fold::idxCentral, true, true, bSaveCMat ? &sModelBase : nullptr );
    {
      const std::string sDescription{ ftC->Description() };
      if( !sDescription.empty() )
        std::cout << sDescription << Common::NewLine;
    }
    std::cout << "Tolerance " << Tolerance << ". Using uncorrelated fit as guess for each replica.\n";
    if( mp.NumScalars( Param::Type::Fixed ) )
    {
      std::cout << " Fixed parameters:\n";
      mp.Dump( std::cout, Guess, Param::Type::Fixed );
    }
    std::cout << " Initial guess:\n";
    mp.Dump( std::cout, Guess, Param::Type::Variable );

    // For correlated fits, perform an uncorrelated fit on the central replica to use as the guess
    if( bCorrelated && !UserGuess )
      Guess = ftC->UncorrelatedFit();

    // Perform all the fits in parallel on separate OpenMP threads
    // As long as exceptions are caught on the same thread that threw them all is well
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
          // By cloning the central replica thread, we can avoid rebuilding central replica correlation matrix
          std::unique_ptr<FitterThread> ft{ ftC->Clone() };
          FitterThread &fitThread( *ft.get() );
#ifndef DEBUG_DISABLE_OMP
          #pragma omp for schedule(dynamic)
#endif
          for( int idx = Fold::idxCentral; idx < ds.NSamples; ++idx )
          {
            try
            {
              bool WasAbort;
#ifndef DEBUG_DISABLE_OMP
              #pragma omp atomic read
#endif
                WasAbort = Abort;
              if( !WasAbort )
              {
                fitThread.SetReplica( idx );
                scalar z{ fitThread.FitOne() };
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
  std::vector<Common::ValWithEr<scalar>> Results( mp.NumScalars( Param::Type::All ) );
  std::vector<scalar> data( NumSummaries );
  for( std::size_t p = 0; p < Results.size(); p++ ) {
    double Central = OutputModel[Fold::idxCentral][p];
    std::size_t Count{ 0 };
    for( int j = 0; j < NumSummaries; ++j )
    {
      double d = OutputModel[j][p];
      if( std::isfinite( d ) )
        data[Count++] = d;
    }
    Results[p].Get( Central, data, Count );
  }
  return Results;
}
