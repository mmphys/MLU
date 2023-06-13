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

// Uncomment this next line to debug without OpenMP
//#define DEBUG_DISABLE_OMP
#ifndef DEBUG_DISABLE_OMP
#include <omp.h>
#endif

// Uncomment the next line to disable Minuit2 (for testing)
//#undef HAVE_MINUIT2

Fitter::Fitter( Model::CreateParams &mcp, DataSet &ds_,
                std::vector<Model::Args> &&ModelArgs_, CovarParams &&cp_, bool bFitCorr_ )
  : bAnalyticDerivatives{ mcp.cl.GotSwitch("analytic") },
    bTestRun{ mcp.cl.GotSwitch("testrun") },
    bCentralGuess{ !mcp.cl.GotSwitch("central") },
    bOverwrite{ mcp.cl.GotSwitch("overwrite") },
    HotellingCutoff{ mcp.cl.SwitchValue<double>( "Hotelling" ) },
    ChiSqDofCutoff{ mcp.cl.SwitchValue<double>( "chisqdof" ) },
    RelEnergySep{ mcp.cl.SwitchValue<double>("sep") },
    MinDof{ mcp.cl.SwitchValue<int>("mindof") },
    Retry{ mcp.cl.SwitchValue<int>("retry") },
    MaxIt{ mcp.cl.SwitchValue<int>("iter") },
    Tolerance{ mcp.cl.SwitchValue<double>("tol") },
    SummaryLevel{ mcp.cl.SwitchValue<int>("summary") },
    bSaveCMat{ mcp.cl.GotSwitch("savecmat") },
    Verbosity{ mcp.cl.SwitchValue<int>("v") },
    UserGuess{ mcp.cl.GotSwitch( "guess" ) },
    Strictness{ mcp.cl.SwitchValue<int>("strict") },
    MonotonicUpperLimit{ mcp.cl.SwitchValue<scalar>("maxE") },
    ErrorDigits{ mcp.cl.SwitchValue<int>("errdig") },
    ds{ ds_ },
    bFitCorr{ bFitCorr_ },
    NumFiles{ static_cast<int>( ModelArgs_.size() ) },
    ModelArgs{ ModelArgs_ }, // Save the original model arguments
    model{ CreateModels( mcp, ModelArgs ) }, // Pass in a copy of model arguments
    NumExponents{ GetNumExponents() },
    mp{ MakeModelParams() },
    bAllParamsKnown{ mp.NumScalars( Param::Type::Variable ) == 0 },
    cp{ std::move( cp_ ) },
    Guess{ mp.NumScalars( Param::Type::All ) },
    OutputModel( ds.NSamples, mp, Common::DefaultModelStats,
                 cp.CovarSampleSize(), cp.bFreeze, cp.Source, cp.RebinSize, cp.CovarNumBoot )
{
  /*assert( NumModelParams == NumFixed + NumVariable && "NumModelParams doesn't match fixed and variable" );
  assert( NumModelParams == NumExponents * NumPerExp + NumOneOff && "NumModelParams doesn't match NumPerExp and NumOneOff" );*/
  if( MinDof < 0 )
    throw std::invalid_argument( "Degrees of freedom (mindof) must be >= 0" );
  if( HotellingCutoff < 0 )
    throw std::invalid_argument( "Hotelling cutoff must be >= 0" );
  if( ChiSqDofCutoff < 0 )
    throw std::invalid_argument( "chi^2 / dof cutoff must be >= 0" );
  if( Retry < 0 )
    throw std::invalid_argument( "Number of fit retries (retry) must be >= 0" );
  if( MaxIt < 0 )
    throw std::invalid_argument( "Maximum fitter iterations (iter) must be >= 0" );
  if( SummaryLevel < ( bAllParamsKnown ? 1 : 0 ) || SummaryLevel > 2 )
    throw std::invalid_argument( "--summary " + std::to_string( SummaryLevel ) + " invalid" );
  // If we've been given a guess, initialise variable parameters with the user-supplied guess
  if( UserGuess )
  {
    const std::vector<scalar> vGuess{Common::ArrayFromString<scalar>(mcp.cl.SwitchValue<std::string>("guess"))};
    if( vGuess.size() != mp.NumScalars( Param::Type::Variable ) )
      throw std::invalid_argument( "Guess contains " + std::to_string( vGuess.size() )
                                 + " parameters, but there are "
                                 + std::to_string( mp.NumScalars( Param::Type::Variable ) )
                                 + " variable parameters" );
    mp.Import<scalar>( Guess, vGuess, Param::Type::Variable, false );
  }
  // Load the constants into the fixed portion of the guess
  ds.GetFixed( Fold::idxCentral, Guess, ParamFixed );
  // Save file attributes
  for( std::size_t i = 0; i < model.size(); ++i )
    OutputModel.FileList.emplace_back( ds(i,bFitCorr).Name_.Filename );
  for( const std::string &f : ds.GetModelFilenames( bFitCorr ? 0 : model.size() ) )
    OutputModel.FileList.emplace_back( f );
  OutputModel.CopyAttributes( ds(0, bFitCorr) );
  OutputModel.NumExponents = NumExponents;
  OutputModel.ModelType = GetModelTypes();
  OutputModel.ModelArgs = GetModelArgs();
}

// Create all the models and return the number of exponents in the fit
std::vector<ModelPtr> Fitter::CreateModels( Model::CreateParams &mcp,
                                            std::vector<Model::Args> ModelArgs )
{
  // Say what we're doing
  {
    std::cout << "Making models";
    const std::string sDescription{ mcp.Description() };
    if( !sDescription.empty() )
      std::cout << " (" << sDescription << ")";
    std::cout << Common::NewLine;
  }
  // Error check
  const int NumModels{ static_cast<int>( ModelArgs.size() ) };
  if( NumModels == 0 )
    throw std::runtime_error( "Fitter::CreateModels() No model arguments" );
  if( bFitCorr && NumModels != ds.corr.size() )
    throw std::runtime_error( "Fitter::CreateModels() " + std::to_string( NumModels )
                             + " models, but " + std::to_string( ds.corr.size() ) + " correlators" );
  else if( !bFitCorr && NumModels > ds.constFile.size() )
    throw std::runtime_error( "Fitter::CreateModels() " + std::to_string( NumModels )
                  + " models, but " + std::to_string( ds.constFile.size() ) + " model files" );
  // Create the models
  std::vector<ModelPtr> model;
  model.reserve( NumModels );
  for( int i = 0; i < ModelArgs.size(); ++i )
  {
    mcp.pCorr = &ds(i, bFitCorr);
    model.emplace_back( Model::MakeModel( mcp, ModelArgs[i] ) );
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
    if( !bFitCorr )
      model[i]->DefineXVector( ds, i );
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
    {
      mp.clear();
      ParamFixed.clear();
    }
    // Ask all the models what parameters they need
    for( ModelPtr &m : model )
    {
      m->param.clear();
      m->AddParameters( mp );
    }
    // Loop through parameters - see if they are available as constants
    for( const Params::value_type &it : mp )
    {
      const Param::Key &pk{ it.first };
      const Param &p{ it.second };
      DataSet::ConstMap::const_iterator cit{ ds.constMap.find( pk ) };
      bool bSwapSourceSink{ false };
      if( cit == ds.constMap.end() && pk.Object.size() == 2 )
      {
        // Try swapping the meson at source and sink
        std::vector<std::string> ObjectReverse( 2 );
        ObjectReverse[0] = pk.Object[1];
        ObjectReverse[1] = pk.Object[0];
        cit = ds.constMap.find( Param::Key( std::move( ObjectReverse ), pk.Name ) );
        bSwapSourceSink = true;
      }
      if( cit != ds.constMap.end() )
      {
        // This is available as a constant
        const ConstantSource &cs{ cit->second };
        const Param &cParam{ ds.GetConstantParam( cs ) };
        if( cParam.size < p.size )
        {
          std::ostringstream os;
          os << "Fitter::MakeModelParams " << pk << "[" << p.size << "] has only "
             << cParam.size << " constants available";
          throw std::runtime_error( os.str().c_str() );
        }
        // Change this parameter to constant
        mp.MakeFixed( pk, bSwapSourceSink );
      }
    }
    // At this point we have a definitive list of parameters - models can save offsets
    mp.AssignOffsets();
    for( ModelPtr &m : model )
      m->SaveParameters( mp );
    // Make list of constants - need to wait until after AssignOffsets()
    for( const Params::value_type &it : mp )
    {
      const Param::Key &pk{ it.first };
      const Param &p{ it.second };
      if( p.type == Param::Type::Fixed )
      {
        DataSet::ConstMap::const_iterator cit{ ds.constMap.find( pk ) };
        assert( cit != ds.constMap.cend() && "Constant unavailable (bug)" );
        ParamFixed.emplace_back( cit->second, p.size, static_cast<int>( p() ) );
      }
    }
    if( mp.NumScalars( Param::Type::Variable ) == 0 )
      bSoluble = true; // Trivial case - no variable parameters
    else
    {
      // Ask the models whether they can guess parameters
      Common::ParamsPairs PP( mp );
      ParamsPairs::StateSize ss{ PP.GetStateSize() };
      ParamsPairs::StateSize ssLast;
      do
      {
        ssLast = ss;
        for( std::size_t i = 0; i < model.size(); ++i )
        {
          model[i]->Guessable( PP );
          if( Verbosity )
            std::cout << "Model " << i << Common::NewLine << PP;
        }
        ss = PP.GetStateSize();
      }
      while( !ss && ss != ssLast );
      bSoluble = ss.Unknown == 0; // Ok if some have ambiguous sign or we only know product
      // If we have unknown product pairs, try ask all the models to minimise number of parameters
      if( bSoluble )
      {
        Common::SignChoice sc( PP, Verbosity );
        mp.SignList = sc;
        std::cout << "Sign groups: " << sc << Common::NewLine;
        mp.DumpSignList( std::cout );
      }
      else if( pass == 0 )
      {
        if( PP.HasUnknownProducts() )
          for( ModelPtr &m : model )
            m->ReduceUnknown( PP );
        else
          pass = 1;
      }
    }
  }
  if( !mp.dispMap.empty() )
    std::cout << mp.dispMap << Common::NewLine;
  return mp;
}

void Fitter::MakeGuess()
{
  // Create list of parameters we know - starting from constants
  const std::size_t NumParams{ mp.NumScalars( Param::Type::All ) };
  if( Guess.size != NumParams )
    throw std::runtime_error( "Fitter::MakeGuess() guess is wrong size" );
  std::vector<bool> ParamKnown( NumParams, false );
  for( const typename Params::value_type &it : mp )
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
    VectorView Data( ds.mFitData.GetCentral() );
    for( std::size_t i = 0; i < model.size(); ++i )
    {
      Data.size( ds.FitTimes[i].size() );
      std::size_t z{ model[i]->Guess( Guess, ParamKnown, mp, Data, ds.FitTimes[i], false ) };
      Data += Data.size();
      if( z )
      {
        NumUnknown += z;
        NumWithUnknowns++;
      }
      if( Verbosity )
      {
        std::cout << " Guess after model " << i << Common::Space << model[i]->Description() << Common::NewLine;
        mp.Dump<scalar>( std::cout, Guess, Param::Type::Variable, nullptr, &ParamKnown );
      }
    }
  }
  while( NumUnknown && NumUnknown != LastUnknown );
  // Since the model is soluble, but we've not guessed all parameters, give the models one last chance
  if( NumUnknown )
  {
    for( std::size_t i = 0; i < model.size(); ++i )
      model[i]->Guess( Guess, ParamKnown, mp, ds.mFitData.GetCentral(), ds.FitTimes[i], true );
    for( bool bKnown : ParamKnown )
      if( !bKnown )
      {
        mp.Dump<scalar>( std::cerr, Guess, Param::Type::Variable, nullptr, &ParamKnown );
        throw std::runtime_error( "Bad model: unable to guess all parameters" );
      }
  }
}

// This should be the only place which knows about different fitters

Fitter * MakeFitterGSL( const std::string &FitterArgs, Model::CreateParams &mcp,
                        DataSet &ds, std::vector<Model::Args> &&ModelArgs,
                        CovarParams &&cp, bool bFitCorr );
#ifdef HAVE_MINUIT2
Fitter * MakeFitterMinuit2( const std::string &FitterArgs, Model::CreateParams &mcp,
                            DataSet &ds, std::vector<Model::Args> &&ModelArgs,
                            CovarParams &&cp, bool bFitCorr );
#endif

Fitter * Fitter::Make( Model::CreateParams &&mcp, DataSet &ds,
                       std::vector<Model::Args> &&ModelArgs, CovarParams &&cp, bool bFitCorr )
{
  Fitter * f;
  std::string FitterArgs{ mcp.cl.SwitchValue<std::string>( "fitter" ) };
  std::string FitterType{ Common::ExtractToSeparator( FitterArgs ) };
  if( Common::EqualIgnoreCase( FitterType, "GSL" ) )
    f = MakeFitterGSL( FitterArgs, mcp, ds, std::move( ModelArgs ), std::move( cp ), bFitCorr );
#ifdef HAVE_MINUIT2
  else if( Common::EqualIgnoreCase( FitterType, "Minuit2" ) )
    f = MakeFitterMinuit2(FitterArgs, mcp, ds, std::move( ModelArgs ), std::move( cp ), bFitCorr);
#endif
  else
    throw std::runtime_error( "Unrecognised fitter: " + FitterType );
  if( !f )
    throw std::runtime_error( "Unrecognised " + FitterType + " fitter options: " + FitterArgs );
  return f;
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
bool Fitter::PerformFit( bool Bcorrelated, double &ChiSq, int &dof_, const std::string &OutBaseName,
                    const std::string &ModelSuffix )
{
  const Common::SeedType Seed{ ds.Seed() };
  bCorrelated = Bcorrelated;
  dof = ds.Extent - static_cast<int>( mp.NumScalars( Param::Type::Variable ) );
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
  if( bTestRun )
  {
    ChiSq = 0;
    return true;
  }

  // See whether this fit already exists
  bool bPerformFit{ true };
  const std::string sModelBase{ OutBaseName + ModelSuffix };
  const std::string ModelFileName{ Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ) };
  if( !bOverwrite && Common::FileExists( ModelFileName ) )
  {
    ModelFile PreBuilt;
    PreBuilt.Read( ModelFileName, "Pre-built: " );
    // TODO: This comparison is incomplete ... but matches previous version
    bool bOK{ dof == PreBuilt.dof };
    bOK = bOK && ds.Extent == PreBuilt.GetExtent();
    bOK = bOK && ds.FitTimes == PreBuilt.FitTimes;
    bOK = bOK && OutputModel.GetColumnNames() == PreBuilt.GetColumnNames();
    if( !bOK )
      throw std::runtime_error( "Pre-built fit not compatible with parameters from this run" );
    bPerformFit = PreBuilt.NewParamsMorePrecise( cp.bFreeze, ds.NSamples );
    if( !bPerformFit )
      ChiSq = dof * PreBuilt.SummaryData( Common::sChiSqPerDof ).Central;
    else
      std::cout << "Overwriting\n";
  }

  bool bOK = true;
  if( bPerformFit )
  {
    // Make somewhere to store the results of the fit for each bootstrap sample
    OutputModel.StdErrorMean.resize( ds.NSamples, ds.Extent ); // Each thread needs its own replica
    OutputModel.ModelPrediction.resize( ds.NSamples, ds.Extent );
    OutputModel.ErrorScaled.resize( ds.NSamples, ds.Extent );
    OutputModel.Name_.Seed = Seed; // TODO: Required?
    OutputModel.FitTimes = ds.FitTimes;
    OutputModel.dof = dof;

    // If there's no user-supplied guess, guess based on the data we're fitting
    if( !UserGuess && !bAllParamsKnown )
      MakeGuess();

    // Fit the central replica, then use this thread as template for all others
    std::unique_ptr<FitterThread> ftC{ MakeThread( bCorrelated, OutputModel ) };
    ftC->Initialise( bSaveCMat ? &sModelBase : nullptr );
    {
      const std::string sDescription{ ftC->Description() };
      if( !sDescription.empty() )
        std::cout << sDescription << Common::NewLine;
    }
    std::cout << "Tolerance " << Tolerance << Common::NewLine;
    if( !bAllParamsKnown )
    {
      std::cout << "Using ";
      if( UserGuess )
        std::cout << "user supplied";
      else if( bCentralGuess )
        std::cout << "central replica fit as";
      else if( bCorrelated )
        std::cout << "uncorrelated fit as";
      else
        std::cout << "simple";
      std::cout << " guess for each replica.\n";
    }
    if( mp.NumScalars( Param::Type::Derived ) )
    {
      std::cout << " Derived parameters: ";
      bool bFirst{ true };
      for( Common::Params::value_type it : mp )
      {
        Param &p{ it.second };
        if( p.type == Param::Type::Derived )
        {
          if( bFirst )
            bFirst = false;
          else
            std::cout << ", ";
          std::cout << it.first;
          if( p.size != 1 )
            std::cout << '[' << p.size << ']';
        }
      }
      std::cout << Common::NewLine;
    }
    if( mp.NumScalars( Param::Type::Fixed ) )
    {
      std::cout << " Fixed parameters:\n";
      mp.Dump( std::cout, Guess, Param::Type::Fixed );
    }
    if( !bAllParamsKnown )
    {
      std::cout << Common::Space << ( UserGuess ? "User supplied" : "Initial" ) << " guess:\n";
      mp.Dump( std::cout, Guess, Param::Type::Variable );
      // For correlated fits, perform an uncorrelated fit on the central replica to use as the guess
      if( bCorrelated && !UserGuess )
        Guess = ftC->UncorrelatedFit();
    }

    // Unless I need the central fit as the guess for each replica, I can do it in parallel below
    if( bCentralGuess )
    {
      // Perform the fit on the central replica. Use that as guess for all other replicas
      ChiSq = ftC->FitOne();
      Guess=ftC->ModelParams;
    }

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
          for( int idx = bCentralGuess ? 0 : Fold::idxCentral; idx < ds.NSamples; ++idx )
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
                fitThread.SwitchReplica( idx );
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
    OutputModel.Guess = Guess;
    OutputModel.FitInput.GetCentral() = ds.mFitData.GetCentral();
    OutputModel.FitInput.Replica = ds.mFitData.Replica;
    OutputModel.SetSummaryNames( Common::sParams );
    OutputModel.MakeCorrSummary();
    OutputModel.CheckParameters( Strictness, MonotonicUpperLimit );
    {
      // Show parameters - in neat columns
      const Common::ValWithEr<scalar> * const v{ &OutputModel.SummaryData() };
      const std::vector<std::string> &Cols{ OutputModel.GetColumnNames() };
      const int maxLen{ static_cast<int>( std::max_element( Cols.begin(), Cols.end(),
                                          [](const std::string &a, const std::string &b)
                                          { return a.length() < b.length(); } )->length() ) };
      const int NumWidth{ static_cast<int>( std::cout.precision() ) + 7 };
      std::string Preamble;
      {
        std::ostringstream os;
        os << OutputModel.GetSummaryNames()[0] << " after " << OutputModel.NumSamples()
           << Common::Space << OutputModel.getData().SeedTypeString() << " replicas\n";
        Preamble = os.str();
      }
      std::cout << Preamble;
      for( int i = 0; i < OutputModel.Nt(); ++i )
      {
        const char cOK{ v[i].Check == 0 ? 'x' : ' ' };
        std::cout << cOK << std::setw(maxLen) << Cols[i] << std::left
                  << Common::Space << std::setw(NumWidth) << v[i].Central
                  << " +"    << std::setw(NumWidth) << (v[i].High - v[i].Central)
                  << " -"    << std::setw(NumWidth) << (v[i].Central - v[i].Low)
                  << " Max " << std::setw(NumWidth) << v[i].Max
                  << " Min "                        << v[i].Min
                  << Common::NewLine << std::right; // Right-aligned is the default
      }
      // Show parameters again - but this time in simplified format
      std::cout << Preamble;
      for( int i = 0; i < OutputModel.Nt(); ++i )
      {
        const char cOK{ v[i].Check == 0 ? 'x' : ' ' };
        if( v[i].Check == 0 )
          bOK = false;
        std::cout << cOK << std::setw(maxLen) << Cols[i] << std::left
                  << Common::Space << v[i].to_string( ErrorDigits )
                  << Common::NewLine << std::right; // Right-aligned is the default
      }
    }
    // Save the file
    if( !bAllParamsKnown )
      OutputModel.Write( ModelFileName );
    if( SummaryLevel >= 1 )
    {
      const bool bVerbose{ true };
      if( bFitCorr )
        OutputModel.WriteSummaryTD( ds,
                    Common::MakeFilename( sModelBase, Common::sModel + "_td", Seed, TEXT_EXT ),
                    bVerbose );
      else
        WriteSummaryTD( Common::ReplaceExtension( ModelFileName, DAT_EXT ), bVerbose );
    }
    if( SummaryLevel >= 2 )
      OutputModel.WriteSummary( Common::ReplaceExtension( ModelFileName, TEXT_EXT ) );
  }
  return bOK;
}

std::vector<std::string> Fitter::GetModelTypes() const
{
  std::vector<std::string> vs;
  for( const ModelPtr &p : model )
  {
    std::ostringstream os;
    os << p->Type();
    vs.push_back( os.str() );
  }
  return vs;
}

std::vector<std::string> Fitter::GetModelArgs() const
{
  std::vector<std::string> vs;
  for( const Model::Args &a : ModelArgs )
    vs.push_back( a.ToString() );
  while( vs.size() && vs.back().empty() )
    vs.resize( vs.size() - 1 );
  return vs;
}

void Fitter::WriteSummaryTD( const std::string &sOutFileName, bool bVerboseSummary )
{
  const ValWithEr veZero( 0, 0, 0, 0, 0, 0 );
  //using namespace Common::CorrSumm;
  assert( std::isnan( Common::NaN ) && "Compiler does not support quiet NaNs" );
  std::ofstream ofs( sOutFileName );
  Common::SummaryHeader<scalar>( ofs, sOutFileName );
  OutputModel.SummaryComments( ofs, bVerboseSummary );
  // Write the list of ensembles
  {
    std::string *pLast{ nullptr };
    ofs << "# Ensembles:";
    for( std::size_t i = 0; i < model.size(); ++i )
      if( i == 0 || Common::CompareIgnoreCase( ds.constFile[i]->Ensemble, *pLast ) )
      {
        pLast = &ds.constFile[i]->Ensemble;
        ofs << Common::Space << *pLast;
      }
    ofs << "\n# Primary key: ensemble " << model[0]->XVectorKeyName() << Common::NewLine;
  }
  // Write column names
  static constexpr int idxData{ 0 };
  //static constexpr int idxTheory{ 1 };
  ofs << "# Global fit results\n";
  ofs << "model ensemble " << model[0]->XVectorKeyName() << Common::Space;
  ValWithEr::Header( "theory", ofs );
  ofs << Common::Space;
  ValWithEr::Header( "data", ofs );
  for( const Param::Key &k : model[0]->XVectorKeyNames() )
  {
    const Param &p{ mp.Find( k, "WriteSummaryTD()" )->second };
    const std::size_t idx{ p() };
    for( std::size_t i = 0; i < p.size; ++i )
    {
      ofs << Common::Space;
      ValWithEr::Header( k.ShortName( i, p.size ), ofs );
    }
  }
  ofs << Common::NewLine;
  // Grab the model data and theory with bootstrapped errors
  const int Extent{ OutputModel.GetExtent() };
  if( !Extent )
    return;
  std::array<VectorView, 2> vv; // View into data and theory replica data. Used to populate Value
  std::vector<scalar> Buffer; // Scratch buffer for ValWithEr<T>
  std::array<std::vector<ValWithEr>,2> Value; // Data and theory values with errors
  for( std::size_t i = 0; i < Value.size(); ++i )
  {
    JackBoot &ThD{ i == idxData ? OutputModel.FitInput : OutputModel.ModelPrediction };
    ThD.MakeStatistics( Value[i] );
  }
  // Write theory and data CORRELATOR values, with each data point in fit on a separate line
  for( std::size_t i = 0; i < OutputModel.FitTimes.size(); ++i )
  {
    if( OutputModel.FitTimes[i].size() != 1 )
      throw std::runtime_error( "Fitter::WriteSummaryTD fit times "
                               + std::to_string( OutputModel.FitTimes[i].size() ) + " != 1" );
    ModelFile &mf{ *ds.constFile[i] };
    Model &m{ *model[i] };
    ofs << i << Common::Space << mf.Ensemble << Common::Space << m.XVectorKey();
    for( std::size_t j = Value.size(); j-- != 0; )
      ofs << Common::Space << Value[j][i];
    for( const Param::Key &k : m.XVectorKeyNames() )
    {
      const Param &p{ mp.Find( k, "WriteSummaryTD()" )->second };
      const std::size_t idx{ p() };
      for( std::size_t j = 0; j < p.size; ++j )
        ofs << Common::Space << OutputModel.SummaryData( 0, idx + j );
    }
    ofs << Common::NewLine;
  }
}

