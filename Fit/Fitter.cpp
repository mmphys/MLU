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

#include "Param.hpp"
#include "Fitter.hpp"
#include "FitterThread.hpp"

// Uncomment this next line to debug without OpenMP
//#define DEBUG_DISABLE_OMP

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
