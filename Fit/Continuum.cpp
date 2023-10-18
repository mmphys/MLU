/*************************************************************************************
 
 Chiral continuum fit
 
 Source file: Continuum.cpp

 Copyright (C) 2023
 
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

#include "Continuum.hpp"
#include "ModelContinuum.hpp"
#include <MLU/DebugInfo.hpp>
#include <chrono>

// Indices for operators in correlator names
const char * pSrcSnk[2] = { "src", "snk" };

void EnsembleMapT::Loaded()
{
  // Assign an index for each Ensemble
  unsigned int idx{};
  for( value_type &it : *this )
    it.second.idx = idx++;
}

CreateParams::CreateParams( const std::vector<std::string> &O, const Common::CommandLine &C,
                            const ContinuumFit &parent )
: Model::CreateParams( O, C ), Parent{ parent }
{
}

const std::string ContinuumFit::sPDG{ "PDG" };
const std::string &ContinuumFit::FieldQSq{ ModelContinuum::FieldQSq };
const std::string &ContinuumFit::FieldEL{ ModelContinuum::FieldEL };
constexpr std::size_t ContinuumFit::idxUnused;

ContinuumFit::ContinuumFit( Common::CommandLine &cl_ )
: cl{cl_},
  NumD{ cl.SwitchValue<unsigned int>("d") },
  NumE{ cl.SwitchValue<unsigned int>("e") },
  NumSamples{cl.SwitchValue<int>("n")},
  doCorr{ !cl.GotSwitch( "uncorr" ) },
  CovarBlock{ cl.GotSwitch( "block" ) },
  ffDefault{Common::FromString<Common::FormFactor>(Common::TrimAt(cl.SwitchValue<std::string>("f")))},
  inBase{ cl.SwitchValue<std::string>("i") },
  outBaseFileName{ cl.SwitchValue<std::string>("o") },
  ds{NumSamples},
  uiFF{0}
{
  GetEnabled( cl.SwitchValue<std::string>("f") );
  Common::MakeAncestorDirs( outBaseFileName );
  PoleMass[idxFF0] = cl.GotSwitch( "poles" ) ? cl.SwitchValue<scalar>("poles") : 0;
  PoleMass[idxFFPlus] = cl.GotSwitch( "polev" ) ? cl.SwitchValue<scalar>("polev") : 0;
}

void ContinuumFit::ParamsAdjust( Common::Params &mp, const Fitter &f )
{
  // Global parameters
  kfPi.Name = "fPi";
  mp.Add( kfPi );
  kmPDGPi.Name = "PDGPi";
  mp.Add( kmPDGPi );
  kPDGH.Name = sPDG + Meson[idxSrc];
  mp.Add( kPDGH );
  kPDGL.Name = sPDG + Meson[idxSnk];
  mp.Add( kPDGL );

  // Per form factor parameters
  for( int idxFF = 0; idxFF < NumFF; ++idxFF )
  {
    if( uiFF & ffMaskFromIndex( idxFF ) )
    {
      const Common::FormFactor ff{ ffIndexReverse( idxFF ) };
      const std::string &sFF{ GetFormFactorString( ff ) };
      if( PoleMass[idxFF] == 0 )
      {
        kPDGDStar[idxFF].Name = GetPoleMassName( ff, f.ds.constFile[0]->Name_ );
        mp.Add( kPDGDStar[idxFF], 1 );
      }
      kDelta[idxFF].Object.push_back( sFF );
      kDelta[idxFF].Name = "Delta";
      mp.Add( kDelta[idxFF], 1, false, Param::Type::Derived );
      // Add constants
      if( c0Enabled[idxFF] )
      {
        kC0[idxFF].Object.push_back( sFF );
        kC0[idxFF].Name = "c0";
        mp.Add( kC0[idxFF] );
      }
      if( c1Enabled[idxFF] )
      {
        kC1[idxFF].Object.push_back( sFF );
        kC1[idxFF].Name = "c1";
        mp.Add( kC1[idxFF] );
      }
      for( int i = 0; i < NumD; ++i )
        if( dEnabled[idxFF][i] )
        {
          kD[idxFF][i].Object.push_back( sFF );
          kD[idxFF][i].Name = "d" + std::to_string( i );
          mp.Add( kD[idxFF][i] );
        }
      for( int i = 0; i < NumE; ++i )
        if( eEnabled[idxFF][i] )
        {
          kE[idxFF][i].Object.push_back( sFF );
          kE[idxFF][i].Name = "e" + std::to_string( i );
          mp.Add( kE[idxFF][i] );
        }
    }
  }
  // If I'm fitting both form factors, I need to impose a constraint
  if( bDoConstraint && ( uiFF & uiFF0 ) && ( uiFF & uiFFPlus ) )
  {
    // Impose the constraint on a fit constant that is enabled
    if( c0Enabled[idxFFPlus] )
      WhichConstraint = -1;
    else
    {
      WhichConstraint = 0;
      while( !eEnabled[idxFFPlus][WhichConstraint] )
        ++WhichConstraint;
      if( WhichConstraint >= NumE )
        throw std::runtime_error( "Can't impose constraint when c0 and all energies disabled on f+" );
    }
    // The constraint becomes a derived parameter
    mp.SetType( ConstraintKey(), Param::Type::Derived );
  }

  // Per ensemble parameters
  idxaInv.resize( EnsembleMap.size() );
  idxmPi.resize( EnsembleMap.size() );
  kaInv.resize( EnsembleMap.size() );
  kmPi.resize( EnsembleMap.size() );
  if( NeedFV() )
  {
    idxFVSim.resize( EnsembleMap.size() );
    idxFVPhys.resize( EnsembleMap.size() );
    kFVSim.resize( EnsembleMap.size() );
    kFVPhys.resize( EnsembleMap.size() );
  }
  if( NeedChiral() )
  {
    idxChiSim.resize( EnsembleMap.size() );
    kChiSim.resize( EnsembleMap.size() );
    kChiPhys.Name = "ChiPhys";
    mp.Add( kChiPhys, 1, false, Param::Type::Derived );
  }
  if( NeedFV() || NeedChiral() )
  {
    idxChiFV.resize( EnsembleMap.size() );
    kChiFV.resize( EnsembleMap.size() );
  }
  if( c1Enabled[idxFF0] || c1Enabled[idxFFPlus] )
  {
    idxDeltaMPi.resize( EnsembleMap.size() );
    kDeltaMPi.resize( EnsembleMap.size() );
  }
  for( const typename EnsembleMapT::value_type &it : EnsembleMap )
  {
    const std::string &sEnsemble{ it.first };
    const std::size_t i{ it.second.idx };
    kaInv[i].Object.push_back( sEnsemble );
    kaInv[i].Name = "aInv";
    mp.Add( kaInv[i] );
    kmPi[i].Object.push_back( sEnsemble );
    kmPi[i].Name = "mPi";
    mp.Add( kmPi[i] );
    if( NeedFV() )
    {
      kFVSim[i].Object.push_back( sEnsemble );
      kFVSim[i].Name = "FVSim";
      mp.Add( kFVSim[i], 1, false, Param::Type::Derived );
      kFVPhys[i].Object.push_back( sEnsemble );
      kFVPhys[i].Name = "FVPhys";
      mp.Add( kFVPhys[i], 1, false, Param::Type::Derived );
    }
    if( NeedChiral() )
    {
      kChiSim[i].Object.push_back( sEnsemble );
      kChiSim[i].Name = "ChiSim";
      mp.Add( kChiSim[i], 1, false, Param::Type::Derived );
    }
    if( NeedFV() || NeedChiral() )
    {
      kChiFV[i].Object.push_back( sEnsemble );
      kChiFV[i].Name = "ChiFV";
      mp.Add( kChiFV[i], 1, false, Param::Type::Derived );
    }
    if( c1Enabled[idxFF0] || c1Enabled[idxFFPlus] )
    {
      kDeltaMPi[i].Object.push_back( sEnsemble );
      kDeltaMPi[i].Name = "DeltaMPi";
      mp.Add( kDeltaMPi[i], 1, false, Param::Type::Derived );
    }
  }
}

void ContinuumFit::SaveParameters( Common::Params &mp, const Fitter &f )
{
  // PDG masses for heavy and light, plus the Delta in this case
  static const std::string sFindError{ "ContinuumFit::SaveParameters() finding model constants" };
  idxfPi = mp.at( kfPi )();
  idxmPDGPi = mp.at( kmPDGPi )();
  idxPDGH = mp.at( kPDGH )();
  idxPDGL = mp.at( kPDGL )();
  if( NeedChiral() )
    idxChiralPhys = mp.at( kChiPhys )();

  // Per form factor parameter offsets
  for( int idxFF = 0; idxFF < NumFF; ++idxFF )
  {
    if( uiFF & ffMaskFromIndex( idxFF ) )
    {
      if( PoleMass[idxFF] == 0 )
        idxPDGDStar[idxFF] = mp.at( kPDGDStar[idxFF] )();
      idxDelta[idxFF] = mp.at( kDelta[idxFF] )();
      if( c0Enabled[idxFF] )
        idxC0[idxFF] = mp.at( kC0[idxFF] )();
      if( c1Enabled[idxFF] )
        idxC1[idxFF] = mp.at( kC1[idxFF] )();
      for( int i = 0; i < NumD; ++i )
        if( dEnabled[idxFF][i] )
          idxD[idxFF][i] = mp.at( kD[idxFF][i] )();
      for( int i = 0; i < NumE; ++i )
        if( eEnabled[idxFF][i] )
          idxE[idxFF][i] = mp.at( kE[idxFF][i] )();
    }
  }
  
  // Per ensemble parameter offsets
  for( std::size_t i = 0; i < EnsembleMap.size(); ++i )
  {
    idxaInv[i] = mp.at( kaInv[i] )();
    idxmPi[i] = mp.at( kmPi[i] )();
    if( NeedFV() )
    {
      idxFVSim[i] = mp.at( kFVSim[i] )();
      idxFVPhys[i] = mp.at( kFVPhys[i] )();
    }
    if( NeedChiral() )
      idxChiSim[i] = mp.at( kChiSim[i] )();
    if( NeedFV() || NeedChiral() )
      idxChiFV[i] = mp.at( kChiFV[i] )();
    if( c1Enabled[idxFF0] || c1Enabled[idxFFPlus] )
      idxDeltaMPi[i] = mp.at( kDeltaMPi[i] )();
  }

  // If I'm imposing a constraint, remember what on!
  if( bDoConstraint && ( uiFF & uiFF0 ) && ( uiFF & uiFFPlus ) )
  {
    idxConstraint = mp.at( ConstraintKey() )();
  }
}

void ContinuumFit::SetReplica( Vector &ModelParams ) const
{
  const scalar mPDGL{ ModelParams[idxPDGL] };
  const scalar mPDGH{ ModelParams[idxPDGH] };
  // Compute delta for each form factor
  for( int idxFF = 0; idxFF < NumFF; ++idxFF )
  {
    if( uiFF & ffMaskFromIndex( idxFF ) )
    {
      if( PoleMass[idxFF] == 0 )
      {
        const scalar PoleMass{ ModelParams[idxPDGDStar[idxFF]] };
        ModelParams[idxDelta[idxFF]] = 0.5 * ( ( PoleMass*PoleMass - mPDGL*mPDGL ) / mPDGH - mPDGH );
      }
      else
        ModelParams[idxDelta[idxFF]] = PoleMass[idxFF];
    }
  }
  // Compute chiral and finite volume corrections for each ensemble
  if( NeedFV() || NeedChiral() )
  {
    const scalar FourPifPi{ FourPi * ModelParams[idxfPi] };
    const scalar Denom{ 1. / ( FourPifPi * FourPifPi ) };
    scalar sChiPhys = 0;
    if( NeedChiral() )
    {
      sChiPhys = Common::DeltaF::ChiralLog( ModelParams[idxmPDGPi], LambdaInv );
      ModelParams[idxChiralPhys] = sChiPhys;
    }
    for( const typename EnsembleMapT::value_type &ei : EnsembleMap )
    {
      const std::size_t i{ ei.second.idx };
      const unsigned int aInv_L{ ei.second.aInv_L };
      const scalar L{ aInv_L / ModelParams[idxaInv[i]] };
      scalar sFVSim = 0, sFVPhys = 0, sChiSim = 0;
      if( NeedFV() )
      {
        sFVSim = Common::DeltaF::FiniteVol( ModelParams[idxmPi[i]], L, aInv_L );
        sFVPhys = Common::DeltaF::FiniteVol( ModelParams[idxmPDGPi], L, aInv_L );
        ModelParams[idxFVSim[i]] = sFVSim;
        ModelParams[idxFVPhys[i]] = sFVPhys;
      }
      if( NeedChiral() )
      {
        sChiSim = Common::DeltaF::ChiralLog( ModelParams[idxmPi[i]], LambdaInv );
        ModelParams[idxChiSim[i]] = sChiSim;
      }
      ModelParams[idxChiFV[i]] = ( DeltaF(sChiSim,sFVSim ) - DeltaF(sChiPhys,sFVPhys) ) * Denom;
    }
  }
  if( c1Enabled[idxFF0] || c1Enabled[idxFFPlus] )
  {
    const scalar mPDGPi{ ModelParams[idxmPDGPi] };
    const scalar mPDGPiSq{ mPDGPi * mPDGPi };
    for( const typename EnsembleMapT::value_type &ei : EnsembleMap )
    {
      const std::size_t i{ ei.second.idx };
      const scalar mPi{ ModelParams[idxmPi[i]] };
      ModelParams[idxDeltaMPi[i]] = ( mPi * mPi - mPDGPiSq ) * LambdaInvSq;
    }
  }
}

void ContinuumFit::ComputeDerived( Vector &ModelParams ) const
{
  if( bDoConstraint && ( uiFF & uiFF0 ) && ( uiFF & uiFFPlus ) )
  {
    // I'm fitting both form factors - impose a constraint on fplus-c0
    const scalar mPDGL{ ModelParams[idxPDGL] };
    const scalar mPDGH{ ModelParams[idxPDGH] };
    const scalar E0{ EOfQSq( mPDGH, mPDGL, 0 ) };
    const scalar EOnLambda{ E0 * LambdaInv };
    scalar ESum[NumFF];
    ESum[idxFF0] = c0Enabled[idxFF0] ? ModelParams[idxC0[idxFF0]] : 0;
    ESum[idxFFPlus] = 0;
    scalar ETerm = EOnLambda;
    for( unsigned int i = 0; i < NumE; ++i )
    {
      for( int idxFF = 0; idxFF < NumFF; ++idxFF )
        if( eEnabled[idxFF][i] && ( idxFF==idxFF0 || WhichConstraint>NumE || WhichConstraint<i ) )
          ESum[idxFF] += ModelParams[idxE[idxFF][i]] * ETerm;
      ETerm *= EOnLambda;
    }
    const scalar Prefactor{ ( E0 + ModelParams[idxDelta[idxFFPlus]] )
                          / ( E0 + ModelParams[idxDelta[idxFF0]] ) };
    scalar Constraint = Prefactor * ESum[idxFF0] - ESum[idxFFPlus];
    // If I'm not using c0 for the plug, I need to scale appropriately
    if( WhichConstraint < NumE )
    {
      const scalar LambdaOnE{ 1. / EOnLambda };
      for( unsigned int i = 0; i <= WhichConstraint; ++i )
        Constraint *= LambdaOnE;
    }
    ModelParams[idxConstraint] = Constraint;
  }
}

void ContinuumFit::ParamCovarList( Common::Params &paramsCovar ) const
{
  if( ( uiFF & uiFF0 ) && ( uiFF & uiFFPlus ) )
  {
    const Common::Param::Key &k{ ConstraintKey() };
    paramsCovar.Add( k, 1, false, Common::Param::Type::Derived );
  }
}

const std::string &ContinuumFit::GetPoleMassName( Common::FormFactor ff,
                                                  const Common::FileNameAtt &fna )
{
  static const std::string PDGDStar{ "PDGDStar" };
  static const std::string PDGD0Star{ "PDGD0Star" };
  static const std::string PDGDsStar{ "PDGDsStar" };
  static const std::string PDGDs0Star{ "PDGDs0Star" };
  const bool bVector{ ff == Common::FormFactor::fplus || ff == Common::FormFactor::fperp };
  // Pole-mass determined by final-state quark (not spectator)
  if( std::toupper( fna.Quark[idxSnk][0] ) == 'L' )
    return bVector ? PDGDStar : PDGD0Star;
  return bVector ? PDGDsStar : PDGDs0Star;
}

void ContinuumFit::AddEnsemble( const std::string &Ensemble )
{
  if( Ensemble.empty() )
    throw( std::runtime_error( "AddEnsemble() Ensemble empty" ) );
  EnsembleMapT::iterator it{ EnsembleMap.find( Ensemble ) };
  if( it == EnsembleMap.end() )
  {
    EnsembleInfo ei;
    const int Prefix{ std::toupper( Ensemble[0] ) };
    switch( Prefix )
    {
      case 'C':
        ei.aInv_L = 24;
        ei.aInv_T = 64;
        break;
      case 'M':
        ei.aInv_L = 32;
        ei.aInv_T = 64;
        break;
      case 'F':
        ei.aInv_L = 48;
        ei.aInv_T = 96;
        break;
      default:
        throw( std::runtime_error( "AddEnsemble() Unrecognised ensemble " + Ensemble ) );
    }
    EnsembleMap.insert( typename EnsembleMapT::value_type( Ensemble, std::move( ei ) ) );
  }
}

void ContinuumFit::GetEnabled( std::string sOptions )
{
  // By default, everything is enabled unless mentioned
  for( int idxFF = 0; idxFF < NumFF; ++idxFF )
  {
    c0Enabled[idxFF] = true;
    c1Enabled[idxFF] = true;
    ChiEnabled[idxFF] = true;
    FVEnabled[idxFF] = true;
    dEnabled[idxFF].resize( NumD, true );
    idxD[idxFF].resize( NumD, idxUnused );
    kD[idxFF].resize( NumD );
    eEnabled[idxFF].resize( NumE, true );
    idxE[idxFF].resize( NumE, idxUnused );
    kE[idxFF].resize( NumE );
  }
  std::vector<Common::FormFactor> ffRepeats;
  while( Common::Trim( sOptions ) )
  {
    // Get form factor and check its not a repeat
    const std::string sThisFF{ Common::ExtractToSeparator( sOptions ) };
    Common::FormFactor ThisFF{ ValidateFF( Common::FromString<Common::FormFactor>( sThisFF ) ) };
    for( Common::FormFactor &ff : ffRepeats )
      if( ff == ThisFF )
      {
        std::ostringstream os;
        os << "ContinuumFit::GetEnabled() form factor " << sThisFF << " repeated";
        throw std::runtime_error( os.str().c_str() );
      }
    ffRepeats.push_back( ThisFF );
    // Now get the options and check they contain no repeats
    const std::string sDisabled{ Common::ExtractToSeparator( sOptions ) };
    for( std::size_t i = 1; i < sDisabled.length(); ++i )
    {
      const int c{ std::toupper( sDisabled[i] ) };
      for( std::size_t j = 0; j < i; ++j )
        if( std::toupper( sDisabled[j] ) == c )
        {
          std::ostringstream os;
          os << "ContinuumFit::GetEnabled() form factor " << sThisFF
             << " repeated options " << sDisabled;
          throw std::runtime_error( os.str().c_str() );
        }
    }
    // Now process options
    const int idxFF{ ffIndex( ThisFF ) };
    for( char c : sDisabled )
    {
      if( c == '0' )
      {
        c0Enabled[idxFF] = false;
        // I don't need either chiral or finite volume corrections either
        ChiEnabled[idxFF] = false;
        FVEnabled[idxFF] = false;
      }
      else if( c == '1' )
        c1Enabled[idxFF] = false;
      else if( std::toupper( c ) == 'X' )
        ChiEnabled[idxFF] = false;
      else if( std::toupper( c ) == 'V' )
        FVEnabled[idxFF] = false;
      else
      {
        unsigned int Num = static_cast<unsigned int>( c - '2' );
        if( Num < NumE )
          eEnabled[idxFF][Num] = false;
        else
        {
          Num = static_cast<unsigned int>( std::toupper( c ) - 'A' );
          if( Num < NumD )
            dEnabled[idxFF][Num] = false;
          else
          {
            std::ostringstream os;
            os << "ContinuumFit::GetEnabled() form factor " << sThisFF
               << " disable " << sDisabled << " unrecognised character " << c;
            throw std::runtime_error( os.str().c_str() );
          }
        }
      }
    }
  }
}

void ContinuumFit::LoadModels()
{
  std::cout << std::setprecision( std::numeric_limits<double>::max_digits10+2 )
            << "Performing chiral continuum fit from ensemble results" << std::endl;
  // Walk the list of parameters on the command-line, loading correlators and making models
  bool bFirst{ true };
  for( std::size_t ArgNum = 0; ArgNum < cl.Args.size(); ++ArgNum )
  {
    // First parameter is the filename we're looking for, optionally followed by comma and arguments
    std::string FileToGlob{ Common::ExtractToSeparator( cl.Args[ArgNum] ) };
    std::vector<std::string> Filenames{ Common::glob( &FileToGlob, &FileToGlob + 1,
                                                      inBase.c_str() ) };
    bool bGlobEmpty{ true };
    // Anything after the comma is a list of arguments
    Model::Args vArgs;
    vArgs.FromString( cl.Args[ArgNum], true );
    // If Ensemble is specified, I will override what I loaded in the file
    std::string Ensemble{ vArgs.Remove( Common::sEnsemble ) };
    // Make sure the default form factor is included as an argument to model (if not user-specified)
    Common::FormFactor thisFF;
    typename Model::Args::iterator itFF{ vArgs.find( sFF ) };
    if( itFF == vArgs.end() )
    {
      thisFF = ffDefault;
      vArgs.emplace( sFF, Common::GetFormFactorString( thisFF ) );
    }
    else
    {
      thisFF = Common::FromString<Common::FormFactor>( itFF->second );
    }
    // Load this batch of files
    for( const std::string &sFileName : Filenames )
    {
      // If there are any arguments, show them as I load each file
      std::string PrintPrefix( 2, ' ' );
      if( !cl.Args[ArgNum].empty() )
      {
        PrintPrefix.append( cl.Args[ArgNum] );
        PrintPrefix.append( 1, ' ' );
      }
      Common::FileNameAtt fna{ sFileName, &OpName };
      if( ( thisFF == Common::FormFactor::fplus || thisFF == Common::FormFactor::fperp )
         && !fna.HasNonZeroMomentum() )
        std::cout << PrintPrefix << "Ignoring zero momentum file " << sFileName << Common::NewLine;
      else
      {
        // Output meson name comes from the first file I load
        if( bFirst )
        {
          bFirst = false;
          if( fna.BaseShortParts.size() < 3 || fna.Spectator.empty() )
            throw std::runtime_error( "Meson names not included in filename " + sFileName );
          Meson[idxSnk] = Common::MesonName( fna.BaseShortParts[1], fna.Spectator );
          Meson[idxSrc] = Common::MesonName( fna.BaseShortParts[2], fna.Spectator );
          if( Common::EqualIgnoreCase( Meson[idxSnk], Meson[idxSrc] ) )
            throw std::runtime_error( "Snk and Src are both " + Meson[idxSnk] );
          bDoConstraint = thisFF == Common::FormFactor::f0 || thisFF == Common::FormFactor::fplus;
        }
        if( thisFF == Common::FormFactor::f0 || thisFF == Common::FormFactor::fpar )
          uiFF |= uiFF0;
        else
          uiFF |= uiFFPlus;
        bGlobEmpty = false;
        ds.LoadModel( std::move( fna ), PrintPrefix.c_str(),
                     Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT
                     | Common::COMPAT_DISABLE_ENSEMBLE );
        ModelArgs.emplace_back( vArgs );
        if( !Ensemble.empty() )
          ds.constFile.back()->Ensemble = std::move( Ensemble );
        AddEnsemble( ds.constFile.back()->Ensemble );
      }
    }
    if( bGlobEmpty )
      throw std::runtime_error( "No files matched " + FileToGlob );
    EnsembleMap.Loaded();
  }
}

void ContinuumFit::GetEnsembleStats()
{
  for( const ModelFilePtr &mp : ds.constFile )
  {
    typename EnsembleStatMap::iterator it{ EnsembleStats.find( mp->Ensemble ) };
    if( it == EnsembleStats.end() )
    {
      EnsembleStats.emplace( std::make_pair( mp->Ensemble, EnsembleStat{ 1, mp->Seed() } ) );
      std::cout << "SampleSize[" << mp->Ensemble << "] = " << mp->SampleSize << Common::NewLine;
    }
    else
    {
      EnsembleStat &es{ it->second };
      if( mp->Seed() != es.Seed )
        throw std::runtime_error( "GetEnsembleStats() seeds don't match on ensemble "
                                 + mp->Ensemble + " file " + std::to_string( es.Num )
                                 + Common::Space + mp->Name_.Filename );
      ++es.Num;
    }
  }
}

void ContinuumFit::SetEnsembleStats()
{
  if( EnsembleStats.size() > 1 )
  {
    std::string &Ens{ f->OutputModel.Ensemble };
    std::vector<Common::ConfigCount> &vcc{ f->OutputModel.ConfigCount };
    vcc.resize( EnsembleStats.size() );
    f->OutputModel.SampleSize = static_cast<int>( EnsembleStats.size() );
    Ens.clear();
    int i = 0;
    for( EnsembleStatMap::value_type v : EnsembleStats )
    {
      // Make sure binSize is correct
      if( i == 0 )
        f->OutputModel.binSize = v.second.Num;
      else
      {
        if( v.second.Num != f->OutputModel.binSize )
          f->OutputModel.binSize = 0;
      }
      // Add this ensemble name
      if( !Ens.empty() )
        Ens.append( 1, ' ' );
      Ens.append( v.first );
      // Update corresponding config count
      vcc[i].Config = i;
      vcc[i].Count = v.second.Num;
      ++i;
    }
  }
}

void ContinuumFit::SetEnsembleStatSeed()
{
  // Clear seed unless they all match
  for( EnsembleStatMap::value_type v : EnsembleStats )
  {
    if( v.second.Seed != f->OutputModel.Seed() )
    {
      f->OutputModel.SetSeed( 0 );
      break;
    }
  }
}

void ContinuumFit::SortModels()
{
  assert( ds.constFile.size() == ModelArgs.size() && "Should be one ModelArgs per ds.constFile" );
  if( ds.constFile.empty() )
    throw std::runtime_error( "Nothing to fit" );

  using FF = Common::FormFactor;

  struct ModelWithArgs
  {
    ModelFilePtr m;
    Model::Args Args;
    FF ff;
  };

  struct ModelWithArgsLess
  {
    const bool CovarBlock;
    ModelWithArgsLess( bool covarBlock ) : CovarBlock{covarBlock} {}
    bool operator()( const ModelWithArgs &lhs, const ModelWithArgs &rhs )
    {
      int i{ Common::CompareIgnoreCase( lhs.m->Ensemble, rhs.m->Ensemble ) };
      if( i )
        return i < 0;
      const int lp2{ lhs.m->Name_.GetFirstNonZeroMomentum().second.p2() };
      const int rp2{ rhs.m->Name_.GetFirstNonZeroMomentum().second.p2() };
      /*if( !CovarBlock )
      {
        if( lp2 != rp2 )
          return lp2 < rp2;
        return lhs.ff < rhs.ff;
      }*/
      if( lhs.ff != rhs.ff )
        return lhs.ff < rhs.ff;
      return lp2 < rp2;
    }
  };

  // Take the models and corresponding arguments out of the data set and put in my list
  std::vector<ModelWithArgs> ei( ds.constFile.size() );
  for( std::size_t i = 0; i < ds.constFile.size(); ++i )
  {
    ei[i].m.reset( ds.constFile[i].release() );
    ei[i].Args = std::move( ModelArgs[i] );
    ei[i].ff = Common::FromString<FF>( ei[i].Args.find( sFF )->second );
  }
  // Sort the models and corresponding arguments by ensemble
  std::sort( ei.begin(), ei.end(), ModelWithArgsLess( CovarBlock ) );
  // Put them back
  for( std::size_t i = 0; i < ds.constFile.size(); ++i )
  {
    ds.constFile[i].reset( ei[i].m.release() );
    ModelArgs[i] = std::move( ei[i].Args );
  }

  // Describe the number of replicas
  std::cout << "Using ";
  if( ds.NSamples == ds.MaxSamples )
    std::cout << "all ";
  else
    std::cout << "first " << ds.NSamples << " of ";
  std::cout << ds.MaxSamples << Common::Space << ds.front().getData().SeedTypeString()
            << " replicas";
  if( ds.NSamples != ds.MaxSamples )
    std::cout << " (all " << ds.MaxSamples << " for var/covar)";
  std::cout << std::endl;
  // Debug
  //for( const ModelFilePtr &m : ds.constFile )
    //std::cout << "  " << m->Name_.Filename << std::endl;
}

void ContinuumFit::LoadExtra()
{
  const std::vector<std::string> &vExtra{ cl.SwitchStrings( "model" ) };
  const std::string PrintPrefix( 2, ' ' );
  const unsigned int CompareFlags{ Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT
                                 | Common::COMPAT_DISABLE_ENSEMBLE };
  for( const std::string &sFileName : vExtra )
  {
    Common::FileNameAtt fna;
    fna.Filename = sFileName;
    fna.Type = "Model constants";
    if( ds.empty() )
    {
      fna.Ext = "h5";
      fna.bSeedNum = true;
      fna.Seed = 0;
    }
    else
    {
      const Common::FileNameAtt &S{ ds.front().Name_ };
      fna.Ext = S.Ext;
      fna.bSeedNum = S.bSeedNum;
      fna.Seed = S.Seed;
    }
    ds.LoadModel( std::move( fna ), "", PrintPrefix.c_str(), CompareFlags );
  }
}

void ContinuumFit::NormaliseData()
{
  static int NumCalled{};
  if( NumCalled++ )
    throw std::runtime_error( "NormaliseData() repeated" );
  using FF = Common::FormFactor;
  for( std::size_t i = 0; i < ModelArgs.size(); ++i )
  {
    const FF ff{ Common::FromString<FF>( ModelArgs[i].find( sFF )->second ) };
    if( ff == FF::fpar || ff ==FF::fperp )
    {
      static const std::string sErrorMsg{ "NormaliseData() find column" };
      // Remove scale in fpar and fperp
      ModelFile &mf{ *ds.constFile[i] };
      const Common::Param::Key kFF{ Common::GetFormFactorString( ff ) };
      const std::size_t idx{ mf.params.Find( kFF, sErrorMsg )->second() };
      const Common::Param::Key kaInv{ Common::Param::Key( mf.Ensemble, "aInv" ) };
      const JackBootColumn aInv{ ds.GetConstant( kaInv ) };
      for( std::size_t rep = JackBoot::idxCentral; rep != mf.NumSamples(); ++rep )
      {
        // For some reason there's a scale of GeV in my data ...
        scalar Factor = std::sqrt( aInv[rep] * 1e-9 );
        if( ff == FF::fperp )
          Factor = 1. / Factor;
        mf(rep,idx) *= Factor;
      }
      // Now update the stats
      ValWithEr &ve{ mf.SummaryData( static_cast<int>( idx ) ) };
      mf.getData().MakeStatistics( ve, idx );
    }
  }
}

void ContinuumFit::MakeOutputFilename()
{
  // Sort operator names (also renumbers references to these names)
  ds.SortOpNames( OpName );
  // Make sorted, concatenated list of operators in the fit for filenames
  sOpNameConcat = OpName[0];
  for( std::size_t i = 1; i < OpName.size(); i++ )
  {
    sOpNameConcat.append( 1, '_' );
    sOpNameConcat.append( OpName[i] );
  }
}

// Make base filename for output
std::string ContinuumFit::GetOutputFilename( unsigned int uiFF )
{
  std::string s{ outBaseFileName };
  // If the output base is given, but doesn't end in '/', then it already contains the correct name
  if( outBaseFileName.empty() || outBaseFileName.back() == '/' )
  {
    // Simplistic, but better is hard to automate
    const Common::FileNameAtt &fna{ ds.front().Name_ };
    if( fna.BaseShortParts.size() < 3 || fna.Spectator.empty() )
      throw std::runtime_error( "Please include output filename in -o" );
    s.append( fna.BaseShortParts[0] );
    s.append( 1, '_' );
    s.append( Meson[idxSnk] );
    s.append( 1, '_' );
    s.append( Meson[idxSrc] );
    for( std::size_t i = 3; i < fna.BaseShortParts.size(); ++i )
    {
      s.append( 1, '_' );
      s.append( fna.BaseShortParts[i] );
    }
  }
  s.append( 1, '.' );
  s.append( doCorr ? "corr" : "uncorr" );
  if( uiFF & uiFF0 )
  {
    s.append( 1, '_' );
    s.append( Common::GetFormFactorString( Common::FormFactor::f0 ) );
  }
  if( uiFF & uiFFPlus )
  {
    s.append( 1, '_' );
    s.append( Common::GetFormFactorString( Common::FormFactor::fplus ) );
  }
  s.append( 1, '.' );
  return s;
}

void ContinuumFit::ShowCovar()
{
  std::cout << "Ensemble    i    j rows cols Condition number\n";
  for( std::size_t i = 0; i < ds.mCovar.size1; )
  {
    // Find the width of this sub-matrix on the same Ensemble
    std::size_t width = 1;
    while( i + width < ds.mCovar.size1
          && Common::EqualIgnoreCase( ds.constFile[i]->Ensemble, ds.constFile[i+width]->Ensemble ) )
      ++width;
    // Show the condition number
    Vector S;
    const scalar Cond{ ds.mCovar.SubMatrix( i, i, width, width ).Cholesky( S ).CholeskyRCond() };
    std::cout << std::setw(8) << ds.constFile[i]->Ensemble
              << std::setw(5) << i
              << std::setw(5) << i
              << std::setw(5) << width
              << std::setw(5) << width
              << Common::Space << Cond
              << Common::NewLine;
    i += width;
  }
}

void ContinuumFit::MakeCovarBlock()
{
  std::cout << "Block Covariance\n";
  for( std::size_t i = 0; i < ds.mCovar.size1; ++i )
  {
    std::cout << Common::Space;
    const ModelContinuum &Mi{ * dynamic_cast<const ModelContinuum *>( f->model[i].get() ) };
    for( std::size_t j = 0; j <= i; ++j )
    {
      const ModelContinuum &Mj{ * dynamic_cast<const ModelContinuum *>( f->model[j].get() ) };
      const bool SameFF{ CovarBlock ? ( Mi.ff == Mj.ff ) : true };
      const bool SameEnsemble{ Common::EqualIgnoreCase( ds.constFile[i]->Ensemble,
                                                        ds.constFile[j]->Ensemble ) };
      const bool bZero{ i != j && !( SameEnsemble && SameFF ) };
      if( bZero )
      {
        ds.mCovar(i,j) = 0;
        ds.mCovar(j,i) = 0;
      }
      std::cout << ( bZero ? '0' : '1' );
    }
    std::cout << Common::NewLine;
  }
}

// Now get the min / max of the specified field
void ContinuumFit::GetMinMax( Common::FormFactor ff, scalar &Min, scalar &Max, int Loop,
                              const std::string &Field ) const
{
  // See whether Min / Max specified in the environment
  const std::string sEnvPrefix{ "MLU" + Common::GetFormFactorString( ff ) + Field };
  const std::string sEnvMin{ sEnvPrefix + "Min" };
  const std::string sEnvMax{ sEnvPrefix + "Max" };
  const char *pMin{ std::getenv( sEnvMin.c_str() ) };
  const char *pMax{ std::getenv( sEnvMax.c_str() ) };
  if( pMin || pMax )
  {
    static const char szEnv[] = "Environment variable ";
    static const char szMissing[] = " missing / empty";
    static const char szInvalid[] = " invalid: ";
    if( !pMin || !*pMin )
      throw std::runtime_error( szEnv + sEnvMin + szMissing );
    if( !pMax || !*pMax )
      throw std::runtime_error( szEnv + sEnvMax + szMissing );
    try
    {
      Min = Common::FromString<scalar>( pMin );
    }
    catch( const std::runtime_error &e )
    {
      throw std::runtime_error( szEnv + sEnvMin + szInvalid + e.what() );
    }
    try
    {
      Max = Common::FromString<scalar>( pMax );
    }
    catch( const std::runtime_error &e )
    {
      throw std::runtime_error( szEnv + sEnvMax + szInvalid + e.what() );
    }
    std::cout << Field
              << Common::CommaSpace << sEnvMin << "=" << pMin
              << Common::CommaSpace << sEnvMax << "=" << pMax
              << Common::NewLine;
    return;
  }
  const ModelFile &om{ f->OutputModel };
  bool bFirst{ true };
  for( std::size_t i = 0; i < f->model.size(); ++i )
  {
    const ModelContinuum &m{ * dynamic_cast<const ModelContinuum *>( f->model[i].get() ) };
    if( ValidateFF( m.ff ) == ff )
    {
      const ValWithEr &ve{ om.SummaryData( static_cast<int>( ( Loop == 0 ? m.qSq : m.EL ).idx ) ) };
      if( bFirst )
      {
        bFirst = false;
        if( !ve.Check )
        {
          std::ostringstream os;
          os << "ContinuumFit::GetMinMax() ValWithEr Check 0 for " << ff << Common::Space << Field;
          throw std::runtime_error( os.str().c_str() );
        }
        Min = ve.Central;
        Max = ve.Central;
      }
      else
      {
        if( Min > ve.Central )
          Min = ve.Central;
        if( Max < ve.Central )
          Max = ve.Central;
      }
    }
  }
  //if( f->Verbosity )
  {
    std::ostringstream os;
    os << std::scientific << std::setprecision( std::numeric_limits<scalar>::max_digits10 + 2 )
       << Min << " <= " << ff << '-' << Field << " <= " << Max << Common::NewLine;
    std::cout << os.str();
  }
}

std::ofstream ContinuumFit::WriteHeader(const std::string &sPrefix, const std::string &FileType) const
{
  const Common::FileNameAtt &fna{ f->OutputModel.Name_ };
  std::string NewType{ fna.Type };
  if( !FileType.empty() )
  {
    NewType.append( 1, '_' );
    NewType.append( FileType );
  }
  Common::FileNameAtt fnaNew;
  fnaNew.Parse( sPrefix, NewType, fna.bSeedNum, fna.Seed, TEXT_EXT );
  std::cout << "Make " << FileType << Common::Space << fnaNew.Filename << Common::NewLine;
  std::ofstream os( fnaNew.Filename );
  Common::SummaryHeader<scalar>( os, fnaNew.Filename );
  f->OutputModel.Sample::SummaryComments( os, true, false );
  os << "# Primary key: ensemble " << f->model[0]->XVectorKeyName() << Common::NewLine;
  os << "# x-units: eV\n";
  return os;
}

void ContinuumFit::WriteFieldName( std::ofstream &os, const std::string &FieldName ) const
{
  os << Common::Space << FieldName << "_witherr ";
  ValWithEr::Header( FieldName, os, Common::Space );
}

struct BootDat
{
  const int ErrorDigits;
  scalar Central;
  Common::Vector<scalar> Buffer;
  ValWithEr ve;
  BootDat( int errorDigits, int Size ) : ErrorDigits{errorDigits}, Buffer( Size ) {}
  inline scalar &operator[]( int idx )
  {
    if( idx == ModelFile::idxCentral )
      return Central;
    return Buffer[idx];
  }
  ValWithEr &Get( std::vector<scalar> &ScratchBuffer )
  {
    ve.Get( Central, Buffer, ScratchBuffer );
    return ve;
  }
};

std::ostream & operator<<( std::ostream &os, BootDat &bd )
{
  return os << Common::Space << bd.ve.to_string( bd.ErrorDigits ) << Common::Space << bd.ve;
}

void ContinuumFit::WriteFitQSq( Common::FormFactor ff, const std::string &sPrefix ) const
{
  const int idxff{ ffIndex( ff ) };

  // Make output file
  const ModelFile &om{ f->OutputModel };
  std::ofstream os{ WriteHeader( sPrefix, "fit" ) };
  
  // I'll need the central value of q^2_max
  const scalar PDGHCentral{ om(ModelFile::idxCentral,idxPDGH) };
  const scalar PDGLCentral{ om(ModelFile::idxCentral,idxPDGL) };
  scalar qSqMaxCentral{ PDGHCentral - PDGLCentral };
  qSqMaxCentral *= qSqMaxCentral;

  // Now write each data row
  ValWithEr ve[3];
  Common::Vector<scalar> Buffer( om.NumSamples() );
  Common::Vector<scalar> PoleBuffer( om.NumSamples() );
  Common::Vector<scalar> yNoPoleBuffer( om.NumSamples() );
  std::vector<scalar> ScratchBuffer( om.NumSamples() );

  scalar x = -999.999;
  for( int Loop = 0; Loop < 2; ++Loop )
  {
    const std::string &FieldName{ Loop == 0 ? FieldQSq : FieldEL };
    if( Loop )
      os << "\n\n";

    // Write header row - i.e. column names
    os << "field x";
    WriteFieldName( os, "y" );
    WriteFieldName( os, "Pole" );
    WriteFieldName( os, "yNoPole" );
    os << Common::NewLine;

    scalar Min, Max;
    GetMinMax( ff, Min, Max, Loop, FieldName );
    const scalar Tick{ ( Max - Min ) / NumTicks };
    scalar Central = -999.999; // value unused
    scalar CentPole = 0, CentYNoPole = 0; // value unused
  for( int nTick = -2; nTick <= NumTicks; ++nTick, x += Tick )
  {
    if( nTick == -2 )
      x = ( Loop == 0 ) ? 0 : EOfQSq(ModelFile::idxCentral, 0);
    else if( nTick == -1 )
      x = ( Loop == 0 ) ? qSqMaxCentral : PDGLCentral;
    if( nTick == 0 )
      x = Min;
    for( int rep = ModelFile::idxCentral; rep < om.NumSamples(); ++rep )
    {
      const scalar Delta{ om(rep,idxDelta[idxff]) };
      const scalar E{ Loop == 0 ? EOfQSq(rep, x) : x };
      const scalar EOnLambda{ E * LambdaInv };
      scalar yNoPole = c0Enabled[idxff] ? om(rep,idxC0[idxff]) : 0;
      scalar Factor = EOnLambda;
      for( unsigned int i = 0; i < NumE; ++i )
      {
        if( eEnabled[idxff][i] )
          yNoPole += om(rep,idxE[idxff][i]) * Factor;
        Factor *= EOnLambda;
      }
      const scalar Pole{ Lambda / ( E + Delta ) };
      const scalar Y{ yNoPole * Pole };
      if( rep == ModelFile::idxCentral )
      {
        Central = Y;
        CentYNoPole = yNoPole;
        CentPole = Pole;
      }
      else
      {
        Buffer[rep] = Y;
        yNoPoleBuffer[rep] = yNoPole;
        PoleBuffer[rep] = Pole;
      }
    }
    ve[0].Get( Central, Buffer, ScratchBuffer );
    ve[1].Get( CentPole, PoleBuffer, ScratchBuffer );
    ve[2].Get( CentYNoPole, yNoPoleBuffer, ScratchBuffer );
    if( nTick < 0 )
      os << "# ";
    os << FieldName << Common::Space << x;
    for( int i = 0; i < sizeof(ve) / sizeof(ve[0]); ++i )
      os << Common::Space << ve[i].to_string( f->ErrorDigits ) << Common::Space << ve[i];
    if( nTick == -2 )
      os << " # q_0";
    else if( nTick == -1 )
      os << " # q^2_max";
    os << Common::NewLine;
  }
  }
}

void ContinuumFit::WriteAdjustedQSq( Common::FormFactor ff, const std::string &sPrefix ) const
{
  const ModelFile &om{ f->OutputModel };
  std::ofstream os{ WriteHeader( sPrefix, "" ) };
  // Write header row - i.e. column names
  os << "model ensemble nSq";
  WriteFieldName( os, FieldQSq );
  WriteFieldName( os, FieldEL );
  WriteFieldName( os, "data" );
  WriteFieldName( os, "adjusted" );
  WriteFieldName( os, "dataNoPole" );
  WriteFieldName( os, "adjustedNoPole" );
  WriteFieldName( os, "Pole" );
  os << Common::NewLine;

  // Now write each data row
  std::vector<scalar> ScratchBuffer( om.NumSamples() );
  BootDat Adjusted( f->ErrorDigits, om.NumSamples() );
  BootDat DataNoPole( f->ErrorDigits, om.NumSamples() );
  BootDat AdjNoPole( f->ErrorDigits, om.NumSamples() );
  BootDat Pole( f->ErrorDigits, om.NumSamples() );
  std::size_t RowNum{};
  for( std::size_t i = 0; i < om.FitTimes.size(); ++i )
  {
    ModelFile &mf{ *ds.constFile[i] };
    const ModelContinuum &m{ * dynamic_cast<const ModelContinuum *>( f->model[i].get() ) };
    if( ValidateFF( m.ff ) == ff )
    {
      const int idxFF{ ffIndex( ff ) };
    // Adjust the q^2 value on each replica
    for( int rep = ModelFile::idxCentral; rep < om.NumSamples(); ++rep )
    {
      const scalar rmPi{ om(rep,idxmPi[m.ei.idx]) };
      const scalar rmPDGPi{ om(rep,idxmPDGPi) };
      const scalar aInv{ om(rep,idxaInv[m.ei.idx]) };
      // Compute the Chiral term
      const scalar TermChiral = ( c0Enabled[idxFF] && ( NeedFV() || NeedChiral() ) )
                                ? om(rep,idxC0[idxFF]) * om(rep,idxChiFV[m.ei.idx]) : 0;
      // Compute the MPi term
      const scalar TermMPi = c1Enabled[idxFF]
              ? om(rep,idxC1[idxFF]) * ( rmPi * rmPi - rmPDGPi * rmPDGPi ) * LambdaInvSq : 0;
      // Compute the Discretisation term
      scalar TermDisc = 0;
      const scalar aLambda{ Lambda / aInv };
      const scalar aLambdaSq{ aLambda * aLambda };
      scalar Factor = aLambdaSq;
      for( unsigned int j = 0; j < NumD; ++j )
      {
        if( dEnabled[idxFF][j] )
          TermDisc += om(rep,idxD[idxFF][j]) * Factor;
        Factor *= aLambdaSq;
      }
      // Get adjusted form factor
      const scalar PoleTerm{ Lambda / ( om(rep,m.EL.idx) + om(rep,idxDelta[idxFF]) ) };
      const scalar Adjust = PoleTerm * ( TermChiral + TermMPi + TermDisc );
      const scalar FFAdjust = om.FitInput(rep,i) - Adjust;
      const scalar InvPoleTerm{ 1. / PoleTerm };
      Adjusted[rep] = FFAdjust;
      Pole[rep] = PoleTerm;
      DataNoPole[rep] = om.FitInput(rep,i) * InvPoleTerm;
      AdjNoPole[rep] = Adjusted[rep] * InvPoleTerm;
    }
    Adjusted.Get( ScratchBuffer );
    DataNoPole.Get( ScratchBuffer );
    AdjNoPole.Get( ScratchBuffer );
    Pole.Get( ScratchBuffer );
    const ValWithEr &veQ{ om.SummaryData( static_cast<int>( m.qSq.idx ) ) };
    const ValWithEr &veEL{ om.SummaryData( static_cast<int>( m.EL.idx ) ) };
    const ValWithEr &F{ mf.SummaryData( static_cast<int>( m.idxFitColumn ) ) };
    os << RowNum++ << Common::Space << mf.Ensemble << Common::Space << m.XVectorKey()
       << Common::Space << veQ .to_string( f->ErrorDigits ) << Common::Space << veQ
       << Common::Space << veEL.to_string( f->ErrorDigits ) << Common::Space << veEL
       << Common::Space << F.to_string( f->ErrorDigits ) << Common::Space << F
       << Adjusted << DataNoPole << AdjNoPole << Pole
       << Common::NewLine;
  }
  }
}

void ContinuumFit::WriteSynthetic()
{
  for( int idxFF = 0; idxFF < NumFF; ++idxFF )
  {
    const unsigned int uiMyFF{ ffMaskFromIndex( idxFF ) };
    if( uiFF & uiMyFF )
    {
      const Common::FormFactor ff{ ffIndexReverse( idxFF ) };
      const std::string sPrefix{ GetOutputFilename( uiMyFF ) + sOpNameConcat };
      WriteFitQSq( ff, sPrefix );
      WriteAdjustedQSq( ff, sPrefix );
    }
  }
}

int ContinuumFit::DoFit()
{
  int iReturn{ EXIT_SUCCESS };

  // Set fit ranges
  std::vector<std::vector<int>> fitTimes( ModelArgs.size() );
  for( int i = 0; i < f->ModelArgs.size(); ++i )
    fitTimes[i] = { static_cast<int>( f->model[i]->GetFitColumn() ) };
  ds.SetFitTimes( fitTimes, false );
  ShowCovar();
  MakeCovarBlock();

  try
  {
    {
      std::stringstream ss;
      if( f->bTestRun )
        ss << "Test run of ";
      ss << (doCorr ? "C" : "Unc") << "orrelated " << f->Type() << " chiral continuum fit";
      Fitter::SayNumThreads( ss );
      const std::string &sMsg{ ss.str() };
      if( !f->bTestRun )
        std::cout << std::string( sMsg.length(), '=' ) << std::endl;
      std::cout << sMsg << std::endl;
    }
    double ChiSq;
    int dof;
    if( !f->PerformFit( doCorr, ChiSq, dof, GetOutputFilename( uiFF ), sOpNameConcat ) )
      iReturn = 3; // Not all parameters resolved
    SetEnsembleStatSeed();
    WriteSynthetic();
  }
  catch(const std::exception &e)
  {
    iReturn = 2;
    std::cout << "Fitter error: " << e.what() << std::endl;
  }
  return iReturn;
}

int ContinuumFit::Run()
{
  LoadModels();
  GetEnsembleStats();
  SortModels();
  MakeOutputFilename();
  LoadExtra();
  NormaliseData();
  CovarParams cp{ cl, ds };
  std::cout << cp << std::endl;
  std::cout << "Covariance matrix entries across form factor (same ensemble) "
            << ( CovarBlock ? "retained" : "set to zero" ) << std::endl;

  // Make the fitter
  f.reset( Fitter::Make( CreateParams{ OpName, cl, *this },
                         ds, std::move(ModelArgs), std::move(cp), false, *this ) );
  SetEnsembleStats();
  // Now do the fit
  return DoFit();
}

int main(int argc, const char *argv[])
{
#ifndef HAVE_MINUIT2
  std::ios_base::sync_with_stdio( false );
#endif
  static const char DefaultNumEnergies[] = "3";
  static const char DefaultNumDiscret[] = "1";
  static const char DefaultFormFactor[] = "f0";
  static const char DefaultEnergySep[] = "0";
  static const char DefaultHotelling[] = "0.05";
  static const char DefaultErrDig[] = "2";
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"d", CL::SwitchType::Single, DefaultNumDiscret},
      {"e", CL::SwitchType::Single, DefaultNumEnergies},
      {"f", CL::SwitchType::Single, DefaultFormFactor},
      {"overwrite", CL::SwitchType::Flag, nullptr},
      {"poles", CL::SwitchType::Multiple, nullptr},
      {"polev", CL::SwitchType::Multiple, nullptr},
      {"model", CL::SwitchType::Multiple, nullptr},
      {"Hotelling", CL::SwitchType::Single, DefaultHotelling},
      {"chisqdof", CL::SwitchType::Single, "0"},
      {"sep", CL::SwitchType::Single, DefaultEnergySep},
      {"mindof", CL::SwitchType::Single, "1"},
      {"retry", CL::SwitchType::Single, "0"},
      {"iter", CL::SwitchType::Single, "0"},
      {"tol", CL::SwitchType::Single, "1e-7"},
      {"block", CL::SwitchType::Flag, nullptr},
      {"summary", CL::SwitchType::Single, "1"},
      {"nopolap", CL::SwitchType::Single, ""},
      {"v", CL::SwitchType::Single, "0"},
      {"strict",  CL::SwitchType::Single, "0"},
      {"maxE",  CL::SwitchType::Single, "10"},
      {"errdig", CL::SwitchType::Single, DefaultErrDig},
      // Covariance parameters
      {"covsrc", CL::SwitchType::Single, nullptr},
      {"covboot", CL::SwitchType::Single, nullptr},
      {"freeze", CL::SwitchType::Flag, nullptr},
      // Other params
      {"fitter", CL::SwitchType::Single, "GSL"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"debug-signals", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      bShowUsage = false;
      if( cl.GotSwitch( "debug-signals" ) )
        Common::Grid_debug_handler_init();
      const auto Start{ std::chrono::steady_clock::now() };
      ContinuumFit cf{ cl };
      iReturn = cf.Run();
      const auto Stop{ std::chrono::steady_clock::now() };
      // Say whether this worked, what the time is and how long it took
      switch( iReturn )
      {
        case 0:
          std::cout << "OK";
          break;
        case 3:
          std::cout << "Warning: not all parameters resolved";
          break;
        default:
          std::cout << "Error " << iReturn;
          break;
      }
      const auto mS{ std::chrono::duration_cast<std::chrono::milliseconds>( Stop - Start ) };
      std::cout << ". Fit duration ";
      if( mS.count() < 1100 )
        std::cout << mS.count() << " milliseconds";
      else
      {
        const auto S{ std::chrono::duration_cast<std::chrono::duration<double>>( Stop - Start ) };
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1) << S.count() << " seconds";
        std::cout << ss.str();
      }
      std::time_t ttNow;
      std::time( &ttNow );
      std::string sNow{ std::ctime( &ttNow ) };
      while( sNow.length() && sNow.back() == '\n' )
        sNow.resize( sNow.length() - 1 );
      std::cout << ". " << sNow << std::endl;
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
    " [Options] Model[,params]...\n"
    "Perform a chiral continuum fit of the per Ensemble data\n"
    "Options:\n"
    "-d       Number of (a Lambda)^n terms ( default: " << DefaultNumDiscret << ")\n"
    "-e       Number of (E_L/Lambda)^n terms ( default: " << DefaultNumEnergies << ")\n"
    "-f       ff,constants[,ff,constants] (eg f0,3,fplus,31)\n"
    "         Disable constants: 0=c0; 1=c1; 2,3,...=e0,e1,...; a,b,...=d0,d1,...;\n"
    "                            X=chiral correction; V=finite volume correction\n"
    "         1st form factor becomes default (default: " << DefaultFormFactor << ")\n"
    "--poles  Scalar pole mass [eV]\n"
    "--polev  Vector pole mass [eV]\n"
    "--model  An additional model file, e.g. to load lattice spacings\n"
    "--Hotelling Minimum Hotelling Q-value on central replica (default " << DefaultHotelling << ")\n"
    "--chisqdof  Maximum chi^2 / dof on central replica\n"
    "--sep    Minimum relative separation between energies (default " << DefaultEnergySep << ")\n"
    "--retry  Maximum number of times to retry fits (default Minuit2=10, GSL=0)\n"
    "--iter   Max iteration count, 0 (default) = unlimited\n"
    "--tol    Tolerance of required fits (default 1e-7)\n"
    "--mindof Minimum degrees of freedom (default 1)\n"
    "--fitter (GSL|Minuit2)[,options] fitter (default GSL,Levenberg-Marquardt)\n"
    "         GSL options: lm, lmaccel, dogleg, ddogleg, subspace2D\n"
    "--strict Mask. Check field shows whether parameters are non-zero and unique\n"
    "         Strictness bits: Off=+/- 1 sigma; On=every replica (default 0)\n"
    "         Bit 0: Difference from 0\n"
    "         Bit 1: Difference from other parameters in same series\n"
    "--maxE   Maximum energy (default 10 - decays so fast effectively undetermined)\n"
    "--errdig Number of significant figures in error (default " << DefaultErrDig << ")\n"
    "--covsrc source[,options[,...]] build (co)variance from source, i.e. one of\n"
    "         Binned    Use the already binned data\n"
    "         Bootstrap Use bootstrap replicas\n"
    "         Raw       Use the raw (unbinned) data\n"
    "         Rebin n[,n[...] Rebin raw data using bin size(s) specified, 0=auto\n"
    "         Reboot[NumBoot] Re-bootstrap binned data, optional NumBoot samples\n"
    "                         Can be used with rebin\n"
    "         H5,f[,g],d Load INVERSE covariance from .h5 file f, group g, dataset d\n"
    "         Default: Binned if available, otherwise Bootstrap\n"
    "Flags:\n"
    "--block   Zero covariance entries for different form factors on each ensemble\n"
    "--covboot Build covariance using secondary bootstrap & this num replicas, 0=all\n"
    "--summary 0 no summaries; 1 model_td.seq.txt only; 2 model_td and model.seq.txt\n"
    "--nopolap List of overlap coefficients which don't depend on momenta\n"
    "          TODO: This is ignored for now. Should be moved into MultiFit only\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-v     Verbosity, 0 (default)=central, 1=detail, 2=all, 3=all detail\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--debug-signals Trap signals (code courtesy of Grid)\n"
    "--help     This message\n";
  }
  Common::Grid_exit_handler_disable = true;
  return iReturn;
}
