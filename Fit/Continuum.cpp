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

ContinuumFit::ContinuumFit( Common::CommandLine &cl_ )
: cl{cl_},
  NumSamples{cl.SwitchValue<int>("n")},
  doCorr{ !cl.GotSwitch( "uncorr" ) },
  CovarBlock{ cl.GotSwitch( "block" ) },
  ffDefault{ ValidateFF( Common::FromString<Common::FormFactor>( Common::TrimAt( cl.SwitchValue<std::string>("f") ) ) ) },
  cEnabled{ GetEnabled( cl.SwitchValue<std::string>("f") ) },
  inBase{ cl.SwitchValue<std::string>("i") },
  outBaseFileName{ cl.SwitchValue<std::string>("o") },
  ds{NumSamples},
  uiFF{0}
{
  Common::MakeAncestorDirs( outBaseFileName );
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
  for( int idxFF = 0; idxFF < idxC.size(); ++idxFF )
  {
    if( uiFF & ffMaskFromIndex( idxFF ) )
    {
      const Common::FormFactor ff{ ffIndexReverse( idxFF ) };
      kPDGDStar[idxFF].Name = GetPoleMassName( ff, f.ds.constFile[0]->Name_ );
      mp.Add( kPDGDStar[idxFF], 1 );
      kDelta[idxFF].Object.push_back( GetFormFactorString( ff ) );
      kDelta[idxFF].Name = "Delta";
      mp.Add( kDelta[idxFF], 1, false, Param::Type::Derived );
      // Add constants
      for( int i = 0; i < idxC[idxFF].size(); ++i )
      {
        if( CNeeded( idxFF, i ) )
        {
          kC[idxFF][i].Object.push_back( GetFormFactorString( ff ) );
          kC[idxFF][i].Name = "c" + std::to_string( i );
          mp.Add( kC[idxFF][i] );
        }
      }
    }
  }
  if( ( uiFF & uiFF0 ) && ( uiFF & uiFFPlus ) )
  {
    // If I'm fitting both form factors, I can impose a constraint on fplus-c0
    mp.SetType( kC[idxFFPlus][CChiral], Param::Type::Derived );
  }

  // Per ensemble parameters
  const bool bChiral{ ( uiFF & uiFF0 && CEnabled( idxFF0, CChiral ) )
                   || ( uiFF & uiFFPlus && CEnabled( idxFFPlus, CChiral ) ) };
  aInv.resize( EnsembleMap.size() );
  mPi.resize( EnsembleMap.size() );
  kaInv.resize( EnsembleMap.size() );
  kmPi.resize( EnsembleMap.size() );
  if( bChiral )
  {
    FVSim.resize( EnsembleMap.size() );
    FVPhys.resize( EnsembleMap.size() );
    kFVSim.resize( EnsembleMap.size() );
    kFVPhys.resize( EnsembleMap.size() );
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
    if( bChiral )
    {
      kFVSim[i].Object.push_back( sEnsemble );
      kFVSim[i].Name = "FVSim";
      mp.Add( kFVSim[i], 1, false, Param::Type::Derived );
      kFVPhys[i].Object.push_back( sEnsemble );
      kFVPhys[i].Name = "FVPhys";
      mp.Add( kFVPhys[i], 1, false, Param::Type::Derived );
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

  // Per form factor parameter offsets
  for( int idxFF = 0; idxFF < idxC.size(); ++idxFF )
  {
    idxDelta[idxFF] = idxCUnused;
    idxPDGDStar[idxFF] = idxCUnused;
    for( int i = 0; i < idxC[idxFF].size(); ++i )
      idxC[idxFF][i] = idxCUnused;
    if( uiFF & ffMaskFromIndex( idxFF ) )
    {
      idxPDGDStar[idxFF] = mp.at( kPDGDStar[idxFF] )();
      idxDelta[idxFF] = mp.at( kDelta[idxFF] )();
      for( int i = 0; i < idxC[idxFF].size(); ++i )
        if( CNeeded( idxFF, i ) )
          idxC[idxFF][i] = mp.at( kC[idxFF][i] )();
    }
  }

  // Per ensemble parameter offsets
  const bool bChiral{ ( uiFF & uiFF0 && CEnabled( idxFF0, CChiral ) )
                   || ( uiFF & uiFFPlus && CEnabled( idxFFPlus, CChiral ) ) };
  for( std::size_t i = 0; i < kaInv.size(); ++i )
  {
    aInv[i] = mp.at( kaInv[i] )();
    mPi[i] = mp.at( kmPi[i] )();
    if( bChiral )
    {
      FVSim[i] = mp.at( kFVSim[i] )();
      FVPhys[i] = mp.at( kFVPhys[i] )();
    }
  }
}

void ContinuumFit::SetReplica( Vector &ModelParams ) const
{
  const scalar mPDGL{ ModelParams[idxPDGL] };
  const scalar mPDGH{ ModelParams[idxPDGH] };
  // Compute delta for each form factor
  for( int idxFF = 0; idxFF < idxC.size(); ++idxFF )
  {
    if( uiFF & ffMaskFromIndex( idxFF ) )
    {
      const scalar PoleMass{ ModelParams[idxPDGDStar[idxFF]] };
      ModelParams[idxDelta[idxFF]] = 0.5 * ( ( PoleMass*PoleMass - mPDGL*mPDGL ) / mPDGH - mPDGH );
    }
  }
  // Compute finite volume corrections for each ensemble
  if( !FVSim.empty() )
  {
    for( const typename EnsembleMapT::value_type &ei : EnsembleMap )
    {
      const std::size_t i{ ei.second.idx };
      const unsigned int aInv_L{ ei.second.aInv_L };
      const scalar L{ aInv_L / ModelParams[aInv[i]] };
      ModelParams[FVSim[i]] = Common::DeltaF::FiniteVol( ModelParams[mPi[i]], L, aInv_L );
      ModelParams[FVPhys[i]] = Common::DeltaF::FiniteVol( ModelParams[idxmPDGPi], L, aInv_L );
    }
  }
}

void ContinuumFit::ComputeDerived( Vector &ModelParams ) const
{
  const scalar mPDGL{ ModelParams[idxPDGL] };
  const scalar mPDGH{ ModelParams[idxPDGH] };
  if( ( uiFF & uiFF0 ) && ( uiFF & uiFFPlus ) )
  {
    // If I'm fitting both form factors, I can impose a constraint on fplus-c0
    const scalar E0{ EOfQSq( mPDGH, mPDGL, 0 ) };
    const scalar EOnLambda{ E0 * LambdaInv };
    const scalar Prefactor{ ( E0 + ModelParams[idxDelta[idxFFPlus]] )
                          / ( E0 + ModelParams[idxDelta[idxFF0]] ) };
    scalar z = ModelParams[idxC[idxFF0][CChiral]];
    if( idxC[idxFF0][CEOnL] != idxCUnused )
      z += ModelParams[idxC[idxFF0][CEOnL]] * EOnLambda;
    if( idxC[idxFF0][CEOnL2] != idxCUnused )
      z += ModelParams[idxC[idxFF0][CEOnL2]] * EOnLambda * EOnLambda;
    if( idxC[idxFF0][CEOnL3] != idxCUnused )
      z += ModelParams[idxC[idxFF0][CEOnL3]] * EOnLambda * EOnLambda * EOnLambda;
    z *= Prefactor;
    if( idxC[idxFFPlus][CEOnL] != idxCUnused )
      z -= ModelParams[idxC[idxFFPlus][CEOnL]] * EOnLambda;
    if( idxC[idxFFPlus][CEOnL2] != idxCUnused )
      z -= ModelParams[idxC[idxFFPlus][CEOnL2]] * EOnLambda * EOnLambda;
    if( idxC[idxFFPlus][CEOnL3] != idxCUnused )
      z -= ModelParams[idxC[idxFFPlus][CEOnL3]] * EOnLambda * EOnLambda * EOnLambda;
    ModelParams[idxC[idxFFPlus][CChiral]] = z;
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
  if( std::toupper( fna.Quark[1][0] ) == 'L' )
    return bVector ? PDGDStar : PDGD0Star;
  return bVector ? PDGDsStar : PDGDs0Star;
}

Common::FormFactor ContinuumFit::ValidateFF( Common::FormFactor ff )
{
  if( ff != Common::FormFactor::fplus && ff != Common::FormFactor::f0 )
  {
    std::ostringstream os;
    os << "Unsupported form factor " << ff;
    throw std::runtime_error( os.str().c_str() );
  }
  return ff;
}

void ContinuumFit::AddEnsemble( const std::string &Ensemble, Common::FormFactor thisFF )
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
  EnsembleFFs.insert( EnsembleFF( Ensemble, thisFF ) );
}

std::array<std::array<bool, ContinuumFit::NumConst>, ContinuumFit::NumFF>
ContinuumFit::GetEnabled( std::string sOptions )
{
  std::array<std::array<bool, NumConst>, NumFF> cEnabled;
  // By default, everything is enabled unless mentioned
  for( std::size_t i = 0; i < NumFF; ++i )
    for( std::size_t j = 0; j < NumConst; ++j )
      cEnabled[i][j] = true;
  while( Common::Trim( sOptions ) )
  {
    const std::string sThisFF{ Common::ExtractToSeparator( sOptions ) };
    Common::FormFactor ThisFF{ ValidateFF( Common::FromString<Common::FormFactor>( sThisFF ) ) };
    const int idxFF{ ffIndex( ThisFF ) };
    for( std::size_t j = 0; j < NumConst; ++j )
      if( !cEnabled[idxFF][j] )
      {
        std::ostringstream os;
        os << "ContinuumFit::GetEnabled() form factor " << sThisFF << " repeated";
        throw std::runtime_error( os.str().c_str() );
      }
    const std::string sDisabled{ Common::ExtractToSeparator( sOptions ) };
    for( std::size_t i = 0; i < sDisabled.length(); ++i )
    {
      const int Num{ sDisabled[i] - '0' };
      if( Num < 0 || Num >= NumConst )
      {
        std::ostringstream os;
        os << "ContinuumFit::GetEnabled() form factor " << sThisFF << " can't disable C" << Num;
        throw std::runtime_error( os.str().c_str() );
      }
      if( !cEnabled[idxFF][Num] )
      {
        std::ostringstream os;
        os << "ContinuumFit::GetEnabled() form factor " << sThisFF << " repeat disable C" << Num;
        throw std::runtime_error( os.str().c_str() );
      }
      cEnabled[idxFF][Num] = false;
    }
  }
  return cEnabled;
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
      vArgs.emplace( sFF, Common::GetFormFactorString( ffDefault ) );
    }
    else
    {
      thisFF = ValidateFF( Common::FromString<Common::FormFactor>( itFF->second ) );
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
        }
        if( thisFF == Common::FormFactor::f0 )
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
        AddEnsemble( ds.constFile.back()->Ensemble, thisFF );
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
      ++it->second.Num;
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
    return;
  }
  const ModelFile &om{ f->OutputModel };
  bool bFirst{ true };
  for( std::size_t i = 0; i < f->model.size(); ++i )
  {
    const ModelContinuum &m{ * dynamic_cast<const ModelContinuum *>( f->model[i].get() ) };
    if( m.ff == ff )
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
  ValWithEr ve;
  Common::Vector<scalar> Buffer( om.NumSamples() );
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
    os << Common::NewLine;

    scalar Min, Max;
    GetMinMax( ff, Min, Max, Loop, FieldName );
    const scalar Tick{ ( Max - Min ) / NumTicks };
    scalar Central = -999.999; // value unused
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
      scalar z = om(rep,idxC[idxff][CChiral]); // "unused" means it's not chiral log adjusted
      if( idxC[idxff][CEOnL] != idxCUnused )
        z += om(rep,idxC[idxff][CEOnL]) * EOnLambda;
      if( idxC[idxff][CEOnL2] != idxCUnused )
        z += om(rep,idxC[idxff][CEOnL2]) * EOnLambda * EOnLambda;
      if( idxC[idxff][CEOnL3] != idxCUnused )
        z += om(rep,idxC[idxff][CEOnL3]) * EOnLambda * EOnLambda * EOnLambda;
      z *= Lambda / ( E + Delta );
      if( rep == ModelFile::idxCentral )
        Central = z;
      else
        Buffer[rep] = z;
    }
    ve.Get( Central, Buffer, ScratchBuffer );
    if( nTick < 0 )
      os << "# ";
    os << FieldName << Common::Space << x << Common::Space
       << ve.to_string( f->ErrorDigits ) << Common::Space << ve;
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
  os << Common::NewLine;

  std::vector<ValWithEr> F; // form factor data values with errors
  om.FitInput.MakeStatistics( F );

  // Now write each data row
  scalar Central = 0;
  ValWithEr ve;
  Common::Vector<scalar> Buffer( om.NumSamples() );
  std::vector<scalar> ScratchBuffer( om.NumSamples() );
  std::size_t RowNum{};
  for( std::size_t i = 0; i < om.FitTimes.size(); ++i )
  {
    ModelFile &mf{ *ds.constFile[i] };
    const ModelContinuum &m{ * dynamic_cast<const ModelContinuum *>( f->model[i].get() ) };
    if( m.ff == ff )
    {
      const int idxFF{ ffIndex( ff ) };
    // Adjust the q^2 value on each replica
    for( int rep = ModelFile::idxCentral; rep < om.NumSamples(); ++rep )
    {
      const scalar rmPi{ om(rep,mPi[m.ei.idx]) };
      const scalar rmPDGPi{ om(rep,idxmPDGPi) };
      // Compute the Chiral term
      const scalar sFVSim{ CEnabled( idxFF, CChiral ) ? om(rep,FVSim[idxFF]) : 0 };
      const scalar sFVPhys{ CEnabled( idxFF, CChiral ) ? om(rep,FVPhys[idxFF]) : 0 };
      const scalar Num = m.DeltaF( rmPi, sFVSim ) - m.DeltaF( rmPDGPi, sFVPhys );
      const scalar Denom = FourPi * om(rep,idxfPi);
      const scalar TermChiral = om(rep,idxC[idxFF][CChiral]) * ( Num / ( Denom * Denom ) );
      // Compute the MPi term
      scalar TermMPi = 0;
      if( CEnabled( idxFF, CMPi ) )
        TermMPi = om(rep,idxC[idxFF][CMPi]) * ( rmPi * rmPi - rmPDGPi * rmPDGPi ) * LambdaInvSq;
      // Compute the c4 term
      scalar TermDiscret = 0;
      if( CEnabled( idxFF, CDiscret ) )
      {
        const scalar aLambda{ Lambda / om(rep,aInv[m.ei.idx]) };
        TermDiscret = om(rep,idxC[idxFF][CDiscret]) * aLambda * aLambda;
      }
      // Get adjusted form factor
      const scalar PoleTerm{ Lambda / ( om(rep,m.EL.idx) + om(rep,idxDelta[idxFF]) ) };
      const scalar Adjust = PoleTerm * ( TermChiral + TermMPi + TermDiscret );
      const scalar FFAdjust = om.FitInput(rep,i) - Adjust;
      if( rep == ModelFile::idxCentral )
        Central = FFAdjust;
      else
        Buffer[rep] = FFAdjust;
    }
    ve.Get( Central, Buffer, ScratchBuffer );
    const ValWithEr &veQ{ om.SummaryData( static_cast<int>( m.qSq.idx ) ) };
    const ValWithEr &veEL{ om.SummaryData( static_cast<int>( m.EL.idx ) ) };
    os << RowNum++ << Common::Space << mf.Ensemble << Common::Space << m.XVectorKey()
       << Common::Space << veQ .to_string( f->ErrorDigits ) << Common::Space << veQ
       << Common::Space << veEL.to_string( f->ErrorDigits ) << Common::Space << veEL
       << Common::Space << F[i].to_string( f->ErrorDigits ) << Common::Space << F[i]
       << Common::Space << ve  .to_string( f->ErrorDigits ) << Common::Space << ve
       << Common::NewLine;
  }
  }
}

void ContinuumFit::WriteSynthetic()
{
  for( int i = 0; i < idxC.size(); ++i )
  {
    const unsigned int uiMyFF{ ffMaskFromIndex( i ) };
    if( uiFF & uiMyFF )
    {
      const Common::FormFactor ff{ ffIndexReverse( i ) };
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
      {"f", CL::SwitchType::Single, DefaultFormFactor},
      {"overwrite", CL::SwitchType::Flag, nullptr},
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
    "-f       ff,constants[,ff,constants] (eg f0,3,fplus,31)\n"
    "         Disable constants (0-4, 0=FV corrections) for each form factor\n"
    "         1st form factor becomes default (default: " << DefaultFormFactor << ")\n"
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
