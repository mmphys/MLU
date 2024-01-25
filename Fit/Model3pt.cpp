/*************************************************************************************
 
 Models for 3pt functions

 Source file: Model3pt.cpp
 
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

#include "Model3pt.hpp"

Model3pt::Model3pt( const Model::CreateParams &cp, Model::Args &Args, int NumExponents_,
                    std::vector<std::string> &&objectID, std::vector<std::string> &&opNames )
: ModelOverlap( cp, Args, NumExponents_, NumExponents_ > 1 ? 2 : 1,
                std::move( objectID ), std::move( opNames ) ),
N{ dynamic_cast<const MultiFitCreateParams &>( cp ).N },
dispType{ dynamic_cast<const MultiFitCreateParams &>( cp ).dispType }
{
  if( !cp.pCorr->Name_.GotDeltaT() )
    throw std::runtime_error( "DeltaT not available in " + cp.pCorr->Name_.Filename );
  DeltaT = cp.pCorr->Name_.DeltaT[0];
  if( NumExponents > 3 )
    throw std::runtime_error( "Model3pt supports a maximum of 3 exponentials: "
                              "Gnd-Gnd, Gnd-Ex (and ^\\dag), Ex-Ex" );
  std::string ESrc = Args.Remove( "ESrc", Common::ModelBase::EnergyPrefix );
  std::string ESnk = Args.Remove( "ESnk", Common::ModelBase::EnergyPrefix );
  const bool bESame{ objectID.size() == 1 && Common::EqualIgnoreCase( ESrc, ESnk ) };
  E.resize( bESame ? 1 : 2 );
  E[idxSrc].Key.Object = { ObjectID( idxSrc ) };
  E[idxSrc].Key.Name = std::move( ESrc );
  std::string sEDiff{ Args.Remove( Common::ModelBase::EDiffPrefix, Common::ModelBase::EDiffPrefix ) };
  if( !bESame )
  {
    E[idxSnk].Key.Object = { ObjectID( idxSnk ) };
    E[idxSnk].Key.Name = std::move( ESnk );
    // Different source and sink energies - save the difference
    EDiff.Key.Object = objectID;
    EDiff.Key.Name = sEDiff;
  }
  // Now create the matrix element
  MEL.Key.Object = objectID;
  static const std::string sMEL{ "MEL" };
  if( cp.pCorr->Name_.Gamma.size() > 1 )
    throw std::runtime_error( "More than one gamma unsupported" );
  Common::Gamma::Algebra MyGamma{ cp.pCorr->Name_.Gamma.empty() ? Common::Gamma::Algebra::Unknown
                                                                : cp.pCorr->Name_.Gamma[0] };
  MEL.Key.Name = Args.Remove( sMEL, sMEL );
  if( MyGamma != Common::Gamma::Algebra::Unknown )
    MEL.Key.Name.append( Common::Gamma::NameShort( MyGamma ) );
}

void Model3pt::AddParameters( Params &mp )
{
  for( ModelParam &p : E )
    AddEnergy( mp, p, NumOverlapExp, N, dispType );
  AddParam( mp, MEL, NumExponents, false );
  ModelOverlap::AddParameters( mp );
  if( !EDiff.Key.Object.empty() )
    AddParam( mp, EDiff, 1, false, Param::Type::Derived );
}

void Model3pt::SaveParameters( const Params &mp )
{
  for( ModelParam &p : E )
    p.idx = mp.at( p.Key )();
  MEL.idx = mp.at( MEL.Key )( 0, 0 );
  ModelOverlap::SaveParameters( mp );
  if( !EDiff.Key.Object.empty() )
    EDiff.idx = mp.at( EDiff.Key )();
}

// Get a descriptive string for the model
std::string Model3pt::Description() const
{
  std::string s{ "C3(" };
  s.append( std::to_string( NumExponents ) );
  s.append( "-exp,DeltaT=" );
  s.append( std::to_string( DeltaT ) );
  s.append( 1, ',' );
  s.append( MEL.Key.Name );
  if( objectID.size() > 1 )
  {
    s.append( 1, ',' );
    s.append( ObjectID( idxSnk ) );
  }
  s.append( 1, ',' );
  s.append( ObjectID( idxSrc ) );
  s.append( ModelOverlap::Description() );
  s.append( 1, ')' );
  return s;
}

// For now, all we can guess are matrix elements - and we assume they are 1
void Model3pt::Guessable( ParamsPairs &PP ) const
{
  // I can guess the matrix element
  PP.SetState( ParamsPairs::State::Known, MEL.Key, MEL.param->size );
  if( !EDiff.Key.Object.empty()
     && PP.keystate[E[idxSrc].idx] == ParamsPairs::State::Known
     && PP.keystate[E[idxSnk].idx] == ParamsPairs::State::Known )
  {
    // I know the energies, so I know their difference
    PP.SetState( ParamsPairs::State::Known, EDiff.Key, EDiff.param->size );
  }
  // TODO: Implement guessing
}

std::size_t Model3pt::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                             const VectorView &FitData, std::vector<int> FitTimes,
                             bool bLastChance ) const
{
  // I can guess the matrix element
  for( std::size_t i = 0; i < MEL.param->size; ++i )
  {
    if( !bKnown[MEL.idx + i] )
    {
      bKnown[MEL.idx + i] = true;
      Guess[MEL.idx + i] = 1;
    }
  }
  return 0; // TODO: Implement guessing
}

// TODO: Make sure we add partial derivatives for all models
double Model3pt::Derivative( int t, int p ) const
{
  return 0;
}

ModelType Model3pt::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::ThreePoint );
  return m;
}

scalar Model3pt::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  scalar z = 0;
  for( int eSnk = 0; eSnk < NumOverlapExp; ++eSnk )
  {
    const scalar &ESnk{ ModelParams[E.back().idx + eSnk] };
    const scalar ASnk = ModelParams[Overlap( eSnk, idxSnk )];
    const scalar expSnk{ - ESnk * ( DeltaT - t ) };
    for( int eSrc = 0; eSrc < NumOverlapExp; ++eSrc )
    {
      if( !eSnk || !eSrc || NumExponents >= 3 )
      {
        const scalar &ESrc{ ModelParams[E[idxSrc].idx + eSrc] };
        const scalar ASrc = ModelParams[Overlap( eSrc, idxSrc )];
        const scalar expSrc{ - ESrc * t };
        const scalar &thisMEL{ ModelParams[ ( *MEL.param )( eSnk, eSrc ) ] };
        const scalar OverProd{ ASnk * ASrc };
        scalar d{ OverProd * thisMEL * std::exp( expSnk + expSrc ) };
        if( !bOverlapAltNorm )
        {
          // If the factors of 1/2E have been absorbed into overlap coefficients, then no need to do this
          d /= 4 * ESnk * ESrc;
        }
        z += d;
        if( !eSnk && !eSrc )
        {
          if( !EDiff.Key.Object.empty() )
            ModelParams[EDiff.idx] = ESrc - ESnk;
        }
      }
    }
  }
  return z;
}

/*std::size_t Model3pt::ParamIndex( std::size_t idxSnk, std::size_t idxSrc ) const
{
  if( idxSnk > 1 || idxSrc > 1 )
    throw std::runtime_error( "Model3pt::ParamIndex index out of bounds" );
  std::size_t idx{ idxSnk + idxSrc };
  if( idxSnk && E.size() > 1 )
    idx++;
  return idx;
}

// TODO: Not sure I use this. Inconsistent definitions of how many parameters
std::size_t Model3pt::NumUnknown( std::vector<bool> &bKnown ) const
{
  std::size_t UnKnown{ ModelOverlap::NumUnknown( bKnown ) };
  for( const ModelParam &p : E )
  {
    for( std::size_t i = 0; i < p.param->size; ++i )
    {
      if( !bKnown[p.idx + i] )
        ++UnKnown;
    }
  }
  return UnKnown;
}*/
