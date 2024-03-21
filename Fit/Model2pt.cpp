/*************************************************************************************
 
 Models for 2pt functions

 Source file: Model2pt.cpp
 
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

#include "Model2pt.hpp"

Model2pt::Model2pt( const Model::CreateParams &cp, Model::Args &Args, int NumExp,
                    std::vector<std::string> &&objectID, std::vector<std::string> &&opNames )
: ModelOverlap( cp, Args, NumExp, static_cast<std::size_t>( NumExp ),
                std::move( objectID ), std::move( opNames ) ),
  N{ dynamic_cast<const MultiFitCreateParams &>( cp ).N },
dispType{ dynamic_cast<const MultiFitCreateParams &>( cp ).dispType }
{
  E.Key.Object = { ObjectID( idxSrc ) };
  E.Key.Name = Args.Remove( "Energy", MLU::ModelBase::EnergyPrefix );
}

void Model2pt::AddParameters( Params &mp )
{
  AddEnergy( mp, E, NumOverlapExp, N, dispType );
  ModelOverlap::AddParameters( mp );
}

void Model2pt::SaveParameters( const Params &mp )
{
  E.idx = mp.at( E.Key )();
  ModelOverlap::SaveParameters( mp );
}

// Get a descriptive string for the model
std::string Model2pt::Description() const
{
  std::string s( 1, '(' );
  s.append( std::to_string( NumExponents ) );
  s.append( "-exp," );
  s.append( ObjectID( idxSrc ) );
  s.append( ModelOverlap::Description() );
  s.append( 1, ')' );
  return s;
}

void Model2pt::Guessable( ParamsPairs &PP ) const
{
  // I can guess the energy
  PP.SetState( ParamsPairs::State::Known, E.Key, NumOverlapExp );
  ModelOverlap::Guessable( PP );
}

// Take a guess as to my parameters
std::size_t Model2pt::Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                             const VectorView &FitData, std::vector<int> FitTimes,
                             bool bLastChance ) const
{
  // I need a minimum of two timeslices per energy level to guess
  static constexpr std::size_t DataPointsPerGuess{ 2 };
  const std::size_t NumICanGuess{ FitData.size() / DataPointsPerGuess };
  bool bOK{ NumICanGuess != 0 };
  const double Stride{ NumICanGuess < NumOverlapExp ? DataPointsPerGuess
                     : NumOverlapExp <= 1 ? 0
                     : static_cast<double>( FitData.size() - DataPointsPerGuess ) / ( NumOverlapExp - 1 ) };
  for( std::size_t i = 0; bOK && i < NumOverlapExp; ++i )
  {
    const std::size_t tGuess{ NumOverlapExp <= 1 ? FitData.size() - DataPointsPerGuess
                              : static_cast<std::size_t>( Stride * ( NumOverlapExp - 1 - i ) + 0.5 ) };
    // Guess the energy
    if( !bKnown[E.idx + i] )
    {
      scalar ThisEnergy;
      if( i < NumICanGuess )
      {
        // Assume each new pair of data points can explain one more excited-state energy
        const scalar C1{ Estimate( Guess, FitData, FitTimes, i, tGuess ) };
        const scalar C2{ Estimate( Guess, FitData, FitTimes, i, tGuess + 1 ) };
        const scalar Log{ std::log( std::abs( C1 / C2 ) ) };
        const int DeltaT{ FitTimes[tGuess + 1] - FitTimes[tGuess] };
        ThisEnergy = Log / DeltaT;
      }
      else if( i == 0 )
      {
        // I just don't have enough data points to fit anything
        break;
      }
      else
      {
        // Very crude guess, but assume each energy level is half the previous Delta E
        ThisEnergy = Guess[E.idx + i - 1];
        if( i == 1 )
          ThisEnergy *= 1.5;
        else
          ThisEnergy += 0.5 * ( Guess[E.idx + i - 1] - Guess[E.idx + i - 2] );
      }
      mp.GuessEnergy( Guess, bKnown, E.Key, i, ThisEnergy );
    }
    // Do we need to guess overlap coefficients?
    const std::size_t OverlapSize{ OverlapCount( i ) };
    const bool bKnow0{ bKnown[Overlap( i, idxSrc )] };
    const bool bKnow1{ bKnown[Overlap( i, idxSnk )] };
    if( !bKnow0 || !bKnow1 )
    {
      const bool Make1From0{  bKnow0 && !bKnow1 };
      const bool Make0From1{ !bKnow0 &&  bKnow1 };
      if( ( OverlapSize == 1 && !bKnow0 )
       || ( OverlapSize == 2 && ( bLastChance || Make1From0 || Make0From1 || ( i && !bKnow0 && !bKnow1 ) ) ) )
      {
        if( i < NumICanGuess )
        {
          const scalar Residual{ Estimate( Guess, FitData, FitTimes, i, tGuess ) };
          const scalar Product{ Residual * std::exp( Guess[E.idx + i] * FitTimes[tGuess] ) };
          if( OverlapSize == 1 )
          {
            if( Product > 0 )
              Guess[Overlap( i, idxSrc )] = std::sqrt( Product );
            else
            {
              scalar Factor;
              if( i > 1 )
                Factor = Guess[Overlap( i - 1, idxSrc )] / Guess[Overlap( i - 2, idxSrc )];
              else
                Factor = 0.5;
              Guess[Overlap( i, idxSrc )] = Guess[Overlap( i - 1, idxSrc )] * Factor;
            }
          }
          else if( Make0From1 )
            Guess[Overlap( i, idxSrc )] = Product / Guess[Overlap( i, idxSnk )];
          else if( Make1From0 )
            Guess[Overlap( i, idxSnk )] = Product / Guess[Overlap( i, idxSrc )];
          else //if( bLastChance )
          {
            Guess[Overlap( i, idxSrc )] = std::sqrt( std::abs( Product ) );
            Guess[Overlap( i, idxSnk )] = Guess[Overlap( i, idxSrc )];
            if( Product < 0 )
              Guess[Overlap( i, idxSnk )] = -Guess[Overlap( i, idxSnk )];
          }
        }
        else if( i == 0 )
        {
          // I just don't have enough data points to fit anything
          break;
        }
        else
        {
          // Very crude guess, but assume each energy level is half the previous Delta E
          Guess[Overlap( i, idxSrc )] = std::sqrt( 2 ) * Guess[Overlap( i - 1, idxSrc )];
          if( OverlapSize == 2 )
            Guess[Overlap( i, idxSnk )] = std::sqrt( 2 ) * Guess[Overlap( i, idxSnk - 1 )];
          else if( i && OverlapCount( i - 1 ) == 2 ) // Combined, but previous level wasn't
            Guess[Overlap( i, idxSrc )] *= std::sqrt( 2 ) * Guess[Overlap( i - 1, idxSnk )];
        }
        if( !bKnow0 )
          bKnown[Overlap( i, idxSrc )] = true;
        if( OverlapSize > 1 && !bKnow1 )
          bKnown[Overlap( i, idxSnk )] = true;
      }
      else
      {
        // I just don't have enough data points to fit anything
        break;
      }
    }
  }
  std::size_t NumUnknown{ 0 };
  for( std::size_t i = 0; i < NumOverlapExp; ++i )
  {
    if( !bKnown[E.idx + i] )
      NumUnknown++;
    for( int j = 0; j < OverlapCount( i ); ++j )
      if( !bKnown[Overlap( i, j )] )
        NumUnknown++;
  }
  return NumUnknown;
}

// TODO: Make sure we add partial derivatives for all models
double Model2pt::Derivative( int t, int p ) const
{
  double d = 0;
  // TODO re-implement
  /*for( int e = 0; e < NumExponents; ++e )
  {
    if( Energy[e] == p )
    {
      // Derivative wrt Energy
      d = t * std::exp( - Params[Energy[e]] * t );
      if( parity == MLU::Parity::Even || parity == MLU::Parity::Odd )
      {
        int NtMinusT = parity == MLU::Parity::Even ? Nt - t : t - Nt;
        d += NtMinusT * std::exp( - Params[Energy[e]] * ( Nt - t ) );
        d *= 0.5;
      }
      d *= - Params[src[e]] * Params[snk[e]];
      break;
    }
    // Perhaps this is a derivative wrt overlap coefficient
    // Work out whether this operator is source or sink or both
    double OtherOverlap = 0;
    int Factor{ 0 };
    if( snk[e] == p )
    {
      OtherOverlap = Params[src[e]];
      ++Factor;
    }
    if( src[e] == p )
    {
      OtherOverlap = Params[snk[e]]; // Doesn't matter if we do this twice ... it's the same operator!
      ++Factor;
    }
    if( Factor )
    {
      switch( parity )
      {
        case MLU::Parity::Even:
          d = SinhCoshAdjust[e] * std::cosh( - Params[Energy[e]] * ( t - HalfNt ) );
          break;
        case MLU::Parity::Odd:
          d = SinhCoshAdjust[e] * std::sinh( - Params[Energy[e]] * ( t - HalfNt ) );
          break;
        default:
          d = std::exp( - Params[Energy[e]] * t );
          break;
      }
      d *= Factor * OtherOverlap;
      break;
    }
  }*/
  return d;
}

// This is used by the guesser to take an estimate of the excited state
scalar Model2pt::Estimate( Vector &Guess, const VectorView &FitData, std::vector<int> FitTimes,
                           std::size_t NumExp, std::size_t Timeslice ) const
{
  scalar Theory{ 0 };
  for( std::size_t i = 0; i < NumExp; ++i )
  {
    Theory += Guess[Overlap( i, idxSrc )] * Guess[Overlap( i, idxSnk )]
              * std::exp( - Guess[E.idx + i] * FitTimes[Timeslice] );
  }
  return FitData[Timeslice] - Theory;
}

ModelType ModelExp::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::Exp );
  return m;
}

scalar ModelExp::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  double z = 0;
  for( int e = 0; e < NumOverlapExp; ++e )
  {
    double d = std::exp( - ModelParams[E.idx + e] * t );
    d *= ModelParams[Overlap( e, idxSrc )] * ModelParams[Overlap( e, idxSnk )];
    if( !bOverlapAltNorm )
      d /= 2 * ModelParams[E.idx + e];
    z += d;
  }
  return z;
}

std::string ModelExp::Description() const
{
  std::string s{ "C2Exp" };
  s.append( Model2pt::Description() );
  return s;
}

ModelType ModelCosh::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::Cosh );
  return m;
}

scalar ModelCosh::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  double z = 0;
  for( int e = 0; e < NumOverlapExp; ++e )
  {
    double d = std::exp( - ModelParams[E.idx + e] * t );
    d += std::exp( - ModelParams[E.idx + e] * ( Nt - t ) );
    d *= ModelParams[Overlap( e, idxSrc )] * ModelParams[Overlap( e, idxSnk )];
    if( !bOverlapAltNorm )
      d /= 2 * ModelParams[E.idx + e];
    z += d;
  }
  return z;
}

std::string ModelCosh::Description() const
{
  std::string s{ "C2Cosh" };
  s.append( Model2pt::Description() );
  return s;
}

ModelType ModelSinh::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::Sinh );
  return m;
}

scalar ModelSinh::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  double z = 0;
  for( int e = 0; e < NumOverlapExp; ++e )
  {
    double d = std::exp( - ModelParams[E.idx + e] * t );
    d -= std::exp( - ModelParams[E.idx + e] * ( Nt - t ) );
    d *= ModelParams[Overlap( e, idxSrc )] * ModelParams[Overlap( e, idxSnk )];
    if( !bOverlapAltNorm )
      d /= 2 * ModelParams[E.idx + e];
    z += d;
  }
  return z;
}

std::string ModelSinh::Description() const
{
  std::string s{ "C2Sinh" };
  s.append( Model2pt::Description() );
  return s;
}
