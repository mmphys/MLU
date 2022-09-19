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

Model2pt::Model2pt( const Model::CreateParams &cp, Model::Args &Args )
: ModelOverlap( cp, Args, GetObjectNameSingle( cp, Args ), Args.Remove( "e", cp.NumExponents, true ) )
{
  E.Key.Object = { ObjectID( idxSrc ) };
  E.Key.Name = Args.Remove( "energy", ::E );
}

void Model2pt::AddParameters( Params &mp )
{
  AddParam( mp, E, NumOverlapExp, true );
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
  std::string s{ "C2(" };
  s.append( ObjectID( idxSrc ) );
  s.append( ModelOverlap::Description() );
  s.append( 1, ')' );
  return s;
}

std::size_t Model2pt::Guessable( std::vector<bool> &bKnown, bool bLastChance ) const
{
  // I can guess the energy
  for( std::size_t i = 0; i < E.param->size; ++i )
    bKnown[E.idx + i] = true;
  return ModelOverlap::Guessable( bKnown, bLastChance );
}

// Take a guess as to my parameters
std::size_t Model2pt::Guess( Vector &Guess, std::vector<bool> &bKnown,
                       const VectorView &FitData, std::vector<int> FitTimes, bool bLastChance ) const
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
      if( i < NumICanGuess )
      {
        // Assume each new pair of data points can explain one more excited-state energy
        const scalar C1{ Estimate( Guess, FitData, FitTimes, i, tGuess ) };
        const scalar C2{ Estimate( Guess, FitData, FitTimes, i, tGuess + 1 ) };
        const scalar Log{ std::log( std::abs( C1 / C2 ) ) };
        const int DeltaT{ FitTimes[tGuess + 1] - FitTimes[tGuess] };
        Guess[E.idx + i] = Log / DeltaT;
        // Enforce monotonic
        if( i > 0 && Guess[E.idx + i] < Guess[E.idx + i - 1] )
          Guess[E.idx + i] = Guess[E.idx + i - 1];
      }
      else if( i == 0 )
      {
        // I just don't have enough data points to fit anything
        bOK = false;
        break;
      }
      else
      {
        // Very crude guess, but assume each energy level is half the previous Delta E
        if( i == 1 )
          Guess[E.idx + i] = Guess[E.idx + i - 1] * 1.5;
        else
          Guess[E.idx + i] = Guess[E.idx + i - 1] + 0.5 * ( Guess[E.idx + i - 1] - Guess[E.idx + i - 2] );
      }
      bKnown[E.idx + i] = true;
    }
    // Do we need to guess overlap coefficients?
    const bool bKnow0{ bKnown[Overlap[0].idx + i] };
    const bool bKnow1{ bKnown[Overlap.back().idx + i] };
    if( !bKnow0 || !bKnow1 )
    {
      const bool Make1From0{  bKnow0 && !bKnow1 };
      const bool Make0From1{ !bKnow0 &&  bKnow1 };
      if( ( Overlap.size() == 1 && !bKnow0 )
       || ( Overlap.size() == 2 && ( bLastChance || Make1From0 || Make0From1 || ( i && !bKnow0 && !bKnow1 ) ) ) )
      {
        if( i < NumICanGuess )
        {
          const scalar Residual{ Estimate( Guess, FitData, FitTimes, i, tGuess ) };
          const scalar Product{ Residual * std::exp( Guess[E.idx + i] * FitTimes[tGuess] ) };
          if( Overlap.size() == 1 )
          {
            if( Product > 0 )
              Guess[Overlap[0].idx + i] = std::sqrt( Product );
            else
            {
              scalar Factor;
              if( i > 1 )
                Factor = Guess[Overlap[0].idx + i - 1] / Guess[Overlap[0].idx + i - 2];
              else
                Factor = 0.5;
              Guess[Overlap[0].idx + i] = Guess[Overlap[0].idx + i - 1] * Factor;
            }
          }
          else if( Make0From1 )
            Guess[Overlap[0].idx + i] = Product / Guess[Overlap[1].idx + i];
          else if( Make1From0 )
            Guess[Overlap[1].idx + i] = Product / Guess[Overlap[0].idx + i];
          else //if( bLastChance )
          {
            Guess[Overlap[0].idx + i] = std::sqrt( std::abs( Product ) );
            Guess[Overlap[1].idx + i] = Guess[Overlap[0].idx + i];
            if( Product < 0 )
              Guess[Overlap[1].idx + i] = -Guess[Overlap[1].idx + i];
          }
        }
        else if( i == 0 )
        {
          // I just don't have enough data points to fit anything
          bOK = false;
          break;
        }
        else
        {
          // Very crude guess, but assume each energy level is half the previous Delta E
          Guess[Overlap[0].idx + i] = std::sqrt( 2 ) * Guess[Overlap[0].idx + i - 1];
          if( Overlap.size() == 2 )
            Guess[Overlap[1].idx + i] = std::sqrt( 2 ) * Guess[Overlap[1].idx + i - 1];
        }
        if( !bKnow0 )
          bKnown[Overlap[0].idx + i] = true;
        if( Overlap.size() > 1 && !bKnow1 )
          bKnown[Overlap[1].idx + i] = true;
      }
      else
      {
        // I just don't have enough data points to fit anything
        bOK = false;
        break;
      }
    }
  }
  std::size_t NumUnknown{ 0 };
  for( std::size_t i = 0; i < NumOverlapExp; ++i )
  {
    if( !bKnown[E.idx + i] )
      NumUnknown++;
    for( std::size_t j = 0; j < Overlap.size(); ++j )
      if( !bKnown[Overlap[j].idx + i] )
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
      if( parity == Common::Parity::Even || parity == Common::Parity::Odd )
      {
        int NtMinusT = parity == Common::Parity::Even ? Nt - t : t - Nt;
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
        case Common::Parity::Even:
          d = SinhCoshAdjust[e] * std::cosh( - Params[Energy[e]] * ( t - HalfNt ) );
          break;
        case Common::Parity::Odd:
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
    Theory += Guess[Overlap[0].idx + i] * Guess[Overlap.back().idx + i]
              * std::exp( - Guess[E.idx + i] * FitTimes[Timeslice] );
  }
  return FitData[Timeslice] - Theory;
}

scalar ModelExp::operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const
{
  double z = 0;
  for( int e = 0; e < NumOverlapExp; ++e )
  {
    double d = std::exp( - ModelParams[E.idx + e] * t );
    d *= ModelParams[Overlap[0].idx + e] * ModelParams[Overlap.back().idx + e];
    if( bNormalisationByEnergy )
      d /= 2 * ModelParams[E.idx + e];
    z += d;
  }
  return z;
}

scalar ModelCosh::operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const
{
  double z = 0;
  for( int e = 0; e < NumOverlapExp; ++e )
  {
    double d = std::exp( - ModelParams[E.idx + e] * t );
    d += std::exp( - ModelParams[E.idx + e] * ( Nt - t ) );
    d *= ModelParams[Overlap[0].idx + e] * ModelParams[Overlap.back().idx + e];
    if( bNormalisationByEnergy )
      d /= 2 * ModelParams[E.idx + e];
    z += d;
  }
  return z;
}

scalar ModelSinh::operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const
{
  double z = 0;
  for( int e = 0; e < NumOverlapExp; ++e )
  {
    double d = std::exp( - ModelParams[E.idx + e] * t );
    d -= std::exp( - ModelParams[E.idx + e] * ( Nt - t ) );
    d *= ModelParams[Overlap[0].idx + e] * ModelParams[Overlap.back().idx + e];
    if( bNormalisationByEnergy )
      d /= 2 * ModelParams[E.idx + e];
    z += d;
  }
  return z;
}