/**
 
 Use Minuit2 as fitting engine (OPTIONAL)
 
 Source file: FitMinuit2.cpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "FitMinuit2.hpp"

const ROOT::Minuit2::MnStrategy FitterThreadMinuit2::Strategy( FitterThreadMinuit2::StrategyLevel );

void FitterThreadMinuit2::DumpParamsFitter( std::ostream &os ) const
{
  if( state.bValid )
    os << Minuit2State;
  else
    parent.mp.Dump( os, ModelParams, Param::Type::Variable, state.bValid ? &state.ModelErrors : nullptr );
}

void FitterThreadMinuit2::ReplicaMessage( std::ostream &os ) const
{
  if( state.bValid )
    os << "edm " << Minuit2State.Edm() << ", ";
}

std::string FitterThreadMinuit2::DescriptionImpl() const
{
  return "Minuit2 fitter";
}

FitterThread * FitterThreadMinuit2::Clone() const
{
  return new FitterThreadMinuit2( *this );
}

// Painful, but Minuit2 calls operator() const
double FitterThreadMinuit2::operator()( const std::vector<double> & par ) const
{
  FitterThreadMinuit2 * pMe{ const_cast<FitterThreadMinuit2 *>( this ) };
  if( !pMe->SaveError( Error, &par[0], par.size() ) )
    return MLU::NaN;
  double chi2;
  if( bCorrelated )
    chi2 = Error.Dot( Error );
  else
  {
    chi2 = 0;
    for( int i = 0; i < Extent; ++i )
      chi2 += Error[i] * Error[i];
  }
  if( chi2 < 0 )
    return MLU::NaN;
  return chi2;
}

void FitterThreadMinuit2::Minimise( int iNumGuesses )
{
  ROOT::Minuit2::MnUserParameters Par;
  if( !state.bValid )
  {
    for( const Params::value_type &it : parent.mp )
    {
      const Param &p{ it.second };
      if( p.type == Param::Type::Variable )
      {
        for( std::size_t i = 0; i < p.size; ++i )
          Par.Add( it.first.FullName( i, p.size ), FitterParams[p.GetOffset(i,Param::Type::Variable)], 0 );
      }
    }
    state.NumCalls = 0;
  }
  ROOT::Minuit2::FunctionMinimum min = Minimiser.Minimize( *this, state.bValid ? Minuit2State : Par,
                                                           Strategy, parent.MaxIt, parent.Tolerance * 1000 );
  const ROOT::Minuit2::MnUserParameterState &M2State{ min.UserState() };
  if( !M2State.IsValid() )
    throw std::runtime_error( ReplicaString( iNumGuesses ) + " did not converge" );
  Minuit2State = M2State;
  state.TestStat = M2State.Fval();
  state.NumCalls += M2State.NFcn();
  state.bValid = true;
  const std::vector<double> M2Params{ M2State.Parameters().Params() };
  const std::vector<double> M2Errors{ M2State.Parameters().Errors() };
  if( M2Params.size() != FitterParams.size )
    throw std::runtime_error( "FitterThreadMinuit2::Minimise incorrect parameter size "
                             + std::to_string( M2Params.size() ) );
  if( M2Errors.size() != state.FitterErrors.size || M2Errors.size() != M2Params.size() )
    throw std::runtime_error( "FitterThreadMinuit2::Minimise incorrect error size "
                             + std::to_string( M2Errors.size() ) );
  // Copy the error the fitter computed for each parameter
  for( std::size_t i = 0; i < M2Errors.size(); ++i )
  {
    FitterParams[i] = M2Params[i];
    state.FitterErrors[i] = M2Errors[i];
  }
  // Save the result
  ( *this )( M2Params );
}

const std::string & FitterMinuit2::Type() const
{
  static const std::string MyType{ "Minuit2" };
  return MyType;
}

Fitter * MakeFitterMinuit2( const std::string &FitterArgs, Model::CreateParams &mcp,
                            DataSet &ds, std::vector<Model::Args> &&ModelArgs,
                            CovarParams &&cp, bool bFitCorr, FitController &fitController )
{
  if( !FitterArgs.empty() )
    return nullptr;
  return new FitterMinuit2(mcp, ds, std::move(ModelArgs), std::move(cp), bFitCorr, fitController);
}
