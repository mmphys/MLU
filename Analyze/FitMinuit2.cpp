/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: FitMinuit2.cpp
 
 Copyright (C) 2021
 
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

#include "FitMinuit2.hpp"

// Fill a Parameters structure from Minuit2 State
void ParamStateMinuit2::CopyParameter( const ROOT::Minuit2::MnUserParameterState &State, Parameters & Par )
{
  // if( State.IsValid() ) // Do I need this
  Par.clear();
  for( const ROOT::Minuit2::MinuitParameter &mp : State.Parameters().Parameters() )
    Par.Add( mp.GetName(), mp.Value(), mp.Error() );
}

// Make a Parameters structure from Minuit2 State
Parameters ParamStateMinuit2::MakeParameter( const ROOT::Minuit2::MnUserParameterState &State )
{
  Parameters Par;
  CopyParameter( State, Par );
  return Par;
}

void ParamStateMinuit2::StandardOut( std::ostream &os ) const
{
  if( bGotMinuit2State )
    os << Minuit2State;
}

void ParamStateMinuit2::ReplicaMessage( std::ostream &os ) const
{
  if( bGotMinuit2State )
    os << "edm " << Edm() << ", ";
}

ParamState * FitterThreadMinuit2::MakeParamState( const Parameters &Params )
{
  return new ParamStateMinuit2( Params );
}

const ROOT::Minuit2::MnStrategy FitterThreadMinuit2::Strategy( FitterThreadMinuit2::StrategyLevel );

// Painful, but Minuit2 calls operator() const
double FitterThreadMinuit2::operator()( const std::vector<double> & par ) const
{
  scalar localBuffer[Extent];
  Vector Error;
  Error.MapView( localBuffer, Extent );
  FitterThreadMinuit2 * pMe{ const_cast<FitterThreadMinuit2 *>( this ) };
  if( !pMe->SaveError( Error, &par[0], par.size() ) )
    return std::numeric_limits<double>::max();
  double chi2;
  if( bCorrelated )
    chi2 = Error.Dot( Cholesky.CholeskySolve( Error ) );
  else
  {
    chi2 = 0;
    for( int i = 0; i < Extent; ++i )
      chi2 += Error[i] * Error[i];
  }
  if( chi2 < 0 )
    throw std::runtime_error( "Chi^2 < 0 on replica " + std::to_string( idx ) );
  return chi2;
}

void FitterThreadMinuit2::Minimise( ParamState &Guess_, int iNumGuesses )
{
  ParamStateMinuit2 &Guess{ *dynamic_cast<ParamStateMinuit2*>( &Guess_ ) };
  ROOT::Minuit2::MnUserParameters Minuit2Par;
  if( !Guess.bGotMinuit2State )
    for( const Parameters::Parameter & p : Guess.parameters )
      Minuit2Par.Add( p.Name, p.Value, p.Error );
  ROOT::Minuit2::FunctionMinimum min = Minimiser.Minimize( *this, Guess.bGotMinuit2State ? Guess.Minuit2State : Minuit2Par,
                                                           Strategy, parent.MaxIt, parent.Tolerance * 1000 );
  const ROOT::Minuit2::MnUserParameterState &state{ min.UserState() };
  if( !state.IsValid() )
    throw std::runtime_error( ReplicaString( iNumGuesses ) + " did not converge" );
  Guess = state;
}

const std::string & FitterMinuit2::Type() const
{
  static const std::string MyType{ "Minuit2" };
  return MyType;
}

Fitter * MakeFitterMinuit2(const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                           const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames)
{
  if( !FitterArgs.empty() )
    return nullptr;
  return new FitterMinuit2( cl, ds, ModelArgs, opNames );
}
