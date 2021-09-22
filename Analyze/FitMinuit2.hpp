/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: FitMinuit2.hpp
 
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

#ifndef FitMinuit2_hpp
#define FitMinuit2_hpp

#include <Minuit2/FCNBase.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/VariableMetricMinimizer.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MnPrint.h>

#include "MultiFit.hpp"

struct ParamStateMinuit2: public ParamState
{
  scalar edm; // Protected by bValid
  bool bGotMinuit2State;
  ROOT::Minuit2::MnUserParameterState Minuit2State;
protected:
  scalar Edm() const { return bValid ? edm : 0; }
        void CopyParameter( const ROOT::Minuit2::MnUserParameterState &State, Parameters & Par );
  Parameters MakeParameter( const ROOT::Minuit2::MnUserParameterState &State );
public:
  // Initialising constructors
  ParamStateMinuit2( const Parameters &par, bool bValid_ = false,
                     scalar TestStat_ = 0, unsigned int NumCalls_ = 0, scalar edm_ = 0 )
  : ParamState( par, bValid_, TestStat_, NumCalls_ ), edm{edm_}, bGotMinuit2State{false}
  {
  }
  ParamStateMinuit2( Parameters &&par, bool bValid_ = false,
                     scalar TestStat_ = 0, unsigned int NumCalls_ = 0, scalar edm_ = 0 )
  : ParamState( std::move(par), bValid_, TestStat_, NumCalls_ ), edm{edm_}, bGotMinuit2State{false}
  {
  }
  ParamStateMinuit2( const ROOT::Minuit2::MnUserParameterState &State )
  : ParamState( MakeParameter(State), State.IsValid(), State.IsValid() ? State.Fval() : 0,
                State.IsValid() ? State.NFcn() : 0 ),
    edm{State.IsValid() ? State.Edm() : 0}, bGotMinuit2State{true}, Minuit2State{State}
  {
  }
  ParamStateMinuit2 &operator=( const ParamStateMinuit2 &rhs )
  {
    ParamState::operator=( rhs );
    edm = rhs.edm;
    bGotMinuit2State = rhs.bGotMinuit2State;
    if( bGotMinuit2State )
      Minuit2State = rhs.Minuit2State;
    return *this;
  }
  void StandardOut( std::ostream &os ) const override; //TODO: Deuglify
  void ReplicaMessage( std::ostream &os ) const override; //TODO: Deuglify
};

// Several of these will be running at the same time on different threads during a fit
class FitterThreadMinuit2 : public FitterThread, ROOT::Minuit2::FCNBase
{
protected:
  static constexpr unsigned int StrategyLevel{ 1 }; // for parameter ERRORS (see MnStrategy) 0=low, 1=medium, 2=high
  static const ROOT::Minuit2::MnStrategy Strategy;
  ROOT::Minuit2::VariableMetricMinimizer Minimiser;
  // Helper functions
  ParamState * MakeParamState( const Parameters &Params ) override;
public:
  FitterThreadMinuit2( const Fitter &fitter_, bool bCorrelated_, ModelFile &OutputModel, vCorrelator &CorrSynthetic_ )
  : FitterThread( fitter_, bCorrelated_, OutputModel, CorrSynthetic_ ) {}
  virtual ~FitterThreadMinuit2() {}
  // These are part of the FCNBase interface
  double Up() const override { return 1.; }
  double operator()( const std::vector<double> &ModelParameters ) const override;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  void Minimise( ParamState &Guess, int iNumGuesses ) override;
  int NumRetriesGuess() const override { return parent.Retry ? parent.Retry + 10 : 20; };
  int NumRetriesFit() const override { return parent.Retry ? parent.Retry : 10; };
};

class FitterMinuit2 : public Fitter
{
  using Fitter::Fitter; // Inherit constructor (with same visibility)
protected:
  const std::string &Type() const override;
  FitterThread * MakeThread( bool bCorrelated, ModelFile &OutputModel, vCorrelator &CorrSynthetic ) override
  {
    return new FitterThreadMinuit2( *this, bCorrelated, OutputModel, CorrSynthetic );
  }
};

Fitter * MakeFitterMinuit2(const Common::CommandLine &cl, const DataSet &ds,
                           const std::vector<std::string> &ModelArgs, const std::vector<std::string> &opNames);

#endif // FitMinuit2_hpp
