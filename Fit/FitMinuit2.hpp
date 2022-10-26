/*************************************************************************************
 
 Use Minuit2 as fitting exgine (OPTIONAL)

 Source file: FitMinuit2.hpp
 
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

#ifndef FitMinuit2_hpp
#define FitMinuit2_hpp

#include "Fitter.hpp"
#include "FitterThread.hpp"
#include <Minuit2/FCNBase.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <Minuit2/VariableMetricMinimizer.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MinimumParameters.h>
#include <Minuit2/MinimumState.h>
#include <Minuit2/MnPrint.h>

// Several of these will be running at the same time on different threads during a fit
class FitterThreadMinuit2 : public FitterThread, ROOT::Minuit2::FCNBase
{
protected:
  static constexpr unsigned int StrategyLevel{ 1 }; // for parameter ERRORS (see MnStrategy) 0=low, 1=medium, 2=high
  static const ROOT::Minuit2::MnStrategy Strategy;
  ROOT::Minuit2::VariableMetricMinimizer Minimiser;
  // Fitter state
  ROOT::Minuit2::MnUserParameterState Minuit2State;
  void DumpParamsFitter( std::ostream &os ) const override; //TODO: Deuglify
  void ReplicaMessage( std::ostream &os ) const override; //TODO: Deuglify
  std::string DescriptionImpl() const override;
public:
  FitterThreadMinuit2( const Fitter &fitter_, bool bCorrelated_, ModelFile &OutputModel )
  : FitterThread( fitter_, bCorrelated_, OutputModel ) {}
  virtual ~FitterThreadMinuit2() {}
  FitterThread * Clone() const override;
  // These are part of the FCNBase interface
  double Up() const override { return 1.; }
  double operator()( const std::vector<double> &ModelParameters ) const override;
  //virtual void SetErrorDef(double def) {theErrorDef = def;}
  void Minimise( int iNumGuesses ) override;
  int NumRetriesGuess() const override { return parent.Retry ? parent.Retry + 10 : 20; };
  int NumRetriesFit() const override { return parent.Retry ? parent.Retry : 10; };
};

class FitterMinuit2 : public Fitter
{
  using Fitter::Fitter; // Inherit constructor (with same visibility)
protected:
  const std::string &Type() const override;
  FitterThread * MakeThread( bool bCorrelated, ModelFile &OutputModel ) override
  {
    return new FitterThreadMinuit2( *this, bCorrelated, OutputModel );
  }
};

#endif // FitMinuit2_hpp
