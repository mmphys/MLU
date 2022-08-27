/*************************************************************************************
 
 Base for all models
 
 Source file: ModelBase.hpp
 
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

#ifndef ModelBase_hpp
#define ModelBase_hpp

#include "MultiFit.hpp"

// This is the name of an energy level
extern const std::string E;

enum class ModelType{ Unknown, Exp, Cosh, Sinh, ThreePoint, Constant };
std::ostream & operator<<( std::ostream &os, const ModelType m );
std::istream & operator>>( std::istream &is, ModelType &m );

// Default parameters for model creation
struct ModelDefaultParams
{
  const bool bForceSrcSnkDifferent;
  const int NumExponents;
  //int NumOps;
  ModelDefaultParams( bool bForceSrcSnkDifferent_, const Common::CommandLine &cl )
  : bForceSrcSnkDifferent{ bForceSrcSnkDifferent_ }, NumExponents{ cl.SwitchValue<int>( "e" ) } {}
};

// This represents the model I'm fitting to
class Model;
using ModelPtr = std::unique_ptr<Model>;

class Model
{
  friend class ModelSet;
public:
  int Nt;
  int NumExponents;
  vString ParamNames;
  vString ParamNamesPerExp;                     // Be careful: these names EXCLUDE the energy
  vInt    ParamIdx;
  std::vector<std::vector<int>> ParamIdxPerExp; // Be careful: Zero'th element of this IS energy
  //protected:
  Vector Params;
protected:
  inline void SetParams( int Nt_, int NumExponents_ ) { Nt=Nt_; NumExponents=NumExponents_; }
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName )
  { ParamNamesPerExp.push_back( E ); } // Make sure to call base if you override (because everything uses Energy)
public:
  virtual ~Model() {}
  // This is how to make a model
  static ModelPtr MakeModel( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  // NOT COUNTING ENERGY, put names of parameters that can be determined in list, return number remaining to determine
  // For now, synonymous with overlap coefficients
  virtual std::size_t UnknownParameterCount( UniqueNames &Names ) const { return 0; }
  // There are more unknowns than correlators - combine overlap coefficients into single operator
  virtual void ReduceUnknown( const UniqueNames &Names ) {}
  // Take a guess as to E0 (guesses will be averaged, and guesses for higher energies based on E0)
  virtual bool GuessE0( scalar &E0, const Vector &Corr ) const { return false; }
  // Take a guess for all parameters other than energies (guaranteed to be valid on entry). Return true for another pass
  virtual bool Guess( Vector &Parameters, std::vector<bool> &bKnown, int pass, const Vector &Corr ) const { return false; }
  // How big a scratchpad is required for each model?
  virtual int GetScratchPadSize() const { return 0; }
  // Cache values based solely on the model parameters (to speed up computation)
  virtual void ModelParamsChanged( Vector &ScratchPad, const Vector &ModelParams ) const {};
  // This is where the actual computation is performed
  virtual scalar operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const = 0;
  // Partial derivative of the model function on a timeslice with respect to each parameter
  virtual double Derivative( int t, int p ) const { return 0; }
};

#endif // ModelBase_hpp
