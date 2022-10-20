/*************************************************************************************
 
 Abstract base for all models. i.e. interface consumed by fitter
 
 Source file: Model.hpp
 
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

#ifndef Model_hpp
#define Model_hpp

#include "MultiFit.hpp"

// This is the name of an energy level
extern const std::string E;

enum class ModelType{ Unknown, Exp, Cosh, Sinh, ThreePoint, Constant, R3 };
std::ostream & operator<<( std::ostream &os, const ModelType m );
std::istream & operator>>( std::istream &is, ModelType &m );

// Models use this structure to refer to the parameters they are interested in
struct ModelParam
{
  ModelParam() = default;
  Param::Key Key;
  Param *param;
  std::size_t idx;
};

// This represents the model I'm fitting to
class Model;
using ModelPtr = std::unique_ptr<Model>;

struct Model
{
  struct Args : public std::map<std::string, std::string, Common::LessCaseInsensitive>
  {
    void FromString( const std::string &s, bool bOptionalValue );
    std::string Remove( const std::string &key, bool * Removed = nullptr );
    template <typename T> T Remove( const std::string &key, T Default, bool bPeek = false )
    {
      typename Model::Args::iterator it = find( key );
      if( it == end() )
        return Default;
      T t{ Common::FromString<T>( it->second ) };
      if( !bPeek )
        erase( it );
      return t;
    }
  };

  // Default parameters for model creation - these apply to all models
  struct CreateParams
  {
    // These apply to all models
    const std::vector<std::string> &OpNames;
    const Common::CommandLine &cl;
    const int NumExponents; // Default if not specified per model
    const bool bOverlapAltNorm;
    // These will change per model
    const Fold *pCorr;
    CreateParams( const std::vector<std::string> &OpNames_, const Common::CommandLine &cl_ );
  };
  std::vector<Params::value_type *> param;
  const int Nt;
  const int NumExponents;
protected:
  Model( const CreateParams &cp, Args &args );
  void AddParam( Params &mp, ModelParam &ModPar, std::size_t NumExp = 1, bool bMonotonic = false,
                 Param::Type Type = Param::Type::Variable );
public:
  virtual ~Model() {}
  // This is how to make a model
  static ModelPtr MakeModel( const CreateParams &cp, Model::Args &Args );
  // These must be implemented by the model
  // Say which parameters this model wants
  virtual void AddParameters( Params &mp ) = 0;
  // Now that a complete model has been agreed, save the parameters this model will use
  virtual void SaveParameters( const Params &mp ) = 0;
  // Get a descriptive string for the model
  virtual std::string Description() const = 0;
  // Mark parameters guessable. Return number of unknown parameters remaining
  // NB: Ignore excited states - i.e. each parameter counts as 1 regardless of size
  virtual std::size_t Guessable( std::vector<bool> &bKnown, bool bLastChance ) const = 0;
  // Guess unknown parameters. LastChance true -> fit about to be abandoned, so guess products of params
  // Return number of unknown parameters remaining. NB: Count excited states - i.e. return sum of param.size
  virtual std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown,
                       const VectorView &FitData, std::vector<int> FitTimes, bool bLastChance ) const = 0;
  // Get model type
  virtual ModelType Type() const = 0;
  // This is where the actual computation is performed
  virtual scalar operator()( int t, Vector &ScratchPad, const Vector &ModelParams ) const = 0;
  // Partial derivative of the model function on a timeslice with respect to each parameter
  virtual double Derivative( int t, int p ) const { return 0; }
  // Model insoluble - please reduce number of parameters if possible?
  virtual void ReduceUnknown() {} // Most models won't be able to do this
  // How big a scratchpad is required for each model?
  virtual int GetScratchPadSize() const { return 0; }
  // Cache values based solely on the model parameters (to speed up computation)
  virtual void ModelParamsChanged( Vector &ScratchPad, const Vector &ModelParams ) const {};
};

#endif // Model_hpp
