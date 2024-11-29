/**
 
 Abstract base for all models. i.e. interface consumed by fitter
 
 Source file: Model.hpp
 
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

#ifndef Model_hpp
#define Model_hpp

#include "MultiFit.hpp"

struct ModelType { int t; }; // An opaque type provided by Model Glue
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
struct Model;
using ModelPtr = std::unique_ptr<Model>;

struct Model
{
  struct Args : public std::map<std::string, std::string, MLU::LessCaseInsensitive>
  {
    void FromString( const std::string &s, bool bOptionalValue );
    std::string ToString() const;
    std::string Remove( const std::string &key, bool * Removed = nullptr );
    template <typename T> T Remove( const std::string &key, T Default, bool bPeek = false )
    {
      typename Model::Args::iterator it = find( key );
      if( it == end() )
        return Default;
      T t{ MLU::FromString<T>( it->second ) };
      if( !bPeek )
        erase( it );
      return t;
    }
  };

  // Default parameters for model creation - feel free to pass a derived type to models
  struct CreateParams
  {
    // These apply to all models
  protected:
    const std::vector<std::string> &OpNames;
    const std::vector<bool> bOpMomIndependent;
    std::vector<bool> CreateMomList( const std::string &MomIndependentOpList );
  public:
    const MLU::CommandLine &cl;
    const bool bOverlapAltNorm;
    // These will change per model
    const Sample *pCorr; // Actually either a Fold or a Model
    CreateParams( const std::vector<std::string> &OpNames_, const MLU::CommandLine &cl_ );
    virtual ~CreateParams(){}
    std::string GetOpName( int Idx ) const { return OpNames[pCorr->Name_.op[Idx]]; }
    bool OpMomIndependent( int Idx ) const { return bOpMomIndependent[pCorr->Name_.op[Idx]]; }
    virtual std::string Description() const;
  };
  std::vector<Params::value_type *> param;
  const int Nt;
  const int NumExponents;
protected:
  Model( const CreateParams &cp, int NumExponents_ )
  : Nt{cp.pCorr->NtUnfolded()}, NumExponents{NumExponents_} {}
  void AddParam( Params &mp, ModelParam &ModPar, std::size_t NumExp = 1, bool bMonotonic = false,
                 Param::Type Type = Param::Type::Variable );
  /**
   Add energy. If N != 0 then use the lattice dispersion relation for this parameter
   */
  void AddEnergy( Params &mp, ModelParam &ModPar, std::size_t NumExp, int N,
                  MLU::DispersionType dispType );
public:
  virtual ~Model() {}
  // This is how to make a model
  static ModelPtr MakeModel( int ModelNum, const CreateParams &mcp, Model::Args &Args );
  // These must be implemented by the model
  /// Models which are a function of a vector (eg continuum fit) must define X Vector, i.e. constants used
  virtual void DefineXVector( DataSet &ds, int i ) {}
  virtual const std::string &XVectorKeyName() const;
  virtual std::string XVectorKey() const;
  virtual std::vector<Param::Key> XVectorKeyNames() const;
  /// Models are expected to create parameters and constants they use
  virtual void AddParameters( Params &mp ) = 0;
  /// Models which fit to Models (not correlators) must return the column number they are fitting
  virtual std::size_t GetFitColumn() const;
  /// A complete model has been agreed, Models must save the parameters (ie offsets) they use
  virtual void SaveParameters( const Params &mp ) = 0;
  // Get a descriptive string for the model
  virtual std::string Description() const = 0;
  // Mark parameters guessable. Return number of unknown parameters remaining
  // NB: Ignore excited states - i.e. each parameter counts as 1 regardless of size
  virtual void Guessable( ParamsPairs &PP ) const = 0;
  // Guess unknown parameters. LastChance true -> fit about to be abandoned, so guess products of params
  // Return number of unknown parameters remaining. NB: Count excited states - i.e. return sum of param.size
  virtual std::size_t Guess( Vector &Guess, std::vector<bool> &bKnown, const Params &mp,
                             const VectorView &FitData, std::vector<int> FitTimes,
                             bool bLastChance ) const = 0;
  // Get model type
  virtual ModelType Type() const = 0;
  /**
   This is where the actual computation is performed
   
   - Warning: Don't change any ModelParams other than Derived
   */
  virtual scalar operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const = 0;
  // Partial derivative of the model function on a timeslice with respect to each parameter
  virtual double Derivative( int t, int p ) const { return 0; }
  // Model insoluble - please reduce number of parameters if possible?
  virtual void ReduceUnknown( const ParamsPairs &PP ) {} // Most models won't be able to do this
  // How big a scratchpad is required for each model?
  virtual int GetScratchPadSize() const { return 0; }
  /// This is called once at the start of each replica
  virtual void SetReplica( Vector &ScratchPad, Vector &ModelParams ) const {}
  /**
   Cache values based solely on the model parameters (to speed up computation)
   Called once per guess (not every timeslice)
   **/
  virtual void ModelParamsChanged( Vector &ScratchPad, const Vector &ModelParams ) const {}
};

// Models which describe a single object. Use with multiple inheritance
struct Object
{
  const std::string &ObjectID( int SnkSrc ) const
  { return objectID[SnkSrc==idxSrc || objectID.size() < 2 ? 0 : 1]; }
protected:
  Object( std::vector<std::string> &&ObjectID ) : objectID{ObjectID} {}
  std::vector<std::string> objectID; // Name of the object propagating
  // Use this to initialise where there's only one objectID, e.g. two-point functions
  static std::vector<std::string> GetObjectNameSingle( const Model::CreateParams &cp, Model::Args &Args );
  // Use this to initialise where there's a sink and source objectID, e.g. three-point functions
  static std::vector<std::string> GetObjectNameSnkSrc( const Model::CreateParams &cp, Model::Args &Args );
};

#endif // Model_hpp
