/*************************************************************************************
 
 Common components for model implementations
 
 Source file: ModelCommon.hpp
 
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

#ifndef ModelCommon_hpp
#define ModelCommon_hpp

#include "Model.hpp"

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

/**
 A model with two overlap coefficients

 Overlap factor normalisation
 
 Call the overlap factors A = <n| O^\dag | \Omega>
 we have two choices of normalisation.
 
 - Parameters:
  - `bOverlapAltNorm = false` : **Default Normalisation**
 
      C2(t) = \sum_{j=1}^n \frac{A}{2 E_j} \exp{- E_j t}

 - `bOverlapAltNorm = true` : Alternate Normalisation
 
      C2(t) = \sum_{j=1}^n        A        \exp{- E_j t}
 
  - Warning: Do not mix and match!
 The same normalisation must be used everywhere.
 I.e. you **cannot** perform a two point fit with one normalisation...
 ...then feed that into a second fit using a different normalisation.
 
 - Warning: **Do not use the alternate normalisation unless you know what you are doing**
 */
struct ModelOverlap : Model, Object
{
  ModelOverlap( const Model::CreateParams &cp, Model::Args &Args,
                std::vector<std::string> &&ObjectID, std::size_t NumOverlapExp );
  void AddParameters( Params &mp ) override;
  void SaveParameters( const Params &mp ) override;
  std::string Description() const override;
  void Guessable( ParamsPairs &PP ) const override;
  void ReduceUnknown( const ParamsPairs &PP ) override;
// TODO: These really shouldn't be accessed directly - otherwise we might forget to normalise
//protected:
  //std::size_t NumUnknown( std::vector<bool> &bKnown ) const;
  const bool bOverlapAltNorm;
  const std::size_t NumOverlapExp;
  void OverlapValidateExp( std::size_t Exp ) const
  {
    if( Exp >= NumOverlapExp )
      throw std::domain_error( "OverlapIdx() Overlap exponent out of bounds" );
  }
  void OverlapValidate( std::size_t Exp, int SourceSink ) const
  {
    OverlapValidateExp( Exp );
    if( SourceSink < idxSrc || SourceSink > idxSnk )
      throw std::domain_error( "OverlapIdx() SourceSink out of bounds" );
  }
  std::size_t OverlapCount( std::size_t Exp ) const
  {
    OverlapValidateExp( Exp );
    return Exp < NumOverlapExpDual ? 2 : 1;
  }
  std::size_t Overlap( std::size_t Exp, int SourceSink ) const
  {
    OverlapValidate( Exp, SourceSink );
    if( Exp < NumOverlapExpDual )
      return vOverlap[SourceSink].idx + Exp;
    return vOverlap[NumOverlapExpDual ? 2 : 0].idx + Exp - NumOverlapExpDual;
  }
private:
  // NumOverlapExpDual is the number of dual overlap factors, e.g. P * W (optional)
  // If present, these are stored in Overlap[0]=Source and Overlap[1]=Sink
  // Followed by NumOverlapExp - NumOverlapExpDual single overlap factors, e.g. P^2
  // in Overlap[2], or Overlap[0] if NumOverlapExpDual==0
  std::size_t NumOverlapExpDual;
  std::vector<ModelParam> vOverlap;
};

#endif // ModelCommon_hpp
