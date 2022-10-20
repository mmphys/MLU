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
  std::size_t Guessable( std::vector<bool> &bKnown, bool bLastChance ) const override;
  void ReduceUnknown() override;
// TODO: These really shouldn't be accessed directly - otherwise we might forget to normalise
//protected:
  //std::size_t NumUnknown( std::vector<bool> &bKnown ) const;
  const bool bOverlapAltNorm;
  std::size_t NumOverlapExp;
  std::vector<ModelParam> Overlap;
};

#endif // ModelCommon_hpp
