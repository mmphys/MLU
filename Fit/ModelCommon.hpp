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

// A model with two overlap coefficients
struct ModelOverlap : Model, Object
{
  ModelOverlap( const Model::CreateParams &cp, Model::Args &Args,
                std::vector<std::string> &&ObjectID, std::size_t NumOverlapExp );
  void AddParameters( struct Params &mp ) override;
  void SaveParameters( const struct Params &mp ) override;
  std::string Description() const override;
  std::size_t Guessable( std::vector<bool> &bKnown, bool bLastChance ) const override;
  void ReduceUnknown() override;
protected:
  std::size_t NumUnknown( std::vector<bool> &bKnown ) const;
  const bool bNormalisationByEnergy;
  std::size_t NumOverlapExp;
  std::vector<ModelParam> Overlap;
};

#endif // ModelCommon_hpp
