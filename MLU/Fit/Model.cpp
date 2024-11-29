/**
 
 Abstract base for all models. i.e. interface consumed by fitter

 Source file: Model.cpp
 
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

// Keep my models separate from my fitting code

#include <MLUconfig.h>
#include "Model.hpp"

std::string Model::Args::Remove( const std::string &key, bool * Removed )
{
  std::string s;
  typename Model::Args::iterator it = find( key );
  const bool bFound{ it != end() };
  if( bFound )
  {
    s = std::move( it->second );
    erase( it );
  }
  if( Removed )
    *Removed = bFound;
  return s;
}

// Read key-Value map. In future this should work for any type - but for now, it's string only
void Model::Args::FromString( const std::string &s, bool bOptionalValue )
{
  for( std::size_t iStart = 0; ( iStart = s.find_first_not_of( MLU::WhiteSpace, iStart ) ) != std::string::npos; )
  {
    // Find where this field stops, trimming trailing whitespace
    const std::size_t NextComma = s.find_first_of( MLU::Comma, iStart );
    std::size_t iStop = ( NextComma == std::string::npos ? s.length() : NextComma );
    while( iStop - iStart && std::isspace( s[iStop - 1] ) )
      iStop--;
    std::string sKey{ s.substr( iStart, iStop - iStart ) };
    iStart = NextComma == std::string::npos ? NextComma : NextComma + 1; // Next loop starts after comma
    // Is there an equal sign?
    std::string Value;
    std::size_t iEquals = sKey.find_first_of( '=' );
    if( iEquals != std::string::npos )
    {
      // We've found an equal sign. Extract the value
      const std::size_t ValStart = sKey.find_first_not_of( MLU::WhiteSpace, iEquals + 1 );
      if( ValStart != std::string::npos )
        Value = sKey.substr( ValStart );
      // Now trim the Key
      while( iEquals && std::isspace( sKey[iEquals - 1] ) )
        iEquals--;
      sKey.resize( iEquals );
    }
    if( !bOptionalValue && Value.empty() )
      throw std::runtime_error( "No value for Key \"" + sKey );
    if( !emplace( std::pair<std::string, std::string>( sKey, Value ) ).second )
      throw std::runtime_error( "Repeated Key \"" + sKey + "\"" );
  }
}

std::string Model::Args::ToString() const
{
  std::string s;
  for( const value_type &v : *this )
  {
    if( !s.empty() )
      s.append( 1, ',' );
    s.append( v.first );
    if( !v.second.empty() )
    {
      s.append( 1, '=' );
      s.append( v.second );
    }
  }
  return s;
}

std::vector<bool> Model::CreateParams::CreateMomList( const std::string &MomIndependentOpList )
{
  std::vector<bool> v( OpNames.size(), false );
  for( const std::string &s : MLU::ArrayFromString( MomIndependentOpList ) )
  {
    int i{ MLU::IndexIgnoreCase( OpNames, s ) };
    if( i < 0 || i >= OpNames.size() )
      std::cout << "WARNING: Operator " << s << " not referred to in input files.\n";
    else
      v[i] = true;
  }
  return v;
}

Model::CreateParams::CreateParams( const std::vector<std::string> &o, const MLU::CommandLine &c )
: OpNames{ o },
  bOpMomIndependent( CreateMomList( c.SwitchValue<std::string>( "nopolap" ) ) ),
  cl{ c },
  bOverlapAltNorm{ cl.GotSwitch( MLU::sOverlapAltNorm.c_str() ) }
{
  if( bOverlapAltNorm )
    std::cout << "WARNING: Use of --" << MLU::sOverlapAltNorm << " is deprecated.\n";
}

std::string Model::CreateParams::Description() const
{
  std::string s;
  if( bOverlapAltNorm )
    s.append( "overlap alternate normalisation (deprecated)" );
  return s;
}

void Model::AddParam( Params &mp, ModelParam &ModPar, std::size_t NumExp, bool bMonotonic, Param::Type Type )
{
  param.emplace_back( &*mp.Add( ModPar.Key, NumExp, bMonotonic, Type ) );
  ModPar.param = &param.back()->second;
}

void Model::AddEnergy( Params &mp, ModelParam &ModPar, std::size_t NumExp, int N,
                       MLU::DispersionType dispType )
{
  param.emplace_back( &*mp.AddEnergy( ModPar.Key, NumExp, N, dispType ) );
  ModPar.param = &param.back()->second;
}

const std::string &Model::XVectorKeyName() const
{
  throw std::runtime_error( "Model " + Description() + " supports correlator fits only" );
}

std::string Model::XVectorKey() const
{
  throw std::runtime_error( "Model " + Description() + " supports correlator fits only" );
}

std::vector<Param::Key> Model::XVectorKeyNames() const
{
  throw std::runtime_error( "Model " + Description() + " supports correlator fits only" );
}

std::size_t Model::GetFitColumn() const
{
  throw std::runtime_error( "Model " + Description() + " supports correlator fits only" );
}

std::vector<std::string> Object::GetObjectNameSingle( const Model::CreateParams &cp, Model::Args &Args )
{
  bool bManualObjectID;
  std::string ObjectID{ Args.Remove( "ObjectID", &bManualObjectID ) };
  if( !bManualObjectID )
  {
    if( !cp.pCorr->Name_.MesonMom.empty() )
      ObjectID = cp.pCorr->Name_.MesonMom[0];
    else if( !cp.pCorr->Name_.Meson.empty() )
      ObjectID = cp.pCorr->Name_.Meson[0];
    else
    {
      ObjectID = cp.pCorr->Name_.Base;
      std::size_t Len{ ObjectID.find( '.' ) };
      if( Len != std::string::npos )
        ObjectID.resize( Len );
    }
  }
  return { ObjectID };
}

std::vector<std::string> Object::GetObjectNameSnkSrc( const Model::CreateParams &cp, Model::Args &Args )
{
  // Now remove object names
  bool bManualObjectID[2];
  std::vector<std::string> ObjectID( 2 );
  ObjectID[0] = Args.Remove( "Src", &bManualObjectID[0] );
  ObjectID[1] = Args.Remove( "Snk", &bManualObjectID[1] );
  const MLU::FileNameAtt &fna{ cp.pCorr->Name_ };
  if( !bManualObjectID[0] )
  {
    if( !fna.MesonMom.empty() )
      ObjectID[0] = fna.MesonMom[0];
    else if( fna.BaseShortParts.size() >= 2
            && MLU::EqualIgnoreCase( fna.BaseShortParts[0], "ZV" ) )
    {
      ObjectID[0] = fna.MakeMesonName( fna.BaseShortParts[1] );
      ObjectID.resize( 1 );
    }
    if( ObjectID[0].empty() )
      throw std::runtime_error( "GetObjectNameSnkSrc(): Src unavailable - specify manually" );
  }
  if( ObjectID.size() > 1 )
  {
    if( !bManualObjectID[1] )
    {
      if( fna.MesonMom.size() < 2 )
        throw std::runtime_error( "GetObjectNameSnkSrc(): Snk unavailable - specify manually" );
      ObjectID[1] = fna.MesonMom[1];
    }
    if( MLU::EqualIgnoreCase( ObjectID[0], ObjectID[1] ) )
      ObjectID.resize( 1 );
  }
  return ObjectID;
}
