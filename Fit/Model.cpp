/*************************************************************************************
 
 Abstract base for all models. i.e. interface consumed by fitter

 Source file: Model.cpp
 
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

// Keep my models separate from my fitting code

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
  for( std::size_t iStart = 0; ( iStart = s.find_first_not_of( Common::WhiteSpace, iStart ) ) != std::string::npos; )
  {
    // Find where this field stops, trimming trailing whitespace
    const std::size_t NextComma = s.find_first_of( Common::Comma, iStart );
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
      const std::size_t ValStart = sKey.find_first_not_of( Common::WhiteSpace, iEquals + 1 );
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

Model::CreateParams::CreateParams( const std::vector<std::string> &o, const Common::CommandLine &c )
: OpNames{ o }, cl{ c },
  bOverlapAltNorm{ cl.GotSwitch( Common::sOverlapAltNorm.c_str() ) }
{
  if( bOverlapAltNorm )
    std::cout << "WARNING: Use of --" << Common::sOverlapAltNorm << " is deprecated.\n";
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

void Model::AddEnergy( Params &mp, ModelParam &ModPar, std::size_t NumExp, int N, bool bEnablePHat )
{
  param.emplace_back( &*mp.AddEnergy( ModPar.Key, NumExp, N, bEnablePHat ) );
  ModPar.param = &param.back()->second;
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
    ObjectID = cp.pCorr->Name_.Base;
    std::size_t Len{ ObjectID.find( '.' ) };
    if( Len != std::string::npos )
      ObjectID.resize( Len );
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
  if( !bManualObjectID[0] )
  {
    if( cp.pCorr->Name_.MesonMom.empty() )
      throw std::runtime_error( "GetObjectNameSnkSrc(): Src unavailable - specify manually" );
    ObjectID[0] = cp.pCorr->Name_.MesonMom[0];
  }
  if( !bManualObjectID[1] )
  {
    if( cp.pCorr->Name_.MesonMom.size() < 2 )
      throw std::runtime_error( "GetObjectNameSnkSrc(): Snk unavailable - specify manually" );
    ObjectID[1] = cp.pCorr->Name_.MesonMom[1];
  }
  if( Common::EqualIgnoreCase( ObjectID[0], ObjectID[1] ) )
    ObjectID.resize( 1 );
  return { ObjectID };
}
