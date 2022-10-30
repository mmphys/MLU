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

// This is the glue between abstract models and specific implementations
// Other than here, these must be kept strictly separate

#include "Model2pt.hpp"
#include "Model3pt.hpp"
#include "ModelConstant.hpp"
#include "ModelRatio.hpp"

const std::string E{ "E" };
const std::string EDiff{ "EDiff" };

static const std::string sModelTypeUnknown{ "Unknown" };
static const std::string sModelTypeExp{ "Exp" };
static const std::string sModelTypeCosh{ "Cosh" };
static const std::string sModelTypeSinh{ "Sinh" };
static const std::string sModelTypeThreePoint{ "3pt" };
static const std::string sModelTypeConstant{ "Const" };
static const std::string sModelTypeRatio{ "R3" };

std::ostream & operator<<( std::ostream &os, const ModelType m )
{
  switch( m )
  {
    case ModelType::Exp:
      os << sModelTypeExp;
      break;
    case ModelType::Cosh:
      os << sModelTypeCosh;
      break;
    case ModelType::Sinh:
      os << sModelTypeSinh;
      break;
    case ModelType::ThreePoint:
      os << sModelTypeThreePoint;
      break;
    case ModelType::Constant:
      os << sModelTypeConstant;
      break;
    case ModelType::R3:
      os << sModelTypeRatio;
      break;
    default:
      os << sModelTypeUnknown;
      break;
  }
  return os;
}

std::istream & operator>>( std::istream &is, ModelType &m )
{
  std::string s;
  if( is >> s )
  {
    if( Common::EqualIgnoreCase( s, sModelTypeExp ) )
      m = ModelType::Exp;
    else if( Common::EqualIgnoreCase( s, sModelTypeCosh ) )
      m = ModelType::Cosh;
    else if( Common::EqualIgnoreCase( s, sModelTypeSinh ) )
      m = ModelType::Sinh;
    else if( Common::EqualIgnoreCase( s, sModelTypeThreePoint ) )
      m = ModelType::ThreePoint;
    else if( Common::EqualIgnoreCase( s, sModelTypeConstant ) )
      m = ModelType::Constant;
    else if( Common::EqualIgnoreCase( s, sModelTypeRatio ) )
      m = ModelType::R3;
    else
    {
      m = ModelType::Unknown;
      is.setstate( std::ios_base::failbit );
    }
  }
  return is;
}

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

Model::CreateParams::CreateParams( const std::vector<std::string> &o, const Common::CommandLine &c )
: OpNames{ o }, cl{ c },
  NumExponents{ cl.SwitchValue<int>( "e" ) },
  bOverlapAltNorm{ cl.GotSwitch( Common::sOverlapAltNorm.c_str() ) }
{
  if( bOverlapAltNorm )
    std::cout << "WARNING: Use of --" << Common::sOverlapAltNorm << " is deprecated.\n";
}

// Model constructor
Model::Model( const CreateParams &cp, Model::Args &Args )
: Nt{cp.pCorr->NtUnfolded}, NumExponents{ Args.Remove( "e", cp.NumExponents ) }
{
}

void Model::AddParam( Params &mp, ModelParam &ModPar, std::size_t NumExp, bool bMonotonic, Param::Type Type )
{
  param.emplace_back( &*mp.Add( ModPar.Key, NumExp, bMonotonic, Type ) );
  ModPar.param = &param.back()->second;
}

// Create a model of the appropriate type - this is the only place with knowledge of this mapping
ModelPtr Model::MakeModel( const Model::CreateParams &cp, Model::Args &Args )
{
  // Now work out what type of model we are creating
  ModelType modelType{ ModelType::Unknown };
  bool bGotModelType{ false };
  {
    std::string s{ Args.Remove( "Model", &bGotModelType ) };
    if( bGotModelType )
    {
      std::istringstream ss( s );
      if( !( ss >> modelType ) || !Common::StreamEmpty( ss ) )
        throw std::runtime_error( "Unknown ModelType \"" + s + "\"" );
    }
  }
  if( !bGotModelType )
  {
    // We haven't been told which model to use. Choose a suitable default
    const bool b3pt{ cp.pCorr->Name_.bGotDeltaT && cp.pCorr->Name_.Gamma.size() == 1 };
    if( b3pt )
    {
      if( !cp.pCorr->Name_.BaseShortParts.empty()
         && Common::EqualIgnoreCase( cp.pCorr->Name_.BaseShortParts[0], "R3" ) )
        modelType = ModelType::R3;
      else
        modelType = ModelType::ThreePoint;
    }
    else switch( cp.pCorr->parity )
    {
      case Common::Parity::Even:
        modelType = ModelType::Cosh;
        break;
      case Common::Parity::Odd:
        modelType = ModelType::Sinh;
        break;
      default:
        modelType = ModelType::Exp;
        break;
    }
  }
  // Make the model
  std::cout << "  " << modelType;
  for( typename Model::Args::value_type v : Args )
  {
    std::cout << Common::CommaSpace << v.first;
    if( !v.second.empty() )
      std::cout << '=' << v.second;
  }
  std::cout << Common::NewLine;
  ModelPtr model;
  switch( modelType )
  {
    case ModelType::Exp:
      model.reset( new ModelExp( cp, Args ) );
      break;
    case ModelType::Cosh:
      model.reset( new ModelCosh( cp, Args ) );
      break;
    case ModelType::Sinh:
      model.reset( new ModelSinh( cp, Args ) );
      break;
    case ModelType::ThreePoint:
      model.reset( new Model3pt( cp, Args ) );
      break;
    case ModelType::Constant:
      model.reset( new ModelConstant( cp, Args ) );
      break;
    case ModelType::R3:
      model.reset( new ModelRatio( cp, Args ) );
      break;
    default:
      throw std::runtime_error( "Unrecognised ModelType for " + cp.pCorr->Name_.Filename );
  }
  // 2 part construction - now that virtual functions in place
  return model;
}
