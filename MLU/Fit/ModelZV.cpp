/*************************************************************************************
 
 Models for ZV extraction

 Source file: ModelZV.cpp

 Copyright (C) 2019-2023
 
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

#include <MLUconfig.h>
#include "ModelZV.hpp"

//const std::string ModelZV::sRaw{ "Raw" };
const std::string ModelZV::sZV{ "ZV" };
const std::string ModelZV::sChildModel{ "C2Model" };

ModelZV::ModelZV( const Model::CreateParams &cp, Model::Args &Args, int NumExp,
                        std::vector<std::string> &&objectID, std::vector<std::string> &&opNames )
: Model3pt(cp, Args, NumExp, std::vector<std::string>(objectID), std::vector<std::string>(opNames)),
  ChildType{ Args.Remove<eModelType>( sChildModel, eModelType::Exp ) }
{
  if( E.size() != 1 )
    throw std::runtime_error( sZV + " only supports 1 energy, not " + std::to_string( E.size() ) );
  // Number of 2-pt exponentials defaults to same as 3-pt
  int NumC2Exp{ static_cast<int>( NumOverlapExp ) };
  {
    // "C2e" argument specifies number of exponentials at source and sink
    bool b;
    std::string s{ Args.Remove( "C2e", &b ) };
    if( b )
      NumC2Exp = MLU::FromString<int>( s );
  }
  if( ChildType != eModelType::Exp && ChildType != eModelType::Cosh && ChildType != eModelType::Sinh )
  {
    std::ostringstream os;
    os << sChildModel << MLU::Space << ChildType << " invalid ZV numerator type";
    throw std::runtime_error( os.str().c_str() );
  }
  {
    std::vector<std::string> oid{ objectID[0] };
    std::vector<std::string> on{ opNames[0], opNames[1] };
    switch( ChildType )
    {
      case eModelType::Exp:
        C2.reset(new ModelExp(cp, Args, NumC2Exp, std::move(oid), std::move(on)));
        break;
      case eModelType::Cosh:
        C2.reset(new ModelCosh(cp, Args, NumC2Exp, std::move(oid), std::move(on)));
        break;
      default:
        C2.reset(new ModelSinh(cp, Args, NumC2Exp, std::move(oid), std::move(on)));
        break;
    }
  }
  ZV.Key.Object.emplace_back( cp.pCorr->Name_.BaseShortParts[1] );
  ZV.Key.Name = sZV;
}

ModelZV::ModelZV( const Model::CreateParams &cp, Model::Args &Args, int NumExp )
: ModelZV(cp, Args, NumExp, GetObjectNameZV(cp, Args), ModelOverlap::GetOpNames(cp, Args))
{
}

void ModelZV::AddParameters( Params &mp )
{
  Model3pt::AddParameters( mp );
  C2->AddParameters( mp );
  AddParam( mp, ZV, 1, false, Param::Type::Derived );
}

void ModelZV::SaveParameters( const Params &mp )
{
  Model3pt::SaveParameters( mp );
  C2->SaveParameters( mp );
  ZV.idx = mp.at( ZV.Key )();
}

// Get a descriptive string for the model
std::string ModelZV::Description() const
{
  std::string s{ sZV + "(" };
  s.append( Model3pt::Description() );
  s.append( 1, ',' );
  s.append( C2->Description() );
  s.append( 1, ')' );
  return s;
}

// For now, all we can guess are matrix elements - and we assume they are 1
void ModelZV::Guessable( ParamsPairs &PP ) const
{
  Model3pt::Guessable( PP ); // Let 3pt guess MEL & EDiff
  // I can guess the matrix element
  PP.SetState( ParamsPairs::State::Known, ZV.Key, ZV.param->size );
  // TODO: Implement guessing
}

ModelType ModelZV::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::ZV );
  return m;
}

// TODO: Make sure we add partial derivatives for all models
double ModelZV::Derivative( int t, int p ) const
{
  return 0;
}

scalar ModelZV::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  const scalar sC3{ Model3pt::operator()( t, ScratchPad, ModelParams ) };
  const scalar sC2{ (*C2)( DeltaT, ScratchPad, ModelParams ) };
  const scalar RawRatio{ sC2 / sC3 };
  ModelParams[ZV.idx] = 2 * ModelParams[E[idxSrc].idx] / ModelParams[ ( *MEL.param )( 0, 0 ) ];
  return RawRatio;
}

std::vector<std::string> ModelZV::GetObjectNameZV( const Model::CreateParams &cp, Model::Args &Args )
{
  bool bManualObjectID;
  std::string ObjectID{ Args.Remove( "ObjectID", &bManualObjectID ) };
  if( !bManualObjectID )
  {
    const MLU::FileNameAtt &fna{ cp.pCorr->Name_ };
    if( fna.BaseShortParts.size() < 2 )
      throw std::runtime_error( "ObjectID unknown - specify manually" );
    ObjectID = fna.MakeMesonName( fna.BaseShortParts[1] );
    MLU::Momentum p0( 0 );
    ObjectID.append( p0.FileString( MLU::Momentum::DefaultPrefix ) );
  }
  return { ObjectID };
}

