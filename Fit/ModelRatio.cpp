/*************************************************************************************
 
 Models for Ratios

 Source file: ModelRatio.cpp

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

#include "ModelRatio.hpp"

const std::string ModelRatio::sRaw{ "Raw" };

ModelRatio::ModelRatio( const Model::CreateParams &cp, Model::Args &Args )
: Model3pt( cp, Args )
{
  static const std::string ChildModel{ "C2Model" };
  ModelType ChildType = Args.Remove<ModelType>( ChildModel, ModelType::Exp );
  if( ChildType != ModelType::Exp && ChildType != ModelType::Cosh && ChildType != ModelType::Sinh )
  {
    std::ostringstream os;
    os << ChildModel << Common::Space << ChildType << " invalid for ratio children";
    throw std::runtime_error( os.str().c_str() );
  }
  C2.reserve( E.size() );
  for( std::size_t i = 0 ; i < E.size(); ++i )
  {
    std::vector<std::string> oid( 1 );
    oid[0] = objectID[i];
    switch( ChildType )
    {
      case ModelType::Exp:
        C2.emplace_back( new ModelExp( cp, Args, std::move( oid ), NumOverlapExp ) );
        break;
      case ModelType::Cosh:
        C2.emplace_back( new ModelCosh( cp, Args, std::move( oid ), NumOverlapExp ) );
        break;
      default:
        C2.emplace_back( new ModelSinh( cp, Args, std::move( oid ), NumOverlapExp ) );
        break;
    }
  }
  // The raw option means the ratio was constructed without the overlap coefficients
  Args.Remove( sRaw, &bRaw );
}

void ModelRatio::AddParameters( Params &mp )
{
  Model3pt::AddParameters( mp );
  for( std::unique_ptr<Model2pt> &m : C2 )
    m->AddParameters( mp );
}

void ModelRatio::SaveParameters( const Params &mp )
{
  Model3pt::SaveParameters( mp );
  for( std::unique_ptr<Model2pt> &m : C2 )
    m->SaveParameters( mp );
}

// Get a descriptive string for the model
std::string ModelRatio::Description() const
{
  std::string s{ "R3(" };
  s.append( ObjectID( idxSrc ) );
  if( objectID.size() > 1 )
  {
    s.append( 1, ',' );
    s.append( ObjectID( idxSnk ) );
  }
  s.append( ModelOverlap::Description() );
  s.append( 1, ')' );
  return s;
}

// TODO: Make sure we add partial derivatives for all models
double ModelRatio::Derivative( int t, int p ) const
{
  return 0;
}

scalar ModelRatio::operator()( int t, Vector &ScratchPad, Vector &ModelParams ) const
{
  const scalar C3{ Model3pt::operator()( t, ScratchPad, ModelParams ) };
  const scalar CSource{ (*C2[idxSrc])( t, ScratchPad, ModelParams ) };
  const scalar CSink{ (*C2.back())( DeltaT - t, ScratchPad, ModelParams ) };
  const scalar CProd{ CSource * CSink };
  const scalar RawRatio{ C3 / CProd };
  // Assumption here is that the overlap factors at source and sink are the same
  assert( C2[idxSrc]->Overlap( 0, idxSrc ) == C2[idxSrc]->Overlap( 0, idxSnk ) );
  assert( C2.back()->Overlap( 0, idxSrc ) == C2.back()->Overlap( 0, idxSnk ) );
  const scalar ZSource{ ModelParams[C2[idxSrc]->Overlap( 0, idxSrc )] };
  const scalar ZSink{ ModelParams[C2.back()->Overlap( 0, idxSrc )] };
  scalar ZProd{ ZSource * ZSink };
  if( bOverlapAltNorm ) // DEPRECATED
  {
    // DEPRECATED: If the overlap factors absorbed 1/2E, we need to get rid of that
    if( E.size() < 2 )
      ZProd *= 2 * ModelParams[E[idxSrc].idx];
    else
      ZProd *= 2 * std::sqrt( ModelParams[E[idxSrc].idx] * ModelParams[E[idxSnk].idx] );
  }
  scalar Result{ RawRatio };
  if( bRaw )
    ModelParams[R3Raw.idx] = ModelParams[ ( *MEL.param )( 0, 0 ) ] / ZProd;
  else
    Result *= ZProd;
  return Result;
}
