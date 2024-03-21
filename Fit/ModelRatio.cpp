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
const std::string ModelRatio::sR3{ "R3" };
const std::string ModelRatio::sChildModel{ "C2Model" };

ModelRatio::ModelRatio( const Model::CreateParams &cp, Model::Args &Args, int NumExp,
                        std::vector<std::string> &&objectID, std::vector<std::string> &&opNames )
: Model3pt(cp, Args, NumExp, std::vector<std::string>(objectID), std::vector<std::string>(opNames)),
  ChildType{ Args.Remove<eModelType>( sChildModel, eModelType::Exp ) }
{
  // Number of 2-pt exponentials defaults to same as 3-pt
  std::vector<int> NumC2Exp( E.size(), static_cast<int>( NumOverlapExp ) );
  {
    bool b[3];
    std::array<std::string, 3> s;
    // "C2e" argument specifies number of exponentials at source and sink
    s[2] = Args.Remove( "C2e", &b[2] );
    if( b[2] )
    {
      const int nE{ MLU::FromString<int>( s[2] ) };
      for( std::size_t i = 0; i < NumC2Exp.size(); ++i )
        NumC2Exp[i] = nE;
    }
    if( NumC2Exp.size() > 1 )
    {
      // If source and sink are different, num exponents can (optionally) be specified individually
      s[idxSrc] = Args.Remove( "C2eSrc", &b[idxSrc] );
      s[idxSnk] = Args.Remove( "C2eSnk", &b[idxSnk] );
      if( b[2] && b[idxSrc] && b[idxSnk] )
        throw std::runtime_error( "C2e, C2eSrc and C2eSnk cannot all be specified" );
      for( int i = 0; i < 2; ++i )
        if( b[i] )
          NumC2Exp[i] = MLU::FromString<int>( s[i] );
    }
  }
  if( ChildType != eModelType::Exp && ChildType != eModelType::Cosh && ChildType != eModelType::Sinh )
  {
    std::ostringstream os;
    os << sChildModel << MLU::Space << ChildType << " invalid for ratio children";
    throw std::runtime_error( os.str().c_str() );
  }
  C2.reserve( E.size() );
  for( std::size_t i = 0 ; i < E.size(); ++i )
  {
    std::vector<std::string> oid{ objectID[i] };
    std::vector<std::string> on{ opNames[i], opNames[i] };
    switch( ChildType )
    {
      case eModelType::Exp:
        C2.emplace_back(new ModelExp(cp, Args, NumC2Exp[i], std::move(oid), std::move(on)));
        break;
      case eModelType::Cosh:
        C2.emplace_back(new ModelCosh(cp, Args, NumC2Exp[i], std::move(oid), std::move(on)));
        break;
      default:
        C2.emplace_back(new ModelSinh(cp, Args, NumC2Exp[i], std::move(oid), std::move(on)));
        break;
    }
  }
  if( cp.pCorr->Name_.Gamma.size() != 1 )
    throw std::runtime_error( "More than one gamma unsupported" );
  R3Raw.Key.Object = objectID;
  R3Raw.Key.Name = sR3 + MLU::Gamma::NameShort( cp.pCorr->Name_.Gamma[0] ) + sRaw;
  // The raw option means the ratio was constructed without the overlap coefficients
  Args.Remove( sRaw, &bRaw );
}

void ModelRatio::AddParameters( Params &mp )
{
  Model3pt::AddParameters( mp );
  for( std::unique_ptr<Model2pt> &m : C2 )
    m->AddParameters( mp );
  AddParam( mp, R3Raw, 1, false, Param::Type::Derived );
}

void ModelRatio::SaveParameters( const Params &mp )
{
  Model3pt::SaveParameters( mp );
  for( std::unique_ptr<Model2pt> &m : C2 )
    m->SaveParameters( mp );
  R3Raw.idx = mp.at( R3Raw.Key )();
}

// Get a descriptive string for the model
std::string ModelRatio::Description() const
{
  std::string s{ "R3(" };
  s.append( Model3pt::Description() );
  for( const std::unique_ptr<Model2pt> & p : C2 )
  {
    s.append( 1, ',' );
    s.append( p->Description() );
  }
  s.append( 1, ')' );
  return s;
}

// For now, all we can guess are matrix elements - and we assume they are 1
void ModelRatio::Guessable( ParamsPairs &PP ) const
{
  Model3pt::Guessable( PP ); // Let 3pt guess MEL & EDiff
  // I can guess the matrix element
  PP.SetState( ParamsPairs::State::Known, R3Raw.Key, R3Raw.param->size );
  // TODO: Implement guessing
}

ModelType ModelRatio::Type() const
{
  ModelType m;
  m.t = static_cast<int>( eModelType::R3 );
  return m;
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
  // Assumption here is that the overlap factors at source and sink of 2pt functions are the same
  for( const std::unique_ptr<Model2pt> &C : C2 )
    assert( C->OverlapCount( 0 ) == 1 );
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
  ModelParams[R3Raw.idx] = ModelParams[ ( *MEL.param )( 0, 0 ) ] / ZProd;
  return bRaw ? RawRatio : RawRatio * ZProd;
}
