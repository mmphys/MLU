/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: FitModel.cpp
 
 Copyright (C) 2020
 
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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

#include "MultiFit.hpp"

static const std::string sModelTypeUnknown{ "Unknown" };
static const std::string sModelTypeExp{ "Exp" };
static const std::string sModelTypeCosh{ "Cosh" };
static const std::string sModelTypeSinh{ "Sinh" };
static const std::string sModelTypeThreePoint{ "3pt" };
static const std::string sModelTypeConstant{ "Const" };

inline std::ostream & operator<<( std::ostream &os, const ModelType m )
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
    else
    {
      m = ModelType::Unknown;
      is.setstate( std::ios_base::failbit );
    }
  }
  return is;
}

// ParamNamesPerExp are all overlap coefficients
class ModelOverlap : public Model
{
public:
  int NumOverlap;
protected:
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
public:
  virtual std::size_t UnknownParameterCount( UniqueNames &Names ) const;
  virtual void ReduceUnknown( const UniqueNames &Names );
  virtual bool Guess( Vector &Parameters, std::vector<bool> &bKnown, int pass, const Vector &Corr ) const;
};

class Model2pt : public ModelOverlap
{
public:
  //Model2pt( const Fold &f, int NumExponents, int NumOps );
  //Model2pt( const Fold &f, const Model2pt &Other );
  virtual inline double Derivative( int t, int p ) const;
  virtual bool GuessE0( scalar &E0, const Vector &Corr ) const;
  virtual int GetScratchPadSize() const { return NumExponents; }
  virtual void UpdateScratchPad( Vector &ScratchPad, const Vector &ModelParams ) const;
};

class ModelExp : public Model2pt
{
public:
  virtual scalar operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const;
};

class ModelSinh : public Model2pt
{
public:
  virtual scalar operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const;
};

class ModelCosh : public Model2pt
{
public:
  virtual scalar operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const;
};

class ModelThreePoint : public ModelOverlap
{
protected:
  inline void SetParams( int Nt_, int HalfNt_, int NumExponents_ ) { ModelOverlap::SetParams( Nt_, HalfNt_, 1 ); }
public:
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  virtual scalar operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const;
};

class ModelConstant : public Model
{
public:
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  virtual scalar operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const;
};

// 2nd phase of construction (now that virtual functions in place)
void ModelOverlap::Construct(vString &Params,const ModelDefaultParams &Default,const Fold &Corr,const vString &OpName)
{
  Model::Construct( Params, Default, Corr, OpName );
  assert( ParamNamesPerExp.size() == 1 && "This needs to be fixed. I assumed only energy prior" );
  NumOverlap = 0;
  // Have we been given source and sink names?
  std::string Snk_Src;
  if( !Params.empty() )
  {
    Snk_Src = std::move( Params[0] );
    Params.erase( Params.begin() );
  }
  if( Snk_Src.empty() )
  {
    // Take sink and source from filename
    assert( std::min( idxSrc, idxSnk ) == 0 );
    assert( std::max( idxSrc, idxSnk ) == 1 );
    for( int i = 0; i <= 1; ++i )
    {
      std::string ParamName{ OpName[Corr.Name_.op[i]] };
      if( Default.bFactor )
      {
        if( i == 0 || !Common::EqualIgnoreCase( ParamName, ParamNamesPerExp[1] ) )
        {
          ParamNamesPerExp.push_back( ParamName );
          NumOverlap++;
        }
      }
      else
      {
        ParamName.append( 1, '_' );
        ParamName.append( pSrcSnk[i] );
        ParamNamesPerExp.push_back( ParamName );
        NumOverlap++;
      }
    }
  }
  else
  {
    // Use the sink and source names we've been given
    vString OpNames{ Common::ArrayFromString( Params[0], Common::Underscore ) };
    if(OpNames.size() < 1 || OpNames.size() > 2 || OpNames[0].empty() || ( OpNames.size() >= 2 && OpNames[1].empty() ))
      throw std::runtime_error( "Invalid Snk_Src operator \"" + Params[0] + "\"" );
    ParamNamesPerExp.push_back( OpNames[0] );
    NumOverlap++;
    if( OpNames.size() >= 2 && !Common::EqualIgnoreCase( OpNames[0], OpNames[1] ) )
    {
      ParamNamesPerExp.push_back( OpNames[1] );
      NumOverlap++;
    }
  }
}

// NOT COUNTING ENERGY, put names of parameters that can be determined in list, return number remaining to determine
// For now, synonymous with overlap coefficients
std::size_t ModelOverlap::UnknownParameterCount( UniqueNames &Names ) const
{
  assert( NumOverlap >= 1 && "Should be at least one overlap coefficient" );
  assert( NumOverlap <= 2 && "Update code to support more than two overlap coefficients" );
  std::size_t NumUnknown{ 0 };
  if( NumOverlap == 1 )
    Names[ParamNamesPerExp[0]];
  else
  {
    bool bKnow0{ Names.find( ParamNamesPerExp[0] ) != Names.end() };
    bool bKnow1{ Names.find( ParamNamesPerExp[1] ) != Names.end() };
    if( !bKnow0 && !bKnow1 )
      NumUnknown = 2;
    else
    {
      // I know at least one of them, so I in fact know both
      Names[ParamNamesPerExp[0]];
      Names[ParamNamesPerExp[1]];
    }
  }
  return NumUnknown;
}

void ModelOverlap::ReduceUnknown( const UniqueNames &Names )
{
  if( NumOverlap == 2 && !Names.count( ParamNamesPerExp[0] ) && !Names.count( ParamNamesPerExp[1] ) )
  {
    ParamNamesPerExp[0] = ParamNamesPerExp[idxSnk] + ParamNamesPerExp[idxSrc];
    ParamNamesPerExp.resize( 1 );
    NumOverlap = 1;
  }
}

// Take a guess as to my parameters
bool ModelOverlap::Guess( Vector &Parameters, std::vector<bool> &bKnown, int pass, const Vector &Corr ) const
{
  static const scalar MELFactor{ std::sqrt( 2 ) };
  const int tGuess{ HalfNt / 2 };
  scalar E0{ Parameters[ParamIdxPerExp[0][0]] };
  bool bNeedAnotherPass{ false };
  bool bGuessSame{ false };
  if( NumOverlap == 1 )
    bGuessSame = true;
  else
  {
    // Two overlap coefficients
    int idx1{ ParamIdxPerExp[0][1] };
    int idx2{ ParamIdxPerExp[0][2] };
    bool bKnow1{ bKnown[idx1] };
    bool bKnow2{ bKnown[idx2] };
    if( bKnow1 != bKnow2 )
    {
      // We know one, but not the other
      scalar Product{ std::abs( Corr[tGuess] ) * std::exp( E0 * tGuess ) };
      int idxKnown  { bKnow1 ? idx1 : idx2 };
      int idxUnknown{ bKnow1 ? idx2 : idx1 };
      scalar MELGuess = Product / Parameters[idxKnown];
      bKnown[idxUnknown] = true;
      Parameters[idxUnknown] = MELGuess;
      for( int e = 1; e < NumExponents; ++e )
      {
        for( int op = 1; op <= NumOverlap; ++op )
        {
          const int idx{ ParamIdxPerExp[e][op] };
          if( !bKnown[idx] )
          {
            bKnown[idx] = true;
            Parameters[idx] = Parameters[ParamIdxPerExp[e - 1][op]] * MELFactor;
          }
        }
      }
    }
    else if( !bKnow1 && !bKnow2 )
    {
      if( pass == 0 )
        bNeedAnotherPass = true;
      else
        bGuessSame = true; // We should solve simultaneous equations to get the overlap coefficients
    }
  }
  if( bGuessSame )
  {
    scalar MELGuess{ std::sqrt( std::abs( Corr[tGuess] ) ) * std::exp( E0 * tGuess / 2 ) };
    for( int e = 0; e < NumExponents; ++e )
    {
      if( e )
        MELGuess *= MELFactor;
      for( int op = 1; op <= NumOverlap; ++op )
      {
        const int idx{ ParamIdxPerExp[e][op] };
        if( !bKnown[idx] )
        {
          bKnown[idx] = true;
          Parameters[idx] = MELGuess;
        }
      }
    }
  }
  return bNeedAnotherPass;
}

/*Model::Model( const Fold &f, int NumExponents_, int NumOps_ )
: Nt{ f.Nt() },
  HalfNt{ f.NtUnfolded / 2 }, // NB: This is half of the ORIGINAL correlator time dimension
  NumExponents{ NumExponents_ },
  NumOps{ NumOps_ }
{
}

Model2pt::Model2pt( const Fold &f, int NumExponents_, int NumOps_ )
: Model( f, NumExponents_, NumOps_ ),
  snk( NumExponents ),
  src( NumExponents ),
  SinhCoshAdjust( NumExponents )
{
  parity = f.parity;
}

Model2pt::Model2pt( const Fold &f, const Model2pt &Other ) : Model2pt( f, Other.NumExponents, Other.NumOps )
{
}

void ThreadModel::Init( const scalar * ModelParameters, std::size_t NumParams, std::size_t Stride )
{
  Vector NewParams;
  NewParams.MapView( const_cast<scalar *>( ModelParameters ), NumParams, Stride );
  if( NewParams != Params )
  {
    Params = NewParams;
    ParamsChanged();
  }
}*/

void ModelThreePoint::Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr,
                                 const vString &OpName )
{
  throw std::runtime_error( "Implement ModelThreePoint" );
}

scalar ModelThreePoint::operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const
{
  return 0;
}

void ModelConstant::Construct(vString &Params,const ModelDefaultParams &Default,const Fold &Corr,const vString &OpName)
{
  throw std::runtime_error( "Implement ModelConstant" );
}

scalar ModelConstant::operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const
{
  return 0;
}

bool Model2pt::GuessE0( scalar &E0, const Vector &Corr ) const
{
  const int tGuess{ HalfNt / 2 };
  E0 = std::log( Corr[tGuess] / Corr[tGuess + 1] );
  return true;
}

void Model2pt::UpdateScratchPad( Vector &ScratchPad, const Vector &ModelParams ) const
{
  for( int e = 0; e < NumExponents; ++e )
    ScratchPad[e] = std::exp( - ModelParams[ParamIdxPerExp[e][0]] * HalfNt );
}

scalar ModelExp::operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const
{
  double z = 0;
  for( int e = 0; e < NumExponents; ++e )
  {
    double d = std::exp( - ModelParams[ParamIdxPerExp[e][0]] * t );
    z += d * ModelParams[ParamIdxPerExp[e][1]] * ModelParams[ParamIdxPerExp[e][NumOverlap > 1 ? 2 : 1]];
  }
  return z;
}

scalar ModelCosh::operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const
{
  double z = 0;
  for( int e = 0; e < NumExponents; ++e )
  {
    double d = ScratchPad[e] * std::cosh( - ModelParams[ParamIdxPerExp[e][0]] * ( t - HalfNt ) );
    z += d * ModelParams[ParamIdxPerExp[e][1]] * ModelParams[ParamIdxPerExp[e][NumOverlap > 1 ? 2 : 1]];
  }
  return z;
}

scalar ModelSinh::operator()( int t, const Vector &ModelParams, Vector &ScratchPad ) const
{
  double z = 0;
  for( int e = 0; e < NumExponents; ++e )
  {
    double d = ScratchPad[e] * std::sinh( - ModelParams[ParamIdxPerExp[e][0]] * ( t - HalfNt ) );
    z += d * ModelParams[ParamIdxPerExp[e][1]] * ModelParams[ParamIdxPerExp[e][NumOverlap > 1 ? 2 : 1]];
  }
  return z;
}

// TODO: Make sure we add partial derivatives for all models
double Model2pt::Derivative( int t, int p ) const
{
  double d = 0;
  // TODO re-implement
  /*for( int e = 0; e < NumExponents; ++e )
  {
    if( Energy[e] == p )
    {
      // Derivative wrt Energy
      d = t * std::exp( - Params[Energy[e]] * t );
      if( parity == Common::Parity::Even || parity == Common::Parity::Odd )
      {
        int NtMinusT = parity == Common::Parity::Even ? Nt - t : t - Nt;
        d += NtMinusT * std::exp( - Params[Energy[e]] * ( Nt - t ) );
        d *= 0.5;
      }
      d *= - Params[src[e]] * Params[snk[e]];
      break;
    }
    // Perhaps this is a derivative wrt overlap coefficient
    // Work out whether this operator is source or sink or both
    double OtherOverlap = 0;
    int Factor{ 0 };
    if( snk[e] == p )
    {
      OtherOverlap = Params[src[e]];
      ++Factor;
    }
    if( src[e] == p )
    {
      OtherOverlap = Params[snk[e]]; // Doesn't matter if we do this twice ... it's the same operator!
      ++Factor;
    }
    if( Factor )
    {
      switch( parity )
      {
        case Common::Parity::Even:
          d = SinhCoshAdjust[e] * std::cosh( - Params[Energy[e]] * ( t - HalfNt ) );
          break;
        case Common::Parity::Odd:
          d = SinhCoshAdjust[e] * std::sinh( - Params[Energy[e]] * ( t - HalfNt ) );
          break;
        default:
          d = std::exp( - Params[Energy[e]] * t );
          break;
      }
      d *= Factor * OtherOverlap;
      break;
    }
  }*/
  return d;
}

// Create a model of the appropriate type - this is the only place with knowledge of this mapping
ModelPtr Model::MakeModel(vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &FileOps)
{
  // Now work out what type of model we are creating
  ModelType modelType{ ModelType::Unknown };
  bool bGotModelType{ false };
  // Use the model we have been told to use
  if( !Params.empty() )
  {
    if( !Params[0].empty() )
    {
      std::istringstream ss( Params[0] );
      if( !( ss >> modelType ) || !Common::StreamEmpty( ss ) )
        throw std::runtime_error( "Unknown ModelType \"" + Params[0] + "\"" );
      bGotModelType = true;
    }
    Params.erase( Params.begin() );
  }
  if( !bGotModelType )
  {
    // We haven't been told which model to use. Choose a suitable default
    const bool b3pt{ Corr.Name_.bGotDeltaT && Corr.Name_.Gamma.size() == 1 };
    if( b3pt )
      modelType = ModelType::ThreePoint;
    else switch( Corr.parity )
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
  ModelPtr model;
  switch( modelType )
  {
    case ModelType::Exp:
      model.reset( new ModelExp );
      break;
    case ModelType::Cosh:
      model.reset( new ModelCosh );
      break;
    case ModelType::Sinh:
      model.reset( new ModelSinh );
      break;
    case ModelType::ThreePoint:
      model.reset( new ModelThreePoint );
      break;
    case ModelType::Constant:
      model.reset( new ModelConstant );
      break;
    default:
      throw std::runtime_error( "Unrecognised ModelType for " + Corr.Name_.Filename );
  }
  // 2 part construction - now that virtual functions in place
  model->SetParams( Corr.Nt(), Corr.NtUnfolded, Default.NumExponents );
  model->Construct( Params, Default, Corr, FileOps );
  return model;
}
