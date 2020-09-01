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
protected:
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
public:
  virtual std::size_t UnknownParameterCount( UniqueNames &Names ) const;
  virtual void ReduceUnknown( const UniqueNames &Names );
};

class Model2pt : public ModelOverlap
{
public:
//protected:
  Common::Parity parity;
  std::vector<int> Energy;
  std::vector<int> snk;
  std::vector<int> src;
  std::vector<double> SinhCoshAdjust;
public:
  //Model2pt( const Fold &f, int NumExponents, int NumOps );
  //Model2pt( const Fold &f, const Model2pt &Other );
  virtual inline double operator()( int t ) const;
  virtual inline double Derivative( int t, int p ) const;
  virtual void ParamsChanged();
};

class ModelExp : public Model2pt
{
};

class ModelSinh : public Model2pt
{
};

class ModelCosh : public Model2pt
{
};

class ModelThreePoint : public ModelOverlap
{
protected:
  inline void SetParams( int Nt_, int HalfNt_, int NumExponents_ ) { ModelOverlap::SetParams( Nt_, HalfNt_, 1 ); }
public:
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  virtual double operator()( int t ) const;
  virtual void ParamsChanged();
};

class ModelConstant : public Model
{
public:
  virtual void Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr, const vString &OpName );
  virtual double operator()( int t ) const;
  virtual void ParamsChanged();
};

std::size_t ModelOverlap::UnknownParameterCount( UniqueNames &Names ) const
{
  constexpr int idxOverlap{ 1 };
  const int NumOverlap{ static_cast<int>( ParamNamesPerExp.size() - idxOverlap ) };
  assert( NumOverlap >= 1 && "Should be at least one overlap coefficient" );
  assert( NumOverlap <= 2 && "Update code to support more than two overlap coefficients" );
  // I will be able to solve for energies
  std::size_t NumUnknown{ 0 };
  if( NumOverlap == 1 )
    Names[ParamNamesPerExp[idxOverlap]];
  else
  {
    bool bKnow0{ Names.find( ParamNamesPerExp[idxOverlap] ) != Names.end() };
    bool bKnow1{ Names.find( ParamNamesPerExp[idxOverlap + 1] ) != Names.end() };
    if( !bKnow0 && !bKnow1 )
      NumUnknown = 2;
    else
    {
      // I know at least one of them, so I in fact know both
      Names[ParamNamesPerExp[1]];
      Names[ParamNamesPerExp[2]];
    }
  }
  return NumUnknown;
}

void ModelOverlap::ReduceUnknown( const UniqueNames &Names )
{
  if( ParamNamesPerExp.size() == 2 && !Names.count( ParamNamesPerExp[0] ) && !Names.count( ParamNamesPerExp[1] ) )
  {
    ParamNamesPerExp[0] = ParamNamesPerExp[idxSnk] + ParamNamesPerExp[idxSrc];
    ParamNamesPerExp.resize( 1 );
  }
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

// 2nd phase of construction (now that virtual functions in place)
void ModelOverlap::Construct(vString &Params,const ModelDefaultParams &Default,const Fold &Corr,const vString &OpName)
{
  Model::Construct( Params, Default, Corr, OpName );
  const std::size_t p{ ParamNamesPerExp.size() };
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
        if( i == 0 || !Common::EqualIgnoreCase( ParamName, ParamNamesPerExp[p] ) )
          ParamNamesPerExp.push_back( ParamName );
      }
      else
      {
        ParamName.append( 1, '_' );
        ParamName.append( pSrcSnk[i] );
        ParamNamesPerExp.push_back( ParamName );
      }
    }
  }
  else
  {
    // Use the sink and source names we've been given
    vString OpNames{ Common::ArrayFromString( Params[0], Common::Underscore ) };
    if( OpNames.size() < 1 || OpNames.size() > 2 || OpNames[0].empty() || ( OpNames.size() >= 2 && OpNames[1].empty() ) )
      throw std::runtime_error( "Invalid Snk_Src operator \"" + Params[0] + "\"" );
    ParamNamesPerExp.push_back( OpNames[0] );
    if( OpNames.size() >= 2 && !Common::EqualIgnoreCase( OpNames[0], OpNames[1] ) )
      ParamNamesPerExp.push_back( OpNames[1] );
  }
}

void ModelThreePoint::Construct( vString &Params, const ModelDefaultParams &Default, const Fold &Corr,
                                 const vString &OpName )
{
  throw std::runtime_error( "Implement ModelThreePoint" );
}

double ModelThreePoint::operator()( int t ) const
{
  return 0;
}

void ModelThreePoint::ParamsChanged()
{
}

void ModelConstant::Construct(vString &Params,const ModelDefaultParams &Default,const Fold &Corr,const vString &OpName)
{
  throw std::runtime_error( "Implement ModelConstant" );
}

double ModelConstant::operator()( int t ) const
{
  return 0;
}

void ModelConstant::ParamsChanged()
{
}

void Model2pt::ParamsChanged()
{
  for( int e = 0; e < NumExponents; ++e )
    SinhCoshAdjust[e] = std::exp( - Energy[e] * HalfNt );
}

double Model2pt::operator()( int t ) const
{
  double z = 0;
  for( int e = 0; e < NumExponents; ++e )
  {
    double d;
    switch( parity )
    {
      case Common::Parity::Even:
        d = SinhCoshAdjust[e] * std::cosh( - Params[Energy[e]] * ( t - HalfNt ) );
        break;
      case Common::Parity::Odd:
        d = SinhCoshAdjust[e] * std::sinh( - Params[Energy[e]] * ( t - HalfNt ) );
        break;
      default:
        d = std::exp( - Params[Energy[e]] * ( t ) );
        break;
    }
    z += d * Params[src[e]] * Params[snk[e]];
  }
  /*if( !std::isfinite( z ) )
  {
    std::cout << "Energy[0]=" << Energy[0] << ", Energy[1]=" << Energy[1] << "\n";
    std::cout << "Coeff[0]=" << Coeff[0] << ", Coeff[1]=" << Coeff[1]
              << "Coeff[2]=" << Coeff[2] << ", Coeff[3]=" << Coeff[3] << "\n";
    std::cout << "Debug here\n";
  }*/
  return z;
}

// TODO: Make sure we add partial derivatives for all models
double Model2pt::Derivative( int t, int p ) const
{
  double d = 0;
  for( int e = 0; e < NumExponents; ++e )
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
  }
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
