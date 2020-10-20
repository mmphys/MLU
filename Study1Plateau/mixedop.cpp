/*************************************************************************************
 
 Construct mixed operator from individual correlators and fit
 
 Source file: mixedop.cpp
 
 Copyright (C) 2019-2020
 
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

#include <cmath>
//#include <iomanip>
#include <stdio.h>

#include "../Analyze/Common.hpp"

using scalar = double;
using Model  = Common::Model <scalar>;
using Fold   = Common::Fold  <scalar>;

static const std::string Sep{ "_" };
static const std::string Energy{ "E" };

// I refer to these as point and wall, but this works for any 2-operator basis rotation
static constexpr int idxPoint{ 0 };
static constexpr int idxWall{ 1 };
static constexpr int ModelNumIndices{ 2 }; // Model should include operators with these indices

static constexpr int idxSrc{ 0 };
static constexpr int idxSnk{ 1 };
static constexpr int NumIndices{ 2 }; // Number of sink and source indices

template <class scalar_ = double>
class Normaliser
{
public:
  using Scalar = scalar_;
protected:
  bool   bGotNormalisation;
  Scalar normalisation; // Multiplicative constant
public:
  Normaliser()                       : bGotNormalisation{ false } {}
  Normaliser( Scalar Normalisation ) : bGotNormalisation{ true }, normalisation{ Normalisation } {}
  inline Scalar operator()( Scalar z ) const { return bGotNormalisation ? normalisation * z : z; }
  inline operator bool() const { return bGotNormalisation; }
  void clear() { bGotNormalisation = false; }
};

struct Parameters
{
  std::string MixedOpName;
  int Exponent;
  bool bSaveCorr;
  bool bAllReplicas;
  bool bTryConjugate;
  std::string InBase;
  std::string OutBase;
  std::regex  RegExExt;
  bool        RegExSwap;
  int tmin;
  int tmax;
  int step;
  bool bNormalise; // Experimenting with normalisations
};

class MixedOp
{
public:
  const std::string &Description;
  const Parameters &Par;
protected:
  struct ModelInfo
  {
    Model       model;
    std::string MixedOpName;
    std::array<int        , ModelNumIndices> OpIdx;
    std::array<std::string, ModelNumIndices> OpName;
  };
  std::vector<ModelInfo> MI;
  Common::SeedType Seed; // Set when first model loaded. Everything else needs to be compatible
  int ModelSnk, ModelSrc; // Based on filename, which model should be used at sink and source?
  int Exponent;
  std::vector<std::string> FileOpNames; // Operator names that came from the filenames ... somewhat redundant
  Common::FileNameAtt fnaName;
  Fold CorrMixed;
  std::vector<std::vector<Fold>> CorrIn;

  // Set by LoadCorrelators()
  int NumSamples; // maximum of a) user parameters b) raw correlators c) model (if doing all replicas)
  int Nt;
  Normaliser<scalar> Normalisation;

  void LoadModel( Common::FileNameAtt &&fna, const std::vector<std::vector<std::string>> &Words );
  void LoadCorrelator( int idxSnk, int idxSrc, const std::string &InFileName );
  void DoOneAngle( int degrees, std::string &Out, std::size_t OutLen, bool bSaveCorr, bool bPrint = true );
  inline bool IsSource() const { return CorrIn[0].size() > 1; }
  inline bool IsSink()   const { return CorrIn   .size() > 1; }
  virtual void LoadCorrelators( const std::string &InBase ) = 0;
  virtual void MixingAngle(double costheta, double sintheta) = 0;
  virtual void SetNormalisation() { Normalisation.clear(); }
public:
  MixedOp( const std::string &description_, const Parameters &par ) : Description{ description_ }, Par{ par } {}
  virtual ~MixedOp() {}
  bool LoadModels( const std::string &FileName, const std::string &Args );
  void Make( const std::string &FileName );
};

class MixedOp_Norm : public MixedOp
{
protected:
  std::array<int, NumIndices> idxNorm;
public:
  using MixedOp::MixedOp;
  void SetNormalisation( int idx );
};

// Mixed operator at source only
class MixedOp_Src : public MixedOp_Norm
{
protected:
  int iSnk; // Index of the sink operator (which is fixed and comes from the input file name)
  static const std::string MyDescription;
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
  virtual void SetNormalisation() { MixedOp_Norm::SetNormalisation( ModelSrc ); }
public:
  MixedOp_Src( const Parameters &par ) : MixedOp_Norm( MyDescription, par )
  {
    CorrIn.resize( 1 );
    CorrIn[0].resize( ModelNumIndices );
  }
};

// Mixed operator at sink only
class MixedOp_Snk : public MixedOp_Norm
{
protected:
  int iSrc; // Index of the source operator (which is fixed and comes from the input file name)
  static const std::string MyDescription;
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
  virtual void SetNormalisation() { MixedOp_Norm::SetNormalisation( ModelSnk ); }
public:
  MixedOp_Snk( const Parameters &par ) : MixedOp_Norm( MyDescription, par )
  {
    CorrIn.resize( ModelNumIndices );
    for( int i = 0; i < ModelNumIndices; ++i )
      CorrIn[i].resize( 1 );
  }
};

// Mixed operator at sink and source
class MixedOp_SnkSrc : public MixedOp_Norm
{
protected:
  static const std::string MyDescription;
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
  virtual void SetNormalisation() { MixedOp_Norm::SetNormalisation( ModelSrc ); }
public:
  MixedOp_SnkSrc( const Parameters &par ) : MixedOp_Norm( MyDescription, par )
  {
    CorrIn.resize( ModelNumIndices );
    for( int i = 0; i < ModelNumIndices; ++i )
      CorrIn[i].resize( ModelNumIndices );
  }
};

const std::string MixedOp_Src   ::MyDescription{ "source" };
const std::string MixedOp_Snk   ::MyDescription{ "sink" };
const std::string MixedOp_SnkSrc::MyDescription{ "source & sink" };

// Load one model given overlap coefficient names that have already been validated
void MixedOp::LoadModel( Common::FileNameAtt &&fna, const std::vector<std::vector<std::string>> &Words )
{
  // Now try to load the specified Model(s)
  const bool bFirst{ MI.empty() };
  const int idx{ static_cast<int>( MI.size() ) };
  const std::string idxString{ std::to_string( idx ) };
  MI.emplace_back();
  ModelInfo &mi{ MI.back() };
  Model &m{ mi.model };
  fna.ParseExtra();
  m.SetName( std::move( fna ) );
  m.Read( ( "Read model " + idxString + Common::Space ).c_str() );
  if( bFirst )
    Seed = m.Seed_;
  else
    MI[0].model.IsCompatible( m, nullptr, false );
  // Check that the exponent we want is available
  if( bFirst )
    Exponent = ( Par.Exponent >= 0 ) ? Par.Exponent : Par.Exponent + m.NumExponents;
  if( Exponent < 0 || Exponent >= m.NumExponents )
    throw std::runtime_error( "Exponent " + std::to_string( Par.Exponent ) + " not available in model" );
  // Display info about the model
  std::cout << "  Base: " << m.Name_.Base << Common::NewLine;
  for( std::size_t i = m.Name_.Extra.size(); i > 0; i-- )
    std::cout << "  Extra[" << std::to_string(i-1) << "]: " << m.Name_.Extra[i-1] << Common::NewLine;
  std::cout << "  Fit: " << ( m.Factorised ? "F" : "Unf" ) << "actorised, ti="
            << m.ti << ", tf=" << m.tf << Common::NewLine;
  for( std::size_t i = 0; i < m.OpNames.size(); i++ )
    std::cout << "  Model op[" << i << "]: " << m.OpNames[i] << Exponent << Common::NewLine;
  // Now see which operators to use
  if( Words.empty() )
  {
    // Get operators from file
    if( m.OpNames.size() != NumIndices )
      throw std::runtime_error( "Ambiguous - model file has " + std::to_string( m.OpNames.size() )
                               + " operators. " + std::to_string( NumIndices ) + " expected." );
    for( int i = 0; i < ModelNumIndices; ++i )
    {
      mi.OpName[i] = m.OpNames[i];
      mi.OpIdx[i]  = m.GetColumnIndex( mi.OpName[i], Exponent );
    }
  }
  else
  {
    int i = 0;
    for( const std::vector<std::string> &Word : Words )
    {
      mi.OpName[i] = Word[0];
      mi.OpIdx[i]  = m.GetColumnIndex( Word[Word.size() == 2 ? 1 : 0], Exponent );
      ++i;
    }
  }
  // Now save the name of the mixed operator
  if( Par.MixedOpName.empty() )
  {
    for( const std::string &op : mi.OpName )
      mi.MixedOpName.append( op );
  }
  else
    mi.MixedOpName = Par.MixedOpName;
}

bool MixedOp::LoadModels( const std::string &FileName, const std::string &Args )
{
  // First validate the arguments
  std::vector<std::vector<std::string>> Words;
  if( Args.find_first_not_of( Common::WhiteSpace ) != std::string::npos )
  {
    std::vector<std::string> OpNames = Common::ArrayFromString( Args );
    bool bOK{ ( OpNames.size() == ModelNumIndices ) };
    for( int i = 0; bOK && i < ModelNumIndices; ++i )
    {
       if( OpNames[i].empty() )
         bOK = false;
    }
    if( !bOK )
      throw std::runtime_error( "Either provide all " + std::to_string( ModelNumIndices ) + " operators or none" );
    Words.reserve( ModelNumIndices );
    for( int i = 0; i < ModelNumIndices; ++i )
    {
      Words.push_back( Common::ArrayFromString( OpNames[i], "=" ) );
      if( Words[i].empty() || Words[i][0].empty() || Words[i].size() > 2 || ( Words[i].size() > 1 && Words[i][1].empty() ) )
        throw std::runtime_error( "Cannot interpret model parameter " + std::to_string( i ) + " \"" + OpNames[i] + "\"" );
    }
  }
  MI.clear(); // Replace models already loaded
  bool bIsModel{ false };
  Common::FileNameAtt fna{ FileName };
  if( Common::EqualIgnoreCase( fna.Ext, TEXT_EXT ) )
  {
    // Entries in the key file are relative to the file
    std::string Base{ std::move( fna.Dir ) };
    Common::AppendSlash( Base );
    const std::size_t BaseLen{ Base.length() };
    // Treat this as a list of fit files to read
    std::map<int, std::string> FitList;
    {
      std::ifstream s( FileName );
      if( !Common::FileExists( FileName ) || s.bad() )
        throw std::runtime_error( "Error reading \"" + FileName + "\"" );
      FitList = Common::KeyValReader<int, std::string>::Read( s );
    }
    std::map<int, std::string>::iterator it;
    MI.reserve( FitList.size() );
    for( int iNum = 0; ( it = FitList.find( iNum ) ) != FitList.end() ; ++iNum )
    {
      Base.resize( BaseLen );
      Base.append( it->second );
      fna.Parse( Base );
      LoadModel( std::move( fna ), Words );
      bIsModel = true;
    }
  }
  else if( Common::EqualIgnoreCase( fna.Type, Common::sModel ) )
  {
    LoadModel( std::move( fna ), Words );
    bIsModel = true;
  }
  return bIsModel;
}

void MixedOp::DoOneAngle( int degrees, std::string &Out, std::size_t OutLen, bool bSaveCorr, bool bPrint )
{
  if( bPrint )
  {
    std::cout << " " << std::to_string( degrees );
    std::cout.flush();
  }
  double costheta, sintheta;
  const int degrees2 = degrees < 0 ? ( 360 - (-degrees) % 360 ) : degrees % 360;
  switch( degrees2 )
  {
    case 0:
      costheta = 1;
      sintheta = 0;
      break;
    case 90:
      costheta = 0;
      sintheta = 1;
      break;
    case 180:
      costheta = -1;
      sintheta = 0;
      break;
    case 270:
      costheta = 0;
      sintheta = -1;
      break;
    default:
    {
      const double theta{ M_PI * degrees2 / 180 };
      costheta = cos( theta );
      sintheta = sin( theta );
    }
  }
  MixingAngle( costheta, sintheta );
  Out.resize( OutLen );
  Out.append( std::to_string( degrees ) );
  if( Normalisation )
    Out.append( "_N" );
  Out.append( 1, '_' );
  Out.append( IsSink()   ? MI[ModelSnk].MixedOpName : FileOpNames[CorrIn[0][0].Name_.op[idxSnk]] );
  Out.append( 1, '_' );
  Out.append( IsSource() ? MI[ModelSrc].MixedOpName : FileOpNames[CorrIn[0][0].Name_.op[idxSrc]] );
  CorrMixed.MakeCorrSummary( nullptr );
  if( bSaveCorr )
    CorrMixed.Write( Common::MakeFilename( Out, CorrIn[0][0].Name_.Type, Seed, DEF_FMT ) );
  CorrMixed.WriteSummary( Common::MakeFilename( Out, CorrIn[0][0].Name_.Type, Seed, TEXT_EXT ) );
}

// Given a filename, work out which models apply, then rotate source and/or sink
void MixedOp::Make( const std::string &FileName )
{
  assert( ( IsSource() || IsSink() ) && "Bug: model must have at least one of sink and source" );
  FileOpNames.clear();
  fnaName.Parse( FileName, &FileOpNames );
  if( fnaName.op.size() != 2 )
    throw std::runtime_error( "Could not extract sink and source operator names from " + FileName );
  if( MI.size() <= 1 )
  {
    ModelSnk = 0;
    ModelSrc = 0;
  }
  else
  {
    std::smatch base_match;
    if( !std::regex_search( FileName, base_match, Par.RegExExt ) || base_match.size() != 3 )
      throw std::runtime_error( "Can't extract sink/source from " + FileName );
    ModelSnk = Common::FromString<int>( base_match[Par.RegExSwap ? 2 : 1] );
    ModelSrc = Common::FromString<int>( base_match[Par.RegExSwap ? 1 : 2] );
  }
  //if( Par.bNormalise )
    SetNormalisation();
  //else
    //Normalisation.clear();

  //Input base is everything except the operator names
  std::string InBase{ fnaName.Dir };
  Common::AppendSlash( InBase );
  InBase.append( fnaName.GetBaseExtra() );
  // Output base name for the output files is the input + fit times + theta
  std::string OutBase{ fnaName.GetBaseExtra() };
  if( IsSink() )
    MI[ModelSnk].model.Name_.AppendExtra( OutBase, 1 );
  if( IsSource() )
    MI[ModelSrc].model.Name_.AppendExtra( OutBase, 1 );
  std::cout << "Rotate " << OutBase << "\n";
  CorrMixed.FileList.clear();
  LoadCorrelators( InBase );
  CorrMixed.resize( NumSamples, Nt );
  CorrMixed.CopyAttributes( CorrIn[0][0] );

  std::string Out{ Par.OutBase + OutBase };
  std::cout << "    Writing " << Description << " mixed operator to:\n    " << Out << Common::NewLine << "   ";
  Out.append( ".theta_" );
  const std::size_t OutLen{ Out.length() };
  bool bDo0  = true; //!Par.bSaveCorr; // If we're not saving the correlators, always do 0 and 90 degrees
  bool bDo90 = bDo0;
  for( int degrees = Par.tmin; ; degrees+=Par.step )
  {
    DoOneAngle( degrees, Out, OutLen, Par.bSaveCorr, degrees == Par.tmin || !( degrees % 10 ) );
    switch( degrees )
    {
      case 0:
        bDo0 = false;
        break;
      case 90:
        bDo90 = false;
        break;
    }
    if( Par.step == 0 || degrees >= Par.tmax )
      break;
  }
  if( bDo0 )
    DoOneAngle( 0, Out, OutLen, false ); // We're outside the range requested, no need to save correlators
  if( bDo90 )
    DoOneAngle( 90, Out, OutLen, false );
  std::cout << Common::NewLine;
}

// Set sink / source normalisation
void MixedOp_Norm::SetNormalisation( int idx )
{
  // This was the old crappy attempt to normalise to unit. Not really needed
  /*const ModelInfo & mi{ MI[idx] };
  if( mi.model.NumExponents < 2 )
    Normalisation.clear();
  else
  {
    std::array<std::array<scalar, 2>, ModelNumIndices> MEL; // Two matrix elements for each operator
    for( int op = 0; op < ModelNumIndices; ++op )
    {
      for( int e = 0; e < 2; ++e )
        MEL[op][e] = mi.model.GetColumnIndex( mi.OpName[op], e );
    }
    scalar z =   MEL[idxPoint][1] * MEL[idxWall][1]
             / ( MEL[idxPoint][0] * MEL[idxWall][1] - MEL[idxPoint][1] * MEL[idxWall][0] );
    Normalisation = z * z;
  }*/
  if( Par.bNormalise )
  {
    for( int idx = 0; idx < NumIndices; ++idx )
    {
      idxNorm[idx] = MI[idx].model.GetColumnIndex( Energy, Exponent );
    }
  }
}

// Keep track of info on correlators as they are loaded
void MixedOp::LoadCorrelator( int idxSnk, int idxSrc, const std::string &InFileName )
{
  Fold &Corr{ CorrIn[idxSnk][idxSrc] };
  Corr.Read( InFileName, "  ", &FileOpNames );
  CorrMixed.FileList.push_back( InFileName );
  if( idxSnk == 0 && idxSrc == 0 )
  {
    Nt = Corr.Nt();
    if( NumSamples == 0 || NumSamples > Corr.NumSamples() )
      NumSamples = Corr.NumSamples();
    CorrMixed.NtUnfolded = Corr.NtUnfolded;
    CorrMixed.parity = Corr.parity;
    CorrMixed.reality = Corr.reality;
    CorrMixed.sign = Corr.sign;
    CorrMixed.Conjugated = false;
    CorrMixed.t0Negated = false;
  }
  else
  {
    CorrIn[0][0].IsCompatible( Corr, &NumSamples );
    if( CorrMixed.parity != Corr.parity )
      CorrMixed.parity = Common::Parity::Unknown;
    if( CorrMixed.reality != Corr.reality )
      CorrMixed.reality = Common::Reality::Unknown;
    if( CorrMixed.sign != Corr.sign )
      CorrMixed.sign = Common::Sign::Unknown;
  }
}

void MixedOp_Src::LoadCorrelators( const std::string &InBase )
{
  const std::string &SinkName{ FileOpNames[fnaName.op[idxSnk]] };
  iSnk = Common::IndexIgnoreCase( MI[ModelSnk].OpName, SinkName );
  if( iSnk == FileOpNames.size() )
    throw std::runtime_error( "Sink operator \"" + SinkName + "\" not member of rotation basis" );
  // If we're using all model replicas then this limits the maximum number of samples
  NumSamples = Par.bAllReplicas ? MI[ModelSrc].model.NumSamples() : 0;
  // Load all of the correlators
  std::string InFile{ InBase };
  InFile.append( 1, '_' );
  const std::size_t InFileLen1{ InFile.length() };
  InFile.append( SinkName ); // Use the sink from the original file
  InFile.append( 1, '_' );
  const std::size_t InFileLen2{ InFile.length() };
  for( int iSrc = 0; iSrc < ModelNumIndices; iSrc++ )
  {
    InFile.resize( InFileLen2 );
    const std::string &SourceName{ MI[ModelSrc].OpName[iSrc] };
    InFile.append( SourceName );
    std::string InFileName{ Common::MakeFilename( InFile, Common::sFold, Seed, DEF_FMT ) };
    if(!Common::FileExists( InFileName ) && Par.bTryConjugate && !Common::EqualIgnoreCase( SinkName, SourceName ) )
    {
      std::cout << "Warning: loading conjugate of " << InFileName << Common::NewLine;
      InFile.resize( InFileLen1 );
      InFile.append( SourceName );
      InFile.append( 1, '_' );
      InFile.append( SinkName );
      InFileName = Common::MakeFilename( InFile, Common::sFold, Seed, DEF_FMT );
    }
    LoadCorrelator( 0, iSrc, InFileName );
  }
}

void MixedOp_Snk::LoadCorrelators( const std::string &InBase )
{
  const std::string &SourceName{ FileOpNames[fnaName.op[idxSrc]] };
  iSrc = Common::IndexIgnoreCase( MI[ModelSrc].OpName, SourceName );
  if( iSrc == FileOpNames.size() )
    throw std::runtime_error( "Source operator \"" + SourceName + "\" not member of rotation basis" );
  // If we're using all model replicas then this limits the maximum number of samples
  NumSamples = Par.bAllReplicas ? MI[ModelSrc].model.NumSamples() : 0;
  // Load all of the correlators
  std::string InFile{ InBase };
  InFile.append( 1, '_' );
  const std::size_t InFileLen1{ InFile.length() };
  for( int iSnk = 0; iSnk < ModelNumIndices; iSnk++ )
  {
    InFile.resize( InFileLen1 );
    const std::string &SinkName{ MI[ModelSrc].OpName[iSnk] };
    InFile.append( SinkName );
    InFile.append( 1, '_' );
    InFile.append( SourceName ); // Use the source from the original file
    std::string InFileName{ Common::MakeFilename( InFile, Common::sFold, Seed, DEF_FMT ) };
    if(!Common::FileExists( InFileName ) && Par.bTryConjugate && !Common::EqualIgnoreCase( SinkName, SourceName ) )
    {
      std::cout << "Warning: loading conjugate of " << InFileName << Common::NewLine;
      InFile.resize( InFileLen1 );
      InFile.append( SourceName );
      InFile.append( 1, '_' );
      InFile.append( SinkName );
      InFileName = Common::MakeFilename( InFile, Common::sFold, Seed, DEF_FMT );
    }
    LoadCorrelator( iSnk, 0, InFileName );
  }
}

void MixedOp_SnkSrc::LoadCorrelators( const std::string &InBase )
{
  // If we're using all model replicas then this limits the maximum number of samples
  NumSamples = Par.bAllReplicas ? std::min( MI[ModelSnk].model.NumSamples(), MI[ModelSrc].model.NumSamples() ) : 0;
  // Load all of the correlators
  std::string InFile{ InBase };
  InFile.append( 1, '_' );
  const std::size_t InFileLen1{ InFile.length() };
  for( int iSnk = 0; iSnk < ModelNumIndices; iSnk++ )
  {
    InFile.resize( InFileLen1 );
    const std::string &SinkName{ MI[ModelSrc].OpName[iSnk] };
    InFile.append( SinkName ); // Use the sink from the original file
    InFile.append( 1, '_' );
    const std::size_t InFileLen2{ InFile.length() };
    for( int iSrc = 0; iSrc < ModelNumIndices; iSrc++ )
    {
      InFile.resize( InFileLen2 );
      const std::string &SourceName{ MI[ModelSrc].OpName[iSrc] };
      InFile.append( SourceName );
      std::string InFileName{ Common::MakeFilename( InFile, Common::sFold, Seed, DEF_FMT ) };
      if(!Common::FileExists( InFileName ) && Par.bTryConjugate && !Common::EqualIgnoreCase( SinkName, SourceName ) )
      {
        std::cout << "Warning: loading conjugate of " << InFileName << Common::NewLine;
        InFile.resize( InFileLen1 );
        InFile.append( SourceName );
        InFile.append( 1, '_' );
        InFile.append( SinkName );
        InFileName = Common::MakeFilename( InFile, Common::sFold, Seed, DEF_FMT );
      }
      LoadCorrelator( iSnk, iSrc, InFileName );
    }
  }
}

//#define NORMALISE( x ) std::sqrt( 2 * x )
#define NORMALISE( x ) std::sqrt( x )

void MixedOp_Src::MixingAngle(double costheta, double sintheta)
{
  const double * pModelSrc{ MI[ModelSrc].model[Model::idxCentral] };
  const double * pModelSnk{ MI[ModelSnk].model[Model::idxCentral] };
  const double * pDataP{ CorrIn[0][idxPoint][Fold::idxCentral] };
  const double * pDataW{ CorrIn[0][idxWall ][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_Snk{ 0 };
  double A_P{ 0 };
  double A_W{ 0 };
  double Op_P{ 0 };
  double Op_W{ 0 };
  //double E_Norm{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_Snk = pModelSnk[MI[ModelSnk].OpIdx[iSnk]];
      A_P   = pModelSrc[MI[ModelSrc].OpIdx[idxPoint]];
      A_W   = pModelSrc[MI[ModelSrc].OpIdx[idxWall]];
      if( Par.bNormalise )
      {
        const double E_Snk{ pModelSnk[idxNorm[idxSnk]] };
        const double E_Src{ pModelSrc[idxNorm[idxSrc]] };
        //E_Norm = 1 / ( 4 * E_Snk * E_Src );
        const scalar n{ NORMALISE( E_Src ) };
        A_Snk *= NORMALISE( E_Snk );
        A_P   *= n;
        A_W   *= n;
      }
      Op_P = costheta / ( A_Snk * A_P );
      Op_W = sintheta / ( A_Snk * A_W );
      pModelSrc += MI[ModelSrc].model.Nt();
      pModelSnk += MI[ModelSnk].model.Nt();
    }
    for( int t = 0; t < Nt; t++ )
    {
      scalar z{ Op_P * *pDataP++ + Op_W * *pDataW++ };
      // if( Par.bNormalise )
        // z *= E_Norm;
      *pDst++ = z;
    }
  }
}

void MixedOp_Snk::MixingAngle(double costheta, double sintheta)
{
  const double * pModelSrc{ MI[ModelSrc].model[Model::idxCentral] };
  const double * pModelSnk{ MI[ModelSnk].model[Model::idxCentral] };
  const double * pDataP{ CorrIn[idxPoint][0][Fold::idxCentral] };
  const double * pDataW{ CorrIn[idxWall ][0][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_Src{ 0 };
  double A_P{ 0 };
  double A_W{ 0 };
  double Op_P{ 0 };
  double Op_W{ 0 };
  //double E_Norm{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_Src = pModelSrc[MI[ModelSrc].OpIdx[iSrc]];
      A_P   = pModelSnk[MI[ModelSnk].OpIdx[idxPoint]];
      A_W   = pModelSnk[MI[ModelSnk].OpIdx[idxWall]];
      if( Par.bNormalise )
      {
        const double E_Snk{ pModelSnk[idxNorm[idxSnk]] };
        const double E_Src{ pModelSrc[idxNorm[idxSrc]] };
        //E_Norm = 1 / ( 4 * E_Snk * E_Src );
        const scalar n{ NORMALISE( E_Snk ) };
        A_Src *= NORMALISE( E_Src );
        A_P   *= n;
        A_W   *= n;
      }
      Op_P = costheta / ( A_Src * A_P );
      Op_W = sintheta / ( A_Src * A_W );
      pModelSrc += MI[ModelSrc].model.Nt();
      pModelSnk += MI[ModelSnk].model.Nt();
    }
    for( int t = 0; t < Nt; t++ )
    {
      scalar z{ Op_P * *pDataP++ + Op_W * *pDataW++ };
      // if( Par.bNormalise )
        // z *= E_Norm;
      *pDst++ = z;
    }
  }
}

void MixedOp_SnkSrc::MixingAngle(double costheta, double sintheta)
{
  const double cos_sq_theta{ costheta * costheta };
  const double sin_sq_theta{ sintheta * sintheta };
  const double cos_sin_theta{ costheta * sintheta };
  const double * pModelSrc{ MI[ModelSrc].model[Model::idxCentral] };
  const double * pModelSnk{ MI[ModelSnk].model[Model::idxCentral] };
  const double * pPP{ CorrIn[idxPoint][idxPoint][Fold::idxCentral] };
  const double * pPW{ CorrIn[idxPoint][idxWall ][Fold::idxCentral] };
  const double * pWP{ CorrIn[idxWall ][idxPoint][Fold::idxCentral] };
  const double * pWW{ CorrIn[idxWall ][idxWall ][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_P_snk{ 0 };
  double A_P_src{ 0 };
  double A_W_snk{ 0 };
  double A_W_src{ 0 };
  double Op_PP{ 0 };
  double Op_WP{ 0 };
  double Op_PW{ 0 };
  double Op_WW{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_P_snk = pModelSnk[MI[ModelSnk].OpIdx[idxPoint]];
      A_P_src = pModelSrc[MI[ModelSrc].OpIdx[idxPoint]];
      A_W_snk = pModelSnk[MI[ModelSnk].OpIdx[idxWall ]];
      A_W_src = pModelSrc[MI[ModelSrc].OpIdx[idxWall ]];
      if( Par.bNormalise )
      {
        const double E_Snk{ pModelSnk[idxNorm[idxSnk]] };
        const double E_Src{ pModelSrc[idxNorm[idxSrc]] };
        const scalar nSnk{ NORMALISE( E_Snk ) };
        const scalar nSrc{ NORMALISE( E_Src ) };
        // const scalar nSnk{ 2 * E_Snk };
        // const scalar nSrc{ 2 * E_Src };
        A_P_snk *= nSnk;
        A_P_src *= nSrc;
        A_W_snk *= nSnk;
        A_W_src *= nSrc;
      }
      Op_PP = cos_sq_theta  / ( A_P_snk * A_P_src );
      Op_WP = cos_sin_theta / ( A_W_snk * A_P_src );
      Op_PW = cos_sin_theta / ( A_P_snk * A_W_src );
      Op_WW = sin_sq_theta  / ( A_W_snk * A_W_src );
      pModelSrc += MI[ModelSrc].model.Nt();
      pModelSnk += MI[ModelSnk].model.Nt();
    }
    for( int t = 0; t < Nt; t++ )
      *pDst++ = Normalisation( Op_PP * *pPP++ + Op_WW * *pWW++ +     Op_WP * *pWP++ + Op_PW * *pPW++ );
  }
}

// Compare normalisation of PW vs WP
bool Debug()
{
  static const char szReading[] = ( "Reading " );
  static const std::string InDir{ "/Users/s1786208/data/Study2/C1/PW/3pt_s/" };
  static const std::string InPW{ InDir + "quark_h0_h0_gT_dt_20_p_0_0_0_g5P_g5W.fold.1835672416.h5" };
  static const std::string InWP{ InDir + "quark_h0_h0_gT_dt_20_p_0_0_0_g5W_g5P.fold.1835672416.h5" };
  static const std::string OutDir{ "/Users/s1786208/data/Study2/C1/PW/Debug/" };
  std::vector<std::string> OpNames;
  Fold fPW( InPW, szReading, &OpNames );
  Fold fWP( InWP, szReading, &OpNames );
  int NumSamples = 0;
  fPW.IsCompatible( fWP, &NumSamples );
  if( !fPW.Name_.bGotDeltaT )
    throw std::runtime_error( "Could not etract DeltaT from" + InPW );
  const int DeltaT{ fPW.Name_.DeltaT };
  for( int type = 0; type < 2; ++type )
  {
    const int Nt{ type == 0 ? fPW.Nt() : DeltaT + 1 };
    Fold fOut( NumSamples, Nt );
    fOut.FileList.push_back( InPW );
    fOut.FileList.push_back( InWP );
    fOut.CopyAttributes( fPW );
    fOut.NtUnfolded = fPW.NtUnfolded;
    fOut.parity = fPW.parity;
    fOut.reality = fPW.reality;
    fOut.sign = fPW.sign;
    fOut.Conjugated = false;
    fOut.t0Negated = false;
    const scalar * pPW{ fPW[Fold::idxCentral] };
    const scalar * pWP{ fWP[Fold::idxCentral] };
    scalar * pOut{ fOut[Fold::idxCentral] };
    for( int i = Fold::idxCentral; i < NumSamples; ++i )
    {
      for( int t = 0; t < Nt; ++t )
      {
        if( type == 0 )
          *pOut++ = pWP[t] / pPW[t];
        else
          *pOut++ = pWP[t] / pPW[DeltaT - t];
      }
      pWP += fWP.Nt();
      pPW += fPW.Nt();
    }
    static const std::array<std::string,2> Sufii{ "_Ratio", "_Ratio_reversed" };
    const std::string Out{ OutDir + fPW.Name_.GetBaseExtra() + Sufii[type] };
    static const Common::SeedType Seed{ fPW.Name_.Seed };
    std::cout << "Writing to " << Out << Common::NewLine;
    fOut.MakeCorrSummary( "Ratio" );
    fOut.Write( Common::MakeFilename( Out, Common::sBootstrap, Seed, DEF_FMT ) );
    fOut.WriteSummary( Common::MakeFilename( Out, Common::sBootstrap, Seed, TEXT_EXT ), true );
  }
  return true;
}

int main( int argc, const char *argv[] )
{
  static const char defaultModel[] = "1";
  static const char defaultExcited[] = "-1";
  static const char defaultTMin[] = "-90";
  static const char defaultTMax[] = "90";
  static const char defaultStep[] = "1";
  static const char DefaultERE[]{ R"([_[:alpha:]]*([[:digit:]]+)_[[:alpha:]]*([[:digit:]]+))" };
  std::ios_base::sync_with_stdio( false );
  //if( Debug() ) return 0;
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"a", CL::SwitchType::Single, defaultModel},
      {"e", CL::SwitchType::Single, defaultExcited},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"m", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "" },
      {"r", CL::SwitchType::Single, DefaultERE },
      {"tmin", CL::SwitchType::Single, defaultTMin},
      {"tmax", CL::SwitchType::Single, defaultTMax},
      {"step", CL::SwitchType::Single, defaultStep},
      {"c", CL::SwitchType::Flag, nullptr},
      {"s", CL::SwitchType::Flag, nullptr},
      {"rep", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"normalise", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      Parameters Par;
      Par.Exponent = cl.SwitchValue<int>("e");
      Par.InBase = Common::AppendSlash( cl.SwitchValue<std::string>("i") );
      Par.OutBase = Common::AppendSlash( cl.SwitchValue<std::string>("o") );
      const std::string modelBase{ cl.SwitchValue<std::string>("m") };
      Par.MixedOpName = cl.SwitchValue<std::string>("n");
      Par.tmin = cl.SwitchValue<int>("tmin");
      Par.tmax = cl.SwitchValue<int>("tmax");
      Par.step = std::abs( cl.SwitchValue<int>("step") );
      Par.bTryConjugate = cl.GotSwitch("c");
      Par.RegExExt = std::regex( cl.SwitchValue<std::string>("r"),std::regex::extended|std::regex::icase);
      Par.RegExSwap = cl.GotSwitch("s");
      Par.bAllReplicas = cl.GotSwitch("rep");
      Par.bSaveCorr = cl.GotSwitch("savecorr");
      Par.bNormalise = cl.GotSwitch("normalise");
      if( Par.tmin > Par.tmax )
        Par.step *= -1;
      std::unique_ptr<MixedOp> Op;
      switch( cl.SwitchValue<int>("a") )
      {
        case 1:
          Op.reset( new MixedOp_Src( Par ) );
          break;
        case 2:
          Op.reset( new MixedOp_Snk( Par ) );
          break;
        case 3:
          Op.reset( new MixedOp_SnkSrc( Par ) );
          break;
        default:
          throw std::runtime_error( "Invalid source/sink type " + cl.SwitchValue<std::string>("a") );
      }
      enum State : int { Empty, ModelLoaded, Rotation };
      State state{ State::Empty };
      for( const std::string &FileName : cl.Args )
      {
        // See whether this is a model
        bool bIsModel{ false };
        {
          std::string Args{ FileName };
          std::string ModelFileName{ Common::ExtractToSeparator( Args ) };
          std::vector<std::string> ModelList{ Common::glob( &ModelFileName, &ModelFileName + 1, modelBase.c_str() ) };
          if( ModelList.size() == 1 && Common::FileExists( ModelList[0] ) )
          {
            if( state == State::ModelLoaded )
              throw std::runtime_error( "Reading consecutive models, i.e. no processing done with prior model" );
            if( Op->LoadModels( ModelList[0], Args ) )
            {
              bIsModel = true;
              state = State::ModelLoaded;
            }
          }
          if( !bIsModel )
          {
            if( state == State::Empty )
              throw std::runtime_error( ModelFileName + " is not a model file" );
            if( !Args.empty() )
              throw std::runtime_error( "Argument \"" + Args + "\" specified for " + ModelFileName );
          }
        }
        if( !bIsModel )
        {
          // This is not a model - see whether it's a number of rotations to perform
          if( state == State::Empty )
            throw std::runtime_error( "A model must be loaded before rotations performed" );
          for( const std::string &BootFile : Common::glob( &FileName, &FileName + 1, Par.InBase.c_str() ) )
          {
            Op->Make( BootFile );
            bShowUsage = false;
          }
        }
      }
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << cl.Name <<
    " <options> Model Bootstrap_wildcard [Model Bootstrap_wildcard [...]]\n"
    "Create a mixed operator from fit parameters and bootstrap replicas\n"
    "<options>\n"
    "-a     Model type: 1=source only; 2=sink only; 3=source/sink; default " << defaultModel << "\n"
    "-e     Which excited state to use in model (default " << defaultExcited << "), -ve counts from highest\n"
    "-i     Input path for folded bootstrap replicas\n"
    "-o     Output path\n"
    "-m     Input path for model files\n"
    "-n     Mixed operator name\n"
    "-r     Extended regex for sink/source type, default\n       " << DefaultERE << "\n"
    "       http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap09.html\n"
    "--tmin Minimum theta, default " << defaultTMin << "\n"
    "--tmax Maximum theta, default " << defaultTMax << "\n"
    "--step Steps between theta, default " << defaultStep << "\n"
    "Flags:\n"
    "-c         Try loading Conjugate operator if bootstrap sample missing\n"
    "-s         Swap source / sink order in regex\n"
    "--rep      Use per replica values of overlap constants in construction of model\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--normalise Normalise mixed correlator to unity (experimental)\n"
    "--help     This message\n";
  }
  return iReturn;
}
