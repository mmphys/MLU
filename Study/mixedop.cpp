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

#include <MLU/Common.hpp>

using scalar = double; // Minuit2 uses double, so haven't tested any other type
using degrees= int;    // Have tested float, but still need to fix file naming when number contains period
using Model  = Common::Model <scalar>;
using Fold   = Common::Fold  <scalar>;
using SSS    = Common::StartStopStep<degrees>;

static const std::string Sep{ "_" };
static const std::string Energy{ "E" };

// I refer to these as point and wall, but this works for any 2-operator basis rotation
static constexpr int idxPoint{ 0 };
static constexpr int idxWall{ 1 };
static constexpr int ModelNumIndices{ 2 }; // Model should include operators with these indices

static constexpr int idxSrc{ 0 };
static constexpr int idxSnk{ 1 };
static constexpr int NumSnkSrcIndices{ 2 }; // Number of sink and source indices

/*template <class scalar_ = scalar>
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
};*/

struct Parameters
{
  bool bSingleModel;
  std::string MixedOpName;
  int Exponent;
  bool bSaveCorr;
  bool bAllReplicas;
  bool bTryConjugate;
  std::string InBase;
  std::string OutBase;
  std::regex  RegExExt1;
  std::regex  RegExExt2;
  bool        RegExSwap;
  //int tmin;
  //int tmax;
  //int step;
  bool bNormalise; // Experimenting with normalisations
  bool bTheta;
  SSS Theta;
  bool bPhi;
  bool bPhiEqualsTheta;
  SSS Phi;
  std::array<std::array<scalar,ModelNumIndices>, NumSnkSrcIndices> OpNorm;
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
    int         idxNorm;
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
  //Normaliser<scalar> Normalisation;

  std::string NormString();
  void LoadModel( Common::FileNameAtt &&fna, const std::vector<std::vector<std::string>> &Words );
  void LoadCorrelator( int idxSnk, int idxSrc, const std::string &InFileName );
  void DoOneAngle( degrees Phi, degrees Theta, std::string &Out, std::size_t OutLen, bool bSaveCorr );
  inline bool IsSource() const { return CorrIn[0].size() > 1; }
  inline bool IsSink()   const { return CorrIn   .size() > 1; }
  virtual void LoadCorrelators( const std::string &InBase ) = 0;
  virtual void MixingAngle( degrees Phi, degrees Theta ) = 0;
  //virtual void SetNormalisation() { Normalisation.clear(); }
public:
  MixedOp( const std::string &description_, const Parameters &par ) : Description{ description_ }, Par{ par } {}
  virtual ~MixedOp() {}
  bool LoadModels( Common::FileNameAtt &&fna, const std::string &Args, bool * pSingle = nullptr );
  void Make( const std::string &FileName );
};

//#define NORMALISE( x ) std::sqrt( 2 * x )
#define NORMALISE( x ) std::sqrt( x )

std::string MixedOp::NormString()
{
  std::string s;
  s.append( 1, 'A' );
  s.append( std::to_string( Exponent ) );
  if( Par.bNormalise )
    s.append( 1, 'E' );
  return s;
}

class MixedOp_Norm : public MixedOp
{
//protected:
  //std::array<int, NumIndices> idxNorm;
public:
  using MixedOp::MixedOp;
  //void SetNormalisation( int idx );
};

// Mixed operator at source only
class MixedOp_Src : public MixedOp_Norm
{
protected:
  int iSnk; // Index of the sink operator (which is fixed and comes from the input file name)
  static const std::string MyDescription;
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle( degrees Phi, degrees Theta );
  //virtual void SetNormalisation() { MixedOp_Norm::SetNormalisation( ModelSrc ); }
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
  virtual void MixingAngle( degrees Phi, degrees Theta );
  //virtual void SetNormalisation() { MixedOp_Norm::SetNormalisation( ModelSnk ); }
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
  virtual void MixingAngle( degrees Phi, degrees Theta );
  //virtual void SetNormalisation() { MixedOp_Norm::SetNormalisation( ModelSrc ); }
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
  m.SetName( std::move( fna ) );
  m.Read( ( "Read model " + idxString + Common::Space ).c_str() );
  if( bFirst )
    Seed = m.Seed();
  else
    MI[0].model.IsCompatible( m, nullptr, Common::COMPAT_DISABLE_BASE );
  // Check that the exponent we want is available
  /*if( bFirst )
    Exponent = ( Par.Exponent >= 0 ) ? Par.Exponent : Par.Exponent + m.NumExponents;
  if( Exponent < 0 || Exponent >= m.NumExponents )
    throw std::runtime_error( "Exponent " + std::to_string( Par.Exponent ) + " not available in model" );
  // Display info about the model
  std::cout << "  Base: " << m.Name_.Base << Common::NewLine;
  for( std::size_t i = m.Name_.Extra.size(); i > 0; i-- )
    std::cout << "  Extra[" << std::to_string(i-1) << "]: " << m.Name_.Extra[i-1] << Common::NewLine;
  std::cout << "  Fit: ti=" << m.ti << ", tf=" << m.tf << Common::NewLine;
  for( std::size_t i = 0; i < m.OpNames.size(); i++ )
    std::cout << "  Model op[" << i << "]: " << m.OpNames[i] << Common::NewLine;
  // Now see which operators to use
  if( Words.empty() )
  {
    // Get operators from file
    if( m.OpNames.size() != ModelNumIndices )
      throw std::runtime_error( "Ambiguous - model file has " + std::to_string( m.OpNames.size() )
                               + " operators. " + std::to_string( ModelNumIndices ) + " expected." );
    for( int i = 0; i < ModelNumIndices; ++i )
    {
      mi.OpName[i] = m.OpNames[i];
      mi.OpIdx[i]  = m.GetColumnIndex( mi.OpName[i], Exponent );
    }
  }
  else
  {
    assert( Words.size() == ModelNumIndices && "Bug: Words.size() != ModelNumIndices" );
    int i = 0;
    for( const std::vector<std::string> &Word : Words )
    {
      mi.OpName[i] = Word[0];
      mi.OpIdx[i]  = m.GetColumnIndex( Word[Word.size() == 2 ? 1 : 0], Exponent );
      ++i;
    }
  }*/
  // Now say which operators we are using for normalisation
  for( int i = 0; i < ModelNumIndices; ++i )
    std::cout << "  normalisation[" << i << "]: " << mi.OpName[i] << Exponent << " = Model[" << mi.OpIdx[i] << "]" << Common::NewLine;
  // Now save the name of the mixed operator
  if( Par.MixedOpName.empty() )
  {
    for( const std::string &op : mi.OpName )
      mi.MixedOpName.append( op );
  }
  else
    mi.MixedOpName = Par.MixedOpName;
  // Now extract the index of the normalisation factor I'll need
  if( Par.bNormalise )
  {
    mi.idxNorm = m.GetColumnIndex( Energy, Exponent );
    std::cout << "  energy normalisation: " << Energy << Exponent << " = Model[" << mi.idxNorm << "]" << Common::NewLine;
  }
}

bool MixedOp::LoadModels( Common::FileNameAtt &&fna, const std::string &Args, bool * pSingle )
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
  if( Common::EqualIgnoreCase( fna.Ext, TEXT_EXT ) )
  {
    // Entries in the key file are relative to the file
    std::string Base{ std::move( fna.Dir ) };
    Common::AppendSlash( Base );
    const std::size_t BaseLen{ Base.length() };
    // Treat this as a list of fit files to read
    std::map<int, std::string> FitList;
    {
      std::ifstream s( fna.Filename );
      if( !Common::FileExists( fna.Filename ) || s.bad() )
        throw std::runtime_error( "Error reading \"" + fna.Filename + "\"" );
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
      if( pSingle )
        *pSingle = false;
    }
  }
  else if( Common::EqualIgnoreCase( fna.Type, Common::sModel ) )
  {
    LoadModel( std::move( fna ), Words );
    bIsModel = true;
    if( pSingle )
      *pSingle = true;
  }
  return bIsModel;
}

struct SinCos
{
  scalar cos, sin;
  SinCos( degrees Theta )
  {
    // ensure 0 <= Theta < 360
    bool bNeg{ Theta < 0 };
    if( bNeg )
      Theta = -Theta;
    if( Theta > 360 )
      Theta -= static_cast<int>( Theta / 360 ) * 360;
    if( bNeg )
      Theta = 360 - Theta;
    if( Theta == 0 )
    {
      cos = 1;
      sin = 0;
    }
    else if( Theta == 90 )
    {
      cos = 0;
      sin = 1;
    }
    else if( Theta == 180 )
    {
      cos = -1;
      sin = 0;
    }
    else if( Theta == 270 )
    {
      cos = 0;
      sin = -1;
    }
    else
    {
      const scalar ThetaRad{ M_PI * Theta / 180 };
      cos = std::cos( ThetaRad );
      sin = std::sin( ThetaRad );
    }
  }
};

void MixedOp::DoOneAngle( degrees Phi, degrees Theta, std::string &Out, std::size_t OutLen, bool bSaveCorr )
{
  //SinCos Phi( phi_ );
  //SinCos Theta( theta_ );
  //const SinCos &angle{ IsSource() ? Theta : Phi };
  MixingAngle( Phi, Theta );
  Out.resize( OutLen );
  Out.append( std::to_string( IsSink() ? Phi : Theta ) );
  if( IsSink() && IsSource() )
  {
    Out.append( "_theta_" );
    Out.append( std::to_string( Theta ) );
  }
  //if( Normalisation )
    //Out.append( "_N" );
  Out.append( 1, '.' );
  Out.append( IsSink()   ? MI[ModelSnk].MixedOpName : FileOpNames[CorrIn[0][0].Name_.op[idxSnk]] );
  Out.append( 1, '_' );
  Out.append( IsSource() ? MI[ModelSrc].MixedOpName : FileOpNames[CorrIn[0][0].Name_.op[idxSrc]] );
  CorrMixed.MakeCorrSummary();
  if( bSaveCorr )
    CorrMixed.Write( Common::MakeFilename( Out, CorrIn[0][0].Name_.Type, Seed, DEF_FMT ) );
  CorrMixed.WriteSummary( Common::MakeFilename( Out, CorrIn[0][0].Name_.Type, Seed, TEXT_EXT ) );
}

// Given a filename, work out which models apply, then rotate source and/or sink
void MixedOp::Make( const std::string &FileName )
{
  assert( ( IsSource() || IsSink() ) && "Bug: model must have at least one of sink and source" );
  assert( MI.size() && "No models loaded" );
  FileOpNames.clear();
  fnaName.Parse( FileName, &FileOpNames );
  if( fnaName.op.size() != 2 )
    throw std::runtime_error( "Could not extract sink and source operator names from " + FileName );
  if( Par.bSingleModel )
  {
    ModelSnk = 0;
    ModelSrc = 0;
  }
  else
  {
    std::smatch base_match;
    if( std::regex_search( fnaName.Base, base_match, Par.RegExExt2 ) && base_match.size() == 3 )
    {
      ModelSnk = Common::FromString<int>( base_match[Par.RegExSwap ? 2 : 1] );
      ModelSrc = Common::FromString<int>( base_match[Par.RegExSwap ? 1 : 2] );
    }
    else if( std::regex_search( fnaName.Base, base_match, Par.RegExExt1 ) && base_match.size() == 2 )
    {
      ModelSnk = Common::FromString<int>( base_match[1] );
      ModelSrc = ModelSnk;
    }
    else
      throw std::runtime_error( "Can't extract sink/source from " + FileName );
  }
  //if( Par.bNormalise )
    //SetNormalisation();
  //else
    //Normalisation.clear();

  //Input base is everything except the operator names
  std::string InBase{ fnaName.Dir };
  //Common::AppendSlash( InBase );
  InBase.append( fnaName.GetBaseExtra() );
  // Output base name for the output files is the input + fit times + theta
  std::string OutBase{ fnaName.GetBaseExtra() };
  if( IsSink() )
    MI[ModelSnk].model.Name_.AppendExtra( OutBase );
  if( IsSource() )
    MI[ModelSrc].model.Name_.AppendExtra( OutBase );
  // Add normalisation info to name
  {
    std::string Norm{ NormString() };
    if( !Norm.empty() )
    {
      OutBase.append( ".N_" );
      OutBase.append( Norm );
    }
  }
  std::cout << "Rotate " << OutBase << "\n";
  CorrMixed.FileList.clear();
  LoadCorrelators( InBase );
  CorrMixed.resize( NumSamples, Nt );
  CorrMixed.CopyAttributes( CorrIn[0][0] );

  std::string Out{ Par.OutBase + OutBase };
  std::cout << "    Writing " << Description << " mixed operator to:\n    " << Out << Common::NewLine << "   ";
  Out.append( IsSink() ? ".phi_" : ".theta_" );
  const std::size_t OutLen{ Out.length() };
  bool bFirstSnk{ true };
  for( degrees Phi : Par.Phi )
  {
    if( Par.bPhi && !Par.bPhiEqualsTheta && ( bFirstSnk || !( static_cast<int>( Phi ) % 10 ) ) )
    {
      std::cout << " " << std::to_string( Phi );
      std::cout.flush();
    }
    bFirstSnk = false;
    bool bFirstSrc{ true };
    for( degrees Theta : Par.Theta )
    {
      if( ( !Par.bPhi || Par.bPhiEqualsTheta ) && ( bFirstSrc || !( static_cast<int>( Theta ) % 10 ) ) )
      {
        std::cout << " " << std::to_string( Theta );
        std::cout.flush();
      }
      bFirstSrc = false;
      DoOneAngle( Par.bPhiEqualsTheta ? Theta : Phi, Theta, Out, OutLen, Par.bSaveCorr );
    }
  }
  std::cout << Common::NewLine;
}

// Set sink / source normalisation
/*void MixedOp_Norm::SetNormalisation( int idx )
{
  // This was the old crappy attempt to normalise to unit. Not really needed
  const ModelInfo & mi{ MI[idx] };
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
  }
}*/

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
    CorrMixed.NtUnfolded_ = Corr.NtUnfolded();
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

void MixedOp_Src::MixingAngle( degrees phi_, degrees theta_ )
{
  const SinCos Theta( theta_ );
  ModelInfo &miSrc{ MI[ModelSrc] };
  double Op_P{ 0 };
  double Op_W{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      double A_P { miSrc.model(i, miSrc.OpIdx[idxPoint]) * Par.OpNorm[idxSrc][idxPoint] };
      double A_W { miSrc.model(i, miSrc.OpIdx[idxWall ]) * Par.OpNorm[idxSrc][idxWall ] };
      if( Par.bNormalise )
      {
        const double Norm_E_Src{ NORMALISE( miSrc.model(i, miSrc.idxNorm) ) };
        A_P   *= Norm_E_Src;
        A_W   *= Norm_E_Src;
      }
      Op_P = Theta.cos / A_P;
      Op_W = Theta.sin / A_W;
    }
    for( int t = 0; t < Nt; t++ )
      CorrMixed(i,t) = Op_P * CorrIn[0][idxPoint](i,t) + Op_W * CorrIn[0][idxWall](i,t);
  }
}

void MixedOp_Snk::MixingAngle( degrees phi_, degrees theta_ )
{
  const SinCos Phi( phi_ );
  ModelInfo &miSnk{ MI[ModelSnk] };
  double Op_P{ 0 };
  double Op_W{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      double A_P { miSnk.model(i,miSnk.OpIdx[idxPoint]) * Par.OpNorm[idxSnk][idxPoint] };
      double A_W { miSnk.model(i,miSnk.OpIdx[idxWall ]) * Par.OpNorm[idxSnk][idxWall ] };
      if( Par.bNormalise )
      {
        const double Norm_E_Snk{ NORMALISE( miSnk.model(i,miSnk.idxNorm) ) };
        A_P   *= Norm_E_Snk;
        A_W   *= Norm_E_Snk;
      }
      Op_P = Phi.cos / A_P;
      Op_W = Phi.sin / A_W;
    }
    for( int t = 0; t < Nt; t++ )
      CorrMixed(i,t) = Op_P * CorrIn[idxPoint][0](i,t) + Op_W * CorrIn[idxWall ][0](i,t);
  }
}

void MixedOp_SnkSrc::MixingAngle( degrees phi_, degrees theta_ )
{
  const SinCos Phi( phi_ );
  const SinCos Theta( theta_ );
  const double AnglePP{ Phi.cos * Theta.cos };
  const double AnglePW{ Phi.cos * Theta.sin };
  const double AngleWP{ Phi.sin * Theta.cos };
  const double AngleWW{ Phi.sin * Theta.sin };
  ModelInfo &miSrc{ MI[ModelSrc] };
  ModelInfo &miSnk{ MI[ModelSnk] };
  const Fold &fPP{ CorrIn[idxPoint][idxPoint] };
  const Fold &fPW{ CorrIn[idxPoint][idxWall ] };
  const Fold &fWP{ CorrIn[idxWall ][idxPoint] };
  const Fold &fWW{ CorrIn[idxWall ][idxWall ] };
  double Op_PP{ 0 };
  double Op_WP{ 0 };
  double Op_PW{ 0 };
  double Op_WW{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      double A_P_snk { miSnk.model(i,miSnk.OpIdx[idxPoint]) * Par.OpNorm[idxSnk][idxPoint] };
      double A_W_snk { miSnk.model(i,miSnk.OpIdx[idxWall ]) * Par.OpNorm[idxSnk][idxWall ] };
      double A_P_src { miSrc.model(i,miSrc.OpIdx[idxPoint]) * Par.OpNorm[idxSrc][idxPoint] };
      double A_W_src { miSrc.model(i,miSrc.OpIdx[idxWall ]) * Par.OpNorm[idxSrc][idxWall ] };
      if( Par.bNormalise )
      {
        const scalar Norm_E_Snk{ NORMALISE( miSnk.model(i,miSnk.idxNorm) ) };
        const scalar Norm_E_Src{ NORMALISE( miSrc.model(i,miSrc.idxNorm) ) };
        A_P_snk *= Norm_E_Snk;
        A_W_snk *= Norm_E_Snk;
        A_P_src *= Norm_E_Src;
        A_W_src *= Norm_E_Src;
      }
      Op_PP = AnglePP / ( A_P_snk * A_P_src );
      Op_PW = AnglePW / ( A_P_snk * A_W_src );
      Op_WP = AngleWP / ( A_W_snk * A_P_src );
      Op_WW = AngleWW / ( A_W_snk * A_W_src );
      //std::cout << "\nCos theta=" << Theta.cos << ", Sin theta=" << Theta.sin << ", Cos phi=" << Phi.cos << ", Sin phi=" << Phi.sin << Common::NewLine;
      //std::cout << "A_P_src=" << A_P_src << ", A_W_src=" << A_W_src << "A_P_snk=" << A_P_snk << ", A_W_snk=" << A_W_src << Common::NewLine;
      //std::cout << "PP=" << fPP(i,0) << ", PW=" << fPW(i,0) << ", WP=" << fWP(i,0) << ", WW=" << fWW(i,0) << Common::NewLine;
      //std::cout << "M-M=" << ( Op_PP * fPP(i,0) + Op_WW * fWW(i,0) + Op_WP * fWP(i,0) + Op_PW * fPW(i,0) ) << Common::NewLine;
    }
    for( int t = 0; t < Nt; t++ )
      CorrMixed(i,t) = Op_PP * fPP(i,t) + Op_WW * fWW(i,t) + Op_WP * fWP(i,t) + Op_PW * fPW(i,t);
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
    fOut.NtUnfolded_ = fPW.NtUnfolded();
    fOut.parity = fPW.parity;
    fOut.reality = fPW.reality;
    fOut.sign = fPW.sign;
    fOut.Conjugated = false;
    fOut.t0Negated = false;
    for( int i = Fold::idxCentral; i < NumSamples; ++i )
    {
      for( int t = 0; t < Nt; ++t )
      {
        if( type == 0 )
          fOut(i,t) = fWP(i,t) / fPW(i,t);
        else
          fOut(i,t) = fWP(i,t) / fPW(i,DeltaT - t);
      }
    }
    static const std::array<std::string,2> Sufii{ "_Ratio", "_Ratio_reversed" };
    const std::string Out{ OutDir + fPW.Name_.GetBaseExtra() + Sufii[type] };
    static const Common::SeedType Seed{ fPW.Name_.Seed };
    std::cout << "Writing to " << Out << Common::NewLine;
    fOut.SetSummaryNames( "Ratio" );
    fOut.MakeCorrSummary();
    fOut.Write( Common::MakeFilename( Out, Common::sBootstrap, Seed, DEF_FMT ) );
    fOut.WriteSummary( Common::MakeFilename( Out, Common::sBootstrap, Seed, TEXT_EXT ), true );
  }
  return true;
}

std::string DoubleERE( const std::string &SingleERE )
{
  return SingleERE + "_" + SingleERE;
}

int main( int argc, const char *argv[] )
{
  static const char defaultExcited[] = "-1";
  static const std::string SingleERE{ R"([[:alpha:]]*([[:digit:]]+))" };
  static const char default_srcn0[] = "1.0";
  static const char default_srcn1[] = "1.0";
  static const char default_snkn0[] = "1.0";
  static const char default_snkn1[] = "1.0";
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"e", CL::SwitchType::Single, defaultExcited},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"m", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "" },
      {"r", CL::SwitchType::Single, SingleERE.c_str() },
      {"r2", CL::SwitchType::Single, nullptr },
      {"theta", CL::SwitchType::Single, nullptr},
      {"phi", CL::SwitchType::Single, nullptr},
      {"c", CL::SwitchType::Flag, nullptr},
      {"s", CL::SwitchType::Flag, nullptr},
      {"rep", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"enorm", CL::SwitchType::Flag, nullptr},
      {"srcn0", CL::SwitchType::Single, default_srcn0},
      {"srcn1", CL::SwitchType::Single, default_srcn1},
      {"snkn0", CL::SwitchType::Single, default_snkn0},
      {"snkn1", CL::SwitchType::Single, default_snkn1},
      {"not090", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      Parameters Par;
      Par.Exponent = cl.SwitchValue<int>("e");
      Par.InBase = cl.SwitchValue<std::string>("i");
      Par.OutBase = cl.SwitchValue<std::string>("o");
      Common::MakeAncestorDirs( Par.OutBase );
      const std::string modelBase{ cl.SwitchValue<std::string>("m") };
      Par.MixedOpName = cl.SwitchValue<std::string>("n");
      Par.bTryConjugate = cl.GotSwitch("c");
      {
        const std::string Single{ cl.SwitchValue<std::string>("r") };
        const std::string Double{ cl.GotSwitch( "r2" ) ? cl.SwitchValue<std::string>("r2") : DoubleERE( Single ) };
        Par.RegExExt1 = std::regex( Single, std::regex::extended | std::regex::icase );
        Par.RegExExt2 = std::regex( Double, std::regex::extended | std::regex::icase );
        if( Par.RegExExt1.mark_count() != 1 )
          throw std::runtime_error( "There should only be 1 marked sub-expression in \"" + Single + "\"" );
        if( Par.RegExExt2.mark_count() != 2 )
          throw std::runtime_error( "There should be 2 marked sub-expressions in \"" + Double + "\"" );
      }
      Par.RegExSwap = cl.GotSwitch("s");
      Par.bAllReplicas = cl.GotSwitch("rep");
      Par.bSaveCorr = cl.GotSwitch("savecorr");
      Par.bNormalise = cl.GotSwitch("enorm");
      Par.OpNorm[idxSrc][idxPoint] = cl.SwitchValue<scalar>("srcn0");
      Par.OpNorm[idxSrc][idxWall] = cl.SwitchValue<scalar>("srcn1");
      Par.OpNorm[idxSnk][idxPoint] = cl.SwitchValue<scalar>("snkn0");
      Par.OpNorm[idxSnk][idxWall] = cl.SwitchValue<scalar>("snkn1");
      const bool Do_0_90{ !cl.GotSwitch("not090") };
      Par.bTheta = cl.GotSwitch("theta");
      if( Par.bTheta )
      {
        Par.Theta = cl.SwitchValue<SSS>("theta");
        if( Do_0_90 )
          Par.Theta.AlwaysIterate = { 0, 90 };
      }
      //else
        //Par.Theta.SetSingle( 0 );
      Par.bPhi = cl.GotSwitch("phi");
      Par.bPhiEqualsTheta = false;
      //Par.Phi.SetSingle( 0 );
      if( Par.bPhi )
      {
        if( Common::EqualIgnoreCase( cl.SwitchValue<std::string>("phi"), "theta" ) )
          Par.bPhiEqualsTheta = true;
        else
        {
          Par.Phi = cl.SwitchValue<SSS>("phi");
          if( Do_0_90 )
            Par.Phi.AlwaysIterate = { 0, 90 };
        }
      }
      std::unique_ptr<MixedOp> Op;
      switch( ( Par.bTheta ? 1 : 0 ) + ( Par.bPhi ? 2 : 0 ) )
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
          throw std::runtime_error( "One or more of --theta and --phi must be specified" );
      }
      enum State : int { Empty, ModelLoaded, Rotation };
      State state{ State::Empty };
      std::size_t Count{ 0 };
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
            Common::FileNameAtt fna( ModelList[0] );
            if( Common::EqualIgnoreCase( fna.Ext, TEXT_EXT ) || Common::EqualIgnoreCase( fna.Type, Common::sModel ) )
            {
              if( state == State::ModelLoaded )
                throw std::runtime_error( "Reading consecutive models, i.e. no processing done with prior model" );
              if( Op->LoadModels( std::move( fna ), Args, &Par.bSingleModel ) )
              {
                bIsModel = true;
                state = State::ModelLoaded;
              }
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
            ++Count;
            bShowUsage = false;
          }
        }
      }
      std::cout << Count << " rotations performed" << Common::NewLine;
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
    "-e     Which excited state to use in model (default " << defaultExcited << "), -ve counts from highest\n"
    "-i     Input path for folded bootstrap replicas\n"
    "-o     Output path\n"
    "-m     Input path for model files\n"
    "-n     Mixed operator name\n"
    "-r     Extended regex to extract either sink or source type, default\n       " << SingleERE << "\n"
    "       http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap09.html\n"
    "--r2   Extended regex for both sink and source, default \"-r_-r\", i.e.\n       " << DoubleERE( SingleERE ) << "\n"
    "--theta Source rotation angles Min:Max[:Step]\n"
    "--phi   Sink   rotation angles Min:Max[:Step] | theta\n"
    "Flags:\n"
    "-c         Try loading Conjugate operator if bootstrap sample missing\n"
    "-s         Swap source / sink order in regex\n"
    "--rep      Use per replica values of overlap constants in construction of model\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--enorm    Energy normalisation, i.e. overlap coefficients *= sqrt(energy)\n"
    "--srcn0    Normalisation at source, operator 0, default " << default_srcn0 << "\n"
    "--srcn1    Normalisation at source, operator 1, default " << default_srcn1 << "\n"
    "--snkn0    Normalisation at sink,   operator 0, default " << default_snkn0 << "\n"
    "--snkn1    Normalisation at sink,   operator 1, default " << default_snkn1 << "\n"
    "--not090   Don't force 0 and 90 degrees\n"
    "--help     This message\n";
  }
  return iReturn;
}
