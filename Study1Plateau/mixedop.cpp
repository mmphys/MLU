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
static const std::string NewLine{ "\n" };

// The order here must match the algebra array for the model code to work
static constexpr int idxg5{ 0 };
static constexpr int idxgT5{ 1 };
static constexpr int ModelNumIndices{ 2 }; // Model should include operators with these indices

static constexpr int idxSrc{ 0 };
static constexpr int idxSnk{ 1 };
static constexpr int NumIndices{ 2 }; // Number of sink and source indices

struct SinkSource
{
  int Sink;
  int Source;
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

  void LoadModel( Common::FileNameAtt &&fna, const std::vector<std::vector<std::string>> &Words );
  void LoadCorrelator( int idxSnk, int idxSrc, const std::string &InFileName );
  void DoOneAngle( int degrees, std::string &Out, std::size_t OutLen, bool bSaveCorr, bool bPrint = true );
  virtual bool IsSource() { return true; }
  virtual bool IsSink() { return false; }
  virtual void LoadCorrelators( const std::string &InBase ) = 0;
  virtual void MixingAngle(double costheta, double sintheta) = 0;
public:
  MixedOp( const std::string &description_, const Parameters &par ) : Description{ description_ }, Par{ par } {}
  virtual ~MixedOp() {}
  bool LoadModels( const std::string &FileName, const std::string &Args );
  void Make( const std::string &FileName );
};

// Mixed operator at source only
class MixedOp_S : public MixedOp
{
protected:
  static const std::string MyDescription;
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
public:
  MixedOp_S( const Parameters &par ) : MixedOp( MyDescription, par )
  {
    CorrIn.resize( 1 );
    CorrIn[0].resize( ModelNumIndices );
  }
};

// Mixed operator at sink and source
class MixedOp_SS : public MixedOp
{
protected:
  static const std::string MyDescription;
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
  virtual bool IsSink() { return true; }
public:
  MixedOp_SS( const Parameters &par ) : MixedOp( MyDescription, par )
  {
    CorrIn.resize( ModelNumIndices );
    for( int i = 0; i < ModelNumIndices; ++i )
      CorrIn[i].resize( ModelNumIndices );
  }
};

const std::string MixedOp_S ::MyDescription{ "source" };
const std::string MixedOp_SS::MyDescription{ "source & sink" };

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
  // Display info about the model
  std::cout << "  Base: " << m.Name_.Base << NewLine;
  for( std::size_t i = m.Name_.Extra.size(); i > 0; i-- )
    std::cout << "  Extra[" << std::to_string(i-1) << "]: " << m.Name_.Extra[i-1] << NewLine;
  std::cout << "  Fit: " << ( m.Factorised ? "F" : "Unf" ) << "actorised, ti="
            << m.ti << ", tf=" << m.tf << NewLine;
  for( std::size_t i = 0; i < m.OpNames.size(); i++ )
    std::cout << "  Model op[" << i << "]: " << m.OpNames[i] << Exponent << NewLine;
  // Check that the exponent we want is available
  if( bFirst )
    Exponent = ( Par.Exponent >= 0 ) ? Par.Exponent : Par.Exponent + m.NumExponents;
  if( Exponent < 0 || Exponent >= m.NumExponents )
    throw std::runtime_error( "Exponent " + std::to_string( Par.Exponent ) + " not available in model" );
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
  Out.append( 1, '_' );
  Out.append( CorrIn   .size() == 1 ? FileOpNames[CorrIn[0][0].Name_.op[idxSnk]] : MI[ModelSrc].MixedOpName );
  Out.append( 1, '_' );
  Out.append( CorrIn[0].size() == 1 ? FileOpNames[CorrIn[0][0].Name_.op[idxSrc]] : MI[ModelSrc].MixedOpName );
  CorrMixed.MakeCorrSummary( nullptr );
  if( bSaveCorr )
    CorrMixed.Write( Common::MakeFilename( Out, Common::sBootstrap, Seed, DEF_FMT ) );
  CorrMixed.WriteSummary( Common::MakeFilename( Out, Common::sBootstrap, Seed, TEXT_EXT ) );
}

// Given a filename, work out which models apply, then rotate source and/or sink
void MixedOp::Make( const std::string &FileName )
{
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
  std::cout << NewLine;
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

void MixedOp_S::LoadCorrelators( const std::string &InBase )
{
  // If we're using all model replicas then this limits the maximum number of samples
  NumSamples = Par.bAllReplicas ? MI[ModelSrc].model.Nt() : 0;
  // Load all of the correlators
  std::string InFile{ InBase };
  InFile.append( 1, '_' );
  const std::size_t InFileLen1{ InFile.length() };
  const std::string &SinkName{ FileOpNames[fnaName.op[idxSnk]] };
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

void MixedOp_SS::LoadCorrelators( const std::string &InBase )
{
  // If we're using all model replicas then this limits the maximum number of samples
  NumSamples = Par.bAllReplicas ? std::min( MI[ModelSnk].model.Nt(), MI[ModelSrc].model.Nt() ) : 0;
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

void MixedOp_S::MixingAngle(double costheta, double sintheta)
{
  const ModelInfo &mi{ MI[ModelSrc] };
  const double * pCoeff{ mi.model[Model::idxCentral] };
  const double * pPP{ CorrIn[0][idxSrc][Fold::idxCentral] };
  const double * pPW{ CorrIn[0][idxSnk][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_PP{ 0 };
  double A_PW{ 0 };
  double Op_PP{ 0 };
  double Op_PW{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_PP = pCoeff[mi.OpIdx[idxSrc]];
      A_PW = pCoeff[mi.OpIdx[idxSnk]];
      Op_PP = costheta / ( A_PP * A_PP );
      Op_PW = sintheta / ( A_PP * A_PW );
      pCoeff += mi.model.Nt();
    }
    for( int t = 0; t < Nt; t++ )
      *pDst++ = Op_PP * *pPP++ + Op_PW * *pPW++;
  }
}

void MixedOp_SS::MixingAngle(double costheta, double sintheta)
{
  const double cos_sq_theta{ costheta * costheta };
  const double sin_sq_theta{ sintheta * sintheta };
  const double cos_sin_theta{ costheta * sintheta };
  const double * pCoeffSrc{ MI[ModelSrc].model[Model::idxCentral] };
  const double * pCoeffSnk{ MI[ModelSnk].model[Model::idxCentral] };
  const double * pPP{ CorrIn[0][0][Fold::idxCentral] };
  const double * pPA{ CorrIn[0][1][Fold::idxCentral] };
  const double * pAP{ CorrIn[1][0][Fold::idxCentral] };
  const double * pAA{ CorrIn[1][1][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_P1_snk{ 0 };
  double A_P1_src{ 0 };
  double A_A1_snk{ 0 };
  double A_A1_src{ 0 };
  double Op_PP{ 0 };
  double Op_AP{ 0 };
  double Op_PA{ 0 };
  double Op_AA{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_P1_snk = pCoeffSnk[MI[ModelSnk].OpIdx[idxSrc]];
      A_P1_src = pCoeffSrc[MI[ModelSrc].OpIdx[idxSrc]];
      A_A1_snk = pCoeffSnk[MI[ModelSnk].OpIdx[idxSnk]];
      A_A1_src = pCoeffSrc[MI[ModelSrc].OpIdx[idxSnk]];
      Op_PP = cos_sq_theta  / ( A_P1_snk * A_P1_src );
      Op_AP = cos_sin_theta / ( A_A1_snk * A_P1_src );
      Op_PA = cos_sin_theta / ( A_P1_snk * A_A1_src );
      Op_AA = sin_sq_theta  / ( A_A1_snk * A_A1_src );
      pCoeffSrc += MI[ModelSrc].model.Nt();
      pCoeffSnk += MI[ModelSnk].model.Nt();
    }
    for( int t = 0; t < Nt; t++ )
      *pDst++ = Op_PP * *pPP++ + Op_AA * *pAA++ +     Op_AP * *pAP++ + Op_PA * *pPA++;
  }
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
      if( Par.tmin > Par.tmax )
        Par.step *= -1;
      std::unique_ptr<MixedOp> Op;
      switch( cl.SwitchValue<int>("a") )
      {
        case 1:
          Op.reset( new MixedOp_S( Par ) );
          break;
        case 3:
          Op.reset( new MixedOp_SS( Par ) );
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
    "-a     Model type: 1=sink only; 2=source only; 3=source/sink; default " << defaultModel << "\n"
    "-e     Which excited state to use in model (default " << defaultExcited << "), -ve counts from highest state\n"
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
    "--help     This message\n";
  }
  return iReturn;
}
