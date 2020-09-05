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
#include <ios>
#include <iostream>
#include <stdio.h>

#include "../Analyze/Common.hpp"

using scalar = double;
using Model  = Common::Model <scalar>;
using Fold   = Common::Fold  <scalar>;

static const std::string Sep{ "_" };
static const std::string NewLine{ "\n" };

/*struct OperatorTrait
{
  //int Imag;
  //int Parity;
  Common::Gamma::Algebra Alg;
  //OperatorTrait( int Imag_, int Parity_, Common::Gamma::Algebra Alg_ ) : Imag{Imag_}, Parity{Parity_}, Alg{Alg_} {}
  OperatorTrait( int, int, Common::Gamma::Algebra Alg_ ) : Alg{Alg_} {}
};

static constexpr int NumMixed{ 5 };
static const std::array<OperatorTrait, NumMixed> opTraits
{
  OperatorTrait( 0, 1, Common::Gamma::Algebra::Gamma5 ),
  OperatorTrait( 0, -1, Common::Gamma::Algebra::GammaTGamma5 ),
  OperatorTrait( 1, 1, Common::Gamma::Algebra::GammaXGamma5 ),
  OperatorTrait( 1, 1, Common::Gamma::Algebra::GammaYGamma5 ),
  OperatorTrait( 1, 1, Common::Gamma::Algebra::GammaZGamma5 ),
};*/
// The order here must match the algebra array for the model code to work
static constexpr int idxg5{ 0 };
static constexpr int idxgT5{ 1 };
//static constexpr int idxgX5{ 2 };
//static constexpr int idxgY5{ 3 };
//static constexpr int idxgZ5{ 4 };

static constexpr int ModelNumIndices{ 2 }; // Model should include operators with these indices

struct SinkSource
{
  int Sink;
  int Source;
};

enum ModelType : int { DiagSinkSource, DiagSourceOnly };

struct Parameters
{
  ModelType modelType;
  std::string modelArgs;
  int NumSamples;
  int Exponent;
  bool bSaveCorr;
  bool bAllReplicas;
  bool bTryConjugate;
  std::string InBase;
  std::string OutBase;
  std::vector<Common::Momentum> Momenta;
  int tmin;
  int tmax;
  int step;
};

/*int ModelOps::Index( const std::string &Op ) const
{
  int j = 0;
  while( j < Names.size() && !Common::EqualIgnoreCase( Names[j], Op ) )
    j++;
  if( j >= Names.size() )
    throw std::runtime_error( "Model doesn't contain operator " + Op );
  return j;
}*/

class MixedOp
{
protected:
  const Parameters & Par;
  const Model &model;
  const int NumParams;
  const int Exponent;
  Fold CorrMixed;
  std::vector<std::string> OpNames;
  std::vector<std::string> FileOpNames; // Operator names that came from the filenames ... somewhat redundant
  // Set by LoadCorrelators()
  int NumSamples; // maximum of a) user parameters b) raw correlators c) model (if doing all replicas)
  int Nt;

  MixedOp( const Model &model_, const Parameters &par );
  virtual const std::string & Description() = 0;
  void LoadCorrelator( Fold &Corr, const std::string &InFileName );
  void LoadCorrelatorComplete();
  virtual void ValidateOpNames() = 0;
  virtual void LoadCorrelators( const std::string &InBase ) = 0;
  virtual void MixingAngle(double costheta, double sintheta) = 0;
  virtual const Fold & Corr0() const = 0;
  virtual void AppendModelOps( std::string &s ) const = 0;
  void DoOneAngle( int degrees, std::string &Out, std::size_t OutLen, bool bSaveCorr, bool bPrint = true );
public:
  virtual ~MixedOp() = default;
  static void Make( const Model &model_, const Parameters & Par );
};

MixedOp::MixedOp( const Model &model_, const Parameters &par ) : Par{ par }, model{ model_ },
  NumParams{ model.Nt() },
  Exponent{ Par.Exponent >= 0 ? Par.Exponent : Par.Exponent + model.NumExponents }
{
  if( model.NumExponents < 2 )
    throw std::runtime_error( "NumExponents=" + std::to_string( model.NumExponents ) + ", at least 2 required" );
  if( Exponent < 0 || Exponent >= model.NumExponents )
    throw std::runtime_error( "NumExponents=" + std::to_string( model.NumExponents ) + ", excited state " + std::to_string( Par.Exponent ) + " not available for model" );
  const int FileOps{ static_cast<int>( model.OpNames.size() ) };
  const int NtExpected{ model.NumExponents * ( FileOps + 1 ) + 1 };
  if( NumParams < NtExpected - 1 || NumParams > NtExpected )
  throw std::runtime_error("Number of parameters=" + std::to_string( NumParams ) + ", "
                           + std::to_string( NtExpected ) + " expected");
}

// Keep track of info on correlators as they are loaded
void MixedOp::LoadCorrelator( Fold &Corr, const std::string &InFileName )
{
  Corr.Read( InFileName, "    ", &FileOpNames );
  CorrMixed.FileList.push_back( InFileName );
  if( &Corr == &Corr0() )
  {
    Nt = Corr.Nt();
    NumSamples = ( Par.NumSamples > 0 && Par.NumSamples < Corr.NumSamples() )
                 ? Par.NumSamples : Corr.NumSamples();
    if( Par.bAllReplicas && NumSamples > model.NumSamples() )
      NumSamples = model.NumSamples();
    CorrMixed.NtUnfolded = Corr.NtUnfolded;
    CorrMixed.parity = Common::Parity::Unknown;
    CorrMixed.reality = Common::Reality::Real;
    CorrMixed.sign = Common::Sign::Positive;
    CorrMixed.Conjugated = false;
    CorrMixed.t0Negated = false;
  }
  else
  {
    Corr0().IsCompatible( Corr, &NumSamples );
  }
}

void MixedOp::LoadCorrelatorComplete()
{
  CorrMixed.resize( NumSamples, Nt );
  CorrMixed.CopyAttributes( Corr0() );
}

// Mixed operator at sink and source
class MixedOp_SS : public MixedOp
{
  friend void MixedOp::Make( const Model &model_, const Parameters & Par );
protected:
  std::array<SinkSource, ModelNumIndices> ss;
  std::array<std::array<Fold, ModelNumIndices>, ModelNumIndices> Corr;
  virtual const Fold & Corr0() const { return Corr[0][0]; }
  virtual void AppendModelOps( std::string &s ) const;
  MixedOp_SS( const Model &model_, const Parameters &par ) : MixedOp( model_, par )
  {
    if( !par.modelArgs.empty() )
      throw std::runtime_error( "Unexpected model argument \"" + par.modelArgs + "\"" );
  }
  virtual ~MixedOp_SS() = default;
  virtual const std::string & Description();
  virtual void ValidateOpNames();
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
};

void MixedOp_SS::AppendModelOps( std::string &s ) const
{
  assert( 0 && "Check naming for source and sink-mixed op" );
}

const std::string & MixedOp_SS::Description()
{
  static const std::string s{ "source & sink" };
  return s;
}

// Find the indices in the file for the operators I need
void MixedOp_SS::ValidateOpNames()
{
  // Do we have the right number of operators
  const int NumOps{ ModelNumIndices * ( model.Factorised ? 1 : 2 ) };
  if( model.OpNames.size() != NumOps )
    throw std::runtime_error( "Expecting " + std::to_string( NumOps ) + " operators, but model has " + std::to_string( model.OpNames.size() ) );
  // Are they the right operators
  if( model.Factorised )
  {
    OpNames = model.OpNames;
    for( int i = 0; i < NumOps; i++ )
    {
      ss[i].Source = i;
      ss[i].Sink = i;
    }
  }
  else
  {
    SinkSource Count[ModelNumIndices];
    for( int i = 0; i < ModelNumIndices; i++ )
    {
      Count[i].Sink = 0;
      Count[i].Source = 0;
    }
    for( int i = 0; i < NumOps; i++ )
    {
      std::string sOp{ model.OpNames[i] };
      bool bOK{ false };
      bool bSource{ false };
      std::size_t p = sOp.length();
      if( p > 4 )
      {
        std::string Suffix{ sOp.substr( p - 4, 4 ) };
        if( Common::EqualIgnoreCase( Suffix, "_src" ) )
        {
          bOK = true;
          bSource = true;
        }
        else if( Common::EqualIgnoreCase( Suffix, "_snk" ) )
          bOK = true;
      }
      if( !bOK )
        throw new std::runtime_error( "Operator " + sOp + " does not end in _src or _snk" );
      sOp.resize( p - 4 );
      int idx{ 0 };
      while( idx < OpNames.size() && !Common::EqualIgnoreCase( sOp, OpNames[i] ) )
        ++idx;
      if( idx == OpNames.size() )
      {
        if( idx == ModelNumIndices )
          throw new std::runtime_error( "Too many operators in model" );
        OpNames.push_back( sOp );
      }
      if( bSource )
      {
        ss[idx].Source = i;
        Count[idx].Source++;
      }
      else
      {
        ss[idx].Sink = i;
        Count[idx].Sink++;
      }
    }
    for( int i = 0; i < ModelNumIndices; i++ )
      if( Count[i].Sink != 1 || Count[i].Source != 1 )
        throw new std::runtime_error( "Expected 1 _src and 1 _snk for each operator" );
  }
}

void MixedOp_SS::LoadCorrelators( const std::string &InBase )
{
  // Load all of the correlators
  for( int iSnk = 0; iSnk < ModelNumIndices ; iSnk++ ) // upper limit NumMixed?
  {
    std::string InFile{ InBase };
    InFile.append( 1, '_' );
    InFile.append( OpNames[iSnk] );
    const std::size_t InFileLen{ InFile.length() };
    for( int iSrc = 0; iSrc <= ( model.Factorised ? iSnk : ModelNumIndices - 1 ); iSrc++ )
    {
      InFile.resize( InFileLen );
      InFile.append( 1, '_' );
      InFile.append( OpNames[iSrc] );
      std::string InFileName{ Common::MakeFilename(InFile,Common::sFold,model.Name_.Seed,DEF_FMT)};
      if( model.Factorised && iSrc != iSnk && !Common::FileExists( InFileName ) )
      {
        InFileName = InBase;
        InFileName.append( 1, '_' );
        InFileName.append( OpNames[iSrc] );
        InFileName.append( 1, '_' );
        InFileName.append( OpNames[iSnk] );
        InFileName = Common::MakeFilename( InFileName, Common::sFold, model.Name_.Seed, DEF_FMT );
      }
      LoadCorrelator( Corr[iSnk][iSrc], InFileName );
    }
  }
}

void MixedOp_SS::MixingAngle(double costheta, double sintheta)
{
  const double cos_sq_theta{ costheta * costheta };
  const double sin_sq_theta{ sintheta * sintheta };
  const double cos_sin_theta{ costheta * sintheta };
  const double * pCoeff{ model[Model::idxCentral] + model.NumExponents + Exponent * model.OpNames.size() };
  const double * pPP{ Corr[idxg5 ][idxg5 ][Fold::idxCentral] };
  const double * pAP{ Corr[idxgT5][idxg5 ][Fold::idxCentral] };
  const double * pPA{ nullptr };
  const double * pAA{ Corr[idxgT5][idxgT5][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_P1_snk{ 0 };
  double A_P1_src{ 0 };
  double A_A1_snk{ 0 };
  double A_A1_src{ 0 };
  double Op_PP{ 0 };
  double Op_AP{ 0 };
  double Op_PA{ 0 };
  double Op_AA{ 0 };
  if( !model.Factorised )
                 pPA= Corr[idxg5 ][idxgT5][Fold::idxCentral];
  for( int i = Fold::idxCentral; i < NumSamples; i++ )//, pCoeff += model.Nt() )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_P1_snk = pCoeff[ss[idxg5 ].Sink];
      A_P1_src = pCoeff[ss[idxg5 ].Source];
      A_A1_snk = pCoeff[ss[idxgT5].Sink];
      A_A1_src = pCoeff[ss[idxgT5].Source];
      Op_PP = cos_sq_theta  / ( A_P1_snk * A_P1_src );
      Op_AP = cos_sin_theta / ( A_A1_snk * A_P1_src );
      Op_PA = cos_sin_theta / ( A_P1_snk * A_A1_src );
      Op_AA = sin_sq_theta  / ( A_A1_snk * A_A1_src );
      pCoeff += NumParams;
    }
    for( int t = 0; t < Nt; t++ )
      if( model.Factorised )
        *pDst++ = Op_PP * *pPP++ + Op_AA * *pAA++ + 2 * Op_AP * *pAP++;
      else
        *pDst++ = Op_PP * *pPP++ + Op_AA * *pAA++ +     Op_AP * *pAP++ + Op_PA * *pPA++;
  }
}

class MixedOp_S : public MixedOp
{
  friend void MixedOp::Make( const Model &model_, const Parameters & Par );
protected:
  std::string sSink;
  std::array<Fold, ModelNumIndices> Corr;
  virtual const Fold & Corr0() const { return Corr[0]; }
  virtual void AppendModelOps( std::string &s ) const;
  MixedOp_S( const Model &model_, const Parameters &par )
  : MixedOp( model_, par ), sSink{ Par.modelArgs.empty() ? model_.OpNames[0] : Par.modelArgs } {}
  virtual ~MixedOp_S() = default;
  virtual const std::string &Description();
  virtual void ValidateOpNames();
  virtual void LoadCorrelators( const std::string &InBase );
  virtual void MixingAngle(double costheta, double sintheta);
};

void MixedOp_S::AppendModelOps( std::string &s ) const
{
  s.append( sSink );
  s.append( 1, '_' );
  s.append( OpNames[0] );
  s.append( OpNames[1] );
}

const std::string &MixedOp_S::Description()
{
  static const std::string s{ "source" };
  return s;
}

// Find the indices in the file for the operators I need
void MixedOp_S::ValidateOpNames()
{
  // Do we have the right number of operators
  if( !model.Factorised )
    throw std::runtime_error( "Model must be factorised" );
  if( model.OpNames.size() != ModelNumIndices )
    throw std::runtime_error( "Expecting " + std::to_string( ModelNumIndices )
                             + " operators, but model has " + std::to_string( model.OpNames.size() ) );
  OpNames = model.OpNames;
}

void MixedOp_S::LoadCorrelators( const std::string &InBase )
{
  // Load all of the correlators
  std::string InFile{ InBase };
  InFile.append( 1, '_' );
  InFile.append( sSink );
  InFile.append( 1, '_' );
  const std::size_t InFileLen{ InFile.length() };
  for( int iSrc = 0; iSrc < ModelNumIndices; iSrc++ )
  {
    InFile.resize( InFileLen );
    InFile.append( OpNames[iSrc] );
    std::string InFileName{ Common::MakeFilename( InFile, Common::sFold, model.Name_.Seed, DEF_FMT ) };
    if(!Common::FileExists( InFileName ) && Par.bTryConjugate && !Common::EqualIgnoreCase( sSink, OpNames[iSrc]))
    {
      std::cout << "Warning: loading conjugate of " << InFileName << Common::NewLine;
      std::string InFileName{ InBase };
      InFileName.append( 1, '_' );
      InFileName.append( OpNames[iSrc] );
      InFileName.append( 1, '_' );
      InFileName.append( sSink );
      InFileName = Common::MakeFilename( InFileName, Common::sFold, model.Name_.Seed, DEF_FMT );
    }
    LoadCorrelator( Corr[iSrc], InFileName );
  }
}

void MixedOp_S::MixingAngle(double costheta, double sintheta)
{
  const double * pCoeff{ model[Model::idxCentral] + model.NumExponents + Exponent * model.OpNames.size() };
  const double * pPP{ Corr[0][Fold::idxCentral] };
  const double * pPW{ Corr[1][Fold::idxCentral] };
  scalar * pDst = CorrMixed[Fold::idxCentral];
  double A_PP{ 0 };
  double A_PW{ 0 };
  double Op_PP{ 0 };
  double Op_PW{ 0 };
  for( int i = Fold::idxCentral; i < NumSamples; i++ )
  {
    if( Par.bAllReplicas || i == Fold::idxCentral )
    {
      A_PP = pCoeff[0];
      A_PW = pCoeff[1];
      Op_PP = costheta / ( A_PP * A_PP );
      Op_PW = sintheta / ( A_PP * A_PW );
      pCoeff += NumParams;
    }
    for( int t = 0; t < Nt; t++ )
      *pDst++ = Op_PP * *pPP++ + Op_PW * *pPW++;
  }
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
  AppendModelOps( Out );
  CorrMixed.MakeCorrSummary( nullptr );
  const Common::SeedType Seed{ model.Name_.Seed };
  if( bSaveCorr )
    CorrMixed.Write( Common::MakeFilename( Out, Common::sBootstrap, Seed, DEF_FMT ) );
  CorrMixed.WriteSummary( Common::MakeFilename( Out, Common::sBootstrap, Seed, TEXT_EXT ) );
}

void MixedOp::Make( const Model &model_, const Parameters & Par )
{
  // Extract momentum from model name to get our base name.
  // Make sure momentum specified either as command-line option, or in filename
  std::string Base{ model_.Name_.Base };
  std::vector<Common::Momentum> fileMom{ 1 };
  const bool bGotFileMom{ fileMom[0].Extract( Base ) };
  if( Par.Momenta.empty() && !bGotFileMom )
    throw std::runtime_error( "Momentum needs to be specified" );
  const std::vector<Common::Momentum> &MyMomenta{ Par.Momenta.size() ? Par.Momenta : fileMom };

  // Display info about the model
  std::cout << "  Base: " << Base << NewLine;
  for( std::size_t i = model_.Name_.Extra.size(); i > 0; i-- )
    std::cout << "  Extra[" << std::to_string(i-1) << "]: " << model_.Name_.Extra[i-1] << NewLine;
  std::cout << "  Fit: " << ( model_.Factorised ? "F" : "Unf" ) << "actorised, ti="
            << model_.ti << ", tf=" << model_.tf << NewLine;
  if( bGotFileMom )
    std::cout << "  Model Momentum: " << fileMom[0] << NewLine;
  for( std::size_t i = 0; i < model_.OpNames.size(); i++ )
    std::cout << "  Model op[" << i << "]: " << model_.OpNames[i] << NewLine;

  // Now make our mixed operator
  std::unique_ptr<MixedOp> mixed;
  {
    switch( Par.modelType )
    {
      case DiagSinkSource:
        mixed.reset( new MixedOp_SS( model_, Par ) );
        break;
      case DiagSourceOnly:
        mixed.reset( new MixedOp_S( model_, Par ) );
        break;
      default:
        throw std::runtime_error( "Model type " + std::to_string( Par.modelType ) + " not supported" );
    }
  }
  // Validate the model
  mixed->ValidateOpNames();

  // Construct optimised model for each momentum
  Base.append( "_p_" );
  const std::size_t BaseLen{ Base.length() };
  for( const Common::Momentum &p : MyMomenta )
  {
    Base.resize( BaseLen );
    Base.append( p.to_string( Sep ) );
    std::string BaseFit{ Base };
    for( std::size_t i = model_.Name_.Extra.size(); i > 1; i-- )
    {
      BaseFit.append( 1, '.' );
      BaseFit.append( model_.Name_.Extra[i-1] );
    }
    const std::string InBase{ Par.InBase + Base };
    std::cout << "  " << BaseFit << "\n";
    mixed->LoadCorrelators( InBase );
    mixed->LoadCorrelatorComplete();

    std::string Out{ Par.OutBase + BaseFit };
    std::cout << "    Writing " << mixed->NumSamples << " samples of " << mixed->Description()
              << " mixed operator to:\n    " << Out << NewLine << "   ";
    Out.append( ".theta_" );
    const std::size_t OutLen{ Out.length() };
    bool bDo0  = true; //!Par.bSaveCorr; // If we're not saving the correlators, always do 0 and 90 degrees
    bool bDo90 = bDo0;
    for( int degrees = Par.tmin; ; degrees+=Par.step )
    {
      mixed->DoOneAngle( degrees, Out, OutLen, Par.bSaveCorr, degrees == Par.tmin || !( degrees % 10 ) );
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
      mixed->DoOneAngle( 0, Out, OutLen, false ); // We're outside the range requested, no need to save correlators
    if( bDo90 )
      mixed->DoOneAngle( 90, Out, OutLen, false );
    std::cout << NewLine;
  }
}

int main( int argc, const char *argv[] )
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"a", CL::SwitchType::Single, "1"},
      {"p", CL::SwitchType::Single, nullptr},
      {"i", CL::SwitchType::Single, "" },
      {"m", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"e", CL::SwitchType::Single, "-1"},
      {"tmin", CL::SwitchType::Single, "-90"},
      {"tmax", CL::SwitchType::Single, "90"},
      {"step", CL::SwitchType::Single, "1"},
      {"c", CL::SwitchType::Flag, nullptr},
      {"rep", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      Parameters Par;
      {
        const std::string sMT{ cl.SwitchValue<std::string>("a") };
        std::istringstream s{ sMT };
        int i;
        if( ! ( s >> i ) )
          throw std::runtime_error( "Model type \"" + sMT + "\" must start with an integer" );
        Par.modelType = static_cast<ModelType>( i );
        // optional comma, followed by model options
        if( !Common::StreamEmpty( s ) )
        {
          // Read an optional comma, tiptoeing around whitespace
          if( s.peek() == ',' )
          {
            s.get();
            if( Common::StreamEmpty( s ) )
              throw std::runtime_error( "Model parameter empty (trailing comma)" );
          }
          getline( s, Par.modelArgs );
        }
      }
      const std::string modelBase{ cl.SwitchValue<std::string>("m") };
      Par.InBase = Common::AppendSlash( cl.SwitchValue<std::string>("i") );
      Par.OutBase = Common::AppendSlash( cl.SwitchValue<std::string>("o") );
      Par.NumSamples = cl.SwitchValue<int>("n");
      Par.Exponent = cl.SwitchValue<int>("e");
      Par.bAllReplicas = cl.GotSwitch("rep");
      Par.bTryConjugate = cl.GotSwitch("c");
      Par.bSaveCorr = cl.GotSwitch("savecorr");
      if( cl.GotSwitch( "p" ) )
        Par.Momenta = Common::ArrayFromString<Common::Momentum>( cl.SwitchValue<std::string>("p") );
      Par.tmin = cl.SwitchValue<int>("tmin");
      Par.tmax = cl.SwitchValue<int>("tmax");
      Par.step = std::abs( cl.SwitchValue<int>("step") );
      if( Par.NumSamples < 0 )
        throw std::runtime_error( "NumSamples must be >= 0" );
      if( Par.tmin > Par.tmax )
        Par.step *= -1;
      bShowUsage = false;
      for(const std::string &ModelFile : Common::glob(cl.Args.begin(),cl.Args.end(),modelBase.c_str()))
      {
        Model model;
        model.Read( ModelFile, "Read model " );
        MixedOp::Make( model, Par );
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> Model1 [Model2 ...]\n"
    "Create a mixed operator from fit parameters and bootstrap replicas\n"
    "<options>\n"
    "-a     Model type 0=source/sink, 1=sink only (default)\n"
    "-p     Comma separated list of momenta (default: same as model file)\n"
    "-i     Input path for folded bootstrap replicas\n"
    "-m     Input path for model files\n"
    "-o     Output path\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-e     Which excited state to use in model (default -1), -ve counts from highest state\n"
    "--tmin Minimum theta (default -90)\n"
    "--tmax Maximum theta (default  90)\n"
    "--step Steps between theta (default 1)\n"
    "Flags:\n"
    "-c         Try loading Conjugate operator if bootstrap sample missing"
    "--rep      Use per replica values of overlap constants in construction of model\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--help     This message\n";
  }
  return iReturn;
}
