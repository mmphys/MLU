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
using Sample = Common::Sample<scalar>;

static const std::string Sep{ "_" };
static const std::string NewLine{ "\n" };

struct OperatorTrait
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
};
// The order here must match the algebra array for the model code to work
static constexpr int idxg5{ 0 };
static constexpr int idxgT5{ 1 };
static constexpr int idxgX5{ 2 };
static constexpr int idxgY5{ 3 };
static constexpr int idxgZ5{ 4 };

static constexpr int ModelNumIndices{ 2 }; // Model should include operators with these indices

struct SinkSource
{
  int Sink;
  int Source;
};

struct Parameters
{
  int NumSamples;
  bool bSaveCorr;
  std::string InBase;
  std::string OutBase;
  std::vector<Common::Momentum> Momenta;
};

void MixingAngle( const Model &model, const std::array<SinkSource, ModelNumIndices> &ss,
                  const std::array<std::array<Fold, NumMixed>, NumMixed> &Corr, Sample &CorrMixed,
                  int degrees )
{
  const int NumParams{ model.Nt() };
  const int NumSamples{ CorrMixed.NumSamples() };
  const int Nt{ CorrMixed.Nt() };
  double costheta, sintheta;
  degrees = ( ( degrees ) + 360 ) % 360;
  switch( degrees )
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
      const double theta{ M_PI * degrees / 180 };
      costheta = cos( theta );
      sintheta = sin( theta );
    }
  }
  const double cos_sq_theta{ costheta * costheta };
  const double sin_sq_theta{ sintheta * sintheta };
  const double cos_sin_theta{ costheta * sintheta };
  const double * pCoeff{ model[Model::idxCentral] + ( model.NumExponents + model.OpNames.size() ) };
  const double A_P1_snk{ pCoeff[ss[idxg5 ].Sink] };
  const double A_P1_src{ pCoeff[ss[idxg5 ].Source] };
  const double A_A1_snk{ pCoeff[ss[idxgT5].Sink] };
  const double A_A1_src{ pCoeff[ss[idxgT5].Source] };
  const double Op_PP{ cos_sq_theta  / ( A_P1_snk * A_P1_src ) };
  const double Op_AP{ cos_sin_theta / ( A_A1_snk * A_P1_src ) };
  const double Op_PA{ cos_sin_theta / ( A_P1_snk * A_A1_src ) };
  const double Op_AA{ sin_sq_theta  / ( A_A1_snk * A_A1_src ) };
  const double * pPP{ Corr[idxg5 ][idxg5 ][Sample::idxCentral] };
  const double * pAP{ Corr[idxgT5][idxg5 ][Sample::idxCentral] };
  const double * pPA{ nullptr };
  const double * pAA{ Corr[idxgT5][idxgT5][Sample::idxCentral] };
  scalar * pDst = CorrMixed[Sample::idxCentral];
  if( !model.Factorised )
                 pPA= Corr[idxg5 ][idxgT5][Sample::idxCentral];
  for( int i = Sample::idxCentral; i < NumSamples; i++ )//, pCoeff += model.Nt() )
  {
    for( int t = 0; t < Nt; t++ )
      if( model.Factorised )
        *pDst++ = Op_PP * *pPP++ + Op_AA * *pAA++ + 2 * Op_AP * *pAP++;
      else
        *pDst++ = Op_PP * *pPP++ + Op_AA * *pAA++ +     Op_AP * *pAP++ + Op_PA * *pPA++;
  }
}

int OpIndex( const Model &model, const std::string &Op )
{
  int j = 0;
  while( j < model.OpNames.size() && !Common::EqualIgnoreCase( model.OpNames[j], Op ) )
    j++;
  if( j >= model.OpNames.size() )
    throw std::runtime_error( "Model doesn't contain operator " + Op );
  return j;
}

// Find the indices in the file for the operators I need
void ValidateModelOps( const Model &model, std::array<SinkSource, ModelNumIndices> &ss )
{
  // Do we have the right number of operators
  const int NumOps{ ModelNumIndices * ( model.Factorised ? 1 : 2 ) };
  if( model.OpNames.size() != NumOps )
    throw std::runtime_error( "Expecting " + std::to_string( NumOps ) + " operators, but model has " + std::to_string( model.OpNames.size() ) );
  // Are they the right operators
  for( int i = 0; i < ModelNumIndices; i++ )
  {
    const std::string &sOp{ Common::Gamma::nameShort[ static_cast<int>( opTraits[i].Alg ) ] };
    if( model.Factorised )
    {
      ss[i].Source = OpIndex( model, sOp );
      ss[i].Sink = ss[i].Source;
    }
    else
    {
      ss[i].Source = OpIndex( model, sOp + "_src" );
      ss[i].Sink = OpIndex( model, sOp + "_snk" );
    }
  }
}

void MakeModel( const std::string & ModelFile, const Parameters & Par )
{
  // Load Model
  Model model;
  std::string GroupName;
  model.Read( ModelFile, GroupName, "Read model " );

  // Validate the model
  std::array<SinkSource, ModelNumIndices> idxOp;
  ValidateModelOps( model, idxOp );
  if( model.NumExponents < 2 )
    throw std::runtime_error( "NumExponents=" + std::to_string( model.NumExponents ) + ", at least 2 required" );
  const int FileOps{ static_cast<int>( model.OpNames.size() ) };
  {
    const int NtExpected{ model.NumExponents * ( FileOps + 1 ) + 1 };
    if( model.Nt() < NtExpected - 1 || model.Nt() > NtExpected )
    throw std::runtime_error( "Number of parameters=" + std::to_string( model.Nt() ) + ", " + std::to_string( NtExpected ) + " expected" );
  }

  // Base name doesn't have momentum
  std::string Base{ model.Name_.Base };
  Common::Momentum fileMom;
  const bool bGotFileMom{ fileMom.Extract( Base ) };

#ifdef DEBUG
  std::cout << "  Base: " << model.Name_.Base << NewLine;
  for( std::size_t i = model.Name_.Extra.size(); i > 0; i-- )
    std::cout << "  Extra[" << std::to_string(i-1) << "]: " << model.Name_.Extra[i-1] << NewLine;
#endif
  std::cout << "  Fit: " << ( model.Factorised ? "F" : "Unf" ) << "actorised, ti=" << model.ti << ", tf=" << model.tf << NewLine;
  if( bGotFileMom )
    std::cout << "  File Momentum: " << fileMom << NewLine;
  for( std::size_t i = 0; i < model.OpNames.size(); i++ )
    std::cout << "  File op[" << i << "]: " << model.OpNames[i] << NewLine;
  // NumSamples is the maximum of user selected and what's in the file
  int NumSamples{ Par.NumSamples > 0 && Par.NumSamples < model.NumSamples() ? Par.NumSamples : model.NumSamples() };
  std::cout << "  Count " << NumSamples << NewLine;

  // Construct optimised model for each momentum
  Sample CorrMixed;
  std::array<std::array<Fold, NumMixed>, NumMixed> Corr;
  Base.append( "_p_" );
  const std::size_t BaseLen{ Base.length() };
  for( const Common::Momentum &p : Par.Momenta )
  {
    Base.resize( BaseLen );
    Base.append( p.to_string( Sep ) );
    std::string BaseFit{ Base };
    for( std::size_t i = model.Name_.Extra.size(); i > 1; i-- )
    {
      BaseFit.append( 1, '.' );
      BaseFit.append( model.Name_.Extra[i-1] );
    }
    const std::string InBase{ Par.InBase + Base };
    std::cout << "  " << BaseFit << "\n";
    int Nt       = -777;
    // Load all of the correlators
    for( int iSnk = 0; iSnk <= idxgT5 ; iSnk++ ) // upper limit NumMixed?
    {
      std::string InFile{ InBase };
      InFile.append( Common::Gamma::NameShort( opTraits[iSnk].Alg, Common::Underscore.c_str() ) );
      const std::size_t InFileLen{ InFile.length() };
      for( int iSrc = 0; iSrc <= ( model.Factorised ? iSnk : idxgT5 ); iSrc++ ) // upper limit NumMixed?
      {
        InFile.resize( InFileLen );
        InFile.append( Common::Gamma::NameShort( opTraits[iSrc].Alg, Common::Underscore.c_str() ) );
        std::string InFileName{ Common::MakeFilename(InFile,Common::sFold,model.Name_.Seed,DEF_FMT)};
        if( model.Factorised && iSrc != iSnk && !Common::FileExists( InFileName ) )
        {
          InFileName = InBase;
          InFileName.append( Common::Gamma::NameShort(opTraits[iSrc].Alg,Common::Underscore.c_str()) );
          InFileName.append( Common::Gamma::NameShort(opTraits[iSnk].Alg,Common::Underscore.c_str()) );
          InFileName = Common::MakeFilename( InFileName, Common::sFold, model.Name_.Seed, DEF_FMT );
        }
        std::string GroupName;
        Corr[iSnk][iSrc].Read( InFileName, GroupName, "    " );
        if( iSnk == 0 && iSrc == 0 )
          Nt = Corr[iSnk][iSrc].Nt();
        else
        {
          if( Nt != Corr[iSnk][iSrc].Nt() )
            throw std::runtime_error( "Nt " + std::to_string( Corr[iSnk][iSrc].Nt() ) + "!=" + std::to_string( Nt ) );
        }
        if( NumSamples > Corr[iSnk][iSrc].NumSamples() )
          NumSamples = Corr[iSnk][iSrc].NumSamples();
      }
    }
    std::string Out{ Par.OutBase + BaseFit };
    std::cout << "    Writing " << NumSamples << " samples to " << Out << NewLine << "   ";
    Out.append( ".theta_" );
    const std::size_t OutLen{ Out.length() };
    for( int degrees = -90; degrees <= 90; degrees ++ )
    {
      if( !( degrees % 10 ) )
      {
        std::cout << " " << std::to_string( degrees );
        std::cout.flush();
      }
      CorrMixed.resize( NumSamples, Nt );
      MixingAngle( model, idxOp, Corr, CorrMixed, degrees );
      Out.resize( OutLen );
      Out.append( std::to_string( degrees ) );
      if( Par.bSaveCorr )
        CorrMixed.Write( Common::MakeFilename(Out, Common::sBootstrap, model.Name_.Seed, DEF_FMT));
      CorrMixed.WriteSummary(Common::MakeFilename(Out,Common::sBootstrap,model.Name_.Seed,TEXT_EXT));
    }
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
      {"p", CL::SwitchType::Single, "0_0_0"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      Parameters Par;
      Par.InBase = Common::AppendSlash( cl.SwitchValue<std::string>("i") );
      Par.OutBase = Common::AppendSlash( cl.SwitchValue<std::string>("o") );
      Par.NumSamples = cl.SwitchValue<int>("n");
      Par.bSaveCorr = cl.GotSwitch("savecorr");
      Par.Momenta = Common::ArrayFromString<Common::Momentum>( cl.SwitchValue<std::string>("p") );
      if( Par.NumSamples < 0 )
        throw std::runtime_error( "NumSamples must be >= 0" );
      bShowUsage = false;
      for( const std::string & ModelFile : cl.Args )
        MakeModel( ModelFile, Par );
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
    "Create a mixed operator from fit parameters and bootstrap replicas, where <options> are:\n"
    "-p     Comma separated list of momenta, default 0_0_0\n"
    "-i     Input path for folded bootstrap replicas\n"
    "-o     Output path\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "Flags:\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--help     This message\n";
  }
  return iReturn;
}
