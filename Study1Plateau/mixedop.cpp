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

using SampleD = Common::SampleD;
using SampleC = Common::SampleC;

static const std::string Sep{ "_" };
static const std::string NewLine{ "\n" };

static constexpr int NumMixed{ 5 };
static const std::array<Common::Gamma::Algebra, NumMixed> MixedAlg
{
  Common::Gamma::Algebra::Gamma5,
  Common::Gamma::Algebra::GammaXGamma5,
  Common::Gamma::Algebra::GammaYGamma5,
  Common::Gamma::Algebra::GammaZGamma5,
  Common::Gamma::Algebra::GammaTGamma5,
};
// The order here must match the algebra array for the model code to work
static constexpr int idxg5{ 0 };
static constexpr int idxgX5{ 1 };
static constexpr int idxgY5{ 2 };
static constexpr int idxgZ5{ 3 };
static constexpr int idxgT5{ 4 };

struct Parameters
{
  int shift;
  int fold;
  int NumExponents;
  int NumSamples;
  bool bSaveCorr;
  std::string InBase;
  std::string OutBase;
  std::vector<Common::Momentum> Momenta;
};

void MixingAngle( Common::SampleC &CorrMixed,
                  const std::array<std::array<Common::SampleC, NumMixed>, NumMixed> &Corr,
                  const Common::SampleD &Model,
                  int degrees, int NumSamples, int Nt )
{
  double costheta, sintheta;
  switch( ( degrees + 360 ) % 360 )
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
  const double A_P1{ Model[SampleD::idxCentral][4] };
  const double A_A1{ Model[SampleD::idxCentral][5] };
  const double Op_PP{ cos_sq_theta / ( A_P1 * A_P1 ) };
  const double Op_AP{ cos_sin_theta / ( A_P1 * A_A1 ) };
  const double Op_AA{ sin_sq_theta / ( A_A1 * A_A1 ) };
  const std::complex<double> * pPP = Corr[idxg5 ][idxg5 ][SampleC::idxCentral];
  const std::complex<double> * pPA = Corr[idxg5 ][idxgT5][SampleC::idxCentral];
  const std::complex<double> * pAP = Corr[idxgT5][idxg5 ][SampleC::idxCentral];
  const std::complex<double> * pAA = Corr[idxgT5][idxgT5][SampleC::idxCentral];
  CorrMixed.resize( NumSamples, Nt );
  std::complex<double> * pDst = CorrMixed[SampleC::idxCentral];
  for( int i = SampleC::idxCentral; i < NumSamples; i++ )
  {
    for( int t = 0; t < Nt; t++ )
      *pDst++ = Op_PP * (*pPP++).real() + Op_AA * (*pAA++).real() + Op_AP * ( ( *pAP++ + *pPA++ ).real() );
  }
}

void MakeModel( const std::string & ModelFile, const Parameters & Par )
{
  // Load Model
  std::cout << ModelFile << NewLine;
  std::string GroupName;
  SampleD Model{ ModelFile, GroupName };
  std::cout << "  Group " << GroupName << "\n";
  int NumSamples{ Model.NumSamples() };
  std::cout << "  Count " << NumSamples;
  if( Par.NumSamples > 0 && NumSamples > Par.NumSamples )
  {
    NumSamples = Par.NumSamples;
    std::cout << " in file, reduced to " << NumSamples;
  }
  std::cout << NewLine;

  // Extract and validate parameters
  Common::FileNameAtt attModel{ ModelFile };
  if( !attModel.bSeedNum )
    throw std::runtime_error( "Seed missing" );
  std::vector<std::string> OpNames;
  if( !Common::ExtractSuffixSplit( attModel.Base, OpNames ) || OpNames.size() != 2 )
    throw std::runtime_error( "Model name has " + std::to_string( OpNames.size() ) + " operators" );
  std::string FitName;
  if( !Common::ExtractSuffix( attModel.Base, FitName, Common::Period.c_str() ) )
    throw std::runtime_error( "Fit parameters missing" );
  std::cout << "  Fit\t" << FitName;
  Common::Momentum fileMom;
  if( fileMom.Extract( attModel.Base ) )
    std::cout << ", p=" << fileMom;
  for( std::size_t i = 0; i < OpNames.size(); i++ )
    std::cout << ( i == 0 ? ", Ops " : ", " ) << OpNames[i];
  std::cout << NewLine;
  const int FileOps{ static_cast<int>( OpNames.size() ) };
  const int NumExponents{ Model.Nt() / ( FileOps + 1 ) };
  std::cout << "  Model " << FileOps << " operators x " << NumExponents << " exponents = "
            << ( ( FileOps + 1 ) * NumExponents ) << " parameters\n";
  if( ( FileOps + 1 ) * NumExponents != Model.Nt() )
    throw std::runtime_error( "Model contains incorrect number of parameters" );

  // Construct optimised model for each momentum
  Common::SampleC CorrMixed;
  std::array<std::array<Common::SampleC, NumMixed>, NumMixed> Corr;
  std::cout << "  Base\t" << attModel.Base << "\n";
  for( const Common::Momentum &p : Par.Momenta )
  {
    const std::string Base{ attModel.Base + "_p_" + p.to_string( Sep ) };
    const std::string BaseFit{ Base + Common::Period + FitName };
    const std::string InBase{ Par.InBase + Base };
    std::cout << "  " << BaseFit << "\n";
    int Nt         = -777;
    // Load all of the correlators
    for( int iSnk = 0; iSnk < NumMixed; iSnk++ )
    {
      std::string InFile{ InBase };
      InFile.append( Common::Gamma::NameShort( MixedAlg[iSnk], Common::Underscore.c_str() ) );
      const std::size_t InFileLen{ InFile.length() };
      for( int iSrc = 0; iSrc < NumMixed; iSrc++ )
      {
        InFile.resize( InFileLen );
        InFile.append( Common::Gamma::NameShort( MixedAlg[iSrc], Common::Underscore.c_str() ) );
        const std::string InFileName{ Common::MakeFilename( InFile, Common::sBootstrap,
                                                           attModel.Seed, DEF_FMT ) };
        std::string GroupName;
        Corr[iSnk][iSrc].Read( InFileName, GroupName );
        std::cout << "    " << InFileName << " (" << GroupName << ")" << "\n";
        if( iSnk == 0 && iSrc == 0 )
        {
          Nt = Corr[iSnk][iSrc].Nt();
          NumSamples = Corr[iSnk][iSrc].NumSamples();
        }
        else
        {
          if( Nt != Corr[iSnk][iSrc].Nt() )
            throw std::runtime_error( "Nt!=" + std::to_string( Nt ) );
          if( NumSamples > Corr[iSnk][iSrc].NumSamples() )
            NumSamples = Corr[iSnk][iSrc].NumSamples();
        }
      }
    }
    std::string OutBase{ Par.OutBase + BaseFit };
    std::cout << "    Writing " << NumSamples << " samples to " << OutBase << NewLine << "   ";
    OutBase.append( ".theta_" );
    const std::size_t OutLen{ OutBase.length() };
    for( int degrees = -90; degrees <= 90; degrees ++ )
    {
      if( !( degrees % 10 ) )
        std::cout << " " << std::to_string( degrees );
      MixingAngle( CorrMixed, Corr, Model, degrees, NumSamples, Nt );
      OutBase.resize( OutLen );
      OutBase.append( std::to_string( degrees ) );
      if( Par.bSaveCorr )
        CorrMixed.Write( Common::MakeFilename(OutBase, Common::sBootstrap, attModel.Seed, DEF_FMT ) );
      Common::SummariseBootstrapCorr( CorrMixed, OutBase, attModel.Seed );
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
      {"s", CL::SwitchType::Single, "0"},
      {"f", CL::SwitchType::Single, "0"},
      {"p", CL::SwitchType::Single, "0_0_0"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"e", CL::SwitchType::Single, "0"},
      {"n", CL::SwitchType::Single, "0"},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      Parameters Par;
      Par.shift = cl.SwitchValue<int>("s");
      Par.fold = cl.SwitchValue<int>("f");
      Par.InBase = Common::AppendSlash( cl.SwitchValue<std::string>("i") );
      Par.OutBase = Common::AppendSlash( cl.SwitchValue<std::string>("o") );
      Par.NumExponents = cl.SwitchValue<int>( "e" );
      Par.NumSamples = cl.SwitchValue<int>("n");
      if( Par.NumSamples < 0 )
        throw std::runtime_error( "NumSamples must be >= 0" );
      Par.bSaveCorr = cl.GotSwitch("savecorr");
      Par.Momenta = Common::ArrayFromString<Common::Momentum>( cl.SwitchValue<std::string>("p") );
      if( Par.shift != 0 )
        throw std::runtime_error( "Shift option not supported yet" );
      if( Par.fold != 0 )
        throw std::runtime_error( "Fold option not supported yet" );
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
    "-s     time Shift (default 0)\n"
    "-f     Fold the correlator, 0=don't fold (default), 1=+parity, -1=-parity\n"
    "-p     Comma separated list of momenta, default 0_0_0"
    "-i     Input path for bootstrap replicas\n"
    "-o     Output path\n"
    "-e     number of Exponents (default same as number of operators)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "Flags:\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--help     This message\n";
  }
  return iReturn;
}
