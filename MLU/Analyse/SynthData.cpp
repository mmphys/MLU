/**

 Mike's lattice QCD utilities: Synthetic data generator (e.g. for fit test)
 
 Source file: SynthData.cpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "SynthData.hpp"

std::ostream & operator<<( std::ostream &os, const MeanSigma &ms )
{
  return os << ms.Mean << "±" << ms.Sigma;
}

std::istream & operator>>( std::istream &is, MeanSigma &ms )
{
  Scalar m, s;
  if( is >> m && MLU::NextCharIs( is, ',' ) && is >> s )
  {
    ms.Mean = m;
    ms.Sigma = s;
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

inline void Append( std::string &String, const std::string &Appendage )
{
  if( !Appendage.empty() )
  {
    String.append( 1, '_' );
    String.append( Appendage );
  }
}

const std::array<std::string, Synth::NumOps> Synth::opNames{ "A", "B" };

/*std::array<std::array<MeanSigma, Synth::NumParams>, Synth::MaxExp> Synth::ParamMS{{
  { MeanSigma( 0.8, 0.05 ), MeanSigma( 50, 0.5 ), MeanSigma( 100, .3 ) },
  { MeanSigma( 1.2, 0.1  ), MeanSigma( 25, 1),  MeanSigma( -200, 5 ) },
}};*/

std::string Synth::GetParamName( int e, int p ) const
{
  std::string s{ p == 0 ? "E" : Synth::opNames[p - 1] };
  s.append( std::to_string( e ) );
  return s;
}

std::ostream & Synth::DumpParams( std::ostream &os ) const
{
  for( int e = 0; e < NumExp; ++e )
    for( int p = 0; p < NumParams; ++p )
      os << GetParamName( e, p ) << MLU::Space << ParamMS[e][p] << MLU::NewLine;
  return os;
}

Synth::Synth( const MLU::CommandLine &cl, const std::string MachineNameActual )
: NumExp{ cl.SwitchValue<int>( "e" ) },
  CorrelatorJitter{ cl.SwitchValue<Scalar>( "j" ) },
  MachineName{ cl.GotSwitch( "m" ) ? cl.SwitchValue<std::string>( "m" ) : MachineNameActual },
  nSample{ cl.SwitchValue<int>( "s" ) },
  nBootSample{ cl.SwitchValue<int>( "n" ) },
  outStem{ cl.SwitchValue<std::string>( "o" ) },
  seed{ MLU::RandomCache::DefaultSeed() },
  bVerbose{ cl.GotSwitch( "v" ) },
  bWarnIfExists{ cl.GotSwitch( "w" ) }
{
  if( NumExp < 1 || NumExp > 2 )
    throw std::runtime_error( "Number of exponentials (" + std::to_string( NumExp ) + ") must be 1 or 2" );
  if( MachineName.empty() )
    throw std::invalid_argument( "Machine name can't be empty" );
  static const char pszGt0[] = ") must be greater than 0";
  if( nSample < 1 )
    throw std::runtime_error( "Number of samples (" + std::to_string( nSample ) + pszGt0 );
  if( nBootSample < 1 )
    throw std::runtime_error( "Number of bootstrap samples (" + std::to_string( nBootSample ) + pszGt0 );
  // Read parameters
  for( int e = 0; e < NumExp; ++e )
    for( int p = 0; p < NumParams; ++p )
      ParamMS[e][p] = cl.SwitchValue<MeanSigma>( GetParamName( e, p ) );
}

void Synth::Make( std::string Basename ) const
{
  // Make filenames and check whether they exist
  std::array<std::string,NumOps> aFilename;
  std::string MyBase{ outStem + Basename };
  Append( MyBase, opNames[0] );
  const std::size_t Len{ MyBase.length() };
  for( int i = 0; i < NumOps; ++i )
  {
    MyBase.resize( Len );
    Append( MyBase, opNames[i] );
    aFilename[i] = MLU::MakeFilename( MyBase, MLU::sFold, seed, DEF_FMT );
    if( MLU::FileExists( aFilename[i] ) )
    {
      if( !bWarnIfExists )
        throw std::runtime_error( aFilename[i] + " exists" );
      std::cout << "Overwriting ";
    }
    else
      std::cout << "Creating ";
    std::cout << aFilename[i] << MLU::NewLine;
  }

  // Now make our model
  const MLU::Parity MyParity{ MLU::Parity::Even };
  const int Nt{ 64 };
  const int NtHalf{ Nt / 2 + ( MyParity == MLU::Parity::Odd ? 0 : 1 ) };
  std::array<Fold, NumOps> out;
  std::array<Matrix *, NumOps> CRaw;
  std::array<Matrix *, NumOps> CBinned;
  for( std::size_t f = 0; f < aFilename.size(); ++f )
  {
    out[f].resize( nBootSample, NtHalf );
    out[f].resizeRaw( nSample );
    out[f].resizeBinned( nSample );
    out[f].NtUnfolded_ = Nt;
    out[f].parity = MyParity;
    out[f].reality = MLU::Reality::Real;
    out[f].sign = MLU::Sign::Positive;
    out[f].SetSeed( seed );
    out[f].SeedMachine_ = MachineName;
    out[f].ConfigCount.reserve( nSample );
    CRaw[f] = &out[f].getRaw();
    CBinned[f] = &out[f].getBinned();
  }
  MLU::ConfigCount CC( 0, 1 );
  std::array<std::array<Scalar, NumParams>, MaxExp> Values;
  std::mt19937                     engine( seed );
  std::normal_distribution<Scalar> random;

  if( bVerbose )
  {
    std::cout << "Idx:";
    for( int e = 0; e < NumExp; ++e )
      for( int p = 0; p < NumParams; ++p )
        std::cout << MLU::Space << GetParamName( e, p );
    std::cout << MLU::NewLine;
  }
  for( int idx = 0; idx < nSample; ++idx )
  {
    // Get random values with the correct distribution
    if( bVerbose )
      std::cout << idx << MLU::Colon;
    for( int e = 0; e < NumExp; ++e )
      for( int p = 0; p < NumParams; ++p )
      {
        Values[e][p] = random( engine ) * ParamMS[e][p].Sigma + ParamMS[e][p].Mean;
        if( bVerbose )
          std::cout << '\t' << Values[e][p];
      }
    if( bVerbose )
      std::cout << MLU::NewLine;
    
    for( std::size_t f = 0; f < aFilename.size(); ++f )
    {
      // Now build our model
      for( int t = 0; t < NtHalf; ++t )
      {
        (*CRaw[f])(idx,t) = 0;
        for( int e = 0; e < NumExp; ++e )
        {
          Scalar z = std::exp(-Values[e][ParamE] * t);
          z += std::exp(-Values[e][ParamE] * (Nt - t));
          z *= Values[e][ParamA] * Values[e][ParamA + f];
          (*CRaw[f])(idx,t) += z;
        }
        (*CRaw[f])(idx,t) *= random( engine ) * CorrelatorJitter + 1.;
        (*CBinned[f])(idx,t) = (*CRaw[f])(idx,t);
      }
      out[f].ConfigCount.push_back( CC );
    }
    CC.Config++;
  }
  // Now write the files
  for( std::size_t f = 0; f < aFilename.size(); ++f )
  {
    out[f].Resample();
    out[f].MakeCorrSummary();
    MLU::MakeAncestorDirs( aFilename[f] );
    out[f].Write( aFilename[f], MLU::sFold.c_str() );
    std::string SummaryName{ aFilename[f] };
    const std::size_t Len{ SummaryName.find_last_of( '.' ) };
    if( Len != std::string::npos )
      SummaryName.resize( Len );
    out[f].WriteSummary( SummaryName + "." + TEXT_EXT );
  }
}

/*****************************************************************

 Create synthetic data, overwriting the file specified on the command-line

*****************************************************************/

int main(const int argc, const char *argv[])
{
  static const char pszDefaultExp[] = "2";
  static const char pszDefaultJitter[] = "0.001";
  static const char pszDefaultSamples[] = "100";
  static const char pszDefaultE0[] = "0.8, 0.05";
  static const char pszDefaultA0[] =  "50, 0.5";
  static const char pszDefaultB0[] = "100, 0.3";
  static const char pszDefaultE1[] = "1.2, 0.1";
  static const char pszDefaultA1[] =  "25, 1";
  static const char pszDefaultB1[] ="-200, 5";
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  const std::string MachineName{ MLU::GetHostName() };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"e", CL::SwitchType::Single, pszDefaultExp},
      {"j", CL::SwitchType::Single, pszDefaultJitter},
      {"m", CL::SwitchType::Single, nullptr},
      {"s", CL::SwitchType::Single, pszDefaultSamples},
      {"n", CL::SwitchType::Single, DEF_NSAMPLE},
      {"o", CL::SwitchType::Single, "" },
      {"E0", CL::SwitchType::Single, pszDefaultE0},
      {"A0", CL::SwitchType::Single, pszDefaultA0},
      {"B0", CL::SwitchType::Single, pszDefaultB0},
      {"E1", CL::SwitchType::Single, pszDefaultE1},
      {"A1", CL::SwitchType::Single, pszDefaultA1},
      {"B1", CL::SwitchType::Single, pszDefaultB1},
      {"v", CL::SwitchType::Flag, nullptr},
      {"w", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      bShowUsage = false;
      std::size_t Made{};
      Synth Sy( cl, MachineName );
      Sy.DumpParams( std::cout );
      for( const std::string &Basename : cl.Args )
      {
        try
        {
          Sy.Make( Basename );
          ++Made;
        }
        catch(const std::exception &e)
        {
          std::cerr << "Error: " << e.what() << std::endl;
          iReturn = EXIT_FAILURE;
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> ContractionFile1 [ContractionFile2 ...]\n"
    "Perform a bootstrap of the specified files, where <options> are:\n"
    "-e     Number of exponentials (default: " << pszDefaultExp << ")\n"
    "-j     Correlator jitter, i.e. relative error (default: " << pszDefaultJitter << ")\n"
    "-m     Machine name (default: " << MachineName << ")\n"
    "-s     Number of raw/binned samples (" << pszDefaultSamples << ")\n"
    "-n     Number of bootstrap samples (" DEF_NSAMPLE ")\n"
    "-o     Output prefix\n"
    "--E0   Mean, error (sigma) for E0 (default: " << pszDefaultE0 << ")\n"
    "--A0   Mean, error (sigma) for A0 (default: " << pszDefaultA0 << ")\n"
    "--B0   Mean, error (sigma) for B0 (default: " << pszDefaultB0 << ")\n"
    "--E1   Mean, error (sigma) for E0 (default: " << pszDefaultE1 << ")\n"
    "--A1   Mean, error (sigma) for A0 (default: " << pszDefaultA1 << ")\n"
    "--B1   Mean, error (sigma) for B0 (default: " << pszDefaultB1 << ")\n"
    "Flags:\n"
    "-v     Verbose\n"
    "-w     Warn only if file exists. Default=error\n"
    "--help This message\n";
  }
  return iReturn;
}
