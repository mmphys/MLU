/**
 
 For now this divides arg 1 by arg 2. Might become generic mathematical manipulation

 Source file: Divide.cpp

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
#include "Divide.hpp"

template <typename T>
inline void Divide( MLU::Matrix<T> &Num, const MLU::Matrix<T> Den )
{
  for( std::size_t i = 0; i != Num.size1; ++i )
    for( std::size_t j = 0; j < Num.size2; ++j )
      Num(i,j) = MLU::ComponentDivide( Num(i,j), Den(i,j) );
}

template <typename T>
inline void Divide( MLU::JackBoot<T> &Num, const MLU::JackBoot<T> Den )
{
  for( std::size_t i = Num.idxReplicaMean; i != Num.NumReplicas(); ++i )
    for( std::size_t j = 0; j < Num.extent(); ++j )
      Num(i,j) = MLU::ComponentDivide( Num(i,j), Den(i,j) );
}

Params::Params( const MLU::CommandLine &cl )
: bForceOverwrite{ cl.GotSwitch( "f" ) },
  bSaveHdf5{ cl.GotSwitch( "h" ) },
  InBase { cl.SwitchValue<std::string>("i") },
  OutBase{ cl.SwitchValue<std::string>("o") }
{
  if( cl.GotSwitch( "v" ) )
  {
    std::vector<std::string> v{ MLU::Split( cl.SwitchValue<std::string>( "v" ), MLU::Comma.c_str() ) };
    for( std::string &s : v )
    {
      bool bError{ true };
      std::size_t pos{ s.find_first_of( '=' ) };
      if( pos != std::string::npos && pos && ( s.length() - pos ) > 1 )
      {
        double d{ MLU::FromString<double>( s.substr( pos + 1, s.length() - pos - 1 ) ) };
        s.resize( pos );
        MLU::Trim( s );
        if( !s.empty() )
        {
          VolFactor.emplace( std::make_pair( s, d ) ); // Spatial volume factor is cubed
          bError = false;
        }
      }
      if( bError )
        throw std::runtime_error( "Volume factor \"" + s + "\" invalid" );
    }
    bool bFirst{ true };
    for( VFMap::value_type &pair : VolFactor )
    {
      if( bFirst )
      {
        std::cout << "Volume factors: ";
        bFirst = false;
      }
      else
        std::cout << MLU::CommaSpace;
      std::cout << pair.first << "=" << pair.second;
      pair.second = 1 / std::sqrt( pair.second * pair.second * pair.second );
    }
    if( !bFirst )
      std::cout << MLU::NewLine;
  }
  MLU::MakeAncestorDirs( OutBase );
}

template <typename T>
void Params::Normalise( MLU::JackBoot<T> &jb, const std::vector<int> &opNum,
                   const std::vector<std::string> &opNames )
{
  using Real = typename MLU::JackBoot<T>::Real;
  static constexpr Real One{ static_cast<Real>( 1 ) };
  // First work out what the volume factor is
  Real dFactor{ One };
  for( int i : opNum )
  {
    VFMap::iterator it{ VolFactor.find( opNames[i] ) };
    if( it != VolFactor.end() )
      dFactor *= it->second;
  }
  if( dFactor != One )
  {
    std::cout << "  Volume factor " << dFactor << MLU::NewLine;
    for( std::size_t i = 0; i < jb.NumReplicas(); ++i )
      for( std::size_t j = 0; j < jb.extent(); ++j )
        jb(i,j) *= dFactor;
  }
}

template <typename FoldBoot>
void Params::DivideBoot( FoldBoot &Numerator, const FoldBoot &Denominator )
{
  int NumSamples = Numerator.NumSamples();
  Numerator.IsCompatible( Denominator, &NumSamples, MLU::COMPAT_DISABLE_BASE );
  // Make sure output doesn't exist
  const std::string OutBaseName{ OutBase + Numerator.Name_.NameNoExt + "." };
  const std::string OutFileName{ OutBaseName + DEF_FMT };
  const std::string SummaryName{ OutBaseName + TEXT_EXT };
  if( !bForceOverwrite && ( MLU::FileExists( OutFileName ) || MLU::FileExists( SummaryName ) ) )
    std::cout << " Skip existing " << OutFileName << MLU::NewLine;
  else
  {
    // Divide the Numerator by the denominator
    Divide( Numerator.getData(), Denominator.getData() );
    // If the binned data match, save the original and ratios
    if(   Numerator.NumSamplesBinned()
       && Numerator.NumSamplesBinned() == Denominator.NumSamplesBinned() )
    {
      Divide( Numerator.getBinned(), Denominator.getBinned() );
      if( Numerator.binSize != Denominator.binSize )
        Numerator.binSize = 1;
    }
    else
    {
      Numerator.resizeBinned( 0 );
      Numerator.binSize = 1;
    }
    // If the raw data match, save the original and ratios
    if( Numerator.NumSamplesRaw() && Numerator.NumSamplesRaw() == Denominator.NumSamplesRaw() )
      Divide( Numerator.getRaw(), Denominator.getRaw() );
    else
      Numerator.resizeRaw( 0 );
    // Merge the file lists
    if( Numerator.FileList.size() != Denominator.FileList.size() )
      std::cout << "Warning: sample FileList lengths unequal, appending" << std::endl;
    Numerator.FileList = MLU::ZipperMerge( Numerator.FileList, Denominator.FileList );
    // Update the summaries
    Numerator.MakeCorrSummary();
    if( bSaveHdf5 )
      Numerator.Write( OutFileName, Numerator.Name_.Type.c_str() );
    Numerator.WriteSummary( SummaryName );
  }
}

void Params::DivideFold( Fold &Numerator, const Fold &Denominator )
{
  // Merge the Bootstrap lists
  if( Numerator.BootstrapList.size() != Denominator.BootstrapList.size() )
    std::cout << "Warning: sample BootstrapList lengths unequal, appending" << std::endl;
  Numerator.BootstrapList = MLU::ZipperMerge( Numerator.BootstrapList, Denominator.BootstrapList );
  DivideBoot( Numerator, Denominator );
}

void Params::ReadNumerator( const std::string &FileName )
{
  static const char sNumerator[] = "Numerator ";
  std::vector<std::string> opNames;
  if( bFold )
  {
    fNumerator.Read( FileName, sNumerator, &opNames );
    Normalise( fNumerator.getData(), fNumerator.Name_.op, opNames );
    DivideFold( fNumerator, fDenominator );
  }
  else
  {
    bNumerator.Read( FileName, sNumerator, &opNames );
    Normalise( bNumerator.getData(), bNumerator.Name_.op, opNames );
    DivideBoot( bNumerator, bDenominator );
  }
}

void Params::Run( const MLU::CommandLine &cl, bool &bShowUsage )
{
  std::vector<std::string> Filename{ MLU::glob( cl.Args.begin(), cl.Args.end(), InBase.c_str() ) };
  if( Filename.empty() )
    throw std::runtime_error( "No input files found" );
  bShowUsage = false;
  bool bDoDenominator{ true };
  for( std::size_t FileNum = 0; FileNum < Filename.size(); ++FileNum )
  {
    // Read the first entry, either as a folded correlator, or as a bootstrap
    if( FileNum == 0 )
    {
      std::string GroupName;
      std::vector<std::string> opNames;
      try
      {
        fDenominator.Read( Filename[FileNum], nullptr, &opNames, &GroupName );
        Normalise( fDenominator.getData(), fDenominator.Name_.op, opNames );
        bFold = true;
      } catch (const std::exception &e) {
      }
      if( !bFold )
      {
        GroupName.clear();
        opNames.clear();
        bDenominator.Read( Filename[FileNum], nullptr, &opNames, &GroupName );
        Normalise( bDenominator.getData(), bDenominator.Name_.op, opNames );
      }
      std::cout << "Denominator " << GroupName << MLU::Space << Filename[FileNum] << MLU::NewLine;
    }
    else
    {
      if( !Filename[FileNum].compare( Filename[0] ) )
        bDoDenominator = false;
      ReadNumerator( Filename[FileNum] );
    }
  }
  if( bDoDenominator )
    ReadNumerator( Filename[0] );
}

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"v", CL::SwitchType::Single, nullptr },
      {"f", CL::SwitchType::Flag, nullptr},
      {"h", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && !cl.Args.empty() )
    {
      Params par{ cl };
      par.Run( cl, bShowUsage );
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name << " <options> Denominator Numerator1 [Numerator2 [...]]\n"
    "Write each Numerator / Denominator to Output.\n"
    "If Denominator is not in the list of numerators, it is also written\n"
    "allows for example: " << cl.Name << " File3.h5 '*'.h5\n"
    "<options> are:\n"
    "-i   Input  filename prefix\n"
    "-o   Output filename prefix\n"
    "-v   Volume factors, typically spatial lattice extent, e.g. 'g5P=1,g5W=32'\n"
    "Flags:\n"
    "-f     Force overwrite if output already exists (default skip)\n"
    "-h     Save .h5 output as well (default .txt only)\n"
    "--help This message\n";
  }
  return iReturn;
}
