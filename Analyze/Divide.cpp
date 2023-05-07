/*************************************************************************************
 
 For now this divides arg 1 by arg 2. Might become generic mathematical manipulation
 Source file: Divide.cpp
 Copyright (C) 2021
 Author: Michael Marshall<Michael.Marshall@ed.ac.uk>
 
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

#include "Divide.hpp"

template <typename T> typename std::enable_if<std::is_floating_point<T>::value>::type
Divide( T * pNum, const T * pDen, int NumSamples, int Nt )
{
  for( int Sample = 0; Sample < NumSamples; ++Sample )
    for( int t = 0; t < Nt; ++t )
      *pNum++ /= *pDen++;
}

template <typename T> typename std::enable_if<std::is_floating_point<T>::value>::type
Divide( std::complex<T> * pNum, const std::complex<T> * pDen, int NumSamples, int Nt )
{
  Divide( reinterpret_cast<T *>( pNum ), reinterpret_cast<const T *>( pDen ), NumSamples, Nt * 2 );
}

Params::Params( const Common::CommandLine &cl )
: bForceOverwrite{ cl.GotSwitch( "f" ) },
  bSaveHdf5{ cl.GotSwitch( "h" ) },
  InBase { cl.SwitchValue<std::string>("i") },
  OutBase{ cl.SwitchValue<std::string>("o") }
{
  if( cl.GotSwitch( "v" ) )
  {
    std::vector<std::string> v{ Common::Split( cl.SwitchValue<std::string>( "v" ), Common::Comma.c_str() ) };
    for( std::string &s : v )
    {
      bool bError{ true };
      std::size_t pos{ s.find_first_of( '=' ) };
      if( pos != std::string::npos && pos && ( s.length() - pos ) > 1 )
      {
        double d{ Common::FromString<double>( s.substr( pos + 1, s.length() - pos - 1 ) ) };
        s.resize( pos );
        Common::Trim( s );
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
        std::cout << Common::CommaSpace;
      std::cout << pair.first << "=" << pair.second;
      pair.second = 1 / std::sqrt( pair.second * pair.second * pair.second );
    }
    if( !bFirst )
      std::cout << Common::NewLine;
  }
  Common::MakeAncestorDirs( OutBase );
}

template <typename T>
void Params::Normalise( Common::JackBoot<T> &jb, const std::vector<int> &opNum,
                   const std::vector<std::string> &opNames )
{
  using Real = typename Common::JackBoot<T>::Real;
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
    std::cout << "  Volume factor " << dFactor << Common::NewLine;
    for( std::size_t i = 0; i < jb.NumReplicas(); ++i )
      for( std::size_t j = 0; j < jb.extent(); ++j )
        jb(i,j) *= dFactor;
  }
}

template <typename FoldBoot>
void Params::DivideBoot( FoldBoot &Numerator, const FoldBoot &Denominator )
{
  int NumSamples = Numerator.NumSamples();
  Numerator.IsCompatible( Denominator, &NumSamples, Common::COMPAT_DISABLE_BASE );
  // Make sure output doesn't exist
  const std::string OutBaseName{ OutBase + Numerator.Name_.NameNoExt + "." };
  const std::string OutFileName{ OutBaseName + DEF_FMT };
  const std::string SummaryName{ OutBaseName + TEXT_EXT };
  if( !bForceOverwrite && ( Common::FileExists( OutFileName ) || Common::FileExists( SummaryName ) ) )
    std::cout << " Skip existing " << OutFileName << Common::NewLine;
  else
  {
    // Divide the Numerator by the denominator
    const int Nt{ Numerator.Nt() };
    Divide( Numerator[Fold::idxCentral], Denominator[Fold::idxCentral], NumSamples + 1, Nt );
    // If the binned data match, save the original and ratios
    if( Numerator.NumSamplesBinned() && Numerator.NumSamplesBinned() == Denominator.NumSamplesBinned() )
    {
      Divide( Numerator.getBinned(), Denominator.getBinned(),
             Numerator.NumSamplesBinned(), Nt );
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
    {
      Divide( Numerator.getRaw(), Denominator.getRaw(),
             Numerator.NumSamplesRaw(), Nt );
    }
    else
      Numerator.resizeRaw( 0 );
    // Merge the file lists
    if( Numerator.FileList.size() != Denominator.FileList.size() )
      std::cout << "Warning: sample FileList lengths unequal, appending" << std::endl;
    Numerator.FileList = Common::ZipperMerge( Numerator.FileList, Denominator.FileList );
    // Update the summaries
    Numerator.MakeCorrSummary( nullptr );
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
  Numerator.BootstrapList = Common::ZipperMerge( Numerator.BootstrapList, Denominator.BootstrapList );
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

void Params::Run( const Common::CommandLine &cl, bool &bShowUsage )
{
  std::vector<std::string> Filename{ Common::glob( cl.Args.begin(), cl.Args.end(), InBase.c_str() ) };
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
      std::cout << "Denominator " << GroupName << Common::Space << Filename[FileNum] << Common::NewLine;
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
  using CL = Common::CommandLine;
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
