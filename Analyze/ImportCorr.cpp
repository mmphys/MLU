/**

 Mike's lattice QCD utilities: Import correlators
 
 Source file: ImportCorr.cpp
 
 Copyright (C) 2022
 
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
**/

#include "ImportCorr.hpp"

// The correlators I want to import - in order

const std::vector<CorrInfo> Importer::corrInfo{
  { "pp_SSLL", "pLS", "pLS", Common::Parity::Even },
  { "pp_SLLL", "pLS", "pLL", Common::Parity::Even },
  { "ap_SSLL", "aLS", "pLS", Common::Parity::Odd },
  { "ap_SLLL", "aLS", "pLL", Common::Parity::Odd },
  { "aa_SSLL", "aLS", "aLS", Common::Parity::Even },
  { "aa_SLLL", "aLS", "aLL", Common::Parity::Even },
};


/*****************************************************************

 Constructor - read options / initialise

*****************************************************************/

Importer::Importer( const Common::CommandLine &cl )
: Nt{ cl.SwitchValue<int>( "nt" ) },
  NtOut{ Nt / 2 + 1 },
  Seed{ cl.SwitchValue<Common::SeedType>( "r" ) },
  outStem{ cl.SwitchValue<std::string>( "o" ) },
  GroupName{ cl.SwitchValue<std::string>( "g" ) }
{
}

/*****************************************************************

 Read in an array / vector of stuff

*****************************************************************/

template <typename T>
void Importer::ReadCorrVector( ::H5::Group &g, const std::string &Name, std::vector<T> &v )
{
  std::cout << "Reading vector from " << Name << Common::NewLine;
  ::H5::Group gName{ g.openGroup( Name ) };
  v.reserve( corrInfo.size() );
  for( const CorrInfo &i : corrInfo )
  {
    v.emplace_back();
    std::cout << "  Reading element from " << i.id << Common::NewLine;
    Common::H5::ReadVector( gName, i.id, v.back() );
  }
}

/*****************************************************************

 Import data from HDF5

*****************************************************************/

void Importer::Import( const std::string &Filename )
{
  // Read source data
  {
    ::H5::H5File f;
    ::H5::Group gRoot;
    Common::H5::OpenFileGroup( f, gRoot, Filename, "Importing ", &GroupName );
    ReadCorrVector( gRoot, "rawdata", vRawData );
    ::H5::Group gFit{ gRoot.openGroup( "fitinput" ) };
    ReadCorrVector( gFit, "fitranges", vFitRanges );
    static const std::string sCentral{ "central" };
    static const std::string sBootstraps{ "Bootstraps" };
    ::H5::Group g{ gFit.openGroup( "binneddata" ) };
    Common::H5::ReadVector( g, sCentral, vBinnedCentral );
    Common::H5::ReadMatrix( g, sBootstraps, mBinnedData );
    g.close();
    g = gFit.openGroup( "unbinnneddata" ); // Spelling error "nnn" is in the source data
    Common::H5::ReadVector( g, sCentral, vUnbinnedCentral );
    Common::H5::ReadMatrix( g, sBootstraps, mUnbinnedData );
  }
  std::cout << "Timeslices" << Common::NewLine;
  std::size_t TotalTimeslices{ 0 };
  bool bBadTimeslice{ false };
  for( std::size_t i = 0; i < corrInfo.size(); ++i )
  {
    TotalTimeslices += vFitRanges[i].size();
    std::cout << Common::Space << corrInfo[i].id << " (" << corrInfo[i].opSrc << ", "
              << corrInfo[i].opSnk << "):";
    for( int t : vFitRanges[i] )
    {
      std::cout << Common::Space << t;
      if( t < 0 || t > Nt / 2 )
        bBadTimeslice = true;
    }
    std::cout << Common::NewLine;
    Common::NoDuplicates( vFitRanges[i], "Timeslice", 1 );
  }
  static const std::string Suffix{ " timeslices != " + std::to_string( TotalTimeslices ) };
  if( TotalTimeslices != vBinnedCentral.size )
    throw std::runtime_error( "Binned central values " + std::to_string( vBinnedCentral.size ) + Suffix );
  if( TotalTimeslices != mBinnedData.size2 )
    throw std::runtime_error( "Binned data " + std::to_string( mBinnedData.size2 ) + Suffix );
  if( TotalTimeslices != vUnbinnedCentral.size )
    throw std::runtime_error( "Unbinned central values " + std::to_string( vUnbinnedCentral.size ) + Suffix );
  if( TotalTimeslices != mUnbinnedData.size2 )
    throw std::runtime_error( "Unbinned data " + std::to_string( mUnbinnedData.size2 ) + Suffix );
  if( bBadTimeslice )
    throw std::runtime_error( "Timeslice < 0 || > " + std::to_string( Nt / 2 ) );
}

/*****************************************************************

 Take single row of input data and spread over correlators

*****************************************************************/

void Importer::SpreadData( int idx, std::vector<Fold> &out, Scalar * &pSource )
{
  for( int i = 0; i < out.size(); ++i )
  {
    Scalar * C{ out[i][idx] };
    for( int t = 0; t < out[i].Nt(); ++t )
      C[t] = 0;
    for( int j = 0; j < vFitRanges[i].size(); ++j )
      C[vFitRanges[i][j]] = *pSource++;
  }
}

void Importer::SpreadDataRaw( std::vector<Fold> &out, const Matrix *pRawData )
{
  if( pRawData )
  {
    const Scalar *pSource{ pRawData->data };
    std::vector<Scalar *> C( out.size() );
    const int NumRawSamples{ static_cast<int>( pRawData->size1 ) };
    for( int f = 0; f < out.size(); ++f )
    {
      C[f] = out[f].resizeRaw( NumRawSamples );
      out[f].bRawBootstrap = true;
    }
    for( int s = 0; s < NumRawSamples; ++s )
      for( int f = 0; f < out.size(); ++f )
      {
        for( int t = 0; t < out[f].Nt(); ++t )
          C[f][t] = 0;
        for( int j = 0; j < vFitRanges[f].size(); ++j )
          C[f][vFitRanges[f][j]] = *pSource++;
         C[f] += out[f].Nt();
      }
  }
}

/*****************************************************************

 Import data from HDF5

*****************************************************************/

void Importer::Write( const std::string &Base, Matrix &Data, Vector &Central, const Matrix *pRawData )
{
  const int NumSamples{ static_cast<int>( Data.size1 ) };
  std::cout << "Writing " << Base << Common::NewLine;
  std::vector<Fold> out( corrInfo.size() );
  for( std::size_t f = 0; f < corrInfo.size(); ++f )
  {
    out[f].resize( NumSamples, NtOut );
    out[f].NtUnfolded = Nt;
    out[f].parity = corrInfo[f].parity;
    out[f].sign = Common::Sign::Positive;
  }
  SpreadData( Fold::idxCentral, out, Central.data );
  for( int s = 0; s < NumSamples; ++s )
    SpreadData( s, out, Data.data );
  SpreadDataRaw( out, pRawData );
  for( std::size_t f = 0; f < corrInfo.size(); ++f )
  {
    out[f].MakeCorrSummary( nullptr );
    std::string Filename{ Base };
    Filename.append( 1, '_' );
    Filename.append( corrInfo[f].opSnk );
    Filename.append( 1, '_' );
    Filename.append( corrInfo[f].opSrc );
    std::string FullName{ Common::MakeFilename( Filename, Common::sFold, Seed, DEF_FMT ) };
    Common::MakeAncestorDirs( FullName );
    out[f].Write( FullName );
    out[f].WriteSummary( Common::MakeFilename( Filename, Common::sFold, Seed, TEXT_EXT ) );
  }
}

/*****************************************************************

 Import data from HDF5

*****************************************************************/

void Importer::Run( const std::string &Filename )
{
  Import( Filename );

  // Get the base of the filename
  std::string Base{ outStem };
  Base.append( Filename );
  const std::size_t LastSlash{ Base.find_last_of( '/' ) };
  std::size_t pos{ Base.find_last_of( '.' ) };
  if( pos != std::string::npos && ( LastSlash == std::string::npos || pos > LastSlash) )
    Base.resize( pos );
  Base.append( 1, '/' );

  // Now write out the binned, then unbinned data
  pos = Base.length();
  Base.append( "binned" );
  Write( Base, mBinnedData, vBinnedCentral, &mUnbinnedData );
  Base.resize( pos );
  Base.append( "unbinned" );
  Write( Base, mUnbinnedData, vUnbinnedCentral, nullptr );
}

/*****************************************************************

 Import correlators

*****************************************************************/

int main(const int argc, const char *argv[])
{
  const char pszDefaultGroupName[] = "/C0/sh/0.51";
  const char pszDefaultNt[] = "96";
  const char pszDefaultSeed[] = "1";
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"g", CL::SwitchType::Single, pszDefaultGroupName},
      {"o", CL::SwitchType::Single, ""},
      {"r", CL::SwitchType::Single, pszDefaultSeed},
      {"nt", CL::SwitchType::Single, pszDefaultNt},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() == 1 )
    {
      bShowUsage = false;
      Importer I( cl );
      I.Run( cl.Args[0] );
      bShowUsage = false;
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) <<
    "Import correlators.\n"
    "usage: " << cl.Name << " <options> Import_Filename\n"
    "Options:\n"
    "-g      Group name, default " << pszDefaultGroupName << "\n"
    "-o      Output file base\n"
    "-r      Random number seed, default " << pszDefaultSeed << "\n"
    "--nt    Timeslices on ensemble, default" << pszDefaultNt << "\n"
    "Flags:\n"
    "--help  This message\n";
  }
  return iReturn;
}
