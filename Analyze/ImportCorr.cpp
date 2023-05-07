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
  { "ap_SLLL", "aLL", "pLS", Common::Parity::Odd },
  { "aa_SSLL", "aLS", "aLS", Common::Parity::Even },
  { "aa_SLLL", "aLS", "aLL", Common::Parity::Even },
};


/*****************************************************************

 Constructor - read options / initialise

*****************************************************************/

Importer::Importer( const Common::CommandLine &cl )
: Seed{ cl.SwitchValue<Common::SeedType>( "r" ) },
  outStem{ cl.SwitchValue<std::string>( "o" ) },
  bDebug{ cl.GotSwitch( "debug" ) },
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
  v.clear();
  v.reserve( corrInfo.size() );
  for( const CorrInfo &i : corrInfo )
  {
    v.emplace_back();
    std::cout << "  Reading element from " << i.id << Common::NewLine;
    Common::H5::ReadVector( gName, i.id, v.back() );
  }
}

/*****************************************************************

 Dump part of a matrix for debugging

*****************************************************************/

void Importer::DebugDump( const std::vector<Matrix> &v ) const
{
  if( bDebug )
  {
    const std::size_t confs{ v.size() >= 3 ? 3 : v.size() };
    for( int c = 0; c < confs; ++c )
    {
      const Matrix &m{ v[c] };
      const std::size_t rows{ m.size1 >= 2 ? 2 : m.size1 };
      const std::size_t cols{ m.size2 >= 5 ? 5 : m.size2 };
      for( int i = 0; i < rows; ++i )
      {
        if( i == 0 )
          std::cout << "Conf " << c << ":";
        else
          std::cout << std::string( 7, ' ' );
        for( int j = 0; j < cols; ++j )
          std::cout << Common::Space << m( i, j );
        std::cout << Common::NewLine;
      }
    }
  }
}

void Importer::DebugDump( const std::vector<std::vector<Matrix>> &v ) const
{
  if( bDebug )
  {
    std::cout << std::setprecision( std::numeric_limits<Scalar>::max_digits10 );
    for( std::size_t f = 0; f < v.size(); ++f )
    {
      std::cout << Common::NewLine << corrInfo[f].id << Common::NewLine;
      DebugDump( v[f] );
    }
  }
}

/*****************************************************************

 Import data from HDF5

*****************************************************************/

void Importer::ReadInput( const std::string &Filename )
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
  DebugDump( vRawData );
  if( vRawData[0][0].size2 > std::numeric_limits<int>::max() || vRawData[0][0].size2 & 1 )
    throw std::runtime_error( "Raw data Nt=" + std::to_string( vRawData[0][0].size2 ) + " invalid" );
  Nt = static_cast<int>( vRawData[0][0].size2 );
  NtOut = Nt / 2 + 1;
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
  const std::string Suffix{ " timeslices != " + std::to_string( TotalTimeslices ) };
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
  const std::string Suffix2{ " rows >= int_max " + std::to_string( std::numeric_limits<int>::max() ) };
  if( mBinnedData.size1 > std::numeric_limits<int>::max() )
    throw std::runtime_error( "Binned data " + std::to_string( mBinnedData.size1 ) + Suffix2 );
  if( mUnbinnedData.size1 > std::numeric_limits<int>::max() )
    throw std::runtime_error( "Binned data " + std::to_string( mUnbinnedData.size1 ) + Suffix2 );
}

/*****************************************************************

 Take single row of input data and spread over correlators

*****************************************************************/

void Importer::SpreadData( std::vector<Fold> &out, const Scalar *pSource, int idx, int NumSamples )
{
  std::vector<Scalar *> C( out.size() );
  for( int s = 0; s < NumSamples; ++s, ++idx )
  {
    for( int f = 0; f < out.size(); ++f )
    {
      for( int t = 0; t < out[f].Nt(); ++t )
        out[f](idx,t) = 0;
      for( int j = 0; j < vFitRanges[f].size(); ++j )
        out[f](idx,vFitRanges[f][j]) = *pSource++;
    }
  }
}

void Importer::SpreadDataRawBoot( std::vector<Fold> &out, const Matrix &BinnedData )
{
  std::vector<Matrix *> C( out.size() );
  const int NumBinnedSamples{ static_cast<int>( BinnedData.size1 ) };
  for( int f = 0; f < out.size(); ++f )
  {
    C[f] = &out[f].resizeRaw( NumBinnedSamples );
    if( !f )
      throw std::runtime_error( "raw bootstrap support removed May 2023. Fix next line" );
    // TODO: out[f].bRawBootstrap = true;
  }
  for( int s = 0; s < NumBinnedSamples; ++s )
    for( int f = 0; f < out.size(); ++f )
    {
      for( int t = 0; t < out[f].Nt(); ++t )
        (*C[f])(s,t) = 0;
      for( int j = 0; j < vFitRanges[f].size(); ++j )
        (*C[f])(s,vFitRanges[f][j]) = BinnedData(s,j);
    }
}

// This folds the raw data, making positive or preserving sign
void Importer::SaveRawData( std::vector<Fold> &out, bool bPreserveSign )
{
  for( std::size_t f = 0; f < out.size(); ++f )
  {
    const Scalar Multiplier{ corrInfo[f].parity == Common::Parity::Odd ? -1. : 1. };
    std::string sName{ corrInfo[f].id };
    sName.append( "_config_" );
    const std::size_t SizeConf{ sName.length() };
    static const std::string sErr{ "Raw data error: " };
    if( out[f].Nt() != NtOut )
      throw std::runtime_error( "Bug: bad NtOut" );
    if( vRawData[f].size() < 1 || vRawData[f].size() > std::numeric_limits<int>::max() )
      throw std::runtime_error( sErr + std::to_string( vRawData[f].size() ) + " configurations" );
    if( vRawData[f][0].size1 < 1 || vRawData[f][0].size1 > std::numeric_limits<int>::max() )
      throw std::runtime_error( sErr + std::to_string( vRawData[f][0].size1 ) + " timeslices" );
    if( vRawData[f][0].size2 != Nt )
      throw std::runtime_error( sErr + "NtIn " + std::to_string( vRawData[f][0].size2 ) );
    const int NumConfigs   { static_cast<int>( vRawData[f].size() ) };
    const int NumTimeslices{ static_cast<int>( vRawData[f][0].size1 ) };
    const int NtIn         { static_cast<int>( vRawData[f][0].size2 ) };
    const int NumSamplesRaw{ NumConfigs * NumTimeslices };
    if( NumSamplesRaw / NumConfigs != NumTimeslices )
      throw std::runtime_error( sErr + std::to_string( NumConfigs ) + " configs and "
                               + std::to_string( NumTimeslices ) + " timeslices" );
    const bool bNegate{ !bPreserveSign && vRawData[f][0]( 0, Nt >> 2 ) < 0 };
    out[f].binSize = NumTimeslices;
    out[f].sign = bNegate ? Common::Sign::Negative : Common::Sign::Positive;
    out[f].ConfigCount.resize( NumConfigs );
    out[f].SampleSize = NumConfigs;
    out[f].binSize = NumTimeslices;
    out[f].FileList.reserve( NumSamplesRaw );
    Matrix &mRawData{ out[f].resizeRaw( NumSamplesRaw ) };
    for( int c = 0; c < NumConfigs; ++c )
    {
      sName.resize( SizeConf );
      sName.append( std::to_string( c ) );
      sName.append( "_t_" );
      const std::size_t SizeT{ sName.length() };
      const Matrix &m{ vRawData[f][c] };
      if( m.size1 != NumTimeslices || m.size2 != NtIn )
        throw std::runtime_error( "Bug: ragged matrices - not possible if read from hdf5 cube" );
      out[f].ConfigCount[c].Config = c;
      out[f].ConfigCount[c].Count = NumTimeslices;
      for( int tSlice = 0; tSlice < NumTimeslices; ++tSlice )
      {
        sName.resize( SizeT );
        sName.append( std::to_string( tSlice ) );
        out[f].FileList.push_back( sName );
        for( int t = 0; t < NtOut; ++t )
        {
          Scalar z;
          if( t == 0 || t == NtOut - 1 )
            z = m( tSlice, t );
          else
            z = 0.5 * ( m( tSlice, t ) + Multiplier * m( tSlice, Nt - t ) );
          if( bNegate )
            z = -z;
          mRawData(c,t) = z;
        }
      }
    }
  }
}

/*****************************************************************

 Import data from HDF5

*****************************************************************/

void Importer::Write( const std::string &Base, bool bPreserveSign )
{
  const int NumSamples{ static_cast<int>( mBinnedData.size1 ) };
  std::cout << "Writing " << Base << Common::NewLine;
  Common::MakeAncestorDirs( Base );
  std::vector<Fold> out( corrInfo.size() );
  for( std::size_t f = 0; f < corrInfo.size(); ++f )
  {
    out[f].resize( NumSamples, NtOut );
    out[f].NtUnfolded = Nt;
    out[f].parity = corrInfo[f].parity;
  }
  SpreadData( out, vBinnedCentral.data, Fold::idxCentral, 1 );
  SpreadData( out, mBinnedData.data, 0, NumSamples );
  SaveRawData( out, bPreserveSign );
  for( std::size_t f = 0; f < corrInfo.size(); ++f )
  {
    out[f].BinAuto();
    out[f].MakeCorrSummary();
  }
  std::string Filename{ Base };
  const std::size_t BaseLen{ Filename.length() };
  for( int i = 0; i < 3; ++i )
  {
    Filename.resize( BaseLen );
    if( i == 0 )
      Filename.append( "binned" );
    else if( i == 1 )
    {
      SpreadDataRawBoot( out, mUnbinnedData );
      Filename.append( "binned_rawboot" );
    }
    else
    {
      SpreadData( out, vUnbinnedCentral.data, Fold::idxCentral, 1 );
      SpreadData( out, mUnbinnedData.data, 0, NumSamples );
      SaveRawData( out, bPreserveSign );
      for( std::size_t f = 0; f < corrInfo.size(); ++f )
      {
        out[f].BinFixed( 1 );
        out[f].ConfigCount.resize( out[f].NumSamplesRaw() );
        for( int j = 0; j < out[f].NumSamplesRaw(); ++j )
        {
          out[f].ConfigCount[j].Config = j;
          out[f].ConfigCount[j].Count = 1;
        }
        out[f].resizeRaw( 0 );
        out[f].MakeCorrSummary();
      }
      Filename.append( "unbinned" );
    }
    if( bPreserveSign )
      Filename.append( 1, 's' );
    Filename.append( 1, '_' );
    const std::size_t NextLen{ Filename.length() };
    for( std::size_t f = 0; f < corrInfo.size(); ++f )
    {
      Filename.resize( NextLen );
      Filename.append( corrInfo[f].opSnk );
      Filename.append( 1, '_' );
      Filename.append( corrInfo[f].opSrc );
      std::string FullName{ Common::MakeFilename( Filename, Common::sFold, Seed, DEF_FMT ) };
      std::cout << Common::Space << out[f].NumSamplesRaw() << Common::Space << FullName << Common::NewLine;
      out[f].Write( FullName );
      out[f].WriteSummary( Common::MakeFilename( Filename, Common::sFold, Seed, TEXT_EXT ) );
    }
  }
}

/*****************************************************************

 Import data from HDF5

*****************************************************************/

void Importer::Import( const std::string &Filename, bool bPreserveSign )
{
  ReadInput( Filename );

  // Get the base of the filename
  std::string Base{ outStem };
  Base.append( Filename );
  const std::size_t LastSlash{ Base.find_last_of( '/' ) };
  std::size_t pos{ Base.find_last_of( '.' ) };
  if( pos != std::string::npos && ( LastSlash == std::string::npos || pos > LastSlash) )
    Base.resize( pos );
  Base.append( 1, '/' );
  Write( Base, bPreserveSign );
}

/*****************************************************************

 Import correlators

*****************************************************************/

int main(const int argc, const char *argv[])
{
  const char pszDefaultGroupName[] = "/C0/sh/0.51";
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
      {"i", CL::SwitchType::Single, ""},
      {"o", CL::SwitchType::Single, ""},
      {"r", CL::SwitchType::Single, pszDefaultSeed},
      {"debug", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      bShowUsage = false;
      Importer I( cl );
      std::size_t Count{ 0 };
      for( std::string &Args : cl.Args )
      {
        bool bOK{ true };
        bool bPreserveSign{ false };
        std::string Files{ Common::ExtractToSeparator( Args ) };
        if( Args.length() )
        {
          if( Args.length() == 1 && std::toupper( Args[0] ) == 'S' )
            bPreserveSign = true;
          else
          {
            iReturn = EXIT_FAILURE;
            std::cerr << "Error: Bad argument \"" << Args << "\". Ignoring " << Files << Common::NewLine;
            bOK = false;
          }
        }
        if( bOK )
        {
          for( auto &File : Common::glob( &Files, &Files + 1, cl.SwitchValue<std::string>("i").c_str() ) )
          {
            try
            {
              I.Import( File, bPreserveSign );
              Count++;
            }
            catch(const std::exception &e)
            {
              std::cerr << "Error: " << e.what() << std::endl;
              iReturn = EXIT_FAILURE;
            }
          }
        }
      }
      std::cout << Count << " files imported.\n";
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
    "usage: " << cl.Name << " <options> File[,Arg]...\n"
    "Args:\n"
    "s       Preserve the sign of sinh correlators"
    "Options:\n"
    "-g      Group name, default " << pszDefaultGroupName << "\n"
    "-i      Input file base\n"
    "-o      Output file base\n"
    "-r      Random number seed, default " << pszDefaultSeed << "\n"
    "Flags:\n"
    "--debug Dump a portion of the raw data\n";
    "--help  This message\n";
  }
  return iReturn;
}
