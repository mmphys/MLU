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

void Divide( Scalar * pNum, const Scalar * pDen, int NumSamples, int Nt )
{
  for( int Sample = 0; Sample < NumSamples; ++Sample )
    for( int t = 0; t < Nt; ++t )
      *pNum++ /= *pDen++;
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
      //{"i1", CL::SwitchType::Single, "" },
      //{"i2", CL::SwitchType::Single, "" },
      //{"o", CL::SwitchType::Single, "" },
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() == 3 )
    {
      bShowUsage = false;
      // Read the list of fits I've chosen to use
      //const std::string InBase1{ cl.SwitchValue<std::string>("i1") };
      //const std::string InBase2{ cl.SwitchValue<std::string>("i2") };
      //const std::string OutBase{ cl.SwitchValue<std::string>("o") };
      Fold Numerator, Denominator;
      Numerator.Read( cl.Args[0], "Numerator " );
      Denominator.Read( cl.Args[1], "Denominator " );
      int NumSamples = Numerator.NumSamples();
      Numerator.IsCompatible( Denominator, &NumSamples, Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_CONFIG_COUNT );
      // Make sure output doesn't exist
      const std::string &OutBaseName{ cl.Args[2] };
      const std::string OutFileName{ Common::MakeFilename( OutBaseName, Common::sFold, Numerator.Name_.Seed, DEF_FMT ) };
      if( Common::FileExists( OutFileName ) )
        throw std::runtime_error( "Output exists: " + OutFileName );
      // Divide the Numerator by the denominator
      const int Nt{ Numerator.Nt() };
      Divide( Numerator[Fold::idxCentral], Denominator[Fold::idxCentral], NumSamples + 1, Nt );
      // If the binned data match, save the original and ratios
      if( Numerator.NumSamplesBinned() && Numerator.NumSamplesBinned() == Denominator.NumSamplesBinned() )
      {
        // Copy the binned data from numerator and denominator into the raw data area of the file to save
        const std::size_t BlockSize{ static_cast<std::size_t>( Numerator.NumSamplesBinned() ) * Nt };
              Scalar * pDst{ Numerator.resizeRaw( 2 * Numerator.NumSamplesBinned() ) };
        const Scalar * pSrc{ Numerator.getBinned() };
        for( std::size_t i = 0; i < BlockSize; ++i )
          *pDst++ = *pSrc++;
        pSrc = Denominator.getBinned();
        for( std::size_t i = 0; i < BlockSize; ++i )
          *pDst++ = *pSrc++;
        Divide( Numerator.getBinned(), Denominator.getBinned(), Denominator.NumSamplesBinned(), Nt );
      }
      else
      {
        Numerator.resizeRaw( 0 );
        Numerator.resizeBinned( 0 );
      }
      if( Numerator.binSize != Denominator.binSize )
        Numerator.binSize = 1;
      // Merge the file lists
      if( Numerator.FileList.size() != Denominator.FileList.size() )
        std::cout << "Warning: sample FileList lengths unequal, appending" << std::endl;
      Numerator.FileList = Common::ZipperMerge( Numerator.FileList, Denominator.FileList );
      // Merge the Bootstrap lists
      if( Numerator.BootstrapList.size() != Denominator.BootstrapList.size() )
        std::cout << "Warning: sample BootstrapList lengths unequal, appending" << std::endl;
      Numerator.BootstrapList = Common::ZipperMerge( Numerator.BootstrapList, Denominator.BootstrapList );
      // Update the summaries
      Numerator.MakeCorrSummary( nullptr );
      Numerator.Write( OutFileName, Common::sFold.c_str() );
      Numerator.WriteSummary( Common::MakeFilename( OutBaseName, Common::sFold, Numerator.Name_.Seed, TEXT_EXT ) );
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name << " <options> Numerator Denominator Output\n"
    "Write Numerator / Denominator to Output. <options> are:\n"
    //"--i1   Input1 filename prefix\n"
    //"--i2   Input2 filename prefix\n"
    //"-o     Output filename prefix\n"
    "Flags:\n"
    "--help This message\n";
  }
  return iReturn;
}
