/*************************************************************************************
 
 Utility for creating statistics for bootstrap sample - ready for GNUplot
 
 Source file: bootstat.cpp
 
 Copyright (C) 2019
 
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

#include <stdio.h>
#include <typeinfo>
#include "Common.hpp"

using C = std::complex<double>;

bool debug( int NumSamples, int Nt )
{
  Common::Sample<C> a( NumSamples, Nt );
  int z = 0;
  for( int x = -1; x < NumSamples; x++ )
    for( int y = 0; y < Nt; y++ )
      a[x][y] = z++;
  C * p = a[0];
  for( int i = 0; i < ( NumSamples + 1 ) * Nt; i++ )
    std::cout << i << " = " << *p++ << "\n";
  return true;
}

int main(int argc, const char *argv[])
{
  // We're not using C-style I/O
  std::ios_base::sync_with_stdio(false);
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"o", CL::SwitchType::Single, "" },
      {"a", CL::SwitchType::Single, ""},
      {"x", CL::SwitchType::Multiple, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    std::string outStem{ cl.SwitchValue<std::string>( "o" ) };
    std::vector<Common::Gamma::Algebra> Alg = Common::ArrayFromString<Common::Gamma::Algebra>( cl.SwitchValue<std::string>( "a" ) );
    std::vector<std::string> Ignore{ cl.SwitchStrings( "x" ) };
    // If there are files specified on the command line, process each file
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      bShowUsage = false;
      std::size_t Count = 0;
      for( const std::string &Filename : cl.Args )
      {
        std::size_t iIgnore = 0;
        while( iIgnore < Ignore.size() && Ignore[iIgnore].compare(Filename) )
          iIgnore++;
        if( iIgnore < Ignore.size() )
          std::cout << "Ignoring " << Filename << std::endl;
        else if( !Common::FileExists(Filename))
          std::cout << "Error: " << Filename << " doesn't exist" << std::endl;
        else
        {
          std::string GroupName;
          std::vector<Common::Gamma::Algebra> ThisAlg{ Alg };
          Common::CorrelatorFileC InFile( Filename, GroupName, ThisAlg );
          InFile.WriteSummary(outStem, ThisAlg );
          Count++;
        }
      }
      std::cout << "Summaries written for " << Count << " correlators" << std::endl;
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> ContractionFile1 [ContractionFile2 ...]\n"
    "Save correlator in format ready for GNUPlot, where <options> are:\n"
    "-o     Output prefix\n"
    "-a     list of gamma algebras we're interested in\n"
    "-x     eXclude file (may be repeated)\n"
    "--help This message\n";
  }
  return iReturn;
}
