/*************************************************************************************
 
 Compare multiple gauge fields
 Source file: GaugeCmp.cpp
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

#include <stdio.h>
//#include <MLU/Common.hpp>
#include <omp.h>
#include <Grid/Grid.h>

using namespace Grid;

//using Grid::operator<<;
//using Grid::operator>>;

// Compare multiple Gauge files to the first
void GaugeCompare( const std::vector<std::string> &File )
{
  // Double precision grids
  GridCartesian         * UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                   GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                   GridDefaultMpi());
  //GridRedBlackCartesian * UrbGrid = SpaceTimeGrid::makeFourDimRedBlackGrid(UGrid);

  std::vector<LatticeGaugeField> U( 3, UGrid );
  for( std::size_t f = 0; f < File.size(); ++f )
  {
    struct stat buf;
    if( stat(File[f].c_str(), &buf) == -1 )
    {
      std::stringstream ss;
      ss << "Configuration " << (f+1) << " " << File[f] << " doesn't exist";
      if( f == 0 )
        throw std::runtime_error( ss.str() );
      std::cout << GridLogError << ss.str() << std::endl;
    }
    else
    {
      std::cout<<GridLogMessage<<"Loading configuration "<<(f+1)<<" "<<File[f]<<std::endl;
      FieldMetaData header;
      NerscIO::readConfiguration(U[f?1:0], header, File[f]);
      if( f )
      {
        std::cout<<GridLogMessage<<"Computing norm"<<std::endl;
        U[2] = U[0] - U[1];
        std::cout << GridLogMessage << "Norm2( " << File[0] << " - " << File[f] << " ) = "
                  << norm2( U[2] ) << std::endl;
      }
    }
  }
}

int main(int argc, char *argv[])
{
  //std::ios_base::sync_with_stdio( false );
  bool bShowHelp{ true };
  int iReturn = EXIT_SUCCESS;
  std::vector<std::string> Args;
  for( std::size_t i = 1; i < argc && argv[i][0] != '-'; ++i )
    Args.emplace_back( argv[i] );
  try
  {
    if( Args.size() >= 2 )
    {
      Grid::Grid_init( &argc, &argv );
      try
      {
        // std::cout << Grid::GridLogMessage << MLUVersionInfoHuman() << std::endl;
        GaugeCompare( Args );
        std::cout << GridLogMessage << "Gauge comparison complete" << std::endl;
        bShowHelp = false;
      } catch( ... ) {
        Grid::Grid_finalize();
        throw;
      }
      Grid::Grid_finalize();
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
  if( bShowHelp )
  {
    std::cout << "Compare multiple Grid gauge fields, report norm2(fn - f1)\n"
                 "Arguments are the gauge-field filenames, all of which are compared to first" << std::endl;
  }
  return iReturn;
}
