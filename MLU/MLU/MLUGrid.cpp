/*************************************************************************************
 
 MLU utilities (with dependencies on Grid)
 
 Source file: MLUGrid.cpp
 
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

// NB This is not part of MLU (which has no dependency on Grid).
// Old code that I can't quite bring myself to delete

#include <MLUconfig.h>
#include "MLUGrid.hpp"

BEGIN_MLU_NAMESPACE

// Load all combinations of gamma_snk and gamma_src from an HDF5 file into the array of correlators
// If the correlators are empty, resize to fit
// Otherwise, if the sizes don't match, throw an error
void ReadComplexArray(std::vector<std::vector<std::complex<double>>> &buffer, const std::vector<Grid::Gamma::Algebra> &alg, const std::string &FileName, unsigned int tOffset)
{
  const unsigned int nGamma{ static_cast<unsigned int>( alg.size() )};
  assert( nGamma > 0 && "Expected one or more gammas" );
  const unsigned int nGammaSq{ nGamma * nGamma };
  assert( buffer.size() == nGammaSq && "Bad number of correlators" );
  unsigned int Nt{ static_cast<unsigned int>( buffer[0].size() ) };
  using Meson = Grid::Hadrons::MContraction::Meson::Result;
  std::cout << "Reading " << FileName << std::endl;
  Grid::Hdf5Reader r( FileName );
  std::vector<Meson> m;
  Grid::read(r, "meson", m);
  std::cout << "Read " << m.size() << " mesons" << std::endl;
  std::vector<int> iCount( nGammaSq, 0 ); // Keep track of how many read
  for( int i = 0; i < m.size(); ++i ) {
    unsigned int iSrc = 0;
    while( iSrc < nGamma && alg[iSrc] != m[i].gamma_src )
      ++iSrc;
    if( iSrc < nGamma ) {
      unsigned int iSnk = 0;
      while( iSnk < nGamma && alg[iSnk] != m[i].gamma_snk )
        ++iSnk;
      if( iSnk < nGamma ) {
        const unsigned int iIndex{ iSnk * nGamma + iSrc };
        if( ++iCount[iIndex] != 1 )
          throw std::runtime_error( "File repeats gamma indices" );
        //std::cout << "meson[" << i << "].gamma_snk=" << m[i].gamma_snk << std::endl;
        //std::cout << "meson[" << i << "].gamma_src=" << m[i].gamma_src << std::endl;
        const unsigned int cSize{ static_cast<unsigned int>( m[i].corr.size() ) };
        if( Nt == 0 ) {
          if( cSize == 0 )
            throw std::runtime_error( "Empty correlator" );
          Nt = cSize;
          for( unsigned int j = 0; j < nGammaSq; j++ )
            buffer[j].resize( Nt );
        } else if( Nt != cSize )
          throw std::runtime_error( "Error: Nt = " + std::to_string( Nt ) + ", but just read correlator of length " + std::to_string( cSize ) );
        for( int t = 0; t < Nt; ++t ) {
          std::complex<double> &z = m[i].corr[ ( t + tOffset ) % Nt ];
          if( !std::isfinite( z.real() ) || !std::isfinite( z.imag() ) )
            throw std::runtime_error( "Correlator contains NaNs" );
          buffer[iIndex][t] = z;
        }
      }
    }
  }
  for( unsigned int i = 0; i < nGammaSq; ++i )
    if( iCount[i] != 1 )
      throw std::runtime_error( "Specified gamma structures not present in HDF5 file" );
}

END_MLU_NAMESPACE
