/*************************************************************************************
 
 MLU utilities (with dependencies on Grid)
 
 Source file: MLUGrid.hpp
 
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

#ifndef MLUGrid_hpp
#define MLUGrid_hpp

#include <MLU/MLU.hpp>
#include <Grid/GridCore.h>
#include <Grid/Eigen/Dense>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

BEGIN_MLU_NAMESPACE

// Load all combinations of gamma_snk and gamma_src from an HDF5 file into the array of correlators
// If the correlators are empty, resize to fit
// Otherwise, if the sizes don't match, throw an error
void ReadComplexArray(std::vector<std::vector<std::complex<double>>> &buffer, const std::vector<Grid::Gamma::Algebra> &alg, const std::string &FileName, unsigned int tOffset);

// Are all the floating point numbers in this Eigen::matrix finite
template <typename T> inline bool IsFinite( const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & m, bool bDiagonalsOnly = false )
{
  for( Eigen::Index row = 0; row < m.rows(); ++row )
    for( Eigen::Index col = 0; col < m.cols(); ++col )
      if( ( !bDiagonalsOnly || row == col ) && !IsFinite( m( row, col ) ) )
        return false;
  return true;
}

END_MLU_NAMESPACE
#endif // MLUGrid_hpp
