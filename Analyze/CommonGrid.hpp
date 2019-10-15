/*************************************************************************************
 
 Common utilities (with dependencies on Grid)
 
 Source file: CommonGrid.hpp
 
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

#ifndef CommonGrid_hpp
#define CommonGrid_hpp

#include "Common.hpp"
#include <Grid/GridCore.h>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

BEGIN_COMMON_NAMESPACE

// Load all combinations of gamma_snk and gamma_src from an HDF5 file into the array of correlators
// If the correlators are empty, resize to fit
// Otherwise, if the sizes don't match, throw an error
void ReadComplexArray(std::vector<std::vector<std::complex<double>>> &buffer, const std::vector<Grid::Gamma::Algebra> &alg, const std::string &FileName, unsigned int tOffset);

END_COMMON_NAMESPACE
#endif // CommonGrid_hpp
