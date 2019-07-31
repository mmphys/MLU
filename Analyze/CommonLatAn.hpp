/*************************************************************************************
 
 Common utilities (with dependencies on c++ stdlib and LatAnalyze)
 
 Source file: CommonLatAn.hpp
 
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

#ifndef CommonLatAn_hpp
#define CommonLatAn_hpp

#include "Common.hpp"
#include <LatAnalyze/Statistics/Dataset.hpp>

BEGIN_COMMON_NAMESPACE

// Are all the floating point numbers in this matrix finite
template <typename T> inline bool IsFinite( const Latan::Mat<T> & m, bool bDiagonalsOnly = false )
{
  for( Latan::Index row = 0; row < m.rows(); ++row )
    for( Latan::Index col = 0; col < m.cols(); ++col )
      if( ( !bDiagonalsOnly || row == col ) && !std::isfinite( m( row, col ) ) )
        return false;
  return true;
}

// Read a list of bootstrapped correlators into a single correlator
Latan::DMatSample ReadBootstrapCorrs( const std::vector<std::string> & FileName, int Fold, int Shift, int NumOps );

END_COMMON_NAMESPACE
#endif // CommonLatAn_hpp
