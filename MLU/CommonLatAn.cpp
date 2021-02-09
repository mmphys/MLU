/*************************************************************************************
 
 Common utilities (with dependencies on c++ stdlib and LatAnalyze)
 
 Source file: CommonLatAn.cpp
 
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

#include "CommonLatAn.hpp"

BEGIN_COMMON_NAMESPACE

// Make summary files of a bootstrap of a correlator
/*void SummariseBootstrap(const Latan::DMatSample &out, const std::string & sOutFileBase,
     const std::string & sType, Latan::SeedType Seed, const std::vector<std::string> ParamNames )
{
  assert( std::isnan( Common::NaN ) && "Compiler does not support quiet NaNs" );
  const int numParams{ static_cast<int>( ParamNames.size() ) };
  assert( numParams == out[Latan::central].rows() && "Parameter names missing" );
  const std::size_t   nSample{ static_cast<std::size_t>( out.size() )};
  std::vector<double> Data( nSample );
  std::string sOutFileName{ MakeFilename(sOutFileBase, SummaryNames[f], Seed, TEXT_EXT) };
  std::ofstream s( sOutFileName );
  s << Comment << sType << NewLine << Comment << sOutFileBase << "\n# Seed " << Seed
  << NewLine << FieldNames << ( ( f == 0 ) ? FieldNames2 : "" );
  s << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
  for(int t = 0; t < nt; t++)
  {
  }
}*/

END_COMMON_NAMESPACE
