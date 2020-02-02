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

// Copy a correlator to a Latan::DMat
/*inline void CopyCorrelator( Latan::DMat &dst, const Correlator &src, int iOffset = 0, bool bSwapRealImag = false )
{
  const std::size_t Nt{ src.size() };
  if( Nt == 0 )
    throw std::runtime_error( "Can't copy an uninitialised Correlator" );
  if( Nt > INT_MAX )
    throw std::runtime_error( "Correlator size > INT_MAX" );
  dst.resize( Nt, 2 ); // LatAnalyze prefers nt * 2 real matrix
  const std::size_t dt{ ( iOffset < 0 ) ? Nt - ( -iOffset % Nt ) : iOffset % Nt };
  for( std::size_t t = 0; t < Nt; t++ )
  {
    const std::complex<double> & z{ src[( t + dt ) % Nt] };
    dst( t, 0 ) = bSwapRealImag ? z.imag() : z.real();
    dst( t, 1 ) = bSwapRealImag ? z.real() : z.imag();
  }
}*/

// Copy a correlator to a Latan::DMat
inline void CopyRealCorrelator( Latan::DMat &dst, const std::vector<double> &src, int iOffset = 0, bool bSwapRealImag = false )
{
  const std::size_t Nt{ src.size() };
  if( Nt == 0 )
    throw std::runtime_error( "Can't copy an uninitialised Correlator" );
  if( Nt > INT_MAX )
    throw std::runtime_error( "Correlator size > INT_MAX" );
  dst.resize( Nt, 2 ); // LatAnalyze prefers nt * 2 real matrix
  const std::size_t dt{ ( iOffset < 0 ) ? Nt - ( -iOffset % Nt ) : iOffset % Nt };
  for( std::size_t t = 0; t < Nt; t++ )
  {
    const std::complex<double> & z{ src[( t + dt ) % Nt], 0 };
    dst( t, 0 ) = bSwapRealImag ? z.imag() : z.real();
    dst( t, 1 ) = bSwapRealImag ? z.real() : z.imag();
  }
}

// Make a copy of str, in which token has been replaced by x
template <typename T> inline std::string tokenReplaceCopy(const std::string &str, const std::string token, const T &x )
{
  std::string sCopy{str};
  Latan::tokenReplace(sCopy, token, x );
  return sCopy;
}

END_COMMON_NAMESPACE
#endif // CommonLatAn_hpp
