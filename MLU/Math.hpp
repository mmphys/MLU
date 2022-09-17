/*****************************************************************************
 
 Basic math support

 Source file: GSLVecMat.hpp
 
 Copyright (C) 2019-2022
 
 Author: Michael Marshall <Mike@lqcd.me>

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
 ****************************************************************************/
/*  END LEGAL */

#ifndef MLU_Math_hpp
#define MLU_Math_hpp namespace Common {
#define MLU_Math_hpp_end };

#include <complex>
#include <type_traits>

// Eigen dense matrices
#include <Grid/Eigen/Dense>

MLU_Math_hpp

using fint = std::uint_fast32_t; // type for random numbers

extern const double NaN;

// test whether a type is complex
template<typename T> struct is_complex                  : public std::false_type {};
template<typename T> struct is_complex<std::complex<T>> : public std::true_type {};

// Allow real and complex numbers to be squared and square rooted
template <typename T> typename std::enable_if< is_complex<T>::value, T>::type
Squared( T z ) { return z * std::conj( z ); }
template <typename T> typename std::enable_if<!is_complex<T>::value, T>::type
Squared( T z ) { return z * z; }

// Component-wise absolute value for complex types
template<typename T> typename std::enable_if<is_complex<T>::value, T>::type
ComponentAbs( T c ) { return { std::abs( c.real() ), std::abs( c.imag() ) }; }

// Component-wise absolute value for scalars
template<typename T> typename std::enable_if<!(is_complex<T>::value), T>::type
ComponentAbs( T r ) { return std::abs( r ); }

// IsFinite() for floats and complex types
template<typename T> inline typename std::enable_if<is_complex<T>::value, bool>::type
IsFinite( const T &c ) { return std::isfinite( c.real() ) && std::isfinite( c.imag() ); }
template<typename T> inline typename std::enable_if<std::is_floating_point<T>::value, bool>::type
IsFinite( const T &c ) { return std::isfinite( c ); }

// Are all the floating point numbers pointed to finite
template <typename T, typename I> inline bool IsFinite( const T * d, I n )
{
  while( n-- )
    if( !IsFinite( *d++ ) )
      return false;
  return true;
}

// Are all the floating point numbers in this vector finite
template <typename T> inline bool IsFinite( const std::vector<T> & v )
{
  for( const T &n : v )
    if( !IsFinite( n ) )
      return false;
  return true;
}

// Are all the floating point numbers in this Eigen::matrix finite
template <typename T> inline bool IsFinite( const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & m, bool bDiagonalsOnly = false )
{
  for( Eigen::Index row = 0; row < m.rows(); ++row )
    for( Eigen::Index col = 0; col < m.cols(); ++col )
      if( ( !bDiagonalsOnly || row == col ) && !IsFinite( m( row, col ) ) )
        return false;
  return true;
}

MLU_Math_hpp_end
#endif
