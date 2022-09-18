/*****************************************************************************
 
 Basic math support

 Source file: Math.hpp
 
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
#include <vector>

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

// A q-value is an integral of the distribution from x to infinity
// i.e. the probability of obtaining a worse test statistic that is still explained by the model
// (I think) this is what people normally discuss when they say "p-value"
// Don't rely on this - the correct distribution (when covariance comes from the data) is Hotelling )below)
template <typename T> T qValueChiSq( T ChiSquared, unsigned int dof );

// Hotelling's T^2 test statistic is computed exactly as one would for Chi Squared
// However, since the covariance matrix comes from the data (rather than being known a priori)
// We find that the test statistic (t^2) has an T^2-distribution with p and m dof
//    p = Covariance matrix (i.e. fit) degrees of freedom: num data points - num model parameters
//    m = Statistically independent samples used to build ovariance matrix: usually number of configs - 1
// See my year 4 notes and
// https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/ckelly/Gparity/hotelling_v10.pdf
// https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution

struct HotellingDist
{
  const unsigned int p; // Dimension of covariance matrix (degrees of freedom of fit)
  const unsigned int m; // Degrees of freedom of covariance matrix (num samples to build covar - 1)
  const unsigned int Nu;
  const double Factor;
  HotellingDist( unsigned int p, unsigned int m );
  double operator()( double TestStatistic ) const;
  inline static bool Usable( unsigned int p, unsigned int m ) { return p <= m; }
  template <typename T> T static qValue( T TestStatistic, unsigned int p, unsigned int m );
};

// This is prior version. Not used except for reading in old files
template <typename T> struct ValWithErOldV1
{
  T Central;
  T Low;
  T High;
  T Check;
};

template <typename T = double>
struct ValWithEr
{
  T Min;
  T Low;
  T Central;
  T High;
  T Max;
  T Check;
  ValWithEr() = default;
  template <typename U=T>
  ValWithEr( typename std::enable_if<!is_complex<U>::value, T>::type dCentral,
             std::vector<T> &Data, std::size_t Count )
  { Get( dCentral, Data, Count ); }
  ValWithEr( T Min, T Low, T Central, T High, T Max, T Check = 1 );
  static void Header( const std::string &FieldName, std::ostream &os, const std::string &Sep = " " );
  ValWithEr<T>& operator = ( const T Scalar );
  ValWithEr<T>& operator = ( const ValWithEr<T> &Other );
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>&>::type
  operator *=( const T Scalar );
  template <typename U=T> typename std::enable_if< is_complex<U>::value, ValWithEr<T>&>::type
  operator *=( const T Scalar );
  ValWithEr<T>  operator * ( const T Scalar ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>>::type
  qValueChiSq( unsigned int dof ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>>::type
  qValueHotelling( unsigned int p, unsigned int m ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value>::type
  Get( T dCentral, std::vector<T> &Data, std::size_t Count );
};

template <typename T>
std::ostream & operator<<( std::ostream &os, const ValWithEr<T> &v );

MLU_Math_hpp_end
#endif
