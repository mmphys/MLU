/*****************************************************************************
 
 Wrapper for GSL vectors and matrices

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

#ifndef MLU_GSLVecMat_hpp
#define MLU_GSLVecMat_hpp

#include <MLU/MLUFirst.hpp>

#include <cassert>
#include <complex>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <math.h>
#include <vector>

// GSL
#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <gsl/gsl_sf_bessel.h>

BEGIN_MLU_NAMESPACE

extern const double NaN;

using fint = std::uint_fast32_t; // type for random numbers
using SeedType = fint; // TODO: Merge all occurences into a single type

/// Test whether a type is complex. Provide underlying Real type
template<typename T> struct is_complex                     : public std::false_type{ using Real=T; };
template<typename T> struct is_complex<std::complex<T>>    : public std::true_type { using Real=T; };

// Allow real and complex numbers to be squared and square rooted
template <typename T> typename std::enable_if< is_complex<T>::value, T>::type
Squared( T z ) { return z * std::conj( z ); }
template <typename T> typename std::enable_if<!is_complex<T>::value, T>::type
Squared( T z ) { return z * z; }

// Allow real and complex numbers to be conjugated
template <typename T> typename std::enable_if< is_complex<T>::value, T>::type
Conjugate( T z ) { return std::conj( z ); }
template <typename T> typename std::enable_if<!is_complex<T>::value, T>::type
Conjugate( T z ) { return z; }

// Allow real and complex numbers to be squared and square rooted
template <typename T> typename std::enable_if< is_complex<T>::value, T>::type
ComponentSquared( T z ) { return { z.real() * z.real(), z.imag() * z.imag() }; }
template <typename T> typename std::enable_if<!is_complex<T>::value, T>::type
ComponentSquared( T z ) { return z * z; }

// Allow real and complex numbers to be squared and square rooted
template <typename T> typename std::enable_if< is_complex<T>::value, T>::type
ComponentDivide( T Num, T Den ) { return { Num.real() / Den.real(), Num.imag() / Den.imag() }; }
template <typename T> typename std::enable_if<!is_complex<T>::value, T>::type
ComponentDivide( T Num, T Den ) { return Num / Den; }

// Component-wise absolute value for complex types
template<typename T> typename std::enable_if<is_complex<T>::value, T>::type
ComponentAbs( T c ) { return { std::abs( c.real() ), std::abs( c.imag() ) }; }

// Component-wise absolute value for scalars
template<typename T> typename std::enable_if<!(is_complex<T>::value), T>::type
ComponentAbs( T r ) { return std::abs( r ); }

// Component-wise minimum for complex types
template<typename T> typename std::enable_if<is_complex<T>::value, T>::type
ComponentMin( const T &a, const T &b ) { return T( std::min( a.real(), b.real() ),
                                                   std::min( a.imag(), b.imag() ) ); }

// Component-wise minimum for scalars
template<typename T> typename std::enable_if<!(is_complex<T>::value), T>::type
ComponentMin( const T &a, const T &b ) { return std::min( a, b ); }

// Component-wise maximum for complex types
template<typename T> typename std::enable_if<is_complex<T>::value, T>::type
ComponentMax( const T &a, const T &b ) { return T( std::max( a.real(), b.real() ),
                                                   std::max( a.imag(), b.imag() ) ); }

// Component-wise maximum for scalars
template<typename T> typename std::enable_if<!(is_complex<T>::value), T>::type
ComponentMax( const T &a, const T &b ) { return std::max( a, b ); }

// Component-wise relative difference for complex types
template<typename T> typename std::enable_if<is_complex<T>::value, typename is_complex<T>::Real>::type
ComponentRelDif( const T &a, const T &b )
{
  const typename is_complex<T>::Real re{ std::abs( a.real() - b.real() ) / a.real() };
  const typename is_complex<T>::Real im{ std::abs( a.imag() - b.imag() ) / a.imag() };
  return std::max( re, im );
}

// Component-wise relative difference for scalars
template<typename T> typename std::enable_if<!is_complex<T>::value, typename is_complex<T>::Real>::type
ComponentRelDif( const T &a, const T &b )
{
  return std::abs( a - b ) / a;
}

// IsFinite() for floats and complex types
template<typename T> inline typename std::enable_if<is_complex<T>::value, bool>::type
IsFinite( const T &c ) { return std::isfinite( c.real() ) && std::isfinite( c.imag() ); }
template<typename T> inline typename std::enable_if<!is_complex<T>::value, bool>::type
IsFinite( const T &c ) { return std::isfinite( c ); }

// Are all the floating point numbers pointed to finite
template <typename T> inline bool IsFinite( const T * d, std::size_t n )
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

/**
 GSL error handling switched off when program starts
 */
class GSLLibraryGlobal
{
private:
  gsl_error_handler_t * pOldErrorHandler;

public:
  GSLLibraryGlobal();
  ~GSLLibraryGlobal();
  // Use this to throw a std::runtime_error describing the GSL error
  static void Error( const char * reason, const char * file, int line, int gsl_errno );
};

// test whether a type is complex
template<typename T> struct is_gsl_scalar             : public std::false_type {};
template<> struct is_gsl_scalar<float>                : public std::true_type {};
template<> struct is_gsl_scalar<double>               : public std::true_type {};
template<> struct is_gsl_scalar<std::complex<float>>  : public std::true_type {};
template<> struct is_gsl_scalar<std::complex<double>> : public std::true_type {};
template<> struct is_gsl_scalar<fint>                 : public std::true_type {};

// My support for vectors and matrices - wrapper for GSL
template <typename T> struct GSLTraits;
template <typename T> struct Vector;
template <typename T> struct Matrix;

template<> struct GSLTraits<double>
{
  using Scalar     = double;
  using Real       = double;
  using GSLScalar  = double;
  using GSLBlockType  = gsl_block;
  using GSLVectorType = gsl_vector;
  using GSLMatrixType = gsl_matrix;
};

template<> struct GSLTraits<float>
{
  using Scalar     = float;
  using Real       = float;
  using GSLScalar  = float;
  using GSLBlockType  = gsl_block_float;
  using GSLVectorType = gsl_vector_float;
  using GSLMatrixType = gsl_matrix_float;
};

template<> struct GSLTraits<std::complex<double>>
{
  using Scalar     = std::complex<double>;
  using Real       = double;
  using GSLScalar  = gsl_complex;
  using GSLBlockType  = gsl_block_complex;
  using GSLVectorType = gsl_vector_complex;
  using GSLMatrixType = gsl_matrix_complex;
};

template<> struct GSLTraits<std::complex<float>>
{
  using Scalar     = std::complex<float>;
  using Real       = float;
  using GSLScalar  = gsl_complex_float;
  using GSLBlockType  = gsl_block_complex_float;
  using GSLVectorType = gsl_vector_complex_float;
  using GSLMatrixType = gsl_matrix_complex_float;
};

struct MLU_block_fint
{
  size_t size;
  fint *data;
};

struct MLU_vector_fint
{
  size_t size;
  size_t stride;
  fint *data;
  MLU_block_fint *block;
  int owner;
};

struct MLU_matrix_fint
{
  size_t size1;
  size_t size2;
  size_t tda;
  fint * data;
  MLU_block_fint * block;
  int owner;
};

MLU_block_fint *MLU_block_fint_alloc( const std::size_t n );
void MLU_block_fint_free( MLU_block_fint *b );
int MLU_vector_fint_memcpy( MLU_vector_fint *dest, const MLU_vector_fint *src );
int MLU_vector_fint_equal( const MLU_vector_fint *u, const MLU_vector_fint *v );
int MLU_matrix_fint_memcpy( MLU_matrix_fint *dest, const MLU_matrix_fint *src );
int MLU_matrix_fint_equal( const MLU_matrix_fint *u, const MLU_matrix_fint *v );

template<> struct GSLTraits<fint>
{
  using Scalar     = fint;
  using Real       = fint;
  using GSLScalar  = fint;
  using GSLBlockType  = MLU_block_fint;
  using GSLVectorType = MLU_vector_fint;
  using GSLMatrixType = MLU_matrix_fint;
};

//#define EXPAND(...) __VA_ARGS__

#define COMMON_GSL_TYPE double
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _ ## func
#define COMMON_GSL_OPTIONAL
#define COMMON_GSL_DOUBLE
#define MLU_CBLAS( x ) cblas_d ## x
#define MLU_CBLAS_RET_REAL( x ) cblas_d ## x
#define MLU_CBLAS_SCALAR_PTR double *
#include "GSLVecMatImp.hpp"
#define COMMON_GSL_TYPE float
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _float ## _ ## func
#define MLU_CBLAS( x ) cblas_s ## x
#define MLU_CBLAS_RET_REAL( x ) cblas_s ## x
#define MLU_CBLAS_SCALAR_PTR float *
#include "GSLVecMatImp.hpp"
#define COMMON_GSL_TYPE std::complex<double>
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _complex ## _ ## func
#define COMMON_GSL_IS_COMPLEX
#define COMMON_GSL_OPTIONAL
#define MLU_CBLAS( x ) cblas_z ## x
#define MLU_CBLAS_RET_REAL( x ) cblas_dz ## x
#define MLU_CBLAS_SCALAR_PTR void *
#include "GSLVecMatImp.hpp"
#define COMMON_GSL_TYPE std::complex<float>
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _complex_float ## _ ## func
#define COMMON_GSL_IS_COMPLEX
#define MLU_CBLAS( x ) cblas_c ## x
#define MLU_CBLAS_RET_REAL( x ) cblas_sc ## x
#define MLU_CBLAS_SCALAR_PTR void *
#include "GSLVecMatImp.hpp"
#define COMMON_GSL_TYPE fint
#define COMMON_GSL_FUNC( x, func ) MLU_ ## x ## _fint_ ## func
#include "GSLVecMatImp.hpp"

template <typename T> inline bool IsFinite( const Vector<T> &v )
{
  return v.IsFinite();
}

template <typename T> inline bool IsFinite( const Matrix<T> &m, bool bDiagonalsOnly = false )
{
  return m.IsFinite( bDiagonalsOnly );
}

// GSL objects don't work well with const objects

template <typename T> class MatrixView
{
  using Scalar = T;
protected:
  const T * data_;
  std::size_t size1_;
  std::size_t size2_;
  std::size_t tda_;
public:
  const T * data() const { return data_; }
  std::size_t size() const { return size1_ * size2_; }
  std::size_t size1() const { return size1_; }
  std::size_t size2() const { return size2_; }
  std::size_t tda() const { return tda_; }
  /*const T * begin() const { return data_; }
  const T * end() const { return data_ + size_; }*/
  void clear() { data_ = nullptr; size1_ = 0; size2_ = 0; tda_ = 0; }
  inline const Scalar & operator()( std::size_t i, std::size_t j ) const
  {
    if( i >= size1_ || j >= size2_ )
    {
      std::ostringstream ss;
      ss << "MatrixView( " << i << ", " << j << " ) out of bounds (" << size1_ << ", " << size2_ << ")";
      throw std::runtime_error( ss.str().c_str() );
    }
    return data_[i * tda_ + j];
  }
  inline void Map( const T * Data, std::size_t Size1, std::size_t Size2, std::size_t Tda )
  { data_ = Data; size1_ = Size1; size2_ = Size2; tda_ = Tda; }
  inline void Map( const T * Data, std::size_t Size1, std::size_t Size2 ) { Map( Data, Size1, Size2, Size2 ); }
  template <typename U=T> typename std::enable_if<is_gsl_scalar<U>::value, MatrixView<T> &>::type
  inline operator=( const Matrix<T> &m )
  {
    Map( reinterpret_cast<T *>( m.data ), m.size1, m.size2, m.tda );
    return *this;
  }
  inline MatrixView<T> &operator=( const MatrixView<T> &m )
  {
    Map( m.data(), m.size1(), m.size2(), m.tda() );
    return *this;
  }
  template <typename U=T> typename std::enable_if<is_gsl_scalar<U>::value>::type
  MapRow( Vector<T> &v, std::size_t Row );
  template <typename U=T> typename std::enable_if<is_gsl_scalar<U>::value>::type
  MapColumn( Vector<T> &v, std::size_t Column );
  // Constructors
  MatrixView() { clear(); }
  MatrixView( const Matrix<T> &m ) { *this = m; }
  MatrixView( const MatrixView<T> &m ) { *this = m; }
};

template <typename T> class VectorView
{
  using Scalar = T;
protected:
  const T * data_;
  std::size_t size_;
  std::size_t stride_;
public:
  const T * data() const { return data_; }
  std::size_t size() const { return size_; }
  void size( std::size_t NewSize ) { size_ = NewSize; }
  std::size_t stride() const { return stride_; }
  const T * begin() const { return data_; }
  const T * end() const { return data_ + size_ * stride_; }
  void clear() { data_ = nullptr; size_ = 0; stride_ = 1; }
  inline const Scalar & operator[]( std::size_t i ) const
  {
    if( i >= size_ )
      throw std::runtime_error("VectorView index " + std::to_string(i) + " >= size " + std::to_string(size_));
    return data_[i * stride_];
  }
  inline void Map( const T * Data, std::size_t Size, std::size_t Stride = 1 )
  { data_ = Data; size_ = Size; stride_ = Stride; }
  template <typename U=T> typename std::enable_if<is_gsl_scalar<U>::value, VectorView<T> &>::type
  inline operator=( const Vector<T> &v )
  {
    Map( reinterpret_cast<T *>( v.data ), v.size, v.stride );
    return *this;
  }
  inline VectorView<T> &operator=( const VectorView<T> &v )
  {
    Map( v.data(), v.size(), v.stride() );
    return *this;
  }
  inline VectorView<T> &operator=( const std::vector<T> &v )
  {
    Map( v.data(), v.size() );
    return *this;
  }
  inline VectorView<T> operator+=( std::size_t Delta ) { data_ += Delta * stride_; return *this; };
  template <typename U=T> typename std::enable_if<is_gsl_scalar<U>::value>::type
    MapRow( const Matrix<T> &m, std::size_t Row );
  void MapRow( const MatrixView<T> &m, std::size_t Row );
  template <typename U=T> typename std::enable_if<is_gsl_scalar<U>::value>::type
    MapColumn( const Matrix<T> &m, std::size_t Column );
  void MapColumn( const MatrixView<T> &m, std::size_t Column );
  // Constructors
  VectorView() { clear(); }
  VectorView( const T * Data, std::size_t Size, std::size_t Stride = 1 ) { Map( Data, Size, Stride ); }
  VectorView( const Vector<T> &m ) { *this = m; }
  VectorView( const VectorView<T> &m ) { *this = m; }
  VectorView( const std::vector<T> &m ) { *this = m; }
};

inline double BesselK1( double x, double * pError = nullptr )
{
  gsl_sf_result gsfRes;
  int gsl_e{ gsl_sf_bessel_K1_e( x, &gsfRes ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "gsl_sf_bessel_K1_e()", __FILE__, __LINE__, gsl_e );
  if( pError )
    *pError = gsfRes.err;
  return gsfRes.val;
}

END_MLU_NAMESPACE
#endif
