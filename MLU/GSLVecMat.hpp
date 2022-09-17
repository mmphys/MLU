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
#define MLU_GSLVecMat_hpp namespace Common {
#define MLU_GSLVecMat_hpp_end };

#include <MLU/Math.hpp>

#include <vector>

// GSL
#define HAVE_INLINE
#define GSL_RANGE_CHECK_OFF
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

MLU_GSLVecMat_hpp

// GSL error handling switched off when program starts
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

//#define EXPAND(...) __VA_ARGS__

#define COMMON_GSL_TYPE double
#define COMMON_GSL_BLAS( x ) gsl_blas_d ## x
#define COMMON_GSL_BLAS_REAL( x ) gsl_blas_d ## x
#define COMMON_GSL_BLAS_CPLX( x ) gsl_blas_d ## x
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _ ## func
#define COMMON_GSL_OPTIONAL
#define COMMON_GSL_DOUBLE
#include "GSLVecMatImp.hpp"
#undef COMMON_GSL_DOUBLE
#undef COMMON_GSL_OPTIONAL
#define COMMON_GSL_TYPE float
#define COMMON_GSL_BLAS( x ) gsl_blas_s ## x
#define COMMON_GSL_BLAS_REAL( x ) gsl_blas_s ## x
#define COMMON_GSL_BLAS_CPLX( x ) gsl_blas_s ## x
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _float ## _ ## func
#include "GSLVecMatImp.hpp"
#define COMMON_GSL_TYPE std::complex<double>
#define COMMON_GSL_BLAS( x ) gsl_blas_z ## x
#define COMMON_GSL_BLAS_REAL( x ) gsl_blas_dz ## x
#define COMMON_GSL_BLAS_CPLX( x ) gsl_blas_z ## x ## u
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _complex ## _ ## func
#define COMMON_GSL_OPTIONAL
#include "GSLVecMatImp.hpp"
#undef COMMON_GSL_OPTIONAL
#define COMMON_GSL_TYPE std::complex<float>
#define COMMON_GSL_BLAS( x ) gsl_blas_c ## x
#define COMMON_GSL_BLAS_REAL( x ) gsl_blas_sc ## x
#define COMMON_GSL_BLAS_CPLX( x ) gsl_blas_c ## x ## c
#define COMMON_GSL_FUNC( x, func ) gsl_ ## x ## _complex_float ## _ ## func
#include "GSLVecMatImp.hpp"
#undef COMMON_GSL_TYPE
#undef COMMON_GSL_BLAS
#undef COMMON_GSL_BLAS_REAL
#undef COMMON_GSL_BLAS_CPLX
#undef COMMON_GSL_FUNC

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

MLU_GSLVecMat_hpp_end
#endif
