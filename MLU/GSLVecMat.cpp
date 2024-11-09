/*****************************************************************************
 
 Wrapper for GSL vectors and matrices

 Source file: GSLVecMat.cpp
 
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

#include <MLUconfig.h>
#include "GSLVecMat.hpp"

BEGIN_MLU_NAMESPACE

/*****************************************************************************
 
 I will use the GNU Scientific Library, but I will have GSL errors throw exceptions

*****************************************************************************/

GSLLibraryGlobal gslLibraryGlobal;

GSLLibraryGlobal::GSLLibraryGlobal()
{
  pOldErrorHandler = gsl_set_error_handler_off();// gsl_set_error_handler( Error );
}

GSLLibraryGlobal::~GSLLibraryGlobal()
{
  gsl_set_error_handler( pOldErrorHandler );
}

void GSLLibraryGlobal::Error( const char * reason, const char * file, int line, int gsl_errno )
{
  if( gsl_errno )
  {
    std::stringstream ss;
    ss << "gsl: " << file << ":" << line << ": errno " << gsl_errno << ": " << gsl_strerror( gsl_errno );
    if( reason )
      ss << ": " << reason;
    throw std::runtime_error( ss.str().c_str() );
  }
}

/*****************************************************************************
 
Support my (cursory) implementation of vectors and matrices of fint

*****************************************************************************/

static constexpr std::size_t MLU_fint_size{ sizeof( fint ) };
static constexpr std::size_t MLU_block_fint_size{ ( sizeof( MLU_block_fint ) + MLU_fint_size - 1 )
                                                  & ~( MLU_fint_size - 1 ) };

MLU_block_fint *MLU_block_fint_alloc( const std::size_t n )
{
  void * pv = malloc( MLU_block_fint_size + MLU_fint_size * n );
  MLU_block_fint * p = static_cast<MLU_block_fint *>( pv );
  if( p )
  {
    p->size = n;
    p->data = reinterpret_cast<fint *>( static_cast<char *>( pv ) + MLU_block_fint_size );
  }
  return p;
}

void MLU_block_fint_free( MLU_block_fint *b )
{
  if( b )
    free( b );
}

int MLU_vector_fint_memcpy( MLU_vector_fint *dest, const MLU_vector_fint *src )
{
  if( !dest || !src )
    throw std::runtime_error( "MLU_vector_fint_memcpy() null src or dest" );
  if( dest->size != src->size )
    throw std::runtime_error( "MLU_vector_fint_memcpy() size mismatch" );
        Vector<fint> &d{ *reinterpret_cast<      Vector<fint> *>( dest ) };
  const Vector<fint> &s{ *reinterpret_cast<const Vector<fint> *>( src ) };
  for( std::size_t i = 0; i < dest->size; ++i )
    d[i] = s[i];
  return 0;
}

int MLU_vector_fint_equal( const MLU_vector_fint *u, const MLU_vector_fint *v )
{
  if( !u || !v )
    throw std::runtime_error( "MLU_vector_fint_equal() null u or v" );
  if( u->size != v->size )
    return 0;
  const Vector<fint> &l{ *reinterpret_cast<const Vector<fint> *>( u ) };
  const Vector<fint> &r{ *reinterpret_cast<const Vector<fint> *>( v ) };
  for( std::size_t i = 0; i < u->size; ++i )
    if( l[i] != r[i] )
      return 0;
  return 1;
}

int MLU_matrix_fint_memcpy( MLU_matrix_fint *dest, const MLU_matrix_fint *src )
{
  if( !dest || !src )
    throw std::runtime_error( "MLU_matrix_fint_memcpy() null src or dest" );
  if( dest->size1 != src->size1 || dest->size2 != src->size2 )
    throw std::runtime_error( "MLU_matrix_fint_memcpy() size mismatch" );
        Matrix<fint> &d{ *reinterpret_cast<      Matrix<fint> *>( dest ) };
  const Matrix<fint> &s{ *reinterpret_cast<const Matrix<fint> *>( src ) };
  for( std::size_t i = 0; i < dest->size1; ++i )
    for( std::size_t j = 0; j < dest->size2; ++j )
      d(i,j) = s(i,j);
  return 0;
}

int MLU_matrix_fint_equal( const MLU_matrix_fint *u, const MLU_matrix_fint *v )
{
  if( !u || !v )
    throw std::runtime_error( "MLU_matrix_fint_equal() null u or v" );
  if( u->size1 != v->size1 || u->size2 != v->size2 )
    return 0;
  const Matrix<fint> &l{ *reinterpret_cast<const Matrix<fint> *>( u ) };
  const Matrix<fint> &r{ *reinterpret_cast<const Matrix<fint> *>( v ) };
  for( std::size_t i = 0; i < u->size1; ++i )
    for( std::size_t j = 0; j < u->size2; ++j )
      if( l(i,j) != r(i,j) )
        return 0;
  return 1;
}

/*****************************************************************************
 
 Views of GSL Matrix and Vector
 ... because GSL doesn't play nicely with const objects

*****************************************************************************/

template <typename T>
template <typename U> typename std::enable_if<is_gsl_scalar<U>::value>::type
MatrixView<T>::MapRow( Vector<T> &v, std::size_t Row )
{
  if( size1_ == 0 || size2_ == 0 )
    throw std::runtime_error( "MatrixView::MapRow() empty matrix" );
  if( Row >= size1_ )
    throw std::runtime_error( "Row " + std::to_string( Row ) + " >= " + std::to_string( size1_ ) );
  v.MapView( data_ + tda_ * Row, size2_, 1 );
}

template <typename T>
template <typename U> typename std::enable_if<is_gsl_scalar<U>::value>::type
MatrixView<T>::MapColumn( Vector<T> &v, std::size_t Column )
{
  if( size1_ == 0 || size2_ == 0 )
    throw std::runtime_error( "MatrixView::MapColumn() empty matrix" );
  if( Column >= size2_ )
    throw std::runtime_error( "Column " + std::to_string(Column) + " >= " + std::to_string(size2_) );
  MapView( data_ + Column, size1_, tda_ );
}

template class MatrixView<fint>;
template class MatrixView<double>;
template class MatrixView<float>;
template class MatrixView<std::complex<double>>;
template class MatrixView<std::complex<float>>;

template <typename T>
template <typename U> typename std::enable_if<is_gsl_scalar<U>::value>::type
VectorView<T>::MapRow( const Matrix<T> &m, std::size_t Row )
{
  if( m.size1 == 0 || m.size2 == 0 )
    throw std::runtime_error( "VectorView::MapRow() empty matrix" );
  if( Row >= m.size1 )
    throw std::runtime_error( "Row " + std::to_string( Row ) + " > " + std::to_string( m.size1 ) );
  Map( reinterpret_cast<Scalar *>( m.data ) + m.tda * Row, m.size2, 1 );
}

template <typename T>
void VectorView<T>::MapRow( const MatrixView<T> &m, std::size_t Row )
{
  if( m.size1() == 0 || m.size2() == 0 )
    throw std::runtime_error( "VectorView::MapRow() empty matrix" );
  if( Row >= m.size1() )
    throw std::runtime_error( "Row " + std::to_string( Row ) + " > " + std::to_string( m.size1() ) );
  Map( m.data() + m.tda() * Row, m.size2(), 1 );
}

template <typename T>
template <typename U> typename std::enable_if<is_gsl_scalar<U>::value>::type
VectorView<T>::MapColumn( const Matrix<T> &m, std::size_t Column )
{
  if( m.size1 == 0 || m.size2 == 0 )
    throw std::runtime_error( "VectorView::MapColumn() empty matrix" );
  if( Column >= m.size2 )
    throw std::runtime_error( "Column " + std::to_string( Column ) + " > " + std::to_string( m.size2 ) );
  Map( reinterpret_cast<Scalar *>( m.data ) + Column, m.size1, m.tda );
}

template <typename T>
void VectorView<T>::MapColumn( const MatrixView<T> &m, std::size_t Column )
{
  if( m.size1() == 0 || m.size2() == 0 )
    throw std::runtime_error( "VectorView::MapColumn() empty matrix" );
  if( Column >= m.size2() )
    throw std::runtime_error( "Column " + std::to_string( Column ) + " > " + std::to_string( m.size2() ) );
  Map( m.data() + Column, m.size1(), m.tda() );
}

template class VectorView<fint>;
template class VectorView<double>;
template class VectorView<float>;
template class VectorView<std::complex<double>>;
template class VectorView<std::complex<float>>;

template void VectorView<float>::MapRow( const Matrix<float> &, std::size_t );
template void VectorView<double>::MapRow( const Matrix<double> &, std::size_t Row );
template void VectorView<std::complex<float>>::MapRow( const Matrix<std::complex<float>> &, std::size_t );
template void VectorView<std::complex<double>>::MapRow( const Matrix<std::complex<double>> &, std::size_t );
template void VectorView<float>::MapColumn( const Matrix<float> &, std::size_t );
template void VectorView<double>::MapColumn( const Matrix<double> &, std::size_t Row );
template void VectorView<std::complex<float>>::MapColumn( const Matrix<std::complex<float>> &, std::size_t );
template void VectorView<std::complex<double>>::MapColumn( const Matrix<std::complex<double>> &, std::size_t );

END_MLU_NAMESPACE
