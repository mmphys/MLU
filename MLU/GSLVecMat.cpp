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

#include "GSLVecMat.hpp"

MLU_GSLVecMat_hpp

/*****************************************************************************
 
 I will use the GNU Scientific Library, but I will have GSL errors throw exceptions

*****************************************************************************/

class GSLLibraryGlobal
{
private:
  gsl_error_handler_t * pOldErrorHandler;
  static void GSLErrorHandler(const char * reason, const char * file, int line, int gsl_errno );

public:
  GSLLibraryGlobal();
  ~GSLLibraryGlobal();
};

GSLLibraryGlobal gslLibraryGlobal;

GSLLibraryGlobal::GSLLibraryGlobal()
{
  pOldErrorHandler = gsl_set_error_handler( GSLErrorHandler );
}

GSLLibraryGlobal::~GSLLibraryGlobal()
{
  gsl_set_error_handler( pOldErrorHandler );
}

void GSLLibraryGlobal::GSLErrorHandler(const char * reason, const char * file, int line, int gsl_errno )
{
  std::stringstream ss;
  ss << "gsl: " << file << ":" << line << ": errno " << gsl_errno << ": " << reason;
  throw std::runtime_error( ss.str().c_str() );
}

/*****************************************************************************
 
 Views of GSL Matrix and Vector
 ... because GSL doesn't play nicely with const objects

*****************************************************************************/

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

MLU_GSLVecMat_hpp_end
