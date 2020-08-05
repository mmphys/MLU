/**

 Mike's lattice QCD utilities
 
 Source file: CommonGSL.hpp
 
 Copyright (C) 2020
 
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
**/

// Wrapper for GSL

#ifdef COMMON_GSL_TYPE

template <> struct Vector<COMMON_GSL_TYPE> : public GSLTraits<COMMON_GSL_TYPE>::GSLVectorType
{
  using Traits    = GSLTraits<COMMON_GSL_TYPE>;
  using Scalar    = Traits::Scalar;
  using Real      = Traits::Real;
  using GSLScalar = Traits::GSLScalar;
  using GSLBlock  = Traits::GSLBlockType;
  using GSLVector = Traits::GSLVectorType;
  using GSLMatrix = Traits::GSLMatrixType;
  using MyVector  = Vector<COMMON_GSL_TYPE>;
  using MyMatrix  = Matrix<COMMON_GSL_TYPE>;
  inline Vector( std::size_t size );
  inline Vector() : Vector( 0 ) {};
  inline Vector( const MyVector &o ) : Vector() { *this = o; };
  inline Vector( MyVector &&o ) : Vector() { *this = o; };
  inline MyVector & operator=( const Vector &o );
  inline MyVector & operator=( Vector &&o );
  inline ~Vector();
  inline void resize( std::size_t size );
  inline void clear();
  inline void MapView( Scalar * data_, std::size_t size_, std::size_t stride_ = 1 );
  inline void MapView( std::vector<Scalar> &v ) { MapView( v.data(), v.size() ); }
  inline const Scalar & operator[]( std::size_t i ) const;
  inline Scalar & operator[]( std::size_t i );
  inline bool IsFinite() const;
  inline Real norm2() const { return COMMON_GSL_BLAS_REAL( nrm2 )( this ); }
  inline Real norm() const { return norm2(); }
  inline void blas_trmv( CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const GSLMatrix &A );
};

Vector<COMMON_GSL_TYPE>::Vector( std::size_t size_ )
{
  size = 0;
  block = nullptr;
  data = nullptr;
  owner = false;
  resize( size_ );
}

Vector<COMMON_GSL_TYPE>::~Vector()
{
  clear();
};

void Vector<COMMON_GSL_TYPE>::resize( std::size_t size_ )
{
  if( size_ != size )
    clear();
  if( size_ )
  {
    block = COMMON_GSL_FUNC( block, alloc )( size_ );
    data = block->data;
    size = size_;
    owner = true;
    stride = 1;
  }
}

void Vector<COMMON_GSL_TYPE>::clear()
{
  if( owner )
  {
    assert( block && "How did I get to be the owner of a null block?" );
    COMMON_GSL_FUNC( block, free )( block );
    owner = false;
  }
  block = nullptr;
  data = nullptr;
  size = 0;
}

void Vector<COMMON_GSL_TYPE>::MapView( Scalar * data_, std::size_t size_, std::size_t stride_ )
{
  clear();
  size = size_;
  data = reinterpret_cast<Real *>( data_ );
  stride = stride_;
}

Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator=( const Vector &o )
{
  resize( o.size );
  COMMON_GSL_FUNC( vector, memcpy )( this, &o );
  return *this;
}

Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator=( Vector &&o )
{
  size = o.size;
  stride = o.stride;
  data = o.data;
  block = o.block;
  owner = o.owner;
  o.owner = false;
  o.clear();
  return *this;
}

const Vector<COMMON_GSL_TYPE>::Scalar & Vector<COMMON_GSL_TYPE>::operator[]( std::size_t i ) const
{
  assert( data && i < size && "Can't access empty vector" );
  return reinterpret_cast<const Scalar *>( data )[i*stride];
}

Vector<COMMON_GSL_TYPE>::Scalar & Vector<COMMON_GSL_TYPE>::operator[]( std::size_t i )
{
  assert( data && i < size && "Can't access empty vector" );
  return reinterpret_cast<Scalar *>( data )[i*stride];
}

bool Vector<COMMON_GSL_TYPE>::IsFinite() const
{
  for( std::size_t i = 0; i < size; ++i )
    if( !::Common::IsFinite( (*this)[i] ) )
      return false;
  return true;
};

inline void Vector<COMMON_GSL_TYPE>::blas_trmv( CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                                                const Vector<COMMON_GSL_TYPE>::GSLMatrix &A )
{
  COMMON_GSL_BLAS( trmv )( Uplo, TransA, Diag, &A, reinterpret_cast<GSLVector *>( this ) );
}

template <> struct Matrix<COMMON_GSL_TYPE> : public GSLTraits<COMMON_GSL_TYPE>::GSLMatrixType
{
  using Traits    = GSLTraits<COMMON_GSL_TYPE>;
  using Scalar    = Traits::Scalar;
  using Real      = Traits::Real;
  using GSLScalar = Traits::GSLScalar;
  using GSLBlock  = Traits::GSLBlockType;
  using GSLVector = Traits::GSLVectorType;
  using GSLMatrix = Traits::GSLMatrixType;
  using MyVector  = Vector<COMMON_GSL_TYPE>;
  using MyMatrix  = Matrix<COMMON_GSL_TYPE>;
  inline Matrix() : Matrix( 0, 0 ) {};
  inline Matrix( std::size_t size1, std::size_t size2 );
  inline Matrix( const MyMatrix &o ) : Matrix() { *this = o; };
  inline Matrix( MyMatrix &&o ) : Matrix() { *this = o; };
  inline MyMatrix & operator=( const Matrix &o );
  inline MyMatrix & operator=( Matrix &&o );
  inline ~Matrix();
  inline std::size_t rows() const { return size1; }
  inline std::size_t cols() const { return size2; }
  inline void resize( std::size_t size1, std::size_t size2 );
  inline void clear();
  inline const Scalar & operator()( std::size_t i, std::size_t j ) const;
  inline Scalar & operator()( std::size_t i, std::size_t j );
  inline bool IsFinite( bool bDiagonalsOnly = false ) const;
  inline Real norm2() const;
  inline Real norm() const { return norm2(); }
  inline void blas_trmm( CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         Scalar alpha, const MyMatrix &A );
#ifdef COMMON_GSL_OPTIONAL
  inline MyMatrix Cholesky() const;
  inline MyMatrix inverse() const;
#endif
};

Matrix<COMMON_GSL_TYPE>::Matrix( std::size_t size1_, std::size_t size2_ )
{
  size1 = 0;
  size2 = 0;
  tda = 0;
  block = nullptr;
  data = nullptr;
  owner = false;
  resize( size1_, size2_ );
}

Matrix<COMMON_GSL_TYPE>::~Matrix()
{
  clear();
};

void Matrix<COMMON_GSL_TYPE>::clear()
{
  if( owner )
  {
    assert( block && "How did I get to be the owner of a null block?" );
    COMMON_GSL_FUNC( block, free )( block );
    owner = false;
  }
  block = nullptr;
  data = nullptr;
  size1 = 0;
  size2 = 0;
  tda = 0;
}

void Matrix<COMMON_GSL_TYPE>::resize( std::size_t size1_, std::size_t size2_ )
{
  if( size1_ != size1 || size2_ != size2 )
    clear();
  if( size1_ && size2_ )
  {
    std::size_t size = size1_ * size2_;
    assert( size / size1_ == size2_ && "Overflow allocating GSL matrix" );
    block = COMMON_GSL_FUNC( block, alloc )( size );
    data = block->data;
    size1 = size1_;
    size2 = size2_;
    tda = size1_;
    owner = true;
  }
}

Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator=( const Matrix &o )
{
  resize( o.size1, o.size2 );
  COMMON_GSL_FUNC( matrix, memcpy )( this, &o );
  return *this;
}

Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator=( Matrix &&o )
{
  size1 = o.size1;
  size2 = o.size2;
  tda = o.tda;
  data = o.data;
  block = o.block;
  owner = o.owner;
  o.owner = false;
  o.clear();
  return *this;
}

const Matrix<COMMON_GSL_TYPE>::Scalar & Matrix<COMMON_GSL_TYPE>::operator()( std::size_t i, std::size_t j ) const
{
  assert( data && i < size1 && j < size2 && "Error accessing matrix" );
  return reinterpret_cast<const Scalar *>( data )[i*tda+j];
}

Matrix<COMMON_GSL_TYPE>::Scalar & Matrix<COMMON_GSL_TYPE>::operator()( std::size_t i, std::size_t j )
{
  assert( data && i < size1 && j < size2 && "Error accessing matrix" );
  return reinterpret_cast<Scalar *>( data )[i*tda+j];
}

bool Matrix<COMMON_GSL_TYPE>::IsFinite( bool bDiagonalsOnly ) const
{
  for( std::size_t row = 0; row < size1; ++row )
    for( std::size_t col = bDiagonalsOnly ? row : 0 ; col < ( bDiagonalsOnly ? row + 1 : size2 ); ++col )
      if( !::Common::IsFinite( (*this)( row, col ) ) )
        return false;
  return true;
};

Matrix<COMMON_GSL_TYPE>::Real Matrix<COMMON_GSL_TYPE>::norm2() const
{
  assert( size1 && size2 && size1 == tda && "Still need to finish implementing norm2" );
  Vector<COMMON_GSL_TYPE> v;
  v.size = size1 * size2;
  v.data = data;
  v.stride = 1;
  return v.norm2();
}

void Matrix<COMMON_GSL_TYPE>::blas_trmm( CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                                         Scalar alpha, const Matrix<COMMON_GSL_TYPE> &A )
{
  COMMON_GSL_BLAS( trmm )( Side, Uplo, TransA, Diag, * reinterpret_cast<GSLScalar*>( &alpha ),
                           reinterpret_cast<const GSLMatrix *>( &A ), reinterpret_cast<GSLMatrix *>( this ) );
}

#ifdef COMMON_GSL_OPTIONAL
inline Matrix<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::Cholesky() const
{
  assert( size1 && size1 == size2 && "Cholesky decomposition of a non-square matrix" );
  Matrix<COMMON_GSL_TYPE> a( *this );
  //Vector<COMMON_GSL_TYPE> x( size1 * size2 );
  COMMON_GSL_FUNC( linalg, cholesky_decomp )( &a );
  return a;
}

Matrix<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::inverse() const
{
  Matrix<COMMON_GSL_TYPE> a( this->Cholesky() );
  COMMON_GSL_FUNC( linalg, cholesky_invert )( &a );
  return a;
}
#endif

inline std::ostream & operator<<( std::ostream &os, const Matrix<COMMON_GSL_TYPE> &m )
{
  os << "Matrix { " << m.size1 << " x " << m.size2;
  if( m.size1 != m.tda )
    os << " (tda=" << m.tda << ")";
  if( !m.owner )
    os << " not";
  os << " owner";
  const int xmax{ m.size1 >= 10 ? 10 : static_cast<int>( m.size1 ) };
  const int ymax{ m.size2 >= 10 ? 10 : static_cast<int>( m.size2 ) };
  for( int x = 0; x < xmax; ++x )
    for( int y = 0; y < ymax; ++y )
      os << ( y ? Space : NewLine ) << m(x,y);
  return os << " }";
}


#endif // COMMON_GSL_TYPE
