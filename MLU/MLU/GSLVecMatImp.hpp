/**

 Mike's lattice QCD utilities
 
 Source file: GSLVecMatImp.hpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

// Wrapper for GSL

#ifdef COMMON_GSL_TYPE

template <> struct Vector<COMMON_GSL_TYPE> : public GSLTraits<COMMON_GSL_TYPE>::GSLVectorType
{
  using Traits    = GSLTraits<COMMON_GSL_TYPE>;
  using Scalar    = Traits::Scalar;
  using value_type= Scalar;
  using Real      = Traits::Real;
  using GSLScalar = Traits::GSLScalar;
  using GSLBlock  = Traits::GSLBlockType;
  using GSLVector = Traits::GSLVectorType;
  using GSLMatrix = Traits::GSLMatrixType;
  using MyVector  = Vector<COMMON_GSL_TYPE>;
  using MyMatrix  = Matrix<COMMON_GSL_TYPE>;
  inline void NotEmpty() const { if( size == 0 ) throw std::runtime_error( "Vector empty" ); }
  inline void Compatible( const MyVector &right ) const
  {
    NotEmpty();
    if( size != right.size )
      throw std::runtime_error( "Vectors incompatible" );
  }
  inline Vector( std::size_t size );
  inline Vector() : Vector( 0 ) {};
  inline Vector( const MyVector &o );
  inline Vector( MyVector &&o ) : Vector() { *this = std::move( o ); };
  inline ~Vector();
  inline MyVector & operator=( Scalar c ); // Set every entry to the same value
  inline MyVector & operator=( const Vector &o );
  inline MyVector & operator=( Vector &&o );
  inline bool operator==( const Vector &o ) const;
  inline bool operator!=( const Vector &o ) const { return !operator==( o ); }
  template <typename T> inline MyVector &operator+=( T c );
  template <typename T> inline MyVector &operator-=( T c );
  template <typename T> inline MyVector &operator*=( T c );
  template <typename T> inline MyVector &operator/=( T c );
  inline void resize( std::size_t size );
  inline void clear();
  template<typename T=Scalar> inline typename std::enable_if<!is_complex<T>::value, T *>::type
    Data() { return data; }
  template<typename T=Scalar> inline typename std::enable_if< is_complex<T>::value, T *>::type
    Data() { return reinterpret_cast<T *>( data ); }
  template<typename T=Scalar> inline typename std::enable_if<!is_complex<T>::value, const T *>::type
    Data() const { return data; }
  template<typename T=Scalar> inline typename std::enable_if< is_complex<T>::value, const T *>::type
    Data() const { return reinterpret_cast<T *>( data ); }
  inline void MapView( Scalar * data_, std::size_t size_, std::size_t stride_ = 1 );
  inline void MapView( std::vector<Scalar> &v ) { MapView( v.data(), v.size() ); }
  inline void MapView( MyVector &v ) { MapView( v.Data(), v.size, v.stride ); }
  inline void MapRow( MyMatrix &m, std::size_t Row );
  inline void MapColumn( MyMatrix &m, std::size_t Column );
  inline const Scalar & operator[]( std::size_t i ) const;
  inline Scalar & operator[]( std::size_t i );
  inline std::size_t MinMax( Scalar &Min, Scalar &Max, bool bIgnoreNaN = false ) const;
#ifdef MLU_CBLAS
  const MLU_CBLAS_SCALAR_PTR DataCBLAS()const{return static_cast<const MLU_CBLAS_SCALAR_PTR>(data);}
        MLU_CBLAS_SCALAR_PTR DataCBLAS()     {return static_cast<      MLU_CBLAS_SCALAR_PTR>(data);}
  inline bool IsFinite() const;
  inline Real norm2() const
  {
    return MLU_CBLAS_RET_REAL( nrm2 )( static_cast<int>( size ), DataCBLAS(),
                                       static_cast<int>( stride ) );
  }
  inline Real norm() const { return norm2(); }
  inline Scalar Dot( const MyVector &right ) const;
  inline Scalar DotConj( const MyVector &right ) const; // conjugate( this ) . right
  inline void blas_trmv( CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag, const GSLMatrix &A );
#endif
};

Vector<COMMON_GSL_TYPE>::Vector( std::size_t size_ )
{
  size = 0;
  block = nullptr;
  data = nullptr;
  owner = false;
  if( size_ )
    resize( size_ );
}

Vector<COMMON_GSL_TYPE>::~Vector()
{
  clear();
};

// Copy constructor
Vector<COMMON_GSL_TYPE>::Vector( const Vector &o ) : Vector( 0 )
{
  if( o.owner )
  {
    // The other vector owns its memory - copy it
    *this = o;
  }
  else
  {
    // The original vector is a view - I'll be a view too
    size = o.size;
    stride = o.stride;
    data = o.data;
  }
}

Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator=( Scalar c )
{
  if( size == 0 )
    throw std::runtime_error( "Can't assign scalar to empty vector" );
  Scalar *p{ Data() };
  if( stride == 0 )
    * p = c;
  else
    for( std::size_t i = 0; i < size; ++i )
      p[i * stride] = c;
  return *this;
}

Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator=( const Vector &o )
{
  if( o.size == 0 )
    clear();
  else if( o.stride == 0 )
  {
    // Copy vector where all members are the same constant
    clear();
    resize( 1 );
    * Data() = * o.Data();
    size = o.size;
    stride = 0;
  }
  else
  {
    resize( o.size );
    COMMON_GSL_FUNC( vector, memcpy )( this, &o );
  }
  return *this;
}

Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator=( Vector &&o )
{
  clear();
  owner = o.owner;
  size = o.size;
  stride = o.stride;
  data = o.data;
  block = o.block;
  o.owner = false;
  o.clear();
  return *this;
}

template <typename T> Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator+=( T c )
{
  for( std::size_t i = 0; i < size; ++i )
    (*this)[i] += c;
  return *this;
}

template <typename T> Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator-=( T c )
{
  
  for( std::size_t i = 0; i < size; ++i )
    (*this)[i] -= c;
  return *this;
}

template <typename T> Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator*=( T c )
{
  for( std::size_t i = 0; i < size; ++i )
    (*this)[i] *= c;
  return *this;
}

template <typename T> Vector<COMMON_GSL_TYPE>& Vector<COMMON_GSL_TYPE>::operator/=( T c )
{
  for( std::size_t i = 0; i < size; ++i )
    (*this)[i] /= c;
  return *this;
}

bool Vector<COMMON_GSL_TYPE>::operator==( const Vector &o ) const
{
  return size == o.size && COMMON_GSL_FUNC( vector, equal )( this, &o );
}

void Vector<COMMON_GSL_TYPE>::resize( std::size_t size_ )
{
  // Same size or smaller and I can keep the existing buffer (and stride) - which might be a mapped view
  if( size_ && size_ <= size )
    size = size_;
  else if( size_ && owner && block->size >= size_ )
  {
    // Re-use existing buffer if it's big enough and I own it
    size = size_;
    stride = 1;
  }
  else
  {
    clear();
    if( size_ )
    {
      block = COMMON_GSL_FUNC( block, alloc )( size_ );
      if( !block )
        throw std::bad_alloc();
      data = block->data;
      size = size_;
      owner = true;
      stride = 1;
    }
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

std::size_t Vector<COMMON_GSL_TYPE>::MinMax( Scalar &Min, Scalar &Max, bool bIgnoreNaN ) const
{
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  std::size_t Count{};
  for( std::size_t i = 0; i < size; ++i )
  {
    const Scalar &z{ (*this)[i] };
    if( !bIgnoreNaN || ::MLU::IsFinite( z ) )
    {
      if( Count++ == 0 )
      {
        Min = z;
        Max = z;
      }
      else
      {
        Min = ComponentMin( Min, z );
        Max = ComponentMax( Max, z );
      }
    }
  }
  if( Count == 0 )
  {
    Min = NaN;
    Max = NaN;
  }
  return Count;
}

#ifdef MLU_CBLAS
bool Vector<COMMON_GSL_TYPE>::IsFinite() const
{
  for( std::size_t i = 0; i < size; ++i )
    if( !::MLU::IsFinite( (*this)[i] ) )
      return false;
  return true;
};

inline Vector<COMMON_GSL_TYPE>::Scalar Vector<COMMON_GSL_TYPE>::Dot( const MyVector &right ) const
{
  Compatible( right );
#ifdef COMMON_GSL_IS_COMPLEX
    Scalar z;
    MLU_CBLAS( dotu_sub )( static_cast<int>( size ), DataCBLAS(), static_cast<int>( stride ),
                           right.DataCBLAS(), static_cast<int>( right.stride ),
                           static_cast<void *>( &z ) );
    return z;
#else
  return MLU_CBLAS( dot )( static_cast<int>( size ), DataCBLAS(), static_cast<int>( stride ),
                           right.DataCBLAS(), static_cast<int>( right.stride ) );
#endif
}

inline Vector<COMMON_GSL_TYPE>::Scalar Vector<COMMON_GSL_TYPE>::DotConj( const MyVector &right ) const
{
  Compatible( right );
#ifdef COMMON_GSL_IS_COMPLEX
    Scalar z;
    MLU_CBLAS( dotc_sub )( static_cast<int>( size ), DataCBLAS(), static_cast<int>( stride ),
                           right.DataCBLAS(), static_cast<int>( right.stride ),
                           static_cast<void *>( &z ) );
    return z;
#else
  return MLU_CBLAS( dot )( static_cast<int>( size ), DataCBLAS(), static_cast<int>( stride ),
                           right.DataCBLAS(), static_cast<int>( right.stride ) );
#endif
}

inline void Vector<COMMON_GSL_TYPE>::blas_trmv( CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                                                const Vector<COMMON_GSL_TYPE>::GSLMatrix &A )
{
  if( size == 0 || A.size2 != size || A.size1 != A.size2 )
    throw std::runtime_error( "MLU::blas_trmv incompatible matrix and vector" );
  const int N{ static_cast<int>( size ) };
  const int ldA{ static_cast<int>( A.tda ) };
  const int incX{ static_cast<int>( stride ) };
  MLU_CBLAS( trmv )( CblasRowMajor, Uplo, TransA, Diag,
                     N, reinterpret_cast<const MLU_CBLAS_SCALAR_PTR>( A.data ), ldA,
                     reinterpret_cast<MLU_CBLAS_SCALAR_PTR>( data ), incX );
}
#endif

inline std::ostream & operator<<( std::ostream &os, const Vector<COMMON_GSL_TYPE> &v )
{
  static constexpr int MaxCols = 8;
  os << "Vector { " << v.size;
  if( v.stride != 1 )
    os << " (stride=" << v.stride << ")";
  if( !v.owner )
    os << " not";
  os << " owner";
  const int imax{ v.size >= MaxCols ? MaxCols : static_cast<int>( v.size ) };
  for( int i = 0; i < imax; ++i )
      os << " " << v[i];
  return os << " }";
}

template <> struct Matrix<COMMON_GSL_TYPE> : public GSLTraits<COMMON_GSL_TYPE>::GSLMatrixType
{
  using Traits    = GSLTraits<COMMON_GSL_TYPE>;
  using Scalar    = Traits::Scalar;
  using value_type= Scalar;
  using Real      = Traits::Real;
  using GSLScalar = Traits::GSLScalar;
  using GSLBlock  = Traits::GSLBlockType;
  using GSLVector = Traits::GSLVectorType;
  using GSLMatrix = Traits::GSLMatrixType;
  using MyVector  = Vector<COMMON_GSL_TYPE>;
  using MyMatrix  = Matrix<COMMON_GSL_TYPE>;
  inline bool Empty() const { return size1 == 0 || size2 == 0; }
  inline void NotEmpty() const { if( Empty() ) throw std::runtime_error( "Matrix empty" ); }
  inline void Square() const { NotEmpty(); if( size1 != size2 ) throw std::runtime_error( "Matrix not square");}
  inline Matrix() : Matrix( 0, 0 ) {};
  inline Matrix( std::size_t size1, std::size_t size2 );
  inline Matrix( const MyMatrix &o );
  inline Matrix( MyMatrix &&o ) : Matrix() { *this = std::move( o ); };
  inline ~Matrix();
  inline MyMatrix & operator=( Scalar c );
  inline MyMatrix & operator=( const Matrix &o );
  inline MyMatrix & operator=( Matrix &&o );
  template <typename T> inline MyMatrix &operator+=( T c );
  template <typename T> inline MyMatrix &operator-=( T c );
  template <typename T> inline MyMatrix &operator*=( T c );
  template <typename T> inline MyMatrix &operator/=( T c );
  inline bool operator==( const Matrix &o ) const;
  inline bool operator!=( const Matrix &o ) const { return !operator==( o ); }
  inline std::size_t rows() const { return size1; }
  inline std::size_t cols() const { return size2; }
  inline void resize( std::size_t size1, std::size_t size2 );
  inline void clear();
  template<typename T=Scalar> inline typename std::enable_if<!is_complex<T>::value, T *>::type
    Data() { return data; }
  template<typename T=Scalar> inline typename std::enable_if< is_complex<T>::value, T *>::type
    Data() { return reinterpret_cast<T *>( data ); }
  template<typename T=Scalar> inline typename std::enable_if<!is_complex<T>::value, const T *>::type
    Data() const { return data; }
  template<typename T=Scalar> inline typename std::enable_if< is_complex<T>::value, const T *>::type
    Data() const { return reinterpret_cast<T *>( data ); }
  inline void MapView( Scalar * data_, std::size_t size1_, std::size_t size2_, std::size_t tda_ );
  inline void MapView( Scalar * data_, std::size_t size1_, std::size_t size2_ )
  { MapView( data_, size1_, size2_, size2_ ); }
  inline void MapView( MyMatrix &o )
  { MapView( o.Data(), o.size1, o.size2, o.tda ); };
  inline void MapView( MyVector &o )
  {
    if( o.size == 0 )
      throw std::runtime_error( "Cannot map matrix view of empty vector" );
    if( o.stride != 1 )
      throw std::runtime_error( "Cannot map matrix view of vector with stride "
                               + std::to_string( o.stride ) );
    MapView( o.Data(), 1, o.size );
  };
  inline MyMatrix SubMatrix( std::size_t i, std::size_t j, std::size_t rows, std::size_t cols );
  inline const Scalar & operator()( std::size_t i, std::size_t j ) const;
  inline Scalar & operator()( std::size_t i, std::size_t j );
  inline void Row( std::size_t idx, MyVector &v );
  inline void Column( std::size_t idx, MyVector &v );
  inline void SetDiagonalOne();
  inline void ZeroUpperTriangle();
#ifdef MLU_CBLAS
  const MLU_CBLAS_SCALAR_PTR DataCBLAS()const{return static_cast<const MLU_CBLAS_SCALAR_PTR>(Data());}
        MLU_CBLAS_SCALAR_PTR DataCBLAS()     {return static_cast<      MLU_CBLAS_SCALAR_PTR>(Data());}
  inline bool IsFinite( bool bDiagonalsOnly = false ) const;
  inline Real norm2() const;
  inline Real norm() const { return norm2(); }
  inline void blas_gemm( CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, Scalar alpha,
                        const MyMatrix &A, const MyMatrix &B, Scalar beta );
  inline void blas_symm( CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, Scalar alpha, const MyMatrix &A,
                         const MyMatrix &B, Scalar beta );
  inline void blas_trmm( CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                         Scalar alpha, const MyMatrix &A );
#endif
#ifdef COMMON_GSL_DOUBLE
  inline MyMatrix Cholesky( MyVector &S ) const;
  inline MyVector CholeskySolve( const MyVector &b ) const;
  inline MyVector CholeskySolve( const MyVector &b, const MyVector &S ) const;
  inline Scalar CholeskyRCond() const;
  inline MyVector CholeskyScale() const;
  inline MyVector CholeskyExtract( bool bSetDiagonalToOne = true );
  inline MyVector GetEigenValues( MyMatrix *pEigenVectors = nullptr ) const;
#endif
  inline void CholeskyScaleApply( const MyVector &S, bool bSetDiagonalToOne = false );
#ifdef COMMON_GSL_OPTIONAL
  inline MyMatrix Inverse() const;
  inline bool Cholesky( bool bZeroUpperTriangle );
  inline void CholeskyInvert();
  //inline void TriangularInvert( CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag );
#endif
  inline void SaveSquare( const std::string &Filename, const std::string &Type,
                 std::vector<std::string> &ColNames, const char *pGnuplotExtra = nullptr ) const;
};

void Vector<COMMON_GSL_TYPE>::MapRow( MyMatrix &m, std::size_t Row )
{
  if( m.size1 == 0 || m.size2 == 0 )
    throw std::runtime_error( "Vector::MapRow() empty matrix" );
  if( Row >= m.size1 )
    throw std::runtime_error( "Row " + std::to_string( Row ) + " >= " + std::to_string( m.size1 ) );
  MapView( reinterpret_cast<Scalar *>( m.data ) + m.tda * Row, m.size2, 1 );
}

void Vector<COMMON_GSL_TYPE>::MapColumn( MyMatrix &m, std::size_t Column )
{
  if( m.size1 == 0 || m.size2 == 0 )
    throw std::runtime_error( "Vector::MapColumn() empty matrix" );
  if( Column >= m.size2 )
    throw std::runtime_error( "Column " + std::to_string( Column ) + " > " + std::to_string( m.size2 ) );
  MapView( reinterpret_cast<Scalar *>( m.data ) + Column, m.size1, m.tda );
}

Matrix<COMMON_GSL_TYPE>::Matrix( std::size_t size1_, std::size_t size2_ )
{
  size1 = 0;
  size2 = 0;
  tda = 0;
  block = nullptr;
  data = nullptr;
  owner = false;
  if( size1_ * size2_ )
    resize( size1_, size2_ );
}

Matrix<COMMON_GSL_TYPE>::~Matrix()
{
  clear();
};

// Copy constructor
Matrix<COMMON_GSL_TYPE>::Matrix( const Matrix &o ) : Matrix( 0, 0 )
{
  if( o.owner )
  {
    // The other matrix owns its memory - copy it
    *this = o;
  }
  else
  {
    // The original matrix is a view - I'll be a view too
    size1 = o.size1;
    size2 = o.size2;
    tda = o.tda;
    data = o.data;
  }
}

Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator=( Scalar c )
{
  if( size1 == 0 || size1 != size2 )
    throw std::runtime_error( "Can't assign scalar to empty matrix" );
  for( int i = 0; i < size1; ++i )
    for( int j = 0; j < size2; ++j )
      (*this)( i, j ) = ( i == j ) ? c : 0;
  return *this;
}

Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator=( const Matrix &o )
{
  resize( o.size1, o.size2 );
  if( size1 && size2 )
    COMMON_GSL_FUNC( matrix, memcpy )( this, &o );
  return *this;
}

Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator=( Matrix &&o )
{
  clear();
  owner = o.owner;
  size1 = o.size1;
  size2 = o.size2;
  tda = o.tda;
  data = o.data;
  block = o.block;
  o.owner = false;
  o.clear();
  return *this;
}

template <typename T> Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator+=( T c )
{
  for( std::size_t i = 0; i < size1; ++i )
    for( std::size_t j = 0; j < size2; ++j )
      (*this)(i,j) += c;
  return *this;
}

template <typename T> Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator-=( T c )
{
  for( std::size_t i = 0; i < size1; ++i )
    for( std::size_t j = 0; j < size2; ++j )
      (*this)(i,j) -= c;
  return *this;
}

template <typename T> Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator*=( T c )
{
  for( std::size_t i = 0; i < size1; ++i )
    for( std::size_t j = 0; j < size2; ++j )
      (*this)(i,j) *= c;
  return *this;
}

template <typename T> Matrix<COMMON_GSL_TYPE>& Matrix<COMMON_GSL_TYPE>::operator/=( T c )
{
  for( std::size_t i = 0; i < size1; ++i )
    for( std::size_t j = 0; j < size2; ++j )
      (*this)(i,j) /= c;
  return *this;
}

bool Matrix<COMMON_GSL_TYPE>::operator==( const Matrix &o ) const
{
  return size1 == o.size1 && size2 == o.size2 && COMMON_GSL_FUNC( matrix, equal )( this, &o );
}

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

void Matrix<COMMON_GSL_TYPE>::MapView(Scalar * data_, std::size_t size1_, std::size_t size2_, std::size_t tda_)
{
  clear();
  size1 = size1_;
  size2 = size2_;
  tda = tda_;
  data = reinterpret_cast<Real *>( data_ );
}


void Matrix<COMMON_GSL_TYPE>::resize( std::size_t size1_, std::size_t size2_ )
{
  std::size_t size = size1_ * size2_;
  assert( size == size1_ * size2_ && "Overflow allocating GSL matrix" );
  // Re-use existing buffer if it's big enough and I own it
  if( size && owner && block->size >= size )
  {
    size1 = size1_;
    size2 = size2_;
    tda = size2_;
  }
  else
  {
    clear();
    if( size )
    {
      block = COMMON_GSL_FUNC( block, alloc )( size );
      if( !block )
        throw std::bad_alloc();
      data = block->data;
      size1 = size1_;
      size2 = size2_;
      tda = size2_;
      owner = true;
    }
  }
}

Matrix<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::SubMatrix( std::size_t i, std::size_t j,
                                                            std::size_t rows, std::size_t cols )
{
  if( i >= size1 )
    throw std::runtime_error( "Row " + std::to_string( i ) + " >= " + std::to_string( size1 ) );
  if( j >= size2 )
    throw std::runtime_error( "Column " + std::to_string( j ) + " >= " + std::to_string( size2 ) );
  if( i + rows > size1 )
    throw std::runtime_error( "Rows " + std::to_string( rows ) + " too high" );
  if( j + cols > size2 )
    throw std::runtime_error( "Columns " + std::to_string( cols ) + " too wide" );
  MyMatrix o;
  o.MapView( Data() + i * tda + j, rows, cols, tda );
  return o;
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

void Matrix<COMMON_GSL_TYPE>::Row( std::size_t idx, MyVector &v )
{
  if( idx >= size1 )
    throw std::runtime_error( "Row " + std::to_string( idx ) + " >= " + std::to_string( size1 ) );
  v.MapView( reinterpret_cast<Scalar *>( data ) + tda * idx, size2, 1 );
}

void Matrix<COMMON_GSL_TYPE>::Column( std::size_t idx, MyVector &v )
{
  if( idx >= size2 )
    throw std::runtime_error( "Column " + std::to_string( idx ) + " > " + std::to_string( size2 ) );
  v.MapView( reinterpret_cast<Scalar *>( data ) + idx, size1, tda );
}

inline void Matrix<COMMON_GSL_TYPE>::SetDiagonalOne()
{
  Square();
  for( std::size_t i = 0; i < size1; ++i )
    ( *this ) ( i, i ) = 1;
}

inline void Matrix<COMMON_GSL_TYPE>::ZeroUpperTriangle()
{
  Square();
  for( std::size_t i = 0; i < size1; ++i )
    for( std::size_t j = i + 1; j < size2; ++j )
      ( *this ) ( i, j ) = 0;
}

#ifdef MLU_CBLAS
bool Matrix<COMMON_GSL_TYPE>::IsFinite( bool bDiagonalsOnly ) const
{
  for( std::size_t row = 0; row < size1; ++row )
    for( std::size_t col = bDiagonalsOnly ? row : 0 ; col < ( bDiagonalsOnly ? row + 1 : size2 ); ++col )
      if( !::MLU::IsFinite( (*this)( row, col ) ) )
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

// https://developer.apple.com/documentation/accelerate/1513282-cblas_dgemm?language=objc
// this = C

void Matrix<COMMON_GSL_TYPE>::blas_gemm( CBLAS_TRANSPOSE_t TransA, CBLAS_TRANSPOSE_t TransB, Scalar alpha,
                                         const Matrix<COMMON_GSL_TYPE> &A, const Matrix<COMMON_GSL_TYPE> &B,
                                         Scalar beta )
{
  NotEmpty();
  if( size1 != A.size1 )
    throw std::runtime_error( "blas_gemm() A doesn't have " + std::to_string( size1 ) + " rows" );
  if( size2 != B.size2 )
    throw std::runtime_error( "blas_gemm() B doesn't have " + std::to_string( size2 ) + " columns" );
  if( A.size2 != B.size1 )
    throw std::runtime_error( "blas_gemm() A and B not correct shape for multiply" );
  const int M{ static_cast<int>( size1 ) };
  const int N{ static_cast<int>( size2 ) };
  const int K{ static_cast<int>( A.size2 ) };
  const int ldA{ static_cast<int>( A.tda ) };
  const int ldB{ static_cast<int>( B.tda ) };
  const int ldC{ static_cast<int>( tda ) };
#ifdef COMMON_GSL_IS_COMPLEX
  MLU_CBLAS( gemm )( CblasRowMajor, TransA, TransB, M, N, K,
                     reinterpret_cast<MLU_CBLAS_SCALAR_PTR>(&alpha), A.DataCBLAS(), ldA,
                     B.DataCBLAS(), ldB, reinterpret_cast<MLU_CBLAS_SCALAR_PTR>(&beta),
                     DataCBLAS(), ldC );
#else
  MLU_CBLAS( gemm )( CblasRowMajor, TransA, TransB, M, N, K,
                     alpha, A.DataCBLAS(), ldA,
                     B.DataCBLAS(), ldB, beta,
                     DataCBLAS(), ldC );
#endif
}

// https://www.ibm.com/docs/en/essl/6.2?topic=mos-ssymm-dsymm-csymm-zsymm-chemm-zhemm-matrix-matrix-product-where-one-matrix-is-real-complex-symmetric-complex-hermitian
// this = C, A is symmetric, C = alpha AB + beta C, or C = alpha BA + beta C

void Matrix<COMMON_GSL_TYPE>::blas_symm( CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo,
                                         Scalar alpha, const Matrix<COMMON_GSL_TYPE> &A,
                                         const Matrix<COMMON_GSL_TYPE> &B, Scalar beta )
{
  NotEmpty();
  if( Side != CblasLeft && Side != CblasRight )
    throw std::runtime_error( "blas_symm() Side invalid" );
  if( Uplo != CblasUpper && Uplo != CblasLower )
    throw std::runtime_error( "blas_symm() Uplo invalid" );
  if( A.size1 == 0 || A.size1 != A.size2 )
    throw std::runtime_error( "blas_symm() A " + std::to_string( size1 ) + " x "
                             + std::to_string( size2 ) + " isn't symmetric" );
  if( Side == CblasLeft )
  {
    // AB
    if( size1 != A.size1 )
      throw std::runtime_error( "blas_symm() A doesn't have " + std::to_string( size1 ) + " rows" );
    if( size2 != B.size2 )
      throw std::runtime_error( "blas_symm() B doesn't have " + std::to_string( size2 ) + " columns" );
    if( A.size2 != B.size1 )
      throw std::runtime_error( "blas_symm() A and B not correct shape for multiply" );
  }
  else
  {
    // BA
    if( size1 != B.size1 )
      throw std::runtime_error( "blas_symm() B doesn't have " + std::to_string( size1 ) + " rows" );
    if( size2 != A.size2 )
      throw std::runtime_error( "blas_symm() A doesn't have " + std::to_string( size2 ) + " columns" );
    if( B.size2 != A.size1 )
      throw std::runtime_error( "blas_symm() A and B not correct shape for multiply" );
  }
  const int M{ static_cast<int>( size1 ) };
  const int N{ static_cast<int>( size2 ) };
  const int ldA{ static_cast<int>( A.tda ) };
  const int ldB{ static_cast<int>( B.tda ) };
  const int ldC{ static_cast<int>( tda ) };
#ifdef COMMON_GSL_IS_COMPLEX
  MLU_CBLAS( symm )( CblasRowMajor, Side, Uplo, M, N,
                     reinterpret_cast<MLU_CBLAS_SCALAR_PTR>(&alpha), A.DataCBLAS(), ldA,
                     B.DataCBLAS(), ldB, reinterpret_cast<MLU_CBLAS_SCALAR_PTR>(&beta),
                     DataCBLAS(), ldC );
#else
  MLU_CBLAS( symm )( CblasRowMajor, Side, Uplo, M, N,
                     alpha, A.DataCBLAS(), ldA,
                     B.DataCBLAS(), ldB, beta,
                     DataCBLAS(), ldC );
#endif
}

// https://developer.apple.com/documentation/accelerate/1513037-cblas_strmm?language=objc
// this = B, A is triangular, B = alpha AB, or B = alpha BA

void Matrix<COMMON_GSL_TYPE>::blas_trmm( CBLAS_SIDE_t Side, CBLAS_UPLO_t Uplo, CBLAS_TRANSPOSE_t TransA, CBLAS_DIAG_t Diag,
                                         Scalar alpha, const Matrix<COMMON_GSL_TYPE> &A )
{
  NotEmpty();
  if( Side != CblasLeft && Side != CblasRight )
    throw std::runtime_error( "blas_trmm() Side invalid" );
  if( Uplo != CblasUpper && Uplo != CblasLower )
    throw std::runtime_error( "blas_trmm() Uplo invalid" );
  if( TransA != CblasNoTrans && TransA != CblasTrans && TransA != CblasConjTrans )
    throw std::runtime_error( "blas_trmm() TransA invalid" );
  if( Diag != CblasUnit && Diag != CblasNonUnit )
    throw std::runtime_error( "blas_trmm() Diag invalid" );
  if( A.size1 == 0 || A.size1 != A.size2 )
    throw std::runtime_error( "blas_trmm() A " + std::to_string( size1 ) + " x "
                             + std::to_string( size2 ) + " isn't triangular" );
  if( Side == CblasLeft )
  {
    // AB
    if( A.size2 != size1 )
      throw std::runtime_error( "blas_trmm() A doesn't have " + std::to_string( size1 ) + " columns");
  }
  else
  {
    // BA
    if( size2 != A.size1 )
      throw std::runtime_error( "blas_trmm() A doesn't have " + std::to_string( size2 ) + " rows" );
  }
  const int M{ static_cast<int>( size1 ) };
  const int N{ static_cast<int>( size2 ) };
  const int ldA{ static_cast<int>( A.tda ) };
  const int ldB{ static_cast<int>( tda ) };
#ifdef COMMON_GSL_IS_COMPLEX
  MLU_CBLAS( trmm )( CblasRowMajor, Side, Uplo, TransA, Diag, M, N,
                     reinterpret_cast<MLU_CBLAS_SCALAR_PTR>(&alpha), A.DataCBLAS(), ldA,
                     DataCBLAS(), ldB );
#else
  MLU_CBLAS( trmm )( CblasRowMajor, Side, Uplo, TransA, Diag, M, N,
                     alpha, A.DataCBLAS(), ldA,
                     DataCBLAS(), ldB );
#endif
}
#endif

/*
 Calling this matrix "A" and returned matrix "R" ...
 Input:   lower left and diagonals of A contain matrix to Cholesky decompose (assumed symmetric)
 Output:  diag(S) = scaling matrix so that (S A S) has 1 on diagonals, i.e. S_i = 1/sqrt{A_ii}
          R upper right not diagonal = S A S. This is symmetric, and diagonal is assumed to be 1
          R lower left  and diagonal = Cholesky decomp of S A S (triangular, so upper right assumed to be 0)
 Example: If A is a variance matrix on input, then R contains the correlation matrix (upper right)
          and the Cholesky decomposition of the correlation matrix (lower left and diagonals),
          while S contains the inverse of the errors
 */

#ifdef COMMON_GSL_DOUBLE
inline Matrix<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::Cholesky( MyVector &S ) const
{
  Square();
  MyMatrix a;
  a = *this; // TODO: If I do this in the constructor and I'm a view, a won't be a copy
  S.resize( size1 );
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_decomp2 )( &a, &S ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to Cholesky decompose2 matrix", __FILE__, __LINE__, gsl_e );
  return a;
}

inline Vector<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::CholeskySolve( const MyVector &b ) const
{
  MyVector x( size1 );
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_solve )( this, &b, &x ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to solve Ax=b using Cholesky decomposition", __FILE__, __LINE__, gsl_e );
  return x;
}

inline Vector<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::CholeskySolve( const MyVector &b, const MyVector &S ) const
{
  MyVector x( size1 );
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_solve2 )( this, &S, &b, &x ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to solve2 Ax=b using Cholesky decomposition", __FILE__, __LINE__, gsl_e );
  return x;
}

inline COMMON_GSL_TYPE Matrix<COMMON_GSL_TYPE>::CholeskyRCond() const
{
  Square();
  MyVector v( size1 * 3 );
  double rcond;
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_rcond )( this, &rcond, &v ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to estimate reciprocal condition number", __FILE__, __LINE__, gsl_e );
  return rcond;
}

inline Vector<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::CholeskyScale() const
{
  MyVector S( size1 );
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_scale )( this, &S ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to Cholesky scale", __FILE__, __LINE__, gsl_e );
  return S;
}

inline void Matrix<COMMON_GSL_TYPE>::CholeskyScaleApply( const MyVector &S, bool bSetDiagonalToOne )
{
  Square();
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_scale_apply )( this, &S ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to apply Cholesky scale", __FILE__, __LINE__, gsl_e );
  for( int i = 0; i < size1; ++i )
    for( int j = 0; j < i; ++j )
      (*this)( j, i ) = (*this)( i, j );
  if( bSetDiagonalToOne )
    SetDiagonalOne();
}

inline Vector<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::CholeskyExtract( bool bSetDiagonalToOne )
{
  Square();
  MyVector S{ CholeskyScale() };
  CholeskyScaleApply( S );
  if( bSetDiagonalToOne )
    SetDiagonalOne();
  return S;
}

inline Vector<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::GetEigenValues( MyMatrix *pEigenVectors ) const
{
  Square();
  MyVector EigenValues{ size1 };
  MyMatrix WorkingCopy( *this );
  int gsl_e;
  if( pEigenVectors )
  {
    pEigenVectors->resize( size1, size2 );
    gsl_eigen_symmv_workspace *ws{ gsl_eigen_symmv_alloc( size1 ) };
    if( !ws )
      throw std::bad_alloc();
    gsl_e = gsl_eigen_symmv( &WorkingCopy, &EigenValues, pEigenVectors, ws );
    gsl_eigen_symmv_free( ws );
  }
  else
  {
    gsl_eigen_symm_workspace *ws{ gsl_eigen_symm_alloc( size1 ) };
    if( !ws )
      throw std::bad_alloc();
    gsl_e = gsl_eigen_symm( &WorkingCopy, &EigenValues, ws );
    gsl_eigen_symm_free( ws );
  }
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to GetEigenValues()", __FILE__, __LINE__, gsl_e );
  return EigenValues;
}
#else
inline void Matrix<COMMON_GSL_TYPE>::CholeskyScaleApply( const MyVector &S, bool bSetDiagonalToOne )
{
  throw std::runtime_error( "CholeskyScaleApply() undefined except for double" );
}
#endif // COMMON_GSL_DOUBLE

#ifdef COMMON_GSL_OPTIONAL
// Dump debugging info for real-symmetric matrix: eigenvalues and pivoted Cholesky decomposition diagonal
Matrix<COMMON_GSL_TYPE> Matrix<COMMON_GSL_TYPE>::Inverse() const
{
  MyMatrix a( *this );
  a.Cholesky( false );
  a.CholeskyInvert();
  return a;
}

inline bool Matrix<COMMON_GSL_TYPE>::Cholesky( bool bZeroUpperTriangle )
{
  Square();
  const bool bOK{ COMMON_GSL_FUNC( linalg, cholesky_decomp )( this ) == 0 };
  if( bOK && bZeroUpperTriangle )
    ZeroUpperTriangle();
  return bOK;
}

void Matrix<COMMON_GSL_TYPE>::CholeskyInvert()
{
  int gsl_e{ COMMON_GSL_FUNC( linalg, cholesky_invert )( this ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to Cholesky invert matrix", __FILE__, __LINE__, gsl_e );
}

/*void Matrix<COMMON_GSL_TYPE>::TriangularInvert( CBLAS_UPLO_t Uplo, CBLAS_DIAG_t Diag )
{
  int gsl_e{ COMMON_GSL_FUNC( linalg, tri_invert )( Uplo, Diag, this ) };
  if( gsl_e )
    GSLLibraryGlobal::Error( "Unable to triangular invert matrix", __FILE__, __LINE__, gsl_e );
  }*/

#endif // COMMON_GSL_OPTIONAL

// Write square matrix to file
void Matrix<COMMON_GSL_TYPE>::SaveSquare( const std::string &Filename, const std::string &Type,
                                          std::vector<std::string> &ColNames,
                                          const char *pGnuplotExtra ) const
{
  if( size1 == 0 )
    throw std::runtime_error( "SaveSquare() empty " + Type + " matrix " + Filename );
  if( size1 != size2 )
    throw std::runtime_error( "SaveSquare() Non-square " + Type + " matrix " + Filename );
  // Header describing what the covariance matrix is for and how to plot it with gnuplot
  std::ofstream s{ Filename };
  s << "# Matrix: " << Filename << "\n# Type: " << Type << "\n";
  // Save a command which can be used to plot this file
  s << "# gnuplot: set xtics rotate noenhanced font 'Arial,4'; set ytics noenhanced font 'Arial,4';"
  << " set title noenhanced '" << Filename << "'; unset key; ";
  if( pGnuplotExtra && *pGnuplotExtra )
    s << pGnuplotExtra << "; ";
  s << "plot '" << Filename << "' matrix columnheaders rowheaders with image pixels\n";
  // Now save the column names. First column is matrix extent
  s << "# The next row contains column headers, starting with matrix dimension\n";
  s << size1;
  for( std::size_t i = 0; i < size1; ++i )
  {
    s << " ";
    if( i < ColNames.size() )
      s << ColNames[i];
    else
      s << "c" << i;
  }
  s << "\n" << std::setprecision( std::numeric_limits<Scalar>::max_digits10 );
  // Now print the actual covariance matrix
  for( std::size_t i = 0; i < size1; ++i )
  {
    if( i < ColNames.size() )
      s << ColNames[i];
    else
      s << "c" << i;
    for( int j = 0; j < size2; ++j )
      s << " " << (*this)( i, j );
    s << "\n";
  }
}

inline std::ostream & operator<<( std::ostream &os, const Matrix<COMMON_GSL_TYPE> &m )
{
  static constexpr int MaxRows = 8;
  static constexpr int MaxCols = 8;
  os << "Matrix { " << m.size1 << " x " << m.size2;
  if( m.size2 != m.tda )
    os << " (tda=" << m.tda << ")";
  if( !m.owner )
    os << " not";
  os << " owner";
  const int xmax{ m.size1 >= MaxRows ? MaxRows : static_cast<int>( m.size1 ) };
  const int ymax{ m.size2 >= MaxCols ? MaxCols : static_cast<int>( m.size2 ) };
  for( int x = 0; x < xmax; ++x )
    for( int y = 0; y < ymax; ++y )
      os << ( y ? " " : "\n" ) << m(x,y);
  return os << " }";
}

#undef COMMON_GSL_TYPE
#undef COMMON_GSL_FUNC
#undef COMMON_GSL_DOUBLE
#undef COMMON_GSL_IS_COMPLEX
#undef COMMON_GSL_OPTIONAL
#undef MLU_CBLAS
#undef MLU_CBLAS_RET_REAL
#undef MLU_CBLAS_SCALAR_PTR

#endif // COMMON_GSL_TYPE
