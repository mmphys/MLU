/*****************************************************************************
 
 HDF5 support

 Source file: HDF5.hpp
 
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

#ifndef MLU_HDF5_hpp
#define MLU_HDF5_hpp namespace Common {
#define MLU_HDF5_hpp_end };

#include <MLUconfig.h>
#include <MLU/GSLVecMat.hpp>
#include <MLU/Math.hpp>
#include <MLU/Physics.hpp>
#include <MLU/String.hpp>

// HDF5 Library
#include <H5Cpp.h>
#include <H5CompType.h>
#include <H5public.h>

// Default output file extension for binary data
#ifndef DEF_FMT
#define DEF_FMT "h5"
#endif

MLU_HDF5_hpp

// My implementation of H5File - adds a definition of complex type
namespace H5 {
  template<typename T> struct is_AttOrDS        : public std::false_type {};
  template<> struct is_AttOrDS<::H5::Attribute> : public std::true_type {};
  template<> struct is_AttOrDS<::H5::DataSet>   : public std::true_type {};

  template <typename T> struct Equiv;
  template<> struct Equiv<float>                        { static const ::H5::PredType& Type; };
  template<> struct Equiv<double>                       { static const ::H5::PredType& Type; };
  template<> struct Equiv<long double>                  { static const ::H5::PredType& Type; };
  template<> struct Equiv<int>                          { static const ::H5::PredType& Type; };
  template<> struct Equiv<std::size_t>                  { static const ::H5::PredType& Type; };
  template<> struct Equiv<std::string>                  { static const ::H5::StrType   Type; };
  template<> struct Equiv<char *>                       { static const ::H5::StrType&  Type; };
#ifndef HAVE_FAST32_IS_SIZE_T
  template<> struct Equiv<std::uint_fast32_t>           { static const ::H5::PredType& Type; };
#endif
  template<> struct Equiv<std::complex<float>>          { static const ::H5::CompType  Type; };
  template<> struct Equiv<std::complex<double>>         { static const ::H5::CompType  Type; };
  template<> struct Equiv<std::complex<long double>>    { static const ::H5::CompType  Type; };
  template<> struct Equiv<ValWithEr<float>>             { static const ::H5::CompType  Type; };
  template<> struct Equiv<ValWithEr<double>>            { static const ::H5::CompType  Type; };
  template<> struct Equiv<ValWithEr<long double>>       { static const ::H5::CompType  Type; };
  template<> struct Equiv<ValWithErOldV1<float>>        { static const ::H5::CompType  Type; };
  template<> struct Equiv<ValWithErOldV1<double>>       { static const ::H5::CompType  Type; };
  template<> struct Equiv<ValWithErOldV1<long double>>  { static const ::H5::CompType  Type; };
  template<> struct Equiv<ConfigCount>                  { static const ::H5::CompType  Type; };

  /**
   Open the specified HDF5File
   Read from the specified group name (if pGroupName is non-null and *pGroupName not empty)
   Otherwise read from the first group in the file, returning the group name if pGroupName not null
   */
  void OpenFileGroup( ::H5::H5File &f, ::H5::Group &g,
                      const std::string &FileName, const char *PrintPrefix = nullptr,
                      std::string * pGroupName = nullptr, unsigned int flags = H5F_ACC_RDONLY );

  // Get first groupname from specified group
  std::string GetFirstGroupName( ::H5::Group & g );

  bool OpenOptional( ::H5::Attribute &a, ::H5::Group &g, const std::string Name );
  bool OpenOptional( ::H5::DataSet  &ds, ::H5::Group &g, const std::string Name );
  bool OpenOptional( ::H5::Group  &gNew, ::H5::Group &g, const std::string Name );

  // Read the gamma algebra attribute string and make sure it's valid
  Gamma::Algebra ReadGammaAttribute( ::H5::Group &g, const char * pAttName );

  template<typename AttOrDS> std::vector<std::string> ReadStrings( const AttOrDS &a );

  // Make a multi-dimensional string attribute
  void WriteAttribute( ::H5::Group &g, const std::string &AttName, const std::vector<std::string> &vs );
  // Make a multi-dimensional string dataset
  void WriteStringData( ::H5::Group &g, const std::string &DSName, const std::vector<std::string> &vs );
  // Read a vector from a dataset
  // TODO: Rationalise
  inline void ReadAttOrDS( ::H5::Attribute &a, void * p, const ::H5::DataType &dt )
  { a.read( dt, p ); }
  inline void ReadAttOrDS( ::H5::DataSet  &ds, void * p, const ::H5::DataType &dt )
  { ds.read( p, dt ); }
  template<typename T, typename AttOrDS> typename std::enable_if<is_AttOrDS<AttOrDS>::value>::type
  ReadVectorHelper( T &v, AttOrDS &a, const std::string &DSName );
  template<typename T>
  void ReadVector( ::H5::Group &g, const std::string &DSName, T &v );
  template<typename T>
  void ReadVector( ::H5::Group &g, const std::string &DSName, std::vector<Common::Matrix<T>> &v );
  // Write a vector to a dataset
  template<typename T>
  bool WriteVector( ::H5::Group &g, const std::string &DSName, const Common::Vector<T> &v );
  template<typename T>
  bool WriteVector( ::H5::Group &g, const std::string &DSName, const std::vector<T> &v );
  // Read a matrix from a dataset
  template<typename T>
  void ReadMatrix( ::H5::Group &g, const std::string &DSName, Common::Matrix<T> &m );
  // Write a matrix to a dataset
  template<typename T>
  bool WriteMatrix( ::H5::Group &g, const std::string &DSName, const Common::Matrix<T> &m );

  // I would like these to be private members of a class full of static functions
  // BUT explicit specialisations must be at class scope (until C++ 17)
  template<typename T>
  void ReadStringsHelper( const T &a, const ::H5::StrType &aType, char * * MDString );

  template<typename T> inline typename T::value_type * GetDataHelper( T & t ) { return t.data(); }
  template<typename T> inline typename std::enable_if<!is_complex<T>::value, T *>::type
    GetDataHelper( Common::Vector<T> & m ) { return m.data; }
  template<typename T> inline typename std::enable_if< is_complex<T>::value, T *>::type
    GetDataHelper( Common::Vector<T> & m ) { return reinterpret_cast<T *>( m.data ); }
  // template<typename T> inline T * GetDataHelper( Common::Matrix<T> & m ) { return m.data; }
  std::string GetErrorClearStack( const ::H5::Exception &e );
};

template<typename T, typename AttOrDS> typename std::enable_if<H5::is_AttOrDS<AttOrDS>::value>::type
H5::ReadVectorHelper( T &v, AttOrDS &ds, const std::string &DSName )
{
  bool bError{ true };
  ::H5::DataSpace dsp = ds.getSpace();
  const int nDims{ dsp.getSimpleExtentNdims() };
  if( nDims == 1 )
  {
    hsize_t Dim;
    dsp.getSimpleExtentDims( &Dim );
    if( Dim > 0 && Dim <= std::numeric_limits<std::size_t>::max() )
    {
      v.resize( static_cast<std::size_t>( Dim ) );
      ReadAttOrDS( ds, GetDataHelper( v ), H5::Equiv<typename T::value_type>::Type );
      bError = false;
    }
  }
  if( bError )
    throw std::runtime_error( "Dimensions (" + std::to_string( nDims ) + ") invalid reading " + DSName );
}

template<typename T>
void H5::ReadVector( ::H5::Group &g, const std::string &DSName, T &v )
{
  ::H5::DataSet ds = g.openDataSet( DSName );
  ReadVectorHelper( v, ds, DSName );
}

// Read a vector of matrices from a dataset
template<typename T>
void H5::ReadVector( ::H5::Group &g, const std::string &DSName, std::vector<Common::Matrix<T>> &v )
{
  ::H5::DataSet ds = g.openDataSet( DSName );
  ::H5::DataSpace dsp = ds.getSpace();
  bool bError{ true };
  const int nDims{ dsp.getSimpleExtentNdims() };
  if( nDims == 3 )
  {
    hsize_t hDims[3];
    dsp.getSimpleExtentDims( hDims );
    std::size_t Dims[3];
    int i;
    for( i = 0; i < 3 && hDims[i] > 0 && hDims[i] <= std::numeric_limits<std::size_t>::max(); ++i )
      Dims[i] = static_cast<std::size_t>( hDims[i] );
    const std::size_t MatrixSize{ Dims[1] * Dims[2] };
    if( i == 3 && MatrixSize / Dims[1] == Dims[2] )
    {
      ::H5::DataSpace dspMemory( 2, &hDims[1] );
      hDims[0] = 1; // This becomes the count
      hsize_t hOffsets[3] = { 0, 0, 0 };
      v.reserve( Dims[0] );
      for( std::size_t i = 0; i < Dims[0]; ++i )
      {
        v.emplace_back( Dims[1], Dims[2] );
        hOffsets[0] = i;
        dsp.selectHyperslab( H5S_SELECT_SET, hDims, hOffsets );
        ds.read( v.back().data, H5::Equiv<T>::Type, dspMemory, dsp );
      }
      bError = false;
    }
  }
  if( bError )
    throw std::runtime_error( "Dimensions (" + std::to_string( nDims ) + ") invalid reading " + DSName );
}

// Write a vector to a dataset
template<typename T>
bool H5::WriteVector( ::H5::Group &g, const std::string &DSName, const Common::Vector<T> &v )
{
  const bool bHasData{ v.size != 0 };
  if( bHasData )
  {
    // Make sure the vector's memory is contiguous (copy it if not)
    Common::Vector<T> vContiguous;
    const Common::Vector<T> * pContiguous = &v;
    if( v.stride != 1 )
    {
      vContiguous = v;
      pContiguous = &vContiguous;
    }
    const hsize_t Dim{ pContiguous->size };
    ::H5::DataSpace dsp( 1, &Dim );
    ::H5::DataSet ds = g.createDataSet( DSName, Equiv<T>::Type, dsp );
    ds.write( pContiguous->data, Equiv<T>::Type );
    ds.close();
    dsp.close();
  }
  return bHasData;
}

template<typename T>
bool H5::WriteVector( ::H5::Group &g, const std::string &DSName, const std::vector<T> &v )
{
  const bool bHasData{ v.size != 0 };
  if( bHasData )
  {
    // Make sure the vector's memory is contiguous (copy it if not)
    const hsize_t Dim{ v.size() };
    ::H5::DataSpace dsp( 1, &Dim );
    ::H5::DataSet ds = g.createDataSet( DSName, Equiv<T>::Type, dsp );
    ds.write( v.data(), Equiv<T>::Type );
    ds.close();
    dsp.close();
  }
  return bHasData;
}

// Read a matrix from a dataset
template<typename T>
void H5::ReadMatrix( ::H5::Group &g, const std::string &DSName, Common::Matrix<T> &m )
{
  ::H5::DataSet ds = g.openDataSet( DSName );
  ::H5::DataSpace dsp = ds.getSpace();
  bool bError{ true };
  if( dsp.getSimpleExtentNdims() == 2 )
  {
    hsize_t Dims[2];
    dsp.getSimpleExtentDims( Dims );
    if( Dims[0] > 0 && Dims[0] < std::numeric_limits<int>::max()
     && Dims[1] > 0 && Dims[1] < std::numeric_limits<int>::max() )
    {
      m.resize( static_cast<int>( Dims[0] ), static_cast<int>( Dims[1] ) );
      ds.read( m.data, H5::Equiv<T>::Type );
      bError = false;
    }
  }
  if( bError )
    throw std::runtime_error( "Hdf5 object dimensions invalid: " + DSName );
}

// Write a matrix to a dataset
template<typename T>
bool H5::WriteMatrix( ::H5::Group &g, const std::string &DSName, const Common::Matrix<T> &m )
{
  const bool bHasData{ m.size1 && m.size2 };
  if( bHasData )
  {
    // Make sure the vector's memory is contiguous (copy it if not)
    Common::Matrix<T> mContiguous;
    const Common::Matrix<T> * pContiguous = &m;
    if( m.tda != m.size2 )
    {
      mContiguous = m;
      pContiguous = &mContiguous;
    }
    const hsize_t Dims[2]{ pContiguous->size1, pContiguous->size2 };
    ::H5::DataSpace dsp( 2, Dims );
    ::H5::DataSet ds = g.createDataSet( DSName, Equiv<T>::Type, dsp );
    ds.write( pContiguous->data, Equiv<T>::Type );
    ds.close();
    dsp.close();
  }
  return bHasData;
}

MLU_HDF5_hpp_end
#endif
