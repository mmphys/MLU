/*************************************************************************************
 
 HDF5 support

 Source file: HDF5.cpp
 
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
 *************************************************************************************/
/*  END LEGAL */

#include <MLUconfig.h>
#include "HDF5.hpp"
#include <iostream>

MLU_HDF5_hpp

// Make the same HDF5 complex type Grid uses
template<typename T> ::H5::CompType MakeComplex()
{
  ::H5::CompType myComplex( sizeof( std::complex<T> ) );
  myComplex.insertMember("re", 0 * sizeof(T), H5::Equiv<T>::Type);
  myComplex.insertMember("im", 1 * sizeof(T), H5::Equiv<T>::Type);
  return myComplex;
}

// Make an HDF5 type representing a value with an error
template<typename T> static ::H5::CompType MakeValWithEr()
{
  ::H5::CompType myType( sizeof( ValWithEr<T> ) );
  myType.insertMember("Min",     offsetof(ValWithEr<T>, Min    ), H5::Equiv<T>::Type);
  myType.insertMember("Low",     offsetof(ValWithEr<T>, Low    ), H5::Equiv<T>::Type);
  myType.insertMember("Central", offsetof(ValWithEr<T>, Central), H5::Equiv<T>::Type);
  myType.insertMember("High",    offsetof(ValWithEr<T>, High   ), H5::Equiv<T>::Type);
  myType.insertMember("Max",     offsetof(ValWithEr<T>, Max    ), H5::Equiv<T>::Type);
  myType.insertMember("Check",   offsetof(ValWithEr<T>, Check  ), H5::Equiv<T>::Type);
  return myType;
}

// Make an HDF5 type representing a value with an error
template<typename T> static ::H5::CompType MakeValWithErOldV1()
{
  ::H5::CompType myType( sizeof( ValWithErOldV1<T> ) );
  myType.insertMember("Central", offsetof(ValWithErOldV1<T>, Central), H5::Equiv<T>::Type);
  myType.insertMember("Low",     offsetof(ValWithErOldV1<T>, Low    ), H5::Equiv<T>::Type);
  myType.insertMember("High",    offsetof(ValWithErOldV1<T>, High   ), H5::Equiv<T>::Type);
  myType.insertMember("Check",   offsetof(ValWithErOldV1<T>, Check  ), H5::Equiv<T>::Type);
  return myType;
}

// Make an HDF5 type representing a config number and count
static ::H5::CompType MakeConfigCount()
{
  ::H5::CompType myType( sizeof( ConfigCount ) );
  myType.insertMember("Config", offsetof(ConfigCount, Config), ::H5::PredType::NATIVE_UINT32);
  myType.insertMember("Count",  offsetof(ConfigCount, Count ), ::H5::PredType::NATIVE_UINT32);
  return myType;
}

const ::H5::PredType& H5::Equiv<float>                      ::Type{ ::H5::PredType::NATIVE_FLOAT };
const ::H5::PredType& H5::Equiv<double>                     ::Type{ ::H5::PredType::NATIVE_DOUBLE };
const ::H5::PredType& H5::Equiv<long double>                ::Type{ ::H5::PredType::NATIVE_LDOUBLE };
const ::H5::PredType& H5::Equiv<int>                        ::Type{ ::H5::PredType::NATIVE_INT };
const ::H5::PredType& H5::Equiv<std::size_t>                ::Type{  sizeof( std::size_t ) == 4
                                                                    ? ::H5::PredType::NATIVE_UINT32
                                                                    : ::H5::PredType::NATIVE_UINT64 };
const ::H5::StrType   H5::Equiv<std::string>                ::Type{ ::H5::PredType::C_S1, H5T_VARIABLE };
const ::H5::StrType&  H5::Equiv<char *>                     ::Type{   H5::Equiv<std::string>::Type };
#ifndef HAVE_FAST32_IS_SIZE_T
const ::H5::PredType& H5::Equiv<std::uint_fast32_t>         ::Type{ sizeof( std::uint_fast32_t ) == 4
                                                              ? ::H5::PredType::STD_U32LE
                                                              : ::H5::PredType::STD_U64LE };
#endif
const ::H5::CompType  H5::Equiv<std::complex<float>>        ::Type{ MakeComplex<float>() };
const ::H5::CompType  H5::Equiv<std::complex<double>>       ::Type{ MakeComplex<double>() };
const ::H5::CompType  H5::Equiv<std::complex<long double>>  ::Type{ MakeComplex<long double>() };
const ::H5::CompType  H5::Equiv<ValWithEr<float>>           ::Type{ MakeValWithEr<float>() };
const ::H5::CompType  H5::Equiv<ValWithEr<double>>          ::Type{ MakeValWithEr<double>() };
const ::H5::CompType  H5::Equiv<ValWithEr<long double>>     ::Type{ MakeValWithEr<long double>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<float>>      ::Type{ MakeValWithErOldV1<float>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<double>>     ::Type{ MakeValWithErOldV1<double>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<long double>>::Type{ MakeValWithErOldV1<long double>() };
const ::H5::CompType  H5::Equiv<ConfigCount>                ::Type{ MakeConfigCount() };

// Open the specified HDF5File and group
void H5::OpenFileGroup( ::H5::H5File &f, ::H5::Group &g,
                        const std::string &FileName, const char *PrintPrefix,
                        std::string * pGroupName, unsigned int flags )
{
  const bool bFindGroupName{ !pGroupName || pGroupName->empty() };
  std::string localGroupName;
  if( !pGroupName )
    pGroupName = &localGroupName;
  f.openFile( FileName, flags );
  g = f.openGroup( bFindGroupName ? std::string("/") : *pGroupName );
  if( bFindGroupName ) {
    *pGroupName = GetFirstGroupName( g );
    g = g.openGroup( *pGroupName );
  }
  if( PrintPrefix )
    std::cout << PrintPrefix << FileName << " (" << *pGroupName << ")\n";
}

// Get first groupname from specified group
std::string H5::GetFirstGroupName( ::H5::Group & g )
{
  hsize_t n = g.getNumObjs();
  for( hsize_t i = 0; i < n; ++i ) {
    H5G_obj_t t = g.getObjTypeByIdx( i );
    if( t == H5G_GROUP )
      return g.getObjnameByIdx( i );
  }
  return std::string();
}

bool H5::OpenOptional( ::H5::Attribute &a, ::H5::Group &g, const std::string Name )
{
  bool bOpen{ false };
  if( g.attrExists( Name ) )
  {
    H5E_auto2_t h5at;
    void      * f5at_p;
    ::H5::Exception::getAutoPrint( h5at, &f5at_p );
    ::H5::Exception::dontPrint();
    try
    {
      a = g.openAttribute( Name );
      bOpen = true;
    }
    catch( const ::H5::Exception &e )
    {
      ::H5::Exception::clearErrorStack();
    }
    catch(...)
    {
      ::H5::Exception::setAutoPrint( h5at, f5at_p );
      throw;
    }
    ::H5::Exception::setAutoPrint( h5at, f5at_p );
  }
  return bOpen;
}

bool H5::OpenOptional( ::H5::DataSet &ds, ::H5::Group &g, const std::string Name )
{
  bool bOpen{ false };
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint( h5at, &f5at_p );
  ::H5::Exception::dontPrint();
  try
  {
    ds = g.openDataSet( Name );
    bOpen = true;
  }
  catch( const ::H5::Exception &e )
  {
    ::H5::Exception::clearErrorStack();
  }
  catch(...)
  {
    ::H5::Exception::setAutoPrint( h5at, f5at_p );
    throw;
  }
  ::H5::Exception::setAutoPrint( h5at, f5at_p );
  return bOpen;
}

bool H5::OpenOptional( ::H5::Group &gNew, ::H5::Group &gOld, const std::string Name )
{
  bool bOpen{ false };
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint( h5at, &f5at_p );
  ::H5::Exception::dontPrint();
  try
  {
    gNew = gOld.openGroup( Name );
    bOpen = true;
  }
  catch( const ::H5::Exception &e )
  {
    ::H5::Exception::clearErrorStack();
  }
  catch(...)
  {
    ::H5::Exception::setAutoPrint( h5at, f5at_p );
    throw;
  }
  ::H5::Exception::setAutoPrint( h5at, f5at_p );
  return bOpen;
}

// Read the gamma algebra attribute string and make sure it's valid
Gamma::Algebra H5::ReadGammaAttribute( ::H5::Group &g, const char * pAttName )
{
  std::string sGamma;
  ::H5::Attribute a = g.openAttribute( pAttName );
  ::H5::StrType s = a.getStrType();
  a.read( s, sGamma );
  a.close();
  for( int idxGamma = 0; idxGamma < Gamma::nGamma; idxGamma++ )
    if( EqualIgnoreCase( sGamma, Gamma::name[idxGamma] ) )
      return static_cast<Gamma::Algebra>( idxGamma );
  throw ::H5::Exception( "Common::ReadGammaAttribute", "Invalid gamma algebra string" );
}

template<> void H5::ReadStringsHelper<::H5::Attribute>( const ::H5::Attribute &a, const ::H5::StrType &aType,
                                                        char * * MDString )
{
  a.read( aType, ( void * ) MDString );
}
template<> void H5::ReadStringsHelper<::H5::DataSet>( const ::H5::DataSet &ds, const ::H5::StrType &aType,
                                                      char * * MDString )
{
  ds.read( ( void * ) MDString, aType );
}

template<typename AttOrDS> std::vector<std::string> H5::ReadStrings( const AttOrDS &a )
{
  ::H5::DataSpace dsp = a.getSpace();
  const ::H5::StrType aType = a.getStrType();
  const int rank{ dsp.getSimpleExtentNdims() };
  if( rank != 1 )
  {
    std::size_t len = H5Iget_name( a.getId(), NULL, 0 );
    char buffer[len+1];
    H5Iget_name( a.getId(), buffer, len+1 );
    std::string n = buffer;
    throw std::runtime_error( n + " number of dimensions=" + std::to_string( rank ) + ", expecting 1" );
  }
  hsize_t NumOps;
  dsp.getSimpleExtentDims( &NumOps );
  std::vector<std::string> vs( NumOps );
  char * * MDString = new char *[NumOps];
  for( std::size_t i = 0; i < NumOps; i++ )
    MDString[i] = nullptr;
  try
  {
    ReadStringsHelper( a, aType, MDString );
  }
  catch(...)
  {
    for( std::size_t i = 0; i < NumOps; i++ )
      if( MDString[i] )
        delete [] MDString[i];
    delete [] MDString;
    throw;
  }
  for( std::size_t i = 0; i < NumOps; i++ )
  {
    vs[i] = MDString[i];
    delete [] MDString[i];
  }
  delete [] MDString;
  return vs;
}

template std::vector<std::string> H5::ReadStrings<::H5::Attribute>( const ::H5::Attribute &a  );
template std::vector<std::string> H5::ReadStrings<::H5::DataSet>  ( const ::H5::DataSet   &ds );

// Make a multi-dimensional string attribute
void H5::WriteAttribute( ::H5::Group &g, const std::string &AttName, const std::vector<std::string> &vs)
{
  std::unique_ptr<char *[]> RawArray( new char * [vs.size()] );
  const char * * const MDString{ const_cast<const char * *>( RawArray.get() ) };
  for( std::size_t i = 0; i < vs.size(); i++ )
    MDString[i] = vs[i].c_str();
  const hsize_t NDimension{ vs.size() };
  ::H5::DataSpace dsN( 1, &NDimension );
  ::H5::Attribute a{ g.createAttribute( AttName, Equiv<std::string>::Type, dsN ) };
  a.write( Equiv<std::string>::Type, MDString );
  a.close();
}

// Make a multi-dimensional string attribute
void H5::WriteStringData( ::H5::Group &g, const std::string &DSName, const std::vector<std::string> &vs)
{
  std::unique_ptr<char *[]> RawArray( new char * [vs.size()] );
  const char * * const MDString{ const_cast<const char * *>( RawArray.get() ) };
  for( std::size_t i = 0; i < vs.size(); i++ )
    MDString[i] = vs[i].c_str();
  const hsize_t NDimension{ vs.size() };
  ::H5::DataSpace dsN( 1, &NDimension );
  ::H5::DataSet ds{ g.createDataSet( DSName, Equiv<std::string>::Type, dsN ) };
  ds.write( MDString, Equiv<std::string>::Type );
  ds.close();
}

std::string H5::GetErrorClearStack( const ::H5::Exception &e )
{
  std::string s{ "HDF5 error in " + e.getFuncName() + ": " + e.getCDetailMsg() };
  ::H5::Exception::clearErrorStack();
  return s;
}

MLU_HDF5_hpp_end
