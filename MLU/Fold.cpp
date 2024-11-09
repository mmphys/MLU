/**

 Mike's lattice QCD utilities
 
 Source file: Fold.cpp
 
 Copyright (C) 2019 - 2023
 
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

#include <MLUconfig.h>
#include "Fold.hpp"

BEGIN_MLU_NAMESPACE

template <typename T>
const std::string &Fold<T>::DefaultGroupName()
{
  return sFold;
}

template <typename T>
bool Fold<T>::bFolded()
{
  return true;
}

template <typename T>
void Fold<T>::SummaryComments( std::ostream &s, bool bVerboseSummary, bool bShowColumnNames ) const
{
  Base::SummaryComments( s, bVerboseSummary, bShowColumnNames );
  if( NtUnfolded() ) s << "# NtUnfolded: " << NtUnfolded() << NewLine;
  if( reality != Reality::Unknown ) s << "# Reality: " << reality << NewLine;
  if( parity != Parity::Unknown ) s << "# Parity: " << parity << NewLine;
  if( sign != Sign::Unknown ) s << "# Sign: " << sign << NewLine;
  if( t0Negated ) s << "# timeslice 0 negated: true" << NewLine;
  if( Conjugated ) s << "# conjugate operator average: true" << NewLine;
}

template <typename T>
void Fold<T>::ReadAttributes( ::H5::Group &g )
{
  Base::ReadAttributes( g );
  ::H5::Attribute a;
  NtUnfolded_ = 0;
  try
  {
    a = g.openAttribute(sNtUnfolded);
    a.read( ::H5::PredType::NATIVE_INT, &NtUnfolded_ );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  reality = Reality::Unknown;
  try
  {
    a = g.openAttribute(sReality);
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    if( EqualIgnoreCase( s, Reality_TextEquiv_Real ) )
      reality = Reality::Real;
    else if( EqualIgnoreCase( s, Reality_TextEquiv_Imaginary ) )
      reality = Reality::Imag;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  parity = Parity::Unknown;
  try
  {
    a = g.openAttribute(sParity);
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    if( EqualIgnoreCase( s, Parity_TextEquiv_Even ) )
      parity = Parity::Even;
    else if( EqualIgnoreCase( s, Parity_TextEquiv_Odd ) )
      parity = Parity::Odd;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  sign = Sign::Unknown;
  try
  {
    a = g.openAttribute(sSign);
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    if( EqualIgnoreCase( s, Sign_TextEquiv_Positive ) )
      sign = Sign::Positive;
    else if( EqualIgnoreCase( s, Sign_TextEquiv_Negative ) )
      sign = Sign::Negative;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  t0Negated = false;
  std::int8_t i8;
  try
  {
    a = g.openAttribute(st0Negated);
    a.read( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    if( i8 )
      t0Negated = true;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  Conjugated = false;
  try
  {
    a = g.openAttribute(sConjugated);
    a.read( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    if( i8 )
      Conjugated = true;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  BootstrapList.clear();
  try
  {
    a = g.openAttribute(sBootstrapList);
    BootstrapList = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
}

template <typename T>
int Fold<T>::WriteAttributes( ::H5::Group &g ) const
{
  int iReturn = Sample<T>::WriteAttributes( g );
  const hsize_t OneDimension{ 1 };
  ::H5::DataSpace ds1( 1, &OneDimension );
  ::H5::Attribute a;
  if( NtUnfolded() )
  {
    a = g.createAttribute( sNtUnfolded, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &NtUnfolded_ );
    a.close();
    iReturn++;
  }
  if( reality == Reality::Real || reality == Reality::Imag )
  {
    const std::string &s{reality==Reality::Real?Reality_TextEquiv_Real:Reality_TextEquiv_Imaginary};
    a = g.createAttribute( sReality, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, s );
    a.close();
    iReturn++;
  }
  if( parity == Parity::Even || parity == Parity::Odd )
  {
    const std::string &s{ parity == Parity::Even ? Parity_TextEquiv_Even : Parity_TextEquiv_Odd };
    a = g.createAttribute( sParity, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, s );
    a.close();
    iReturn++;
  }
  if( sign == Sign::Positive || sign == Sign::Negative )
  {
    const std::string &s{ sign == Sign::Positive ? Sign_TextEquiv_Positive : Sign_TextEquiv_Negative};
    a = g.createAttribute( sSign, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, s );
    a.close();
    iReturn++;
  }
  const std::int8_t i8{ 1 };
  if( t0Negated )
  {
    a = g.createAttribute( st0Negated, ::H5::PredType::STD_U8LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    iReturn++;
  }
  if( Conjugated )
  {
    a = g.createAttribute( sConjugated, ::H5::PredType::STD_U8LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT8, &i8 );
    a.close();
    iReturn++;
  }
  if( BootstrapList.size() )
  {
    H5::WriteAttribute( g, sBootstrapList, BootstrapList );
    iReturn++;
  }
  return iReturn;
}

template class Fold<double>;
template class Fold<float>;
template class Fold<std::complex<double>>;
template class Fold<std::complex<float>>;

END_MLU_NAMESPACE
