/*****************************************************************************
 
 Basic math support

 Source file: Math.cpp

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

#include "Math.hpp"
#include <gsl/gsl_cdf.h>

MLU_Math_hpp

/*****************************************************************************
 
 Chi^2 distribution

 ****************************************************************************/

template <typename T> T qValueChiSq( T ChiSquared, unsigned int dof )
{
  return gsl_cdf_chisq_Q( ChiSquared, dof );
}

template float  qValueChiSq( float  ChiSquared, unsigned int dof );
template double qValueChiSq( double ChiSquared, unsigned int dof );
template long double qValueChiSq( long double ChiSquared, unsigned int dof );

/*****************************************************************************
 
 Hotelling's T^2 distribution

 ****************************************************************************/

HotellingDist::HotellingDist( unsigned int p_, unsigned int m_ )
: p{p_}, m{m_}, Nu{m - p + 1}, Factor{ Nu / ( static_cast<double>( p ) * m ) }
{
  if( !Usable( p, m ) )
    throw std::runtime_error( "T^2 distribution m (degrees of freedom of covariance matrix) "
                             + std::to_string( m ) + " < p (fit degrees of freedom) " + std::to_string( p ) );
}

double HotellingDist::operator()( double TestStatistic ) const
{
  const double ModifiedStatistic{ TestStatistic * Factor };
  const double FDistQ{ gsl_cdf_fdist_Q( ModifiedStatistic, p, Nu ) };
#ifdef _DEBUG
  const double FDistP{ gsl_cdf_fdist_P( ModifiedStatistic, p, Nu ) };
  std::cout << "t^2=" << TestStatistic << ", Factor=" << Factor
            << ", t^2 * Factor=" << ModifiedStatistic << ", p [dof]=" << p
            << ", m [SampleSize-1]=" << m << ", Nu [m - p + 1]=" << Nu << Common::NewLine
            << "gsl_cdf_fdist_Q( t^2 * Factor, p, Nu )=" << FDistQ << Common::NewLine
            << "gsl_cdf_fdist_P( t^2 * Factor, p, Nu )=" << FDistP << Common::NewLine
            << "gsl_cdf_fdist_Q + gsl_cdf_fdist_P = " << ( FDistQ + FDistP ) << Common::NewLine;
#endif
  return FDistQ;
}

template <typename T> T HotellingDist::qValue( T TestStatistic, unsigned int p, unsigned int m )
{
  HotellingDist TSq( p, m );
  return TSq( TestStatistic );
}

template float  HotellingDist::qValue( float  TestStatistic, unsigned int p, unsigned int m );
template double HotellingDist::qValue( double TestStatistic, unsigned int p, unsigned int m );
template long double HotellingDist::qValue( long double TestStatistic, unsigned int p, unsigned int m );

/*****************************************************************************
 
 A value with error bars

 ****************************************************************************/

template <typename T> void ValWithEr<T>::Header( const std::string &FieldName, std::ostream &os, const std::string &Sep )
{
  os << FieldName << "_min" << Sep << FieldName << "_low" << Sep << FieldName << Sep
     << FieldName << "_high" << Sep << FieldName << "_max" << Sep << FieldName << "_check";
}

template <typename T>
ValWithEr<T>::ValWithEr( T min_, T low_, T central_, T high_, T max_, ValWithEr<T>::Scalar check_ )
: Min{min_}, Low{low_}, Central{central_}, High{high_}, Max{max_}, Check{check_} {}

template <typename T> ValWithEr<T>& ValWithEr<T>::operator=( const T Scalar )
{
  Min     = Scalar;
  Low     = Scalar;
  Central = Scalar;
  High    = Scalar;
  Max     = Scalar;
  Check   = 1;
  return *this;
}

template <typename T> ValWithEr<T>& ValWithEr<T>::operator=( const ValWithEr<T> &Other )
{
  Min     = Other.Min;
  Low     = Other.Low;
  Central = Other.Central;
  High    = Other.High;
  Max     = Other.Max;
  Check   = Other.Check;
  return *this;
}

template <typename T> template <typename U> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>&>::type
ValWithEr<T>::operator *=( const T Scalar )
{
  Min     *= Scalar;
  Low     *= Scalar;
  Central *= Scalar;
  High    *= Scalar;
  Max     *= Scalar;
  if( Scalar < 0 )
  {
    std::swap( Min, Max );
    std::swap( Low, High );
  }
  return *this;
}

template <typename T> template <typename U> typename std::enable_if<is_complex<U>::value, ValWithEr<T>&>::type
ValWithEr<T>::operator *=( const T Scalar )
{
  Min     *= Scalar;
  Low     *= Scalar;
  Central *= Scalar;
  High    *= Scalar;
  Max     *= Scalar;
  return *this;
}

template <typename T> ValWithEr<T> ValWithEr<T>::operator * ( const T Scalar ) const
{
  ValWithEr<T> Product( *this );
  Product *= Scalar;
  return Product;
}

template <typename T> template <typename U> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>>::type
ValWithEr<T>::qValueChiSq( unsigned int dof ) const
{
  // There is a deliberate inversion here: high statistic => low q-value
  return ValWithEr<T>( gsl_cdf_chisq_Q( Max,     dof ),
                       gsl_cdf_chisq_Q( High,    dof ),
                       gsl_cdf_chisq_Q( Central, dof ),
                       gsl_cdf_chisq_Q( Low,     dof ),
                       gsl_cdf_chisq_Q( Min,     dof ),
                       Check );
}

template <typename T> template <typename U> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>>::type
ValWithEr<T>::qValueHotelling( unsigned int p, unsigned int m ) const
{
  HotellingDist TSq( p, m );
  // There is a deliberate inversion here: high statistic => low q-value
  return ValWithEr<T>( TSq(Max), TSq(High), TSq(Central), TSq(Low), TSq(Min), Check );
}

// Sort the list of values, then extract the lower and upper 68th percentile error bars
template <typename T> template <typename U> typename std::enable_if<!is_complex<U>::value>::type
ValWithEr<T>::Get( T Central_, std::vector<T> &Data, std::size_t Count )
{
  assert( Data.size() >= Count && "ValWithErr<T>::Get() data too small" );
  Central = Central_;
  if( Count == 0 )
  {
    if( Data.size() == 0 )
    {
      // I wasn't trying to take an estimate over many samples
      Min = Central;
      Max = Central;
      Low = Central;
      High  = Central;
      Check = 1;
    }
    else
    {
      // Tried to take an estimate over many samples, but all values NaN
      Min = NaN;
      Max = NaN;
      Low = NaN;
      High  = NaN;
      Check = 0;
    }
    return;
  }
  const typename std::vector<T>::iterator itStart{ Data.begin() };
  const typename std::vector<T>::iterator itEnd  { itStart + Count };
  std::sort( itStart, itEnd );
  const std::size_t Index{ OneSigmaIndex( Count ) };
  Min  = Data[0];
  Max  = Data[Count - 1];
  Low  = Data[Index];
  High = Data[Count - 1 - Index];
  Check = static_cast<T>( static_cast<double>( Count ) / Data.size() );
}

// Sort the list of values, then extract the lower and upper 68th percentile error bars
// Sort real and imaginary parts separately
template <typename T> template <typename U> typename std::enable_if<is_complex<U>::value>::type
ValWithEr<T>::Get( T Central_, const std::vector<T> &Data, std::size_t Count )
{
  assert( Data.size() >= Count && "ValWithErr<T>::Get() data too small" );
  Central = Central_;
  if( Count == 0 )
  {
    if( Data.size() == 0 )
    {
      // I wasn't trying to take an estimate over many samples
      Min = Central;
      Max = Central;
      Low = Central;
      High  = Central;
      Check = 1;
    }
    else
    {
      // Tried to take an estimate over many samples, but all values NaN
      Min = NaN;
      Max = NaN;
      Low = NaN;
      High  = NaN;
      Check = 0;
    }
    return;
  }
  using Scalar = typename T::value_type;
  const std::size_t Index{ OneSigmaIndex( Count ) };
  std::vector<Scalar> Buffer( Count );
  for( int reim = 0; reim < 2; ++reim )
  {
    for( std::size_t i = 0; i < Count; ++i )
      Buffer[i] = reim ? Data[i].imag() : Data[i].real();
    const typename std::vector<Scalar>::iterator itStart{ Buffer.begin() };
    std::sort( itStart, itStart + Count );
    if( reim )
    {
       Min.imag( Buffer[0] );
       Max.imag( Buffer[Count - 1] );
       Low.imag( Buffer[Index] );
      High.imag( Buffer[Count - 1 - Index] );
    }
    else
    {
      Min.real( Buffer[0] );
      Max.real( Buffer[Count - 1] );
      Low.real( Buffer[Index] );
     High.real( Buffer[Count - 1 - Index] );
    }
  }
  Check = static_cast<Scalar>( Count ) / Data.size();
}

template <typename T>
void ValWithEr<T>::Get( T dCentral, const VectorView<T> &Source, std::vector<T> &ScratchBuffer )
{
  if( ScratchBuffer.size() < Source.size() )
    ScratchBuffer.resize( Source.size() );
  std::size_t Count{ 0 };
  for( std::size_t i = 0; i < Source.size(); ++i )
  {
    if( IsFinite( Source[i] ) )
      ScratchBuffer[Count++] = Source[i];
  }
  Get( dCentral, ScratchBuffer, Count );
}

template class ValWithEr<float>;
template class ValWithEr<double>;
template class ValWithEr<long double>;
template class ValWithEr<std::complex<float>>;
template class ValWithEr<std::complex<double>>;

template typename std::enable_if<!is_complex<float>::value, ValWithEr<float>>::type
ValWithEr<float>::qValueChiSq( unsigned int dof ) const;
template typename std::enable_if<!is_complex<double>::value, ValWithEr<double>>::type
ValWithEr<double>::qValueChiSq( unsigned int dof ) const;
template typename std::enable_if<!is_complex<long double>::value, ValWithEr<long double>>::type
ValWithEr<long double>::qValueChiSq( unsigned int dof ) const;

template typename std::enable_if<!is_complex<float>::value, ValWithEr<float>>::type
ValWithEr<float>::qValueHotelling( unsigned int p, unsigned int m ) const;
template typename std::enable_if<!is_complex<double>::value, ValWithEr<double>>::type
ValWithEr<double>::qValueHotelling( unsigned int p, unsigned int m ) const;
template typename std::enable_if<!is_complex<long double>::value, ValWithEr<long double>>::type
ValWithEr<long double>::qValueHotelling( unsigned int p, unsigned int m ) const;

template typename std::enable_if<!is_complex<float>::value>::type
ValWithEr<float>::Get( float Central_, std::vector<float> &Data, std::size_t Count );
template typename std::enable_if<!is_complex<double>::value>::type
ValWithEr<double>::Get( double Central_, std::vector<double> &Data, std::size_t Count );
template typename std::enable_if<!is_complex<long double>::value>::type
ValWithEr<long double>::Get( long double Central_, std::vector<long double> &Data, std::size_t Count );

template <typename T>
std::ostream & operator<<( std::ostream &os, const ValWithEr<T> &v )
{
  static const char *Space{ " " };
  return os << v.Min << Space << v.Low << Space << v.Central << Space << v.High << Space << v.Max << Space << v.Check;
}

template std::ostream & operator<<( std::ostream &os, const ValWithEr<float> &v );
template std::ostream & operator<<( std::ostream &os, const ValWithEr<double> &v );
template std::ostream & operator<<( std::ostream &os, const ValWithEr<std::complex<float>> &v );
template std::ostream & operator<<( std::ostream &os, const ValWithEr<std::complex<double>> &v );

MLU_Math_hpp_end
