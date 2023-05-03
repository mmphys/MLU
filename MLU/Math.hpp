/*****************************************************************************
 
 Basic math support

 Source file: Math.hpp
 
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

#ifndef MLU_Math_hpp
#define MLU_Math_hpp namespace Common {
#define MLU_Math_hpp_end };

#include <MLU/GSLVecMat.hpp>

#include <algorithm>
#include <limits>
#include <type_traits>

MLU_Math_hpp

// https://en.wikipedia.org/wiki/68–95–99.7_rule
constexpr double OneSigma{ 0.682689492137086 };
// How far from in each end of a sorted bootstrap ensemble we should take 1 sigma errors
constexpr double OneSigmaFromEnd{ ( 1. - OneSigma ) * 0.5 };
// Get the index into a sorted bootstrap ensemble to take 1 sigma errors
inline std::size_t OneSigmaIndex( std::size_t Count )
{
  return static_cast<std::size_t>( OneSigmaFromEnd * Count + 0.5 );
}

/**
 q-value according to chi^2 distribution, i.e. integral of chi^2( dof ) from ChiSquared to infinity
 
 - Parameter ChiSquared: This is the test statistic (dot product of scaled error vector with itself)
 - Parameter dof: Number of degrees of freedom (num data points - num model parameters)

 This is the probability of obtaining a worse test statistic that is still explained by the model.
 This is what people normally mean when they say "p-value".

 - Warning: When covariance comes from the data, this is not the correct distribution.
 Use the Hotelling distribution instead (below)
 */
template <typename T> T qValueChiSq( T ChiSquared, unsigned int dof );

/**
 q-value according to Hotelling's T^2 distribution, i.e. integral of Hotelling( p, m ) from TestStatistic to infinity
 
 - Parameter p: Covariance matrix (i.e. fit) degrees of freedom: num data points - num model parameters
 - Parameter m: Statistically independent samples used to build ovariance matrix: usually number of configs - 1

 When the covariance matrix comes from the data (rather than being known a priori),
 Hotelling found that the test statistic (t^2) has an T^2-distribution with p and m degrees of freedom.
 
 See: My year 4 notes;
 [Chris Kelly's Hotelling Study](https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/ckelly/Gparity/hotelling_v10.pdf);
 [Wikipedia T^2 distribution](https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution)
 
 - Warning: This is the correct distribution when covariance comes from the data (not chi^2)
 - Throws: std::runtime_error when p <= m. In order to avoid this, check `Usable()` returns `true` before calling `qValue()`
 */
struct HotellingDist
{
  const unsigned int p; // Dimension of covariance matrix (degrees of freedom of fit)
  const unsigned int m; // Degrees of freedom of covariance matrix (num samples to build covar - 1)
  const unsigned int Nu;
  const double Factor;
  HotellingDist( unsigned int p, unsigned int m );
  double operator()( double TestStatistic ) const;
  inline static bool Usable( unsigned int p, unsigned int m ) { return p <= m; }
  template <typename T> T static qValue( T TestStatistic, unsigned int p, unsigned int m );
};

// This is prior version. Not used except for reading in old files
template <typename T> struct ValWithErOldV1
{
  T Central;
  T Low;
  T High;
  T Check;
};

template <typename T> struct ValSigFig
{
  using value_type = T;
  T Value;
  unsigned char SigFig;
  int Exponent;
protected:
  int GetExponent() const;
public:
  ValSigFig( T Value, unsigned char SigFig );
  ValSigFig() : ValSigFig( 0, 2u ) {}
  ValSigFig( T Value, const ValSigFig<T> &Error );
  bool Zero() const { return Value == 0.; }
  bool Negative() const { return std::signbit( Value ); }
  int Class() const { return std::fpclassify( Value ); }
  bool ScientificHigh() const { return Exponent > std::max( 4, SigFig + 2 ); }
  bool ScientificLow() const { return Exponent < -4; }
  bool Scientific() const { return ScientificHigh() || ScientificLow(); }
  T GetUnit() const { return std::pow( 10., Exponent - ( SigFig - 1 ) ); }
  void AdjustExp( int Adjust );
  void Truncate();
  operator std::string() const;
};

/**
 Value with error
 
 The size of the Data buffer vector is the maximum number of samples.
 NaNs will be discarded and the proportion of non-NaN values recorded in Count as a fraction [0,1]
 
 - Warning: Non-`const std::vector<T> &Data` versions of `Get()` sort the data vector in-place
 */
template <typename T = double> struct ValWithEr
{
  using value_type = T;
  using Scalar = typename is_complex<T>::Scalar;
  using ValSigFigT = ValSigFig<T>;

  T Min;
  T Low;
  T Central;
  T High;
  T Max;
  Scalar Check;

  static void Header( const std::string &FieldName, std::ostream &os, const std::string &Sep = " " );

  ValWithEr() = default;
  template <typename U=T>
  ValWithEr( T dCentral, std::vector<T> &Data, std::size_t Count ) { Get( dCentral, Data, Count ); }
  ValWithEr( T Min, T Low, T Central, T High, T Max, Scalar Check = 1 );

  ValWithEr<T>& operator = ( const T Scalar );
  ValWithEr<T>& operator = ( const ValWithEr<T> &Other );
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>&>::type
  operator *=( const T Scalar );
  template <typename U=T> typename std::enable_if< is_complex<U>::value, ValWithEr<T>&>::type
  operator *=( const T Scalar );
  ValWithEr<T>  operator * ( const T Scalar ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>>::type
  qValueChiSq( unsigned int dof ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, ValWithEr<T>>::type
  qValueHotelling( unsigned int p, unsigned int m ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value>::type
  Get( T dCentral, std::vector<T> &Data, std::size_t Count );
  template <typename U=T> typename std::enable_if< is_complex<U>::value>::type
  Get( T dCentral, const std::vector<T> &Data, std::size_t Count );
  void Get( T dCentral, const VectorView<T> &Source, std::vector<T> &ScratchBuffer = std::vector<T>() );
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, std::string>::type
  to_string( unsigned char SigFigError ) const;
  template <typename U=T> typename std::enable_if<!is_complex<U>::value, std::string>::type
  to_string() const { return to_string( 2u ); }
};

template <typename T>
std::ostream & operator<<( std::ostream &os, const ValWithEr<T> &v );

MLU_Math_hpp_end
#endif
