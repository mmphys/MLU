/*************************************************************************************
 
 Physics primitives

 Source file: Physics.cpp

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

#include "Common.hpp"

MLU_Physics_hpp

namespace Gamma
{
  const std::array<const std::string, nGamma> name{ // Long name, per Grid
    "MinusGamma5",
    "Gamma5",
    "MinusGammaT",
    "GammaT",
    "MinusGammaTGamma5",
    "GammaTGamma5",
    "MinusGammaX",
    "GammaX",
    "MinusGammaXGamma5",
    "GammaXGamma5",
    "MinusGammaY",
    "GammaY",
    "MinusGammaYGamma5",
    "GammaYGamma5",
    "MinusGammaZ",
    "GammaZ",
    "MinusGammaZGamma5",
    "GammaZGamma5",
    "MinusIdentity",
    "Identity",
    "MinusSigmaXT",
    "SigmaXT",
    "MinusSigmaXY",
    "SigmaXY",
    "MinusSigmaXZ",
    "SigmaXZ",
    "MinusSigmaYT",
    "SigmaYT",
    "MinusSigmaYZ",
    "SigmaYZ",
    "MinusSigmaZT",
    "SigmaZT",
  };
  const std::array<const std::string, nGamma> nameShort{ // My abbreviations
    "g-5",
    "g5",
    "g-T",
    "gT",
    "g-T5",
    "gT5",
    "g-X",
    "gX",
    "g-X5",
    "gX5",
    "g-Y",
    "gY",
    "g-Y5",
    "gY5",
    "g-Z",
    "gZ",
    "g-Z5",
    "gZ5",
    "g-I",
    "gI",
    "g-sXT",
    "gsXT",
    "g-sXY",
    "gsXY",
    "g-sXZ",
    "gsXZ",
    "g-sYT",
    "gsYT",
    "g-sYZ",
    "gsYZ",
    "g-sZT",
    "gsZT",
  };
  const std::string sGammaSpatial{ "gXYZ" };

  std::string NameShort( Algebra alg, const char * pszPrefix, const char * pszOpSuffix )
  {
    std::string sName;
    const bool bGotSuffix{ pszOpSuffix && *pszOpSuffix };
    if( alg != Algebra::Unknown || bGotSuffix )
    {
      if( pszPrefix && *pszPrefix )
        sName.append( pszPrefix );
      if( alg == Algebra::Spatial )
        sName.append( sGammaSpatial );
      else if( alg != Algebra::Unknown )
      {
        int idx = static_cast<int>( alg );
        sName.append( ( idx >= 0 && idx < nGamma ) ? nameShort[idx] : "?" );
      }
      if( bGotSuffix )
        sName.append( pszOpSuffix );
    }
    return sName;
  }

  std::ostream& operator<<(std::ostream& os, const Algebra &a)
  {
    return os << NameShort( a );
  }

  std::istream& operator>>(std::istream& is, Algebra &a)
  {
    std::string s;
    if( is >> s )
    {
      int i;
      for( i = 0; i < nGamma && !EqualIgnoreCase( s, nameShort[i] ); i++ )
        ;
      if( i < nGamma )
        a = static_cast<Algebra>( i );
      else if( EqualIgnoreCase( s, sGammaSpatial ) )
        a = Algebra::Spatial;
      else
        is.setstate( std::ios_base::failbit );
    }
    return is;
  }
};

std::ostream & operator<<( std::ostream &os, const ConfigCount &cc )
{
  return os << "{" << cc.Config << "," << cc.Count << "}";
}

// For my current project, I always have the heavier quark at the sink
std::string ConcatHeavyFirstNoExcept( const std::string &Quark1, const std::string &Quark2 )
{
  std::string s;
  if( QuarkWeightNoExcept( Quark1 ) < QuarkWeightNoExcept( Quark2 ) )
  {
    s.append( Quark2 );
    s.append( Underscore );
    s.append( Quark1 );
  }
  else
  {
    s.append( Quark1 );
    s.append( Underscore );
    s.append( Quark2 );
  }
  return s;
}

std::string MesonName( const std::string &q1, const std::string &q2 )
{
  std::string s;
  const bool bSwap{ QuarkWeight( q1 ) < QuarkWeight( q2 ) };
  const std::string &sHeavy{ bSwap ? q2 : q1 };
  const std::string &sLight{ bSwap ? q1 : q2 };
  const int Heavy{ std::toupper( sHeavy[0] ) };
  const int Light{ std::toupper( sLight[0] ) };
  switch( Heavy )
  {
    case 'H':
      if( Light == 'S' )
        s = "Ds";
      else if( Light == 'L' )
        s = "D";
      break;
    case 'S':
      if( Light == 'L' )
        s = "K";
      break;
    case 'L':
      if( Light == 'L' )
        s = "Pi";
      break;
  }
  return s;
}

const Momentum p0(0,0,0);
const std::string Momentum::DefaultPrefix{ "p" };
const std::string Momentum::SinkPrefix{ DefaultPrefix + "s" };
const std::string Momentum::SquaredSuffix{ "2" };
const std::string Momentum::DefaultPrefixSquared{ DefaultPrefix + SquaredSuffix };
const std::string Momentum::SinkPrefixSquared{ SinkPrefix + SquaredSuffix };
const std::string Momentum::DefaultRegexName{ "[pP][[:alnum:]]*" };

Momentum& Momentum::operator=( const int p2Original )
{
  int p2{ p2Original };
  int p[3];
  const bool bNegative{ p2 < 0 };
  if( bNegative )
  {
    p2 = - p2;
    throw std::domain_error( "Negative Momentum^2 " + std::to_string( p2Original ) );
  }
  for( int i = 0; i < 3; ++i )
    p[i] = 0;
  for( int i = 0; i < 3 && p2; ++i )
  {
    int root = 1;
    while( ( root + 1 ) * ( root + 1 ) <= p2 )
      ++root;
    p2 -= root * root;
    p[i] = root;
  }
  if( p2 )
    throw std::domain_error( std::to_string( p2Original )
                             + " can't be expressed as a three-momentum" );
  x = p[0];
  y = p[1];
  z = p[2];
  bp2 = true;
  if( bNegative )
    x = -x;
  return *this;
}

std::string Momentum::p2_string( const std::string &Separator, const std::string &Name ) const
{
  std::string s{ Separator };
  s.append( Name );
  if( !s.empty() )
  {
    s.append( SquaredSuffix );
    s.append( Separator.empty() ? Underscore : Separator );
  }
  s.append( std::to_string( p2() ) );
  return s;
}

std::string Momentum::to_string3d( const std::string &separator, bool bNegative ) const
{
  return std::to_string( bNegative ? -x : x ) + separator
  + std::to_string( bNegative ? -y : y ) + separator
  + std::to_string( bNegative ? -z : z );
}

std::string Momentum::to_string4d( const std::string &separator, bool bNegative ) const
{
  return to_string3d( separator, bNegative ) + separator + "0";
}

std::string Momentum::FileString( const std::string &Name, const std::string &Separator ) const
{
  if( bp2 )
    return p2_string( Separator, Name );
  std::string s{ Separator };
  s.append( Name );
  const std::string &Sep{ Separator.empty() ? Underscore : Separator };
  if( !s.empty() )
    s.append( Sep );
  s.append( to_string3d( Sep ) );
  return s;
}

std::regex Momentum::MakeRegex( bool bp2, const std::string &Name )
{
  static const char szDigits[] = "_(-?[0-9]+)";
  std::string sPattern( "_(" );
  sPattern.append( Name );
  sPattern.append( ")" );
  if( bp2 )
  {
    sPattern.append( Momentum::SquaredSuffix );
    sPattern.append( szDigits );
  }
  else
  {
    for( int i = 0; i < 3; ++i )
      sPattern.append( szDigits );
  }
  return std::regex( sPattern );
}

// Strip out momentum info from string if present
bool Momentum::Extract(std::string &Prefix, const std::string &MomName, bool IgnoreSubsequentZeroNeg)
{
  Momentum Other;
  bool bGotMomentum = false;
  std::smatch match;
  for( int Pass = 0; Pass < 2; ++Pass )
  {
    const std::regex Pattern{ MakeRegex( Pass, MomName ) };
    while( std::regex_search( Prefix, match, Pattern ) )
    {
      if( Pass )
        Other = std::stoi( match[2] );
      else
      {
        Other.bp2 = false;
        Other.x = std::stoi( match[2] );
        Other.y = std::stoi( match[3] );
        Other.z = std::stoi( match[4] );
      }
      if( !bGotMomentum )
      {
        bGotMomentum = true;
        *this = Other;
      }
      else if( *this != Other )
      {
        if( IgnoreSubsequentZeroNeg && !*this ) // i.e. previous momentum was zero
          *this = Other;
        else if( IgnoreSubsequentZeroNeg && ( !Other || EqualsNeg( Other ) ) )
          ;
        else
          throw std::runtime_error( "Multiple momenta " + MomName + Space + FileString( "", "" )
                                   + " != " + Other.FileString( "", "" ) );
      }
      const std::string sSuffix{ match.suffix() };
      Prefix = match.prefix();
      Prefix.append( sSuffix );
    }
  }
  return bGotMomentum;
}

void Momentum::Replace( std::string &s, const std::string &MomName, bool bNegative ) const
{
  const std::string Replace( FileString( MomName ) );
  std::smatch match;
  bool bFound{ false };
  for( int Pass = 0; Pass < 2; ++Pass )
  {
    const std::regex Pattern{ MakeRegex( Pass, MomName ) };
    std::string Search{ std::move( s ) };
    s.clear();
    while( std::regex_search( Search, match, Pattern ) )
    {
      s.append( match.prefix() );
      s.append( Replace );
      Search = match.suffix();
      bFound = true;
    }
    if( !bFound )
    {
      std::ostringstream ss;
      ss << "momentum " << MomName << " not found in " << s;
      throw std::runtime_error( ss.str() );
    }
    s.append( Search );
  }
}

// Use Lattice Dispersion relation and N=L/a to boost am to aE(p)
double Momentum::LatticeDispersion( double am, unsigned int N, bool bEnablePHat,
                                    bool bGetGroundFromExcited ) const
{
  if( ! ( *this ) )
    return am;
  double Val = std::sinh( 0.5 * am );
  Val *= Val;
  const double NInv{ 1. / N };
  for( int i = 0; i < 3; ++i )
  {
    if( (*this)[i] )
    {
      double w = M_PI * NInv * (*this)[i]; // w = a p / 2
      if( bEnablePHat )
        w = std::sin( w );                 // w = a p_hat / 2
      w *= w;
      if( bGetGroundFromExcited )
        w = -w;
      Val += w;
    }
  }
  return ( Val <= 0 ) ? 0 : 2 * std::asinh( std::sqrt( Val ) );
}

std::ostream& operator<<( std::ostream& os, const Momentum &p )
{
  return os << p.to_string3d( Underscore );
}

std::istream& operator>>( std::istream& is, Momentum &p )
{
  char c = 0;
  if(!(is>>p.x && is.get(c) && c==Underscore[0] && is>>p.y && is.get(c) && c==Underscore[0] && is>>p.z))
    is.setstate( std::ios_base::failbit );
  return is;
}

const std::array<std::string, 4> FormFactorString{ "f0", "fplus", "fpar", "fperp" };

const std::string &GetFormFactorString( FormFactor ff )
{
  using T = typename std::underlying_type<FormFactor>::type;
  T i{ static_cast<T>( ff ) };
  return i >= 0 && i < FormFactorString.size() ? FormFactorString[i] : sUnknown;
}

std::ostream& operator<<( std::ostream& os, FormFactor ff )
{
  return os << GetFormFactorString( ff );
}

std::istream& operator>>( std::istream& is, FormFactor &ff )
{
  std::string s;
  if( is >> s )
  {
    using T = typename std::underlying_type<FormFactor>::type;
    for( T i = 0; i < FormFactorString.size(); ++i )
    {
      if( EqualIgnoreCase( s, FormFactorString[i] ) )
      {
        ff = static_cast<FormFactor>( i );
        return is;
      }
    }
  }
  is.setstate( std::ios_base::failbit );
  return is;
}

MLU_Physics_hpp_end
