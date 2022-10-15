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

const Momentum p0(0,0,0);
const std::string Momentum::DefaultPrefix{ "p" };
const std::string Momentum::SinkPrefix{ DefaultPrefix + "s" };
const std::string Momentum::SquaredSuffix{ "2" };

std::string Momentum::p2_string( const std::string &separator, const std::string Prefix ) const
{
  std::string s{ separator };
  s.append( Prefix );
  s.append( SquaredSuffix );
  s.append( separator );
  s.append( std::to_string( p2() ) );
  return s;
}

std::string Momentum::to_string( const std::string &separator, bool bNegative ) const
{
  return std::to_string( bNegative ? -x : x ) + separator
  + std::to_string( bNegative ? -y : y ) + separator
  + std::to_string( bNegative ? -z : z );
}

std::string Momentum::to_string4d( const std::string &separator, bool bNegative ) const
{
  return to_string( separator, bNegative ) + separator + "0";
}

std::regex Momentum::MakeRegex( const std::string &MomName )
{
  std::string sPattern( Common::Underscore );
  sPattern.append( MomName );
  sPattern.append( "_(-?[0-9]+)_(-?[0-9]+)_(-?[0-9]+)" );
  return std::regex( sPattern );
}

// Strip out momentum info from string if present
bool Momentum::Extract( std::string &Prefix, const std::string &MomName, bool IgnoreSubsequentZeroNeg )
{
  bool bGotMomentum = false;
  std::smatch match;
  const std::regex Pattern{ MakeRegex( MomName ) };
  while( std::regex_search( Prefix, match, Pattern ) )
  {
    int px{ std::stoi( match[1] ) };
    int py{ std::stoi( match[2] ) };
    int pz{ std::stoi( match[3] ) };
    if( !bGotMomentum )
    {
      bGotMomentum = true;
      x = px;
      y = py;
      z = pz;
    }
    else if( x != px || y != py || z != pz )
    {
      if( IgnoreSubsequentZeroNeg && !(*this) ) // i.e. previous momentum was zero
      {
        x = px;
        y = py;
        z = pz;
      }
      else if(IgnoreSubsequentZeroNeg && ((px==0 && py==0 && pz==0) || (px==-x && py==-y && pz==-z)))
        ;
      else
      {
        static const std::string Sep{ "," };
        static const std::string m1{ to_string( Sep ) };
        x = px;
        y = py;
        z = pz;
        throw std::runtime_error( "Multiple momenta: " + m1 + " != " + to_string( Sep ) );
      }
    }
    const std::string sSuffix{ match.suffix() };
    Prefix = match.prefix();
    Prefix.append( sSuffix );
  }
  return bGotMomentum;
}

void Momentum::Replace( std::string &s, const std::string &MomName, bool bNegative ) const
{
  std::string Replace( Common::Underscore );
  Replace.append( MomName );
  Replace.append( Common::Underscore );
  Replace.append( to_string( Common::Underscore, bNegative ) );
  std::smatch match;
  const std::regex Pattern{ MakeRegex( MomName ) };
  bool bFound{ false };
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

// Use Lattice Dispersion relation and N=L/a to boost am to aE(p)
double Momentum::LatticeDispersion( double am, unsigned int N ) const
{
  if( ! ( *this ) )
    return am;
  double Val = std::sinh( 0.5 * am );
  Val *= Val;
  const double NInv{ 1. / N };
  double w;
  if( x )
  {
    w = std::sin( M_PI * x * NInv );
    Val += w * w;
  }
  if( y )
  {
    w = std::sin( M_PI * y * NInv );
    Val += w * w;
  }
  if( z )
  {
    w = std::sin( M_PI * z * NInv );
    Val += w * w;
  }
  return 2 * std::asinh( std::sqrt( Val ) );
}

// Take a squared momentum and make an approximate 3d momentum from it
Momentum FileNameMomentum::FromSquared( const int p2Original )
{
  int p2{ p2Original };
  Momentum p;
  const bool bNegative{ p2 < 0 };
  if( bNegative )
    p2 = - p2;
  for( int i = 0; i < 3 && p2; ++i )
  {
    int root = 1;
    while( ( root + 1 ) * ( root + 1 ) <= p2 )
      root++;
    p2 -= root * root;
    switch( i )
    {
      case 0:
        p.x = root;
        break;
      case 1:
        p.y = root;
        break;
      case 2:
        p.z = root;
        break;
    }
  }
  if( p2 )
    throw std::runtime_error( std::to_string( p2Original ) + " can't be expressed as a three-momentum" );
  if( bNegative )
    p.x = -p.x;
  return p;
}

std::string FileNameMomentum::FileString( const std::string &separator ) const
{
  if( bp2 )
    return p2_string( separator, Name );
  std::string s{ separator };
  s.append( Name );
  s.append( separator );
  s.append( to_string( separator ) );
  return s;
}

std::ostream& operator<<( std::ostream& os, const Momentum &p )
{
  return os << p.to_string( Underscore );
}

std::istream& operator>>( std::istream& is, Momentum &p )
{
  char c = 0;
  if(!(is>>p.x && is.get(c) && c==Underscore[0] && is>>p.y && is.get(c) && c==Underscore[0] && is>>p.z))
    is.setstate( std::ios_base::failbit );
  return is;
}

MLU_Physics_hpp_end
