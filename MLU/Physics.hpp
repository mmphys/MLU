/*****************************************************************************
 
 Physics primitives

 Source file: Physics.hpp
 
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

#ifndef MLU_Physics_hpp
#define MLU_Physics_hpp namespace Common {
#define MLU_Physics_hpp_end };

#include <array>
#include <ostream>
#include <regex>
#include <string>

MLU_Physics_hpp

extern const std::string Underscore;

namespace Gamma
{
  enum class Algebra
  {
    Spatial = -2,        // Gamma_spatial, i.e. X Y Z
    Unknown = -1,
    MinusGamma5,        // 0
    Gamma5,             // 1
    MinusGammaT,        // 2
    GammaT,             // 3
    MinusGammaTGamma5,  // 4
    GammaTGamma5,       // 5
    MinusGammaX,        // 6
    GammaX,             // 7
    MinusGammaXGamma5,  // 8
    GammaXGamma5,       // 9
    MinusGammaY,        // 10
    GammaY,             // 11
    MinusGammaYGamma5,  // 12
    GammaYGamma5,       // 13
    MinusGammaZ,        // 14
    GammaZ,             // 15
    MinusGammaZGamma5,  // 16
    GammaZGamma5,       // 17
    MinusIdentity,      // 18
    Identity,           // 19
    MinusSigmaXT,       // 20
    SigmaXT,            // 21
    MinusSigmaXY,       // 22
    SigmaXY,            // 23
    MinusSigmaXZ,       // 24
    SigmaXZ,            // 25
    MinusSigmaYT,       // 26
    SigmaYT,            // 27
    MinusSigmaYZ,       // 28
    SigmaYZ,            // 29
    MinusSigmaZT,       // 30
    SigmaZT,            // 31
  };
  static constexpr unsigned int                nGamma = 32;
  extern const std::array<const std::string, nGamma> name;      // Long name, per Grid
  extern const std::array<const std::string, nGamma> nameShort; // My abbreviations
  std::string NameShort( Algebra alg, const char * pszPrefix = nullptr, const char * pszOpSuffix = nullptr );
  std::ostream& operator<<(std::ostream& os, const Algebra &a);
  std::istream& operator>>(std::istream& is, Algebra &a);
};

struct ConfigCount {
  int         Config;
  int         Count = 0;
  ConfigCount( int Config_, int Count_ ) : Config{ Config_ }, Count{ Count_ } {}
  ConfigCount() = default;
  ConfigCount( const ConfigCount & ) = default;
  ConfigCount( ConfigCount && ) = default;
  ConfigCount& operator=( const ConfigCount & ) = default;
  ConfigCount& operator=( ConfigCount && ) = default;
};

std::ostream & operator<<( std::ostream &os, const ConfigCount &cc );

/**
 Specific to my project
 
 The absolute weights are meaningless. Relative masses are important.
 */
inline int QuarkWeight( const std::string &Quark )
{
  if( Quark.empty() )
    throw std::runtime_error( "QuarkWeight() called for empty Quark name" );
  int Weight;
  switch( std::toupper( Quark[0] ) )
  {
    case 'H':
      Weight = 1;
      break;
    case 'S':
      Weight = 0;
      break;
    case 'L':
      Weight = -1;
      break;
    default:
      throw std::runtime_error( "Relative weight of " + Quark + " quark unknown" );
      break;
  }
  return Weight;
}

std::string MakeMesonName( const std::string &Quark, const std::string &Spec );

// Generic representation of momentum
struct Momentum
{
  static const std::string DefaultPrefix;
  static const std::string SinkPrefix;
  static const std::string SquaredSuffix;
  static const std::string DefaultPrefixSquared;
  static const std::string SinkPrefixSquared;
  static const std::string DefaultRegexName;
  int x;
  int y;
  int z;
  bool bp2;
  Momentum( int _x, int _y, int _z, bool _bp2 ) : x(_x), y(_y), z(_z), bp2{_bp2} {}
  Momentum( int _x, int _y, int _z ) : Momentum( _x, _y , _z, false ) {}
  Momentum( bool _bp2 ) : Momentum( 0, 0 , 0, _bp2 ) {}
  Momentum() : Momentum( 0, 0, 0, false ) {}
  Momentum( const Momentum &o ) : Momentum( o.x, o.y, o.z, o.bp2 ) {}
  Momentum( int p2 ) : x(0), y(0), z(0), bp2{true} { *this = p2; }

  inline int p2() const { return x * x + y * y + z * z; }

  int &operator[]( int i )
  {
    if( i == 0 )
      return x;
    else if( i == 1 )
      return y;
    else if( i == 2 )
      return z;
    throw std::domain_error( "Bad Momentum index " + std::to_string( i ) );
  }

  const int &operator[]( int i ) const
  {
    if( i == 0 )
      return x;
    else if( i == 1 )
      return y;
    else if( i == 2 )
      return z;
    throw std::domain_error( "Bad Momentum index " + std::to_string( i ) );
  }

  /// Set momentum to p^2 (throw error if not possible)
  Momentum& operator=( const int p2 );
  Momentum& operator=( const Momentum& rhs )
  {
    if( this != &rhs )
    {
      x = rhs.x;
      y = rhs.y;
      z = rhs.z;
      bp2 = rhs.bp2;
    }
    return *this;
  }
  inline explicit operator bool() const { return x!=0 || y!=0 || z!=0; }
  inline Momentum operator-() const { return Momentum( -x, -y, -z, bp2 ); }
  inline Momentum abs() const { return Momentum( x<0?-x:x, y<0?-y:y, z<0?-z:z, bp2 ); }
  inline bool operator==( const Momentum &rhs ) const
  {
    if( bp2 || rhs.bp2 )
      return p2() == rhs.p2();
    return x == rhs.x && y == rhs.y && z == rhs.z;
  }
  inline bool operator!=( const Momentum &rhs ) const { return ! ( *this == rhs ); }
  Momentum& operator+=(const Momentum& rhs)
  {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
  }
  Momentum& operator-=(const Momentum& rhs)
  {
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
    return *this;
  }
  // friends defined inside class body are inline and are hidden from non-ADL lookup
  friend Momentum operator+( Momentum lhs, const Momentum &rhs ) // passing lhs by value helps optimize chained a+b+c
  {
    lhs += rhs; // reuse compound assignment
    return lhs; // return the result by value (uses move constructor)
  }
  friend Momentum operator-( Momentum lhs, const Momentum &rhs ) // passing lhs by value helps optimize chained a-b-c
  {
    lhs -= rhs; // reuse compound assignment
    return lhs; // return the result by value (uses move constructor)
  }

  inline bool IsNeg() const { return x<0 || ( x==0 && ( y<0 || ( y==0 && z < 0 ))); }
  inline bool EqualsNeg(const Momentum &m) const { return x==-m.x || y==-m.y || z==-m.z; }
  /// Convert to p^2 string: _Name_p^2
  std::string p2_string  ( const std::string &separator, const std::string &Name=DefaultPrefix )const;
  /// Convert to x-y-z components
  std::string to_string3d( const std::string &separator, bool bNegative = false ) const;
  /// Convert to x-y-z-t components (assume t=0)
  std::string to_string4d( const std::string &separator, bool bNegative = false ) const;
  /// Convert momentum to original format, with name
  std::string FileString( const std::string &Name, const std::string &Separator = Underscore ) const;
  /// Make a regular expression suitable for parsing momenta 3d or squared
  static std::regex MakeRegex( bool bp2, const std::string &Name = DefaultRegexName );
  // Strip out momentum info from string if present
  bool Extract( std::string &s, const std::string &MomName, bool IgnoreSubsequentZeroNeg = false );
  bool Extract( std::string &s ) { return Extract( s, DefaultPrefix, true ); }
  void Replace( std::string &s, const std::string &MomName, bool bNegative = false ) const;
  /**
   Use Lattice Dispersion relation and N=L/a to boost am to aE(p)
   
   - Parameters:
      - am: input mass (or energy - see `bGetGroundFromExcited`) in lattice units
      - N: integer spatial extent of lattice (i.e. number of sites)
      - bEnablePHat: `true` **default** use p_hat in dispersion relation (per the definition), `false` just use p
      - bGetGroundFromExcited: `false` **default** get aE(p) given am, `true` get am given aE(p)
   
   PhDYear3Diary.pdf eq (22)
   */
  double LatticeDispersion( double am, unsigned int N, bool bEnablePHat = true,
                            bool bGetGroundFromExcited = false ) const;
};

extern const Momentum p0;

std::ostream& operator<<( std::ostream& os, const Momentum &p );
std::istream& operator>>( std::istream& is, Momentum &p );

inline bool operator<( const Momentum &l, const Momentum &r )
{
  if( l == r )
    return false;
  const int lp2{ l.p2() };
  const int rp2{ r.p2() };
  if( lp2 != rp2 )
    return lp2 < rp2;
  if( l.x != r.x )
    return l.x < r.x;
  if( l.y != r.y )
    return l.y < r.y;
  return l.z < r.z;
}

MLU_Physics_hpp_end
#endif
