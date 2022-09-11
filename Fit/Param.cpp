/*************************************************************************************
 
 Manage model parameters

 Source file: Param.cpp
 
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

#include "Param.hpp"

static const std::array<std::string,3> ParamTypeHuman{ "All", "Variable", "Fixed" };

std::ostream &operator<<( std::ostream &os, const ParamBase::Type &type )
{
  std::size_t i{ static_cast<std::size_t>( type ) };
  if( i < ParamTypeHuman.size() )
    os << ParamTypeHuman[i];
  else
    os << Common::sUnknown << "(" << i << ")";
  return os;
}

std::istream &operator>>( std::istream &is, ParamBase::Type &Type )
{
  std::string s;
  if( is >> s )
  {
    const int Idx{ Common::IndexIgnoreCase( ParamTypeHuman, s ) };
    if( Idx != ParamTypeHuman.size() )
    {
      Type = static_cast<Param::Type>( Idx );
      return is;
    }
  }
  throw std::runtime_error( "Param::Type \"" + s + "\" unrecognised" );
}

void ParamBase::Validate( ParamBase::Type Type, const char * pErrorIfAll )
{
  if( Type == Param::Type::All && pErrorIfAll )
    throw std::runtime_error( pErrorIfAll );
  std::size_t i{ static_cast<std::size_t>( Type ) };
  if( i >= ParamTypeHuman.size() )
    throw std::runtime_error( "Param::Type \"" + std::to_string( i ) + "\" unrecognised" );
}

std::string Param::Key::FullName( std::size_t idx, std::size_t Size ) const
{
  std::ostringstream os;
  os << *this;
  if( Size > 1 )
    os << idx;
  return os.str();
}

bool Param::Key::operator==( const Key &rhs ) const
{
  return !Common::CompareIgnoreCase( Object, rhs.Object )
      && !Common::CompareIgnoreCase( Name, rhs.Name );
}

bool Param::Key::operator!=( const Key &rhs ) const
{
  return Common::CompareIgnoreCase( Object, rhs.Object )
      || Common::CompareIgnoreCase( Name, rhs.Name );
}

bool Param::Key::Less::operator()( const Key &lhs, const Key &rhs ) const
{
  int i = Common::CompareIgnoreCase( lhs.Object, rhs.Object );
  if( i == 0 )
    i = Common::CompareIgnoreCase( lhs.Name, rhs.Name );
  return i < 0;
}

std::ostream &operator<<( std::ostream &os, const Param::Key &key )
{
  if( !key.Object.empty() )
    os << key.Object << '-';
  return os << key.Name;
}

std::ostream &operator<<( std::ostream &os, const Params::iterator &it )
{
  const Param &p{ it->second };
  os << it->first << '[' << p.size;
  if( p.type != Param::Type::Variable )
    os << ',' << p.type;
  if( p.bMonotonic )
    os << 'M';
  return os << ']';
}

std::size_t Param::operator()( std::size_t Idx, Param::Type ListType ) const
{
  Param::Validate( ListType );
  if( Idx >= size )
  {
    std::ostringstream s;
    s << "Param::operator() Index " << Idx << " >= size " << size;
    throw std::runtime_error( s.str().c_str() );
  }
  if( ListType != type && ListType != Param::Type::All )
  {
    std::ostringstream s;
    s << "Param::operator() can't read " << type << " parameter from " << ListType << " list";
    throw std::runtime_error( s.str().c_str() );
  }
  return ( ListType == Param::Type::All ? OffsetAll : OffsetMyType ) + Idx;
}

template <typename T>
T &Param::operator()( std::vector<T> &v, std::size_t Idx, Param::Type ListType ) const
{
  std::size_t Offset{ (*this)( Idx, ListType ) };
  if( Offset >= v.size() )
  {
    std::ostringstream s;
    s << "Param::operator() Index " << Idx << " >= size " << v.size() << " on std::vector of " << ListType;
    throw std::runtime_error( s.str().c_str() );
  }
  return v[Offset];
}

/*Param::Flags operator~( Param::Flags f )
{
  return static_cast<Param::Flags>( ~ static_cast<Param::FlagsU>( f ) );
}

Param::Flags operator|( Param::Flags l, Param::Flags r )
{
  return static_cast<Param::Flags>( static_cast<Param::FlagsU>( l ) | static_cast<Param::FlagsU>( r ) );
}

Param::Flags operator&( Param::Flags l, Param::Flags r )
{
  return static_cast<Param::Flags>( static_cast<Param::FlagsU>( l ) & static_cast<Param::FlagsU>( r ) );
}

Param::Flags &operator|=( Param::Flags &l, const Param::Flags &r )
{
  l = l | r;
  return l;
}

Param::Flags &operator&=( Param::Flags &l, const Param::Flags &r )
{
  l = l & r;
  return l;
}*/

Params::iterator Params::Add( const Param::Key &key, std::size_t size, bool bMonotonic, Param::Type type_ )
{
  if( !size )
  {
    std::ostringstream ss;
    ss << "Param::Add key " << key << " has zero members";
    throw std::runtime_error( ss.str().c_str() );
  }
  iterator it{ find( key ) };
  if( it == end() )
  {
    // Doesn't exist - Add it
    try
    {
      it = emplace( std::make_pair( key, Param( size, bMonotonic, type_ ) ) ).first;
    }
    catch( std::exception &e )
    {
      std::ostringstream ss;
      ss << "Param::Add cannot add key " << key;
      throw std::runtime_error( ss.str().c_str() );
    }
  }
  else
  {
    // Exists - see whether it needs to be altered
    Param &p{ it->second };
    if( p.type != type_ )
    {
      // Type is changing.
      const bool bCantShrink{ type_ == Param::Type::Fixed && size < p.size };
      const bool bCantGrow{ type_ == Param::Type::Variable && size > p.size };
      if( bCantShrink || bCantGrow )
      {
        std::ostringstream s;
        s << "Param::Add key " << key << " can't " << ( bCantGrow ? "grow" : "shrink" )
          << " from " << p.type << "[" << p.size << ""
          << "] to " << type_ << "[" << size << Common::Space << "] parameters";
        throw std::runtime_error( s.str().c_str() );
      }
      // Can't become variable ... but can go the other way
      if( type_ != Param::Type::Variable )
      {
        p.type = type_;
        p.bMonotonic = bMonotonic; // Doesn't really do anything unless Variable
      }
    }
    // Has there been a request to change Monotonic?
    if( ( p.bMonotonic && !bMonotonic ) || ( !p.bMonotonic && bMonotonic ) )
    {
      // Silently ignore requests to remove monotonic from Variable parameters
      if( !( !bMonotonic && p.type == Param::Type::Variable ) )
        p.bMonotonic = bMonotonic;
    }
    // Allow growth
    if( p.size < size )
      p.size = size;
  }
  return it;
}

void Params::AssignOffsets()
{
  NumFixed = 0;
  NumVariable = 0;
  MaxLen = 0;
  for( value_type &it : *this )
  {
    const Param::Key &k{ it.first };
    Param &p{ it.second };
    // Validate the parameter type
    if( p.type == Param::Type::All )
    {
      std::ostringstream s;
      s << "Params::AssignOffsets() key " << k << " invalid type " << p.type;
      throw std::runtime_error( s.str().c_str() );
    }
    Param::Validate( p.type );
    // Assign positions of these parameters in tables
    p.OffsetAll = NumScalars( Param::Type::All );
    std::size_t &ThisNum{ p.type == Param::Type::Variable ? NumVariable : NumFixed };
    p.OffsetMyType = ThisNum;
    ThisNum += p.size;
    // To simplify printing, save the maximum length of this field name
    p.FieldLen = k.Len();
    if( p.size > 1 )
    {
      // When there's more than one parameter, it is followed by a decimal number
      std::ostringstream o;
      o << ( p.size - 1 );
      p.FieldLen += o.str().size();
    }
    if( MaxLen < p.FieldLen )
      MaxLen = p.FieldLen;
  }
}

void Params::Export( Vector &vType, const Vector &All, Param::Type type ) const
{
  if( All.size != NumScalars( Param::Type::All ) )
  {
    std::ostringstream es;
    es << "Params::Export() In[" << Param::Type::All << "] size=" << All.size
       << " but should be " << NumScalars( Param::Type::All );
    throw( es.str().c_str() );
  }
  if( vType.size != NumScalars( type ) )
  {
    std::ostringstream es;
    es << "Params::Export() Out[" << type << "] size=" << All.size
       << " but should be " << NumScalars( Param::Type::All );
    throw( es.str().c_str() );
  }
  for( const value_type &it : *this )
  {
    const Param &p{ it.second };
    if( p.type == type )
    {
      const std::size_t &Offset{ type == Param::Type::All ? p.OffsetAll : p.OffsetMyType };
      for( std::size_t i = 0; i < p.size; ++i )
      {
        if( i == 0 || !p.bMonotonic || type != Param::Type::Variable )
          vType[Offset + i] = All[p.OffsetAll + i];
        else
        {
          // Monotonically increasing set of values. Implemented as a sum of squares
          const scalar diff{ All[p.OffsetAll + i] - All[p.OffsetAll + i - 1] };
          if( diff < 0 )
          {
            std::ostringstream es;
            es << "Params::Export Monotonic invalid key " << it.first
               << "[" << i << "] " << All[p.OffsetAll + i]
               << " < [" << ( i - 1 ) << "] " << All[p.OffsetAll + i - 1];
            throw std::runtime_error( es.str().c_str() );
          }
          vType[Offset + i] = std::sqrt( diff );
        }
      }
    }
  }
}

void Params::Import( Vector &All, const VectorView &vType, Param::Type type ) const
{
  if( All.size < NumScalars( Param::Type::All ) )
    throw( "Params::Import() All is too short" );
  if( vType.size() < NumScalars( type ) )
  {
    std::ostringstream es;
    es << "Params::Import() " << type << " is too short";
    throw( es.str().c_str() );
  }
  for( const value_type &it : *this )
  {
    const Param &p{ it.second };
    if( p.type == type )
    {
      const std::size_t &Offset{ type == Param::Type::All ? p.OffsetAll : p.OffsetMyType };
      for( std::size_t i = 0; i < p.size; ++i )
      {
        if( i == 0 || !p.bMonotonic || type != Param::Type::Variable )
          All[p.OffsetAll + i] = vType[Offset + i];
        else
        {
          // Monotonically increasing set of values. Implemented as a sum of squares
          All[p.OffsetAll + i] = All[p.OffsetAll + i - 1] + vType[Offset + i] * vType[Offset + i];
        }
      }
    }
  }
}

void Params::Dump( std::ostream &os, const Vector &Values, Param::Type ShowType,
                   const Vector *pErrors, const std::vector<bool> *pbKnown ) const
{
  const Param::Type VectorType{ Param::Type::All };
  if( VectorType != Param::Type::All && VectorType != ShowType )
  {
    std::ostringstream es;
    es << "Params::Dump() Can't show " << ShowType << " from a vector containing " << VectorType;
    throw( es.str().c_str() );
  }
  if( Values.size < NumScalars( VectorType ) || ( pErrors && pErrors->size < NumScalars( VectorType ) ) )
  {
    std::ostringstream es;
    es << "Params::Dump() " << VectorType << " data too short";
    throw( es.str().c_str() );
  }
  for( value_type it : *this )
  {
    Param &p{ it.second };
    if( ShowType == Param::Type::All || p.type == ShowType )
    {
      const std::size_t &Offset{ VectorType == Param::Type::All ? p.OffsetAll : p.OffsetMyType };
      for( std::size_t i = 0; i < p.size; ++i )
      {
        if( !pbKnown || (*pbKnown)[Offset + i] )
        {
          os << std::string( MaxLen - p.FieldLen + 2, ' ' ) << it.first;
          if( p.size )
            os << i;
          os << Common::Space << Values[Offset + i];
          if( pErrors && p.type == Param::Type::Variable )
            os << "\t+/- " << (*pErrors)[p.OffsetMyType + i];
          os << Common::NewLine;
        }
      }
    }
  }
}
