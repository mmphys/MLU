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

MLU_Param_hpp

static const std::array<std::string,4> ParamTypeHuman{ "All", "Variable", "Fixed", "Derived" };

std::size_t Param::Key::Len() const
{
  std::size_t Size{ Name.length() };
  for( const std::string &s : Object )
    Size += s.length() + 1;
  return Size;
}

std::string Param::Key::FullName( std::size_t idx, std::size_t Size ) const
{
  std::ostringstream os;
  os << *this;
  if( Size > 1 )
    os << idx;
  return os.str();
}

std::string Param::Key::ShortName( std::size_t idx, std::size_t Size ) const
{
  std::string s{ Name };
  if( Size > 1 )
    s.append( std::to_string( idx ) );
  return s;
}

void Param::Key::ValidateKey()
{
  for( std::size_t i = 0; i < Object.size(); ++i )
    if( Object[i].empty() )
    {
      std::ostringstream es;
      es << "Object ID " << i << " of " << Object.size() << " empty";
      throw std::runtime_error( es.str().c_str() );
    }
}

bool Param::Key::SameObject( const Key &rhs ) const
{
  if( Object.size() != rhs.Object.size() )
    return false;
  for( std::size_t i = 0; i < Object.size(); ++i )
    if( !EqualIgnoreCase( Object[i], rhs.Object[i] ) )
      return false;
  return true;
}

bool Param::Key::Less::operator()( const Key &lhs, const Key &rhs ) const
{
  if( lhs.Object.size() != rhs.Object.size() )
    return lhs.Object.size() < rhs.Object.size();
  for( std::size_t i = 0; i < lhs.Object.size(); ++i )
  {
    int c{ CompareIgnoreCase( lhs.Object[i], rhs.Object[i] ) };
    if( c )
      return c < 0;
  }
  return CompareIgnoreCase( lhs.Name, rhs.Name ) < 0;
}

// Can this parameter be read from ListType?
void Param::Validate( Param::Type ListType ) const
{
  std::size_t i{ static_cast<std::size_t>( ListType ) };
  if( i >= ParamTypeHuman.size() )
  {
    std::ostringstream s;
    s << "Param::Type " << i << " invalid";
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
  std::size_t j{ static_cast<std::size_t>( type ) };
  if( j >= ParamTypeHuman.size() || type == Type::All )
  {
    std::ostringstream s;
    s << "Param::Type " << type << " invalid";
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
  if( type != ListType && ListType != Param::Type::All )
  {
    std::ostringstream s;
    s << "Param::Validate() can't read " << type << " parameter from " << ListType << " list";
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
}

// Underlying implementation for all operator()
std::size_t Param::GetOffset( std::size_t Idx, Type ListType ) const
{
  Validate( ListType );
  if( Idx >= size )
  {
    std::ostringstream s;
    s << "Param::GetOffset index " << Idx << " >= size " << size;
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
  return ( ListType == Param::Type::All ? OffsetAll : OffsetMyType ) + Idx;
}

std::size_t Param::operator()( std::size_t Idx, Param::Type ListType ) const
{
  if( !pKey ) //|| pKey->Object.size() != 1 ) // TODO: Do I need to support Matrix Element src/snk?
  {
    std::ostringstream s;
    s << "Param::operator() not 1-dimensional object";
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
  return GetOffset( Idx, ListType );
}

std::size_t Param::operator()( std::size_t idxSnk, std::size_t idxSrc, Type ListType ) const
{
  if( !pKey || pKey->Object.size() < 1 || pKey->Object.size() > 2 )
  {
    std::ostringstream s;
    s << "Param::operator() not 2-dimensional object";
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
  if( idxSnk > 1 || idxSrc > 1 )
    throw std::runtime_error( "Model3pt::ParamIndex index out of bounds" );
  if( bSwapSourceSink )
    std::swap( idxSnk, idxSrc );
  std::size_t idx{ idxSnk + idxSrc };
  if( idxSnk && pKey->Object.size() > 1 )
    idx++;
  return GetOffset( idx, ListType );
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

std::ostream &operator<<( std::ostream &os, const Param::Type &type )
{
  std::size_t i{ static_cast<std::size_t>( type ) };
  if( i < ParamTypeHuman.size() )
    os << ParamTypeHuman[i];
  else
    os << "Param::Type(" << i << ")";
  return os;
}

std::istream &operator>>( std::istream &is, Param::Type &Type )
{
  std::string s;
  if( is >> s )
  {
    const int Idx{ IndexIgnoreCase( ParamTypeHuman, s ) };
    if( Idx != ParamTypeHuman.size() )
    {
      Type = static_cast<Param::Type>( Idx );
      return is;
    }
  }
  throw std::runtime_error( "Param::Type \"" + s + "\" unrecognised" );
}

std::ostream &operator<<( std::ostream &os, const Param::Key &key )
{
  for( const std::string &s : key.Object )
    os << s << '-';
  return os << key.Name;
}

std::istream &operator>>( std::istream &is, Param::Key &key )
{
  std::string s;
  if( !( is >> s ) )
    throw std::runtime_error( "Param::Key empty" );
  key.Object = ArrayFromString( s, "-" );
  key.Name = std::move( key.Object.back() );
  key.Object.pop_back();
  return is;
}

std::ostream &operator<<( std::ostream &os, const Params::DispMap::value_type &dv )
{
  const Common::Param::Key &dk{ dv.first };
  const Common::DispEntry &de{ dv.second };
  return os << dk << " derived from " << de.ParentKey << " N=" << de.N << ", p="
            << de.p.to_string3d( Common::Comma );
}

std::ostream &operator<<( std::ostream &os, const Params::DispMap &dm )
{
  for( Params::DispMap::const_iterator it = dm.cbegin(); it != dm.cend() ; )
  {
    os << *it;
    if( ++it != dm.cend() )
      os << Common::NewLine;
  }
  return os;
}

Params::Params( const std::vector<std::string> &ParamNames )
{
  Param::Key k;
  for( const std::string &s : ParamNames )
  {
    k.Name = s;
    Add( k );
  }
  AssignOffsets();
}

Params::Params( const std::vector<std::string> &ParamNames, std::vector<std::size_t> &vIdx )
: Params( ParamNames )
{
  Common::Param::Key k;
  vIdx.resize( ParamNames.size() );
  for( std::size_t i = 0; i < ParamNames.size(); ++i )
  {
    k.Name = ParamNames[i];
    vIdx[i] = at( k )();
  }
}

Params::Params( const std::vector<Param::Key> &Keys )
{
  for( const Param::Key &k : Keys )
    Add( k );
  AssignOffsets();
}

Params::iterator Params::Find( const Param::Key &key, const std::string &ErrorPrefix )
{
  iterator it{ find( key ) };
  if( it == end() )
  {
    // Doesn't exist
    std::ostringstream ss;
    ss << ErrorPrefix << Space << key << " not found";
    throw std::runtime_error( ss.str().c_str() );
  }
  return it;
}

Params::const_iterator Params::Find( const Param::Key &key, const std::string &ErrorPrefix ) const
{
  const_iterator it{ find( key ) };
  if( it == cend() )
  {
    // Doesn't exist
    std::ostringstream ss;
    ss << ErrorPrefix << Space << key << " not found";
    throw std::runtime_error( ss.str().c_str() );
  }
  return it;
}

bool Params::operator==( const Params &rhs ) const
{
  if( size() != rhs.size() || NumScalars( Param::Type::All ) != rhs.NumScalars( Param::Type::All ) )
    return false;
  for( const value_type &it : *this )
  {
    const Param::Key &k{ it.first };
    const Param &p{ it.second };
    const_iterator it2{ rhs.find( k ) };
    if( it2 == rhs.end() || p.size != it2->second.size )
      return false;
  }
  return true;
}

bool Params::operator<=( const Params &rhs ) const
{
  for( const value_type &it : *this )
  {
    const Param::Key &k{ it.first };
    const Param &p{ it.second };
    const_iterator it2{ rhs.find( k ) };
    if( it2 == rhs.end() || p.size > it2->second.size )
      return false;
  }
  if( NumScalars( Param::Type::All ) != rhs.NumScalars( Param::Type::All ) )
    return false;
  return true;
}

bool Params::operator<( const Params &rhs ) const
{
  return NumScalars( Param::Type::All ) < rhs.NumScalars( Param::Type::All ) && *this <= rhs;
}

void Params::clear() noexcept
{
  Base::clear();
  dispMap.clear();
}

Params::iterator Params::Add( const Param::Key &key, std::size_t NumExp, bool bMonotonic, Param::Type type_ )
{
  // Work out how many elements are needed
  std::size_t size;
  if( key.Object.size() <= 1 ) // TODO: This won't work for matrix elements with src=snk
    size = NumExp;
  else if( key.Object.size() == 2 )
  {
    if( NumExp < 1 || NumExp > 3 )
      throw std::runtime_error( "Params::Add 3-point parameters support a maximum of 3 exponentials: "
                                "1) Gnd-Gnd, 2) Gnd-Ex (and ^\\dag), 3) Ex-Ex" );
    size = NumExp;
    if( NumExp >= 2 && key.Object.size() == 2 )
      ++size;
  }
  else
  {
    std::ostringstream es;
    es << "Params::Add key " << key << ": " << key.Object.size() << " members unsupported";
    throw std::runtime_error( es.str().c_str() );
  }
  if( !size )
  {
    std::ostringstream ss;
    ss << "Params::Add key " << key << " has zero members";
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
      ss << "Params::Add cannot add key " << key;
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
      const bool bCantGrow  { type_ != Param::Type::Fixed && size > p.size };
      if( bCantShrink || bCantGrow )
      {
        std::ostringstream s;
        s << "Param::Add key " << key << " can't " << ( bCantGrow ? "grow" : "shrink" )
          << " from " << p.type << "[" << p.size << ""
          << "] to " << type_ << "[" << size << "] parameters";
        throw std::runtime_error( s.str().c_str() );
      }
      // Can become fixed ... but can't go the other way
      if( type_ == Param::Type::Fixed )
      {
        p.type = type_;
        p.bMonotonic = bMonotonic; // Doesn't really do anything unless Variable
      }
    }
    // Has there been a request to change Monotonic?
    if( ( p.bMonotonic && !bMonotonic ) || ( !p.bMonotonic && bMonotonic ) )
    {
      // Silently ignore requests to remove monotonic from Variable parameters
      if( !( !bMonotonic && p.type != Param::Type::Fixed ) )
        p.bMonotonic = bMonotonic;
    }
    // Allow growth
    if( p.size < size )
      p.size = size;
  }
  return it;
}

Params::iterator Params::AddEnergy( const Param::Key &key, std::size_t NumExp, int N, bool bEnablePHat )
{
  bool bMonotonic{ true };
  Param::Type pt{ Param::Type::Variable };
  if( N )
  {
    // I would like to use the dispersion relation
    if( key.Object.size() != 1 )
    {
      std::ostringstream os;
      os << "Params::AddEnergy unexpected composite key " << key;
      throw std::runtime_error( os.str().c_str() );
    }
    std::string Key0{ key.Object[0] };
    Common::Momentum p;
    if( p.Extract( Key0 ) && p )
    {
      // I'm adding an energy parameter for non-zero momentum
      Common::Momentum p0( p.bp2 );
      Key0.append( p0.FileString( Common::Momentum::DefaultPrefix ) );
      Param::Key k0( std::move( Key0 ), key.Name );
      Add( k0, NumExp, true, Param::Type::Variable );
      dispMap.emplace( key, DispEntry( k0, N, p, bEnablePHat ) );
      bMonotonic = false;
      pt = Param::Type::Derived;
    }
  }
  return Add( key, NumExp, bMonotonic, pt );
}

template <typename T>
void Params::GuessEnergy( Vector<T> &Guess, std::vector<bool> &bKnown,
                          const Param::Key &key, std::size_t idx, T Energy ) const
{
  const_iterator it{ Find( key, "Params::GuessEnergy()" ) };
  const Param &p{ it->second };
  std::size_t idxEnergy{ p( idx ) };
  if( bKnown[idxEnergy] )
    std::cout << "  Ignoring guess " << key << "[" << idx << "] = " << Energy
              << ". Keeping previous guess " << Guess[idxEnergy] << NewLine;
  else
  {
    // Enforce monotonic
    if( p.bMonotonic && idx > 0 )
    {
      if( !bKnown[idxEnergy - 1] )
      {
        std::ostringstream os;
        os << "Can't guess " << key << "[" << idx << "] = " << Energy
           << " when " << key << "[" << (idx - 1) << "] unknown" << NewLine;
        throw std::runtime_error( os.str().c_str() );
      }
      if( Energy < Guess[idxEnergy - 1] )
      {
        std::cout << "  Adjusting guess " << key << "[" << idx << "] = " << Energy
                  << " to " << Guess[idxEnergy - 1] << NewLine;
        Energy = Guess[idxEnergy - 1];
      }
    }
    Guess[idxEnergy] = Energy;
    bKnown[idxEnergy] = true;
    // Now see whether we are using the dispersion relation for this
    DispMap::const_iterator itd{ dispMap.find( key ) };
    if( itd != dispMap.cend() )
    {
      // Propagate this up to parent using dispersion relation
      const DispEntry &de{ itd->second };
      const_iterator itParent{ Find( de.ParentKey, "Params::GuessEnergy() parent" ) };
      const Param &pParent{ itParent->second };
      std::size_t idxParent{ pParent() }; // Zero'th entry of the parent
      if( bKnown[idxParent] )
        std::cout << "  Ignoring parent guess " << de.ParentKey << "[0] = " << Energy
                  << ". Keeping previous guess " << Guess[idxParent] << NewLine;
      else
      {
        Guess[idxParent] = de.p.LatticeDispersion( Energy, de.N, de.bEnablePHat, true );
        bKnown[idxParent] = true;
      }
    }
    // Now propagate guesses down through all children
    PropagateEnergy( Guess, &bKnown );
  }
}

template void Params::GuessEnergy<float>( Vector<float> &Guess, std::vector<bool> &bKnown,
                                  const Param::Key &key, std::size_t idx, float Energy ) const;
template void Params::GuessEnergy<double>( Vector<double> &Guess, std::vector<bool> &bKnown,
                                  const Param::Key &key, std::size_t idx, double Energy ) const;

template <typename T>
void Params::PropagateEnergy( Vector<T> &Guess, std::vector<bool> *bKnown ) const
{
  bool bChanged;
  do
  {
    bChanged = false;
    // Loop through every entry in the dispersion table
    for( const Params::DispMap::value_type &v : dispMap )
    {
      static const std::string ErrorMsg{ "ParamsPairs::PropagateEnergy() Bug - dispersion" };
      // Get (derived) parameter in dispersion table
      const Param::Key &key{ v.first };
      const DispEntry &de{ v.second };
      Params::const_iterator it{ Find( key, ErrorMsg ) };
      const Param &p{ it->second };
      std::size_t idx{ p() };
      // Get parent parameter
      Params::const_iterator itParent{ Find( de.ParentKey, ErrorMsg ) };
      const Param &pParent{ itParent->second };
      std::size_t idxParent{ pParent() };
      if( pParent.size < p.size )
        throw std::runtime_error( "ParamsPairs::PropagateEnergy() parent size "
                + std::to_string( pParent.size ) + " < dispersion " + std::to_string( p.size ) );
      // Propagate each parameter
      for( std::size_t i = 0; i < p.size; ++i )
      {
        if( !bKnown || ( !(*bKnown)[idx + i] && (*bKnown)[idxParent + i] ) )
        {
          Guess[idx + i] = de.p.LatticeDispersion( Guess[idxParent + i], de.N, de.bEnablePHat );
          if( bKnown )
          {
            (*bKnown)[idx + i] = true;
            bChanged = true;
          }
        }
      }
    }
  }
  while( bChanged );
}

template void Params::PropagateEnergy<float>( Vector<float> &Guess, std::vector<bool> *bKnown ) const;
template void Params::PropagateEnergy<double>(Vector<double> &Guess, std::vector<bool> *bKnown) const;

// Make a parameter fixed. Must exist
Param &Params::SetType( const Param::Key &key, Param::Type type )
{
  static const std::string sErrPrefix{ "Params::SetType()" };
  if( type == Param::Type::All || type == Param::Type::Variable )
  {
    std::ostringstream os;
    os << sErrPrefix << " can't set type to " << type;
    throw os.str().c_str();
  }
  iterator it{ Find( key, sErrPrefix ) };
  Param &p{ it->second };
  p.type = type;
  p.bMonotonic = false;
  dispMap.erase( key ); // Also remove it from dispersion list (if present)
  return p;
}

// Make a parameter fixed. Must exist
Param &Params::SetType( const Param::Key &key, Param::Type type, bool bSwapSourceSink )
{
  Param &p{ SetType( key, type ) };
  p.bSwapSourceSink = bSwapSourceSink;
  return p;
}

Params::const_iterator Params::FindPromiscuous( const Param::Key &key ) const
{
  if( empty() )
    return end();
  if( bSingleObject )
    return find( Param::Key( begin()->first.Object[0], key.Name ) );
  return find( key );
}

void Params::AssignOffsets()
{
  const Param::Key * pFirstKey{ nullptr };
  bSingleObject = true;
  NumFixed = 0;
  NumVariable = 0;
  NumDerived = 0;
  MaxLen = 0;
  for( value_type &it : *this )
  {
    const Param::Key &k{ it.first };
    Param &p{ it.second };
    // See whether this is a multi-object set of parameters
    if( !pFirstKey )
      pFirstKey = &k;
    else if( bSingleObject )
      bSingleObject = pFirstKey->SameObject( k );
    // Save the pointer to the key
    p.pKey = &k;
    // Validate the parameter type
    p.Validate();
    // Assign positions of these parameters in tables
    p.OffsetAll = NumScalars( Param::Type::All );
    std::size_t &ThisNum{ SizeType( p.type) };
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

std::size_t Params::MaxExponents() const
{
  std::size_t Max{ 0 };
  for( const value_type &it : *this )
  {
    const Param::Key &k{ it.first };
    const Param &p{ it.second };
    std::size_t ThisSize{ p.size };
    if( k.Object.size() > 1 && ThisSize > 2 )
      ThisSize--;
    if( Max < ThisSize )
      Max = ThisSize;
  }
  return Max;
}

template <typename T>
void Params::Export( Vector<T> &vType, const Vector<T> &All, Param::Type type ) const
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
        if( i == 0 || !p.bMonotonic || type == Param::Type::Fixed )
          vType[Offset + i] = All[p.OffsetAll + i];
        else
        {
          // Monotonically increasing set of values. Implemented as a sum of squares
          const T diff{ All[p.OffsetAll + i] - All[p.OffsetAll + i - 1] };
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

template void Params::Export<float>( Vector<float> &vType,
                const Vector<float> &All, Param::Type type ) const;
template void Params::Export<double>( Vector<double> &vType,
                const Vector<double> &All, Param::Type type ) const;

// Import values
// If Ref is given, then we are importing errors, which we add in quadrature
template <typename T>
void Params::Import( Vector<T> &All, const VectorView<T> &Source, Param::Type SourceType,
                     bool bSourceMonotonic, const Vector<T> * pRef ) const
{
  if( All.size < NumScalars( Param::Type::All ) )
    throw( "Params::Import() All is too short" );
  if( Source.size() < NumScalars( SourceType ) )
  {
    std::ostringstream es;
    es << "Params::Import() " << SourceType << " is too short";
    throw( es.str().c_str() );
  }
  for( const value_type &it : *this )
  {
    const Param &p{ it.second };
    if( p.type == SourceType )
    {
      const std::size_t &Offset{ SourceType == Param::Type::All ? p.OffsetAll : p.OffsetMyType };
      for( std::size_t i = 0; i < p.size; ++i )
      {
        if( i == 0 || !p.bMonotonic || !bSourceMonotonic || SourceType == Param::Type::Fixed )
          All[p.OffsetAll + i] = Source[Offset + i];
        else
        {
          // Monotonically increasing set of values. Implemented as a sum of squares
          T Previous{ All[p.OffsetAll + i - 1] };
          T Current{ Source[Offset + i] };
          if( pRef )
          {
            const T PrevVal{ (*pRef)[p.OffsetAll + i - 1] };
            const T CurrVal{ (*pRef)[p.OffsetAll + i] };
            const T RelE1{ Previous / PrevVal };
            const T RelE2{ Current / CurrVal };
            All[p.OffsetAll + i] = std::sqrt( RelE1 * RelE1 + RelE2 * RelE2 );
          }
          else
            All[p.OffsetAll + i] = Previous + Current * Current;
        }
      }
    }
  }
}

template void Params::Import<float>( Vector<float> &All, const VectorView<float> &vType,
              Param::Type type, bool bSourceMonotonic, const Vector<float> * pRef ) const;
template void Params::Import<double>( Vector<double> &All, const VectorView<double> &vType,
              Param::Type type, bool bSourceMonotonic, const Vector<double> * pRef ) const;

template <typename T>
void Params::Dump( std::ostream &os, const Vector<T> &Values, Param::Type ShowType,
                   const Vector<T> *pErrors, const std::vector<bool> *pbKnown ) const
{
  const Param::Type VectorType{ Param::Type::All };
  if( VectorType != Param::Type::All && VectorType != ShowType )
  {
    std::ostringstream es;
    es << "Params::Dump() Can't show " << ShowType << " from a vector containing " << VectorType;
    throw( es.str().c_str() );
  }
  if( Values.size < NumScalars( VectorType ) || ( pErrors && pErrors->size < NumScalars( Param::Type::Variable ) ) )
  {
    std::ostringstream es;
    es << "Params::Dump() " << VectorType << " data too short";
    throw( std::runtime_error( es.str().c_str() ) );
  }
  for( value_type it : *this )
  {
    Param &p{ it.second };
    if( ShowType == Param::Type::All || p.type == ShowType )
    {
      const std::size_t &Offset{ VectorType == Param::Type::All ? p.OffsetAll : p.OffsetMyType };
      for( std::size_t i = 0; i < p.size; ++i )
      {
        os << std::string( MaxLen - p.FieldLen + 2, ' ' ) << it.first;
        if( p.size > 1 )
          os << i;
        os << " ";
        if( !pbKnown || (*pbKnown)[Offset + i] )
        {
          os << Values[Offset + i];
          if( pErrors && p.type == Param::Type::Variable )
            os << "\t+/- " << (*pErrors)[p.OffsetMyType + i];
        }
        else
          os << "unknown";
        os << "\n";
      }
    }
  }
}

template void Params::Dump<float>( std::ostream &os, const Vector<float> &Values, Param::Type ShowType,
                                   const Vector<float> *pErrors, const std::vector<bool> *pbKnown ) const;
template void Params::Dump<double>( std::ostream &os, const Vector<double> &Values, Param::Type ShowType,
                                    const Vector<double> *pErrors, const std::vector<bool> *pbKnown ) const;

// Get the name of the parameters. Short names don't include object name - but has to be only 1 object
std::vector<std::string> Params::GetNames( Param::Type type, bool bLongNames ) const
{
  if( !bLongNames )
  {
    // We must also show long names if any object names differ
    const Param::Key *FirstKey = nullptr;
    for( const value_type &it : *this )
    {
      const Param &p{ it.second };
      if( type == Param::Type::All || p.type == type )
      {
        const Param::Key &k{ it.first };
        if( !FirstKey )
          FirstKey = &k;
        else if( !FirstKey->SameObject( k ) )
        {
          // Object names differ, shjow long names
          bLongNames = true;
          break;
        }
      }
    }
  }
  std::vector<std::string> Names( NumScalars( type ) );
  for( const value_type &it : *this )
  {
    const Param &p{ it.second };
    if( type == Param::Type::All || p.type == type )
    {
      const std::size_t &Offset{ type == Param::Type::All ? p.OffsetAll : p.OffsetMyType };
      const Param::Key &k{ it.first };
      for( std::size_t i = 0; i < p.size; ++i )
        Names[Offset + i] = bLongNames ? k.FullName( i, p.size ) : k.ShortName( i, p.size );
    }
  }
  return Names;
}

void Params::WriteNames( std::ostream &os ) const
{
  // Write parameter names
  bool bFirst{ true };
  for( const Params::value_type &it : *this )
  {
    const Param &p{ it.second };
    for( std::size_t i = 0; i < p.size; ++i )
    {
      if( bFirst )
        bFirst = false;
      else
        os << Space;
      os << GetName( it, i );
    }
  }
}

void Params::WriteNamesValWithEr( std::ostream &os ) const
{
  // Write parameter names
  bool bFirst{ true };
  for( const Params::value_type &it : *this )
  {
    const Param &p{ it.second };
    for( std::size_t i = 0; i < p.size; ++i )
    {
      if( bFirst )
        bFirst = false;
      else
        os << Space;
      ValWithEr<double>::Header( GetName( it, i ), os, Space ); // Header is type independent
    }
  }
}

void Params::ReadH5 ( ::H5::Group gParent, const std::string GroupName )
{
  clear();
  std::string sError;
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  try
  {
    // Open the group
    ::H5::Group g{ gParent.openGroup( GroupName ) };
    // Get the object IDs TODO: add support for renaming the objects
    std::vector<std::string> vObjects;
    try
    {
      vObjects = H5::ReadStrings( g.openAttribute( sObjectNames ) );
    }
    catch( const ::H5::Exception &e )
    {
      H5::GetErrorClearStack( e );
    }
    // Get number of params
    int iNumParams;
    ::H5::Attribute a{ g.openAttribute( sCount ) };
    a.read( ::H5::PredType::NATIVE_INT, &iNumParams );
    a.close();
    // Load each param
    std::string SubGroup{ GroupName };
    const std::size_t SubGroupLen{ SubGroup.length() };
    for( std::size_t i = 0; i < iNumParams; ++i )
    {
      Param::Key k;
      // Open sub-group for each parameter
      SubGroup.resize( SubGroupLen );
      SubGroup.append( std::to_string( i ) );
      ::H5::Group gSub{ g.openGroup( SubGroup ) };
      // Read in the Object IDs
      try
      {
        std::vector<int> vOID;
        a = gSub.openAttribute( sObjectNames );
        H5::ReadVectorHelper( vOID, a, sObjectNames );
        k.Object.reserve( vOID.size() );
        for( int i : vOID )
          k.Object.push_back( vObjects[i] );
        a.close();
      }
      catch( const ::H5::Exception &e )
      {
        H5::GetErrorClearStack( e );
      }
      // Parameter name
      a = gSub.openAttribute( sName );
      a.read( a.getStrType(), k.Name );
      a.close();
      // Parameter size
      a = gSub.openAttribute( sSize );
      int iSize;
      a.read( ::H5::PredType::NATIVE_INT, &iSize );
      a.close();
      // Get the type (default Variable)
      Param::Type type{ Param::Type::Variable };
      if( H5::OpenOptional( a, gSub, sTypeName ) )
      {
        std::string s;
        a.read( a.getStrType(), s );
        a.close();
        type = FromString<Param::Type>( s );
      }
      // Read monotonic
      bool bMonotonic{ false };
      if( H5::OpenOptional( a, gSub, sMonotonic ) )
      {
        std::int8_t i8;
        a.read( ::H5::PredType::NATIVE_INT8, &i8 );
        a.close();
        if( i8 )
          bMonotonic = true;
      }
      // Now create this parameter
      insert( std::make_pair( std::move( k ), Param( iSize, bMonotonic, type ) ) );
    }
  }
  catch( const ::H5::Exception &e )
  {
    clear();
    sError = H5::GetErrorClearStack( e );
  }
  catch(...)
  {
    clear();
    ::H5::Exception::setAutoPrint( h5at, f5at_p );
    throw;
  }
  ::H5::Exception::setAutoPrint( h5at, f5at_p );
  if( sError.size() )
    throw std::runtime_error( sError );
  AssignOffsets();
}

void Params::WriteH5( ::H5::Group gParent, const std::string GroupName ) const
{
  ::H5::Group g{ gParent.createGroup( GroupName ) };
  if( empty() )
    return;
  // Get list of object IDs
  UniqueNames OID;
  for( value_type it : *this )
  {
    const Param::Key &k{ it.first };
    for( std::size_t i = 0; i < k.Object.size(); ++i )
      OID[k.Object[i]];
  }
  int iNumObjects = 0;
  for( typename UniqueNames::value_type &it : OID )
    it.second = iNumObjects++;
  if( iNumObjects )
  {
    // Write all the object names
    std::vector<std::string> vObjects( iNumObjects );
    int i{ 0 };
    for( typename UniqueNames::value_type &it : OID )
      vObjects[i++] = it.first;
    H5::WriteAttribute( g, sObjectNames, vObjects );
  }
  hsize_t Dims[1];
  Dims[0] = 1;
  ::H5::DataSpace ds1( 1, Dims );
  {
    int iNumParams{ static_cast<int>( size() ) };
    ::H5::Attribute a = g.createAttribute( sCount, ::H5::PredType::STD_U16LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &iNumParams );
  }
  // Now loop through and write each parameter
  {
    int i{ 0 };
    std::string SubGroup{ GroupName };
    const std::size_t SubGroupLen{ SubGroup.length() };
    for( value_type it : *this )
    {
      const Param::Key &k{ it.first };
      const Param &p{ it.second };
      // Create a sub-group for each parameter
      SubGroup.resize( SubGroupLen );
      SubGroup.append( std::to_string( i++ ) );
      ::H5::Group gSub{ g.createGroup( SubGroup ) };
      // Write all the object IDS for this parameter
      if( ! k.Object.empty() )
      {
        std::vector<int> vOIDs( k.Object.size() );
        for( std::size_t j = 0; j < k.Object.size(); ++j )
          vOIDs[j] = OID[k.Object[j]];
        Dims[0] = vOIDs.size();
        ::H5::DataSpace dsOID( 1, Dims );
        ::H5::Attribute a = gSub.createAttribute( sObjectNames, ::H5::PredType::STD_U16LE, dsOID );
        a.write( ::H5::PredType::NATIVE_INT, &vOIDs[0] );
      }
      {
        // Parameter name
        ::H5::Attribute a = gSub.createAttribute( sName, H5::Equiv<std::string>::Type, ds1 );
        a.write( H5::Equiv<std::string>::Type, k.Name );
      }
      // Parameter Type - default is variable
      if( p.type != Param::Type::Variable )
      {
        std::ostringstream ss;
        ss << p.type;
        ::H5::Attribute a = gSub.createAttribute( sTypeName, H5::Equiv<std::string>::Type, ds1 );
        a.write( H5::Equiv<std::string>::Type, ss.str() );
      }
      {
        // Size
        int iSize{ static_cast<int>( p.size ) };
        ::H5::Attribute a = gSub.createAttribute( sSize, ::H5::PredType::STD_U16LE, ds1 );
        a.write( ::H5::PredType::NATIVE_INT, &iSize );
      }
      const std::int8_t i8{ 1 };
      if( p.bMonotonic )
      {
        ::H5::Attribute a = gSub.createAttribute( sMonotonic, ::H5::PredType::STD_U8LE, ds1 );
        a.write( ::H5::PredType::NATIVE_INT8, &i8 );
      }
    }
  }
}

std::string Params::GetName( const value_type &param, std::size_t idx ) const
{
  const Param::Key &k{ param.first };
  const Param &p{ param.second };
  if( SingleObject() )
    return k.ShortName( idx, p.size );
  return k.FullName( idx, p.size );
}

void Params::KeepCommon( const Params &Other )
{
  for( iterator it = begin(); it != end(); )
  {
    const Param::Key &k{ it->first };
    const_iterator it2{ Other.find( k ) };
    if( it2 == Other.end() )
    {
      // Not in the other list - delete my copy
      it = erase( it );
    }
    else
    {
      // In the other list - shrink if other is smaller
      Param &p{ it->second };
      if( p.size > it2->second.size )
        p.size = it2->second.size;
      if( p.size )
        ++it;
      else
        it = erase( it );
    }
  }
}

void Params::Merge( const Params &Other )
{
  for( const value_type &it : Other )
  {
    const Param::Key &k{ it.first };
    iterator itMe{ find( k ) };
    if( itMe == end() )
    {
      // Not in my list - add it
      insert( it );
    }
    else
    {
      // In the other list - grow if other is larger
      Param &p{ itMe->second };
      if( p.size < it.second.size )
        p.size = it.second.size;
    }
  }
}

template <typename T> void Params::AdjustSigns( Vector<T> &data ) const
{
  for( std::size_t i = 0; i < SignList.size(); ++i )
  {
    if( data[SignList[i][0]] < 0 )
    {
      for( std::size_t j = 0; j < SignList[i].size(); ++j )
      {
        data[SignList[i][j]] = -data[SignList[i][j]];
      }
    }
  }
}

void Params::DumpSignList( std::ostream &os ) const
{
  os << "Sign indices: ";
  for( std::size_t i = 0; i < SignList.size(); ++i )
  {
    if( i )
      os << "; ";
    for( std::size_t j = 0; j < SignList[i].size(); ++j )
    {
      if( j )
        os << Comma;
      os << SignList[i][j];
    }
  }
  os << NewLine;
}

template void Params::AdjustSigns<float>( Vector<float> &data ) const;
template void Params::AdjustSigns<double>( Vector<double> &data ) const;

const std::string Params::sCount{ "Count" };
const std::string Params::sObjectNames{ "Object" };
const std::string Params::sName{ "Name" };
const std::string Params::sTypeName{ "Type" };
const std::string Params::sSize{ "Size" };
const std::string Params::sMonotonic{ "Monotonic" };

std::ostream &operator<<( std::ostream &os, const Params::value_type &param )
{
  const Param &p{ param.second };
  os << param.first << '[' << p.size;
  if( p.type != Param::Type::Variable )
    os << ',' << p.type;
  if( p.bMonotonic )
    os << 'M';
  return os << ']';
}

ParamsPairs::StateSize ParamsPairs::GetStateSize() const
{
  StateSize ss;
  for( const State &s : keystate )
  {
    switch( s )
    {
      case State::Known:
        ++ss.Known;
        break;
      case State::Unknown:
        ++ss.Unknown;
        break;
      default:
        ++ss.Other;
    }
  }
  for( const Pair &pair : pairs )
  {
    std::array<const State *, Pair::size> aState{ GetPairState( pair, true ) };
    if( *aState[0] == State::ProductOnly && *aState[1] == State::ProductOnly )
    {
      ss.Other -= 2;
      ss.Unknown += 2;
    }
  }
  return ss;
}

void ParamsPairs::clear()
{
  keystate.resize( params.NumScalars( Param::Type::All ) );
  for( State &s : keystate )
    s = State::Unknown;
  pairs.clear();
  // Initialise state of Fixed=Known and Derived/Variable=Unknown parameters
  for( const Params::value_type &it : params )
  {
    const Param &p{ it.second };
    if( p.type == Param::Type::Fixed )
    {
      for( std::size_t i = 0; i < p.size; ++i )
        keystate[p(i)] = State::Known;
    }
  }
  PropagateDisp(); // Propagate anything known down through dispersion relations
}

void ParamsPairs::SetState( State NewState, const Param::Key &key, std::size_t Size, std::size_t Idx )
{
  Params::const_iterator it{ params.find( key ) };
  if( it == params.cend() )
  {
    std::ostringstream os;
    os << "ParamsPairs::SetState() key not found " << key << "[" << Idx << "]x" << Size
       << " setting new state " << NewState;
    throw std::runtime_error( os.str().c_str() );
  }
  SetState( NewState, it, Size, Idx );
}

void ParamsPairs::KnowProduct( const Key &key0, const Key &key1 )
{
  // Find each key's state
  bool bChanged{ false };
  Pair pProd( key0, key1 );
  std::array<State *, Pair::size> aState{ GetPairState( pProd, true ) };
  pairs.insert( std::move( pProd ) ); // It's ok if they are in the pair list already
  const State newState{ MaxState( *aState[0], State::ProductOnly ) };
  if( *aState[0] != newState )
  {
    *aState[0] = newState;
    bChanged = true;
  }
  if( PromoteState( aState ) )
    bChanged = true;
  if( bChanged )
    PropagateKnown();
}

void ParamsPairs::KnowProduct( const Param::Key &key0, const Param::Key &key1, std::size_t Size )
{
  Key k0( key0, 0 );
  Key k1( key1, 0 );
  for( k0.Index = 0; k0.Index < Size; ++k0.Index )
  {
    k1.Index = k0.Index;
    KnowProduct( k0, k1 );
  }
}

bool ParamsPairs::HasUnknownProducts() const
{
  for( const PairSet::value_type &pair : pairs )
  {
    std::array<const State *, Pair::size> aState{ GetPairState( pair, true ) };
    if( *aState[0] == State::ProductOnly && *aState[1] == State::ProductOnly )
      return true;
  }
  return false;
}

std::size_t ParamsPairs::NumKnownProducts( const Param::Key &k0, const Param::Key &k1,
                                           std::size_t Size ) const
{
  Pair pair( Key( k0, Size ), Key( k1, Size ) );
  while( pair[0].Index-- )
  {
    pair[1].Index = pair[0].Index;
    std::array<const ParamsPairs::State *, ParamsPairs::Pair::size> aState = GetPairState( pair );
    if( *aState[0] != State::ProductOnly || *aState[1] != State::ProductOnly )
      return pair[0].Index + 1;
  }
  return 0;
}

void ParamsPairs::GetPairState( std::array<Params::const_iterator, ParamsPairs::Pair::size> &it,
                                const Pair &pair, bool bUnknownOk ) const
{
  for( std::size_t i = 0; i < it.size(); ++i )
  {
    it[i] = params.Find( pair[i], "ParamsPairs::GetPairState()" );
    const Param &p{ it[i]->second };
    State state = keystate[ p( pair[i].Index ) ];
    if( !bUnknownOk && state == State::Unknown ) // Shouldn't happen (this is debugging)
    {
      std::ostringstream os;
      os << "ParamsPairs::GetPairState() key " << pair[i] << " state " << state;
      throw std::runtime_error( os.str().c_str() );
    }
  }
}

std::array<const ParamsPairs::State *, ParamsPairs::Pair::size>
ParamsPairs::GetPairState( const Pair &pair, bool bUnknownOk ) const
{
  std::array<Params::const_iterator, ParamsPairs::Pair::size> it;
  GetPairState( it, pair, bUnknownOk );
  std::array<const State *, Pair::size> aState;
  for( std::size_t i = 0; i < aState.size(); ++i )
  {
    const Param &p{ it[i]->second };
    aState[i] = &keystate[ p( pair[i].Index ) ];
  }
  return aState;
}

std::array<ParamsPairs::State *, ParamsPairs::Pair::size>
ParamsPairs::GetPairState( const Pair &pair, bool bUnknownOk )
{
  std::array<Params::const_iterator, ParamsPairs::Pair::size> it;
  GetPairState( it, pair, bUnknownOk );
  std::array<State *, Pair::size> aState;
  for( std::size_t i = 0; i < aState.size(); ++i )
  {
    const Param &p{ it[i]->second };
    aState[i] = &keystate[ p( pair[i].Index ) ];
  }
  return aState;
}

// Promote other member of the pair to the highest state. Return true if something changed
bool ParamsPairs::PromoteState( std::array<State *, Pair::size> &aState )
{
  bool bChanged{ *aState[0] != *aState[1] };
  if( bChanged )
  {
    const State maxState{ MaxState( *aState[0], *aState[1] ) };
    *aState[0] = maxState;
    *aState[1] = maxState;
  }
  return bChanged;
}

// Promote other member of the pair to the highest state. Return true if something changed
bool ParamsPairs::PromoteState( const Pair &pair )
{
  std::array<State *, Pair::size> aState{ GetPairState( pair, true ) };
  return PromoteState( aState );
}

void ParamsPairs::SetState( State NewState, Params::const_iterator &it, std::size_t Size,
                            std::size_t Index )
{
  const Param::Key &key{ it->first };
  if( NewState != State::Known && NewState != State::AmbiguousSign )
  {
    std::ostringstream os;
    os << "ParamsPairs::SetState() key " << key << "[" << Index << "]x" << Size
       << "=" << NewState << " invalid";
    throw std::runtime_error( os.str().c_str() );
  }
  const Param &p{ it->second };
  bool bNeedPropagate{ false };
  for( std::size_t i = 0; i < Size; ++i )
  {
    State &state{ keystate[ p( Index + i ) ] };
    if( NewState == state || ( NewState == State::AmbiguousSign && state == State::Known ) )
    {
      // Silently ignore when state doesn't change, or we already know for sure
    }
    else
    {
      // Transition to the new state
      state = NewState;
      bNeedPropagate = true;
    }
  }
  if( bNeedPropagate )
    PropagateKnown();
  PropagateDisp(); // Propagate anything known down through dispersion relations
}

void ParamsPairs::PropagateKnown()
{
  bool bChanged;
  do
  {
    bChanged = false;
    for( const Pair &pair : pairs )
    {
      if( PromoteState( pair ) )
        bChanged = true;
    }
  }
  while( bChanged );
}

void ParamsPairs::PropagateDisp()
{
  bool bChanged;
  do
  {
    bChanged = false;
    // Loop through every entry in the dispersion table
    for( const Params::DispMap::value_type &v : params.dispMap )
    {
      static const std::string ErrorMsg{ "ParamsPairs::PropagateDisp() Bug - dispersion" };
      // Get (derived) parameter in dispersion table
      const Param::Key &key{ v.first };
      const DispEntry &de{ v.second };
      Params::const_iterator it{ params.Find( key, ErrorMsg ) };
      const Param &p{ it->second };
      std::size_t idx{ p() };
      // Get parent parameter
      Params::const_iterator itParent{ params.Find( de.ParentKey, ErrorMsg ) };
      const Param &pParent{ itParent->second };
      std::size_t idxParent{ pParent() };
      if( pParent.size < p.size )
        throw std::runtime_error( "ParamsPairs::PropagateDisp() parent size "
                + std::to_string( pParent.size ) + " < dispersion " + std::to_string( p.size ) );
      // Propagate each parameter
      for( std::size_t i = 0; i < p.size; ++i )
      {
        if( keystate[idx + i] != keystate[idxParent + i] )
        {
          if(     keystate[idxParent+i] == State::Known
             || ( keystate[idxParent+i] == State::AmbiguousSign && keystate[idx+i] != State::Known ) )
          {
            // Propagate downwards from parent
            keystate[idx+i] = keystate[idxParent+i];
            bChanged = true;
          }
          else if( keystate[idx+i] == State::Known
              || ( keystate[idx+i] == State::AmbiguousSign && keystate[idxParent+i] != State::Known))
          {
            // Propagate upwards to parent
            keystate[idxParent + i] = keystate[idx + i];
            bChanged = true;
          }
        }
      }
    }
  }
  while( bChanged );
}

bool ParamsPairs::Key::operator==( const Key &rhs ) const
{
  const Param::Key &pkRhs{ rhs };
  const Param::Key &pkLhs{ *this };
  return pkLhs == pkRhs && Index == rhs.Index;
}

bool ParamsPairs::Key::Less::operator()( const Key &lhs, const Key &rhs ) const
{
  const Param::Key &pkLhs{ lhs };
  const Param::Key &pkRhs{ rhs };
  if( pkLhs != pkRhs )
    return Param::Key::Less()( pkLhs, pkRhs );
  return lhs.Index < rhs.Index;
}

ParamsPairs::Pair::Pair( const Key &key0, const Key &key1 ) : std::array<Key, 2>{ key0, key1 }
{
  if( key0 == key1 )
  {
    std::ostringstream os;
    os << "ParamsPairs::Pair::Pair() repeated key " << key0;
    throw std::runtime_error( os.str().c_str() );
  }
  if( Param::Key::Less()( key1, key0 ) )
    std::swap( (*this)[0], (*this)[1] );
}

bool ParamsPairs::Pair::Less::operator()( const Pair &lhs, const Pair &rhs ) const
{
  if( lhs[0] != rhs[0] )
    return ParamsPairs::Key::Less()( lhs[0], rhs[0] );
  return ParamsPairs::Key::Less()( lhs[1], rhs[1] );
}

std::ostream &operator<<( std::ostream &os, const ParamsPairs::State &state )
{
  switch( state )
  {
    case ParamsPairs::State::Unknown:
      os << "unknown";
      break;
    case ParamsPairs::State::Known:
      os << "known";
      break;
    case ParamsPairs::State::AmbiguousSign:
      os << "ambiguousSign";
      break;
    case ParamsPairs::State::ProductOnly:
      os << "productOnly";
      break;
    default:
      os << "undefinedState"
         << static_cast<typename std::underlying_type<ParamsPairs::State>::type>( state );
      break;
  }
  return os;
}

std::ostream &operator<<( std::ostream &os, const ParamsPairs::Key &key )
{
  const Param::Key &base{ key };
  os << base << "[" << key.Index << "]";
  return os;
}

std::ostream &operator<<( std::ostream &os, const ParamsPairs::Pair &pair )
{
  for( std::size_t i = 0; i < pair.size; ++i )
  {
    if( i )
      os << Comma;
    os << pair[i];
  }
  return os;
}

std::ostream &operator<<( std::ostream &os, const ParamsPairs &PP )
{
  const char Indent[] = "  ";
  for( const Params::value_type &vt : PP.params )
  {
    const Param::Key &key{ vt.first };
    const Param &p{ vt.second };
    const std::size_t Index{ p() };
    os << Indent << key << Space;
    for( std::size_t i = 0; i < p.size; ++i )
    {
      if( i )
        os << CommaSpace;
      os << "[" << i << "]=" << PP.keystate[Index + i];
    }
    os << NewLine;
  }
  // Now show pairs
  if( !PP.pairs.empty() )
  {
    os << Indent << "Pairs:";
    for( const ParamsPairs::PairSet::value_type &pair : PP.pairs )
      os << Space << pair;
    os << NewLine;
  }
  return os;
}

SignChoice::SignChoice( const ParamsPairs &PP_, bool bShowSigns_ )
: PP{ PP_ }, bShowSigns{bShowSigns_}
{
  // Add all the pairs
  for( const ParamsPairs::Pair &pair : PP.pairs )
    Add( pair );
  // Add individual AmbiguousSign objects
  for( const Params::value_type &it : PP.params )
  {
    const Param p{ it.second };
    Key key{ it.first, 0 };
    for( key.Index = 0; key.Index < p.size; ++key.Index )
    {
      std::size_t Index{ p( key.Index ) };
      const ParamsPairs::State &state{ PP.keystate[Index] };
      SignChoice::List::iterator lit;
      SignChoice::Set::iterator sit;
      if( state == ParamsPairs::State::AmbiguousSign && !find( key, lit, sit ) )
        list.emplace_back( Set( { key } ) );
    }
  }
}

void SignChoice::Add( const Pair &pair )
{
  std::array<bool, Pair::size> bFound;
  std::array<List::iterator, Pair::size> lit;
  std::array<Set::iterator, Pair::size> sit;
  find( pair, bFound, lit, sit );
  const int NumFound{ ( bFound[0] ? 1 : 0 ) + ( bFound[1] ? 1 : 0 ) };
  if( NumFound == 2 )
  {
    // Found both. Merge if each is in a different set
    if( lit[0] != lit[1] )
    {
      for( Set::iterator it = lit[1]->begin(); it != lit[1]->end(); )
      {
        if( !lit[0]->insert( *it ).second ) // insert to new
        {
          // Couldn't insert
          std::ostringstream os;
          os << "SignChoice::Add() couldn't move key " << (*it);
          throw std::runtime_error( os.str().c_str() );
        }
        it = lit[1]->erase( it ); // remove from old
      }
      list.erase( lit[1] );
    }
  }
  else if( NumFound == 1 )
  {
    // Add the new member to same set as first
    int kExist{ bFound[0] ? 0 : 1 };
    int kNew{ bFound[0] ? 1 : 0 };
    if( !lit[kExist]->insert( pair[kNew] ).second )
    {
      // Couldn't insert
      std::ostringstream os;
      os << "SignChoice::Add() couldn't insert key " << pair[kNew];
      throw std::runtime_error( os.str().c_str() );
    }
  }
  else
  {
    // Couldn't find either - they go in a new list
    lit[0] = list.insert( list.end(), Set( { pair[0], pair[1] } ) );
  }
}

bool SignChoice::find( const Key &key, List::iterator &lit, Set::iterator &sit )
{
  for( lit = list.begin(); lit != list.end(); ++lit )
  {
    sit = lit->find( key );
    if( sit != lit->end() )
      return true;
  }
  return false;
}

SignChoice::operator std::vector<std::vector<std::size_t>>() const
{
  // Convert SignChoice into list of groups of parameters to sign-flip collectively
  std::vector<std::vector<std::size_t>> SignList( list.size() );
  std::size_t i{ 0 };
  static const std::string ErrorMsg{ "SignChoice::operator std::vector<std::vector<std::size_t>>()" };
  for( const SignChoice::Set &s : list )
  {
    std::vector<std::size_t> &l{ SignList[i++] };
    l.reserve( s.size() );
    // Start with an AmbiguousSign
    SignChoice::Set::iterator itAmbig{ s.size() > 1 ? s.cbegin() : s.cend() };
    for( ; itAmbig != s.cend(); ++itAmbig )
    {
      const Key &key{ *itAmbig };
      Params::const_iterator pit{ PP.params.Find( key, ErrorMsg ) };
      const Param &p{ pit->second };
      std::size_t Index{ p( key.Index ) };
      const ParamsPairs::State &state{ PP.keystate[Index] };
      if( state == ParamsPairs::State::AmbiguousSign )
      {
        l.emplace_back( Index );
        break;
      }
    }
    // Now add all the other entries
    for( SignChoice::Set::iterator it = s.cbegin(); it != s.cend(); ++it )
    {
      if( it != itAmbig )
      {
        const Key &key{ *it };
        Params::const_iterator pit{ PP.params.Find( key, ErrorMsg ) };
        const Param &p{ pit->second };
        l.emplace_back( p( key.Index ) );
      }
    }
  }
  return SignList;
}

void SignChoice::find( const Pair &pair, std::array<bool, Pair::size> &bFound,
                       std::array<List::iterator, Pair::size> &lit,
                       std::array<Set::iterator, Pair::size> &sit )
{
  for( std::size_t i = 0; i < Pair::size; ++i )
    bFound[i] = find( pair[i], lit[i], sit[i] );
}

std::ostream &operator<<( std::ostream &os, const SignChoice &sc )
{
  bool bFirst{ true };
  for( const SignChoice::List::value_type &l : sc.list )
  {
    if( bFirst )
      bFirst = false;
    else
      os << "; ";
    bool bOne{ true };
    for( const SignChoice::Set::value_type &s : l )
    {
      if( bOne )
        bOne = false;
      else
        os << CommaSpace;
      os << s;
    }
  }
  return os;
}

MLU_Param_hpp_end
