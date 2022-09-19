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

static const std::array<std::string,3> ParamTypeHuman{ "All", "Variable", "Fixed" };

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

bool Param::Key::operator==( const Key &rhs ) const
{
  if( Object.size() != rhs.Object.size() )
    return false;
  for( std::size_t i = 0; i < Object.size(); ++i )
   if( CompareIgnoreCase( Object[i], rhs.Object[i] ) )
      return false;
  return EqualIgnoreCase( Name, rhs.Name );
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
    s << "Params::Type " << i << " invalid";
    if( pKey )
      s << ", Key=" << *pKey;
    throw std::runtime_error( s.str().c_str() );
  }
  std::size_t j{ static_cast<std::size_t>( type ) };
  if( j >= ParamTypeHuman.size() || type == Type::All )
  {
    std::ostringstream s;
    s << "Params::Type " << type << " invalid";
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
  if( !pKey || pKey->Object.size() != 1 )
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

Params::iterator Params::Add( const Param::Key &key, std::size_t NumExp, bool bMonotonic, Param::Type type_ )
{
  // Work out how many elements are needed
  std::size_t size;
  if( key.Object.size() == 1 )
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
      const bool bCantGrow{ type_ == Param::Type::Variable && size > p.size };
      if( bCantShrink || bCantGrow )
      {
        std::ostringstream s;
        s << "Param::Add key " << key << " can't " << ( bCantGrow ? "grow" : "shrink" )
          << " from " << p.type << "[" << p.size << ""
          << "] to " << type_ << "[" << size << "] parameters";
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
    // Save the pointer to the key
    p.pKey = &k;
    // Validate the parameter type
    p.Validate();
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
        if( i == 0 || !p.bMonotonic || type != Param::Type::Variable )
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
                const Vector<float> &All, Param::Type type = Param::Type::Variable ) const;
template void Params::Export<double>( Vector<double> &vType,
                const Vector<double> &All, Param::Type type = Param::Type::Variable ) const;

// Import values
// If Ref is given, then we are importing errors, which we add in quadrature
template <typename T>
void Params::Import( Vector<T> &All, const VectorView<T> &vType, Param::Type type, const Vector<T> *pRef ) const
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
          T Previous{ All[p.OffsetAll + i - 1] };
          T Current{ vType[Offset + i] };
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
                                     Param::Type type, const Vector<float> *pRef ) const;
template void Params::Import<double>( Vector<double> &All, const VectorView<double> &vType,
                                      Param::Type type, const Vector<double> *pRef ) const;

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
        if( !pbKnown || (*pbKnown)[Offset + i] )
        {
          os << std::string( MaxLen - p.FieldLen + 2, ' ' ) << it.first;
          if( p.size )
            os << i;
          os << " " << Values[Offset + i];
          if( pErrors && p.type == Param::Type::Variable )
            os << "\t+/- " << (*pErrors)[p.OffsetMyType + i];
          os << "\n";
        }
      }
    }
  }
}

template void Params::Dump<float>( std::ostream &os, const Vector<float> &Values, Param::Type ShowType,
                                   const Vector<float> *pErrors, const std::vector<bool> *pbKnown ) const;
template void Params::Dump<double>( std::ostream &os, const Vector<double> &Values, Param::Type ShowType,
                                    const Vector<double> *pErrors, const std::vector<bool> *pbKnown ) const;

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
    ::H5::Attribute a{ g.openAttribute( sObjectNames ) };
    std::vector<std::string> vObjects{ H5::ReadStrings( a ) };
    a.close();
    // Get number of params
    int iNumParams;
    a = g.openAttribute( sCount );
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
      {
        std::vector<int> vOID;
        a = gSub.openAttribute( sObjectNames );
        H5::ReadVectorHelper( vOID, a, sObjectNames );
        k.Object.reserve( vOID.size() );
        for( int i : vOID )
          k.Object.push_back( vObjects[i] );
        a.close();
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
      std::vector<int> vOIDs( k.Object.size() );
      for( std::size_t j = 0; j < k.Object.size(); ++j )
        vOIDs[j] = OID[k.Object[j]];
      Dims[0] = vOIDs.size();
      {
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

const std::string Params::sCount{ "Count" };
const std::string Params::sObjectNames{ "Object" };
const std::string Params::sName{ "Name" };
const std::string Params::sTypeName{ "Type" };
const std::string Params::sSize{ "Size" };
const std::string Params::sMonotonic{ "Monotonic" };

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

MLU_Param_hpp_end