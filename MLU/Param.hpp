/*****************************************************************************
 
 Manage model parameters

 Source file: Param.hpp
 
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

#ifndef MLU_Param_hpp
#define MLU_Param_hpp namespace Common {
#define MLU_Param_hpp_end };

#include <MLU/HDF5.hpp>

MLU_Param_hpp

struct Param
{
  /**
  Parameter types
   
   Variable: These are optimised by the fitter
   
   Fixed: Loaded from previous fit results
   
   Derived: Computed by models, but not optimised by the fitter
   */
  enum class Type{ All, Variable, Fixed, Derived };

  std::size_t size;   // How many values (e.g. energy levels) this parameter has
  bool bMonotonic;    // strictly increasing - implemented as p_n = p_{n-1} + a^2
  bool bSwapSourceSink = false;
  Type type;
  // The second member of a pair points to its mate
  Param * pProductWith = nullptr;

  struct Key
  {
    Key() = default;
    std::vector<std::string> Object; // Name of the object parameters describe, e.g. D_s meson, p^2=1
    std::string Name;   // Name of the parameter
    std::size_t Len() const;
    bool empty() const { return Object.empty() && Name.empty(); }
    std::string FullName( std::size_t idx=0, std::size_t Size=std::numeric_limits<std::size_t>::max() ) const;
    std::string ShortName( std::size_t idx=0, std::size_t Size=std::numeric_limits<std::size_t>::max() ) const;
    bool SameObject( const Key &rhs ) const;
    bool operator==( const Key &rhs ) const { return SameObject( rhs ) && EqualIgnoreCase( Name, rhs.Name ); }
    bool operator!=( const Key &rhs ) const { return !operator==( rhs ); }
    struct Less { bool operator()( const Key &lhs, const Key &rhs ) const; };
    // Constructors
    Key( const std::string  &object, const std::string  &name )
    : Object{{object}}, Name{name} { ValidateKey(); }
    Key( const std::string  &object,       std::string &&name )
    : Object{{object}}, Name{std::move(name)} { ValidateKey(); }
    Key(       std::string &&object, const std::string  &name )
    : Object{{std::move(object)}}, Name{name} { ValidateKey(); }
    Key(       std::string &&object,       std::string &&name )
    : Object{{std::move(object)}}, Name{std::move(name)} { ValidateKey(); }
    Key( const std::vector<std::string>  &object, const std::string  &name )
    : Object{object}, Name{name} { ValidateKey(); }
    Key( const std::vector<std::string>  &object,       std::string &&name )
    : Object{object}, Name{std::move(name)} { ValidateKey(); }
    Key(       std::vector<std::string> &&object, const std::string  &name )
    : Object{std::move(object)}, Name{name} { ValidateKey(); }
    Key(       std::vector<std::string> &&object,       std::string &&name )
    : Object{std::move(object)}, Name{std::move(name)} { ValidateKey(); }
  protected:
    void ValidateKey();
  };

  Param( std::size_t Size, bool bMonotonic_, Type type_ = Type::Variable )
  : size{Size}, bMonotonic{bMonotonic_}, type{type_} {}
  void Validate( Type ListType = Type::All ) const;
  // Underlying implementation for all operator(). Don't call this direct
  std::size_t GetOffset( std::size_t Idx, Type ListType ) const;
  std::size_t operator()( std::size_t Idx = 0, Type ListType = Type::All ) const;
  std::size_t operator()( std::size_t idxSnk, std::size_t idxSrc, Type ListType = Type::All ) const;
  template <typename T>
  T &operator()( std::vector<T> &v, std::size_t Idx = 0, Type ListType = Type::All ) const;

protected:
  friend class Params;
  const Key *pKey = nullptr;// Only valid after Params::AssignOffsets()
  std::size_t OffsetAll;    // Offset into a structure containing all parameters
  std::size_t OffsetMyType; // Offset into a structure containing parameters of selected type only
  std::size_t FieldLen;     // How long is this name including digits
};

std::ostream &operator<<( std::ostream &os, const Param::Type &type );
std::istream &operator>>( std::istream &is, Param::Type &type );
std::ostream &operator<<( std::ostream &os, const Param::Key &key );
std::istream &operator>>( std::istream &is, Param::Key &key );

/**
 List of parameters
 
 - Warning: sizes (i.e. protected members) only valid after AssignOffsets()
 */
struct Params : std::map<Param::Key, Param, Param::Key::Less>
{
  //using MapT = std::map<Param::Key, Param, Param::Key::Less>;
  //using MapPairT = std::pair<iterator, bool>;
  //iterator Add( const Param::Key &key, std::size_t Size = 1, bool bMonotonic = false,
                //Param::Type Type = Param::Type::Variable );
  Params::iterator Add( const Param::Key &key, std::size_t NumExp = 1, bool bMonotonic = false,
                        Param::Type Type = Param::Type::Variable );
  Params::iterator MakeFixed( const Param::Key &key, bool bSwapSourceSink );
  void AssignOffsets();
  std::size_t NumScalars( Param::Type Type ) const
  {
    switch( Type )
    {
      case Param::Type::All:
        return NumVariable + NumFixed + NumDerived;
      case Param::Type::Variable:
        return NumVariable;
      case Param::Type::Fixed:
        return NumFixed;
      case Param::Type::Derived:
        return NumDerived;
    }
    throw std::runtime_error( "Unknown Param::Type " + std::to_string( static_cast<int>( Type ) ) );
  }
  std::size_t MaxExponents() const;
  template <typename T>
  void Export( Vector<T> &vType, const Vector<T> &All, Param::Type type = Param::Type::Variable ) const;
  template <typename T>
  void Import( Vector<T> &All, const VectorView<T> &vType, Param::Type type = Param::Type::Variable,
               const Vector<T> * pRef = nullptr ) const;
  template <typename T>
  void Dump( std::ostream &os, const Vector<T> &Values, Param::Type type = Param::Type::Variable,
             const Vector<T> *pErrors = nullptr, const std::vector<bool> *pbKnown = nullptr ) const;
  std::vector<std::string> GetNames( Param::Type type, bool bLongNames ) const;
  void ReadH5 ( ::H5::Group gParent, const std::string GroupName );
  void WriteH5( ::H5::Group gParent, const std::string GroupName ) const;
  bool SingleObject() const { return bSingleObject; }
  /**
   I might only know the product of some parameters.
   If there's a choice, make the second of each pair negative
   - Parameters:
    - Products: even list of parameter pairs
   */
  void SetProducts( const std::string &sProducts );
protected:
  // Only valid after AssignOffsets()
  bool bSingleObject;
  std::size_t NumFixed;
  std::size_t NumVariable;
  std::size_t NumDerived;
  std::size_t MaxLen;
  static const std::string sCount;
  static const std::string sObjectNames;
  static const std::string sName;
  static const std::string sTypeName;
  static const std::string sSize;
  static const std::string sMonotonic;
  std::size_t &SizeType( Param::Type Type )
  {
    switch( Type )
    {
      case Param::Type::Variable:
        return NumVariable;
      case Param::Type::Fixed:
        return NumFixed;
      case Param::Type::Derived:
        return NumDerived;
      case Param::Type::All:
        break;
    }
    throw std::runtime_error( "Unknown Param::Type " + std::to_string( static_cast<int>( Type ) ) );
  }
};

std::ostream &operator<<( std::ostream &os, const Params::value_type &param );

MLU_Param_hpp_end
#endif
