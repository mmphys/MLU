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
  enum class Type{ All, Variable, Fixed };

  std::size_t size;   // How many values (e.g. energy levels) this parameter has
  bool bMonotonic;    // strictly increasing - implemented as p_n = p_{n-1} + a^2
  Type type;

  struct Key
  {
    Key() = default;
    std::vector<std::string> Object; // Name of the object these parameters describe, e.g. D_s meson, p^2=1
    std::string Name;   // Name of the parameter
    std::size_t Len() const;
    std::string FullName( std::size_t idx=0, std::size_t Size=std::numeric_limits<std::size_t>::max() ) const;
    bool operator==( const Key &rhs ) const;
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

struct Params : std::map<Param::Key, Param, Param::Key::Less>
{
  //using MapT = std::map<Param::Key, Param, Param::Key::Less>;
  //using MapPairT = std::pair<iterator, bool>;
  //iterator Add( const Param::Key &key, std::size_t Size = 1, bool bMonotonic = false,
                //Param::Type Type = Param::Type::Variable );
  Params::iterator Add( const Param::Key &key, std::size_t NumExp = 1, bool bMonotonic = false,
                        Param::Type Type = Param::Type::Variable );
  void AssignOffsets();
  std::size_t NumScalars( Param::Type Type ) const
  {
    return Type == Param::Type::All      ? NumVariable + NumFixed
         : Type == Param::Type::Variable ? NumVariable : NumFixed;
  }
  template <typename T>
  void Export( Vector<T> &vType, const Vector<T> &All, Param::Type type = Param::Type::Variable ) const;
  template <typename T>
  void Import( Vector<T> &All, const VectorView<T> &vType, Param::Type type = Param::Type::Variable,
               const Vector<T> * pRef = nullptr ) const;
  template <typename T>
  void Dump( std::ostream &os, const Vector<T> &Values, Param::Type type = Param::Type::Variable,
             const Vector<T> *pErrors = nullptr, const std::vector<bool> *pbKnown = nullptr ) const;
  void ReadH5 ( ::H5::Group gParent, const std::string GroupName );
  void WriteH5( ::H5::Group gParent, const std::string GroupName ) const;
protected:
  std::size_t NumFixed;
  std::size_t NumVariable;
  std::size_t MaxLen;
  static const std::string sCount;
  static const std::string sObjectNames;
  static const std::string sName;
  static const std::string sTypeName;
  static const std::string sSize;
  static const std::string sMonotonic;
};

std::ostream &operator<<( std::ostream &os, const Params::iterator &it );

MLU_Param_hpp_end
#endif
