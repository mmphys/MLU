/*************************************************************************************
 
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
 *************************************************************************************/
/*  END LEGAL */

#ifndef Param_hpp
#define Param_hpp

#include "MultiFit.hpp"

struct ParamBase
{
  enum class Type{ All, Variable, Fixed };
  static void Validate( ParamBase::Type Type, const char * pErrorIfAll = nullptr );

  std::size_t size;   // How many values (e.g. energy levels) this parameter has
  bool bMonotonic;    // strictly increasing - implemented as p_n = p_{n-1} + a^2
  Type type;

  ParamBase( std::size_t Size, bool bMonotonic_, Type type_ = Type::Variable )
  : size{Size}, bMonotonic{bMonotonic_}, type{type_} {}
protected:
  friend class Params;
  std::size_t OffsetAll;    // Offset into a structure containing all parameters
  std::size_t OffsetMyType; // Offset into a structure containing parameters of selected type only
  std::size_t FieldLen;     // How long is this name including digits
};

std::ostream &operator<<( std::ostream &os, const ParamBase::Type &type );
std::istream &operator>>( std::istream &is, ParamBase::Type &type );

struct Param : public ParamBase
{
  using Type = enum ParamBase::Type;

  struct Key
  {
    Key() = default;
    std::string Object; // Name of the object these parameters describe, e.g. D_s meson, p^2=1
    std::string Name;   // Name of the parameter
    std::size_t Len() const { return Object.length() + ( Object.empty() ? 0 : 1 ) + Name.length(); }
    std::string FullName( std::size_t idx=0, std::size_t Size=std::numeric_limits<std::size_t>::max() ) const;
    Key( std::string object, std::string name ) : Object{object}, Name{name} {}
    bool operator==( const Key &rhs ) const;
    bool operator!=( const Key &rhs ) const;
    struct Less { bool operator()( const Key &lhs, const Key &rhs ) const; };
  };

  Param( std::size_t Size, bool bMonotonic_, Type type_ = Type::Variable )
  : ParamBase( Size, bMonotonic_, type_ ) {}
  std::size_t operator()( std::size_t Idx = 0, Param::Type ListType = Param::Type::All ) const;
  template <typename T>
  T &operator()( std::vector<T> &v, std::size_t Idx = 0, Param::Type ListType = Param::Type::All ) const;
};

struct Params : std::map<Param::Key, Param, Param::Key::Less>
{
  //using MapT = std::map<Param::Key, Param, Param::Key::Less>;
  //using MapPairT = std::pair<iterator, bool>;
  iterator Add( const Param::Key &key, std::size_t Size = 1, bool bMonotonic = false,
                Param::Type = Param::Type::Variable );
  void AssignOffsets();
  std::size_t NumScalars( Param::Type Type ) const
  {
    return Type == Param::Type::All      ? NumVariable + NumFixed
         : Type == Param::Type::Variable ? NumVariable : NumFixed;
  }
  void Export( Vector &vType, const Vector &All, Param::Type type = Param::Type::Variable ) const;
  void Import( Vector &All, const VectorView &vType, Param::Type type = Param::Type::Variable ) const;
  void Dump( std::ostream &os, const Vector &Values, Param::Type type = Param::Type::Variable,
             const Vector *pErrors = nullptr ) const;
protected:
  std::size_t NumFixed;
  std::size_t NumVariable;
  std::size_t MaxLen;
};

std::ostream &operator<<( std::ostream &os, const Param::Key &key );
std::ostream &operator<<( std::ostream &os, const Params::iterator &it );

/*Param::Flags operator~( Param::Flags f );
Param::Flags operator|( Param::Flags l, Param::Flags r );
Param::Flags operator&( Param::Flags l, Param::Flags r );
Param::Flags &operator|=( Param::Flags &l, const Param::Flags &r );
Param::Flags &operator&=( Param::Flags &l, const Param::Flags &r );*/

struct ModelParam
{
  ModelParam  () = default;
  Param::Key Key;
  std::size_t idx;
  typename Params::const_iterator it;
};

#endif // Param_hpp
