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

struct Parameters
{
  struct Parameter
  {
    std::string Name;
    scalar      Value;
    scalar      Error;
    Parameter() {} // Random Value and Error is ok - would like this to show up
    Parameter( const std::string &name_, scalar value_, scalar error_ ) : Name{name_}, Value{value_}, Error{error_} {}
  };
  using iterator = std::vector<Parameter>::iterator;
  using const_iterator = std::vector<Parameter>::const_iterator;
protected:
  std::size_t maxLen = 0;
  std::vector<Parameter> Params;
public:
  // Helpers
  inline void clear() { maxLen = 0; Params.clear(); }
  inline std::size_t size() const { return Params.size(); }
  inline std::size_t MaxLen() const { return maxLen; }
  //inline       Parameter & operator[]( std::size_t i )       { return Params[i]; }
  inline const Parameter & operator[]( std::size_t i ) const { return Params[i]; }
  inline void Add( const std::string &Name, scalar Value, scalar Error )
  {
    Params.emplace_back( Name, Value, Error );
    std::size_t Len = Name.length();
    if( maxLen < Len )
      maxLen = Len;
  }
  inline iterator begin() { return Params.begin(); }
  inline iterator end()   { return Params.end(); }
  inline const_iterator begin() const { return Params.begin(); }
  inline const_iterator end()   const { return Params.end(); }
};

std::ostream & operator<<( std::ostream &os, const Parameters &Params );

struct ParamState
{
  Parameters parameters;
  bool bValid;
  scalar TestStat;
  unsigned int NumCalls;
        scalar getTestStat() const { return bValid ? TestStat : 0; }
  unsigned int getNumCalls() const { return bValid ? NumCalls : 0; }
  ParamState( const Parameters &parameters_, bool bValid_=false, scalar TestStat_=0, unsigned int NumCalls_=0 )
  : parameters{parameters_}, bValid{bValid_}, TestStat{TestStat_}, NumCalls{NumCalls_} {}
  ParamState( Parameters &&parameters_, bool bValid_=false, scalar TestStat_=0, unsigned int NumCalls_=0 )
  : parameters{std::move( parameters_ )}, bValid{bValid_}, TestStat{TestStat_}, NumCalls{NumCalls_} {}
  virtual ~ParamState(){}
  virtual void StandardOut( std::ostream &os ) const = 0; //TODO: Deuglify
  virtual void ReplicaMessage( std::ostream &os ) const = 0; //TODO: Deuglify
};

std::ostream & operator<<( std::ostream &os, const ParamState &State );

#endif // Param_hpp
