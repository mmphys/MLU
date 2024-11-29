/**
 
 Make Ensemble info
 
 Source file: MakeEnsembleInfo.hpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLU/MLU.hpp>

#include <stdio.h>
#include <ostream>

using SeedType = std::uint_fast32_t;

using Scalar = double;
using Model = MLU::Model<Scalar>;

// A value with a standard deviation
struct ValStddev
{
  Scalar Value;
  Scalar Stddev;
  ValStddev Average( const ValStddev &o ) const;
  ValStddev AddQuadrature( unsigned int lMul, const ValStddev &r, unsigned int rMul ) const;
  inline std::string to_string( unsigned char SigFigError ) const
  {
    return MLU::ValSigFig<Scalar>::Show( Value, Stddev, SigFigError );
  }
};

std::ostream &operator<<( std::ostream &os, const ValStddev &v );

struct GlobalInfo
{
  std::string Name;
  SeedType  Seed;
  ValStddev Value;
};

static constexpr int MaxNumParams{ 3 };
static constexpr int idxaInv{ 0 };
static constexpr int idxmPi{ 1 };
static constexpr int idxZV{ 2 };

struct EnsembleInfo
{
  struct ValT : public ValStddev
  {
    SeedType Seed;
    ValT( Scalar value, Scalar stddev, SeedType seed ) : ValStddev{value, stddev}, Seed{seed} {}
  };
  std::string Ensemble;
  int         L;
  int         T;
  std::array<ValT, MaxNumParams> Value;
};

static constexpr int NumEnsembles{ 6 };

inline SeedType RandomNumber()
{
  return std::random_device{}(); // Hardware generated
}

class Maker
{
protected:
  using EnsMPiMapT = std::map<std::string, std::string, MLU::LessCaseInsensitive>;
  using EnsMPiReaderT = MLU::KeyValReader<std::string, std::string, MLU::LessCaseInsensitive>;
  const MLU::Params params;
  const EnsMPiMapT MPiMap;
  MLU::Params MakeParams();
  EnsMPiMapT MakeMPiMap( const char *mPiList );
  void MakeGaussian( Model &m, MLU::Param::Key k, const ValStddev &v, SeedType Seed ) const;
  //void MakeEnsembleInfo( std::string sFileName ) const;
  // Read Raj's Z_{V,mixed} data
  void ReadZV();
public:
  Maker( const char *mPiList ) : params{MakeParams()}, MPiMap{MakeMPiMap( mPiList )} { ReadZV(); }
  void Run( std::string sFileName ) const;
protected:
  int NumParams;
};

struct ZVInfo
{
  std::vector<Scalar> ZV;
  std::vector<Scalar> Err;
  std::vector<Scalar> Mu;
  void Hydrate( ::H5::Group &g );
protected:
  void Hydrate( ::H5::Group &g, std::vector<Scalar> &v, const std::string &sName );
};
