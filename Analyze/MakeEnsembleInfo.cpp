/*************************************************************************************
 
 Make Ensemble info
 
 Source file: MakeEnsembleInfo.cpp
 
 Copyright (C) 2023
 
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
#include <MLU/Common.hpp>

#include <stdio.h>
#include <ostream>

using SeedType = std::uint_fast32_t;

using Scalar = double;
using Model = Common::Model<Scalar>;

// A value with a standard deviation
struct ValStddev
{
  Scalar Value;
  Scalar Stddev;
  ValStddev Average( const ValStddev &o ) const;
  ValStddev AddQuadrature( unsigned int lMul, const ValStddev &r, unsigned int rMul ) const;
};

ValStddev ValStddev::Average( const ValStddev &o ) const
{
  ValStddev Avg;
  Avg.Value = 0.5 * ( Value + o.Value );
  if( Stddev == o.Stddev )
    Avg.Stddev = Stddev;
  else
  {
    // Add in quadrature
    const Scalar A{   Stddev /   Value };
    const Scalar B{ o.Stddev / o.Value };
    Avg.Stddev = std::sqrt( A*A + B*B ) * Avg.Value;
  }
  return Avg;
}

ValStddev ValStddev::AddQuadrature( unsigned int lMul, const ValStddev &r, unsigned int rMul ) const
{
  if( lMul == rMul && Stddev == r.Stddev )
    return Average( r );
  const Scalar Norm{ static_cast<Scalar>( 1. / ( lMul + rMul ) ) };
  ValStddev v;
  v.Value = Norm * ( lMul * Value + rMul * r.Value );
  // Add errors in quadrature
  const Scalar L{   Stddev /   Value };
  const Scalar R{ r.Stddev / r.Value };
  v.Stddev = std::sqrt( Norm * ( lMul * L * L + rMul * R * R ) ) * v.Value;
  return v;
}

std::ostream &operator<<( std::ostream &os, const ValStddev &v )
{
  return os << v.Value << "\t[" << ( v.Value - v.Stddev ) << ",\t" << ( v.Value + v.Stddev ) << "]";
}

struct GlobalInfo
{
  std::string Name;
  SeedType  Seed;
  ValStddev Value;
};

// PDG meson masses from https://pdglive.lbl.gov Units eV
const ValStddev DStar2010{ 2.01026e9, 0.05e6 }; // D^*(2010)^\pm
const ValStddev DStar2007{ 2.00685e9, 0.05e6 }; // D^*(2007)^0
const ValStddev PDGDPM{ 1.86966e9, 0.05e6 }; // D^\pm
const ValStddev PDGD0{ 1.86484e9, 0.05e6 }; // D^0
const ValStddev PDGKPM{ 493.677e6, 0.016e6 }; // K^\pm
const ValStddev PDGK0{ 497.611e6, 0.013e6 }; // K^0
const ValStddev PDGPIPM{ 139.57039e6, 0.00018e6 }; // \pi^\pm
const ValStddev PDGPI0{ 134.9768e6, 0.0005e6 }; // \pi^0

const std::array<GlobalInfo, 9> globalInfo{{
  // Gradient flow scale
  // Data from https://arxiv.org/pdf/1411.7017.pdf \cite{Blum:2016} table 1, page 6
  { "w0", 4195934094, { 0.8742e-9, 4.6e-12 } }, // Units eV^{-1}
  { "fPI", 2055826088, {130.19e6, 0.89e6} }, // Units eV
  { "fK", 528365601, {155.51e6, 0.83e6} }, // Units eV
  // PDG meson masses from https://pdglive.lbl.gov Units eV
  { "PDGD0*", 1270112283, { 2.343e9, 10e6 } },
  { "PDGD*", 1358725311, DStar2010.Average( DStar2007 ) },
  { "PDGD", 1081938155, PDGDPM.Average( PDGD0 ) },
  { "PDGDs", 1794043036, { 1.96835e9, 0.07e6 } }, // PDG D_s^\pm
  { "PDGK", 3899441903, PDGKPM.Average( PDGK0 ) },
  { "PDGPI", 4013980664, PDGPIPM.AddQuadrature( 2, PDGPI0, 1 ) },
}};

struct EnsembleInfo
{
  std::string Ensemble;
  SeedType    Seed;
  int         L;
  int         T;
  ValStddev   aInv;
};

static constexpr int NumEnsembles{ 6 };

// Data from https://arxiv.org/pdf/1812.08791.pdf \cite{Boyle:2018knm}
// Inverse lattice spacings in eV
static const std::array<EnsembleInfo, NumEnsembles> EnsembleArray{{
  { "C1",  2124653957, 24, 64, { 1.7848e9, 5e6 } },
  { "C2",  2924767711, 24, 64, { 1.7848e9, 5e6 } },
  { "F1M", 3490653410, 48, 96, { 2.7080e9, 1e7 } },
  { "M1",  3144081093, 32, 64, { 2.3833e9, 8.6e6 } },
  { "M2",   277434998, 32, 64, { 2.3833e9, 8.6e6 } },
  { "M3",   492385535, 32, 64, { 2.3833e9, 8.6e6 } }
}};

static constexpr int NumParams{ 1 };

static const std::array<std::string, NumParams> ParamNames{{ "aInv" }};

inline SeedType RandomNumber()
{
  return std::random_device{}(); // Hardware generated
}

void MakeSeeds()
{
  std::random_device rd; // Hardware generated
  for( const EnsembleInfo &ei : EnsembleArray )
    std::cout << ei.Ensemble << '\t' << rd() << Common::NewLine;
  std::cout << "w0" << '\t' << rd() << Common::NewLine;
}

class Maker
{
protected:
  const Common::Params params;
  Common::Params MakeParams();
  void MakeGaussian( Model &m, Common::Param::Key k, const ValStddev &v, SeedType Seed ) const;
  void MakeEnsembleInfo( std::string sFileName ) const;
public:
  Maker() : params{ MakeParams() } {}
  void Run( std::string sFileName ) const;
};

Common::Params Maker::MakeParams()
{
  Common::Params params;
  Common::Param::Key k;
  for( const GlobalInfo & gi : globalInfo )
  {
    k.Name = gi.Name;
    params.Add( k );
  }
  k.Object.resize( 1 );
  for( const EnsembleInfo &ei : EnsembleArray )
  {
    k.Object[0] = ei.Ensemble;
    for( const std::string &ParamName : ParamNames )
    {
      k.Name = ParamName;
      params.Add( k );
    }
  }
  params.AssignOffsets();
  return params;
}

void Maker::MakeGaussian( Model &m, Common::Param::Key k, const ValStddev &v, SeedType Seed ) const
{
  if( Seed == 0 )
    Seed = RandomNumber();
  const std::size_t Column{ params.Find( k, "MakeGaussian() Bug" )->second() };
  std::cout << k << '\t' << v << '\t' << Seed << Common::NewLine;
  std::mt19937                     engine( Seed );
  std::normal_distribution<Scalar> random( v.Value, v.Stddev );
  m( Model::idxCentral, Column ) = v.Value;
  for( std::size_t i = 0; i < m.NumSamples(); ++i )
    m( i, Column ) = random( engine );
}

void Maker::Run( std::string sFileName ) const
{
  std::cout << "Making " << sFileName << Common::NewLine;
  Model m{ static_cast<int>( Common::RandomCache::DefaultNumReplicas() ), params, {} };
  Common::Param::Key k;
  for( const GlobalInfo & gi : globalInfo )
  {
    k.Name = gi.Name;
    MakeGaussian( m, k, gi.Value, gi.Seed );
  }
  k.Name = ParamNames[0];
  k.Object.resize( 1 );
  for( const EnsembleInfo &ei : EnsembleArray )
  {
    k.Object[0] = ei.Ensemble;
    MakeGaussian( m, k, ei.aInv, ei.Seed );
  }
  m.SetSummaryNames( "Central" );
  m.MakeCorrSummary();
  Common::MakeAncestorDirs( sFileName );
  m.Write( sFileName );
}

int main(int argc, char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  //MakeSeeds();
  int iReturn = EXIT_SUCCESS;
  std::string sFileName{ argc < 2 ? "EnsembleInfo.h5" : argv[1] };
  try
  {
    if( Common::FileExists( sFileName ) )
      throw std::runtime_error( sFileName + " exists" );
    Maker{}.Run( sFileName );
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  return iReturn;
}
