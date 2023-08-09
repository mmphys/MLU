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

const std::array<GlobalInfo, 11> globalInfo{{
  // Gradient flow scale
  // Data from https://arxiv.org/pdf/1411.7017.pdf \cite{Blum:2016} table 1, page 6
  { "w0", 4195934094, { 0.8742e-9, 4.6e-12 } }, // Units eV^{-1}
  { "fPi", 2055826088, {130.19e6, 0.89e6} }, // Units eV
  { "fK", 528365601, {155.51e6, 0.83e6} }, // Units eV
  // PDG meson masses from https://pdglive.lbl.gov Units eV
  { "PDGD0Star", 1270112283, { 2.343e9, 10e6 } },
  { "PDGDStar", 1358725311, DStar2010.Average( DStar2007 ) },
  { "PDGD", 1081938155, PDGDPM.Average( PDGD0 ) },
  { "PDGDs", 1794043036, { 1.96835e9, 0.07e6 } }, // PDG D_s^\pm
  { "PDGDsStar", 1247053099, { 2.1122e9, 0.4e6 } }, // PDG D_s^{*\pm}
  { "PDGDs0Star", 4246507467, { 2.3178e9, 0.5e6 } }, // PDG D_{s0}^*(2317)^\pm
  { "PDGK", 3899441903, PDGKPM.Average( PDGK0 ) },
  { "PDGPi", 4013980664, PDGPIPM.AddQuadrature( 2, PDGPI0, 1 ) },
}};

static constexpr int NumParams{ 2 };
static constexpr int idxaInv{ 0 };
static constexpr int idxmPi{ 1 };
static const std::array<std::string, NumParams> ParamNames{{ "aInv", "mPi" }};

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
  std::array<ValT, NumParams> Value;
};

static constexpr int NumEnsembles{ 6 };

// Data from https://arxiv.org/pdf/1812.08791.pdf \cite{Boyle:2018knm}, Table 1, pg 6
// Inverse lattice spacings in eV
static const std::array<EnsembleInfo, NumEnsembles> EnsembleArray{{
  { "C1",  24, 64, {{ { 1.7848e9, 5e6,  2124653957 }, { 339.76e6, 1.22e6, 3524755446 } }} },
  { "C2",  24, 64, {{ { 1.7848e9, 5e6,  2924767711 }, { 430.63e6, 1.38e6, 2427481952 } }} },
  { "F1M", 48, 96, {{ { 2.7080e9, 1e7,  3490653410 }, { 232.01e6, 1.01e6, 1818009048 } }} },
  { "M1",  32, 64, {{ { 2.3833e9, 8.6e6,3144081093 }, { 303.56e6, 1.38e6,  233350358 } }} },
  { "M2",  32, 64, {{ { 2.3833e9, 8.6e6, 277434998 }, { 360.71e6, 1.58e6,  184399682 } }} },
  { "M3",  32, 64, {{ { 2.3833e9, 8.6e6, 492385535 }, { 410.76e6, 1.74e6,  585791459 } }} },
}};

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
  using EnsMPiMapT = std::map<std::string, std::string, Common::LessCaseInsensitive>;
  using EnsMPiReaderT = Common::KeyValReader<std::string, std::string, Common::LessCaseInsensitive>;
  const Common::Params params;
  const EnsMPiMapT MPiMap;
  Common::Params MakeParams();
  EnsMPiMapT MakeMPiMap( const char *mPiList );
  void MakeGaussian( Model &m, Common::Param::Key k, const ValStddev &v, SeedType Seed ) const;
  void MakeEnsembleInfo( std::string sFileName ) const;
public:
  Maker( const char *mPiList ) : params{MakeParams()}, MPiMap{MakeMPiMap( mPiList )} {}
  void Run( std::string sFileName ) const;
};

Maker::EnsMPiMapT Maker::MakeMPiMap( const char *mPiList )
{
  EnsMPiMapT MPiMap;
  if( !mPiList || !*mPiList )
    return EnsMPiMapT();

  std::ifstream s( mPiList );
  if( !Common::FileExists( mPiList ) || s.bad() )
    throw std::runtime_error( std::string( "Error reading m_pi list \"" ) + mPiList + "\"" );
  return EnsMPiReaderT::Read( s );
}

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
  k.Object.resize( 1 );
  for( const EnsembleInfo &ei : EnsembleArray )
  {
    k.Object[0] = ei.Ensemble;
    for( int i = 0; i < NumParams; ++i )
    {
      k.Name = ParamNames[i];
      if( MPiMap.empty() || i != idxmPi )
        MakeGaussian( m, k, ei.Value[i], ei.Value[i].Seed );
      else
      {
        static const std::string sErrorMsg{ "Loading m_pi" };
        // Get my offset and aInv
        const std::size_t cDst{ params.Find( k, sErrorMsg )->second() };
        k.Name = ParamNames[idxaInv];
        const std::size_t caInv{ params.Find( k, sErrorMsg )->second() };
        // Load a mPi from a file
        EnsMPiMapT::const_iterator it = MPiMap.find( ei.Ensemble );
        if( it == MPiMap.cend() )
          throw std::runtime_error( "No m_pi model file for ensemble " + ei.Ensemble );
        Model mPi;
        mPi.Read( it->second, "  " );
        m.FileList.push_back( it->second );
        Common::Param::Key eKey( "l_l_p2_0", Common::ModelBase::EnergyPrefix );
        Common::Params::const_iterator pit{ mPi.params.FindPromiscuous( eKey ) };
        if( pit == mPi.params.cend() )
        {
          std::ostringstream os;
          os << "Ensemble " << ei.Ensemble << " key " << eKey << " not found in " << it->second;
          throw std::runtime_error( os.str().c_str() );
        }
        const std::size_t cE{ pit->second() };
        // save a mPi * a^{-1} = mPi
        for( int i = Model::idxCentral; i < m.NumSamples(); ++i )
          m( i, cDst ) = mPi( i, cE ) * m( i, caInv );
      }
    }
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
    Maker( argc < 3 ? nullptr : argv[2] ).Run( sFileName );
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
