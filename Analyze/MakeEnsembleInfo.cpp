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
};

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
// Inverse lattice spacings in GeV
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

static std::array<std::array<std::size_t, NumParams>, NumEnsembles> ParamIndex;

void MakeSeeds()
{
  std::random_device rd; // Hardware generated
  for( const EnsembleInfo &ei : EnsembleArray )
  {
    std::cout << ei.Ensemble << Common::Space << rd() << Common::NewLine;
  }
}

Common::Params MakeParams()
{
  Common::Params params;
  Common::Param::Key k;
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
  for( std::size_t i = 0; i < EnsembleArray.size(); ++i )
  {
    k.Object[0] = EnsembleArray[i].Ensemble;
    for( std::size_t j = 0; j < ParamNames.size(); ++j )
    {
      k.Name = ParamNames[j];
      Common::Params::iterator it{  };
      ParamIndex[i][j] = params.Find( k, "Bug looking for key I just inserted" )->second();
    }
  }
  return params;
}

void MakeGaussianAInv( Model &m, std::size_t idxEnsemble )
{
  std::size_t Column{ ParamIndex[idxEnsemble][0] };
  const EnsembleInfo &ei{ EnsembleArray[idxEnsemble] };
  std::cout << ei.Ensemble << '\t'
            << ei.aInv.Value << "\t["
            << ( ei.aInv.Value - ei.aInv.Stddev ) << ",\t"
            << ( ei.aInv.Value + ei.aInv.Stddev ) << "]\n";
  std::mt19937                     engine( ei.Seed );
  std::normal_distribution<Scalar> random( ei.aInv.Value, ei.aInv.Stddev );
  m( Model::idxCentral, Column ) = ei.aInv.Value;
  for( std::size_t i = 0; i < m.NumSamples(); ++i )
    m( i, Column ) = random( engine );
}

void MakeEnsembleInfo( std::string sFileName )
{
  std::cout << "Making " << sFileName << Common::NewLine;
  Common::Params params{ MakeParams() };
  Model m{ static_cast<int>( Common::RandomCache::DefaultNumReplicas() ), params, {} };
  for( std::size_t i = 0; i < EnsembleArray.size(); ++i )
    MakeGaussianAInv( m, i );
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
    MakeEnsembleInfo( sFileName );
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
