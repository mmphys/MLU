/*************************************************************************************

Comparison of Grid serialisation of Eigen tensor vs NamedTensor

Copyright (C) 2021

Author: Michael Marshall <michael.marshall@ed.ac.uk>

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

#include <stdio.h>
//#include <MLU/Common.hpp>
#include <omp.h>
#include <Grid/Grid.h>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// The only way I could get these iterators to work is to put the begin() and end() functions in the Eigen namespace
// So if Eigen ever defines these, we'll have a conflict and have to change this
namespace Eigen {
  template <typename ET>
  inline typename std::enable_if<Grid::EigenIO::is_tensor<ET>::value, typename Grid::EigenIO::Traits<ET>::scalar_type *>::type
  begin( ET & et ) { return reinterpret_cast<typename Grid::EigenIO::Traits<ET>::scalar_type *>(et.data()); }
  template <typename ET>
  inline typename std::enable_if<Grid::EigenIO::is_tensor<ET>::value, typename Grid::EigenIO::Traits<ET>::scalar_type *>::type
  end( ET & et ) { return begin(et) + et.size() * Grid::EigenIO::Traits<ET>::count; }
}

using ESO = Eigen::StorageOptions;
using TensorRank4 = Eigen::Tensor<Complex, 4, ESO::RowMajor>;

class NoiseTensorReadable : Grid::Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(NoiseTensorReadable,
                                  TensorRank4,                tensor,
                                  std::vector<std::string>,   IndexNames,
                                  NamedTensorDefaultMetadata, MetaData );
};

// Demonstrate the flexibility of Eigen serialisation cf NamedTensor
bool Debug()
{
  static const std::string FileName{ "TensorDemo.h5" };
  static constexpr int iLong{ 7 };
  static constexpr int iShort{ 2 };
  using ESO = Eigen::StorageOptions;
  using TensorRank2 = Eigen::Tensor<int, 2, ESO::RowMajor>;
  using NoiseTensor = Grid::Hadrons::MDistil::NoiseTensor;
  int Counter {};
  TensorRank2 t2Tall( iLong, iShort );
  TensorRank2 t2Wide( iShort, iLong );
  NoiseTensor tNoise( 2, 3, 5, 4 ); // nNoise, nT, nVec, nS
  tNoise.MetaData.Version = 777;
  for( auto &s : t2Tall ) s = Counter++;
  for( auto &s : t2Wide ) s = Counter++;
  for( auto &s : tNoise.tensor ) s = Counter++;
  const std::string st2Tall{ "Tall2d" };
  const std::string st2Wide{ "Wide2d" };
  const std::string st2Noise{ "Noise4d" };
  std::cout << st2Tall << ": " << t2Tall << std::endl;
  std::cout << st2Wide << ": " << t2Wide << std::endl;
  {
    std::cout << "Writing " << FileName << std::endl;
    Grid::Hdf5Writer w( FileName );
    write( w, st2Tall, t2Tall );
    write( w, st2Wide, t2Wide );
    write( w, st2Noise, tNoise );
  }
  std::cout << "Reading back " << FileName << " exchanging tall/wide" << std::endl;
  Grid::Hdf5Reader r( FileName );
  NoiseTensor tGoodNoise( 2, 3, 5, 4 ); // nNoise, nT, nVec, nS
  NoiseTensor tBadNoise( 2, 3, 4, 4 ); // nNoise, nT, nVec, nS
  read( r, st2Wide, t2Tall );
  std::cout << st2Tall << ": " << t2Tall << std::endl;
  read( r, st2Tall, t2Wide );
  std::cout << st2Wide << ": " << t2Wide << std::endl;
  NoiseTensorReadable rNoise;
  std::cout << "rNoise (before read): IndexNames " << rNoise.IndexNames << ", MetaData.Version=" << rNoise.MetaData.Version << ", tensor:" << rNoise.tensor << std::endl;
  read( r, st2Noise, rNoise );
  std::cout << "rNoise: IndexNames " << rNoise.IndexNames << ", MetaData.Version=" << rNoise.MetaData.Version << ", tensor:" << rNoise.tensor << std::endl;
  read( r, st2Noise, tGoodNoise );
  std::cout << "tGoodNoise: IndexNames " << tGoodNoise.IndexNames << ", MetaData.Version=" << tGoodNoise.MetaData.Version << ", tensor:" << tGoodNoise.tensor << std::endl;
  read( r, st2Noise, tBadNoise );
  std::cout << "tBadNoise: IndexNames " << tBadNoise.IndexNames << ", MetaData.Version=" << tBadNoise.MetaData.Version << ", tensor:" << tBadNoise.tensor << std::endl;
  return true;
}

int main(int argc, char *argv[])
{
  //std::ios_base::sync_with_stdio( false );
  Grid::Grid_init(&argc,&argv);
  //std::cout << Grid::GridLogMessage << MLUVersionInfoHuman() << std::endl;
  int iReturn = EXIT_SUCCESS;
  try
  {
    //#pragma omp parallel
    if( !Debug() )
      iReturn = EXIT_FAILURE;
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  Grid::Grid_finalize();
  return iReturn;
}
