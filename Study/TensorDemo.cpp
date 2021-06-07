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

//using namespace Grid;
//using namespace Hadrons;
using Grid::operator<<;
using Grid::operator>>;
using Grid::Complex;
using Grid::iMatrix;

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

using   GridScalar = float;
using   GridTensor = iMatrix<iMatrix<GridScalar, 3>, 2>;
using  vGridTensor = std::vector<GridTensor>;
using vvGridTensor = std::vector<vGridTensor>;

class NoiseTensorReadable : Grid::Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(NoiseTensorReadable,
                                  TensorRank4,                tensor,
                                  std::vector<std::string>,   IndexNames,
                                  Grid::Hadrons::NamedTensorDefaultMetadata, MetaData );
};

void ShowNoiseTensor( const std::string &Name, const TensorRank4 &tensor,
                      const std::vector<std::string> &IndexNames,
                      const Grid::Hadrons::NamedTensorDefaultMetadata &md )
{
  std::cout << Name << ": IndexNames " << IndexNames << ", MetaData.Version=" << md.Version;
#if ( (!defined(GRID_SYCL)) && (!defined(GRID_CUDA)) && (!defined(GRID_HIP)) )
  // This causes problems (under CUDA at least) because Grid::Complex is not std::complex
  std::cout << ", noise_tensor:" << tensor;
#endif
  std::cout << std::endl;
}

template <typename T> void ShowEigenTensor( const std::string &Name, const T& t )
{
  auto Num = t.size();
  auto Size = sizeof(typename T::Scalar);
  std::cout << Name << " (size " << Num << " * " << Size << " = " << Num*Size << ") : ";
  std::cout << t;
  std::cout << std::endl;
}

template <typename T> void ShowEigenTensorSize( const std::string &Name, const T& t )
{
  auto Num = t.size();
  auto Size = sizeof(typename T::Scalar);
  std::cout << Name << " (size " << Num << " * " << Size << " = " << Num*Size << ") : ";
  std::cout << std::endl;
}

// Demonstrate the flexibility of Eigen serialisation cf NamedTensor
bool Debug()
{
  static const std::string FileName{ "TensorDemo.h5" };
  // Create tensors
  std::vector<std::string> vString{ {"Do"},{"androids"},{"dream"},{"of"},{"electric"},{"sheep?"} };
  using SingleGridTensor = Grid::iMatrix<Grid::iVector<Grid::iMatrix<Grid::iVector<Grid::Complex, 5>, 4>, 3>, 2>;
  static constexpr int iLong{ 7 };
  static constexpr int iShort{ 2 };
  using ESO = Eigen::StorageOptions;
  using TensorRank2 = Eigen::Tensor<int, 2, ESO::RowMajor>;
  using TensorGrid = Eigen::Tensor<iMatrix<Complex, 2>, 2, ESO::RowMajor>;
  using NoiseTensor = Grid::Hadrons::MDistil::NoiseTensor;
  int Counter {};
  SingleGridTensor GT1;
  TensorRank2 t2Tall( iLong-1, iShort );
  TensorRank2 t2Wide( iShort, iLong+1 );
  TensorGrid  t2Grid( iShort, iLong );
  NoiseTensor tNoise( 2, 3, 5, 4 ); // nNoise, nT, nVec, nS
  tNoise.MetaData.Version = 777;
  vGridTensor tGrid( 4 );
  vvGridTensor tGridReg( 2 );
  vvGridTensor tGridRag( 3 );
  const std::string svString{ "PhilipKindredDick" };
  const std::string sGT1{ "SingleGridT" };
  const std::string st2Tall{ "Tall2d" };
  const std::string st2Wide{ "Wide2d" };
  const std::string st2Grid{ "EigenGrid" };
  const std::string st2Noise{ "Noise4d" };
  const std::string stGrid{ "Grid" };
  const std::string stGridReg{ "GridReg" };
  const std::string stGridRag{ "GridRag" };
  // Fill tensors with something interesting (sequential numbers) so we can check I/O
  for( auto &s : GT1    ) s = Counter++;
  for( auto &s : t2Tall ) s = Counter++;
  for( auto &s : t2Wide ) s = Counter++;
  for( auto &s : t2Grid ) s = Counter++;
  for( auto &s : tNoise.tensor ) s = Counter++;
  for( auto &s : tGrid )
    for( auto &x : s )
      x = Counter++;
  for( auto &v1 : tGridReg )
  {
    v1.resize(3);
    for( auto &s : v1 )
    {
      for( auto &x : s )
        x = Counter++;
    }
  }
  for( std::size_t i = 0; i < tGridRag.size(); ++i )
  {
    tGridRag[i].resize(i + 1);
    for( auto &s : tGridRag[i] )
    {
      for( auto &x : s )
        x = Counter++;
    }
  }
  // Show tensors
  std::cout << sGT1      << ": " << GT1      << std::endl;
  ShowEigenTensor( st2Tall, t2Tall );
  ShowEigenTensor( st2Wide, t2Wide );
  ShowEigenTensorSize( st2Grid, t2Grid );
  std::cout << stGrid  << ": " << tGrid  << std::endl;
  std::cout << stGridReg << ": " << tGridReg << std::endl;
  std::cout << stGridRag << ": " << tGridRag << std::endl;
  // Write tensors
  {
    std::cout << "Writing " << FileName << std::endl;
    Grid::Hdf5Writer w( FileName );
    write( w, svString, vString );
    write( w, sGT1,     GT1 );
    write( w, st2Tall,  t2Tall );
    write( w, st2Wide,  t2Wide );
    write( w, st2Grid,  t2Grid );
    write( w, st2Noise, tNoise );
    write( w, stGrid,   tGrid );
    write( w, stGridReg,tGridReg );
    write( w, stGridRag,tGridRag );
  }
  // Read back tensors - deliberately exchange tall / wide
  std::cout << "Reading back " << FileName << " exchanging tall/wide" << std::endl;
  Grid::Hdf5Reader r( FileName );
  std::vector<std::string> Novel;
  read( r, svString, Novel );
  for( std::size_t i = 0; i < Novel.size(); ++i )
  {
    if( i ) std::cout << " ";
    std::cout << Novel[i];
  }
  std::cout << std::endl;
  SingleGridTensor rGT1;
  read( r, sGT1, rGT1 );
  std::cout << sGT1      << ": " << rGT1      << std::endl;
  NoiseTensor tGoodNoise( 2, 3, 5, 4 ); // nNoise, nT, nVec, nS
  NoiseTensor tBadNoise( 2, 3, 4, 4 ); // nNoise, nT, nVec, nS
  std::cout << "Grid::EigenIO::EigenResizeCounter=" << Grid::EigenIO::EigenResizeCounter << std::endl;
  read( r, st2Wide, t2Tall );
  ShowEigenTensor( st2Tall, t2Tall );
  std::cout << "Grid::EigenIO::EigenResizeCounter=" << Grid::EigenIO::EigenResizeCounter << std::endl;
  read( r, st2Tall, t2Wide );
  ShowEigenTensor( st2Wide, t2Wide );
  std::cout << "Grid::EigenIO::EigenResizeCounter=" << Grid::EigenIO::EigenResizeCounter << std::endl;
  TensorGrid  r2Grid;
  ShowEigenTensorSize( st2Grid + " before read", r2Grid );
  read( r, st2Grid, r2Grid );
  ShowEigenTensorSize( st2Grid, r2Grid );
  std::cout << "Grid::EigenIO::EigenResizeCounter=" << Grid::EigenIO::EigenResizeCounter << std::endl;
  NoiseTensorReadable rNoise;
  ShowNoiseTensor( "rNoise (before read)", rNoise.tensor, rNoise.IndexNames, rNoise.MetaData );
  read( r, st2Noise, rNoise );
  ShowNoiseTensor( "rNoise (after  read)", rNoise.tensor, rNoise.IndexNames, rNoise.MetaData );
  read( r, st2Noise, tGoodNoise );
  ShowNoiseTensor( "rNoise (before read)", tGoodNoise.tensor, tGoodNoise.IndexNames, tGoodNoise.MetaData );
  //read( r, st2Noise, tBadNoise );
  ShowNoiseTensor( "rNoise (before read)", tBadNoise.tensor, tBadNoise.IndexNames, tBadNoise.MetaData );
  vGridTensor rGrid;
  read( r, stGrid, rGrid );
  std::cout << stGrid << ": " << rGrid << std::endl;
  vvGridTensor rGridReg;
  read( r, stGridReg, rGridReg );
  std::cout << stGridReg << ": " << rGridReg << std::endl;
  vvGridTensor rGridRag;
  read( r, stGridRag, rGridRag );
  std::cout << stGridRag << ": " << rGridRag << std::endl;
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
