//
//  Debug.cpp
//  A place to play with debug code
//
//  Created by Michael Marshall on 09/05/2020.
//  Copyright © 2020 sopa. All rights reserved.
//

#include <stdio.h>
//#include <MLU/Common.hpp>
#include <omp.h>
#include <Grid/Grid.h>

template<typename T>
std::ostream & operator<<( std::ostream &o, const std::vector<T> &vT )
{
  static constexpr int ShowNum{ 5 };
  std::cout << "[";
  for( std::size_t i = 0; i < vT.size(); ++i )
  {
    if( i < ShowNum || i > vT.size() - ShowNum - 1 )
      std::cout << " " << vT[i];
    else if( i == ShowNum )
      std::cout << " ...";
  }
  std::cout << " ]";
  return o;
}

bool Debug()
{
  static constexpr int RaggedRows{ 8 };
  static constexpr int MinCols{ 1 };
  static constexpr int FixedCols{ MinCols + 1 };
  static constexpr int FixedRows{ 1000 };
  // Build the ragged vector
  std::vector<std::vector<int>> vRagged{ RaggedRows };
  int NumCols{ MinCols + RaggedRows - 1 };
  int Counter = 0;
  for( int r = 0; r < RaggedRows; ++r, --NumCols )
  {
    //NumCols = r % 3 + 1;
    vRagged[r].resize( NumCols );
    for( int c = 0; c < NumCols; ++c )
      vRagged[r][c] = Counter++;
  }
  // Build the uniform vector
  std::vector<std::vector<int>> vUniform{ FixedRows };
  Counter = 0;
  for( int r = 0; r < FixedRows; ++r )
  {
    vUniform[r].resize( FixedCols );
    for( int c = 0; c < FixedCols; ++c )
      vUniform[r][c] = Counter++;
  }
  std::cout << "Uniform vector:" << vUniform << std::endl;
  std::cout << "Ragged  vector:" << vRagged << std::endl;
#ifdef HAVE_HDF5
  const std::string FileName{ "VectorDebug.h5" };
  const std::string tagU{ "UniformVector" };
  const std::string tagR{ "RaggedVector" };
  {
    Grid::Hdf5Writer writer(FileName);
    write(writer, tagU , vUniform);
    write(writer, tagR , vRagged);
  }
  std::cout << "Reading back ..." << std::endl;
  Grid::Hdf5Reader reader(FileName);
  std::vector<std::vector<int>> vReadU;
  std::vector<std::vector<int>> vReadR;
  read( reader, tagU, vReadU );
  read( reader, tagR, vReadR );
  std::cout << "Uniform vector:" << vReadU << std::endl;
  std::cout << "Ragged  vector:" << vReadR << std::endl;
#endif
  return true;
}

int main(int argc, char *argv[])
{
  //std::ios_base::sync_with_stdio( false );
  Grid::Grid_init(&argc,&argv);
  //std::cout << Grid::GridLogMessage << MLUVersionInfoHuman() << std::endl;
  std::cout << Grid::GridLogMessage << "Hello" << std::endl;
  int iReturn = EXIT_SUCCESS;
  try
  {
    //#pragma omp parallel
    std::cout << "Hello world!\n";
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
