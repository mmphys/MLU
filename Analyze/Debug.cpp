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
  std::cout << "[";
  for( const T& t : vT )
    std::cout << " " << t;
  std::cout << " ]";
  return o;
}

bool Debug()
{
  static constexpr int NumRows{ 5 };
  static constexpr int MinCols{ 2 };
  std::vector<std::vector<int>> vUniform{ NumRows };
  std::vector<std::vector<int>> vRagged{ NumRows };
  int NumCols{ MinCols + NumRows - 1 };
  int Increment{ 0 };
  for( int r = 0; r < NumRows; ++r, --NumCols )
  {
    vUniform[r].resize( MinCols );
    for( int c = 0; c < MinCols; ++c )
      vUniform[r][c] = Increment++;
    vRagged[r].resize( NumCols );
    for( int c = 0; c < NumCols; ++c )
      vRagged[r][c] = c;
  }
  std::cout << "Uniform vector:" << vUniform << std::endl;
  std::cout << "Ragged  vector:" << vRagged << std::endl;
#ifdef HAVE_HDF5
  Grid::Hdf5Writer writer("VectorDebug.h5");
  write(writer, "UniformVector" , vUniform);
  write(writer, "RaggedVector" , vRagged);
#endif
  return true;
}

int main(int argc, char *argv[])
{
  //std::ios_base::sync_with_stdio( false );
  Grid::Grid_init(&argc,&argv);
  //std::cout << Grid::GridLogMessage << MLUVersionInfoHuman() << std::endl;
  std::cout << Grid::GridLogMessage << "Hello" << std::endl
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
