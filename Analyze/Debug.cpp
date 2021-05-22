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

using  vvi = std::vector<std::vector<int>>;
using vvvi = std::vector<vvi>;

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

static constexpr int MinCols{ 1 };

// Fill a 2-d vector with a ragged shape
void FillRagged( vvi &v, int &Counter )
{
  const std::size_t Rows{ v.size() };
  std::size_t Cols{ MinCols + Rows - 1 };
  for( std::size_t r = 0; r < Rows; ++r, --Cols )
  {
    if( Rows > 10 )
      Cols = r % 3 + MinCols;
    v[r].resize( Cols );
    for( int c = 0; c < Cols; ++c )
      v[r][c] = Counter++;
  }
}

// Fill a 2-d vector with a regular shape
void FillRegular( vvi &v, int &Counter )
{
  const std::size_t Rows{ v.size() };
  static constexpr int Cols{ MinCols + 1 };
  for( std::size_t r = 0; r < Rows; ++r )
  {
    v[r].resize( Cols );
    for( int c = 0; c < Cols; ++c )
      v[r][c] = Counter++;
  }
}

bool Debug()
{
  // Make 2-d data sets
  static constexpr int Small{ 5 };
  static constexpr int Large{ 100 };
  vvi vRagSmall{ Small };
  vvi vRegSmall{ Small };
  vvi vRagLarge{ Large };
  vvi vRegLarge{ Large };
  int Counter;
  Counter = 0; FillRagged(  vRagSmall, Counter );
  Counter = 0; FillRegular( vRegSmall, Counter );
  Counter = 0; FillRagged(  vRagLarge, Counter );
  Counter = 0; FillRegular( vRegLarge, Counter );

  // Make 3-d data sets
  static constexpr int OuterDim{ 4 };
  vvvi vRagRag{ OuterDim };
  vvvi vRagReg{ OuterDim };
  vvvi vRegReg{ OuterDim };
  int ctrRagRag{ 0 };
  int ctrRagReg{ 0 };
  int ctrRegReg{ 0 };
  for( int i = 0; i < OuterDim; ++i )
  {
    vRagRag[i].resize( i + MinCols );
    vRagReg[i].resize( i + MinCols );
    vRegReg[i].resize( Small );
    FillRagged ( vRagRag[i], ctrRagRag );
    FillRegular( vRagReg[i], ctrRagReg );
    FillRegular( vRegReg[i], ctrRegReg );
  }

  const std::string sRagSmall{ "RaggedSmall" };
  const std::string sRegSmall{ "RegularSmall" };
  const std::string sRagLarge{ "RaggedLarge" };
  const std::string sRegLarge{ "RegularLarge" };
  const std::string sRagRag  { "RaggedRagged" };
  const std::string sRagReg  { "RaggedRegular" };
  const std::string sRegReg  { "RegularRegular" };
  std::cout << sRagSmall << ": " << vRagSmall << std::endl;
  std::cout << sRegSmall << ": " << vRegSmall << std::endl;
  std::cout << sRagLarge << ": " << vRagLarge << std::endl;
  std::cout << sRegLarge << ": " << vRegLarge << std::endl;
  std::cout << sRagRag   << ": " << vRagRag   << std::endl;
  std::cout << sRagReg   << ": " << vRagReg   << std::endl;
  std::cout << sRegReg   << ": " << vRegReg   << std::endl;
  const std::string FileName{ "VectorDebug.h5" };
  {
    std::cout << "Writing " << FileName << std::endl;
    Grid::Hdf5Writer writer(FileName);
    write(writer, sRagSmall, vRagSmall);
    write(writer, sRegSmall, vRegSmall);
    write(writer, sRagLarge, vRagLarge);
    write(writer, sRegLarge, vRegLarge);
    write(writer, sRagRag,   vRagRag);
    write(writer, sRagReg,   vRagReg);
    write(writer, sRegReg,   vRegReg);
  }
  std::cout << "Reading back " << FileName << std::endl;
  Grid::Hdf5Reader reader(FileName);
  vvi rRegSmall;
  vvi rRagSmall;
  vvi rRegLarge;
  vvi rRagLarge;
  vvvi rRagRag;
  vvvi rRagReg;
  vvvi rRegReg;
  read( reader, sRegSmall, rRegSmall );
  read( reader, sRagSmall, rRagSmall );
  read( reader, sRegLarge, rRegLarge );
  read( reader, sRagLarge, rRagLarge );
  read( reader, sRagRag  , rRagRag   );
  read( reader, sRagReg  , rRagReg   );
  read( reader, sRegReg  , rRegReg   );
  std::cout << sRagSmall << ": " << rRagSmall << std::endl;
  std::cout << sRegSmall << ": " << rRegSmall << std::endl;
  std::cout << sRagLarge << ": " << rRagLarge << std::endl;
  std::cout << sRegLarge << ": " << rRegLarge << std::endl;
  std::cout << sRagRag   << ": " << rRagRag   << std::endl;
  std::cout << sRagReg   << ": " << rRagReg   << std::endl;
  std::cout << sRegReg   << ": " << rRegReg   << std::endl;
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
