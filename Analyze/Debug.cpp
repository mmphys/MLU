//
//  Debug.cpp
//  A place to play with debug code
//
//  Created by Michael Marshall on 09/05/2020.
//  Copyright © 2020 sopa. All rights reserved.
//

#include <stdio.h>
#include <MLU/Common.hpp>
#include <omp.h>
#include <Grid/Grid.h>

#include <boost/spirit/home/x3/version.hpp>

using Grid::operator<<;
using Grid::operator>>;

using GridScalar = float;
using GridTensor = Grid::iVector<Grid::iVector<GridScalar, 3>, 2>;
using vGridTensor = std::vector<GridTensor>;

using  vvi = std::vector<std::vector<int>>;
using vvvi = std::vector<vvi>;

/*template<typename T>
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
}*/

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
  vvvi vTrickyReg(2);
  Counter = 0;
  for( int i = 0; i < vTrickyReg.size(); ++i )
  {
    vTrickyReg[i].resize(2 - i);
    for( int j = 0; j < vTrickyReg[i].size(); ++j )
    {
      vTrickyReg[i][j].resize(2);
      for( int k = 0; k < vTrickyReg[i][j].size(); ++k )
        vTrickyReg[i][j][k] = Counter++;
    }
  }
  vGridTensor tGrid( 4 );
  for( auto &v : tGrid )
    for( auto &s : v )
      s = Counter++;

  const std::string sRagSmall{ "RaggedSmall" };
  const std::string sRegSmall{ "RegularSmall" };
  const std::string sRagLarge{ "RaggedLarge" };
  const std::string sRegLarge{ "RegularLarge" };
  const std::string sRagRag  { "RaggedRagged" };
  const std::string sRagReg  { "RaggedRegular" };
  const std::string sRegReg  { "RegularRegular" };
  const std::string sTrickyReg{ "TrickyRegular" };
  const std::string sInt  { "Int1d" };
  const std::string stGrid{ "Grid4d" };
  std::vector<int> vFibo{ 1, 1, 2, 3, 5, 8, 13, 21 };
  std::cout << sRagSmall << ": " << vRagSmall << std::endl;
  std::cout << sRegSmall << ": " << vRegSmall << std::endl;
  std::cout << sRagLarge << ": " << vRagLarge << std::endl;
  std::cout << sRegLarge << ": " << vRegLarge << std::endl;
  std::cout << sRagRag   << ": " << vRagRag   << std::endl;
  std::cout << sRagReg   << ": " << vRagReg   << std::endl;
  std::cout << sRegReg   << ": " << vRegReg   << std::endl;
  std::cout << sTrickyReg<< ": " << vTrickyReg<< std::endl;
  std::cout << sInt      << ": " << vFibo     << std::endl;
  std::cout << stGrid    << ": " << tGrid     << std::endl;
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
    write(writer, sTrickyReg,vTrickyReg);
    write(writer, sInt,      vFibo);
    write(writer, stGrid,    tGrid );
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
  vvvi rTrickyReg;
  std::vector<int> rFibo;
  read( reader, sRegSmall, rRegSmall );
  try{
    read( reader, sRagSmall, rRagSmall );
  } catch(...) {
    std::cout << "Exception reading " << sRagSmall << std::endl;
  }
  read( reader, sRegLarge, rRegLarge );
  try{
    read( reader, sRagLarge, rRagLarge );
  } catch(...) {
    std::cout << "Exception reading " << sRagLarge << std::endl;
  }
  try{
    read( reader, sRagRag  , rRagRag   );
  } catch(...) {
    std::cout << "Exception reading " << sRagRag << std::endl;
  }
  try{
    read( reader, sRagReg  , rRagReg   );
  } catch(...) {
    std::cout << "Exception reading " << sRagReg << std::endl;
  }
  read( reader, sRegReg  , rRegReg   );
  read( reader, sTrickyReg,rTrickyReg);
  read( reader, sInt     , rFibo     );
  vGridTensor rGrid;
  read( reader, stGrid, rGrid );
  std::cout << sRagSmall << ": " << rRagSmall << std::endl;
  std::cout << sRegSmall << ": " << rRegSmall << std::endl;
  std::cout << sRagLarge << ": " << rRagLarge << std::endl;
  std::cout << sRegLarge << ": " << rRegLarge << std::endl;
  std::cout << sRagRag   << ": " << rRagRag   << std::endl;
  std::cout << sRagReg   << ": " << rRagReg   << std::endl;
  std::cout << sRegReg   << ": " << rRegReg   << std::endl;
  std::cout << sTrickyReg<< ": " << rTrickyReg<< std::endl;
  std::cout << sInt      << ": " << rFibo     << std::endl;
  std::cout << stGrid    << ": " << rGrid     << std::endl;
  return true;
}

std::string MakeSeed( int timeslice, int hit = 1 )
{
  std::string s{"Study2-Source_Z2_p_0_0_0_t_"};
  s.append( std::to_string( timeslice ) );
  if( hit != 1 )
  {
    s.append( "_hit_" );
    s.append( std::to_string( hit ) );
  }
    static const std::string Suffix{"-880"};
  s.append( "-880" );
  return s;
}

class MesonFile: Grid::Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(MesonFile, std::vector<std::vector<Grid::Complex> >, data);
};

void Convert( const std::string &FileName )
{
  static constexpr int nchannel{ 4 };
  Common::Gamma::Algebra Gammas[nchannel][2] = {
    {Common::Gamma::Algebra::Gamma5      ,Common::Gamma::Algebra::Gamma5},
    {Common::Gamma::Algebra::GammaTGamma5,Common::Gamma::Algebra::GammaTGamma5},
    {Common::Gamma::Algebra::GammaTGamma5,Common::Gamma::Algebra::Gamma5},
    {Common::Gamma::Algebra::Gamma5      ,Common::Gamma::Algebra::GammaTGamma5}
  };
  std::vector<Common::Gamma::Algebra> Alg;
  using Corr = Common::CorrelatorFileC;
  Corr corr( FileName, Alg, Alg, nullptr, "Converting " );
  MesonFile mf;
  mf.data.resize( nchannel );
  for( int channel = 0; channel < nchannel; ++channel )
  {
    const std::complex<double> * pData = corr( Gammas[channel][1], Gammas[channel][0] );
    mf.data[channel].resize( corr.Nt() );
    for( int i = 0; i < corr.Nt(); ++i )
      mf.data[channel][i] = *pData++;
  }
  Grid::XmlWriter w( corr.Name_.NameNoExt + ".xml" );
  Grid::write( w, "MesonFile", mf );
}

int main(int argc, char *argv[])
{
  //std::ios_base::sync_with_stdio( false );
  Grid::Grid_init(&argc,&argv);
  std::cout << Grid::GridLogMessage << MLUVersionInfoHuman() << std::endl;
  std::cout << Grid::GridLogMessage << "Hello" << std::endl;
  std::cout << Grid::GridLogMessage << "Boost Spirit X3 version 0x"
            << std::hex << SPIRIT_X3_VERSION << std::dec << std::endl;
  std::cout << Grid::GridLogMessage << MakeSeed( 0, 1 ) << std::endl;
  std::cout << Grid::GridLogMessage << MakeSeed( 5, 8 ) << std::endl;
  int iReturn = EXIT_SUCCESS;
  try
  {
    //#pragma omp parallel
    std::cout << "Hello world!\n";
    if( argc <= 1 )
    {
      if( !Debug() )
        iReturn = EXIT_FAILURE;
    }
    else
    {
      for( const std::string &f : Common::glob(&argv[1], &argv[argc]))
        Convert( f );
    }
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
