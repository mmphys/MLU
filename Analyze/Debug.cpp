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
  bool Load( const std::string &FileName );
};

std::ostream & operator<<( std::ostream &os, const Common::FileNameAtt &fna )
{
  os << fna.NameNoExt;
  if( !fna.Ext.empty() )
    os << Common::Period << fna.Ext;
  return os;
}

static constexpr int nchannel{ 4 };
static constexpr int idxSrc{ 0 };
static constexpr int idxSnk{ 1 };

using Algebra = Common::Gamma::Algebra;
static const Algebra Gammas2pt[nchannel][2] = {
  { Algebra::Gamma5      , Algebra::Gamma5 },
  { Algebra::GammaTGamma5, Algebra::GammaTGamma5 },
  { Algebra::GammaTGamma5, Algebra::Gamma5 },
  { Algebra::Gamma5      , Algebra::GammaTGamma5 }
};
static const Algebra Gammas3pt[nchannel][2] = {
  { Algebra::Gamma5      , Algebra::GammaX },
  { Algebra::Gamma5      , Algebra::GammaY },
  { Algebra::Gamma5      , Algebra::GammaZ },
  { Algebra::Gamma5      , Algebra::GammaT }
};

// Load correlator - .h5 or .xml
bool MesonFile::Load( const std::string &FileName )
{
  bool b3pt{ false };
  std::size_t pos{ FileName.find_last_of( '.' ) };
  if( pos != std::string::npos && Common::EqualIgnoreCase( FileName.substr( pos + 1 ), "h5" ) )
  {
    {
      Common::FileNameAtt fna( FileName );
      if( fna.bGotDeltaT && fna.Gamma.size() )
        b3pt = true;
    }
    const Algebra (&Gammas)[nchannel][2]{ b3pt ? Gammas3pt : Gammas2pt };
    std::vector<Common::Gamma::Algebra> Alg;
    using Corr = Common::CorrelatorFileC;
    Corr corr( FileName, Alg, Alg, nullptr, "  Reading hdf5 " );
    data.resize( nchannel );
    for( int channel = 0; channel < nchannel; ++channel )
    {
      const std::complex<double> * pData = corr( Gammas[channel][idxSnk], Gammas[channel][idxSrc] );
      data[channel].resize( corr.Nt() );
      for( int i = 0; i < corr.Nt(); ++i )
        data[channel][i] = *pData++;
    }
  }
  else
  {
    std::cout << "  Reading xml " << FileName << std::endl;
    Grid::XmlReader r( FileName );
    Grid::read( r, "MesonFile", *this );
    if( data.size() != nchannel )
      throw std::runtime_error( "Expected " + std::to_string( nchannel )
                              + " channels, but read " + std::to_string( data.size() ) );
    if( data[0].size() == 0 )
      throw std::runtime_error( "Nt=0" );
    for( std::size_t i = 1; i < nchannel; ++i )
      if( data[i].size() != data[0].size() )
        throw std::runtime_error( "Nt[" + std::to_string( i ) + "] != " + std::to_string( data[0].size() ) );
  }
  return b3pt;
}

// Convert Correlator file from .h5 to .xml

void Convert( const std::string &FileName, std::string OutBase )
{
  Common::FileNameAtt fna( FileName );
  if( Common::EqualIgnoreCase( fna.Ext, "xml" ) )
    std::cout << "Doing nothing: " << fna.NameNoExt << fna.Ext << " is xml" << std::endl;
  else
  {
    MesonFile mf;
    mf.Load( FileName );
    OutBase.append( fna.NameNoExt );
    OutBase.append( ".xml" );
    Grid::XmlWriter w( OutBase );
    Grid::write( w, "MesonFile", mf );
  }
}

// Compare two Correlator files (.h5 and .xml)

struct Comparator
{
public:
  const double tol;
  const double precision;
  const unsigned int NumDigits;
  // These values are the result of a comparison
  bool bSame;
  bool bSignFlip;
  bool bAltOK; // True if this passed a relative test after failing an absolute test
  double Value;
  std::string StringValue;
protected:
  std::stringstream ss;
  unsigned int GetDigits( double Tol, double p )
  {
    if( p == 0 )
      p = Tol;
    unsigned int Digits = 0;
    if( p < 1. )
      Digits = static_cast<unsigned int>( 0.999 - std::log10( p ) );
    return Digits;
  }
  void MakePos( double &Num )
  {
    if( Num < 0 )
    {
      bSignFlip = !bSignFlip;
      Num = -Num;
    }
  }
public:
  Comparator( double Tol, double Precision )
  : tol{Tol}, precision{Precision}, NumDigits{GetDigits(Tol,Precision)} {}
  std::string CompareType() const { return std::string( 1, precision ? '-' : '/' ); }
  bool Compare( double Num1, double Num2 )
  {
    bSame = false;
    bSignFlip = false;
    bAltOK = false;
    MakePos( Num1 );
    MakePos( Num2 );
    if( precision )
    {
      Value = std::abs( Num1 - Num2 );
      unsigned long long iNum = static_cast<unsigned long long >( Value / precision + 0.5 );
      Value = static_cast<double>( iNum * precision );
      ss << std::setprecision( NumDigits ) << Value;
      StringValue = ss.str();
      bSame = StringValue.size() == 1 && StringValue[0] == '0';
    }
    if( !bSame && tol )
    {
      double tmpValue = Num1 / Num2;
      ss.str(std::string());
      ss << std::setprecision( NumDigits + ( tmpValue >= 1 ? 1 : 0 ) ) << tmpValue;
      std::string tmpStringValue = ss.str();
      bSame = tmpStringValue.size() == 1 && tmpStringValue[0] == '1';
      if( bSame || !precision )
      {
        Value = tmpValue;
        StringValue = tmpStringValue;
        if( precision )
          bAltOK = true;
      }
    }
    ss.str(std::string());
    return bSame;
  }
};

int Compare( const std::vector<std::string> &FileName, std::string OutBase, Comparator &c )
{
  bool b3pt{ false };
  const std::size_t OutLen{ OutBase.size() };
  std::vector<MesonFile> in ( FileName.size() );
  int Count{ 0 };
  int SignFlip{ 0 };
  for( std::size_t f = 0; f < FileName.size(); ++f )
  {
    if( in[f].Load( FileName[f] ) )
      b3pt = true;
    const Algebra (&Gammas)[nchannel][2]{ b3pt ? Gammas3pt : Gammas2pt };
    const std::size_t Nt{ in[f].data[0].size() };
    if( f )
    {
      if( Nt != in[0].data[0].size() )
        throw std::runtime_error( "Nt mismatch" );
      std::ofstream out; // Will be closed automatically if there's an exception
      {
        std::size_t pStart{ FileName[f].find_last_of( '/' ) };
        if( pStart == std::string::npos )
          pStart = 0;
        else
          pStart++;
        std::size_t pEnd{ FileName[f].find_last_of( '.' ) };
        if( pEnd == std::string::npos || pEnd < pStart )
          pEnd = FileName[f].size();
        OutBase.resize( OutLen );
        OutBase.append( FileName[f].substr( pStart, pEnd - pStart ) );
        OutBase.append( ".txt" );
        std::cout << "  Creating " << OutBase << std::endl;
        out.open( OutBase );
      }
      static const std::string Comment{ "# " };
      const std::string HeadReal{ "Snk Src t Re(f1) Re(f2) Re(f2)" + c.CompareType() + "Re(f1)" };
      const std::string HeadImag{ "Snk Src t Im(f1) Im(f2) Im(f2)" + c.CompareType() + "Im(f1)" };
      using TReIm = double (Grid::Complex::*)() const;
      TReIm ReIm = b3pt ? static_cast<TReIm>(&Grid::Complex::imag) : static_cast<TReIm>(&Grid::Complex::real);
      out << Comment << "File 1 " << FileName[0] << Common::NewLine
          << Comment << "File 2 " << FileName[f] << Common::NewLine
          << ( b3pt ? HeadImag : HeadReal )
          << std::endl;
      for( std::size_t channel = 0; channel < nchannel; ++channel )
      {
        if( channel == nchannel - 1 && b3pt )
        {
          ReIm = &Grid::Complex::real;
          out << "\n\n" << HeadReal << std::endl;
        }
        for( std::size_t t = 0; t < Nt; ++t )
        {
          const double Num1{ (in[0].data[channel][t].*ReIm)() };
          const double Num2{ (in[f].data[channel][t].*ReIm)() };
          if( !c.Compare( Num1, Num2 ) )
            Count++;
          if( c.bSignFlip )
            SignFlip++;
          out << Common::Gamma::NameShort( Gammas[channel][idxSnk] ) << Common::Space
              << Common::Gamma::NameShort( Gammas[channel][idxSrc] ) << Common::Space
              << t << Common::Space
              << std::scientific << std::setprecision( std::numeric_limits<double>::digits10 + 1 )
              << Num1 << Common::Space
              << Num2 << Common::Space
              << c.StringValue;
          if( !c.bSame )
            out << " different";
          if( c.bAltOK )
            out << " rel";
          if( c.bSignFlip )
            out << " -";
          out << std::endl;
        }
      }
    }
  }
  std::stringstream s;
  s << "  ";
  if( Count == 0 )
    s << "Same";
  else
    s << Count << " differences";
  if( SignFlip )
    s << " and " << SignFlip << " sign flips";
  s << " at" << std::scientific << std::setprecision(0);
  if( c.precision )
    s << " absolute dufference " << c.precision;
  if( c.tol )
    s << " relative error " << c.tol;
  std::cout << s.str() << std::endl;
  return Count;
}

void GridDebug(int argc, char *argv[])
{
  //std::ios_base::sync_with_stdio( false );
  Grid::Grid_init( &argc, &argv );
  try
  {
    std::cout << Grid::GridLogMessage << MLUVersionInfoHuman() << std::endl;
    std::cout << Grid::GridLogMessage << "Hello" << std::endl;
    std::cout << Grid::GridLogMessage << "Boost Spirit X3 version 0x"
              << std::hex << SPIRIT_X3_VERSION << std::dec << std::endl;
    std::cout << Grid::GridLogMessage << MakeSeed( 0, 1 ) << std::endl;
    std::cout << Grid::GridLogMessage << MakeSeed( 5, 8 ) << std::endl;
    Debug();
  } catch( ... ) {
    Grid::Grid_finalize();
    throw;
  }
  Grid::Grid_finalize();
}

// Grab the string following this parameter from the command-line
const char * NextString( int argc, char *argv[], int &i )
{
  if( argv[i][2] )
    throw new std::runtime_error( std::string( argv[i] ) + " should be single-character" );
  if( i >= argc - 1 )
    throw new std::runtime_error( std::string( argv[i] ) + " should be followed by an argument" );
  return argv[++i];
}

int main(int argc, char *argv[])
{
  int iReturn = EXIT_SUCCESS;
  try
  {
    if( argc <= 1 )
    {
      GridDebug( argc, argv );
    }
    else
    {
      bool bConvertXml{ false };
      std::vector<const char *> argsRaw;
      const char * InBase = "", * OutBase = "";
      static constexpr double tolDefault{ 1e-5 };
      double tol{ 0 };
      double precision{ 0 };
      for( int i = 1; i < argc; ++i )
      {
        if( argv[i][0] == '-' )
        {
          switch( argv[i][1] )
          {
            case 'i':
              InBase = NextString( argc, argv, i );
              break;
            case 'o':
              OutBase = NextString( argc, argv, i );
              break;
            case 'd':
              precision = Common::FromString<double>( NextString( argc, argv, i ) );
              break;
            case 'p':
              tol = Common::FromString<double>( NextString( argc, argv, i ) );
              break;
            case 'x':
              bConvertXml = true;
              break;
            case '?':
            case 'h':
              std::cout << "-i Input prefix, prepended to each input file\n"
                           "-o Output prefix\n"
                           "-d compare using absolute Difference\n"
                           "-p compare using relative Precision (error)\n"
                           "-x Convert input files to .xml. Otherwise compares files\n"
                           "-h This help message\n"
                           "If neither -d nor -p specified, use -p " << tolDefault << std::endl;
              break;
            default:
              std::cout << "Warning: assuming " << argv[i] << " is a Grid option" << std::endl;
              i = argc - 1; // Stop - We've hit Grid options
              break;
          }
        }
        else
          argsRaw.emplace_back( argv[i] );
      }
      // Default arguments
      if( !tol && !precision )
        tol = tolDefault;
      // Parse arguments
      std::vector<std::string> args{ Common::glob( argsRaw.begin(), argsRaw.end(), InBase ) };
      if( !bConvertXml )
      {
        if( args.size() < 2 )
          std::cout << args.size() << " files to compare" << std::endl;
        else
        {
          Comparator c( tol, precision );
          if( Compare( args, OutBase, c ) )
            iReturn = EXIT_FAILURE;
        }
      }
      else
        for( const std::string &f : args )
          Convert( f, OutBase );
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
  return iReturn;
}
