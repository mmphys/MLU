//
//  Debug.cpp
//  A place to play with debug code
//
//  Created by Michael Marshall on 09/05/2020.
//  Copyright © 2020 sopa. All rights reserved.
//

#include <functional>
#include <stdio.h>
#include <ostream>
#include <random>
#include <MLU/Common.hpp>
#include <omp.h>
#include <Grid/Grid.h>

//#include <boost/spirit/home/x3/version.hpp>

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
  const double relError;
  const double absError;
  const unsigned int relDigits;
  const unsigned int absDigits;
  // These values are the result of a comparison
  bool bSame;
  bool bSignFlip;
  bool bAltOK; // True if this passed a relative test after failing an absolute test
  double Value;
  std::string StringValue;
protected:
  std::stringstream ss;
  unsigned int GetDigits( double p )
  {
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
  Comparator( double RelError, double AbsError )
  : relError{RelError}, absError{AbsError}, relDigits{GetDigits(RelError)}, absDigits{GetDigits(absError)} {}
  std::string CompareType() const { return std::string( 1, absError ? '-' : '/' ); }
  bool Compare( double Num1, double Num2 )
  {
    bSame = false;
    bSignFlip = false;
    bAltOK = false;
    MakePos( Num1 );
    MakePos( Num2 );
    if( relError )
    {
      Value = Num1 / Num2;
      ss << std::setprecision( relDigits + ( Value >= 1 ? 1 : 0 ) ) << Value;
      StringValue = ss.str();
      bSame = StringValue.size() == 1 && StringValue[0] == '1';
    }
    if( !bSame && absError )
    {
      double tmpValue = std::abs( Num1 - Num2 );
      unsigned long long iNum = static_cast<unsigned long long >( tmpValue / absError + 0.5 );
      tmpValue = static_cast<double>( iNum * absError );
      ss.str(std::string());
      ss << std::setprecision( absDigits ) << tmpValue;
      std::string tmpStringValue = ss.str();
      bSame = tmpStringValue.size() == 1 && tmpStringValue[0] == '0';
      if( bSame || !relError )
      {
        Value = tmpValue;
        StringValue = tmpStringValue;
        if( relError )
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
          << Comment << "File 2 " << FileName[f] << Common::NewLine << std::defaultfloat;
      if( c.absError )
        out << Comment << "Absolute error: Num2 - Num1 <= " << c.absError << Common::NewLine;
      if( c.relError )
        out << Comment << "Relative error: Num2 / Num1 <= " << c.relError << Common::NewLine;
      out << ( b3pt ? HeadImag : HeadReal ) << std::endl;
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
          static constexpr int PrintPrecision{ std::numeric_limits<double>::digits10 + 1 };
          out << Common::Gamma::NameShort( Gammas[channel][idxSnk] ) << Common::Space
              << Common::Gamma::NameShort( Gammas[channel][idxSrc] ) << Common::Space
              << std::setw(3) << t << Common::Space
              << std::scientific << std::setprecision( PrintPrecision )
              << std::setw( PrintPrecision + 7 ) << Num1 << Common::Space
              << std::setw( PrintPrecision + 7 ) << Num2 << Common::Space
              << c.StringValue;
          if( !c.bSame )
            out << " different";
          if( c.bAltOK )
            out << " abs";
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
  if( c.absError )
    s << " absolute error " << c.absError;
  if( c.relError )
    s << " relative error " << c.relError;
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
    //std::cout << Grid::GridLogMessage << "Boost Spirit X3 version 0x"
    //          << std::hex << SPIRIT_X3_VERSION << std::dec << std::endl;
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

bool FitRangeTest( int argc, char *argv[] )
{
  bool bReturn{ true };
  try
  {
    if( argc <= 1 )
      std::cout << "Please specify one or more fit ranges" << std::endl;
    else
    {
      // Make FitRanges
      std::vector<std::string> vs;
      vs.reserve( argc - 1 );
      for (int i = 1; i < argc; ++i)
        vs.push_back( argv[i] );
      Common::FitRanges fr( vs, 1 );
      // Show FitRanges
      std::cout << fr << std::endl;
      // Iterate FitRanges
      std::size_t Count{};
      for( Common::FitRangesIterator it = fr.begin(); !it.PastEnd(); ++it )
      {
        std::cout << Count++ << ": " << it.AbbrevString( "-", ", " )
                  << " --- " << it.AbbrevString()
                  << " --- " << it
                  << std::endl;
      }
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    bReturn = false;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    bReturn = false;
  }
  return bReturn;
}

// Data for reproduction of tables in section 8, pg 12
// https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/ckelly/Gparity/hotelling_v10.pdf
struct CKelly
{
  struct Entry {
    double ChiSq;
    double ChiSqDof;
    double pHotelling;
    double pChi2;
  };

  using Table = std::vector<Entry>;

  const Common::HotellingDist TSq;
  const Table &t;

  CKelly( unsigned int Dof, unsigned int n, const Table &t_ ) : TSq(Dof, n - 1), t{t_} {}
  inline double RelativeError( double New, double Original ) const { return (New - Original) / Original; };
  bool Test() const;
  void Test( const Entry &e ) const;
};

std::ostream & operator<<( std::ostream &os, const CKelly::Entry &e )
{
  return os << "chisq=" << e.ChiSq << " chisq/dof=" << e.ChiSqDof
            << " p(T2)=" << e.pHotelling << " p(chi2)=" << e.pChi2;
}

void CKelly::Test( const Entry &e ) const
{
  const double ChisqQ{ Common::qValueChiSq( e.ChiSq, TSq.p ) };
  const double FDistQ{ TSq( e.ChiSq ) };
  std::cout << "CHISQ=" << e.ChiSq << " CHISQ/dof=" << e.ChiSq / TSq.p
            << " P(T2)=" << FDistQ << " P(CHI2)=" << ChisqQ
            << "\n" << std::string( 67, ' ' )
            << " RelError(T2)="    << std::setw(12) << RelativeError( FDistQ, e.pHotelling )
            << " RelError(CHI2)="  << std::setw(12) << RelativeError( ChisqQ, e.pChi2 )
            << std::endl;
}

bool CKelly::Test() const
{
  std::cout << "\nTEST FOR k=" << TSq.p << ", N=" << (TSq.m+1) << " lower/UPPER = original/REPRODUCTION" << std::endl;
  for( const Entry &e : t )
  {
    std::cout << e << std::endl;
    Test( e );
  }
  return true;
}

const CKelly::Table CKellyTable1{
  { 0, 0, 1, 1 },
  { 1.86667, 0.0666667, 1, 1 },
  { 3.73333, 0.133333, 1, 1 },
  { 5.6, 0.2, 1, 0.999998 },
  { 7.46667, 0.266667, 0.999996, 0.999963 },
  { 9.33333, 0.333333, 0.999954, 0.99964 },
  { 11.2, 0.4, 0.999732, 0.998019 },
  { 13.0667, 0.466667, 0.998924, 0.992597 },
  { 14.9333, 0.533333, 0.996705, 0.97913 },
  { 16.8, 0.6, 0.991781, 0.952436 },
  { 18.6667, 0.666667, 0.982501, 0.908104 },
  { 20.5333, 0.733333, 0.967142, 0.844255 },
  { 22.4, 0.8, 0.944241, 0.762445 },
  { 24.2667, 0.866667, 0.912903, 0.667385 },
  { 26.1333, 0.933333, 0.872966, 0.565716 },
  { 28, 1, 0.825027, 0.464448 },
  { 29.8667, 1.06667, 0.770332, 0.36962 },
  { 31.7333, 1.13333, 0.710579, 0.285498 },
  { 33.6, 1.2, 0.647696, 0.214354 },
  { 35.4667, 1.26667, 0.583629, 0.156687 },
  { 37.3333, 1.33333, 0.520173, 0.11169 },
  { 39.2, 1.4, 0.458856, 0.0777593 },
  { 41.0667, 1.46667, 0.400877, 0.0529546 },
  { 42.9333, 1.53333, 0.347096, 0.0353251 },
  { 44.8, 1.6, 0.298046, 0.0231134 },
  { 46.6667, 1.66667, 0.253978, 0.0148518 },
  { 48.5333, 1.73333, 0.214913, 0.00938263 },
  { 50.4, 1.8, 0.180691, 0.0058338 },
  { 52.2667, 1.86667, 0.151031, 0.0035734 },
  { 54.1333, 1.93333, 0.125566, 0.00215824 },
};

CKelly CK1( 28, 100, CKellyTable1 );

const CKelly::Table CKellyTable2{
  { 0, 0, 1, 1 },
  { 1.86667, 0.0666667, 1, 1 },
  { 3.73333, 0.133333, 1, 1 },
  { 5.6, 0.2, 0.999999, 0.999998 },
  { 7.46667, 0.266667, 0.999969, 0.999963 },
  { 9.33333, 0.333333, 0.999701, 0.99964 },
  { 11.2, 0.4, 0.998341, 0.998019 },
  { 13.0667, 0.466667, 0.993751, 0.992597 },
  { 14.9333, 0.533333, 0.982224, 0.97913 },
  { 16.8, 0.6, 0.959083, 0.952436 },
  { 18.6667, 0.666667, 0.920095, 0.908104 },
  { 20.5333, 0.733333, 0.863031, 0.844255 },
  { 22.4, 0.8, 0.788605, 0.762445 },
  { 24.2667, 0.866667, 0.700429, 0.667385 },
  { 26.1333, 0.933333, 0.604111, 0.565716 },
  { 28, 1, 0.505964, 0.464448 },
  { 29.8667, 1.06667, 0.411788, 0.36962 },
  { 31.7333, 1.13333, 0.326039, 0.285498 },
  { 33.6, 1.2, 0.251481, 0.214354 },
  { 35.4667, 1.26667, 0.189248, 0.156687 },
  { 37.3333, 1.33333, 0.139159, 0.11169 },
  { 39.2, 1.4, 0.100138, 0.0777593 },
  { 41.0667, 1.46667, 0.0706187, 0.0529546 },
  { 42.9333, 1.53333, 0.048873, 0.0353251 },
  { 44.8, 1.6, 0.0332355, 0.0231134 },
  { 46.6667, 1.66667, 0.022235, 0.0148518 },
  { 48.5333, 1.73333, 0.0146506, 0.00938263 },
  { 50.4, 1.8, 0.00951694, 0.0058338 },
  { 52.2667, 1.86667, 0.00610064, 0.0035734 },
  { 54.1333, 1.93333, 0.00386248, 0.00215824 },
};

CKelly CK2( 28, 1000, CKellyTable2 );

const CKelly::Table CKellyTable3{
  { 0, 0, 1, 1 },
  { 1.86667, 0.0666667, 1, 1 },
  { 3.73333, 0.133333, 1, 1 },
  { 5.6, 0.2, 0.999998, 0.999998 },
  { 7.46667, 0.266667, 0.999964, 0.999963 },
  { 9.33333, 0.333333, 0.999647, 0.99964 },
  { 11.2, 0.4, 0.998053, 0.998019 },
  { 13.0667, 0.466667, 0.99272, 0.992597 },
  { 14.9333, 0.533333, 0.979458, 0.97913 },
  { 16.8, 0.6, 0.953136, 0.952436 },
  { 18.6667, 0.666667, 0.909359, 0.908104 },
  { 20.5333, 0.733333, 0.846205, 0.844255 },
  { 22.4, 0.8, 0.765141, 0.762445 },
  { 24.2667, 0.866667, 0.67076, 0.667385 },
  { 26.1333, 0.933333, 0.569599, 0.565716 },
  { 28, 1, 0.468603, 0.464448 },
  { 29.8667, 1.06667, 0.373793, 0.36962 },
  { 31.7333, 1.13333, 0.289461, 0.285498 },
  { 33.6, 1.2, 0.217937, 0.214354 },
  { 35.4667, 1.26667, 0.159787, 0.156687 },
  { 37.3333, 1.33333, 0.114267, 0.11169 },
  { 39.2, 1.4, 0.0798266, 0.0777593 },
  { 41.0667, 1.46667, 0.05456, 0.0529546 },
  { 42.9333, 1.53333, 0.0365354, 0.0353251 },
  { 44.8, 1.6, 0.0240016, 0.0231134 },
  { 46.6667, 1.66667, 0.0154875, 0.0148518 },
  { 48.5333, 1.73333, 0.00982729, 0.00938263 },
  { 50.4, 1.8, 0.00613831, 0.0058338 },
  { 52.2667, 1.86667, 0.00377785, 0.0035734 },
  { 54.1333, 1.93333, 0.00229303, 0.00215824 },
};

CKelly CK3( 28, 10000, CKellyTable3 );

const CKelly::Table CKellyTable4{
  { 0, 0, 1, 1 },
  { 0.333333, 0.0666667, 0.99716, 0.996969 },
  { 0.666667, 0.133333, 0.985721, 0.984748 },
  { 1, 0.2, 0.964968, 0.962566 },
  { 1.33333, 0.266667, 0.935871, 0.931465 },
  { 1.66667, 0.333333, 0.899932, 0.893072 },
  { 2, 0.4, 0.858766, 0.849145 },
  { 2.33333, 0.466667, 0.813911, 0.801359 },
  { 2.66667, 0.533333, 0.766735, 0.751212 },
  { 3, 0.6, 0.718414, 0.699986 },
  { 3.33333, 0.666667, 0.669918, 0.648742 },
  { 3.66667, 0.733333, 0.622029, 0.598332 },
  { 4, 0.8, 0.575358, 0.549416 },
  { 4.33333, 0.866667, 0.530366, 0.502488 },
  { 4.66667, 0.933333, 0.487386, 0.457898 },
  { 5, 1, 0.446646, 0.41588 },
  { 5.33333, 1.06667, 0.408284, 0.376568 },
  { 5.66667, 1.13333, 0.372369, 0.340016 },
  { 6, 1.2, 0.338912, 0.306219 },
  { 6.33333, 1.26667, 0.307884, 0.275122 },
  { 6.66667, 1.33333, 0.279218, 0.246634 },
  { 7, 1.4, 0.252827, 0.22064 },
  { 7.33333, 1.46667, 0.228605, 0.197007 },
  { 7.66667, 1.53333, 0.206435, 0.175588 },
  { 8, 1.6, 0.186193, 0.156236 },
  { 8.33333, 1.66667, 0.167753, 0.138797 },
  { 8.66667, 1.73333, 0.150989, 0.123121 },
  { 9, 1.8, 0.135775, 0.109064 },
  { 9.33333, 1.86667, 0.121992, 0.0964848 },
  { 9.66667, 1.93333, 0.109523, 0.0852502 },
};

CKelly CK4( 5, 100, CKellyTable4 );

const CKelly::Table MikeTable1{
  { 234.8155764934, 2.257842081667, -1, -2 },
};

CKelly CompareMikeFit( 104, 128, MikeTable1 );

class Callable
{
protected:
  int Magic;
public:
  Callable( int Magic_ ) : Magic{Magic_} {}
  void DoMagic( int Plus );
};

void Callable::DoMagic( int Plus )
{
  std::cout << "Magic " << Magic << " + " << Plus << " = ";
  Magic += Plus;
  std::cout << Magic << std::endl;
}

static const std::function<void(int)> NullFunction;

void DoCaller( int i, const std::function<void(int)> &f = NullFunction )
{
  for( int j = 0; j < i; ++j )
  {
    std::cout << j << ": ";
    if( f )
      f( j );
  }
  std::cout << "Done\n";
}

bool TryCallback()
{
  Callable mA( 7 ), mB( 42 );
  std::function<void(int)> CallBack{ std::bind( &Callable::DoMagic, mB, std::placeholders::_1 ) };
  DoCaller( 3, CallBack );
  DoCaller( 4 );
  return false;
}

int main(int argc, char *argv[])
{
  std::cout << "Your random number is " << std::random_device()() << std::endl;
  //if( !TryCallback() ) return EXIT_SUCCESS;
  //if( CK1.Test() && CK2.Test() && CK3.Test() && CK4.Test() && CompareMikeFit.Test() ) return EXIT_SUCCESS;
  //if( FitRangeTest( argc, argv ) ) return EXIT_SUCCESS;
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
      static constexpr double relErrorDefault{ 1e-5 };
      double relError{ 0 };
      double absError{ 0 };
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
              absError = Common::FromString<double>( NextString( argc, argv, i ) );
              break;
            case 'p':
              relError = Common::FromString<double>( NextString( argc, argv, i ) );
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
                           "If neither -d nor -p specified, use -p " << relErrorDefault << std::endl;
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
      if( !relError && !absError )
        relError = relErrorDefault;
      // Parse arguments
      std::vector<std::string> args{ Common::glob( argsRaw.begin(), argsRaw.end(), InBase ) };
      if( !bConvertXml )
      {
        if( args.size() < 2 )
          std::cout << args.size() << " files to compare" << std::endl;
        else
        {
          Comparator c( relError, absError );
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
