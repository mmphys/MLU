/**

 Mike's lattice QCD utilities
 
 Source file: Common.cpp
 
 Copyright (C) 2019
 
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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
**/

// Common utilities (no dependencies other than c++ stdlib)

#include "Common.hpp"

#include <cstddef>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <sys/stat.h>
#include <H5CompType.h>

// This should be the only use of the AutoConf header
// Any useful values should be exposed as functions

#include <MLUconfig.h>

extern "C" const char * MLUVersionInfoHuman()
{
  static const char PackageString[] = MLU_PACKAGE_STRING ", " MLU_GIT_SUMMARY;
  return PackageString;
}

BEGIN_COMMON_NAMESPACE

// I will use the GNU Scientific Library, but I will have GSL errors throw exceptions
class GSLLibraryGlobal
{
private:
  gsl_error_handler_t * pOldErrorHandler;
  static void GSLErrorHandler(const char * reason, const char * file, int line, int gsl_errno );

public:
  GSLLibraryGlobal();
  ~GSLLibraryGlobal();
};

GSLLibraryGlobal gslLibraryGlobal;

GSLLibraryGlobal::GSLLibraryGlobal()
{
  pOldErrorHandler = gsl_set_error_handler( GSLErrorHandler );
}

GSLLibraryGlobal::~GSLLibraryGlobal()
{
  gsl_set_error_handler( pOldErrorHandler );
}

void GSLLibraryGlobal::GSLErrorHandler(const char * reason, const char * file, int line, int gsl_errno )
{
  std::stringstream ss;
  ss << "gsl: " << file << ":" << line << ": errno " << gsl_errno << ": " << reason;
  throw std::runtime_error( ss.str() );
}

// Text required for summaries of correlators
namespace CorrSumm {
  const char sep[] = " ";
  const char Comment[] = "# ";
  const char * FieldNames[NumFields] = { "corr", "exp", "cosh" };
};

extern const std::string Empty{ "" };
extern const std::string Space{ " " };
extern const std::string WhiteSpace{ " \t\n\r\f\v" };
extern const std::string Underscore{ "_" };
extern const std::string Period{ "." };
extern const std::string Colon{ ":" };
extern const std::string NewLine{ "\n" };
extern const std::string Comma{ "," };
extern const std::string CommaSpace{ ", " };
extern const std::string Hash{ "#" };
extern const std::string EqualSign{ "=" };
const std::string sBootstrap{ "bootstrap" };
const std::string sFold{ "fold" };
const std::string sModel{ "model" };
const std::string sParams{ "params" };
const std::string sCovmat{ "covmat" };
const std::string sCovmatIn{ "covmat_in" };
const std::string sCovmatInv{ "covmat_inv" };
const std::string sCovmatInvCholesky{ sCovmatInv + "chol" };
const std::string sCormat{ "cormat" };
const std::string sCormatCholesky{ sCormat + "_chol" };
const std::string sCormatInvCholesky{ sCormat + "_invchol" };
const std::string sNtUnfolded{ "NtUnfolded" };
const std::string st0Negated{ "t0Negated" };
const std::string sConjugated{ "Conjugated" };
const std::string sRawBootstrap{ "RawBootstrap" };
const std::string sTI{ "TI" };
const std::string sTF{ "TF" };
const std::string sDoF{ "DoF" };
const std::string sNumExponents{ "NumExponents" };
const std::string sNumFiles{ "NumFiles" };
const std::string sFactorised{ "Factorised" };
const std::string sCovarFrozen{ "CovarFrozen" };
const std::string s_C{ "_C" };
const std::string sStdErrorMean{ "StdErrorMean" };
const std::string sErrorScaled{ "ErrorScaled" };
const std::string sCovariance{ "Covariance" };
const std::string sCovarianceIn{ sCovariance + "In" };
const std::string sCovarianceInv{ sCovariance + "Inv" };
const std::string sCovarianceInvCholesky{ sCovarianceInv + "Cholesky" };
const std::string sCorrelation{ "Correlation" };
const std::string sCorrelationCholesky{ sCorrelation + "Cholesky" };
const std::string sCorrelationInvCholesky{ sCorrelation + "InvCholesky" };
const std::string sFitInput{ "FitInput" };
const std::string sModelPrediction{ "ModelPrediction" };
const std::string sOperators{ "Operators" };
const std::string sAuxNames{ "AuxNames" };
const std::string sSummaryNames{ "SummaryNames" };
const std::string sColumnNames{ "ColumnNames" };
const std::string sSeed{ "Seed" };
const std::string sSeedMachine{ "SeedMachine" };
const std::string sRandom{ "Random" };
const std::string sSampleSize{ "SampleSize" };
const std::string sConfigCount{ "ConfigCount" };
const std::string sFileList{ "FileList" };
const std::string sBootstrapList{ "BootstrapList" };
const std::string sBinSize{ "BinSize" };
const std::string sCovarSource{ "CovarSource" };
const std::string sCovarRebin{ "CovarRebin" };
const std::string sCovarSampleSize{ "CovarSampleSize" };
const std::string sCovarNumBoot{ "CovarNumBoot" };
const std::string sNE{ " != " };
const std::vector<std::string> sCorrSummaryNames{ "corr", "bias", "exp", "cosh" };
const double NaN{ std::nan( "" ) };

namespace Gamma
{
  const std::array<const std::string, nGamma> name{ // Long name, per Grid
    "MinusGamma5",
    "Gamma5",
    "MinusGammaT",
    "GammaT",
    "MinusGammaTGamma5",
    "GammaTGamma5",
    "MinusGammaX",
    "GammaX",
    "MinusGammaXGamma5",
    "GammaXGamma5",
    "MinusGammaY",
    "GammaY",
    "MinusGammaYGamma5",
    "GammaYGamma5",
    "MinusGammaZ",
    "GammaZ",
    "MinusGammaZGamma5",
    "GammaZGamma5",
    "MinusIdentity",
    "Identity",
    "MinusSigmaXT",
    "SigmaXT",
    "MinusSigmaXY",
    "SigmaXY",
    "MinusSigmaXZ",
    "SigmaXZ",
    "MinusSigmaYT",
    "SigmaYT",
    "MinusSigmaYZ",
    "SigmaYZ",
    "MinusSigmaZT",
    "SigmaZT",
  };
  const std::array<const std::string, nGamma> nameShort{ // My abbreviations
    "g-5",
    "g5",
    "g-T",
    "gT",
    "g-T5",
    "gT5",
    "g-X",
    "gX",
    "g-X5",
    "gX5",
    "g-Y",
    "gY",
    "g-Y5",
    "gY5",
    "g-Z",
    "gZ",
    "g-Z5",
    "gZ5",
    "g-I",
    "gI",
    "g-sXT",
    "gsXT",
    "g-sXY",
    "gsXY",
    "g-sXZ",
    "gsXZ",
    "g-sYT",
    "gsYT",
    "g-sYZ",
    "gsYZ",
    "g-sZT",
    "gsZT",
  };
  const std::string sGammaSpatial{ "gXYZ" };

  std::string NameShort( Algebra alg, const char * pszPrefix, const char * pszOpSuffix )
  {
    std::string sName;
    const bool bGotSuffix{ pszOpSuffix && *pszOpSuffix };
    if( alg != Algebra::Unknown || bGotSuffix )
    {
      if( pszPrefix && *pszPrefix )
        sName.append( pszPrefix );
      if( alg == Algebra::Spatial )
        sName.append( sGammaSpatial );
      else if( alg != Algebra::Unknown )
      {
        int idx = static_cast<int>( alg );
        sName.append( ( idx >= 0 && idx < nGamma ) ? nameShort[idx] : "?" );
      }
      if( bGotSuffix )
        sName.append( pszOpSuffix );
    }
    return sName;
  }

  std::ostream& operator<<(std::ostream& os, const Algebra &a)
  {
    return os << NameShort( a );
  }

  std::istream& operator>>(std::istream& is, Algebra &a)
  {
    std::string s;
    if( is >> s )
    {
      int i;
      for( i = 0; i < nGamma && !EqualIgnoreCase( s, nameShort[i] ); i++ )
        ;
      if( i < nGamma )
        a = static_cast<Algebra>( i );
      else if( EqualIgnoreCase( s, sGammaSpatial ) )
        a = Algebra::Spatial;
      else
        is.setstate( std::ios_base::failbit );
    }
    return is;
  }
};

// Default delimeters for the next couple of functions
extern const char szDefaultDelimeters[] = " \t,";

// Remove anything past the last delimeter from string, returning the removed part in suffix
// Return success / fail
bool ExtractSuffix( std::string &String, std::string &Suffix, const char * pszDelimeters )
{
  if( !pszDelimeters || !*pszDelimeters )
    pszDelimeters = szDefaultDelimeters;
  std::size_t NumDelims{ 0 };
  while( pszDelimeters[NumDelims] )
    NumDelims++;
  std::size_t Len{ String.length() };
  bool bFoundDelim{ false };
  while( !bFoundDelim && Len )
  {
    const char c{ String[--Len] };
    for( std::size_t i = 0; !bFoundDelim && i < NumDelims; i++ )
    {
      bFoundDelim = ( pszDelimeters[i] == c );
      if( bFoundDelim )
      {
        Suffix = String.substr( Len + 1 );
        Trim( Suffix );
        // Skip past multiple delimeters if they are all whitespace
        if( c == ' ' || c == '\t' || c == '\r' || c == '\n' )
        {
          while( Len && ( String[Len - 1] == ' ' || String[Len - 1] == '\t'
                         || String[Len - 1] == '\r' || String[Len - 1] == '\n' ) )
            --Len;
        }
        String.resize( Len );
      }
    }
  }
  return bFoundDelim;
}

// Remove the directory from the start of FileName (leave the trailing '/' in place)
std::string ExtractDirPrefix( std::string &FileName )
{
  std::string Dir;
  std::size_t pos{ FileName.find_last_of( '/' ) };
  if( pos != std::string::npos )
  {
    Dir = FileName.substr( 0, pos + 1 );
    FileName.erase( 0, pos + 1 );
  }
  return Dir;
}

// Split String into an array using specified delimeters

std::vector<std::string> Split( const std::string &String, const char * pszDelimeters )
{
  std::vector<std::string> a;
  if( !pszDelimeters || !*pszDelimeters )
    pszDelimeters = szDefaultDelimeters;
  std::size_t NumDelims{ 0 };
  while( pszDelimeters && pszDelimeters[NumDelims] )
    NumDelims++;
  const std::size_t Len{ String.length() };
  std::size_t Start{ 0 };
  while( Start < Len )
  {
    // Look for the next delimeter
    std::size_t Pos{ Start };
    bool bFoundDelim{ false };
    char c = 0;
    while( Pos < Len && !bFoundDelim )
    {
      c = String[Pos];
      for( std::size_t i = 0; !bFoundDelim && i < NumDelims; i++ )
        bFoundDelim = ( pszDelimeters[i] == c );
      if( !bFoundDelim )
        Pos++;
    }
    // Append this substring to list of items to return
    a.push_back( String.substr( Start, Pos - Start ) );
    // Skip past this delimeter
    Start = Pos + 1;
    // Skip past multiple delimeters if they are all whitespace
    if( c == ' ' || c == '\t' || c == '\r' || c == '\n' )
    {
      while( Start < Len && ( String[Start] == ' ' || String[Start] == '\t'
                             || String[Start] == '\r' || String[Start] == '\n' ) )
        ++Start;
    }
  }
  return a;
}

// Extract suffix, then split strings. Default delimeters '.' and '_' respectively
bool ExtractSuffixSplit( std::string &String, std::vector<std::string> &Suffii,
                        const char * pszStringDelim, const char * pszSuffixDelim )
{
  if( !pszStringDelim || !*pszStringDelim )
    pszStringDelim = Period.c_str();
  if( !pszSuffixDelim || !*pszSuffixDelim )
    pszSuffixDelim = Underscore.c_str();
  std::string Suffix;
  const bool bExtracted{ ExtractSuffix( String, Suffix, pszStringDelim ) };
  if( bExtracted )
    Suffii = Split( Suffix, pszSuffixDelim );
  return bExtracted;
}

// Zipper merge v1 and v2 if same size (otherwise just append)
std::vector<std::string> ZipperMerge( const std::vector<std::string> &v1, const std::vector<std::string> &v2 )
{
  std::vector<std::string> myFileList;
  // Merge the file lists
  if( v1.size() == v2.size() )
  {
    // Same length - merge the lists like a zipper
    myFileList.reserve( v1.size() + v2.size() );
    auto p2 = v2.begin();
    for( auto p1 = v1.begin(); p1 != v1.end(); ++p1, ++p2 )
    {
      myFileList.emplace_back( *p1 );
      myFileList.emplace_back( *p2 );
    }
  }
  else
  {
    myFileList = v1;
    myFileList.insert( myFileList.end(), v2.begin(), v2.end() );
  }
  return myFileList;
}

// Dump the environment to stdout, prefixed by optional message
void DumpEnv(int argc, const char * const *argv, const char * pStr )
{
  static const char sIndent2[]{ "    " };
  static const char * sIndent1{ sIndent2 + 2 };
  static const char sQuote[]{ "\"" };
  if( pStr )
    std::cout << pStr << std::endl;
  std::cout << sIndent1 << argc << " arguments:" << std::endl;
  for( int i = 0; i < argc; i++ )
  {
    std::cout << sIndent2 << "argv[" << i << "] = ";
    if( argv[i] )
      std::cout << sQuote << argv[i] << sQuote;
    else
      std::cout << "nullptr";
    std::cout << std::endl;
  }
}

// Does the specified file exist?
bool FileExists( const std::string& Filename )
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

// Make the ancestor directories leading up to last element
// NB: Don't test whether this worked, so that it if this is done simultaneously
// by many threads, it will still work
void MakeAncestorDirs( const std::string& Filename )
{
  bool WasSlash{ true };
  std::size_t SegmentStart = 0;
  for( std::size_t i = 0; i < Filename.length(); ++i )
  {
    const bool IsSlash{ Filename[i] == '/' };
    if( IsSlash )
    {
      if( !WasSlash )
      {
        // First slash after non-slash - try to make this bit
        WasSlash = true;
        const std::size_t SegmentLen{ i - SegmentStart };
        // Don't bother with '.' or '..'
        if( SegmentLen > 2 || Filename[i - SegmentLen] != '.'
           || ( SegmentLen == 2 && Filename[i - SegmentLen + 1] != '.' ) )
        {
          // Make this part of the name
          mkdir( Filename.substr( 0, i ).c_str(), 0777 );
        }
      }
    }
    else
    {
      if( WasSlash )
      {
        // First non-slash after slash - remember where this segment starts
        WasSlash = false;
        SegmentStart = i;
      }
    }
  }
}

// Wrapper for posix gethostname()
std::string GetHostName()
{
  char Buffer[256];
  const int BufLen{ sizeof( Buffer ) - 1 };
  if( gethostname( Buffer, BufLen ) )
    throw std::runtime_error( "gethostname() returned error " + std::to_string( errno ) );
  Buffer[BufLen] = 0;
  return std::string( Buffer );
}

template class MatrixView<double>;
template class MatrixView<float>;
template class MatrixView<std::complex<double>>;
template class MatrixView<std::complex<float>>;

template <typename T>
void VectorView<T>::MapRow( const Matrix<T> &m, std::size_t Row )
{
  if( m.size1 == 0 || m.size2 == 0 )
    throw std::runtime_error( "VectorView::MapRow() empty matrix" );
  if( Row >= m.size1 )
    throw std::runtime_error( "Row " + std::to_string( Row ) + " > " + std::to_string( m.size1 ) );
  Map( reinterpret_cast<Scalar *>( m.data ) + m.tda * Row, m.size2, 1 );
}

template <typename T>
void VectorView<T>::MapRow( const MatrixView<T> &m, std::size_t Row )
{
  if( m.size1() == 0 || m.size2() == 0 )
    throw std::runtime_error( "VectorView::MapRow() empty matrix" );
  if( Row >= m.size1() )
    throw std::runtime_error( "Row " + std::to_string( Row ) + " > " + std::to_string( m.size1() ) );
  Map( m.data() + m.tda() * Row, m.size2(), 1 );
}

template <typename T>
void VectorView<T>::MapColumn( const Matrix<T> &m, std::size_t Column )
{
  if( m.size1 == 0 || m.size2 == 0 )
    throw std::runtime_error( "VectorView::MapColumn() empty matrix" );
  if( Column >= m.size2 )
    throw std::runtime_error( "Column " + std::to_string( Column ) + " > " + std::to_string( m.size2 ) );
  Map( reinterpret_cast<Scalar *>( m.data ) + Column, m.size1, m.tda );
}

template <typename T>
void VectorView<T>::MapColumn( const MatrixView<T> &m, std::size_t Column )
{
  if( m.size1() == 0 || m.size2() == 0 )
    throw std::runtime_error( "VectorView::MapColumn() empty matrix" );
  if( Column >= m.size2() )
    throw std::runtime_error( "Column " + std::to_string( Column ) + " > " + std::to_string( m.size2() ) );
  Map( m.data() + Column, m.size1(), m.tda() );
}

template class VectorView<double>;
template class VectorView<float>;
template class VectorView<std::complex<double>>;
template class VectorView<std::complex<float>>;

// Read a bootstrap replica from an HDF5 group
template <typename T>
void BootRep<T>::Read( ::H5::Group &g, const std::string &Name )
{
  H5::ReadVector( g, Name + s_C, Central );
  H5::ReadMatrix( g, Name, Replica );
}

// Write a bootstrap replica to an HDF5 group
template <typename T>
void BootRep<T>::Write( ::H5::Group &g, const std::string &Name ) const
{
  H5::WriteVector( g, Name + s_C, Central );
  H5::WriteMatrix( g, Name, Replica );
}

HotellingDist::HotellingDist( unsigned int p_, unsigned int m_ )
: p{p_}, m{m_}, Nu{m - p + 1}, Factor{ Nu / ( static_cast<double>( p ) * m ) }
{
  if( !Usable( p, m ) )
    throw std::runtime_error( "T^2 distribution m (degrees of freedom of covariance matrix) "
                             + std::to_string( m ) + " < p (fit degrees of freedom) " + std::to_string( p ) );
}

double HotellingDist::operator()( double TestStatistic ) const
{
  const double ModifiedStatistic{ TestStatistic * Factor };
  const double FDistQ{ gsl_cdf_fdist_Q( ModifiedStatistic, p, Nu ) };
#ifdef _DEBUG
  const double FDistP{ gsl_cdf_fdist_P( ModifiedStatistic, p, Nu ) };
  std::cout << "t^2=" << TestStatistic << ", Factor=" << Factor
            << ", t^2 * Factor=" << ModifiedStatistic << ", p [dof]=" << p
            << ", m [SampleSize-1]=" << m << ", Nu [m - p + 1]=" << Nu << Common::NewLine
            << "gsl_cdf_fdist_Q( t^2 * Factor, p, Nu )=" << FDistQ << Common::NewLine
            << "gsl_cdf_fdist_P( t^2 * Factor, p, Nu )=" << FDistP << Common::NewLine
            << "gsl_cdf_fdist_Q + gsl_cdf_fdist_P = " << ( FDistQ + FDistP ) << Common::NewLine;
#endif
  return FDistQ;
}

// Sort the list of values, then extract the lower and upper 68th percentile error bars
template<typename T>
void ValWithEr<T>::Get( T Central_, std::vector<T> &Data, std::size_t Count )
{
  assert( Data.size() >= Count && "ValWithErr<T>::Get() data too small" );
  Central = Central_;
  if( Count == 0 )
  {
    if( Data.size() == 0 )
    {
      // I wasn't trying to take an estimate over many samples
      Min = Central;
      Max = Central;
      Low = Central;
      High  = Central;
      Check = 1;
    }
    else
    {
      // Tried to take an estimate over many samples, but all values NaN
      Min = NaN;
      Max = NaN;
      Low = NaN;
      High  = NaN;
      Check = 0;
    }
    return;
  }
  const typename std::vector<T>::iterator itStart{ Data.begin() };
  const typename std::vector<T>::iterator itEnd  { itStart + Count };
  std::sort( itStart, itEnd );
  const double OneSigma{ ( 1. - 0.682689492137086 ) * 0.5 }; // https://en.wikipedia.org/wiki/68–95–99.7_rule
  std::size_t Index = OneSigma * Count + 0.5;
  if( Index >= Count )
    Index = Count - 1;
  Min  = Data[0];
  Max  = Data[Count - 1];
  Low  = Data[Index];
  High = Data[Count - 1 - Index];
  Check = static_cast<T>( Count ) / Data.size();
}

template class ValWithEr<float>;
template class ValWithEr<double>;
template class ValWithEr<long double>;

const Momentum p0(0,0,0);
const std::string Momentum::DefaultPrefix{ "p" };
const std::string Momentum::SinkPrefix{ DefaultPrefix + "s" };
const std::string Momentum::SquaredSuffix{ "2" };

std::string Momentum::p2_string( const std::string &separator, const std::string Prefix ) const
{
  std::string s{ separator };
  s.append( Prefix );
  s.append( SquaredSuffix );
  s.append( separator );
  s.append( std::to_string( p2() ) );
  return s;
}

std::string Momentum::to_string( const std::string &separator, bool bNegative ) const
{
  return std::to_string( bNegative ? -x : x ) + separator
  + std::to_string( bNegative ? -y : y ) + separator
  + std::to_string( bNegative ? -z : z );
}

std::string Momentum::to_string4d( const std::string &separator, bool bNegative ) const
{
  return to_string( separator, bNegative ) + separator + "0";
}

std::regex Momentum::MakeRegex( const std::string &MomName )
{
  std::string sPattern( Common::Underscore );
  sPattern.append( MomName );
  sPattern.append( "_(-?[0-9]+)_(-?[0-9]+)_(-?[0-9]+)" );
  return std::regex( sPattern );
}

// Strip out momentum info from string if present
bool Momentum::Extract( std::string &Prefix, const std::string &MomName, bool IgnoreSubsequentZeroNeg )
{
  bool bGotMomentum = false;
  std::smatch match;
  const std::regex Pattern{ MakeRegex( MomName ) };
  while( std::regex_search( Prefix, match, Pattern ) )
  {
    int px{ std::stoi( match[1] ) };
    int py{ std::stoi( match[2] ) };
    int pz{ std::stoi( match[3] ) };
    if( !bGotMomentum )
    {
      bGotMomentum = true;
      x = px;
      y = py;
      z = pz;
    }
    else if( x != px || y != py || z != pz )
    {
      if( IgnoreSubsequentZeroNeg && !(*this) ) // i.e. previous momentum was zero
      {
        x = px;
        y = py;
        z = pz;
      }
      else if(IgnoreSubsequentZeroNeg && ((px==0 && py==0 && pz==0) || (px==-x && py==-y && pz==-z)))
        ;
      else
      {
        static const std::string Sep{ "," };
        static const std::string m1{ to_string( Sep ) };
        x = px;
        y = py;
        z = pz;
        throw std::runtime_error( "Multiple momenta: " + m1 + " != " + to_string( Sep ) );
      }
    }
    const std::string sSuffix{ match.suffix() };
    Prefix = match.prefix();
    Prefix.append( sSuffix );
  }
  return bGotMomentum;
}

void Momentum::Replace( std::string &s, const std::string &MomName, bool bNegative ) const
{
  std::string Replace( Common::Underscore );
  Replace.append( MomName );
  Replace.append( Common::Underscore );
  Replace.append( to_string( Common::Underscore, bNegative ) );
  std::smatch match;
  const std::regex Pattern{ MakeRegex( MomName ) };
  bool bFound{ false };
  std::string Search{ std::move( s ) };
  while( std::regex_search( Search, match, Pattern ) )
  {
    s.append( match.prefix() );
    s.append( Replace );
    Search = match.suffix();
    bFound = true;
  }
  if( !bFound )
  {
    std::ostringstream ss;
    ss << "momentum " << MomName << " not found in " << s;
    throw std::runtime_error( ss.str() );
  }
  s.append( Search );
}

// Take a squared momentum and make an approximate 3d momentum from it
Momentum FileNameMomentum::FromSquared( const int p2Original )
{
  int p2{ p2Original };
  Momentum p;
  const bool bNegative{ p2 < 0 };
  if( bNegative )
    p2 = - p2;
  for( int i = 0; i < 3 && p2; ++i )
  {
    int root = 1;
    while( ( root + 1 ) * ( root + 1 ) <= p2 )
      root++;
    p2 -= root * root;
    switch( i )
    {
      case 0:
        p.x = root;
        break;
      case 1:
        p.y = root;
        break;
      case 2:
        p.z = root;
        break;
    }
  }
  if( p2 )
    throw std::runtime_error( std::to_string( p2Original ) + " can't be expressed as a three-momentum" );
  if( bNegative )
    p.x = -p.x;
  return p;
}

std::string FileNameMomentum::FileString( const std::string &separator ) const
{
  if( bp2 )
    return p2_string( separator, Name );
  std::string s{ separator };
  s.append( Name );
  s.append( separator );
  s.append( to_string( separator ) );
  return s;
}

std::ostream& operator<<( std::ostream& os, const Momentum &p )
{
  return os << p.to_string( Underscore );
}

std::istream& operator>>( std::istream& is, Momentum &p )
{
  char c = 0;
  if(!(is>>p.x && is.get(c) && c==Underscore[0] && is>>p.y && is.get(c) && c==Underscore[0] && is>>p.z))
    is.setstate( std::ios_base::failbit );
  return is;
}

// These are the attributes I like to use in my filenames
void FileNameAtt::Parse( const std::string &Filename_, std::vector<std::string> * pOpNames,
                         const std::vector<std::string> * pIgnoreMomenta,
                         const std::vector<std::string> * pIgnoreRegEx )
{
  clear();
  Filename = Filename_;
  std::size_t pos = Filename.find_last_of( '/' );
  if( pos == std::string::npos )
    pos = 0;
  else
    Dir = Filename.substr( 0, ++pos );
  Base = Filename.substr( pos );
  NameNoExt = Base;
  int i = 0;
  while( i < 3 && ( pos = Base.find_last_of( '.' ) ) != std::string::npos )
  {
    switch( i )
    {
      case 0:
        Ext = Base.substr( pos + 1 );
        break;
      case 1:
        SeedString = Base.substr( pos + 1 );
        try {
          Seed = FromString<unsigned int>( SeedString );
          bSeedNum = true;
        } catch(...) {
          std::cout << "Ignoring invalid seed in " << Filename << std::endl;
        }
        break;
      case 2:
        Type = Base.substr( pos + 1 );
        break;
    }
    Base.resize(pos);
    if( i == 0 )
      NameNoExt = Base;
    i++;
  }
  // Now see whether we can extract operator names
  if( pOpNames )
    ParseOpNames( *pOpNames, 2, pIgnoreMomenta, pIgnoreRegEx );
  else
    ParseShort( pIgnoreMomenta, pIgnoreRegEx );
}

void FileNameAtt::ParseShort( const std::vector<std::string> * pIgnoreMomenta,
                              const std::vector<std::string> * pIgnoreRegEx )
{
  // Remove momenta we are going to ignore in the filename
  Momentum pIgnore;
  if( pIgnoreMomenta )
  {
    for( const std::string &s : *pIgnoreMomenta )
    {
      while( pIgnore.Extract( Base, s ) ) // There might be more than one copy
        ; //momIgnore.emplace_back( s, pIgnore );
    }
  }
  if( pIgnoreRegEx )
  {
    for( const std::string &s : *pIgnoreRegEx )
    {
      // Remove strings from the filename. No need to expose these strings (for now)
      while( ExtractToken( Base, s ) )
        ;
    }
  }
  BaseShort = Base;
  // Now get the remaining momenta and save them
  std::smatch match;
  std::regex Pattern{ Momentum::MakeRegex( "([pP][[:alnum:]]*)" ) };
  for( int pLoop = 0; pLoop < 2; ++pLoop )
  {
    while( std::regex_search( BaseShort, match, Pattern ) )
    {
      // Extract the momentum
      const std::string sMom{ match[1] };
      std::vector<FileNameMomentum> fnp;
      fnp.reserve( 1 );
      if( pLoop == 0 )
      {
        fnp.emplace_back( sMom, std::stoi( match[2] ), std::stoi( match[3] ), std::stoi( match[4] ) );
      }
      else
      {
        fnp.emplace_back( sMom, std::stoi( match[2] ) );
      }
      const std::string sSuffix{ match.suffix() };
      BaseShort = match.prefix();
      BaseShort.append( sSuffix );
      // Now see whether we already have it
      auto it = p.find( sMom );
      if( it == p.end() )
        p.emplace( std::make_pair( sMom, fnp[0] ) );
      else if( ! ( it->second == fnp[0] ) )
      {
        std::stringstream ss;
        ss << "Repeated momentum " << sMom << CommaSpace << it->second << " != " << pIgnore;
        throw std::runtime_error( ss.str() );
      }
    }
    // Second time through the loop, we're looking for p^2
    if( pLoop == 0 )
    {
      const std::string sPattern{ "_([pP][[:alnum:]]*)" + Momentum::SquaredSuffix + "_(-?[0-9]+)" };
      Pattern = sPattern;
    }
  }
  // Extract other attributes from filename
  ExtractTimeslice( BaseShort, bGotTimeslice, Timeslice );
  ExtractDeltaT( BaseShort, bGotDeltaT, DeltaT );
  Gamma.clear();
  Gamma = ExtractGamma( BaseShort );
}

// Parse operator names from the end of Base, building up a list of all the operator names
// NB: Because I parse from the end, op[0] is the last operator, op[1] second last, etc
std::vector<std::string> FileNameAtt::ParseOpNames( int NumOps,
                                                    const std::vector<std::string> * pIgnoreMomenta,
                                                    const std::vector<std::string> * pIgnoreRegEx )
{
  static const char Sep[] = "_.";
  constexpr std::size_t NumSeps{ sizeof( Sep ) / sizeof( Sep[0] ) - 1 };
  std::vector<std::string> o;
  o.reserve( NumOps );
  std::size_t LastPos = Base.length();
  int i = 0;
  for( ; LastPos && i < NumOps; ++i )
  {
    // Search for the separators - on the last one, check for a period as well as underscore
    LastPos--; // Skip past the last separator
    std::size_t pos = Base.find_last_of( Sep, LastPos, NumSeps );
    if( pos == std::string::npos )
    {
      o.push_back( Base.substr( 0, LastPos + 1 ) );
      LastPos = 0;
    }
    else if( i != NumOps - 1 && Base[pos] == '.' )
    {
      break;
    }
    else
    {
      o.push_back( Base.substr( pos + 1, LastPos - pos ) );
      LastPos = pos;
    }
  }
  if( i != NumOps )
    throw std::runtime_error( "Can't extract " + std::to_string( NumOps ) + " operators from " + Base );
  if( i )
  {
    Base.resize( LastPos ); // Shorten the base
  }
  ParseShort( pIgnoreMomenta, pIgnoreRegEx );
  return o;
}

void FileNameAtt::ParseOpNames( std::vector<std::string> &OpNames, int NumOps,
                                const std::vector<std::string> * pIgnoreMomenta,
                                const std::vector<std::string> * pIgnoreRegEx )
{
  std::vector<std::string> Ops{ ParseOpNames( NumOps, pIgnoreMomenta, pIgnoreRegEx ) };
  op.clear();
  op.reserve( NumOps );
  for( std::string &s : Ops )
  {
    int iOp = 0;
    while( iOp < OpNames.size() && !EqualIgnoreCase( s, OpNames[iOp] ) )
      iOp++;
    if( iOp == OpNames.size() )
      OpNames.emplace_back( std::move( s ) );
    op.push_back( iOp );
  }
}

void FileNameAtt::ParseExtra( unsigned int MaxElements, const std::vector<std::string> * pIgnoreMomenta,
                              const std::vector<std::string> * pIgnoreRegEx )
{
  std::size_t pos;
  while( MaxElements-- && ( pos = Base.find_last_of( '.' ) ) != std::string::npos )
  {
    Extra.push_back( Base.substr( pos + 1 ) );
    Base.resize( pos );
  }
  ParseShort( pIgnoreMomenta, pIgnoreRegEx );
}

// Append the extra info to the string
void FileNameAtt::AppendExtra( std::string &s, int Last, int First ) const
{
  int Num{ static_cast<int>( Extra.size() ) };
  if( First >= 0 && First < Num )
    Num = First;
  while( Num-- > Last )
  {
    s.append( 1, '.' );
    s.append( Extra[Num] );
  }
}

std::string FileNameAtt::GetBaseExtra( int Last, int First ) const
{
  std::string s{ Base };
  AppendExtra( s, Last, First );
  return s;
}

/*std::string FileNameAtt::GetBaseShortExtra( int Last, int First ) const
{
  std::string s{ BaseShort };
  AppendExtra( s, Last, First );
  return s;
}*/

// Make a new name based on this one, overriding specified elements
std::string FileNameAtt::DerivedName( const std::string &Suffix, const std::string &Snk, const std::string &Src,
                                      const std::string &Ext ) const
{
  std::string s{ Base };
  s.append( 1, '_' );
  s.append( Snk );
  s.append( 1, '_' );
  s.append( Src );
  AppendExtra( s );
  s.append( 1, '.' );
  s.append( Type );
  s.append( 1, '.' );
  s.append( SeedString );
  s.append( 1, '.' );
  s.append( Ext );
  return s;
}

const FileNameMomentum &FileNameAtt::GetMomentum( const std::string &Name ) const
{
  auto it = p.find( Name );
  if( it == p.end() )
    it = p.find( Name + Momentum::SquaredSuffix );
  if( it == p.end() )
    throw std::runtime_error( "Momentum " + Name + " not available" );
  return it->second;
}

const FileNameMomentum &FileNameAtt::GetFirstNonZeroMomentum() const
{
  if( p.empty() )
    throw std::runtime_error( "No momenta available" );
  for( auto it = p.begin(); it != p.end(); ++it )
    if( it->second )
      return it->second;
  return p.begin()->second;
}

void FileNameAtt::AppendMomentum( std::string &s, const FileNameMomentum &fnp, const std::string &Name ) const
{
  if( fnp.bp2 )
    s.append( fnp.p2_string( Underscore, Name ) );
  else
  {
    s.append( Underscore );
    s.append( Name );
    s.append( fnp.to_string( Underscore ) );
  }
}

// Make a filename "Base.Type.seed.Ext"
std::string MakeFilename(const std::string &Base, const std::string &Type, SeedType Seed, const std::string &Ext)
{
  const char Sep = '.';
  std::string s{ Base };
  s.append( 1, Sep );
  s.append( Type );
  s.append( 1, Sep );
  s.append( std::to_string( Seed ) );
  s.append( 1, Sep );
  s.append( Ext );
  return s;
}

// If present, remove Token from a string. Return true if removed
bool ExtractToken( std::string &Prefix, const std::string &Token )
{
  bool bExtracted{ false };
  std::smatch match;
  const std::regex pattern{ "(^|_)" + Token + "(_|$)" };
  while( std::regex_search( Prefix, match, pattern  ) )
  {
    if( bExtracted )
      throw std::runtime_error( "Multiple " + Token + " tokens in " + Prefix );
    bExtracted = true;
    std::string s{ match.prefix() };
    if( match[1].length() && match[2].length() )
      s.append( match[2] );
    s.append( match.suffix() );
    Prefix = s;
  }
  return bExtracted;
}

// If present, remove integer preceded by Token from a string
void ExtractInteger( std::string &Prefix, bool &bHasValue, int &Value, const std::string Token )
{
  std::smatch match;
  const std::regex pattern{ "_" + Token + "_?([0-9]+)" };
  while( std::regex_search( Prefix, match, pattern  ) )
  {
    int ThisValue = std::stoi( match[1] );
    if( !bHasValue )
    {
      Value = ThisValue;
      bHasValue = true;
    }
    else if( ThisValue != Value )
      throw std::runtime_error( "Multiple regex " + Token + ": " + std::to_string(Value) + " and " + std::to_string(ThisValue) );
    const std::string sSuffix{ match.suffix() };
    Prefix = match.prefix();
    Prefix.append( sSuffix );
  }
}

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &Prefix, bool &bHasTimeslice, int &Timeslice )
{
  ExtractInteger( Prefix, bHasTimeslice, Timeslice, "[tT]" );
}

// Strip out DeltaT from a string if present
void ExtractDeltaT( std::string &Prefix, bool &bHasDeltaT, int &DeltaT )
{
  ExtractInteger( Prefix, bHasDeltaT, DeltaT, "[dD][tT]" );
}

// Append DeltaT to string
void AppendDeltaT( std::string &s, int DeltaT )
{
  s.append( "_dt_" );
  s.append( std::to_string( DeltaT ) );
}

// Remove any gammas from Prefix
std::vector<Gamma::Algebra> ExtractGamma( std::string &Prefix )
{
  static const std::regex GammaPattern{ "_(g[^_]+)" };

  std::vector<Gamma::Algebra> v;
  std::smatch match;
  std::string Search{ std::move( Prefix ) };
  while( std::regex_search( Search, match, GammaPattern ) )
  {
    Prefix.append( match.prefix() );
    Gamma::Algebra a;
    std::stringstream ss( match[1] );
    if( ss >> a && ( ss.eof() || ( ss >> std::ws && ss.eof() ) ) )
      v.push_back( a );
    else
    {
      Prefix.append( 1, '_' );
      Prefix.append( match[1] );
    }
    Search = match.suffix();
  }
  Prefix.append( Search );
  return v;
}

// Append Gamma to string
void AppendGamma( std::string &s, Gamma::Algebra g, const std::string &Sep )
{
  std::ostringstream os;
  os << Sep << g;
  s.append( os.str() );
}

// Remove any gammas from Prefix
void ReplaceGamma( std::string &Prefix, Gamma::Algebra gFrom, Gamma::Algebra gTo )
{
  std::regex GammaPattern;
  {
    std::ostringstream ss;
    ss << Common::Underscore << gFrom << "(_|$)";
    GammaPattern = ss.str();
  }
  std::string GammaReplace;
  {
    std::ostringstream ss;
    ss << Common::Underscore << gTo;
    GammaReplace = ss.str();
  }
  std::smatch match;
  std::string Search{ std::move( Prefix ) };
  while( std::regex_search( Search, match, GammaPattern ) )
  {
    Prefix.append( match.prefix() );
    Prefix.append( GammaReplace );
    Prefix.append( match[1] );
    Search = match.suffix();
  }
  Prefix.append( Search );
}

// Make the same HDF5 complex type Grid uses
template<typename T> static ::H5::CompType MakeComplex()
{
  ::H5::CompType myComplex( sizeof( std::complex<T> ) );
  myComplex.insertMember("re", 0 * sizeof(T), H5::Equiv<T>::Type);
  myComplex.insertMember("im", 1 * sizeof(T), H5::Equiv<T>::Type);
  return myComplex;
}

// Make an HDF5 type representing a value with an error
template<typename T> static ::H5::CompType MakeValWithErOldV1()
{
  ::H5::CompType myType( sizeof( ValWithErOldV1<T> ) );
  myType.insertMember("Central", offsetof(ValWithErOldV1<T>, Central), H5::Equiv<T>::Type);
  myType.insertMember("Low",     offsetof(ValWithErOldV1<T>, Low    ), H5::Equiv<T>::Type);
  myType.insertMember("High",    offsetof(ValWithErOldV1<T>, High   ), H5::Equiv<T>::Type);
  myType.insertMember("Check",   offsetof(ValWithErOldV1<T>, Check  ), H5::Equiv<T>::Type);
  return myType;
}

// Make an HDF5 type representing a value with an error
template<typename T> static ::H5::CompType MakeValWithEr()
{
  ::H5::CompType myType( sizeof( ValWithEr<T> ) );
  myType.insertMember("Min",     offsetof(ValWithEr<T>, Min    ), H5::Equiv<T>::Type);
  myType.insertMember("Low",     offsetof(ValWithEr<T>, Low    ), H5::Equiv<T>::Type);
  myType.insertMember("Central", offsetof(ValWithEr<T>, Central), H5::Equiv<T>::Type);
  myType.insertMember("High",    offsetof(ValWithEr<T>, High   ), H5::Equiv<T>::Type);
  myType.insertMember("Max",     offsetof(ValWithEr<T>, Max    ), H5::Equiv<T>::Type);
  myType.insertMember("Check",   offsetof(ValWithEr<T>, Check  ), H5::Equiv<T>::Type);
  return myType;
}

// Make an HDF5 type representing a config number and count
static ::H5::CompType MakeConfigCount()
{
  ::H5::CompType myType( sizeof( ConfigCount ) );
  myType.insertMember("Config", offsetof(ConfigCount, Config), ::H5::PredType::NATIVE_UINT32);
  myType.insertMember("Count",  offsetof(ConfigCount, Count ), ::H5::PredType::NATIVE_UINT32);
  return myType;
}

const ::H5::PredType& H5::Equiv<float>      ::Type{ ::H5::PredType::NATIVE_FLOAT };
const ::H5::PredType& H5::Equiv<double>     ::Type{ ::H5::PredType::NATIVE_DOUBLE };
const ::H5::PredType& H5::Equiv<long double>::Type{ ::H5::PredType::NATIVE_LDOUBLE };
const ::H5::PredType& H5::Equiv<int>        ::Type{ ::H5::PredType::NATIVE_INT };
const ::H5::StrType   H5::Equiv<std::string>::Type{ ::H5::PredType::C_S1, H5T_VARIABLE };
const ::H5::StrType&  H5::Equiv<char *>     ::Type{   H5::Equiv<std::string>::Type };
const ::H5::PredType& H5::Equiv<std::uint_fast32_t>::Type{ sizeof( std::uint_fast32_t ) == 4
                                                    ? ::H5::PredType::STD_U32LE : ::H5::PredType::STD_U64LE };
const ::H5::CompType  H5::Equiv<std::complex<float>>      ::Type{ MakeComplex<float>() };
const ::H5::CompType  H5::Equiv<std::complex<double>>     ::Type{ MakeComplex<double>() };
const ::H5::CompType  H5::Equiv<std::complex<long double>>::Type{ MakeComplex<long double>() };
const ::H5::CompType  H5::Equiv<ValWithEr<float>>      ::Type{ MakeValWithEr<float>() };
const ::H5::CompType  H5::Equiv<ValWithEr<double>>     ::Type{ MakeValWithEr<double>() };
const ::H5::CompType  H5::Equiv<ValWithEr<long double>>::Type{ MakeValWithEr<long double>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<float>>      ::Type{ MakeValWithErOldV1<float>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<double>>     ::Type{ MakeValWithErOldV1<double>() };
const ::H5::CompType  H5::Equiv<ValWithErOldV1<long double>>::Type{ MakeValWithErOldV1<long double>() };
const ::H5::CompType  H5::Equiv<ConfigCount> ::Type{ MakeConfigCount() };

// Open the specified HDF5File and group
void H5::OpenFileGroup(::H5::H5File &f, ::H5::Group &g, const std::string &FileName, const char *PrintPrefix,
                       std::string * pGroupName, unsigned int flags)
{
  const bool bFindGroupName{ !pGroupName || pGroupName->empty() };
  std::string localGroupName;
  if( !pGroupName )
    pGroupName = &localGroupName;
  f.openFile( FileName, flags );
  g = f.openGroup( bFindGroupName ? std::string("/") : *pGroupName );
  if( bFindGroupName ) {
    *pGroupName = GetFirstGroupName( g );
    g = g.openGroup( *pGroupName );
  }
  if( PrintPrefix )
    std::cout << PrintPrefix << FileName << " (" << *pGroupName << ")\n";
}

// Get first groupname from specified group
std::string H5::GetFirstGroupName( ::H5::Group & g )
{
  hsize_t n = g.getNumObjs();
  for( hsize_t i = 0; i < n; ++i ) {
    H5G_obj_t t = g.getObjTypeByIdx( i );
    if( t == H5G_GROUP )
      return g.getObjnameByIdx( i );
  }
  return std::string();
}

// Read the gamma algebra attribute string and make sure it's valid
Gamma::Algebra H5::ReadGammaAttribute( ::H5::Group &g, const char * pAttName )
{
  std::string sGamma;
  ::H5::Attribute a = g.openAttribute( pAttName );
  ::H5::StrType s = a.getStrType();
  a.read( s, sGamma );
  a.close();
  for( int idxGamma = 0; idxGamma < Gamma::nGamma; idxGamma++ )
    if( EqualIgnoreCase( sGamma, Gamma::name[idxGamma] ) )
      return static_cast<Gamma::Algebra>( idxGamma );
  throw ::H5::Exception( "Common::ReadGammaAttribute", "Invalid gamma algebra string" );
}

template<> void H5::ReadStringsHelper<::H5::Attribute>( const ::H5::Attribute &a, const ::H5::StrType &aType, char * * MDString )
{
  a.read( aType, ( void * ) MDString );
}
template<> void H5::ReadStringsHelper<::H5::DataSet>( const ::H5::DataSet &ds, const ::H5::StrType &aType, char * * MDString )
{
  ds.read( ( void * ) MDString, aType );
}

// Make a multi-dimensional string attribute
void H5::WriteAttribute( ::H5::Group &g, const std::string &AttName, const std::vector<std::string> &vs)
{
  std::unique_ptr<char *[]> RawArray( new char * [vs.size()] );
  const char * * const MDString{ const_cast<const char * *>( RawArray.get() ) };
  for( std::size_t i = 0; i < vs.size(); i++ )
    MDString[i] = vs[i].c_str();
  const hsize_t NDimension{ vs.size() };
  ::H5::DataSpace dsN( 1, &NDimension );
  ::H5::Attribute a{ g.createAttribute( AttName, Equiv<std::string>::Type, dsN ) };
  a.write( Equiv<std::string>::Type, MDString );
  a.close();
}

// Make a multi-dimensional string attribute
void H5::WriteStringData( ::H5::Group &g, const std::string &DSName, const std::vector<std::string> &vs)
{
  std::unique_ptr<char *[]> RawArray( new char * [vs.size()] );
  const char * * const MDString{ const_cast<const char * *>( RawArray.get() ) };
  for( std::size_t i = 0; i < vs.size(); i++ )
    MDString[i] = vs[i].c_str();
  const hsize_t NDimension{ vs.size() };
  ::H5::DataSpace dsN( 1, &NDimension );
  ::H5::DataSet ds{ g.createDataSet( DSName, Equiv<std::string>::Type, dsN ) };
  ds.write( MDString, Equiv<std::string>::Type );
  ds.close();
}

const std::string sUnknown{ "unknown" };

const std::string sReality{ "Reality" };
const std::string Reality_TextEquiv_Real{ "real" };
const std::string Reality_TextEquiv_Imaginary{ "imaginary" };

std::ostream& operator<<(std::ostream& os, const Reality &reality)
{
  switch( reality )
  {
    case Reality::Real:
      os << Reality_TextEquiv_Real;
      break;
    case Reality::Imag:
      os << Reality_TextEquiv_Imaginary;
      break;
    default:
      os << sUnknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Reality &reality)
{
  reality = Reality::Unknown;
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Reality_TextEquiv_Real ) )
      reality = Reality::Real;
    else if( EqualIgnoreCase( s, Reality_TextEquiv_Imaginary ) )
      reality = Reality::Imag;
    else if( !EqualIgnoreCase( s, sUnknown ) )
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

const std::string sParity{ "Parity" };
const std::string Parity_TextEquiv_Even{ "even" };
const std::string Parity_TextEquiv_Odd{ "odd" };

std::ostream& operator<<(std::ostream& os, const Parity &parity)
{
  switch( parity )
  {
    case Parity::Even:
      os << Parity_TextEquiv_Even;
      break;
    case Parity::Odd:
      os << Parity_TextEquiv_Odd;
      break;
    default:
      os << sUnknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Parity &parity)
{
  parity = Parity::Unknown;
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Parity_TextEquiv_Even ) )
      parity = Parity::Even;
    else if( EqualIgnoreCase( s, Parity_TextEquiv_Odd ) )
      parity = Parity::Odd;
    else if( !EqualIgnoreCase( s, sUnknown ) )
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

const std::string sSign{ "Sign" };
const std::string Sign_TextEquiv_Positive{ "positive" };
const std::string Sign_TextEquiv_Negative{ "negative" };

std::ostream& operator<<(std::ostream& os, const Sign &sign)
{
  switch( sign )
  {
    case Sign::Positive:
      os << Sign_TextEquiv_Positive;
      break;
    case Sign::Negative:
      os << Sign_TextEquiv_Negative;
      break;
    default:
      os << sUnknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Sign &sign)
{
  sign = Sign::Unknown;
  std::string s;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Sign_TextEquiv_Positive ) )
      sign = Sign::Positive;
    else if( EqualIgnoreCase( s, Sign_TextEquiv_Negative ) )
      sign = Sign::Negative;
    else if( !EqualIgnoreCase( s, sUnknown ) )
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

void GenerateRandom( std::vector<fint> &Random, SeedType Seed, std::size_t NumBoot, std::size_t NumSamples )
{
  const std::size_t RandomLen{ NumSamples * NumBoot };
  if( RandomLen / NumSamples != NumBoot || NumSamples > std::numeric_limits<fint>::max() )
    throw std::runtime_error( "Too many bootstrap replicas " + std::to_string( NumBoot )
                             + " (with " + std::to_string( NumSamples ) + " per replica)" );
  Random.resize( RandomLen );
  std::mt19937                        engine( Seed );
  std::uniform_int_distribution<fint> random( 0, static_cast<fint>( NumSamples - 1 ) );
  for( std::size_t i = 0; i < RandomLen; ++i )
    Random[i] = random( engine );
}

const std::array<std::string, 3> aSampleSource{ "Binned", "Raw", "Bootstrap" };

std::ostream& operator<<( std::ostream& os, SampleSource sampleSource )
{
  const int enumIdx{ static_cast<int>( sampleSource ) };
  if( enumIdx >= 0 && enumIdx < aSampleSource.size() )
    return os << aSampleSource[enumIdx];
  return os << "Unknown" << std::to_string( enumIdx );
}

std::istream& operator>>( std::istream& is, SampleSource &sampleSource )
{
  std::string sEnum;
  if( is >> sEnum )
  {
    const int enumIdx{ IndexIgnoreCase( aSampleSource, sEnum ) };
    if( enumIdx != aSampleSource.size() )
    {
      sampleSource = static_cast<SampleSource>( enumIdx );
      return is;
    }
  }
  throw std::runtime_error( "SampleSource \"" + sEnum + "\" unrecognised" );
}

template <typename T>
void Model<T>::ReadAttributes( ::H5::Group &g )
{
  Base::ReadAttributes( g );
  ::H5::Attribute a;
  a = g.openAttribute(sNumExponents);
  a.read( ::H5::PredType::NATIVE_INT, &NumExponents );
  a.close();
  a = g.openAttribute(sNumFiles);
  a.read( ::H5::PredType::NATIVE_INT, &NumFiles );
  a.close();
  a = g.openAttribute(sTI);
  a.read( ::H5::PredType::NATIVE_INT, &ti );
  a.close();
  a = g.openAttribute(sTF);
  a.read( ::H5::PredType::NATIVE_INT, &tf );
  a.close();
  a = g.openAttribute(sDoF);
  a.read( ::H5::PredType::NATIVE_INT, &dof );
  a.close();
  a = g.openAttribute(sOperators);
  OpNames = H5::ReadStrings( a );
  a.close();
  std::int8_t i8;
  a = g.openAttribute(sFactorised);
  a.read( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  Factorised = ( i8 != 0 );
  a = g.openAttribute(sCovarFrozen);
  a.read( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  CovarFrozen = ( i8 != 0 );
  try
  {
    H5::ReadMatrix( g, sCovarianceIn+s_C, CovarIn );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovariance+s_C, Covar );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelation+s_C, Correl );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationCholesky+s_C, CorrelCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceInv+s_C, CovarInv );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationInvCholesky+s_C, CorrelInvCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceInvCholesky+s_C, CovarInvCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    StdErrorMean.Read( g, sStdErrorMean );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    FitInput.Read( g, sFitInput );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    ModelPrediction.Read( g, sModelPrediction );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    ErrorScaled.Read( g, sErrorScaled );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute( sCovarSource );
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    std::stringstream ss( s );
    ss >> CovarSource;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
    CovarSource = SS::Binned; // Default if missing
  }
  CovarRebin.clear();
  try
  {
    a = g.openAttribute( sCovarRebin );
    ::H5::DataSpace dsp = a.getSpace();
    const int rank{ dsp.getSimpleExtentNdims() };
    if( rank != 1 )
      throw std::runtime_error( sCovarRebin + " dimensions " + std::to_string( rank ) + ", expecting 1" );
    hsize_t Num;
    dsp.getSimpleExtentDims( &Num );
    if( Num > std::numeric_limits<int>::max() )
      throw std::runtime_error( sCovarRebin + " too many items " + std::to_string( Num ) );
    std::vector<int> Buffer( Num );
    a.read( ::H5::PredType::NATIVE_INT, &Buffer[0] );
    a.close();
    CovarRebin = std::move( Buffer );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute(sCovarSampleSize);
    a.read( ::H5::PredType::NATIVE_INT, &CovarSampleSize );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    CovarSampleSize = 0;
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute(sCovarNumBoot);
    a.read( ::H5::PredType::NATIVE_INT, &CovarNumBoot );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    CovarNumBoot = 0;
    ::H5::Exception::clearErrorStack();
  }
}

template <typename T>
void Model<T>::ValidateAttributes()
{
  Base::ValidateAttributes();
  const int NumOps{ static_cast<int>( OpNames.size() ) };
  const int NumExpected{ NumExponents * ( NumOps + 1 ) + 1 };
  if( Base::Nt_ != NumExpected )
  {
    std::ostringstream s;
    s << "Have " << Base::Nt_ << " parameters, but expected " << NumExpected << ", i.e. "
      << NumExponents << " exponents * ( " << NumOps << " operators + 1 energy ) + chi squared per degree of freedom";
    throw std::runtime_error( s.str().c_str() );
  }
}

template <typename T>
int Model<T>::WriteAttributes( ::H5::Group &g )
{
  int iReturn = Base::WriteAttributes( g ) + 8;
  const hsize_t OneDimension{ 1 };
  ::H5::DataSpace ds1( 1, &OneDimension );
  ::H5::Attribute a;
  a = g.createAttribute( sNumExponents, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &NumExponents );
  a.close();
  a = g.createAttribute( sNumFiles, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &NumFiles );
  a.close();
  a = g.createAttribute( sTI, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &ti );
  a.close();
  a = g.createAttribute( sTF, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &tf );
  a.close();
  a = g.createAttribute( sDoF, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &dof );
  a.close();
  std::int8_t i8{ static_cast<std::int8_t>( Factorised ? 1 : 0 ) };
  a = g.createAttribute( sFactorised, ::H5::PredType::STD_U8LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  i8 = static_cast<std::int8_t>( CovarFrozen ? 1 : 0 );
  a = g.createAttribute( sCovarFrozen, ::H5::PredType::STD_U8LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  H5::WriteMatrix( g, sCovarianceIn+s_C, CovarIn );
  H5::WriteMatrix( g, sCovariance+s_C, Covar );
  H5::WriteMatrix( g, sCorrelation+s_C, Correl );
  H5::WriteMatrix( g, sCorrelationCholesky+s_C, CorrelCholesky );
  H5::WriteMatrix( g, sCovarianceInv+s_C, CovarInv );
  H5::WriteMatrix( g, sCorrelationInvCholesky+s_C, CorrelInvCholesky );
  H5::WriteMatrix( g, sCovarianceInvCholesky+s_C, CovarInvCholesky );
  StdErrorMean.Write( g, sStdErrorMean );
  FitInput.Write( g, sFitInput );
  ModelPrediction.Write( g, sModelPrediction );
  ErrorScaled.Write( g, sErrorScaled );
  {
    std::ostringstream os;
    os << CovarSource;
    a = g.createAttribute( sCovarSource, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, os.str() );
    a.close();
    iReturn++;
  }
  if( !CovarRebin.empty() )
  {
    hsize_t Dims[1] = { CovarRebin.size() };
    ::H5::DataSpace dsp( 1, Dims );
    a = g.createAttribute( sCovarRebin, ::H5::PredType::NATIVE_INT, dsp );
    a.write( ::H5::PredType::NATIVE_INT, &CovarRebin[0] );
    a.close();
    dsp.close();
    iReturn++;
  }
  if( CovarSampleSize )
  {
    a = g.createAttribute( sCovarSampleSize, ::H5::PredType::STD_U32LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &CovarSampleSize );
    a.close();
    iReturn++;
  }
  if( CovarNumBoot )
  {
    a = g.createAttribute( sCovarNumBoot, ::H5::PredType::STD_U32LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &CovarNumBoot );
    a.close();
    iReturn++;
  }
  H5::WriteAttribute( g, sOperators, OpNames );
  return iReturn;
}

template <typename T>
void Model<T>::SummaryComments( std::ostream & s, bool bVerboseSummary ) const
{
  Base::SummaryComments( s, bVerboseSummary );
  s << "# NumExponents: " << NumExponents << NewLine
    << "# NumFiles: " << NumFiles << NewLine
    << "# Factorised: " << Factorised << NewLine
    << "# Frozen covariance: " << CovarFrozen << NewLine
    << "# Covariance source: " << CovarSource << NewLine;
  if( !CovarRebin.empty() )
  {
    s << "# Covariance rebin:";
    for( auto i : CovarRebin )
      s << Space << i;
    s << NewLine;
  }
}

template class Model<double>;
template class Model<float>;
template class Model<std::complex<double>>;
template class Model<std::complex<float>>;

// Make non-random numbers 0 ... m.size1 - 1 to simplify central replica code
// If the sample is already bootstrapped, then we won't need these numbers
template <typename T>
void DataSet<T>::MakeCentralReplicaNonRandom()
{
  for( int i = 0; i < 2; ++i )
  {
    std::vector<fint> &c{ RandomCentralBuf[i] };
    const SS ss{ i == 0 ? SS::Raw : SS::Binned };
    const int Num{ corr.empty() || ( ss == SS::Raw && corr[0].bRawBootstrap ) ? 0 : corr[0].NumSamples( ss ) };
    if( !Num )
      c.clear();
    else
    {
      c.resize( Num );
      for( fint i = 0; i < c.size(); ++i )
        c[i] = i;
    }
  }
}

template <typename T>
void DataSet<T>::clear()
{
  NSamples = 0;
  Extent = 0;
  MinExponents = 0;
  MaxExponents = 0;
  corr.clear();
  FitTimes.clear();
  constFile.clear();
  ConstantNames.clear();
  ConstantNamesPerExp.clear();
  constMap.clear();
}

// Specify which times I'm fitting to, as a list of timeslices for each correlator
template <typename T>
void DataSet<T>::SetFitTimes( const std::vector<std::vector<int>> &fitTimes_ )
{
  if( fitTimes_.size() != corr.size() )
    throw std::runtime_error( std::to_string( fitTimes_.size() ) + " FitTimes but "
                             + std::to_string( corr.size() ) + " correlators" );
  std::vector<std::vector<int>> ft{ fitTimes_ };
  std::size_t extent_ = 0;
  for( int i = 0; i < ft.size(); ++i )
  {
    if( !ft[i].empty() )
    {
      std::sort( ft[i].begin(), ft[i].end() );
      if( ft[i][0] < 0 || ft[i].back() >= corr[i].Nt()
         || std::adjacent_find( ft[i].begin(), ft[i].end() ) != ft[i].end() )
        throw std::runtime_error( "FitTimes[" + std::to_string( i ) + "]=[" + std::to_string( ft[i][0] )
                                 + "..." + std::to_string( ft[i].back() ) + "] invalid" );
      extent_ += ft[i].size();
    }
  }
  SetValidatedFitTimes( extent_, std::move( ft ) );
}

template <typename T>
void DataSet<T>::SetFitTimes( int tMin, int tMax )
{
  const int extent_{ tMax - tMin + 1 };
  std::vector<std::vector<int>> ft{ corr.size() };
  for( int i = 0; i < corr.size(); ++i )
  {
    if( tMin < 0 || extent_ < 1 || tMax >= corr[i].Nt() )
      throw std::runtime_error("Fit range ["+std::to_string(tMin)+", "+std::to_string(tMax)+"] invalid");
    ft[i].reserve( extent_ );
    for( int j = tMin; j <= tMax; ++j )
      ft[i].push_back( j );
  }
  SetValidatedFitTimes( corr.size() * extent_, std::move( ft ) );
}

// Cache the data for this data set
template <typename T>
void DataSet<T>::SetValidatedFitTimes( std::size_t Extent_, std::vector<std::vector<int>> &&FitTimes_ )
{
  if( Extent_ == 0 )
    throw std::runtime_error( "Fit range empty" );
  if( Extent_ > std::numeric_limits<int>::max() )
    throw std::runtime_error( "Fit range stupidly big" );
  Extent = static_cast<int>( Extent_ );
  FitTimes = std::move( FitTimes_ );
  // Get the central replica
  vCentral.resize( Extent );
  {
    int dst{ 0 };
    for( int f = 0; f < corr.size(); ++f )
    {
      const T * pSrc{ corr[f][Fold<T>::idxCentral] };
      for( int t : FitTimes[f] )
        vCentral[dst++] = pSrc[t];
    }
  }
  // Cache the raw data
  for( int i = 0; i < 3; ++i )
  {
    SampleSource ss{ i == 0 ? SampleSource::Raw : i == 1 ? SampleSource::Binned : SampleSource::Bootstrap };
    Matrix<T> &m{ i == 0 ? mRaw : i == 1 ? mBinned : mBoot };
    const int iCount{ corr[0].NumSamples( ss ) };
    if( !iCount )
      m.clear();
    else
    {
      m.resize( iCount, Extent );
      for( int idx = 0; idx < iCount; ++idx )
      {
        int dst{ 0 };
        for( int f = 0; f < corr.size(); ++f )
        {
          const T * pSrc{ corr[f].get( ss, idx ) };
          for( int t : FitTimes[f] )
            m( idx, dst++ ) = pSrc[t];
        }
      }
    }
  }
}

// Get the constants from the appropriate timeslice
template <typename T>
void DataSet<T>::GetFixed( int idx, Vector<T> &vResult, const std::vector<FixedParam> &Params ) const
{
  std::vector<const T *> Src( constFile.size() );
  for( int f = 0; f < constFile.size(); ++f )
    Src[f] = constFile[f][idx];
  for( const FixedParam &p : Params )
    vResult[p.idx] = Src[p.src.File][p.src.idx];
}

// Initialise random numbers I will need if I'm to get (co)variance on non-central replicas
template <typename T>
bool DataSet<T>::InitRandomNumbers( SS ss )
{
  bool bMade{ false };
  if( ss == SS::Binned || ss == SS::Raw )
  {
    const int i{ ss == SS::Raw ? 0 : 1 };
    MatrixView<fint> &r{ RandomViews[i] };
    std::vector<fint> &Buffer{ RandomBuffer[i] };
    const int Count{corr.empty() || ( ss == SS::Raw && corr[0].bRawBootstrap ) ? 0 : corr[0].NumSamples( ss )};
    if( !Count )
    {
      r.clear();
      Buffer.clear();
    }
    else
    {
      // Make Random numbers
      if(   corr[0].RandNum() // Do I have original bootstrap random numbers
         && corr[0].SampleSize == Count ) // Are the random numbers over same range (0...SampleSize-1)
      {
        // Re-use existing random numbers
        // I assume enough replicas are available because this is checked on creation and load
        Buffer.clear();
        r.Map( corr[0].RandNum(), MaxSamples, Count );
      }
      else
      {
        // Generate random numbers
        GenerateRandom( Buffer, corr[0].Name_.Seed, MaxSamples, Count );
        r.Map( Buffer.data(), MaxSamples, Count );
        bMade = true;
      }
    }
  }
  return bMade;
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// From what is already a bootstrap replica, and therefore only the central replica
template <typename T>
void DataSet<T>::MakeCovarFromBootstrap( SS ss, Matrix<T> &Covar ) const
{
  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, replica );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - vCentral[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += z[i] * z[j];
  }
  const T Norm{ static_cast<T>( 1 ) / static_cast<T>( m.size1 ) };
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * Norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

template <typename T>
void DataSet<T>::MakeVarFromBootstrap( SS ss, Vector<T> &Var ) const
{
  Var.resize( Extent );
  for( int i = 0; i < Extent; ++i )
    Var[i] = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, replica );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - vCentral[i];
    for( int i = 0; i < Extent; ++i )
      Var[i] += Squared( z[i] );
  }
  const T Norm{ static_cast<T>( 1 ) / static_cast<T>( m.size1 ) };
  for( int i = 0; i < Extent; ++i )
    Var[i] *= Norm;
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// From the underlying raw/(re)binned data (without bootstrapping)
template <typename T>
void DataSet<T>::MakeCovarFromNonBootstrap( int idx, SS ss, Matrix<T> &Covar ) const
{
  VectorView<T> Mean;
  VectorView<fint> Replace;
  if( idx == Fold<T>::idxCentral )
  {
    // Central replica - no replacement
    Mean = vCentral;
    Replace = RandomCentral( ss );
  }
  else
  {
    // bootstrap replica
    Mean.MapRow( mBoot, idx );
    Replace.MapRow( RandomNumbers( ss ), idx );
  }

  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, Replace[replica] );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - Mean[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += z[i] * z[j];
  }
  const T Norm{ static_cast<T>( 1 ) / ( static_cast<T>( m.size1 - 1 ) * static_cast<T>( m.size1 ) ) };
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * Norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

template <typename T>
void DataSet<T>::MakeVarFromNonBootstrap( int idx, SS ss, Vector<T> &Var ) const
{
  VectorView<T> Mean;
  VectorView<fint> Replace;
  if( idx == Fold<T>::idxCentral )
  {
    // Central replica - no replacement
    Mean = vCentral;
    Replace = RandomCentral( ss );
  }
  else
  {
    // bootstrap replica
    Mean.MapRow( mBoot, idx );
    Replace.MapRow( RandomNumbers( ss ), idx );
  }

  Var.resize( Extent );
  for( int i = 0; i < Extent; ++i )
    Var[i] = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const Matrix<T> &m{ Cache( ss ) };
  for( int replica = 0; replica < m.size1; ++replica )
  {
    Data.MapRow( m, Replace[replica] );
    for( int i = 0; i < Extent; ++i )
      z[i] = Data[i] - Mean[i];
    for( int i = 0; i < Extent; ++i )
      Var[i] += Squared( z[i] );
  }
  const T Norm{ static_cast<T>( 1 ) / ( static_cast<T>( m.size1 - 1 ) * static_cast<T>( m.size1 ) ) };
  for( int i = 0; i < Extent; ++i )
    Var[i] *= Norm;
}

// Make a covariance matrix estimate of \Sigma_{\bar{\vb{x}}}, i.e. error of the mean
// Call the appropriate one of the previous two functions, depending on what the data are
template <typename T>
void DataSet<T>::MakeCovariance( int idx, SS ss, Matrix<T> &Covar ) const
{
  if( ss == SampleSource::Bootstrap || ( ss == SampleSource::Raw && corr[0].bRawBootstrap ) )
  {
    if( idx != Fold<T>::idxCentral )
      throw std::runtime_error( "Can only make covariance from bootstrap on central replica" );
    MakeCovarFromBootstrap( ss, Covar );
  }
  else
    MakeCovarFromNonBootstrap( idx, ss, Covar );
}

template <typename T>
void DataSet<T>::MakeVariance( int idx, SS ss, Vector<T> &Var ) const
{
  if( ss == SampleSource::Bootstrap || ( ss == SampleSource::Raw && corr[0].bRawBootstrap ) )
  {
    if( idx != Fold<T>::idxCentral )
      throw std::runtime_error( "Can only make variance from bootstrap on central replica" );
    MakeVarFromBootstrap( ss, Var );
  }
  else
    MakeVarFromNonBootstrap( idx, ss, Var );
}

// Make covariance matrix using a secondary bootstrap - Random numbers must be provided by caller
// I.e. unfrozen covariance matrix
template <typename T> void DataSet<T>::
MakeCovariance( int idx, SS ss, Matrix<T> &Covar, const MatrixView<fint> &Random ) const
{
  if( ss == SampleSource::Bootstrap || ( ss == SampleSource::Raw && corr[0].bRawBootstrap ) )
    throw std::runtime_error( "Can't rebootstrap a bootstrap" );
  const Matrix<T> &m{ Cache( ss ) };
  if( m.size1 != Random.size2() )
  {
    std::ostringstream os;
    os << ss << " data have " << m.size1 << " samples, but random numbers expect "
       << Random.size2() << " samples";
    throw std::runtime_error( os.str().c_str() );
  }
  VectorView<T> Mean;
  VectorView<fint> Replace;
  if( idx == Fold<T>::idxCentral )
  {
    // Central replica - no replacement
    Mean = vCentral;
    Replace = RandomCentral( ss );
  }
  else
  {
    // bootstrap replica
    Mean.MapRow( mBoot, idx );
    Replace.MapRow( RandomNumbers( ss ), idx );
  }

  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  VectorView<T> Data;
  Vector<T> z( Extent );
  const T InvBootLen{ static_cast<T>( 1 ) / static_cast<T>( m.size1 ) };
  for( int replica = 0; replica < Random.size1(); ++replica )
  {
    // Perform the covariance (inner) bootstrap
    for( int i = 0; i < Extent; ++i )
      z[i] = 0;
    for( int inner = 0; inner < m.size1; ++inner )
    {
      Data.MapRow( m, Replace[ Random( replica, inner ) ] );
      for( int i = 0; i < Extent; ++i )
        z[i] += 0;
    }
    // Now contribute this to the covariance
    for( int i = 0; i < Extent; ++i )
      z[i] = z[i] * InvBootLen - Mean[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += z[i] * z[j];
  }
  const T Norm{ static_cast<T>( 1 ) / static_cast<T>( Random.size1() ) };
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * Norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

// Write covariance matrix to file
template <typename T>
void DataSet<T>::SaveMatrixFile( const Matrix<T> &m, const std::string &Type, const std::string &Filename,
                                 std::vector<std::string> &Abbreviations,
                                 const std::vector<std::string> *FileComments, const char *pGnuplotExtra ) const
{
  // For now, all the matrices I write are square and match dataset, but this restriction can be safely relaxed
  if( m.size1 != Extent || m.size1 != m.size2 )
    throw std::runtime_error( Type + " matrix dimensions (" + std::to_string( m.size1 ) + Common::CommaSpace
                             + std::to_string( m.size2 ) + ") don't match dataset" );
  // Header describing what the covariance matrix is for and how to plot it with gnuplot
  std::ofstream s{ Filename };
  s << "# Matrix: " << Filename << "\n# Type: " << Type << "\n# Files: " << corr.size() << Common::NewLine;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    s << "# File" << f << ": " << corr[f].Name_.NameNoExt << "\n# Times" << f << ":";
    for( int t : FitTimes[f] )
      s << Common::Space << t;
    s << Common::NewLine;
    // Say which operators are in each file
    if( FileComments && !(*FileComments)[f].empty() )
      s << (*FileComments)[f]; // These are expected to have a trailing NewLine
    // Make default abbreviation - a letter for each file
    if( Abbreviations.size() < f )
      Abbreviations.emplace_back( 1, 'A' + f );
  }
  // Save a command which can be used to plot this file
  s << "# gnuplot: set xtics rotate noenhanced; set ytics noenhanced; set title noenhanced '" << Filename
    << "'; set key font 'Arial,8' top right noenhanced; ";
  if( pGnuplotExtra && *pGnuplotExtra )
    s << pGnuplotExtra << "; ";
  s << "plot '" << Filename << "' matrix columnheaders rowheaders with image pixels\n";
  // Now save the column names. First column is matrix extent
  s << "# The next row contains column headers, starting with matrix extent\n" << Extent;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    for( int t : FitTimes[f] )
      s << Common::Space << Abbreviations[f] << t;
  }
  s << Common::NewLine << std::setprecision( std::numeric_limits<T>::max_digits10 );
  // Now print the actual covariance matrix
  int i{ 0 };
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    for( int t : FitTimes[f] )
    {
      s << Abbreviations[f] << t;
      for( int j = 0; j < Extent; ++j )
        s << Common::Space << m( i, j );
      s << Common::NewLine;
      ++i;
    }
  }
}

// Add a constant to my list of known constants - make sure it isn't already there
template <typename T>
void DataSet<T>::AddConstant( const std::string &Name, std::size_t File, std::size_t idx )
{
  if( constMap.find( Name ) != constMap.end() )
    throw std::runtime_error( "Constant \"" + Name + "\" loaded from multiple model files" );
  constMap.insert( { Name, ConstantSource( File, idx ) } );
}

template <typename T>
void DataSet<T>::AddConstant( const std::string &Name, std::size_t File, std::size_t idx, int e )
{
  std::string s{ Name };
  s.append( std::to_string( e ) );
  AddConstant( s, File, idx );
}

template <typename T>
int DataSet<T>::LoadCorrelator( Common::FileNameAtt &&FileAtt, unsigned int CompareFlags,
                                const char * PrintPrefix )
{
  if( corr.size() >= std::numeric_limits<int>::max() )
    throw std::runtime_error( "More than an integer's worth of correlators loaded!!??" );
  const int i{ static_cast<int>( corr.size() ) };
  corr.emplace_back();
  corr[i].SetName( std::move( FileAtt ) );
  corr[i].Read( PrintPrefix );
  // See whether this correlator is compatible with prior correlators
  if( i )
  {
    corr[0].IsCompatible( corr[i], &NSamples, CompareFlags );
    if( MaxSamples > corr[i].NumSamples() )
      MaxSamples = corr[i].NumSamples();
  }
  else
  {
    MaxSamples = corr[i].NumSamples();
    if( NSamples == 0 || NSamples > MaxSamples )
      NSamples = MaxSamples;
    OriginalBinSize = corr[i].binSize; // Because this will be overwritten
    MakeCentralReplicaNonRandom();
  }
  return i;
}

// Load a model file
template <typename T>
void DataSet<T>::LoadModel( Common::FileNameAtt &&FileAtt, const std::string &Args )
{
  static const std::string EnergyPrefix{ "E" };
  // Break trailing arguments for this file into an array of strings
  std::vector<std::string> vThisArg{ Common::ArrayFromString( Args ) };
  // This is a pre-built model (i.e. the result of a previous fit)
  const std::size_t i{ constFile.size() };
  constFile.emplace_back();
  constFile[i].SetName( std::move( FileAtt ) );
  constFile[i].Read( "  " );
  // Keep track of minimum number of replicas across all files
  if( NSamples == 0 )
    NSamples = constFile[i].NumSamples();
  else if( NSamples > constFile[i].NumSamples() )
    NSamples = constFile[i].NumSamples();
  // Keep track of minimum and maximum number of exponents across all files
  if( i )
  {
    if( MinExponents > constFile[i].NumExponents )
      MinExponents = constFile[i].NumExponents;
    if( MaxExponents < constFile[i].NumExponents )
      MinExponents = constFile[i].NumExponents;
  }
  else
  {
    MinExponents = constFile[i].NumExponents;
    MaxExponents = constFile[i].NumExponents;
  }
  // Now see which parameters we want to read from the model
  const std::vector<std::string> &ColumnNames{ constFile[i].GetColumnNames() };
  const std::vector<std::string> &OpNames{ constFile[i].OpNames };
  if( vThisArg.empty() )
  {
    // Load every constant in this file
    for( int j = 0; j < ColumnNames.size(); ++j )
    {
      // Add a map back to this specific parameter (making sure not already present)
      AddConstant( ColumnNames[j], i, j );
      // Now take best stab at whether this should be per exponent or individual
      std::string s{ ColumnNames[j] };
      int Exp;
      if( Common::ExtractTrailing( s, Exp )
         && ( Common::EqualIgnoreCase( EnergyPrefix, s ) || Common::IndexIgnoreCase( OpNames, s ) != OpNames.size() ) )
      {
        ConstantNamesPerExp[s];
      }
      else
        ConstantNames[ColumnNames[j]];
    }
  }
  else
  {
    // Load only those constants specifically asked for
    for( const std::string &ThisArg : vThisArg )
    {
      std::vector<std::string> vPar{ Common::ArrayFromString( ThisArg, "=" ) };
      if( vPar.empty() || vPar[0].empty() || vPar.size() > 2 || ( vPar.size() > 1 && vPar[1].empty() ) )
        throw std::runtime_error( "Cannot interpret model parameter string \"" + Args + "\"" );
      const std::string &vLookFor{ vPar.back() };
      const std::string &vLoadAs{ vPar[0] };
      // Have we asked for a per-exponent constant (which includes energies)
      if( Common::EqualIgnoreCase( EnergyPrefix, vLookFor ) || Common::IndexIgnoreCase( OpNames, vLookFor )!=OpNames.size() )
      {
        for( int e = 0; e < constFile[i].NumExponents; ++e )
          AddConstant( vLoadAs, i, constFile[i].GetColumnIndex( vLookFor, e ), e );
        ConstantNamesPerExp[vLoadAs];
      }
      else
      {
        AddConstant( vLoadAs, i, constFile[i].GetColumnIndex( vLookFor ) );
        ConstantNames[vLoadAs];
      }
    }
  }
}

// Sort the operator names and renumber all the loaded correlators referring to them
template <typename T>
void DataSet<T>::SortOpNames( std::vector<std::string> &OpNames )
{
  int NumOps{ static_cast<int>( OpNames.size() ) };
  if( OpNames.size() > 1 )
  {
    // Sort the names
    UniqueNames OpSorted;
    for( int i = 0; i < NumOps; ++i )
      OpSorted.emplace( std::move( OpNames[i] ), i );
    // Extract the sorted names and indices (to renumber operators in correlator names)
    std::vector<std::string> SortedNames;
    std::vector<int> SortIndex( NumOps );
    SortedNames.reserve( NumOps );
    int idx{ 0 };
    for( UniqueNames::iterator it = OpSorted.begin(); it != OpSorted.end(); ++it )
    {
      SortedNames.emplace_back( it->first );
      SortIndex[it->second] = idx++;
    }
    // Renumber the operators and save the sorted operator names
    for( auto &f : corr )
      for( int i = 0; i < f.Name_.op.size(); ++i )
        f.Name_.op[i] = SortIndex[f.Name_.op[i]];
    OpNames = SortedNames;
  }
}

template <typename T>
void DataSet<T>::Rebin( const std::vector<int> &NewSize )
{
  RebinSize.clear();
  RebinSize.reserve( corr.size() );
  for( std::size_t i = 0; i < corr.size(); ++i )
  {
    if( !corr[i].NumSamplesRaw() )
    {
      std::ostringstream os;
      os << "Raw samples unavailable rebinning corr " << i << CommaSpace << corr[i].Name_.Filename;
      throw std::runtime_error( os.str().c_str() );
    }
    // Rebin to the size specified (last bin size repeats to end).
    // bin size 0 => auto (i.e. all measurements on same config binned together)
    RebinSize.push_back( NewSize.empty() ? 0 : NewSize[i < NewSize.size() ? i : NewSize.size() - 1] );
    if( RebinSize.back() )
      corr[i].Bin( RebinSize.back(), SS::Raw );
    else
      corr[i].Bin( SS::Raw );
    if( corr[i].NumSamplesRaw() != corr[0].NumSamplesRaw() )
    {
      std::ostringstream os;
      os << "Rebinned corr " << i << " has " << corr[i].NumSamplesRaw()
         << " samples, others have " << corr[i].NumSamplesRaw();
      throw std::runtime_error( os.str().c_str() );
    }
  }
  MakeCentralReplicaNonRandom();
}

template class DataSet<double>;
template class DataSet<float>;
template class DataSet<std::complex<double>>;
template class DataSet<std::complex<float>>;

// Read an array (real or complex) from an HDF5 file
#define template_ReadArrayArgs( T ) \
std::vector<T> &buffer, const std::string &FileName, const std::string &ObjectName, \
const char *PrintPrefix, std::string * pGroupName

template<typename T> void ReadArray( template_ReadArrayArgs( T ) )
{
  ::H5::H5File f;
  ::H5::Group  g;
  H5::OpenFileGroup( f, g, FileName, PrintPrefix, pGroupName );
  ::H5::DataSet ds = g.openDataSet(ObjectName);
  ::H5::DataSpace dsp = ds.getSpace();
  const int nDims{dsp.getSimpleExtentNdims()};
  if( nDims != 1 )
    throw std::runtime_error("Object " + ObjectName + " in " + FileName + " has " + std::to_string( nDims ) + " dimensions" );
  hsize_t Nt;
  dsp.getSimpleExtentDims( &Nt );
  hsize_t BufferSize{ static_cast<hsize_t>( buffer.size() ) };
  if( BufferSize == 0 )
    buffer.resize( Nt );
  else if( BufferSize != Nt )
    throw std::runtime_error("Object " + ObjectName + " in " + FileName + " has " + std::to_string( Nt ) + " entries, doesn't match Nt=" + std::to_string( BufferSize ) );
  ds.read( &buffer[0], H5::Equiv<T>::Type );
  if( !IsFinite( buffer ) )
     throw std::runtime_error( "Error: Infinite/NaN values in " + FileName );
}

// Make sure all the specialisations I support are instantiated
template void ReadArray<float>( template_ReadArrayArgs( float ) );
template void ReadArray<double>( template_ReadArrayArgs( double ) );
template void ReadArray<long double>( template_ReadArrayArgs( long double ) );
template void ReadArray<std::complex<float>>( template_ReadArrayArgs( std::complex<float> ) );
template void ReadArray<std::complex<double>>( template_ReadArrayArgs( std::complex<double> ) );
template void ReadArray<std::complex<long double>>( template_ReadArrayArgs(std::complex<long double>));

std::ostream& operator<<( std::ostream& os, const CommandLine &cl)
{
  os << "Command-line has " << cl.Args.size() << " arguments and " << cl.Switches.size() << " switches";
  for( int i = 0; i < cl.Args.size(); i++ )
    os << "\nArg[" << i << "]=\"" << cl.Args[i] << "\"";
  int i = 0;
  //for( const CommandLine::SwitchPair &p : cl.Switches ) {
  for( CommandLine::SwitchMap::const_iterator p = cl.Switches.begin(); p != cl.Switches.end(); ++p ) {
    const std::string &Switch{ p->first };
    const std::vector<std::string> &Values{ p->second };
    os << "\nSwitch[" << i << "]=\"" << Switch << "\", was specified";
    const std::size_t sz{ Values.size() };
    if( sz ) {
      os << " " << Values.size() << " time";
      if( sz > 1 )
        os << "s";
    }
    for( int j = 0; j < Values.size(); j++ )
      os << "\n      [" << i << "][" << j << "]=\"" << Values[j] << "\"";
    ++i;
  }
  return os;
}

bool CommandLine::IsValuePresent( const char * & p )
{
  const char * const pOriginal{ p };
  while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
    p++;
  if( *p == '=' ) {
    ++p;
    while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
      p++;
  }
  // There's a value present if there are characters in the string
  // ... or if we skipped past whitespace and/or equal sign to get to end of string
  return *p || pOriginal != p;
}

void CommandLine::Parse( int argc, const char *argv[], const std::vector<SwitchDef> &defs )
{
  // Save the name of the executable without any path
  Name = argv[0];
  auto pos = Name.find_last_of('/');
  if( pos != std::string::npos )
    Name = Name.substr( pos + 1 );
  // Now parse the command-line
  bool bError{ false };
  int SwitchNo = -1; // Not waiting for switch value
  bool bInMultiChar{ false };
  for( int i = 1; i < argc; i++ ) {
    const char *p = argv[i];
    if( *p == '-' ) { // Switch
      if( SwitchNo != -1 ) { // We were waiting on a switch that requires a value
        Switches[defs[SwitchNo].Switch].push_back( "" );
        SwitchNo = -1;
      }
      std::string SwitchName;
      if( *++p == '-' ) {
        // multi-char switch.
        const char * pEnd = ++p;
        if( *p == 0 ) {
          // This is the special switch "--"
          p -= 2;
          bInMultiChar = true;
        } else {
          while( *pEnd && *pEnd != '=' && *pEnd != ' ' && *pEnd != '\t' && *pEnd != '\r' && *pEnd != '\n' )
            pEnd++;
        }
        SwitchName = std::string( p, pEnd - p );
        p = pEnd;
      } else {
        // Single character switch
        if( *p == 0 ) {
          std::cerr << "Unrecognised switch \"-\"" << std::endl;
          bError = true;
          continue;
        }
        SwitchName.resize( 1 );
        SwitchName[0] = *p++;
      }
      // I've just decoded a short or long switch name - see whether it's valid
      SwitchNo = 0;
      while( SwitchNo < defs.size() && SwitchName.compare( defs[SwitchNo].Switch ) )
        SwitchNo++;
      if( SwitchNo >= defs.size() ) {
        std::cerr << "Unrecognised switch \"" << argv[i] << "\"" << std::endl;
        bError = true;
        SwitchNo = -1;
        continue;
      }
      if( defs[SwitchNo].Type == SwitchType::Flag ) {
        // This is just a switch - it should not have a value
        // Swallow any trailing whitespace
        while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
          p++;
        if( *p ) {
          std::cerr << "Switch \"" << SwitchName << "\" should not have a value \"" << p << "\"" << std::endl;
          bError = true;
        } else if( GotSwitch( SwitchName ) ) {
          std::cerr << "Switch \"" << SwitchName << "\" should not be repeated" << std::endl;
          bError = true;
        } else
          Switches[SwitchName]; // First time we've seen this flag
        SwitchNo = -1;
        continue;
      } else if( defs[SwitchNo].Type == SwitchType::Single && GotSwitch( SwitchName ) ) {
        std::cerr << "Switch \"" << SwitchName << "\" should not be repeated" << std::endl;
        bError = true;
        SwitchNo = -1;
        continue;
      }
      // Use the remainder of this switch as a value ... or wait for next param if empty
      if( !IsValuePresent( p ) )
        continue;
    }
    if( SwitchNo != -1 ) {
      // Get the value of the current switch
      Switches[defs[SwitchNo].Switch].push_back( p );
      if( !bInMultiChar )
        SwitchNo = -1;
    }
    else // Argument
      Args.push_back( p );
  }
  if( SwitchNo != -1 && !bInMultiChar ) {
    std::cerr << "Still processing switch \"" << defs[SwitchNo].Switch << "\"" << std::endl;
    bError = true;
  }
  if( bError ) {
    // std::cerr << (*this) << std::endl;
    throw std::runtime_error( "Invalid command-line" );
  }
  // Now put defaults in for any missing switches
  for( int i = 0; i < defs.size(); i++ ) {
    if( defs[i].Type == Single && defs[i].Default && !GotSwitch( defs[i].Switch ) )
      Switches[defs[i].Switch].push_back( defs[i].Default );
  }
}
END_COMMON_NAMESPACE
