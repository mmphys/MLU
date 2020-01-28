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

#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <sys/stat.h>
#include <H5CompType.h>

BEGIN_COMMON_NAMESPACE

// Text required for summaries of correlators
namespace CorrSumm {
  const char sep[] = " ";
  const char Comment[] = "# ";
  const char NewLine[] = "\n";

  const char * SummaryNames[NumSummaries] = { "corr", "mass", "cosh" };

  const char FieldNames[] = "t y y_low y_high y_check";
  const char FieldNames2[] = " im im_low im_high im_check";
  const char * SummaryHeader[NumSummaries] =
  {
    "correlator",
    "mass",
    "cosh mass",
  };
};

extern const std::string Underscore{ "_" };
extern const std::string Period{ "." };
const std::string sBootstrap{ "bootstrap" };
const std::string sModel{ "model" };
const double NaN{ std::nan( "" ) };

namespace Gamma
{
  const std::array<std::string, nGamma> name{ // Long name, per Grid
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
  const std::array<std::string, nGamma> nameShort{ // My abbreviations
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

  std::string NameShort( Algebra alg, const char * pszPrefix )
  {
    std::string sName;
    if( alg != Algebra::Unknown )
    {
      if( pszPrefix )
        sName.append( pszPrefix );
      int idx = static_cast<int>( alg );
      sName.append( ( idx >= 0 && idx < nGamma ) ? nameShort[idx] : "?" );
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
  while( pszDelimeters && pszDelimeters[NumDelims] )
    NumDelims++;
  const std::size_t Len{ String.length() };
  std::size_t Start{ Len };
  bool bFoundDelim{ false };
  char c = 0;
  while( !bFoundDelim && Start )
  {
    c = String[--Start];
    for( std::size_t i = 0; !bFoundDelim && i < NumDelims; i++ )
    {
      bFoundDelim = ( pszDelimeters[i] == c );
      if( bFoundDelim )
      {
        Suffix = String.substr( Start + 1 );
        // Skip past multiple delimeters if they are all whitespace
        if( c == ' ' || c == '\t' || c == '\r' || c == '\n' )
        {
          while( Start && ( String[Start - 1] == ' ' || String[Start - 1] == '\t'
                         || String[Start - 1] == '\r' || String[Start - 1] == '\n' ) )
            --Start;
        }
        String.resize( Start );
      }
    }
  }
  return bFoundDelim;
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

// Does the specified file exist?
bool FileExists( const std::string& Filename )
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

// Sort the list of values, then extract the lower and upper 68th percentile error bars
void ValWithEr::Get( double Central_, std::vector<double> &Data, std::size_t Count )
{
  assert( Data.size() >= Count && "ValWithErr<T>::Get: data too small" );
  Central = Central_;
  if( Count == 0 ) {
    ErLow = NaN;
    ErHigh  = NaN;
    Check = 0;
    return;
  }
  const typename std::vector<double>::iterator itStart{ Data.begin() };
  const typename std::vector<double>::iterator itEnd  { itStart + Count };
  std::sort( itStart, itEnd );
  std::size_t Index = 0.16 * Count + 0.5;
  if( Index >= Count )
    Index = Count - 1;
  ErLow  = Central - Data[Index];
  ErHigh = Data[Count - 1 - Index] - Central;
  Check = static_cast<double>( Count ) / Data.size();
}

std::string Momentum::p2_string( const std::string &separator ) const
{
  std::string s{ separator };
  s.append( "p2" );
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

// Strip out momentum info from string if present
bool Momentum::Extract( std::string &Prefix, bool IgnoreSubsequentZeroNeg )
{
  bool bGotMomentum = false;
  std::smatch match;
  static const std::regex pattern{ R"(_[pP]_(-?[0-9]+)_(-?[0-9]+)_(-?[0-9]+))" };
  while( std::regex_search( Prefix, match, pattern  ) )
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
      if( IgnoreSubsequentZeroNeg && !(*this) )
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
    Prefix = match.prefix();
    Prefix.append( match.suffix() );
  }
  return bGotMomentum;
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
void FileNameAtt::Parse( const std::string &Filename_ )
{
  clear();
  Filename = Filename_;
  std::size_t pos = Filename.find_last_of( '/' );
  if( pos == std::string::npos )
    pos = 0;
  else
    Dir = Filename.substr( 0, ++pos );
  Base = Filename.substr( pos );
  int i = 0;
  for( ; i < 3 && ( pos = Base.find_last_of( '.' ) ) != std::string::npos ; i++ ) {
    switch( i ) {
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
  }
  /*if( i < 3 ) {
    std::cout << "Warning: Missing type ";
    if( i < 2 ) {
      std::cout << "+ Seq ";
      if( i < 1 ) {
        std::cout << "+ extension ";
      }
    }
    std::cout << "in " << Filename << std::endl;
  }*/
}

// The base should end with an operator in my list
FileNameAtt::FileNameAtt( const std::string &Filename, std::vector<std::string> &OpNames )
{
  Parse( Filename );
  char Sep = '_';
  std::size_t pos = Base.find_last_of( Sep );
  if( pos != std::string::npos ) {
    std::string sOp{ Base.substr( pos + 1 ) }; // Operator name
    int iOp1 = 0;
    while( iOp1 < OpNames.size() && !EqualIgnoreCase( sOp, OpNames[iOp1] ) )
      iOp1++;
    if( iOp1 == OpNames.size() )
      OpNames.emplace_back( sOp );
    std::string sTmp{ Base.substr( 0, pos ) };  // Truncated string
    pos = sTmp.find_last_of( Sep );
    if( pos != std::string::npos ) {
      sOp = sTmp.substr( pos + 1 ); // Operator name
      int iOp2 = 0;
      while( iOp2 < OpNames.size() && !EqualIgnoreCase( sOp, OpNames[iOp2] ) )
        iOp2++;
      if( iOp2 == OpNames.size() )
        OpNames.emplace_back( sOp );
      // Got valid sink and source operators
      op.resize( 2 );
      op[0] = iOp1;
      op[1] = iOp2;
      Base.resize( pos );
      return;
    }
  }
  throw std::runtime_error( "Invalid operator names at end of " + Base );
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

// Strip out timeslice info from a string if present
void ExtractTimeslice( std::string &Prefix, bool &bHasTimeslice, int &Timeslice )
{
  std::smatch match;
  static const std::regex pattern{ R"(_[tT]_?([0-9]+))" };
  while( std::regex_search( Prefix, match, pattern  ) )
  {
    int ThisTimeslice = std::stoi( match[1] );
    if( !bHasTimeslice )
    {
      Timeslice = std::stoi( match[1] );
      bHasTimeslice = true;
    }
    else if( ThisTimeslice != Timeslice )
      throw std::runtime_error( "Multiple momenta: " + std::to_string(Timeslice) + " and " + std::to_string(ThisTimeslice) );
    Prefix = match.prefix();
    Prefix.append( match.suffix() );
  }
}

const H5::PredType& H5EquivType<float>      ::predType{ H5::PredType::NATIVE_FLOAT };
const H5::PredType& H5EquivType<double>     ::predType{ H5::PredType::NATIVE_DOUBLE };
const H5::PredType& H5EquivType<long double>::predType{ H5::PredType::NATIVE_LDOUBLE };

// Return the same HDF5 complex type Grid uses
H5::CompType & H5File::ComplexType( int fpsize )
{
  static H5::CompType m_Complexf(sizeof(std::complex<float>));
  static H5::CompType m_Complexd(sizeof(std::complex<double>));
  static H5::CompType m_Complexl(sizeof(std::complex<long double>));
  {
#ifndef __INTEL_COMPILER
    // mutex not available for this compiler, so initialisation not thread-safe
    static std::mutex ComplexTypeSync;
    std::lock_guard<std::mutex> guard( ComplexTypeSync );
#endif
    static bool bInitialised = false;
    if( !bInitialised ) {
      m_Complexf.insertMember("re", 0 * sizeof(float), H5EquivType<float>::predType);
      m_Complexf.insertMember("im", 1 * sizeof(float), H5EquivType<float>::predType);
      m_Complexd.insertMember("re", 0 * sizeof(double), H5EquivType<double>::predType);
      m_Complexd.insertMember("im", 1 * sizeof(double), H5EquivType<double>::predType);
      m_Complexl.insertMember("re", 0 * sizeof(long double), H5EquivType<long double>::predType);
      m_Complexl.insertMember("im", 1 * sizeof(long double), H5EquivType<long double>::predType);
      bInitialised = true;
    }
  }
  if( fpsize == sizeof( float ) )
    return m_Complexf;
  if( fpsize == sizeof( double ) )
    return m_Complexd;
  assert( fpsize == sizeof( long double ) );
  return m_Complexl;
}

// Get first groupname from specified group
std::string GetFirstGroupName( H5::Group & g )
{
  hsize_t n = g.getNumObjs();
  for( hsize_t i = 0; i < n; ++i ) {
    H5G_obj_t t = g.getObjTypeByIdx( i );
    if( t == H5G_GROUP )
      return g.getObjnameByIdx( i );
  }
  return std::string();
}

// Open the specified HDF5File and group
void OpenHdf5FileGroup(H5::H5File &f, H5::Group &g, const std::string &FileName, std::string &GroupName, const char *PrintPrefix, unsigned int flags )
{
  const bool bFindGroupName{ GroupName.empty() };
  f.openFile( FileName, flags );
  g = f.openGroup( bFindGroupName ? std::string("/") : GroupName );
  if( bFindGroupName ) {
    GroupName = GetFirstGroupName( g );
    g = g.openGroup( GroupName );
  }
  if( PrintPrefix )
    std::cout << PrintPrefix << FileName << " (" << GroupName << ")\n";
}

// Read the gamma algebra attribute string and make sure it's valid
Gamma::Algebra ReadGammaAttribute( H5::Group &g, const char * pAttName )
{
  std::string sGamma;
  H5::Attribute a = g.openAttribute( pAttName );
  H5::StrType s = a.getStrType();
  a.read( s, sGamma );
  a.close();
  for( int idxGamma = 0; idxGamma < Gamma::nGamma; idxGamma++ )
    if( EqualIgnoreCase( sGamma, Gamma::name[idxGamma] ) )
      return static_cast<Gamma::Algebra>( idxGamma );
  throw H5::Exception( "Common::ReadGammaAttribute", "Invalid gamma algebra string" );
}

static const std::string Signal_TextEquiv_Real{ "real" };
static const std::string Signal_TextEquiv_Imaginary{ "imaginary" };
static const std::string Signal_TextEquiv_Unknown{ "unknown" };

std::ostream& operator<<(std::ostream& os, const Signal &sig)
{
  switch( sig )
  {
    case Signal::Real:
      os << Signal_TextEquiv_Real;
      break;
    case Signal::Imag:
      os << Signal_TextEquiv_Imaginary;
      break;
    default:
      os << Signal_TextEquiv_Unknown;
      break;
  }
  return os;
}

std::istream& operator>>(std::istream& is, Signal &sig)
{
  std::string s;
  sig = Signal::Unknown;
  if( is >> s )
  {
    if( EqualIgnoreCase( s, Signal_TextEquiv_Real ) )
      sig = Signal::Real;
    else if( EqualIgnoreCase( s, Signal_TextEquiv_Imaginary ) )
      sig = Signal::Imag;
    else if( !EqualIgnoreCase( s, Signal_TextEquiv_Unknown ) )
      is.setstate( std::ios_base::failbit );
  }
  return is;
}

// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      std::string &GroupName, const std::string &ObjectName )
{
  H5::H5File f;
  H5::Group  g;
  OpenHdf5FileGroup( f, g, FileName, GroupName );
  H5::DataSet ds = g.openDataSet(ObjectName);
  H5::DataSpace dsp = ds.getSpace();
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
  ds.read( &buffer[0], Common::H5File::ComplexType() );
  for( hsize_t t = 0; t < Nt; ++t )
    if( !std::isfinite( buffer[t].real() ) || !std::isfinite( buffer[t].imag() ) )
       throw std::runtime_error( "Error: Infinite/NaN values in " + FileName );
}

// Read a real array from an HDF5 file
void ReadRealArray(std::vector<double> &buffer, const std::string &FileName,
                      std::string &GroupName, const std::string &ObjectName )
{
  H5::H5File f;
  H5::Group  g;
  OpenHdf5FileGroup( f, g, FileName, GroupName );
  H5::DataSet ds = g.openDataSet(ObjectName);
  H5::DataSpace dsp = ds.getSpace();
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
  ds.read( &buffer[0], H5::PredType::NATIVE_DOUBLE );
  for( hsize_t t = 0; t < Nt; ++t )
    if( !std::isfinite( buffer[t] ) )
       throw std::runtime_error( "Error: Infinite/NaN values in " + FileName );
}

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

void CommandLine::SkipPastSep( const char * & p )
{
  while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
    p++;
  if( *p == '=' ) {
    ++p;
    while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
      p++;
  }
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
      // Swallow any trailing whitespace
      while( *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' )
        p++;
      if( defs[SwitchNo].Type == SwitchType::Flag ) {
        // This is just a switch - it should not have a value
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
      SkipPastSep( p );
      if( *p == 0 )
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

// Make summary files of a bootstrap of a correlator
void SummariseBootstrapCorr(const Common::SampleC &out, const std::string & sOutFileBase, SeedType Seed )//, int momentum_squared)
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  const int nt{ out.Nt() };
  const int nSample{ out.NumSamples() };
  std::vector<double> Data( nSample );
  std::vector<double> DataImag( nSample );
  using C = std::complex<double>;
  using type_ReIm = double (*)( const C & );
  type_ReIm myReal, myImag;
  if( out.Signal_ == Signal::Imag )
  {
    myReal = &std::imag;
    myImag = &std::real;
  }
  else
  {
    myReal = &std::real;
    myImag = &std::imag;
  }
  for(int f = 0; f < NumSummaries; f++)
  {
    std::string sOutFileName{ MakeFilename(sOutFileBase, SummaryNames[f], Seed, TEXT_EXT) };
    std::ofstream s( sOutFileName );
    if( !s )
      throw std::runtime_error( "Unable to create " + sOutFileName );
    s << Comment << SummaryHeader[f] << NewLine << Comment << sOutFileBase << "\n# Seed " << Seed
      << "\n# Signal " << out.Signal_ << NewLine << FieldNames << ( ( f == 0 ) ? FieldNames2 : "" )
      << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
    for(int t = 0; t < nt; t++)
    {
      std::size_t Count = 0;
      if( f == 0 )
      {
        const std::complex<double> *pc = out[0] + t;
        for( int i = 0; i < nSample; i++, pc += nt )
        {
          double re = myReal( * pc );
          double im = myImag( * pc );
          if( std::isfinite( re ) && std::isfinite( im ) )
          {
            Data[Count] = re;
            DataImag[Count] = im;
            Count++;
          }
        }
        pc = out[Common::SampleC::idxCentral] + t;
        Common::ValWithEr Re( myReal( * pc ), Data, Count );
        Common::ValWithEr Im( myImag( * pc ), DataImag, Count );
        s << t << sep << Re.Central << sep << Re.ErLow << sep << Re.ErHigh
               << sep << ( static_cast<double>( Count ) / nSample )
               << sep << Im.Central << sep << Im.ErLow << sep << Im.ErHigh
               << sep << ( static_cast<double>( Count ) / nSample ) << std::endl;
      }
      else
      {
        double dCentral = 0;
        for( int i = Common::SampleC::idxCentral; i < nSample; i++ )
        {
          double DThis;
          switch(f)
          {
            case 1: // mass
              DThis = std::log( myReal( out[i][t] ) / myReal( out[i][(t + 1 + nt) % nt] ) );
              break;
            case 2: // cosh mass
              DThis = std::acosh((myReal( out[i][(t - 1 + nt) % nt] ) + myReal( out[i][(t + 1) % nt] ) ) / (2 * myReal( out[i][t] )));
              break;
            default:
              DThis = 0;
          }
          if( i == Common::SampleC::idxCentral )
            dCentral = DThis;
          else if( std::isfinite( DThis ) )
            Data[Count++] = DThis;
        }
        Common::ValWithEr Re( dCentral, Data, Count );
        s << t << sep << Re.Central << sep << Re.ErLow << sep << Re.ErHigh << sep
          << ( static_cast<double>( Count ) / nSample ) << std::endl;
      }
    }
  }
}
END_COMMON_NAMESPACE
