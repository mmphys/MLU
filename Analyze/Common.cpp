/*************************************************************************************
 
 Common utilities (no dependencies other than c++ stdlib)
 
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
 *************************************************************************************/
/*  END LEGAL */

#include "Common.hpp"

#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <sys/stat.h>
#include <H5CompType.h>

BEGIN_COMMON_NAMESPACE

const std::string sBootstrap{ "bootstrap" };
const std::string sModel{ "model" };
const double NaN{ std::nan( "" ) };

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

// These are the attributes I like to use in my filenames
FileNameAtt::FileNameAtt( const std::string &Filename ) : Seed( 0 )
{
  std::size_t pos = Filename.find_last_of( '/' );
  if( pos == std::string::npos )
    pos = 0;
  else
    pos++;
  Base = Filename.substr( pos );
  int i = 0;
  for( ; i < 3 && ( pos = Base.find_last_of( '.' ) ) != std::string::npos ; i++ ) {
    switch( i ) {
      case 0:
        Ext = Base.substr( pos + 1 );
        break;
      case 1:
        try {
          Seed = FromString<unsigned int>( Base.substr( pos + 1 ) );
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
  if( i < 3 ) {
    std::cout << "Warning: Missing type ";
    if( i < 2 ) {
      std::cout << "+ Seq ";
      if( i < 1 ) {
        std::cout << "+ extension ";
      }
    }
    std::cout << "in " << Filename << std::endl;
  }
}

// The base should end with an operator in my list
FileNameAtt::FileNameAtt( const std::string &Filename, std::vector<std::string> &OpNames )
  : FileNameAtt( Filename )
{
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

// Return the same HDF5 complex type Grid uses
#ifndef __INTEL_COMPILER
// mutex not available for this compiler, so initialisation not thread-safe
static std::mutex ComplexTypeSync;
#endif

H5::CompType & H5File::ComplexType(void)
{
  static bool bInitialised = false;
  static H5::CompType m_Complex(sizeof(std::complex<double>));
  {
#ifndef __INTEL_COMPILER
    std::lock_guard<std::mutex> guard( ComplexTypeSync );
#endif
    if( !bInitialised ) {
      m_Complex.insertMember("re", 0 * sizeof(double), H5::PredType::NATIVE_DOUBLE);
      m_Complex.insertMember("im", 1 * sizeof(double), H5::PredType::NATIVE_DOUBLE);
      bInitialised = true;
    }
  }
  return m_Complex;
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

// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      const std::string &GroupName, const std::string &ObjectName )
{
  std::cout << "Reading complex array from " << FileName;
  H5::H5File f(FileName, H5F_ACC_RDONLY);
  const bool bFindGroupName{ GroupName.empty() };
  H5::Group g = f.openGroup( bFindGroupName ? std::string("/") : GroupName );
  std::cout << ", group ";
  if( bFindGroupName ) {
    std::string FirstGroupName = GetFirstGroupName( g );
    std::cout << FirstGroupName << "\n";
    g = g.openGroup(FirstGroupName);
  }
  else
    std::cout << GroupName << "\n";
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
                      const std::string &GroupName, const std::string &ObjectName )
{
  std::cout << "Reading real array from " << FileName;
  H5::H5File f(FileName, H5F_ACC_RDONLY);
  const bool bFindGroupName{ GroupName.empty() };
  H5::Group g = f.openGroup( bFindGroupName ? std::string("/") : GroupName );
  std::cout << ", group ";
  if( bFindGroupName ) {
    std::string FirstGroupName = GetFirstGroupName( g );
    std::cout << FirstGroupName << "\n";
    g = g.openGroup(FirstGroupName);
  }
  else
    std::cout << GroupName << "\n";
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

enum ExtractFilenameReturn {Good, Bad, No_trajectory};

// Extract the contraction name and trajectory number from filename
static ExtractFilenameReturn ExtractFilenameParts(const std::string &Filename, std::string &Contraction, int &traj)
{
  ExtractFilenameReturn r{Bad};
  Contraction.clear();
  traj = 0;
  auto pos = Filename.find_last_of('/');
  if( pos == std::string::npos )
    pos = 0;
  else
    pos++;
  auto delim = Filename.find_first_of('.', pos);
  if( delim != std::string::npos )
  {
    Contraction = Filename.substr(pos, delim - pos);
    pos = delim + 1;
    delim = Filename.find_first_of('.', pos);
    if( delim != std::string::npos && delim != pos )
    {
      r = Good;
      while( r == Good && pos < delim )
      {
        auto c = Filename[pos++] - '0';
        if( c >=0 && c < 10 )
          traj = traj * 10 + c;
        else
          r = No_trajectory;
      }
    }
  }
  return r;
}

Manifest::Manifest(const std::vector<std::string> &Args, const std::string &sIgnore)
{
  // Get the list of files to ignore
  std::vector<std::string> Ignore;
  for( std::size_t start = 0; start < sIgnore.length() ; ) {
    auto end = sIgnore.find(':', start);
    std::size_t SepLen;
    if( end == std::string::npos ) {
      SepLen = 0;
      end = sIgnore.length();
    }
    else
      SepLen = 1;
    if( end > start )
      Ignore.push_back( sIgnore.substr( start, end - start ) );
    start = end + SepLen;
  }
  // Now walk the list of arguments.
  // Any file that's not in the ignore list gets added to the manifest
  if( Args.size() == 0 )
    return;
  bool parsed = true;
  std::map<std::string, TrajList> & Contractions = (* this);
  for( const std::string &Filename : Args )
  {
    // See whether this file is in the ignore list
    std::size_t iIgnore = 0;
    while( iIgnore < Ignore.size() && Ignore[iIgnore].compare(Filename) )
      iIgnore++;
    if( iIgnore < Ignore.size() )
      std::cout << "Ignoring " << Filename << std::endl;
    else if( !FileExists(Filename))
    {
      parsed = false;
      std::cout << "Error: " << Filename << " doesn't exist" << std::endl;
    }
    else
    {
      std::string Contraction;
      int         traj;
      switch( ExtractFilenameParts( Filename, Contraction, traj ) )
      {
        case Good:
        {
          auto itc = Contractions.find( Contraction );
          if( itc == Contractions.end() )
            itc = Contractions.emplace( Contraction, TrajList( Contraction ) ).first;
          TrajList & cl{ itc->second };
          if( cl.TrajFile.size() == 0 )
            cl.Name = Contraction;
          auto it = cl.TrajFile.find( traj );
          if( it == cl.TrajFile.end() ) {
            cl.TrajFile.emplace( traj, new TrajFile( Filename ));
          }
          else
          {
            const Common::TrajFile &tf{ *it->second };
            if( !Filename.compare( tf.Filename ) )
              std::cout << "Ignoring repetition of " << Filename << std::endl;
            else
            {
              parsed = false;
              std::cout << "Error " << Filename << " conflicts with " << tf.Filename << std::endl;
            }
          }
        }
          break;
          
        case No_trajectory:
          std::cout << "Ignoring non-numeric trajectory " << Filename << std::endl;
          break;
          
        default:
          parsed = false;
          std::cout << "Error: " << Filename << " is not a contraction file" << std::endl;
          break;
      }
    }
  }
  if( !parsed )
    throw std::runtime_error( "Error parsing command-line arguments" );
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

CommandLine::CommandLine( int argc, char *argv[], const std::vector<SwitchDef> &defs )
{
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
    std::cout << (*this) << std::endl;
    exit( EXIT_FAILURE );
  }
  // Now put defaults in for any missing switches
  for( int i = 0; i < defs.size(); i++ ) {
    if( defs[i].Type == Single && defs[i].Default && !GotSwitch( defs[i].Switch ) )
      Switches[defs[i].Switch].push_back( defs[i].Default );
  }
}

END_COMMON_NAMESPACE
