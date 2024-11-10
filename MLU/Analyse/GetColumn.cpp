/*************************************************************************************
 
 Extract a column from an HDF5 file, e.g. a parameter from a fit (Model<>) file

 Source file: GetColumn.cpp
 
 Copyright (C) 2019-2022
 
 Author: Michael Marshall <Mike@lqcd.me>

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

#include <MLUconfig.h>
#include <MLU/MLU.hpp>

struct Extractor
{
  using Scalar = double;
  using ValErType = MLU::ValWithEr<Scalar>;
  static const ValErType veUnknown;

  struct NameValue
  {
    std::vector<std::string> Name;
    std::vector<ValErType> Value;
  };
  const std::string inBase;
  const std::vector<std::string> ExactNames;
  const std::vector<std::string> PartialNames;
  const std::string DefaultGroup;
  const int ErrDig;
  Extractor( MLU::CommandLine &cl );
  bool Run( const std::vector<std::string> &Args );
protected:
  NameValue ReadFile( const std::string &Filename );
  bool Show( const NameValue &NV, const std::vector<std::string> Selection, bool bExact ) const;
  void ShowOne( const std::string &Name, const ValErType &Value ) const;
};

const Extractor::ValErType Extractor::veUnknown( 0, 0, 0, 0, 0, 0 );

Extractor::Extractor( MLU::CommandLine &cl )
: inBase{ cl.SwitchValue<std::string>("i") },
  ExactNames{ MLU::ArrayFromString( cl.SwitchValue<std::string>( "exact" ) ) },
  PartialNames{ MLU::ArrayFromString( cl.SwitchValue<std::string>( "partial" ) ) },
  DefaultGroup{ cl.SwitchValue<std::string>( "group" ) },
  ErrDig{ cl.SwitchValue<int>( "errdig" ) }
{
}

bool Extractor::Run( const std::vector<std::string> &Args )
{
  bool bOK{ true };
  for( const std::string &s : MLU::glob( Args.begin(), Args.end(), inBase.c_str() ) )
  {
    try
    {
      NameValue NV{ ReadFile( s ) };
      if( ExactNames.empty() && PartialNames.empty() )
        for( std::size_t i = 0; i < NV.Name.size(); ++i )
          ShowOne( NV.Name[i], NV.Value[i] );
      else
      {
        if( !ExactNames.empty() && !Show( NV, ExactNames, true ) )
          bOK = false;
        if( !PartialNames.empty() && !Show( NV, PartialNames, false ) )
          bOK = false;
      }
    }
    catch(const std::exception &e)
    {
      bOK = false;
      std::cerr << "Error: " << e.what() << ". File: " << s << std::endl;
    }
  }
  return bOK;
}

// Show exact or partial matches
bool Extractor::Show( const NameValue &NV, const std::vector<std::string> Selection,
                      bool bExact ) const
{
  bool bOK{ Selection.empty() };
  for( const std::string &s : Selection )
  {
    std::size_t Count{};
    const bool bLookingForDigit{ !s.empty() && std::isdigit( s.back() ) };
    // If I ask for xxx0, also show xxx
    const bool bLookingForZero{ !s.empty() && s.back() == '0' };
    const std::string sWithoutZero{ bLookingForZero ? s.substr( 0, s.size() - 1 ) : "" };
    for( std::size_t i = 0; i < NV.Name.size(); ++i )
    {
      const std::string &ni{ NV.Name[i] };
      const std::string NameWildcard{ !bLookingForDigit && std::isdigit( ni.back() )
                                      ? ni.substr( 0, ni.size() - 1 ) : "" };
      bool bFound{ false };
      if( bExact )
      {
        // The whole string must match - with or without trailing zero
        bFound = MLU::EqualIgnoreCase( s, ni );
        if( !bFound )
        {
          if( !sWithoutZero.empty() )
            bFound = MLU::EqualIgnoreCase( sWithoutZero, ni );
          else if( !NameWildcard.empty() )
            bFound = MLU::EqualIgnoreCase( s, NameWildcard );
        }
      }
      else
      {
        // Either I find what I'm looking for anywhere in the string
        bFound = ni.find( s ) != std::string::npos;
        // ... or the end matches without the trailing zero
        if( !bFound && bLookingForZero && ni.size() >= sWithoutZero.size() )
          bFound = MLU::EqualIgnoreCase( sWithoutZero,
                                            ni.substr( ni.size() - sWithoutZero.size() ) );
      }
      if( bFound )
      {
        ShowOne( ni, NV.Value[i] );
        ++Count;
        bOK = true;
      }
    }
    if( !Count )
      ShowOne( s, veUnknown );
  }
  return bOK;
}

void Extractor::ShowOne( const std::string &Name, const ValErType &Value ) const
{
  std::cout << Name << MLU::Space << Value.to_string( ErrDig ) << MLU::Space
            << Value << MLU::NewLine;
}

Extractor::NameValue Extractor::ReadFile( const std::string &Filename )
{
  if( !MLU::FileExists( Filename ) )
    throw std::runtime_error( Filename + " doesn't exist" );
  ::H5::H5File f;
  ::H5::Group g;
  std::string GroupName{ DefaultGroup };
  MLU::H5::OpenFileGroup( f, g, Filename, nullptr, &GroupName );
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  std::vector<std::string> ColumnNames;
  std::vector<ValErType> ValEr;
  try
  {
    try
    {
      ::H5::Attribute a{ g.openAttribute( MLU::sColumnNames ) };
      ColumnNames = MLU::H5::ReadStrings( a );
    }
    catch(const ::H5::Exception &)
    {
      throw std::runtime_error( MLU::sColumnNames + " not available in " + GroupName );
    }
    bool bOK{ false };
    try
    {
      ::H5::DataSet ds{ g.openDataSet( MLU::sSummaryDSName ) };
      ::H5::DataSpace dsp{ ds.getSpace() };
      int nDims{ dsp.getSimpleExtentNdims() };
      if( nDims == 2 )
      {
        hsize_t Dim[2];
        dsp.getSimpleExtentDims( Dim );
        const std::size_t ValSize{ static_cast<std::size_t>( Dim[0] * Dim[1] ) };
        if( Dim[1] == ColumnNames.size() && ValSize == Dim[0] * Dim[1] )
        {
          ValEr.resize( ValSize );
          ds.read( ValEr.data(), MLU::H5::Equiv<ValErType>::Type );
          ValEr.resize( ColumnNames.size() );
          bOK = true;
        }
      }
    }
    catch(const ::H5::Exception &)
    {
    }
    if( !bOK )
      throw std::runtime_error( MLU::sSummaryDSName + " not available in " + GroupName );
  }
  catch(...)
  {
    ::H5::Exception::clearErrorStack();
    ::H5::Exception::setAutoPrint(h5at, f5at_p);
    throw;
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  
  return NameValue{ ColumnNames, ValEr };
}

int main(int argc, const char *argv[])
{
  static const char DefaultErrDig[] = "2";
  std::ios_base::sync_with_stdio( false );
  int iReturn = EXIT_SUCCESS;
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      // Fitter parameters
      {"i", CL::SwitchType::Single, "" },
      {"exact", CL::SwitchType::Single, ""},
      {"partial", CL::SwitchType::Single, ""},
      {"group", CL::SwitchType::Single, ""},
      {"errdig", CL::SwitchType::Single, DefaultErrDig},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    Extractor Fan{ cl };
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      bShowUsage = false;
      if( !Fan.Run( cl.Args ) )
        iReturn = EXIT_FAILURE;
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
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " [Options] files\n"
    "Extract data from specified column names from files\n"
    "Options:\n"
    "-i        Input prefix\n"
    "--exact   Comma separated list of column names which must match exactly\n"
    "--partial Comma separated list of column names which must partially match\n"
    "--group   HDF5 group to read from (default: first)\n"
    "--errdig  Number of significant figures in error (default " << DefaultErrDig << ")\n"
    "Flags:\n"
    "--help    This message\n";
  }
  return iReturn;
}
