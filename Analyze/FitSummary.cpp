/*************************************************************************************
 
 Summarise fits
 
 Source file: FitSummary.cpp
 
 Copyright (C) 2020
 
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

//#include <ios>

using scalar = double;
using Model = Common::Model<scalar>;

using FitTimes = std::pair<int, int>;

using Detail = std::map<FitTimes, std::string>;

using Summary = std::map<std::string, Detail>;

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      const std::string inBase{ Common::AppendSlash( cl.SwitchValue<std::string>("i") ) };
      std::string outBaseFileName{ Common::AppendSlash( cl.SwitchValue<std::string>("o") ) };

      bShowUsage = false;
      Summary Sum;
      for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
      {
        Common::FileNameAtt n{ sFileName, };
        n.ParseExtra();
        std::size_t NumExtra{ n.Extra.size() };
        std::string sFitType;
        bool bOK = false;
        int t[2];
        if( !Common::FileExists( sFileName ) )
          std::cout << sFileName << " doesn't exist\n";
        else if( NumExtra )
        {
          sFitType = n.Extra[NumExtra - 1];
          bOK = Common::FileExists( sFileName );
          for( int i = 0; bOK && i < 2 ; i++ )
          {
            bOK = false;
            std::size_t pos = sFitType.find_last_of( '_' );
            if( pos != std::string::npos )
            {
              std::stringstream ss{ sFitType.substr( pos + 1 ) };
              if( ss >> t[i] )
              {
                sFitType.resize( pos );
                bOK = true;
              }
            }
          }
          if( !bOK )
            std::cout << "No timestamp info in " << sFileName << "\n";
        }
        if( bOK )
        {
          std::string sOutFile{ outBaseFileName };
          sOutFile.append( n.Base );
          sOutFile.append( 1, '.' );
          sOutFile.append( sFitType );
          for( std::size_t i = NumExtra - 1; i-- > 0; )
          {
            sOutFile.append( 1, '.' );
            sOutFile.append( n.Extra[i] );
          }
          sOutFile.append( 1, '.' );
          sOutFile.append( Common::sParams );
          sOutFile.append( 1, '.' );
          sOutFile.append( n.SeedString );
          sOutFile.append( 1, '.' );
          sOutFile.append( TEXT_EXT );

          // Add this filename to the correct list
          //std::cout << sOutFile << " " << t[1] << "-" << t[0] << "\n";
          Sum[sOutFile].emplace( std::make_pair( std::make_pair( t[0], t[1] ), sFileName ) );
        }
      }
      // Now make each summary file
      static const std::string Sep{ " " };
      static const std::string sIndent{ "  " };
      for( Summary::iterator it = Sum.begin(); it != Sum.end(); ++it )
      {
        const std::string &sSummaryName{ it->first };
        Detail &Det{ it->second };
        std::cout << "Making " << sSummaryName << "\n";
        std::ofstream s;
        Model m;
        std::size_t tf_last = std::numeric_limits<int>::min();
        for( Detail::iterator dt = Det.begin(); dt != Det.end(); ++dt )
        {
          const FitTimes &Fit{ dt->first };
          const std::string &sModelName{ dt->second };
          const int ti{ Fit.second };
          const int tf{ Fit.first };
          std::cout << sIndent << ti << "-" << tf;
          bool bOK{ false };
          try
          {
            m.Read( sModelName, " " );
            bOK = true;
          }
          catch (const std::exception &e)
          {
            std::cout << sIndent << "Error: " << e.what() << "\n";
          }
          catch (...)
          {
            std::cout << sIndent << "Not a model - ignoring\n";
          }
          if( bOK )
          {
            const std::vector<std::string> &ParamNames{ m.GetColumnNames() };
            const int NumParams{ static_cast<int>( ParamNames.size() ) };
            // Write header before each block
            if( tf_last != tf )
            {
              tf_last = tf;
              if( !s.is_open() )
              {
                s.open( sSummaryName );
                Common::SummaryHeader<scalar>( s, sSummaryName );
                m.SummaryComments( s );
              }
              else
              {
                // two blank lines at start of new data block
                s << "\n" << std::endl;
              }
              // Name the data series
              s << "# [tf=" << tf << "]" << std::endl;
              // Column names, with the series value embedded in the column header (best I can do atm)
              s << "tf=" << tf << Sep << "ti";
              for( int p = 0; p < NumParams; p++ )
                s << Sep << ParamNames[p] << Sep << ParamNames[p] << "_low" << Sep << ParamNames[p]
                  << "_high" << Sep << ParamNames[p] << "_check";
              s << " ChiSq Dof" << std::endl;
            }
            // Now write this row
            const Common::ValWithEr<scalar> * pSum{ m.getSummaryData() };
            s << tf << Sep << ti;
            for( int p = 0; p < NumParams; p++ )
              s << Sep << *pSum++;
            s << Sep << ( pSum[-1].Central * m.dof ) << Sep << m.dof << std::endl;
          }
        }
      }
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
    " <options> Model1 [Model2 ...]\n"
    "Summarise the fits contained in each model (ready for plotting), where <options> are:\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "Flags:\n"
    "--help     This message\n";
  }
  return iReturn;
}
