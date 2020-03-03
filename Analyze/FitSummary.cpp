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

#include <set>

using scalar = double;
using Model = Common::Model<scalar>;

using FitTimes = std::pair<int, int>;
struct FitData
{
  scalar ChisqPerDof;
  int ti;
  int tf;
  int dof;
  std::string Parameters;
  int ChiSeq = 0;
  FitData() = default;
  FitData( scalar ChisqPerDof_, int ti_, int tf_, int dof_, const std::string &Parameters_ )
  : ChisqPerDof{ ChisqPerDof_ }, ti{ ti_ }, tf{ tf_ }, dof{ dof_ }, Parameters{ Parameters_ } {}
};

struct Detail
{
  std::map<FitTimes, FitData> Fits;
  int Nt;
  std::string ColumnNames;
  std::string Comments;
  Common::SeedType Seed;
  std::string SeedMachine;
};

using Summary = std::map<std::string, Detail>;

using chisqdof_ti_tf = std::tuple<scalar, int, int>;
using SortByDof = std::map<chisqdof_ti_tf, std::string>;

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
      // Read every model, saving all the data I'll need to reconstruct the summary
      for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
      {
        Common::FileNameAtt n{ sFileName, };
        n.ParseExtra( 2 );
        std::size_t NumExtra{ n.Extra.size() };
        std::string sFitType;
        try
        {
          if( !Common::FileExists( sFileName ) )
            throw std::runtime_error( "doesn't exist" );
          else if( !NumExtra )
            throw std::runtime_error( "doesn't include fit type in name" );
          sFitType = n.Extra[NumExtra - 1];
          for( int i = 0; i < 2 ; i++ )
          {
            std::size_t pos = sFitType.find_last_of( '_' );
            if( pos == std::string::npos )
              break;
            sFitType.resize( pos );
          }
          // The master key for this record is the base filename of the output
          std::string sOutFile{ n.Base };
          sOutFile.append( 1, '.' );
          sOutFile.append( sFitType );
          for( std::size_t i = NumExtra - 1; i-- > 0; )
          {
            sOutFile.append( 1, '.' );
            sOutFile.append( n.Extra[i] );
          }
          Detail & Det{ Sum[sOutFile] };
          // Check the model characteristics
          Model m;
          m.Read( sFileName, " " );
          std::ostringstream ss;
          m.WriteColumnNames( ss );
          if( Det.Fits.empty() )
          {
            Det.Nt = m.Nt();
            Det.ColumnNames = ss.str();
            ss.str("");
            m.SummaryComments( ss );
            ss << "# column(1) is the row order when sorted by ChiSqPerDof\n";
            Det.Comments = ss.str();
            Det.Seed = m.Seed_;
            Det.SeedMachine = m.SeedMachine_;
          }
          else
          {
            if( Det.Nt != m.Nt() || !Common::EqualIgnoreCase( Det.ColumnNames, ss.str() ) )
              throw std::runtime_error( "Fit column names don't match" );
            if( Det.Seed != m.Seed_ || !Common::EqualIgnoreCase( Det.SeedMachine, m.SeedMachine_ ) )
              throw std::runtime_error( "Fit seeds don't match" );
          }
          // Save the summary of the parameters for this file
          ss.str("");
          m.WriteSummaryData( ss );
          Det.Fits.emplace( std::make_pair(std::make_pair( m.ti, m.tf ),
                                           FitData(m.getSummaryData()[Det.Nt-1].Central,
                                                   m.ti, m.tf, m.dof, ss.str())));
        }
        catch( const std::exception &e )
        {
          std::cout << sFileName << " Error: " << e.what() << std::endl;
        }
      }
      // Now process the list of summary files we need to make
      static const std::string Sep{ " " };
      static const std::string sIndent{ "  " };
      for( Summary::iterator it = Sum.begin(); it != Sum.end(); ++it )
      {
        const std::string &sSummaryName{ outBaseFileName + it->first };
        std::cout << sSummaryName << "\n";
        Detail &Det{ it->second };
        // Sort the fits by chi-squared
        std::set<chisqdof_ti_tf> vChiSort;
        for( auto dt = Det.Fits.begin(); dt != Det.Fits.end(); ++dt )
        {
          const FitData &d{ dt->second };
          vChiSort.emplace( d.ChisqPerDof, d.ti, d.tf );
        }
        // Update the fits with sorted sequence number for chi-squared and save sorted file
        std::string sFileName=Common::MakeFilename( sSummaryName, Common::sParams + "_sort", Det.Seed, TEXT_EXT );
        std::ofstream s( sFileName );
        Common::SummaryHeader<scalar>( s, sFileName );
        s << Det.Comments;
        s << "ChiSeq ti tf dof " << Det.ColumnNames << "\n";
        int ChiSeq{ 0 };
        for( auto dt = vChiSort.begin(); dt != vChiSort.end(); ++dt )
        {
          const chisqdof_ti_tf &z{ *dt };
          FitData & d{ Det.Fits[std::make_pair( std::get<1>( z ), std::get<2>( z ) )] };
          d.ChiSeq = ChiSeq++;
          s << d.ChiSeq << Sep << d.ti << Sep << d.tf << Sep << d.dof << Sep << d.Parameters << "\n";
        }
        s.close();
        // Now make the output file
        sFileName=Common::MakeFilename( sSummaryName, Common::sParams, Det.Seed, TEXT_EXT );
        std::size_t t_last = std::numeric_limits<int>::min();
        for( auto dt = Det.Fits.begin(); dt != Det.Fits.end(); ++dt )
        {
          const FitData &d{ dt->second };
          // Write header before each block
          if( t_last != d.ti )
          {
            t_last = d.ti;
            if( !s.is_open() )
            {
              s.open( sFileName );
              Common::SummaryHeader<scalar>( s, sFileName );
              s << Det.Comments;
              s << "# columnheader(1) is the block (index) title\n";
            }
            else
            {
              // two blank lines at start of new data block
              s << "\n\n";
            }
            // Name the data series
            s << "# [ti=" << d.ti << "]\n";
            // Column names, with the series value embedded in the column header (best I can do atm)
            s << "ti=" << d.ti << " ti tf dof " << Det.ColumnNames << "\n";
          }
          // Now write this row
          s << d.ChiSeq << Sep << d.ti << Sep << d.tf << Sep << d.dof << Sep << d.Parameters << "\n";
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
