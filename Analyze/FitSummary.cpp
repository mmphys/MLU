/*************************************************************************************
 
 Summarise fits
 
 Source file: FitSummary.cpp
 
 Copyright (C) 2020-2022
 
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

#include "FitSummary.hpp"

bool bReverseSort{ false };

const std::string &Sep{ Common::Space };
const std::string &NL{ Common::NewLine };
const std::string sIndent{ "  " };
const std::string sError{ "Error: " };

std::string FitTimes::tfLabel() const
{
  std::string s;
  for( std::size_t i = 1; i < time.size(); ++i )
  {
    if( i != 1 )
      s.append( 1, '_' );
    s.append( std::to_string( time[i] ) );
  }
  return s;
}

// The old file format didn't have any information containing the timeslices in the fit
// In this case I guess based on the filename and the number of files
unsigned int FitTimes::GuessOldNumDataPoints( int NumFiles ) const
{
  unsigned int Num{ 0 };
  for( std::size_t i = 0; i < time.size(); i += 2 )
    Num += time[i+1] - time[i] + 1;
  if( time.size() == 2 )
    Num *= NumFiles;
  return Num;
}

bool FitTimes::Parse( std::string Times )
{
  for( auto &c : Times )
    if( c == '_' )
      c = ' ';
  try { time = Common::ArrayFromString<int>( Times ); } catch(...) {}
  if( time.size() & 1 )
    time.clear();
  return !time.empty();
}

std::ostream &operator<<( std::ostream &s, const FitTimes &ft )
{
  for( auto t : ft.time )
    s << Sep << t;
  return s << Sep << ft.tfLabel();
}

bool operator<( const FitTimes &lhs, const FitTimes &rhs )
{
  if( lhs.time.empty() || lhs.time.size() != rhs.time.size() )
    throw std::runtime_error( "Mismatched numbers of fit ranges" );
  for( std::size_t i = 0; i < lhs.time.size(); ++i )
    if( lhs.time[i] != rhs.time[i] )
      return lhs.time[i] < rhs.time[i];
  return false;
}

std::ostream &operator<<( std::ostream &s, const FitData &d )
{
  return s << d.Seq << d.ft << Sep << d.NumDataPoints << Sep << d.dof << Sep << d.SampleSize << Sep << d.Parameters << NL;
}

bool operator<( const TestStatKey &lhs, const TestStatKey &rhs )
{
  if( lhs.TestStatistic != rhs.TestStatistic )
  {
    if( bReverseSort )
      return lhs.TestStatistic > rhs.TestStatistic;
    return lhs.TestStatistic < rhs.TestStatistic;
  }
  return lhs.ft < rhs.ft;
}

// Parse the command line, splitting into separate lists for each base
Summariser::Summariser( const Common::CommandLine &cl )
: inBase{ Common::AppendSlash( cl.SwitchValue<std::string>("i") ) },
  outBaseFileName{ Common::AppendSlash( cl.SwitchValue<std::string>("o") ) },
  StatisticName{ cl.SwitchValue<std::string>( "stat" ) }
{
  bReverseSort = !cl.GotSwitch( "inc" );
  Common::MakeAncestorDirs( outBaseFileName );
  for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
  {
    Common::FileNameAtt n{ sFileName, };
    std::size_t NumExtra{ n.Extra.size() };
    if( !Common::FileExists( sFileName ) )
      throw std::runtime_error( sFileName + " doesn't exist" );
    if( NumExtra == 0 )
      throw std::runtime_error( "No fit type in " + sFileName );
    std::string sFitType{ n.Extra[0] };
    const std::size_t pos_ti = sFitType.find_first_of( '_' );
    FitTimes ft;
    if( pos_ti == std::string::npos || !ft.Parse( sFitType.substr( pos_ti + 1 ) ) )
      throw std::runtime_error( "Error extracting fit ranges from " + sFileName );
    sFitType.resize( pos_ti );
    // The master key for this record is the base filename of the output
    std::string sOutFile{ n.Base };
    sOutFile.append( 1, '.' );
    sOutFile.append( sFitType );
    n.AppendOps( sOutFile );
    lBase[sOutFile].emplace_back( FileInfo( ft, sFileName ) );
  }
}

// Process every base name and it's list of models separately
void Summariser::Run()
{
  for( BaseList::iterator it = lBase.begin(); it != lBase.end(); ++it )
  {
    const std::string &sOutFile{ it->first };
    const std::string &sSummaryName{ outBaseFileName + sOutFile };
    std::cout << sSummaryName << NL;
    // Read every model, saving all the data I'll need to reconstruct the summary
    int NumSamples{ 0 };
    std::string ColumnNames;
    std::string Comments;
    std::vector<std::string> FileNameOps; // For now I get these, but don't check them
    Common::SeedType Seed{ 0 };
    int ModelNum = 0;
    std::array<Model, 2> Models;
    FitMap Fits;
    for( FileInfo &ThisFile : it->second )
    {
      try
      {
        // Check the model characteristics
        Model &m{ Models[ModelNum] };
        m.SetName( ThisFile.FileName, &FileNameOps );
        m.Read( "" );
        if( m.dof < 0 )
          throw std::runtime_error("dof=" + std::to_string( m.dof ) + " (<0 invalid)");
        std::ostringstream ss;
        m.SummaryColumnNames( ss );
        if( ModelNum == 0 )
        {
          ModelNum++;
          ColumnNames = ss.str();
          ss.str("");
          m.SummaryComments( ss );
          ss << "# column(1) is the row order when sorted by pValueH (Hotelling p-value)\n";
          Comments = ss.str();
          Seed = m.Name_.Seed;
          NumSamples = m.NumSamples();
        }
        else
        {
          Models[0].IsCompatible( m, &NumSamples );
          if( !Common::EqualIgnoreCase( ColumnNames, ss.str() ) )
            throw std::runtime_error( "Fit column names don't match" );
        }
        // Save the summary of the parameters for this file
        ss.str("");
        m.SummaryContents( ss );
        const Common::ValWithEr<scalar> &pValueH{ m.getSummaryData( StatisticName ) };
        scalar TestStat = pValueH.Central;
        // save the model in my list
        int NumDataPoints{ 0 };
        for( const std::vector<int> &v : m.FitTimes )
          NumDataPoints += static_cast<int>( v.size() );
        FitData fd(TestStat,ThisFile.ft,NumDataPoints,m.dof,m.SampleSize,ss.str(),Models.size());
        Fits.emplace( std::make_pair( ThisFile.ft, fd ) );
      }
      catch( const std::exception &e )
      {
        //std::cout << sFileName << Sep << sError << e.what() << NL;
        std::cout << sIndent << sError << e.what() << NL;
      }
    }
    // Update the fits with sorted sequence number for test statistic and save sorted file
    {
      // Sort the fits by selected test statistic
      std::set<TestStatKey> SortSet;
      for( FitMap::value_type &dt : Fits )
        SortSet.emplace( TestStatKey{ dt.second.Stat, dt.first } );
      // Now write this out to file
      std::string sFileName=Common::MakeFilename(sSummaryName,Common::sParams+"_sort",Seed,TEXT_EXT);
      std::ofstream s( sFileName );
      Common::SummaryHeader<scalar>( s, sFileName );
      s << Comments;
      s << "Seq " << ColumnNames << NL;
      int SortSeq{ 0 };
      for( const TestStatKey &z : SortSet )
      {
        FitData & d{ Fits[z.ft] };
        d.Seq = SortSeq++;
        s << d.Seq << Common::Space << d.Parameters;
      }
    }
    // Now make the output file with blocks for each ti, sorted by tf
    {
      std::ofstream s;
      const std::string sFileName{Common::MakeFilename(sSummaryName,Common::sParams,Seed,TEXT_EXT)};
      std::size_t t_last = std::numeric_limits<int>::min();
      for( auto dt = Fits.begin(); dt != Fits.end(); ++dt )
      {
        const FitData &d{ dt->second };
        // Write header before each block
        if( t_last != d.ft.ti() )
        {
          t_last = d.ft.ti();
          if( !s.is_open() )
          {
            s.open( sFileName );
            Common::SummaryHeader<scalar>( s, sFileName );
            s << Comments;
            s << "# columnheader(1) is the block (index) title\n";
          }
          else
          {
            // two blank lines at start of new data block
            s << "\n\n";
          }
          // Name the data series
          s << "# [ti=" << d.ft.ti() << "]\n";
          // Column names, with the series value embedded in the column header (best I can do atm)
          s << "ti=" << d.ft.ti() << Common::Space << ColumnNames << NL;
        }
        // Now write this row
        s << d.Seq << Common::Space << d.Parameters;
      }
    }
  }
}

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
      {"stat", CL::SwitchType::Single, Common::sPValueH.c_str() },
      {"inc",  CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      Summariser Sum( cl );
      bShowUsage = false;
      Sum.Run();
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << sError << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << sError << "Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> Model1 [Model2 ...]\n"
    "Summarise the fits contained in each model (ready for plotting), where <options> are:\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "--stat Statistic (default: " << Common::sPValueH << ")\n"
    "Flags:\n"
    "--inc  Sort increasing (default sort test stat decreasing)\n"
    "--help This message\n";
  }
  return iReturn;
}
