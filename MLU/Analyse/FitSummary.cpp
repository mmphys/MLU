/**
 
 Summarise fits
 
 Source file: FitSummary.cpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "FitSummary.hpp"

const std::string &Sep{ MLU::Space };
const std::string &NL{ MLU::NewLine };
const std::string sIndent{ "  " };
const std::string sError{ "Error: " };

void FitTimes::tiLabelSuffix( std::string &s ) const
{
  for( std::size_t i = 2; i < time.size(); ++i )
  {
    s.append( 1, '_' );
    s.append( std::to_string( time[i] ) );
  }
}

std::string FitTimes::tiLabel() const
{
  std::string s;
  if( time.empty() )
    s = std::string( 1, '0' );
  else
  {
    s = std::to_string( time[0] );
    tiLabelSuffix( s );
  }
  return s;
}

std::string FitTimes::tfLabel() const
{
  std::string s;
  if( time.size() < 1 )
    s = std::string( 1, '0' );
  else
  {
    s = std::to_string( time[1] );
    tiLabelSuffix( s );
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

bool FitTimes::ParseTimes( std::string Times )
{
  for( auto &c : Times )
    if( c == '_' )
      c = ' ';
  try { time = MLU::ArrayFromString<int>( Times ); } catch(...) {}
  if( time.size() & 1 )
    time.clear();
  return !time.empty();
}

bool FitTimes::ParseTypeTimes( std::string &TypeTimes )
{
  bool bOK{ false };
  const std::size_t pos{ TypeTimes.find_first_of( '_' ) };
  if( pos != std::string::npos && pos < TypeTimes.length() - 1 )
  {
    bOK = ParseTimes( TypeTimes.substr( pos + 1 ) );
    if( bOK )
      TypeTimes.resize( pos );
  }
  return bOK;
}

bool FitTimes::operator<( const FitTimes &rhs ) const
{
  // First sort by the values of the times in common
  for( std::size_t i = 0; i < std::min( time.size(), rhs.time.size() ); ++i )
    if( time[i] != rhs.time[i] )
      return time[i] < rhs.time[i];
  // All other things equal, the shorter sequence sorts first
  return time.size() < rhs.time.size();
}

/*std::ostream &operator<<( std::ostream &s, const FitTimes &ft )
{
  for( auto t : ft.time )
    s << Sep << t;
  return s << Sep << ft.tfLabel();
}*/

/*std::ostream &operator<<( std::ostream &s, const FitData &d )
{
  return s << d.Seq << d.ft << Sep << d.NumDataPoints << Sep << d.dof << Sep << d.SampleSize << Sep << d.Parameters << NL;
}*/

void FitData::WriteLongFormat( std::ostream &os, const FitTimes &ft ) const
{
  os << Seq << MLU::Space << idxLoadOrder << MLU::Space;
  if( ft.NumFitTimes() == 0 )
    os << "0 0 0 0";
  else
  {
    os << ft.ti( 0 ) << MLU::Space << ft.tf( 0 )
       << MLU::Space << ft.tiLabel() << MLU::Space << ft.tfLabel();
  }
  os << MLU::Space << NumDataPoints;
  os << MLU::Space << dof;
  os << MLU::Space << SampleSize;
  os << MLU::Space << CovarSampleSize;
  for( const veScalar &v : Parameters )
    os << MLU::Space << v;
  for( int i = 1; i < ft.NumFitTimes(); ++i )
    os << MLU::Space << ft.ti( i ) << MLU::Space << ft.tf( i );
  os << MLU::NewLine;
}

void FitData::WriteTableFormat( std::ostream &os, const FitTimes &ft, std::size_t StatIndex,
                                unsigned char ErrorDigits ) const
{
  for( int i = 0; i < ft.NumFitTimes(); ++i )
  {
    if( i )
      os << MLU::CommaSpace;
    os << "$ \\left[ " << ft.ti( i );
    if( ft.tf( i ) != ft.ti( i ) )
      os << "\\text{-}" << ft.tf( i );
    os << " \\right] $";
  }
  for( std::size_t i = 0; i < Parameters.size(); ++i )
  {
    const veScalar &v{ Parameters[i] };
    os << " & ";
    if( v.Check )
    {
      std::string s{ v.to_string( ErrorDigits ) };
      if( i >= StatIndex )
      {
        // Don't show the error bar on the statistic - meaningless
        const std::size_t pos{ s.find( '(' ) };
        if( pos != std::string::npos )
          s.resize( pos );
      }
      os << s;
    }
  }
  os << "\\\\\n";
}

bool TestStatKey::bReverseSort{ false };

bool TestStatKey::operator<( const TestStatKey &rhs ) const
{
  if( TestStatistic != rhs.TestStatistic )
  {
    if( bReverseSort )
      return TestStatistic > rhs.TestStatistic;
    return TestStatistic < rhs.TestStatistic;
  }
  return ft < rhs.ft;
}

// Parse the command line, splitting into separate lists for each base
Summariser::Summariser( const MLU::CommandLine &cl )
: inBase{ MLU::AppendSlash( cl.SwitchValue<std::string>("i") ) },
  outBaseFileName{ MLU::AppendSlash( cl.SwitchValue<std::string>("o") ) },
  StatisticName{ cl.SwitchValue<std::string>( "stat" ) },
  Strictness{ cl.SwitchValue<int>( "strict" ) },
  MonotonicUpperLimit{ cl.SwitchValue<scalar>( "maxE" ) },
  bAll{ cl.GotSwitch( "all" ) },
  bFast{ cl.GotSwitch( "fast" ) },
  bTableN{ cl.GotSwitch( "tabn" ) },
  ErrorDigits{ static_cast<unsigned char>( cl.SwitchValue<unsigned int>( "errdig" ) ) },
  TableRowsPerPage{ cl.SwitchValue<std::size_t>( "tablen" ) },
  vIgnoreMomenta{ MLU::ArrayFromString( cl.SwitchValue<std::string>("pignore") ) },
  ParamSubset{ MLU::ArrayFromString( cl.SwitchValue<std::string>("params") ) }
{
  if( !ParamSubset.empty() && bAll )
    throw std::runtime_error( "--all and --params can't be used together" );
  MLU::NoDuplicates( vIgnoreMomenta, "Ignored momenta", 0 );
  if( Strictness < -1 || Strictness > 3 )
    throw std::runtime_error( "--strict " + std::to_string( Strictness ) + " invalid" );
  if( bAll && bFast )
    throw std::runtime_error( "--fast implies all models have same parameters. --all incompatible" );
  TestStatKey::bReverseSort = !cl.GotSwitch( "inc" );
  MLU::MakeAncestorDirs( outBaseFileName );
  int SeqNum{ 0 };
  for( const std::string &sFileName : MLU::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
  {
    MLU::FileNameAtt n{ sFileName, nullptr, &vIgnoreMomenta };
    if( !MLU::FileExists( sFileName ) )
      throw std::runtime_error( sFileName + " doesn't exist" );
    FitTimes ft;
    std::string sFitType;
    bool bIsFit{ false };
    if( n.Extra.size() )
    {
      sFitType = n.Extra[0];
      bIsFit = ft.ParseTypeTimes( sFitType );
    }
    if( !bIsFit )
    {
      // If there are no fit times, give every item a unique number so it can be sorted
      ft.time.resize( 2 );
      ft.time[0] = 0;
      ft.time[1] = ++SeqNum;
    }
    // The master key for this record is the base filename of the output
    std::string sOutFile{ n.Base };
    if( bIsFit && sFitType.size() )
    {
      sOutFile.append( 1, '.' );
      sOutFile.append( sFitType );
    }
    n.AppendOps( sOutFile, MLU::Period );
    std::cout << sOutFile << MLU::Space << sFileName << MLU::NewLine;
    // Have we seen this base before?
    BaseList::iterator it{ lBase.find( sOutFile ) };
    if( it == lBase.end() )
    {
      const MLU::Momentum &np{ n.p.empty() ? MLU::p0 : n.GetFirstNonZeroMomentum().second };
      it = lBase.emplace( sOutFile, BaseInfo( np ) ).first;
    }
    it->second.vFI.emplace_back( ft, sFileName, bIsFit );
  }
}

bool Summariser::ReadModel( Model &m, FileInfoIterator &it, std::vector<FileInfo> &Files,
                            std::vector<std::string> &FileNameOps, bool bShow )
{
  m.SetName( MLU::FileNameAtt{ it->FileName, &FileNameOps, &vIgnoreMomenta } );
  try
  {
    m.Read( bShow ? " " : nullptr );
    // Check the model characteristics
    if( m.dof < 0 )
      throw std::runtime_error("dof=" + std::to_string( m.dof ) + " (<0 invalid)");
    if( !m.CheckParameters( Strictness, MonotonicUpperLimit ) )
      throw std::runtime_error( "Parameter(s) compatible with 0 and/or not unique (strictness "
                               + std::to_string( Strictness ) + ")");
  }
  catch( const std::exception &e )
  {
    std::cout << sIndent << sError << e.what() << NL;
    it = Files.erase( it );
    return false;
  }
  ++it;
  return true;
}

bool Summariser::GetCommonParameters( std::vector<FileInfo> &Files, bool bMaximum )
{
  bool bFirst{ true };
  for( FileInfoIterator itFI = Files.begin(); itFI != Files.end(); )
  {
    if( bFirst )
    {
      bFirst = false;
      ++itFI;
    }
    else if( ReadModel( Models[1], itFI, Files, FileNameOps, true ) )
    {
      // Subsequent model
      if( MaxFitTimes < Models[1].FitTimes.size() )
        MaxFitTimes = Models[1].FitTimes.size();
      if( bAll )
        Params.Merge( Models[1].params );
      else
        Params.KeepCommon( Models[1].params );
      // Remove statistics columns not present in this model
      using UNS = MLU::UniqueNameSet;
      const UNS mStat{ Models[1].GetStatColumnNames() };
      if( bAll )
      {
        for( const std::string &s : mStat )
          StatColumnNames.insert( s ); // I'm merging - add this to my list of names
      }
      else
      {
        for( UNS::iterator it = StatColumnNames.begin(); it != StatColumnNames.end(); )
        {
          if( mStat.find( *it ) == mStat.end() )
            it = StatColumnNames.erase( it ); // Not in the other list - delete my copy
          else
            ++it;
        }
      }
    }
  }
  // Not finding any common parameters on the first pass is an error
  if( Params.empty() )
    std::cout << sIndent << sError << "No common parameters" << NL;
  else
    Params.AssignOffsets(); // Because we might have deleted a few
  return !Params.empty();
}

void Summariser::BuildFitMap( std::vector<FileInfo> &Files )
{
  const unsigned int CompatFlags{ bDoPassOne ? MLU::COMPAT_DISABLE_NT
                                             : MLU::COMPAT_DEFAULT };
  int ModelNum = 0;
  std::size_t idxLoadOrder{ 0 };
  Fits.clear();
  for( FileInfoIterator itFI = Files.begin(); itFI != Files.end(); )
  {
    Model &m{ Models[ModelNum] };
    FileInfo &ThisFile{ *itFI };
    if( ModelNum == 0 )
      ++itFI;
    if( ModelNum == 0 || ReadModel( m, itFI, Files, FileNameOps, !bDoPassOne ) )
    {
      try
      {
        if( ModelNum == 0 )
        {
          ++ModelNum;
          std::ostringstream ss;
          m.SummaryComments( ss );
          ss << "# column(1) is the row order when sorted by " << StatisticName << MLU::Space
             << ( TestStatKey::bReverseSort ? "de" : "a" ) << "scending\n";
          Comments = ss.str();
          Seed = m.Name_.Seed;
          NumSamples = m.NumSamples();
        }
        else
        {
          Models[0].IsCompatible( m, &NumSamples, CompatFlags, false );
          if( !bDoPassOne && ( Params != m.params || StatColumnNames != m.GetStatColumnNames() ) )
            throw std::runtime_error( "Fit column names don't match" );
        }
        // Get the test statistic
        const int idxpValue{ m.GetColumnIndexNoThrow( StatisticName ) };
        if( idxpValue < 0 && ThisFile.bIsFit )
          throw std::runtime_error( StatisticName + " unavailable - required for fit" );
        scalar TestStat = idxpValue < 0 ? 1 : m.SummaryData(idxpValue).Central;
        // save the model in my list
        const int NumDataPoints{ static_cast<int>( MLU::GetExtent( m.FitTimes ) ) };
        FitData fd( TestStat, NumDataPoints, m.dof, m.SampleSize, m.CovarSampleSize,
                    idxLoadOrder++, m.GetValWithEr( Params, StatColumnNames ) );
        Fits.emplace( std::make_pair( ThisFile.ft, fd ) );
      }
      catch( const std::exception &e )
      {
        std::cout << sIndent << sError << e.what() << NL;
      }
    }
  }
}

void Summariser::SummaryColumnNames( std::ostream &os, std::size_t NumFitTimes,
              const MLU::Params &ParamNames, const MLU::UniqueNameSet &StatNames ) const
{
  os << "Load " << Model::SummaryColumnPrefix << MLU::Space;
  ParamNames.WriteNamesValWithEr( os );
  for( const std::string &Name : StatNames )
  {
    os << MLU::Space;
    veScalar::Header( Name, os, MLU::Space );
  }
  // Write additional fit time pair names
  for( std::size_t i = 1; i < NumFitTimes; ++i )
    os << MLU::Space << "ti" << i << MLU::Space << "tf" << i;
}

// Write sorted list. Must be done first because this fills in the sort sequence number
void Summariser::WriteSorted( const std::string &sFileName, const std::set<TestStatKey> &SortSet,
                              bool SetSeq )
{
  // Create file and write header
  std::ofstream s( sFileName );
  MLU::SummaryHeader<scalar>( s, sFileName );
  s << Comments << "Seq ";
  SummaryColumnNames( s, MaxFitTimes, Params, StatColumnNames );
  s << NL;
  // Update the fits with sorted sequence number and write each row
  int SortSeq{ 0 };
  for( const TestStatKey &z : SortSet )
  {
    FitData & d{ Fits[z.ft] };
    if( SetSeq )
      d.Seq = SortSeq++;
    d.WriteLongFormat( s, z.ft );
  }
}

// Now make the output file with blocks for each ti, sorted by tf
void Summariser::WriteUnsorted( const std::string &sFileName ) const
{
  std::ofstream s;
  std::size_t t_last = std::numeric_limits<std::size_t>::max();
  for( const FitMap::value_type &dt : Fits )
  {
    const FitTimes &ft{ dt.first };
    const FitData &d{ dt.second };
    // Write header before each block
    if( t_last != ft.ti() )
    {
      t_last = ft.ti();
      if( !s.is_open() )
      {
        s.open( sFileName );
        MLU::SummaryHeader<scalar>( s, sFileName );
        s << Comments;
        s << "# columnheader(1) is the block (index) title\n";
      }
      else
      {
        // two blank lines at start of new data block
        s << "\n\n";
      }
      // Name the data series
      s << "# [ti=" << ft.ti() << "]\n";
      // Column names, with the series value embedded in the column header (best I can do atm)
      s << "ti=" << ft.ti() << MLU::Space;
      SummaryColumnNames( s, MaxFitTimes, Params, StatColumnNames );
      s << NL;
    }
    // Now write this row
    d.WriteLongFormat( s, ft );
  }
}

void Summariser::WriteTableHeader( std::ofstream &os ) const
{
  os << "\\begin{tabular}{|c|";
  if( bTableN )
    os << "c|";
  for( std::size_t i = 0; i < NumParams() + NumStats(); i++ )
    os << "c|";
  os << "}\n\\hline\n";
  // Write column names
  if( bTableN )
    os << "$n^2$ & ";
  os << "Fit";
  for( const MLU::Params::value_type &it : Params )
  {
    const MLU::Param &p{ it.second };
    for( std::size_t i = 0; i < p.size; ++i )
      os << " & " << Params.GetName( it, i );
  }
  // Write stat column names
  for( const std::string &s : StatColumnNames )
  {
    os << " & ";
    if( MLU::EqualIgnoreCase( s, MLU::sChiSqPerDof ) )
      os << "$\\flatfrac{\\chi^2}{\\textrm{dof}}$";
    else
      os << s;
  }
  os << "\\\\\n\\hline\n";
}

void Summariser::WriteTableFooter( std::ofstream &os ) const
{
  os << "\\hline\n\\end{tabular}\n";
}

void Summariser::WriteTabular( const std::string &sFileName, const BaseInfo &bi ) const
{
  const std::string sMom{ std::to_string( bi.p.p2() ) };
  // Write header
  std::ofstream os( sFileName );
  os << "{\\tiny\n";
  std::size_t RowNumber{ 0 };
  bool bFirst{ true };
  // Write each row
  for( const FitMap::value_type &dt : Fits )
  {
    if( !RowNumber++ )
    {
      if( bFirst )
        bFirst = false;
      else
        os << MLU::NewLine;
      WriteTableHeader( os );
    }
    const FitTimes &ft{ dt.first };
    const FitData &d{ dt.second };
    // Now write this row
    if( bTableN )
      os << sMom << " & ";
    d.WriteTableFormat( os, ft, NumParams(), ErrorDigits );
    if( RowNumber == TableRowsPerPage )
    {
      WriteTableFooter( os );
      RowNumber = 0;
    }
  }
  if( RowNumber )
    WriteTableFooter( os );
  os << "}\n";
}

// Process every base name and it's list of models separately
void Summariser::Run()
{
  for( BaseList::iterator it = lBase.begin(); it != lBase.end(); ++it )
  {
    const std::string &sSummaryName{ outBaseFileName + it->first };
    std::cout << sSummaryName << NL;
    BaseInfo &bi{ it->second };
    std::vector<FileInfo> &Files{ bi.vFI };
    FileNameOps.clear(); // For now I get these, but don't check them
    // Load the first Model
    {
      bool Model0Loaded{ false };
      for( FileInfoIterator itFI = Files.begin(); !Model0Loaded && itFI != Files.end(); )
        Model0Loaded = ReadModel( Models[0], itFI, Files, FileNameOps, true );
      if( !Model0Loaded )
        continue;
    }
    MaxFitTimes = Models[0].FitTimes.size(); // Maximum number of TI-TF pairs
    if( ParamSubset.empty() )
      Params = Models[0].params;
    else
    {
      std::vector<std::size_t> vIdx;
      Params = MLU::Params( ParamSubset, vIdx );
    }
    StatColumnNames = Models[0].GetStatColumnNames();
    // Optional first pass - build a list of common parameters and statistics columns
    bDoPassOne = !bFast && Files.size() > 1;
    if( bDoPassOne && !GetCommonParameters( Files, false ) )
      continue;
    // Make sure the requested statistic is included
    if( StatColumnNames.find( StatisticName ) == StatColumnNames.end() )
      StatColumnNames.insert( StatisticName );
    // Read every model, saving all the data I'll need to reconstruct the summary
    BuildFitMap( Files );
    if( !Fits.empty() )
    {
      {
        // Sort the fits by selected test statistic
        std::set<TestStatKey> SortSet;
        for( FitMap::value_type &dt : Fits )
          SortSet.emplace( TestStatKey{ dt.second.Stat, dt.first } );
        WriteSorted( MLU::MakeFilename(sSummaryName,MLU::sParams+"_sort",Seed,TEXT_EXT),
                     SortSet, true );
      }
      WriteUnsorted( MLU::MakeFilename(sSummaryName,MLU::sParams,Seed,TEXT_EXT) );
      WriteTabular( MLU::MakeFilename(sSummaryName,MLU::sParams+"_table",Seed,TEXT_EXT), bi );
      Fits.clear();
    }
  }
}

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  static const char DefaultErrDig[] = "2";
  static const char DefaultTableLen[] = "66";
  static const char DefaultIgnoreMomenta[]{ "" };
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"stat", CL::SwitchType::Single, MLU::sPValueH.c_str() },
      {"inc",  CL::SwitchType::Flag, nullptr},
      {"strict",  CL::SwitchType::Single, "-1"},
      {"maxE",  CL::SwitchType::Single, "10"},
      {"fast",  CL::SwitchType::Flag, nullptr},
      {"all",  CL::SwitchType::Flag, nullptr},
      {"errdig", CL::SwitchType::Single, DefaultErrDig},
      {"tabn", CL::SwitchType::Flag, nullptr},
      {"tablen", CL::SwitchType::Single, DefaultTableLen},
      {"pignore",CL::SwitchType::Single, DefaultIgnoreMomenta},
      {"params",CL::SwitchType::Single, ""},
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
    "-i       Input  filename prefix\n"
    "-o       Output filename prefix\n"
    "--stat   Statistic (default: " << MLU::sPValueH << ")\n"
    "--strict Mask. Only include models where all parameters are non-zero and unique\n"
    "         Strictness bits: Off=+/- 1 sigma; On=every replica\n"
    "         Bit 0: Difference from 0\n"
    "         Bit 1: Difference from other parameters in same series\n"
    "--maxE   Maximum energy (default 10 - decays so fast effectively undetermined)\n"
    "--errdig Number of significant figures in error (default " << DefaultErrDig << ")\n"
    "--tablen Number of rows per page in table (default " << DefaultTableLen << ")\n"
    "Flags:\n"
    "--inc    Sort increasing (default sort test stat decreasing)\n"
    "--fast   Skip first pass reading models to find common parameters\n"
    "--all    Save all parameters. Default: only save common parameters\n"
    "--tabn   Include n^2 in the parameters table\n"
    "--pignore List of momenta to ignore (default: '" << DefaultIgnoreMomenta << "')\n"
    "--params Subset of parameters to save\n"
    "--help   This message\n";
  }
  return iReturn;
}
