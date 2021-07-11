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

#include <MLU/Common.hpp>

#include <set>
#include <gsl/gsl_randist.h>

using scalar = double;
using Model = Common::Model<scalar>;

const std::string &Sep{ Common::Space };
const std::string &NL{ Common::NewLine };
const std::string sIndent{ "  " };
const std::string sError{ "Error: " };

struct FitTimes
{
  int ti;
  int tf;
  std::string Extra;
  unsigned int NumDataPoints;
};

std::ostream &operator<<( std::ostream &s, const FitTimes &ft )
{
  s << ft.ti << ':' << ft.tf << "(";
  if( !ft.Extra.empty() )
    s << ft.Extra << ",";
  s << ft.NumDataPoints << ")";
  return s;
}

bool operator<( const FitTimes &lhs, const FitTimes &rhs )
{
  if( lhs.ti != rhs.ti )
    return lhs.ti < rhs.ti;
  if( lhs.tf != rhs.tf )
    return lhs.tf < rhs.tf;
  return Common::CompareIgnoreCase( lhs.Extra, rhs.Extra ) < 0;
}

struct FitData
{
  scalar SortBy;
  scalar Stat;
  FitTimes ft;
  int dof;
  std::string Parameters;
  std::size_t idxModel;
  int Seq = 0;
  FitData() = default;
  FitData(scalar sortBy_,scalar stat_,FitTimes ft_,int dof_,const std::string &Parameters_,std::size_t idx_Model)
  : SortBy{sortBy_},Stat{stat_},ft{ft_},dof{dof_},Parameters{Parameters_},idxModel{idx_Model}{}
  void Write( std::ostream &s, unsigned int NumDataPoints ) const
  {
    std::string tfLabel{ std::to_string( ft.tf ) };
    if( !ft.Extra.empty() )
    {
      tfLabel.append( 1, '_' );
      tfLabel.append( ft.Extra );
    }
    s << Seq << Sep << ft.ti << Sep << ft.tf << Sep << tfLabel << Sep << NumDataPoints << Sep << dof
      << Sep << Parameters << NL;
  }
};

struct TestStatKey
{
  scalar TestStatistic;
  FitTimes ft;
};

bool operator<( const TestStatKey &lhs, const TestStatKey &rhs )
{
  if( lhs.TestStatistic != rhs.TestStatistic )
    return lhs.TestStatistic < rhs.TestStatistic;
  return lhs.ft < rhs.ft;
}

// When I first parse file names, I split them into separate lists for each base
struct FileInfo
{
  FitTimes ft;
  std::string FileName;
  FileInfo( const FitTimes ft_, const std::string &FileName_ ) : ft{ft_}, FileName{FileName_} {}
};

struct Summariser
{
  const std::string inBase;
  const std::string outBaseFileName;
  const scalar Threshold;
  using BaseList = std::map<std::string, std::vector<FileInfo>>;
  BaseList lBase;
protected:
  BaseList MakeBaseList( const Common::CommandLine &cl );
  void SaveCMat(const std::vector<Model> &Corr, const int NumBoot, const std::string &sFileName,
                const std::vector<int> &OpIndices, bool bInvertNeg, bool bSingleModel );
public:
  void Run();
  Summariser( const Common::CommandLine &cl );
};

void Summariser::SaveCMat( const std::vector<Model> &Corr, const int NumBoot, const std::string &sFileName,
                           const std::vector<int> &OpIndices, bool bInvertNeg, bool bSingleModel )
{
  using Matrix = Eigen::MatrixXd; // dynamic sized matrix of complex double
  // Make covariance
  const int idx{ Model::idxCentral };
  const int NtCorr{ static_cast<int>( OpIndices.size() ) };
  const int NumFiles{ static_cast<int>( Corr.size() ) };
  const int Extent{ NtCorr * NumFiles };
  std::vector<scalar> VarianceInv( Extent );
  Matrix Covar( Extent, Extent );
  for( int x = 0; x < Extent; x++ )
  {
    const int ix{ OpIndices[x / NumFiles] };
    const Model &CorrX{ Corr[x % NumFiles] };
    const scalar * CentralX = CorrX[idx];
    const scalar * ReplicaX = CorrX[ 0 ];
    const int NtX { CorrX.Nt() };
    for( int y = 0; y <= x; y++ )
    {
      const int iy{ OpIndices[y / NumFiles] };
      const Model &CorrY{ Corr[y % NumFiles] };
      const scalar * CentralY = CorrY[idx];
      const scalar * ReplicaY = CorrY[ 0 ];
      const int NtY { CorrY.Nt() };
      double z = 0;
      const scalar * DataX{ ReplicaX };
      const scalar * DataY{ ReplicaY };
      for( int i = 0; i < NumBoot; i++ )
      {
        scalar zThis = ( DataX[ix] - CentralX[ix] ) * ( DataY[iy] - CentralY[iy] );
        if( bInvertNeg && (   ( CentralX[ix] <  0 && CentralY[iy] >= 0 )
                           || ( CentralX[ix] >= 0 && CentralY[iy] <  0 ) ) )
          zThis = -zThis;
        z += zThis;
        DataX += NtX;
        DataY += NtY;
      }
      z /= NumBoot;
      Covar( x, y ) = z;
      if( x != y )
        Covar( y, x ) = z;
      if( y == x )
        VarianceInv[x] = 1. / sqrt( z );
    }
  }
  // Turn covariance matrix into correlation matrix
  for( int x = 0; x < Extent; x++ )
    for( int y = 0; y <= x; y++ )
    {
      if( y == x )
        Covar( x, y ) = 1.;
      else
      {
        double z = Covar( x, y ) * VarianceInv[x] * VarianceInv[y];
        Covar( x, y ) = z;
        Covar( y, x ) = z;
      }
    }
  if( !Common::IsFinite( Covar ) )
    throw std::runtime_error( "Covariance matrix isn't finite" );
  // Now write file
  std::ofstream s{ sFileName };
  s << "# Correlation matrix\n# Files: " << NumFiles << Common::NewLine;
  for( int i = 0; i < NumFiles; i++ )
    s << "# File" << ( i + 1 ) << ": " << Corr[i].Name_.Filename << Common::NewLine;
  s << "# Operators: " << NtCorr;
  for( int o = 0; o < NtCorr; o++ )
    s << " " << Corr[0].GetColumnNames()[OpIndices[o]];
  s << "\n# gnuplot: plot '" << sFileName
    << "' matrix columnheaders rowheaders with image pixels" << Common::NewLine << Extent;
  std::vector<std::string> sOpNames;
  for( int o = 0; o < NtCorr; o++ )
  {
    const int iColumn{ OpIndices[o] };
    for( int f = 0; f < NumFiles; f++ )
    {
      const bool bThisNeg{ Corr[f][idx][iColumn] < 0 };
      std::string sOp{ "\"" + Corr[f].GetColumnNames()[iColumn] };
      if( bThisNeg )
        sOp.append( 1, '-' );
      if( !bSingleModel )
      {
        if( !bThisNeg )
          sOp.append( 1, ' ' );
        sOp.append( std::to_string( Corr[f].ti ) );
        sOp.append( 1, ',' );
        sOp.append( std::to_string( Corr[f].tf ) );
      }
      sOp.append( 1, '"' );
      s << Common::Space << sOp;
      sOpNames.push_back( sOp );
    }
  }
  s << Common::NewLine;
  for( int o = 0; o < NtCorr; o++ )
    for( int f = 0; f < NumFiles; f++ )
    {
      const int i{ o * NumFiles + f };
      s << sOpNames[i];
      for( int j = 0; j < Extent; j++ )
        s << Common::Space << Covar( i, j );
      s << Common::NewLine;
    }
}

// Parse the command line, splitting into separate lists for each base
Summariser::Summariser( const Common::CommandLine &cl )
: inBase{ Common::AppendSlash( cl.SwitchValue<std::string>("i") ) },
  outBaseFileName{ Common::AppendSlash( cl.SwitchValue<std::string>("o") ) },
  Threshold{ cl.SwitchValue<scalar>("p") }
{
  for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
  {
    Common::FileNameAtt n{ sFileName, };
    n.ParseExtra( 2 );
    std::size_t NumExtra{ n.Extra.size() };
    if( !Common::FileExists( sFileName ) )
      throw std::runtime_error( sFileName + " doesn't exist" );
    if( !NumExtra )
      throw std::runtime_error( "No fit type in " + sFileName );
    std::string sFitType{ n.Extra[NumExtra - 1] };
    const std::size_t pos_ti = sFitType.find_first_of( '_' );
    bool bOK{ false };
    FitTimes ft;
    if( pos_ti != std::string::npos )
    {
      const std::size_t pos_tf{ sFitType.find_first_of( '_', pos_ti + 1 ) };
      if( pos_tf != std::string::npos )
      {
        const std::size_t pos_extra{ sFitType.find_first_of( '_', pos_tf + 1 ) };
        if( pos_extra != std::string::npos )
          ft.Extra = sFitType.substr( pos_extra + 1 );
        for( std::size_t i = pos_ti + 1; i < sFitType.length(); ++i )
          if( sFitType[i] == '_' )
            sFitType[i] = ' ';
        std::vector<int> Ranges;
        try{
          Ranges = Common::ArrayFromString<int>( sFitType.substr( pos_ti + 1 ) );
        } catch(...) {}
        if( !Ranges.empty() && !( Ranges.size() % 2 ) ) // Must contain an even number
        {
          ft.ti = Ranges[0];
          ft.tf = Ranges[1];
          ft.NumDataPoints = 0;
          for( std::size_t i = 0; i < Ranges.size(); i += 2 )
            ft.NumDataPoints += Ranges[i + 1] - Ranges[i] + 1;
          bOK = true;
        }
      }
    }
    if( !bOK )
      throw std::runtime_error( "Error extracting fit ranges from " + sFileName );
    sFitType.resize( pos_ti );
    // The master key for this record is the base filename of the output
    std::string sOutFile{ n.Base };
    sOutFile.append( 1, '.' );
    sOutFile.append( sFitType );
    if( NumExtra == 2 )
    {
      sOutFile.append( 1, '.' );
      sOutFile.append( n.Extra[0] );
    }
    lBase[sOutFile].emplace_back( FileInfo( ft, sFileName ) );
  }
}

// Process every base name and it's list of models separately
void Summariser::Run()
{
  for( BaseList::iterator it = lBase.begin(); it != lBase.end(); ++it )
  {
    const std::string &sOutFile{ it->first };
    std::vector<FileInfo> FileList{ it->second };
    const std::string &sSummaryName{ outBaseFileName + sOutFile };
    std::cout << sSummaryName << NL;
    std::vector<Model> Models;
    Models.reserve( FileList.size() );
    // Read every model, saving all the data I'll need to reconstruct the summary
    int NumSamples{ 0 };
    std::string ColumnNames;
    std::string Comments;
    Common::SeedType Seed{ 0 };
    using FitMap = std::map<FitTimes, FitData>;
    FitMap Fits;
    for( FileInfo &ThisFile : FileList )
    {
      try
      {
        // Check the model characteristics
        Model m;
        m.Read( ThisFile.FileName, "" );
        if( m.dof < 1 )
          throw std::runtime_error("dof=" + std::to_string( m.dof ) + " i.e. extrapolation, not fit");
        std::ostringstream ss;
        m.WriteColumnNames( ss );
        ss << Sep;
        using vEr = Common::ValWithEr<scalar>;
        vEr::Header( "pvalue", ss );
        ss << Sep;
        vEr::Header( "pvalueH", ss );
        if( Models.empty() )
        {
          ColumnNames = ss.str();
          ss.str("");
          m.SummaryComments( ss );
          ss << "# column(1) is the row order when sorted by ChiSqPerDof\n";
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
        m.WriteSummaryData( ss );
        const unsigned int NumDataPoints{ ThisFile.ft.Extra.empty() ? (m.tf - m.ti + 1) * m.NumFiles
                                        : ThisFile.ft.NumDataPoints};
        ThisFile.ft.NumDataPoints = NumDataPoints;
        const vEr & ChisqDof{ m.getSummaryData()[m.Nt() - 1] };
        const vEr ChiSq{ ChisqDof * m.dof }; // Really wish I'd saved the test statistic, not reduced test statistic
        // Chi squared statistic and Q-value (probability of a worse statistic)
        // There is a deliberate inversion here: high statistic => low q-value
        const vEr qValueChiSq{ ChiSq.qValueChiSq( m.dof ) };
        ss << Sep << qValueChiSq;
        // Hotelling t statistic and Q-value
        const vEr Hotelling{ ChiSq.qValueHotelling( m.dof, NumDataPoints ) };
        ss << Sep << Hotelling;
        // save the model in my list
        Fits.emplace( std::make_pair( ThisFile.ft,
                      FitData(1-Hotelling.Central,Hotelling.Central,ThisFile.ft,m.dof,ss.str(),Models.size())));
        Models.emplace_back( std::move( m ) );
      }
      catch( const std::exception &e )
      {
        //std::cout << sFileName << Sep << sError << e.what() << NL;
        std::cout << sIndent << sError << e.what() << NL;
      }
    }
    // Update the fits with sorted sequence number for test statistic and save sorted file
    static const std::string InitialColumnNames{ " ti tf tfLabel NumDataPoints dof " };
    {
      // Sort the fits by chi-squared
      std::set<TestStatKey> SortSet;
      for( auto dt = Fits.begin(); dt != Fits.end(); ++dt )
      {
        const FitData &d{ dt->second };
        SortSet.emplace( TestStatKey{ d.SortBy, dt->first } );
      }
      std::string sFileName=Common::MakeFilename(sSummaryName,Common::sParams+"_sort",Seed,TEXT_EXT);
      std::ofstream s( sFileName );
      Common::SummaryHeader<scalar>( s, sFileName );
      s << Comments;
      s << "Seq" << InitialColumnNames << ColumnNames << NL;
      int SortSeq{ 0 };
      for( const TestStatKey &z : SortSet )
      {
        FitData & d{ Fits[z.ft] };
        d.Seq = SortSeq++;
        d.Write( s, z.ft.NumDataPoints );
      }
    }
    // Now make the output file with blocks for each ti, sorted by tf
    const bool bSingleModel{ Models.size() == 1 };
    std::vector<Model> CorrModels;
    CorrModels.reserve( Models.size() );
    {
      std::ofstream s;
      const std::string sFileName{Common::MakeFilename(sSummaryName,Common::sParams,Seed,TEXT_EXT)};
      std::size_t t_last = std::numeric_limits<int>::min();
      for( auto dt = Fits.begin(); dt != Fits.end(); ++dt )
      {
        const FitData &d{ dt->second };
        // Write header before each block
        if( t_last != d.ft.ti )
        {
          t_last = d.ft.ti;
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
          s << "# [ti=" << d.ft.ti << "]\n";
          // Column names, with the series value embedded in the column header (best I can do atm)
          s << "ti=" << d.ft.ti << InitialColumnNames << ColumnNames << NL;
        }
        // Now write this row
        d.Write( s, d.ft.NumDataPoints );
        // If the statistic is above threshold, include it in correlation matrix
        if( d.Stat >= Threshold )
          CorrModels.emplace_back( std::move( Models[d.idxModel] ) );
      }
    }
    // Now plot the correlation matrix for selected models
    Models.clear();
    try
    {
      if( !CorrModels.empty() )
      {
        const Model &m{ CorrModels[0] };
        std::string BaseName;
        if( bSingleModel )
        {
          BaseName = m.Name_.Base;
          for( std::size_t i = m.Name_.Extra.size(); i > 0; )
          {
            BaseName.append( 1, '.' );
            BaseName.append( m.Name_.Extra[--i] );
          }
        }
        else
          BaseName = sSummaryName;
        for( int i = 0; i < 2; i++ )
        {
          const int NumOperators{ static_cast<int>( m.OpNames.size() ) };
          std::vector<int> OpIndices{};
          OpIndices.clear();
          OpIndices.reserve( m.NumExponents * NumOperators );
          for( int o = 0; o < NumOperators; o++ )
            for( int e = 0; e < m.NumExponents; e++ )
              OpIndices.push_back( m.GetColumnIndex( m.OpNames[o], e ) );
          std::string sFileName=Common::MakeFilename(BaseName+".small",Common::sCormat,Seed,TEXT_EXT);
          std::cout << "Making " << sFileName << NL;
          SaveCMat( CorrModels, NumSamples, sFileName, OpIndices, i, bSingleModel );
          if( m.NumExponents > 1 )
          {
            OpIndices.clear();
            for( int e = 0; e < m.NumExponents; e++ )
            {
              OpIndices.push_back( e );
              for( int o = 0; o < NumOperators; o++ )
                OpIndices.push_back( m.GetColumnIndex( m.OpNames[o], e ) );
            }
            sFileName=Common::MakeFilename(BaseName, Common::sCormat, Seed, TEXT_EXT );
            std::cout << "Making " << sFileName << NL;
            SaveCMat( CorrModels, NumSamples, sFileName, OpIndices, i, bSingleModel );
          }
          BaseName.append( ".cpos" );
        }
      }
    }
    catch( const std::exception &e )
    {
      std::cout << sIndent << sError << e.what() << NL;
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
      {"p", CL::SwitchType::Single, "0.05" },
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
    "-p     p-value threshold for covariance matrix\n"
    "Flags:\n"
    "--help This message\n";
  }
  return iReturn;
}
