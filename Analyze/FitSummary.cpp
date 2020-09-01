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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

using scalar = double;
using Model = Common::Model<scalar>;

using FitTimes = std::pair<int, int>;
struct FitData
{
  scalar SortBy;
  scalar Stat;
  int ti;
  int tf;
  int dof;
  std::string Parameters;
  std::size_t idxModel;
  int Seq = 0;
  FitData() = default;
  FitData(scalar sortBy_,scalar stat_,int ti_,int tf_,int dof_,const std::string &Parameters_,std::size_t idx_Model)
  : SortBy{sortBy_},Stat{stat_},ti{ti_},tf{tf_},dof{dof_},Parameters{Parameters_},idxModel{idx_Model}{}
};

/*struct Detail
{
  std::map<FitTimes, FitData> Fits;
  int Nt;
  std::string ColumnNames;
  std::string Comments;
  Common::SeedType Seed;
  std::string SeedMachine;
};

using Summary = std::map<std::string, Detail>;*/
using FitMap = std::map<FitTimes, FitData>;

using chisqdof_ti_tf = std::tuple<scalar, int, int>;
using SortByDof = std::map<chisqdof_ti_tf, std::string>;

// When I first parse file names, I split them into separate lists for each base
using BaseList = std::map<std::string, std::vector<std::string>>;
using Matrix = Eigen::MatrixXd; // dynamic sized matrix of complex double

void SaveCMat(const std::vector<Model> &Corr, const int NumBoot, const std::string &sFileName,
              const std::vector<int> &OpIndices, bool bInvertNeg, bool bSingleModel )
{
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

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  const std::string &Sep{ Common::Space };
  const std::string &NL{ Common::NewLine };
  static const std::string sIndent{ "  " };
  static const std::string sError{ "Error: " };
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
    BaseList lBase;
    const std::string outBaseFileName{ Common::AppendSlash( cl.SwitchValue<std::string>("o") ) };
    const scalar Threshold{ cl.SwitchValue<scalar>("p") };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      const std::string inBase{ Common::AppendSlash( cl.SwitchValue<std::string>("i") ) };
      bShowUsage = false;
      // Parse the command line, splitting into separate lists for each base
      for( const std::string &sFileName : Common::glob( cl.Args.begin(), cl.Args.end(), inBase.c_str()))
      {
        Common::FileNameAtt n{ sFileName, };
        n.ParseExtra( 2 );
        std::size_t NumExtra{ n.Extra.size() };
        std::string sFitType;
        if( !Common::FileExists( sFileName ) )
          std::cout << sFileName << Sep << sError << "doesn't exist" << NL;
        else if( !NumExtra )
          std::cout << sFileName << Sep << sError << "doesn't include fit type in name" << NL;
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
        lBase[sOutFile].emplace_back( std::move( sFileName ) );
      }
    }
    // Process every base name and it's list of models separately
    for( BaseList::iterator it = lBase.begin(); it != lBase.end(); ++it )
    {
      const std::string &sOutFile{ it->first };
      const std::vector<std::string> ModelFileName{ it->second };
      const std::string &sSummaryName{ outBaseFileName + sOutFile };
      std::cout << sSummaryName << NL;
      std::vector<Model> Models;
      Models.reserve( ModelFileName.size() );
      // Read every model, saving all the data I'll need to reconstruct the summary
      int NumSamples{ 0 };
      std::string ColumnNames;
      std::string Comments;
      Common::SeedType Seed{ 0 };
      FitMap Fits;
      for( const std::string &sFileName : ModelFileName )
      {
        try
        {
          // Check the model characteristics
          Model m;
          m.Read( sFileName, "" );
          if( m.dof < 1 )
            throw std::runtime_error("dof=" + std::to_string( m.dof ) + " i.e. extrapolation, not fit");
          std::ostringstream ss;
          m.WriteColumnNames( ss );
          ss << " pvalue pvalue_low pvalue_high pvalue_check pvalueH pvalueH_low pvalueH_high pvalueH_check";
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
          using vEr = Common::ValWithEr<scalar>;
          const vEr * const pChisqDof{ m.getSummaryData() + m.Nt()-1 };
          // Chi squared statistic and Q-value (probability of a worse statistic)
          // There is a deliberate inversion here: high statistic => low q-value
          scalar QValue = gsl_cdf_chisq_Q( pChisqDof->Central * m.dof, m.dof );
          ss << Sep << QValue
             << Sep << gsl_cdf_chisq_Q( pChisqDof->High * m.dof, m.dof )
             << Sep << gsl_cdf_chisq_Q( pChisqDof->Low * m.dof, m.dof )
             << Sep << pChisqDof->Check;
          // Hotelling t statistic and Q-value
          const int Extent{(m.tf-m.ti+1)*m.NumFiles}; // Extent of covariance matrix, i.e. # data points in fit
          const int p{ Extent - m.dof }; // Number of params in fit
          scalar tFactor{ static_cast<scalar>( m.dof * m.dof ) / ( p * ( Extent - 1 ) ) };
          const vEr Hotelling{pChisqDof->Central * tFactor, pChisqDof->High * tFactor, pChisqDof->Low * tFactor};
          /*std::cout << "Debug:\n"
          << "\n\tm.dof = " << m.dof << "\tExtent = " << Extent << Common::NewLine
          << "\tstatic_cast<scalar>( m.dof * m.dof ) = " << static_cast<scalar>( m.dof * m.dof )
          << "\n\tp * ( Extent - 1 ) = " << ( p * ( Extent - 1 ) )
          << "\n\ttFactor = " << tFactor
          << "\n\t*pChisqDof = " << *pChisqDof
          << "\n\tHotelling = " << Hotelling
          << Common::NewLine;*/
          ss << Sep << gsl_cdf_fdist_Q( Hotelling.Central, p, m.dof )
             << Sep << gsl_cdf_fdist_Q( Hotelling.Low, p, m.dof )
             << Sep << gsl_cdf_fdist_Q( Hotelling.High, p, m.dof )
             << Sep << pChisqDof->Check;
          // save the model in my list
          Fits.emplace( std::make_pair(std::make_pair( m.ti, m.tf ),
                                  FitData(1-QValue,QValue,m.ti,m.tf,m.dof,ss.str(),Models.size())));
          Models.emplace_back( std::move( m ) );
        }
        catch( const std::exception &e )
        {
          //std::cout << sFileName << Sep << sError << e.what() << NL;
          std::cout << sIndent << sError << e.what() << NL;
        }
      }
      // Update the fits with sorted sequence number for chi-squared and save sorted file
      {
        // Sort the fits by chi-squared
        std::set<chisqdof_ti_tf> vChiSort;
        for( auto dt = Fits.begin(); dt != Fits.end(); ++dt )
        {
          const FitData &d{ dt->second };
          vChiSort.emplace( d.SortBy, d.ti, d.tf );
        }
        std::string sFileName=Common::MakeFilename(sSummaryName,Common::sParams+"_sort",Seed,TEXT_EXT);
        std::ofstream s( sFileName );
        Common::SummaryHeader<scalar>( s, sFileName );
        s << Comments;
        s << "Seq ti tf dof " << ColumnNames << NL;
        int ChiSeq{ 0 };
        for( const chisqdof_ti_tf &z : vChiSort )
        {
          FitData & d{ Fits[std::make_pair( std::get<1>( z ), std::get<2>( z ) )] };
          d.Seq = ChiSeq++;
          s << d.Seq << Sep << d.ti << Sep << d.tf << Sep << d.dof << Sep << d.Parameters << NL;
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
          if( t_last != d.ti )
          {
            t_last = d.ti;
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
            s << "# [ti=" << d.ti << "]\n";
            // Column names, with the series value embedded in the column header (best I can do atm)
            s << "ti=" << d.ti << " ti tf dof " << ColumnNames << NL;
          }
          // Now write this row
          s << d.Seq << Sep << d.ti << Sep << d.tf << Sep << d.dof << Sep << d.Parameters << NL;
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
