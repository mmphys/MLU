/**
 
 Summarise fits
 
 Source file: FitSummary.hpp
 
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

#ifndef Fit_Summary_hpp
#define Fit_Summary_hpp

#include <MLU/MLU.hpp>

#include <set>
#include <gsl/gsl_randist.h>

using scalar = double;
using Model = MLU::Model<scalar>;
using veScalar = MLU::ValWithEr<typename Model::scalar_type>;

struct FitTimes
{
  std::vector<int> time;
  int NumFitTimes() const { return static_cast<int>( time.size() / 2 ); }
  int ti( std::size_t idx = 0 ) const { return time[2 * idx]; };
  int tf( std::size_t idx = 0 ) const { return time[2 * idx + 1]; };
  std::string tiLabel() const;
  std::string tfLabel() const;
  unsigned int GuessOldNumDataPoints( int NumFiles ) const;
  bool ParseTimes( std::string Times );
  bool ParseTypeTimes( std::string &TypeTimes );
  bool operator<( const FitTimes &rhs ) const;
protected:
  void tiLabelSuffix( std::string &s ) const;
};

struct FitData
{
  scalar Stat;
  int NumDataPoints;
  int dof;
  int SampleSize;
  int CovarSampleSize;
  std::size_t idxLoadOrder;
  std::vector<veScalar> Parameters;
  int Seq = 0;
  FitData() = default;
  FitData( scalar stat_, int NumDataPoints_, int dof_, int SampleSize_, int CovarSampleSize_,
           std::size_t idxLoadOrder_, std::vector<veScalar> &&Parameters_ )
  : Stat{stat_}, NumDataPoints{NumDataPoints_}, dof{dof_}, SampleSize{SampleSize_},
    CovarSampleSize{CovarSampleSize_}, idxLoadOrder{idxLoadOrder_}, Parameters{Parameters_} {}
  void WriteLongFormat( std::ostream &os, const FitTimes &ft ) const;
  void WriteTableFormat( std::ostream &os, const FitTimes &ft, std::size_t StatIndex,
                         unsigned char ErrorDigits ) const;
};

struct TestStatKey
{
  static bool bReverseSort;
  scalar TestStatistic;
  FitTimes ft;
  bool operator<( const TestStatKey &rhs ) const;
};

// When I first parse file names, I split them into separate lists for each base
struct FileInfo
{
  FitTimes ft;
  std::string FileName;
  bool bIsFit;
  FileInfo( const FitTimes ft_, const std::string &FileName_, bool bIsFit_ )
  : ft{ft_}, FileName{FileName_}, bIsFit{bIsFit_} {}
};

struct BaseInfo
{
  const MLU::Momentum p;
  std::vector<FileInfo> vFI;
  BaseInfo( const MLU::Momentum &p_ ) : p{p_} {}
};

struct Summariser
{
  const std::string inBase;
  const std::string outBaseFileName;
  const std::string StatisticName;
  const int Strictness;
  const scalar MonotonicUpperLimit;
  const bool bAll;
  const bool bFast;
  const bool bTableN;
  const unsigned char ErrorDigits;
  const std::size_t TableRowsPerPage;
  const std::vector<std::string> vIgnoreMomenta;
  const std::vector<std::string> ParamSubset;
  using FitMap = std::map<FitTimes, FitData>;
  using BaseList = std::map<std::string, BaseInfo>;
  BaseList lBase;
  void Run();
  Summariser( const MLU::CommandLine &cl );
protected:
  using FileInfoIterator = std::vector<FileInfo>::iterator;
  bool bDoPassOne;
  std::array<Model, 2> Models;
  std::vector<std::string> FileNameOps;
  std::size_t MaxFitTimes; // Maximum number of TI-TF pairs
  MLU::Params Params;
  std::size_t NumParams() const { return Params.NumScalars( MLU::Param::Type::All ); }
  MLU::UniqueNameSet StatColumnNames;
  std::size_t NumStats() const { return StatColumnNames.size(); }
  FitMap Fits;
  // Saved about first model by BuildFitMap()
  std::string Comments;
  MLU::SeedType Seed{ 0 };
  int NumSamples{ 0 };
  bool ReadModel( Model &m, FileInfoIterator &it, std::vector<FileInfo> &Files,
                  std::vector<std::string> &FileNameOps, bool bShow );
  bool GetCommonParameters( std::vector<FileInfo> &Files, bool bMaximum );
  void BuildFitMap( std::vector<FileInfo> &Files );
  void SummaryColumnNames( std::ostream &os, std::size_t NumFitTimes,
              const MLU::Params &ParamNames, const MLU::UniqueNameSet &StatNames ) const;
  void WriteSorted( const std::string &sFileName, const std::set<TestStatKey> &SortSet, bool SetSeq );
  void WriteUnsorted( const std::string &sFileName ) const;
  void WriteTableHeader( std::ofstream &os ) const;
  void WriteTableFooter( std::ofstream &os ) const;
  void WriteTabular( const std::string &sFileName, const BaseInfo &bi ) const;
};

#endif //Fit_Summary_hpp
