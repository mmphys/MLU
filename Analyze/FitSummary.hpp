/*************************************************************************************
 
 Summarise fits
 
 Source file: FitSummary.hpp
 
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

#ifndef Fit_Summary_hpp
#define Fit_Summary_hpp

#include <MLU/Common.hpp>

#include <set>
#include <gsl/gsl_randist.h>

using scalar = double;
using Model = Common::Model<scalar>;

extern bool bReverseSort; // TODO: shouldn't be global

struct FitTimes
{
  std::vector<int> time;
  int ti() const { return time[0]; };
  int tf() const { return time[1]; };
  std::string tfLabel() const;
  unsigned int GuessOldNumDataPoints( int NumFiles ) const;
  bool ParseTimes( std::string Times );
  bool ParseTypeTimes( std::string &TypeTimes );
};

struct FitData
{
  scalar Stat;
  FitTimes ft;
  int NumDataPoints;
  int dof;
  int SampleSize;
  std::string Parameters;
  std::size_t idxModel;
  int Seq = 0;
  FitData() = default;
  FitData( scalar stat_, FitTimes ft_, int NumDataPoints_, int dof_, int SampleSize_,
           const std::string &Parameters_,std::size_t idx_Model)
  : Stat{stat_}, ft{ft_}, NumDataPoints{NumDataPoints_}, dof{dof_}, SampleSize{SampleSize_},
    Parameters{Parameters_}, idxModel{idx_Model} {}
};

struct TestStatKey
{
  scalar TestStatistic;
  FitTimes ft;
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

struct Summariser
{
  const std::string inBase;
  const std::string outBaseFileName;
  const std::string StatisticName;
  const int Strictness;
  const scalar MonotonicUpperLimit;
  using FitMap = std::map<FitTimes, FitData>;
  using BaseList = std::map<std::string, std::vector<FileInfo>>;
  BaseList lBase;
  void Run();
  Summariser( const Common::CommandLine &cl );
};

#endif //Fit_Summary_hpp
