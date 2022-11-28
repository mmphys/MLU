/**
 
 Mike's lattice QCD utilities: Bootstrapper

 Source file: bootstrap.hpp
 
 Copyright (C) 2019-2021
 
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
**/


#ifndef bootstrap_hpp
#define bootstrap_hpp

#include <MLU/Common.hpp>
using MomentumMap = Common::FileNameAtt::MomentumMap;
using MomentumMapValue = typename MomentumMap::value_type;
using Algebra = Common::Gamma::Algebra;
//using namespace Common;

//#include "CommonGrid.hpp"
//#include "CommonLatAn.hpp"

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <set>
//#include <filesystem> // C++17
//#include <LatAnalyze/Core/OptParser.hpp>
//#include <LatAnalyze/Statistics/Dataset.hpp>
//#include <LatAnalyze/Io/Io.hpp>
//#include <LatAnalyze/Io/Hdf5File.hpp>

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif

enum StudySubject{ Z2=1, GFWW, GFPW };
enum GroupMomenta{ None, Squared, Abs };
enum BinOrder{ Auto, Old, VeryOld };

struct OperatorAttributes
{
  Common::Gamma::Algebra alg;
  Common::RPS rps;
};

struct TrajFile
{
  bool bHasTimeslice = false;
  int  Timeslice;
  MomentumMap p;
  int Config;
  Common::Momentum pFirstNonZero;
  bool bTimeRev;
  int DeltaT;
  TrajFile(bool bHasTimeslice_, int Timeslice_, const MomentumMap &p_, int Config_, bool bTimeRev_, int DeltaT_)
  :bHasTimeslice{bHasTimeslice_},Timeslice{Timeslice_},p{p_},Config{Config_},bTimeRev{bTimeRev_},DeltaT{DeltaT_}
  {
    for( const MomentumMapValue & mmv : p )
    {
      if( mmv.second )
      {
        pFirstNonZero = mmv.second;
        break;
      }
    }
  }
  void Reverse( Common::CorrelatorFileC &File ) const;
};

// This describes one contraction and each of its trajectory files
struct TrajList
{
  const std::string Name;                     // name of the contraction
  const std::string sShortPrefix;
  const std::string sShortSuffix;
  const std::string OpSuffixSnk;
  const std::string OpSuffixSrc;
  const bool b3pt;
  // These members only present if b3pt
  bool bRev; // If set, the gamma in filename is quark 1(rhs/reversed), otherwise quark 2(lhs)
  Common::Gamma::Algebra Alg3pt;
  // Filenames with corresponding timeslice info
  std::map<std::string, TrajFile> FileInfo;
  TrajList(const std::string &Name_, const std::string &sShortPrefix_, const std::string &sShortSuffix_,
           const std::string &opSuffixSnk_, const std::string &opSuffixSrc_,
           bool b3pt_, bool brev_, Common::Gamma::Algebra Alg3pt_ )
  : Name{Name_}, sShortPrefix{sShortPrefix_}, sShortSuffix{ sShortSuffix_}, OpSuffixSnk{opSuffixSnk_},
    OpSuffixSrc{opSuffixSrc_}, b3pt{b3pt_}, bRev{brev_}, Alg3pt{Alg3pt_} {}
  // 2pt Trajectory list
  TrajList(const std::string &Name_, const std::string &sShortPrefix_, const std::string &sShortSuffix_,
           const std::string &opSuffixSnk_, const std::string &opSuffixSrc_ )
  : TrajList( Name_, sShortPrefix_, sShortSuffix_, opSuffixSnk_, opSuffixSrc_,
              false, false, Common::Gamma::Algebra::MinusSigmaZT ) {}
  bool OpSuffiiSame() const { return Common::EqualIgnoreCase( OpSuffixSnk, OpSuffixSrc ); }
};

class BootstrapParams
{
public:
  const bool b2ptSymOp;
  const bool b2ptSortZeroMom;
  const bool bWarnIfExists;
  const bool bVerboseSummaries;
  const int TimesliceDetail;
  const int nSample;
  const int binSize;
  const bool binAuto;
  const BinOrder binOrder;
  const Common::SeedType seed;
  const std::string outStem;
  const std::string MachineName;

  using CorrFile = Common::CorrelatorFileC;
  using vCorrFile = std::vector<CorrFile>;
  using Iter = typename vCorrFile::iterator;
protected:
  Common::SampleC RandomSample;
  Common::SeedType GetSeedType( const Common::CommandLine &cl );
  static std::vector<Common::ConfigCount> CountConfigs( const Iter &first, const Iter &last, bool bGammai );
  bool GatherInput( Common::SampleC &out, const Iter &first, const Iter &last, const TrajList &Traj,
                    Algebra Snk, Algebra Src, bool bAlignTimeslices ) const;
  int PerformBootstrap( const Iter &first, const Iter &last, const TrajList &Traj, const std::string &Suffix,
                        bool bAlignTimeslices, bool bSaveBootstrap, bool bSaveSummaries,
                        const std::vector<Algebra> &SinkAlgebra, const std::vector<Algebra> &SourceAlgebra ) const;
public:
  BootstrapParams( const Common::CommandLine &cl, const std::string MachineNameActual );
  int PerformBootstrap( vCorrFile &f, const TrajList &Traj, const std::vector<Algebra> &SinkAlgebra,
                        const std::vector<Algebra> &SourceAlgebra ) const;
  void Study1Bootstrap( StudySubject Study, const std::string &StudyPath, const Common::Momentum &mom,
                        std::vector<Common::Gamma::Algebra> Alg,
                        const std::string &Heavy, const std::string &Light, bool Factorised) const;
};

// This is a list of all the contractions we've been asked to process
struct Manifest : public std::map<std::string, TrajList>
{
  const BootstrapParams &par;
  const bool bShowOnly;
  const bool bTimeRev3pt;
  const bool bRevRev3pt;
  //const bool bEnable2ptSort; // Enables 2pt prefix sorting.
  const GroupMomenta GroupP;
  const std::string InStem;
  const std::string DefaultGroup;
  const std::string DefaultDataSet;
  const std::vector<std::string> vIgnoreMomenta;
  const std::vector<std::string> vIgnoreRegEx;
  std::vector<Algebra> AlgSource;  // Also the sink algebra for two-point functions. Empty=wildcard, loaded from first file
protected:
  std::vector<Algebra> AlgCurrentLoad; // These are the current algebra to load. Gammai removed and Gamma X, Y & Z added
  std::vector<Common::NegateStar> AlgCurrentLoadNeg; // Whether to negate currents on load
public:
  const std::vector<Algebra> AlgCurrent; // Current algebra we wish to save for 3pt. May include Gammai. Empty = wildcard
protected:
  std::unique_ptr<std::regex> SSRegEx;
  std::vector<Algebra> GetCurrentAlgebra( const Common::CommandLine &cl );
  static GroupMomenta GetGroupP( const Common::CommandLine &cl );
  int QuarkWeight( const char q ) const;
  bool NeedsTimeReverse( std::string &Contraction, MomentumMap &p, bool &bRev,
                         std::string &OpSuffixSnk, std::string &OpSuffixSrc ) const;
public:
  // Process list of files on the command-line, breaking them up into individual trajectories
  Manifest( const Common::CommandLine &cl, const BootstrapParams &par );
  /*Manifest(const std::vector<std::string> &Files, const std::vector<std::string> &Ignore,
           bool bSwapQuarks, const GroupMomenta GroupP, const std::vector<std::string> &vIgnoreMomenta,
           std::regex *SSRegEx, bool bSwapSnkSrcRegEx = false);*/
  void MakeStudies( const std::vector<std::string> &Args );
  void BuildManifest( const std::vector<std::string> &Args, const std::vector<std::string> &Ignore );
  bool RunManifest();
};

#endif // bootstrap_hpp
