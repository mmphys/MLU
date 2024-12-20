/**
 
 Mike's lattice QCD utilities: Bootstrapper

 Source file: bootstrap.hpp
 
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

#ifndef bootstrap_hpp
#define bootstrap_hpp

#include <MLU/MLU.hpp>
using MomentumMap = MLU::MomentumMap;
using MomentumMapValue = typename MomentumMap::value_type;
using Algebra = MLU::Gamma::Algebra;

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <set>

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif

enum StudySubject{ Z2=1, GFWW, GFPW };
enum GroupMomenta{ None, Squared, Abs };
enum BinOrder{ Auto, Old, VeryOld };

struct OperatorAttributes
{
  MLU::Gamma::Algebra alg;
  MLU::RPS rps;
};

struct TrajFile
{
  bool bHasTimeslice = false;
  int  Timeslice;
  MomentumMap p;
  int Config;
  MLU::Momentum pFirstNonZero;
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
  void Reverse( MLU::CorrelatorFileC &File ) const;
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
  MLU::Gamma::Algebra Alg3pt;
  // Filenames with corresponding timeslice info
  std::map<std::string, TrajFile> FileInfo;
  TrajList(const std::string &Name_, const std::string &sShortPrefix_, const std::string &sShortSuffix_,
           const std::string &opSuffixSnk_, const std::string &opSuffixSrc_,
           bool b3pt_, bool brev_, MLU::Gamma::Algebra Alg3pt_ )
  : Name{Name_}, sShortPrefix{sShortPrefix_}, sShortSuffix{ sShortSuffix_}, OpSuffixSnk{opSuffixSnk_},
    OpSuffixSrc{opSuffixSrc_}, b3pt{b3pt_}, bRev{brev_}, Alg3pt{Alg3pt_} {}
  // 2pt Trajectory list
  TrajList(const std::string &Name_, const std::string &sShortPrefix_, const std::string &sShortSuffix_,
           const std::string &opSuffixSnk_, const std::string &opSuffixSrc_ )
  : TrajList( Name_, sShortPrefix_, sShortSuffix_, opSuffixSnk_, opSuffixSrc_,
              false, false, MLU::Gamma::Algebra::MinusSigmaZT ) {}
  bool OpSuffiiSame() const { return MLU::EqualIgnoreCase( OpSuffixSnk, OpSuffixSrc ); }
};

class BootstrapParams
{
public:
  const bool b2ptSymOp;
  const bool b2ptSortZeroMom;
  const bool bOverwrite;
  const bool bVerboseSummaries;
  const int TimesliceDetail;
  const int nSample;
  const int binSize;
  const bool binAuto;
  const BinOrder binOrder;
  const MLU::SeedType seed;
  const std::string outStem;
  const std::string Ensemble;

  using CorrFile = MLU::CorrelatorFileC;
  using vCorrFile = std::vector<CorrFile>;
  using Iter = typename vCorrFile::iterator;
protected:
  MLU::SeedType GetSeedType( const MLU::CommandLine &cl );
  static std::vector<MLU::ConfigCount> CountConfigs( const Iter &first, const Iter &last, bool bGammai );
  bool GatherInput( MLU::SampleC &out, const Iter &first, const Iter &last, const TrajList &Traj,
                    Algebra Snk, Algebra Src, bool bAlignTimeslices ) const;
  int PerformBootstrap( const Iter &first, const Iter &last, const TrajList &Traj, const std::string &Suffix,
                        bool bAlignTimeslices, bool bSaveBootstrap, bool bSaveSummaries,
                        const std::vector<Algebra> &SinkAlgebra, const std::vector<Algebra> &SourceAlgebra ) const;
public:
  BootstrapParams( const MLU::CommandLine &cl );
  int PerformBootstrap( vCorrFile &f, const TrajList &Traj, const std::vector<Algebra> &SinkAlgebra,
                        const std::vector<Algebra> &SourceAlgebra ) const;
  void Study1Bootstrap( StudySubject Study, const std::string &StudyPath, const MLU::Momentum &mom,
                        std::vector<MLU::Gamma::Algebra> Alg,
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
  std::vector<MLU::NegateStar> AlgCurrentLoadNeg; // Whether to negate currents on load
public:
  const std::vector<Algebra> AlgCurrent; // Current algebra we wish to save for 3pt. May include Gammai. Empty = wildcard
protected:
  std::unique_ptr<std::regex> SSRegEx;
  std::vector<Algebra> GetCurrentAlgebra( const MLU::CommandLine &cl );
  static GroupMomenta GetGroupP( const MLU::CommandLine &cl );
  int QuarkWeight( const char q ) const;
  bool NeedsTimeReverse( std::string &Contraction, MomentumMap &p, bool &bRev,
                         std::string &OpSuffixSnk, std::string &OpSuffixSrc ) const;
public:
  // Process list of files on the command-line, breaking them up into individual trajectories
  Manifest( const MLU::CommandLine &cl, const BootstrapParams &par );
  /*Manifest(const std::vector<std::string> &Files, const std::vector<std::string> &Ignore,
           bool bSwapQuarks, const GroupMomenta GroupP, const std::vector<std::string> &vIgnoreMomenta,
           std::regex *SSRegEx, bool bSwapSnkSrcRegEx = false);*/
  void MakeStudies( const std::vector<std::string> &Args );
  void BuildManifest( const std::vector<std::string> &Args, const std::vector<std::string> &Ignore );
  bool RunManifest();
};

#endif // bootstrap_hpp
