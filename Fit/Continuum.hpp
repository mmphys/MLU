/*************************************************************************************
 
 Chiral continuum fit
 
 Source file: Continuum.hpp
 
 Copyright (C) 2023
 
 Author: Michael Marshall <Mike@lqcd.me>
 
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

#ifndef Continuum_hpp
#define Continuum_hpp

#include "MultiFit.hpp"
#include "Fitter.hpp"
#include "Model.hpp"
#include "ModelContinuum.hpp"

// Stats on the ensembles we loaded from files
struct EnsembleStat
{
  int Num;
  Common::SeedType Seed;
};
using EnsembleStatMap = std::map<std::string, EnsembleStat, Common::LessCaseInsensitive>;
using ESPair = std::pair<typename EnsembleStatMap::iterator, bool>;

struct EnsembleInfo
{
  unsigned int aInv_L;
  unsigned int aInv_T;
};
using EnsembleMapT = std::map<std::string, EnsembleInfo>;

struct ContinuumFit;

struct CreateParams : public Model::CreateParams
{
  const ContinuumFit &Parent;
  CreateParams( const std::vector<std::string> &OpNames, const Common::CommandLine &cl,
                const ContinuumFit &Parent );
};

struct ContinuumFit : public FitController
{
  static const std::string &FieldQSq;
  static const std::string &FieldEL;
protected:
  static constexpr scalar InvGeV{ 1e-9 };
  static constexpr scalar InvGeVSq{ InvGeV * InvGeV };
  static constexpr scalar Lambda{ 1e9 };
  static constexpr scalar InvLambda{ 1. / Lambda };
  static constexpr int NumTicks{ 200 };
  static constexpr int NumFF{ 2 };
  static constexpr int NumConst{ 5 };
  static const std::string sPDG;
  Common::CommandLine &cl;
public:
  EnsembleMapT EnsembleMap;
protected:

  struct EnsembleFF
  {
    std::string Ensemble;
    Common::FormFactor ff;
    EnsembleFF( const std::string &ensemble, Common::FormFactor FF ) : Ensemble{ensemble}, ff{FF} {}
  };
  struct EnsembleFFLess
  {
    bool operator()( const EnsembleFF &lhs, const EnsembleFF &rhs ) const
    {
      int i = Common::CompareIgnoreCase( lhs.Ensemble, rhs.Ensemble );
      if( i )
        return i < 0;
      return lhs.ff < rhs.ff;
    }
  };
  using EnsembleFFSetT = std::set<EnsembleFF, EnsembleFFLess>;

public:
  const int NumSamples;
  const bool doCorr;
  const bool CovarBlock;
  const Common::FormFactor ffDefault;
  const std::array<std::array<bool, NumConst>, NumFF> cEnabled;
  const std::string inBase;
  const std::string outBaseFileName;
  DataSet ds;
  std::vector<std::string> OpName;
  std::vector<Model::Args> ModelArgs;
  std::unique_ptr<Fitter> f;
  explicit ContinuumFit( Common::CommandLine &cl );
  virtual ~ContinuumFit() {}
  void ParamsAdjust( Common::Params &mp, const Fitter &f ) const override;
  void SaveParameters( Common::Params &mp, const Fitter &f ) override;
  void SetReplica( Vector &ModelParams ) const override;
  void ComputeDerived( Vector &ModelParams ) const override;
  int Run();
  inline bool CEnabled( int idxFF, int C ) const { return cEnabled[ idxFF ][C]; }
  inline bool CEnabled( Common::FormFactor ff, int C ) const { return CEnabled( ffIndex( ff ), C ); }
  inline bool CNeeded( int idxFF, int C ) const { return C == 0 || CEnabled( idxFF, C ); }
  inline bool CNeeded( Common::FormFactor ff, int C ) const { return CNeeded( ffIndex( ff ), C ); }
protected:
  EnsembleStatMap EnsembleStats;
  EnsembleFFSetT  EnsembleFFs;
  // Which form factor am I fitting
  static constexpr unsigned int uiFF0 = 1;
  static constexpr unsigned int uiFFPlus = 2;
  unsigned int uiFF;
  inline unsigned int ffMaskFromIndex( int idx ) const { return idx ? uiFFPlus : uiFF0; }
  inline unsigned int ffMask( Common::FormFactor ff ) const
  {
    if( ff == Common::FormFactor::f0 )
      return uiFF0;
    if( ff == Common::FormFactor::fplus )
      return uiFFPlus;
    std::ostringstream os;
    os << "ContinuumFit::ffMask() " << ff << " invalid";
    throw std::runtime_error( os.str().c_str() );
  }
  // Parameters from the model just created
  static constexpr int idxFF0 = 0;
  static constexpr int idxFFPlus = 1;
  inline Common::FormFactor ffIndexReverse( int idx ) const
  { return idx ? Common::FormFactor::fplus : Common::FormFactor::f0; }
public:
  inline int ffIndex( Common::FormFactor ff ) const
  {
    if( ff == Common::FormFactor::f0 )
      return idxFF0;
    if( ff == Common::FormFactor::fplus )
      return idxFFPlus;
    std::ostringstream os;
    os << "ContinuumFit::ffIndex() " << ff << " invalid";
    throw std::runtime_error( os.str().c_str() );
  }
protected:
  static constexpr std::size_t idxCUnused{ std::numeric_limits<std::size_t>::max() };
  std::array<std::array<std::size_t, NumConst>, NumFF> idxC;
  std::array<std::size_t, NumFF> idxDelta;
  std::array<std::size_t, NumFF> idxPDGDStar;
  std::size_t idxPDGH;
  std::size_t idxPDGL;
  std::size_t idxmPDGPi;
  std::vector<std::size_t> aInv, mPi, FVSim, FVPhys;
  inline static scalar EOfQSq( scalar PDGH, scalar PDGL, scalar qSq )
  {
    const scalar E{ ( PDGH * PDGH + PDGL * PDGL - qSq ) / ( 2 * PDGH ) };
    return E;
  }
  inline scalar EOfQSq( int rep, scalar qSq ) const
  {
    return EOfQSq( f->OutputModel(rep,idxPDGH), f->OutputModel(rep,idxPDGL), qSq );
  }
  std::string sOpNameConcat; // Sorted, concatenated list of operators in the fit for filenames
  std::array<std::string, 2> Meson;
  static const std::string &GetPoleMassName( Common::FormFactor ff, const Common::FileNameAtt &fna );
  static Common::FormFactor ValidateFF( Common::FormFactor ff );
  void AddEnsemble( const std::string &Ensemble, Common::FormFactor thisFF );
  std::array<std::array<bool, NumConst>, NumFF> GetEnabled( std::string sOptions );
  void LoadModels();
  void GetEnsembleStats();
  void SetEnsembleStats();
  void SetEnsembleStatSeed();
  void SortModels();
  void LoadExtra();
  void MakeOutputFilename();
  std::string GetOutputFilename( unsigned int uiFF );
  void ShowCovar();
  void MakeCovarBlock();
  void GetMinMax( Common::FormFactor ff, scalar &Min, scalar &Max,
                  ModelParam ModelContinuum::* mp, const std::string &Field ) const;
  std::ofstream WriteHeader( const std::string &sPrefix, const std::string &FileType ) const;
  void WriteFieldName( std::ofstream &os, const std::string &FieldName ) const;
  void WriteFitQSq( Common::FormFactor ff, const std::string &sPrefix ) const;
  void WriteAdjustedQSq( Common::FormFactor ff, const std::string &sPrefix ) const;
  void WriteSynthetic();
  int DoFit();
};

#endif // Continuum_hpp
