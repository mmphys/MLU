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

struct EnsembleInfo
{
  unsigned int idx;
  unsigned int aInv_L;
  unsigned int aInv_T;
};
struct EnsembleMapT : public std::map<std::string, EnsembleInfo, Common::LessCaseInsensitive>
{
  void Loaded(); // Call this once all ensembles loaded
};

struct ContinuumFit;

struct CreateParams : public Model::CreateParams
{
  const ContinuumFit &Parent;
  CreateParams( const std::vector<std::string> &OpNames, const Common::CommandLine &cl,
                const ContinuumFit &Parent );
};

struct ContinuumFit : public FitController
{
  static const std::string sPDG;
  static const std::string &FieldQSq;
  static const std::string &FieldEL;
  static constexpr scalar FourPi{ 4. * M_PI };
  static constexpr scalar InvGeV{ 1e-9 };
  static constexpr scalar InvGeVSq{ InvGeV * InvGeV };
  static constexpr scalar Lambda{ 1e9 };
  static constexpr scalar LambdaInv{ 1. / Lambda };
  static constexpr scalar LambdaInvSq{ LambdaInv * LambdaInv };
  static constexpr int NumTicks{ 200 };
  static constexpr int NumFF{ 2 };
  static constexpr const std::size_t idxUnused{ std::numeric_limits<std::size_t>::max() };

  Common::CommandLine &cl;
  const unsigned int NumD, NumE;
  const int NumSamples;
  const bool doCorr;
  const bool CovarBlock;
  const Common::FormFactor ffDefault;
  const std::string inBase;
  const std::string outBaseFileName;
  DataSet ds;
  std::array<std::string, 2> Meson;
  std::vector<std::string> OpName;
  std::vector<Model::Args> ModelArgs;
  std::unique_ptr<Fitter> f;
  EnsembleMapT EnsembleMap;

  std::array<bool, NumFF> c0Enabled, c1Enabled, ChiEnabled, FVEnabled;
  std::array<std::vector<bool>, NumFF> dEnabled, eEnabled;

  // Global parameters
  std::size_t idxfPi;
  std::size_t idxmPDGPi;
  std::size_t idxPDGH;
  std::size_t idxPDGL;
  std::size_t idxChiralPhys;
  // Per form factor parameters
  std::array<std::size_t, NumFF> idxC0, idxC1;
  std::array<std::vector<std::size_t>, NumFF> idxD, idxE; // Discretisation / energy
  std::array<std::size_t, NumFF> idxDelta;
  std::array<std::size_t, NumFF> idxPDGDStar;
  // Per ensemble parameters
  std::vector<std::size_t> idxaInv, idxmPi, idxFVSim, idxFVPhys, idxChiSim, idxChiFV;

  // Global keys
  Common::Param::Key kfPi;
  Common::Param::Key kmPDGPi;
  Common::Param::Key kPDGH;
  Common::Param::Key kPDGL;
  Common::Param::Key kChiPhys;
  // Per form factor keys
  std::array<Common::Param::Key, NumFF> kC0, kC1;
  std::array<std::vector<Common::Param::Key>, NumFF> kD, kE;
  std::array<Common::Param::Key, NumFF> kDelta;
  std::array<Common::Param::Key, NumFF> kPDGDStar;
  // Per ensemble keys
  std::vector<Common::Param::Key> kaInv, kmPi, kFVSim, kFVPhys, kChiSim, kChiFV;

protected:
  // The global constraint
  std::size_t idxConstraint;
  unsigned int WhichConstraint;
  inline const Common::Param::Key &ConstraintKey() const
  {
    if( WhichConstraint < eEnabled[idxFFPlus].size() )
      return kE[idxFFPlus][WhichConstraint];
    return kC0[idxFFPlus];
  }

  /// Do I need to compute Finite Volume corrections on each ensemble
  inline bool NeedFV() const { return ( uiFF & uiFF0 && FVEnabled[idxFF0] )
                                   || ( uiFF & uiFFPlus && FVEnabled[idxFFPlus] ); }
  /// Do I need to compute chiral corrections on each ensemble
  inline bool NeedChiral() const { return ( uiFF & uiFF0 && ChiEnabled[idxFF0] )
                                       || ( uiFF & uiFFPlus && ChiEnabled[idxFFPlus] ); }

  inline scalar DeltaF( scalar Chiral, scalar FV ) const
  {
    return -0.75 * ( Chiral + FV );
  }

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
  EnsembleFFSetT  EnsembleFFs;

  // Stats on the ensembles we loaded from files
  struct EnsembleStat
  {
    int Num;
    Common::SeedType Seed;
  };
  using EnsembleStatMap = std::map<std::string, EnsembleStat, Common::LessCaseInsensitive>;
  EnsembleStatMap EnsembleStats;
  //using ESPair = std::pair<typename EnsembleStatMap::iterator, bool>;

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

  explicit ContinuumFit( Common::CommandLine &cl );
  virtual ~ContinuumFit() {}
  void ParamsAdjust( Common::Params &mp, const Fitter &f ) override;
  void SaveParameters( Common::Params &mp, const Fitter &f ) override;
  void SetReplica( Vector &ModelParams ) const override;
  void ComputeDerived( Vector &ModelParams ) const override;
  int Run();
protected:
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
  static const std::string &GetPoleMassName( Common::FormFactor ff, const Common::FileNameAtt &fna );
  static Common::FormFactor ValidateFF( Common::FormFactor ff );
  void AddEnsemble( const std::string &Ensemble, Common::FormFactor thisFF );
  void GetEnabled( std::string sOptions );
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
  void GetMinMax( Common::FormFactor ff, scalar &Min, scalar &Max, int Loop,
                  const std::string &Field ) const;
  std::ofstream WriteHeader( const std::string &sPrefix, const std::string &FileType ) const;
  void WriteFieldName( std::ofstream &os, const std::string &FieldName ) const;
  void WriteFitQSq( Common::FormFactor ff, const std::string &sPrefix ) const;
  void WriteAdjustedQSq( Common::FormFactor ff, const std::string &sPrefix ) const;
  void WriteSynthetic();
  int DoFit();
};

#endif // Continuum_hpp
