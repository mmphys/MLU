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

struct CreateParams : public Model::CreateParams
{
 const EnsembleMapT EnsembleMap;
  CreateParams( const std::vector<std::string> &OpNames, const Common::CommandLine &cl );
};

struct ContinuumFit
{
  static const std::string &FieldQSq;
  static const std::string &FieldEL;
protected:
  static constexpr scalar InvGeV{ 1e-9 };
  static constexpr scalar InvGeVSq{ InvGeV * InvGeV };
  static constexpr scalar Lambda{ 1e9 };
  static constexpr scalar InvLambda{ 1. / Lambda };
  static constexpr int NumTicks{ 200 };
  static const std::string sPDG;
  Common::CommandLine &cl;
public:
  const int NumSamples;
  const bool doCorr;
  const bool CovarBlock;
  const std::string sFFValue;
  const Common::FormFactor ff;
  const std::string inBase;
  std::string outBaseFileName;
  DataSet ds;
  std::vector<std::string> OpName;
  std::vector<Model::Args> ModelArgs;
  std::unique_ptr<Fitter> f;
  explicit ContinuumFit( Common::CommandLine &cl );
  int Run();
protected:
  EnsembleStatMap EnsembleStats;
  // Parameters from the model just created. Set by GetIndices()
  static constexpr int NumConst{ 5 };
  static constexpr std::size_t idxCUnused{ std::numeric_limits<std::size_t>::max() };
  std::array<std::size_t, NumConst> idxC;
  std::size_t idxPDGH;
  std::size_t idxPDGL;
  inline scalar EOfQSq( int rep, scalar qSq ) const
  {
    const scalar PDGH{ f->OutputModel(rep,idxPDGH) };
    const scalar PDGL{ f->OutputModel(rep,idxPDGL) };
    const scalar E{ ( PDGH * PDGH + PDGL * PDGL - qSq ) / ( 2 * PDGH ) };
    return E;
  }
  std::size_t idxDelta;
  std::string sOpNameConcat; // Sorted, concatenated list of operators in the fit for filenames
  std::array<std::string, 2> Meson;
  void LoadModels();
  void GetEnsembleStats();
  void SetEnsembleStats();
  void SetEnsembleStatSeed();
  void SortModels();
  void LoadExtra();
  void MakeOutputFilename();
  void MakeCovarBlock();
  // Get the indices of parameters from the model just created
  void GetIndices();
  void GetMinMax( scalar &Min, scalar &Max, const std::string &Field ) const;
  std::ofstream WriteHeader( const std::string &FileType ) const;
  void WriteFieldName( std::ofstream &os, const std::string &FieldName ) const;
  void WriteFitQSq() const;
  void WriteAdjustedQSq() const;
  void WriteSynthetic() const;
  int DoFit();
};

#endif // Continuum_hpp
