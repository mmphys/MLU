/**

 Mike's lattice QCD utilities
 
 Source file: DataSet.hpp
 
 Copyright (C) 2019 - 2023
 
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

// Common utilities (no dependencies other than c++ stdlib)

#ifndef MLU_DataSet_hpp
#define MLU_DataSet_hpp

#include <MLU/Fold.hpp>
#include <MLU/Model.hpp>

BEGIN_COMMON_NAMESPACE

struct ConstantSource
{
  std::size_t File; // Index of the constant file in the DataSet this comes from
  const Param::Key pKey;
  ConstantSource( std::size_t File_, const Param::Key &pKey_ ) : File{File_}, pKey{pKey_} {}
};

template <typename T>
struct DataSet
{
  using SS = SampleSource;
  using ConstMap = std::map<Param::Key, ConstantSource, Param::Key::Less>;
  struct FixedParam
  {
    ConstantSource  src;   // Where to get values from
    std::size_t     Count; // How many to get
    int             idx;   // Where to put them (destination)
    FixedParam(const ConstantSource &src_, std::size_t Count_, int idx_):src{src_}, Count{Count_}, idx{idx_}{}
  };
  int NSamples;   // Number of samples we are using. These are guaranteed to exist
  int MaxSamples; // Maximum number of samples available. Guaranteed to exist. >= NSamples.
  int Extent = 0; // Number of data points in our fit (i.e. total number of elements in FitTimes)
  std::vector<Fold<T>>          corr;     // Correlator files
  std::vector<std::vector<int>> FitTimes; // The actual timeslices we are fitting to in each correlator
  inline SeedType Seed() const
  { return corr.empty() ? RandomCache::DefaultSeed() : corr[0].Seed(); }
  ConstMap constMap;
  std::vector<int> RebinSize; // Bin sizes if I've rebinned the raw data
  // Cached data for selected FitTimes
  JackBoot<T> mFitData;
  // This is where the covariance matrix comes from
  // mBinned  mCovar  Meaning
  // data     empty   Fully unfrozen. Covariance matrix built from binned data on each replica
  // data     data    Semi-frozen. mCovar scaled by variance on each replica
  // empty    data    Frozen. The same covariance matrix is used on every replica
  Matrix<T> mBinned;
  Matrix<T> mCovar;
  Vector<T> mCovarCentralMean;     // This is the mean on the central replica
  SampleSource CovarSource = SS::Bootstrap;
  int idxJackBootCovarSource = 0;
  bool bFrozenCovarSource = false;
  void SetCovarSource( SampleSource ss, int idxJackBoot, bool bFrozen )
  {
    CovarSource = ss;
    idxJackBootCovarSource = idxJackBoot;
    bFrozenCovarSource = bFrozen;
  }
protected:
  std::array<std::vector<fint>, 2> RandomCentralBuf; // No replacement. Saves having to rebuild continually
  std::array< MatrixView<fint>, 2> RandomViews;
  std::array<std::vector<fint>, 2> RandomBuffer;
  std::vector<Model<T>> constFile;// Each of the constant files (i.e. results from previous fits) I've loaded
  void AddConstant( const Param::Key &Key, std::size_t File );
  void SetValidatedFitTimes( std::vector<std::vector<int>> &&FitTimes );
  /// Cache the fit data, i.e. JackBoot[0]. Parameters describe the covariance source
  void CacheRawData();
public:
  explicit DataSet( int nSamples = 0 ) : NSamples{ nSamples } {}
  inline bool empty() const { return corr.empty() && constFile.empty(); }
  void clear();
  const Param &GetConstantParam( const ConstantSource &cs ) const;
  int  LoadCorrelator( Common::FileNameAtt &&FileAtt, unsigned int CompareFlags = COMPAT_DEFAULT,
                       const char * PrintPrefix = "  " );
  void LoadModel     ( Common::FileNameAtt &&FileAtt, const std::string &Args );
  std::vector<std::string> GetModelFilenames() const;
  void SortOpNames( std::vector<std::string> &OpNames );
  void Rebin( const std::vector<int> &RebinSize );
  /// Minimum number of samples supported by all correlators
  int NumSamples( SampleSource ss, int idxJackBoot = 0 ) const;
  inline int NumSamplesBinned() const { return NumSamples( SS::Binned, 0 ); }
  void SetFitTimes( const std::vector<std::vector<int>> &FitTimes ); // A list of all the timeslices to include
  void SetFitTimes( int tMin, int tMax ); // All fit ranges are the same
        JackBoot<T> &GetData()       { return mFitData; }
  const JackBoot<T> &GetData() const { return mFitData; }
  void GetFixed( int idx, Vector<T> &vResult, const std::vector<FixedParam> &Params ) const;
  void SaveMatrixFile( const Matrix<T> &m, const std::string &Type, const std::string &FileName,
                       std::vector<std::string> &Abbreviations,
                       const std::vector<std::string> *FileComments = nullptr,
                       const char *pGnuplotExtra = nullptr ) const;
};

END_COMMON_NAMESPACE
#endif // MLU_DataSet
