/**

 Mike's lattice QCD utilities
 
 Source file: Model.hpp
 
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

#ifndef MLU_Model_hpp
#define MLU_Model_hpp

#include <MLU/Sample.hpp>

BEGIN_COMMON_NAMESPACE

struct ModelBase
{
  static const std::string EnergyPrefix;
  static const std::string EDiffPrefix;
  static const std::string ConstantPrefix;
  static const char SummaryColumnPrefix[];
};

template <typename T> struct DataSet;

template <typename T>
struct Model : public Sample<T>, public ModelBase
{
  using Base = Sample<T>;
  using Traits = typename Base::Traits;
  using scalar_type = typename Traits::scalar_type;
  using value_type = typename Traits::value_type;
  int NumExponents = 0;
  int dof = 0;
  bool CovarFrozen = false;
  using SS = SampleSource;
  SS CovarSource = SS::Binned;  // Where am I building the covariance from
  std::vector<int> CovarRebin;  // Did I rebin the data before computing covariances
  // 0 = build correlation from CovarSource, then scale by the variance of the data on each bootstrap replica
  // else this is the number of bootstrap replicas to use when estimating covariance
  int CovarNumBoot = 0;
  int CovarSampleSize = 0;      // What is the appropriate parameter for m to use in the T^2 distribution
  Vector<T> Guess;
  Matrix<T> CovarIn;      // Optional. From correlation source
  Matrix<T> Covar;        // As used in the fit
  Matrix<T> Correl;       // As used in the fit
  Matrix<T> CorrelCholesky;
  Matrix<T> CovarInv;
  Matrix<T> CorrelInv;
  Matrix<T> CorrelInvCholesky;
  Matrix<T> CovarInvCholesky;
  Matrix<T> CorrelParam;    // Correlation matrix of the fit parameters
  std::vector<std::string> CorrelParamNames;
  JackBoot<T> StdErrorMean; // From all samples of binned data
  JackBoot<T> FitInput;
  JackBoot<T> ModelPrediction;
  JackBoot<T> ErrorScaled;  // (Theory - Data) * CholeskyScale

  // These variables were added in current format
  Params params;
  std::vector<std::vector<int>> FitTimes;
  std::vector<std::string> ModelType;
  std::vector<std::string> ModelArgs;

  // Helper functions
  int GetExtent() const { return static_cast<int>( ::Common::GetExtent( FitTimes ) ); };
  int NumFitParams() const { return static_cast<int>( params.NumScalars( Param::Type::Variable ) ); };
  int NumParams() const { return static_cast<int>( params.NumScalars( Param::Type::All ) ); };
  JackBootColumn<T> Column( const Param::Key &k, std::size_t Index = 0 );
  JackBootColumn<T> ColumnFrozen( const Param::Key &k, std::size_t Index = 0 );
  int NumStatColumns() const { return Base::ColumnNames.size() <= NumParams() ? 0
                                : static_cast<int>( Base::ColumnNames.size() - NumParams() ); }
  UniqueNameSet GetStatColumnNames() const;

  Model() : Base::Sample{} {}
  Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns );
  Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns,
         int CovarSampleSize_, bool CovarFrozen_, SampleSource CovarSource_,
         std::vector<int> CovarRebin_, int CovarNumBoot_ );
  const std::string & DefaultGroupName() override { return sModel; }
  inline bool NewParamsMorePrecise( bool covarFrozen_, int NumSamples ) const
  {
    // Unfreezing the covariance matrix makes a big difference!
    if( ( CovarFrozen && !covarFrozen_ ) || ( !CovarFrozen && covarFrozen_ ) )
      return CovarFrozen;
    return NumSamples > Base::NumSamples();
  }
  using Base::Read;
  void Read( const char *PrintPrefix = nullptr, std::string * pGroupName = nullptr ) override;
  void ReadAttributes( ::H5::Group &g ) override;
  void ValidateAttributes() override;
  int WriteAttributes( ::H5::Group &g ) const override;
  /**
   Make Check field in summary reflect whether each var parameter is different from zero and unique
   
   - Parameters:
     - Strictness: -1 to disable, otherwise bitmask: bit 0=When comparing to zero; bit 1=When comparing to other parameters of same set. Strictness bits are: true to make sure every single replica has same sign; False for +/- 1 sigma
   */
  bool CheckParameters( int Strictness = -1,
                        scalar_type MonotonicUpperLimit = std::numeric_limits<scalar_type>::max() );
  void SummaryComments( std::ostream & s, bool bVerboseSummary = false,
                        bool bShowColumnNames = true ) const override;
  void SummaryColumnNames( std::ostream &os ) const override;
  void SummaryContents( std::ostream &os ) const override;
  std::vector<ValWithEr<scalar_type>> GetValWithEr( const Params &ParamNames,
                                                    const UniqueNameSet &StatNames ) const;
  void WriteSummaryTD( const DataSet<T> &ds, const std::string &sOutFileName,
                       bool bVerboseSummary = false );
protected:
  int OldFormatNumExponents;
  std::vector<std::string> OldFormatOpNames;
  void ReorderOldFormat( int NumOps, int NumExponents, Vector<T> &m );
  void ReorderOldFormat( int NumOps, int NumExponents, Matrix<T> &m );
  void CommonConstruct( const std::vector<std::string> &ExtraColumns );
  void SummaryContentsPrefix( std::ostream &os ) const;
  void SummaryContentsSuffix( std::ostream &os ) const;
};

END_COMMON_NAMESPACE
#endif // MLU_Model
