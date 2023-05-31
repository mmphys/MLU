/**

 Mike's lattice QCD utilities
 
 Source file: Sample.hpp
 
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

#ifndef MLU_Sample_hpp
#define MLU_Sample_hpp

#include <MLU/Utility.hpp>
#include <MLU/JackBoot.hpp>

BEGIN_COMMON_NAMESPACE

/**
 This is for a sample of anything, but usually a sample of correlators.
 There is always a central replica (specified by user, not calculated as average) in addition to NumSamples copies
 If AuxNames are specified, then there are 3 * AuxNames.count() auxiliarry records (value, low and high)
 Fields are in memory in reverse order, so that index -1 is always the central value
*/
template <typename T>
class Sample
{
public:
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  using value_type = typename Traits::value_type;
  using ValEr = ValWithEr<scalar_type>;
  static constexpr bool is_complex { Traits::is_complex };
  static constexpr int scalar_count { Traits::scalar_count };
  static constexpr int idxCentral{ -1 };
  static constexpr int NumJackBoot{ 2 }; // Number of resamples. The first is primary
protected:
  int Nt_ = 0;
//public: // Hope I don't regret this - DO NOT RESIZE THESE directly
  Matrix<T> Raw;
  std::array<  Matrix<T>, NumJackBoot> Binned;
  std::array<JackBoot<T>, NumJackBoot> Data;
protected:
  static constexpr int NumExtraSamples{ 1 }; // central replica
  std::vector<std::string> SummaryNames;
  std::vector<ValEr> m_SummaryData;
  std::vector<std::string> ColumnNames;
public:
  FileNameAtt Name_;
  std::string Ensemble; // Name of the Ensemble the data belong to
  std::string SeedMachine_; // name of the machine that ran the bootstrap
  int binSize = 1;
  int SampleSize = 0; // Number of samples (after binning) used to create bootstrap (set during bootstrap)
  std::vector<Common::ConfigCount> ConfigCount; // Info on every config in the bootstrap in order
  std::vector<std::string> FileList; // Info on every config in the bootstrap in order
protected:
  inline void AllocSummaryBuffer() { m_SummaryData.resize( SummaryNames.size() * Nt_ ); }
  inline void ValidateJackBoot( int idxJackBoot ) const
  {
    if( idxJackBoot < 0 || idxJackBoot >= NumJackBoot )
      throw std::runtime_error( "JackBoot index " + std::to_string( idxJackBoot ) + " invalid" );
  }
public:
  virtual ~Sample() {}
  explicit Sample( const std::vector<std::string> &SummaryNames_, int NumSamples = 0, int Nt = 0 )
  {
    resize( NumSamples, Nt );
    SetSummaryNames( SummaryNames_ );
  }
  explicit Sample( int NumSamples = 0, int Nt = 0 ) : Sample( sCorrSummaryNames, NumSamples, Nt ) {}
  explicit Sample( const std::string &FileName, const char *PrintPrefix = nullptr,
          std::vector<std::string> * pOpNames = nullptr, std::string * pGroupName = nullptr )
  : Sample( 0, 0 )
  {
    Read( FileName, PrintPrefix, pOpNames, pGroupName );
  }
  inline int Nt() const { return Nt_; }
  virtual int NtUnfolded() const { return Nt_; }
  void resize( int NumReplicas, int Nt );
  virtual void clear( bool bClearName = true );
  inline SeedType Seed( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return Data[idxJackBoot].Seed;
  }
  inline void SetSeed( SeedType Seed, int idxJackBoot = 0 )
  {
    ValidateJackBoot( idxJackBoot );
    Data[idxJackBoot].Seed = Seed;
  }
  inline std::string SeedString( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return Data[idxJackBoot].SeedString();
  }
  inline typename JackBootBase::Norm Norm( SampleSource ss ) const
  {
    if( ss == SampleSource::Raw || ss == SampleSource::Binned )
      return JackBootBase::Norm::RawBinned;
    if( Seed() == SeedWildcard )
      return JackBootBase::Norm::Jackknife;
    return JackBootBase::Norm::Bootstrap;
  }
  inline int NumSamples( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    if( Data[idxJackBoot].Replica.size1 > std::numeric_limits<int>::max() )
      throw std::runtime_error( "Too many resampled data rows "
                               + std::to_string( Data[idxJackBoot].Replica.size1 ) );
    return static_cast<int>( Data[idxJackBoot].Replica.size1 );
  }
  inline int NumSamplesRaw() const
  {
    if( Raw.size1 > std::numeric_limits<int>::max() )
      throw std::runtime_error( "Too many raw data rows " + std::to_string( Raw.size1 ) );
    return static_cast<int>( Raw.size1 );
  }
  inline int NumSamplesBinned( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    if( Binned[idxJackBoot].size1 > std::numeric_limits<int>::max() )
      throw std::runtime_error("Too many raw data rows " + std::to_string(Binned[idxJackBoot].size1));
    return static_cast<int>( Binned[idxJackBoot].size1 );
  }
  inline int NumSamples( SampleSource ss, int idxJackBoot = 0 ) const
  {
    if( ss == SampleSource::Bootstrap )
      return NumSamples( idxJackBoot );
    if( ss == SampleSource::Binned )
      return NumSamplesBinned( idxJackBoot );
    if( idxJackBoot )
      throw std::runtime_error( "Invalid raw idxJackBoot " + std::to_string( idxJackBoot ) );
    return NumSamplesRaw();
  }
  // Return the summary data for the nth type (0=central, 1=bias, etc)
  inline       ValEr &SummaryData( int Row, int Column )
  {
    if( Row < 0 || Row >= SummaryNames.size() )
      throw std::runtime_error( "SummaryData row " + std::to_string( Row ) + " of " + std::to_string( SummaryNames.size() ) + " invalid" );
    if( Column < 0 || Column >= Nt_ )
      throw std::runtime_error( "SummaryData column " + std::to_string( Column ) + " of " + std::to_string( Nt_ ) + " invalid" );
    return m_SummaryData[ Row * Nt_ + Column ];
  }
  inline const ValEr &SummaryData( int Row, int Column ) const
  {
    if( Row < 0 || Row >= SummaryNames.size() )
      throw std::runtime_error( "SummaryData row " + std::to_string( Row ) + " of " + std::to_string( SummaryNames.size() ) + " invalid" );
    if( Column < 0 || Column >= Nt_ )
      throw std::runtime_error( "SummaryData column " + std::to_string( Column ) + " of " + std::to_string( Nt_ ) + " invalid" );
    return m_SummaryData[ Row * Nt_ + Column ];
  }
  inline       ValEr &SummaryData( int Column = 0 )       { return SummaryData( 0, Column ); }
  inline const ValEr &SummaryData( int Column = 0 ) const { return SummaryData( 0, Column ); }
  // Return the central replica summary data for the nth column
  inline       ValEr &SummaryData( int Row, const std::string &ColumnName )
  {
    int idxField{ IndexIgnoreCase( GetColumnNames(), ColumnName ) };
    if( idxField >= GetColumnNames().size() )
      throw std::runtime_error( "SummaryData column " + ColumnName + " not available" );
    return SummaryData( Row, idxField );
  }
  inline const ValEr &SummaryData( int Row, const std::string &ColumnName ) const
  {
    int idxField{ IndexIgnoreCase( GetColumnNames(), ColumnName ) };
    if( idxField >= GetColumnNames().size() )
      throw std::runtime_error( "SummaryData column " + ColumnName + " not available" );
    return SummaryData( Row, idxField );
  }
  inline       ValEr &SummaryData( const std::string &ColumnName )
  { return SummaryData(0, ColumnName ); }
  inline const ValEr &SummaryData( const std::string &ColumnName ) const
  { return SummaryData(0, ColumnName ); }
  inline void WriteSummaryData( std::ostream &s, int idx = 0 ) const
  {
    s << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
    for( int t = 0; t < Nt_; ++t )
      s << ( t == 0 ? "" : " " ) << SummaryData( idx, t );
  }
  inline const std::vector<std::string> &GetSummaryNames() const { return SummaryNames; };
  inline void SetSummaryNames( const std::vector<std::string> &summaryNames_ )
  {
    SummaryNames = summaryNames_;
    if( is_complex )
      for( const std::string & s : summaryNames_ )
        SummaryNames.push_back( s + "_im" );
    AllocSummaryBuffer();
  }
  inline void SetSummaryNames( const std::string &sAvgName )
  {
    std::string sBias{ "bias" };
    std::vector<std::string> v{ sAvgName, sBias };
    SetSummaryNames( v );
  }
  inline void SetSummaryNames( const char *pAvgName )
  {
    if( pAvgName )
      SetSummaryNames( std::string( pAvgName ) );
    else
      SetSummaryNames( sCorrSummaryNames );
  }
  inline void SetColumnNames( const std::vector<std::string> &columnNames_ )
  {
    if( columnNames_.size() != Nt_ )
      throw std::runtime_error( "There should be names for " + std::to_string( Nt_ ) + " columns" );
    ColumnNames = columnNames_;
  }
  inline const std::vector<std::string> & GetColumnNames() const { return ColumnNames; }
  inline int GetColumnIndexNoThrow( const std::string & ColumnName ) const
  {
    int idxField{ IndexIgnoreCase( GetColumnNames(), ColumnName ) };
    if( idxField >= GetColumnNames().size() )
      idxField = -1;
    return idxField;
  }
  inline int GetColumnIndexNoThrow( const std::string & ColumnName, int idx ) const
  {
    return GetColumnIndexNoThrow( ColumnName + std::to_string( idx ) );
  }
  inline int GetColumnIndex( const std::string & ColumnName ) const
  {
    int i = GetColumnIndexNoThrow( ColumnName );
    if( i < 0 )
      throw std::runtime_error( "Column " + ColumnName + " not found" );
    return i;
  }
  inline int GetColumnIndex( const std::string & ColumnName, int idx ) const
  {
    return GetColumnIndex( ColumnName + std::to_string( idx ) );
  }
  inline bool CanRehydrate() const { return !Binned[0].Empty(); }
  JackBootColumn<T> Column( std::size_t column, int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return Data[idxJackBoot].Column( column );
  }
  inline JackBootColumn<T> Column( int column, int idxJackBoot = 0 ) const
  {
    return Column( static_cast<std::size_t>( column ), idxJackBoot );
  }
  inline JackBootColumn<T> ColumnFrozen( std::size_t column, int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return JackBootColumn<T>( Data[idxJackBoot]( JackBoot<T>::idxCentral, column ) );
  }
  inline JackBootColumn<T> ColumnFrozen( int column, int idxJackBoot = 0 ) const
  {
    return ColumnFrozen( static_cast<std::size_t>( column ), idxJackBoot );
  }
  inline       T & operator()( std::size_t i, std::size_t j )       { return Data[0]( i, j ); }
  inline const T & operator()( std::size_t i, std::size_t j ) const { return Data[0]( i, j ); }
  inline       T & operator()(         int i, std::size_t j )
  { return Data[0]( static_cast<std::size_t>( i ), j ); }
  inline const T & operator()(         int i, std::size_t j ) const
  { return Data[0]( static_cast<std::size_t>( i ), j ); }
  inline       T & operator()(         int i,         int j )
  { return Data[0]( static_cast<std::size_t>( i ), static_cast<std::size_t>( j ) ); }
  inline const T & operator()(         int i,         int j ) const
  { return Data[0]( static_cast<std::size_t>( i ), static_cast<std::size_t>( j ) ); }
  inline const JackBoot<T> &getData( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return Data[idxJackBoot];
  }
  inline       JackBoot<T> &getData( int idxJackBoot = 0 )
  {
    ValidateJackBoot( idxJackBoot );
    return Data[idxJackBoot];
  }
  inline       Matrix<T> &getRaw()       { return Raw; }
  inline const Matrix<T> &getRaw() const { return Raw; }
  inline       Matrix<T> &resizeRaw( int NumSamplesRaw )
  {
    if( !NumSamplesRaw )
      Raw.clear();
    else if( !Nt_ )
      throw std::runtime_error( "Cannot resizeRaw() when Nt==0" );
    else
      Raw.resize( NumSamplesRaw, Nt_ );
    return getRaw();
  }
  inline       Matrix<T> &getBinned( int idxJackBoot = 0 )
  {
    ValidateJackBoot( idxJackBoot );
    return Binned[idxJackBoot];
  }
  inline const Matrix<T> &getBinned( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return Binned[idxJackBoot];
  }
  inline       Matrix<T> &get( SampleSource ss, int idxJackBoot = 0 )
  {
    ValidateJackBoot( idxJackBoot );
    if( ss == SampleSource::Bootstrap )
      return Data[idxJackBoot].Replica;
    if( ss == SampleSource::Binned )
      return Binned[idxJackBoot];
    if( idxJackBoot )
      throw std::runtime_error( "Invalid raw idxJackBoot " + std::to_string( idxJackBoot ) );
    return Raw;
  }
  inline const Matrix<T> &get( SampleSource ss, int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    if( ss == SampleSource::Bootstrap )
      return Data[idxJackBoot].Replica;
    if( ss == SampleSource::Binned )
      return Binned[idxJackBoot];
    if( idxJackBoot )
      throw std::runtime_error( "Invalid raw idxJackBoot " + std::to_string( idxJackBoot ) );
    return Raw;
  }
  inline       Matrix<T> &resizeBinned( int NumSamplesBinned, int idxJackBoot = 0 )
  {
    ValidateJackBoot( idxJackBoot );
    if( !NumSamplesBinned )
      Binned[idxJackBoot].clear();
    else if( !Nt_ )
      throw std::runtime_error( "Cannot resizeBinned() when Nt==0" );
    else
      Binned[idxJackBoot].resize( NumSamplesBinned, Nt_ );
    return getBinned( idxJackBoot );
  }
  inline bool IsFinite( int idxJackBoot = 0 ) const
  {
    ValidateJackBoot( idxJackBoot );
    return ::Common::IsFinite( Data[idxJackBoot] );
  }
  void WriteColumnNames( std::ostream &s ) const;
  /// Initialise *pNumSamples either to 0, or to the size of the first Sample before first call
  template <typename U>
  void IsCompatible( const Sample<U> &o, int * pNumSamples = nullptr, unsigned int CompareFlags = COMPAT_DEFAULT, bool bSuffix = true ) const;
  template <typename U> void CopyAttributes( const Sample<U> &o );
  void BinAuto( int JackBootDest = 0 ); // Auto bin based on config count
  void BinFixed( int binSize_, int JackBootDest = 0 );
  void Resample( int JackBootDest = 0,
                  int NumReplicas = static_cast<int>( RandomCache::DefaultNumReplicas() ),
                  int BinnedSource = 0, SeedType Seed = RandomCache::DefaultSeed() );
  /// Is it possible to regenerate primary bootstrap replica for any seed?
  inline bool CanResample() const { return !Binned[0].Empty(); }
  virtual void SetName( const std::string &FileName, std::vector<std::string> * pOpNames = nullptr )
  {
    Name_.Parse( FileName, pOpNames );
  }
  inline void SetName( FileNameAtt &&FileName ) { Name_ = std::move( FileName ); }
  void Read( const char *PrintPrefix = nullptr, std::string * pGroupName = nullptr );
  inline void Read( const std::string &FileName, const char *PrintPrefix = nullptr,
                    std::vector<std::string> * pOpNames=nullptr, std::string * pGroupName=nullptr )
  {
    SetName( FileName, pOpNames );
    Read( PrintPrefix, pGroupName );
  }
  void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void MakeCorrSummary();
  void WriteSummary( const std::string &sOutFileName, bool bVerboseSummary = false );
public: // Override these for specialisations
  virtual const std::string & DefaultGroupName() { return sBootstrap; }
  virtual bool bFolded() { return false; }
  // Descendants should call base first
  virtual void SummaryComments( std::ostream & s, bool bVerboseSummary = false ) const;
  virtual void SummaryColumnNames( std::ostream &os ) const;
  /**
   Write a summary of this sample, showing central value, +/- errors and abs min/max

   If there are no ColumnNames, write one row per timeslice.
   Otherwise write all columns on a single row
   
   - Important: There is no final NewLine written
   */
  virtual void SummaryContents( std::ostream &os ) const;
  virtual void ReadAttributes( ::H5::Group &g ) {}
  virtual void ValidateAttributes() {} // Called once data read to validate attributes against data
  virtual int WriteAttributes( ::H5::Group &g ) const { return 0; }
};

using SampleC = Sample<std::complex<double>>;
using SampleD = Sample<double>;

END_COMMON_NAMESPACE
#endif // MLU_Sample_hpp
