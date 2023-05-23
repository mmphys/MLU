/**

 Mike's lattice QCD utilities
 
 Source file: DataSet.cpp
 
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

#include "DataSet.hpp"

BEGIN_COMMON_NAMESPACE

template <typename T>
void DataSet<T>::clear()
{
  NSamples = 0;
  Extent = 0;
  corr.clear();
  FitTimes.clear();
  constFile.clear();
  constMap.clear();
}

// Specify which times I'm fitting to, as a list of timeslices for each correlator
template <typename T>
void DataSet<T>::SetFitTimes( const std::vector<std::vector<int>> &fitTimes_ )
{
  if( fitTimes_.size() != corr.size() )
    throw std::runtime_error( std::to_string( fitTimes_.size() ) + " FitTimes but "
                             + std::to_string( corr.size() ) + " correlators" );
  std::vector<std::vector<int>> ft{ fitTimes_ };
  for( int i = 0; i < ft.size(); ++i )
  {
    if( !ft[i].empty() )
    {
      std::sort( ft[i].begin(), ft[i].end() );
      if( ft[i][0] < 0 || ft[i].back() >= corr[i].Nt()
         || std::adjacent_find( ft[i].begin(), ft[i].end() ) != ft[i].end() )
        throw std::runtime_error( "FitTimes[" + std::to_string( i ) + "]=[" + std::to_string( ft[i][0] )
                                 + "..." + std::to_string( ft[i].back() ) + "] invalid" );
    }
  }
  SetValidatedFitTimes( std::move( ft ) );
}

template <typename T>
void DataSet<T>::SetFitTimes( int tMin, int tMax )
{
  const int extent_{ tMax - tMin + 1 };
  std::vector<std::vector<int>> ft{ corr.size() };
  for( int i = 0; i < corr.size(); ++i )
  {
    if( tMin < 0 || extent_ < 1 || tMax >= corr[i].Nt() )
      throw std::runtime_error("Fit range ["+std::to_string(tMin)+", "+std::to_string(tMax)+"] invalid");
    ft[i].reserve( extent_ );
    for( int j = tMin; j <= tMax; ++j )
      ft[i].push_back( j );
  }
  SetValidatedFitTimes( std::move( ft ) );
}

// Cache the data for this data set
template <typename T>
void DataSet<T>::SetValidatedFitTimes( std::vector<std::vector<int>> &&FitTimes_ )
{
  std::size_t Extent_{ GetExtent( FitTimes_ ) };
  if( Extent_ == 0 )
    throw std::runtime_error( "Fit range empty" );
  if( Extent_ > std::numeric_limits<int>::max() )
    throw std::runtime_error( "Fit range stupidly big" );
  Extent = static_cast<int>( Extent_ );
  FitTimes = std::move( FitTimes_ );
  CacheRawData();
}

template <typename T>
void DataSet<T>::CacheRawData()
{
  SampleSource ss{ CovarSource };
  int idxJackBoot{ idxJackBootCovarSource };
  bool bFrozen{ bFrozenCovarSource };
  // How many data rows are available?
  const int NumJackBoot{ NumSamples( SS::Bootstrap, 0 ) };
  if( NumJackBoot < 1 )
    throw std::runtime_error( "DataSet<T>::CacheRawData JackBoot data unavailable" );
  // Is the requested source available
  const int NumRequestedSource{ NumSamples( ss, idxJackBoot ) };
  if( NumRequestedSource < 1 )
  {
    std::ostringstream os;
    os << "DataSet<T>::CacheRawData " << ss << '[' << idxJackBoot << "] unavailable";
    throw std::runtime_error( os.str().c_str() );
  }
  // Binned data must be available if unfrozen
  const int NumBinned{ NumSamples( SS::Binned, 0 ) };
  if( !bFrozen && NumBinned < 1 )
  {
    std::ostringstream os;
    os << "DataSet<T>::CacheRawData " << ss << '[' << idxJackBoot << "] unfrozen impossible (no binned data)";
    throw std::runtime_error( os.str().c_str() );
  }
  // Work out what to copy
  const bool FullyUnfrozen{ idxJackBoot == 0 && ss == SS::Binned && !bFrozen };
  const bool bCovarSrcData{ idxJackBoot == 0 && ss == SS::Bootstrap };
  const bool CopyCovarSrc{ !FullyUnfrozen && !bCovarSrcData };
  mFitData.resize( NumJackBoot, Extent );
  if( FullyUnfrozen )
  {
    mCovar.clear();
    mCovarCentralMean.clear();
  }
  else
  {
    mCovar.resize( Extent, Extent );
    mCovarCentralMean.resize( Extent );
  }
  if( bFrozen )
    mBinned.clear();
  else
    mBinned.resize( NumBinned, Extent );
  Matrix<T> CovarBuffer;
  if( CopyCovarSrc )
    CovarBuffer.resize( NumRequestedSource, Extent );
  // Copy the resampled data
  for( int idx = static_cast<int>( JackBootBase::idxReplicaMean ); idx < NumJackBoot; ++idx )
  {
    int dst{ 0 };
    for( int f = 0; f < corr.size(); ++f )
    {
      const JackBoot<T> &mSrc{ corr[f].getData() };
      for( int t : FitTimes[f] )
        mFitData( idx, dst++ ) = mSrc( idx, t );
    }
  }
  // Copy the binned data
  if( !bFrozen )
  {
    for( int idx = 0; idx < NumBinned; ++idx )
    {
      int dst{ 0 };
      for( int f = 0; f < corr.size(); ++f )
      {
        const Matrix<T> &mSrc{ corr[f].getBinned( 0 ) };
        for( int t : FitTimes[f] )
          mBinned( idx, dst++ ) = mSrc( idx, t );
      }
    }
  }
  // Copy the covariance source
  if( CopyCovarSrc )
  {
    for( int idx = 0; idx < NumRequestedSource; ++idx )
    {
      int dst{ 0 };
      for( int f = 0; f < corr.size(); ++f )
      {
        const Matrix<T> &mSrc{ corr[f].get( ss, idxJackBoot ) };
        for( int t : FitTimes[f] )
          CovarBuffer( idx, dst++ ) = mSrc( idx, t );
      }
    }
    // Make correlation matrix from this data
    if( ss == SampleSource::Bootstrap )
    {
      // For Bootstrap, we check environment for whether to use mean of source or mean of resamples
      const std::size_t idxMean{ JackBootBase::getCovarMeanIdx() };
      int dst{ 0 };
      for( int f = 0; f < corr.size(); ++f )
      {
        const JackBoot<T> &mSrc{ corr[f].getData( idxJackBoot ) };
        for( int t : FitTimes[f] )
          mCovarCentralMean[dst++] = mSrc( idxMean, t );
      }
    }
    else
    {
      // For binned/raw data we create the mean from the data
      JackBoot<T>::MakeMean( mCovarCentralMean, CovarBuffer );
    }
    JackBoot<T>::MakeCovar( mCovar, mCovarCentralMean, CovarBuffer, corr[0].Norm( ss ) );
    CovarBuffer.clear();
  }
  else if( bCovarSrcData )
  {
    // Make correlation matrix from this data
    mCovarCentralMean = mFitData.GetCovarMean();
    JackBoot<T>::MakeCovar( mCovar, mCovarCentralMean, mFitData.Replica, corr[0].Norm( ss ) );
  }
}

// Get the constants from the appropriate timeslice
template <typename T>
void DataSet<T>::GetFixed( int idx, Vector<T> &vResult, const std::vector<FixedParam> &Params ) const
{
  // Now copy the data into vResult
  for( const FixedParam &p : Params )
  {
    const Param &param{ GetConstantParam( p.src ) };
    std::size_t SrcIdx{ param.GetOffset( 0, Param::Type::All ) };
    for( std::size_t i = 0; i < p.Count; ++i )
      vResult[p.idx + i] = constFile[p.src.File](static_cast<std::size_t>( idx ), SrcIdx + i);
  }
}

// Write covariance matrix to file
template <typename T>
void DataSet<T>::SaveMatrixFile( const Matrix<T> &m, const std::string &Type, const std::string &Filename,
                                 std::vector<std::string> &Abbreviations,
                                 const std::vector<std::string> *FileComments, const char *pGnuplotExtra ) const
{
  // For now, all the matrices I write are square and match dataset, but this restriction can be safely relaxed
  if( m.size1 != Extent || m.size1 != m.size2 )
    throw std::runtime_error( Type + " matrix dimensions (" + std::to_string( m.size1 ) + Common::CommaSpace
                             + std::to_string( m.size2 ) + ") don't match dataset" );
  // Header describing what the covariance matrix is for and how to plot it with gnuplot
  std::ofstream s{ Filename };
  s << "# Matrix: " << Filename << "\n# Type: " << Type << "\n# Files: " << corr.size() << Common::NewLine;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    s << "# File" << f << ": " << corr[f].Name_.NameNoExt << "\n# Times" << f << ":";
    for( int t : FitTimes[f] )
      s << Common::Space << t;
    s << Common::NewLine;
    // Say which operators are in each file
    if( FileComments && !(*FileComments)[f].empty() )
      s << (*FileComments)[f]; // These are expected to have a trailing NewLine
    // Make default abbreviation - a letter for each file
    if( Abbreviations.size() < f )
      Abbreviations.emplace_back( 1, 'A' + f );
  }
  // Save a command which can be used to plot this file
  s << "# gnuplot: set xtics rotate noenhanced font 'Arial,4'; set ytics noenhanced font 'Arial,4';"
    << " set title noenhanced '" << Filename << "'; unset key; ";
  if( pGnuplotExtra && *pGnuplotExtra )
    s << pGnuplotExtra << "; ";
  s << "plot '" << Filename << "' matrix columnheaders rowheaders with image pixels\n";
  // Now save the column names. First column is matrix extent
  s << "# The next row contains column headers, starting with matrix extent\n" << Extent;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    for( int t : FitTimes[f] )
      s << Common::Space << Abbreviations[f] << t;
  }
  s << Common::NewLine << std::setprecision( std::numeric_limits<T>::max_digits10 );
  // Now print the actual covariance matrix
  int i{ 0 };
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    for( int t : FitTimes[f] )
    {
      s << Abbreviations[f] << t;
      for( int j = 0; j < Extent; ++j )
        s << Common::Space << m( i, j );
      s << Common::NewLine;
      ++i;
    }
  }
}

// Add a constant to my list of known constants - make sure it isn't already there
template <typename T>
void DataSet<T>::AddConstant( const Param::Key &Key, std::size_t File )
{
  static const char Invalid[] = " invalid";
  std::ostringstream os;
  os << "DataSet::AddConstant " << Key << Space;
  if( constMap.find( Key ) != constMap.end() )
  {
    os << "loaded from multiple model files";
    throw std::runtime_error( os.str().c_str() );
  }
  if( File >= constFile.size() )
  {
    os << "file " << File << Invalid;
    throw std::runtime_error( os.str().c_str() );
  }
  constMap.insert( { Key, ConstantSource( File, Key ) } );
}

template <typename T>
const Param &DataSet<T>::GetConstantParam( const ConstantSource &cs ) const
{
  Params::const_iterator it{ constFile[cs.File].params.find( cs.pKey ) };
  if( it == constFile[cs.File].params.cend() )
  {
    std::ostringstream os;
    os << "Constant not found: File " << cs.File << CommaSpace << cs.pKey;
    throw std::runtime_error( os.str().c_str() );
  }
  return it->second;
}

template <typename T>
int DataSet<T>::LoadCorrelator( Common::FileNameAtt &&FileAtt, unsigned int CompareFlags,
                                const char * PrintPrefix )
{
  if( corr.size() >= std::numeric_limits<int>::max() )
    throw std::runtime_error( "More than an integer's worth of correlators loaded!!??" );
  const int i{ static_cast<int>( corr.size() ) };
  corr.emplace_back();
  corr[i].SetName( std::move( FileAtt ) );
  corr[i].Read( PrintPrefix );
  // See whether this correlator is compatible with prior correlators
  if( i )
  {
    corr[0].IsCompatible( corr[i], &NSamples, CompareFlags );
    if( MaxSamples > corr[i].NumSamples() )
      MaxSamples = corr[i].NumSamples();
  }
  else
  {
    MaxSamples = corr[i].NumSamples();
    if( NSamples == 0 || NSamples > MaxSamples )
      NSamples = MaxSamples;
  }
  return i;
}

// Load a model file
template <typename T>
void DataSet<T>::LoadModel( Common::FileNameAtt &&FileAtt, const std::string &Args )
{
  // This is a pre-built model (i.e. the result of a previous fit)
  const std::size_t i{ constFile.size() };
  constFile.emplace_back();
  constFile[i].SetName( std::move( FileAtt ) );
  constFile[i].Read( "  " );
  // Keep track of minimum number of replicas across all files
  if( NSamples == 0 )
    NSamples = constFile[i].NumSamples();
  else if( NSamples > constFile[i].NumSamples() )
    NSamples = constFile[i].NumSamples();
  // Now see which parameters we want to read from the model
  if( Args.find_first_not_of( Common::WhiteSpace ) == std::string::npos )
  {
    // Load every constant in this file
    for( const Params::value_type &it : constFile[i].params )
      AddConstant( it.first, i );
  }
  else
  {
    // Load only those constants specifically asked for
    using ReadMap = std::map<Param::Key, Param::Key, Param::Key::Less>;
    using KVR = KeyValReader<Param::Key, Param::Key, Param::Key::Less>;
    ReadMap vThisArg{ KVR::Read( Common::ArrayFromString( Args ), &EqualSign, true ) };
    for( const ReadMap::value_type &it : vThisArg )
    {
      const Param::Key OldKey{ it.first };
      const Param::Key NewKey{ it.second };
      Params::iterator itP{ constFile[i].params.find( OldKey ) };
      if( itP == constFile[i].params.end() )
      {
        std::ostringstream es;
        es << "Parameter " << OldKey << " not found";
        throw std::runtime_error( es.str().c_str() );
      }
      AddConstant( NewKey.empty() ? OldKey : NewKey, i );
    }
  }
}

template <typename T>
std::vector<std::string> DataSet<T>::GetModelFilenames() const
{
  std::vector<std::string> v;
  for( const Model<T> &m : constFile )
    v.emplace_back( m.Name_.Filename );
  return v;
}

// Sort the operator names and renumber all the loaded correlators referring to them
template <typename T>
void DataSet<T>::SortOpNames( std::vector<std::string> &OpNames )
{
  int NumOps{ static_cast<int>( OpNames.size() ) };
  if( OpNames.size() > 1 )
  {
    // Sort the names
    UniqueNames OpSorted;
    for( int i = 0; i < NumOps; ++i )
      OpSorted.emplace( std::move( OpNames[i] ), i );
    // Extract the sorted names and indices (to renumber operators in correlator names)
    std::vector<std::string> SortedNames;
    std::vector<int> SortIndex( NumOps );
    SortedNames.reserve( NumOps );
    int idx{ 0 };
    for( UniqueNames::iterator it = OpSorted.begin(); it != OpSorted.end(); ++it )
    {
      SortedNames.emplace_back( it->first );
      SortIndex[it->second] = idx++;
    }
    // Renumber the operators and save the sorted operator names
    for( auto &f : corr )
      for( int i = 0; i < f.Name_.op.size(); ++i )
        f.Name_.op[i] = SortIndex[f.Name_.op[i]];
    OpNames = SortedNames;
  }
}

template <typename T>
int DataSet<T>::NumSamples( SampleSource ss, int idxJackBoot ) const
{
  int NSB{ corr.empty() ? 0 : corr[0].NumSamples( ss, idxJackBoot ) };
  for( std::size_t i = 1; NSB && i < corr.size(); ++i )
  {
    const int ThisNSB{ corr[i].NumSamples( ss, idxJackBoot ) };
    if( NSB > ThisNSB )
      NSB = ThisNSB;
  }
  return NSB;
}

template <typename T>
void DataSet<T>::Rebin( const std::vector<int> &NewSize )
{
  constexpr int DestJackBoot{ 1 };
  RebinSize.clear();
  RebinSize.reserve( corr.size() );
  for( std::size_t i = 0; i < corr.size(); ++i )
  {
    if( !corr[i].NumSamplesRaw() )
    {
      std::ostringstream os;
      os << "Raw samples unavailable rebinning corr " << i << CommaSpace << corr[i].Name_.Filename;
      throw std::runtime_error( os.str().c_str() );
    }
    // Rebin to the size specified (last bin size repeats to end).
    // bin size 0 => auto (i.e. all measurements on same config binned together)
    RebinSize.push_back( NewSize.empty() ? 0 : NewSize[i < NewSize.size() ? i : NewSize.size() - 1] );
    if( RebinSize.back() )
      corr[i].BinFixed( RebinSize.back(), DestJackBoot );
    else
      corr[i].BinAuto( DestJackBoot );
    if( corr[i].NumSamplesBinned( DestJackBoot ) != corr[0].NumSamplesBinned( DestJackBoot ) )
    {
      std::ostringstream os;
      os << "Rebinned corr " << i << " has " << corr[i].NumSamplesBinned( DestJackBoot )
         << " samples, others have " << corr[0].NumSamplesBinned( DestJackBoot );
      throw std::runtime_error( os.str().c_str() );
    }
  }
}

template class DataSet<double>;
template class DataSet<float>;
template class DataSet<std::complex<double>>;
template class DataSet<std::complex<float>>;

END_COMMON_NAMESPACE
