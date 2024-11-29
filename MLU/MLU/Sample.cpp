/**

 Mike's lattice QCD utilities
 
 Source file: Sample.cpp
 
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

#include <MLUconfig.h>
#include "Sample.hpp"

BEGIN_MLU_NAMESPACE

template <typename T>
void Sample<T>::resize( int NumReplicas, int Nt )
{
  const bool NtChanged{ Nt_ != Nt };
  if( NtChanged || NumSamples() != NumReplicas )
  {
    if( NtChanged )
    {
      Nt_ = Nt;
      AllocSummaryBuffer();
      ColumnNames.clear();
      Raw.clear();
      Binned[0].clear();
      Data[0].clear();
    }
    Data[0].resize( NumReplicas, Nt );
    for( int i = 1; i < NumJackBoot; ++i )
    {
      Binned[i].clear();
      Data[i].clear();
    }
  }
}

template <typename T>
void Sample<T>::clear( bool bClearName )
{
  Nt_ = 0;
  Raw.clear();
  for( int i = 0; i < NumJackBoot; ++i )
  {
    Binned[i].clear();
    Data[i].clear();
  }
  SummaryNames.clear();
  m_SummaryData.clear();
  ColumnNames.clear();
  SeedMachine_.clear();
  binSize = 1;
  SampleSize = 0;
  ConfigCount.clear();
  FileList.clear();
  if( bClearName )
    Name_.clear();
}

template <typename T>
void Sample<T>::WriteColumnNames( std::ostream &s ) const
{
  const std::string Sep{ " " };
  if( Nt_ == 0 )
    throw std::runtime_error( "Can't write header - Nt_ = 0" );
  std::string Buffer{ "t" };
  const std::size_t BufferLen{ Buffer.size() };
  const std::string * ps{ &Buffer };
  for( std::size_t t = 0; t < Nt_; t++ )
  {
    if( ColumnNames.empty() )
    {
      Buffer.resize( BufferLen );
      Buffer.append( std::to_string( t ) );
    }
    else
      ps = &ColumnNames[t];
    if( t )
      s << Sep;
    MLU::ValWithEr<T>::Header( *ps, s, Sep );
  }
}

// Initialise *pNumSamples either to 0, or to the size of the first Sample before first call
template <typename T> template <typename U>
void Sample<T>::IsCompatible( const Sample<U> &o, int * pNumSamples, unsigned int CompareFlags,
                              bool bSuffix ) const
{
  static const std::string sPrefix{ "Incompatible " + sBootstrap + " samples - " };
  const std::string sSuffix{ bSuffix ? ":\n  " + Name_.Filename + "\nv " + o.Name_.Filename : "" };
  if( !( CompareFlags & COMPAT_DISABLE_BASE ) && !MLU::EqualIgnoreCase( o.Name_.Base, Name_.Base ) )
    throw std::runtime_error( sPrefix + "base " + o.Name_.Base + sNE + Name_.Base + sSuffix );
  if( !( CompareFlags & COMPAT_DISABLE_TYPE ) && !MLU::EqualIgnoreCase( o.Name_.Type, Name_.Type ) )
    throw std::runtime_error( sPrefix + "type " + o.Name_.Type + sNE + Name_.Type + sSuffix );
  if( !MLU::EqualIgnoreCase( o.Name_.Ext, Name_.Ext ) )
    throw std::runtime_error( sPrefix + "extension " + o.Name_.Ext + sNE + Name_.Ext + sSuffix );
  if( pNumSamples )
  {
    if( *pNumSamples == 0 || *pNumSamples > NumSamples() )
      *pNumSamples = NumSamples();
    if( *pNumSamples > o.NumSamples() )
      *pNumSamples = o.NumSamples();
  }
  else if( o.NumSamples() != NumSamples() )
    throw std::runtime_error( sPrefix + "NumSamples " + std::to_string(o.NumSamples()) +
                             sNE + std::to_string(NumSamples()) + sSuffix );
  if( !( CompareFlags & COMPAT_DISABLE_NT ) && o.Nt() != Nt_ )
    throw std::runtime_error( sPrefix + "Nt " + std::to_string(o.Nt()) +
                             sNE + std::to_string(Nt_) + sSuffix );
  if( ! ( CompareFlags & COMPAT_DISABLE_ENSEMBLE ) )
  {
    if( o.SampleSize && SampleSize && o.SampleSize != SampleSize )
      throw std::runtime_error( sPrefix + "SampleSize " + std::to_string(o.SampleSize) +
                               sNE + std::to_string(SampleSize) + sSuffix );
    // NB: If the random numbers were different, the random number cache would have caught this
    if( o.Seed() != Seed() )
      throw std::runtime_error( "Seed " + o.SeedString() + sNE + SeedString() );
    const std::size_t CSize {   ConfigCount.size() };
    const std::size_t CSizeO{ o.ConfigCount.size() };
    if( CSize && CSizeO )
    {
      if( CSizeO != CSize )
        throw std::runtime_error( sPrefix + "Number of configs " +
                                 std::to_string(CSizeO) + sNE + std::to_string(CSize) + sSuffix );
      for( std::size_t i = 0; i < CSize; i++ )
      {
        const MLU::ConfigCount &l{   ConfigCount[i] };
        const MLU::ConfigCount &r{ o.ConfigCount[i] };
        if( r.Config != l.Config )
          throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                   sNE + std::to_string(l.Config) + sSuffix );
      }
    }
    if( !Ensemble.empty() && !o.Ensemble.empty() && !EqualIgnoreCase( Ensemble, o.Ensemble ) )
      throw std::runtime_error( "Ensemble " + Ensemble + " != " + o.Ensemble );
  }
}

template void Sample<float>::IsCompatible( const Sample<float> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<float>::IsCompatible( const Sample<double> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<float>::IsCompatible( const Sample<std::complex<float>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<float>::IsCompatible( const Sample<std::complex<double>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<double>::IsCompatible( const Sample<float> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<double>::IsCompatible( const Sample<double> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<double>::IsCompatible( const Sample<std::complex<float>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<double>::IsCompatible( const Sample<std::complex<double>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<float>>::IsCompatible( const Sample<float> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<float>>::IsCompatible( const Sample<double> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<float>>::IsCompatible( const Sample<std::complex<float>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<float>>::IsCompatible( const Sample<std::complex<double>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<double>>::IsCompatible( const Sample<float> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<double>>::IsCompatible( const Sample<double> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<double>>::IsCompatible( const Sample<std::complex<float>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;
template void Sample<std::complex<double>>::IsCompatible( const Sample<std::complex<double>> &o,
                              int * pNumSamples, unsigned int CompareFlags, bool bSuffix ) const;

// Take ownership of the FileList
template <typename T>
template <typename U> void Sample<T>::CopyAttributes( const Sample<U> &in )
{
  binSize = in.binSize;
  SampleSize = in.SampleSize;
  ConfigCount = in.ConfigCount;
  Ensemble = in.Ensemble;
  if( FileList.empty() )
    FileList = in.FileList;
  // Copy the random numbers for this sample
  SeedMachine_ = in.SeedMachine_;
  SetSeed( in.Seed() );
}

template void Sample<float>::CopyAttributes( const Sample<float> &in );
template void Sample<float>::CopyAttributes( const Sample<double> &in );
template void Sample<float>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<float>::CopyAttributes( const Sample<std::complex<double>> &in );
template void Sample<double>::CopyAttributes( const Sample<float> &in );
template void Sample<double>::CopyAttributes( const Sample<double> &in );
template void Sample<double>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<double>::CopyAttributes( const Sample<std::complex<double>> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<float> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<double> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<std::complex<float>>::CopyAttributes( const Sample<std::complex<double>> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<float> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<double> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<std::complex<float>> &in );
template void Sample<std::complex<double>>::CopyAttributes( const Sample<std::complex<double>> &in );

// Auto bin: all measurements on each config binned into a single measurement
template <typename T> void Sample<T>::BinAuto( int JackBootDest )
{
  if( !NumSamplesRaw() )
    throw std::runtime_error( "No raw data to bin" );
  // Work out how many bins there are and whether they are all the same size
  int NumSamplesRaw{ 0 };
  int NewBinSize{ ConfigCount[0].Count };
  for( const MLU::ConfigCount &cc : ConfigCount )
  {
    if( cc.Count < 1 )
      throw std::runtime_error( "ConfigCount invalid" );
    NumSamplesRaw += cc.Count;
    if( NewBinSize && NewBinSize != cc.Count )
      NewBinSize = 0; // Indicates the bin size varies per ConfigCount
  }
  if( NumSamplesRaw != this->NumSamplesRaw() )
    throw std::runtime_error( "ConfigCount doesn't match raw data" );
  // Rebin
  Matrix<T> &mSrc = getRaw();
  Matrix<T> &mDst = resizeBinned( static_cast<int>( ConfigCount.size() ), JackBootDest );
  if( JackBootDest == 0 )
    binSize = NewBinSize; // File attribute is based on the primary sample
  Vector<T> vSrc, vDst;
  std::size_t rowSrc = 0;
  for( std::size_t Bin = 0; Bin < ConfigCount.size(); ++Bin )
  {
    for( std::size_t i = 0; i < ConfigCount[Bin].Count; ++i )
    {
      for( std::size_t t = 0; t < Nt_; ++t )
      {
        if( i )
          mDst( Bin, t ) += mSrc( rowSrc, t );
        else
          mDst( Bin, t ) = mSrc( rowSrc, t );
      }
      ++rowSrc;
    }
    for( std::size_t t = 0; t < Nt_; ++t )
      mDst( Bin, t ) /= ConfigCount[Bin].Count;
  }
}

template <typename T> void Sample<T>::BinFixed( int binSize_, int JackBootDest )
{
  if( binSize_ < 1 )
    throw std::runtime_error( "BinSize " + std::to_string( binSize_ ) + " invalid" );
  if( !NumSamplesRaw() )
    throw std::runtime_error( "No raw data to bin" );
  // Work out how many samples there are and check they match ConfigCount
  int NumSamplesRaw{ 0 };
  for( const MLU::ConfigCount &cc : ConfigCount )
  {
    if( cc.Count != ConfigCount[0].Count )
      throw std::runtime_error( "Can't manually bin - raw sample counts differ per config" );
    NumSamplesRaw += cc.Count;
  }
  if( NumSamplesRaw != this->NumSamplesRaw() )
    throw std::runtime_error( "ConfigCount doesn't match raw data" );
  const bool PartialLastBin = NumSamplesRaw % binSize_;
  const int NewNumSamples{ static_cast<int>( NumSamplesRaw / binSize_ + ( PartialLastBin ? 1 : 0 ) ) };
  // Rebin
  Matrix<T> &mSrc = getRaw();
  Matrix<T> &mDst = resizeBinned( NewNumSamples, JackBootDest );
  if( JackBootDest == 0 )
    binSize = binSize_; // File attribute is based on the primary sample
  Vector<T> vSrc, vDst;
  std::size_t rowSrc = 0;
  int ThisBinSize{ binSize_ };
  for( std::size_t Bin = 0; Bin < NewNumSamples; ++Bin )
  {
    if( PartialLastBin && Bin == NewNumSamples - 1 )
      ThisBinSize = NumSamplesRaw % binSize_;
    for( std::size_t i = 0; i < ThisBinSize; ++i )
    {
      for( std::size_t t = 0; t < Nt_; ++t )
      {
        if( i )
          mDst( Bin, t ) += mSrc( rowSrc, t );
        else
          mDst( Bin, t ) = mSrc( rowSrc, t );
      }
      ++rowSrc;
    }
    for( std::size_t t = 0; t < Nt_; ++t )
      mDst( Bin, t ) /= ThisBinSize;
  }
}

/// Perform bootstrap or jackknife, depending on Seed
template <typename T>
void Sample<T>::Resample( int JackBootDest, int NumReplicas, int BinnedSource, SeedType Seed )
{
  ValidateJackBoot( JackBootDest );
  ValidateJackBoot( BinnedSource );
  Data[JackBootDest].Seed = Seed;
  Data[JackBootDest].Resample( Binned[BinnedSource], NumReplicas );
  if( JackBootDest == 0 )
  {
    SampleSize = NumSamplesBinned( BinnedSource );
    SeedMachine_ = GetHostName();
  }
}

template<typename T> inline typename std::enable_if<!SampleTraits<T>::is_complex>::type
CopyOldFormat( T &Dest, const double &Real, const double &Imag )
{
  if( Imag )
    throw std::runtime_error( "Complex sample has imaginary component" );
  using Scalar = typename SampleTraits<T>::scalar_type;
  Dest = static_cast<Scalar>( Real );
}

template<typename T> inline typename std::enable_if<SampleTraits<T>::is_complex>::type
CopyOldFormat( T &Dest, const double &Real, const double &Imag )
{
  using Scalar = typename SampleTraits<T>::scalar_type;
  Dest.real( static_cast<Scalar>( Real ) );
  Dest.imag( static_cast<Scalar>( Imag ) );
}

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void Sample<T>::Read( bool bSetSeed, SeedType NewSeed,
                      const char *PrintPrefix, std::string *pGroupName )
{
  /*if( !Name_.bSeedNum )
    throw std::runtime_error( "Seed missing from " + Name_.Filename );
  if( Name_.Type.empty() )
    throw std::runtime_error( "Type missing from " + Name_.Filename );*/
  ::H5::H5File f;
  ::H5::Group  g;
  H5::OpenFileGroup( f, g, Name_.Filename, PrintPrefix, pGroupName );
  bool bOK{ false };
  bool bNoReturn{ false };
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  clear( false );
  try // to load from LatAnalyze format
  {
    SetSeed( Name_.Seed ); // Seed comes from filename if not loaded from hdf5
    unsigned short att_Type;
    ::H5::Attribute a;
    a = g.openAttribute("type");
    a.read( ::H5::PredType::NATIVE_USHORT, &att_Type );
    a.close();
    if( att_Type == 2 && ( !bSetSeed || Name_.Seed == NewSeed ) )
    {
      unsigned long  att_nSample;
      a = g.openAttribute("nSample");
      a.read( ::H5::PredType::NATIVE_ULONG, &att_nSample );
      a.close();
      std::vector<double> buffer;
      for( int i = idxCentral; i < NumSamples(); i++ )
      {
        ::H5::DataSet ds = g.openDataSet( "data_" + ( i == idxCentral ? "C" : "S_" + std::to_string( i ) ) );
        ::H5::DataSpace dsp = ds.getSpace();
        if( dsp.getSimpleExtentNdims() == 2 )
        {
          hsize_t Dim[2];
          dsp.getSimpleExtentDims( Dim );
          if( Dim[1] != 2 || Dim[0] * att_nSample > std::numeric_limits<int>::max() )
            break;
          if( i == idxCentral )
          {
            resize( static_cast<int>( att_nSample ), static_cast<int>( Dim[0] ) );
            buffer.resize( 2 * Nt_ );
          }
          else if( Dim[0] != static_cast<unsigned int>( Nt_ ) )
            break;
          ds.read( buffer.data(), ::H5::PredType::NATIVE_DOUBLE );
          if( !MLU::IsFinite( buffer ) )
            throw std::runtime_error( "NANs at row " + std::to_string(i) + "in " + Name_.Filename );
          for( int t = 0; t < Nt_; t++ )
            CopyOldFormat( Data[0]( i, t ), buffer[t], buffer[t + Nt_] );
          // Check whether this is the end
          if( i == 0 )
            bNoReturn = true;
          if( i == NumSamples() - 1 )
            bOK = true;
        }
      }
    }
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
    ::H5::Exception::clearErrorStack();
  }
  catch(...)
  {
    ::H5::Exception::setAutoPrint(h5at, f5at_p);
    clear( false );
    throw;
  }
  if( !bOK && !bNoReturn )
  {
    clear( false );
    try // to load from my format
    {
      unsigned short att_nAux;
      ::H5::Attribute a;
      // Auxilliary rows removed 30 Apr 2023. Was never used
      a = g.openAttribute("nAux");
      a.read( ::H5::PredType::NATIVE_USHORT, &att_nAux );
      a.close();
      if( att_nAux )
        throw std::runtime_error( std::to_string( att_nAux ) + " unexpected auxilliary rows in "
                                 + Name_.Filename );
      unsigned int att_nSample;
      a = g.openAttribute("nSample");
      a.read( ::H5::PredType::NATIVE_UINT, &att_nSample );
      a.close();
      if( att_nSample > std::numeric_limits<int>::max() )
        throw std::runtime_error( "nSample " + std::to_string( att_nSample ) + " invalid" );
      try
      {
        SeedType tmp;
        a = g.openAttribute( RandomCache::sSeed );
        a.read( H5::Equiv<SeedType>::Type, &tmp );
        a.close();
        SetSeed( tmp );
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        SetSeed( 0 ); // In this format, not present means 0
      }
      // If I'm not trying to load a specific seed, just use the one from the file
      if( !bSetSeed )
        NewSeed = Seed();
      try
      {
        a = g.openAttribute( sEnsemble );
        a.read( a.getStrType(), Ensemble );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        Ensemble.clear();
      }
      try
      {
        a = g.openAttribute( RandomCache::sSeedMachine );
        a.read( a.getStrType(), SeedMachine_ );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        SeedMachine_.clear();
      }
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sSummaryNames);
        SummaryNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        // Auxiliary names should match the number of auxiliary records
        a = g.openAttribute(sColumnNames);
        ColumnNames = H5::ReadStrings( a );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        int tmp;
        a = g.openAttribute( sSampleSize );
        a.read( ::H5::PredType::NATIVE_INT, &tmp );
        a.close();
        SampleSize = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        int tmp;
        a = g.openAttribute( sBinSize );
        a.read( ::H5::PredType::NATIVE_INT, &tmp );
        a.close();
        binSize = tmp;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      // Try to load the FileList from dataSet first, falling back to attribute
      bool bGotFileList{ false };
      try
      {
        ::H5::DataSet ds = g.openDataSet( sFileList );
        FileList = H5::ReadStrings( ds );
        ds.close();
        bGotFileList = true;
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      if( !bGotFileList )
      {
        try
        {
          a = g.openAttribute(sFileList);
          FileList = H5::ReadStrings( a );
          a.close();
          bGotFileList = true;
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
        }
      }
      try
      {
        a = g.openAttribute(sConfigCount);
        ::H5::DataSpace dsp = a.getSpace();
        const int rank{ dsp.getSimpleExtentNdims() };
        if( rank != 1 )
          throw std::runtime_error( sConfigCount + " number of dimensions=" + std::to_string( rank ) + ", expecting 1" );
        hsize_t NumConfig;
        dsp.getSimpleExtentDims( &NumConfig );
        if( NumConfig > std::numeric_limits<int>::max() )
          throw std::runtime_error( sConfigCount + " too many items in ConfigCount: " + std::to_string( NumConfig ) );
        ConfigCount.resize( NumConfig );
        a.read( H5::Equiv<MLU::ConfigCount>::Type, &ConfigCount[0] );
        a.close();
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
      }
      try
      {
        ReadAttributes( g );
      }
      catch(const std::exception &e)
      {
        throw;
      }
      catch(...)
      {
        throw std::runtime_error( "Extended attributes for " + DefaultGroupName() + " missing" );
      }
      // Load Binned data
      try
      {
        H5::ReadMatrix( g, "samples_Binned", Binned[0] );
        if( Binned[0].size1 != SampleSize )
          throw std::runtime_error( "SampleSize = " + std::to_string( SampleSize ) + " but "
                      + std::to_string( Binned[0].size1 ) + " binned replicas read" );
        if( Nt_ == 0 )
          Nt_ = static_cast<int>( Binned[0].size2 );
        else if( Binned[0].size2 != Nt_ )
          throw std::runtime_error( "nT = " + std::to_string( Nt_ ) + " but "
                      + std::to_string( Binned[0].size2 ) + " binned columns read" );
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        Binned[0].clear();
      }
      // Load raw data
      try
      {
        H5::ReadMatrix( g, "samples_Raw", Raw );
        if( Nt_ == 0 )
          Nt_ = static_cast<int>( Raw.size2 );
        else if( Raw.size2 != Nt_ )
          throw std::runtime_error( "nT = " + std::to_string( Nt_ ) + " but "
                      + std::to_string( Raw.size2 ) + " raw columns read" );
        if( FileList.size() && Raw.size1 != FileList.size() )
          std::cout << "* Warning: FileList has " << FileList.size() << " entries but "
                    << Raw.size1 << " raw replicas read\n";
      }
      catch(const ::H5::Exception &)
      {
        ::H5::Exception::clearErrorStack();
        Raw.clear();
      }
      // If random numbers present, put them in cache (which checks compatibility)
      if( Seed() == NewSeed )
      {
        try
        {
          Matrix<fint> Random;
          H5::ReadMatrix( g, RandomCache::sRandom, Random );
          if( Random.size2 != SampleSize )
            throw std::runtime_error( "SampleSize = " + std::to_string( SampleSize ) + " but "
                                     + std::to_string( Random.size2 ) + " random columns read" );
          if( Random.size1 != att_nSample )
            throw std::runtime_error( "nSample = " + std::to_string( att_nSample ) + " but "
                                     + std::to_string( Random.size1 ) + " random replicas read" );
          // Put random numbers in cache - will throw if they don't match previous entries
          RandomCache::Global.Put( Seed(), Random );
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
        }
      }
      const int NewNumReplicas{ NewSeed == SeedWildcard ? SampleSize // Jackknife
                                      : static_cast<int>( RandomCache::DefaultNumReplicas() ) };
      // If I don't need to resample, load optional items
      bool bNeedResample{ Seed() != NewSeed || att_nSample < NewNumReplicas };
      if( !bNeedResample )
      {
        try
        {
          // Load resampled data
          Data[0].Read( g, "data" );
          if( Data[0].NumReplicas() != att_nSample )
            throw std::runtime_error( "nSample = " + std::to_string( att_nSample ) + " but "
                                     + std::to_string( Data[0].NumReplicas() ) + " replicas read" );
          if( Nt_ == 0 )
            Nt_ = static_cast<int>( Data[0].extent() );
          else if( Data[0].extent() != Nt_ )
            throw std::runtime_error( "nT = " + std::to_string( Nt_ ) + " but "
                        + std::to_string( Data[0].extent() ) + " columns read" );
          if( !SummaryNames.empty() )
          {
            AllocSummaryBuffer();
            try // to load the summaries
            {
              ::H5::DataSet ds = g.openDataSet( sSummaryDSName );
              ::H5::DataSpace dsp = ds.getSpace();
              if( dsp.getSimpleExtentNdims() != 2 )
                throw std::runtime_error( "Dimension error reading " + sSummaryDSName );
              hsize_t Dim[2];
              dsp.getSimpleExtentDims( Dim );
              if( Dim[0] != SummaryNames.size() || Dim[1] != static_cast<unsigned int>( Nt_ ) )
                throw std::runtime_error( "Bad size reading " + sSummaryDSName );
              const ::H5::DataType ErrType{ ds.getDataType() };
              if( ErrType == H5::Equiv<ValWithErOldV1<scalar_type>>::Type )
              {
                const std::size_t Count{ SummaryNames.size() * Nt_ };
                std::vector<ValWithErOldV1<scalar_type>> tmp( Count );
                ds.read( &tmp[0], ErrType );
                for( std::size_t i = 0; i < Count; ++i )
                {
                  m_SummaryData[i].Central = tmp[i].Central;
                  m_SummaryData[i].Low     = tmp[i].Low;
                  m_SummaryData[i].High    = tmp[i].High;
                  m_SummaryData[i].Check   = tmp[i].Check;
                  m_SummaryData[i].Min     = 0;
                  m_SummaryData[i].Max     = 0;
                }
              }
              else
                ds.read( &m_SummaryData[0], H5::Equiv<ValWithEr<scalar_type>>::Type );
            }
            catch(const ::H5::Exception &)
            {
              ::H5::Exception::clearErrorStack();
              MakeCorrSummary(); // Rebuild summaries if I can't load them
            }
          }
        }
        catch(const ::H5::Exception &)
        {
          // I couldn't load the data
          ::H5::Exception::clearErrorStack();
          bNeedResample = true;
        }
      }
      // Either resampled or binned data must be available
      if( Data[0].NumReplicas() == 0 && Binned[0].Empty() )
        throw std::runtime_error( "Neither binned nor resampled data available" );
      if( bNeedResample )
      {
        if( !CanResample() )
        {
          static const char szTo[]{ " -> " };
          std::ostringstream os;
          os << "Can't resample";
          if( att_nSample < NewNumReplicas )
            os << Space << att_nSample << szTo << NewNumReplicas << " replicas";
          if( Seed() != Name_.Seed)
            os << " seed " << SeedString() << szTo
            << RandomCache::SeedString( Name_.Seed );
          throw std::runtime_error( os.str().c_str() );
        }
        Data[0].Seed = NewSeed;
        Data[0].Resample( Binned[0], NewNumReplicas );
        MakeCorrSummary();
      }
      // If we get here then we've loaded OK
      ValidateAttributes();
      bOK = true;
    }
    catch(const ::H5::Exception &)
    {
      bOK = false;
      ::H5::Exception::clearErrorStack();
    }
    catch(...)
    {
      ::H5::Exception::setAutoPrint(h5at, f5at_p);
      clear( false );
      throw;
    }
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  if( !bOK )
  {
    clear( false );
    throw std::runtime_error( "Unable to read sample from " + Name_.Filename );
  }
}

template <typename T>
void Sample<T>::Read( const char *PrintPrefix, std::string *pGroupName )
{
  Read( true, RandomCache::DefaultSeed(), PrintPrefix, pGroupName );
}

template <typename T>
void Sample<T>::Write( const char * pszGroupName, unsigned int Flags )
{
  const std::string GroupName{ pszGroupName == nullptr || *pszGroupName == 0 ? DefaultGroupName() : pszGroupName };
  bool bOK = false;
  try // to write in my format
  {
    ::H5::H5File f( Name_.Filename, Flags );
    hsize_t Dims[2];
    Dims[0] = 1;
    ::H5::DataSpace ds1( 1, Dims );
    ::H5::Group g = f.createGroup( GroupName );
    ::H5::Attribute a = g.createAttribute( "nAux", ::H5::PredType::STD_U16LE, ds1 );
    int tmp{ NumExtraSamples - 1 };
    a.write( ::H5::PredType::NATIVE_INT, &tmp );
    a.close();
    a = g.createAttribute( "nSample", ::H5::PredType::STD_U32LE, ds1 );
    tmp = NumSamples();
    a.write( ::H5::PredType::NATIVE_INT, &tmp );
    a.close();
    int NumAttributes = 2;
    if( Seed() )
    {
      const SeedType tmp{ Seed() };
      a = g.createAttribute( RandomCache::sSeed, ::H5::PredType::STD_U32LE, ds1 );
      a.write( H5::Equiv<SeedType>::Type, &tmp );
      a.close();
      NumAttributes++;
    }
    if( !Ensemble.empty() )
    {
      a = g.createAttribute( sEnsemble, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, Ensemble );
      a.close();
      NumAttributes++;
    }
    if( !SeedMachine_.empty() )
    {
      a = g.createAttribute( RandomCache::sSeedMachine, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, SeedMachine_ );
      a.close();
      NumAttributes++;
    }
    if( SummaryNames.size() )
    {
      H5::WriteAttribute( g, sSummaryNames, SummaryNames );
      NumAttributes++;
    }
    if( ColumnNames.size() )
    {
      H5::WriteAttribute( g, sColumnNames, ColumnNames );
      NumAttributes++;
    }
    if( SampleSize )
    {
      a = g.createAttribute( sSampleSize, ::H5::PredType::STD_U32LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &SampleSize );
      a.close();
      NumAttributes++;
    }
    if( binSize != 1 )
    {
      a = g.createAttribute( sBinSize, ::H5::PredType::STD_I32LE, ds1 );
      a.write( ::H5::PredType::NATIVE_INT, &binSize );
      a.close();
      NumAttributes++;
    }
    ::H5::DataSpace dsp;
    if( ConfigCount.size() )
    {
      Dims[0] = ConfigCount.size();
      dsp = ::H5::DataSpace( 1, Dims );
      a = g.createAttribute( sConfigCount, H5::Equiv<MLU::ConfigCount>::Type, dsp );
      a.write( H5::Equiv<MLU::ConfigCount>::Type, &ConfigCount[0] );
      a.close();
      dsp.close();
      NumAttributes++;
    }
    NumAttributes += WriteAttributes( g );
    if( FileList.size() )
      H5::WriteStringData( g, sFileList, FileList );
    // If I don't have binned data, I must save the bootstrap (because it can't be recreated)
    if( Binned[0].Empty() || RandomCache::SaveFatFiles() )
      Data[0].Write( g, "data" );
    // Only save random numbers if asked, and then only for bootstrap (jackknife don't need random)
    if( RandomCache::SaveFatFiles() && Seed() != SeedWildcard && SampleSize )
    {
      const fint &cSeed{ Seed() };
      Matrix<fint> mRnd{ RandomCache::Global.Get( cSeed, SampleSize, Data[0].extent() ) };
      Dims[0] = mRnd.size1;
      Dims[1] = mRnd.size2;
      dsp = ::H5::DataSpace( 2, Dims );
      // The values in this table go from 0 ... NumSamples_ - 1, so choose a space-minimising size in the file
      ::H5::DataSet ds = g.createDataSet( RandomCache::sRandom, SampleSize<=static_cast<int>(std::numeric_limits<std::uint8_t>::max())
        ? ::H5::PredType::STD_U8LE : SampleSize<=static_cast<int>(std::numeric_limits<std::uint16_t>::max())
        ? ::H5::PredType::STD_U16LE: ::H5::PredType::STD_U32LE, dsp );
      ds.write( mRnd.Data(), H5::Equiv<std::uint_fast32_t>::Type );
      ds.close();
      dsp.close();
    }
    if( !Raw.Empty() )
      H5::WriteMatrix( g, "samples_Raw", Raw );
    if( !Binned[0].Empty() )
      H5::WriteMatrix( g, "samples_Binned", Binned[0] );
    if( SummaryNames.size() )
    {
      Dims[0] = SummaryNames.size();
      Dims[1] = Nt_;
      dsp = ::H5::DataSpace( 2, Dims );
      ::H5::DataSet ds = g.createDataSet( sSummaryDSName, H5::Equiv<ValWithEr<scalar_type>>::Type, dsp );
      ds.write( &SummaryData(0,0), H5::Equiv<ValWithEr<scalar_type>>::Type );
      ds.close();
      dsp.close();
    }
    g = f.openGroup( "/" );
    static const std::string sGridDatasetThreshold{ "_Grid_dataset_threshold" };
    if( H5::ExistsAttribute( g, sGridDatasetThreshold ) )
    {
      int PrevNumAttributes;
      a = g.openAttribute( sGridDatasetThreshold );
      a.read( ::H5::PredType::NATIVE_INT, &PrevNumAttributes );
      if( NumAttributes > PrevNumAttributes )
        a.write( ::H5::PredType::NATIVE_INT, &NumAttributes );
    }
    else
    {
      a = g.createAttribute( sGridDatasetThreshold, ::H5::PredType::STD_U32LE, ds1);
      a.write( ::H5::PredType::NATIVE_INT, &NumAttributes );
    }
    a.close();
    ds1.close();
    g.close();
    bOK = true;
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
  }
  if( !bOK )
    throw std::runtime_error( "Unable to write sample to "+Name_.Filename+", group "+GroupName );
}

// Make a summary of the data
// If a name is specified, simply save the central replica
// Otherwise, calculate exp and cosh mass as well

template <typename T>
void Sample<T>::MakeCorrSummary()
{
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  // Now perform summaries
  AllocSummaryBuffer();
  const int tMid{ bFolded() ? Nt_ : Nt_ / 2 };
  JackBoot<scalar_type> Buffer( NumSamples(), Nt_ );
  std::vector<ValEr> veBuf( Nt_ );
  const int NumRealFields{ static_cast<int>( SummaryNames.size() ) / scalar_count };
  for( int i = 0; i < Traits::scalar_count; i++ ) // real or imaginary
  {
    for( int f = 0; f < std::min( NumRealFields, 4 ); ++f ) // each field
    {
      const int SummaryRow{ f + i * NumRealFields };
      // Fill in the JackBoot sample for this field
      if( f == 1 )
      {
        // Same as previous row - except average is taken from average of all samples
        JackBoot<scalar_type>::MakeMean( Buffer.GetCentral(), Buffer.Replica );
      }
      else
      {
        for( int n = idxCentral; n < NumSamples(); ++n )
        {
          for( int t = 0; t < Nt_; ++t )
          {
            // Index of previous and next timeslice relative to current timeslice
            const int tNext{ ( t == Nt_ - 1 ? 0       : t + 1 ) };
            const int tPrev{ ( t == 0       ? Nt_ - 1 : t - 1 ) };
            const scalar_type z{ Traits::RealImag( (*this)(n,t), i ) };
            const scalar_type zNext{ Traits::RealImag( (*this)(n,tNext), i ) };
            const scalar_type zPrev{ Traits::RealImag( (*this)(n,tPrev), i ) };
            scalar_type d;
            switch( f )
            {
              case 0: // central value
                d = z;
                break;
              case 2: // exponential mass
                if( t == 0 )
                  d = NaN;
                else if( t <= tMid )
                  d = std::log( zPrev / z );
                else
                  d = std::log( zNext / z );
                break;
              case 3: // cosh mass
                d = std::acosh( ( zPrev + zNext ) / ( z * 2 ) );
                break;
            }
            Buffer(n,t) = d;
          }
        }
      }
      // Now get the averages
      Buffer.MakeStatistics( veBuf );
      for( int t = 0; t < Nt_; ++t )
        SummaryData( SummaryRow, t ) = veBuf[t];
    }
  }
}

template <typename T>
void Sample<T>::WriteSummary( const std::string &sOutFileName, bool bVerboseSummary )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  if( SummaryNames.empty() )
    throw std::runtime_error( "Summaries can't be written because they've not been created" );
  std::ofstream s( sOutFileName );
  SummaryHeader<T>( s, sOutFileName );
  SummaryComments( s, bVerboseSummary );
  SummaryColumnNames( s );
  s << NewLine;
  SummaryContents( s );
  s << NewLine;
}

template <typename T>
void Sample<T>::SummaryComments( std::ostream & s, bool bVerboseSummary,
                                 bool bShowColumnNames ) const
{
  s << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
  if( !Ensemble.empty() )
    s << "# Ensemble: " << Ensemble << NewLine;
  s << "# Seed: " << SeedString() << NewLine;
  if( !SeedMachine_.empty() ) s << "# Seed machine: " << SeedMachine_ << NewLine;
  s << "# Bootstrap: " << NumSamples() << NewLine << "# SampleSize: " << SampleSize << NewLine;
  if( binSize != 1 ) s << "# BinSize: " << binSize << NewLine;
  if( !ConfigCount.empty() )
  {
    s << "# Configs: " << ConfigCount.size();
    for( const MLU::ConfigCount &cc : ConfigCount )
      s << CommaSpace << cc;
    s << NewLine;
  }
  s << "# NumColumns: " << ( ColumnNames.size() ? ColumnNames.size() : Nt() );
  if( ColumnNames.empty() )
    s << " (rows t0 ... t" << ( Nt() - 1 ) << ")";
  else if( bShowColumnNames )
  {
    s << NewLine << "# " << sColumnNames << ": ";
    for( std::size_t i = 0; i < ColumnNames.size(); ++i )
    {
      if( i )
        s << CommaSpace;
      s << ColumnNames[i];
    }
  }
  s << NewLine;
  if( bVerboseSummary )
  {
    s << "# FileCount: " << FileList.size() << NewLine;
    std::size_t i = 0;
    for( const std::string &f : FileList )
      s << "# File " << i++ << ": " << f << NewLine;
  }
}

template <typename T>
void Sample<T>::SummaryColumnNames( std::ostream &os ) const
{
  if( ColumnNames.empty() )
  {
    // First column is timeslice, then each of the summary columns per timeslice
    os << "t";
    for( const std::string &Name : SummaryNames )
    {
      os << Space;
      ValWithEr<T>::Header( Name, os, Space );
    }
  }
  else
  {
    // Show the bootstrap central replica for every column on one row
    bool bFirst{ true };
    for( const std::string &Name : ColumnNames )
    {
      if( bFirst )
        bFirst = false;
      else
        os << Space;
      ValWithEr<T>::Header( Name, os, Space );
    }
  }
}

template <typename T>
void Sample<T>::SummaryContents( std::ostream &os ) const
{
  os << std::setprecision(std::numeric_limits<scalar_type>::digits10+2) << std::boolalpha;
  if( ColumnNames.empty() )
  {
    // One row per timeslice, showing each of the summary values
    for( int t = 0; t < Nt_; ++t )
    {
      os << t;
      for( int f = 0; f < SummaryNames.size(); f++ )
        os << Space << SummaryData( f, t );
      if( t != Nt_ - 1 )
        os << NewLine;
    }
  }
  else
  {
    // A single row, showing the first summary field for each column
    for( int t = 0; t < Nt_; t++ )
    {
      if( t )
        os << Space;
      os << SummaryData( t );
    }
  }
}

template class Sample<double>;
template class Sample<float>;
template class Sample<std::complex<double>>;
template class Sample<std::complex<float>>;

END_MLU_NAMESPACE
