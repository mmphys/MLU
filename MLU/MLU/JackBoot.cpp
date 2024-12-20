/**
 
 Jackknife and Bootstrap support

 Source file: JackBoot.cpp
 
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
#include "JackBoot.hpp"
#include <random>

// Uncomment this next line to debug without OpenMP
//#define DEBUG_DISABLE_OMP
#ifndef DEBUG_DISABLE_OMP
#include <omp.h>
#endif

BEGIN_MLU_NAMESPACE

const std::string s_C{ "_C" };
const std::string JackBootBase::s_S{ "_S" };
const std::string JackBootBase::sBootstrap{ "bootstrap" };
const std::string JackBootBase::sJackknife{ "jackknife" };

/*****************************************************************************

 This is a cache of random numbers used to construct bootstrap samples
 The number of replicas in each sample grows to support the largest request ever made

 *****************************************************************************/

const std::string RandomCache::sSeed{ "Seed" };
const std::string RandomCache::sSeedMachine{ "SeedMachine" };
const std::string RandomCache::sRandom{ "Random" };
bool        RandomCache::CacheDirPrefixOK{};
std::string RandomCache::CacheDirPrefix{ sRandom + "/" };
bool        RandomCache::DefaultSeedOK{};
SeedType    RandomCache::DefaultSeed_{ 1835672416u }; // Mike's default seed
bool        RandomCache::DefaultNumReplicasOK{};
std::size_t RandomCache::DefaultNumReplicas_{ 10000u }; // Mike's default number of replicas

RandomCache RandomCache::Global;

/// Get default seed from envionment variable MLUSeed (only once then cache it)
SeedType RandomCache::DefaultSeed()
{
  if( !DefaultSeedOK )
  {
    DefaultSeedOK = true;
    const char * pDefaultString{ std::getenv( "MLUSeed" ) };
    if( pDefaultString && * pDefaultString )
    {
      if( EqualIgnoreCase( sRandom, pDefaultString ) )
        DefaultSeed_ = std::random_device{}(); // Hardware generated random number
      else
        DefaultSeed( pDefaultString );
    }
  }
  return DefaultSeed_;
}

SeedType RandomCache::DefaultSeed( SeedType NewDefaultSeed )
{
  DefaultSeed_ = NewDefaultSeed;
  DefaultSeedOK = true;
  return DefaultSeed_;
}

void RandomCache::DefaultSeed( const std::string &SeedOrFile_ )
{
  std::string SeedOrFile{ SeedOrFile_ };
  bool bGotSeed{ true };
  SeedType MySeed;
  try // to interpret the switch as the random number
  {
    MySeed = Seed( SeedOrFile );
  }
  catch(const std::exception &e)
  {
    bGotSeed = false;
  }
  if( !bGotSeed )
  {
    Matrix<fint> Random;
    std::string Machine;
    const bool LoadFromOutsideCache{ SeedOrFile[0] == '/' };
    SeedOrFile = PrependCachePath( SeedOrFile );
    bGotSeed = Read( SeedOrFile, Random, MySeed, Machine );
    if( bGotSeed )
    {
      SetHostName( Machine );
      Global.Put( MySeed, Random, LoadFromOutsideCache ); // Only write if loaded from outside cache
      std::cout << "Seed " << MySeed << " generated on " << Machine << "\n";
    }
  }
  if( !bGotSeed )
    throw std::runtime_error( "Seed neither number nor file: " + SeedOrFile );
  DefaultSeed( MySeed );
}

std::size_t RandomCache::DefaultNumReplicas()
{
  if( !DefaultNumReplicasOK )
  {
    const char * pDefaultString{ std::getenv( "MLUNumBoot" ) };
    if( pDefaultString && * pDefaultString )
    {
      try
      {
        DefaultNumReplicas_ = FromString<std::size_t>( pDefaultString );
        if( DefaultNumReplicas_ < 1 )
          DefaultNumReplicas_ = 1;
      }
      catch(...){}
    }
    DefaultNumReplicasOK = true;
  }
  return DefaultNumReplicas_;
}

std::size_t RandomCache::DefaultNumReplicas( std::size_t NewDefaultNumReplicas )
{
  DefaultNumReplicas_ = NewDefaultNumReplicas;
  DefaultNumReplicasOK = true;
  return DefaultNumReplicas_;
}

bool RandomCache::Key::Less::operator()( const Key &lhs, const Key &rhs ) const
{
  if( lhs.Seed != rhs.Seed )
    return lhs.Seed < rhs.Seed;
  return lhs.NumSamples < rhs.NumSamples;
}

void RandomCache::Write( ::H5::Group &g, const Matrix<fint> &Random,
                         SeedType Seed, const std::string &Machine )
{
  hsize_t Dim;
  Dim = 1;
  ::H5::DataSpace ds1( 1, &Dim );
  const ::H5::DataType &dtRandom{ Random.size2 <= 0x100u ? ::H5::PredType::STD_U8LE
    : Random.size2 <= 0x10000u ? ::H5::PredType::STD_U16LE : ::H5::PredType::STD_U32LE };
  H5::WriteMatrix( g, RandomCache::sRandom, Random, dtRandom );
  ::H5::Attribute a = g.createAttribute( RandomCache::sSeed, ::H5::PredType::STD_U32LE, ds1 );
  a.write( H5::Equiv<SeedType>::Type, &Seed );
  a.close();
  a = g.createAttribute( RandomCache::sSeedMachine, H5::Equiv<std::string>::Type, ds1 );
  a.write( H5::Equiv<std::string>::Type, Machine );
  a.close();
}

void RandomCache::Read( ::H5::Group &g, Matrix<fint> &Random,
                        SeedType &Seed, std::string &Machine )
{
  ::H5::Attribute a = g.openAttribute( RandomCache::sSeed );
  a.read( H5::Equiv<SeedType>::Type, &Seed );
  a.close();
  a = g.openAttribute( RandomCache::sSeedMachine );
  a.read( a.getStrType(), Machine );
  a.close();
  H5::ReadMatrix( g, sRandom, Random );
}

bool RandomCache::Read( const std::string &Filename,
                        Matrix<fint> &Random, SeedType &Seed, std::string &Machine )
{
  bool bOK{};
  // Try to get the seed from HDF5 file
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  try
  {
    ::H5::H5File f;
    ::H5::Group  g;
    H5::OpenFileGroup( f, g, Filename, "Random cache load: ", nullptr );
    Read( g, Random, Seed, Machine );
    bOK = true;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  catch(...)
  {
    ::H5::Exception::setAutoPrint(h5at, f5at_p);
    throw;
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  return bOK;
}

void RandomCache::SetCachePrefix( std::string &NewCachePrefix )
{
  CacheDirPrefix = NewCachePrefix;
  CacheDirPrefixOK = true;
}

// True if random numbers should be saved with data files
bool RandomCache::SaveFatFiles()
{
  static bool bInitialised{};
  static bool bFat{};
  if( !bInitialised )
  {
    bInitialised = true;
    if( std::getenv( "MLUFat" ) ) // I.e. set to anything (including empty)
      bFat = true;
  }
  return bFat || GetCachePrefix().empty();
}

const std::string &RandomCache::GetCachePrefix()
{
  if( !CacheDirPrefixOK )
  {
    CacheDirPrefixOK = true;
    const char * pDefaultString{ std::getenv( "MLUCache" ) };
    if( pDefaultString )
    {
      CacheDirPrefix = pDefaultString;
      std::cout << "MLUCache ";
      if( CacheDirPrefix.empty() )
        std::cout << "disabled via MLUCache environment variable";
      else
      {
        std::cout << CacheDirPrefix;
        if( CacheDirPrefix.back() != '/' )
          std::cout << " Warning: doesn't end with /";
      }
      std::cout << "\n";
    }
  }
  return CacheDirPrefix;
}

std::string RandomCache::PrependCachePath( const std::string &sFilename )
{
  std::string s;
  if( !sFilename.empty() && sFilename[0] != '/' )
    s.append( GetCachePrefix() );
  s.append( sFilename );
  return s;
}

std::string RandomCache::GetCachePath( const Key &key )
{
  std::string Filename{ GetCachePrefix() };
  if( !Filename.empty() )
  {
    Filename.append( GetHostName() );
    Filename.append( 1, '.' );
    Filename.append( std::to_string( key.NumSamples ) );
    Filename.append( 1, '.' );
    Filename.append( sRandom );
    Filename.append( 1, '.' );
    Filename.append( SeedString( key.Seed ) );
    Filename.append( "." DEF_FMT );
  }
  return Filename;
}

Matrix<fint> RandomCache::Make( fint &Seed, std::size_t NumSamples, std::size_t NumReplicas )
{
  if( Seed == SeedWildcard )
    throw std::runtime_error( "RandomCache::Make Seed==SeedWildcard (" + std::to_string(Seed) + ")" );
  Matrix<fint> m( NumReplicas, NumSamples );
  std::mt19937                        engine( Seed );
  std::uniform_int_distribution<fint> random( 0, static_cast<fint>( NumSamples - 1 ) );
  for( std::size_t i = 0; i < NumReplicas; ++i )
    for( std::size_t j = 0; j < NumSamples; ++j )
      m(i,j) = random( engine );
  return m;
}

void RandomCache::Compare( const Matrix<fint> &l, const Matrix<fint> &r )
{
  if( !l.Empty() && !r.Empty() && l.size2 != r.size2 )
  {
    std::ostringstream os;
    os << "Random mismatch number of samples " << l.size2 << " != " << r.size2;
    throw std::runtime_error( os.str().c_str() );
  }
  const std::size_t CompareReps{ std::min( l.size1, r.size1 ) };
  for( std::size_t i = 0 ; i < CompareReps; ++i )
    for( std::size_t j = 0 ; j < l.size2; ++j )
      if( l(i,j) != r(i,j) )
      {
        std::ostringstream os;
        os << "Seed mismatch (" << i << ',' << j << ") " << l(i,j) << " != " << r(i,j);
        throw std::runtime_error( os.str().c_str() );
      }
}

void RandomCache::SaveRandom( const Key &key, Value &value )
{
  std::string Filename{ GetCachePath( key ) };
  if( !Filename.empty() )
  {
    bool bOK{};
    try // to write in my format
    {
      if( MLU::FileExists( Filename ) )
      {
        // See whether compatible with existing file ... but only do this once
        if( !value.bOnDisk )
        {
          std::string Machine;
          SeedType LoadSeed;
          Matrix<fint> m;
          if( !Read( Filename, m, LoadSeed, Machine ) )
            throw std::runtime_error( "Random numbers incompatible with cache " + Filename );
          Compare( m, value.m );
          value.bOnDisk = true;
        }
        bOK = true;
      }
      else
      {
        // Didn't exist - write it
        MakeAncestorDirs( Filename );
        ::H5::H5File f( Filename, H5F_ACC_TRUNC );
        ::H5::Group g = f.createGroup( sRandom );
        Write( g, value.m, key.Seed, GetHostName() );
        bOK = true;
      }
    }
    catch(const ::H5::Exception &) {}
    catch(const std::exception &e) {}
    if( !bOK )
      std::cout << "Warning: Unable to write random numbers to " << Filename << "\n";
  }
}

Matrix<fint> RandomCache::Get( fint &Seed, std::size_t NumSamples, std::size_t NumReplicas )
{
  if( Seed == SeedWildcard )
    Seed = DefaultSeed();
  const Key key{ Seed, NumSamples };
  Value &value{ Map[key] };
  Matrix<fint> &m{ value.m };
  assert(((m.size1==0 && m.size2==0)||m.size2==NumSamples)&&"RandomCache::Get() NumSamples mismatch");
  if( m.size1 == 0 )
  {
    // First time I've seen this key - can I load it from cache
    std::string Machine;
    SeedType LoadSeed;
    const std::string Filename{ GetCachePath( key ) };
    value.bOnDisk = Read( Filename, m, LoadSeed, Machine );
  }
  if( m.size1 < NumReplicas )
  {
    // Not enough replicas - make some more
    Matrix<fint> New{ Make( Seed, NumSamples, NumReplicas ) };
    Compare( m, New );
    m = std::move( New );
    SaveRandom( key, value );
  }
  Matrix<fint> Copy;
  Copy.MapView( m );
  return Copy;
}

Matrix<fint> RandomCache::Get( const fint &Seed, std::size_t NumSamples, std::size_t NumReplicas )
{
  if( Seed == SeedWildcard )
    throw std::runtime_error( "Invalid Seed " + std::to_string( Seed ) );
  SeedType SeedCopy{ Seed };
  return Get( SeedCopy, NumSamples, NumReplicas );
}

void RandomCache::Put( fint Seed, Matrix<fint> &New, bool bWrite )
{
  if( Seed == SeedWildcard )
    throw std::runtime_error( "Invalid Seed " + SeedString( Seed ) );
  if( New.Empty() )
    throw std::runtime_error( "Can't put empty random numbers in cache" );
  // Are the random numbers valid?
  const std::size_t NumSamples{ New.size2 };
  for( std::size_t i = 0 ; i < New.size1; ++i )
    for( std::size_t j = 0 ; j < NumSamples; ++j )
      if( New(i,j) >= NumSamples )
      {
        std::ostringstream os;
        os << "NewCache::Put() Seed " << SeedString( Seed )
           << ", replica " << i << ", sample " << j
           << " = " << New(i,j) << " invalid. Should be < " << NumSamples << " samples.";
        throw std::runtime_error( os.str().c_str() );
      }
  // Do we already have this in the cache?
  Key key{ Seed, NumSamples };
  typename MapT::iterator it{ Map.find( key ) };
  if( it == Map.end() )
  {
    // Nothing in cache - move the random numbers into the cache
    Value InsertVal{ std::move( New ), false };
    auto pair{ Map.emplace( key, std::move( InsertVal ) ) };
    assert( pair.second && "Bug: I just searched for this, but it wasn't inserted" );
    // Make the caller's matrix point back to the new entry in the cache
    it = pair.first;
    Value &value{ it->second };
    New.MapView( value.m );
    if( bWrite )
      SaveRandom( key, value );
  }
  else
  {
    // Something in cache - are the numbers compatible?
    Value &value{ it->second };
    Matrix<fint> &Old{ value.m };
    assert( Old.size1 && Old.size2 == NumSamples && "Invalid random numbers in cache" );
    Compare( Old, New );
    if( New.size1 > Old.size1  )
    {
      Old = std::move( New ); // New has more replicas - move this into cache
      if( bWrite )
        SaveRandom( key, value );
    }
    New.MapView( Old ); // Make callers matrix point to cache
  }
}

/*****************************************************************************
 
 true use the central replica when computing (co)variance
 false use mean of replicas when computing (co)variance

 false unless MLUCovarCentral environment variable set

 *****************************************************************************/

bool JackBootBase::UseCentralCovar()
{
  static bool bUseCentral{};
  static bool bInitialised{};
  if( !bInitialised )
  {
    bInitialised = true;
    bUseCentral = std::getenv( "MLUCovarCentral" ); // true if set to anything
  }
  return bUseCentral;
}

/*****************************************************************************
 
 Get random numbers appropriate for given parameters
 NumReplicas should be the total number of replicas required (so that random numbers created once for all replicas)
 
 *****************************************************************************/

void JackBootBase::GetRandom( Vector<fint> &vRandom, fint Seed,
                              std::size_t NumSamples, std::size_t NumReplicas, std::size_t Replica )
{
  if( NumSamples == 1 || ( NumSamples == 2 && Seed == SeedWildcard ) )
  {
    std::ostringstream os;
    os << "JackBootBase::GetRandom() Cannot resample " << NumSamples << " samples";
    throw std::runtime_error( os.str().c_str() );
  }
  if( Seed == SeedWildcard )
    NumReplicas = NumSamples;
  if( Replica == idxCentral || Replica == idxReplicaMean )
  {
    vRandom.resize( NumSamples );
    for( fint i = 0; i < NumSamples; ++i )
      vRandom[i] = i;
  }
  else if( Replica >= NumReplicas )
  {
    std::ostringstream ss;
    ss << "JackBoot<T>::GetRandom Request for replica " << Replica << " of " << NumReplicas
       << " replicas";
    throw std::runtime_error( ss.str().c_str() );
  }
  // The presence or absence of random numbers tells us whether we are doing bootstrap or jackknife
  else if( Seed == SeedWildcard )
  {
    vRandom.resize( NumSamples - 1 );
    for( fint i = 0, r = 0; r < NumSamples; ++r )
      if( r != Replica )
        vRandom[i++] = r;
  }
  else
  {
    // Bootstrap on non-central replica
    Matrix<fint> Random{ RandomCache::Global.Get( Seed, NumSamples, NumReplicas ) };
    // Random numbers come from pre-drawn random numbers for this replica
    vRandom.MapRow( Random, Replica );
  }
}

std::ostream &operator<<( std::ostream &os, JackBootBase::Norm norm )
{
  switch( norm )
  {
    case JackBootBase::Norm::RawBinned:
      os << "RawBinned";
      break;
    case JackBootBase::Norm::Bootstrap:
      os << "Bootstrap";
      break;
    case JackBootBase::Norm::Jackknife:
      os << "Jackknife";
      break;
    default:
      os << "SampleCovar";
      break;
  }
  return os;
}

template <typename T>
void JackBoot<T>::clear()
{
  Central.clear();
  Replica.clear();
  ReplicaMean.clear();
  Seed = RandomCache::DefaultSeed();
}

template <typename T>
typename JackBoot<T>::Real JackBoot<T>::GetNorm( Norm norm, std::size_t NumReplicas )
{
  if( NumReplicas == 0 )
  {
    std::ostringstream ss;
    ss << "JackBoot::GetNorm() Norm undefined for " << norm
       << " with " << NumReplicas << " replicas";
    throw std::runtime_error( ss.str().c_str() );
  }
  double A;
  switch( norm )
  {
    case Norm::Bootstrap:
      // Normalisation for bootstrap (error of mean)
      A = 1.0 / NumReplicas;
      break;
    case Norm::Jackknife:
      // Normalisation for error of mean for Jackknife
      A = static_cast<double>( NumReplicas - 1 ) / NumReplicas;
      break;
    case Norm::RawBinned:
      // Normalisation for error of mean
      A = 1.0 / ( NumReplicas * ( NumReplicas - 1 ) );
      break;
    default:
      // Normalisation for variance of the data NOT error of mean
      A = 1.0 / ( NumReplicas - 1 );
      break;
  }
  return static_cast<Real>( A );
}

/*****************************************************************************

 Make the mean of any matrix
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::MakeMean( Vector<T> &vMean, const MatrixView<T> &Source )
{
  if( Source.size1() == 0 || Source.size2() == 0 )
    throw std::runtime_error( "JackBoot::MakeMean Source matrix empty" );
  std::size_t NumReplicas{ Source.size1() };
  std::size_t Extent{ Source.size2() };
  vMean.resize( Extent );
  vMean = static_cast<T>( 0 );
  for( std::size_t i = 0; i < NumReplicas; ++i )
    for( std::size_t j = 0; j < Extent; ++j )
      vMean[j] += Source(i,j);
  if( NumReplicas > 1 )
    vMean *= static_cast<Real>( 1 ) / NumReplicas;
}

/*****************************************************************************

 Make the mean of a matrix while resampling it - e.g. from raw/binned source data
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::MakeMean( Vector<T> &vMean, Vector<fint> &vRandom, const MatrixView<T> &Source,
                            std::size_t Replica, std::size_t NumReplicas, SeedType Seed )
{
  if( Source.size1() == 0 || Source.size2() == 0 )
    throw std::runtime_error( "JackBoot::MakeMean Source matrix empty" );
  std::size_t NumSamples{ Source.size1() };
  std::size_t Extent{ Source.size2() };
  GetRandom( vRandom, Seed, NumSamples, NumReplicas, Replica );
  vMean.resize( Extent );
  vMean = static_cast<T>( 0 );
  for( std::size_t i = 0; i < vRandom.size; ++i )
    for( std::size_t j = 0; j < Extent; ++j )
      vMean[j] += Source(vRandom[i],j);
  if( vRandom.size > 1 )
    vMean *= static_cast<Real>( 1 ) / vRandom.size;
}

/*****************************************************************************

 Perform jackknife or bootstrap (based on seed) resampling
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::Resample( const MatrixView<T> &Source, std::size_t NumReplicas )
{
  std::size_t Extent{ Source.size2() };
  std::size_t NumSamples{ Source.size1() };
  if( Extent == 0 || NumSamples == 0 )
    throw std::runtime_error( "JackBoot::Resample Source matrix empty" );
  if( Seed == SeedWildcard )
  {
    if( NumSamples < 2 )
      throw std::runtime_error( "JackBoot::Resample Jackknife needs minimum of 2 samples" );
    NumReplicas = NumSamples;
  }
  else if( !NumReplicas )
    throw std::runtime_error( "JackBoot::Resample At least one bootstrap replica required" );
  // Now compute each replica
  MakeMean( Central, Source );
  Replica.resize( NumReplicas, Extent );
  if( Seed == SeedWildcard )
  {
    // Jackknife resample
    const Real AvgNorm{ static_cast<Real>( 1. ) / ( NumSamples - 1u ) };
    Vector<T> nxBar{ Central };
    nxBar *= static_cast<Real>( NumSamples );
#ifndef DEBUG_DISABLE_OMP
    #pragma omp parallel for default(shared) schedule(dynamic)
#endif
    for( std::size_t idx = 0; idx < NumReplicas; ++idx )
      for( std::size_t i = 0; i < Extent; ++i )
        Replica( idx, i ) = ( nxBar[i] - Source( idx, i ) ) * AvgNorm;
  }
  else
  {
    // Bootstrap resample
    const Real AvgNorm{ static_cast<Real>( 1. ) / NumSamples };
    Matrix<fint> Random{ RandomCache::Global.Get( Seed, NumSamples, NumReplicas ) };
#ifndef DEBUG_DISABLE_OMP
    #pragma omp parallel default(shared)
#endif
    {
      VectorView<T> Sample;
#ifndef DEBUG_DISABLE_OMP
    #pragma omp for schedule(dynamic)
#endif
    for( std::size_t idx = 0; idx < NumReplicas; ++idx )
    {
      for( std::size_t Resample = 0; Resample < NumSamples; ++Resample )
      {
        Sample.MapRow( Source, Random( idx, Resample ) );
        for( std::size_t i = 0; i < Extent; ++i )
        {
          if( Resample == 0 )
            Replica( idx, i ) = Sample[i];
          else
            Replica( idx, i ) += Sample[i];
        }
      }
      // Turn the sum into an average
      for( std::size_t i = 0; i < Extent; ++i )
        Replica( idx, i ) *= AvgNorm;
    }
    }
  }
  // Take mean of replicas. Jaccknife should be same as Central, but compute to find errors
  MakeMean( ReplicaMean, Replica );
}

/*****************************************************************************

 Make statistics for this sample

 Needs to work even when empty
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::MakeStatistics( ValWithEr<T> *vStats, std::size_t StartCol,
                                  std::size_t NumCols ) const
{
  if( NumCols )
  {
    if( StartCol + NumCols > extent() )
    {
      std::ostringstream os;
      os << "JackBoot<T>::MakeStatistics() StartCol " << StartCol << " + NumCols " << NumCols
         << " > extent " << extent();
      throw std::runtime_error( os.str().c_str() );
    }

    // Default to showing the central values
    for( std::size_t i = 0; i < NumCols; ++i )
    {
      vStats[i] = Central[i + StartCol];
      vStats[i].Check = NumReplicas() ? 1 : 0;
    }
    if( NumReplicas() > 1 )
    {
      if( Seed == SeedWildcard )
      {
        // Jackknife statistics
        std::vector<std::size_t> Count( NumCols, 0 );
        Vector<T> Var( NumCols );
        Var = static_cast<T>( 0 );
        // Sum the Variance entries on each replica
        for( std::size_t replica = 0; replica < NumReplicas(); ++replica )
        {
          for( std::size_t j = 0; j < NumCols; ++j )
          {
            if( ::MLU::IsFinite( Replica( replica, j + StartCol ) ) )
            {
              Var[j] += Squared( Replica( replica, j + StartCol ) - Central[j + StartCol] );
              ++Count[j];
            }
          }
        }
        // Normalise
        Vector<T> Column;
        for( std::size_t j = 0; j < NumCols; ++j )
        {
          if( Count[j] )
          {
            Var[j] *= GetNorm( Norm::Jackknife, Count[j] );
            Var[j] = std::sqrt( Var[j] );
            vStats[j].Low -= Var[j];
            vStats[j].High += Var[j];
            Column.MapColumn( const_cast<Matrix<T> &>( Replica ), j + StartCol );
            Column.MinMax( vStats[j].Min, vStats[j].Max, true );
          }
          vStats[j].Check = static_cast<Real>( Count[j] ) / NumReplicas();
        }
      }
      else
      {
        // Bootstrap statistics
        std::vector<T> Buffer( NumReplicas() );
        for( std::size_t j = 0; j < NumCols; ++j )
        {
          std::size_t Count{};
          for( std::size_t replica = 0; replica < NumReplicas(); ++replica )
          {
            if( ::MLU::IsFinite( Replica( replica, j + StartCol ) ) )
              Buffer[Count++] = Replica( replica, j + StartCol );
          }
          vStats[j].Get( Central[j + StartCol], Buffer, Count );
        }
      }
    }
  }
}

/*****************************************************************************
 
 Make covariance matrix or variance vector given any matrix (raw/binned, bootstrap or jackknife) and a vector of central values.
 No resampling is possible (i.e. central replica only).
 Caller specifies normalisation required, most of which compute `\Sigma_{\bar{\vb{x}}}`, i.e. error of the mean
 but Norm::SampleCovar can be used to get the sample covariance if required.
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::MakeCovar( Matrix<T> &Covar, const VectorView<T> &Mean,
                             const MatrixView<T> &Source, Norm norm )
{
  CheckVarCovar( "static JackBoot::MakeCovar()", Mean, Source );
  const Real A{ GetNorm( norm, Source.size1() ) };
  MakeCovarHelper( Covar, Source, Mean, A );
}

template <typename T>
void JackBoot<T>::MakeVar( Vector<T> &Var, const VectorView<T> &Mean,
                           const MatrixView<T> &Source, Norm norm )
{
  CheckVarCovar( "static JackBoot::MakeVar()", Mean, Source );
  const Real A{ GetNorm( norm, Source.size1() ) };
  MakeVarHelper( Var, Source, Mean, A );
}

/*****************************************************************************
 
 Make the (co)variance of a matrix while resampling it - e.g. from raw/binned source data

*****************************************************************************/

template <typename T>
void JackBoot<T>::MakeCovar( Matrix<T> &Covar, const MatrixView<T> &Source,
                             std::size_t Replica, std::size_t NumReplicas, SeedType Seed )
{
  Vector<T> Mean;
  Vector<fint> Random;
  MakeMean( Mean, Random, Source, Replica, NumReplicas, Seed );
  const std::size_t Extent{ Source.size2() };
  const std::size_t NumSamples{ Source.size1() };
  // Clear the covariance matrix
  Covar.resize( Extent, Extent );
  Covar = static_cast<T>( 0 );
  VectorView<T> Data;         // This is a view into each replica
  Vector<T> Error( Extent );  // This is the error on each replica
  // Sum the covariance entries on each replica - lower left triangle
  for( std::size_t replica = 0; replica < Random.size; ++replica )
  {
    Data.MapRow( Source, Random[replica] );
    for( std::size_t i = 0; i < Extent; ++i )
    {
      Error[i] = Data[i] - Mean[i];
      for( std::size_t j = 0; j <= i; ++j )
        Covar( i, j ) += Conjugate( Error[i] ) * Error[j];
    }
  }
  // Normalise and fill in upper right triangle
  const Real norm{ GetNorm( Norm::RawBinned, NumSamples ) };
  for( std::size_t i = 0; i < Extent; ++i )
    for( std::size_t j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

template <typename T>
void JackBoot<T>::MakeVar( Vector<T> &Var, const MatrixView<T> &Source,
                           std::size_t Replica, std::size_t NumReplicas, SeedType Seed )
{
  Vector<T> Mean;
  Vector<fint> Random;
  MakeMean( Mean, Random, Source, Replica, NumReplicas, Seed );
  const std::size_t Extent{ Source.size2() };
  const std::size_t NumSamples{ Source.size1() };
  // Clear the covariance matrix
  Var.resize( Extent );
  Var = static_cast<T>( 0 );
  VectorView<T> Data;         // This is a view into each replica
  // Sum the covariance entries on each replica - lower left triangle
  for( std::size_t replica = 0; replica < Random.size; ++replica )
  {
    Data.MapRow( Source, Random[replica] );
    for( std::size_t i = 0; i < Extent; ++i )
      Var[i] += Squared( Data[i] - Mean[i] );
  }
  // Normalise and fill in upper right triangle
  const Real norm{ GetNorm( Norm::RawBinned, NumSamples ) };
  for( std::size_t i = 0; i < Extent; ++i )
    Var[i] *= norm;
}

/*****************************************************************************
 
 Make covariance matrix or variance vector given a bootstrap or jackknife and the source matrix it was created from.
 Resampling is possible for any replica, with the (co)variance coming from the source data, resampled appropriately.
 Normalisation is always `Norm::RawBinned`, i.e. `\Sigma_{\bar{\vb{x}}}`error of the mean
 
 *****************************************************************************/

/*template <typename T>
void JackBoot<T>::MakeCovar( Matrix<T> &Covar, const MatrixView<T> &Source, std::size_t idx ) const
{
  VectorView<T> Mean;
  Vector<fint> Random;
  const Real A{ CheckVarCovar( "JackBoot::MakeCovar()", Mean, Random, Source, idx ) };
  MakeCovarHelper( Covar, Source, Mean, A, &Random );
}

template <typename T>
void JackBoot<T>::MakeVar( Vector<T> &Var, const MatrixView<T> &Source, std::size_t idx ) const
{
  VectorView<T> Mean;
  Vector<fint> Random;
  const Real A{ CheckVarCovar( "JackBoot::MakeVar()", Mean, Random, Source, idx ) };
  MakeVarHelper( Var, Source, Mean, A, &Random );
}*/

/*****************************************************************************
 
 Common validation for covariance and variance

 ****************************************************************************/

/// Common validation for (co)variance given any matrix and a vector of central values
template <typename T>
void JackBoot<T>::CheckVarCovar( const std::string &sPrefix, const VectorView<T> &vMean,
                                 const MatrixView<T> &Source )
{
  // Error check
  const std::size_t Extent{ Source.size2() };
  if( Extent == 0 )
    throw std::runtime_error( sPrefix + " empty extent" );
  if( vMean.size() != Extent )
    throw std::runtime_error( sPrefix + " mismatch: " + std::to_string( vMean.size() )
                             + " central values != extent " + std::to_string( vMean.size() ) );
}

/// Common validation for (co)variance given a bootstrap / jackknife and the matrix it was made from
template <typename T> typename JackBoot<T>::Real
JackBoot<T>::CheckVarCovar( const std::string &sPrefix, VectorView<T> &Mean, Vector<fint> &vRandom,
                            const MatrixView<T> &Source, std::size_t idx ) const
{
  // TODO: Check whether using bootstrap mean is root cause of why my correlated fits go high
  const Vector<T> &CovarianceMean{ GetCovarMean() };
  // Error check
  CheckVarCovar( sPrefix, CovarianceMean, Source );
  const Norm norm{ Norm::RawBinned };
  const std::size_t Extent{ Source.size2() };
  const std::size_t NumSamples{ Source.size1() };
  const std::size_t NumReplicas{ Replica.size1 };
  if( idx == idxCentral )
  {
    // I must not modify this mean
    Mean = CovarianceMean;
    // Don't do any replacements
    vRandom.resize( Extent );
    for( fint i = 0; i < Extent; ++i )
      vRandom[i] = i;
  }
  else
  {
    if( idx >= NumReplicas )
    {
      std::ostringstream ss;
      ss << sPrefix << " " << norm << " replica " << idx << " unavailable >= " << NumReplicas;
      throw std::runtime_error( ss.str().c_str() );
    }
    if( Replica.size2 != Extent )
    {
      std::ostringstream ss;
      ss << sPrefix << " " << norm << " replica " << idx << " extent " << Replica.size2
         << " != Source extent " << Extent;
      throw std::runtime_error( ss.str().c_str() );
    }
    Mean.MapRow( Replica, idx );
    GetRandom( vRandom, Seed, NumSamples, NumReplicas, idx );
  }
  return GetNorm( norm, vRandom.size );
}

/*****************************************************************************
 
 Do the actual work of computing (co)variance. Optionally with replacement

 ****************************************************************************/

template <typename T>
void JackBoot<T>::MakeCovarHelper( Matrix<T> &Covar, const MatrixView<T> &Source,
                                   const VectorView<T> &Mean, Real norm, Vector<fint> * pRandom )
{
  const std::size_t Extent{ Source.size2() };
  const std::size_t NumReplicas{ pRandom ? pRandom->size : Source.size1() };
  // Clear the covariance matrix
  Covar.resize( Extent, Extent );
  Covar = static_cast<T>( 0 );
  VectorView<T> Data;         // This is a view into each replica
  Vector<T> Error( Extent );  // This is the error on each replica
  // Sum the covariance entries on each replica - lower left triangle
  for( std::size_t replica = 0; replica < NumReplicas; ++replica )
  {
    Data.MapRow( Source, pRandom ? (*pRandom)[replica] : replica );
    for( std::size_t i = 0; i < Extent; ++i )
    {
      Error[i] = Data[i] - Mean[i];
      for( std::size_t j = 0; j <= i; ++j )
        Covar( i, j ) += Conjugate( Error[i] ) * Error[j];
    }
  }
  // Normalise and fill in upper right triangle
  for( std::size_t i = 0; i < Extent; ++i )
    for( std::size_t j = 0; j <= i; ++j )
    {
      const T z{ Covar( i, j ) * norm };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

template <typename T>
void JackBoot<T>::MakeVarHelper( Vector<T> &Var, const MatrixView<T> &Source,
                                 const VectorView<T> &Mean, Real norm, Vector<fint> * pRandom )
{
  const std::size_t Extent{ Source.size2() };
  const std::size_t NumReplicas{ pRandom ? pRandom->size : Source.size1() };
  // Clear the covariance matrix
  Var.resize( Extent );
  Var = static_cast<T>( 0 );
  VectorView<T> Data;         // This is a view into each replica
  // Sum the Variance entries on each replica
  for( std::size_t replica = 0; replica < NumReplicas; ++replica )
  {
    Data.MapRow( Source, pRandom ? (*pRandom)[replica] : replica );
    for( std::size_t i = 0; i < Extent; ++i )
      Var[i] += Squared( Data[i] - Mean[i] );
  }
  // Normalise
  for( std::size_t i = 0; i < Extent; ++i )
    Var[i] *= norm;
}

/*****************************************************************************
 
 Jackknife or Bootstrap replicas with central value

 ****************************************************************************/

/// Read a bootstrap replica from an HDF5 group
template <typename T>
void JackBoot<T>::Read( ::H5::Group &g, const std::string &Name )
{
  try
  {
    // Load the data and make sure they're valid
    const std::string sReplica{ EqualIgnoreCase( Name, "data" ) ? Name + s_S : Name };
    H5::ReadMatrix( g, sReplica, Replica );
    if( Replica.size1 == 0 || Replica.size2 == 0 )
      throw std::runtime_error( "Central replica empty reading JackBoot " + sReplica );
    if( !Replica.IsFinite() )
      throw std::runtime_error( "NaNs reading JackBoot " + sReplica );
    // Load the central replica and make sure it's valid
    const std::string sCentral{ Name + s_C };
    H5::ReadVector( g, sCentral, Central );
    if( Central.size != Replica.size2 )
      throw std::runtime_error( "Central replica extent " + std::to_string( Central.size )
                               + " != data extent " + std::to_string( Replica.size2 )
                               + " reading JackBoot " + sCentral );
    if( !Central.IsFinite() )
      throw std::runtime_error( "NaNs reading JackBoot " + sCentral + " central replica" );
    // Load (optional) replica mean and make sure it's valid
    MakeMean( ReplicaMean, Replica );
    if( !ReplicaMean.IsFinite() )
      throw std::runtime_error( "NaNs creating replica mean JackBoot " + Name );
  }
  catch(...)
  {
    clear();
    throw;
  }
}

/// Write a bootstrap replica to an HDF5 group
template <typename T>
void JackBoot<T>::Write( ::H5::Group &g, const std::string &Name ) const
{
  H5::WriteVector( g, Name + s_C, Central );
  H5::WriteMatrix( g, Name.compare( "data" ) ? Name : Name + s_S, Replica );
}

template class JackBoot<double>;
template class JackBoot<float>;
template class JackBoot<std::complex<double>>;
template class JackBoot<std::complex<float>>;

END_MLU_NAMESPACE
