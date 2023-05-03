/*****************************************************************************
 
 Jackknife and Bootstrap support

 Source file: JackBoot.cpp
 
 Copyright (C) 2019-2023

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
 ****************************************************************************/
/*  END LEGAL */

#include "JackBoot.hpp"
#include <random>

MLU_JackBoot_hpp

const std::string s_C{ "_C" };
static const std::string s_S{ "_S" };

// TODO: I need these from Common.hpp
std::string GetHostName();
void SetHostName( const std::string &NewHostName );
void MakeAncestorDirs( const std::string& Filename );

/*****************************************************************************

 This is a cache of random numbers used to construct bootstrap samples
 The number of replicas in each sample grows to support the largest request ever made

 *****************************************************************************/

const std::string RandomCache::sSeed{ "Seed" };
const std::string RandomCache::sSeedMachine{ "SeedMachine" };
const std::string RandomCache::sRandom{ "Random" };
std::string RandomCache::SaveRandomPrefix;
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
    const char * pDefaultString{ std::getenv( "MLUSeed" ) };
    if( pDefaultString && * pDefaultString )
    {
      if( EqualIgnoreCase( pDefaultString, "Jackknife" ) )
        DefaultSeed_ = SeedWildcard;
      else
      {
        try
        {
          DefaultSeed_ = FromString<SeedType>( pDefaultString );
        }
        catch(...){}
      }
    }
    DefaultSeedOK = true;
  }
  return DefaultSeed_;
}

SeedType RandomCache::DefaultSeed( SeedType NewDefaultSeed )
{
  DefaultSeed_ = NewDefaultSeed;
  DefaultSeedOK = true;
  return DefaultSeed_;
}

void RandomCache::DefaultSeed( const std::string &SeedOrFile )
{
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
    // Try to get the seed from HDF5 file
    H5E_auto2_t h5at;
    void      * f5at_p;
    ::H5::Exception::getAutoPrint(h5at, &f5at_p);
    ::H5::Exception::dontPrint();
    try
    {
      ::H5::H5File f;
      ::H5::Group  g;
      H5::OpenFileGroup( f, g, SeedOrFile, "Reading seed from: ", nullptr );
      Matrix<fint> Random;
      H5::ReadMatrix( g, sRandom, Random );
      ::H5::Attribute a = g.openAttribute( RandomCache::sSeed );
      a.read( H5::Equiv<SeedType>::Type, &MySeed );
      a.close();
      std::string SeedMachine;
      a = g.openAttribute( RandomCache::sSeedMachine );
      a.read( a.getStrType(), SeedMachine );
      a.close();
      SetHostName( SeedMachine );
      bGotSeed = true;
      Global.Put( MySeed, Random );
      std::cout << "Seed " << MySeed << " generated on " << SeedMachine << "\n";
    }
    catch(const ::H5::Exception &)
    {
      bGotSeed = false;
      ::H5::Exception::clearErrorStack();
    }
    catch(...)
    {
      ::H5::Exception::setAutoPrint(h5at, f5at_p);
      throw;
    }
    ::H5::Exception::setAutoPrint(h5at, f5at_p);
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

void RandomCache::SaveRandom( const Key &key, const Matrix<fint> &Random )
{
  if( !SaveRandomPrefix.empty() )
  {
    const std::string SeedMachine{ GetHostName() };
    std::string Filename{ SaveRandomPrefix };
    Filename.append( SeedMachine );
    Filename.append( 1, '.' );
    Filename.append( std::to_string( key.NumSamples ) );
    Filename.append( 1, '.' );
    Filename.append( sRandom );
    Filename.append( 1, '.' );
    Filename.append( SeedString( key.Seed ) );
    Filename.append( ".h5" );
    bool bOK{};
    try // to write in my format
    {
      MakeAncestorDirs( Filename );
      ::H5::H5File f( Filename, H5F_ACC_TRUNC );
      hsize_t Dims[2];
      Dims[0] = 1;
      ::H5::DataSpace ds1( 1, Dims );
      ::H5::Group g = f.createGroup( sRandom );
      H5::WriteMatrix( g, RandomCache::sRandom, Random );
      ::H5::Attribute a = g.createAttribute( RandomCache::sSeed, ::H5::PredType::STD_U32LE, ds1 );
      a.write( H5::Equiv<SeedType>::Type, &key.Seed );
      a.close();
      a = g.createAttribute( RandomCache::sSeedMachine, H5::Equiv<std::string>::Type, ds1 );
      a.write( H5::Equiv<std::string>::Type, SeedMachine );
      a.close();
      bOK = true;
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
  Matrix<fint> &m{ Map[key] };
  assert(((m.size1==0 && m.size2==0)||m.size2==NumSamples)&&"RandomCache::Get() NumSamples mismatch");
  if( m.size1 < NumReplicas )
  {
    // Not enough replicas - make some more
    Matrix<fint> New{ Make( Seed, NumSamples, NumReplicas ) };
    // Compare
    for( std::size_t i = 0 ; i < m.size1; ++i )
      for( std::size_t j = 0 ; j < NumSamples; ++j )
        if( New(i,j) != m(i,j) )
        {
          std::ostringstream os;
          os << "RandomCache::Get() Seed " << Seed << " samples " << NumSamples
             << ". Error growing from " << m.size1 << " to " << NumReplicas
             << " replicas. Old(" << i << ',' << j << ") " << m(i,j)
             << " != New(" << i << ',' << j << ") " << New(i,j);
          throw std::runtime_error( os.str().c_str() );
        }
    m = std::move( New );
    SaveRandom( key, m );
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

void RandomCache::Put( fint Seed, Matrix<fint> &New )
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
    auto pair{ Map.emplace( key, std::move( New ) ) };
    assert( pair.second && "Bug: I just searched for this, but it wasn't inserted" );
    // Make the caller's matrix point back to the new entry in the cache
    it = pair.first;
    New.MapView( it->second );
    SaveRandom( key, New );
  }
  else
  {
    // Something in cache - are the numbers compatible?
    Matrix<fint> &Old{ it->second };
    assert( Old.size1 && Old.size2 == NumSamples && "Invalid random numbers in cache" );
    std::size_t CompareReps{ std::min( New.size1, Old.size1 ) };
    for( std::size_t i = 0 ; i < CompareReps; ++i )
      for( std::size_t j = 0 ; j < NumSamples; ++j )
        if( New(i,j) != Old(i,j) )
        {
          std::ostringstream os;
          os << "RandomCache::Put() Seed " << SeedString( Seed )
             << " samples " << NumSamples
             << ". New samples incompatible with cache. Old(" << i << ',' << j << ") " << Old(i,j)
             << " != New(" << i << ',' << j << ") " << New(i,j);
          throw std::runtime_error( os.str().c_str() );
        }
    if( New.size1 > Old.size1  )
    {
      Old = std::move( New ); // New has more replicas - move this into cache
      SaveRandom( key, Old );
    }
    New.MapView( Old ); // Make callers matrix point to cache
  }
}

/*****************************************************************************
 
 Get random numbers appropriate for given parameters
 NumReplicas should be the total number of replicas required (so that random numbers created once for all replicas)
 
 *****************************************************************************/

void JackBootBase::GetRandom( Vector<fint> &vRandom, fint Seed,
                              std::size_t NumSamples, std::size_t NumReplicas, std::size_t Replica )
{
  if( Replica >= NumReplicas )
  {
    std::ostringstream ss;
    ss << "JackBoot<T>::GetRandom Request for replica " << Replica << " of " << NumReplicas
       << " replicas";
    throw std::runtime_error( ss.str().c_str() );
  }
  // The presence or absence of random numbers tells us whether we are doing bootstrap or jackknife
  if( Seed == SeedWildcard )
  {
    // Jackknife on non-central replica
    if( NumReplicas != NumSamples )
    {
      std::ostringstream ss;
      ss << "JackBoot<T>::GetRandom Jackknife NumReplicas " << NumReplicas
         << " != NumSamples " << NumSamples;
      throw std::runtime_error( ss.str().c_str() );
    }
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
void JackBoot<T>::MakeMean( Vector<T> &vCentral, const MatrixView<T> &Source )
{
  if( Source.size1() == 0 || Source.size2() == 0 )
    throw std::runtime_error( "JackBoot::MakeMean Source matrix empty" );
  std::size_t NumReplicas{ Source.size1() };
  std::size_t Extent{ Source.size2() };
  vCentral.resize( Extent );
  vCentral = static_cast<T>( 0 );
  for( std::size_t i = 0; i < NumReplicas; ++i )
    for( std::size_t j = 0; j < Extent; ++j )
      vCentral[j] += Source(i,j);
  for( std::size_t j = 0; j < Extent; ++j )
    vCentral[j] /= NumReplicas;
}

/*****************************************************************************

 Perform jackknife or bootstrap (based on seed) resampling
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::Resample( const MatrixView<T> &Source, std::size_t NumReplicas )
{
  if( Source.size1() == 0 || Source.size2() == 0 )
    throw std::runtime_error( "JackBoot::Resample Source matrix empty" );
  std::size_t Extent{ Source.size2() };
  std::size_t NumSamples{ Source.size1() };
  if( Seed == SeedWildcard )
  {
    if( NumSamples < 2 )
      throw std::runtime_error( "JackBoot::Resample Can't jackknife with "
                    + std::to_string( NumSamples ) + " samples. Minimum of 2 required" );
    NumReplicas = NumSamples;
  }
  else if( NumReplicas < 2 )
    throw std::runtime_error( "JackBoot::Resample Can't bootstrap into "
                  + std::to_string( NumReplicas ) + " bootstrap replicas (minimum 2)" );
  MakeMean( Central, Source );
  Replica.resize( NumReplicas, Extent );
  Vector<fint> Random;
  VectorView<T> Sample;
  for( std::size_t idx = 0; idx < NumReplicas; ++idx )
  {
    GetRandom( Random, Seed, NumSamples, NumReplicas, idx );
    for( std::size_t Resample = 0; Resample < Random.size; ++Resample )
    {
      Sample.MapRow( Source, Random[Resample] );
      for( std::size_t i = 0; i < Extent; ++i )
      {
        if( Resample == 0 )
          Replica( idx, i ) = Sample[i];
        else
          Replica( idx, i ) += Sample[i];
      }
    }
    // Turn the sum into an average
    if( NumSamples > 1 )
      for( std::size_t i = 0; i < Extent; ++i )
        Replica( idx, i ) /= NumSamples;
  }
}

/*****************************************************************************

 Make statistics for this sample

 Needs to work even when empty
 
 *****************************************************************************/

template <typename T>
void JackBoot<T>::MakeStatistics( std::vector<ValWithEr<T>> &vStats ) const
{
  vStats.resize( extent() );
  if( extent() )
  {
    // Default to showing the central values
    for( std::size_t i = 0; i < extent(); ++i )
    {
      vStats[i] = Central[i];
      vStats[i].Check = NumReplicas() ? 1 : 0;
    }
    if( NumReplicas() > 1 )
    {
      if( Seed == SeedWildcard )
      {
        // Jackknife statistics
        std::vector<std::size_t> Count( extent(), 0 );
        Vector<T> Var( extent() );
        Var = static_cast<T>( 0 );
        // Sum the Variance entries on each replica
        for( std::size_t replica = 0; replica < NumReplicas(); ++replica )
        {
          for( std::size_t j = 0; j < extent(); ++j )
          {
            if( IsFinite( Replica( replica, j ) ) )
            {
              Var[j] += Squared( Replica( replica, j ) - Central[j] );
              ++Count[j];
            }
          }
        }
        // Normalise
        Vector<T> Column;
        for( std::size_t j = 0; j < extent(); ++j )
        {
          if( Count[j] )
          {
            Var[j] *= GetNorm( Norm::Jackknife, Count[j] );
            Var[j] = std::sqrt( Var[j] );
            vStats[j].Low -= Var[j];
            vStats[j].High += Var[j];
            Column.MapColumn( const_cast<Matrix<T> &>( Replica ), j );
            Column.MinMax( vStats[j].Min, vStats[j].Max, true );
          }
          vStats[j].Check = static_cast<Real>( Count[j] ) / NumReplicas();
        }
      }
      else
      {
        // Bootstrap statistics
        std::vector<T> Buffer( NumReplicas() );
        for( std::size_t j = 0; j < extent(); ++j )
        {
          std::size_t Count{};
          for( std::size_t replica = 0; replica < NumReplicas(); ++replica )
          {
            if( IsFinite( Replica( replica, j ) ) )
              Buffer[Count++] = Replica( replica, j );
          }
          vStats[j].Get( Central[j], Buffer, Count );
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
void JackBoot<T>::MakeCovar( Matrix<T> &Covar, const MatrixView<T> &Source,
                             const VectorView<T> &Mean, Norm norm )
{
  CheckVarCovar( "static JackBoot::MakeCovar()", Source, Mean );
  const Real A{ GetNorm( norm, Source.size1() ) };
  MakeCovarHelper( Covar, Source, Mean, A );
}

template <typename T>
void JackBoot<T>::MakeVar( Vector<T> &Var, const MatrixView<T> &Source,
                           const VectorView<T> &Mean, Norm norm )
{
  CheckVarCovar( "static JackBoot::MakeVar()", Source, Mean );
  const Real A{ GetNorm( norm, Source.size1() ) };
  MakeVarHelper( Var, Source, Mean, A );
}

/*****************************************************************************
 
 Make covariance matrix or variance vector given a bootstrap or jackknife and the source matrix it was created from.
 Resampling is possible for any replica, with the (co)variance coming from the source data, resampled appropriately.
 Normalisation is always `Norm::RawBinned`, i.e. `\Sigma_{\bar{\vb{x}}}`error of the mean
 
 *****************************************************************************/

template <typename T>
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
}

/*****************************************************************************
 
 Common validation for covariance and variance

 ****************************************************************************/

/// Common validation for (co)variance given any matrix and a vector of central values
template <typename T>
void JackBoot<T>::CheckVarCovar( const std::string &sPrefix,
                                 const MatrixView<T> &Source, const VectorView<T> &vCentral )
{
  // Error check
  const std::size_t Extent{ Source.size2() };
  if( Extent == 0 )
    throw std::runtime_error( sPrefix + " empty extent" );
  if( vCentral.size() != Extent )
    throw std::runtime_error( sPrefix + " mismatch: " + std::to_string( vCentral.size() )
                             + " central values != extent " + std::to_string( vCentral.size() ) );
}

/// Common validation for (co)variance given a bootstrap / jackknife and the matrix it was made from
template <typename T> typename JackBoot<T>::Real
JackBoot<T>::CheckVarCovar( const std::string &sPrefix, VectorView<T> &Mean, Vector<fint> &vRandom,
                            const MatrixView<T> &Source, std::size_t idx ) const
{
  // Error check
  CheckVarCovar( sPrefix, Source, Central );
  const Norm norm{ Norm::RawBinned };
  const std::size_t Extent{ Source.size2() };
  const std::size_t NumSamples{ Source.size1() };
  const std::size_t NumReplicas{ Replica.size1 };
  if( idx == idxCentral )
  {
    // I must not modify this mean
    Mean = Central;
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
  H5::ReadVector( g, Name + s_C, Central );
  H5::ReadMatrix( g, Name.compare( "data" ) ? Name : Name + s_S, Replica );
  if( Central.size != Replica.size2 )
  {
    std::ostringstream os;
    os << "JackBoot sample " << Name << ": central extent " << Central.size
       << " != replica columns " << Replica.size2;
    throw std::runtime_error( os.str().c_str() );
  }
  if( !Central.IsFinite() )
    throw std::runtime_error( "NaNs reading JackBoot " + Name + " central replica" );
  if( !Replica.IsFinite() )
    throw std::runtime_error( "NaNs reading JackBoot " + Name + " replicas" );
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

MLU_JackBoot_hpp_end
