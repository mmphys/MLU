/*****************************************************************************
 
 Jackknife and Bootstrap support

 Source file: JackBoot.hpp
 
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

#ifndef MLU_JackBoot_hpp
#define MLU_JackBoot_hpp

#include <MLU/HDF5.hpp>
#include <MLU/Posix.hpp>

#include <cstdint>

BEGIN_MLU_NAMESPACE

static constexpr SeedType SeedWildcard{};

extern const std::string s_C;

/**
 This is a cache of random numbers used to construct bootstrap samples
 The number of replicas in each sample grows to support the largest request ever made
 */
struct RandomCache
{
  static const std::string sSeed;
  static const std::string sSeedMachine;
  static const std::string sRandom;
  static SeedType DefaultSeed();
  static SeedType DefaultSeed( SeedType NewDefaultSeed );
  /// Set default seed, eg from command-line switch. NB: Can be a random number file
  static void     DefaultSeed( const std::string &SeedOrFile );
  static std::size_t DefaultNumReplicas();
  static std::size_t DefaultNumReplicas( std::size_t NewDefaultNumReplicas );
  static RandomCache Global; // This is the global cache random numbers come from
  inline static SeedType Seed( const std::string &TextSeed )
  {
    if( EqualIgnoreCase( TextSeed, "jackknife" ) )
      return SeedWildcard;
    return FromString<SeedType>( TextSeed );
  }
  inline static std::string SeedString( SeedType Seed )
  {
    if( Seed == SeedWildcard )
      return "jackknife";
    return std::to_string( Seed );
  }
  static void Write( ::H5::Group &g, const Matrix<fint> &Random,
                     SeedType Seed, const std::string &Machine );
  static void Read( ::H5::Group &g, Matrix<fint> &Random, SeedType &Seed, std::string &Machine );
  static bool Read( const std::string &Filename, Matrix<fint> &Random,
                    SeedType &Seed, std::string &Machine );
  static void SetCachePrefix( std::string &NewCachePrefix );
  static bool SaveFatFiles(); // True if random numbers should be saved with data files
protected:
  static const std::string &GetCachePrefix();
  static bool DefaultSeedOK;
  static SeedType DefaultSeed_;
  static bool DefaultNumReplicasOK;
  static std::size_t DefaultNumReplicas_; // Mike's default number of replicas
  struct Key
  {
    fint  Seed; // bootstrap seed for random numbers
    std::size_t NumSamples; // how many raw data samples were bootstrapped
    struct Less { bool operator()( const Key &lhs, const Key &rhs ) const; };
  };
  struct Value
  {
    Matrix<fint> m;
    bool bOnDisk = false;
  };
  using MapT = std::map<Key, Value, Key::Less>;
  MapT Map;
  static bool CacheDirPrefixOK;
  static std::string CacheDirPrefix;
public:
  static std::string PrependCachePath( const std::string &sFilename );
protected:
  static std::string GetCachePath( const Key &key );
  /// Make random numbers for a bootstrap sample
  Matrix<fint> Make( fint &Seed, std::size_t NumSamples, std::size_t NumReplicas
                                                          = DefaultNumReplicas());
  static void Compare( const Matrix<fint> &l, const Matrix<fint> &r );
  static void SaveRandom( const Key &key, Value &value );
  void Put( fint Seed, Matrix<fint> &Random, bool bWrite );
public:
  /** Get random numbers for the given Seed and NumSamples
   
   Might return more replicas than asked for (if previous calls have asked for more replcas)
   
   - Warning: If Seed == SeedWildcard, it will be changed to the default value from MLUSeed environment variable,
   or ``std::mt19937::default_seed`` =  5489u if not present
   */
  Matrix<fint> Get(       fint &Seed, std::size_t NumSamples, std::size_t NumReplicas
                                                                = DefaultNumReplicas() );
  Matrix<fint> Get( const fint &Seed, std::size_t NumSamples, std::size_t NumReplicas
                                                                = DefaultNumReplicas() );
  /// Put random numbers in cache (probably loaded from file). Throw error if incompatible
  void Put( fint Seed, Matrix<fint> &Random ) { Put( Seed, Random, true ); }
  /// Get the random numbers for the specified replica.
};

struct JackBootBase
{
  static const std::string sBootstrap;
  static const std::string sJackknife;
  /// Guaranteed that incrementing idxCentral moves to record zero
  static constexpr std::size_t idxCentral{ std::numeric_limits<std::size_t>::max() };
  static constexpr std::size_t idxReplicaMean{ idxCentral - 1u };
  static const std::string s_S;
  static bool UseCentralCovar(); // false unless MLUCovarCentral environment variable set
  static std::size_t getCovarMeanIdx() { return UseCentralCovar() ? idxCentral : idxReplicaMean; }
  enum class Norm{ RawBinned,   // Sample covariance OF MEAN from raw / binned data
                   Bootstrap,   // Sample covariance OF MEAN from bootstraps
                   Jackknife,   // Sample covariance OF MEAN from Jackknife
                   SampleCovar};// Sample covariance from raw / binned data
  static void GetRandom( Vector<fint> &vRandom, fint Seed,
                         std::size_t NumSamples, std::size_t NumReplicas, std::size_t Replica );
  inline static const std::string &SeedTypeString( SeedType Seed )
  {
    if( Seed == SeedWildcard )
      return sJackknife;
    return sBootstrap;
  }
};

std::ostream &operator<<( std::ostream &os, JackBootBase::Norm norm );

template<typename T> struct JackBoot;

/// This is a single column from a JackBoot ... but can be set to a constant instead
template<typename T> struct JackBootColumn
{
protected:
  union {
    T Constant;
    std::size_t Column;
  };
  const JackBoot<T> *pJB;
public:
  inline bool IsConstant() const { return pJB == nullptr; }
  inline bool IsIdentity() const { return IsConstant() && Constant == static_cast<T>( 1 ); }
  inline const T &operator[]( std::size_t Row ) const { return pJB ? (*pJB)(Row,Column) : Constant; }
  explicit JackBootColumn( const JackBoot<T> &JB, std::size_t column ) : pJB{&JB}, Column{column} {}
  explicit JackBootColumn( T constant ) : Constant{constant}, pJB{nullptr} {}
  void operator=( T constant ) { Constant = constant; pJB = nullptr; }
  JackBootColumn() : JackBootColumn{ static_cast<T>( 0 ) } {}
  JackBootColumn( const JackBootColumn &o )
  {
    pJB = o.pJB;
    if( pJB )
      Column = o.Column;
    else
      Constant = o.Constant;
  }
};

/// Jackknife or Bootstrap replicas with central value
template<typename T> struct JackBoot : public JackBootBase
{
  static constexpr std::size_t idxCentral{ JackBootBase::idxCentral };
  static constexpr std::size_t idxReplicaMean{ JackBootBase::idxReplicaMean };
  using Norm = JackBootBase::Norm;
  using Real = typename is_complex<T>::Real;
protected:
  Vector<T> Central;      // Averages of the source data that were resampled
  Vector<T> ReplicaMean;  // Mean of the replicas. Bootstrap/Jackknife differs from/same as Central
public:
  Matrix<T> Replica;      // Jackknife / bootstrap replicas
  /// Read-only view of a column
  JackBootColumn<T> Column( std::size_t column ) const { return JackBootColumn<T>( *this, column ); }
  // This is the average from the original data
  inline       Vector<T> &GetCentral()       { return Central; }
  inline const Vector<T> &GetCentral() const { return Central; }
  // This is the mean of the replicas. Use this to build (co)variance from replicas
  inline       Vector<T> &GetReplicaMean()       { return ReplicaMean; }
  inline const Vector<T> &GetReplicaMean() const { return ReplicaMean; }
  // Use this as the mean when building (co)variance from replicas - TODO: Definitely not central?
  inline       Vector<T> &GetCovarMean()       { return UseCentralCovar() ? Central : ReplicaMean; }
  inline const Vector<T> &GetCovarMean() const { return UseCentralCovar() ? Central : ReplicaMean; }
  /// bootstrap: random number seed for this sample. SeedWildcard = jackknife
  SeedType Seed = RandomCache::DefaultSeed();
  void clear();
  inline T &operator()( std::size_t Row, std::size_t Column )
  {
    if( Row == idxCentral )
      return Central[Column];
    if( Row == idxReplicaMean )
      return ReplicaMean[Column];
    return Replica( Row, Column );
  }
  inline const T &operator()( std::size_t Row, std::size_t Column ) const
  {
    if( Row == idxCentral )
      return Central[Column];
    if( Row == idxReplicaMean )
      return ReplicaMean[Column];
    return Replica( Row, Column );
  }
  inline std::string SeedString() const { return RandomCache::SeedString( Seed ); }
  inline const std::string &SeedTypeString() const { return JackBootBase::SeedTypeString( Seed ); }
  static Real GetNorm( Norm norm, std::size_t NumReplicas );
  Norm MyNorm() const { return Seed == SeedWildcard ? Norm::Jackknife : Norm::Bootstrap; }
  /// Make the mean of any matrix
  static void MakeMean( Vector<T> &vMean, const MatrixView<T> &Source );
  /// Make the mean of a matrix while resampling it - e.g. from raw/binned source data
  static void MakeMean( Vector<T> &vMean, Vector<fint> &vRandom, const MatrixView<T> &Source,
                        std::size_t Replica,
                        std::size_t NumReplicas = RandomCache::DefaultNumReplicas(),
                        SeedType Seed = RandomCache::DefaultSeed() );
  /// Make the mean of my samples
  void MakeMean() { MakeMean( ReplicaMean, Replica ); }
  /// Perform jackknife or bootstrap (based on seed) resampling
  /// NumBoot is the number of bootstrap replicas (ignored for jackknife)
  void Resample( const MatrixView<T> &Source, std::size_t NumBoot );
protected:
  /// Make statistics for this sample
  void MakeStatistics( ValWithEr<T> *vStats, std::size_t StartCol, std::size_t NumCols ) const;
public:
  /// Make statistics for one column
  inline void MakeStatistics( ValWithEr<T> &ve, std::size_t Column ) const
  {
    MakeStatistics( &ve, Column, 1 );
  }
  /// Make statistics for this sample
  inline void MakeStatistics( std::vector<ValWithEr<T>> &vStats ) const
  {
    vStats.resize( extent() );
    if( extent() )
      MakeStatistics( &vStats[0], 0, extent() );
  }
  /**
   Make covariance matrix or variance vector given any matrix (raw/binned, bootstrap or jackknife) and a vector of central values.
   
   No resampling is possible (i.e. central replica only).
   Caller specifies normalisation required, most of which compute `\Sigma_{\bar{\vb{x}}}`, i.e. error of the mean
   but Norm::SampleCovar can be used to get the sample covariance if required.
   **/
  static void MakeCovar( Matrix<T> &Covar, const VectorView<T> &Mean,
                         const MatrixView<T> &Source, Norm norm );
  static void MakeVar( Vector<T> &Var, const VectorView<T> &Mean,
                       const MatrixView<T> &Source, Norm norm );
  /// Make the (co)variance of a matrix while resampling it - e.g. from raw/binned source data
  static void MakeCovar( Matrix<T> &Covar,
                         const MatrixView<T> &Source, std::size_t Replica,
                         std::size_t NumReplicas = RandomCache::DefaultNumReplicas(),
                         SeedType Seed = RandomCache::DefaultSeed() );
  static void MakeVar( Vector<T> &Var,
                       const MatrixView<T> &Source, std::size_t Replica,
                       std::size_t NumReplicas = RandomCache::DefaultNumReplicas(),
                       SeedType Seed = RandomCache::DefaultSeed() );
  /// Make covariance matrix or variance vector from this JackBoot
  void MakeCovar( Matrix<T> &Covar ) const { MakeCovar( Covar, GetCovarMean(), Replica, MyNorm() ); }
  void MakeVar( Vector<T> &Var ) const { MakeVar( Var, GetCovarMean(), Replica, MyNorm() ); }
  /**
   Make covariance matrix or variance vector given a bootstrap or jackknife and the source matrix it was created from.
   
   Resampling is possible for any replica, with the (co)variance coming from the source data, resampled appropriately.
   Normalisation is always `Norm::RawBinned`, i.e. `\Sigma_{\bar{\vb{x}}}`error of the mean
   */
  /*void MakeCovar( Matrix<T> &Covar, const MatrixView<T> &Source, std::size_t idx ) const;
  void MakeVar( Vector<T> &Var, const MatrixView<T> &Source, std::size_t idx ) const;*/

  inline std::size_t extent() const { return Central.size; } // Num data points (columns) per replica
  inline std::size_t NumReplicas() const { return Replica.size1; }
  inline void resize( std::size_t NumReplicas, std::size_t Extent )
  {
    Central.resize( Extent );
    Replica.resize( NumReplicas, Extent );
    ReplicaMean.resize( Extent );
  }
  JackBoot() = default;
  JackBoot( std::size_t NumReplicas, std::size_t Extent ) { resize( NumReplicas, Extent ); }
  inline bool IsFinite() const
  {
    return Central.IsFinite() && ReplicaMean.IsFinite() && Replica.IsFinite();
  }
  void Read( ::H5::Group &g, const std::string &Name );
  void Write( ::H5::Group &g, const std::string &Name ) const;
  inline void ValidateReplica( std::size_t idx ) const
  {
    if( Central.size == 0 )
      throw std::runtime_error( "Can't access empty JackBoot" );
    if( idx == idxReplicaMean )
    {
      if( ReplicaMean.size != Central.size )
        throw std::runtime_error( "BoootRep ReplicaMean extent mismatch" );
    }
    else if( idx != idxCentral )
    {
      if( Replica.size2 != Central.size )
        throw std::runtime_error( "BoootRep extent mismatch" );
      if( idx >= Replica.size1 )
        throw std::runtime_error( "BoootRep index " + std::to_string( idx ) + " >= NumReplicas "
                                 + std::to_string( NumReplicas() ) );
    }
  }
  inline void MapRow( Vector<T> &v, std::size_t idx )
  {
    ValidateReplica( idx );
    if( idx == idxCentral )
      v.MapView( Central );
    else if( idx == idxReplicaMean )
     v.MapView( ReplicaMean );
    else
      v.MapRow( Replica, idx );
  }
  inline void MapRow( VectorView<T> &v, std::size_t idx ) const
  {
    ValidateReplica( idx );
    if( idx == idxCentral )
      v = Central;
    else if( idx == idxReplicaMean )
      v = ReplicaMean;
    else
      v.MapRow( Replica, idx );
  }
  inline void MapColumn( VectorView<T> &v, std::size_t Column ) const
  {
    ValidateReplica( 0 );
    if( Column >= Replica.size2 )
      throw std::runtime_error( "BoootRep column " + std::to_string( Column ) + " >= "
                               + std::to_string( extent() ) );
    v.MapColumn( Replica, Column );
  }
  inline VectorView<T> operator[]( std::size_t idx )
  {
    VectorView<T> vv;
    MapRow( vv, idx );
    return vv;
  }
protected:
  /// Common validation for (co)variance given any matrix and a vector of central values
  static void CheckVarCovar( const std::string &sPrefix, const VectorView<T> &vMean,
                             const MatrixView<T> &Source );
  /// Common validation for (co)variance given a bootstrap / jackknife and the matrix it was made from
  Real CheckVarCovar( const std::string &sPrefix, VectorView<T> &Mean, Vector<fint> &Random,
                      const MatrixView<T> &Source, std::size_t idx ) const;
  /// Do the actual work of computing (co)variance. Optionally with replacement
  static void MakeCovarHelper( Matrix<T> &Covar, const MatrixView<T> &Source,
                               const VectorView<T> &Mean, Real norm, Vector<fint> * pRandom=nullptr );
  static void MakeVarHelper( Vector<T> &Var, const MatrixView<T> &Source,
                             const VectorView<T> &Mean, Real norm, Vector<fint> * pRandom = nullptr );
};

template <typename T>
inline bool IsFinite( const JackBoot<T> &jb )
{
  return jb.IsFinite();
}

END_MLU_NAMESPACE
#endif
