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
#define MLU_JackBoot_hpp namespace Common {
#define MLU_JackBoot_hpp_end };

#include <MLU/HDF5.hpp>

#include <cstdint>

MLU_JackBoot_hpp

using SeedType = fint; // TODO: Merge all occurences into a single type
static constexpr SeedType SeedWildcard{};

extern const std::string s_C;

/**
 This is a cache of random numbers used to construct bootstrap samples
 The number of replicas in each sample grows to support the largest request ever made
 */
struct RandomCache
{
  static SeedType DefaultSeed();
  static std::size_t DefaultNumReplicas();
  static RandomCache Global; // This is the global cache random numbers come from
protected:
  struct Key
  {
    fint  Seed; // bootstrap seed for random numbers
    std::size_t NumSamples; // how many raw data samples were bootstrapped
    struct Less { bool operator()( const Key &lhs, const Key &rhs ) const; };
  };
  using Value = Matrix<fint>;
  using MapT = std::map<Key, Value, Key::Less>;
  MapT Map;
  Matrix<fint> Make( fint &Seed, std::size_t NumReplicas, std::size_t NumSamples );
public:
  /** Get random numbers for the given Seed and NumSamples
   
   Might return more replicas than asked for (if previous calls have asked for more replcas)
   
   - Warning: If Seed == SeedWildcard, it will be changed to the default value from MLUSeed environment variable,
   or ``std::mt19937::default_seed`` =  5489u if not present
   */
  Matrix<fint> Get(       fint &Seed, std::size_t NumSamples, std::size_t NumReplicas );
  Matrix<fint> Get( const fint &Seed, std::size_t NumSamples, std::size_t NumReplicas );
  /// Put random numbers in cache (probably loaded from file). Throw error if incompatible
  void Put( fint &Seed, Matrix<fint> &Random );
};

struct JackBootBase
{
  enum class Norm{ RawBinned,   // Sample covariance OF MEAN from raw / binned data
                   Bootstrap,   // Sample covariance OF MEAN from bootstraps
                   Jackknife,   // Sample covariance OF MEAN from Jackknife
                   SampleCovar};// Sample covariance from raw / binned data
  static void GetRandom( Vector<fint> &vRandom, fint Seed,
                         std::size_t NumSamples, std::size_t NumReplicas, int Replica );
};

std::ostream &operator<<( std::ostream &os, JackBootBase::Norm norm );

/// Jackknife or Bootstrap replicas with central value
template<typename T> struct JackBoot : public JackBootBase
{
  using Norm = JackBootBase::Norm;
  using Real = typename is_complex<T>::Scalar;
  static constexpr int idxCentral{ -1 };
  Vector<T> Central;
  Matrix<T> Replica;
  /// bootstrap: random number seed for this sample. SeedWildcard = jackknife
  SeedType Seed = RandomCache::DefaultSeed();
  void clear() { Central.clear(); Replica.clear(); Seed = RandomCache::DefaultSeed(); }
  static Real GetNorm( Norm norm, std::size_t NumReplicas );
  /// Make the mean of any matrix
  static void MakeMean( Vector<T> &vCentral, const MatrixView<T> &Source );
  /// Perform jackknife or bootstrap (based on seed) resampling
  /// NumBoot is the number of bootstrap replicas (ignored for jackknife)
  void Resample( const MatrixView<T> &Source, std::size_t NumBoot );
  /**
   Make covariance matrix or variance vector given any matrix (raw/binned, bootstrap or jackknife) and a vector of central values.
   
   No resampling is possible (i.e. central replica only).
   Caller specifies normalisation required, most of which compute `\Sigma_{\bar{\vb{x}}}`, i.e. error of the mean
   but Norm::SampleCovar can be used to get the sample covariance if required.
   **/
  static void MakeCovar( Matrix<T> &Covar, const MatrixView<T> &Source,
                         const VectorView<T> &Mean, Norm norm );
  static void MakeVar( Vector<T> &Var, const MatrixView<T> &Source,
                       const VectorView<T> &Mean, Norm norm );
  /**
   Make covariance matrix or variance vector given a bootstrap or jackknife and the source matrix it was created from.
   
   Resampling is possible for any replica, with the (co)variance coming from the source data, resampled appropriately.
   Normalisation is always `Norm::RawBinned`, i.e. `\Sigma_{\bar{\vb{x}}}`error of the mean
   */
  void MakeCovar( Matrix<T> &Covar, const MatrixView<T> &Source, int idx ) const;
  void MakeVar( Vector<T> &Var, const MatrixView<T> &Source, int idx ) const;

  inline int extent() const { return static_cast<int>( Central.size ); } // Num samples per replica
  inline int NumReplicas() const { return static_cast<int>( Replica.size1 ); }
  inline void resize( int NumReplicas, int Extent )
  {
    Central.resize( Extent );
    Replica.resize( NumReplicas, Extent );
  }
  JackBoot() = default;
  JackBoot( int NumReplicas, int Extent ) { resize( NumReplicas, Extent ); }
  void Read( ::H5::Group &g, const std::string &Name );
  void Write( ::H5::Group &g, const std::string &Name ) const;
  inline T * operator[]( int idx )
  {
    if( Central.size == 0 )
      throw std::runtime_error( "Can't access empty JackBoot" );
    if( idx == -1 )
      return Central.Data();
    if( Replica.size2 != Central.size )
      throw std::runtime_error( "BoootRep extent mismatch" );
    if( idx >= Replica.size1 )
      throw std::runtime_error( "BoootRep index " + std::to_string( idx ) + " >= NumReplicas "
                               + std::to_string( NumReplicas() ) );
    return Replica.Data() + Replica.tda * idx;
  }
  inline void MapRow( Vector<T> &v, int idx )
  {
    if( idx == idxCentral )
      v.MapView( Central );
    else
      v.MapRow( Replica, idx );
  }
  inline void MapRow( VectorView<T> &v, int idx ) const
  {
    if( idx == idxCentral )
      v = Central;
    else
      v.MapRow( Replica, idx );
  }
  inline void MapColumn( VectorView<T> &v, int Column ) const
  {
    v.MapColumn( Replica, Column );
  }
protected:
  /// Common validation for (co)variance given any matrix and a vector of central values
  static void CheckVarCovar( const std::string &sPrefix,
                             const MatrixView<T> &Source, const VectorView<T> &vCentral );
  /// Common validation for (co)variance given a bootstrap / jackknife and the matrix it was made from
  Real CheckVarCovar( const std::string &sPrefix, VectorView<T> &Mean, Vector<fint> &Random,
                      const MatrixView<T> &Source, int idx ) const;
  /// Do the actual work of computing (co)variance. Optionally with replacement
  static void MakeCovarHelper( Matrix<T> &Covar, const MatrixView<T> &Source,
                               const VectorView<T> &Mean, Real norm, Vector<fint> * pRandom=nullptr );
  static void MakeVarHelper( Vector<T> &Var, const MatrixView<T> &Source,
                             const VectorView<T> &Mean, Real norm, Vector<fint> * pRandom = nullptr );
};

MLU_JackBoot_hpp_end
#endif
