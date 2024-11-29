/**
 
 Covariance matrices for fits
 
 Source file: Covar.hpp
 
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

#ifndef Covar_hpp
#define Covar_hpp

#include "MultiFit.hpp"

// Fitter parameters
struct CovarParams
{
  const DataSet &ds;
  // Whether to rebin the data before computing covariances
  std::vector<int> RebinSize;
  // Where am I building the covariance from
  using SS = MLU::SampleSource;
  SS Source;
  int idxJackBoot;

  // 0 = build correlation from CovarSource, scaled by variance of data on each bootstrap replica
  // else estimate covariance by re-bootstrapping each replica using this # of secondary replicas
  int CovarNumBoot = 0;
  // Freeze covariance to the central replica
  bool bFreeze;
  // Manually loaded inverse covariance matrix
  Matrix Covar;
  // Constructor
  CovarParams( const MLU::CommandLine &cl, DataSet &ds, bool bNoInit = false );
  /// Can we compute (co)variance on non-central replicas? Only possible if we have rebinned data
  inline bool SupportsUnfrozen() const { return ds.NumSamplesBinned(); }
  // Source is a bootstrap - only valid to compute variance on central replica
  inline bool SourceIsBootstrap() const { return Source == SS::Bootstrap; }
  // This is the appropriate parameter for me to use in the T^2 distribution
  inline int CovarSampleSize() const
  {
    return SourceIsBootstrap() ? ds.MaxSampleSize
                               : ds.front().NumSamples( Source==SS::Raw ? SS::Raw : SS::Binned );
  }
  /**
   Can I make the covariance matrix in one step for this replica?
   If false, then I make correlation matrix from source, then scale it by variance of data on this replica
   */
  /*inline bool OneStep( int Replica ) const
  {
    // TODO: support covar from replica data
    return CovarNumBoot || (Replica == Fold::idxCentral && (Source == SS::Binned || Source == SS::Bootstrap) );
  }*/
  // How many covariance samples are there
  inline int CovarCount() const { return ds.front().NumSamples( Source ); }
};

std::ostream & operator<<( std::ostream &os, const CovarParams &cp );

#endif // Covar_hpp
