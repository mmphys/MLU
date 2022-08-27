/*************************************************************************************
 
 Covariance matrices for fits
 
 Source file: Covar.hpp
 
 Copyright (C) 2019-2022
 
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
 *************************************************************************************/
/*  END LEGAL */

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
  using SS = Common::SampleSource;
  SS Source;
  // 0 = build correlation from CovarSource, then scale by the variance of the data on each bootstrap replica
  // else this is the number of bootstrap replicas to use when estimating covariance
  int CovarNumBoot;
  // Random numbers to get the covariance from a bootstrap
  using fint = Common::fint;
  Common::MatrixView<fint> CovarRandom;
  std::vector<fint> vCovarRandom; // Contains the random numbers iff I generated them
  // Freeze covariance to the central replica
  bool bFreeze;
  // Manually loaded correlation matrix
  Matrix Covar;
  // Constructor
  CovarParams( const Common::CommandLine &cl, DataSet &ds );
  // Can we compute (co)variance on non-central replicas - i.e. do we have rebinned data
  inline bool SupportsUnfrozen() const { return ds.corr[0].NumSamplesBinned(); }
  // Source is a bootstrap - only valid to compute variance on central replica
  inline bool SourceIsBootstrap() const
  { return Source == SS::Bootstrap || ( Source == SS::Raw && ds.corr[0].bRawBootstrap ); }
  // This is the appropriate parameter for m to use in the T^2 distribution
  inline int CovarSampleSize() const
  {
    return SourceIsBootstrap() ? ds.corr[0].SampleSize
                               : ds.corr[0].NumSamples( Source==SS::Raw ? SS::Raw : SS::Binned );
  }
  // Make correlation matrix, then scale it by variance of data on each replica
  inline bool OneStep( int Replica ) const
  {
    // TODO: support covar from replica data
    return CovarNumBoot || (Replica == Fold::idxCentral && (Source == SS::Binned || Source == SS::Bootstrap) );
  }
  // How many covariance samples are there
  inline int CovarCount() const { return ds.corr[0].NumSamples( Source ); }
};

std::ostream & operator<<( std::ostream &os, const CovarParams &cp );

#endif // Covar_hpp
