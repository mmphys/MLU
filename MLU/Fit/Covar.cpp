/**
 
 Covariance matrices for fits
 
 Source file: Covar.cpp
 
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
#include "Covar.hpp"

// Work out where the covariance matrix comes from
CovarParams::CovarParams( const MLU::CommandLine &cl, DataSet &dsrw, bool bNoInit ) : ds{ dsrw }
{
  if( bNoInit )
    return;
  // Are we using a frozen covariance matrix? Or does it vary for each sample?
  bFreeze = SupportsUnfrozen() ? cl.GotSwitch( "freeze" ) : true;

  idxJackBoot = 0;
  const std::string sCovarSourceSwitch{ "covsrc" };
  bool bLoaded{};
  if( !cl.GotSwitch( sCovarSourceSwitch ) )
  {
    // Covariance source unspecified: Unfrozen - use binned data; Frozen - use bootstrap
    // TODO: Choose between the following two
    // Source = SupportsUnfrozen() ? SS::Binned : SS::Bootstrap;
    Source = SS::Bootstrap;
  }
  else
  {
    // There could be multiple options separated by '.'
    const std::string &sCovOptions{ cl.SwitchValue<std::string>( sCovarSourceSwitch ) };
    std::vector<std::string> vCovarSwitch{ MLU::ArrayFromString( sCovOptions, "." ) };
    static const std::string TooMany{ "Too many " };
    for( std::size_t optNum = 0; optNum < vCovarSwitch.size(); ++optNum )
    {
      // Decode the user's choice of covariance source
      const std::string sCovSrc{ MLU::ExtractToSeparator( vCovarSwitch[optNum] ) };
      if( MLU::EqualIgnoreCase( sCovSrc, "Reboot" ) )
      {
        if( optNum > 1 )
          throw std::runtime_error( "Too many " + sCovarSourceSwitch + " options " + sCovOptions );
        const int NumReplicas{ vCovarSwitch[optNum].empty()
          ? static_cast<int>( MLU::RandomCache::DefaultNumReplicas() )
          : MLU::FromString<int>( vCovarSwitch[optNum] ) };
        std::cout << "Covar Reboot " << NumReplicas << " replicas\n";
        for( std::unique_ptr<Sample> &f : dsrw.corr )
          f->Resample( 1, NumReplicas, idxJackBoot );
        Source = SS::Bootstrap;
        idxJackBoot = 1;
      }
      else if( optNum )
        throw std::runtime_error( "Too many " + sCovarSourceSwitch + " options " + sCovOptions );
      else if( MLU::EqualIgnoreCase( sCovSrc, "Rebin" ) )
      {
        // Rebin the raw data. Overwrites the raw data in the correlators
        dsrw.Rebin( MLU::ArrayFromString<int>( vCovarSwitch[optNum] ) );
        RebinSize = ds.RebinSize;
        std::cout << "Covar Rebin";
        for( int i : RebinSize )
          std::cout << MLU::Space << i;
        std::cout << MLU::NewLine;
        Source = SS::Binned;
        idxJackBoot = 1;
      }
      else if( MLU::EqualIgnoreCase( sCovSrc, "h5" ) )
      {
        // Load inverse covariance matrix from hdf5 - for debugging and comparison with other fit code
        std::vector<std::string> Opts{ MLU::ArrayFromString( vCovarSwitch[optNum] ) };
        if( Opts.size() < 2 || Opts.size() > 3 )
          throw std::runtime_error( "Options should contain file[,group],dataset. Bad options: " + vCovarSwitch[optNum] );
        std::string sRoot( "/" );
        const int iHaveGroupName{ Opts.size() == 2 ? 0 : 1 };
        std::string &GroupName{ iHaveGroupName ? Opts[1] : sRoot };
        std::string &DataSetName{ Opts[1 + iHaveGroupName] };
        ::H5::H5File f;
        ::H5::Group  g;
        MLU::H5::OpenFileGroup( f, g, Opts[0], "Loading covariance matrix from ", &GroupName );
        try
        {
          MLU::H5::ReadMatrix( g, DataSetName, Covar );
        }
        catch(const ::H5::Exception &)
        {
          ::H5::Exception::clearErrorStack();
          throw std::runtime_error( "Unable to load covariance matrix " + vCovarSwitch[optNum] );
        }
        Source = SS::Raw;
        bFreeze = true; // because we only have inv_cov. TODO: invert inv_cov if unfrozen? (Better to load cov)
        bLoaded = true;
      }
      else
      {
        // See what the user has asked for
        Source = MLU::FromString<SS>( sCovSrc );
        if( !vCovarSwitch[optNum].empty() )
          throw std::runtime_error( "Covariance source " + sCovSrc + " unexpected parameters: " + vCovarSwitch[optNum] );
      }
    }
  }
  dsrw.SetCovarSource( Source, idxJackBoot, bFreeze );

  // Check that all the input files match for raw/rebinned (checked for binned data on load)
  if( ( Source == SS::Raw || Source == SS::Binned ) && !bLoaded )
  {
    // Check whether all files have the same number of samples
    for( std::size_t i = 0; i < ds.corr.size(); ++i )
    {
      const int Count{ ds.corr[i]->NumSamples( Source, idxJackBoot ) };
      if( Count == 0 )
      {
        std::ostringstream os;
        os << Source << " samples not available for covariance matrix";
        if( i )
          os << " file " << ds.corr[i]->Name_.Filename;
        throw std::runtime_error( os.str().c_str() );
      }
      if( Count != ds.corr[0]->NumSamples( Source ) )
      {
        std::ostringstream os;
        os << "Can't use " << Source << " samples for covariance:" << MLU::NewLine
           << MLU::Space << ds.corr[0]->NumSamples( Source ) << " samples in "
           << ds.corr[0]->Name_.Filename << MLU::NewLine
           << MLU::Space << ds.corr[i]->NumSamples( Source ) << " samples in "
           << ds.corr[i]->Name_.Filename;
        throw std::runtime_error( os.str().c_str() );
      }
    }
  }

  // Has the user requested we get covariance using a bootstrap
  // TODO: Reinstate
  /*if( !cl.GotSwitch( "covboot" ) )
    CovarNumBoot = 0;
  else
  {
    CovarNumBoot = cl.SwitchValue<int>( "covboot" );
    if( CovarNumBoot < 0 || CovarNumBoot == 1 )
      throw std::runtime_error( "Rebootstrapping with " + std::to_string( CovarNumBoot ) + " replicas" );
    if( SourceIsBootstrap() )
    {
      std::ostringstream os;
      os << "Can't use ";
      if( Source == SS::Raw )
        os << "bootstrapped ";
      os << Source << " samples when bootstrapping covariance matrix";
      throw std::runtime_error( os.str().c_str() );
    }
    // 0 => Use the maximum number of samples available
    if( CovarNumBoot == 0 || CovarNumBoot > ds.MaxSamples )
      CovarNumBoot = ds.MaxSamples;
    // If I'm bootstrapping, I will need some random numbers
    if(   ds.corr[0].RandNum() // Do I have original bootstrap random numbers
       && ds.corr[0].SampleSize == CovarCount() // Are the random numbers over same range (0...SampleSize)
       && CovarNumBoot <= ds.corr[0].NumSamples() ) // Are there enough replicas available
    {
      // Re-use existing random numbers
      vCovarRandom.clear();
      CovarRandom.Map( ds.corr[0].RandNum(), CovarNumBoot, CovarCount() );
    }
    else
    {
      MLU::GenerateRandom( vCovarRandom, ds.corr[0].Seed_, CovarNumBoot, CovarCount() );
      CovarRandom.Map( vCovarRandom.data(), CovarNumBoot, CovarCount() );
    }
  }*/
}

std::ostream & operator<<( std::ostream &os, const CovarParams &cp )
{
  if( cp.bFreeze )
    os << "F";
  else
    os << "Unf";
  os << "rozen covariance matrix (mean from "
     << ( MLU::JackBootBase::UseCentralCovar() ? "binned" : "resampled" )
     << " data) constructed from " << cp.CovarCount() << MLU::Space;
  if( cp.CovarCount() != cp.CovarSampleSize() )
    os << "(" << cp.CovarSampleSize() << " independent) ";
  os << cp.Source << " samples";
  return os;
}
