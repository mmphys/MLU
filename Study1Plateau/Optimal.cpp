//
//  Optimal.cpp
//  Utility
//
//  Created by Michael Marshall on 25/11/2019.
//  Copyright © 2019 sopa. All rights reserved.
//

#include <cmath>
//#include <iomanip>
#include <ios>
#include <iostream>
#include <stdio.h>

#include "../Analyze/Common.hpp"

using SampleD = Common::SampleD;
using SampleC = Common::SampleC;

int main( int argc, const char *argv[] )
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"ti", CL::SwitchType::Single, nullptr},
      {"tf", CL::SwitchType::Single, nullptr},
      {"dti", CL::SwitchType::Single, "1"},
      {"dtf", CL::SwitchType::Single, "1"},
      {"s", CL::SwitchType::Single, "0"},
      {"v", CL::SwitchType::Single, "0"},
      {"m", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"f", CL::SwitchType::Single, "0"},
      {"e", CL::SwitchType::Single, "0"},
      {"n", CL::SwitchType::Single, "0"},
      {"uncorr", CL::SwitchType::Flag, nullptr},
      {"savecorr", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles )
    {
      const int ti{ cl.SwitchValue<int>("ti") };
      const int tf{ cl.SwitchValue<int>("tf") };
      const int dti_max{ cl.SwitchValue<int>("dti") };
      const int dtf_max{ cl.SwitchValue<int>("dtf") };
      const int shift{ cl.SwitchValue<int>("s") };
      const int fold{ cl.SwitchValue<int>("f") };
      int NumExponents{ cl.SwitchValue<int>( "e" ) };
      std::string modelBaseFileName{ cl.SwitchValue<std::string>("m") };
      std::string outBaseFileName{ cl.SwitchValue<std::string>("o") };
      const int Verbosity{ cl.SwitchValue<int>("v") };
      const bool bSaveCorr{ cl.GotSwitch("savecorr") };
      const bool doCorr{ !cl.GotSwitch( "uncorr" ) };

      int NumOps = static_cast<int>( std::sqrt( static_cast<double>( NumFiles ) ) + 0.5 );
      if( NumOps * NumOps != NumFiles )
        throw std::invalid_argument( "Number of files should be a perfect square" );

      std::size_t i = 0;
      Common::SeedType Seed = 0;
      std::vector<std::string> OpNames;
      std::vector<Common::FileNameAtt> FileNames;
      std::string sSummaryBase{};
      for( const std::string &sFileName : cl.Args )
      {
        FileNames.emplace_back( sFileName, OpNames );
        if( i == 0 )
        {
          sSummaryBase = FileNames[0].Base;
          Seed = FileNames[0].Seed;
          outBaseFileName.append( sSummaryBase );
        }
        else
        {
          static const std::string sFile{ "File " };
          static const std::string sBad{ " doesn't match " };
          if( !Common::EqualIgnoreCase( FileNames[i].Base, FileNames[0].Base ) )
            throw std::runtime_error( sFile + std::to_string( i ) + " base " + FileNames[i].Base + sBad + FileNames[0].Base );
          if( !Common::EqualIgnoreCase( FileNames[i].Type, FileNames[0].Type ) )
            throw std::runtime_error( sFile + std::to_string( i ) + " type " + FileNames[i].Type + sBad + FileNames[0].Type );
          if( FileNames[i].Seed != FileNames[0].Seed )
            throw std::runtime_error( sFile + std::to_string( i ) + " seed " + std::to_string( FileNames[i].Seed ) + sBad + std::to_string( FileNames[0].Seed ) );
          if( !Common::EqualIgnoreCase( FileNames[i].Ext, FileNames[0].Ext ) )
            throw std::runtime_error( sFile + std::to_string( i ) + " extension " + FileNames[i].Ext + sBad + FileNames[0].Ext );
        }
        i++;
      }
      if( OpNames.size() != NumOps )
        throw std::runtime_error( std::to_string( OpNames.size() ) + " operators provided, but " + std::to_string( NumOps ) + " expected for " + std::to_string( NumFiles ) + " files" );
      for( int snk = 0; snk < NumOps; ++snk )
        for( int src = 0; src < NumOps; ++src )
        {
          Common::FileNameAtt &f{ FileNames[snk * NumOps + src] };
          if( f.op[0] != src || f.op[1] != snk )
            throw std::runtime_error( "Warning: Operator order should be sink-major, source minor" );
        }

      bShowUsage = false;
      if( NumExponents == 0 )
        NumExponents = NumOps;
      SampleD sd;
      {
        SampleC sc;
        sc.Read( cl.Args, fold, shift, NumOps, "  " );
        const int NumSamples{ sc.NumSamples() };
        const int Nt{ sc.Nt() };
        sd.resize( NumSamples, Nt );
        const std::complex<double> * pSrc = sc[SampleC::idxAux];
        double * pDst = sd[SampleD::idxAux];
        for( int i = SampleC::idxAux; i != NumSamples; i++ )
          for( int t = 0; t != Nt; t++ )
            *pDst++ = (pSrc++)->real();
      }
      sSummaryBase.append( 1, '.' );
      if( doCorr )
        sSummaryBase.append( "corr" );
      else
        sSummaryBase.append( "uncorr" );
      sSummaryBase.append( 1, '_' );
      sSummaryBase.append( std::to_string( ti ) );
      sSummaryBase.append( 1, '_' );
      sSummaryBase.append( std::to_string( tf ) );
      sSummaryBase.append( 1, '.' );
      std::string OutputRunBase{ outBaseFileName + sSummaryBase + "theta_"};
      sSummaryBase.append( OpNames[0] );
      for( int op = 1; op < NumOps; op++ )
      {
        sSummaryBase.append( 1, '_' );
        sSummaryBase.append( OpNames[op] );
      }
      const int NumFiles{ NumOps * NumOps };
      const int NSamples{ cl.SwitchValue<int>("n") <= 0 || cl.SwitchValue<int>("n") > sd.NumSamples() ? sd.NumSamples() : cl.SwitchValue<int>("n") };
      const int Nt      { sd.Nt() / NumFiles };
      const int NumParams{ NumExponents * ( 1 + NumOps ) };
      std::string sModelBase{ sSummaryBase };
      std::string GroupName;
      SampleD ModelParams;
      ModelParams.Read( Common::MakeFilename( sModelBase, Common::sModel, Seed, DEF_FMT ), GroupName );
      const double A_P1{ ModelParams[SampleD::idxCentral][4] };
      const double A_A1{ ModelParams[SampleD::idxCentral][5] };
      ModelParams.resize( 0, 0 );
      if( NumOps == 2)
      {  // I've only coded for two operators
        SampleC OptimalCorr{ NSamples, Nt };
        for( int degrees = -90; degrees <= 90; degrees ++ )
        {
          double costheta, sintheta;
          switch( ( degrees + 360 ) % 360 )
          {
            case 0:
              costheta = 1;
              sintheta = 0;
              break;
            case 90:
              costheta = 0;
              sintheta = 1;
              break;
            case 180:
              costheta = -1;
              sintheta = 0;
              break;
            case 270:
              costheta = 0;
              sintheta = -1;
              break;
            default:
            {
              const double theta{ M_PI * degrees / 180 };
              costheta = cos( theta );
              sintheta = sin( theta );
            }
          }
          const double cos_sq_theta{ costheta * costheta };
          const double sin_sq_theta{ sintheta * sintheta };
          const double cos_sin_theta{ costheta * sintheta };
          const typename SampleD::scalar_type * p = sd[SampleD::idxCentral];
          std::complex<double> * pDst = OptimalCorr[SampleC::idxCentral];
          for( int i = SampleD::idxCentral; i < NSamples; i++ )
          {
            const double Op_PP{ cos_sq_theta / ( A_P1 * A_P1 ) };
            const double Op_AP{ cos_sin_theta / ( A_P1 * A_A1 ) };
            const double Op_AA{ sin_sq_theta / ( A_A1 * A_A1 ) };
            for( int t = 0; t < Nt; t++, p++ )
              *pDst++ = Op_PP * p[0 * Nt] + Op_AA * p[3 * Nt] + Op_AP * ( p[1 * Nt] + p[2 * Nt] );
          }
          std::string sOptimusPrime{ OutputRunBase };
          sOptimusPrime.append( std::to_string(degrees) );
          if( bSaveCorr )
            OptimalCorr.Write(Common::MakeFilename(sOptimusPrime, Common::sBootstrap, Seed, DEF_FMT));
          Common::SummariseBootstrapCorr( OptimalCorr, sOptimusPrime, Seed );
        }
      }
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> Bootstrap1 [Bootstrap2 ...]\n"
    "Perform a multi-exponential fit of the specified bootstrap replicas, where <options> are:\n"
    "--ti   Initial fit time\n"
    "--tf   Final   fit time\n"
    "--dti  Number of initial fit times (default 1)\n"
    "--dti  Number of final   fit times (default 1)\n"
    "-s     time Shift (default 0)\n"
    "-v     Verbosity, 0 = normal (default), 1=debug, 2=save covar, 3=Every iteration\n"
    "-m     Model base filename\n"
    "-o     Output base filename\n"
    "-f     Fold the correlator, 0=don't fold (default), 1=+parity, -1=-parity\n"
    "-e     number of Exponents (default same as number of operators)\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "Flags:\n"
    "--uncorr   Uncorrelated fit (default correlated)\n"
    "--savecorr Save bootstrap replicas of correlators\n"
    "--help     This message\n";
  }
  return iReturn;
}
