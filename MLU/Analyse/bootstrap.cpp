/**

 Mike's lattice QCD utilities: Bootstrapper
 
 Source file: bootstrap.cpp
 
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

// Perform a bootstrap

#include <MLUconfig.h>
#include "bootstrap.hpp"

using Scalar = MLU::SampleC::value_type;
using Vector = MLU::Vector<Scalar>;
using Matrix = MLU::Matrix<Scalar>;

static const char szError[]{ "Error: " };
static const char szWarning[]{ "Warning: " };

// Show me the averages for each timeslice
/*void ShowTimeSliceAvg(const Latan::Dataset<Latan::DMat> &data) {
  const int nFile{static_cast<int>(data.size())};
  if( !nFile )
    std::cout << "No timeslices to average over!" << std::endl;
  else {
    const int nt{static_cast<int>(data[0].rows())};
    std::cout << "Averages of " << nFile << " files over " << nt << " timeslices:" << std::endl;
    std::vector<std::complex<double>> Avg(nt);
    for (unsigned int t = 0; t < nt; ++t) {
      Avg[t] = 0;
      for (unsigned int j = 0; j < nFile; ++j) {
        std::complex<double> d( data[j](t,0), data[j](t,1) );
        //std::cout << "d[" << j << "]=" << d << std::endl;
        Avg[t] += d;
        //std::cout << "a[" << j << "]=" << a << std::endl;
      }
      Avg[t] /= nFile;
      std::cout << "C(" << t << ")=" << Avg[t] << std::endl;
    }
  }
}*/

static const std::string sTimeRev{ "timerev" };

std::istream& operator>>(std::istream& is, BinOrder &binOrder )
{
  std::string s;
  if( is >> s )
  {
    if( MLU::EqualIgnoreCase( s, "Old" ) )
      binOrder = BinOrder::Old;
    else if( MLU::EqualIgnoreCase( s, "VeryOld" ) )
      binOrder = BinOrder::VeryOld;
    else if( MLU::EqualIgnoreCase( s, "Auto" ) )
      binOrder = BinOrder::Auto;
    else
      is.setstate( std::ios_base::failbit );
  }
  else
    is.setstate( std::ios_base::failbit );
  return is;
}

void TrajFile::Reverse( MLU::CorrelatorFileC &File ) const
{
  if( bTimeRev )
  {
    File.Name_.Filename.append( 1, ',' );
    File.Name_.Filename.append( sTimeRev );
    const int Nt{ File.Nt() };
    const int T0{ File.Timeslice() };
    const int DeltaT{ File.Name_.GotDeltaT() ? File.Name_.DeltaT[0] : 0 };
    using ValueT = typename MLU::SampleTraits<std::complex<double>>::value_type;
    std::vector<ValueT> Buffer( Nt );
    for( int row = 0; row < File.NumOps(); ++row )
    {
      ValueT * pCorr{ File[row] };
      for( int t = 0; t < Nt; ++t )
        Buffer[( T0 + t + Nt ) % Nt] = pCorr[( T0 + DeltaT - t + Nt ) % Nt];
      for( int t = 0; t < Nt; ++t )
        pCorr[t] = std::conj( Buffer[t] );
    }
  }
}

BootstrapParams::BootstrapParams( const MLU::CommandLine &cl )
: b2ptSymOp{ cl.GotSwitch( "symop" ) },
  b2ptSortZeroMom{ !cl.GotSwitch( "nosort" ) },
  bOverwrite{ cl.GotSwitch( "overwrite" ) },
  bVerboseSummaries{ !cl.GotSwitch( "terse" ) },
  TimesliceDetail{ cl.SwitchValue<int>( "t" ) },
  nSample{ cl.SwitchValue<int>( "n" ) },
  binSize{ cl.SwitchValue<int>( "b" ) },
  binAuto{ binSize == 0 },
  binOrder{ cl.SwitchValue<BinOrder>( "border" ) },
  seed{ MLU::RandomCache::DefaultSeed() },
  outStem{ cl.SwitchValue<std::string>( "o" ) },
  Ensemble{ cl.GotSwitch("e") ? cl.SwitchValue<std::string>("e") : "" }
{
  if( TimesliceDetail < 0 || TimesliceDetail > 2 )
    throw std::invalid_argument( "Timeslice detail " + std::to_string( TimesliceDetail ) + " invalid" );
  // Binning
  if( binSize < 0 )
    throw std::runtime_error( "Bin size must be positive if specified" );
  if( binAuto && binOrder != BinOrder::Auto )
    throw std::runtime_error( "Auto binning only works with Auto order" );
  if( MLU::GetHostName().empty() )
    throw std::invalid_argument( "Machine name can't be empty" );
  if( Ensemble.empty() )
    throw std::invalid_argument( "Ensemble name can't be empty" );
  MLU::MakeAncestorDirs( outStem );
}

void CopyTimeSlice( std::complex<double> *&pDst, const std::complex<double> *pSrc, int Nt, int TOffset )
{
  for( int t = 0; t < Nt; t++ )
    *pDst++  = pSrc[ ( t + TOffset ) % Nt ];
}

bool BootstrapParams::GatherInput( MLU::SampleC &out, const Iter &first, const Iter &last,
                        const TrajList &Traj, Algebra Snk, Algebra Src, bool bAlignTimeslices ) const
{
  // Factorised 2pt functions with Snk!=Src will have both entries loaded
  const int OpFactor{ ( !Traj.b3pt && b2ptSymOp && Snk != Src && Traj.OpSuffiiSame() ) ? 2 : 1 };
  // Count how many input records there are: in total; and per configuration
  int NumSamplesRaw{ 0 };
  std::map<int, int> ConfigCount;
  for( Iter it = first; it != last; ++it )
  {
    // How many entries do I need for this file
    int iNum{ Snk == Algebra::Spatial ? 0 : OpFactor };
    if( !iNum )
    {
      const MLU::Momentum &p{ it->Name_.GetFirstNonZeroMomentum().second };
      if( p.x )
        iNum++;
      if( p.y )
        iNum++;
      if( p.z )
        iNum++;
    }
    if( iNum )
    {
      const int ConfigNum{ static_cast<int>( it->Name_.Seed ) };
      auto itcc = ConfigCount.find( ConfigNum );
      if( itcc == ConfigCount.end() )
        ConfigCount.insert( { ConfigNum, iNum } );
      else
        itcc->second += iNum;
      NumSamplesRaw += iNum;
    }
  }
  if( !NumSamplesRaw )
    return false;

  out.ConfigCount.clear();
  out.ConfigCount.reserve( ConfigCount.size() );
  for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
    out.ConfigCount.emplace_back( it->first, it->second );
  ConfigCount.clear();
  std::size_t Row{};
  Matrix &mDst{ out.resizeRaw( NumSamplesRaw ) };

  // Now gather the raw data, aligning timeslices
  const int Nt{ out.Nt() };
  out.FileList.clear();
  out.FileList.reserve( NumSamplesRaw );
  for( Iter file = first; file != last; ++file )
  {
    const int TOffset{ bAlignTimeslices ? file->Timeslice() : 0 };
    std::string Filename{ file->Name_.Filename };
    const std::size_t FilenameLen{ Filename.length() };
    if( Snk == Algebra::Spatial )
    {
      const MLU::Momentum &p{ file->Name_.GetFirstNonZeroMomentum().second };
      if( p.x )
      {
        MLU::AppendGamma( Filename, Algebra::GammaX, MLU::Comma );
        out.FileList.emplace_back( Filename );
        const double Scale{ 1. / p.x };
        const std::complex<double> * const pSrc = (*file)( Algebra::GammaX, Src );
        for( int t = 0; t < Nt; t++ )
          mDst(Row,t) = pSrc[ ( t + TOffset ) % Nt ] * Scale;
        ++Row;
      }
      if( p.y )
      {
        Filename.resize( FilenameLen );
        MLU::AppendGamma( Filename, Algebra::GammaY, MLU::Comma );
        out.FileList.emplace_back( Filename );
        const double Scale{ 1. / p.y };
        const std::complex<double> * const pSrc = (*file)( Algebra::GammaY, Src );
        for( int t = 0; t < Nt; t++ )
          mDst(Row,t) = pSrc[ ( t + TOffset ) % Nt ] * Scale;
        ++Row;
      }
      if( p.z )
      {
        Filename.resize( FilenameLen );
        MLU::AppendGamma( Filename, Algebra::GammaZ, MLU::Comma );
        out.FileList.emplace_back( Filename );
        const double Scale{ 1. / p.z };
        const std::complex<double> * const pSrc = (*file)( Algebra::GammaZ, Src );
        for( int t = 0; t < Nt; t++ )
          mDst(Row,t) = pSrc[ ( t + TOffset ) % Nt ] * Scale;
        ++Row;
      }
    }
    else
    {
      for( int o = 0; o < OpFactor; o++ )
      {
        const Algebra ThisSnk( o ? Src : Snk );
        const Algebra ThisSrc( o ? Snk : Src );
        if( OpFactor > 1 )
        {
          Filename.resize( FilenameLen );
          MLU::AppendGamma( Filename, ThisSnk, MLU::Comma );
          MLU::AppendGamma( Filename, ThisSrc, MLU::Comma );
        }
        out.FileList.emplace_back( Filename );
        const std::complex<double> * const pSrc = (*file)( ThisSnk, ThisSrc );
        for( int t = 0; t < Nt; t++ )
          mDst(Row,t) = pSrc[ ( t + TOffset ) % Nt ];
        ++Row;
      }
    }
  }

  // Now Bin the data
  //if( binAuto && out.ConfigCount.size() > 1 ) out.Bin();
  //else out.Bin( ( binAuto ? 1 : binSize ) * OpFactor );
  if( binAuto )
    out.BinAuto();
  else
    out.BinFixed( binSize * OpFactor );
  return true;
}

int BootstrapParams::PerformBootstrap( const Iter &first, const Iter &last, const TrajList &Traj,
                                       const std::string &Suffix, bool bAlignTimeslices, bool bSaveBootstrap,
                                       bool bSaveSummaries, const std::vector<Algebra> &SinkAlgebra,
                                       const std::vector<Algebra> &SourceAlgebra ) const
{
  static constexpr int NumBlanks{2};
  static const std::string sBlanks( NumBlanks, ' ' );
  static const std::string sStar{ std::string( NumBlanks - 1, ' ' ) + "*" };
  const int NumFiles{ static_cast<int>( std::distance( first, last ) ) };
  if( NumFiles == 0 || ( !bSaveSummaries && !bSaveBootstrap ) )
    return 0;
  if( !binAuto && NumFiles % binSize )
    std::cout << szWarning << "last bin partially filled (" << ( NumFiles % binSize ) << " of " << binSize << ")\n";
  // Sort ... unless we want the old sort order (e.g. to check old results can be replicated)
  if( binOrder == BinOrder::Auto )
  {
    // sort all of the files by config, then timeslice
    std::sort( first, last, [&](const CorrFile &l, const CorrFile &r)
    {
      // Sort by trajectory
      if( l.Name_.Seed != r.Name_.Seed )
        return l.Name_.Seed < r.Name_.Seed;
      // Sort by timeslice
      if( l.bHasTimeslice != r.bHasTimeslice )
        return r.bHasTimeslice;
      if( l.bHasTimeslice && r.bHasTimeslice && l.Timeslice_ != r.Timeslice_ )
        return l.Timeslice_ < r.Timeslice_;
      // Sort by base filename
      return l.Name_.Base.compare( r.Name_.Base ) < 0;
    } );
  }
  int iCount{ 0 };
  // Make somewhere to put the raw data
  MLU::SampleC out( nSample, first->Nt() );
  out.Ensemble = Ensemble;
  for( int Snk = 0; Snk < SinkAlgebra.size(); Snk++ )
  {
    const int UpperLimit{ !Traj.b3pt && b2ptSymOp && Traj.OpSuffiiSame()
                          ? Snk + 1 : static_cast<int>( SourceAlgebra.size() ) };
    for( int Src = 0; Src < UpperLimit; Src++ )
    {
      // Now construct the base output filename
      static const char pszSep[] = "_";
      std::string sOutBase{ outStem };
      sOutBase.append( Traj.sShortPrefix );
      if( Traj.b3pt )
      {
        // The current insertion is what was at the contraction sink
        sOutBase.append( MLU::Gamma::NameShort( SinkAlgebra[Snk], pszSep, nullptr ) );
      }
      sOutBase.append( Traj.sShortSuffix );
      sOutBase.append( Suffix );
      sOutBase.append(MLU::Gamma::NameShort(Traj.b3pt ? (Traj.bRev ? SourceAlgebra[Src] : Traj.Alg3pt ):SinkAlgebra[Snk],
                                               pszSep, Traj.OpSuffixSnk.c_str() ));
      sOutBase.append(MLU::Gamma::NameShort(Traj.b3pt && Traj.bRev ? Traj.Alg3pt : SourceAlgebra[Src],
                                               pszSep, Traj.OpSuffixSrc.c_str() ));
      // Skip bootstrap if output exists
      const std::string sOutFile{ MLU::MakeFilename( sOutBase, MLU::sBootstrap, seed, DEF_FMT ) };
      const std::string sSummary{ MLU::MakeFilename( sOutBase, MLU::sBootstrap, seed, TEXT_EXT)};
      const bool bExistsFile{ MLU::FileExists( sOutFile ) };
      const bool bExistsSummary{ MLU::FileExists( sSummary ) };
      const bool bOutFileExists{ bExistsFile || ( !bSaveBootstrap && bExistsSummary ) };
      const std::string &sWarnFile{ bExistsFile ? sOutFile : sSummary };
      static const char szOver[]{ "overwriting " };
      if( bOutFileExists && !bOverwrite )
        std::cout << sStar << szError << szOver << sWarnFile << std::endl;
      else
      {
        if( bOutFileExists )
          std::cout << sStar << szWarning << szOver << sWarnFile << std::endl;
        // Count how many files there are for each config - needs to work regardless of sort
        if( GatherInput( out, first, last, Traj, SinkAlgebra[Snk], SourceAlgebra[Src], bAlignTimeslices ) )
        {
          std::cout << sBlanks << nSample << " samples to " << sOutFile << std::endl;
          out.Resample();
          // Now save the audit data for the bootstrap
          out.MakeCorrSummary();
          try {
          if( bSaveBootstrap )
            out.Write( sOutFile );
          if( bSaveSummaries )
            out.WriteSummary( sSummary, bVerboseSummaries );
          iCount++;
          }
          catch(const std::exception &e)
          {
            std::cout << szWarning << e.what() << std::endl;
          } catch( ... ) {
            std::cout << szError << "Cannot write - unknown exception" << std::endl;
          }
        }
      }
    }
  }
  return iCount;
}

int BootstrapParams::PerformBootstrap( vCorrFile &f, const TrajList &Traj, const std::vector<Algebra> &SinkAlgebra,
                                       const std::vector<Algebra> &SourceAlgebra ) const
{
  if( f.empty() )
    throw std::invalid_argument( "Can't perform a bootstrap without correlators" );
  int iCount{ 0 };
  // Now perform bootstraps for any individual timeslices
  Iter first = f.begin();
  Iter last  = f.end();
  // sort files by timeslice if we are doing timeslice summaries ... or old sort order asked for
  if( binOrder != BinOrder::Auto || TimesliceDetail > 0 )
  {
    std::sort( f.begin(), f.end(), [&](const CorrFile &l, const CorrFile &r)
    {
      // Sort first by timeslice
      if( l.bHasTimeslice != r.bHasTimeslice )
        return r.bHasTimeslice;
      if( l.bHasTimeslice && r.bHasTimeslice && l.Timeslice_ != r.Timeslice_ )
        return l.Timeslice_ < r.Timeslice_;
      // Very old sort order was config then filename
      if( binOrder == BinOrder::VeryOld )
      {
        if( l.Name_.Seed != r.Name_.Seed )
          return l.Name_.Seed < r.Name_.Seed;
        return l.Name_.Base.compare( r.Name_.Base ) < 0;
      }
      // More recent sort order is filename then config
      int iCompare = l.Name_.Base.compare( r.Name_.Base );
      if( iCompare )
        return iCompare < 0;
      return l.Name_.Seed < r.Name_.Seed;
    } );
  }
  if( TimesliceDetail > 0 )
  {
    // If more than one timeslice, save summary info for individual timeslices
    int Timeslice{ f[0].Timeslice() };
    Iter i = std::find_if( first + 1, last, [Timeslice]( const CorrFile &cf ) { return cf.Timeslice() != Timeslice; } );
    if( i != last )
    {
      while( first != last )
      {
        std::string PrefixT{ "_t_" };
        PrefixT.append( std::to_string( Timeslice ) );
        iCount += PerformBootstrap( first, i, Traj, PrefixT, false, TimesliceDetail > 1, true, SinkAlgebra, SourceAlgebra );
        first = i;
        if( first != last )
        {
          Timeslice = first->Timeslice();
          i = std::find_if( first + 1, last, [Timeslice]( const CorrFile &cf ) { return cf.Timeslice() != Timeslice; } );
        }
      }
      first = f.begin();
    }
  }
  // Now perform a single bootstrap, combining all the separate timeslices
  iCount += PerformBootstrap( first, last, Traj, "", true, true, true, SinkAlgebra, SourceAlgebra );
  return iCount;
}

/*struct Momentum {
  int x;
  int y;
  int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
};

struct SubStrings {
  std::string              Name;
  std::vector<std::string> Strings;
};

const std::string & MesonSuffix( const std::string MesonTypeName )
{
  static const std::string axial{ "ax" };
  static const std::string axial_name{ "_ax" };
  static const std::string other_name{};
  return axial.compare( MesonTypeName ) == 0 ? axial_name : other_name;
}*/

// Make a default manifest
// Heavy-light meson decays. 2pt function with current inserted and momentum on light meson

void BootstrapParams::Study1Bootstrap(StudySubject Study, const std::string &StudyPath,
                      const MLU::Momentum &mom, std::vector<Algebra> Alg,
                      const std::string &Heavy, const std::string &Light, bool Factorised ) const
{
  /*
  static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
  static const std::vector<std::string> algNames = { "g5", "gT5" };
  static const int NumAlgebra{ static_cast<int>( alg.size() ) };
  static const int NumCorr{ NumAlgebra * NumAlgebra };*/

  static constexpr unsigned int Nt{64};
  static constexpr unsigned int CStart{3000};
  static constexpr unsigned int CSkip{40};
  const unsigned int CCount{ Study == Z2 ? 10u : 1u };//10};
  const unsigned int CEnd{CStart + CSkip * CCount};
  static const std::string Dot{ "." };    // used inside filenames
  static const std::string H5{ Dot + DEF_FMT };    // used inside filenames
  static const std::string Sep{ "_" };    // used inside filenames
  static const std::string Space{ " " };  // whitespace as field separator / human readable info
  static const std::string Sink{ Sep + "sink" };    // used inside filenames
  static const std::string sProp{ "prop" + Sep };    // used inside filenames
  static const std::string sMom0{ Sep + "p" + Sep + MLU::p0.to_string3d( Sep ) };
  const std::string sMom{ Sep + "p" + Sep + mom.to_string3d( Sep )};
  const std::string sMomNeg{ Sep + "p" + Sep + mom.to_string3d( Sep, true )};
  //assert( alg.size() == algNames.size() && "Abbreviations should match gamma algebra" );
  const std::string CorrPrefix{ Heavy + Sep + Light + sMom };
  /*std::vector<std::string> CorrSuffixes( NumCorr );
  for( int iSrc = 0; iSrc < NumAlgebra; iSrc++ )
    for( int iSnk = 0; iSnk < NumAlgebra; iSnk++ ) {
      const int iCorr{ iSnk * NumAlgebra + iSrc };
      std::stringstream ss;
      ss << Sep << algNames[iSnk] << Sep << algNames[iSrc];
      CorrSuffixes[iCorr] = ss.str();
    }*/
  const int TimeSliceInc{ Study == Z2 ? 1 : 4 };
  const unsigned int NumEnds{ Factorised && !MLU::EqualIgnoreCase( Heavy, Light ) ? 2u : 1u };
  /*const std::size_t NumSamplesT{ NumEnds * CCount };
  const std::size_t NumSamples{ NumSamplesT * ( Nt / TimeSliceInc ) };
  std::vector<Latan::Dataset<Latan::DMat>> bsData( NumCorr );
  std::vector<Latan::Dataset<Latan::DMat>> bsDataT( NumCorr );
  std::vector<MLU::Correlator> buffer( NumCorr );
  for( int iCorr = 0 ; iCorr < NumCorr; iCorr++ ) {
    bsData[ iCorr ].resize( NumSamples );
    bsDataT[ iCorr ].resize( NumSamplesT );
    buffer[iCorr].resize( Nt );
  }
  int CorrIndex{ 0 };
  BootstrapParams par{ bsParams };
  par.bSaveSummaries = true;
  par.bSaveBootstrap = false;*/
  int CorrIndex{ 0 };
  vCorrFile InFiles( NumEnds * CCount * ( Nt / TimeSliceInc ) );
  for( int t = 0; t < Nt; t += TimeSliceInc ) {
    const std::string tPrefix{ Sep + "t" + ( Study == Z2 ? "" : Sep ) + std::to_string( t ) };
    //int CorrIndexT{ 0 };
    for( unsigned int iEnd = 0; iEnd < NumEnds; iEnd++ ) {
      const std::string & Left{ iEnd == 0 ? Heavy : Light };
      const std::string & Right{ iEnd == 0 ? Light : Heavy };
      for( int iConfig = CStart; iConfig < CEnd; iConfig += CSkip )
      {
        std::string sFileName{ StudyPath };
        switch( Study )
        {
          case Z2:
            sFileName.append( sProp + Left + tPrefix + sMom + Sep + sProp + Right + tPrefix + sMom0 + Sink + sMomNeg );
            break;
          default:
            sFileName.append( Left + Sep + Right + sMom + tPrefix );
            break;
        }
        sFileName.append( Dot + std::to_string( iConfig ) + H5 );
        /*MLU::ReadComplexArray(buffer, alg, sFileName, 0);
        for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
        {
          MLU::CopyCorrelator( bsData[iCorr][CorrIndex], buffer[iCorr], t );
          MLU::CopyCorrelator( bsDataT[iCorr][CorrIndexT], buffer[iCorr] );
        }
        CorrIndex++;
        CorrIndexT++;*/
        InFiles[CorrIndex++].Read( sFileName, Alg, Alg, &t );
      }
    }
    // If there's more than one configuration, perform a bootstrap of this timeslice
    //for( int iCorr = 0; CCount > 1 && iCorr < NumCorr; ++iCorr )
      //par.PerformBootstrap( bsDataT[iCorr], CorrPrefix + tPrefix + CorrSuffixes[iCorr] );
  }
  // Now perform a bootstrap overall
  /*par.bSaveBootstrap = true;
  for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
    par.PerformBootstrap( bsData[iCorr], CorrPrefix + CorrSuffixes[iCorr] );*/
  static const std::string opSuffixPoint{ "P" };
  static const std::string opSuffixWall{ "W" };
  const std::string *popSnk, *popSrc;
  switch( Study )
  {
    case GFWW:
      popSrc= &opSuffixWall;
      popSnk= &opSuffixPoint;
      break;
    case GFPW:
      popSrc= &opSuffixWall;
      popSnk= &opSuffixWall;
      break;
    default:
      popSrc= &opSuffixPoint;
      popSnk= &opSuffixPoint;
      break;
  }
  const std::string sShortPrefix{ Heavy + Sep + Light };
  const std::string sShortSuffix{ mom.p2_string( Sep ) };
  PerformBootstrap( InFiles, TrajList( sShortPrefix + sShortSuffix, sShortPrefix, sShortSuffix, *popSnk, *popSrc ),
                    Alg, Alg );
}

GroupMomenta Manifest::GetGroupP( const MLU::CommandLine &cl )
{
  GroupMomenta GroupP{ GroupMomenta::None };
  const bool p2{ cl.GotSwitch( "p2" ) };
  if( cl.GotSwitch( "pa" ) )
  {
    if( p2 )
      throw std::invalid_argument( "Can't group momenta by both p^2 and abs( p )" );
    GroupP = GroupMomenta::Abs;
  }
  else if( p2 )
    GroupP = GroupMomenta::Squared;
  return GroupP;
}

std::vector<Algebra> Manifest::GetCurrentAlgebra( const MLU::CommandLine &cl )
{
  bool bGammaSpatial{ false };
  AlgCurrentLoad.clear();
  AlgCurrentLoadNeg.clear();
  std::vector<Algebra> Alg3pt;
  if( cl.GotSwitch( "c" ) )
  {
    std::vector<MLU::NegateStar> NegRequested;
    Alg3pt = MLU::ArrayFromString<Algebra>( cl.SwitchValue<std::string>( "c" ), &NegRequested );
    MLU::NoDuplicates( Alg3pt, "Current algebra", 1 );
    // Copy algebra to be saved into algebra to be loaded, removing Spatial
    MLU::NegateStar bGammaSpatialNeg{ MLU::NegateStar::None };
    std::vector<bool> AlgGammaXYZFound( 3, false );
    for( int i = 0; i < Alg3pt.size(); ++i )
    {
      bool bCopy;
      switch( Alg3pt[i] )
      {
        case Algebra::Spatial:
          bGammaSpatial = true;
          bGammaSpatialNeg = NegRequested[i];
          bCopy = false;
          break;
        case Algebra::GammaX:
          AlgGammaXYZFound[0] = true;
          bCopy = true;
          break;
        case Algebra::GammaY:
          AlgGammaXYZFound[1] = true;
          bCopy = true;
          break;
        case Algebra::GammaZ:
          AlgGammaXYZFound[2] = true;
          bCopy = true;
          break;
        default:
          bCopy = true;
          break;
      }
      if( bCopy )
      {
        AlgCurrentLoad.push_back( Alg3pt[i] );
        AlgCurrentLoadNeg.push_back( NegRequested[i] );
      }
    }
    if( bGammaSpatial )
    {
      // Gamma Spatial has been asked for. Make sure: Spatial last
      // I.e. make sure the index of anything we're going to save is the same as anything we load
      Alg3pt = AlgCurrentLoad;
      Alg3pt.push_back( Algebra::Spatial );
      // Make sure: GammaX, Y and Z in algebra list to load
      for( int i = 0; i < 3; ++i )
      {
        if( !AlgGammaXYZFound[i] )
        {
          AlgCurrentLoad.push_back( i == 0 ? Algebra::GammaX : i == 1 ? Algebra::GammaY : Algebra::GammaZ );
          AlgCurrentLoadNeg.push_back( bGammaSpatialNeg );
        }
      }
    }
  }
  return Alg3pt;
}

Manifest::Manifest( const MLU::CommandLine &cl, const BootstrapParams &par_ )
/*Manifest::Manifest(const std::vector<std::string> &Args, const std::vector<std::string> &Ignore,
                   bool bSwapQuarks, const GroupMomenta GroupP, const std::vector<std::string> &vIgnoreMomenta,
                   std::regex *SSRegEx, bool bSwapSnkSrcRegEx )*/
: par{ par_ },
  bShowOnly{ cl.GotSwitch( "show" ) },
  bTimeRev3pt{ cl.GotSwitch( sTimeRev ) },
  bRevRev3pt{ cl.GotSwitch( "revrev" ) },
  GroupP{ GetGroupP( cl ) },
  InStem{ cl.SwitchValue<std::string>( "i" ) },
  DefaultGroup{ cl.SwitchValue<std::string>( "g" ) },
  DefaultDataSet{ cl.SwitchValue<std::string>( "d" ) },
  vIgnoreMomenta{ MLU::ArrayFromString( cl.SwitchValue<std::string>("pignore") ) },
  vIgnoreRegEx{ MLU::ArrayFromString( cl.SwitchValue<std::string>("ignore") ) },
  AlgSource{ MLU::ArrayFromString<Algebra>( cl.SwitchValue<std::string>( "a" ) ) },
  AlgCurrent{ GetCurrentAlgebra( cl ) }
{
  MLU::NoDuplicates( vIgnoreMomenta, "Ignored momenta", 0 );
  MLU::NoDuplicates( vIgnoreRegEx, "Ignored regular expressions", 0 );
  MLU::NoDuplicates( AlgSource, "Source algebra", 0 );
  for( Algebra a : AlgCurrent )
  {
    if( a == Algebra::Spatial )
    {
      if( GroupP != GroupMomenta::Squared )
        throw std::runtime_error( "Saving Gamma spatial only makes sense when grouping by p^2" );
      break;
    }
  }
  if( cl.GotSwitch( "ssre" ) )
  {
    const std::string &ssre{ cl.SwitchValue<std::string>( "ssre" ) };
    if( !ssre.empty() )
      SSRegEx.reset( new std::regex( ssre, std::regex::extended | std::regex::icase ) );
  }
  if( bTimeRev3pt && AlgSource.size() != 1 )
    throw std::runtime_error( "When --" + sTimeRev + " specified, only 1 source/sink algebra supported"
                              " and it must match the filename" );
  if( bRevRev3pt && !bTimeRev3pt )
    throw std::runtime_error( "--revrev is only valid with --" + sTimeRev );
}

void Manifest::MakeStudies( const std::vector<std::string> &Args )
{
  static const std::string qHeavy{ "h1" };
  static const std::string qLight{ "l" };
  static const MLU::Momentum StudyMomenta[] = {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 1, 1, 0 },
    { 1, 1, 1 },
    { 2, 0, 0 },
  };
  for( const std::string &s : Args )
  {
    const int iStudy{ MLU::FromString<int>( s ) };
    std::cout << "Study " << iStudy << ": ";
    const StudySubject Study{ static_cast<StudySubject>( iStudy ) };
    switch( Study )
    {
      case Z2:
      case GFWW:
      case GFPW:
        // Heavy-light semi-leptonics
        std::cout << "Making manifest for " << qHeavy << " and " << qLight << std::endl;
        for( const MLU::Momentum &m : StudyMomenta )
        {
          const bool Factorised{ Study != GFPW };
          par.Study1Bootstrap( Study, InStem, m, AlgSource, qHeavy, qLight, Factorised );
          if( !Factorised )
            par.Study1Bootstrap( Study, InStem, m, AlgSource, qLight, qHeavy, false );
          par.Study1Bootstrap( Study, InStem, m, AlgSource, qHeavy, qHeavy, false );
          par.Study1Bootstrap( Study, InStem, m, AlgSource, qLight, qLight, false );
        }
        break;
      default:
        std::cout << "undefined" << std::endl;
    }
  }
}

/*int Manifest::QuarkWeight( const char q ) const
{
  int iWeight{ std::toupper( q ) };
  switch( iWeight )
  {
    case 'H':
      iWeight = std::numeric_limits<int>::max();
      break;
    case 'L':
      iWeight = std::numeric_limits<int>::min();
      break;
    case 'S':
      iWeight = std::numeric_limits<int>::min() + 1;
      break;
  }
  return iWeight;
}*/

bool Manifest::NeedsTimeReverse( std::string &Contraction, MomentumMap &p, bool &bRev,
                                 std::string &OpSuffixSnk, std::string &OpSuffixSrc ) const
{
  bool bNeedsReverse{ false };
  if( p.size() == 2 )
  {
    std::vector<std::string> Parts( MLU::ArrayFromString( Contraction, MLU::Underscore ) );
    if( Parts.size() == 3 && !MLU::EqualIgnoreCase( Parts[1], Parts[2] ) )
    {
      /*const int Weight1{ QuarkWeight( Parts[1][0] ) };
      const int Weight2{ QuarkWeight( Parts[2][0] ) };
      bNeedsReverse = Weight1 < Weight2; // Should be other way around, but preserves GF wall-source*/
      bNeedsReverse = bRevRev3pt ? bRev : !bRev;
      if( bNeedsReverse )
      {
        bRev = !bRev;
        Contraction = Parts[0];
        Contraction.append( 1, '_' );
        Contraction.append( Parts[2] );
        Contraction.append( 1, '_' );
        Contraction.append( Parts[1] );
        MLU::Momentum &pOne{ p.begin()->second };
        MLU::Momentum &pTwo{ (++p.begin())->second };
        std::swap( pOne, pTwo );
        std::swap( OpSuffixSnk, OpSuffixSrc );
      }
    }
  }
  return bNeedsReverse;
}

void Manifest::BuildManifest( const std::vector<std::string> &Args, const std::vector<std::string> &Ignore )
{
  static const std::string Sep{ "_" };
  // Now walk the list of arguments.
  // Any file that's not in the ignore list gets added to the manifest
  bool parsed = true;
  for( const std::string &Filename : MLU::glob( Args.begin(), Args.end(), InStem.c_str() ) )
  {
    // See whether this file is in the ignore list
    std::size_t iIgnore = 0;
    while( iIgnore < Ignore.size() && Ignore[iIgnore].compare(Filename) )
      iIgnore++;
    if( iIgnore < Ignore.size() )
      std::cout << "Ignoring " << Filename << std::endl;
    else if( !MLU::FileExists(Filename))
    {
      parsed = false;
      std::cout << szError << Filename << " doesn't exist" << std::endl;
    }
    else
    {
      // Parse the name. Not expecting a type, so if present, put it back on the end of Base
      MLU::FileNameAtt Name_{ Filename, nullptr, &vIgnoreMomenta, &vIgnoreRegEx, true };
      if( !Name_.bSeedNum )
        throw std::runtime_error( "Contraction files must contain a configuration number" );
      if( Name_.Gamma.size() > 1 )
        throw std::runtime_error( "Multiple gamma insertions unsupported" );
      const bool b3pt{ !Name_.Gamma.empty() && Name_.GotDeltaT() };
      if( b3pt && AlgCurrent.empty() )
        throw std::runtime_error( "3pt functions present, but current insertion not specified" );
      if( !Name_.Gamma.empty() && !Name_.GotDeltaT() )
        std::cout << szWarning << "Gamma insertion without DeltaT. Possible 3pt function treated as 2pt" << std::endl;
      if( Name_.Gamma.empty() && Name_.GotDeltaT() )
        std::cout << szWarning << "DeltaT without gamma insertion. Possible 3pt function treated as 2pt" << std::endl;
      std::string Contraction{ Name_.GetBaseShortExtra() };
      // NB: bRev true means that the gamma in the name belongs with Q1
      bool bRev{ b3pt ? MLU::ExtractToken( Contraction, "[Rr][Ee][Vv]" ) : false };
      bool bOpSuffix{ false }; // Suffix for sink and source operators
      std::string OpSuffixSnk;
      std::string OpSuffixSrc;
      if( SSRegEx )
      {
        std::smatch base_match;
        if( std::regex_search( Contraction, base_match, *SSRegEx ) && base_match.size() == 3 )
        {
          bOpSuffix = true;
          static constexpr bool bSwapSnkSrcRegEx{ false }; // TODO: make this command-line argument
          OpSuffixSnk = base_match[bSwapSnkSrcRegEx ? 2 : 1];
          OpSuffixSrc = base_match[bSwapSnkSrcRegEx ? 1 : 2];
          const std::string Suffix{ base_match.suffix() };
          Contraction = base_match.prefix();
          Contraction.append( Suffix );
        }
      }
      if( !b3pt && par.b2ptSortZeroMom && !Name_.HasNonZeroMomentum() )
      {
        // Sort all the separate, underscore delimited component parts
        std::vector<std::string> vs{ MLU::ArrayFromString( Contraction, MLU::Underscore ) };
        Contraction.clear();
        std::sort(vs.begin(), vs.end());
        for( int i = 0; i < vs.size(); i++ )
        {
          if( i )
            Contraction.append(1, '_');
          Contraction.append( vs[i] );
        }
      }
      const bool bTimeRev{ b3pt && bTimeRev3pt && AlgSource.size() == 1 && Name_.Gamma[0] == AlgSource[0]
                           && NeedsTimeReverse( Contraction, Name_.p, bRev, OpSuffixSnk, OpSuffixSrc ) };
      // Add attributes back into contraction string
      std::string sShortPrefix{ Contraction }; // Prefix is everything before gamma
      std::string sShortSuffix;  // Everything after the gamma
      if( Name_.GotDeltaT() )
      {
        sShortSuffix = "_dt_";
        sShortSuffix.append( std::to_string( Name_.DeltaT[0] ) );
      }
      // Add the momenta into the correlator if they were present
      for( const MomentumMapValue & mmv : Name_.p )
      {
        if( GroupP == GroupMomenta::Squared )
        {
          sShortSuffix.append( mmv.second.p2_string( Sep, mmv.first ) );
        }
        else
        {
          sShortSuffix.append( Sep );
          sShortSuffix.append( mmv.first );
          sShortSuffix.append( Sep );
          if( GroupP == GroupMomenta::Abs )
            sShortSuffix.append( mmv.second.abs().to_string3d( Sep ) );
          else
            sShortSuffix.append( mmv.second.to_string3d( Sep ) );
        }
      }
      // The contraction is everything that makes this unique
      for( Algebra a : Name_.Gamma )
        Contraction.append( MLU::Gamma::NameShort( a, Sep.c_str() ) );
      Contraction.append( sShortSuffix );
      if( bOpSuffix )
      {
        Contraction.append( 1, '_' );
        Contraction.append( OpSuffixSnk );
        Contraction.append( OpSuffixSrc );
      }
      // Look for the contraction list this file belongs to
      auto itc = find( Contraction );
      if( itc == end() )
        itc = emplace( Contraction, TrajList( Contraction, sShortPrefix, sShortSuffix, OpSuffixSnk,
                                              OpSuffixSrc, b3pt, bRev,
                                              b3pt ? Name_.Gamma[0] : Algebra::MinusSigmaZT ) ).first;
      TrajList & cl{ itc->second };
      auto it = cl.FileInfo.find( Name_.Filename );
      if( it == cl.FileInfo.end() )
        cl.FileInfo.emplace( Name_.Filename, TrajFile( Name_.bGotTimeslice, Name_.Timeslice, Name_.p, Name_.Seed, bTimeRev, Name_.GotDeltaT() ? Name_.DeltaT[0] : 0 ) );
      else
        std::cout << "Ignoring repetition of " << Name_.Filename << std::endl;
    }
  }
  if( !parsed )
    throw std::runtime_error( "Remove non-existent files (or use '-x filename' to eXclude)" );
}

bool Manifest::RunManifest()
{
  // Walk the list of contractions, performing a separate bootstrap for each
  int BootstrapCount{ 0 };
  int CorrelatorCount{ 0 };
  if( bShowOnly )
    std::cout << "Contraction, Files, Configs, Config..., Count..." << MLU::NewLine;
  bool bOK{ true };
  for( auto itc = this->begin(); bOK && itc != this->end(); itc++ )
  {
    const std::string &Contraction{itc->first};
    const TrajList &l{itc->second};
    const unsigned int nFile{ static_cast<unsigned int>( l.FileInfo.size() ) };
    if( bShowOnly )
    {
      // Count each configuration
      using CCMap = std::map<int, int>;
      CCMap ConfigCount;
      for( auto it = l.FileInfo.begin(); it != l.FileInfo.end(); ++it )
      {
        const int ConfigNum{ static_cast<int>( it->second.Config ) };
        auto p = ConfigCount.find( ConfigNum );
        if( p == ConfigCount.end() )
          ConfigCount.insert( { ConfigNum, 1 } );
        else
          p->second++;
      }
      //Print a summary of each configuration
      static const std::string SepTab{ ",\t" };
      std::cout << Contraction << (l.bRev ? "[Rev]" : "") << SepTab << nFile << SepTab << ConfigCount.size();
      for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
        std::cout << MLU::Comma << it->first;
      for( auto it = ConfigCount.begin(); it != ConfigCount.end(); ++it )
        std::cout << MLU::Comma << it->second;
      std::cout << MLU::NewLine;
      //for( const auto &v : l.FileInfo )
        //std::cout << "  " << v.first << (v.second.bTimeRev ? " trev" : "") << MLU::NewLine;
    }
    else
    {
      std::cout << "Loading " << nFile << " files for " << Contraction << MLU::NewLine;
      //const bool b3pt{ !Alg3pt.empty() };
      std::vector<MLU::CorrelatorFileC> InFiles( nFile );
      unsigned int j{ 0 };
      for( auto it = l.FileInfo.begin(); it != l.FileInfo.end(); ++j, ++it )
      {
        const std::string &Filename{ it->first };
        const TrajFile &tf{ it->second };
        std::cout << "  t=" << tf.Timeslice << ( ( tf.bHasTimeslice && tf.Timeslice ) ? "->0" : "   " )
                  << '\t' << Filename << std::endl;
        std::string GroupName{ DefaultGroup };
        InFiles[j].Read( Filename, l.b3pt ? AlgCurrentLoad : AlgSource, AlgSource,
                         tf.bHasTimeslice ? &tf.Timeslice : nullptr, nullptr, &GroupName,
                         l.b3pt && tf.pFirstNonZero.IsNeg() && AlgCurrentLoadNeg.size() ? &AlgCurrentLoadNeg : nullptr, DefaultDataSet.c_str() );
        tf.Reverse( InFiles[j] );
      }
      try
      {
        CorrelatorCount += par.PerformBootstrap( InFiles, l, l.b3pt ? AlgCurrent : AlgSource, AlgSource );
        BootstrapCount++;
      }
      catch(const std::exception &e)
      {
        std::cerr << szError << e.what() << std::endl;
        bOK = false;
      }
    }
  }
  std::cout << CorrelatorCount << " bootstrap files written for " << BootstrapCount << " / " << this->size()
            << " correlators" << std::endl;
  return bOK;
}

/*****************************************************************

 Perform a bootstrap of all the files specified on the command line.
 Organise all the files by contraction, and sort them by trajectory.
 Ensure trajectories are processed in the same order so the
 bootstrap replicas for each contraction can be combined later during analysis

*****************************************************************/

int main(const int argc, const char *argv[])
{
  static const char DefaultERE[]{ R"(^([PWpw])([PWpw])_)" };
  static const char DefaultIgnore[]{ "[hH][iI][tT]_[^_]+" };
  static const char DefaultIgnoreMomenta[]{ "pq2" };
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"n", CL::SwitchType::Single, DEF_NSAMPLE},
      {"b", CL::SwitchType::Single, "0"},
      {"e", CL::SwitchType::Single, nullptr},
      {"border", CL::SwitchType::Single, "Auto"},
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"a", CL::SwitchType::Single, ""},
      {"c", CL::SwitchType::Single, nullptr},
      {"g", CL::SwitchType::Single, "" },
      {"d", CL::SwitchType::Single, "" },
      {"t", CL::SwitchType::Single, "0"},
      {"x", CL::SwitchType::Multiple, nullptr},
      {"symop", CL::SwitchType::Flag, nullptr},
      {"s", CL::SwitchType::Flag, nullptr},
      {"overwrite", CL::SwitchType::Flag, nullptr},
      {"p2",CL::SwitchType::Flag, nullptr},
      {"pa",CL::SwitchType::Flag, nullptr},
      {"ignore",CL::SwitchType::Single, DefaultIgnore},
      {"pignore",CL::SwitchType::Single, DefaultIgnoreMomenta},
      {"show", CL::SwitchType::Flag, nullptr},
      {"nosort", CL::SwitchType::Flag, nullptr},
      {sTimeRev.c_str(), CL::SwitchType::Flag, nullptr},
      {"revrev", CL::SwitchType::Flag, nullptr},
      {"ssre", CL::SwitchType::Single, DefaultERE },
      {"terse", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) )
    {
      BootstrapParams par( cl );
      Manifest Man{ cl, par };
      // If there are files specified on the command line,
      // parse the input file names, grouping by correlator, indexed by trajectory.
      if( !cl.Args.size() )
        throw std::runtime_error( "No work to do" );
      bShowUsage = false;
      if( cl.GotSwitch( "s" ) )
      {
        Man.MakeStudies( cl.Args );
      }
      else
      {
        Man.BuildManifest( cl.Args, cl.SwitchStrings( "x" ) );
        Man.RunManifest();
      }
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << szError << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << szError << "Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> ContractionFile1 [ContractionFile2 ...]\n"
    "Perform a bootstrap of the specified files, where <options> are:\n"
    "-n     Number of samples (" DEF_NSAMPLE ")\n"
    "-b     Bin size, or 0 (default)=auto (1 config=no binning, else 1 bin/config)\n"
    "-e     Ensemble name\n"
    "--border Bin Order: `Auto' (default)=config then timeslice then filename\n"
    "        `Old'=timeslice/filename/config, `VeryOld'=timeslice/config/filename\n"
    "-i     Input  prefix\n"
    "-o     Output prefix\n"
    "-a     list of gamma Algebras we're interested in at source (and sink for 2pt)\n"
    "-c     list of gamma algebras for current insertion         (Enable 3-pt mode)\n"
    "       Precede each by - and/or * to negate / conjugate negative momenta\n"
    "-g     Group name to read correlators from\n"
    "-d     DataSet name to read correlators from\n"
    "-t     timeslice detail 0 (none=default), 1 (.txt) or 2 (.txt+.h5)\n"
    "-x     eXclude file (may be repeated)\n"
    "Flags:\n"
    "-s     Perform bootstrap for specified study numbers\n"
    "--overwrite Overwrite existing files (default skip)\n"
    "--p2   group momenta by P^2\n"
    "--pa   group momenta by Abs( p )\n"
    "--pignore List of momenta to ignore (default: " << DefaultIgnoreMomenta << ")\n"
    "--ignore  List of regular expressions to ignore (default: " << DefaultIgnore << ")\n"
    "--nosort  Zero momentum 2pt functions of q1-q2 are grouped together with q2-q1\n"
    "          This option disables this (only affects zero momentum 2pt functions)\n"
    "--" << sTimeRev << " Fold 3pt functions together with time-reversed partner\n"
    "--revrev  When --" << sTimeRev << " in effect, reverses \"_rev\" correlators\n"
    "--show Show how files would be bootstrapped, but don't execute\n"
    "--ssre Sink / Source extended Regular Expression ('' to disable), default:\n       " << DefaultERE << "\n"
    "       http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap09.html\n"
    "--terse Terse summaries (no file list)\n"
    "--symop Symmetric operators e.g. g5-gT5 same as gT5-g5 (2pt only, experimental)\n"
    "--help This message\n"
    "Environment variables\n"
    "MLUCache   Prefix for location to save random numbers, or empty for no cache\n"
    "MLUHost    Host name to use with seed cache\n"
    "MLUNumBoot Default number of bootstrap replicas\n"
    "MLUSeed    Random number seed or text file to load from (default=random)\n";
  }
  return iReturn;
}
