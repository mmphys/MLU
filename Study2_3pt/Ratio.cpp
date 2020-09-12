/*************************************************************************************
 
 Create XML for a 3-pt current insertion study
 Initially use Z2 wall sources.
 Source file: xml3pt.cpp
 Copyright (C) 2020
 Author: Michael Marshall<Michael.Marshall@ed.ac.uk>
 
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

#include "../Analyze/Common.hpp"

static const std::vector<std::string> vHeavy{ "h0", "h1", "h2", "h3" };
static const int NumHeavy{ static_cast<int>( vHeavy.size() ) };
static const std::vector<int> vDeltaT{ 12, 14, 16, 20 };
static const int NumDeltaT{ static_cast<int>( vDeltaT.size() ) };
using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Model = Common::Model<Scalar>;

struct ModelInfo
{
  Model m;
  int idxE0;
};

class RatioMaker
{
public:
  const std::string &inBase;
  const std::string &outBase;
  const std::string &modelBase;
  const std::string gSnk;
  const std::string gSrc;
  const std::string Spectator;
  const Common::Momentum p;
protected:
  std::vector<ModelInfo> model;
  int MaxModelSamples;
  Common::SeedType Seed;
protected:
  std::string Filename3pt(const std::string &HeavySnk,const std::string &HeavySrc,const std::string &Current,int DeltaT)const;
  std::string FilenameRatio(int Ratio, const std::string &HeavySnk, const std::string &HeavySrc, int DeltaT ) const;
  std::string HeavyKey( const std::string &Heavy ) const;
public:
  RatioMaker( const std::string &inBase, const std::string &outBase, const std::string &modelBase,
              const std::string &gSnk, const std::string &gSrc, const std::string &Spectator,
              const Common::Momentum &p, const std::string &FitListName );
  void MakeRatios( int iSnk, int iSrc, int DeltaT );
};

std::string RatioMaker::Filename3pt( const std::string &HeavySnk, const std::string &HeavySrc, const std::string &Current, int DeltaT ) const
{
  std::string s{ inBase };
  s.append( "quark_" );
  s.append( HeavySnk );
  s.append( 1, '_' );
  s.append( HeavySrc );
  s.append( 1, '_' );
  s.append( Current );
  s.append( "_dt_" );
  s.append( std::to_string( DeltaT ) );
  s.append( "_p_" );
  s.append( p.to_string( Common::Underscore ) );
  s.append( 1, '_' );
  s.append( gSnk );
  s.append( 1, '_' );
  s.append( gSrc );
  s.append( ".fold." );
  s.append( std::to_string( Seed ) );
  s.append( 1, '.' );
  s.append( DEF_FMT );
  return s;
}

std::string RatioMaker::FilenameRatio(int Ratio, const std::string &HeavySnk, const std::string &HeavySrc, int DeltaT ) const
{
  std::string s{ outBase };
  s.append( 1, 'R' );
  s.append( std::to_string( Ratio ) );
  s.append( 1, '_' );
  s.append( HeavySnk );
  s.append( "_p_" );
  s.append( p.to_string( Common::Underscore ) );
  s.append( 1, '_' );
  s.append( HeavySrc );
  s.append( "_p_" );
  s.append( p.to_string( Common::Underscore ) );
  s.append( "_dt_" );
  s.append( std::to_string( DeltaT ) );
  s.append( 1, '_' );
  s.append( gSnk );
  s.append( 1, '_' );
  s.append( gSrc );
  s.append( ".fold." );
  s.append( std::to_string( Seed ) );
  s.append( 1, '.' );
  //s.append( DEF_FMT );
  return s;
}

std::string RatioMaker::HeavyKey( const std::string &Heavy ) const
{
  std::string s{ Heavy };
  s.append( 1, '_' );
  s.append( Spectator );
  s.append( "_p_" );
  s.append( p.to_string( Common::Underscore ) );
  return s;
}

RatioMaker::RatioMaker( const std::string &inBase_, const std::string &outBase_, const std::string &modelBase_,
                        const std::string &gSnk_, const std::string &gSrc_, const std::string &Spectator_,
                        const Common::Momentum &p_, const std::string &FitListName )
: inBase{inBase_}, outBase{outBase_}, modelBase{modelBase_}, gSnk{gSnk_}, gSrc{gSrc_}, Spectator{Spectator_}, p{ p_ },
  model( NumHeavy )
{
  // Map from case insensitive string to string
  using StringMapI= std::map<std::string, std::string, Common::LessCaseInsensitive>;
  // Here, SSRI = StringStringReaderInsensitive (not selective seratonin reuptake inhibitor)
  using SSRI= Common::KeyValReader<std::string, std::string, Common::LessCaseInsensitive>;

  if( p.p2() )
    throw std::runtime_error( "Only zero-momentum currently supported" );
  std::cout << "Reading fit list from " << FitListName << Common::NewLine;
  std::ifstream s( FitListName );
  if( !Common::FileExists( FitListName ) || s.bad() )
    throw std::runtime_error( "Fit list \"" + FitListName + "\" not found" );
  StringMapI FitList;
  FitList = SSRI::Read( s );
  for( int i = 0; i < NumHeavy; ++i )
  {
    std::string sH = HeavyKey( vHeavy[i] );
    auto it = FitList.find( sH );
    if( it == FitList.end() )
      throw std::runtime_error( "Fit " + sH + " not specified in FitListName" );
    std::string ModelFileName{ modelBase };
    ModelFileName.append( it->second );
    std::string MsgPrefix( 2, ' ' );
    MsgPrefix.append( sH );
    MsgPrefix.append( Common::Space );
    model[i].m.Read( ModelFileName, MsgPrefix.c_str() );
    model[i].idxE0 = model[i].m.GetColumnIndex( "E0" );
    if( i == 0 || MaxModelSamples > model[i].m.NumSamples() )
      MaxModelSamples = model[i].m.NumSamples();
    if( i )
      model[0].m.IsCompatible( model[i].m, nullptr, false );
    else
    {
      if( !model[0].m.Name_.bSeedNum )
        throw std::runtime_error( "Seed \"" + model[0].m.Name_.SeedString + "\" invalid" );
      Seed = model[0].m.Name_.Seed;
    }
  }
}

void RatioMaker::MakeRatios( int iSnk, int iSrc, int DeltaT )
{
  static const std::string Current{ "gT" };
  std::cout << vHeavy[iSnk] << Common::CommaSpace << gSnk << " <- " << vHeavy[iSrc] << Common::CommaSpace << gSrc
            << " dT=" << DeltaT << " p=" << p.to_string( Common::Space ) << Common::NewLine;
  // Load the correlators we need
  Common::DataSet<Scalar> ds;
  for( int i = 0; i < 2; i++ )
    for( int j = 0; j < 2; j++ )
    {
      Common::FileNameAtt fna{ Filename3pt( vHeavy[i ? iSnk : iSrc], vHeavy[j ? iSnk : iSrc], Current, DeltaT ) };
      ds.LoadCorrelator( std::move( fna ), false );
    }
  // Ensure the correlator includes timeslice DeltaT
  if( DeltaT >= ds.corr[0].Nt() )
    throw std::runtime_error( "DeltaT " + std::to_string( DeltaT ) + " >= Nt " + std::to_string( ds.corr[0].Nt() ) );
  // Make somewhere to put R2
  Common::Fold<Scalar> out( ds.NSamples, DeltaT + 1 );
  for( const Common::Fold<Scalar> &f : ds.corr )
    out.FileList.push_back( f.Name_.Filename );
  out.CopyAttributes( ds.corr[0] );
  out.NtUnfolded = ds.corr[0].Nt();
  Scalar * pDst{ out[Common::Fold<Scalar>::idxCentral] };
  // Now get pointers into the source correlators
  const int NumSrc{ static_cast<int>( ds.corr.size() ) };
  const Scalar * pSrc[NumSrc];
  for( int i = 0; i < NumSrc; ++i )
    pSrc[i] = ds.corr[i][Common::Fold<Scalar>::idxCentral];
  const Scalar * pE[2];
  pE[0] = model[iSrc].m[Common::Fold<Scalar>::idxCentral] + model[iSrc].idxE0;
  pE[1] = model[iSnk].m[Common::Fold<Scalar>::idxCentral] + model[iSnk].idxE0;
  for( int idx = Common::Fold<Scalar>::idxCentral; idx < ds.NSamples; ++idx )
  {
    const double EProd = (*pE[0]) * (*pE[1]);
    for( int t = 0; t <= DeltaT; ++t )
    {
      double z[4];
      for( int i = 0; i < 4; ++i )
        z[i] = *pSrc[i]++;
      const double d = 2 * std::sqrt( EProd * z[1] * z[2] / ( z[0] * z[3] ) );
      if( !std::isfinite( d ) )
        throw std::runtime_error( "R2 Overflow" );
      *pDst++ = d;
      //*pDst++ = 2 * std::sqrt( EProd * ( (*pSrc[1]++) * (*pSrc[2]++) ) / ( (*pSrc[0]++) * (*pSrc[3]++) ) );
    }
    for( int i = 0; i < NumSrc; ++i )
      pSrc[i] += ds.corr[i].Nt() - ( DeltaT + 1 );
    pE[0] += model[iSrc].m.Nt();
    pE[1] += model[iSnk].m.Nt();
  }
  out.MakeCorrSummary( nullptr );
  std::string OutFileName{ FilenameRatio( 2, vHeavy[iSnk], vHeavy[iSrc], DeltaT ) };
  const std::size_t FileNameLen{ OutFileName.length() };
  OutFileName.append( DEF_FMT );
  std::cout << "->" << OutFileName << Common::NewLine;
  out.Write( OutFileName, Common::sFold.c_str() );
  OutFileName.resize( FileNameLen );
  OutFileName.append( TEXT_EXT );
  out.WriteSummary( OutFileName );
}

static const char DefaultSnk[]{ "g5P" };
static const char *DefaultSrc{ DefaultSnk };
static const char DefaultSpectator[]{ "s" };

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"m", CL::SwitchType::Single, "" },
      {"s", CL::SwitchType::Single, DefaultSpectator },
      //{"n", CL::SwitchType::Single, "0"},
      {"snk", CL::SwitchType::Single, DefaultSnk },
      {"src", CL::SwitchType::Single, DefaultSrc },
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() == 1 )
    {
      // Read the list of fits I've chosen to use
      RatioMaker rm( cl.SwitchValue<std::string>("i"), cl.SwitchValue<std::string>("o"), cl.SwitchValue<std::string>("m"),
                     cl.SwitchValue<std::string>("snk"), cl.SwitchValue<std::string>("src"),
                     cl.SwitchValue<std::string>("s"), Common::Momentum( 0, 0, 0 ), cl.Args[0] );
      bShowUsage = false;
      //const int NSamples{ cl.SwitchValue<int>("n") };
      for( int iSnk = 0; iSnk < NumHeavy - 1; ++iSnk )
      {
        for( int iSrc = iSnk + 1; iSrc < NumHeavy; ++iSrc )
        {
          for( int DeltaT : vDeltaT )
          {
            rm.MakeRatios( iSnk, iSrc, DeltaT );
          }
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name << " <options> FitListFile\n"
    "Create 3pt ratios using E0 from fits in FitListFile. <options> are:\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-m     Model filename prefix\n"
    "-s     Spectator (Default: " << DefaultSpectator << ")\n"
    //"-n     Number of samples to fit, 0 (default) = all available from bootstrap\n"
    "Flags:\n"
    "--snk  Sink   operator (Default: " << DefaultSnk << ")\n"
    "--src  Source operator (Default: " << DefaultSrc << ")\n"
    "--help This message\n";
  }
  return iReturn;
}
