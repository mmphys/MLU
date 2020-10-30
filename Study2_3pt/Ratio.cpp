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

//static const std::vector<std::string> vHeavy{ "h0", "h1", "h2", "h3" };
//static const int NumHeavy{ static_cast<int>( vHeavy.size() ) };
//static const std::vector<int> vDeltaT{ 12, 14, 16, 20 };
//static const int NumDeltaT{ static_cast<int>( vDeltaT.size() ) };
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
  //const std::string gSnk;
  //const std::string gSrc;
  const std::string Spectator;
  const Common::Momentum p;
protected:
  std::regex RegExExt;
  const bool RegExSwap;
  std::map<int, ModelInfo> model;
  int MaxModelSamples;
  Common::SeedType Seed;
protected:
  //std::string Filename3pt(const std::string &HeavySnk,const std::string &HeavySrc,const std::string &Current,int DeltaT)const;
  //std::string FilenameRatio(int Ratio, const std::string &HeavySnk, const std::string &HeavySrc, int DeltaT ) const;
  std::string HeavyKey( const std::string &Heavy ) const;
public:
  RatioMaker( const std::string &inBase, const std::string &outBase, const std::string &modelBase,
              //const std::string &gSnk, const std::string &gSrc,
              const std::string &Spectator, const Common::Momentum &p,
              std::regex RegExExt, const bool RegExSwap, const std::string &FitListName );
  //void MakeRatios( int iSnk, int iSrc, int DeltaT );
  void MakeRatios( std::string &sFileName );
};

/*std::string RatioMaker::Filename3pt( const std::string &HeavySnk, const std::string &HeavySrc, const std::string &Current, int DeltaT ) const
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
}*/

RatioMaker::RatioMaker( const std::string &inBase_, const std::string &outBase_, const std::string &modelBase_,
                        //const std::string &gSnk_, const std::string &gSrc_,
                        const std::string &Spectator_, const Common::Momentum &p_,
                        std::regex regExExt_, const bool regExSwap_, const std::string &fitListName_ )
: inBase{inBase_}, outBase{outBase_}, modelBase{modelBase_}, //gSnk{gSnk_}, gSrc{gSrc_},
  Spectator{Spectator_}, p{ p_ }, RegExExt{ regExExt_ }, RegExSwap{ regExSwap_ }
{
  // Here, SSRI = StringStringReaderInsensitive (not selective seratonin reuptake inhibitor)
  // using SSRI= Common::KeyValReader<std::string, std::string, Common::LessCaseInsensitive>;

  if( p.p2() )
    throw std::runtime_error( "Only zero-momentum currently supported" );
  const std::string FitListName{ modelBase + fitListName_ };
  std::cout << "Reading fit list from " << FitListName << Common::NewLine;
  std::ifstream s( FitListName );
  if( !Common::FileExists( FitListName ) || s.bad() )
    throw std::runtime_error( "Fit list \"" + FitListName + "\" not found" );
  std::map<int, std::string> FitList{ Common::KeyValReader<int, std::string>::Read( s ) };
  bool bFirst{ true };
  int idxFirst{ 0 };
  for( std::map<int, std::string>::iterator it = FitList.begin(); it != FitList.end(); ++it )
  {
    const int i{ it->first };
    const std::string sH{ std::to_string( i ) };//HeavyKey( vHeavy[i] ) };
    std::string ModelFileName{ modelBase };
    ModelFileName.append( it->second );
    std::string MsgPrefix( 2, ' ' );
    MsgPrefix.append( sH );
    MsgPrefix.append( Common::Space );
    ModelInfo &mi{ model[i] };
    mi.m.Read( ModelFileName, MsgPrefix.c_str() );
    mi.idxE0 = mi.m.GetColumnIndex( "E0" );
    if( bFirst || MaxModelSamples > mi.m.NumSamples() )
      MaxModelSamples = mi.m.NumSamples();
    if( bFirst )
    {
      if( !mi.m.Name_.bSeedNum )
        throw std::runtime_error( "Seed \"" + mi.m.Name_.SeedString + "\" invalid" );
      Seed = mi.m.Name_.Seed;
      idxFirst = i;
      bFirst = false;
    }
    else
      model[idxFirst].m.IsCompatible( mi.m, nullptr, false );
  }
}

void RatioMaker::MakeRatios( std::string &FileName )
{
  std::string Dir{ Common::ExtractDirPrefix( FileName ) };
  std::cout << "Processing " << Dir << FileName << Common::NewLine;
  std::smatch base_match;
  if( !std::regex_search( FileName, base_match, RegExExt ) || base_match.size() != 5 )
    throw std::runtime_error( "Can't extract sink/source from " + FileName );
  const int ModelSnk{ Common::FromString<int>( base_match[RegExSwap ? 4 : 2] ) };
  const int ModelSrc{ Common::FromString<int>( base_match[RegExSwap ? 2 : 4] ) };
  std::string sSnk{ base_match[RegExSwap ? 3 : 1] };
  sSnk.append( base_match[RegExSwap ? 4 : 2] );
  std::string sSrc{ base_match[RegExSwap ? 1 : 3] };
  sSrc.append( base_match[RegExSwap ? 2 : 4] );
  const ModelInfo &miSnk{ model.at( ModelSnk ) };
  const ModelInfo &miSrc{ model.at( ModelSrc ) };

  //int iSnk, int iSrc, int DeltaT;
  //static const std::string Current{ "gT" };
  //std::cout << vHeavy[iSnk] << Common::CommaSpace << gSnk << " <- " << vHeavy[iSrc] << Common::CommaSpace << gSrc
          // << " dT=" << DeltaT << " p=" << p.to_string( Common::Space ) << Common::NewLine;
  // Load the correlators we need
  Dir.append( base_match.prefix() );
  const std::size_t Len{ Dir.length() };
  int DeltaT{ 0 };
  Common::DataSet<Scalar> ds;
  for( int i = 0; i < 2; i++ )
    for( int j = 0; j < 2; j++ )
    {
      Dir.resize( Len );
      Dir.append( i ? sSnk : sSrc );
      Dir.append( 1, '_' );
      Dir.append( j ? sSnk : sSrc );
      Dir.append( base_match.suffix() );
      //Common::FileNameAtt fna{ Filename3pt( vHeavy[i ? iSnk : iSrc], vHeavy[j ? iSnk : iSrc], Current, DeltaT ) };
      Common::FileNameAtt fna{ Dir };
      if( i == 0 && j == 0 )
      {
        if( !fna.bGotDeltaT )
          throw std::runtime_error( "DeltaT missing" );
        DeltaT = fna.DeltaT;
      }
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
  pE[0] = miSrc.m[Common::Fold<Scalar>::idxCentral] + miSrc.idxE0;
  pE[1] = miSnk.m[Common::Fold<Scalar>::idxCentral] + miSnk.idxE0;
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
    pE[0] += miSrc.m.Nt();
    pE[1] += miSnk.m.Nt();
  }
  out.MakeCorrSummary( nullptr );
  std::string OutFileName{ outBase };
  OutFileName.append( "R2_" );
  OutFileName.append( sSnk );
  OutFileName.append( 1, '_' );
  OutFileName.append( sSrc );
  OutFileName.append( base_match.suffix() );
  std::cout << "->" << OutFileName << Common::NewLine;
  out.Write( OutFileName, Common::sFold.c_str() );
  const std::size_t FileNameLen{ OutFileName.find_last_of( '.' ) };
  if( FileNameLen == std::string::npos )
    OutFileName.append( 1, '_' );
  else
    OutFileName.resize( FileNameLen + 1 );
  OutFileName.append( TEXT_EXT );
  out.WriteSummary( OutFileName );
}

//static const char DefaultSnk[]{ "g5P" };
//static const char *DefaultSrc{ DefaultSnk };
static const char DefaultSpectator[]{ "s" };

int main(int argc, const char *argv[])
{
  static const char DefaultERE[]{ R"(([[:alpha:]]*)([[:digit:]]+)_([[:alpha:]]*)([[:digit:]]+))" };
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
      {"r", CL::SwitchType::Single, DefaultERE },
      //{"n", CL::SwitchType::Single, "0"},
      //{"snk", CL::SwitchType::Single, DefaultSnk },
      //{"src", CL::SwitchType::Single, DefaultSrc },
      {"w", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() > 1 )
    {
      // Read the list of fits I've chosen to use
      const std::string InBase{ cl.SwitchValue<std::string>("i") };
      RatioMaker rm( InBase, cl.SwitchValue<std::string>("o"), cl.SwitchValue<std::string>("m"),
                     //cl.SwitchValue<std::string>("snk"), cl.SwitchValue<std::string>("src"),
                     cl.SwitchValue<std::string>("s"), Common::Momentum( 0, 0, 0 ),
                    std::regex( cl.SwitchValue<std::string>("r"), std::regex::extended | std::regex::icase ),
                    cl.GotSwitch("w"),
                    cl.Args[0] );
      bShowUsage = false;
      std::vector<std::string> FileList{ Common::glob( ++cl.Args.begin(), cl.Args.end(), InBase.c_str() ) };
      std::size_t Count{ 0 };
      for( std::string &sFile : FileList )
      {
        try
        {
          rm.MakeRatios( sFile );
          ++Count;
        }
        catch(const std::exception &e)
        {
          std::cerr << "Error: " << e.what() << std::endl;
        }
        catch( ... )
        {
          std::cerr << "Error: Unknown exception" << std::endl;
        }
      }
      std::cout << Count << " ratios created" << Common::NewLine;
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name << " <options> FitList RatioFile1 [RatioFile2 [...]]\n"
    "Create 3pt ratios using E0 from fits in FitListFile. <options> are:\n"
    "-i     Input  filename prefix\n"
    "-o     Output filename prefix\n"
    "-m     Model filename prefix\n"
    "-s     Spectator (Default: " << DefaultSpectator << ")\n"
    "-r     Extended regex for sink/source type, default\n       " << DefaultERE << "\n"
    "       http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap09.html\n"
    "-w     Swap source / sink order in regex\n"
    //"-n     Number of samples to fit, 0 (default) = all available from bootstrap\n"
    "Flags:\n"
    //"--snk  Sink   operator (Default: " << DefaultSnk << ")\n"
    //"--src  Source operator (Default: " << DefaultSrc << ")\n"
    "--help This message\n";
  }
  return iReturn;
}
