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

#include <MLU/Common.hpp>

//static const std::vector<std::string> vHeavy{ "h0", "h1", "h2", "h3" };
//static const int NumHeavy{ static_cast<int>( vHeavy.size() ) };
//static const std::vector<int> vDeltaT{ 12, 14, 16, 20 };
//static const int NumDeltaT{ static_cast<int>( vDeltaT.size() ) };
using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Model = Common::Model<Scalar>;

using StringMapCI = std::multimap<std::string, std::string, Common::LessCaseInsensitive>;
using StringMapCIReader = Common::KeyValReaderMulti<std::string, std::string, Common::LessCaseInsensitive>;

using QP=std::pair<std::string, Common::Momentum>;
struct LessQP
{
  bool operator()( const QP &lhs, const QP &rhs ) const
  {
    int i = Common::CompareIgnoreCase( lhs.first, rhs.first );
    if( i )
      return i < 0;
    return lhs.second < rhs.second;
  }
};

static const Common::Momentum p0(0, 0, 0);

struct ModelInfo
{
  Model m;
  std::string FileName2pt;
  int idxE0;
};

class RatioMaker
{
public:
  const std::string &inBase;
  const std::string &C2Base;
  const std::string &outBase;
  const std::string &modelBase;
  //const std::string gSnk;
  //const std::string gSrc;
  const std::string Spectator;
  const Common::Momentum p;
protected:
  std::regex RegExExt;
  const bool RegExSwap;
  std::map<QP, ModelInfo, LessQP> model;
  int MaxModelSamples;
  Common::SeedType Seed;
protected:
  //std::string Filename3pt(const std::string &HeavySnk,const std::string &HeavySrc,const std::string &Current,int DeltaT)const;
  //std::string FilenameRatio(int Ratio, const std::string &HeavySnk, const std::string &HeavySrc, int DeltaT ) const;
  std::string HeavyKey( const std::string &Heavy ) const;
public:
  RatioMaker( const std::string &inBase, const std::string &C2Base,
              const std::string &outBase, const std::string &modelBase,
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

RatioMaker::RatioMaker( const std::string &inBase_, const std::string &C2Base_,
                        const std::string &outBase_, const std::string &modelBase_,
                        //const std::string &gSnk_, const std::string &gSrc_,
                        const std::string &Spectator_, const Common::Momentum &p_,
                        std::regex regExExt_, const bool regExSwap_, const std::string &fitListName_ )
: inBase{inBase_}, C2Base{C2Base_}, outBase{outBase_}, modelBase{modelBase_},
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
  StringMapCI FitList{ StringMapCIReader::Read( s ) };
  bool bFirst{ true };
  const ModelInfo *idxFirst{ nullptr };
  for( StringMapCI::iterator it = FitList.begin(); it != FitList.end(); ++it )
  {
    std::string ModelFileName{ modelBase };
    ModelFileName.append( it->second );
    Common::FileNameAtt fna( ModelFileName );
    const std::string &sH{ it->first };
    const QP qp( sH, fna.bGotMomentum ? fna.p : p0 );
    std::string MsgPrefix( 2, ' ' );
    MsgPrefix.append( sH );
    MsgPrefix.append( Common::Space );
    ModelInfo &mi{ model[qp] };
    {
      mi.FileName2pt = C2Base;
      std::vector<std::string> v2pt{ Common::ArrayFromString( fna.Base, Common::Period ) };
      const int Len{ static_cast<int>( v2pt.size() ) };
      if( Len < 3 )
        throw std::runtime_error( "Expected " + fna.Base + " to have at least 3 components" );
      for( int i = 0; i < Len - 2; ++i )
      {
        if( i )
          mi.FileName2pt.append( 1, '.' );
        mi.FileName2pt.append( v2pt[i] );
      }
    }
    mi.m.SetName( std::move( fna ) );
    mi.m.Read( MsgPrefix.c_str() );
    mi.idxE0 = mi.m.GetColumnIndex( "E0" );
    if( bFirst || MaxModelSamples > mi.m.NumSamples() )
      MaxModelSamples = mi.m.NumSamples();
    if( bFirst )
    {
      if( !mi.m.Name_.bSeedNum )
        throw std::runtime_error( "Seed \"" + mi.m.Name_.SeedString + "\" invalid" );
      Seed = mi.m.Name_.Seed;
      idxFirst = &mi;
      bFirst = false;
    }
    else
      idxFirst->m.IsCompatible( mi.m, nullptr, false );
  }
}

void RatioMaker::MakeRatios( std::string &FileName )
{
  std::cout << "Processing " << FileName << Common::NewLine;
  std::string Dir{ Common::ExtractDirPrefix( FileName ) };
  int DeltaT;
  Common::Gamma::Algebra gFrom;
  std::string SinkSourceOp( 1, '_' );
  {
    std::vector<std::string> OpNames;
    Common::FileNameAtt fna{ FileName, &OpNames };
    if( OpNames.empty() )
      throw std::runtime_error( "Sink / source perator names mising" );
    SinkSourceOp.append( OpNames[fna.op[1]] );
    SinkSourceOp.append( 1, '_' );
    SinkSourceOp.append( OpNames[fna.op[0]] );
    if( !fna.bGotDeltaT )
      throw std::runtime_error( "DeltaT missing" );
    DeltaT = fna.DeltaT;
    if( fna.Gamma.size() != 1 )
      throw std::runtime_error( FileName + " has " + std::to_string( fna.Gamma.size() ) + " currents" );
    gFrom = fna.Gamma[0];
  }
  Common::Momentum MomSnk, MomSrc;
  {
    std::string sCopy{ FileName };
    if( !MomSrc.Extract( sCopy ) )
      throw std::runtime_error( "No source momentum in " + FileName );
    if( !MomSnk.Extract( sCopy, Common::Momentum::SinkPrefix ) )
      throw std::runtime_error( "No sink momentum in " + FileName );
  }
  std::smatch base_match;
  if( !std::regex_search( FileName, base_match, RegExExt ) || base_match.size() != 4 )
    throw std::runtime_error( "Can't extract sink/source from " + FileName );
  const std::string sSnk = base_match[RegExSwap ? 3 : 2];
  const std::string sSrc = base_match[RegExSwap ? 2 : 3];
  const std::string FileNameSuffix{ base_match.suffix() };
  const ModelInfo &miSnk{ model.at( QP( sSnk, MomSnk ) ) };
  const ModelInfo &miSrc{ model.at( QP( sSrc, MomSrc ) ) };

  // Now read the two-point functions
  std::string C2NameSnk{ miSnk.FileName2pt };
  C2NameSnk.append( SinkSourceOp );
  static const std::string Spaces( 2, ' ' );
  Fold C2Snk( Common::MakeFilename( C2NameSnk, Common::sFold, Seed, DEF_FMT ), Spaces.c_str() );
  std::string C2NameSrc{ miSrc.FileName2pt };
  C2NameSrc.append( SinkSourceOp );
  Fold C2Src( Common::MakeFilename( C2NameSrc, Common::sFold, Seed, DEF_FMT ), Spaces.c_str() );

  //int iSnk, int iSrc, int DeltaT;
  //static const std::string Current{ "gT" };
  //std::cout << vHeavy[iSnk] << Common::CommaSpace << gSnk << " <- " << vHeavy[iSrc] << Common::CommaSpace << gSrc
          // << " dT=" << DeltaT << " p=" << p.to_string( Common::Space ) << Common::NewLine;
  // Load the correlators we need
  Dir.append( base_match.prefix() );
  Dir.append( base_match[1] );
  const std::size_t Len{ Dir.length() };
  Common::DataSet<Scalar> ds;
  for( int i = 0; i < 2; i++ )
    for( int j = 0; j < 2; j++ )
    {
      Dir.resize( Len );
      Dir.append( i ? sSnk : sSrc );
      Dir.append( 1, '_' );
      Dir.append( j ? sSnk : sSrc );
      Dir.append( FileNameSuffix );
      if( i == 0 )
        MomSrc.Replace( Dir, Common::Momentum::SinkPrefix, true );
      if( j == 1 )
        MomSnk.Replace( Dir, Common::Momentum::DefaultPrefix );
      if( i == j )
        Common::ReplaceGamma( Dir, gFrom, Common::Gamma::Algebra::GammaT );
      //Common::FileNameAtt fna{ Filename3pt( vHeavy[i ? iSnk : iSrc], vHeavy[j ? iSnk : iSrc], Current, DeltaT ) };
      ds.LoadCorrelator( Common::FileNameAtt( Dir ), false );
    }
  // Ensure the correlator includes timeslice DeltaT
  const int NTHalf{ ds.corr[0].Nt() / 2 };
  if( DeltaT >= NTHalf )
    throw std::runtime_error( "DeltaT " + std::to_string( DeltaT ) + " > Nt/2 " + std::to_string( NTHalf ) );
  // Make somewhere to put R2
  static constexpr int NumRatio{ 2 };
  std::vector<Common::Fold<Scalar>> out;
  Scalar * pDst[NumRatio];
  for( int i = 0; i < 2; i++ )
  {
    out.emplace_back( ds.NSamples, DeltaT + 1 );
    for( const Common::Fold<Scalar> &f : ds.corr )
      out[i].FileList.push_back( f.Name_.Filename );
    out[i].CopyAttributes( ds.corr[0] );
    out[i].NtUnfolded = ds.corr[0].Nt();
    pDst[i] = out[i][Common::Fold<Scalar>::idxCentral];
  }
  // Now get pointers into the source correlators
  const int NumSrc{ static_cast<int>( ds.corr.size() ) };
  const Scalar * pSrc[NumSrc];
  for( int i = 0; i < NumSrc; ++i )
    pSrc[i] = ds.corr[i][Common::Fold<Scalar>::idxCentral];
  const Scalar * pE[2];
  pE[0] = miSrc.m[Common::Fold<Scalar>::idxCentral] + miSrc.idxE0;
  pE[1] = miSnk.m[Common::Fold<Scalar>::idxCentral] + miSnk.idxE0;
  static constexpr int NumC2{ 2 };
  const Scalar * pC2[NumC2];
  pC2[0] = C2Src[Common::Fold<Scalar>::idxCentral];
  pC2[1] = C2Snk[Common::Fold<Scalar>::idxCentral];
  for( int idx = Common::Fold<Scalar>::idxCentral; idx < ds.NSamples; ++idx )
  {
    const double EProd = (*pE[0]) * (*pE[1]);
    double C2[NumC2];
    for( int i = 0; i < NumC2; ++i )
      C2[i] = pC2[i][DeltaT] - 0.5 * pC2[i][NTHalf] * std::exp( - pE[i][0] * ( NTHalf - DeltaT ) );
    const double C2Prod = std::abs( C2[0] * C2[1] );
    for( int t = 0; t <= DeltaT; ++t )
    {
      double z[4];
      for( int i = 0; i < 4; ++i )
        z[i] = *pSrc[i]++;
      const double n = std::abs( EProd * z[1] * z[2] );
      if( !std::isfinite( n ) )
        throw std::runtime_error( "R2 Overflow" );
      *pDst[0]++ = 4 * std::sqrt( n / C2Prod ); // R1
      *pDst[1]++ = 2 * std::sqrt( n / std::abs( z[0] * z[3] ) ); // R2
    }
    for( int i = 0; i < NumSrc; ++i )
      pSrc[i] += ds.corr[i].Nt() - ( DeltaT + 1 );
    pE[0] += miSrc.m.Nt();
    pE[1] += miSnk.m.Nt();
    pC2[0] += C2Src.Nt();
    pC2[1] += C2Snk.Nt();
  }
  for( int i = 0; i < 2; i++ )
  {
    out[i].MakeCorrSummary( nullptr );
    std::string OutFileName{ outBase };
    OutFileName.append( 1, 'R' );
    OutFileName.append( 1, '1' + i );
    OutFileName.append( 1, '_' );
    OutFileName.append( sSnk );
    OutFileName.append( 1, '_' );
    OutFileName.append( sSrc );
    OutFileName.append( base_match.suffix() );
    std::cout << "->" << OutFileName << Common::NewLine;
    out[i].Write( OutFileName, Common::sFold.c_str() );
    const std::size_t FileNameLen{ OutFileName.find_last_of( '.' ) };
    if( FileNameLen == std::string::npos )
      OutFileName.append( 1, '_' );
    else
      OutFileName.resize( FileNameLen + 1 );
    OutFileName.append( TEXT_EXT );
    out[i].WriteSummary( OutFileName );
  }
}

//static const char DefaultSnk[]{ "g5P" };
//static const char *DefaultSrc{ DefaultSnk };
static const char DefaultSpectator[]{ "s" };

int main(int argc, const char *argv[])
{
  static const char DefaultERE[]{ R"((quark_|anti_)([[:alpha:]]+[[:digit:]]*)_([[:alpha:]]+[[:digit:]]*))" };
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i", CL::SwitchType::Single, "" },
      {"c", CL::SwitchType::Single, "" },
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
      RatioMaker rm( InBase, cl.SwitchValue<std::string>("c"),
                     cl.SwitchValue<std::string>("o"), cl.SwitchValue<std::string>("m"),
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
    "-i     Input3 filename prefix\n"
    "-c     Input2 filename prefix\n"
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
