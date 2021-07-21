/*************************************************************************************
 
 Create Ratios, e.g. R1, R2 and associated values such as Z_V
 Source file: Ratio.cpp
 Copyright (C) 2020-2021
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

#include "Ratio.hpp"

/*namespace Common {
template<std::size_t Num>
std::istream & operator>>( std::istream &is, std::array<std::string,Num> &c )
{
  for( std::size_t i = 0; i < c.size(); ++i )
    if( ! ( is >> c[i] ) )
      throw std::runtime_error( "Error reading array of " + std::to_string(Num) + " strings" );
  return is;
}
}*/

struct SnkSrcString
{
  int         DeltaT;
  std::string Snk;
  std::string Src;
  std::string s;
};

std::istream & operator>>( std::istream &is, SnkSrcString &sss )
{
  if( is >> sss.DeltaT >> sss.Snk >> sss.Src >> sss.s )
    return is;
  throw std::runtime_error( "Error reading SnkSrcString" );
}

QPMapModelInfo::QPMapModelInfo( const std::string &FitListName, const std::string &modelBase, const std::string &C2Base )
{
  std::map<QP, ModelInfo, LessQP> &model( *this );
  std::cout << "Reading 2pt fit list from " << FitListName << Common::NewLine;
  std::ifstream s( FitListName );
  if( !Common::FileExists( FitListName ) || s.bad() )
    throw std::runtime_error( "2pt fit file \"" + FitListName + "\" not found" );
  using StringMapCI = std::multimap<std::string, std::string, Common::LessCaseInsensitive>;
  using StringMapCIReader = Common::KeyValReaderMulti<std::string, std::string, Common::LessCaseInsensitive>;
  StringMapCI FitList{ StringMapCIReader::Read( s ) };
  const ModelInfo *idxFirst{ nullptr };
  for( StringMapCI::iterator it = FitList.begin(); it != FitList.end(); ++it )
  {
    std::string ModelFileName{ modelBase };
    ModelFileName.append( it->second );
    Common::FileNameAtt fna( ModelFileName );
    const std::string &sH{ it->first };
    const QP qp( sH, fna.p.empty() ? Common::p0 : fna.p.begin()->second );
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
    if( !idxFirst )
    {
      idxFirst = &mi;
      if( !mi.m.Name_.bSeedNum )
        throw std::runtime_error( "Seed \"" + mi.m.Name_.SeedString + "\" invalid" );
      Seed = mi.m.Name_.Seed;
      MaxModelSamples = mi.m.NumSamples();
    }
    else
    {
      idxFirst->m.IsCompatible( mi.m, &MaxModelSamples, Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_CONFIG_COUNT );
    }
  }
}

QDTMapModelInfo::QDTMapModelInfo( const std::string &FitListName, const std::string &modelBase, const std::string &C2Base )
{
  std::map<QDT, ModelInfo, LessQDT> &model( *this );
  std::cout << "Reading Z_V fit list from " << FitListName << Common::NewLine;
  std::ifstream s( FitListName );
  if( !Common::FileExists( FitListName ) || s.bad() )
    throw std::runtime_error( "Z_V fit file \"" + FitListName + "\" not found" );
  using StringMapCI = std::multimap<std::string, SnkSrcString, Common::LessCaseInsensitive>;
  using StringMapCIReader = Common::KeyValReaderMulti<std::string, SnkSrcString, Common::LessCaseInsensitive>;
  StringMapCI FitList{ StringMapCIReader::Read( s ) };
  const ModelInfo *idxFirst{ nullptr };
  for( StringMapCI::iterator it = FitList.begin(); it != FitList.end(); ++it )
  {
    SnkSrcString & sss{ it->second };
    for( const std::string &ModelFileName : Common::glob( &sss.s, &sss.s + 1, modelBase.c_str() ) ) {
    std::vector<std::string> OpNames;
    Common::FileNameAtt fna( ModelFileName, &OpNames );
    //if( !fna.bGotDeltaT )
      //throw std::runtime_error( "DeltaT missing from " + FitListName );
    if( fna.op.size() != 2 )
      throw std::runtime_error( "Expected 2 operator names, but " + std::to_string(fna.op.size()) + " provided" );
    const std::string &sH{ it->first };
    const QDT qdt( sH, sss.DeltaT, sss.Snk, sss.Src );
    std::string MsgPrefix( 2, ' ' );
    MsgPrefix.append( sH );
    MsgPrefix.append( Common::Space );
    ModelInfo &mi{ model[qdt] };
    mi.m.SetName( std::move( fna ) );
    mi.m.Read( MsgPrefix.c_str() );
    mi.idxE0 = mi.m.GetColumnIndex( "E0" );
    if( !idxFirst )
    {
      idxFirst = &mi;
      if( !mi.m.Name_.bSeedNum )
        throw std::runtime_error( "Seed \"" + mi.m.Name_.SeedString + "\" invalid" );
      Seed = mi.m.Name_.Seed;
      MaxModelSamples = mi.m.NumSamples();
    }
    else
    {
      idxFirst->m.IsCompatible( mi.m, &MaxModelSamples, Common::COMPAT_DISABLE_BASE );
    }
  } }
}

Maker * Maker::Make( const std::string &Type, std::string &TypeParams,
                     const std::string &inBase, const std::string &C2Base,
                     const std::string &modelBase,const std::string &outBase,
                     std::regex RegExExt, const bool RegExSwap,
                     const Freeze fEnergy, const Freeze fZV, const bool bSymmetrise, const std::string &FitListName )
{
  Maker * r;
  if( !Common::CompareIgnoreCase( Type, "R1R2" ) )
    r = new R1R2Maker( TypeParams, inBase, C2Base, modelBase, outBase, RegExExt, RegExSwap, fEnergy, fZV,
                       bSymmetrise, FitListName );
  else if( !Common::CompareIgnoreCase( Type, "ZV" ) )
    r = new ZVMaker( TypeParams, inBase, C2Base, modelBase, outBase, RegExExt, RegExSwap, fEnergy, fZV,
                     bSymmetrise, FitListName );
  else
    throw std::runtime_error( "I don't know how to make type " + Type );
  if( !TypeParams.empty() )
    throw std::runtime_error( "Unused parameters: " + TypeParams );
  return r;
}

// Incoming file should be a three-point function with a heavy sink and light(er) source

void Maker::Make( std::string &FileName )
{
  std::cout << "Processing " << FileName << Common::NewLine;
  std::vector<std::string> OpNames;
  Common::FileNameAtt fna{ FileName, &OpNames };
  if( OpNames.empty() )
    throw std::runtime_error( "Sink / source operator names mising" );
  if( !fna.bGotDeltaT )
    throw std::runtime_error( "DeltaT missing" );
  if( fna.Gamma.size() != 1 )
    throw std::runtime_error( FileName + " has " + std::to_string( fna.Gamma.size() ) + " currents" );
  const Common::FileNameMomentum &MomSrc{ fna.GetMomentum() };
  const Common::FileNameMomentum &MomSnk{ fna.GetMomentum( Common::Momentum::SinkPrefix )  };
  std::smatch base_match;
  if( !std::regex_search( fna.BaseShort, base_match, RegExExt ) || base_match.size() != 4 )
    throw std::runtime_error( "Can't extract sink/source from " + fna.BaseShort );
  const std::string qSnk = base_match[RegExSwap ? 3 : 2];
  const std::string qSrc = base_match[RegExSwap ? 2 : 3];
  const std::string FileNameSuffix{ base_match.suffix() };
  std::string FileNamePrefix{ base_match.prefix() };
  FileNamePrefix.append( base_match[1] );
  fna.BaseShort = FileNamePrefix;
  const ModelInfo &miSnk{ model.at( QP( qSnk, MomSnk ) ) };
  const ModelInfo &miSrc{ model.at( QP( qSrc, MomSrc ) ) };
  Make( fna, FileNameSuffix, qSnk, qSrc, miSnk, miSrc, OpNames[fna.op[1]], OpNames[fna.op[0]] );
}

void ZVMaker::Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                    const std::string &qSnk, const std::string &qSrc,
                    const ModelInfo &miSnk, const ModelInfo &miSrc, const std::string &opSnk, const std::string &opSrc )
{
  const Common::FileNameMomentum &MomSrc{ fna.GetMomentum() };
  const Common::FileNameMomentum &MomSnk{ fna.GetMomentum( Common::Momentum::SinkPrefix )  };
  // Complain about errors
  {
    std::ostringstream os;
    if( Common::CompareIgnoreCase( qSnk, qSrc ) )
      os << "qSnk " << qSnk << " != qSrc " << qSrc;
    else if( fna.Gamma[0] != Common::Gamma::Algebra::GammaT )
      os << "Gamma " << fna.Gamma[0] << " != " << Common::Gamma::Algebra::GammaT;
    else if( MomSnk )
      os << "Sink momentum " << MomSnk;
    else if( MomSrc )
      os << "Source momentum " << MomSrc;
    if( os.tellp() != std::streampos( 0 ) )
    {
      os << " file " << fna.Filename;
      throw std::runtime_error( os.str() );
    }
  }

  // Read the three-point correlator
  static const std::string Spaces( 2, ' ' );
  Fold C3;
  C3.Read( fna.Filename, Spaces.c_str() );

  // Ensure the correlator includes timeslice DeltaT
  const int NTHalf{ C3.Nt() / 2 };
  if( fna.DeltaT > NTHalf )
    throw std::runtime_error( "DeltaT " + std::to_string( fna.DeltaT ) + " > Nt/2 " + std::to_string( NTHalf ) );

  // Now read the two-point function
  std::string C2Name{ miSnk.FileName2pt };
  AppendOps( C2Name, opSnk, opSrc );
  C2Name = Common::MakeFilename( C2Name, Common::sFold, model.Seed, DEF_FMT );
  Fold C2;
  C2.Read( C2Name, Spaces.c_str() );
  int NumSamples {};
  C3.IsCompatible( C2, &NumSamples, Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT
                                    | Common::COMPAT_DISABLE_CONFIG_COUNT );
  if( C2.NtUnfolded != C3.Nt() )
  {
    std::ostringstream os;
    os << "C2.NtUnfolded " << C2.NtUnfolded << " != C3.Nt() " << C3.Nt();
    throw std::runtime_error( os.str() );
  }

  // Make somewhere to put Z_V
  Fold out( NumSamples, fna.DeltaT + 1 );
  out.FileList.push_back( C3.Name_.Filename );
  out.FileList.push_back( C2.Name_.Filename );
  out.FileList.push_back( miSrc.m.Name_.Filename );
  out.CopyAttributes( C3 );
  out.NtUnfolded = C3.Nt();
  Scalar * pDst{ out[Fold::idxCentral] };

  const Scalar * pE{ miSnk.m[Fold::idxCentral] + miSnk.idxE0 };
  const Scalar * pC2 { C2[Fold::idxCentral] };
  const Scalar * pSrc{ C3[Fold::idxCentral] };
  for( int idx = Fold::idxCentral; idx < NumSamples; ++idx )
  {
    const double CTilde{ pC2[fna.DeltaT] - 0.5 * pC2[NTHalf] * std::exp( - pE[0] * ( NTHalf - fna.DeltaT ) ) };
    for( int t = 0; t <= fna.DeltaT; ++t )
    {
      const double z{ CTilde / *pSrc++ };
      if( !std::isfinite( z ) )
        throw std::runtime_error( "ZV Overflow" );
      *pDst++ = z;
    }
    pSrc += C3.Nt() - ( fna.DeltaT + 1 );
    if( fEnergy != Freeze::Central )
      pE += miSnk.m.Nt();
    pC2 += C2.Nt();
  }
  // Make output file
  out.MakeCorrSummary( nullptr );
  std::string OutFileName{ outBase };
  OutFileName.append( "ZV_" );
  OutFileName.append( qSnk );
  //OutFileName.append( Suffix );
  Common::AppendDeltaT( OutFileName, fna.DeltaT );
  AppendOps( OutFileName, opSnk, opSrc );
  OutFileName = Common::MakeFilename( OutFileName, Common::sFold, model.Seed, DEF_FMT );
  std::cout << "->" << OutFileName << Common::NewLine;
  out.Write( OutFileName, Common::sFold.c_str() );
  const std::size_t FileNameLen{ OutFileName.find_last_of( '.' ) };
  if( FileNameLen == std::string::npos )
    OutFileName.append( 1, '.' );
  else
    OutFileName.resize( FileNameLen + 1 );
  OutFileName.append( TEXT_EXT );
  out.WriteSummary( OutFileName );
}

void R1R2Maker::Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                    const std::string &qSnk, const std::string &qSrc,
                    const ModelInfo &miSnk, const ModelInfo &miSrc, const std::string &opSnk, const std::string &opSrc )
{
  const Common::FileNameMomentum &MomSrc{ fna.GetMomentum() };
  const Common::FileNameMomentum &MomSnk{ fna.GetMomentum( Common::Momentum::SinkPrefix )  };
  // Make sure I have the model for Z_V already loaded
  const ModelInfo &miZVSnk{ ZVmi.at( QDT( qSnk, fna.DeltaT, opSrc, opSnk ) ) };
  const ModelInfo &miZVSrc{ ZVmi.at( QDT( qSrc, fna.DeltaT, opSnk, opSrc ) ) };
  // Now read the two-point functions
  static const std::string Spaces( 2, ' ' );
  std::vector<Fold> Corr2pt;
  static constexpr int NumC2{ 2 };
  Corr2pt.reserve( NumC2 );
  {
    std::string C2NameSrc{ miSrc.FileName2pt };
    AppendOps( C2NameSrc, opSrc, opSrc );
    C2NameSrc = Common::MakeFilename( C2NameSrc, Common::sFold, model.Seed, DEF_FMT );
    Corr2pt.emplace_back( C2NameSrc, Spaces.c_str() );
    std::string C2NameSnk{ miSnk.FileName2pt };
    AppendOps( C2NameSnk, opSnk, opSnk );
    C2NameSnk = Common::MakeFilename( C2NameSnk, Common::sFold, model.Seed, DEF_FMT );
    Corr2pt.emplace_back( C2NameSnk, Spaces.c_str() );
  }

  // Load the correlators we need
  // The first two files are the numerator
  // Followed by one or two files for the denominator. Both required for R2
  std::string Dir{ fna.Dir };
  Dir.append( fna.BaseShort );
  const std::size_t Len{ Dir.length() };
  Common::DataSet<Scalar> ds;
  static constexpr int iInitial{ 0 };
  static constexpr int iFinal{ 1 };
  const Scalar * pSrc[4];
  int nT[4];
  for( int Snk_Src = 0; Snk_Src < 4; ++Snk_Src )
  {
    int iSnk, iSrc;
    switch( Snk_Src )
    {
      case 0: // This is the way the file specified on command-line
        iSnk = iFinal;
        iSrc = iInitial;
        break;
      case 1: // Reversed compared with command-line
        iSnk = iInitial;
        iSrc = iFinal;
        break;
      case 2:
        iSnk = iInitial;
        iSrc = iInitial;
        break;
      case 3:
        iSnk = iFinal;
        iSrc = iFinal;
        break;
    }
    Dir.resize( Len );
    Dir.append( iSnk == iInitial ? qSrc : qSnk );
    Dir.append( 1, '_' );
    Dir.append( iSrc == iInitial ? qSrc : qSnk );
    Dir.append( fnaSuffix );
    Common::AppendGammaDeltaT( Dir, iSnk == iSrc ? Common::Gamma::Algebra::GammaT : fna.Gamma[0], fna.DeltaT );
    fna.AppendMomentum( Dir, iSrc == iFinal ? MomSnk : MomSrc, MomSrc.Name );
    fna.AppendMomentum( Dir, iSnk == iInitial ? MomSrc : MomSnk, MomSnk.Name );
    Dir.append( Common::Underscore );
    Dir.append( iSnk == iInitial ? opSrc : opSnk );
    Dir.append( Common::Underscore );
    Dir.append( iSrc == iInitial ? opSrc : opSnk );
    Dir = Common::MakeFilename( Dir, fna.Type, fna.Seed, fna.Ext );
    if( Common::FileExists( Dir ) )
    {
      const int cNum{ ds.LoadCorrelator( Common::FileNameAtt( Dir ), Common::COMPAT_DISABLE_BASE ) };
      Fold &f { ds.corr[cNum] };
      pSrc[cNum] = f[Fold::idxCentral];
      nT[cNum] = f.Nt();
    }
    else if( iSnk == iSrc )
    {
      ds.corr.resize( 2 );
      break; // No point loading the other denominator file
    }
    else
      throw std::runtime_error( Dir + " doesn't exist" );
  }
  // Ensure the correlator includes timeslice DeltaT
  const int NTHalf{ ds.corr[0].Nt() / 2 };
  if( fna.DeltaT > NTHalf )
    throw std::runtime_error( "DeltaT " + std::to_string( fna.DeltaT ) + " > Nt/2 " + std::to_string( NTHalf ) );
  // Make somewhere to put R2
  const bool bMakeR2{ ds.corr.size() == 4 };
  const int NumRatio{ bMakeR2 ? 2 : 1 };
  std::vector<Fold> out;
  std::vector<Scalar *> pDst( NumRatio );
  for( int i = 0; i < NumRatio; i++ )
  {
    out.emplace_back( ds.NSamples, fna.DeltaT + 1 );
    for( int j = 0; j < ( i == 0 ? 2 : ds.corr.size() ); ++j )
      out[i].FileList.push_back( ds.corr[j].Name_.Filename );
    if( i == 0 )
    {
      out[i].FileList.push_back( Corr2pt[0].Name_.Filename );
      out[i].FileList.push_back( Corr2pt[1].Name_.Filename );
    }
    out[i].CopyAttributes( ds.corr[0] );
    out[i].NtUnfolded = ds.corr[0].Nt();
    pDst[i] = out[i][Fold::idxCentral];
  }
  const Scalar * pE[NumC2];
  pE[0] = miSrc.m[Fold::idxCentral] + miSrc.idxE0;
  pE[1] = miSnk.m[Fold::idxCentral] + miSnk.idxE0;
  const Scalar * pC2[NumC2];
  pC2[0] = Corr2pt[0][Fold::idxCentral];
  pC2[1] = Corr2pt[1][Fold::idxCentral];
  const Scalar * pZVSrc{ miZVSrc.m[Fold::idxCentral] + miZVSrc.idxE0 };
  const Scalar * pZVSnk{ miZVSnk.m[Fold::idxCentral] + miZVSnk.idxE0 };
  for( int idx = Fold::idxCentral; idx < ds.NSamples; ++idx )
  {
    // Get the energies and compute the Correlator at Delta T with backward propagating wave subtracted
    // Allow energies to be frozen to 1 (which means don't subtract backward propagating wave)
    double  E[NumC2];
    double C2[NumC2];
    for( int i = 0; i < NumC2; ++i )
    {
      C2[i] = pC2[i][fna.DeltaT];
      if( fEnergy == Freeze::One )
      {
        // If we're freezing energy to one, it makes no sense to try to subtract half the midpoint
        E[i] = 1;
        if( NTHalf == fna.DeltaT )
          C2[i] *= 0.5; // In the middle of the lattice, this is needed for consistency
      }
      else
      {
         E[i] = *pE[i];
        C2[i] -= 0.5 * pC2[i][NTHalf] * std::exp( - E[i] * ( NTHalf - fna.DeltaT ) );
      }
    }
    const double EProd = E[0] * E[1];
    const double C2Prod = std::abs( C2[0] * C2[1] );
    // Now compute the ratios
    for( int t = 0; t <= fna.DeltaT; ++t )
    {
      const double n = std::abs( EProd * pSrc[0][t] * pSrc[1][t] ); // fna.DeltaT - t
      if( !std::isfinite( n ) )
        throw std::runtime_error( "Numerator Overflow" );
      const Scalar ZVSrc{ fZV == Freeze::One ? 1 : *pZVSrc };
      const Scalar ZVSnk{ fZV == Freeze::One ? 1 : *pZVSnk };
      const Scalar R1{ 2 * std::sqrt( ( n / C2Prod ) * ZVSrc * ZVSnk ) };
      pDst[0][t] = R1;
      if( bMakeR2 )
      {
        const Scalar R2{ 2 * std::sqrt( n / std::abs( pSrc[2][t] * pSrc[3][t] ) ) };
        pDst[1][t] = R2;
      }
    }
    // Symmetrise the output
    if( bSymmetrise )
    {
      assert( ! ( fna.DeltaT & 1 ) && "DeltaT must be even" );
      for( int i = 0; i < NumRatio; ++i )
      {
        for( int t = 1; t < fna.DeltaT / 2; ++t )
        {
          pDst[i][t] = ( pDst[i][t] + pDst[i][fna.DeltaT - t] ) / 2;
          pDst[i][fna.DeltaT - t] = pDst[i][t];
        }
      }
    }
    for( int i = 0; i < NumRatio; ++i )
      pDst[i] += out[i].Nt();
    for( int i = 0; i < ds.corr.size(); ++i )
      pSrc[i] += nT[i];
    if( fEnergy == Freeze::None )
    {
      pE[0] += miSrc.m.Nt();
      pE[1] += miSnk.m.Nt();
    }
    if( fZV == Freeze::None )
    {
      pZVSrc+= miZVSrc.m.Nt();
      pZVSnk+= miZVSnk.m.Nt();
    }
    pC2[0] += Corr2pt[0].Nt();
    pC2[1] += Corr2pt[1].Nt();
  }
  for( int i = 0; i < NumRatio ; i++ )
  {
    out[i].MakeCorrSummary( nullptr );
    std::string OutFileName{ outBase };
    OutFileName.append( 1, 'R' );
    OutFileName.append( 1, '1' + i );
    OutFileName.append( 1, '_' );
    OutFileName.append( qSnk );
    OutFileName.append( 1, '_' );
    OutFileName.append( qSrc );
    OutFileName.append( fnaSuffix );
    Common::AppendGammaDeltaT( OutFileName, fna.Gamma[0], fna.DeltaT );
    fna.AppendMomentum( OutFileName, MomSrc ? MomSrc : MomSnk, MomSrc.Name );
    AppendOps( OutFileName, opSnk, opSrc );
    OutFileName = Common::MakeFilename( OutFileName, fna.Type, fna.Seed, fna.Ext );
    std::cout << "->" << OutFileName << Common::NewLine;
    out[i].Write( OutFileName, Common::sFold.c_str() );
    const std::size_t FileNameLen{ OutFileName.find_last_of( '.' ) };
    if( FileNameLen == std::string::npos )
      OutFileName.append( 1, '.' );
    else
      OutFileName.resize( FileNameLen + 1 );
    OutFileName.append( TEXT_EXT );
    out[i].WriteSummary( OutFileName );
  }
}

/*void DebugInput( const std::string &Input, std::array<std::string,3> &a )
{
  std::cout << "Input \"" << Input << "\"" << std::endl;
  {
    std::istringstream ss( Input );
    ss >> a;
  }
  for( std::size_t i = 0; i < a.size(); ++i )
    std::cout << "\ta[" << i << "]=\"" << a[i] << "\"" << std::endl;
}

bool Debug()
{
  std::array<std::string,3> a;
  DebugInput( "alpha bravo delta gamma", a );
  DebugInput( "vega upsilon", a );
  return false;
}*/

Freeze GetFreeze( const Common::CommandLine &cl, const std::string &sOne, const std::string &sCentral )
{
  Freeze f;
  const bool bOne{ cl.GotSwitch( sOne ) };
  const bool bCentral{ cl.GotSwitch( sCentral ) };
  if( bOne )
  {
    if( bCentral )
      throw std::invalid_argument( "Cannot specify both " + sOne + " and " + sCentral );
    f = Freeze::One;
  }
  else if( bCentral )
    f = Freeze::Central;
  else
    f = Freeze::None;
  return f;
}

int main(int argc, const char *argv[])
{
  //if(!Debug()) return EXIT_FAILURE;
  static const char DefaultType[]{ "ZV" };
  static const char DefaultERE[]{ R"((quark_|anti_)([[:alpha:]]+[[:digit:]]*)_([[:alpha:]]+[[:digit:]]*))" };
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i2", CL::SwitchType::Single, "" },
      {"i3", CL::SwitchType::Single, "" },
      {"im", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"type", CL::SwitchType::Single, DefaultType },
      {"ssre", CL::SwitchType::Single, DefaultERE },
      {"swap", CL::SwitchType::Flag, nullptr},
      {"eone", CL::SwitchType::Flag, nullptr},
      {"ec", CL::SwitchType::Flag, nullptr},
      {"zvone", CL::SwitchType::Flag, nullptr},
      {"zvc", CL::SwitchType::Flag, nullptr},
      {"nosym", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() > 1 )
    {
      // Read the list of fits I've chosen to use
      std::string TypeParams{ cl.SwitchValue<std::string>("type") };
      const std::string Type{ Common::ExtractToSeparator( TypeParams ) };
      const std::string InBase{ cl.SwitchValue<std::string>("i3") };
      std::unique_ptr<Maker> m( Maker::Make( Type, TypeParams, InBase, cl.SwitchValue<std::string>("i2"),
                                             cl.SwitchValue<std::string>("im"), cl.SwitchValue<std::string>("o"),
                                             std::regex( cl.SwitchValue<std::string>("ssre"),
                                                         std::regex::extended | std::regex::icase ),
                                             cl.GotSwitch("swap"),
                                             GetFreeze( cl, "eone", "ec" ),
                                             GetFreeze( cl, "zvone", "zvc" ),
                                             !cl.GotSwitch("nosym"),
                                             cl.Args[0] ) );
      bShowUsage = false;
      std::vector<std::string> FileList{ Common::glob( ++cl.Args.begin(), cl.Args.end(), InBase.c_str() ) };
      std::size_t Count{ 0 };
      std::size_t Done{ 0 };
      for( std::string &sFile : FileList )
      {
        try
        {
          ++Count;
          m->Make( sFile );
          ++Done;
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
      std::cout << Done << " of " << Count << Common::Space << Type << " ratios created" << Common::NewLine;
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
    "--i2   Input2 filename prefix\n"
    "--i3   Input3 filename prefix\n"
    "--im   Input model filename prefix\n"
    "-o     Output filename prefix\n"
    "--type Ratio type[,parameters[,...]] (default " << DefaultType << ")\n"
    "--ssre Extended regex for sink/source type, default\n       " << DefaultERE << "\n"
    "       http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap09.html\n"
    //"-n     Number of samples to fit, 0 (default) = all available from bootstrap\n"
    "Flags:\n"
    "--swap Swap source / sink order in regex\n"
    "--eone Freeze the energy to 1 (in R1R2 ratios)\n"
    "--ec   Freeze the energy fit to it's central value\n"
    "--zvone Freeze ZV to 1 (in R1R2 ratios)\n"
    "--zvc  Freeze ZV fit to it's central value\n"
    "--nosym Don't symmetrise the waves\n"
    "--help This message\n";
  }
  return iReturn;
}
