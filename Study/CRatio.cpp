/*************************************************************************************
 
 Create Ratios, e.g. R1, R2 and associated values such as Z_V
 Source file: Ratio.cpp
 Copyright (C) 2020-2022
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

#include "CRatio.hpp"

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

static const char EFitSwitch[] = "efit";
const std::string DefaultColumnName{ "E0" };
static const char LoadFilePrefix[] = "  ";

std::istream & operator>>( std::istream &is, QDT &qdt )
{
  if( is >> qdt.q >> qdt.deltaT >> qdt.opSnk >> qdt.opSrc )
    return is;
  throw std::runtime_error( "Error reading QDT" );
}

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

template<typename FileT>
void FileCache<FileT>::clear()
{
  Files.clear();
  NameMap.clear();
  opNames.clear();
}

template<typename FileT>
int FileCache<FileT>::GetIndex( const std::string &Filename )
{
  // If it's already in our list, return it's index
  typename TNameMap::iterator it = NameMap.find( Filename );
  if( it != NameMap.end() )
    return it->second;
  // Needs to be added - check we're not going to overflow
  if( NameMap.size() >= std::numeric_limits<int>::max() )
    throw std::runtime_error( "Too many entries in FileCache - switch to std::size_t" );
  int iIndex = static_cast<int>( NameMap.size() );
  const std::size_t OldOpNameSize{ opNames.size() };
  it = NameMap.emplace( Filename, iIndex ).first;
  try
  {
    Files.resize( iIndex + 1 );
    std::string Fullname{ Base };
    Fullname.append( Filename );
    Files[iIndex].SetName( Fullname, &opNames );
  }
  catch(...)
  {
    // Maintain the invariant - NameMap and Files always match
    opNames.resize( OldOpNameSize );
    Files.resize( iIndex );
    NameMap.erase( it );
    throw;
  }
  return iIndex;
}

template<typename FileT>
FileT & FileCache<FileT>::operator[]( int iIndex )
{
  if( iIndex < 0 || iIndex >= Files.size() )
    throw std::runtime_error( "FileCache index " + std::to_string( iIndex ) + " out of bounds" );
  FileT &f{ Files[iIndex] };
  if( !f.Nt() ) f.Read( PrintPrefix ); // Read file if not in cache
  return f;
}

QP QuarkReader::Convert( const std::string &Filename ) const
{
  Common::FileNameAtt fna( Filename );
  if( fna.p.empty() )
    return QP( *this, Common::p0 );
  return QP( *this, fna.p.begin()->second );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::KeyFileCache( const std::string &modelBase )
: model{modelBase}
{
  clear();
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
void KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::clear()
{
  model.clear();
  KeyMap.clear();
  freeze = Freeze::One;
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
template<typename K, typename R>
typename std::enable_if<std::is_same<K, R>::value, std::map<Key, std::string, LessKey>>::type
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::ReadNameMap( std::ifstream &s )
{
  return Common::KeyValReader<Key, std::string, LessKey>::Read( s );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
template<typename K, typename R>
typename std::enable_if<!std::is_same<K, R>::value, std::map<Key, std::string, LessKey>>::type
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::ReadNameMap( std::ifstream &s )
{
  // Read the file list
  using KeyNameRead = std::multimap<KeyRead, std::string, LessKeyRead>;
  using StringMapCIReader = Common::KeyValReader<KeyRead, std::string, LessKeyRead, KeyNameRead>;
  KeyNameRead ReadList{ StringMapCIReader::Read( s ) };
  // Now make a unique key for each row from name and filename
  KeyReadMapT m;
  for( const auto &k : ReadList )
    m.emplace( k.first.Convert( k.second ), k.second );
  return m;
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
void KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::Read( const std::string &filename, const char *Prefix )
{
  clear();
  model.PrintPrefix = Prefix;
  std::string sParams{ filename };
  std::string Filename{ Common::ExtractToSeparator( sParams ) };
  if( Filename.length() == 1 && Filename[0] == '1' )
  {
    freeze = Freeze::One;
    FrozenOptions( sParams );
    if( !sParams.empty() )
      throw std::runtime_error( "KeyFileCache bad options: " + sParams );
    return;
  }
  else if( sParams.empty() )
    freeze = Freeze::None;
  else if( sParams.length() == 1 && std::toupper( sParams[0] ) == 'C' )
    freeze = Freeze::Central;
  else
    throw std::runtime_error( "Freeze type " + sParams + " unrecognised" );

  // Load the file
  std::string FitListName{ model.Base };
  FitListName.append( Filename );
  std::cout << "Reading 2pt fit list from " << FitListName << Common::NewLine;
  {
    KeyReadMapT m;
    {
      std::ifstream s( FitListName );
      if( !Common::FileExists( FitListName ) || s.bad() )
        throw std::runtime_error( "2pt fit file \"" + FitListName + "\" not found" );
      m = ReadNameMap( s );
    }
    // Put all the files in our cache and save their index
    for( const auto &i : m )
      KeyMap.emplace( i.first, model.GetIndex( i.second ) );
  }
  /*const ModelInfo *idxFirst{ nullptr };
  for( typename FileMap::iterator it = fileMap.begin(); it != fileMap.end(); ++it )
  {
    ModelInfo &mi{ it->second };
    std::string ModelFileName{ ModelBase };
    ModelFileName.append( it->first );
    Common::FileNameAtt fna( ModelFileName );
    std::string MsgPrefix( 2, ' ' );
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
  }*/
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
Vector KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::GetVector( const Key &key, const std::string &Name )
{
  Vector v;
  if( freeze == Freeze::One )
  {
    static Scalar One{ 1 };
    v.MapView( &One, std::numeric_limits<int>::max(), 0 );
  }
  else
  {
    Model &m{ model[KeyMap[key]] };
    Scalar * p{ m[Model::idxCentral] + m.GetColumnIndex( Name ) };
    if( freeze == Freeze::Central )
      v.MapView( p, std::numeric_limits<int>::max(), 0 );
    else
      v.MapView( p, m.NumSamples() - Model::idxCentral, m.Nt() );
  }
  return v;
}

void QPModelMap::clear()
{
  Spectator.clear();
  Base::clear();
}

void QPModelMap::FrozenOptions( std::string &Options )
{
  if( Options.empty() )
    throw std::runtime_error( "QPModelMap spectator must be passed in options when frozen" );
  Spectator = std::move( Options );
}

std::string QPModelMap::Get2ptName( const QP &key )
{
  std::string s;//{ C2Base };
  if( freeze == Freeze::One )
  {
    s.append( key.q );
    Append( s, Spectator );
    Common::FileNameMomentum p( key.p );
    p.Name = Common::Momentum::DefaultPrefix;
    s.append( p.FileString() );
  }
  else
  {
    const Model &m{ model[KeyMap[key]] };
    s = m.Name_.Base;
    if( m.Name_.Extra.size() != 1 )
      throw std::runtime_error( "Expected " + s + " to have model info removed by parser" );
  }
  return s;
}

//template class ModelMap<QDT, LessQDT, QDT, LessQDT>;
//template class ModelMap<QP, LessQP, QuarkReader, Common::LessCaseInsensitive>;

/*QPMapModelInfo::QPMapModelInfo( const std::string &Filename,
                                const std::string &modelBase, const std::string &C2Base )
{
  std::string FitListName{ modelBase };
  FitListName.append( Filename );
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
}*/

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
  const Common::FileNameMomentum &MomSnk{ fna.GetMomentum( Common::Momentum::SinkPrefix ) };
  std::smatch base_match;
  if( !std::regex_search( fna.BaseShort, base_match, RegExExt ) || base_match.size() != 4 )
    throw std::runtime_error( "Can't extract sink/source from " + fna.BaseShort );
  const std::string qSnk = base_match[RegExSwap ? 3 : 2];
  const std::string qSrc = base_match[RegExSwap ? 2 : 3];
  const std::string FileNameSuffix{ base_match.suffix() };
  std::string FileNamePrefix{ base_match.prefix() };
  FileNamePrefix.append( base_match[1] );
  fna.BaseShort = FileNamePrefix;
  QP QPSnk{ QP( qSnk, MomSnk ) };
  QP QPSrc{ QP( qSrc, MomSrc ) };
  Vector ESnk{ model.GetVector( QPSnk ) };
  Vector ESrc{ model.GetVector( QPSrc ) };
  Make( fna, FileNameSuffix, QPSnk, QPSrc, ESnk, ESrc, OpNames[fna.op[1]], OpNames[fna.op[0]] );
}

void ZVMaker::Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                    const QP &QPSnk, const QP &QPSrc,
                    const Vector &ESnk, const Vector &ESrc, const std::string &opSnk, const std::string &opSrc )
{
  const std::string &qSnk{ QPSnk.q };
  const std::string &qSrc{ QPSrc.q };
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
  Fold C3;
  C3.Read( fna.Filename, LoadFilePrefix );

  // Ensure the correlator includes timeslice DeltaT
  const int NTHalf{ C3.Nt() / 2 };
  if( fna.DeltaT > NTHalf )
    throw std::runtime_error( "DeltaT " + std::to_string( fna.DeltaT ) + " > Nt/2 " + std::to_string( NTHalf ) );

  // Now read the two-point function
  std::string C2Name{ model.Get2ptName( QPSnk ) };
  AppendOps( C2Name, opSnk, opSrc );
  C2Name = Common::MakeFilename( C2Name, Common::sFold, fna.Seed, DEF_FMT );
  Fold &C2{ Cache2[C2Name] };
  //C2.Read( C2Name, LoadFilePrefix );
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
  out.FileList.push_back( model[QPSrc].Name_.Filename );
  out.CopyAttributes( C3 );
  out.NtUnfolded = C3.Nt();
  Scalar * pDst{ out[Fold::idxCentral] };

  const Scalar * pC2 { C2[Fold::idxCentral] };
  const Scalar * pSrc{ C3[Fold::idxCentral] };
  for( int idx = Fold::idxCentral; idx < NumSamples; ++idx )
  {
    const double CTilde{ pC2[fna.DeltaT] - 0.5 * pC2[NTHalf] * std::exp( - ESrc[idx - Fold::idxCentral] * ( NTHalf - fna.DeltaT ) ) };
    for( int t = 0; t <= fna.DeltaT; ++t )
    {
      const double z{ CTilde / *pSrc++ };
      if( !std::isfinite( z ) )
        throw std::runtime_error( "ZV Overflow" );
      *pDst++ = z;
    }
    pSrc += C3.Nt() - ( fna.DeltaT + 1 );
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
  OutFileName = Common::MakeFilename( OutFileName, Common::sFold, fna.Seed, DEF_FMT );
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

R1R2Maker::R1R2Maker( std::string &TypeParams, const Common::CommandLine &cl )
: Maker( TypeParams, cl ), ZVmi( modelBase )
{
  // TODO fix ZV freeze
  //GetFreeze( cl, "zvone", "zvc" );
  ZVmi.Read( TypeParams, LoadFilePrefix );
  TypeParams.clear(); // I used this to load my map from
}

void R1R2Maker::Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                    const QP &QPSnk, const QP &QPSrc,
                    const Vector &ESnk, const Vector &ESrc, const std::string &opSnk, const std::string &opSrc )
{
  const Common::FileNameMomentum &MomSrc{ QPSrc.p };
  const Common::FileNameMomentum &MomSnk{ QPSnk.p };
  // Make sure I have the model for Z_V already loaded
  const std::string &qSnk{ QPSnk.q };
  const std::string &qSrc{ QPSrc.q };
  QDT QDTSnk( qSnk, fna.DeltaT, opSrc, opSnk );
  QDT QDTSrc( qSrc, fna.DeltaT, opSnk, opSrc );
  Vector ZVSnk{ ZVmi.GetVector( QDTSnk ) };
  Vector ZVSrc{ ZVmi.GetVector( QDTSrc ) };
  // Now read the two-point functions
  static constexpr int NumC2{ 2 };
  std::array<Fold *, 2> Corr2pt;
  {
    std::string C2NameSrc{ model.Get2ptName( QPSrc ) };
    AppendOps( C2NameSrc, opSrc, opSrc );
    C2NameSrc = Common::MakeFilename( C2NameSrc, Common::sFold, fna.Seed, DEF_FMT );
    int hFile0{ Cache2.GetIndex( C2NameSrc ) };
    std::string C2NameSnk{ model.Get2ptName( QPSnk ) };
    AppendOps( C2NameSnk, opSnk, opSnk );
    C2NameSnk = Common::MakeFilename( C2NameSnk, Common::sFold, fna.Seed, DEF_FMT );
    Corr2pt[1] = &Cache2[C2NameSnk];
    Corr2pt[0] = &Cache2[hFile0];
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
      out[i].FileList.push_back( Corr2pt[0]->Name_.Filename );
      out[i].FileList.push_back( Corr2pt[1]->Name_.Filename );
    }
    out[i].CopyAttributes( ds.corr[0] );
    out[i].NtUnfolded = ds.corr[0].Nt();
    pDst[i] = out[i][Fold::idxCentral];
  }
  const std::array<const Vector *, NumC2> pE{ &ESrc, &ESnk };
  const Scalar * pC2[NumC2];
  pC2[0] = (*Corr2pt[0])[Fold::idxCentral];
  pC2[1] = (*Corr2pt[1])[Fold::idxCentral];
  for( int idx = Fold::idxCentral; idx < ds.NSamples; ++idx )
  {
    // Get the energies and compute the Correlator at Delta T with backward propagating wave subtracted
    // Allow energies to be frozen to 1 (which means don't subtract backward propagating wave)
    double  E[NumC2];
    double C2[NumC2];
    for( int i = 0; i < NumC2; ++i )
    {
      C2[i] = pC2[i][fna.DeltaT];
      if( model.freeze == Freeze::One )
      {
        // If we're freezing energy to one, it makes no sense to try to subtract half the midpoint
        E[i] = 1;
        if( NTHalf == fna.DeltaT )
          C2[i] *= 0.5; // In the middle of the lattice, this is needed for consistency
      }
      else
      {
         E[i] = (*pE[i])[idx - Fold::idxCentral];
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
      const Scalar R1{ 2 * std::sqrt( ( n / C2Prod )
                                     * ZVSrc[idx - Fold::idxCentral] * ZVSnk[idx - Fold::idxCentral] ) };
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
    pC2[0] += Corr2pt[0]->Nt();
    pC2[1] += Corr2pt[1]->Nt();
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

Maker::Maker( std::string &TypeParams, const Common::CommandLine &cl )
: modelBase{cl.SwitchValue<std::string>("im")},
  outBase{cl.SwitchValue<std::string>("o")},
  bSymmetrise{!cl.GotSwitch("nosym")},
  RegExSwap{cl.GotSwitch("swap")},
  RegExExt{std::regex( cl.SwitchValue<std::string>("ssre"), std::regex::extended | std::regex::icase )},
  model{modelBase},
  Cache2{cl.SwitchValue<std::string>("i2"), LoadFilePrefix },
  Cache3{cl.SwitchValue<std::string>("i3"), LoadFilePrefix }
{
  Common::MakeAncestorDirs( outBase );
  model.Read( cl.SwitchValue<std::string>( EFitSwitch ), LoadFilePrefix );
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
      {EFitSwitch, CL::SwitchType::Single, ""},
      {"swap", CL::SwitchType::Flag, nullptr},
      {"zv", CL::SwitchType::Flag, nullptr},
      {"nosym", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() > 1 )
    {
      // Read the list of fits I've chosen to use
      std::string TypeParams{ cl.SwitchValue<std::string>("type") };
      const std::string Type{ Common::ExtractToSeparator( TypeParams ) };
      std::unique_ptr<Maker> m;
      if( Common::EqualIgnoreCase( Type, "R1R2" ) )
        m.reset( new R1R2Maker( TypeParams, cl ) );
      else if( Common::EqualIgnoreCase( Type, "ZV" ) )
        m.reset( new ZVMaker( TypeParams, cl ) );
      else
        throw std::runtime_error( "I don't know how to make type " + Type );
      if( !TypeParams.empty() )
        throw std::runtime_error( "Unused parameters: " + TypeParams );
      bShowUsage = false;
      std::size_t Prefix3ptLen{ m->Cache3.Base.size() };
      std::vector<std::string> FileList{ Common::glob( cl.Args.begin(), cl.Args.end(), m->Cache3.Base.c_str() ) };
      std::size_t Count{ 0 };
      std::size_t Done{ 0 };
      for( std::string &sFile : FileList )
      {
        try
        {
          ++Count;
          //sFile.erase( 0, Prefix3ptLen );
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
