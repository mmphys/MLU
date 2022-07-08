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

Model DummyModel{};
const std::string DefaultColumnName{ "E0" };
const char LoadFilePrefix[] = "  ";
static Scalar ConstantOne{ 1 };

std::istream & operator>>( std::istream &is, QDT &qdt )
{
  if( is >> qdt.q >> qdt.deltaT >> qdt.opSnk >> qdt.opSrc )
    return is;
  throw std::runtime_error( "Error reading QDT" );
}

std::ostream & operator<<( std::ostream &os, const QDT &qdt )
{
  return os << "q=" << qdt.q << ", deltaT=" << qdt.deltaT
            << ", sink=" << qdt.opSnk << ", source=" << qdt.opSrc;
}

std::ostream & operator<<( std::ostream &os, const QP &qp )
{
  return os << "q=" << qp.q << ", " << qp.p.Name << "=" << qp.p;
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
int FileCache<FileT>::GetIndex( const std::string &Filename, const char * printPrefix )
{
  if( !printPrefix )
    printPrefix = PrintPrefix;
  // If it's already in our list, return it's index
  typename TNameMap::iterator it = NameMap.find( Filename );
  if( it != NameMap.end() )
    return it->second;
  // Needs to be added
  if( !Common::FileExists( Base + Filename ) )
  {
    // File doesn't exist
    if( printPrefix )
      std::cout << printPrefix << "- " << Filename << std::endl;
    return BadIndex;
  }
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
FileT & FileCache<FileT>::operator()( int iIndex, const char * printPrefix )
{
  if( iIndex < 0 || iIndex >= Files.size() )
    throw std::runtime_error( "FileCache index " + std::to_string( iIndex ) + " out of bounds" );
  FileT &f{ Files[iIndex] };
  if( !f.Nt() )
    f.Read( printPrefix ); // Read file if not in cache
  else if( printPrefix )
    std::cout << printPrefix << "C " << f.Name_.Filename << std::endl;
  return f;
}

QP QuarkReader::operator()( const std::string &Filename ) const
{
  Common::FileNameAtt fna( Filename );
  if( fna.p.empty() )
    return QP( *this, Common::p0 );
  return QP( *this, fna.p.begin()->second );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
void KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::clear()
{
  model.clear();
  KeyMap.clear();
  freeze = Freeze::Constant;
  FieldMap.clear();
  ConstantMap.clear();
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
template<typename K, typename R>
typename std::enable_if<std::is_same<K, R>::value, std::map<Key, std::string, LessKey>>::type
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::ReadNameMap( const std::string &Filename, const char *p )
{
  return Common::KeyValReader<Key, std::string, LessKey>::Read( Filename, p );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
template<typename K, typename R>
typename std::enable_if<!std::is_same<K, R>::value, std::map<Key, std::string, LessKey>>::type
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::ReadNameMap( const std::string &Filename, const char *p )
{
  // Read the file list
  using KeyNameRead = std::multimap<KeyRead, std::string, LessKeyRead>;
  using StringMapCIReader = Common::KeyValReader<KeyRead, std::string, LessKeyRead, KeyNameRead>;
  KeyNameRead ReadList{ StringMapCIReader::Read( Filename, p ) };
  // Now make a unique key for each row from name and filename
  KeyReadMapT m;
  for( const auto &k : ReadList )
    m.emplace( k.first( k.second ), k.second );
  return m;
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
void KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::Read(const std::string &filename, const char *Prefix)
{
  clear();
  std::string sParams{ filename };
  std::string Filename{ Common::ExtractToSeparator( sParams ) };
  if( Filename.empty() )
  {
    freeze = Freeze::Constant;
    FrozenOptions( sParams );
    ConstantMap = Common::KeyValReader<std::string, Scalar>::Read( Common::ArrayFromString( sParams ) );
    return;
  }

  freeze = Freeze::None;
  {
    // See whether next parameter is 'C'
    std::string RemainingParams{ sParams };
    std::string Options{ Common::ExtractToSeparator( RemainingParams ) };
    if( Options.length() == 1 && std::toupper( sParams[0] ) == 'C' )
    {
      freeze = Freeze::Central;
      sParams = RemainingParams;
    }
  }
  FieldMap = Common::KeyValReader<std::string, std::string>::Read( Common::ArrayFromString( sParams ) );

  // Load the map from the key to the filename
  std::string FitListName{ model.Base };
  FitListName.append( Filename );
  static const char szDescription[]{ "2pt fit list" };
  std::cout << "Reading " << szDescription << " from " << FitListName << Common::NewLine;
  KeyReadMapT m{ ReadNameMap( FitListName, szDescription ) };
  // Put all the filenames in our cache and save their index
  // This delays loading until the files are referenced
  for( const auto &i : m )
  {
    const int Handle{ model.GetIndex( i.second, Prefix ) };
    if( Handle != FileCacheT::BadIndex )
      KeyMap.emplace( i.first, Handle );
  }
}

// Throws informative errors if frozen to constant, or index out of bounds
template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
int KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::GetIndex( const Key &key ) const
{
  if( freeze != Freeze::Constant )
  {
    // If it's already in our list, return it's index
    typename KeyMapT::const_iterator it = KeyMap.find( key );
    if( it != KeyMap.end() )
      return it->second;
    std::ostringstream os;
    os << "Please add entry to cache for: " << key;
    throw std::runtime_error( os.str().c_str() );
  }
  std::ostringstream os;
  os << "KeyFileCache: GetIndex(" << key << ") invalid when frozen";
  throw std::runtime_error( os.str().c_str() );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
M & KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::operator()( const Key &key, const char * PrintPrefix )
{
  if( freeze==Freeze::Constant )
    return DummyModel;
  return model( GetIndex( key ), PrintPrefix );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename M>
Vector KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, M>::GetVector(const Key &key, const std::string &Name)
{
  Vector v;
  if( freeze == Freeze::Constant )
  {
    // Freeze to value specified on command line, or 1 if not specified
    typename ConstantMapT::iterator it = ConstantMap.find( Name );
    v.MapView( it == ConstantMap.end() ? &ConstantOne : &it->second, std::numeric_limits<int>::max(), 0 );
  }
  else
  {
    M &m{ operator[]( key ) };
    typename FieldMapT::iterator it = FieldMap.find( Name );
    Scalar * p{ m[M::idxCentral] + m.GetColumnIndex( it == FieldMap.end() ? Name : it->second ) };
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
  std::string Spec{ Common::ExtractToSeparator( Options ) };
  // Spectator is immediately followed by Key=Value pairs, so check spectator is actually present
  if( Spec.empty() || Spec.find_first_of( Common::EqualSign ) != std::string::npos )
    throw std::runtime_error( "Spectator must be passed in options when energies frozen" );
  Spectator = std::move( Spec );
}

std::string QPModelMap::Get2ptName( const QP &key )
{
  std::string s;//{ C2Base };
  if( freeze == Freeze::Constant )
  {
    s.append( key.q );
    Append( s, Spectator );
    Common::FileNameMomentum p( key.p );
    p.Name = Common::Momentum::DefaultPrefix;
    s.append( p.FileString() );
  }
  else
  {
    const Model &m{ model[GetIndex( key )] };
    s = m.Name_.Base;
    if( m.Name_.Extra.size() != 1 )
      throw std::runtime_error( "Expected " + s + " to have model info removed by parser" );
  }
  return s;
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
  std::smatch base_match;
  if( !std::regex_search( fna.BaseShort, base_match, RegExExt ) || base_match.size() != 4 )
    throw std::runtime_error( "Can't extract sink/source from " + fna.BaseShort );
  const std::string qSnk = base_match[RegExSwap ? 3 : 2];
  const std::string qSrc = base_match[RegExSwap ? 2 : 3];
  const std::string FileNameSuffix{ base_match.suffix() };
  std::string FileNamePrefix{ base_match.prefix() };
  FileNamePrefix.append( base_match[1] );
  fna.BaseShort = FileNamePrefix;
  SSInfo Snk( EFit, OpNames[fna.op[1]], QP( qSnk, fna.GetMomentum( Common::Momentum::SinkPrefix ) ) );
  SSInfo Src( EFit, OpNames[fna.op[0]], QP( qSrc, fna.GetMomentum() ) );
  Make( fna, FileNameSuffix, Snk, Src );
}

ZVMaker::ZVMaker( const std::string &TypeParams, const Common::CommandLine &cl ) : Maker( cl )
{
  if( !TypeParams.empty() )
    throw std::runtime_error( "ZVMaker does not recognise any parameters: " + TypeParams );
}

void ZVMaker::Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                    const SSInfo &Snk, const SSInfo &Src )
{
  const std::string &qSnk{ Snk.qp.q };
  const std::string &qSrc{ Src.qp.q };
  const Common::FileNameMomentum &MomSrc{ Snk.qp.p };
  const Common::FileNameMomentum &MomSnk{ Src.qp.p };
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

  // Now read the two-point function corresponding to the model (Src and Snk are the same)
  std::string C2Name{ EFit.Get2ptName( Snk.qp ) };
  // DON'T DO THIS - IT JUST MAKES THE ERRORS LARGER
  // Wall sink - Point source ==> swap source and sink on the 2-point for accuracy
  //if( std::toupper( opSnk.back() ) == 'W' && std::toupper( opSrc.back() ) == 'P' )
    //AppendOps( C2Name, opSrc, opSnk );
  //else
    AppendOps( C2Name, Snk.op, Src.op );
  C2Name = Common::MakeFilename( C2Name, Common::sFold, fna.Seed, DEF_FMT );
  Fold &C2{ Cache2[C2Name] };
  //C2.Read( C2Name, LoadFilePrefix );
  int NumSamples {};
  C3.IsCompatible( C2, &NumSamples, Common::COMPAT_DISABLE_BASE | Common::COMPAT_DISABLE_NT );
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
  if( EFit.freeze != Freeze::Constant )
    out.FileList.push_back( Src.EModel.Name_.Filename );
  out.CopyAttributes( C3 );
  out.NtUnfolded = C3.Nt();
  Scalar * pDst{ out[Fold::idxCentral] };

  const Scalar * pC2 { C2[Fold::idxCentral] };
  const Scalar * pSrc{ C3[Fold::idxCentral] };
  for( int idx = Fold::idxCentral; idx < NumSamples; ++idx )
  {
    const double CTilde{ pC2[fna.DeltaT] - 0.5 * pC2[NTHalf] * std::exp( - Src.E[idx - Fold::idxCentral] * ( NTHalf - fna.DeltaT ) ) };
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
  AppendOps( OutFileName, Snk.op, Src.op );
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

RMaker::RMaker( const std::string &TypeParams, const Common::CommandLine &cl )
: Maker( cl ), bAltR3{cl.GotSwitch("r3a")}, ZVmi( modelBase )
{
  ZVmi.Read( TypeParams, LoadFilePrefix );
}

void RMaker::Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                   const SSInfo &Snk, const SSInfo &Src )
{
  // Make sure I have the model for Z_V already loaded
  const std::string &qSnk{ Snk.qp.q };
  const std::string &qSrc{ Src.qp.q };
  QDT QDTSnk( qSnk, fna.DeltaT, Src.op, Snk.op );
  QDT QDTSrc( qSrc, fna.DeltaT, Snk.op, Src.op );
  Model &ZVModelSnk{ZVmi( QDTSnk, LoadFilePrefix )};
  Model &ZVModelSrc{ZVmi( QDTSrc, LoadFilePrefix )};
  Vector ZVSnk{ ZVmi.GetVector( QDTSnk ) };
  Vector ZVSrc{ ZVmi.GetVector( QDTSrc ) };

  // Get the names of all 2pt correlators we will need
  std::array<CorrT, 2> Corr2;
  static constexpr int NumC2{ static_cast<int>( Corr2.size() ) };
  Corr2[0].Name = EFit.Get2ptName( Src.qp );
  AppendOps( Corr2[0].Name, Src.op, Src.op );
  Corr2[1].Name = EFit.Get2ptName( Snk.qp );
  AppendOps( Corr2[1].Name, Snk.op, Snk.op );
  // Must put both in the cache before trying to dereference either
  for( int i = 0; i < NumC2; ++i )
  {
    Corr2[i].Name = Common::MakeFilename( Corr2[i].Name, Common::sFold, fna.Seed, DEF_FMT );
    Corr2[i].Handle = Cache2.GetIndex( Corr2[i].Name );
    if( Corr2[i].Handle == CorrCache::BadIndex )
      throw std::runtime_error( "2pt functions are required" );
  }
  const std::array<const Vector *, NumC2> pE{ &Src.E, &Snk.E };

  // Get the names of all 3pt correlators we will need
  // The first two files are the numerator
  // Followed by one or two files for the denominator. Both required for R2
  std::array<CorrT, 4> Corr3;
  static constexpr int NumC3{ static_cast<int>( Corr3.size() ) };
  static constexpr int iInitial{ 0 };
  static constexpr int iFinal{ 1 };
  constexpr int NumRatio{ 4 };
  std::vector<bool> MakeRatio( NumRatio, true );
  for( int Snk_Src = 0; Snk_Src < NumC3; ++Snk_Src )
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
    Corr3[Snk_Src].Name = fna.Dir;
    Corr3[Snk_Src].Name.append( fna.BaseShort );
    Corr3[Snk_Src].Name.append( iSnk == iInitial ? qSrc : qSnk );
    Corr3[Snk_Src].Name.append( 1, '_' );
    Corr3[Snk_Src].Name.append( iSrc == iInitial ? qSrc : qSnk );
    Corr3[Snk_Src].Name.append( fnaSuffix );
    Common::AppendGammaDeltaT( Corr3[Snk_Src].Name, iSnk==iSrc ? Common::Gamma::Algebra::GammaT : fna.Gamma[0],
                               fna.DeltaT );
    fna.AppendMomentum( Corr3[Snk_Src].Name, iSrc == iFinal ? Snk.qp.p : Src.qp.p, Src.qp.p.Name );
    fna.AppendMomentum( Corr3[Snk_Src].Name, iSnk == iInitial ? Src.qp.p : Snk.qp.p, Snk.qp.p.Name );
    Corr3[Snk_Src].Name.append( Common::Underscore );
    Corr3[Snk_Src].Name.append( iSnk == iInitial ? Src.op : Snk.op );
    Corr3[Snk_Src].Name.append( Common::Underscore );
    Corr3[Snk_Src].Name.append( iSrc == iInitial ? Src.op : Snk.op );
    Corr3[Snk_Src].Name = Common::MakeFilename( Corr3[Snk_Src].Name, fna.Type, fna.Seed, fna.Ext );
    Corr3[Snk_Src].Handle = Cache3.GetIndex( Corr3[Snk_Src].Name );
    if( Corr3[Snk_Src].Handle == CorrCache::BadIndex )
    {
      if( Snk_Src < 2 )
        throw std::runtime_error( "Forward and reversed 3pt functions are required" );
      MakeRatio[1] = false;
      Corr3[Snk_Src == 2 ? 3 : 2].Handle = CorrCache::BadIndex; // Other not required
      break;
    }
  }

  // Load 2pt functions
  int NSamples{ MaxSamples };
  const unsigned int CompareFlags{ Common::COMPAT_DISABLE_BASE };
  const Scalar * pC2[NumC2];
  for( int i = 0; i < NumC2; ++i )
  {
    Corr2[i].Corr = &Cache2[Corr2[i].Handle];
    if( i )
      Corr2[0].Corr->IsCompatible( *Corr2[i].Corr, &NSamples, CompareFlags );
    pC2[i] = (*Corr2[i].Corr)[Fold::idxCentral];
  }
  const int nT2{ Corr2[0].Corr->Nt() };

  // Load 3pt functions
  const Scalar * pSrc[NumC3];
  for( int i = 0; i < NumC3; ++i )
  {
    if( Corr3[i].Handle != CorrCache::BadIndex )
    {
      Corr3[i].Corr = &Cache3[Corr3[i].Handle];
      if( i == 0 )
      {
        Corr2[0].Corr->IsCompatible( *Corr3[i].Corr, &NSamples, CompareFlags | Common::COMPAT_DISABLE_NT );
        if( Corr2[0].Corr->NtUnfolded != Corr3[i].Corr->NtUnfolded )
          throw std::runtime_error( "2-pt NtUnfolded " + std::to_string( Corr2[0].Corr->NtUnfolded ) +
                                    " != 3-pt NtUnfolded " + std::to_string( Corr3[i].Corr->NtUnfolded ) );
      }
      else
        Corr3[0].Corr->IsCompatible( *Corr3[i].Corr, &NSamples, CompareFlags );
      pSrc[i] = (*Corr3[i].Corr)[Fold::idxCentral];
    }
  }
  const int nT3{ Corr3[0].Corr->Nt() };

  // Check compatibility with the models
  if( ZVmi.freeze != Freeze::Constant )
  {
    const unsigned int ModelCompareFlags{CompareFlags|Common::COMPAT_DISABLE_TYPE|Common::COMPAT_DISABLE_NT};
    int * pNumSamples = ZVmi.freeze == Freeze::None ? &NSamples : nullptr;
    Corr2[0].Corr->IsCompatible( ZVModelSnk, pNumSamples, ModelCompareFlags );
    Corr2[0].Corr->IsCompatible( ZVModelSrc, pNumSamples, ModelCompareFlags );
  }

  // Get the overlap coefficients we need if we're making R3
  std::array<Vector, 2> OverlapCoeff;
  if( !bAltR3 && ( MakeRatio[2] || MakeRatio[3] ) )
  {
    static const char Zero[] = "0";
    OverlapCoeff[0] = EFit.GetVector( Src.qp, Src.op + Zero );
    OverlapCoeff[1] = EFit.GetVector( Snk.qp, Snk.op + Zero );
  }

  // Ensure 3pt correlator includes timeslice DeltaT
  const int NTHalf{ nT3 / 2 };
  if( fna.DeltaT > NTHalf )
    throw std::runtime_error("DeltaT " + std::to_string( fna.DeltaT ) + " > Nt/2 " + std::to_string( NTHalf ));

  // Make somewhere to put ratios
  std::array<Fold, NumRatio> out;
  static const std::array<std::string, NumRatio> RatioNames{ "R1", "R2", "R3", "R3" }; // R3 isn't symmetric
  static const std::array<bool, NumRatio> NameReverse{ false, false, false, true };
  static const std::array<bool, NumRatio> RatioSymmetric{ true, true, false, false };
  std::vector<Scalar *> pDst( NumRatio, nullptr );
  for( int i = 0; i < NumRatio; i++ )
  {
    if( MakeRatio[i] )
    {
      out[i].resize( NSamples, fna.DeltaT + 1 );
      // Energies
      if( EFit.freeze != Freeze::Constant )
      {
        out[i].FileList.push_back( Snk.EModel.Name_.Filename );
        out[i].FileList.push_back( Src.EModel.Name_.Filename );
      }
      // Three-point names
      const int Num3ptStart{ i == 3 ? 1 : 0 };
      const int Num3pt{ i == 0 ? 2 : i == 1 ? 4 : 1 };
      for( int j = 0; j < Num3pt; ++j )
        out[i].FileList.push_back( Corr3[Num3ptStart+j].Corr->Name_.Filename );
      if( i != 1 )
      {
        // Two-point names
        for( int j = 0; j < NumC2; ++j )
          out[i].FileList.push_back( Corr2[j].Corr->Name_.Filename );
        // Z_V
        if( ZVmi.freeze != Freeze::Constant )
        {
          out[i].FileList.push_back( ZVModelSnk.Name_.Filename );
          out[i].FileList.push_back( ZVModelSrc.Name_.Filename );
        }
      }
      // Now copy the rest of the attributes
      out[i].CopyAttributes( *Corr3[0].Corr );
      out[i].NtUnfolded = nT3;
      pDst[i] = out[i][Fold::idxCentral];
    }
  }

  // Now make the ratios
  for( int idx = Fold::idxCentral; idx < NSamples; ++idx )
  {
    const int idxV{ idx - Fold::idxCentral };
    // Get the energies and compute the Correlator at Delta T with backward propagating wave subtracted
    // Allow energies to be frozen to 1 (which means don't subtract backward propagating wave)
    double  E[NumC2];
    double C2[NumC2];
    for( int i = 0; i < NumC2; ++i )
    {
      E[i] = (*pE[i])[idxV];
      C2[i] = pC2[i][fna.DeltaT];
      if( EFit.freeze == Freeze::Constant )
      {
        // If we're freezing energy to constant, it makes no sense to try to subtract half the midpoint
        if( NTHalf == fna.DeltaT )
          C2[i] *= 0.5; // In the middle of the lattice, this is needed for consistency
      }
      else
      {
        C2[i] -= 0.5 * pC2[i][NTHalf] * std::exp( - E[i] * ( NTHalf - fna.DeltaT ) );
      }
    }
    const double EProd = E[0] * E[1];
    const double C2Prod = std::abs( C2[0] * C2[1] );
    // Now compute the ratios
    for( int t = 0; t <= fna.DeltaT; ++t )
    {
      if( MakeRatio[0] || MakeRatio[1] )
      {
        const double n = std::abs( EProd * pSrc[0][t] * pSrc[1][t] ); // fna.DeltaT - t
        if( !std::isfinite( n ) )
          throw std::runtime_error( "Numerator Overflow" );
        if( MakeRatio[0] )
          pDst[0][t] = 2 * std::sqrt( ( n / C2Prod ) * ZVSrc[idxV] * ZVSnk[idxV] ); // R1
        if( MakeRatio[1] )
          pDst[1][t] = 2 * std::sqrt( n / std::abs( pSrc[2][t] * pSrc[3][t] ) );
      }
      if( bAltR3 )
      {
        // This is the alternate R3 - not the one we expect to use
        if( MakeRatio[2] )
        {
          // R3 - forward  (1st propagator)
          const Scalar R3Exp{ std::exp( ( E[0] - E[1] ) * t + E[1] * fna.DeltaT ) };
          const Scalar R3SqrtNum{ EProd * R3Exp * ZVSrc[idxV] * ZVSnk[idxV] };
          pDst[2][t] = 2 * pSrc[0][t] * std::sqrt( R3SqrtNum / ( pC2[0][t] * pC2[1][fna.DeltaT - t] ) );
        }
        if( MakeRatio[3] )
        {
          // R3 backward (2nd propagator)
          const Scalar R3Exp{ std::exp( ( E[1] - E[0] ) * t + E[0] * fna.DeltaT ) };
          const Scalar R3SqrtNum{ EProd * R3Exp * ZVSrc[idxV] * ZVSnk[idxV] };
          pDst[3][t] = 2 * pSrc[1][t] * std::sqrt( R3SqrtNum / ( pC2[1][t] * pC2[0][fna.DeltaT - t] ) );
        }
      }
      else
      {
        // This is the R3 we expect to use
        const Scalar OverlapNorm{ Overlap2ENorm ? 1 : 0.25 / EProd };
        const Scalar OverProd{ OverlapCoeff[0][idxV] * OverlapCoeff[1][idxV]
                              * std::sqrt( OverlapNorm * ZVSrc[idxV] * ZVSnk[idxV] ) };
        if( MakeRatio[2] ) // R3 - forward  (1st propagator)
          pDst[2][t] = OverProd * pSrc[0][t] / ( pC2[0][t] * pC2[1][fna.DeltaT - t] );
        if( MakeRatio[3] ) // R3 backward (2nd propagator)
          pDst[3][t] = OverProd * pSrc[1][t] / ( pC2[1][t] * pC2[0][fna.DeltaT - t] );
      }
    }
    // Symmetrise the symmetric ratios, i.e. R1 and R2, but not R3
    if( bSymmetrise )
    {
      assert( ! ( fna.DeltaT & 1 ) && "DeltaT must be even" );
      for( int i = 0; i < NumRatio; ++i )
      {
        if( MakeRatio[i] && RatioSymmetric[i] )
        {
          for( int t = 0; t < fna.DeltaT / 2; ++t )
          {
            pDst[i][t] = ( pDst[i][t] + pDst[i][fna.DeltaT - t] ) / 2;
            pDst[i][fna.DeltaT - t] = pDst[i][t];
          }
        }
      }
    }
    for( int i = 0; i < NumRatio; ++i )
      pDst[i] += out[i].Nt();
    for( int i = 0; i < NumC3; ++i )
      pSrc[i] += nT3;
    pC2[0] += nT2;
    pC2[1] += nT2;
  }
  for( int i = 0; i < NumRatio ; i++ )
  {
    if( MakeRatio[i] )
    {
      out[i].MakeCorrSummary( nullptr );
      std::string OutFileName{ outBase };
      OutFileName.append( RatioNames[i] );
      OutFileName.append( 1, '_' );
      OutFileName.append( NameReverse[i] ? qSrc : qSnk );
      OutFileName.append( 1, '_' );
      OutFileName.append( NameReverse[i] ? qSnk : qSrc );
      OutFileName.append( fnaSuffix );
      Common::AppendGammaDeltaT( OutFileName, fna.Gamma[0], fna.DeltaT );
      fna.AppendMomentum( OutFileName, Src.qp.p ? Src.qp.p : Snk.qp.p, Src.qp.p.Name );
      AppendOps( OutFileName, NameReverse[i] ? Src.op : Snk.op, NameReverse[i] ? Snk.op : Src.op );
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
}

Maker::Maker( const Common::CommandLine &cl )
: MaxSamples{cl.SwitchValue<int>("n")},
  modelBase{cl.SwitchValue<std::string>("im")},
  outBase{cl.SwitchValue<std::string>("o")},
  bSymmetrise{!cl.GotSwitch("nosym")},
  Overlap2ENorm{cl.GotSwitch("2e")},
  RegExSwap{cl.GotSwitch("swap")},
  RegExExt{std::regex( cl.SwitchValue<std::string>("ssre"), std::regex::extended | std::regex::icase )},
  EFit{modelBase},
  Cache2{cl.SwitchValue<std::string>("i2"), LoadFilePrefix },
  Cache3{ "", LoadFilePrefix }
{
  Common::MakeAncestorDirs( outBase );
  EFit.Read( cl.SwitchValue<std::string>( "efit" ), LoadFilePrefix );
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
      {"n", CL::SwitchType::Single, "0"},
      {"o", CL::SwitchType::Single, "" },
      {"type", CL::SwitchType::Single, DefaultType },
      {"ssre", CL::SwitchType::Single, DefaultERE },
      {"efit", CL::SwitchType::Single, ""},
      {"2e", CL::SwitchType::Flag, nullptr},
      {"swap", CL::SwitchType::Flag, nullptr},
      {"zv", CL::SwitchType::Flag, nullptr},
      {"nosym", CL::SwitchType::Flag, nullptr},
      {"r3a", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      // Read the list of fits I've chosen to use
      std::string TypeParams{ cl.SwitchValue<std::string>("type") };
      const std::string Type{ Common::ExtractToSeparator( TypeParams ) };
      std::unique_ptr<Maker> m;
      if( Common::EqualIgnoreCase( Type, "R" ) )
        m.reset( new RMaker( TypeParams, cl ) );
      else if( Common::EqualIgnoreCase( Type, "ZV" ) )
        m.reset( new ZVMaker( TypeParams, cl ) );
      else
        throw std::runtime_error( "I don't know how to make type " + Type );
      bShowUsage = false;
      const std::string Prefix3pt{ cl.SwitchValue<std::string>("i3") };
      std::vector<std::string> FileList{ Common::glob( cl.Args.begin(), cl.Args.end(), Prefix3pt.c_str() ) };
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
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) <<
    "usage: " << cl.Name << " <options> RatioFile1 [RatioFile2 [...]]\n"
    "Create 3pt ratios using E0 from fits in FitListFile. <options> are:\n"
    "--i2   Input2 filename prefix\n"
    "--i3   Input3 filename prefix\n"
    "--im   Input model filename prefix\n"
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-o     Output filename prefix\n"
    "--type Ratio type[,options[,...]] (default " << DefaultType << ")\n"
    "       ZV Make ZV (no options)\n"
    "       R  Make R ratios. options as per --efit, except a) Fit_list->ZV_List\n"
    "          b) Spectator not required\n"
    "--efit Fit_list[,options[,...]] List of fit files for energies,\n"
    "         each row has 'Quark file' (one row per momenta).\n"
    "         First option can be 'c' to freeze values to central value,\n"
    "         then comma separated list of Key=Key_in_file pairs to rename keys.\n"
    "       If Fit_list='', then first option=spectator and is required,\n"
    "         then comma separated list of Variable=Value pairs follows.\n"
    "         Any constants not mentioned are frozen to 1\n"
    "--ssre Extended regex for sink/source type, default\n       " << DefaultERE << "\n"
    "       http://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap09.html\n"
    //"-n     Number of samples to fit, 0 (default) = all available from bootstrap\n"
    "Flags:\n"
    "--2e    Overlap coefficients are already normalised by 1/sqrt{2E_0}\n"
    "--swap  Swap source / sink order in regex\n"
    "--nosym Disable Out[t]=(Out[t]+Out[deltaT-t])/2 for symmetric ratios\n"
    "--r3a   Use alternate definition of R3 (i.e. R3a)\n"
    "--help  This message\n";
  }
  return iReturn;
}
