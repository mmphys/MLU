/**
 
 Create Ratios, e.g. R1, R2 and associated values such as Z_V

 Source file: CRatio.cpp

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
#include "CRatio.hpp"
#include <MLU/DebugInfo.hpp>

#define _USE_MATH_DEFINES
#include <cmath>

static const char szReading[] = " Reading ";
MLU::DispersionType dispType;

Model DummyModel{};
const char LoadFilePrefix[] = "  ";
static Scalar ConstantOne{ 1 };

std::ostream &operator<<( std::ostream &os, const MP &mp )
{
  os << "meson=" << mp.Meson << ", " << mp.pName;
  if( mp.p.bp2 )
    os << "^2=" << mp.p.p2();
  else
    os << "=" << mp.p.to_string3d( MLU::Space );
  return os;
}

MP MPReader::operator()( const std::string &Filename ) const
{
  MP mp;
  mp.Meson = *this;
  mp.pName = MLU::Momentum::DefaultPrefix;
  mp.p.Extract( mp.Meson, mp.pName );
  return mp;
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
  std::string Fullname;
  if( Filename.empty() || Filename[0] != '/' )
    Fullname = Base;
  Fullname.append( Filename );
  if( !MLU::FileExists( Fullname ) )
  {
    // File doesn't exist
    if( printPrefix )
      std::cout << printPrefix << "- missing - " << Fullname << std::endl;
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

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename S>
void KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, S>::clear()
{
  model.clear();
  KeyMap.clear();
  freeze = Freeze::Constant;
  FieldMap.clear();
  ConstantMap.clear();
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename S>
template<typename K, typename R>
typename std::enable_if<std::is_same<K, R>::value, std::map<Key, std::string, LessKey>>::type
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, S>::ReadNameMap( const std::string &Filename, const char *p )
{
  return MLU::KeyValReader<Key, std::string, LessKey>::Read( Filename, p );
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename S>
template<typename K, typename R>
typename std::enable_if<!std::is_same<K, R>::value, std::map<Key, std::string, LessKey>>::type
KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, S>::ReadNameMap( const std::string &Filename, const char *p )
{
  // Read the file list
  using KeyNameRead = std::multimap<KeyRead, std::string, LessKeyRead>;
  using StringMapCIReader = MLU::KeyValReader<KeyRead, std::string, LessKeyRead, KeyNameRead>;
  KeyNameRead ReadList{ StringMapCIReader::Read( Filename, p ) };
  // Now make a unique key for each row from name and filename
  KeyReadMapT m;
  for( const auto &k : ReadList )
    m.emplace( k.first( k.second ), k.second );
  return m;
}

template<typename Key, typename LessKey, typename KeyRead, typename LessKeyRead, typename S>
void KeyFileCache<Key, LessKey, KeyRead, LessKeyRead, S>::Read(const std::string &filename, const char *Prefix)
{
  clear();
  std::string sParams{ filename };
  std::string Filename{ MLU::ExtractToSeparator( sParams ) };
  if( Filename.empty() )
  {
    freeze = Freeze::Constant;
    ConstantMap = MLU::KeyValReader<std::string, Scalar>::Read( MLU::ArrayFromString( sParams ) );
    return;
  }

  freeze = Freeze::None;
  {
    // See whether next parameter is 'C'
    std::string RemainingParams{ sParams };
    std::string Options{ MLU::ExtractToSeparator( RemainingParams ) };
    if( Options.length() == 1 && std::toupper( sParams[0] ) == 'C' )
    {
      freeze = Freeze::Central;
      sParams = RemainingParams;
    }
  }
  FieldMap = MLU::KeyValReader<std::string, std::string>::Read( MLU::ArrayFromString( sParams ) );

  static const char szDescription[]{ "2pt fit list" };
  std::cout << "Reading " << szDescription << " from " << Filename << MLU::NewLine;
  KeyReadMapT m{ ReadNameMap( Filename, szDescription ) };
  // All the files in the cache are relative to the directory containing the fit list
  model.SetBase( MLU::GetDirPrefix( Filename ) );
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
template<typename Key, typename L, typename KR, typename LKR, typename S>
int KeyFileCache<Key,L,KR,LKR,S>::GetIndex( const Key &key ) const
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

template<typename Key, typename L, typename KR, typename LKR, typename S>
MLU::Model<S> &KeyFileCache<Key,L,KR,LKR,S>::operator()( const Key &key, const char * PrintPrefix )
{
  if( freeze==Freeze::Constant )
    return DummyModel;
  return model( GetIndex( key ), PrintPrefix );
}

template<typename Key, typename L, typename KR, typename LKR, typename S>
MLU::JackBootColumn<S> KeyFileCache<Key,L,KR,LKR,S>::GetColumn( const Key &key, const MLU::Param::Key &pKey, std::size_t Index )
{
  if( freeze == Freeze::Constant )
  {
    // Freeze to value specified on command line, or 1 if not specified
    typename ConstantMapT::iterator it = ConstantMap.find( pKey.FullName( Index ) );
    return MLU::JackBootColumn<S>( it == ConstantMap.end() ? ConstantOne : it->second );
  }
  else
  {
    M &m{ operator[]( key ) };
    // TODO: Reinstate field mapping (not quite sure how to define it?)
    //typename FieldMapT::iterator it = FieldMap.find( Name );
    //Scalar * p{ m[M::idxCentral] + m.GetColumnIndex( it == FieldMap.end() ? Name : it->second ) };
    if( freeze == Freeze::Central )
      return m.ColumnFrozen( pKey, Index );
    return m.Column( pKey, Index );
  }
}

void Maker::Run( std::size_t &NumOK, std::size_t &Total, std::vector<std::string> &FileList )
{
  Total += FileList.size();
  for( std::string &sFile : FileList )
  {
    try
    {
      Make( sFile );
      ++NumOK;
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
}

// Incoming file should be a three-point function with a heavy sink and light(er) source
void ZVRCommon::Make( std::string &FileName )
{
  std::cout << "Processing " << FileName << MLU::NewLine;
  std::vector<std::string> OpNames;
  MLU::FileNameAtt fna{ FileName, &OpNames };
  if( OpNames.empty() )
    throw std::runtime_error( "Sink / source operator names mising" );
  if( !fna.GotDeltaT() )
    throw std::runtime_error( "DeltaT missing" );
  if( fna.Gamma.size() != 1 )
    throw std::runtime_error( FileName + " has " + std::to_string( fna.Gamma.size() ) + " currents" );
  if( fna.MesonMom.size() < 2 )
    throw std::runtime_error( "Can't extract sink/source mesons/momenta from " + fna.GetBaseShort() );
  const std::string FileNameSuffix{ fna.GetBaseShort( 3 ) };
  SSInfo Snk( EFit, OpNames[fna.op[1]], MP( fna.Meson[1], fna.MesonP[1] ), fna.BaseShortParts[1] );
  SSInfo Src( EFit, OpNames[fna.op[0]], MP( fna.Meson[0], fna.MesonP[0] ), fna.BaseShortParts[2] );
  fna.BaseShortParts.resize( 1 );
  ZVRMake( fna, FileNameSuffix, Snk, Src );
}

ZVMaker::ZVMaker( const std::string &TypeParams, const MLU::CommandLine &cl ) : ZVRCommon( cl )
{
  if( !TypeParams.empty() )
    throw std::runtime_error( "ZVMaker does not recognise any parameters: " + TypeParams );
}

void ZVMaker::ZVRMake( const MLU::FileNameAtt &fna, const std::string &fnaSuffix,
                       const SSInfo &Snk, const SSInfo &Src )
{
  // Complain about errors
  {
    std::ostringstream os;
    if( MLU::CompareIgnoreCase( Snk.q, Src.q ) )
      os << "qSnk " << Snk.q << " != qSrc " << Src.q;
    else if( fna.Gamma[0] != MLU::Gamma::Algebra::GammaT )
      os << "Gamma " << fna.Gamma[0] << " != " << MLU::Gamma::Algebra::GammaT;
    else if( Snk.mp.p )
      os << "Sink momentum " << Snk.mp.p.FileString( Snk.mp.pName );
    else if( Src.mp.p )
      os << "Source momentum " << Src.mp.p.FileString( Src.mp.pName );
    if( os.tellp() != std::streampos( 0 ) )
    {
      os << " file " << fna.Filename;
      throw std::runtime_error( os.str().c_str() );
    }
  }

  // Read the three-point correlator
  Fold C3;
  C3.Read( fna.Filename, LoadFilePrefix );

  // Ensure the correlator includes timeslice DeltaT
  const int NTHalf{ C3.Nt() / 2 };
  if( fna.DeltaT[0] > NTHalf )
    throw std::runtime_error( "DeltaT " + std::to_string( fna.DeltaT[0] ) + " > Nt/2 " + std::to_string( NTHalf ) );
  const bool bMidpoint{ fna.DeltaT[0] == NTHalf };

  // Now read the two-point function corresponding to the model (Src and Snk are the same)
  std::string C2Name{ Src.mp.C2Name() }; // Src and Snk same, so we can choose either

  // DON'T DO THIS - IT JUST MAKES THE ERRORS LARGER
  // Wall sink - Point source ==> swap source and sink on the 2-point for accuracy
  //if( std::toupper( opSnk.back() ) == 'W' && std::toupper( opSrc.back() ) == 'P' )
    //AppendOps( C2Name, opSrc, opSnk );
  //else
    AppendOps( C2Name, Snk.op, Src.op );
  C2Name = MLU::ExistingFilename( C2Name, MLU::sFold, fna.Seed, DEF_FMT );
  Fold &C2{ Cache2[C2Name] };
  //C2.Read( C2Name, LoadFilePrefix );
  int NumSamples {};
  C3.IsCompatible( C2, &NumSamples, MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
  if( C2.NtUnfolded() != C3.Nt() )
  {
    std::ostringstream os;
    os << "C2.NtUnfolded " << C2.NtUnfolded() << " != C3.Nt() " << C3.Nt();
    throw std::runtime_error( os.str() );
  }

  // Make somewhere to put Z_V
  Fold out( NumSamples, fna.DeltaT[0] + 1 );
  out.FileList.push_back( C3.Name_.Filename );
  out.FileList.push_back( C2.Name_.Filename );
  if( EFit.freeze != Freeze::Constant )
    out.FileList.push_back( Src.EModel.Name_.Filename );
  out.CopyAttributes( C3 );
  out.NtUnfolded_ = C3.Nt();
  for( int idx = Fold::idxCentral; idx < NumSamples; ++idx )
  {
    const Scalar CTilde{ bMidpoint ? C2(idx,fna.DeltaT[0]) * 0.5
                        : C2(idx,fna.DeltaT[0]) - Src.A[idx] * Src.A[idx]
                        * std::exp( - Src.E[idx] * ( C3.Nt() - fna.DeltaT[0] ) )};
    for( int t = 0; t <= fna.DeltaT[0]; ++t )
    {
      const double z{ CTilde / C3(idx,t) };
      if( !std::isfinite( z ) )
        throw std::runtime_error( "ZV Overflow" );
      out(idx,t) = z;
    }
  }
  // Make output file
  out.MakeCorrSummary();
  std::string OutFileName{ outBase };
  OutFileName.append( "ZV_" );
  OutFileName.append( Snk.q );
  //OutFileName.append( Suffix );
  MLU::AppendDeltaT( OutFileName, fna.DeltaT[0] );
  AppendOps( OutFileName, Snk.op, Src.op );
  OutFileName = MLU::MakeFilename( OutFileName, MLU::sFold, C3.Seed(), DEF_FMT );
  std::cout << "->" << OutFileName << MLU::NewLine;
  out.Write( OutFileName, MLU::sFold.c_str() );
  MLU::ReplaceExtension( OutFileName, TEXT_EXT );
  out.WriteSummary( OutFileName );
}

RMaker::RMaker( const std::string &TypeParams, const MLU::CommandLine &cl )
: ZVRCommon( cl ), bAltR3{cl.GotSwitch("r3a")}
{
  ZVmi.Read( TypeParams, LoadFilePrefix );
}

void RMaker::ZVRMake( const MLU::FileNameAtt &fna, const std::string &fnaSuffix,
                      const SSInfo &Snk, const SSInfo &Src )
{
  // Make sure I have the model for Z_V already loaded
  Model &ZVModelSnk{ ZVmi( Snk.q, LoadFilePrefix ) };
  Model &ZVModelSrc{ ZVmi( Src.q, LoadFilePrefix ) };
  const std::vector<MLU::Param::Key> vZVKeys{ MLU::Param::Key( Src.q, "ZV" ),
                                                 MLU::Param::Key( Snk.q, "ZV" ) };
  const MLU::Param::Key &kZVSrc( vZVKeys[0] );
  const MLU::Param::Key &kZVSnk{ vZVKeys[1] };
  Column ZVSnk{ ZVmi.GetColumn( Snk.q, kZVSnk ) };
  Column ZVSrc{ ZVmi.GetColumn( Src.q, kZVSrc ) };

  // Get the names of all 2pt correlators we will need
  std::array<CorrT, 2> Corr2;
  static constexpr int NumC2{ static_cast<int>( Corr2.size() ) };
  Corr2[0].Name = Src.mp.C2Name();
  AppendOps( Corr2[0].Name, Src.op, Src.op );
  Corr2[1].Name = Snk.mp.C2Name();
  AppendOps( Corr2[1].Name, Snk.op, Snk.op );
  // Must put both in the cache before trying to dereference either
  for( int i = 0; i < NumC2; ++i )
  {
    Corr2[i].Name = MLU::ExistingFilename( Corr2[i].Name, MLU::sFold, fna.Seed, DEF_FMT );
    Corr2[i].Handle = Cache2.GetIndex( Corr2[i].Name );
    if( Corr2[i].Handle == CorrCache::BadIndex )
      throw std::runtime_error( "2pt functions are required" );
  }
  const std::array<const Column *, NumC2> pE{ &Src.E, &Snk.E };

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
    Corr3[Snk_Src].Name.append( fna.GetBaseShort() );
    MLU::Append( Corr3[Snk_Src].Name, iSnk == iInitial ? Src.q : Snk.q );
    MLU::Append( Corr3[Snk_Src].Name, iSrc == iInitial ? Src.q : Snk.q );
    Corr3[Snk_Src].Name.append( fnaSuffix );
    MLU::AppendGammaDeltaT( Corr3[Snk_Src].Name, iSnk==iSrc ? MLU::Gamma::Algebra::GammaT : fna.Gamma[0],
                               fna.DeltaT[0] );
    fna.AppendMomentum( Corr3[Snk_Src].Name, iSrc == iFinal ? Snk.mp.p : Src.mp.p, Src.mp.pName );
    fna.AppendMomentum( Corr3[Snk_Src].Name, iSnk==iInitial ? Src.mp.p : Snk.mp.p, Snk.mp.pName );
    MLU::Append( Corr3[Snk_Src].Name, iSnk == iInitial ? Src.op : Snk.op );
    MLU::Append( Corr3[Snk_Src].Name, iSrc == iInitial ? Src.op : Snk.op );
    Corr3[Snk_Src].Name = MLU::ExistingFilename(Corr3[Snk_Src].Name, fna.Type, fna.Seed, fna.Ext);
    Corr3[Snk_Src].Handle = Cache3.GetIndex( Corr3[Snk_Src].Name );
    if( Corr3[Snk_Src].Handle == CorrCache::BadIndex )
    {
      if( Snk_Src == 1 && Corr3[0].Handle == CorrCache::BadIndex )
        throw std::runtime_error( "At least one (forward or reversed) 3pt numerator required" );
      MakeRatio[1] = false; // R2 needs every correlator
      if( Snk_Src < 2 )
      {
        // One of our numerators are missing
        MakeRatio[0] = false; // R1 needs both numerators
        MakeRatio[2 + Snk_Src] = false; // Don't make the R3 we don't have
      }
      else
      {
        // One of our 3pt denominators are missing. Ignore the other, since both are required for R2
        Corr3[Snk_Src == 2 ? 3 : 2].Handle = CorrCache::BadIndex; // Other not required
      }
      break;
    }
  }

  // Load 2pt functions
  int NSamples{ MaxSamples };
  const unsigned int CompareFlags{ MLU::COMPAT_DISABLE_BASE };
  const Fold * pC2[NumC2];
  for( int i = 0; i < NumC2; ++i )
  {
    Corr2[i].Corr = &Cache2[Corr2[i].Handle];
    if( i )
      Corr2[0].Corr->IsCompatible( *Corr2[i].Corr, &NSamples, CompareFlags );
    pC2[i] = Corr2[i].Corr;
  }
  const int nT2{ Corr2[0].Corr->Nt() };

  // Load 3pt functions
  const Fold * pSrc[NumC3];
  for( int i = 0; i < NumC3; ++i )
  {
    if( Corr3[i].Handle != CorrCache::BadIndex )
    {
      Corr3[i].Corr = &Cache3[Corr3[i].Handle];
      if( i == 0 )
      {
        Corr2[0].Corr->IsCompatible( *Corr3[i].Corr, &NSamples, CompareFlags | MLU::COMPAT_DISABLE_NT );
        if( Corr2[0].Corr->NtUnfolded() != Corr3[i].Corr->NtUnfolded() )
          throw std::runtime_error( "2-pt NtUnfolded "
                                   + std::to_string( Corr2[0].Corr->NtUnfolded() )
                                   + " != 3-pt NtUnfolded "
                                   + std::to_string( Corr3[i].Corr->NtUnfolded() ) );
      }
      else
        Corr3[0].Corr->IsCompatible( *Corr3[i].Corr, &NSamples, CompareFlags );
      pSrc[i] = Corr3[i].Corr;
    }
  }
  const int nT3{ Corr3[0].Corr->Nt() };

  // Check compatibility with the models
  if( ZVmi.freeze != Freeze::Constant )
  {
    const unsigned int ModelCompareFlags{CompareFlags|MLU::COMPAT_DISABLE_TYPE|MLU::COMPAT_DISABLE_NT};
    int * pNumSamples = ZVmi.freeze == Freeze::None ? &NSamples : nullptr;
    Corr2[0].Corr->IsCompatible( ZVModelSnk, pNumSamples, ModelCompareFlags );
    Corr2[0].Corr->IsCompatible( ZVModelSrc, pNumSamples, ModelCompareFlags );
  }

  // Get the overlap coefficients we need if we're making R3
  std::array<Column, 2> OverlapCoeff;
  if( !bAltR3 && ( MakeRatio[2] || MakeRatio[3] ) )
  {
    OverlapCoeff[0] = EFit.GetColumn( Src.mp, Src.mp.PK( Src.op ) );
    OverlapCoeff[1] = EFit.GetColumn( Snk.mp, Snk.mp.PK( Snk.op ) );
  }

  // Ensure 3pt correlator includes timeslice DeltaT
  const int NTHalf{ nT3 / 2 };
  if( fna.DeltaT[0] > NTHalf )
    throw std::runtime_error("DeltaT " + std::to_string( fna.DeltaT[0] ) + " > Nt/2 " + std::to_string( NTHalf ));

  // Make somewhere to put ratios
  std::array<Fold, NumRatio> out;
  static const std::array<std::string, NumRatio> RatioNames{ "R1", "R2", "R3", "R3" }; // R3 isn't symmetric
  static const std::array<bool, NumRatio> NameReverse{ false, false, false, true };
  static const std::array<bool, NumRatio> RatioSymmetric{ true, true, false, false };
  for( int i = 0; i < NumRatio; i++ )
  {
    if( MakeRatio[i] )
    {
      out[i].resize( NSamples, fna.DeltaT[0] + 1 );
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
      out[i].NtUnfolded_ = nT3;
    }
  }

  // Now make the ratios
  for( int idx = Fold::idxCentral; idx < NSamples; ++idx )
  {
    // Get the energies and compute the Correlator at Delta T with backward propagating wave subtracted
    // Allow energies to be frozen to 1 (which means don't subtract backward propagating wave)
    double  E[NumC2];
    double C2[NumC2];
    for( int i = 0; i < NumC2; ++i )
    {
      E[i] = (*pE[i])[idx];
      C2[i] = (*pC2[i])(idx,fna.DeltaT[0]);
      if( EFit.freeze == Freeze::Constant )
      {
        // If we're freezing energy to constant, it makes no sense to try to subtract half the midpoint
        if( NTHalf == fna.DeltaT[0] )
          C2[i] *= 0.5; // In the middle of the lattice, this is needed for consistency
      }
      else
      {
        C2[i] -= 0.5 * (*pC2[i])(idx,NTHalf) * std::exp( - E[i] * ( NTHalf - fna.DeltaT[0] ) );
      }
    }
    const double EProd = E[0] * E[1];
    const double C2Prod = std::abs( C2[0] * C2[1] );
    // Now compute the ratios
    for( int t = 0; t <= fna.DeltaT[0]; ++t )
    {
      if( MakeRatio[0] || MakeRatio[1] )
      {
        const double n = std::abs( EProd * (*pSrc[0])(idx,t) * (*pSrc[1])(idx,t) ); // fna.DeltaT - t
        if( !std::isfinite( n ) )
          throw std::runtime_error( "Numerator Overflow" );
        if( MakeRatio[0] )
          out[0](idx,t) = 2 * std::sqrt( ( n / C2Prod ) * ZVSrc[idx] * ZVSnk[idx] ); // R1
        if( MakeRatio[1] )
          out[1](idx,t) = 2 * std::sqrt( n / std::abs( (*pSrc[2])(idx,t) * (*pSrc[3])(idx,t) ) );
      }
      if( bAltR3 )
      {
        // This is the alternate R3 - not the one we expect to use
        if( MakeRatio[2] )
        {
          // R3 - forward  (1st propagator)
          const Scalar R3Exp{ std::exp( ( E[0] - E[1] ) * t + E[1] * fna.DeltaT[0] ) };
          const Scalar R3SqrtNum{ EProd * R3Exp * ZVSrc[idx] * ZVSnk[idx] };
          out[2](idx,t) = 2 * (*pSrc[0])(idx,t) * std::sqrt( R3SqrtNum / ( (*pC2[0])(idx,t) * (*pC2[1])(idx,fna.DeltaT[0] - t) ) );
        }
        if( MakeRatio[3] )
        {
          // R3 backward (2nd propagator)
          const Scalar R3Exp{ std::exp( ( E[1] - E[0] ) * t + E[0] * fna.DeltaT[0] ) };
          const Scalar R3SqrtNum{ EProd * R3Exp * ZVSrc[idx] * ZVSnk[idx] };
          out[3](idx,t) = 2 * (*pSrc[1])(idx,t) * std::sqrt( R3SqrtNum / ( (*pC2[1])(idx,t) * (*pC2[0])(idx,fna.DeltaT[0] - t) ) );
        }
      }
      else
      {
        // This is the R3 we expect to use
        const Scalar OverlapNorm{ bOverlapAltNorm ? ( 0.25 / EProd ) : 1 };
        const Scalar OverProd{ OverlapCoeff[0][idx] * OverlapCoeff[1][idx]
                              * std::sqrt( OverlapNorm * ZVSrc[idx] * ZVSnk[idx] ) };
        if( MakeRatio[2] ) // R3 - forward  (1st propagator)
          out[2](idx,t) = OverProd * (*pSrc[0])(idx,t) / ( (*pC2[0])(idx,t) * (*pC2[1])(idx,fna.DeltaT[0] - t) );
        if( MakeRatio[3] ) // R3 backward (2nd propagator)
          out[3](idx,t) = OverProd * (*pSrc[1])(idx,t) / ( (*pC2[1])(idx,t) * (*pC2[0])(idx,fna.DeltaT[0] - t) );
      }
    }
    // Symmetrise the symmetric ratios, i.e. R1 and R2, but not R3
    if( bSymmetrise )
    {
      assert( ! ( fna.DeltaT[0] & 1 ) && "DeltaT must be even" );
      for( int i = 0; i < NumRatio; ++i )
      {
        if( MakeRatio[i] && RatioSymmetric[i] )
        {
          for( int t = 0; t < fna.DeltaT[0] / 2; ++t )
          {
            out[i](idx,t) = ( out[i](idx,t) + out[i](idx,fna.DeltaT[0] - t) ) / 2;
            out[i](idx,fna.DeltaT[0] - t) = out[i](idx,t);
          }
        }
      }
    }
  }
  for( int i = 0; i < NumRatio ; i++ )
  {
    if( MakeRatio[i] )
    {
      out[i].MakeCorrSummary();
      std::string OutFileName{ outBase };
      OutFileName.append( RatioNames[i] );
      OutFileName.append( 1, '_' );
      OutFileName.append( NameReverse[i] ? Src.q : Snk.q );
      OutFileName.append( 1, '_' );
      OutFileName.append( NameReverse[i] ? Snk.q : Src.q );
      OutFileName.append( fnaSuffix );
      MLU::AppendGammaDeltaT( OutFileName, fna.Gamma[0], fna.DeltaT[0] );
      fna.AppendMomentum( OutFileName, Src.mp.p ? Src.mp.p : Snk.mp.p, Src.mp.pName );
      AppendOps( OutFileName, NameReverse[i] ? Src.op : Snk.op, NameReverse[i] ? Snk.op : Src.op );
      OutFileName = MLU::MakeFilename( OutFileName, fna.Type, Corr2[0].Corr->Seed(), fna.Ext );
      std::cout << "->" << OutFileName << MLU::NewLine;
      out[i].Write( OutFileName, MLU::sFold.c_str() );
      // Now save ZV together with the ratio
      if( !ZVSrc.IsIdentity() || !ZVSnk.IsIdentity() )
      {
        static const std::string sError{ "RMaker::ZVRMake() writing ZV" };
        Model m( NSamples, MLU::Params( vZVKeys ), {} );
        m.FileList.push_back( ZVModelSrc.Name_.Filename );
        m.FileList.push_back( ZVModelSnk.Name_.Filename );
        m.CopyAttributes( *Corr3[0].Corr );
        std::size_t iZVSrc{ m.params.Find( kZVSrc, sError )->second() };
        std::size_t iZVSnk{ m.params.Find( kZVSnk, sError )->second() };
        for( int rep = Model::idxCentral; rep < NSamples; ++rep )
        {
          m(rep,iZVSrc) = ZVSrc[rep];
          m(rep,iZVSnk) = ZVSnk[rep];
        }
        m.Write( OutFileName, MLU::sModel.c_str(), H5F_ACC_RDWR );
      }
      MLU::ReplaceExtension( OutFileName, TEXT_EXT );
      out[i].WriteSummary( OutFileName );
    }
  }
}

const std::vector<std::string> FormFactor::ParamNames{ "EL", "mL", "mH", "qSq", "kMu",
  "melV0", "melVi", "fPar", "fPerp", "fPlus", "f0", "ELLat", "qSqLat",
  "tPlus", "tMinus", "t0", "zre", "zim", "ZV", "melV0Raw", "melViRaw" };

FormFactor::FormFactor( std::string TypeParams )
: N{ MLU::FromString<unsigned int>( MLU::ExtractToSeparator( TypeParams ) ) },
  ap{ 2 * M_PI / N },
  apInv{ 1. / ap },
  bAdjustGammaSpatial{ MLU::EqualIgnoreCase( MLU::ExtractToSeparator(TypeParams), "spatial" ) }
{
  // optional ZV: Can be constant (number), or cache of ZV files
  std::string sZV;
  if( !TypeParams.empty() )
    sZV = MLU::ExtractToSeparator( TypeParams );
  if( sZV.empty() )
  {
    ZV[0] = 1.;
    ZV[1] = 1.;
  }
  else if( MLU::FileExists( sZV ) )
  {
    ZVmi.Read( sZV, LoadFilePrefix );
    if( ZVmi.empty() )
      throw std::runtime_error( "No ZV fits in " + sZV );
  }
  else
  {
    const Scalar s{ MLU::FromString<Scalar>( sZV ) };
    ZV[0] = s;
    ZV[1] = s;
  }
  // optional: Filename for Fully non-perturbative ZV correction Z_{V,mixed}
  if( !TypeParams.empty() )
  {
    const std::string sEnsembleInfo{ MLU::ExtractToSeparator( TypeParams ) };
    if( !MLU::FileExists( sEnsembleInfo ) )
      throw std::runtime_error( "Ensemble info file doesn't exist " + sEnsembleInfo );
    try { mEnsembleInfo.SetName( sEnsembleInfo ); } catch ( const std::runtime_error &e ) {}
    try
    {
      mEnsembleInfo.Read( "Z_{V,mixed} from " );
    }
    catch( const ::H5::Exception &e )
    {
      throw std::runtime_error( "HDF5 error reading ensemble info file " + sEnsembleInfo );
    }
  }
  std::cout << ( mEnsembleInfo.Nt() ? "Ful" : "Most" ) << "ly non-perturbative results.\n";
  if( !TypeParams.empty() )
    throw std::runtime_error( "Unrecognised form factor option: " + TypeParams );
}

void FormFactor::Write( std::string &OutFileName, const Model &CopyAttributesFrom,
                        std::vector<std::string> &&SourceFileNames, int NumSamples,
                        const Column &MHeavy, const Column &ELight,
                        const Column &vT, const MLU::Momentum &p,
                        const Column &ZVPrevSrc, const Column &ZVPrevSnk, const Column &ZVMixed,
                        // These only required for non-zero momentum
                        const Column *pMLight, const Column *pvXYZ )
{
  assert( Model::idxCentral == -1 && "Bug: Model::idxCentral != -1" );
  if( p )
  {
    if( !pMLight )
      throw std::invalid_argument( "pMLight must be specified for non-zero momentum" );
    if( !pvXYZ )
      throw std::invalid_argument( "pvXYZ must be specified for non-zero momentum" );
  }
  const Column &MLight{ p ? *pMLight : ELight };
  std::vector<std::size_t> vIdx;
  Model Out( NumSamples, MLU::Params( ParamNames, vIdx ), {} );
  Out.FileList = std::move( SourceFileNames );
  Out.CopyAttributes( CopyAttributesFrom );
  Out.Name_.Seed = CopyAttributesFrom.Name_.Seed;
  Out.binSize = CopyAttributesFrom.binSize;
  for( int idx = Model::idxCentral; idx < NumSamples; ++idx )
  {
    Out(idx,vIdx[qSq]) = MHeavy[idx] * MHeavy[idx] + MLight[idx] * MLight[idx] - 2 * MHeavy[idx] * ELight[idx];
    Out(idx,vIdx[mH]) = MHeavy[idx];
    Out(idx,vIdx[mL]) = MLight[idx];
    Out(idx,vIdx[EL]) = ELight[idx];
    Out(idx,vIdx[ELLat]) = p.LatticeDispersion( MLight[idx], N, dispType );
    Out(idx,vIdx[qSqLat]) = MHeavy[idx] * MHeavy[idx] + MLight[idx] * MLight[idx] - 2 * MHeavy[idx] * Out(idx,vIdx[ELLat]);
    // z-expansion (23) pg 8 https://arxiv.org/pdf/2104.09883.pdf
    Out(idx,vIdx[tPlus ]) = MHeavy[idx] + MLight[idx];
    Out(idx,vIdx[tMinus]) = MHeavy[idx] - MLight[idx];
    Out(idx,vIdx[t0    ]) = 0.5 * Out(idx,vIdx[tMinus]);
    {
      std::complex<Scalar> Left{ Out(idx,vIdx[tPlus]) - Out(idx,vIdx[qSq]) };
      Left = std::sqrt( Left );
      std::complex<Scalar> Right{ Out(idx,vIdx[tPlus]) - Out(idx,vIdx[t0]) };
      Right = std::sqrt( Right );
      const std::complex<Scalar> z{ ( Left - Right ) / ( Left + Right ) };
      Out(idx,vIdx[z_re]) = z.real();
      Out(idx,vIdx[z_im]) = z.imag();
    }
    // Renormalisation
    const Scalar ThisZV{ ( ZV[0][idx] == ZV[1][idx] ? ZV[0][idx] : std::sqrt(ZV[0][idx]*ZV[1][idx]) )
                         * ZVMixed[idx] };
    Out(idx,vIdx[idxZV]) = ThisZV * std::sqrt( ZVPrevSrc[idx] * ZVPrevSnk[idx] );
    Out(idx,vIdx[melV0Raw]) = vT[idx];
    Out(idx,vIdx[melV0]) = vT[idx] * ThisZV;
    Out(idx,vIdx[melViRaw]) = p ? (*pvXYZ)[idx] : 0;
    Out(idx,vIdx[melVi]) = p ? (*pvXYZ)[idx] * ThisZV : 0;
    // Now make form factors
    const Scalar    Root2MH{ std::sqrt( MHeavy[idx] * 2 ) };
    const Scalar InvRoot2MH{ 1. / Root2MH };
    Out(idx,vIdx[fPar]) = Out(idx,vIdx[melV0]) * InvRoot2MH;
    if( p )
    {
      Out(idx,vIdx[kMu]) = ap;
      if( bAdjustGammaSpatial )
      {
        int FirstNonZero{ p.x ? p.x : p.y ? p.y : p.z };
        if( ( p.y && p.y != FirstNonZero ) || ( p.z && p.z != FirstNonZero ) )
          throw std::runtime_error( "Can't adjust gXYZ when momentum components differ" );
        Out(idx,vIdx[melVi]) /= FirstNonZero;
      }
      Out(idx,vIdx[fPerp]) = Out(idx,vIdx[melVi]) * InvRoot2MH * apInv;
      Out(idx,vIdx[fPlus]) = ( MHeavy[idx] - ELight[idx] ) * Out(idx,vIdx[fPerp]);
      Out(idx,vIdx[f0]) = ( ELight[idx] * ELight[idx] - MLight[idx] * MLight[idx] ) * Out(idx,vIdx[fPerp]);
    }
    else
    {
      Out(idx,vIdx[kMu]) = 0;
      Out(idx,vIdx[melVi]) = 0;
      Out(idx,vIdx[fPerp]) = 0;
      Out(idx,vIdx[fPlus]) = 0;
      Out(idx,vIdx[f0]) = 0;
    }
    Out(idx,vIdx[fPlus]) += Out(idx,vIdx[fPar]);
    Out(idx,vIdx[fPlus]) *= InvRoot2MH;
    Out(idx,vIdx[f0]) += ( MHeavy[idx] - ELight[idx] ) * Out(idx,vIdx[fPar]);
    Out(idx,vIdx[f0]) *= Root2MH / ( MHeavy[idx] * MHeavy[idx] - MLight[idx] * MLight[idx] );
  }
  // Write output
  std::cout << " Writing " << OutFileName <<MLU::NewLine;
  Out.SetSummaryNames( MLU::sParams );
  Out.MakeCorrSummary();
  Out.Write( OutFileName );
  OutFileName.resize( OutFileName.length() - Out.Name_.Ext.length() );
  OutFileName.append( TEXT_EXT );
  Out.WriteSummary( OutFileName );
}

bool FMaker::RatioFile::Less::operator()( const RatioFile &lhs, const RatioFile &rhs ) const
{
  /*int i{ MLU::CompareIgnoreCase( lhs.Prefix, rhs.Prefix ) };
  if( i )
    return i < 0;*/
  int i = MLU::CompareIgnoreCase( lhs.Sink, rhs.Sink );
  if( i )
    return i < 0;
  i = MLU::CompareIgnoreCase( lhs.Source, rhs.Source );
  if( i )
    return i < 0;
  return lhs.p < rhs.p;
}

void FMaker::Make( std::string &FileName )
{
  // Parse Filename and check it contains what we need
  //std::cout << "Processing " << FileName << MLU::NewLine;
  std::vector<std::string> OpNames;
  MLU::FileNameAtt fna{ FileName, &OpNames };
  if( OpNames.empty() )
    throw std::runtime_error( "Sink / source operator names mising" );
  if( fna.Gamma.size() != 1 )
    throw std::runtime_error( FileName + " has " + std::to_string( fna.Gamma.size() ) + " currents" );
  const bool bGammaT{ fna.Gamma[0] == MLU::Gamma::Algebra::GammaT };
  if( !bGammaT && fna.Gamma[0] != MLU::Gamma::Algebra::Spatial )
  {
    std::ostringstream ss;
    ss << "Ignoring current " << fna.Gamma[0] << " - should be "
       << MLU::Gamma::Algebra::GammaT << " or " << MLU::Gamma::Algebra::Spatial << " only";
    throw std::runtime_error( ss.str().c_str() );
  }
  if( fna.p.empty() )
    throw std::runtime_error( FileName + " has no momenta" );
  if( fna.Meson.size() != 2 || fna.MesonMom.size() != 2 )
    throw std::runtime_error( FileName + " unable to decode mesons" );
  //std::cout << "  " << fna.MesonMom[0] << " -> " << fna.Gamma[0] << " -> " << fna.MesonMom[1] << MLU::NewLine;
  const MLU::Momentum &p{ fna.p.begin()->second };
  RatioFile rf( p, fna.Meson[0], fna.Meson[1] );
  std::vector<std::string> &vs{ map[bGammaT ? 0 : 1][rf] };
  vs.emplace_back( FileName );
}

void FMaker::Run( std::size_t &NumOK, std::size_t &Total, std::vector<std::string> &FileList )
{
  // Parse each filename and group by decay, momentum and gamma
  for( std::string &f : FileList )
  {
    try
    {
      Make( f ); // Parse
    }
    catch( const std::exception &e )
    {
      std::cerr << "Error: " << e.what() << std::endl;
    }
    catch( ... )
    {
      std::cerr << "Error: Unknown exception" << std::endl;
    }
  }
  // Walk the list of every prefix, decay and momentum
  for( const typename FileMap::value_type &vt : map[0] )
  {
    try
    {
      // Get information about this decay
      const RatioFile &rf{ vt.first };
      const std::vector<std::string> &vFileGammaT{ vt.second };
      std::vector<std::string> OpNames;
      MLU::FileNameAtt fna{ vFileGammaT[0], &OpNames };
      std::string sOpNames;
      fna.AppendOps( sOpNames, MLU::Period, &OpNames );
      const MLU::MomentumPair &p{ fna.GetFirstNonZeroMomentum() };
      std::string Description;
      {
        std::ostringstream os;
        os << rf.Source << " --> " << rf.Sink << p.second.FileString( p.first, MLU::Space );
        Description = os.str();
      }
      std::cout << Description << MLU::NewLine;
      const std::vector<std::string> &Quark{ fna.Quark };
      const int iLight{ MLU::QuarkWeight( Quark[0] ) <= MLU::QuarkWeight( Quark[1] ) ? 0 : 1 };
      const int iHeavy{ 1 - iLight };
      // Get ZV
      if( !ZVmi.empty() )
      {
        ZV[0] = ZVmi.GetColumn( Quark[0], MLU::Param::Key( Quark[0], "ZV" ) );
        ZV[1] = ZVmi.GetColumn( Quark[1], MLU::Param::Key( Quark[1], "ZV" ) );
      }
      // If non-zero momentum, build a list of corresponding spatial correlators
      std::vector<Model> mGammaXYZ;
      Model * pmMLight = nullptr; // Doesn't own what it points to
      Column MLight;
      if( !rf.p )
        mGammaXYZ.resize( 1 ); // A single, dummy model so the loop works
      else
      {
        // Get the zero-momentum version of the light
        const MLU::MomentumPair p0( p.first, p.second.bp2 );
        MP mpMLight( fna.Meson[iLight], p0.second, p0.first );
        pmMLight = &EFit( mpMLight, szReading );
        MLight = EFit.GetColumn( mpMLight, mpMLight.PK() );
        // Read all the gamma spatial
        typename FileMap::iterator it{ map[1].find( rf ) };
        if( it == map[1].end() )
          throw std::runtime_error( "Gamma spatial not found for " + Description );
        std::vector<std::string> SpatialNames{ std::move( it->second ) };
        map[1].erase( it );
        for( std::string &File : SpatialNames )
        {
          try
          {
            Model m;
            m.Read( File, " reading spatial model ", &OpNames );
            mGammaXYZ.emplace_back( std::move( m ) );
          }
          catch( const std::exception &e )
          {
            std::cerr << "Error: " << e.what() << std::endl;
          }
        }
        if( mGammaXYZ.empty() )
          throw std::runtime_error( "Gamma spatial not found for " + Description );
      }
      // For each temporal correlator - process with all corresponding spatial
      int iSeq{ 0 };
      std::string OutFileName{ outBase };
      OutFileName.append( 1, 'F' );
      {
        std::string sShort{ fna.GetBaseShort() };
        if( !sShort.empty() && std::toupper( sShort[0] ) == 'R' )
          sShort.erase( 0, 1 );
        OutFileName.append( sShort );
      }
      OutFileName.append( p.second.FileString( p.first ) );
      const std::size_t OutFileNameLen{ OutFileName.length() };
      for( const std::string &File : vFileGammaT )
      {
        try
        {
          Model mGT;
          mGT.Read( File, " reading temporal model ", &OpNames );
          // Most of the data can come from the temporal file
          MLU::Param::Key k( fna.MesonMom[iLight], mGT.EnergyPrefix );
          Column ELight = mGT.Column( k );
          k.Object[0] = fna.MesonMom[iHeavy];
          Column MHeavy = mGT.Column( k );
          k.Object[0] = fna.MesonMom[0];
          k.Object.emplace_back( fna.MesonMom[1] );
          k.Name = "MELgT";
          Column vT;
          try
          {
            vT = mGT.Column( k );
          } catch ( const std::runtime_error &e )
          {
            k.Name = "MEL"; // Old files just had a matrix element
            vT = mGT.Column( k );
          }
          // Get the renormalisation from the file if present
          Model mRenorm;
          Column ZVPrevSrc, ZVPrevSnk, ZVMixed;
          ZVMixed = 1;
          if( mEnsembleInfo.Nt() )
          {
            MLU::Param::Key k( "ZVmixed" );
            if( !mGT.Ensemble.empty() )
              k.Object.push_back( mGT.Ensemble );
            try
            {
              ZVMixed = mEnsembleInfo.Column( k );
            }
            catch(...)
            {
            }
          }
          bool bGotZVPrev{};
          std::string sPrevZVFilename;
          try
          {
            if( !mGT.FileList.empty() )
            {
              sPrevZVFilename = mGT.FileList[0];
              if( MLU::FileExists( sPrevZVFilename ) )
              {
                std::string GroupName{ MLU::sModel };
                mRenorm.Read( sPrevZVFilename, " Prior renormalisation ", &OpNames, &GroupName );
                MLU::Param::Key k( fna.Quark[0], "ZV" );
                ZVPrevSrc = mRenorm.Column( k );
                k.Object[0] = fna.Quark[1];
                ZVPrevSnk = mRenorm.Column( k );
                bGotZVPrev = true;
              }
            }
          }
          catch( const std::runtime_error &e ) {}
          catch( const ::H5::Exception &e ) {}
          if( !bGotZVPrev )
          {
            ZVPrevSrc = 1.;
            ZVPrevSnk = 1.;
            std::cout << " Previous ZV unavailable " << sPrevZVFilename << MLU::NewLine;
          }
          // NB: There's a single, dummy spatial model in the list for zero momentum
          for( Model &mSpatial : mGammaXYZ )
          {
            ++Total;
            try
            {
              int NumSamples{ mGT.NumSamples() };
              std::vector<std::string> vSourceFiles{ mGT.Name_.Filename };
              Column vXYZ;
              if( rf.p )
              {
                mGT.IsCompatible( mSpatial, &NumSamples,
                                  MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
                mGT.IsCompatible( *pmMLight, &NumSamples,
                                  MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
                vSourceFiles.emplace_back( mSpatial.Name_.Filename );
                vSourceFiles.emplace_back( pmMLight->Name_.Filename );
                MLU::Param::Key kSpatial( k.Object, "MELgXYZ" );
                try
                {
                  vXYZ = mSpatial.Column( kSpatial );
                } catch ( const std::runtime_error &e )
                {
                  kSpatial.Name = "MEL"; // Old files just had a matrix element
                  vXYZ = mSpatial.Column( kSpatial );
                }
              }
              OutFileName.resize( OutFileNameLen );
              if( iSeq )
              {
                OutFileName.append( 1, '.' );
                OutFileName.append( std::to_string( iSeq ) );
              }
              OutFileName.append( sOpNames );
              OutFileName = MLU::MakeFilename( OutFileName, MLU::sModel, mGT.Seed(), DEF_FMT );
              Write( OutFileName, mGT, std::move( vSourceFiles ), NumSamples, MHeavy, ELight,
                     vT, rf.p, ZVPrevSrc, ZVPrevSnk, ZVMixed, &MLight, &vXYZ );
              ++NumOK;
              ++iSeq;
            }
            catch( const std::exception &e )
            {
              std::cerr << "Error: " << e.what() << std::endl;
            }
          }
        }
        catch( const std::exception &e )
        {
          std::cerr << "Error: " << e.what() << std::endl;
        }
      }
    }
    catch( const std::exception &e )
    {
      std::cerr << "Error: " << e.what() << std::endl;
    }
    catch( ... )
    {
      std::cerr << "Error: Unknown exception" << std::endl;
    }
  }
}

int FFitConstMaker::Weight( const std::string &Quark )
{
  int Weight = 0;
  if( !Quark.empty() )
    switch( std::toupper( Quark[0] ) )
    {
      case 'L':
        Weight = 1;
        break;
      case 'S':
        Weight = 2;
        break;
      case 'H':
        Weight = 3;
        break;
    }
  return Weight;
}

// Incoming file should be a three-point function with (generally) a heavy source and light(er) sink
void FFitConstMaker::Make( std::string &FileName )
{
  // Parse Filename and check it contains what we need
  std::cout << "Processing " << FileName << MLU::NewLine;
  std::vector<std::string> OpNames;
  MLU::FileNameAtt fna{ FileName, &OpNames };
  if( OpNames.empty() )
    throw std::runtime_error( "Sink / source operator names mising" );
  if( !fna.GotDeltaT() )
    throw std::runtime_error( "DeltaT missing" );
  if( fna.Gamma.size() != 1 )
    throw std::runtime_error( FileName + " has " + std::to_string( fna.Gamma.size() ) + " currents" );
  if( fna.Gamma[0] != MLU::Gamma::Algebra::GammaT )
  {
    std::ostringstream ss;
    ss << "Ignoring current " << fna.Gamma[0] << " - should be " << MLU::Gamma::Algebra::GammaT << " only";
    throw std::runtime_error( ss.str().c_str() );
  }
  if( fna.p.empty() )
    throw std::runtime_error( FileName + " has no momenta" );
  if( fna.BaseShortParts.size() < 3 || fna.BaseShortParts[0].size() < 2
     || std::toupper( fna.BaseShortParts[0][0] ) != 'R' || fna.Extra.empty() )
    throw std::runtime_error( FileName + " is not a ratio file" );
  RatioNum = MLU::FromString<int>( fna.BaseShortParts[0].substr( 1 ) );
  qSnk = fna.BaseShortParts[1];
  qSrc = fna.BaseShortParts[2];
  const bool bSnkHeavier{ Weight( qSnk ) > Weight( qSrc ) };
  const std::string &MesonHeavy{ fna.Meson[bSnkHeavier ? 1 : 0] };
  const std::string &MesonLight{ fna.Meson[bSnkHeavier ? 0 : 1] };
  {
    const MLU::MomentumPair &np{ fna.GetFirstNonZeroMomentum() };
    p = np.second;
    pName = np.first;
  }
  // Get the fit type and fit ranges
  const std::string &FitTypeOriginal{ fna.Extra[0] };
  {
    std::string FitTypeBuffer{ FitTypeOriginal };
    FitType = MLU::ExtractToSeparator( FitTypeBuffer, MLU::Underscore );
    for( std::size_t pos = 0; ( pos = FitTypeBuffer.find( '_', pos ) ) != std::string::npos; )
      FitTypeBuffer[pos] = ' ';
    FitParts = MLU::ArrayFromString<int>( FitTypeBuffer );
    if( FitParts.empty() || FitParts.size() % 2 )
      throw std::runtime_error( FileName + " fit info " + ( FitParts.empty() ? "missing" : "bad" ) );
  }
  // Get filename prefix
  Prefix = std::to_string( RatioNum );
  Prefix.append( 1, '_' );
  Prefix.append( qSnk );
  Prefix.append( 1, '_' );
  Prefix.append( qSrc );
  // Get filename suffix
  std::string Suffix1;
  MLU::AppendDeltaT( Suffix1, fna.DeltaT[0] );
  Suffix1.append( p.FileString( pName ) );
  MLU::Append( Suffix1, fna.GetBaseShort( 3 ) );
  fna.AppendExtra( Suffix1 );
  std::string Suffix2{ MLU::Period };
  Suffix2.append( OpNames[fna.op[1]] );
  Suffix2.append( 1, '_' );
  Suffix2.append( OpNames[fna.op[0]] );
  Suffix2.append( 1, '.' );
  Suffix2.append( fna.Type );
  Suffix2.append( 1, '.' );
  Suffix2.append( fna.SeedString );
  Suffix2.append( 1, '.' );
  Suffix2.append( fna.Ext );
  Suffix = Suffix1;
  Suffix.append( 1, '.' );
  Suffix.append( FitTypeOriginal );
  Suffix.append( Suffix2 );
  // Debug - show filename
#if DEBUG
  std::cout << RatioNum << MLU::Space << qSnk << MLU::Space << qSrc << MLU::Space << fna.Gamma[0];
  for( std::size_t i = 3; i < fna.BaseShortParts.size(); ++i )
    std::cout << MLU::Space << fna.BaseShortParts[i];
  std::cout << " { " << FitType;
  for( std::size_t i = 0; i < FitParts.size(); i += 2 )
  {
    if( i )
      std::cout << MLU::Comma;
    std::cout << " [" << FitParts[i] << "," << FitParts[i+1] << "]";
  }
  std::cout << " }";
  for( std::size_t i = fna.Extra.size(); i-- > 1; )
    std::cout << " ." << fna.Extra[i];
  std::cout << " Suffix=" << Suffix << MLU::NewLine;
#endif
  // Check whether output exists
  std::string OutFileName{ outBase };
  OutFileName.append( 1, 'F' );
  OutFileName.append( Prefix );
  OutFileName.append( Suffix );
  //if( MLU::FileExists( OutFileName ) )
    //throw std::runtime_error( OutFileName + " exists" );
  // Open files
  Model mT, mXYZ;
  Column vXYZ;
  mT.Read( FileName, szReading );
  assert( Model::idxCentral == -1 && "Bug: Model::idxCentral != -1" );
  int NumSamples {};
  std::string FileNameXYZ;
  if( p )
  {
    std::string GlobName;
    {
      std::ostringstream ss;
      ss << "R" << Prefix << MLU::Underscore << MLU::Gamma::Algebra::Spatial << Suffix1 << ".*" << Suffix2;
      GlobName = ss.str();
    }
    std::vector<std::string> FileList{ MLU::glob( &GlobName, &GlobName + 1, i3Base.c_str() ) };
    if( FileList.empty() )
      throw std::runtime_error( GlobName + " not found" );
    if( FileList.size() != 1 )
      std::cout << MLU::Space << GlobName << " ambiguous" << MLU::NewLine;
    FileNameXYZ = std::move( FileList[0] );
    mXYZ.Read( FileNameXYZ, szReading );
    vXYZ = mXYZ.Sample<Scalar>::Column( 0 );
    mT.IsCompatible( mXYZ, &NumSamples, MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
  }
  // Load energies
  MP mpMHeavy( MesonHeavy, MLU::p0 );
  MP mpMLight( MesonLight, MLU::p0 );
  MP mpELight( MesonLight, p );
  Model &mMHeavy( EFit( mpMHeavy, szReading ) );
  mT.IsCompatible( mMHeavy, &NumSamples, MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
  Model &mMLight( EFit( mpMLight, szReading ) );
  mT.IsCompatible( mMLight, &NumSamples, MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
  Model &mELight( p ? EFit( mpELight, szReading ) : mMLight );
  if( p )
    mT.IsCompatible( mELight, &NumSamples, MLU::COMPAT_DISABLE_BASE | MLU::COMPAT_DISABLE_NT );
  const Column MHeavy{ EFit.GetColumn( mpMHeavy, mpMHeavy.PK() ) };
  const Column MLight{ EFit.GetColumn( mpMLight, mpMLight.PK() ) };
  const Column ELight{ EFit.GetColumn( mpELight, mpELight.PK() ) };

  // Make output
  std::vector<std::string> SourceFileNames{ FileName };
  if( p )
    SourceFileNames.emplace_back( FileNameXYZ );
  SourceFileNames.emplace_back( mMHeavy.Name_.Filename );
  SourceFileNames.emplace_back( mMLight.Name_.Filename );
  if( p )
    SourceFileNames.emplace_back( mELight.Name_.Filename );
  Column PrevOne{ 1. };
  Write( OutFileName, mT, std::move( SourceFileNames ), NumSamples, MHeavy, ELight, mT.Sample<Scalar>::Column(0), p, PrevOne, PrevOne, PrevOne, &MLight, &vXYZ );
}

Maker::Maker( const MLU::CommandLine &cl )
: MaxSamples{cl.SwitchValue<int>("n")},
  outBase{cl.SwitchValue<std::string>("o")},
  bSymmetrise{!cl.GotSwitch("nosym")},
  bOverlapAltNorm{ cl.GotSwitch( MLU::sOverlapAltNorm.c_str() ) },
  Cache2{ LoadFilePrefix },
  Cache3{ LoadFilePrefix }
{
  if( bOverlapAltNorm )
    std::cout << "WARNING: Use of --" << MLU::sOverlapAltNorm << " is deprecated.\n";
  Cache2.SetBase( cl.SwitchValue<std::string>("i2") );
  //Cache3.SetBase( cl.SwitchValue<std::string>("i3") );
  EFit.Read( cl.SwitchValue<std::string>( "efit" ), LoadFilePrefix );
  MLU::MakeAncestorDirs( outBase );
}

int main(int argc, const char *argv[])
{
  //if(!Debug()) return EXIT_FAILURE;
  static const char DefaultType[]{ "ZV" };
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i2", CL::SwitchType::Single, "" },
      {"i3", CL::SwitchType::Single, "" },
      {"n", CL::SwitchType::Single, "0"},
      {"o", CL::SwitchType::Single, "" },
      {"type", CL::SwitchType::Single, DefaultType },
      {"efit", CL::SwitchType::Single, ""},
      {MLU::sOverlapAltNorm.c_str(), CL::SwitchType::Flag, nullptr},
      {"dispersion", CL::SwitchType::Single, nullptr},
      {"nosym", CL::SwitchType::Flag, nullptr},
      {"r3a", CL::SwitchType::Flag, nullptr},
      {"debug-signals", CL::SwitchType::Flag, nullptr},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      if( cl.GotSwitch( "debug-signals" ) )
        MLU::Grid_debug_handler_init();
      dispType = cl.GotSwitch( "dispersion" ) ? cl.SwitchValue<MLU::DispersionType>( "dispersion" )
                                              : MLU::DispersionType::LatFreeScalar;
      // Read the list of fits I've chosen to use
      std::string TypeParams{ cl.SwitchValue<std::string>("type") };
      const std::string Type{ MLU::ExtractToSeparator( TypeParams ) };
      std::unique_ptr<Maker> m;
      if( MLU::EqualIgnoreCase( Type, "R" ) )
        m.reset( new RMaker( TypeParams, cl ) );
      else if( MLU::EqualIgnoreCase( Type, "ZV" ) )
        m.reset( new ZVMaker( TypeParams, cl ) );
      else if( MLU::EqualIgnoreCase( Type, "F" ) )
        m.reset( new FMaker( TypeParams, cl ) );
      else if( MLU::EqualIgnoreCase( Type, "FFC" ) )
        m.reset( new FFitConstMaker( TypeParams, cl ) );
      else
        throw std::runtime_error( "I don't know how to make type " + Type );
      bShowUsage = false;
      const std::string Prefix3pt{ cl.SwitchValue<std::string>("i3") };
      std::vector<std::string> FileList{ MLU::glob( cl.Args.begin(), cl.Args.end(),
                                                       Prefix3pt.c_str() ) };
      std::size_t Count{ 0 };
      std::size_t Done{ 0 };
      try
      {
        m->Run( Done, Count, FileList );
      }
      catch(const std::exception &e)
      {
        std::cerr << "Error: " << e.what() << std::endl;
      }
      catch( ... )
      {
        std::cerr << "Error: Unknown exception" << std::endl;
      }
      std::cout << Done << " of " << Count << MLU::Space << Type << " ratios created\n";
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
    "-n     Number of samples to fit, 0 = all available from bootstrap (default)\n"
    "-o     Output filename prefix\n"
    "--type Ratio type[,Options[,...]] (default " << DefaultType << ")\n"
    "       ZV  Make ZV (no Options)\n"
    "       R   Make R ratios. Options as per --efit, except Fit_list->ZV_List\n"
    "       F   Make form factors from matrix element fits.\n"
    "           Options\n"
    "           L[,spatial[,ZV[,FullyNPR]]]\n"
    "           L=spatial extent of lattice, (e.g. 24, 32, 64, ...)\n"
    "           spatial=scale vector matrix elements by momentum. DEPRECATED.\n"
    "           ZV=Filename containing ZV fits for mostly NPR\n"
    "              or Integer ZV\n"
    "              or empty = no renormalisation\n"
    "           FullyNPR=Filename for Fully non-perturbative ZV correction ZV_mixed\n"
    "       FFC Make form factors from fit to constant. Deprecated. Options\n"
    "           L i.e. spatial extent of lattice, (e.g. 24, 32, 64, ...)\n"
    "           'spatial' flag to divide spatial current gXYZ by integer momentum\n"
    "           (Should not be needed except for very old bootstrap files)\n"
    "--efit Fit_list[,options[,...]] List of fit files for energies,\n"
    "         each row has 'Meson_Momentum file' (one row per momentum).\n"
    "         First option can be 'c' to freeze values to central value,\n"
    "         then comma separated list of Key=Key_in_file pairs to rename keys.\n"
    "       If Fit_list='', comma separated list of Variable=Value pairs follows.\n"
    "         Any constants not mentioned are frozen to 1\n"
    "Flags:\n"
    "--" << MLU::sOverlapAltNorm << " Alternate normalisation for overlap factors. DEPRECATED\n"
    "--dispersion Which dispersion relation to use (default: LatFreeScalar)\n"
    "--nosym   Disable Out[t]=(Out[t]+Out[deltaT-t])/2 for symmetric ratios\n"
    "--r3a     Use alternate definition of R3 (i.e. R3a)\n"
    "--debug-signals Trap signals (code courtesy of Grid)\n"
    "--help    This message\n";
  }
  MLU::Grid_exit_handler_disable = true;
  return iReturn;
}
