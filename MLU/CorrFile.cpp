/**

 Correlator file as saved by Hadrons (Grid)
 
 Source file: CorrFile.cpp
 
 Copyright (C) 2019 - 2023
 
 Author: Michael Marshall <Michael.Marshall@ed.ac.uk>
 
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
**/

#include "CorrFile.hpp"

BEGIN_MLU_NAMESPACE

template <typename T>
bool CorrelatorFile<T>::IsFinite() const
{
  return MLU::IsFinite( m_pData.get(), static_cast<size_t>( NumSnk_ ) * NumSrc_ * Nt_ );
}

template <typename T>
void CorrelatorFile<T>::resize( int NumSnk, int NumSrc, int Nt )
{
  const std::size_t OldMemSize{ static_cast<std::size_t>( NumSnk_ ) * NumSrc_ * Nt_ };
  const std::size_t NewMemSize{ static_cast<std::size_t>( NumSnk  ) * NumSrc  * Nt  };
  if( NumSnk_ != NumSnk || NumSrc_ != NumSrc )
  {
    AlgSnk_.clear();
    AlgSrc_.clear();
    NumSnk_ = NumSnk;
    NumSrc_ = NumSrc;
  }
  Nt_ = Nt;
  if( OldMemSize != NewMemSize )
  {
    if( NewMemSize == 0 )
      m_pData.reset( nullptr );
    else
      m_pData.reset( new T[ NewMemSize ] );
  }
}

template <typename T>
void CorrelatorFile<T>::swap( CorrelatorFile &o )
{
  using std::swap;
  swap( NumSnk_, o.NumSnk_ );
  swap( NumSrc_, o.NumSrc_ );
  swap( Nt_, o.Nt_ );
  swap( m_pData, o.m_pData );
  swap( AlgSnk_, o.AlgSnk_ );
  swap( AlgSrc_, o.AlgSrc_ );
  Name_.swap( o.Name_ );
  swap( bHasTimeslice, o.bHasTimeslice );
  swap( Timeslice_, o.Timeslice_ );
}

// Read from file. If GroupName empty, read from first group and return name in GroupName
template <typename T>
void CorrelatorFile<T>::Read(const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk,
                             std::vector<Gamma::Algebra> &AlgSrc, const int * pTimeslice,
                             const char * PrintPrefix, std::string *pGroupName,
                             const std::vector<NegateStar> *pAlgSnkNeg, const char * pDSName )
{
  const bool bSameAlgebras{&AlgSnk == &AlgSrc};
  if( pAlgSnkNeg )
  {
    assert( pAlgSnkNeg->size() == AlgSnk.size() );
    assert( pAlgSnkNeg->size() && "Algebras unknown. Negating makes no sense" );
    assert( !bSameAlgebras && "Work out whether two negatives make a positive, then implement" );
  }
  // Parse the name. This is PreBootstrap, i.e. no type and not my naming convention
  Name_.Parse( FileName, nullptr, nullptr, nullptr, true );
  if( !Name_.bSeedNum )
    throw std::runtime_error( "Configuration number missing from " + FileName );
  // Strip out timeslice info if present
  bHasTimeslice = pTimeslice;
  if( bHasTimeslice )
    Timeslice_ = * pTimeslice;
  else
  {
    std::string sCopy{ Name_.Base };
    ExtractTimeslice( sCopy, bHasTimeslice, Timeslice_ );
  }
  // Now load the file
  ::H5::H5File f;
  ::H5::Group  g;
  std::string localGroupName;
  if( !pGroupName )
    pGroupName = &localGroupName;
  H5::OpenFileGroup( f, g, FileName, PrintPrefix, pGroupName );
  bool bOK = false;
  H5E_auto2_t h5at;
  void      * f5at_p;
  ::H5::Exception::getAutoPrint(h5at, &f5at_p);
  ::H5::Exception::dontPrint();
  try // to load a single correlator
  {
    if( ( AlgSnk.empty() || ( AlgSnk.size() == 1 && AlgSnk[0] == Gamma::Algebra::Unknown ) )
      && ( AlgSrc.empty() || ( AlgSrc.size() == 1 && AlgSrc[0] == Gamma::Algebra::Unknown ) ) )
    {
      ::H5::DataSet ds = g.openDataSet( ( pDSName && *pDSName ) ? pDSName : "correlator" );
      ::H5::DataSpace dsp = ds.getSpace();
      if( dsp.getSimpleExtentNdims() == 1 )
      {
        hsize_t Dim[1];
        dsp.getSimpleExtentDims( Dim );
        if( Dim[0] <= std::numeric_limits<int>::max() )
        {
          resize( 1, 1, static_cast<int>( Dim[0] ) );
          AlgSnk_.resize( 1 );
          AlgSnk_[0] = Gamma::Algebra::Unknown;
          AlgSrc_.resize( 1 );
          AlgSrc_[0] = Gamma::Algebra::Unknown;
          if( AlgSnk.empty() )
            AlgSnk.push_back( Gamma::Algebra::Unknown );
          if( !bSameAlgebras && AlgSrc.empty() )
            AlgSrc.push_back( Gamma::Algebra::Unknown );
          T * const pData{ (*this)[0] };
          ds.read( pData, H5::Equiv<T>::Type );
          // Negate this correlator if requested
          if( pAlgSnkNeg )
            for( int i = 0; i < Nt_; i++ )
              ApplyNegateStar( pData[i], (*pAlgSnkNeg)[0] );
          bOK = true;
        }
      }
    }
  }
  catch(const ::H5::Exception &)
  {
    bOK = false;
    ::H5::Exception::clearErrorStack();
  }
  if( !bOK )
  {
    try // to load from array of correlators indexed by gamma matrix
    {
      unsigned short NumVec;
      ::H5::Attribute a;
      a = g.openAttribute("_Grid_vector_size");
      a.read( ::H5::PredType::NATIVE_USHORT, &NumVec );
      a.close();
      // Must be a perfect square and have at least as many as entries as requested
      const unsigned short NumFileOps{static_cast<unsigned short>( std::sqrt( NumVec ) + 0.5 )};
      if( NumFileOps * NumFileOps == NumVec && NumFileOps >= AlgSnk.size() && NumFileOps >= AlgSrc.size() )
      {
        bOK = true;
        std::vector<int> count;
        for( unsigned short i = 0; bOK && i < NumVec; i++ ) // Loop through all vectors in file
        {
          bOK = false;
          ::H5::Group gi = g.openGroup( *pGroupName + "_" + std::to_string( i ) );
          ::H5::DataSet ds = gi.openDataSet( "corr" );
          ::H5::DataSpace dsp = ds.getSpace();
          if( dsp.getSimpleExtentNdims() == 1 )
          {
            hsize_t Dim[1];
            dsp.getSimpleExtentDims( Dim );
            if( Dim[0] <= std::numeric_limits<int>::max() )
            {
              const int ThisNt{ static_cast<int>( Dim[0] ) };
              if( i == 0 )
              {
                // First correlator - resize and save operators if known
                resize(AlgSnk.size() ? static_cast<int>( AlgSnk.size() ) : NumFileOps,
                       AlgSrc.size() ? static_cast<int>( AlgSrc.size() ) : NumFileOps,
                       ThisNt);
                count.resize( NumSnk_ * NumSrc_, 0 ); // I want to check each operator combination appears once only
                if( AlgSnk.size() )
                {
                  AlgSnk_.resize( AlgSnk.size() );
                  std::copy( AlgSnk.cbegin(), AlgSnk.cend(), AlgSnk_.begin() );
                }
                else
                {
                  AlgSnk_.clear();
                  AlgSnk_.reserve( NumFileOps );
                }
                if( AlgSrc.size() )
                {
                  AlgSrc_.resize( AlgSrc.size() );
                  std::copy( AlgSrc.cbegin(), AlgSrc.cend(), AlgSrc_.begin() );
                }
                else
                {
                  AlgSrc_.clear();
                  AlgSrc_.reserve( NumFileOps );
                }
              }
              else if( ThisNt != Nt_ )
              {
                break;
              }
              // Read the gamma algebra strings and make sure they are valid
              const Gamma::Algebra gSnk{ H5::ReadGammaAttribute( gi, "gamma_snk" ) };
              int idxSnk;
              for( idxSnk = 0; idxSnk < AlgSnk_.size() && AlgSnk_[idxSnk] != gSnk; idxSnk++ )
                ;
              if( idxSnk == AlgSnk_.size() && AlgSnk_.size() < NumSnk_ )
                AlgSnk_.push_back( gSnk );
              bOK = true; // We can safely ignore gamma structures we're not interested in
              if( idxSnk < AlgSnk_.size() )
              {
                const Gamma::Algebra gSrc{ H5::ReadGammaAttribute( gi, "gamma_src" ) };
                int idxSrc;
                for( idxSrc = 0; idxSrc < AlgSrc_.size() && AlgSrc_[idxSrc] != gSrc; idxSrc++ )
                  ;
                if( idxSrc == AlgSrc_.size() && AlgSrc_.size() < NumSrc_ )
                  AlgSrc_.push_back( gSrc );
                if( idxSrc < AlgSrc_.size() )
                {
                  const int idx{ idxSnk * NumSrc_ + idxSrc };
                  T * const pData{ (*this)[idx] };
                  ds.read( pData, H5::Equiv<T>::Type );
                  count[idx]++;
                  if( pAlgSnkNeg )
                    for( int i = 0; i < Nt_; i++ )
                      ApplyNegateStar( pData[i], (*pAlgSnkNeg)[idxSnk] );
                }
              }
            }
          }
        }
        // Make sure that everything we wanted was loaded once and only once
        for( int i = 0; bOK && i < NumSrc_ * NumSnk_; i++ )
          if( count[i] != 1 )
            bOK = false;
      }
    }
    catch(const ::H5::Exception &)
    {
      bOK = false;
      ::H5::Exception::clearErrorStack();
    }
  }
  ::H5::Exception::setAutoPrint(h5at, f5at_p);
  if( !bOK )
    throw std::runtime_error( "Unable to read sample from " + FileName );
  if( !IsFinite() )
    throw std::runtime_error( "Values read are not all finite" );
  // If I'm discovering which operators are in the file, copy them back to caller
  // Bear in mind that caller may have passed in the same array for each gamma algebra
  bool bCopyBackSrc{ AlgSrc.empty() };
  bool bCopyBackSnk{ !bSameAlgebras && AlgSnk.empty() };
  if( bCopyBackSrc )
  {
    AlgSrc.resize( AlgSrc_.size() );
    std::copy( AlgSrc_.cbegin(), AlgSrc_.cend(), AlgSrc.begin() );
  }
  if( bCopyBackSnk )
  {
    AlgSnk.resize( AlgSnk_.size() );
    std::copy( AlgSnk_.cbegin(), AlgSnk_.cend(), AlgSnk.begin() );
  }
}

template <typename T>
void CorrelatorFile<T>::WriteSummary( const std::string &Prefix, const std::vector<Gamma::Algebra> &AlgSnk, const std::vector<Gamma::Algebra> &AlgSrc )
{
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  const int nt{ Nt() };
  const std::vector<Gamma::Algebra> &MySnk{ AlgSnk.size() ? AlgSnk : AlgSnk_ };
  const std::vector<Gamma::Algebra> &MySrc{ AlgSrc.size() ? AlgSrc : AlgSrc_ };
  const int NumSnk{ static_cast<int>( MySnk.size() ) };
  const int NumSrc{ static_cast<int>( MySrc.size() ) };
  std::string sOutFileName{ Prefix };
  sOutFileName.append( Name_.Base );
  std::size_t Len{ sOutFileName.length() };
  std::string sSuffix( 1, '.' );
  sSuffix.append( Name_.SeedString );
  sSuffix.append( 1, '.' );
  sSuffix.append( TEXT_EXT );
  for( int Snk = 0; Snk < NumSnk; Snk++ )
  {
    static const char pszSep[] = "_";
    sOutFileName.resize( Len );
    sOutFileName.append( MLU::Gamma::NameShort( AlgSnk[Snk], pszSep ) );
    std::size_t Len2{ sOutFileName.length() };
    for( int Src = 0; Src < NumSrc; Src++ )
    {
      sOutFileName.resize( Len2 );
      sOutFileName.append( MLU::Gamma::NameShort( AlgSrc[Src], pszSep ) );
      sOutFileName.append( sSuffix );
      SummaryHelper( sOutFileName, (*this)( AlgSnk[Snk], AlgSrc[Src] ), nt );
    }
  }
}

template class CorrelatorFile<double>;
template class CorrelatorFile<float>;
template class CorrelatorFile<std::complex<double>>;
template class CorrelatorFile<std::complex<float>>;

template <typename T>
void swap( CorrelatorFile<T> &l, CorrelatorFile<T> &r )
{
  l.swap( r );
}

template void swap( CorrelatorFileC &l, CorrelatorFileC &r );
template void swap( CorrelatorFileD &l, CorrelatorFileD &r );
template void swap( CorrelatorFile<float> &l, CorrelatorFile<float> &r );
template void swap( CorrelatorFile<std::complex<float>> &l, CorrelatorFile<std::complex<float>> &r );

END_MLU_NAMESPACE
