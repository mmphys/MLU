/**

 Mike's lattice QCD utilities
 
 Source file: Utility.hpp
 
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

// Common utilities (no dependencies other than c++ stdlib)

#ifndef MLU_CorrFile_hpp
#define MLU_CorrFile_hpp

#include <MLU/Utility.hpp>

BEGIN_COMMON_NAMESPACE

// Correlator file. Could be either single correlator, or multiple gammas

template <typename T>
class CorrelatorFile
{
  // Data members
/*public:
  using Traits = SampleTraits<T>;
  using scalar_type = typename Traits::scalar_type;
  using value_type = typename Traits::value_type;
  static constexpr bool is_complex { Traits::is_complex };*/
private:
  int NumSnk_ = 0;
  int NumSrc_ = 0;
  int Nt_ = 0;
  std::unique_ptr<T[]> m_pData;
  std::vector<Gamma::Algebra> AlgSnk_;
  std::vector<Gamma::Algebra> AlgSrc_;
public:
  FileNameAtt Name_;
  bool bHasTimeslice = false;
  int  Timeslice_ = 0;
  // Member functions
private:
  inline void RangeCheck( int Sample ) const
  {
    if( Sample < 0 || Sample >= NumSnk_ * NumSrc_ )
      throw std::out_of_range( "Sample " + std::to_string( Sample ) );
  }
  inline int SinkIndex( Gamma::Algebra g ) const
  {
    int idx;
    for( idx = 0; idx < AlgSnk_.size() && AlgSnk_[idx] != g; idx++ )
      ;
    return idx;
  }
  inline int SourceIndex( Gamma::Algebra g ) const
  {
    int idx;
    for( idx = 0; idx < AlgSrc_.size() && AlgSrc_[idx] != g; idx++ )
      ;
    return idx;
  }
public:
  // Constructors (copy operations missing for now - add them if they become needed)
  CorrelatorFile() {}
  CorrelatorFile( CorrelatorFile && ) = default; // Move constructor
  CorrelatorFile( const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk,
                  std::vector<Gamma::Algebra> &AlgSrc, const int * pTimeslice = nullptr,
                  const char * PrintPrefix = nullptr, std::string *pGroupName = nullptr )
  {
    Read( FileName, AlgSnk, AlgSrc, pTimeslice, PrintPrefix, pGroupName );
  }
  // Operators
  inline CorrelatorFile& operator=( CorrelatorFile && r ) = default; // Move assignment
  inline       T * operator[]( int Sample )
  {
    RangeCheck( Sample );
    return & m_pData[static_cast<std::size_t>( Sample ) * Nt_];
  }
  inline const T * operator[]( int Sample ) const
  {
    RangeCheck( Sample );
    return & m_pData[static_cast<std::size_t>( Sample ) * Nt_];
  }
  inline       T * operator()( Gamma::Algebra gSink, Gamma::Algebra gSource )
  { return (*this)[ SinkIndex( gSink ) * NumSrc_ + SourceIndex( gSource ) ]; }
  inline const T * operator()( Gamma::Algebra gSink, Gamma::Algebra gSource ) const
  { return (*this)[ SinkIndex( gSink ) * NumSrc_ + SourceIndex( gSource ) ]; }
  inline int NumSnk() const { return NumSnk_; }
  inline int NumSrc() const { return NumSrc_; }
  inline int NumOps() const { return NumSnk() * NumSrc(); }
  inline int Nt() const { return Nt_; }
  inline int Timeslice() const { return bHasTimeslice ? Timeslice_ : 0; }
  inline const std::vector<Gamma::Algebra> &AlgSnk() const { return AlgSnk_; }
  inline const std::vector<Gamma::Algebra> &AlgSrc() const { return AlgSrc_; }
  inline bool IsFinite() const
  { return Common::IsFinite( reinterpret_cast<typename SampleTraits<T>::scalar_type *>( m_pData.get() ),
      static_cast<size_t>( NumSnk_ * NumSrc_ * ( SampleTraits<T>::is_complex ? 2 : 1 ) ) * Nt_ ); }
  void resize( int NumSnk, int NumSrc, int Nt );
  /// Swap function so that this type is sortable
  void swap( CorrelatorFile &o );
  void Read (const std::string &FileName, std::vector<Gamma::Algebra> &AlgSnk, std::vector<Gamma::Algebra> &AlgSrc,
             const int * pTimeslice = nullptr, const char * PrintPrefix = nullptr,
             std::string *pGroupName = nullptr, const std::vector<NegateStar> *pAlgSnkNeg = nullptr,
             const char * pDSName = nullptr );
  //void Write( const std::string &FileName, const char * pszGroupName = nullptr );
  void WriteSummary(const std::string &Prefix, const std::vector<Gamma::Algebra> &AlgSnk,
                    const std::vector<Gamma::Algebra> &AlgSrc);
};

using CorrelatorFileC = CorrelatorFile<std::complex<double>>;
using CorrelatorFileD = CorrelatorFile<double>;

template <typename T>
void swap( CorrelatorFile<T> &l, CorrelatorFile<T> &r );

END_COMMON_NAMESPACE
#endif // MLU_CorrFile
