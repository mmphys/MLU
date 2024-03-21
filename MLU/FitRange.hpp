/**
 
 Mike's lattice QCD utilities: Fit ranges
 
 Source file: FitRange.hpp
 
 Copyright (C) 2022

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

#ifndef FitRange_hpp
#define FitRange_hpp

#include <MLU/MLUFirst.hpp>

#include <array>
#include <algorithm>
#include <cassert>
#include <cctype>
#include <complex>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

BEGIN_MLU_NAMESPACE

struct FitRangesIterator;
struct FitRangesIteratorElement;

/**
 Abstract base clase for a single range of fit times to scan
 
 Can have dependencies on zero or more other FitRange
 
 Specialisations:

 `FitRangeAbsolute`: scans a range of start stop times (optional thinning)

 `FitRangeRelative`: scans start stop times relative to anither FitRange (optional thinning)
 */
struct FitRange
{
  int MinDP; // Minimum number of data points in this range
  static FitRange * Deserialise( const std::string &String, std::size_t MyIndex, int MinDP );
  static std::vector<int> GetColonList( std::istream &is,
                                        std::size_t MaxLen=std::numeric_limits<std::size_t>::max() );
  virtual ~FitRange() {}
  virtual bool Validate( int Nt = std::numeric_limits<int>::max() ) const = 0;
  virtual FitRangesIteratorElement * GetStart( const FitRangesIterator &Parent ) const = 0;
  virtual FitRangesIteratorElement * GetEnd( const FitRangesIterator &Parent ) const = 0;
  virtual const std::vector<std::size_t> &GetDependencies() const { return DependsOn; }
  virtual void Print( std::ostream &os ) const = 0;
protected:
  std::vector<std::size_t> DependsOn;
};

inline std::ostream & operator<<( std::ostream &os, const FitRange &fr )
{
  fr.Print( os );
  return os;
}

/**
 A collection of one or more FitRange

 Can be iterated over, i.e. provides `begin()` and `end()`.

 NB: the fit ranges are evaluated (and therefore sorted in memory) so that each
 FitRange appears before any others it depends on.
 
 `operator[]` is redefined to return each FitRange & in OriginalIndex order
 
 `operator()` can be used to return each FitRange & in MemoryIndex order
 */
struct FitRanges : public std::vector<std::unique_ptr<FitRange>>
{
  using Base = std::vector<std::unique_ptr<FitRange>>;
  using Base::Base; // Import base constructors
  void Deserialise( const std::vector<std::string> &vString, int MinDP );
  FitRanges( std::vector<std::string> vString, int MinDP ) { Deserialise( vString, MinDP ); }
  FitRangesIterator begin() const;
  FitRangesIterator end() const;
  /// Default is to access FitRanges using the original index order
  inline std::size_t MemIndex( std::size_t OriginalIndex ) const { return vIndex[OriginalIndex]; }
  inline       FitRange &operator[]( std::size_t OriginalIndex )
  { return * Base::operator[]( MemIndex( OriginalIndex ) ).get(); }
  inline const FitRange &operator[]( std::size_t OriginalIndex ) const
  { return * Base::operator[]( MemIndex( OriginalIndex ) ).get(); }
  /// FitRanges can also to be accessed using the in-memory (dependency) order
  inline       FitRange &operator()( std::size_t MemoryIndex )
  { return * Base::operator[]( MemoryIndex ).get(); }
  inline const FitRange &operator()( std::size_t MemoryIndex ) const
  { return * Base::operator[]( MemoryIndex ).get(); }
protected:
  /// This is a map from the original order to the order in memory. I.e. the nth entry of this array contains the in-memory order
  std::vector<std::size_t> vIndex;
};

std::ostream & operator<<( std::ostream &os, const FitRanges &fr );

/**
 Abstract base class for each element within a FitRangesIterator
 
 Can be queried for the fit range being scanned.
 */
struct FitRangesIteratorElement
{
  const FitRange &fitRange;
  const FitRangesIterator &Parent;
  FitRangesIteratorElement( const FitRange &fitRange_, const FitRangesIterator &parent )
  : fitRange{fitRange_}, Parent{parent} {}
  virtual ~FitRangesIteratorElement() {}
  virtual FitRangesIteratorElement * Clone() const = 0;
  /// Number of timeslices
  virtual int Extent() const = 0;
  virtual int TI() const = 0;
  virtual int TF() const = 0;
  virtual void SetStart() = 0;
  virtual bool PastEnd() const = 0;
  virtual bool Increment( bool bWrap ) = 0;
  std::string AbbrevString( const std::string &Sep ) const;
  virtual void Print( std::ostream &os ) const { os << AbbrevString( "-" ); }
  const std::vector<int> &GetFitTimes() const { return FitTimes; }
protected:
  std::vector<int> FitTimes;
};

inline std::ostream & operator<<( std::ostream &os, const FitRangesIteratorElement &elem )
{
  elem.Print( os );
  return os;
}

/**
 Iterator over `FitRanges`. A collection of one or more `FitRangesIteratorElement`

 Typical iteration loop:
 
 `for(auto it = fr.begin(); !it.PastEnd(); ++it)`

 `operator[]` is redefined to return each Element & in OriginalIndex order
 
 `operator()` can be used to return each Element & in MemoryIndex order
 */
struct FitRangesIterator : public std::vector<std::unique_ptr<FitRangesIteratorElement>>
{
  const FitRanges &Ranges; // This is the FitRange I'm iterating over
  using Element = FitRangesIteratorElement;
  using Base = std::vector<std::unique_ptr<Element>>;
  using Base::Base;
  FitRangesIterator( const FitRanges &ranges_, bool bEnd );
  FitRangesIterator( const FitRangesIterator &it );
  bool PastEnd() const;
  FitRangesIterator &operator++(); // Prefix increment
  FitRangesIterator operator++(int) { FitRangesIterator old(*this); this->operator++(); return old; }
  std::string AbbrevString( const std::string &Sep1, const std::string &Sep2 ) const;
  std::string AbbrevString( const std::string &Sep = "_" ) const { return AbbrevString( Sep, Sep ); }
  /// Default is to access FitRangesIterator using the original index order
  inline       Element &operator[]( std::size_t OriginalIndex )
  { return * Base::operator[]( Ranges.MemIndex( OriginalIndex ) ).get(); }
  inline const Element &operator[]( std::size_t OriginalIndex ) const
  { return * Base::operator[]( Ranges.MemIndex( OriginalIndex ) ).get(); }
  /// FitRangesIterator can also to be accessed using the in-memory (dependency) order
  inline       Element &operator()( std::size_t MemoryIndex )
  { return * Base::operator[]( MemoryIndex ).get(); }
  inline const Element &operator()( std::size_t MemoryIndex ) const
  { return * Base::operator[]( MemoryIndex ).get(); }
protected:
  inline bool GotMinDP() const;
};

std::ostream & operator<<( std::ostream &os, const FitRangesIterator &it );

END_MLU_NAMESPACE
#endif // FitRange_hpp
