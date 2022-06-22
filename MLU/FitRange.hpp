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
#define FitRange_hpp namespace Common {
#define FitRange_hpp_end };

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

FitRange_hpp

/*
 
 FitRange:  A single range of initial and final fit times to scan
            Can have dependencies on zero or more other FitRange
            Abstract base clase (see below for specialisations)

 FitRanges: A collection of one or more FitRange
            Can be iterated over, i.e. provides begin() and end()
            NB: the fit ranges are evaluated (and therefore sorted in memory) so that each
            FitRange appears before any others it depends on, with vIndex providing a map
            vIndex[OriginalIndex]=MemoryIndex
            operator[] is redefined to return const FitRange & in OriginalIndex order

 FitTime:   A single start-stop pair

 FitRangesIterator:
            A collection of one or more FitTime, iterating over FitRanges
            operator[] is redefined to return const FitTime & in OriginalIndex order

 FitRange specialisations
 
 FitRangeAbsolute:  ti:tf[dti[:dtf]]
                    All combinations of Start [ti,ti+dti) ... [tf,tf+dtf)

 FitRangeRelative:  Rn:ti:tf[dti[:dtf]]
                    All combinations of Range n ti + [ti,ti+dti) ... Range n tf +  [tf,tf+dtf)

 */

// A single pair of start-stop times
struct FitTime
{
  int ti;
  int tf;
};

struct FitRange;
struct FitRangesIterator;

struct FitRanges : public std::vector<std::unique_ptr<FitRange>>
{
  std::vector<std::size_t> vIndex;
  using Base = std::vector<std::unique_ptr<FitRange>>;
  using Base::Base; // Import base constructors
  void Deserialise( const std::vector<std::string> &vString );
  FitRanges( std::vector<std::string> vString ) { Deserialise( vString ); }
  FitRangesIterator begin() const;
  FitRangesIterator end() const;
  const FitRange & operator[]( std::size_t Index ) const { return * Base::operator[]( vIndex[Index] ).get(); }
};

struct FitRange
{
  using vDepend = std::vector<std::size_t>;
  using vFitTime = std::vector<FitTime>;
  virtual bool Validate( int Nt = std::numeric_limits<int>::max() ) const = 0;
  virtual void GetStart( FitTime &ft, const FitRangesIterator &it ) const = 0;
  virtual void GetEnd( FitTime &ft, const FitRangesIterator &it ) const = 0;
  virtual bool PastEnd( const FitTime &ft, const FitRangesIterator &it ) const = 0;
  virtual bool Increment( FitTime &ft, const FitRangesIterator &it, bool bWrap ) const = 0;
  virtual void Print( std::ostream &os ) const = 0;
  virtual const vDepend &GetDependencies() const { return DependsOn; };
  virtual ~FitRange() {};
  static FitRange * Deserialise( const std::string &String );
protected:
  vDepend DependsOn;
};

std::ostream & operator<<( std::ostream &os, const FitRange &fr );
//std::istream & operator>>( std::istream &is, FitRange &fr );

struct FitRangesIterator : public FitRange::vFitTime
{
  const FitRanges &Ranges; // This is the FitRange I'm iterating over
  using Base = FitRange::vFitTime;
  using Base::Base;
  FitRangesIterator( const FitRanges &ranges_, bool bEnd );
  FitRangesIterator( const FitRangesIterator &it ) : Base( it ), Ranges{ it.Ranges }, RangeMemOrder{Ranges}{}
  bool PastEnd() const { return Ranges.back()->PastEnd( this->back(), *this ); }
  FitRangesIterator &operator++(); // Prefix increment
  FitRangesIterator operator++(int) { FitRangesIterator old(*this); this->operator++(); return old; }
  std::string to_string( const std::string &Sep1, const std::string &Sep2 ) const;
  std::string to_string( const std::string &Sep ) const { return to_string( Sep, Sep ); }
  std::string to_string() const;
  const FitTime &operator[]( std::size_t Index ) const { return Base::operator[]( Ranges.vIndex[Index] ); }
protected:
  const FitRanges::Base &RangeMemOrder;
};

struct FitRangeAbsolute : FitRange
{
protected:
  int ti;
  int tf;
  int dti;
  int dtf;
public:
  FitRangeAbsolute( int ti_, int tf_, int dti_, int dtf_ ) : ti{ti_}, tf{tf_}, dti{dti_}, dtf{dtf_} {}
  FitRangeAbsolute( int ti_, int tf_, int dt ) : FitRangeAbsolute( ti_, tf_, dt, dt ) {}
  FitRangeAbsolute( int ti_, int tf_ ) : FitRangeAbsolute( ti_, tf_, 1, 1 ) {}
  FitRangeAbsolute() : FitRangeAbsolute( 0, 0, 1, 1 ) {}
  bool Validate( int Nt = std::numeric_limits<int>::max() ) const override;
  void GetStart( FitTime &ft, const FitRangesIterator &it ) const override;
  void GetEnd( FitTime &ft, const FitRangesIterator &it ) const override;
  bool PastEnd( const FitTime &ft, const FitRangesIterator &it ) const override;
  bool Increment( FitTime &ft, const FitRangesIterator &it, bool bWrap ) const override;
  void Print( std::ostream &os ) const override;
  static FitRange * Deserialise( std::istringstream &is );
};

struct FitRangeRelative : FitRange
{
protected:
  int ti;
  int tf;
  int dti;
  int dtf;
public:
  FitRangeRelative( int Depends, int ti_, int tf_, int dti_, int dtf_ )
  : ti{ti_}, tf{tf_}, dti{dti_}, dtf{dtf_} { DependsOn.push_back( Depends ); }
  FitRangeRelative( int Depends, int ti_, int tf_, int dt ) : FitRangeRelative( Depends, ti_, tf_, dt, dt ) {}
  FitRangeRelative( int Depends, int ti_, int tf_ ) : FitRangeRelative( Depends, ti_, tf_, 1, 1 ) {}
  FitRangeRelative() : FitRangeRelative( 0, 0, 0, 1, 1 ) {}
  bool Validate( int Nt = std::numeric_limits<int>::max() ) const override;
  void GetStart( FitTime &ft, const FitRangesIterator &it ) const override;
  void GetEnd( FitTime &ft, const FitRangesIterator &it ) const override;
  bool PastEnd( const FitTime &ft, const FitRangesIterator &it ) const override;
  bool Increment( FitTime &ft, const FitRangesIterator &it, bool bWrap ) const override;
  void Print( std::ostream &os ) const override;
  static FitRange * Deserialise( std::istringstream &is );
};

FitRange_hpp_end
#endif // FitRange_hpp
