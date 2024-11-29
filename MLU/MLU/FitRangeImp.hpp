/**
 
 Mike's lattice QCD utilities: Fit ranges implementation. Only to be included within MLU itself - not to be exposed to users of MLU
 
 Source file: FitRangeImp.hpp
 
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

#ifndef FitRangeImp_hpp
#define FitRangeImp_hpp

#include "FitRange.hpp"

BEGIN_MLU_NAMESPACE

/*
 FitRange specialisations
 
 FitRangeAbsolute:  ti:tf[:dti[:dtf]]
                    All combinations of Start [ti,ti+dti) ... [tf,tf+dtf)

 FitRangeRelative:  Rn:ti:tf[:dti[:dtf]]
                    All combinations of Range n ti + [ti,ti+dti) ... Range n tf +  [tf,tf+dtf)
 */

struct FitRangeBlock : FitRange
{
  FitRangeBlock( int ti_, int tf_, int dti_, int dtf_ ) : ti{ti_}, tf{tf_}, dti{dti_}, dtf{dtf_} {}
  FitRangeBlock( int ti_, int tf_, int dt ) : FitRangeBlock( ti_, tf_, dt, dt ) {}
  FitRangeBlock( int ti_, int tf_ ) : FitRangeBlock( ti_, tf_, 1, 1 ) {}
  FitRangeBlock() : FitRangeBlock( 0, 0, 1, 1 ) {}
  virtual ~FitRangeBlock() {}
  /**
   Specifies how to thin a range.
   
   Sequence of `Delta`-`Num` pairs, where Delta is the distance between timeslices and Num is the number of timeslices
   Final `Num` is optional and means all remaining timeslices.
   
   `Delta=0` means a gap of `Num` timeslices (not valid as first entry)
   */
  void GetThinning( std::istream &is );
  void ShowThinning( std::ostream &os ) const;
protected:
  int ti;
  int tf;
  int dti;
  int dtf;
  std::vector<int> Thinning;
  friend class FitRangesIteratorBlock;
  friend class FitRangesIteratorAbsolute;
  friend class FitRangesIteratorRelative;
};

// A single pair of start-stop times
struct FitRangesIteratorBlock : FitRangesIteratorElement
{
  int ti;
  int tf;
  using Base = FitRangesIteratorElement;
  using Base::Base;
  FitRangesIteratorBlock( const FitRange &fitRange_, const FitRangesIterator &parent, int ti_,int tf_)
  : FitRangesIteratorElement( fitRange_, parent ), ti{ti_}, tf{tf_} {}
  virtual ~FitRangesIteratorBlock() {}
  int Extent() const override { return tf - ti + 1; } // Number of timeslices
  int TI() const override { return ti; }
  int TF() const override { return tf; }
protected:
  /// Returns true if FitTimes same as last time and therefore we need to keep incrementing
  bool SetFitTimes();
  void Print( std::ostream &os ) const override;
};

struct FitRangeAbsolute : FitRangeBlock
{
  using FitRangeBlock::FitRangeBlock;
  virtual ~FitRangeAbsolute() {}
  bool Validate( int Nt = std::numeric_limits<int>::max() ) const override;
  FitRangesIterator::Element * GetStart( const FitRangesIterator &Parent ) const override;
  FitRangesIterator::Element * GetEnd( const FitRangesIterator &Parent ) const override;
  void Print( std::ostream &os ) const override;
  static FitRange * Deserialise( std::istringstream &is, std::size_t MyIndex );
};

struct FitRangesIteratorAbsolute : FitRangesIteratorBlock
{
  using FitRangesIteratorBlock::FitRangesIteratorBlock;
  virtual ~FitRangesIteratorAbsolute() {}
  FitRangesIteratorElement * Clone() const override
  { return new FitRangesIteratorAbsolute( fitRange, Parent, ti, tf ); }
  void SetStart() override;
  bool PastEnd() const override;
  bool Increment( bool bWrap ) override;
};

struct FitRangeRelative : FitRangeBlock
{
  using FitRangeBlock::FitRangeBlock;
  virtual ~FitRangeRelative() {}
  bool Validate( int Nt = std::numeric_limits<int>::max() ) const override;
  FitRangesIterator::Element * GetStart( const FitRangesIterator &Parent ) const override;
  FitRangesIterator::Element * GetEnd( const FitRangesIterator &Parent ) const override;
  void Print( std::ostream &os ) const override;
  static FitRange * Deserialise( std::istringstream &is, std::size_t MyIndex );
};

struct FitRangesIteratorRelative : FitRangesIteratorBlock
{
  using FitRangesIteratorBlock::FitRangesIteratorBlock;
  virtual ~FitRangesIteratorRelative() {}
  FitRangesIteratorElement * Clone() const override
  { return new FitRangesIteratorRelative( fitRange, Parent, ti, tf ); }
  void SetStart() override;
  bool PastEnd() const override;
  bool Increment( bool bWrap ) override;
};

END_MLU_NAMESPACE
#endif // FitRangeImp_hpp
