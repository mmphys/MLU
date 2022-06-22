/**

 Mike's lattice QCD utilities: Fit ranges
 
 Source file: FitRange.cpp

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

#include "FitRange.hpp"
#include <list>
#include <cstdint>

FitRange_hpp

/*****************************************************************

 FitRanges

*****************************************************************/

// Deserialise a set of fit ranges
void FitRanges::Deserialise( const std::vector<std::string> &vString )
{
  clear();
  reserve( vString.size() );
  vIndex.resize( vString.size() );
  try
  {
    // Deserialise all the fit ranges
    Base New;
    New.reserve( vString.size() );
    for( std::size_t i = 0; i < vString.size(); ++i )
      New.emplace_back( FitRange::Deserialise( vString[i] ) );
    // Now sort them so that fit ranges appear BEFORE the fit ranges they depend on
    // Use Kahn's algorithm, https://en.wikipedia.org/wiki/Topological_sorting
    struct DepTrack
    {
      std::unique_ptr<FitRange> Me;
      std::size_t Index;
      FitRange::vDepend DependsOn;
      DepTrack( std::unique_ptr<FitRange> &&me, std::size_t index ) : Me{ std::move( me ) }, Index{index}
        { DependsOn = Me->GetDependencies(); }
      DepTrack( DepTrack &&o ) : Me{ std::move( o.Me ) }, Index{o.Index}, DependsOn{std::move(o.DependsOn)} {}
    };
    std::list<DepTrack> HasEdges;
    std::list<DepTrack> NoEdges;
    // 1. Split into those with/out edges
    for( std::size_t i = 0; i < New.size(); ++i )
    {
      DepTrack dt( std::move( New[i] ), i );
      if( dt.DependsOn.empty() )
        NoEdges.emplace_back( std::move( dt ) );
      else
        HasEdges.emplace_back( std::move( dt ) );
    }
    // 2. Remove dependencies
    std::list<DepTrack> Sorted;
    while( !NoEdges.empty() )
    {
      // This can go into our sorted list
      Sorted.emplace_front( std::move( NoEdges.front() ) );
      NoEdges.pop_front();
      std::size_t ThisIndex{ Sorted.front().Index };
      for( auto dt = HasEdges.begin(); dt != HasEdges.end(); )
      {
        // Remove any dependencies on pThis
        for( auto i = dt->DependsOn.begin(); i != dt->DependsOn.end(); )
        {
          if( *i == ThisIndex )
            i = dt->DependsOn.erase( i );
          else
            ++i;
        }
        // If no dependencies left, move to NoEdges
        if( dt->DependsOn.empty() )
        {
          NoEdges.emplace_back( std::move( *dt ) );
          dt = HasEdges.erase( dt );
        }
        else
          ++dt;
      }
    }
    if( !HasEdges.empty() )
      throw std::runtime_error( "Fit times have cyclic dependency" );
    std::size_t i = 0;
    for( DepTrack &dt : Sorted )
    {
      push_back( std::move( dt.Me ) );
      vIndex[dt.Index] = i++;
    }
  }
  catch( ... )
  {
    clear();
    throw;
  }
}

/*****************************************************************

 FitRange

*****************************************************************/

// Deserialise a single fit range, by trying all possible deserialisations in turn

FitRange * FitRange::Deserialise( const std::string &String )
{
  using DeserialiseFunc = FitRange * (*)( std::istringstream &is );
  static const std::array<DeserialiseFunc, 2> aDeserialise{ &FitRangeAbsolute::Deserialise,
    &FitRangeRelative::Deserialise };

  for( DeserialiseFunc f : aDeserialise )
  {
    try
    {
      std::istringstream is( String );
      FitRange * p = ( *f )( is );
      if( p )
      {
        if( is.eof() || ( is >> std::ws && is.eof() ) )
          return p;
        delete p;
      }
    }
    catch(...){}
  }
  throw std::runtime_error( "Unrecognised fit range: " + String );
}

std::ostream & operator<<( std::ostream &os, const FitRange &fr )
{
  fr.Print( os );
  return os;
}

/*****************************************************************

 FitRangesIterator

*****************************************************************/

//FitRanges can be iterated by FitRangesIterator
FitRangesIterator FitRanges::begin() const
{
  return FitRangesIterator( *this, false );
}

FitRangesIterator FitRanges::end() const
{
  return FitRangesIterator( *this, true );
}

// Construct from FitRange at either start or end
FitRangesIterator::FitRangesIterator( const FitRanges &ranges_, bool bEnd )
: Base( ranges_.size() ), Ranges{ ranges_ }, RangeMemOrder{ ranges_ }
{
  for( std::size_t i = Ranges.size(); i--; )
  {
    if( bEnd && i == Ranges.size() - 1 )
      RangeMemOrder[i]->GetEnd( Base::operator[]( i ), *this );
    else
      RangeMemOrder[i]->GetStart( Base::operator[]( i ), *this );
  }
}

// Prefix increment
FitRangesIterator &FitRangesIterator::operator++()
{
  // Don't increment if we are already at the end
  if( !PastEnd() )
  {
    std::size_t i = 0;
    while( i < size() && RangeMemOrder[i]->Increment( Base::operator[]( i ), *this, i != size() - 1 ) )
      ++i;
    while( i-- )
    {
      if( !RangeMemOrder[i]->GetDependencies().empty() )
        RangeMemOrder[i]->GetStart( Base::operator[]( i ), *this );
    }
  }
  return *this;
}

// String representation
std::string FitRangesIterator::to_string( const std::string &Sep1, const std::string &Sep2 ) const
{
  std::string s;
  for( std::size_t i = 0; i < size(); ++i )
  {
    if( i )
      s.append( Sep2 );
    s.append( std::to_string( (*this)[i].ti ) );
    s.append( Sep1 );
    s.append( std::to_string( (*this)[i].tf ) );
  }
  return s;
}

std::string FitRangesIterator::to_string() const
{
  std::string s;
  for( std::size_t i = 0; i < size(); ++i )
  {
    if( i )
      s.append( ", " );
    s.append( std::to_string( i ) );
    s.append( " [" );
    s.append( std::to_string( (*this)[i].ti ) );
    s.append( ", " );
    s.append( std::to_string( (*this)[i].tf ) );
    s.append( "]" );
  }
  return s;
}

/*****************************************************************

 FitRangeAbsolute - ti:tf[dti[:dtf]]
 All combinations of Start [ti,ti+dti) ... [tf,tf+dtf)

*****************************************************************/

bool FitRangeAbsolute::Validate( int Nt ) const
{
  return !( ti < 0 || ti >= Nt || tf < ti || tf >= Nt || dti < 1 || dtf < 1 || dti >= Nt || dtf >= Nt
          || ti + dti - 1 >= Nt || tf + dtf - 1 >= Nt );
}

void FitRangeAbsolute::GetStart( FitTime &ft, const FitRangesIterator &it ) const
{
  ft.ti = ti;
  ft.tf = tf;
}

void FitRangeAbsolute::GetEnd( FitTime &ft, const FitRangesIterator &it ) const
{
  ft.ti = ti;
  ft.tf = tf + dtf;
}

bool FitRangeAbsolute::PastEnd( const FitTime &ft, const FitRangesIterator &it ) const
{
  return ft.tf >= tf + dtf;
}

bool FitRangeAbsolute::Increment( FitTime &ft, const FitRangesIterator &it, bool bWrap ) const
{
  bool bOverflow{ false };
  if( ++ft.ti >= ti + dti )
  {
    ft.ti = ti;
    if( ++ft.tf >= tf + dtf )
    {
      if( bWrap )
        ft.tf = tf;
      bOverflow = true;
    }
  }
  return bOverflow;
}

void FitRangeAbsolute::Print( std::ostream &os ) const
{
  char colon{ ':' };
  os << ti << colon << tf;
  if( dti != dtf )
    os << colon << dti << colon << dtf;
  else if( dti != 1 )
    os << colon << dti;
}

FitRange * FitRangeAbsolute::Deserialise( std::istringstream &is )
{
  int Numbers[4];
  int i{ 0 };
  for( bool bMore = true; bMore && i < 4 && is >> Numbers[i]; )
  {
    if( ++i < 4 && !is.eof() && is.peek() == ':' )
      is.get();
    else
      bMore = false;
  }
  FitRangeAbsolute * p = nullptr;
  if( i >= 2 )
  {
    p = new FitRangeAbsolute( Numbers[0], Numbers[1] );
    if( i >= 3 )
    {
      p->dti = Numbers[2];
      p->dtf = ( i == 3 ? Numbers[2] : Numbers[3] );
    }
  }
  return p;
}

/*****************************************************************

 FitRangeRelative - Rn:ti:tf[dti[:dtf]]
 All combinations of Range n ti + [ti,ti+dti) ... Range n tf +  [tf,tf+dtf)

*****************************************************************/

bool FitRangeRelative::Validate( int Nt ) const
{
  return !( std::abs( ti ) >= Nt || std::abs( tf ) >= Nt || dti < 1 || dtf < 1 || dti >= Nt || dtf >= Nt );
}

void FitRangeRelative::GetStart( FitTime &ft, const FitRangesIterator &it ) const
{
  ft.ti = it[DependsOn[0]].ti + ti;
  ft.tf = it[DependsOn[0]].tf + tf;
}

void FitRangeRelative::GetEnd( FitTime &ft, const FitRangesIterator &it ) const
{
  ft.ti = it[DependsOn[0]].ti + ti;
  ft.tf = it[DependsOn[0]].tf + tf + dtf;
}

bool FitRangeRelative::PastEnd( const FitTime &ft, const FitRangesIterator &it ) const
{
  return ft.tf >= it[DependsOn[0]].tf + tf + dtf;
}

bool FitRangeRelative::Increment( FitTime &ft, const FitRangesIterator &it, bool bWrap ) const
{
  bool bOverflow{ false };
  if( ++ft.ti >= it[DependsOn[0]].ti + ti + dti )
  {
    ft.ti = it[DependsOn[0]].ti + ti;
    if( ++ft.tf >= it[DependsOn[0]].tf + tf + dtf )
    {
      if( bWrap )
        ft.tf = it[DependsOn[0]].tf + tf;
      bOverflow = true;
    }
  }
  return bOverflow;
}

void FitRangeRelative::Print( std::ostream &os ) const
{
  char colon{ ':' };
  os << 'R' << DependsOn[0] << colon << ti << colon << tf;
  if( dti != dtf )
    os << colon << dti << colon << dtf;
  else if( dti != 1 )
    os << colon << dti;
}

FitRange * FitRangeRelative::Deserialise( std::istringstream &is )
{
  FitRangeRelative * p = nullptr;
  is >> std::ws;
  char c = is.get();
  if( c != std::istringstream::traits_type::eof() && std::toupper( c ) == 'R' )
  {
    static constexpr int Count{ 5 };
    int Numbers[Count];
    int i{ 0 };
    for( bool bMore = true; bMore && i < Count && is >> Numbers[i]; )
    {
      if( ++i < Count && !is.eof() && is.peek() == ':' )
        is.get();
      else
        bMore = false;
    }
    if( i >= 3 )
    {
      p = new FitRangeRelative( Numbers[0], Numbers[1], Numbers[2] );
      if( i >= 4 )
      {
        p->dti = Numbers[3];
        p->dtf = ( i == 4 ? Numbers[3] : Numbers[4] );
      }
    }
  }
  return p;
}

FitRange_hpp_end
