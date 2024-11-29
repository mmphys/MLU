/**

 Mike's lattice QCD utilities: Fit ranges
 
 Source file: FitRange.cpp

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
#include "FitRangeImp.hpp"
#include <list>
#include <cstdint>

BEGIN_MLU_NAMESPACE

/*****************************************************************

 FitRanges

*****************************************************************/

// Deserialise a set of fit ranges
void FitRanges::Deserialise( const std::vector<std::string> &vString, int MinDP )
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
      New.emplace_back( FitRange::Deserialise( vString[i], i, MinDP ) );
    // Now sort them so that fit ranges appear BEFORE the fit ranges they depend on
    // Use Kahn's algorithm, https://en.wikipedia.org/wiki/Topological_sorting
    struct DepTrack
    {
      std::unique_ptr<FitRange> Me;
      std::size_t Index;
      std::vector<std::size_t> DependsOn;
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

FitRangesIterator FitRanges::begin() const
{
  return FitRangesIterator( *this, false );
}

FitRangesIterator FitRanges::end() const
{
  return FitRangesIterator( *this, true );
}

std::ostream & operator<<( std::ostream &os, const FitRanges &fr )
{
  for( std::size_t i = 0; i < fr.size(); ++i )
  {
    if( i )
      os << ", ";
    os << i << " [" << fr[i] << "]";
  }
  return os;
}

/*****************************************************************

 FitRange

*****************************************************************/

// Deserialise a single fit range, by trying all possible deserialisations in turn

FitRange * FitRange::Deserialise( const std::string &String, std::size_t MyIndex, int MinDP )
{
  using DeserialiseFunc = FitRange * (*)( std::istringstream &, std::size_t );
  static const std::array<DeserialiseFunc, 2> aDeserialise{ &FitRangeAbsolute::Deserialise,
    &FitRangeRelative::Deserialise };
  if( MinDP < 0 )
    throw std::runtime_error( "Minimum data points " + std::to_string( MinDP ) + " < 0" );
  for( DeserialiseFunc f : aDeserialise )
  {
    try
    {
      std::istringstream is( String );
      std::unique_ptr<FitRange> p{ ( *f )( is, MyIndex ) };
      if( p )
      {
        p->MinDP = MinDP;
        int NextChar{ is.get() };
        if( NextChar == std::istringstream::traits_type::eof()
          || ( std::tolower( NextChar ) == 'd'
              && ( is >> p->MinDP )
              && ( is.eof() || ( is >> std::ws && is.eof() ) ) ) )
        {
          if( p->MinDP < 0 )
            throw std::runtime_error( "Minimum data points " + std::to_string( p->MinDP ) + " < 0" );
          return p.release();
        }
      }
    }
    catch(...){}
  }
  throw std::runtime_error( "Unrecognised fit range: " + String );
}

std::vector<int> FitRange::GetColonList( std::istream &is, std::size_t MaxLen )
{
  std::vector<int> v;
  while( !is.eof() && ( is >> std::ws && !is.eof() ) && v.size() < MaxLen )
  {
    if( !v.empty() )
    {
      if( is.peek() != ':' )
        break;
      is.get();
    }
    int i;
    if( is >> i )
      v.push_back( i );
    else
    {
      if( !v.empty() )
        throw std::runtime_error( "FitRange::GetColonList() colon followed by non-integer" );
      break;
    }
  }
  return v;
}

/*****************************************************************

 FitRangesIterator

*****************************************************************/

FitRangesIterator::FitRangesIterator( const FitRangesIterator &it ) : Ranges{ it.Ranges }
{
  reserve( it.size() );
  for( const std::unique_ptr<FitRangesIteratorElement> &elem : it )
    emplace_back( elem->Clone() );
}

// Construct from FitRange at either start or end
FitRangesIterator::FitRangesIterator( const FitRanges &ranges_, bool bEnd )
: Base( ranges_.size() ), Ranges{ ranges_ }
{
  Base &base{ *this };
  for( std::size_t i = Ranges.size(); i--; )
  {
    Element * e;
    if( bEnd && i == Ranges.size() - 1 )
      e = Ranges(i).GetEnd( *this );
    else
      e = Ranges(i).GetStart( *this );
    base[i].reset( e );
  }
  if( !bEnd && !GotMinDP() )
    operator++();
}

bool FitRangesIterator::PastEnd() const
{
  return back()->PastEnd();
}

// Prefix increment
FitRangesIterator &FitRangesIterator::operator++()
{
  // Keep incrementing until minimum data points satisfied or we hit the end
  bool bMinDPOK{ false };
  while( !bMinDPOK && !PastEnd() )
  {
    std::size_t i = 0;
    while( i < size() && (*this)(i).Increment( i != size() - 1 ) )
      ++i;
    while( i-- )
    {
      if( !Ranges(i).GetDependencies().empty() )
        (*this)(i).SetStart();
    }
    // See whether minimum number of data points for each range are satisfied
    bMinDPOK = GotMinDP();
  }
  return *this;
}

// String representation
std::string FitRangesIterator::AbbrevString( const std::string &Sep1, const std::string &Sep2 ) const
{
  std::string s;
  for( std::size_t i = 0; i < size(); ++i )
  {
    if( i )
      s.append( Sep2 );
    s.append( (*this)[i].AbbrevString( Sep1 ) );
  }
  return s;
}

inline bool FitRangesIterator::GotMinDP() const
{
  bool bMinDPOK = true;
  for( std::size_t i = 0; bMinDPOK && i < size(); ++i )
    bMinDPOK = (*this)(i).Extent() >= Ranges(i).MinDP;
  return bMinDPOK;
}

std::string FitRangesIteratorElement::AbbrevString( const std::string &Sep ) const
{
  std::string s;
  int i_, f_;
  if( FitTimes.empty() )
  {
    i_ = 1;
    f_ = 0;
  }
  else
  {
    i_ = FitTimes[0];
    f_ = FitTimes[0];
    for( std::size_t j = 1; j < FitTimes.size(); ++j )
    {
      if( i_ > FitTimes[j] )
        i_ = FitTimes[j];
      if( f_ < FitTimes[j] )
        f_ = FitTimes[j];
    }
  }
  s.append( std::to_string( i_ ) );
  s.append( Sep );
  s.append( std::to_string( f_ ) );
  return s;
}

std::ostream & operator<<( std::ostream &os, const FitRangesIterator &it )
{
  for( std::size_t i = 0; i < it.size(); ++i )
  {
    if( i )
      os<< ", ";
    os << i << " [" << it[i] << "]";
  }
  return os;
}

bool FitRangesIteratorBlock::SetFitTimes()
{
  bool bSameLastTime{ false };
  const int Extent{ tf - ti + 1 };
  const FitRangeBlock &r{ dynamic_cast<const FitRangeBlock &>( fitRange ) };
  if( tf < ti )
    FitTimes.clear();
  else if( r.Thinning.empty() )
  {
    FitTimes.resize( Extent );
    for( int j = 0; j < Extent; ++j )
      FitTimes[j] = ti + j;
  }
  else
  {
    // Non-zero extent + thinning
    std::vector<int> v;
    int t = ti;
    for( std::size_t Thindex=0; t <= tf && Thindex < r.Thinning.size(); Thindex += 2 )
    {
      const bool bLast{ r.Thinning.size() - Thindex < 2 };
      const int Step{ r.Thinning[Thindex] };
      if( Step == 0 )
        t = bLast ? tf + 1 : t + r.Thinning[Thindex + 1]; // Skip past this number of timeslices
      else
      {
        int Num{ Thindex + 1 < r.Thinning.size() ? r.Thinning[Thindex + 1]
                                                 : ( tf - t ) / Step + 1 };
        while( Num-- && t <= tf )
        {
          v.push_back( t );
          t += Step;
        }
      }
    }
    bSameLastTime = FitTimes == v;
    if( !bSameLastTime )
      FitTimes = std::move( v );
  }
  return bSameLastTime;
}

void FitRangesIteratorBlock::Print( std::ostream &os ) const
{
  const FitRangeBlock &r{ dynamic_cast<const FitRangeBlock &>( fitRange ) };
  const bool bPrintThinned{ !r.Thinning.empty() && Extent() != FitTimes.size() };
  if( bPrintThinned )
  {
    for( std::size_t i = 0; i < FitTimes.size(); ++i )
    {
      if( i )
        os << ";";
      os << FitTimes[i];
    }
  }
  else
    Base::Print( os );
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

FitRangesIterator::Element * FitRangeAbsolute::GetStart( const FitRangesIterator &Parent ) const
{
  FitRangesIteratorAbsolute *it = new FitRangesIteratorAbsolute( *this, Parent );
  it->SetStart();
  return it;
}

FitRangesIterator::Element * FitRangeAbsolute::GetEnd( const FitRangesIterator &Parent ) const
{
  FitRangesIteratorAbsolute *it = new FitRangesIteratorAbsolute( *this, Parent );
  it->ti = ti;
  it->tf = tf + dtf;
  return it;
}

void FitRangeAbsolute::Print( std::ostream &os ) const
{
  char colon{ ':' };
  os << ti << colon << tf;
  if( dti != dtf )
    os << colon << dti << colon << dtf;
  else if( dti != 1 )
    os << colon << dti;
  ShowThinning( os );
}

void FitRangeBlock::GetThinning( std::istream &is )
{
  if( !is.eof() && is >> std::ws && !is.eof() && std::tolower( is.peek() ) == 't' )
  {
    is.get();
    Thinning = GetColonList( is );
    if( Thinning.empty() )
      throw std::runtime_error( "Thinning unrecognised" );
    for( std::size_t i = 0; i < Thinning.size(); i += 2 )
    {
      if( Thinning[i] < 0 || ( i == 0 && Thinning[i] == 0 ) ) // First spec can't be a skip
        throw std::runtime_error("Thinning separation " + std::to_string( Thinning[i] ) + " invalid");
      if( i < Thinning.size() - 1 && Thinning[i + 1] <= 0 )
        throw std::runtime_error("Thinning run " + std::to_string( Thinning[i + 1] ) + " invalid");
    }
  }
}

void FitRangeBlock::ShowThinning( std::ostream &os ) const
{
  if( !Thinning.empty() )
  {
    os << 't';
    for( std::size_t i = 0; i < Thinning.size(); ++i )
    {
      if( i )
        os << ':';
      os << Thinning[i];
    }
  }
}

FitRange * FitRangeAbsolute::Deserialise( std::istringstream &is, std::size_t MyIndex )
{
  std::vector<int> Numbers = GetColonList( is, 4 );
  std::unique_ptr<FitRangeAbsolute> p;
  if( Numbers.size() >= 2 )
  {
    p.reset( new FitRangeAbsolute( Numbers[0], Numbers[1] ) );
    if( Numbers.size() >= 3 )
    {
      p->dti = Numbers[2];
      p->dtf = Numbers.back();
    }
    p->GetThinning( is );
  }
  return p.release();
}

/*****************************************************************

 Iterate over FitRangeAbsolute

*****************************************************************/

void FitRangesIteratorAbsolute::SetStart()
{
  const FitRangeAbsolute &r{ dynamic_cast<const FitRangeAbsolute &>( fitRange ) };
  ti = r.ti;
  tf = r.tf;
  SetFitTimes();
}

bool FitRangesIteratorAbsolute::PastEnd() const
{
  const FitRangeAbsolute &r{ dynamic_cast<const FitRangeAbsolute &>( fitRange ) };
  return ti >= r.ti + r.dti;
}

bool FitRangesIteratorAbsolute::Increment( bool bWrap )
{
  const FitRangeAbsolute &r{ dynamic_cast<const FitRangeAbsolute &>( fitRange ) };
  bool bOverflow{ false };
  bool bSameAsLastTime;
  do
  {
    if( ++tf >= r.tf + r.dtf )
    {
      tf = r.tf;
      if( ++ti >= r.ti + r.dti )
      {
        if( !bWrap )
        {
          // I'm at the end
          FitTimes.clear();
          return true;
        }
        ti = r.ti;
        bOverflow = true;
      }
    }
    bSameAsLastTime = SetFitTimes();
  }
  while( !bOverflow && bSameAsLastTime );
  return bOverflow;
}

/*****************************************************************

 FitRangeRelative - Rn:ti:tf[dti[:dtf]]
 All combinations of Range n ti + [ti,ti+dti) ... Range n tf +  [tf,tf+dtf)

*****************************************************************/

bool FitRangeRelative::Validate( int Nt ) const
{
  return !( std::abs( ti ) >= Nt || std::abs( tf ) >= Nt || dti < 1 || dtf < 1 || dti >= Nt || dtf >= Nt );
}

FitRangesIterator::Element * FitRangeRelative::GetStart( const FitRangesIterator &Parent ) const
{
  FitRangesIteratorRelative *it = new FitRangesIteratorRelative( *this, Parent );
  it->SetStart();
  return it;
}

FitRangesIterator::Element * FitRangeRelative::GetEnd( const FitRangesIterator &Parent ) const
{
  if( DependsOn.size() != 1 )
    throw std::runtime_error( "FitRangesIteratorRelative::SetStart() Bug DependsOn.size() = " + std::to_string( DependsOn.size() ) );
  FitRangesIteratorRelative *it = new FitRangesIteratorRelative( *this, Parent );
  it->ti = Parent[DependsOn[0]].TI() + ti;
  it->tf = Parent[DependsOn[0]].TF() + tf + dtf;
  return it;
}

void FitRangeRelative::Print( std::ostream &os ) const
{
  char colon{ ':' };
  os << 'R' << DependsOn[0] << colon << ti << colon << tf;
  if( dti != dtf )
    os << colon << dti << colon << dtf;
  else if( dti != 1 )
    os << colon << dti;
  ShowThinning( os );
}

FitRange * FitRangeRelative::Deserialise( std::istringstream &is, std::size_t MyIndex )
{
  FitRangeRelative * p = nullptr;
  is >> std::ws;
  char c = is.get();
  if( c != std::istringstream::traits_type::eof() && std::toupper( c ) == 'R' )
  {
    std::vector<int> Numbers = GetColonList( is, 5 );
    if( Numbers.size() >= 3 )
    {
      p = new FitRangeRelative( Numbers[1], Numbers[2] );
      p->DependsOn.push_back( MyIndex + static_cast<std::size_t>( Numbers[0] ) );
      if( Numbers.size() >= 4 )
      {
        p->dti = Numbers[3];
        p->dtf = Numbers.back();
      }
    }
  }
  return p;
}

/*****************************************************************

 Iterate over FitRangeAbsolute

*****************************************************************/

void FitRangesIteratorRelative::SetStart()
{
  const std::vector<std::size_t> &DependsOn{ fitRange.GetDependencies() };
  if( DependsOn.size() != 1 )
    throw std::runtime_error( "FitRangesIteratorRelative::SetStart() Bug DependsOn.size() = " + std::to_string( DependsOn.size() ) );
  const FitRangeRelative &r{ dynamic_cast<const FitRangeRelative &>( fitRange ) };
  ti = Parent[DependsOn[0]].TI() + r.ti;
  tf = Parent[DependsOn[0]].TF() + r.tf;
  SetFitTimes();
}

bool FitRangesIteratorRelative::PastEnd() const
{
  const std::vector<std::size_t> &DependsOn{ fitRange.GetDependencies() };
  if( DependsOn.size() != 1 )
    throw std::runtime_error( "FitRangesIteratorRelative::SetStart() Bug DependsOn.size() = " + std::to_string( DependsOn.size() ) );
  const FitRangeRelative &r{ dynamic_cast<const FitRangeRelative &>( fitRange ) };
  return ti >= Parent[DependsOn[0]].TI() + r.ti + r.dti;
}

bool FitRangesIteratorRelative::Increment( bool bWrap )
{
  const std::vector<std::size_t> &DependsOn{ fitRange.GetDependencies() };
  if( DependsOn.size() != 1 )
    throw std::runtime_error( "FitRangesIteratorRelative::SetStart() Bug DependsOn.size() = " + std::to_string( DependsOn.size() ) );
  const FitRangeRelative &r{ dynamic_cast<const FitRangeRelative &>( fitRange ) };
  bool bOverflow{ false };
  bool bSameAsLastTime;
  do
  {
    const int OtherTF{ Parent[DependsOn[0]].TF() };
    if( ++tf >= OtherTF + r.tf + r.dtf )
    {
      tf = OtherTF + r.tf;
      const int OtherTI{ Parent[DependsOn[0]].TI() };
      if( ++ti >= OtherTI + r.ti + r.dti )
      {
        if( !bWrap )
        {
          // I'm at the end
          FitTimes.clear();
          return true;
        }
        ti = OtherTI + r.ti;
        bOverflow = true;
      }
    }
    bSameAsLastTime = SetFitTimes();
  }
  while( !bOverflow && bSameAsLastTime );
  return bOverflow;
}

END_MLU_NAMESPACE
