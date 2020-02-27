/*************************************************************************************
 
 Utility for creating statistics for bootstrap sample - ready for GNUplot
 
 Source file: fold.cpp
 
 Copyright (C) 2020
 
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
 *************************************************************************************/
/*  END LEGAL */

#include <stdio.h>
#include <typeinfo>
#include "Common.hpp"

struct FoldProp
{
  Common::Reality reality = Common::Reality::Real;
  Common::Parity parity = Common::Parity::Even;
  Common::Sign sign = Common::Sign::Positive;
  bool t0Abs = true;
  bool Conjugate = false;
  void Parse( const char * const pc, std::size_t const Len )
  {
    bool bOK = true;
    static constexpr int NumOptions{ 8 };
    static const char Options[] = "RIEOPN0C";
    int iCount[NumOptions];
    for( int i = 0; i < NumOptions; i++ )
      iCount[i] = 0;
    for( int i = 0; bOK && i < Len; i++ )
    {
      char c = pc[i];
      if( c >= 'a' && c <= 'z' )
        c -= 'a' - 'A';
      int idx = 0;
      while( idx < NumOptions && c != Options[idx] )
        idx++;
      if( idx < NumOptions )
        iCount[idx]++;
      else
        bOK = false;
    }
    if( bOK )
    {
      if( iCount[0] == 1 && iCount[1] == 0 )
        reality = Common::Reality::Real;
      else if( iCount[0] == 0 && iCount[1] == 1 )
        reality = Common::Reality::Imag;
      else if( iCount[0] != 0 || iCount[1] != 0 )
        bOK = false;
      if( iCount[2] == 1 && iCount[3] == 0 )
        parity = Common::Parity::Even;
      else if( iCount[2] == 0 && iCount[3] == 1 )
        parity = Common::Parity::Odd;
      else if( iCount[2] != 0 || iCount[3] != 0 )
        bOK = false;
      if( iCount[4] == 1 && iCount[5] == 0 )
        sign = Common::Sign::Positive;
      else if( iCount[4] == 0 && iCount[5] == 1 )
        sign = Common::Sign::Negative;
      else if( iCount[4] != 0 || iCount[5] != 0 )
        bOK = false;
      if( iCount[6] == 1 )
        t0Abs = false;
      else if( iCount[6] != 0 )
        bOK = false;
      if( iCount[7] == 1 )
        Conjugate = true;
      else if( iCount[7] != 0 )
        bOK = false;
    }
    if( !bOK )
    {
      const std::string s{ pc, Len };
      throw std::runtime_error( "Option string '" + s + "' invalid" );
    }
  }
  FoldProp() = default;
  explicit FoldProp( const char * pc, std::size_t Len ) { Parse( pc, Len ); }
};

inline std::ostream & operator<<( std::ostream &os, const FoldProp &f )
{
  return os << f.reality << ", " << f.parity << ", " << f.sign << (f.t0Abs ? ", abs( C(0) )" : "" )
            << ( f.Conjugate ? " and conjugate operators" : "" );
}

int main(int argc, const char *argv[])
{
  using scalar = double;
  using Fold = Common::Fold<scalar>;
  using SC = Common::Sample<std::complex<scalar>>;
  // We're not using C-style I/O
  std::ios_base::sync_with_stdio(false);
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"i", CL::SwitchType::Single, "" },
      {"o", CL::SwitchType::Single, "" },
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const std::string outPrefix{ cl.SwitchValue<std::string>( "o" ) };
    std::string inFileName{ cl.SwitchValue<std::string>( "i" ) };
    // If there are files specified on the command line, process each file
    if( !cl.GotSwitch( "help" ) && cl.Args.size() )
    {
      bShowUsage = false;
      const std::size_t inFileNameLen{ inFileName.length() };
      std::size_t Count = 0;
      Fold out;
      for( const std::string &Arg : cl.Args )
      {
        inFileName.resize( inFileNameLen );
        // Look for a comma
        FoldProp f;
        std::size_t pos = Arg.find_last_of(',');
        if( pos == std::string::npos )
          inFileName.append( Arg );
        else
        {
          f.Parse( &Arg[pos + 1], Arg.length() - pos - 1 );
          inFileName.append( Arg.c_str(), pos );
        }
        std::cout << "FoldProp: " << f << std::endl;
        static const char pIndent[] = "  ";
        std::vector<std::string> FileList{ Common::glob( &inFileName, &inFileName + 1 ) };
        for( const std::string & FileName : FileList )
        {
          if( !Common::FileExists( FileName ) )
            std::cout << pIndent << "Error: " << FileName << " doesn't exist" << std::endl;
          else
          {
            std::vector<std::string> OpNames;
            SC in{ FileName, pIndent, &OpNames };
            const int Nt{ in.Nt() };
            const int NumSamples{ in.NumSamples() };
            if( f.Conjugate )
            {
              if( Common::EqualIgnoreCase( OpNames[0], OpNames[1] ) )
                throw std::runtime_error( "Folding with conjugate specified, but sink and source are both " + OpNames[0] );
              std::string ConjFileName{ in.Name_.Dir };
              ConjFileName.append( in.Name_.Base );
              ConjFileName.append( Common::Underscore );
              ConjFileName.append( OpNames[in.Name_.op[0]] );
              ConjFileName.append( Common::Underscore );
              ConjFileName.append( OpNames[in.Name_.op[1]] );
              ConjFileName.append( Common::Period );
              ConjFileName.append( in.Name_.Type );
              ConjFileName.append( Common::Period );
              ConjFileName.append( in.Name_.SeedString );
              ConjFileName.append( Common::Period );
              ConjFileName.append( in.Name_.Ext );
              SC in2{ ConjFileName, "+ " };
              if( in2.NumSamples() != NumSamples || in2.Nt() != Nt )
                throw std::runtime_error( "Incompatible boostrap samples:\n  "
                                         + FileName + "\n+ " + ConjFileName );
              if( in2.Seed_ != in.Seed_ )
                throw std::runtime_error( "Seed " + std::to_string( in2.Seed_ ) + " doesn't match " + std::to_string( in2.Seed_ ) );
              if( !Common::EqualIgnoreCase( in2.SeedMachine_, in.SeedMachine_ ) )
                throw std::runtime_error( "Machine " + in2.SeedMachine_ + " doesn't match " + in.SeedMachine_ );
              std::complex<scalar> * dst{ in[SC::idxCentral] };
              const std::complex<scalar> * src{ in2[SC::idxCentral] };
              for( int s = SC::idxCentral; s < NumSamples; s++ )
                for( int t = 0; t < Nt; t++, dst++, src++ )
                  *dst = ( *dst + *src ) * 0.5;
            }
            // Now fold the correlator, obeying the fold properties
            out.NtUnfolded = Nt;
            out.parity = f.parity;
            out.reality = f.reality;
            out.sign = f.sign;
            out.t0Negated = false;
            out.Conjugated = f.Conjugate;
            const int NtHalf{ Nt / 2 + 1 };
            out.resize( NumSamples, NtHalf );
            scalar * dst{ out[Fold::idxCentral] };
            const scalar * src{ SC::Traits::ScalarPtr( in[SC::idxCentral] )};
            if( f.reality == Common::Reality::Imag )
              src++;
            for( int s = SC::idxCentral; s < NumSamples; s++ )
            {
              // Timeslice 0 is a little special
              scalar d = *src;
              if( s == SC::idxCentral && f.t0Abs && d < 0 )
                out.t0Negated = true;
              if( out.t0Negated )
                d = std::abs( d );
              *dst++ = d;
              // Now do all the other timeslices
              for( int t = 1; t < NtHalf; t++ )
              {
                if( f.parity == Common::Parity::Odd )
                  d = src[t * SC::scalar_count] - src[(Nt - t) * SC::scalar_count];
                else
                  d = src[t * SC::scalar_count] + src[(Nt - t) * SC::scalar_count];
                *dst++ = d * ( f.sign == Common::Sign::Negative ? -0.5 : 0.5 );
              }
              src += SC::scalar_count * Nt;
            }
            out.Seed_ = in.Seed_;
            out.SeedMachine_ = in.SeedMachine_;
            out.MakeCorrSummary( nullptr );
            // Now save the folded correlator
            std::string OutFileName{ outPrefix };
            OutFileName.append( in.Name_.Base );
            OutFileName.append( Common::Underscore );
            OutFileName.append( OpNames[in.Name_.op[1]] );
            OutFileName.append( Common::Underscore );
            OutFileName.append( OpNames[in.Name_.op[0]] );
            OutFileName.append( Common::Period );
            OutFileName.append( Common::sFold );
            OutFileName.append( Common::Period );
            OutFileName.append( in.Name_.SeedString );
            OutFileName.append( Common::Period );
            std::size_t OutLen{ OutFileName.length() };
            OutFileName.append( in.Name_.Ext );
            std::cout << "->" << OutFileName << std::endl;
            out.Write( OutFileName, Common::sFold.c_str() );
            OutFileName.resize( OutLen );
            OutFileName.append( TEXT_EXT );
            out.WriteSummary( OutFileName );
            Count++;
          }
        }
      }
      std::cout << "Folded " << Count << " bootstrapped correlators" << std::endl;
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> BootstrapFile1,<FoldOptions1> [ContractionFile2,<FoldOptions2> ...]\n"
    "Save correlator in format ready for GNUPlot, where <options> are:\n"
    "-i     Input prefix\n"
    "-o     Output prefix\n"
    "--help This message\n"
    "and <FoldOptions> is a case-insensitive string consisting of zero or more of\n"
    "R/I    Real (default) or imaginary\n"
    "E/O    Even (default) or Odd\n"
    "P/N    Positive (default) or Negative in first half of timeslices\n"
    "0      Disable taking the absolute value of timeslice 0. Default: take abs(c(0))\n"
    "C      Fold the Conjugate operators together with this, i.e. swap source and sink\n";
  }
  return iReturn;
}
