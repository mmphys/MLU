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
      {"a", CL::SwitchType::Single, nullptr },
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
        Common::FoldProp f;
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
            std::vector<std::string> myFileList;
            out.binSize = in.binSize;
            out.BootstrapList.clear();
            out.BootstrapList.emplace_back( FileName );
            if( !f.Conjugate )
              myFileList = in.FileList;
            else
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
              static const std::string sNE{ " != " };
              static const std::string sPrefix{ "Incompatible boostrap samples - " };
              std::string sSuffix{ ":\n  " + FileName + "\n+ " + ConjFileName };
              if( in2.NumSamples() != NumSamples )
                throw std::runtime_error( sPrefix + "NumSamples " + std::to_string(in2.NumSamples()) +
                                         sNE + std::to_string(NumSamples) + sSuffix );
              if( in2.Nt() != Nt )
                throw std::runtime_error( sPrefix + "Nt " + std::to_string(in2.Nt()) +
                                         sNE + std::to_string(Nt) + sSuffix );
              if( in2.Seed_ != in.Seed_ )
                throw std::runtime_error( "Seed " + std::to_string( in2.Seed_ ) + sNE + std::to_string( in2.Seed_ ) );
              if( !Common::EqualIgnoreCase( in2.SeedMachine_, in.SeedMachine_ ) )
                throw std::runtime_error( "Machine " + in2.SeedMachine_ + sNE + in.SeedMachine_ );
              if( in2.SampleSize && in.SampleSize && in2.SampleSize != in.SampleSize )
                throw std::runtime_error( sPrefix + "SampleSize " + std::to_string(in2.SampleSize) +
                                         sNE + std::to_string(in.SampleSize) + sSuffix );
              const std::size_t CSize { in .ConfigCount.size() };
              const std::size_t CSize2{ in2.ConfigCount.size() };
              if( CSize && CSize2 )
              {
                if( CSize2 != CSize )
                  throw std::runtime_error( sPrefix + "Number of configs " +
                                      std::to_string(CSize2) + sNE + std::to_string(CSize) + sSuffix );
                for( std::size_t i = 0; i < CSize; i++ )
                {
                  const Common::ConfigCount &l{ in .ConfigCount[i] };
                  const Common::ConfigCount &r{ in2.ConfigCount[i] };
                  if( r.Config != l.Config )
                    throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                             sNE + std::to_string(l.Config) + sSuffix );
                  if( r.Count != l.Count )
                    throw std::runtime_error( sPrefix + "Config " + std::to_string(r.Config) +
                                             ", NumTimeslices " + std::to_string(r.Count) + sNE +
                                             std::to_string(l.Count) + sSuffix );
                }
              }
              const std::size_t FSize { in .FileList.size() };
              const std::size_t FSize2{ in2.FileList.size() };
              if( FSize2 == FSize )
              {
                // Same length - merge the lists like a zipper
                std::string sAlt;
                sAlt.append( OpNames[in.Name_.op[1]] );
                sAlt.append( Common::Underscore );
                sAlt.append( OpNames[in.Name_.op[0]] );
                sAlt.append( Common::Comma );
                sAlt.append( OpNames[in.Name_.op[0]] );
                sAlt.append( Common::Underscore );
                sAlt.append( OpNames[in.Name_.op[1]] );
                myFileList.reserve( FSize + FSize2 );
                auto p1 = in .FileList.begin();
                auto p2 = in2.FileList.begin();
                for( std::size_t i = 0; i < FSize; ++i, ++p1, ++p2 )
                {
                  bool bSame{ Common::EqualIgnoreCase(*p1, *p2) };
                  myFileList.emplace_back( *p1 );
                  myFileList.emplace_back( bSame ? sAlt : *p2 );
                }
              }
              else
              {
                std::cout << pIndent << "Warning: bootstrap FileLists uneven length, "
                          << std::to_string( FSize2 ) << sNE << std::to_string( FSize )
                          << ". Appending lists" << std::endl;
                myFileList = in.FileList;
                myFileList.insert( myFileList.end(), in2.FileList.begin(), in2.FileList.end() );
              }
              std::complex<scalar> * dst{ in[SC::idxCentral] };
              const std::complex<scalar> * src{ in2[SC::idxCentral] };
              for( int s = SC::idxCentral; s < NumSamples; s++ )
                for( int t = 0; t < Nt; t++, dst++, src++ )
                  *dst = ( *dst + *src ) * 0.5;
              out.binSize *= 2;
              out.BootstrapList.emplace_back( ConjFileName );
            }
            // Now fold the correlator, obeying the fold properties
            out.NtUnfolded = Nt;
            out.parity = f.rps.parity;
            out.reality = f.rps.reality;
            out.sign = f.rps.sign;
            out.t0Negated = false;
            out.Conjugated = f.Conjugate;
            const int NtHalf{ Nt / 2 + 1 };
            out.resize( NumSamples, NtHalf );
            out.SampleSize = in.SampleSize;
            out.ConfigCount = in.ConfigCount;
            out.FileList = std::move( myFileList );
            scalar * dst{ out[Fold::idxCentral] };
            const scalar * src{ SC::Traits::ScalarPtr( in[SC::idxCentral] )};
            if( f.rps.reality == Common::Reality::Imag )
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
                if( f.rps.parity == Common::Parity::Odd )
                  d = src[t * SC::scalar_count] - src[(Nt - t) * SC::scalar_count];
                else
                  d = src[t * SC::scalar_count] + src[(Nt - t) * SC::scalar_count];
                *dst++ = d * ( f.rps.sign == Common::Sign::Negative ? -0.5 : 0.5 );
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
    "-a     Perform Average rather than treat forward/backward/conjugate separately\n"
    "-i     Input prefix\n"
    "-o     Output prefix\n"
    "--help This message\n"
    "and <FoldOptions> is a case-insensitive string consisting of zero or more of\n"
    "r/i    Real (default) or imaginary\n"
    "e/o    Even (default) or Odd\n"
    "p/n0    Positive (default) or Negative in first half of timeslices\n"
    "0      Disable taking the absolute value of timeslice 0. Default: take abs(c(0))\n"
    "c      Fold Conjugate operators (i.e. swap source and sink) together\n";
  }
  return iReturn;
}
