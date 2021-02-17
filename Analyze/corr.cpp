/*************************************************************************************
 
 Choose real or imaginary part of correlator from bootstrap sample
 Create statistics, ready for GNUplot
 
 Source file: corr.cpp
 
 Copyright (C) 2020-21
 
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
#include <MLU/Common.hpp>

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
        Common::FoldProp f;
        std::size_t pos = Arg.find_last_of(',');
        bool bRealImagOnly{ true };
        if( pos == std::string::npos )
          inFileName.append( Arg );
        else
        {
          bRealImagOnly = f.Parse( &Arg[pos + 1], Arg.length() - pos - 1 );
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
            out.BootstrapList.clear();
            out.BootstrapList.emplace_back( FileName );
            int BinSize2{ 0 };
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
              SC in2{ ConjFileName, "+ ", &OpNames };
              BinSize2 = in2.binSize;
              in.IsCompatible( in2 );
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
                          << std::to_string( FSize2 ) << Common::sNE << std::to_string( FSize )
                          << ". Appending lists" << std::endl;
                myFileList = in.FileList;
                myFileList.insert( myFileList.end(), in2.FileList.begin(), in2.FileList.end() );
              }
              std::complex<scalar> * dst{ in[SC::idxCentral] };
              const std::complex<scalar> * src{ in2[SC::idxCentral] };
              for( int s = SC::idxCentral; s < NumSamples; s++ )
                for( int t = 0; t < Nt; t++, dst++, src++ )
                  *dst = ( *dst + *src ) * 0.5;
              out.BootstrapList.emplace_back( ConjFileName );
            }
            // Now fold the correlator, obeying the fold properties
            const int NtHalf{ bRealImagOnly ? Nt : Nt / 2 + ( f.rps.parity == Common::Parity::Odd ? 0 : 1 ) };
            out.resize( NumSamples, NtHalf );
            out.FileList = std::move( myFileList );
            out.CopyAttributes( in );
            out.NtUnfolded = Nt;
            out.t0Negated = false;
            out.Conjugated = f.Conjugate;
            scalar * dst{ out[Fold::idxCentral] };
            const scalar * src{ SC::Traits::ScalarPtr( in[SC::idxCentral] )};
            if( f.rps.reality == Common::Reality::Imag )
              src++;
            else
              f.rps.reality = Common::Reality::Real;
            if( bRealImagOnly )
            {
              std::size_t Len{ ( static_cast<std::size_t>( NumSamples ) + 1 ) * Nt };
              for( std::size_t i = 0; i < Len; ++i )
              {
                *dst++ = *src;
                src += SC::scalar_count;
              }
            }
            else
            {
              out.binSize += BinSize2;
              if( f.rps.parity != Common::Parity::Odd )
                f.rps.parity = Common::Parity::Even;
              if( f.rps.sign != Common::Sign::Negative )
                f.rps.sign = Common::Sign::Positive;
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
            }
            out.reality = f.rps.reality;
            out.parity = f.rps.parity;
            out.sign = f.rps.sign;
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
    "r/i    Real (default) or imaginary\n"
    "e/o    Even (default) or Odd\n"
    "p/n0   Positive (default) or Negative in first half of timeslices\n"
    "0      Disable taking the absolute value of timeslice 0. Default: take abs(c(0))\n"
    "c      Fold Conjugate operators (i.e. swap source and sink) together\n"
    "NB: If <FoldOptions> is missing (or only contains r/i), no folding occurs\n";
  }
  return iReturn;
}
