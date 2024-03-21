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
#include <MLU/MLU.hpp>

using scalar = double;
using Complex = std::complex<scalar>;
using Fold = MLU::Fold<scalar>;
using SC = MLU::Sample<Complex>;
using VScalar = MLU::Vector<scalar>;
using VComplex = MLU::Vector<Complex>;
using MScalar = MLU::Matrix<scalar>;
using MComplex = MLU::Matrix<Complex>;
using JBScalar = MLU::JackBoot<scalar>;
using JBComplex = MLU::JackBoot<Complex>;
using CTraits = MLU::SampleTraits<Complex>;

/** Copy Src to Dst.
 
 `bRealImagOnly == true` just copy real or imaginary component - but don't fold.
 `bRealImagOnly == false` fold forward and backward waves as specified
 */
void CopyFold( MScalar &dst, const MComplex src, int NtHalf,
               MLU::FoldProp &f, bool bRealImagOnly )
{
  std::size_t NumRows{ src.size1 };
  std::size_t Nt{ src.size2 };
  if( dst.size1 != NumRows )
    throw std::runtime_error( "CopyFold() bug: NumRows mismatch" );
  if( dst.size2 != NtHalf )
    throw std::runtime_error( "CopyFold() bug: NtHalf mismatch" );
  const int RI{ f.rps.reality == MLU::Reality::Imag ? 1 : 0 };
  if( bRealImagOnly )
  {
    for( std::size_t i = 0; i < NumRows; ++i )
      for( std::size_t j = 0; j < Nt; ++j )
        dst(i,j) = CTraits::RealImag( src(i,j), RI );
  }
  else
  {
    const scalar Norm{ f.rps.sign == MLU::Sign::Negative ? -0.5 : 0.5 };
    for( std::size_t i = 0; i < NumRows; ++i )
      for( std::size_t t = 0; t < NtHalf; ++t )
      {
        scalar d = CTraits::RealImag( src(i,t), RI );
        if( t == 0 )
        {
          // Timeslice 0 is a little special
          if( f.t0Abs )
            d = std::abs( d );
        }
        else
        {
          // All the other timeslices
          if( f.rps.parity == MLU::Parity::Odd )
            d -= CTraits::RealImag( src( i, Nt - t ), RI );
          else
            d += CTraits::RealImag( src( i, Nt - t ), RI );
          d *= Norm;
        }
        dst(i,t) = d;
      }
  }
}

void CopyFold( VScalar &dst, const VComplex &src, int NtHalf,
               MLU::FoldProp &f, bool bRealImagOnly )
{
  MScalar mDst;
  MComplex mSrc;
  mDst.MapView( dst );
  mSrc.MapView( const_cast<VComplex &>( src ) );
  CopyFold( mDst, mSrc, NtHalf, f, bRealImagOnly );
}

int main(int argc, const char *argv[])
{
  // We're not using C-style I/O
  std::ios_base::sync_with_stdio(false);
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = MLU::CommandLine;
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
    MLU::MakeAncestorDirs( outPrefix );
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
        MLU::FoldProp f;
        std::size_t pos = Arg.find_last_of(',');
        bool bRealImagOnly{ true };
        if( pos == std::string::npos )
          inFileName.append( Arg );
        else
        {
          bRealImagOnly = f.Parse( &Arg[pos + 1], Arg.length() - pos - 1 );
          inFileName.append( Arg.c_str(), pos );
        }
        // I always pick either the real or imaginary component of source
        if( f.rps.reality != MLU::Reality::Imag )
          f.rps.reality = MLU::Reality::Real;
        if( !bRealImagOnly )
        {
          // If I'm folding, then I need to choose odd or even and the sign
          if( f.rps.parity != MLU::Parity::Odd )
            f.rps.parity = MLU::Parity::Even;
          if( f.rps.sign != MLU::Sign::Negative )
            f.rps.sign = MLU::Sign::Positive;
        }
        std::cout << "FoldProp: " << f << std::endl;
        static const char pIndent[] = "  ";
        std::vector<std::string> FileList{ MLU::glob( &inFileName, &inFileName + 1 ) };
        for( const std::string & FileName : FileList )
        {
          if( !MLU::FileExists( FileName ) )
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
            if( !f.Conjugate )
              myFileList = in.FileList;
            else
            {
              // Try to load the file with source and sink swapped
              if( MLU::EqualIgnoreCase( OpNames[0], OpNames[1] ) )
                throw std::runtime_error( "Folding with conjugate specified, but sink and source are both " + OpNames[0] );
              std::string ConjFileName{ in.Name_.Dir };
              ConjFileName.append( in.Name_.Base );
              ConjFileName.append( MLU::Underscore );
              ConjFileName.append( OpNames[in.Name_.op[0]] );
              ConjFileName.append( MLU::Underscore );
              ConjFileName.append( OpNames[in.Name_.op[1]] );
              ConjFileName.append( MLU::Period );
              ConjFileName.append( in.Name_.Type );
              ConjFileName.append( MLU::Period );
              ConjFileName.append( in.Name_.SeedString );
              ConjFileName.append( MLU::Period );
              ConjFileName.append( in.Name_.Ext );
              SC in2{ ConjFileName, "+ ", &OpNames };
              in.IsCompatible( in2 );
              const std::size_t FSize { in .FileList.size() };
              const std::size_t FSize2{ in2.FileList.size() };
              if( FSize2 == FSize )
              {
                // Same length - merge the lists like a zipper
                std::string sAlt;
                sAlt.append( OpNames[in.Name_.op[1]] );
                sAlt.append( MLU::Underscore );
                sAlt.append( OpNames[in.Name_.op[0]] );
                sAlt.append( MLU::Comma );
                sAlt.append( OpNames[in.Name_.op[0]] );
                sAlt.append( MLU::Underscore );
                sAlt.append( OpNames[in.Name_.op[1]] );
                myFileList.reserve( FSize + FSize2 );
                auto p1 = in .FileList.begin();
                auto p2 = in2.FileList.begin();
                for( std::size_t i = 0; i < FSize; ++i, ++p1, ++p2 )
                {
                  bool bSame{ MLU::EqualIgnoreCase(*p1, *p2) };
                  myFileList.emplace_back( *p1 );
                  myFileList.emplace_back( bSame ? sAlt : *p2 );
                }
              }
              else
              {
                std::cout << pIndent << "Warning: bootstrap FileLists uneven length, "
                          << std::to_string( FSize2 ) << MLU::sNE << std::to_string( FSize )
                          << ". Appending lists" << std::endl;
                myFileList = in.FileList;
                myFileList.insert( myFileList.end(), in2.FileList.begin(), in2.FileList.end() );
              }
              for( int s = SC::idxCentral; s < NumSamples; s++ )
                for( int t = 0; t < Nt; ++t )
                  in(s,t) = ( in(s,t) + in2(s,t) ) * 0.5;
              out.BootstrapList.emplace_back( ConjFileName );
              // Merge the raw data (if present and same size)
              if( in.NumSamplesRaw() )
              {
                if( in.NumSamplesRaw() != in2.NumSamplesRaw() )
                  in.resizeRaw( 0 );
                else
                {
                  for( int s = 0; s < in.NumSamplesRaw(); s++ )
                    for( int t = 0; t < Nt; ++t )
                      in.getRaw()(s,t) = ( in.getRaw()(s,t) + in2.getRaw()(s,t) ) * 0.5;
                }
              }
              // Merge the binned data (if present and same size)
              if( in.NumSamplesBinned() )
              {
                if( in.NumSamplesBinned() != in2.NumSamplesBinned() )
                  in.resizeBinned( 0 );
                else
                {
                  for( int s = 0; s < in.NumSamplesBinned(); s++ )
                    for( int t = 0; t < Nt; ++t )
                      in.getBinned()(s,t) = ( in.getBinned()(s,t) + in2.getBinned()(s,t) ) * 0.5;
                }
              }
            }
            // Now fold the correlator, obeying the fold properties
            const int NtHalf{ bRealImagOnly ? Nt : Nt / 2 + ( f.rps.parity == MLU::Parity::Odd ? 0 : 1 ) };
            out.resize( NumSamples, NtHalf );
            out.FileList = std::move( myFileList );
            out.CopyAttributes( in );
            out.NtUnfolded_ = Nt;
            out.t0Negated = false;
            out.Conjugated = f.Conjugate;
            out.reality = f.rps.reality;
            out.parity = f.rps.parity;
            out.sign = f.rps.sign;
            // Copy the resampled data
            JBScalar &jbDst{ out.getData() };
            const JBComplex &jbSrc{ in.getData() };
            CopyFold( jbDst.GetCentral(), jbSrc.GetCentral(), NtHalf, f, bRealImagOnly );
            CopyFold( jbDst.GetReplicaMean(), jbSrc.GetReplicaMean(), NtHalf, f, bRealImagOnly );
            CopyFold( jbDst.Replica, jbSrc.Replica, NtHalf, f, bRealImagOnly );
            if( in.NumSamplesRaw() )
            {
              out.resizeRaw( in.NumSamplesRaw() );
              CopyFold( out.getRaw(), in.getRaw(), NtHalf, f, bRealImagOnly );
            }
            if( in.NumSamplesBinned() )
            {
              out.resizeBinned( in.NumSamplesBinned() );
              CopyFold( out.getBinned(), in.getBinned(), NtHalf, f, bRealImagOnly );
            }
            out.MakeCorrSummary();
            // Now save the folded correlator
            std::string OutFileName{ outPrefix };
            OutFileName.append( in.Name_.Base );
            OutFileName.append( MLU::Underscore );
            OutFileName.append( OpNames[in.Name_.op[1]] );
            OutFileName.append( MLU::Underscore );
            OutFileName.append( OpNames[in.Name_.op[0]] );
            OutFileName.append( MLU::Period );
            OutFileName.append( MLU::sFold );
            OutFileName.append( MLU::Period );
            OutFileName.append( in.Name_.SeedString );
            OutFileName.append( MLU::Period );
            std::size_t OutLen{ OutFileName.length() };
            OutFileName.append( in.Name_.Ext );
            std::cout << "->" << OutFileName << std::endl;
            out.Write( OutFileName, MLU::sFold.c_str() );
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
