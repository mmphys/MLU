/*************************************************************************************
 
 Common utilities (with dependencies on c++ stdlib and LatAnalyze)
 
 Source file: CommonLatAn.cpp
 
 Copyright (C) 2019
 
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

#include "CommonLatAn.hpp"

BEGIN_COMMON_NAMESPACE

// Read a list of bootstrapped correlators into a single correlator
static void ReadBootstrapCorrsHelper( Latan::DMat &m, std::vector<double> ReadBuffer, int iFile, int NtSource, int NtDest, int Fold, int Shift )
{
  for( int k = 0; k < NtDest; ++k ) {
    // NB: This is the wrong order! But reflects a bug in LatAnalyze
    // Presumably it was originally written this way, and now can't change because all data are in this format
    //const double re{ReadBuffer[ k ]};
    //const double im{ReadBuffer[ k + dim[0] * 1 ]};
    double Source1 = ReadBuffer[ ( NtSource + Shift + k ) % NtSource ];
    int Dest    = iFile * NtDest + k;
    if( Fold == 0 ) {
      m( Dest, 0 ) = Source1;
    }
    else
    {
      double Source2 = ReadBuffer[ ( 2 * NtSource + Shift - k ) % NtSource ];
      if( Fold > 0 )
        m( Dest, 0 ) = ( Source1 + Source2 ) * 0.5;
      else
        m( Dest, 0 ) = ( std::abs( Source2 ) + std::abs( Source1 ) ) * 0.5;
        //m( Dest, 0 ) = ( Source2 - Source1 ) * 0.5;
    }
  }
}

Latan::DMatSample ReadBootstrapCorrs( const std::vector<std::string> & FileName, int Fold, int Shift, int NumOps )
{
  assert( abs( Fold ) <= 2 && "Error: Invalid Fold parameter" );
  const bool bAlternateFold{ abs( Fold ) > 1 };
  if( bAlternateFold ) {
    if( Fold > 0 )
      --Fold;
    else
      ++Fold;
  }
  const int FirstFold{ Fold };
  const int NumFiles{ static_cast<int>( FileName.size() ) };
  // These will be initialised when we read in the first correlator
  int NtSource = -999; // Number of timeslices in the source file(s)
  int NtDest   = -888; // Number of timeslices after folding
  int NSamples = -777; // Number of bootstrap samples per correlator
  std::vector<double> ReadBuffer;
  Latan::DMatSample   Corr;
  for( int i = 0; i < NumFiles; ++i ) {
    std::cout << "Reading correlator from " << FileName[i];
    Common::H5File f( FileName[i], H5F_ACC_RDONLY );
    H5::Group gBoot = f.openGroup( std::string("/") );
    std::cout << ", group ";
    std::string GroupName{ Common::GetFirstGroupName( gBoot ) };
    std::cout << GroupName << "\n";
    gBoot = gBoot.openGroup( GroupName );
    {
      // Check the attributes of this bootstrap sample
      short int thisType;
      H5::Attribute a = gBoot.openAttribute("type");
      a.read(H5::PredType::NATIVE_SHORT, &thisType);
      assert( thisType == 2 && "Error: bad type attribute" );
      int thisNSample;
      a = gBoot.openAttribute("nSample");
      a.read(H5::PredType::NATIVE_INT, &thisNSample);
      assert( thisNSample > 0 && "Error: bad number of samples" );
      if( i != 0 )
        assert( thisNSample == NSamples && "Error: inconsistent number of samples" );
      else {
        // Initialise defaults first time around
        NSamples = thisNSample;
        std::cout << NSamples << " bootstrap samples per correlator\n";
      }
    }
    H5::DataSet ds = gBoot.openDataSet("data_C");
    H5::DataSpace dsp = ds.getSpace();
    const int nDims{dsp.getSimpleExtentNdims()};
    assert( nDims == 2 && "Error: Wrong number of dimensions");
    hsize_t dim[2];
    dsp.getSimpleExtentDims( dim );
    assert( dim[1] == 2 && "Correlator should contain real and imaginary parts for each timeslice" );
    if( i != 0 ) {
      assert( NtSource == dim[0] && "Error: inconsistent number of timeslices" );
      if( bAlternateFold ) {
        Fold = FirstFold;
        const int row{ i / NumOps };
        const int col{ i % NumOps };
        if( row & 1 )
          Fold *= -1;
        if( col & 1 )
          Fold *= -1;
      }
    }
    else {
      // Initialise defaults first time around
      assert( dim[0] > 0 && dim[0] <= std::numeric_limits<int>::max() && "Error: invalid number of timeslices" );
      NtSource = static_cast<int>( dim[0] );
      NtDest   = Fold ? NtSource / 2 + 1 : NtSource;
      std::cout << NtSource << " timeslices per correlator";
      if( Fold )
        std::cout << " folded into " << NtDest << " timeslices with "
                  << (bAlternateFold ? "alternating, " : "")
                  << (Fold == 1 ? "positive" : "negative") << " parity";
      std::cout << std::endl;
      ReadBuffer.resize( dim[0] * dim[1] );
      Corr.resize( NSamples );
      Corr.resizeMat( NtDest * NumFiles, 1 );
    }
    ds.read( ReadBuffer.data(), H5::PredType::NATIVE_DOUBLE );
    ReadBootstrapCorrsHelper( Corr[Latan::central], ReadBuffer, i, NtSource, NtDest, Fold, Shift );
    for( int j = 0; j < NSamples; ++j ) {
      std::string dsName{ "data_S_" };
      dsName.append( std::to_string( j ) );
      ds = gBoot.openDataSet( dsName );
      dsp = ds.getSpace();
      assert( dsp.getSimpleExtentNdims() == 2 && "Wrong number of dimensions" );
      dsp.getSimpleExtentDims( dim );
      assert( dim[0] == NtSource && dim[1] == 2 && "Error: Inconsistent Correlator dimensions" );
      ds.read( ReadBuffer.data(), H5::PredType::NATIVE_DOUBLE );
      ReadBootstrapCorrsHelper( Corr[j], ReadBuffer, i, NtSource, NtDest, Fold, Shift );
    }
  }
  return Corr;
}

static const char sep[] = " ";
static const char Comment[] = "# ";
static const char NewLine[] = "\n";

// Make summary files of a bootstrap of a correlator
void SummariseBootstrap(const Latan::DMatSample &out, const std::string & sOutFileBase,
     const std::string & sType, Latan::SeedType Seed, const std::vector<std::string> ParamNames )
{
/*  assert( std::isnan( Common::NaN ) && "Compiler does not support quiet NaNs" );
  const int numParams{ static_cast<int>( ParamNames.size() ) };
  assert( numParams == out[Latan::central].rows() && "Parameter names missing" );
  const std::size_t   nSample{ static_cast<std::size_t>( out.size() )};
  std::vector<double> Data( nSample );
  std::string sOutFileName{ MakeFilename(sOutFileBase, SummaryNames[f], Seed, TEXT_EXT) };
  std::ofstream s( sOutFileName );
  s << Comment << sType << NewLine << Comment << sOutFileBase << "\n# Seed " << Seed
  << NewLine << FieldNames << ( ( f == 0 ) ? FieldNames2 : "" );
  s << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
  for(int t = 0; t < nt; t++)
  {
  }
*/}

const char * SummaryNames[3] = { "corr", "mass", "cosh" };//, "sinh" };

// Make summary files of a bootstrap of a correlator
void SummariseBootstrapCorr(const Latan::DMatSample &out, const std::string & sOutFileBase, Latan::SeedType Seed )//, int momentum_squared)
{
  static const char FieldNames[] = "t y y_low y_high y_check";
  static const char FieldNames2[] = " im im_low im_high im_check";
  static const char * SummaryHeader[] =
  {
    "correlator",
    "mass",
    "cosh mass",
    //"sinh mass",
  };
  const int NumSummaries{ static_cast<int>( sizeof( SummaryNames ) / sizeof( SummaryNames[0] ) ) };
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  const int           nt{static_cast<int>(out[Latan::central].rows())};
  const std::size_t   nSample{ static_cast<std::size_t>( out.size() )};
  std::vector<double> Data( nSample );
  for(int f = 0; f < NumSummaries; f++)
  {
    std::string sOutFileName{ MakeFilename(sOutFileBase, SummaryNames[f], Seed, TEXT_EXT) };
    std::ofstream s( sOutFileName );
    s << Comment << SummaryHeader[f] << NewLine << Comment << sOutFileBase << "\n# Seed " << Seed
      << NewLine << FieldNames << ( ( f == 0 ) ? FieldNames2 : "" )
      << std::setprecision(std::numeric_limits<double>::digits10+2) << std::endl;
    for(int t = 0; t < nt; t++)
    {
      std::size_t Count = 0;
      if( f == 0 )
      {
        for( Latan::Index i = 0; i < nSample; i++ )
          if( std::isfinite( out[i](t, 0) ) )
            Data[Count++] = out[i](t, 0);
        Common::ValWithEr Re( out[Latan::central](t,0), Data, Count );
        std::size_t ImCount = 0;
        for( Latan::Index i = 0; i < nSample; i++ )
          if( std::isfinite( out[i](t, 1) ) )
            Data[ImCount++] = out[i](t, 1);
        Common::ValWithEr Im( out[Latan::central](t,1), Data, ImCount );
        s << t << sep << Re.Central << sep << Re.ErLow << sep << Re.ErHigh
               << sep << ( static_cast<double>( Count ) / nSample )
               << sep << Im.Central << sep << Im.ErLow << sep << Im.ErHigh
               << sep << ( static_cast<double>( ImCount ) / nSample ) << std::endl;
      }
      else
      {
        double dCentral = 0;
        for( Latan::Index i = Latan::central; i < static_cast<Latan::Index>( nSample ); i++ )
        {
          double DThis;
          switch(f)
          {
            case 1: // mass
              DThis = std::log( abs( out[i](t, 0) / out[i]((t + 1 + nt) % nt, 0) ) );
              /*assert( momentum_squared >= 0 && momentum_squared <= 6 && "Unsupported momentum" );
               if( momentum_squared >= 1 && momentum_squared <= 6 )
               {
               // Lattice dispersion relation (not very well implemented)
               static const double ss_half{ sin(0.5) * sin(0.5) };
               static const double ss_one{ sin(1.) * sin(1.) };
               double sum;
               if( momentum_squared <= 3 )
               sum = ss_half * momentum_squared;
               else
               sum = ss_half * ( momentum_squared - 4 ) + ss_one;
               double dSinh_HalfE0{ sinh( DThis / 2 ) };
               DThis = 2 * asinh( sqrt( dSinh_HalfE0 * dSinh_HalfE0 - sum ) );
               }*/
              break;
            case 2: // cosh mass
              DThis = std::acosh((out[i]((t - 1 + nt) % nt, 0) + out[i]((t + 1) % nt, 0)) / (2 * out[i](t, 0)));
              break;
              //case 3: // sinh mass
              //DThis = std::asinh((out[i]((t - 1 + nt) % nt, 0) - out[i]((t + 1) % nt, 0)) / (2 * //out[i](t, 0)));
              //break;
            default:
              DThis = 0;
          }
          if( i == Latan::central )
            dCentral = DThis;
          else if( std::isfinite( DThis ) )
            Data[Count++] = DThis;
        }
        Common::ValWithEr Re( dCentral, Data, Count );
        s << t << sep << Re.Central << sep << Re.ErLow << sep << Re.ErHigh << sep
          << ( static_cast<double>( Count ) / nSample ) << std::endl;
      }
    }
  }
}

END_COMMON_NAMESPACE
