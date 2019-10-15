/*
 * bootstrap.cpp - based on ...
 * resample.cpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2016 Antonin Portelli
 *   Modified           2019 Michael Marshall
 *
 * LatAnalyze 3 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatAnalyze 3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatAnalyze 3.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "CommonGrid.hpp"
#include "CommonLatAn.hpp"

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Io/Hdf5File.hpp>
#include <H5CompType.h> // HDF5 Library

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif
// Default output file format for bootstrap data
#ifndef DEF_FMT
#define DEF_FMT "h5"
#endif
// Default output file format for text based summaries of bootstrap data
#ifndef TEXT_EXT
#define TEXT_EXT "txt"
#endif

// Make a copy of str, in which token has been replaced by x
template <typename T> std::string tokenReplaceCopy(const std::string &str, const std::string token, const T &x )
{
  std::string sCopy{str};
  Latan::tokenReplace(sCopy, token, x );
  return sCopy;
}

// Make a filename of a specific type
std::string MakeFilename(const std::string &Base, const char *Type, Latan::SeedType Seed, const std::string &Ext)
{
  std::string type{Type};
  type.append( 1, '.' );
  type.append( std::to_string( Seed ) );
  return tokenReplaceCopy( Base, "type", type ) + '.' + Ext;
}

static const char * SummaryNames[] = { "corr", "mass", "cosh" };//, "sinh" };

// Make summary files of this data set
void MakeSummaries(const Latan::DMatSample &out, const std::string & sOutFileBase, Latan::SeedType Seed )//, int momentum_squared)
{
  const int           nt{static_cast<int>(out[Latan::central].rows())};
  const Latan::Index  nSample{out.size()};
  // Now make summary data files
  std::vector<double> TmpAvg( nt );
  std::vector<std::vector<double>> TmpSamples( nt );
  for(int t = 0; t < nt; t++)
    TmpSamples[t].resize( nSample );
  static const char sep[] = " ";
  static const char * SummaryHeader[] =
  {
    "# correlator\nt corr_re corr_stddev corr_im corr_im_stddev",
    "# "   "mass\nt mass mass_low mass_high check",
    "# cosh_mass\nt mass mass_low mass_high check",
    //"# sinh_mass\nt mass mass_low mass_high check",
  };
  for(int f = 0; f < sizeof(SummaryNames)/sizeof(SummaryNames[0]); f++)
  {
    std::string sOutFileName{ MakeFilename(sOutFileBase, SummaryNames[f], Seed, TEXT_EXT) };
    std::ofstream s(sOutFileName);
    s << std::setprecision(std::numeric_limits<double>::digits10 + 1) << SummaryHeader[f] << std::endl;
    for(int t = 0; t < nt; t++)
    {
      if( f == 0 )
      {
        double Avg = 0;
        double StdDev = 0;
        double AvgIm = 0;
        double StdDevIm = 0;
        for(int i = 0; i < nSample; i++)
        {
          Avg   += out[i](t,0);
          AvgIm += out[i](t,1);
        }
        Avg   /= nSample;
        AvgIm /= nSample;
        for(int i = 0; i < nSample; i++)
        {
          double d = out[i](t,0) - Avg;
          StdDev += d * d;
          d = out[i](t,1) - AvgIm;
          StdDevIm += d * d;
        }
        StdDev   /= nSample;
        StdDevIm /= nSample;
        StdDev    = std::sqrt( StdDev );
        StdDevIm  = std::sqrt( StdDevIm );
        s << t << sep << Avg << sep << StdDev << sep << AvgIm << sep << StdDevIm << sep << std::endl;
      }
      else
      {
        TmpAvg[t] = 0;
        int Count = 0;
        for(int i = 0; i < nSample; i++)
        {
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
            if( std::isfinite( DThis ) ) {
              TmpSamples[t][Count++] = DThis;
              TmpAvg[t] += DThis;
            }
          }
        }
        static const double NaN{ std::nan( "" ) };
        assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
        double dErrPlus, dErrMinus;
        if( Count ) {
          TmpAvg[t] /= Count;
          const typename std::vector<double>::iterator itStart{ TmpSamples[t].begin() };
          const typename std::vector<double>::iterator itEnd  { itStart + Count };
          std::sort( itStart, itEnd );
          int Index = static_cast<int>( 0.16 * Count + 0.5 );
          dErrMinus = TmpAvg[t] - TmpSamples[t][Index];
          dErrPlus  = TmpSamples[t][Count - Index] - TmpAvg[t];
        } else {
          TmpAvg[t] = NaN;
          dErrPlus  = NaN;
          dErrMinus = NaN;
        }
        s << t << sep << TmpAvg[t] << sep << dErrMinus << sep << dErrPlus << sep << ( static_cast<double>( Count ) / nSample ) << std::endl;
      }
    }
  }
}

// Show me the averages for each timeslice
void ShowTimeSliceAvg(const Latan::Dataset<Latan::DMat> &data) {
  const int nFile{static_cast<int>(data.size())};
  if( !nFile )
    std::cout << "No timeslices to average over!" << std::endl;
  else {
    const int nt{static_cast<int>(data[0].rows())};
    std::cout << "Averages of " << nFile << " files over " << nt << " timeslices:" << std::endl;
    std::vector<std::complex<double>> Avg(nt);
    for (unsigned int t = 0; t < nt; ++t) {
      Avg[t] = 0;
      for (unsigned int j = 0; j < nFile; ++j) {
        std::complex<double> d( data[j](t,0), data[j](t,1) );
        //std::cout << "d[" << j << "]=" << d << std::endl;
        Avg[t] += d;
        //std::cout << "a[" << j << "]=" << a << std::endl;
      }
      Avg[t] /= nFile;
      std::cout << "C(" << t << ")=" << Avg[t] << std::endl;
    }
  }
}

struct BootstrapParams
{
  std::string outStem;
  std::string OutFileExt;
  Latan::SeedType seed;
  Latan::Index nSample;
  Latan::Index binSize;
  bool bShowAverages;
  bool bSaveSummaries;
  bool bSaveBootstrap;
  bool dumpBoot;
};

bool PerformBootstrap( const BootstrapParams &par, Latan::Dataset<Latan::DMat> &in, const std::string &Contraction )//, int momentum_squared )
{
  const bool bPerformBootstrap{ par.bSaveSummaries || par.bSaveBootstrap };
  //if( !bPerformBootstrap && !par.dumpBoot )
    //throw new std::runtime_error( "No bootstrap to perform" );
  // Skip bootstrap if output exists
  static const char * pszBootstrap{"bootstrap"};
  static const std::string sBootstrap{pszBootstrap};
  const std::string sOutFileBase{ tokenReplaceCopy( par.outStem, "corr", Contraction ) };
  std::string sOutFileName{ MakeFilename( sOutFileBase, pszBootstrap, par.seed, par.OutFileExt ) };
  if( Common::FileExists( sOutFileName ) || (!par.bSaveBootstrap &&
      Common::FileExists( MakeFilename(sOutFileBase, SummaryNames[0], par.seed, TEXT_EXT) )) )
  {
    std::cout << "Skipping " << Contraction << " because output already exists" << std::endl;
    return false;
  }
  // Make sure there are at least 10x more bootstrap samples than data points
  std::cout << Contraction << std::endl;
  if( par.bShowAverages )
    ShowTimeSliceAvg( in );
  if( bPerformBootstrap ) {
    if( ( par.nSample / 10 ) < in.size() )
      throw new std::runtime_error( "nSample=" + std::to_string( par.nSample ) + " too small for dataset with " + std::to_string( in.size() ) + " elements" );
    std::cout << "-- resampling (" << par.nSample << " samples)..." << std::endl;
    if( par.binSize != 1 )
      in.bin( par.binSize );
    Latan::DMatSample out = in.bootstrapMean( par.nSample, par.seed );
    if( par.bSaveBootstrap ) {
      std::cout << "Saving sample to " << sOutFileName << std::endl;
      Latan::Io::save<Latan::DMatSample>(out, sOutFileName, Latan::File::Mode::write, sBootstrap);
    }
    if( par.bSaveSummaries )
      MakeSummaries( out, sOutFileBase, par.seed );//, momentum_squared );
  }
  if( par.dumpBoot ) { // Dump the bootstrap sequences
    std::string SeqFileName{ MakeFilename( sOutFileBase, "bootseq", par.seed, TEXT_EXT ) };
    std::cout << "Saving bootstrap sequence to '" << SeqFileName << "' ...";
    std::ofstream file(SeqFileName);
    file << "# bootstrap sequences" << std::endl;
    file << "#      bin size: " << par.binSize << std::endl;
    file << "#          seed: " << std::to_string( par.seed ) << std::endl;
    in.dumpBootstrapSeq(file, par.nSample, par.seed);
    std::cout << std::endl;
  }
  return true;
}

/*struct Momentum {
  int x;
  int y;
  int z;
  Momentum( int _x, int _y, int _z ) : x(_x), y(_y), z(_z) {}
};

struct SubStrings {
  std::string              Name;
  std::vector<std::string> Strings;
};

const std::string & MesonSuffix( const std::string MesonTypeName )
{
  static const std::string axial{ "ax" };
  static const std::string axial_name{ "_ax" };
  static const std::string other_name{};
  return axial.compare( MesonTypeName ) == 0 ? axial_name : other_name;
}*/

// Make a default manifest
// Heavy-light meson decays. 2pt function with current inserted and momentum on light meson

void Z2Bootstrap( const BootstrapParams &bsParams, const Common::Momentum &mom, const std::string &Heavy, const std::string &Light )
{
  static constexpr unsigned int Nt{64};
  static constexpr unsigned int CStart{3000};
  static constexpr unsigned int CSkip{40};
  static constexpr unsigned int CCount{10};
  static constexpr unsigned int CEnd{CStart + CSkip * CCount};
  static const std::string sPath{
#ifdef DEBUG
    "/Users/s1786208"               // Laptop
#else
    "/home/dp008/dp008/dc-mars3"  // Tesseract
#endif
    "/data/201910Plat/mesons/C1/Z2/" };
  static const std::string Dot{ "." };    // used inside filenames
  static const std::string H5{ Dot + "h5" };    // used inside filenames
  static const std::string Sep{ "_" };    // used inside filenames
  static const std::string Space{ " " };  // whitespace as field separator / human readable info
  static const std::string Sink{ Sep + "sink" };    // used inside filenames
  static const std::string sProp{ "prop" + Sep };    // used inside filenames
  static const std::string sMom0{ Sep + "p" + Sep + Common::Momentum(0,0,0).to_string( Sep ) };
  const std::string sMom{ Sep + "p" + Sep + mom.to_string( Sep )};
  const std::string sMomNeg{ Sep + "p" + Sep + mom.to_string( Sep, true )};
  using Algebra = Grid::Gamma::Algebra;
  static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
  static const std::vector<std::string> algNames = { "g5", "gT5" };
  assert( alg.size() == algNames.size() && "Abbreviations should match gamma algebra" );
  static const int NumAlgebra{ static_cast<int>( alg.size() ) };
  static const int NumCorr{ NumAlgebra * NumAlgebra };
  const std::string CorrPrefix{ Heavy + Sep + Light + sMom };
  std::vector<std::string> CorrSuffixes( NumCorr );
  for( int iSrc = 0; iSrc < NumAlgebra; iSrc++ )
    for( int iSnk = 0; iSnk < NumAlgebra; iSnk++ ) {
      const int iCorr{ iSnk * NumAlgebra + iSrc };
      std::stringstream ss;
      ss << Sep << algNames[iSnk] << Sep << algNames[iSrc];
      CorrSuffixes[iCorr] = ss.str();
    }
  const unsigned int NumEnds{ Common::EqualIgnoreCase( Heavy, Light ) ? 1u : 2u };
  const std::size_t NumSamplesT{ NumEnds * CCount };
  const std::size_t NumSamples{ NumSamplesT * Nt };
  std::vector<Latan::Dataset<Latan::DMat>> bsData( NumCorr );
  std::vector<Latan::Dataset<Latan::DMat>> bsDataT( NumCorr );
  std::vector<Common::Correlator> buffer( NumCorr );
  for( int iCorr = 0 ; iCorr < NumCorr; iCorr++ ) {
    bsData[ iCorr ].resize( NumSamples );
    bsDataT[ iCorr ].resize( NumSamplesT );
    buffer[iCorr].resize( Nt );
  }
  int CorrIndex{ 0 };
  BootstrapParams par{ bsParams };
  par.bSaveSummaries = true;
  par.bSaveBootstrap = false;
  for( int t = 0; t < Nt; ++t ) {
    const std::string tPrefix{ Sep + "t" + std::to_string( t ) };
    int CorrIndexT{ 0 };
    for( unsigned int iEnd = 0; iEnd < NumEnds; iEnd++ ) {
      const std::string & Left{ iEnd == 0 ? Heavy : Light };
      const std::string & Right{ iEnd == 0 ? Light : Heavy };
      for( int iConfig = CStart; iConfig < CEnd; iConfig += CSkip ) {
        const std::string sFileName{ sPath + sProp + Left + tPrefix + sMom + Sep + sProp + Right + tPrefix + sMom0 + Sink + sMomNeg + Dot + std::to_string( iConfig ) + H5};
        Common::ReadComplexArray(buffer, alg, sFileName, 0);
        for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
        {
          Common::CopyCorrelator( bsData[iCorr][CorrIndex], buffer[iCorr], t );
          Common::CopyCorrelator( bsDataT[iCorr][CorrIndexT], buffer[iCorr] );
        }
        CorrIndex++;
        CorrIndexT++;
      }
    }
    // Now perform a bootstrap of this timeslice
    for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
      PerformBootstrap( par, bsDataT[iCorr], CorrPrefix + CorrSuffixes[iCorr] + tPrefix );
  }
  // Now perform a bootstrap overall
  par.bSaveBootstrap = true;
  for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
    PerformBootstrap( par, bsData[iCorr], CorrPrefix + CorrSuffixes[iCorr] );
}

/*****************************************************************

 Perform a bootstrap of all the files specified on the command line.
 Organise all the files by contraction, and sort them by trajectory.
 Ensure trajectories are processed in the same order so the
 bootstrap replicas for each contraction can be combined later during analysis

*****************************************************************/

static const char DefaultOutputStem[] = "@corr@.@type@";

static void Usage( const char * pExecutableName, Latan::OptParser &opt )
{
  std::string cmd{ pExecutableName };
  auto pos = cmd.find_last_of('/');
  if( pos != std::string::npos && pos < cmd.length() - 1 )
    cmd = cmd.substr( pos + 1 );
  std::cerr << "usage: " << pExecutableName << " <options> ContractionFile1 [ContractionFile2 ...]"
  << "\nOutput stem (if present) must contain '@corr@' and '@type', e.g. '" << DefaultOutputStem << "'"
  << "\nPossible options:\n" << opt << std::endl;
}

/*int DebugTest1()
{
  static const std::string Filename{ "prop_h1_t0_p_0_0_0_prop_h1_t0_p_0_0_0_sink_p_0_0_0.3000.h5" };
  MesonFile f{ Filename, 0, false, 0, Grid::Gamma::Algebra::GammaZ,
    Grid::Gamma::Algebra::GammaT };
  std::vector<std::complex<double>> buf( 64 );
  f.GetCorellator( buf );
  return 1;
}*/

int main(int argc, char *argv[])
{
  /*int i = DebugTest1();
  if( i )
    return i;*/
  // parse command line/options ////////////////////////////////////////////////////////
  Latan::OptParser       opt;
  opt.addOption("n", "nsample"   , Latan::OptParser::OptType::value,   true,
                "number of samples", DEF_NSAMPLE);
  opt.addOption("b", "bin"       , Latan::OptParser::OptType::value,   true,
                "bin size", "1");
  opt.addOption("r", "seed"      , Latan::OptParser::OptType::value,   true,
                "random generator seed (default: random)");
  opt.addOption("o", "output"    , Latan::OptParser::OptType::value,   true,
                "output file", DefaultOutputStem);
  opt.addOption("f", "format"    , Latan::OptParser::OptType::value,   true,
                "output file format", DEF_FMT);
  opt.addOption("d", "dump-boot" , Latan::OptParser::OptType::trigger, true,
                "dump bootstrap sequence");
  opt.addOption("t", "txt"       , Latan::OptParser::OptType::trigger, true,
                "don't save text versions of correlators (for GnuPlot)");
  opt.addOption("a", "average"   , Latan::OptParser::OptType::trigger, true,
                "show correlator averages");
  opt.addOption("x", "exclude"   , Latan::OptParser::OptType::value,   true,
                "exclude files (colon separated list)", "");
  opt.addOption("" , "help"      , Latan::OptParser::OptType::trigger, true,
                "show this help message and exit");
  bool parsed = opt.parse(argc, argv) and !opt.gotOption("help") /*and opt.getArgs().size() > 0*/;
  // If the output stem was specified, ensure it includes "@corr@"
  BootstrapParams par;
  par.bSaveBootstrap = true;
  if( parsed )
  {
    par.outStem = opt.optionValue("o");
    if (par.outStem.find("@corr@") == std::string::npos or par.outStem.find("@type@") == std::string::npos)
    {
      parsed=false;
      std::cout << "Output stem " << par.outStem << " invalid" << std::endl;
    }
  }
  if (!parsed) {
    Usage( argv[0], opt );
    return EXIT_FAILURE;
  }

  par.nSample = opt.optionValue<Latan::Index>("n");
  par.binSize = opt.optionValue<Latan::Index>("b");
  if (opt.gotOption("r"))
    par.seed = opt.optionValue<Latan::SeedType>("r");
  else {
    std::random_device rd;
    par.seed = rd();
  }
  par.OutFileExt = opt.optionValue("f");
  par.dumpBoot = opt.gotOption("d");
  par.bShowAverages = opt.gotOption("a");
  par.bSaveSummaries = !opt.gotOption("t");
  std::cout << "Creating bootstrap output " << par.outStem << ", seed=" << par.seed << std::endl;

  int iReturn = EXIT_SUCCESS;
  try
  {
    ///////////////////////////
    // If there are files specified on the command line,
    // parse the input file names, grouping by correlator, indexed by trajectory.
    ///////////////////////////
    if( opt.getArgs().size() )
    {
      Common::Manifest Manifest{ opt.getArgs(), opt.optionValue("x") };
      // Walk the list of contractions, performing a separate bootstrap for each
      int BootstrapCount = 0;
      Latan::DMatSample out(par.nSample);
      for( auto itc = Manifest.begin(); itc != Manifest.end(); itc++ )
      {
        const std::string &Contraction{itc->first};
        const Common::TrajList &l{itc->second};
        const unsigned int nFile{ static_cast<unsigned int>(l.TrajFile.size())};
        std::cout << "Loading " << nFile << " files for " << Contraction << std::endl;
        unsigned int j = 0; // which trajectory are we processing
        bool bError = false;
        Common::Correlator buf;
        Latan::Dataset<Latan::DMat> data(nFile);
        for( auto it = l.TrajFile.begin(); it != l.TrajFile.end(); it++, j++ )
        {
          const int &traj{it->first}; // Trajectory number
          const Common::TrajFile &tf{*it->second};
          const std::string &Filename{tf.Filename};
          std::cout << "\t" << traj << "\t" << Filename << std::endl;
          try
          {
            Common::ReadComplexArray( buf, Filename );
            Common::CopyCorrelator( data[j], buf );
            PerformBootstrap( par, data, Contraction );
            BootstrapCount++;
          } catch( const std::exception &e ) {
            std::cerr << "\tError: " << e.what() << std::endl;
          } catch( ... ) {
            std::cerr << "\tError while reading " << Filename << std::endl;
            bError = true;
          }
        }
      }
      std::cout << "Bootstraps written for " << BootstrapCount
      << " / " << Manifest.size() << " correlators" << std::endl;
    } else {
      ///////////////////////////
      // Nothing specified on the command-line
      // Perform default bootstraps
      // Heavy-light semi-leptonics
      ///////////////////////////
      const std::string qHeavy{ "h1" };
      const std::string qLight{ "l" };
      std::cout << "No files specified. Making default manifest for "
      << qHeavy << " and " << qLight << std::endl;
      static const Common::Momentum Z2Momenta[] = {
        { 0, 0, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
        { 1, 1, 1 },
        { 2, 0, 0 },
      };
      for( const Common::Momentum &m : Z2Momenta )
      {
        Z2Bootstrap( par, m, qHeavy, qHeavy );
        Z2Bootstrap( par, m, qHeavy, qLight );
        Z2Bootstrap( par, m, qLight, qLight );
      }
    }
  }
  catch(const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  }
  return iReturn;
}
