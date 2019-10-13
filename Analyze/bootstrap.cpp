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

#include "Common.hpp"

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Io/Hdf5File.hpp>
#include <H5CompType.h>

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
template <typename T> std::string tokenReplaceCopy(const std::string &str, const std::string token, const T &x );
// Show me the averages for each timeslice
void ShowTimeSliceAvg(const Latan::Dataset<Latan::DMat> &data);
// Make summary files of this data set
void MakeSummaries(const Latan::DMatSample &out, const std::string & sOutFileBase, Latan::SeedType Seed, int momentum_squared);
// Make a filename of a specific type
std::string MakeFilename(const std::string &Base, const char *Type, Latan::SeedType Seed, const std::string &Ext);
// Make a default manifest if one not specified
bool ReadCorrelatorsZ2( Common::Manifest & manifest, const std::string &Heavy, const std::string &Light );

/*struct MesonFile : Common::TrajFile
{
  Grid::Gamma::Algebra GammaSnk, GammaSrc;
  MesonFile( const std::string & Filename_, int offset_, bool bImaginary_, int iMultiplier_, Grid::Gamma::Algebra GammaSnk_, Grid::Gamma::Algebra GammaSrc_ )
  : TrajFile( Filename_, offset_, bImaginary_, iMultiplier_ ),
    GammaSnk{ GammaSnk_ }, GammaSrc{ GammaSrc_ } {}
  virtual ~MesonFile() {}
  virtual void GetCorellator( std::vector<std::complex<double>> &buffer );
};*/

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
  std::string outStem;
  if (parsed)
  {
    outStem = opt.optionValue("o");
    if (outStem.find("@corr@") == std::string::npos or outStem.find("@type@") == std::string::npos)
    {
      parsed=false;
      std::cout << "Output stem " << outStem << " invalid" << std::endl;
    }
  }
  if (!parsed) {
    Usage( argv[0], opt );
    return EXIT_FAILURE;
  }

  // Now parse the input file names, grouping by correlator, indexed by trajectory
  Common::Manifest Manifest{opt.getArgs(), opt.optionValue("x")};
  if( Manifest.empty() ) {
    // Heavy-light semi-leptonics
    const std::string qHeavy{ "h1" };
    const std::string qLight{ "l" };
    std::cout << "No files specified. Making default manifest for "
          << qHeavy << " and " << qLight << std::endl;
    //MakeManifestHL( Manifest, qHeavy, qLight );
    bool bOK = ReadCorrelatorsZ2( Manifest, qHeavy, qHeavy )
            && ReadCorrelatorsZ2( Manifest, qHeavy, qLight )
            && ReadCorrelatorsZ2( Manifest, qLight, qLight );
    if( !bOK )
      return EXIT_FAILURE;
    //Usage( argv[0], opt );
    //return EXIT_FAILURE;
  }

  Latan::Index    nSample = opt.optionValue<Latan::Index>("n");
  Latan::Index    binSize = opt.optionValue<Latan::Index>("b");
  Latan::SeedType seed;
  if (opt.gotOption("r"))
    seed = opt.optionValue<Latan::SeedType>("r");
  else {
    std::random_device rd;
    seed = rd();
  }
  const std::string & OutFileExt{opt.optionValue("f")};
  const bool dumpBoot{ opt.gotOption("d") };
  const bool bShowAverages{ opt.gotOption("a") };
  const bool bSaveSummaries{ !opt.gotOption("t") };
  std::cout << "Creating bootstrap output " << outStem << ", seed=" << seed << std::endl;

  // Walk the list of contractions, performing a separate bootstrap for each
  int BootstrapCount = 0;
  Latan::DMatSample out(nSample);
  for( auto itc = Manifest.begin(); itc != Manifest.end(); itc++ )
  {
    static const char * pszBootstrap{"bootstrap"};
    static const std::string sBootstrap{pszBootstrap};
    const std::string &Contraction{itc->first};
    Common::TrajList &l{itc->second};
    const std::string sOutFileBase{tokenReplaceCopy(outStem, "corr", Contraction)};
    std::string sOutFileName{MakeFilename(sOutFileBase, pszBootstrap, seed, OutFileExt)};
    if( Common::FileExists( sOutFileName ) )
      std::cout << "Skipping " << Contraction << " because " << sOutFileName << " already exists" << std::endl;
    else
    {
      // load all the files for this correlator
      std::cout << Contraction << std::endl;
      const unsigned int nFile{ static_cast<unsigned int>( l.data.size() ) };
      const unsigned int nt{ static_cast<unsigned int>( l.data[0].rows() ) };
      out.resizeMat(nt, 2);
      if( bShowAverages )
        ShowTimeSliceAvg(l.data);
      // Resample
      std::cout << "-- resampling (" << nSample << " samples)..." << std::endl;
      if( binSize != 1 )
        l.data.bin(binSize);
      Latan::DMatSample out = l.data.bootstrapMean(nSample, seed);
      std::cout << "Saving sample to '" << sOutFileName << "' ...";
      Latan::Io::save<Latan::DMatSample>(out, sOutFileName, Latan::File::Mode::write, sBootstrap);
      std::cout << " done" << std::endl;
      if( bSaveSummaries )
        MakeSummaries(out, sOutFileBase, seed, l.momentum);
      if( dumpBoot ) { // Dump the bootstrap sequences
        std::string SeqFileName{MakeFilename(sOutFileBase, "bootseq", seed, TEXT_EXT)};
        std::cout << "Saving bootstrap sequence to '" << SeqFileName << "' ...";
        std::ofstream file(SeqFileName);
        file << "# bootstrap sequences" << std::endl;
        file << "#      bin size: " << binSize << std::endl;
        file << "#          seed: " << std::to_string( seed ) << std::endl;
        for( auto it = l.TrajFile.begin(); it != l.TrajFile.end(); it++ )
        {
          const int &traj{it->first};
          //std::string &Filename{it->second};
          const Common::TrajFile &tf{*it->second};
          const std::string &Filename{tf.Filename};
          file << "#   traj / file: " << traj << " " << Filename << std::endl;
        }
        l.data.dumpBootstrapSeq(file, nSample, seed);
        std::cout << std::endl;
      }
      BootstrapCount++;
    }
  }
  std::cout << "Bootstraps written for " << BootstrapCount
            << " / " << Manifest.size() << " correlators" << std::endl;
  return EXIT_SUCCESS;
}

// Make a copy of str, in which token has been replaced by x
template <typename T> std::string tokenReplaceCopy(const std::string &str, const std::string token, const T &x )
{
  std::string sCopy{str};
  Latan::tokenReplace(sCopy, token, x );
  return sCopy;
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

// Make a filename of a specific type
std::string MakeFilename(const std::string &Base, const char *Type, Latan::SeedType Seed, const std::string &Ext)
{
  std::string type{Type};
  type.append( 1, '.' );
  type.append( std::to_string( Seed ) );
  return tokenReplaceCopy( Base, "type", type ) + '.' + Ext;
}

// Make summary files of this data set
void MakeSummaries(const Latan::DMatSample &out, const std::string & sOutFileBase, Latan::SeedType Seed, int momentum_squared)
{
  const int           nt{static_cast<int>(out[Latan::central].rows())};
  const Latan::Index  nSample{out.size()};
  // Now make summary data files
  std::vector<double> TmpAvg( nt );
  std::vector<std::vector<double>> TmpSamples( nt );
  for(int t = 0; t < nt; t++)
    TmpSamples[t].resize( nSample );
  static const char sep[] = " ";
  static const char * SummaryNames[] = { "corr", "mass", "cosh", "sinh" };
  static const char * SummaryHeader[] =
  {
    "# correlator\nt corr_re corr_stddev corr_im corr_im_stddev",
    "# "   "mass\nt mass mass_low mass_high check",
    "# cosh_mass\nt mass mass_low mass_high check",
    "# sinh_mass\nt mass mass_low mass_high check",
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
              case 3: // sinh mass
                DThis = std::asinh((out[i]((t - 1 + nt) % nt, 0) - out[i]((t + 1) % nt, 0)) / (2 * out[i](t, 0)));
                break;
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

struct Momentum {
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
}

// Make a default manifest
// Heavy-light meson decays. 2pt function with current inserted and momentum on light meson

bool ReadCorrelatorsZ2( Common::Manifest & manifest, const std::string &Heavy, const std::string &Light )
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
  static const Common::Momentum Momenta[] = {
    { 0, 0, 0 },
    { 1, 0, 0 },
    { 1, 1, 0 },
    { 1, 1, 1 },
    { 2, 0, 0 },
  };
  static constexpr int NumMomenta{ sizeof( Momenta ) / sizeof( Momenta[0] ) };
  static const std::string sMom0{ Sep + "p" + Sep + Common::Momentum(0,0,0).to_string( Sep ) };
  using Algebra = Grid::Gamma::Algebra;
  static const std::vector<Algebra> alg = { Algebra::Gamma5, Algebra::GammaTGamma5 };
  static const std::vector<std::string> algNames = { "g5", "gT5" };
  assert( alg.size() == algNames.size() && "Abbreviations should match gamma algebra" );
  static const int NumAlgebra{ static_cast<int>( alg.size() ) };
  static const int NumCorr{ NumAlgebra * NumAlgebra };
  std::vector<std::vector<std::complex<double>>> buffer( NumCorr );
  std::vector<std::string> CorrSuffixes( NumCorr );
  std::vector<std::string> CorrNames( NumCorr );
  std::vector<std::string> CorrNamesT( NumCorr );
  std::vector<int>         CorrIndex( NumCorr );
  std::vector<int>         CorrIndexT( NumCorr );
  const unsigned int NumEnds{ Common::EqualIgnoreCase( Heavy, Light ) ? 1u : 2u };
  for( int iSrc = 0; iSrc < NumAlgebra; iSrc++ )
    for( int iSnk = 0; iSnk < NumAlgebra; iSnk++ ) {
      const int iCorr{ iSnk * NumAlgebra + iSrc };
      std::stringstream ss;
      ss << Sep << algNames[iSnk] << Sep << algNames[iSrc];
      CorrSuffixes[iCorr] = ss.str();
      buffer[iCorr].resize( Nt );
    }
  bool bOK{ true };
  for( int iMom = 0; iMom < NumMomenta; iMom++ )
  {
    const std::string sMom{ Sep + "p" + Sep + Momenta[iMom].to_string( Sep )};
    const std::string sMomNeg{ Sep + "p" + Sep + Momenta[iMom].to_string( Sep, true )};
    for( int iCorr = 0 ; iCorr < NumCorr; iCorr++ ) {
      CorrIndex[iCorr] = 0;
      CorrNames[iCorr] = Heavy + Sep + Light + sMom + CorrSuffixes[iCorr];
      manifest.emplace( std::make_pair( CorrNames[iCorr], Common::TrajList( CorrNames[iCorr], iMom, NumEnds * CCount * Nt) ) );
    }
    for( int t = 0; t < Nt; ++t ) {
      const std::string tPrefix{ Sep + "t" + std::to_string( t ) };
      for( int iCorr = 0 ; iCorr < NumCorr; iCorr++ ) {
        CorrIndexT[iCorr] = 0;
        CorrNamesT[iCorr] = CorrNames[iCorr] + tPrefix;
        manifest.emplace( std::make_pair( CorrNamesT[iCorr], Common::TrajList( CorrNamesT[iCorr], iMom, NumEnds * CCount) ) );
      }
      for( unsigned int iEnd = 0; iEnd < NumEnds; iEnd++ ) {
        const std::string & Left{ iEnd == 0 ? Heavy : Light };
        const std::string & Right{ iEnd == 0 ? Light : Heavy };
        for( int iConfig = CStart; iConfig < CEnd; iConfig += CSkip ) {
          const std::string sFileName{ sPath + sProp + Left + tPrefix + sMom + Sep + sProp + Right + tPrefix + sMom0 + Sink + sMomNeg + Dot + std::to_string( iConfig ) + H5};
          try {
            Common::ReadComplexArray(buffer, alg, sFileName, 0);
            for( int iCorr = 0; iCorr < NumCorr; ++iCorr )
            {
              {
                // Put this in the correlator for this timeslice
                Latan::DMat & m{ manifest.at( CorrNamesT[iCorr] ).data[CorrIndexT[iCorr]++] };
                m.resize( Nt, 2 );
                for( unsigned int i = 0; i < Nt; i++ ) {
                  std::complex<double> &d{ buffer[iCorr][i] };
                  m( i, 0 ) = d.real();
                  m( i, 1 ) = d.imag();
                }
              }
              {
                // Put this in the correlator for all timeslices, adjusted for time
                Latan::DMat & m{ manifest.at( CorrNames[iCorr] ).data[CorrIndex[iCorr]++] };
                m.resize( Nt, 2 );
                for( unsigned int i = 0; i < Nt; i++ ) {
                  std::complex<double> &d{ buffer[iCorr][(i + t) % Nt] };
                  m( i, 0 ) = d.real();
                  m( i, 1 ) = d.imag();
                }
              }
            }
          }
          catch(std::runtime_error &e) {
            std::cerr << "\tError " << e.what() << " while reading " << sFileName << std::endl;
            bOK = false;
          } catch(...) {
            std::cerr << "\tError reading " << sFileName << std::endl;
            bOK = false;
          }
        }
      }
    }
  }
  return bOK;
}
