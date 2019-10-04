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
void MakeManifest( Common::Manifest & manifest, const std::string &Heavy, const std::string &Light );

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

int main(int argc, char *argv[])
{
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
    const std::string qHeavy{ "h1" };
    const std::string qLight{ "l" };
    std::cout << "No files specified. Making default manifest for "
              << qHeavy << " and " << qLight << std::endl;
    MakeManifest( Manifest, qHeavy, qLight );
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
      const unsigned int nFile{ static_cast<unsigned int>(l.TrajFile.size())};
      Latan::Dataset<Latan::DMat> data(nFile);
      unsigned int nt = 0;
      unsigned int j = 0; // which trajectory are we processing
      std::vector<std::complex<double>> buf;
      bool bError = false;
      bool bResized = false;
      for( auto it = l.TrajFile.begin(); it != l.TrajFile.end(); it++, j++ )
      {
        const int &traj{it->first};
        //std::string &Filename{it->second};
        const Common::TrajFile &tf{it->second};
        const std::string &Filename{tf.Filename};
        std::cout << "\t" << traj << "\t" << Filename << "\t" << tf.offset << std::endl;
        bool bThisLoad = true;
        try {
          Common::ReadComplexArray(buf, Filename);//, Contraction);
        } catch(...) {
          std::cerr << "\tError reading " << Filename << std::endl;
          bThisLoad = false;
          bError = true;
        }
        if( bThisLoad )
        {
          if( !bResized )
          {
            nt = static_cast<unsigned int>(buf.size());
            out.resizeMat(nt, 2);
            bResized = true;
          }
          Latan::DMat & m{data[j]};
          m.resize(nt, 2); // LatAnalyze prefers nt * 2 real matrix
          bool bNumbersOK = true;
          for (unsigned int t = 0; t < nt; ++t)
          {
            const std::complex<double> &d{ buf[( t + tf.offset ) % nt] };
            m(t, 0) = ( tf.bImaginary ? d.imag() : d.real() );
            if( tf.iMultiplier )
              m(t, 0) *= tf.iMultiplier;
            m(t, 1) = ( tf.bImaginary ? d.real() : d.real() );
            if( !std::isfinite( m(t, 0) ) || !std::isfinite( m(t, 1) ) )
              bNumbersOK = false;
          }
          if( !bNumbersOK ) { // Only print this message once per file
            std::cerr << "\t\tError: Infinite/NaN values in " << Filename << std::endl;
            bError = true;
          }
        }
      }
      if( bError )
        std::cerr << "\tOutput not created for " << Contraction << std::endl;
      else
      {
        if( bShowAverages )
          ShowTimeSliceAvg(data);
        // Resample
        std::cout << "-- resampling (" << nSample << " samples)..." << std::endl;
        if( binSize != 1 )
          data.bin(binSize);
        Latan::DMatSample out = data.bootstrapMean(nSample, seed);
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
            const Common::TrajFile &tf{it->second};
            const std::string &Filename{tf.Filename};
            file << "#   traj / file: " << traj << " " << Filename << std::endl;
          }
          data.dumpBootstrapSeq(file, nSample, seed);
          std::cout << std::endl;
        }
        BootstrapCount++;
      }
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

void MakeManifest( Common::Manifest & manifest, const std::string &qHeavy, const std::string &qLight )
{
  const std::string sep{"_"};
  const std::string sepBig{"."};
  static const int Nt{64};
  for( int p2 = 0; p2 < 5; p2++ ) {
    // All I need are the number of momenta >=0 that match p2
    std::vector<Momentum> Momenta;
    for( int px = -p2; px <= p2; px++ )
      for( int py = -p2; py <= p2; py++ )
        for( int pz = -p2; pz <= p2; pz++ )
          if( px*px + py*py +pz*pz == p2 ) {
            Momenta.push_back( Momentum( px, py, pz ) );
          }
    const unsigned int Np{ static_cast<unsigned int>( Momenta.size() ) };
    //for( Momentum &m : Momenta )
      //std::cout << "    " << m.x << ", " << m.y << ", " << m.z << std::endl;
    //assert( Np == 0 && "Der" );
    std::cout << "p^2=" << p2 << " => Np=" << Np << std::endl;
    const std::array<std::string, 2> Mesons{ "ps", "ax" };
    for( const std::string & MesonSrc : Mesons ) {
      for( const std::string & MesonSnk : Mesons ) {
        const std::vector<std::string> g0{ "g0" };
        const std::vector<std::string> gi{ "gx", "gy", "gz" };
        const SubStrings s0{ "V0", g0 };
        const SubStrings si{ "Vi", gi };
        const std::array<SubStrings, 2> VectorType{ s0, si };
        for( const SubStrings & vt : VectorType ) {
          for( int DeltaT = 12; DeltaT <= 24; DeltaT += 4 ) {
            // I have a correlator
            std::string corr{ vt.Name };
            corr.append( sep );
            corr.append( std::to_string( DeltaT ) );
            corr.append( sep );
            corr.append( "p" );
            corr.append( std::to_string( p2 ) );
            corr.append( sep );
            corr.append( MesonSnk );
            corr.append( sep );
            corr.append( MesonSrc );
            std:: cout << "    " << corr << std::endl;
            // Add this to my list of correlators, and get a reference to the trajectory list
            Common::TrajList & cl{manifest[corr]};
            cl.Name = corr;
            cl.momentum = p2;
            unsigned int traj{0};
            for( int ContractFor = 0; ContractFor < 2; ContractFor++ ) {
              const std::string & qLeft{ ContractFor ? qHeavy : qLight };
              const std::string & qRight{ ContractFor ? qLight : qHeavy };
              for( int p=0; p < Np; p++ ) {
                for( int tsrc=0; tsrc < 64; tsrc += 4 ) {
                  for( const std::string & current : vt.Strings ) {
                    const char cDim{ current[1] };
                    if( cDim == '0' || p2 == 0
                       || ( cDim == 'x' && Momenta[p].x != 0 )
                       || ( cDim == 'y' && Momenta[p].y != 0 )
                       || ( cDim == 'z' && Momenta[p].z != 0 ) )
                    {
                      // Create file name
                      std::string f{ sep };
                      f.append( qLeft );
                      f.append( sep );
                      f.append( qRight );
                      f.append( sep );
                      f.append( "p2" );
                      f.append( sep );
                      f.append( std::to_string( p2 ) );
                      f.append( "/sink" );
                      f.append( MesonSuffix( MesonSnk ) );
                      f.append( sep );
                      switch( ContractFor )
                      {
                        case 0:
                          f.append( "0" );
                          break;
                        case 1:
                          f.append( "p" );
                          f.append( sep );
                          f.append( std::to_string( p ) );
                          break;
                        default:
                          assert( 0 && "Unsupported contraction type" );
                      }
                      f.append( sep );
                      f.append( std::to_string( DeltaT ) );
                      f.append( sep );
                      f.append( qLight );
                      f.append( sep );
                      f.append( std::to_string( tsrc ) );
                      f.append( sepBig );
                      f.append( std::to_string( ( tsrc + DeltaT ) % Nt ) );
                      f.append( sepBig );
                      f.append( "source" );
                      f.append( MesonSuffix( MesonSrc ) );
                      f.append( sep );
                      switch( ContractFor )
                      {
                        case 0:
                          f.append( "p" );
                          f.append( sep );
                          f.append( std::to_string( p ) );
                          break;
                        case 1:
                          f.append( "0" );
                          break;
                        default:
                          assert( 0 && "Unsupported contraction type" );
                      }
                      f.append( sep );
                      f.append( std::to_string( tsrc ) );
                      f.append( sepBig );
                      f.append( std::to_string( tsrc ) );
                      f.append( sepBig );
                      f.append( current );
                      f.append( sep );
                      f.append( std::to_string( DeltaT ) );
                      f.append( sep );
                      f.append( qLeft );
                      f.append( sep );
                      f.append( qRight );
                      f.append( sep );
                      f.append( std::to_string( tsrc ) );
                      f.append( sep );
                      f.append( "p" );
                      f.append( sep );
                      f.append( std::to_string( p ) );
                      f.append( ".3200" );
                      f.append( ".h5" );
                      // For momenta <> 2
                      bool bFound{ false };
                      for( int n = 0; n < 2 && !bFound; n++ )
                      {
                        std::string FileName{ "../contract/" };
                        if( n == 0 ) {
                          switch( ContractFor )
                          {
                            case 0:
                              FileName.append( "c30p" );
                              break;
                            case 1:
                              FileName.append( "c3p0" );
                              break;
                            case 2:
                              FileName.append( "c2p0" );
                              break;
                            default:
                              assert( 0 && "Unsupported contraction type" );
                          }
                        } else {
                          FileName.append( "c30p" );
                          // Try it in
                        }
                        FileName.append( f );
                        bFound = Common::FileExists( FileName );
                        if( bFound )
                          f = std::move( FileName );
                        else
                          std::cout << FileName << " not found" << std::endl;
                      }
                      assert( bFound && "File doesn't exist" );
                      Common::TrajFile & tf{ cl.TrajFile[traj++] };
                      tf.Filename = f;
                      tf.offset = tsrc;
                      tf.bImaginary = ( cDim != '0' );
                      tf.iMultiplier = ( ( cDim == 'x' && Momenta[p].x < 0 )
                                        || ( cDim == 'y' && Momenta[p].y < 0 )
                                        || ( cDim == 'z' && Momenta[p].z < 0 ) ) ? -1 : 0;
                      ;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
