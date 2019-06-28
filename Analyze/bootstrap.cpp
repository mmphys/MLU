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

#include <cmath>
#include <iomanip>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <sys/stat.h>
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Io/Hdf5File.hpp>
#include <H5CompType.h>

// Default number of bootstrap replicas
#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif
// Default output file format
#ifndef DEF_FMT
#define DEF_FMT "h5"
#endif

// This describes one contraction and each of its trajectory files
struct OneContraction
{
  std::string Name;                     // name of the contraction
  std::map<int, std::string> Filename;  // list of all the files, sorted by trajectory number
};

// This is a list of all the contractions we've been asked to process
using ContractList = std::map<std::string, OneContraction>;

// Return the same HDF5 complex type Grid uses
namespace H5CustomTypes { static H5::CompType & Complex(void); };
// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
        const std::string &GroupName, const std::string &ObjectName = std::string( "correlator" ) );
// Make a copy of str, in which token has been replaced by x
template <typename T> std::string tokenReplaceCopy(const std::string &str, const std::string token, const T &x );
// Show me the averages for each timeslice
void ShowTimeSliceAvg(const Latan::Dataset<Latan::DMat> &data);
// Does the specified file exist?
bool FileExists(const std::string& Filename);
// Make summary files of this data set
void MakeSummaries(const Latan::DMatSample &out, const std::string & sOutFileBase);
// Process list of files on the command-line, breaking them up into individual trajectories
bool ParseFileList(ContractList &Contractions, const std::vector<std::string> &Args, const std::string &sIgnore);

/*****************************************************************

 Perform a bootstrap of all the files specified on the command line.
 Organise all the files by contraction, and sort them by trajectory.
 Ensure trajectories are processed in the same order so the
 bootstrap replicas for each contraction can be combined later during analysis

*****************************************************************/

int main(int argc, char *argv[])
{
  // parse command line/options ////////////////////////////////////////////////////////
  Latan::OptParser       opt;
  static const char      DefaultOutputStem[] = "@corr@.@type@";
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
                "save text versions of correlators (for use by GnuPlot)");
  opt.addOption("a", "average"   , Latan::OptParser::OptType::trigger, true,
                "show correlator averages");
  opt.addOption("x", "exclude"   , Latan::OptParser::OptType::value,   true,
                "exclude files (colon separated list)", "");
  opt.addOption("" , "help"      , Latan::OptParser::OptType::trigger, true,
                "show this help message and exit");
  bool parsed = opt.parse(argc, argv) and !opt.gotOption("help") and opt.getArgs().size() > 0;
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
  // Now parse the input file names, grouping by correlator, indexed by trajectory
  ContractList Contractions;
  if( parsed )
  {
    const std::string sIgnore{opt.optionValue("x")};
    parsed = ParseFileList(Contractions, opt.getArgs(), sIgnore);
  }
  if (!parsed)
  {
    std::cerr << "usage: " << argv[0] << " <options> ContractionFile1 [ContractionFile2 ...]"
    << "\nOutput stem (if present) must contain '@corr@' and '@type', e.g. '" << DefaultOutputStem << "'"
    << "\nPossible options:\n" << opt << std::endl;
    return EXIT_FAILURE;
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
  const bool bSaveSummaries{ opt.gotOption("t") };
  std::cout << "Creating bootstrap output " << outStem << ", seed=" << seed << std::endl;

  // Walk the list of contractions, performing a separate bootstrap for each
  int BootstrapCount = 0;
  Latan::DMatSample out(nSample);
  for( auto itc = Contractions.begin(); itc != Contractions.end(); itc++ )
  {
    static const std::string szBootstrap{"bootstrap"};
    const std::string &Contraction{itc->first};
    OneContraction &l{itc->second};
    const std::string sOutFileBase{tokenReplaceCopy(outStem, "corr", Contraction)};
    std::string sOutFileName{tokenReplaceCopy(sOutFileBase, "type", szBootstrap + '.' + std::to_string( seed )) + '.' + OutFileExt};
    if( FileExists( sOutFileName ) )
      std::cout << "Skipping " << Contraction << " because " << sOutFileName << " already exists" << std::endl;
    else
    {
      // load all the files for this correlator
      std::cout << Contraction << std::endl;
      const unsigned int nFile{ static_cast<unsigned int>(l.Filename.size())};
      Latan::Dataset<Latan::DMat> data(nFile);
      unsigned int nt = 0;
      unsigned int j = 0; // which trajectory are we processing
      std::vector<std::complex<double>> buf;
      bool bLoadOK = true;
      for( auto it = l.Filename.begin(); bLoadOK && it != l.Filename.end(); it++, j++ )
      {
        const int &traj{it->first};
        std::string &Filename{it->second};
        std::cout << "\t" << traj << "\t" << Filename << std::endl;
        try {
          ReadComplexArray(buf, Filename, Contraction);
        } catch(...) {
          bLoadOK = false;
        }
        if( !bLoadOK )
          std::cerr << "\tError reading " << Filename << std::endl;
        else {
          if(j == 0)
          {
            nt = static_cast<unsigned int>(buf.size());
            out.resizeMat(nt, 2);
          }
          Latan::DMat & m{data[j]};
          m.resize(nt, 2); // LatAnalyze prefers nt * 2 real matrix
          for (unsigned int t = 0; t < nt; ++t)
          {
            m(t, 0) = buf[t].real();
            m(t, 1) = buf[t].imag();
            if( !std::isfinite( m(t, 0) ) || !std::isfinite( m(t, 1) ) )
              bLoadOK = false;
          }
          if(!bLoadOK )
            std::cerr << "\tError: Infinite/NaN values while reading " << Filename << std::endl;
        }
      }
      if(!bLoadOK )
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
        Latan::Io::save<Latan::DMatSample>(out, sOutFileName, Latan::File::Mode::write, szBootstrap);
        std::cout << " done" << std::endl;
        if( bSaveSummaries )
          MakeSummaries(out, sOutFileBase);
        if( dumpBoot ) { // Dump the bootstrap sequences
          std::string SeqFileName{tokenReplaceCopy(sOutFileBase, "type", std::string("bootseq") + '.' + std::to_string( seed )) + ".txt"};
          std::cout << "Saving bootstrap sequence to '" << SeqFileName << "' ...";
          std::ofstream file(SeqFileName);
          file << "# bootstrap sequences" << std::endl;
          file << "#      bin size: " << binSize << std::endl;
          file << "#          seed: " << std::to_string( seed ) << std::endl;
          for( auto it = l.Filename.begin(); it != l.Filename.end(); it++ )
          {
            const int &traj{it->first};
            std::string &Filename{it->second};
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
            << " / " << Contractions.size() << " correlators" << std::endl;
  return EXIT_SUCCESS;
}

// Return the same HDF5 complex type Grid uses
namespace H5CustomTypes {
  static H5::CompType & Complex(void)
  {
    static bool bInitialised = false;
    static H5::CompType m_Complex(sizeof(std::complex<double>));
    {
#ifndef __INTEL_COMPILER
      // mutex not available for this compiler, so initialisation not thread-safe
      static std::mutex sync;
      std::lock_guard<std::mutex> guard( sync );
#endif
      if( !bInitialised ) {
        m_Complex.insertMember("re", 0 * sizeof(double), H5::PredType::NATIVE_DOUBLE);
        m_Complex.insertMember("im", 1 * sizeof(double), H5::PredType::NATIVE_DOUBLE);
        bInitialised = true;
      }
    }
    return m_Complex;
  }
};

// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      const std::string &GroupName, const std::string &ObjectName )
{
  H5::H5File f(FileName, H5F_ACC_RDONLY);
  H5::Group g = f.openGroup(GroupName);
  H5::DataSet ds = g.openDataSet(ObjectName);
  H5::DataSpace dsp = ds.getSpace();
  const int nDims{dsp.getSimpleExtentNdims()};
  if( nDims != 1 ) {
    std::cerr << "Error: " << FileName << ", " << GroupName << ", " << ObjectName << " has " << nDims << " dimensions" << std::endl;
    assert(0 && "Wrong number of dimensions");
  } else {
    hsize_t dim[1];
    dsp.getSimpleExtentDims(dim);
    buffer.resize(dim[0]);
    ds.read(&buffer[0], H5CustomTypes::Complex());
  }
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

// Does the specified file exist?
bool FileExists(const std::string& Filename)
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

// Make summary files of this data set
void MakeSummaries(const Latan::DMatSample &out, const std::string & sOutFileBase)
{
  const int           nt{static_cast<int>(out[Latan::central].rows())};
  const Latan::Index  nSample{out.size()};
  // Now make summary data files
  std::vector<double> TmpAvg( nt );
  std::vector<std::vector<double>> TmpSamples( nt );
  for(int t = 0; t < nt; t++)
    TmpSamples[t].resize( nSample );
  static const char sep[] = " ";
  static const char * SummaryNames[] = { "corr", "mass3pt", "mass" };
  static const char * SummaryHeader[] =
  {
    "# correlator\nt corr_re corr_stddev corr_im corr_im_stddev",
    "# three-point mass\nt mass mass_low mass_high",
    "# mass\nt mass mass_low mass_high",
  };
  for(int f = 0; f < sizeof(SummaryNames)/sizeof(SummaryNames[0]); f++)
  {
    std::string sOutFileName{ tokenReplaceCopy(sOutFileBase, "type", SummaryNames[f]) + ".dat" };
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
        for(int i = 0; i < nSample; i++)
        {
          {
            double DThis;
            switch(f)
            {
              case 1: // three-point mass
                DThis = ( out[i]((t + nt - 1) % nt, 0) + out[i]((t + 1) % nt, 0) ) / (2 * out[i](t, 0) );
                break;
              case 2: // mass
                DThis = - log( out[i]((t + nt + 1) % nt, 0) / out[i](t, 0) );
                break;
              default:
                DThis = 0;
            }
            TmpSamples[t][i] = DThis;
            TmpAvg[t] += DThis;
          }
        }
        TmpAvg[t] /= nSample;
        std::sort(TmpSamples[t].begin(), TmpSamples[t].end());
        int Index = static_cast<int>( 0.16 * nSample + 0.5 );
        s << t << sep << TmpAvg[t] << sep << (TmpAvg[t] - TmpSamples[t][Index])
        << sep << (TmpSamples[t][nSample - Index] - TmpAvg[t]) << std::endl;
      }
    }
  }
}

enum ExtractFilenameReturn {Good, Bad, No_trajectory};

// Extract the contraction name and trajectory number from filename
ExtractFilenameReturn ExtractFilenameParts(const std::string &Filename, std::string &Contraction, int &traj)
{
  ExtractFilenameReturn r{Bad};
  Contraction.clear();
  traj = 0;
  auto pos = Filename.find_last_of('/');
  if( pos == std::string::npos )
    pos = 0;
  else
    pos++;
  auto delim = Filename.find_first_of('.', pos);
  if( delim != std::string::npos )
  {
    Contraction = Filename.substr(pos, delim - pos);
    pos = delim + 1;
    delim = Filename.find_first_of('.', pos);
    if( delim != std::string::npos && delim != pos )
    {
      r = Good;
      while( r == Good && pos < delim )
      {
        auto c = Filename[pos++] - '0';
        if( c >=0 && c < 10 )
          traj = traj * 10 + c;
        else
          r = No_trajectory;
      }
    }
  }
  return r;
}

bool ParseFileList(ContractList &Contractions, const std::vector<std::string> &Args, const std::string &sIgnore)
{
  std::vector<std::string> Ignore;
  for( std::size_t start = 0; start < sIgnore.length() ; ) {
    auto end = sIgnore.find(':', start);
    std::size_t SepLen;
    if( end == std::string::npos ) {
      SepLen = 0;
      end = sIgnore.length();
    }
    else
      SepLen = 1;
    if( end > start )
      Ignore.push_back( sIgnore.substr( start, end - start ) );
    start = end + SepLen;
  }
  bool parsed = ( Args.size() > 0 );
  for( const std::string &Filename : Args )
  {
    // See whether this file is in the ignore list
    std::size_t iIgnore = 0;
    while( iIgnore < Ignore.size() && Ignore[iIgnore].compare(Filename) )
      iIgnore++;
    if( iIgnore < Ignore.size() )
      std::cout << "Ignoring " << Filename << std::endl;
    else if( !FileExists(Filename))
    {
      parsed = false;
      std::cout << "Error: " << Filename << " doesn't exist" << std::endl;
    }
    else
    {
      std::string Contraction;
      int         traj;
      switch( ExtractFilenameParts( Filename, Contraction, traj ) )
      {
        case Good:
        {
          OneContraction & cl{Contractions[Contraction]};
          if( cl.Filename.size() == 0 )
            cl.Name = Contraction;
          auto it = cl.Filename.find( traj );
          if( it == cl.Filename.end() )
            cl.Filename[traj] = Filename;
          else
          {
            if( !Filename.compare(it->second) )
              std::cout << "Ignoring repetition of " << Filename << std::endl;
            else
            {
              parsed = false;
              std::cout << "Error " << Filename << " conflicts with " << it->second << std::endl;
            }
          }
        }
          break;
          
        case No_trajectory:
          std::cout << "Ignoring non-numeric trajectory " << Filename << std::endl;
          break;
          
        default:
          parsed = false;
          std::cout << "Error: " << Filename << " is not a contraction file" << std::endl;
          break;
      }
    }
  }
  return parsed;
}
