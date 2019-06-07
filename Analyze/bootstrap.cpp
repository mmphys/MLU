/*************************************************************************************
 
 Compare perambulators
 
 Source file: bootstrap.cpp
 
 Copyright (C) 2019
 
 Author: Multiple (this version inherited from Andrew Yong)
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

#include <sys/stat.h>
#include <Grid/Grid.h>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
//#include <LatAnalyze/Io/Hdf5File.hpp>

using namespace std;
using namespace Latan;
using namespace Grid;

namespace Contractor
{
  using namespace Grid;
  class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        std::string, name);
    };

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times,
                                        std::string, translations,
                                        bool, translationAverage);
    };

  class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<Contractor::A2AMatrixPar>,  a2aMatrix,
                                        ProductPar, contraction,
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
    };
}

#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif

template <typename Reader, typename Result>
void readFile(Result &out, const std::string &name, const std::string &filename)
{
  Reader reader(filename);
  read(reader, name, out);
}

template <typename T>
std::string tokenReplaceCopy(std::string &str, const std::string token, const T &x )
{
  std::string sCopy{str};
  tokenReplace(sCopy, token, x );
  return sCopy;
}

#ifdef TURNS_OUT_THIS_WASNT_BUGGY
template <typename T>
void NonBuggyPtVectorMean(T &m, const Dataset<T> &ds, const std::vector<int> &selection)
{
  const int NumReplicas{static_cast<int>(ds.size())};
  if( NumReplicas ) {
    const int rows{static_cast<int>(m.rows())};
    const int cols{static_cast<int>(m.cols())};
    //m.resize(rows, cols); // resizing vector to nt * 2 matrix
    for( unsigned int t = 0; t < rows; t++ )
      for( unsigned int j = 0; j < cols; j++ ) {
        m(t,j) = ds[selection[0]](t,j);
        for( unsigned int i = 1; i < NumReplicas; ++i )
          m(t,j) += ds[selection[i]](t,j);
        m(t,j) /= NumReplicas;
      }
  }
}

template <typename T>
void NonBuggyBootstrapMean(Sample<T> &s, const Dataset<T> &ds, const SeedType seed, unsigned int nt)
{
  const Index                         nSample{s.size()};
  const int                           NumReplicas{static_cast<int>(ds.size())};
  std::mt19937                        gen(seed);
  std::uniform_int_distribution<int>  dis(0, NumReplicas - 1);
  std::vector<int>                    selection(NumReplicas);

  std::cout << "Taking average over " << nSample << " samples" << std::endl;
  for (int j = 0; j < NumReplicas; ++j)
    selection[j] = j;
  NonBuggyPtVectorMean(s[central], ds, selection);
  for (Index i = 0; i < nSample; ++i) {
    for (unsigned int j = 0; j < NumReplicas; ++j)
      selection[j] = dis(gen);
    NonBuggyPtVectorMean(s[i], ds, selection);
  }
}
#endif //TURNS_OUT_THIS_WASNT_BUGGY

// Show me the averages for each timeslice
void ShowTimeSliceAvg(const Dataset<DMat> &data) {
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

// Extract the contraction name and trajectory number from filename
enum ExtractFilenameReturn {Good, Bad, No_trajectory};
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

// Per contraction map from trajectory numbers to filename
struct ContractionList
{
  std::string Contraction;
  std::map<int, std::string> Filename;
};

// Map from contraction names to ContractionList
using ContractList = std::map<std::string, ContractionList>;

bool FileExists(const std::string& Filename)
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

int main(int argc, char *argv[])
{
  // parse command line/options ////////////////////////////////////////////////////////
  OptParser              opt;
  string                 outStem;
  Latan::Index           binSize, nSample;
  random_device          rd;
  SeedType               seed = rd();
  static const char      DefaultOutputStem[] = "@corr@.bootstrap.h5";
  opt.addOption("n", "nsample"   , OptParser::OptType::value,   true,
                "number of samples", DEF_NSAMPLE);
  opt.addOption("b", "bin"       , OptParser::OptType::value,   true,
                "bin size", "1");
  opt.addOption("o", "output"    , OptParser::OptType::value,   true,
                "output file", DefaultOutputStem);
  opt.addOption("r", "seed"      , OptParser::OptType::value,   true,
                "random generator seed (default: random)");
  opt.addOption("v", "verbose"   , OptParser::OptType::trigger, true,
                "show extra detail, e.g. Correlator averages");
  opt.addOption("" , "help"      , OptParser::OptType::trigger, true,
                "show this help message and exit");
  bool parsed = opt.parse(argc, argv) and !opt.gotOption("help") and opt.getArgs().size() > 0;
  if (parsed)
  {
    // If the output stem was specified, ensure it includes "@corr@"
    outStem = opt.optionValue("o");
    if (outStem.find("@corr@") == std::string::npos)
    {
      parsed=false;
      std::cout << "Output stem " << outStem << " invalid" << std::endl;
    }
  }
  // Now parse the input file names, grouping by correlator, indexed by trajectory
  ContractList Contractions;
  if( parsed )
  {
    for( int i = 0; i < opt.getArgs().size(); i++ )
    {
      const std::string &Filename{opt.getArgs()[i]};
      if( !FileExists(Filename))
      {
        parsed = false;
        std::cout << "Error: File " << (i+1) << " nonexistent - " << Filename << std::endl;
      }
      else
      {
        std::string Contraction;
        int         traj;
        switch( ExtractFilenameParts( Filename, Contraction, traj ) )
        {
          case Good:
            {
              ContractionList & cl{Contractions[Contraction]};
              if( cl.Filename.size() == 0 )
                cl.Contraction = Contraction;
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
            std::cout << "Error: File " << (i+1) << " not a contraction file - " << Filename << std::endl;
            break;
        }
      }
    }
  }
  if (!parsed)
  {
    cerr << "usage: " << argv[0] << " <options> ContractionFile1 [ContractionFile2 ...]"
    << "\nOutput stem (if present) must contain '@corr@', e.g. '" << DefaultOutputStem << "'"
    << "\nPossible options:\n" << opt << endl;
    return EXIT_FAILURE;
  }
  
  nSample     = opt.optionValue<Latan::Index>("n");
  binSize     = opt.optionValue<Latan::Index>("b");
  if (opt.gotOption("r"))
  {
    seed = opt.optionValue<SeedType>("r");
  }

  // load data /////////////////////////////////////////////////////////////////
  cout << "Creating bootstrap output " << outStem << endl;

  // Walk the list of contractions, performing a separate bootstrap for each
  int BootstrapCount = 0;
  Latan::DMatSample out(nSample);
  for( auto itc = Contractions.begin(); itc != Contractions.end(); itc++ )
  {
    const std::string &Contraction{itc->first};
    ContractionList &l{itc->second};
    std::string sOutFileName{tokenReplaceCopy(outStem, "corr", Contraction)};
    if( FileExists( sOutFileName ) )
      std::cout << "Skipping " << Contraction << " because " << sOutFileName << " already exists" << std::endl;
    else
    {
      std::cout << Contraction << std::endl;
      const unsigned int nFile{ static_cast<unsigned int>(l.Filename.size())};
      Latan::Dataset<Latan::DMat> data(nFile);
      // this reads the DataSet name of .h5 file into buf's vector<ComplexD>; (for resizing variables)
      unsigned int nt = 0;
      unsigned int j = 0; // which file are we processing
      for( auto it = l.Filename.begin(); it != l.Filename.end(); it++, j++ )
      {
        const int &traj{it->first};
        std::string &Filename{it->second};
        std::cout << "\t" << traj << "\t" << Filename << std::endl;
        Contractor::CorrelatorResult buf;
        readFile<Hdf5Reader>(buf, Contraction, Filename);
        if(j == 0)
        {
          nt = static_cast<unsigned int>(buf.correlator.size());
          for( unsigned int i = 0; i < nSample; i++)
            out[i].resize(nt,2);
          out[central].resize(nt,2);
        }
        Latan::DMat & m{data[j]};
        m.resize(nt, 2); // resizing vector to nt * 2 matrix
        for (unsigned int t = 0; t < nt; ++t)
        {
          m(t, 0) = buf.correlator[t].real();
          m(t, 1) = buf.correlator[t].imag();
        }
      }

      if( opt.gotOption("verbose") )
        ShowTimeSliceAvg(data);

      // Resampling //////////////////////////////////////////////////////////////////
      cout << "-- resampling (" << nSample << " samples)..." << endl;
      if( binSize != 1 )
        data.bin(binSize);
#ifdef TURNS_OUT_THIS_WASNT_BUGGY
      NonBuggyBootstrapMean(out, data, seed, nt);
#else
      DMatSample out = data.bootstrapMean(nSample, seed);
#endif //TURNS_OUT_THIS_WASNT_BUGGY
      cout << "Saving sample to '" << sOutFileName << "' ...";
      Latan::Io::save<DMatSample>(out, sOutFileName, Latan::File::Mode::write, "bootstrap");
      cout << " done" << endl;
      BootstrapCount++;
    }
  }
  std::cout << "Bootstraps written for " << BootstrapCount << " correlators" << endl;
  return EXIT_SUCCESS;
}
