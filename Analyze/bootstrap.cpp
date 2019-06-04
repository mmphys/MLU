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

// Debug - show me the averages for each timeslice
#ifdef DEBUG
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
#endif //DEBUG

int main(int argc, char *argv[])
{
  // parse command line/options ////////////////////////////////////////////////////////
  OptParser              opt;
  bool                   parsed;
  string                 manFilename, outStem;
  Latan::Index           binSize, nSample;
  random_device          rd;
  SeedType               seed = rd();
  
  opt.addOption("n", "nsample"   , OptParser::OptType::value,   true,
                "number of samples", DEF_NSAMPLE);
  opt.addOption("b", "bin"       , OptParser::OptType::value,   true,
                "bin size", "1");
  opt.addOption("r", "seed"      , OptParser::OptType::value,   true,
                "random generator seed (default: random)");
  opt.addOption("" , "help"      , OptParser::OptType::trigger, true,
                "show this help message and exit");
  parsed = opt.parse(argc, argv);
  if (!parsed or (opt.getArgs().size() != 2) or opt.gotOption("help"))
  {
    cerr << "usage: " << argv[0];
    cerr << " <manifest file> <output stem> <options>" << endl;
    cerr << endl << "output stem must contain '@corr@', "
         << "e.g. 'foo_@corr@.h5'" << endl;
    cerr << endl << "Possible options:" << endl << opt << endl;
    
    return EXIT_FAILURE;
  }
  
  nSample     = opt.optionValue<Latan::Index>("n");
  binSize     = opt.optionValue<Latan::Index>("b");
  if (opt.gotOption("r"))
  {
    seed = opt.optionValue<SeedType>("r");
  }
  manFilename = opt.getArgs()[0];
  outStem     = opt.getArgs()[1];

  // load data /////////////////////////////////////////////////////////////////
  cout << "Manifest file " << manFilename << ", output stem " << outStem
    << "\n-- loading data..." << endl;
  
  // Which correlators are we interested in
  //vector<string> corrNames{"phiL_0_rhoL", "phiL_g5_0_rhoL_g5"};
  vector<string> corrNames{"phiL_0_rhoL"};
  //vector<string> corrNames{"phiL_g5_0_rhoL_g5"};
  const unsigned int nCorrs{ static_cast<unsigned int>(corrNames.size())};

  // Read the manifest file, ie the /path/to/.h5_file
  vector<string> inFilename{};
  {
    bool bComments = false;
    vector<string> actualManifest{readManifest(manFilename)};
    for( string &s : actualManifest ) {
      if( (s[0] == '/' && s[1] == '/') || s[0] == '%' ) {
        if( !bComments ) {
          bComments = true;
          std::cout << "Ignoring comment(s) in manifest:" << std::endl;
        }
        std::cout << s << std::endl;
      }
      else
        inFilename.push_back( s );
    }
  }
  const unsigned int nFile{ static_cast<unsigned int>(inFilename.size())};
  vector<Dataset<DMat>> data(nCorrs, Latan::Dataset<Latan::DMat>(nFile)); // resizing the data such that each element contains a Latan::DMat

  // this reads the DataSet name of .h5 file into buf's vector<ComplexD>; (for resizing variables)
  Contractor::CorrelatorResult  buf;
  readFile<Hdf5Reader>(buf, corrNames[0], tokenReplaceCopy(inFilename[0], "corr", corrNames[0]));
  const unsigned int nt{ static_cast<unsigned int>(buf.correlator.size())};
  DMatSample out(nSample);
  for( unsigned int i = 0; i < nSample; i++)
    out[i].resize(nt,2);
  out[central].resize(nt,2);
  // data = Latan::Dataset<Latan::DMat>(nFile);

  for(unsigned int i=0; i < nCorrs; ++i)
  {
    for (unsigned int j = 0; j < nFile; ++j) // for each .h5 file, read real and imag part into data
    {
      std::cout << '\r' << Latan::ProgressBar(j+1, nFile);
      if( i != 0 || j != 0 )
        // maybe set dataset name as a command line argument?
        readFile<Hdf5Reader>(buf, corrNames[i], tokenReplaceCopy(inFilename[j], "corr", corrNames[i]));
      Latan::DMat & m{data[i][j]};
      m.resize(nt, 2); // resizing vector to nt * 2 matrix
      for (unsigned int t = 0; t < nt; ++t)
      {
        m(t, 0) = buf.correlator[t].real(); // seems like each element of data is a 'nt x 2' matrix.
        m(t, 1) = buf.correlator[t].imag();
      }
    }
    std::cout << std::endl;

    // Debug - show me the averages for each timeslice
#ifdef DEBUG
    ShowTimeSliceAvg(data[i]);
#endif

    // Resampling //////////////////////////////////////////////////////////////////
    if(1) {
      cout << "-- resampling (" << nSample << " samples)..." << endl;
      std::string sOutFileName{tokenReplaceCopy(outStem, "corr", corrNames[i])};
      if( binSize != 1 )
        data[i].bin(binSize);
#ifdef TURNS_OUT_THIS_WASNT_BUGGY
      NonBuggyBootstrapMean(out, data[i], seed, nt);
#else
      DMatSample out = data[i].bootstrapMean(nSample, seed);
#endif //TURNS_OUT_THIS_WASNT_BUGGY
      cout << "Saving sample to '" << sOutFileName << "' ...";
      Io::save<DMatSample>(out, sOutFileName);
      cout << " done" << endl;
    }
  }
  cout << "All correlators written" << endl;
  int * a = new int[3];
  for(int i = 0; i<4; i++)
    a[i] = 7;
  delete [] a;
  return EXIT_SUCCESS;
}
