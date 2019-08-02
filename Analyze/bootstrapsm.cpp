#include <Grid/Grid.h>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <Hadrons/Global.hpp>

using namespace std;
using namespace Latan;
using namespace Grid;
using namespace Hadrons;

class WeakMesonDecay: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(WeakMesonDecay,
                                        std::vector<SpinMatrix>, corr);
    };



typedef WeakMesonDecay Result;

template <typename Reader, typename Result>
void readFile(Result &out, const std::string &name, 
              const std::string &filename)
{
  Reader reader(filename);

  read(reader, name, out);
}

#ifndef DEF_NSAMPLE
#define DEF_NSAMPLE "10000"
#endif



int main(int argc, char *argv[])
{

  
  // parse command line/options ////////////////////////////////////////////////////////
  OptParser              opt;
  bool                   parsed;
  string                 manFilename, outStem, groupName;
  Latan::Index           binSize, nSample;
  random_device          rd;
  SeedType               seed = rd();

  opt.addOption("g", "name"   , OptParser::OptType::value,   true,
                "h5 file group name", "");  
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
  groupName   = opt.optionValue("g");

  if (opt.gotOption("r"))
  {
    seed = opt.optionValue<SeedType>("r");
  }
  manFilename = opt.getArgs()[0];
  outStem     = opt.getArgs()[1];

  // load data /////////////////////////////////////////////////////////////////
  cout << "-- loading data..." << endl;
  
  vector<string>         inFilename, outFilename;
  unsigned int           nFile, nt;
  vector<Dataset<DMat>>  data;
  // vector<SpinMatrix>  data;
  Result                 buf;
  ComplexD               result(0,0);

// Reads the manifest file, ie the /path/to/.h5_file
  inFilename = readManifest(manFilename);
  readFile<Hdf5Reader>(buf, groupName, inFilename[0]);// this reads the DataSet name of .h5 file into buf's vector<ComplexD>; (for resizing variables)
  nFile = inFilename.size();
  nt    = buf.corr.size();
  data.resize(1, Latan::Dataset<Latan::DMat>(nFile)); // resizing the data such that each element contains a Latan::DMat

for (unsigned int j = 0; j < nFile; ++j) // for each .h5 file, read real and imag part into data
{
  data[0][j].resize(nt,2);
  for(unsigned int t=0; t < nt; ++t)
  {
    result = ComplexD(0,0);
    std::cout << '\r' << Latan::ProgressBar(j+1, nFile);
    readFile<Hdf5Reader>(buf, groupName, inFilename[j]);
    for(int i=1; i< 4; i++) // iterate through spatial Lorentz index
    {
      result = result + trace(Gamma::gmu[i]*buf.corr[t]*GammaL(Gamma::gmu[0])); // trace needs momentum component
      // cout<< "\nAfter:" << result << result.real() << endl;
    }
    data[0][j](t,0) = result.real();
    data[0][j](t,1) = result.imag();
  }

    
  }


  std::cout << std::endl; 

  // Resampling //////////////////////////////////////////////////////////////////
  cout << "-- resampling (" << nSample << " samples)..." << endl;

  DMatSample out;

  for (unsigned int j = 0; j < data.size(); ++j)
  {
    data[j].bin(binSize);
    out = data[j].bootstrapMean(nSample, seed);
    Io::save<DMatSample>(out, outStem);
    cout << "Sample saved to '" << outStem << "'" << endl;
  }
  
  return EXIT_SUCCESS;
}