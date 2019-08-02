#include <Grid/Grid.h>
#include <LatAnalyze/Core/OptParser.hpp>
#include <LatAnalyze/Functional/CompiledModel.hpp>
#include <LatAnalyze/Io/Io.hpp>
#include <LatAnalyze/Statistics/MatSample.hpp>
#include <LatAnalyze/Core/Math.hpp>
#include <LatAnalyze/Core/Plot.hpp>
#include <LatAnalyze/Statistics/Dataset.hpp>
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>


#include <string>

using namespace std;
using namespace Latan;
using namespace Grid;
using namespace H5;
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

// void saveCorrelator(const WeakMesonDecay &result, const std::string outFilename, const std::string groupName
void saveCorrelator(const WeakMesonDecay &result, const std::string outFilename, const std::string groupName
                    )
{
    std::string filename;
    std::cout << "Saving correlator to '" << outFilename << "'" << std::endl;
    // makeFileDir(dir);
    ResultWriter writer(outFilename);
    write(writer, groupName, result);
}


int main(int argc, char *argv[])
{
    // parse arguments /////////////////////////////////////////////////////////
    OptParser            opt;
    bool                 parsed;
    string               manFilename, outFilename;
    unsigned int         traj,dt;

   
    opt.addOption("" , "traj"       , OptParser::OptType::value  , false,
                  "trajectory number", "1000");
    opt.addOption("" , "dt"       , OptParser::OptType::value  , false,
                  "tl-tH", "");
    opt.addOption("", "help"      , OptParser::OptType::trigger, true,
                  "show this help message and exit");
    opt.addOption("o", "output", OptParser::OptType::value  , true,
                  "output directory", ".");
    parsed = opt.parse(argc, argv);
    if (!parsed or (opt.getArgs().size() != 1) or opt.gotOption("help"))
    {
        cerr << "usage: " << argv[0] << " <options> <manifest file> " << endl;
        cerr << endl << "Possible options:" << endl << opt << endl;
        // cerr << endl << "output stem must contain '@t_pi@', "
        //  << "e.g. 'foo_@traj@.h5'" << endl;

        return EXIT_FAILURE;
    }  
    
    manFilename = opt.getArgs()[0];
    traj        = opt.optionValue<int>("traj");
    dt          = opt.optionValue<int>("dt");
    outFilename = opt.optionValue<string>("o");



    // load data /////////////////////////////////////////////////////////////////
    cout << "-- loading data...." << endl;

    vector<string>    inFilename;
    Result            buf,sm;
    unsigned int      nFile, nt;

    // Reads the manifest file, ie the /path/to/.h5_file
    inFilename = readManifest(manFilename);
    readFile<Hdf5Reader>(buf, "weakdecay", inFilename[0]);// this reads the DataSet name of .h5 file into buf's vector<ComplexD>; (for resizing variables)
    nFile = inFilename.size();
    nt    = buf.corr.size();
    sm.corr.resize(nt);

    const SpinMatrixD N(1./nFile);
    //const iScalar<Grid::iMatrix<Grid::iScalar<std::__1::complex<double> >, 4> > N(1./nFile); // this is the way to implement multiplication by integers/doubles with matrix
    
    for(int t=0; t<nt; t++)
    {
        for(int i=0; i< nFile; i++)
        {
            readFile<Hdf5Reader>(buf, "weakdecay", inFilename[i]);
            sm.corr[t] += buf.corr[(t+i)%nt];
        }
        
        sm.corr[t]*=N;
        
    }

    auto savename = outFilename + "/dt" + std::to_string(dt) + "_" + std::to_string(traj) + ".h5";
    saveCorrelator(sm, savename, "weakdecay");

    return EXIT_SUCCESS;

}

