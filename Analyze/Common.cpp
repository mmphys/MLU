#include "Common.hpp"

#include <iostream>
#include <mutex> // Apparently empty under __INTEL_COMPILER
#include <sys/stat.h>
#include <H5CompType.h>

BEGIN_COMMON_NAMESPACE

// Does the specified file exist?
bool FileExists(const std::string& Filename)
{
  struct stat buf;
  return stat(Filename.c_str(), &buf) != -1;
}

// Return the same HDF5 complex type Grid uses
#ifndef __INTEL_COMPILER
// mutex not available for this compiler, so initialisation not thread-safe
static std::mutex ComplexTypeSync;
#endif

H5::CompType & H5File::ComplexType(void)
{
  static bool bInitialised = false;
  static H5::CompType m_Complex(sizeof(std::complex<double>));
  {
#ifndef __INTEL_COMPILER
    std::lock_guard<std::mutex> guard( ComplexTypeSync );
#endif
    if( !bInitialised ) {
      m_Complex.insertMember("re", 0 * sizeof(double), H5::PredType::NATIVE_DOUBLE);
      m_Complex.insertMember("im", 1 * sizeof(double), H5::PredType::NATIVE_DOUBLE);
      bInitialised = true;
    }
  }
  return m_Complex;
}

// Get first groupname from specified group
std::string GetFirstGroupName( H5::Group & g )
{
  hsize_t n = g.getNumObjs();
  for( hsize_t i = 0; i < n; ++i ) {
    H5G_obj_t t = g.getObjTypeByIdx( i );
    if( t == H5G_GROUP )
      return g.getObjnameByIdx( i );
  }
  return std::string();
}

// Read a complex array from an HDF5 file
void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      const std::string &GroupName, const std::string &ObjectName )
{
  std::cout << "Reading complex array from " << FileName;
  H5::H5File f(FileName, H5F_ACC_RDONLY);
  const bool bFindGroupName{ GroupName.empty() };
  H5::Group g = f.openGroup( bFindGroupName ? std::string("/") : GroupName );
  std::cout << ", group ";
  if( bFindGroupName ) {
    std::string FirstGroupName = GetFirstGroupName( g );
    std::cout << FirstGroupName << "\n";
    g = g.openGroup(FirstGroupName);
  }
  else
    std::cout << GroupName << "\n";
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
    ds.read(&buffer[0], Common::H5File::ComplexType());
  }
}

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
        std::cout << " folded into " << NtDest << " timeslices with " << (Fold == 1 ? "positive" : "negative") << " parity";
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

enum ExtractFilenameReturn {Good, Bad, No_trajectory};

// Extract the contraction name and trajectory number from filename
static ExtractFilenameReturn ExtractFilenameParts(const std::string &Filename, std::string &Contraction, int &traj)
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

Manifest::Manifest(const std::vector<std::string> &Args, const std::string &sIgnore)
{
  std::map<std::string, TrajList> & Contractions = (* this);
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
          TrajList & cl{Contractions[Contraction]};
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
}

END_COMMON_NAMESPACE
