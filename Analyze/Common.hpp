#ifndef Common_hpp
#define Common_hpp

#include <complex>
#include <map>
#include <string>
#include <vector>
#include <H5Cpp.h>
#include <LatAnalyze/Statistics/Dataset.hpp>

#define BEGIN_COMMON_NAMESPACE namespace Common {
#define END_COMMON_NAMESPACE   };

BEGIN_COMMON_NAMESPACE

class H5File : public H5::H5File {
public:
  using H5::H5File::H5File; // Inherit the base class's constructors
  // Return the same HDF5 complex type Grid uses
  static H5::CompType & ComplexType( void );
};

// Get first groupname from specified group
std::string GetFirstGroupName( H5::Group & g );

// Does the specified file exist?
bool FileExists(const std::string& Filename);

// Are all the floating point numbers pointed to finite
template <typename T, typename I> inline bool IsFinite( const T * d, I n )
{
  while( n-- )
    if( !std::isfinite( *d++ ) )
      return false;
  return true;
}

// Are all the floating point numbers in this matrix finite
template <typename T> inline bool IsFinite( const Latan::Mat<T> & m, bool bDiagonalsOnly = false )
{
  for( Latan::Index row = 0; row < m.rows(); ++row )
    for( Latan::Index col = 0; col < m.cols(); ++col )
      if( ( !bDiagonalsOnly || row == col ) && !std::isfinite( m( row, col ) ) )
        return false;
  return true;
}

// Are all the floating point numbers in this vector finite
template <typename T> inline bool IsFinite( const std::vector<T> & v )
{
  for( const T n : v )
    if( !std::isfinite( n ) )
      return false;
  return true;
}

// Read a complex array from an HDF5 file

void ReadComplexArray(std::vector<std::complex<double>> &buffer, const std::string &FileName,
                      const std::string &GroupName  = std::string(),
                      const std::string &ObjectName = std::string( "correlator" ) );

// Read a list of bootstrapped correlators into a single correlator
Latan::DMatSample ReadBootstrapCorrs( const std::vector<std::string> & FileName, int Fold, int Shift, int NumOps );

// This describes one contraction and each of its trajectory files
struct TrajList
{
  std::string Name;                     // name of the contraction
  std::map<int, std::string> Filename;  // list of all the files, sorted by trajectory number
};

// This is a list of all the contractions we've been asked to process
class Manifest : public std::map<std::string, TrajList> {
public:
  // Process list of files on the command-line, breaking them up into individual trajectories
  explicit Manifest(const std::vector<std::string> &Files, const std::string &sIgnore = "");
};

END_COMMON_NAMESPACE
#endif // Common_hpp
