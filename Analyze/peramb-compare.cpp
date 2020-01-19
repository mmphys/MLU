/*************************************************************************************
 
 Compare perambulators
 
 Source file: peramb-compare.cpp
 
 Copyright (C) 2019
 
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

#include <stdio.h>
#include <typeinfo>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "Common.hpp"

using namespace Grid;
using namespace Hadrons;

#define DEFAULT_EPSILON "0.00001"

// Very simple iterators for Eigen tensors
// The only way I could get these iterators to work is to put the begin() and end() functions in the Eigen namespace
// So if Eigen ever defines these, we'll have a conflict and have to change this
namespace Eigen {
  template <typename ET>
  inline typename std::enable_if<EigenIO::is_tensor<ET>::value, typename EigenIO::Traits<ET>::scalar_type *>::type
  begin( ET & et ) { return reinterpret_cast<typename Grid::EigenIO::Traits<ET>::scalar_type *>(et.data()); }
  template <typename ET>
  inline typename std::enable_if<EigenIO::is_tensor<ET>::value, typename EigenIO::Traits<ET>::scalar_type *>::type
  end( ET & et ) { return begin(et) + et.size() * EigenIO::Traits<ET>::count; }
}

template<typename T>
inline bool ApproxEq( std::complex<T> x, std::complex<T> y, T epsilon, T &diffMax) {
  T ReX = x.real();
  T MagX = std::abs(x);
  T ReDiff = std::abs( ( ReX - y.real() ) / ReX );
  T MagDiff = std::abs( ( MagX - std::abs(y) ) / MagX );
  if( diffMax < ReDiff )
    diffMax = ReDiff;
  return ReDiff < epsilon && MagDiff < epsilon;
}

int TestPer(const std::string &pszF1, const std::string &pszF2, const double epsilon) {
  std::cout << "Comparing: " << pszF1 << " with " << pszF2 << std::endl;
  using Per = Grid::Hadrons::MDistil::PerambTensor;
  //std::array<std::string,6> sIndexNames{"Nt", "nvec", "LI", "nnoise", "Nt_inv", "SI"};
  //Per per1(sIndexNames,64,100,100,1,1,4);
  //Per per2(sIndexNames,64,100,100,1,1,4);
  Per per1;
  Per per2;
  per1.read( pszF1 );
  per2.read( pszF2 );

  // Cast to eigen tensor
  using PerT = typename Per::ET;
  PerT & pt1{ per1.tensor };
  PerT & pt2{ per2.tensor };

  //Norm of the difference and compare individual members
  PerT diff = pt1 - pt2;
  //Eigen::Tensor<double, 0, PerT::Options> dResult = (pt1 == pt2).all();
  int i = 0;
  double dSum = 0;
  double diffMax = 0;
  Complex * p2 = begin(pt2);
  for(Complex v1 : pt1) {
    dSum += std::abs( v1 );
    if( !ApproxEq( v1, *p2, epsilon, diffMax ) ) {
      std::cout << "C(" << i << ") is different: " << v1 << " vs " << *p2 << std::endl;
      return EXIT_FAILURE;
    }
    i++;
    p2++;
  }
  double dNorm = std::sqrt(dSum);
  std::cout << i << " entries, Sum=" << dSum << ", Norm=" << dNorm << ", Avg=" << dNorm / i << std::endl;
  std::cout << "All " << i << " Complex elements compare within " << diffMax << std::endl;
  return EXIT_SUCCESS;
}

int main(int argc, const char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  using CL = Common::CommandLine;
  CL cl;
  try
  {
    const std::initializer_list<CL::SwitchDef> list = {
      {"e", CL::SwitchType::Single, DEFAULT_EPSILON},
      {"help", CL::SwitchType::Flag, nullptr},
    };
    cl.Parse( argc, argv, list );
    const int NumFiles{ static_cast<int>( cl.Args.size() ) };
    if( !cl.GotSwitch( "help" ) && NumFiles == 2 )
    {
      bShowUsage = false;
      iReturn = TestPer( cl.Args[0], cl.Args[1], cl.SwitchValue<double>("e") );
    }
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: " << cl.Name <<
    " <options> Perambulator1 Perambulator2\n"
    "Compare two perambulators, where <options> are:\n"
    "-e     epsilon (i.e. comparison tolerance), default=" DEFAULT_EPSILON "\n"
    "Flags:\n"
    "--help     This message\n";
  }
  return iReturn;
}
