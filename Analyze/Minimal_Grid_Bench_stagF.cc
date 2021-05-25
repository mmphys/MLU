    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./benchmarks/Benchmark_staggered.cc

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: paboyle <paboyle@ph.ed.ac.uk>

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

#if 0

#include <Grid/Grid.h>

#else

#define NAMESPACE_BEGIN(x) namespace x {
#define NAMESPACE_END(x) }

#include <sys/signal.h>

#include <complex>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include <Grid/GridCore.h>
#include <Grid/GridQCDcore.h>
#include <Grid/qcd/action/Action.h>
#include <Grid/qcd/action/gauge/Gauge.h>
#include <Grid/qcd/action/fermion/FermionOperatorImpl.h>
#include <Grid/qcd/action/fermion/FermionOperator.h>
#include <Grid/qcd/action/fermion/WilsonKernels.h>
#include <Grid/qcd/action/fermion/Fermion.h>
//#include <Grid/qcd/action/fermion/WilsonImpl.h>
//#include <Grid/qcd/action/pseudofermion/PseudoFermion.h>
//#include <Grid/qcd/QCD.h>
//#include <Grid/qcd/action/fermion/WilsonKernels.h>
//#include <Grid/qcd/action/fermion/Fermion.h>

#include "Config.h"
#include <Grid/threads/Threads.h>
#include <Grid/threads/Accelerator.h>
#include <Grid/util/Coordinate.h>
#include <Grid/util/Init.h>
#include <Grid/perfmon/Timer.h>
#include <Grid/log/Log.h>
#include <Grid/allocator/MemoryManager.h>
#include <Grid/allocator/MemoryStats.h>
#include <Grid/allocator/AlignedAllocator.h>
#include <Grid/simd/Simd.h>
#include <Grid/simd/Grid_vector_types.h>
#include <Grid/communicator/SharedMemory.h>
#include <Grid/communicator/Communicator_base.h>
#include <Grid/cartesian/Cartesian_base.h>
#include <Grid/cartesian/Cartesian_full.h>

//static constexpr int Nd{4};

#endif

using namespace std;
using namespace Grid;

int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  Coordinate latt_size   = GridDefaultLatt();
  Coordinate simd_layout = GridDefaultSimd(Nd,vComplexF::Nsimd());
  Coordinate mpi_layout  = GridDefaultMpi();
  GridCartesian Grid(latt_size,simd_layout,mpi_layout);

  typedef typename ImprovedStaggeredFermionF::FermionField FermionField;
  FermionField result(&Grid); result=Zero();
  double t0=usecond();
  std::cout<<GridLogMessage << "Before dummy norm"<<std::endl;
  auto n2r = norm2(result);
  std::cout<<GridLogMessage << "Dummy norm result "<< n2r <<std::endl;
  std::cout<<GridLogMessage << "Calling Ds"<<std::endl;
  double t1=usecond();

  std::cout<<GridLogMessage << "Before final norm"<<std::endl;
  n2r = norm2(result);
  std::cout<<GridLogMessage << "Final norm result "<< n2r <<std::endl;

  std::cout<<GridLogMessage << "Called Ds"<<std::endl;
  std::cout<<GridLogMessage << "norm result "<< n2r<<std::endl;

  Grid_finalize();
}
