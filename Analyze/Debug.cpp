//
//  Debug.cpp
//  A place to play with debug code
//
//  Created by Michael Marshall on 09/05/2020.
//  Copyright © 2020 sopa. All rights reserved.
//

#include "Debug.hpp"

#ifdef TEST_DEBUG_NEW
static thread_local bool bUseDevice{ false };

struct DebugNewHeader
{
  std::size_t szBufSize;
  bool bOnDevice;
};
static constexpr std::size_t DebugNewHeaderSize = ( sizeof( DebugNewHeader ) + 15 ) & (~15);

std::ostream& operator<<(std::ostream& os, const DebugNewHeader &h)
{
  return os << &h << ": " << h.szBufSize << " bytes of " << ( h.bOnDevice ? "device" : "host" ) << " memory";
}

void * DebugNew( std::size_t size )
{
  DebugNewHeader * pHeader = static_cast<DebugNewHeader *>( std::malloc( size + DebugNewHeaderSize ) );
  pHeader->szBufSize = size + DebugNewHeaderSize;
  pHeader->bOnDevice = bUseDevice;
  std::cout << "DebugNew( " << size << " ) = " << *pHeader << std::endl;
  return reinterpret_cast<char *>( pHeader ) + DebugNewHeaderSize;
}

void DebugDelete( void * mem )
{
  DebugNewHeader * pHeader = reinterpret_cast<DebugNewHeader *>( static_cast<char *>( mem ) - DebugNewHeaderSize );
  std::cout << "DebugDelete( " << mem << " ) = " << *pHeader << std::endl;
  std::free( pHeader );
}

void * operator new( std::size_t size ) { return DebugNew( size ); }
void operator delete( void * mem ) { DebugDelete( mem ); }

class CTest
{
public:
  using Scalar = Grid::ComplexD;
  using ET = Eigen::Tensor<Scalar, 6, Eigen::RowMajor>;
  std::vector<Scalar> v;
  //std::vector<Scalar,Grid::alignedAllocator<Scalar>> vGrid;
  std::vector<Scalar> vGrid;
  ET tensor;
  Eigen::TensorMap<ET> tMap;
  CTest();
  static void * operator new( std::size_t size ) { return DebugNew( size ); }
  static void operator delete( void * mem ) { DebugDelete( mem ); }
};

CTest::CTest() : v(2), vGrid(8*1*2*3*4*5), tensor(8,1,2,3,4,5), tMap(&vGrid[0],8,1,2,3,4,5)
{
  std::cout << "CTest::CTest" << std::endl;
}

bool TestDebugNew()
{
  CTest t1;
  bUseDevice = true;
  std::unique_ptr<CTest> pBingo{ new CTest };
  t1.tensor.resize(1,2,1,2,1,2);
  bUseDevice = false;
  t1.tensor.resize(2,2,2,2,2,2);
  return true;
}
#endif

int main(int argc, char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn = EXIT_SUCCESS;
  try
  {
    #ifdef TEST_DEBUG_NEW
    std::cout << "main() :: before Debug()" << std::endl;
    if( TestDebugNew() ) return iReturn;
    #endif
    #pragma omp parallel
    std::cout << "Hello world!\n";
  }
  catch(const std::exception &e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  } catch( ... ) {
    std::cerr << "Error: Unknown exception" << std::endl;
    iReturn = EXIT_FAILURE;
  }
  return iReturn;
}
