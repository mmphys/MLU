/**
 
 Debugging info / signal handling - Courtesy of Grid (util/Init.h)

 Source file: DebugInfo.cpp

 Copyright (C) 2015-2024
 
 Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 Author: Michael Marshall (minor edits to integrate with MLU 2022)

 This file is part of Meson Lattice Utilities (MLU).
 
 MLU is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
 
 MLU is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License along
 with MLU. If not, see <https://www.gnu.org/licenses/>

**/

#include <MLUconfig.h>
#include "DebugInfo.hpp"

#include <cstdio>
#include <cxxabi.h>
#include <fenv.h>
#include <memory>
#include <string>
#include <signal.h>

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#endif

#ifdef __APPLE__

/**
 Provide feenableexcept for Mac OS
 
 - Author: Azusa Yamaguchi <ayamaguc@staffmail.ed.ac.uk>
 Peter Boyle <paboyle@ph.ed.ac.uk>
 Peter Boyle <peterboyle@MacBook-Pro.local>
 paboyle <paboyle@ph.ed.ac.uk>

*/
static int feenableexcept (unsigned int excepts)
{
#if 0
  // Fails on Apple M1
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
  unsigned int old_excepts;  // previous masks
  int iold_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  iold_excepts  = (int) old_excepts;
  return ( fesetenv (&fenv) ? -1 : iold_excepts );
#endif
  return 0;
}
#endif

BEGIN_MLU_NAMESPACE

std::string demangle(const char* name) {
    
  int status = -4; // some arbitrary value to eliminate the compiler warning
    
  // enable c++11 by passing the flag -std=c++11 to g++
  std::unique_ptr<char, void(*)(void*)> res {
    abi::__cxa_demangle(name, NULL, NULL, &status),
      std::free
      };
    
  return (status==0) ? res.get() : name ;
}

#define _NBACKTRACE (256)
void * Grid_backtrace_buffer[_NBACKTRACE];

#define BACKTRACEFILE() {            \
    char string[20];              \
    std::sprintf(string,"MLUbacktrace"); \
    std::FILE * fp = std::fopen(string,"w");        \
    BACKTRACEFP(fp)              \
      std::fclose(fp);              \
  }

#ifdef HAVE_EXECINFO_H
#define BACKTRACEFP(fp) {            \
    int symbols    = backtrace        (Grid_backtrace_buffer,_NBACKTRACE); \
    char **strings = backtrace_symbols(Grid_backtrace_buffer,symbols);  \
    for (int i = 0; i < symbols; i++){          \
      std::fprintf (fp,"BackTrace Strings: %d %s\n",i, demangle(strings[i]).c_str()); std::fflush(fp); \
    }                  \
  }
#else
#define BACKTRACEFP(fp) {            \
    std::fprintf (fp,"BT %d %lx\n",0, __builtin_return_address(0)); std::fflush(fp); \
    std::fprintf (fp,"BT %d %lx\n",1, __builtin_return_address(1)); std::fflush(fp); \
    std::fprintf (fp,"BT %d %lx\n",2, __builtin_return_address(2)); std::fflush(fp); \
    std::fprintf (fp,"BT %d %lx\n",3, __builtin_return_address(3)); std::fflush(fp); \
  }
#endif

void Grid_sa_signal_handler(int sig,siginfo_t *si,void * ptr)
{
  fprintf(stderr,"Caught signal %d\n",si->si_signo);
  fprintf(stderr,"  mem address %llx\n",(unsigned long long)si->si_addr);
  fprintf(stderr,"         code %d\n",si->si_code);
  // Linux/Posix
#ifdef __linux__
  // And x86 64bit
#ifdef __x86_64__
  ucontext_t * uc= (ucontext_t *)ptr;
  struct sigcontext *sc = (struct sigcontext *)&uc->uc_mcontext;
  fprintf(stderr,"  instruction %llx\n",(unsigned long long)sc->rip);
#define REG(A)  printf("  %s %lx\n",#A,sc-> A);
  REG(rdi);
  REG(rsi);
  REG(rbp);
  REG(rbx);
  REG(rdx);
  REG(rax);
  REG(rcx);
  REG(rsp);
  REG(rip);


  REG(r8);
  REG(r9);
  REG(r10);
  REG(r11);
  REG(r12);
  REG(r13);
  REG(r14);
  REG(r15);
#endif
#endif
  fflush(stderr);
  BACKTRACEFP(stderr);
  fprintf(stderr,"Called backtrace\n");
  fflush(stdout);
  fflush(stderr);
  exit(0);
  return;
};

bool Grid_exit_handler_disable = false;

void Grid_exit_handler(void)
{
  if( !Grid_exit_handler_disable )
  {
    BACKTRACEFP(stdout);
    fflush(stdout);
  }
}
void Grid_debug_handler_init(void)
{
  struct sigaction sa;
  sigemptyset (&sa.sa_mask);
  sa.sa_sigaction= Grid_sa_signal_handler;
  sa.sa_flags    = SA_SIGINFO;
  sigaction(SIGSEGV,&sa,NULL);
  sigaction(SIGTRAP,&sa,NULL);
  sigaction(SIGBUS,&sa,NULL);
  sigaction(SIGUSR2,&sa,NULL);

  // In my case I let floating point exceptions occur and user code looks for NaN
  /*feenableexcept( FE_INVALID|FE_OVERFLOW|FE_DIVBYZERO);
  sigaction(SIGFPE,&sa,NULL);*/
  sigaction(SIGKILL,&sa,NULL);
  sigaction(SIGILL,&sa,NULL);

  atexit(Grid_exit_handler);
}

END_MLU_NAMESPACE
