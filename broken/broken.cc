//
//  main.cpp
//  broken
//
//  Created by s1786208 on 12/07/2019.
//  Copyright © 2019 sopa. All rights reserved.
//  From code submitted by Peter Boyle as g++ bug report
//  https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66153

#include <vector>
#include <complex>
#include <type_traits>
#include <iostream>

//#define USE_CONSTEXPR

typedef std::complex<double> ComplexD;

template <class T> class TypeMapper {
public:
#ifdef USE_CONSTEXPR
  static constexpr int NestLevel{ T::NestLevel };
#else
  enum { NestLevel = T::NestLevel };
#endif
};

template<> class TypeMapper<ComplexD> {
public:
#ifdef USE_CONSTEXPR
  static constexpr int NestLevel{ 0 };
#else
  enum { NestLevel = 0 };
#endif
};

template<class obj> class Container {
public:
  std::vector<obj> data;
  Container(int size) : data (size){};
};

template<class obj> class Recursive {
public:
#ifdef USE_CONSTEXPR
  static constexpr int NestLevel{ TypeMapper<obj>::NestLevel + 1 };
#else
  enum { NestLevel = TypeMapper<obj>::NestLevel + 1};
#endif
  obj internal;
};

template<int N,class obj,typename std::enable_if<N==obj::NestLevel >::type * = nullptr > auto function(const obj &arg)-> obj
{
  std::cout<<"Leaf "<<obj::NestLevel<<std::endl;
  return arg;
  }
  template<int N,class obj,typename std::enable_if<N!=obj::NestLevel >::type * = nullptr > auto function(const obj &arg)-> obj
  {
    std::cout<<"Node "<<obj::NestLevel<<std::endl;
    obj ret;
    ret.internal=function<N>(arg.internal);
    return ret;
  }
  
  template<int N,class obj> auto function(const Container<obj> & arg)-> Container<decltype(function<N>(arg.data[0]))>
  {
    Container<decltype(function<N>(arg.data[0]))> ret(arg.data.size());
    for(int ss=0;ss<arg.data.size();ss++){
      ret.data[ss] = function<N>(arg.data[ss]);
    }
    return ret;
  }
  
  
  int main(int argc,char **argv)
  {
    Container<Recursive<Recursive<ComplexD> > > array(10);
    Container<Recursive<Recursive<ComplexD> > > ret(10);
    
    ret = function<1>(array);
  }
