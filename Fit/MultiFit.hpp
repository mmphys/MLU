/*************************************************************************************
 
 Fast (using OpenMP) multi-model fits to lattice QCD correlators
 
 Source file: MultiFit.hpp
 
 Copyright (C) 2019-2022
 
 Author: Michael Marshall <Mike@lqcd.me>
 
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

#ifndef MultiFit_hpp
#define MultiFit_hpp

#include <MLU/Common.hpp>

// Uncomment the next line if your cmath doesn't define M_PI etc by default
//#define _USE_MATH_DEFINES
#include <cmath>
#include <set>
//#include <sys/stat.h>

using scalar = double;
using Matrix = Common::Matrix<scalar>;
using Vector = Common::Vector<scalar>;
using MatrixView = Common::MatrixView<scalar>;
using VectorView = Common::VectorView<scalar>;
using Sample = Common::Sample<scalar>;
using SamplePtr = std::unique_ptr<Sample>;
using Fold = Common::Fold<scalar>;
using FoldPtr = std::unique_ptr<Fold>;
using vCorrelator = std::vector<Fold>;
using ModelFile = Common::Model<scalar>;
using ModelFilePtr = std::unique_ptr<ModelFile>;
using DataSet = Common::DataSet<scalar>;
using JackBoot = Common::JackBoot<scalar>;
using ValWithEr = Common::ValWithEr<scalar>;
using ConstantSource = Common::ConstantSource;
using vString = std::vector<std::string>;
using vInt = std::vector<int>;
using UniqueNames = Common::UniqueNames;
using Param = Common::Param;
using Params = Common::Params;
using ParamsPairs = Common::ParamsPairs;
using GSLLG = Common::GSLLibraryGlobal;

// Indices for operators in correlator names
constexpr int idxSrc{ 0 };
constexpr int idxSnk{ 1 };
extern const char * pSrcSnk[2];

#endif // MultiFit_hpp
