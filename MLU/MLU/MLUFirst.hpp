/**

 Mike's lattice QCD utilities
 
 Source file: MLUFirst.hpp
 
 Copyright (C) 2019 - 2024
 
 Author: Michael Marshall
 
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

// This must be the first header included in any MLU program (because it defines our namespace)

#ifndef MLUFirst_hpp
#define MLUFirst_hpp

#define MLU_NAMESPACE       MLU
#define BEGIN_MLU_NAMESPACE namespace MLU_NAMESPACE {
#define END_MLU_NAMESPACE   };

// Undefine to run without Hadrons guessers (very old code)
#define MLU_HADRONS_HAS_GUESSERS

#endif // MLUFirst_hpp
