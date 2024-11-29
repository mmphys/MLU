/**
 
 Debugging info / signal handling - Courtesy of Grid (util/Init.h)

 Source file: DebugInfo.hpp
 
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

#ifndef MLU_DebugInfo_hpp
#define MLU_DebugInfo_hpp

#include <type_traits>

#include <MLU/MLUFirst.hpp>

BEGIN_MLU_NAMESPACE

extern bool Grid_exit_handler_disable; // Set to true just before exit so nothing printed
void Grid_debug_handler_init(void);

END_MLU_NAMESPACE
#endif
