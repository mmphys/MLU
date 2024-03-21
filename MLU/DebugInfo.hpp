/*****************************************************************************
 
 Debugging info / signal handling - Courtesy of Grid (util/Init.h)

 Source file: DebugInfo.hpp
 
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
 ****************************************************************************/
/*  END LEGAL */

#ifndef MLU_DebugInfo_hpp
#define MLU_DebugInfo_hpp

#include <type_traits>

#include <MLU/MLUFirst.hpp>

BEGIN_MLU_NAMESPACE

extern bool Grid_exit_handler_disable; // Set to true just before exit so nothing printed
void Grid_debug_handler_init(void);

END_MLU_NAMESPACE
#endif
