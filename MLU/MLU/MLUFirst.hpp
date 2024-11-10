/**

 Mike's lattice QCD utilities
 
 Source file: MLUFirst.hpp
 
 Copyright (C) 2019 - 2023
 
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
