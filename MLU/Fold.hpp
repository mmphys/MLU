/**

 Mike's lattice QCD utilities
 
 Source file: Fold.hpp
 
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

// Common utilities (no dependencies other than c++ stdlib)

#ifndef MLU_Fold_hpp
#define MLU_Fold_hpp

#include <MLU/Sample.hpp>

BEGIN_COMMON_NAMESPACE

template <typename T>
class Fold : public Sample<T>
{
public:
  int NtUnfolded = 0;
  Reality reality = Reality::Unknown;
  Parity parity = Parity::Unknown;
  Sign sign = Sign::Unknown;
  bool t0Negated = false;
  bool Conjugated = false;
  std::vector<std::string> BootstrapList;
  using Base = Sample<T>;
  using Base::Base;
  const std::string &DefaultGroupName() override;
  bool bFolded() override;
  void SummaryComments( std::ostream &os, bool bVerboseSummary = false ) const override;
  void ReadAttributes( ::H5::Group &g ) override;
  int WriteAttributes( ::H5::Group &g ) const override;
};

END_COMMON_NAMESPACE
#endif // MLU_Fold_hpp
