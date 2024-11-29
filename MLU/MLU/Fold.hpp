/**

 Mike's lattice QCD utilities
 
 Source file: Fold.hpp
 
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

#ifndef MLU_Fold_hpp
#define MLU_Fold_hpp

#include <MLU/Sample.hpp>

BEGIN_MLU_NAMESPACE

template <typename T>
class Fold : public Sample<T>
{
public:
  int NtUnfolded_ = 0;
  Reality reality = Reality::Unknown;
  Parity parity = Parity::Unknown;
  Sign sign = Sign::Unknown;
  bool t0Negated = false;
  bool Conjugated = false;
  std::vector<std::string> BootstrapList;
  using Base = Sample<T>;
  using Base::Base;
  int NtUnfolded() const override { return NtUnfolded_; }
  const std::string &DefaultGroupName() override;
  bool bFolded() override;
  void SummaryComments( std::ostream &os, bool bVerboseSummary = false,
                        bool bShowColumnNames = true ) const override;
  void ReadAttributes( ::H5::Group &g ) override;
  int WriteAttributes( ::H5::Group &g ) const override;
};

END_MLU_NAMESPACE
#endif // MLU_Fold_hpp
