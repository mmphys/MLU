/*************************************************************************************
 
 Create Ratios, e.g. R1, R2 and associated values such as Z_V
 Source file: Ratio.hpp
 Copyright (C) 2020-2021
 Author: Michael Marshall<Michael.Marshall@ed.ac.uk>
 
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

#pragma once
#include <MLU/Common.hpp>

//static const std::vector<std::string> vHeavy{ "h0", "h1", "h2", "h3" };
//static const int NumHeavy{ static_cast<int>( vHeavy.size() ) };
//static const std::vector<int> vDeltaT{ 12, 14, 16, 20 };
//static const int NumDeltaT{ static_cast<int>( vDeltaT.size() ) };
using Scalar = double;
using Fold = Common::Fold<Scalar>;
using Model = Common::Model<Scalar>;

enum class Freeze
{
  None,
  Central,
  One
};

// Info about a model
struct ModelInfo
{
  Model m;
  int idxE0;  // Parameter Index of E0
  std::string FileName2pt; // Filename of the two-point function containing the raw data
};

// A quark with a sepecified momentum
using QP=std::pair<std::string, Common::Momentum>;

// Case insensitive compare of QP
struct LessQP
{
  bool operator()( const QP &lhs, const QP &rhs ) const
  {
    int i = Common::CompareIgnoreCase( lhs.first, rhs.first );
    if( i )
      return i < 0;
    return lhs.second < rhs.second;
  }
};

// Map from QP -> ModelInfo
struct QPMapModelInfo : public std::map<QP, ModelInfo, LessQP>
{
  int MaxModelSamples;
  Common::SeedType Seed;
  QPMapModelInfo( const std::string &FitListName, const std::string &modelBase, const std::string &C2Base );
};

// A quark with a sepecified DeltaT
struct QDT
{
  std::string q;
  int deltaT;
  std::string opSnk;
  std::string opSrc;
  QDT( std::string q_, int deltaT_, const std::string &opSnk_, const std::string &opSrc_ )
  : q{q_}, deltaT{deltaT_}, opSnk{opSnk_}, opSrc{opSrc_} {}
};

// Case insensitive compare of QP
struct LessQDT
{
  bool operator()( const QDT &lhs, const QDT &rhs ) const
  {
    int i = Common::CompareIgnoreCase( lhs.q, rhs.q );
    if( i )
      return i < 0;
    if( lhs.deltaT != rhs.deltaT )
      return lhs.deltaT < rhs.deltaT;
    i = Common::CompareIgnoreCase( lhs.opSnk, rhs.opSnk );
    if( i )
      return i < 0;
    i = Common::CompareIgnoreCase( lhs.opSrc, rhs.opSrc );
    return i < 0;
  }
};

// Map from a DeltaT -> ModelInfo
struct QDTMapModelInfo : public std::map<QDT, ModelInfo, LessQDT>
{
  int MaxModelSamples;
  Common::SeedType Seed;
  QDTMapModelInfo( const std::string &FitListName, const std::string &modelBase, const std::string &C2Base );
};


// Base class for ratio construction
class Maker
{
public:
  const std::string &inBase;
  const std::string &C2Base;
  const std::string &modelBase;
  const std::string &outBase;
  const Common::Momentum p;
  const Freeze fEnergy;
  const Freeze fZV;
protected:
  std::regex RegExExt;
  const bool RegExSwap;
  QPMapModelInfo model;
protected:
  void AppendOp( std::string &s, const std::string &Op )
  {
    s.append( Common::Underscore );
    s.append( Op );
  }
  void AppendOps( std::string &s, const std::string &Snk, const std::string &Src)
  {
    AppendOp( s, Snk );
    AppendOp( s, Src );
  }
  std::string HeavyKey( const std::string &Heavy ) const;
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const std::string &qSnk, const std::string &qSrc,
                     const ModelInfo &miSnk, const ModelInfo &miSrc,
                     const std::string &opSnk, const std::string &opSrc ) = 0;
public:
  Maker( const std::string &inBase_, const std::string &C2Base_,
         const std::string &modelBase_,const std::string &outBase_,
         std::regex RegExExt_, const bool RegExSwap_,
        const Freeze fEnergy_, const Freeze fZV_, const std::string &FitListName )
  : inBase{inBase_}, C2Base{C2Base_}, modelBase{modelBase_}, outBase{outBase_},
    fEnergy{fEnergy_}, fZV{fZV_}, RegExExt{RegExExt_}, RegExSwap{RegExSwap_},
    model(modelBase_ + FitListName, modelBase_, C2Base_) {}
  virtual ~Maker() {}
  static Maker * Make( const std::string &Type, std::string &TypeParams,
                       const std::string &inBase, const std::string &C2Base,
                       const std::string &modelBase,const std::string &outBase,
                       std::regex RegExExt, const bool RegExSwap,
                       const Freeze fEnergy, const Freeze fZV, const std::string &FitListName );
  void Make( std::string &sFileName );
};

// Make Z_V: eq 2.7 pg 3 https://arxiv.org/pdf/1305.7217.pdf
class ZVMaker : public Maker
{
protected:
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const std::string &qSnk, const std::string &qSrc,
                     const ModelInfo &miSnk, const ModelInfo &miSrc, const std::string &opSnk, const std::string &opSrc );
public:
  ZVMaker( std::string &TypeParams, const std::string &inBase, const std::string &C2Base,
           const std::string &modelBase, const std::string &outBase,
           std::regex RegExExt, const bool RegExSwap,
          const Freeze fEnergy, const Freeze fZV, const std::string &FitListName )
  : Maker( inBase, C2Base, modelBase, outBase, RegExExt, RegExSwap, fEnergy, fZV, FitListName ) {}
};

// Make R1 and R2: eq 2.8 pg 4 https://arxiv.org/pdf/1305.7217.pdf
class R1R2Maker : public Maker
{
protected:
  QDTMapModelInfo ZVmi;
  virtual void Make( const Common::FileNameAtt &fna, const std::string &fnaSuffix,
                     const std::string &qSnk, const std::string &qSrc,
                     const ModelInfo &miSnk, const ModelInfo &miSrc, const std::string &opSnk, const std::string &opSrc );
public:
  R1R2Maker( std::string &TypeParams, const std::string &inBase, const std::string &C2Base,
             const std::string &modelBase, const std::string &outBase,
             std::regex RegExExt, const bool RegExSwap,
             const Freeze fEnergy, const Freeze fZV, const std::string &FitListName )
  : Maker( inBase, C2Base, modelBase, outBase, RegExExt, RegExSwap, fEnergy, fZV, FitListName ),
    ZVmi( modelBase + TypeParams, modelBase, C2Base )
  {
    TypeParams.clear(); // I used this to load my map from
  }
};
