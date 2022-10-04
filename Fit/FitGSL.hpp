/*************************************************************************************
 
 Use GSL as fitting exgine (OPTIONAL)

 Source file: FitGSL.hpp
 
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

#ifndef FitGSL_hpp
#define FitGSL_hpp

#include "Fitter.hpp"
#include "FitterThread.hpp"
#include <gsl/gsl_multifit_nlinear.h>

// Several of these will be running at the same time on different threads during a fit
class FitterThreadGSL : public FitterThread
{
protected:
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_workspace * ws;
  int  f( const Vector &FitParams, Vector &Errors );
  int df( const Vector &x, Matrix &J );
  static inline FitterThreadGSL* ftGSL( void * data ) { return reinterpret_cast<FitterThreadGSL*>(data); }
  static int  sf( const gsl_vector * x, void *data, gsl_vector * f_ )
  { return ftGSL(data)-> f( *reinterpret_cast<const Vector*>(x), *reinterpret_cast<Vector*>(f_) ); }
  static int sdf( const gsl_vector * x, void *data, gsl_matrix * J  )
  { return ftGSL(data)->df( *reinterpret_cast<const Vector*>(x), *reinterpret_cast<Matrix*>(J) ); }
  // Fitter state
  int ConvergeReason;
  std::size_t nevalf;
  std::size_t nevaldf;
  void DumpParamsFitter( std::ostream &os ) const override; //TODO: Deuglify
  void ReplicaMessage( std::ostream &os ) const override; //TODO: Deuglify
  void InitialiseGSL();
  static void SayConvergeReason( std::ostream &os, int ConvergeReason );
public:
  FitterThreadGSL( const Fitter &Fitter, bool bCorrelated, ModelFile &OutputModel );
  FitterThreadGSL( const FitterThreadGSL &ftGSL );
  FitterThreadGSL( FitterThreadGSL &&ftGSL ) = delete;
  FitterThreadGSL &operator=( FitterThreadGSL &&ftGSL ) = delete;
  const FitterThreadGSL &operator=( const FitterThreadGSL &ftGSL ) = delete;
  virtual ~FitterThreadGSL();
  FitterThread * Clone() const override;
  void Minimise( int iNumGuesses ) override;
  int NumRetriesGuess() const override { return parent.Retry; };
  int NumRetriesFit() const override { return parent.Retry; };
  std::string Description() const override;
};

struct FitterGSL : public Fitter
{
  enum class TRS{ lm, lmaccel, dogleg, ddogleg, subspace2D };
  const TRS trs;
  explicit FitterGSL( TRS trs_, const Common::CommandLine &cl, const DataSet &ds,
                      std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                      CovarParams &&cp )
  : Fitter( cl, ds, ModelArgs, opNames, std::move( cp ) ), trs{trs_} {}
  static FitterGSL * Make( const std::string &FitterArgs, const Common::CommandLine &cl, const DataSet &ds,
                           std::vector<Model::Args> &ModelArgs, const std::vector<std::string> &opNames,
                           CovarParams &&cp );
protected:
  const std::string &Type() const override;
  FitterThread * MakeThread( bool bCorrelated, ModelFile &OutputModel ) override
  {
    return new FitterThreadGSL( *this, bCorrelated, OutputModel );
  }
};

#endif // FitGSL_hpp
