/*************************************************************************************
 
 Create a Hadrons Application, top-down building all required dependencies
 Source file: HadronsApp.hpp
 Copyright (C) 2021
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

#include "HadronsApp.hpp"

/**************************
 Base class for my application maker
**************************/

// Map from (case insensitive) quark names to Quarks
using QuarkList    = std::map<std::string, Quark  , Common::LessCaseInsensitive>;
using QuarkNameList = std::set<std::string, Common::LessCaseInsensitive>;

class AppMaker
{
protected:
  const AppParams &appPar;
  const std::string sXmlTagName;
  Application application;
  HModList l;

  QuarkList Q;

  void WriteTaxonomy( std::ostringstream &s, const std::vector<Taxonomy> &Taxa ) const;
  void ReadQuarks( XmlReader &r );
  std::vector<Quark *> ValidateQuarkList( const std::string &sList, const std::string &sName, std::size_t Min );
public:
  std::string JobFilePrefix;
  bool        Run;

  static int MakeThreePoint( int argc, char *argv[], const std::string &sXmlFilename, const std::string &sRunSuffix );
  explicit AppMaker( const AppParams &appPar_, const std::string &sXmlTagName_ )
  : appPar{appPar_}, sXmlTagName{ sXmlTagName_ }, l( application, appPar_ ) {}
  virtual ~AppMaker() {}
  void SetupBase( XmlReader &r, const std::string &sRunSuffix, Grid::GridCartesian * grid );
  virtual std::string RunID() const = 0;
protected:
  virtual void Setup( XmlReader &r ) = 0;
  virtual void Make() = 0;
};

/**************************
 Study 2: generic study for heavy-light meson decays
**************************/

class Study2 : public AppMaker
{
public:
  struct MakePar: Serializable {
      GRID_SERIALIZABLE_CLASS_MEMBERS(MakePar,
                                      bool,         TwoPoint,
                                      bool,         HeavyQuark,
                                      bool,         HeavyAnti,
            Grid::Hadrons::Application::TrajRange,  Timeslices,
                                      std::string,  Taxa,
                                      std::string,  Momenta,
                                      bool,         DoNegativeMomenta,
                                      std::string,  deltaT,
                                      std::string,  gamma,
                                      std::string,  Heavy,
                                      std::string,  Spectator,
                                      std::string,  SpectatorExtra2pt)
  };
protected:
  MakePar makePar;
  bool ThreePoint;
  std::vector<Taxonomy> Taxa;
  std::vector<Common::Momentum> Momenta;
  std::vector<int> deltaT;
  std::vector<Gamma::Algebra> gamma;
  std::vector<Quark *> HeavyQuarks;
  std::vector<Quark *> SpectatorQuarks;
  std::vector<Quark *> SpectatorQuarks2pt;

public:
  explicit Study2( const AppParams &appPar, const std::string &sXmlTagName ) : AppMaker{appPar,sXmlTagName} {}
  virtual std::string RunID() const;
protected:
  virtual void Setup( XmlReader &r );
  virtual void Make();
};

/**************************
 Study 3: compute R2 and R1 ratios for meson decays
**************************/

class Study3 : public AppMaker
{
public:
  struct HeavyMomenta: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(HeavyMomenta,
                                  std::string,  qHeavy,
                                  std::string,  Momenta )
  };

  struct Decay: Serializable {
  GRID_SERIALIZABLE_CLASS_MEMBERS(Decay,
                                  std::string,  name,
                                  std::string,  qLight,
                                  std::string,  qSpectator,
                    std::vector<HeavyMomenta>,  HeavyMom)
  };

  struct MakePar: Serializable {
      GRID_SERIALIZABLE_CLASS_MEMBERS(MakePar,
                                      bool,         TwoPoint,
                                      std::string,  NumHits,
                                      bool,         HeavyQuark,
                                      bool,         HeavyAnti,
                                      bool,         R1Term1,
                                      bool,         R1Term1Backwards,
                                      bool,         R1Term2,
                                      bool,         R2Terms,
            Grid::Hadrons::Application::TrajRange,  Timeslices,
                                      std::string,  Taxa,
                                      bool,         DoNegativeMomenta,
                                      std::string,  deltaT,
                                      std::string,  gamma,
                               std::vector<Decay>,  Decays)
  };
protected:
  MakePar makePar;
  bool ThreePoint;
  std::vector<Taxonomy> Taxa;
  std::vector<int> deltaT;
  std::vector<Gamma::Algebra> gamma;
  QuarkNameList HeavyQuarks;
  int CountHeavyMomenta;
  std::set<int> UniqueP2;

  // Check whether 0-momentum in the list. Add the negative of each momentum if requested
  bool Check0Negate( std::vector<Common::Momentum> &Momenta, bool bNegate );
  void Contract2pt( const Taxonomy &tax, const Quark &q1, const Quark &q2, int idxMomentum, int t,
                    const std::vector<Common::Momentum> &Momenta, bool bGotp0 );

public:
  explicit Study3( const AppParams &appPar, const std::string &sXmlTagName ) : AppMaker{appPar,sXmlTagName} {}
  virtual std::string RunID() const;
protected:
  virtual void Setup( XmlReader &r );
  virtual void Make();
  void MakeStudy3( unsigned int t, const Decay &d );
};

