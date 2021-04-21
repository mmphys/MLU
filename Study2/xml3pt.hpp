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
using QuarkList = std::map<std::string, Quark, Common::LessCaseInsensitive>;

class AppMaker
{
protected:
  const AppParams &appPar;
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
  explicit AppMaker( const AppParams &appPar_ ) : appPar{appPar_}, l( application, appPar_ ) {}
  virtual ~AppMaker() {}
  void Setup( XmlReader &r, const std::string &sRunSuffix );
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
  static const std::string sXmlTagName;
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
  explicit Study2( const AppParams &appPar_ ) : AppMaker{appPar_} {}
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
  static const std::string sXmlTagName;
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
                                      bool,         HeavyQuark,
                                      bool,         HeavyAnti,
                                      bool,         R2Ratio,
                                      bool,         HeavyLightBackwards,
            Grid::Hadrons::Application::TrajRange,  Timeslices,
                                      std::string,  Taxa,
                                      bool,         DoNegativeMomenta,
                                      std::string,  deltaT,
                                      std::string,  gamma,
                                      std::string,  Heavy,
                               std::vector<Decay>,  Decays)
  };
protected:
  MakePar makePar;
  bool ThreePoint;
  std::vector<Taxonomy> Taxa;
  std::vector<int> deltaT;
  std::vector<Gamma::Algebra> gamma;
  std::vector<Quark *> HeavyQuarks;
  int CountHeavyMomenta;

public:
  explicit Study3( const AppParams &appPar_ ) : AppMaker{appPar_} {}
  virtual std::string RunID() const;
protected:
  virtual void Setup( XmlReader &r );
  virtual void Make();
  void MakeStudy3( const Decay &d );
};

