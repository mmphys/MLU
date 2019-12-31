/*************************************************************************************
 
 Create XML for heavy-light semi-leptonics
 
 Source file: hlxml.cpp
 
 Copyright (C) 2019
 
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
 *************************************************************************************/
/*  END LEGAL */

#include "Common.hpp"
#include <Hadrons/Utilities/Contractor.hpp>
//#include <Hadrons/Application.hpp>
//#include <Hadrons/Modules.hpp>
//#include <LatAnalyze/Core/OptParser.hpp>

//using namespace Grid;
//using namespace Hadrons;

int main(int argc, const char *argv[])
{
  int iReturn{ EXIT_SUCCESS };
  try
  {
  using Common::CommandLine;
  static const char longSwitch[]{ "balderdash" };
  const std::initializer_list<CommandLine::SwitchDef> l = {
    {"b", CommandLine::SwitchType::Multiple, nullptr},
    {longSwitch, CommandLine::SwitchType::Single, "Hoping this not set!"},
    {"p", CommandLine::SwitchType::Single, "Porcupine" },
    {"Y", CommandLine::SwitchType::Flag, nullptr}
  };
  CommandLine cl( argc, argv, l );
  std::cout << cl << std::endl;
  std::cout << "b[1]=" << cl.SwitchValue<double>("b", 1) << std::endl;
  std::cout << "b[2]=" << cl.SwitchValue<int>   ("b", 2) << std::endl;
  std::cout << "b[3]=" << cl.SwitchValue<int>   ("b", 3) << std::endl;
  std::cout << "b[1]=" << cl.SwitchValue<std::complex<double>>("b", 1) << std::endl;
  std::cout << longSwitch << "=" << cl.SwitchValue<std::string>(longSwitch) << std::endl;
  std::cout << longSwitch << "[1]=" << cl.SwitchValue<std::string>(longSwitch, 1) << std::endl;

  // Decode command-line
  const std::string sOutFileName{ argc == 2 && argv[1] && argv[1][0] ? argv[1] : "hl_test.xml" };
  std::cout << "Building: " << sOutFileName << std::endl;

  // Test
  Contractor::ContractorPar par;
  par.global.nt = 64;
  par.global.trajCounter.start = 3200;
  par.global.trajCounter.step = 40;
  par.global.trajCounter.end = par.global.trajCounter.start + par.global.trajCounter.step;
  par.global.output = ".";
  par.global.diskVectorDir = "tmp_";
  par.global.diskVectorDir += "boinkly";
  //std::cout << par.global << std::endl;

  // Write output file
  std::cout << "Writing: " << sOutFileName << std::endl;
  Grid::XmlWriter w{ sOutFileName };
  write(w, "global",    par.global);
  write(w, "a2aMatrix", par.a2aMatrix);
  write(w, "product",   par.product);
  } catch (const std::exception &e) {
    std::cerr << "Exception: " << e.what() << std::endl;
    iReturn = EXIT_FAILURE;
  }
  return iReturn;
}
