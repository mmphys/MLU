/*************************************************************************************
 
 Utility for making repairs
 
 Source file: repair.cpp
 
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

#include <stdio.h>
#include <typeinfo>
#include <MLU/Common.hpp>

int main(int argc, char *argv[])
{
  // We're not using C-style I/O
  std::ios_base::sync_with_stdio(false);
  if( argc != 2 ) {
    std::cout << "usage: repair filename\n";
    exit(1);
  }
  const std::string FileName{ argv[1] };
  std::cout << "Repairing " << FileName << "\n";
  static const std::string dsThresholdName{ "_Grid_dataset_threshold" };
  try {
    unsigned int dsThresholdVal{ 6 };
    H5::H5File f{ FileName, H5F_ACC_RDWR };
    H5::Group g{ f.openGroup( "/" ) };
    H5::Attribute a;
    bool bOK = true;
    H5E_auto2_t h5at;
    void      * f5at_p;
    H5::Exception::getAutoPrint(h5at, &f5at_p);
    H5::Exception::dontPrint();
    try {
      a = g.openAttribute(dsThresholdName);
    } catch(const H5::Exception &) {
      bOK = false;
      H5::Exception::clearErrorStack();
    }
    H5::Exception::setAutoPrint(h5at, f5at_p);
    if( bOK )
    {
      a.read( H5::PredType::NATIVE_UINT, &dsThresholdVal );
      std::cout << "Already exists: " << dsThresholdName << "=" << dsThresholdVal << "\n";
      return EXIT_SUCCESS;
    }
    // Create the attribute as a single 32 bit unsigned int
    std::cout << "Creating attribute: " << dsThresholdName << "=" << dsThresholdVal << "\n";
    hsize_t dims[1] = { 1 };
    H5::DataSpace dsAttr = H5::DataSpace( 1, dims );
    a = g.createAttribute(dsThresholdName, H5::PredType::STD_U32LE, dsAttr );
    a.write( H5::PredType::NATIVE_UINT, &dsThresholdVal );
  }
  catch(const H5::Exception &)
  {
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
