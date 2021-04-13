/*************************************************************************************
 
 Create XML for a 3-pt current insertion study
 Initially use Z2 wall sources.
 Source file: xml3pt.cpp
 Copyright (C) 2020
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

#include <MLU/Common.hpp>
#include <stdexcept>
#include <Grid/Grid.h>
#include <Hadrons/Modules.hpp>

using namespace Grid;

using GImpl = Grid::WilsonImplD;
//using GaugeField = typename GImpl::GaugeField;
//using GaugeMat = typename GImpl::GaugeLinkField;

Grid::Hadrons::MGauge::GaugeFix::Par GFPar;
FieldMetaData fmdBuffer;

template<typename Impl> class GaugeFixedStatistics
{
public:
  void operator()(Lattice<vLorentzColourMatrixD> & data,FieldMetaData &header)
  {
    header.link_trace=WilsonLoops<Impl>::linkTrace(data);
    header.plaquette =WilsonLoops<Impl>::avgPlaquette(data);
    header.hdr_version = fmdBuffer.hdr_version;
    header.storage_format = fmdBuffer.storage_format;
    {
      std::ostringstream ss;
      ss << fmdBuffer.ensemble_id << " (" << GFPar.gaugeFix << " gauge fixed)";
      header.ensemble_id = ss.str();
    }
    header.ensemble_label = fmdBuffer.ensemble_label;
    header.sequence_number = fmdBuffer.sequence_number;
    header.creation_date = fmdBuffer.creation_date;
  }
};
typedef GaugeFixedStatistics<PeriodicGimplD> PeriodicGaugeFixedStatistics;
typedef GaugeFixedStatistics<ConjugateGimplD> ConjugateGaugeFixedStatistics;

void GaugeFix( const std::string &ParamsXml, const std::string &InFileName, const std::string &OutFileName )
{
  std::cout<<GridLogMessage << "Reading gauge fixing parameters from " << ParamsXml << Common::NewLine;
  static const std::string sXmlTopLevel{ "GaugeFix" };
  static const std::string sXmlTagName{ "GaugeFix" };
  bool bEnableFix{ true };
  {
    XmlReader r( ParamsXml, false, sXmlTopLevel );
    read( r, sXmlTagName, GFPar );
    // Read optional boolean tag "Disable"
    std::string s;
    read( r, "EnableFix", s );
    std::istringstream ss( s );
    ss >> std::boolalpha >> bEnableFix; // Doesn't matter if this fails
  }

  std::cout<<GridLogMessage << "Reading NERSC Gauge file " << InFileName << Common::NewLine;
  GridCartesian * UGrid = Grid::SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                               GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                               GridDefaultMpi());
  LatticeGaugeField U( UGrid );
  NerscIO::readConfiguration( U, fmdBuffer, InFileName );
  if( !bEnableFix )
  {
    std::cout<<GridLogMessage << "Skipping gauge fixing" << Common::NewLine;
  }
  else
  {
    std::cout<<GridLogMessage << "Gauge fixing using " << GFPar.gaugeFix << " gauge" << Common::NewLine;
    throw std::runtime_error( "Debug: fail on attempt to gauge fix" );
    //GaugeMat M( UGrid );
    FourierAcceleratedGaugeFixer<PeriodicGaugeImpl<GImpl>>::SteepestDescentGaugeFix( U, GFPar.alpha, GFPar.maxiter, GFPar.Omega_tol, GFPar.Phi_tol, GFPar.Fourier, GFPar.gaugeFix );
  }
  std::cout<<GridLogMessage << "Writing NERSC Gauge file " << OutFileName << Common::NewLine;
  NerscIO::writeConfiguration<PeriodicGaugeFixedStatistics>( U, OutFileName, 0, 0 );
}

int main(int argc, char *argv[])
{
  std::ios_base::sync_with_stdio( false );
  int iReturn{ EXIT_SUCCESS };
  bool bShowUsage{ true };
  bool bInitDone{ false };
  if( argc > 3 )
  {
    try
    {
      const std::string ParamsXml{ argv[1] };
      const std::string InFileName{ argv[2] };
      const std::string OutFileName{ argv[3] };
      if( !Common::FileExists( ParamsXml ) )
      {
        throw std::invalid_argument( "Xml parameters file \"" + ParamsXml + "\" doesn't exist" );
      }
      if( !Common::FileExists( InFileName ) )
      {
        throw std::invalid_argument( "NERSC Gauge file \"" + InFileName + "\" doesn't exist" );
      }
      bShowUsage = false;
      if( Common::FileExists( OutFileName ) )
      {
        std::cout << "NERSC Gauge file \"" + OutFileName + "\" already exists - skipping" << std::endl;
      }
      else
      {
        bInitDone = true;
        Grid::Grid_init( &argc, &argv );
        GaugeFix( ParamsXml, InFileName, OutFileName );
        std::cout<<GridLogMessage<< "NERSC Gauge file " << OutFileName << " created" << Common::NewLine;
      }
    }
    catch(const std::exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      iReturn = EXIT_FAILURE;
    } catch( ... ) {
      std::cerr << "Error: Unknown exception" << std::endl;
      iReturn = EXIT_FAILURE;
    }
  }
  if( bInitDone )
    Grid::Grid_finalize();
  if( bShowUsage )
  {
    ( iReturn == EXIT_SUCCESS ? std::cout : std::cerr ) << "usage: GaugeFix XmlParams InGauge OutGauge <GridOptions>\n"
    "Gauge fix InGauge using XmlParams writing result to OutGauge.\n";
  }
  return iReturn;
}
