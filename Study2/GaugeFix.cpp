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

using GImpl = PeriodicGimplD;
using GaugeMat = typename GImpl::GaugeLinkField;
//using GaugeLorentz = typename GImpl::GaugeField;

Grid::Hadrons::MGauge::GaugeFix::Par GFPar;

template <class fobj, class sobj>
struct XformSimpleUnmunger
{
  void operator()(sobj &in, fobj &out)
  {
    for (int i = 0; i < Nc; i++) {
      for (int j = 0; j < Nc; j++) {
        out()()(i, j) = in()()(i, j);
      }
    }
  }
};

// Copies header fields from the gauge field I read into the modified one I'm writing
class NerscIOExtra : public NerscIO
{
public:
  struct NerscCache
  {
    FieldMetaData fmd;
    std::string   Message;
  };
  static NerscCache Cache;

  template<typename Impl> class CachedGaugeFixedStatistics
  {
  public:
    void operator()(Lattice<vLorentzColourMatrixD> & data,FieldMetaData &header)
    {
      header.link_trace=WilsonLoops<Impl>::linkTrace(data);
      header.plaquette =WilsonLoops<Impl>::avgPlaquette(data);
      header.hdr_version = Cache.fmd.hdr_version;
      header.storage_format = Cache.fmd.storage_format;
      {
        std::ostringstream ss;
        ss << Cache.fmd.ensemble_id << " (gauge fixed: " << Cache.Message << ")";
        header.ensemble_id = ss.str();
      }
      header.ensemble_label = Cache.fmd.ensemble_label;
      header.sequence_number = Cache.fmd.sequence_number;
    }
  };

  typedef CachedGaugeFixedStatistics<PeriodicGimplD> PeriodicCachedGaugeFixedStatistics;

  template<typename Impl> class CachedGaugeXformStatistics
  {
  public:
    void operator()(Lattice<vColourMatrixD> & data,FieldMetaData &header)
    {
      header.link_trace = 0;
      header.plaquette  = 0;
      header.hdr_version = Cache.fmd.hdr_version;
      header.storage_format = Cache.fmd.storage_format;
      {
        std::ostringstream ss;
        ss << Cache.fmd.ensemble_id << " (gauge xform: " << Cache.Message << ")";
        header.ensemble_id = ss.str();
      }
      header.ensemble_label = Cache.fmd.ensemble_label;
      header.sequence_number = Cache.fmd.sequence_number;
    }
  };
  typedef CachedGaugeXformStatistics<PeriodicGimplD> CachedPeriodicGaugeXformStatistics;
  
  template<class XformStats=CachedPeriodicGaugeXformStatistics>
  static inline void writeGaugeXform(Lattice<vColourMatrixD > &xform,
                                        std::string file,
                                        int, // two_row,
                                        int) // bits32)
  {
    typedef vColourMatrixD vobj;
    typedef typename vobj::scalar_object sobj;
    
    FieldMetaData header;
    ///////////////////////////////////////////
    // Following should become arguments
    ///////////////////////////////////////////
    header.sequence_number = 1;
    header.ensemble_id     = "UKQCD";
    header.ensemble_label  = "DWF";
    
    using fobj3D = ColourMatrixD;
    
    GridBase *grid = xform.Grid();
    
    GridMetaData(grid,header);
    assert(header.nd==4);
    XformStats Stats; Stats(xform,header);
    MachineCharacteristics(header);
    
    uint64_t offset;
    
    // Sod it -- always write 3x3 double
    header.floating_point = std::string("IEEE64BIG");
    header.data_type      = std::string("1D_SU3_GAUGE_3x3");
    XformSimpleUnmunger<fobj3D,sobj> munge;
    if ( grid->IsBoss() ) {
      NerscIO::truncate(file);
      offset = NerscIO::writeHeader(header,file);
    }
    grid->Broadcast(0,(void *)&offset,sizeof(offset));
    
    uint32_t nersc_csum,scidac_csuma,scidac_csumb;
    BinaryIO::writeLatticeObject<vobj,fobj3D>(xform,file,munge,offset,header.floating_point,
                                              nersc_csum,scidac_csuma,scidac_csumb);
    header.checksum = nersc_csum;
    if ( grid->IsBoss() ) {
      NerscIO::writeHeader(header,file);
    }
    
    std::cout<<GridLogMessage <<"Written NERSC xform on "<< file << " checksum "
    <<std::hex<<header.checksum<<std::endl;
  }
};

NerscIOExtra::NerscCache NerscIOExtra::Cache;

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
    std::istringstream is( s );
    is >> std::boolalpha >> bEnableFix; // Doesn't matter if this fails
  }
  // Cache the parameters I'm using
  if( bEnableFix )
  {
    std::ostringstream os;
    os << "alpha=" << GFPar.alpha
       << ", maxiter=" << GFPar.maxiter
       << ", Omega_tol=" << GFPar.Omega_tol
       << ", Phi_tol=" << GFPar.Phi_tol
       << ", gaugeFix=" << GFPar.gaugeFix
       << std::boolalpha << ", Fourier=" << GFPar.Fourier;
    NerscIOExtra::Cache.Message = os.str();
  }
  else
    NerscIOExtra::Cache.Message.clear();

  std::cout<<GridLogMessage << "Reading NERSC Gauge file " << InFileName << Common::NewLine;
  GridCartesian * UGrid = Grid::SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                               GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                               GridDefaultMpi());
  LatticeGaugeField U( UGrid );
  NerscIO::readConfiguration( U, NerscIOExtra::Cache.fmd, InFileName );
  if( !bEnableFix )
  {
    std::cout<<GridLogMessage << "Skipping gauge fixing" << Common::NewLine;
  }
  else
  {
    std::cout<<GridLogMessage << "Gauge fixing: " << NerscIOExtra::Cache.Message << Common::NewLine;
    GaugeMat xform( UGrid );
    FourierAcceleratedGaugeFixer<GImpl>::SteepestDescentGaugeFix( U, xform, GFPar.alpha, GFPar.maxiter, GFPar.Omega_tol, GFPar.Phi_tol, GFPar.Fourier, GFPar.gaugeFix );
    std::cout<<GridLogMessage << "Writing NERSC Gauge xform file " << OutFileName << Common::NewLine;
    NerscIOExtra::writeGaugeXform<NerscIOExtra::CachedPeriodicGaugeXformStatistics>( xform, OutFileName + ".xform", 0, 0 );
  }
  std::cout<<GridLogMessage << "Writing NERSC Gauge file " << OutFileName << Common::NewLine;
  NerscIO::writeConfiguration<NerscIOExtra::PeriodicCachedGaugeFixedStatistics>( U, OutFileName, 0, 0 );
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
