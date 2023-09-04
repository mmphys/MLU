/**

 Mike's lattice QCD utilities
 
 Source file: Model.cpp
 
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

#include "Model.hpp"
#include "DataSet.hpp"

BEGIN_COMMON_NAMESPACE

const std::string ModelBase::EnergyPrefix{ "E" };
const std::string ModelBase::EDiffPrefix{ "EDiff" };
const std::string ModelBase::ConstantPrefix{ "K" };

template <typename T>
JackBootColumn<T> Model<T>::Column( const Param::Key &k, std::size_t Index )
{
  Params::const_iterator it = params.FindPromiscuous( k );
  if( it == params.cend() )
  {
    std::ostringstream os;
    os << k.FullName( Index, std::numeric_limits<std::size_t>::max() ) << " not found";
    throw std::runtime_error( os.str().c_str() );
  }
  const Param &p{ it->second };
  if( Index >= p.size )
  {
    std::ostringstream os;
    os << k.FullName( Index, p.size ) << " not found - only " << p.size << " elements";
    throw std::runtime_error( os.str().c_str() );
  }
  return Base::Column( p.GetOffset( Index, Param::Type::All ) );
}

template <typename T>
JackBootColumn<T> Model<T>::ColumnFrozen( const Param::Key &k, std::size_t Index )
{
  Params::const_iterator it = params.FindPromiscuous( k );
  if( it == params.cend() )
  {
    std::ostringstream os;
    os << k.FullName( Index, std::numeric_limits<std::size_t>::max() ) << " not found";
    throw std::runtime_error( os.str().c_str() );
  }
  const Param &p{ it->second };
  if( Index >= p.size )
  {
    std::ostringstream os;
    os << k.FullName( Index, p.size ) << " not found - only " << p.size << " elements";
    throw std::runtime_error( os.str().c_str() );
  }
  return Base::ColumnFrozen( p.GetOffset( Index, Param::Type::All ) );
}

template <typename T>
UniqueNameSet Model<T>::GetStatColumnNames() const
{
  UniqueNameSet Names;
  const std::size_t NumParams{ params.NumScalars( Common::Param::Type::All ) };
  const std::size_t NumColumns{ Base::ColumnNames.size() };
  for( std::size_t i = NumParams; i < NumColumns; ++i )
    Names.insert( Base::ColumnNames[i] );
  return Names;
}

template <typename T>
Model<T>::Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns )
: Base::Sample( NumSamples, static_cast<int>( Params_.NumScalars( Param::Type::All ) + ExtraColumns.size() ) ),
  params{Params_}
{
  CommonConstruct( ExtraColumns );
}

template <typename T>
Model<T>::Model( int NumSamples, Params Params_, const std::vector<std::string> &ExtraColumns,
                 int CovarSampleSize_, bool CovarFrozen_, SampleSource CovarSource_, std::vector<int> CovarRebin_, int CovarNumBoot_ )
: Base::Sample( NumSamples, static_cast<int>( Params_.NumScalars( Param::Type::All ) + ExtraColumns.size() ) ),
  params{Params_}, CovarSampleSize{CovarSampleSize_}, CovarFrozen{CovarFrozen_},
  CovarSource{CovarSource_}, CovarRebin{CovarRebin_}, CovarNumBoot{CovarNumBoot_}
{
  CommonConstruct( ExtraColumns );
}

template <typename T>
void Model<T>::Read( const char *PrintPrefix, std::string *pGroupName )
{
  // Just read the data using the Seed from the file 
  Base::Read( false, Common::SeedWildcard, PrintPrefix, pGroupName );
}

template <typename T>
void Model<T>::ReadAttributes( ::H5::Group &g )
{
  Base::ReadAttributes( g );
  ::H5::Attribute a;
  a = g.openAttribute(sDoF);
  a.read( ::H5::PredType::NATIVE_INT, &dof );
  a.close();
  std::int8_t i8;
  a = g.openAttribute(sCovarFrozen);
  a.read( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  CovarFrozen = ( i8 != 0 );
  try
  {
    a = g.openAttribute( sNumExponents );
    a.read( ::H5::PredType::NATIVE_INT, &NumExponents );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadVector( g, sGuess+s_C, Guess );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceIn+s_C, CovarIn );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovariance+s_C, Covar );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelation+s_C, Correl );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationCholesky+s_C, CorrelCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceInv+s_C, CovarInv );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationInv+s_C, CorrelInv );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationInvCholesky+s_C, CorrelInvCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCovarianceInvCholesky+s_C, CovarInvCholesky );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    H5::ReadMatrix( g, sCorrelationParam+s_C, CorrelParam );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute( sCorrelationParamNames );
    CorrelParamNames = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    StdErrorMean.Read( g, sStdErrorMean );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    FitInput.Read( g, sFitInput );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    ModelPrediction.Read( g, sModelPrediction );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    ErrorScaled.Read( g, sErrorScaled );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute( sCovarSource );
    std::string s{};
    a.read( a.getStrType(), s );
    a.close();
    std::stringstream ss( s );
    ss >> CovarSource;
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
    CovarSource = SS::Binned; // Default if missing
  }
  CovarRebin.clear();
  try
  {
    a = g.openAttribute( sCovarRebin );
    ::H5::DataSpace dsp = a.getSpace();
    const int rank{ dsp.getSimpleExtentNdims() };
    if( rank != 1 )
      throw std::runtime_error( sCovarRebin + " dimensions " + std::to_string( rank ) + ", expecting 1" );
    hsize_t Num;
    dsp.getSimpleExtentDims( &Num );
    if( Num > std::numeric_limits<int>::max() )
      throw std::runtime_error( sCovarRebin + " too many items " + std::to_string( Num ) );
    std::vector<int> Buffer( Num );
    a.read( ::H5::PredType::NATIVE_INT, &Buffer[0] );
    a.close();
    CovarRebin = std::move( Buffer );
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute(sCovarSampleSize);
    a.read( ::H5::PredType::NATIVE_INT, &CovarSampleSize );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    CovarSampleSize = this->SampleSize;
    ::H5::Exception::clearErrorStack();
  }
  try
  {
    a = g.openAttribute(sCovarNumBoot);
    a.read( ::H5::PredType::NATIVE_INT, &CovarNumBoot );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    CovarNumBoot = 0;
    ::H5::Exception::clearErrorStack();
  }
  ModelType.clear();
  try
  {
    a = g.openAttribute( sModelTypeList );
    ModelType = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  ModelArgs.clear();
  try
  {
    a = g.openAttribute( sModelArgsList );
    ModelArgs = H5::ReadStrings( a );
    a.close();
  }
  catch(const ::H5::Exception &)
  {
    ::H5::Exception::clearErrorStack();
  }
  ::H5::Group gParam;
  if( H5::OpenOptional( gParam, g, sParam ) )
  {
    OldFormatNumExponents = std::numeric_limits<int>::min();
    // Read new format containing parameters
    gParam.close();
    params.ReadH5( g, sParam );
    // Get the fit times
    std::size_t NumFitTimes;
    std::string sAttName{ sFitTime };
    const std::size_t Len{ sAttName.length() };
    a = g.openAttribute( sAttName );
    a.read( H5::Equiv<std::size_t>::Type, &NumFitTimes );
    a.close();
    FitTimes.resize( NumFitTimes );
    for( std::size_t i = 0; i < NumFitTimes; ++i )
    {
      sAttName.resize( Len );
      sAttName.append( std::to_string( i ) );
      a = g.openAttribute( sAttName );
      ::H5::DataSpace dsp = a.getSpace();
      const int nDims{ dsp.getSimpleExtentNdims() };
      bool bError{ true };
      if( nDims == 1 )
      {
        hsize_t Dim;
        dsp.getSimpleExtentDims( &Dim );
        if( Dim > 0 && Dim <= std::numeric_limits<std::size_t>::max() )
        {
          FitTimes[i].resize( static_cast<std::size_t>( Dim ) );
          a.read( ::H5::PredType::NATIVE_INT, &FitTimes[i][0] );
          bError = false;
        }
      }
      if( bError )
        throw std::runtime_error( "Error reading attribute " + sAttName  );
      a.close();
    }
  }
  else
  {
    // Read old format
    a = g.openAttribute(sNumExponents);
    a.read( ::H5::PredType::NATIVE_INT, &OldFormatNumExponents );
    a.close();
    int NumFiles;
    a = g.openAttribute(sNumFiles);
    a.read( ::H5::PredType::NATIVE_INT, &NumFiles );
    a.close();
    int ti;
    a = g.openAttribute(sTI);
    a.read( ::H5::PredType::NATIVE_INT, &ti );
    a.close();
    int tf;
    a = g.openAttribute(sTF);
    a.read( ::H5::PredType::NATIVE_INT, &tf );
    a.close();
    a = g.openAttribute(sOperators);
    OldFormatOpNames = H5::ReadStrings( a );
    a.close();
    // Recreate FitTimes
    const int Extent{ tf - ti + 1 };
    if( Extent <= 0 )
      throw std::runtime_error( "TI " + std::to_string( ti ) + " > TF " + std::to_string( tf ) );
    FitTimes.resize( NumFiles );
    for( std::vector<int> &v : FitTimes )
    {
      v.resize( Extent );
      for( int i = 0; i < Extent; ++i )
        v[i] = ti + i;
    }
  }
}

template <typename T>
void Model<T>::ValidateAttributes()
{
  if( OldFormatNumExponents != std::numeric_limits<int>::min() )
  {
    // Validate old format
    const int NumOps{ static_cast<int>( OldFormatOpNames.size() ) + 1 };
    const int NumExpected{ OldFormatNumExponents * NumOps + 1 };
    if( Base::Nt_ != NumExpected )
    {
      std::ostringstream s;
      s << "Have " << Base::Nt_ << " parameters, but expected " << NumExpected << ", i.e. " << OldFormatNumExponents
      << " exponents * ( " << ( NumOps - 1 ) << " operators + 1 energy ) + chi squared per degree of freedom";
      throw std::runtime_error( s.str().c_str() );
    }
    // Now build parameters
    params.clear();
    Param::Key k;
    k.Object.push_back( this->Name_.Base );
    std::size_t Len{ k.Object[0].find_first_of( '.' ) };
    if( Len != std::string::npos )
      k.Object[0].resize( Len );
    k.Name = EnergyPrefix;
    params.Add( k, OldFormatNumExponents, true, Param::Type::Variable );
    for( const std::string &s : OldFormatOpNames )
    {
      k.Name = s;
      params.Add( k, OldFormatNumExponents, false, Param::Type::Variable );
    }
    params.AssignOffsets();
    OldFormatOpNames.clear();
    // If there's more than one exponent, data need reordering
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Data[0].GetCentral() );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Data[0].GetReplicaMean() );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Data[0].Replica );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Binned[0] );
    ReorderOldFormat( NumOps, OldFormatNumExponents, Base::Raw );
  }
  Base::ValidateAttributes();
  if( Base::Nt_ < NumParams() )
  {
    std::ostringstream s;
    s << "Model has " << NumParams() << " parameters, but " << Base::Nt_ << " columns";
    throw std::runtime_error( s.str().c_str() );
  }
}

template <typename T>
int Model<T>::WriteAttributes( ::H5::Group &g ) const
{
  int iReturn = Base::WriteAttributes( g ) + 4;
  const hsize_t OneDimension{ 1 };
  ::H5::DataSpace ds1( 1, &OneDimension );
  ::H5::Attribute a;
  a = g.createAttribute( sNumExponents, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &NumExponents );
  a.close();
  a = g.createAttribute( sDoF, ::H5::PredType::STD_U16LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &dof );
  a.close();
  std::int8_t i8{ static_cast<std::int8_t>( CovarFrozen ? 1 : 0 ) };
  a = g.createAttribute( sCovarFrozen, ::H5::PredType::STD_U8LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT8, &i8 );
  a.close();
  H5::WriteVector( g, sGuess+s_C, Guess );
  H5::WriteMatrix( g, sCovarianceIn+s_C, CovarIn );
  H5::WriteMatrix( g, sCovariance+s_C, Covar );
  H5::WriteMatrix( g, sCorrelation+s_C, Correl );
  H5::WriteMatrix( g, sCorrelationCholesky+s_C, CorrelCholesky );
  H5::WriteMatrix( g, sCovarianceInv+s_C, CovarInv );
  H5::WriteMatrix( g, sCorrelationInv+s_C, CorrelInv );
  H5::WriteMatrix( g, sCorrelationInvCholesky+s_C, CorrelInvCholesky );
  H5::WriteMatrix( g, sCovarianceInvCholesky+s_C, CovarInvCholesky );
  H5::WriteMatrix( g, sCorrelationParam, CorrelParam );
  if( CorrelParamNames.size() )
  {
    H5::WriteAttribute( g, sCorrelationParamNames, CorrelParamNames );
    iReturn++;
  }
  StdErrorMean.Write( g, sStdErrorMean );
  FitInput.Write( g, sFitInput );
  ModelPrediction.Write( g, sModelPrediction );
  ErrorScaled.Write( g, sErrorScaled );
  {
    std::ostringstream os;
    os << CovarSource;
    a = g.createAttribute( sCovarSource, H5::Equiv<std::string>::Type, ds1 );
    a.write( H5::Equiv<std::string>::Type, os.str() );
    a.close();
    iReturn++;
  }
  if( !CovarRebin.empty() )
  {
    hsize_t Dims[1] = { CovarRebin.size() };
    ::H5::DataSpace dsp( 1, Dims );
    a = g.createAttribute( sCovarRebin, ::H5::PredType::NATIVE_INT, dsp );
    a.write( ::H5::PredType::NATIVE_INT, &CovarRebin[0] );
    a.close();
    dsp.close();
    iReturn++;
  }
  a = g.createAttribute( sCovarSampleSize, ::H5::PredType::STD_U32LE, ds1 );
  a.write( ::H5::PredType::NATIVE_INT, &CovarSampleSize );
  a.close();
  if( CovarNumBoot )
  {
    a = g.createAttribute( sCovarNumBoot, ::H5::PredType::STD_U32LE, ds1 );
    a.write( ::H5::PredType::NATIVE_INT, &CovarNumBoot );
    a.close();
    iReturn++;
  }
  if( ModelType.size() )
    H5::WriteAttribute( g, sModelTypeList, ModelType );
  if( ModelArgs.size() )
    H5::WriteAttribute( g, sModelArgsList, ModelArgs );
  params.WriteH5( g, sParam );
  // Write the fit times
  std::size_t NumFitTimes{ FitTimes.size() };
  std::string sAttName{ sFitTime };
  const std::size_t Len{ sAttName.length() };
  a = g.createAttribute( sAttName, ::H5::PredType::STD_U32LE, ds1 );
  a.write( H5::Equiv<std::size_t>::Type, &NumFitTimes );
  a.close();
  for( std::size_t i = 0; i < NumFitTimes; ++i )
  {
    sAttName.resize( Len );
    sAttName.append( std::to_string( i ) );
    hsize_t Dims[1] = { FitTimes[i].size() };
    ::H5::DataSpace dsp( 1, Dims );
    a = g.createAttribute( sAttName, ::H5::PredType::STD_U16LE, dsp );
    a.write( ::H5::PredType::NATIVE_INT, &FitTimes[i][0] );
    a.close();
    dsp.close();
    iReturn++;
  }
  return iReturn;
}

template <typename T>
bool Model<T>::CheckParameters( int Strictness, scalar_type MonotonicUpperLimit )
{
  if( Strictness < 0 )
    return true;
  bool bResult{ !this->SummaryNames.empty() };
  if( bResult && MonotonicUpperLimit < std::numeric_limits<scalar_type>::max() )
  {
    const bool bVeryStrictNonZero{ ( Strictness & 1 ) != 0 };
    const bool bVeryStrictDifferent{ ( Strictness & 2 ) != 0 };
    using TypeVE = Common::ValWithEr<scalar_type>;
    scalar_type TypeVE::* const pLowZ { bVeryStrictNonZero ? &TypeVE::Min : &TypeVE::Low };
    scalar_type TypeVE::* const pHighZ{ bVeryStrictNonZero ? &TypeVE::Max : &TypeVE::High };
    scalar_type TypeVE::* const pLowD { bVeryStrictDifferent ? &TypeVE::Min : &TypeVE::Low };
    scalar_type TypeVE::* const pHighD{ bVeryStrictDifferent ? &TypeVE::Max : &TypeVE::High };
    for( const Params::value_type &it : params )
    {
      //const Param::Key &pk{ it.first };
      const Param &p{ it.second };
      if( p.type == Param::Type::Variable )
      {
        const std::size_t o{ p.GetOffset( 0, Param::Type::All ) };
        for( std::size_t i = 0; i < p.size; ++i )
        {
          TypeVE &VEoi{ Base::SummaryData( static_cast<int>( o + i ) ) };
          // All elements must be statistically different from zero
          bool bOK{  ( VEoi.*pLowZ > 0 && VEoi.*pHighZ > 0 )
                  || ( VEoi.*pLowZ < 0 && VEoi.*pHighZ < 0 ) };
          // Monotonic parameters (energies) must be < limit
          if( bOK && p.bMonotonic && VEoi.Central > MonotonicUpperLimit )
            bOK = false;
          // All elements must be different from each other
          for( std::size_t j = 0; bOK && j < i; ++j )
          {
            TypeVE &VEoj{ Base::SummaryData( static_cast<int>( o + j ) ) };
            bOK = VEoi.*pLowD > VEoj.*pHighD || VEoi.*pHighD < VEoj.*pLowD;
          }
          // Save results
          VEoi.Check = bOK ? 1 : 0;
          if( !bOK )
            bResult = false;
        }
      }
    }
  }
  return bResult;
}

template <typename T>
void Model<T>::SummaryComments( std::ostream & s, bool bVerboseSummary,
                               bool bShowColumnNames ) const
{
  Base::SummaryComments( s, true, bShowColumnNames );
  for( std::size_t i = 0; i < FitTimes.size(); ++i )
  {
    s << "# Fit Times " << i << Colon;
    for( const int t : FitTimes[i] )
      s << Space << t;
    s << NewLine;
  }
  for( std::size_t i = 0; i < ModelType.size(); ++i )
    s << "# ModelType " << i << ": " << ModelType[i] << NewLine;
  for( std::size_t i = 0; i < ModelArgs.size(); ++i )
    s << "# ModelArgs " << i << ": " << ModelArgs[i] << NewLine;
  s << "# NumExponents: " << NumExponents << NewLine;
  if( !params.empty() )
  {
    s << "# Params: ";
    bool bFirst{ true };
    for( const Params::value_type &p : params )
    {
      if( bFirst )
        bFirst = false;
      else
        s << CommaSpace;
      s << p;
    }
    s << NewLine;
  }
  s << "# Frozen covariance: " << CovarFrozen << NewLine
    << "# Covariance source: " << CovarSource << NewLine;
  if( !CovarRebin.empty() )
  {
    s << "# Covariance rebin:";
    for( auto i : CovarRebin )
      s << Space << i;
    s << NewLine;
  }
  // Now write info about statistics
  if( dof < 1 )
    s << "# Extrapolation : All p-values=1\n";
  // Chi^2 per dof
  int iCol{ this->GetColumnIndexNoThrow( sChiSqPerDof ) };
  if( iCol >= 0 )
    s << "# " << sChiSqPerDof << " : " << (*this)(Model<T>::idxCentral,iCol)
      << " (dof=" << dof << ")" << NewLine;
  // Hotelling pValue
  s << "# " << sPValueH << " : Hotelling (p=" << dof << ", m-1=" << ( CovarSampleSize - 1 ) << ")";
  bool bHotellingUsable{ Common::HotellingDist::Usable( dof, CovarSampleSize - 1 ) };
  if( !bHotellingUsable )
    s << ". Not usable - set to ChiSquared p-value instead";
  else
  {
    iCol = this->GetColumnIndexNoThrow( sPValueH );
    if( iCol >= 0 )
      s << " = " << (*this)(Model<T>::idxCentral,iCol);
  }
  s << NewLine;
  // Chi^2 pValue
  iCol = this->GetColumnIndexNoThrow( sPValue );
  if( iCol >= 0 )
    s << "# " << sPValue << " : Chi^2 = " << (*this)(Model<T>::idxCentral,iCol) << NewLine;
}

const char ModelBase::SummaryColumnPrefix[] = "ti tf tiLabel tfLabel NumDataPoints dof SampleSize CovarSampleSize";

template <typename T>
void Model<T>::SummaryColumnNames( std::ostream &os ) const
{
  os << SummaryColumnPrefix << Common::Space;
  Base::SummaryColumnNames( os );
  for( std::size_t i = 1; i < FitTimes.size(); ++i )
    os << Space << "ti" << i << Space << "tf" << i;
}

template <typename T>
void Model<T>::SummaryContentsPrefix( std::ostream &os ) const
{
  if( FitTimes.size() )
  {
    os << FitTimes[0][0] << Space << FitTimes[0].back();
    // Write tiLabel
    os << Space << FitTimes[0][0];
    for( std::size_t i = 1; i < FitTimes.size(); ++i )
      os << Underscore << FitTimes[i][0] << Underscore << FitTimes[i].back();
    // Write tfLabel
    os << Space << FitTimes[0].back();
    for( std::size_t i = 1; i < FitTimes.size(); ++i )
      os << Underscore << FitTimes[i][0] << Underscore << FitTimes[i].back();
  }
  else
    os << "0 0 0 0";
  os << Space << GetExtent();
  os << Space << dof;
  os << Space << this->SampleSize;
  os << Space << CovarSampleSize;
}

template <typename T>
void Model<T>::SummaryContentsSuffix( std::ostream &os ) const
{
  // Write all the fit times
  for( std::size_t i = 1; i < FitTimes.size(); ++i )
    os << Space << FitTimes[i][0] << Space << FitTimes[i].back();
  os << NewLine;
}

template <typename T>
void Model<T>::SummaryContents( std::ostream &os ) const
{
  SummaryContentsPrefix( os );
  os << Space;
  Base::SummaryContents( os );
  SummaryContentsSuffix( os );
}

template <typename T> std::vector<ValWithEr<typename Model<T>::scalar_type>>
Model<T>::GetValWithEr( const Params &ParamNames, const UniqueNameSet &StatNames ) const
{
  using Scalar = typename Model<T>::scalar_type;
  const ValWithEr<Scalar> Zero(0,0,0,0,0,0);
  const std::size_t NumScalars{ ParamNames.NumScalars( Param::Type::All ) };
  const std::size_t NumStats{ StatNames.size() };
  std::vector<ValWithEr<Scalar>> v;
  v.reserve( NumScalars + NumStats );
  for( const Params::value_type &it : ParamNames )
  {
    const Param::Key &k{ it.first };
    const Param &p{ it.second };
    const std::size_t NumToWrite{ p.size };
    Params::const_iterator itMe{ params.find( k ) };
    const std::size_t NumHave{ itMe == params.cend() ? 0 : itMe->second.size };
    // Write out the parameters I have
    if( NumHave )
    {
      const Param &pMe{ itMe->second };
      std::size_t MyOffset{ pMe.GetOffset( 0, Param::Type::All ) };
      for( std::size_t i = 0; i < std::min( NumHave, NumToWrite ); ++i )
        v.emplace_back( Base::SummaryData( static_cast<int>( MyOffset + i ) ) );
    }
    // Now write out dummy values for those I don't have
    for( std::size_t i = NumHave; i < NumToWrite; ++i )
      v.emplace_back( Zero );
  }
  // Now write out all the statistics columns
  for( const std::string &StatName : StatNames )
  {
    const int idx{ Base::GetColumnIndexNoThrow( StatName ) };
    v.emplace_back( idx < 0 ? Zero : Base::SummaryData( idx ) );
  }
  return v;
}

template <typename T>
void Model<T>::WriteSummaryTD( const DataSet<T> &ds, const std::string &sOutFileName,
                               bool bVerboseSummary )
{
  const ValWithEr<T> veZero( 0, 0, 0, 0, 0, 0 );
  using Scalar = typename is_complex<T>::Real;
  using namespace CorrSumm;
  assert( std::isnan( NaN ) && "Compiler does not support quiet NaNs" );
  std::ofstream ofs( sOutFileName );
  SummaryHeader<T>( ofs, sOutFileName );
  SummaryComments( ofs, bVerboseSummary );
  // Write column names
  static constexpr int idxData{ 0 };
  //static constexpr int idxTheory{ 1 };
  ofs << "# Correlator data used in fit\n";
  ofs << "field seq fitseq model t tfit ";
  ValWithEr<T>::Header( "data", ofs );
  ofs << Space;
  ValWithEr<T>::Header( "theory", ofs );
  ofs << NewLine;
  // Grab the model data and theory with bootstrapped errors
  const int Extent{ GetExtent() };
  if( !Extent )
    return;
  std::array<VectorView<T>, 2> vv; // View into data and theory replica data. Used to populate Value
  std::vector<T> Buffer; // Scratch buffer for ValWithEr<T>
  std::array<std::vector<ValWithEr<T>>,2> Value; // Data and theory values with errors
  for( std::size_t i = 0; i < Value.size(); ++i )
  {
    JackBoot<T> &ThD{ i == idxData ? FitInput : ModelPrediction };
    ThD.MakeStatistics( Value[i] );
  }
  // Write theory and data CORRELATOR values, with each data point in fit on a separate line
  int Seq = 0;
  int FitSeq = 0;
  std::vector<std::vector<int>> SortedTimes{ FitTimes };
  std::ostringstream osBuffer;
  osBuffer << std::boolalpha << std::setprecision(std::numeric_limits<scalar_type>::digits10+2);
  for( std::size_t m = 0; m < FitTimes.size(); ++m )
  {
    std::sort( SortedTimes[m].begin(), SortedTimes[m].end() );
    for( int t = 0, idxSort = 0; t < ds(m, true).Nt(); ++t )
    {
      const bool bInFit{ t == SortedTimes[m][idxSort] };
      std::ostream &os{ bInFit ? dynamic_cast<std::ostream &>( ofs ) : osBuffer };
      os << "corr " << Seq << Space << ( bInFit ? FitSeq : -1 ) << Space << m
      << Space << t << Space << ( bInFit ? t : -1 );
      if( bInFit )
      {
        for( std::size_t i = 0; i < Value.size(); ++i )
          os << Space << Value[i][FitSeq];
        ++idxSort;
        ++FitSeq;
      }
      else
        os << Space << ds(m, true).SummaryData( t ) << Space << veZero;
      os << NewLine;
      ++Seq;
    }
  }
  if( !osBuffer.str().empty() )
  {
    ofs << "\n\n# Correlator data not used in fit\n" << osBuffer.str();
    osBuffer.str( "" );
  }
  // Write theory and data EFFECTIVE MASS, with each data point in fit on a separate line
  ofs << "\n\n# Effective masses used in fit\n";
  Seq = 0;
  FitSeq = 0;
  ValWithEr<T> v;
  for( std::size_t m = 0; m < FitTimes.size(); ++m )
  {
    ++Seq;    // Skip past t=0
    ++FitSeq; // Skip past t=0
    for( int t = 1, idxSort = 1; t < ds(m, true).Nt(); ++t )
    {
      // Effective masses are for the point half-way between the measurements
      // Keeping track of 2 * the effective timeslice allows us to use integer arithmetic
      const int TwoT{ 2 * t - 1 }; // Average of this and previous timeslice
      const int TwoTFit{ idxSort >= FitTimes[m].size() ? std::numeric_limits<int>::max()
        : FitTimes[m][idxSort] + FitTimes[m][idxSort - 1] };
      const bool bShowFit{ TwoTFit <= TwoT };
      if( bShowFit )
      {
        // Show the fit
        const double dblT{ 0.5 * TwoTFit };
        const int DeltaT{ FitTimes[m][idxSort] - FitTimes[m][idxSort - 1] };
        const Scalar DeltaTInv{ static_cast<Scalar>( DeltaT == 1 ? 1. : ( 1. / DeltaT ) ) };
        ofs << "log " << Seq << Space << FitSeq << Space << m << Space << dblT << Space << dblT;
        for( std::size_t i = 0; i < Value.size(); ++i )
        {
          JackBoot<T> &ThD{ i == idxData ? FitInput : ModelPrediction };
          Buffer.resize( ThD.NumReplicas() ); // Shouldn't vary, but must be right size for check
          int bootCount{ 0 };
          for( int bootrep = 0; bootrep < ThD.NumReplicas(); ++bootrep )
          {
            Buffer[bootCount] = std::log( ThD(bootrep,FitSeq-1) / ThD(bootrep,FitSeq) ) * DeltaTInv;
            if( IsFinite( Buffer[bootCount] ) )
              ++bootCount;
          }
          v.Get( std::log( ThD(JackBoot<T>::idxCentral,FitSeq - 1)
                          / ThD(JackBoot<T>::idxCentral,FitSeq) ) * DeltaTInv, Buffer, bootCount,
                ThD.Seed == SeedWildcard );
          ofs << Space << v;
        }
        ofs << NewLine;
        ++idxSort;
        ++FitSeq;
      }
      if( TwoT != TwoTFit )
      {
        // Show the raw data
        const double dblT{ 0.5 * TwoT };
        osBuffer << "log " << Seq;
        if( bShowFit )
          osBuffer << ".5"; // Fit and data points out by 0.5 - make sure seq is different
        osBuffer << Space << -1 << Space << m << Space << dblT << Space << -1
                 << Space << ds.corr[m]->SummaryData( 2, t ) << Space << veZero << NewLine;
      }
      ++Seq;
    }
  }
  if( !osBuffer.str().empty() )
  {
    ofs << "\n\n# Effective masses not used in fit\n" << osBuffer.str();
    osBuffer.str( "" );
  }
}

template <typename T>
void Model<T>::ReorderOldFormat( int NumOps, int NumExponents, Vector<T> &v )
{
  if( NumExponents > 1 && NumOps > 1 )
  {
    std::vector<T> Buffer( NumOps * NumExponents );
    if( v.size < Buffer.size() )
      throw std::runtime_error( "Model<T>::ReorderOldFormat( Vector<T> & )" );
    for( int e = 0; e < NumExponents; ++e )
      for( int o = 0; o < NumOps; ++o )
        Buffer[o * NumExponents + e] = v[e * NumOps + o];
    for( std::size_t j = 0; j < Buffer.size(); ++j )
      v[j] = Buffer[j];
  }
}

template <typename T>
void Model<T>::ReorderOldFormat( int NumOps, int NumExponents, Matrix<T> &m )
{
  if( NumExponents > 1 && NumOps > 1 )
  {
    std::vector<T> Buffer( NumOps * NumExponents );
    if( m.size2 < Buffer.size() )
      throw std::runtime_error( "Model<T>::ReorderOldFormat( Matrix<T> & )" );
    for( std::size_t i = 0; i < m.size1; ++i )
    {
      for( int e = 0; e < NumExponents; ++e )
        for( int o = 0; o < NumOps; ++o )
          Buffer[o * NumExponents + e] = m(i, e * NumOps + o);
      for( std::size_t j = 0; j < Buffer.size(); ++j )
        m(i,j) = Buffer[j];
    }
  }
}

template <typename T> void Model<T>::CommonConstruct( const std::vector<std::string> &ExtraColumns )
{
  std::vector<std::string> Cols{ params.GetNames( Param::Type::All, false ) };
  Cols.reserve( Cols.size() + ExtraColumns.size() );
  for( const std::string &s : ExtraColumns )
    Cols.push_back( s );
  this->SetColumnNames( Cols );
}

template class Model<double>;
template class Model<float>;
template class Model<std::complex<double>>;
template class Model<std::complex<float>>;


END_COMMON_NAMESPACE
