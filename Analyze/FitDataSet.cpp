/*************************************************************************************
 
 Multi-exponential fits
 
 Source file: FitDataSet.cpp
 
 Copyright (C) 2020
 
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

// Everything to do with data in files that I'm fitting to

#include "MultiFit.hpp"

void DataSet::clear()
{
  NSamples = 0;
  Extent = 0;
  MinExponents = 0;
  MaxExponents = 0;
  corr.clear();
  FitTimes.clear();
  constFile.clear();
  ConstantNames.clear();
  ConstantNamesPerExp.clear();
  constMap.clear();
}

// Specify which times I'm fitting to, as a list of timeslices for each correlator
void DataSet::SetFitTimes( const std::vector<std::vector<int>> &fitTimes_ )
{
  if( fitTimes_.size() != corr.size() )
    throw std::runtime_error( std::to_string( fitTimes_.size() ) + " FitTimes but "
                             + std::to_string( corr.size() ) + " correlators" );
  std::vector<std::vector<int>> ft{ fitTimes_ };
  std::size_t extent_ = 0;
  for( int i = 0; i < ft.size(); ++i )
  {
    if( !ft[i].empty() )
    {
      std::sort( ft[i].begin(), ft[i].end() );
      if( ft[i][0]<0 || ft[i].back()>corr[i].Nt() || std::adjacent_find( ft[i].begin(), ft[i].end() )!=ft[i].end() )
        throw std::runtime_error( "FitTimes[" + std::to_string( i ) + "]=[" + std::to_string( ft[i][0] ) + "..."
                                 + std::to_string( ft[i].back() ) + "] invalid" );
      extent_ += ft[i].size();
    }
  }
  if( extent_ == 0 )
    throw std::runtime_error( "Fit range empty" );
  if( extent_ > std::numeric_limits<int>::max() )
    throw std::runtime_error( "Fit range stupidly big" );
  Extent = static_cast<int>( extent_ );
  FitTimes = std::move( ft );
}

void DataSet::SetFitTimes( int tMin, int tMax )
{
  const int extent_{ tMax - tMin + 1 };
  if( tMin < 0 || extent_ < 1 )
    throw std::runtime_error( "Fit range [" + std::to_string( tMin ) + ", " + std::to_string( tMax ) + "] invalid" );
  std::vector<std::vector<int>> ft{ corr.size() };
  for( int i = 0; i < corr.size(); ++i )
  {
    if( tMax >= corr[i].Nt() )
      throw std::runtime_error( "Fit range [" + std::to_string( tMin ) + ", " + std::to_string( tMax ) + "] invalid" );
    ft[i].reserve( extent_ );
    for( int j = tMin; j <= tMax; ++j )
      ft[i].push_back( j );
  }
  Extent = static_cast<int>( corr.size() * extent_ );
  FitTimes = std::move( ft );
}

// Get the data for the timeslices I'm fitting to
void DataSet::GetData( int idx, Vector &vResult ) const
{
  int dst{ 0 };
  vResult.resize( Extent );
  for( int f = 0; f < corr.size(); ++f )
  {
    const scalar * pSrc{ corr[f][idx] };
    for( int t : FitTimes[f] )
      vResult[dst++] = pSrc[t];
  }
}

// Get the constants from the appropriate timeslice
void DataSet::GetFixed( int idx, Vector &vResult, const std::vector<FixedParam> &Params ) const
{
  std::vector<const scalar *> Src( constFile.size() );
  for( int f = 0; f < constFile.size(); ++f )
    Src[f] = constFile[f][idx];
  for( const FixedParam &p : Params )
    vResult[p.idx] = Src[p.src.File][p.src.idx];
}

// Make the inverse of the error (i.e. inverse of square root of variance)
void DataSet::MakeInvErr( int idx, Vector &Var ) const
{
  Var.resize( Extent );
  for( int i = 0; i < Extent; ++i )
    Var[i] = 0;
  Vector Data( Extent );
  Vector Mean( Extent );
  GetData( idx, Mean );
  for( int replica = 0; replica < NSamples; ++replica )
  {
    GetData( replica, Data );
    for( int i = 0; i < Extent; ++i )
    {
      double z{ ( Data[i] - Mean[i] ) };
      Var[i] += z * z;
    }
  }
  for( int i = 0; i < Extent; ++i )
    Var[i] = std::sqrt( static_cast<scalar>( NSamples ) / Var[i] );
}

// Make covariance using mean from sample idx
void DataSet::MakeCovariance( int idx, Matrix &Covar ) const
{
  Covar.resize( Extent, Extent );
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
      Covar( i, j ) = 0;
  Vector Data( Extent );
  Vector Mean( Extent );
  GetData( idx, Mean );
  for( int replica = 0; replica < MaxSamples; ++replica )
  {
    GetData( replica, Data );
    for( int i = 0; i < Extent; ++i )
      Data[i] -= Mean[i];
    for( int i = 0; i < Extent; ++i )
      for( int j = 0; j <= i; ++j )
        Covar( i, j ) += Data[i] * Data[j];
  }
  for( int i = 0; i < Extent; ++i )
    for( int j = 0; j <= i; ++j )
    {
      const scalar z{ Covar( i, j ) / MaxSamples };
      Covar( i, j ) = z;
      if( i != j )
        Covar( j, i ) = z;
    }
}

// Write covariance matrix to file
void DataSet::SaveCovariance(const std::string FileName, const Matrix &Covar, const std::vector<ModelPtr> *pModels)const
{
  // Header describing what the covariance matrix is for and how to plot it with gnuplot
  assert( Covar.size1 == Extent && Covar.size1 == Covar.size2 && "Build covar before saving it" );
  std::ofstream s{ FileName };
  s << "# Correlation matrix\n# Files: " << corr.size() << Common::NewLine;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    s << "# File" << f << ": " << corr[f].Name_.NameNoExt;
    // Say which operators are in each file
    if( pModels )
    {
      if( (*pModels).size() != corr.size() )
        throw std::runtime_error( "ModelSet doesn't match DataSet" );
      for( const std::string &op : (*pModels)[f]->ParamNames )
        s << ", " << op;
      if( !(*pModels)[f]->ParamNamesPerExp.empty() )
      {
        s << " Per E:";
        for( std::size_t i = 0; i < (*pModels)[f]->ParamNamesPerExp.size(); ++i )
          s << ( i ? ", " : Common::Space ) << (*pModels)[f]->ParamNamesPerExp[i];
      }
    }
    s << Common::NewLine;
    s << "# Times" << f << ":";
    for( int t : FitTimes[f] )
      s << Common::Space << t;
    s << Common::NewLine;
  }
  s << "# gnuplot: plot '" << FileName << "' matrix columnheaders rowheaders with image pixels\n" << Extent;
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    const std::string FileLetter( 1, 'A' + f );
    for( int t : FitTimes[f] )
    {
      s << Common::Space;
      if( pModels )
      {
        for( std::size_t op = 1; op < (*pModels)[f]->ParamNamesPerExp.size(); ++op )
          s << (*pModels)[f]->ParamNamesPerExp[op];
      }
      else
        s << FileLetter;
      s << t;
    }
  }
  s << Common::NewLine;
  // Now print the actual covariance matrix
  int i{ 0 };
  for( std::size_t f = 0; f < corr.size(); ++f )
  {
    const std::string FileLetter( 1, 'A' + f );
    for( int t : FitTimes[f] )
    {
      if( pModels )
      {
        for( std::size_t op = 1; op < (*pModels)[f]->ParamNamesPerExp.size(); ++op )
          s << (*pModels)[f]->ParamNamesPerExp[op];
      }
      else
        s << FileLetter;
      s << t;
      for( int j = 0; j < Extent; ++j )
        s << Common::Space << Covar( i, j );
      s << Common::NewLine;
      ++i;
    }
  }
}

// Add a constant to my list of known constants - make sure it isn't already there
void DataSet::AddConstant( const std::string &Name, std::size_t File, std::size_t idx )
{
  if( constMap.find( Name ) != constMap.end() )
    throw std::runtime_error( "Constant \"" + Name + "\" loaded from multiple model files" );
  constMap.insert( { Name, ConstantSource( File, idx ) } );
}

void DataSet::AddConstant( const std::string &Name, std::size_t File, std::size_t idx, int e )
{
  std::string s{ Name };
  s.append( std::to_string( e ) );
  AddConstant( s, File, idx );
}

// Load a correlator or a model a file. If it's a correlator, there may be leftover arguments
void DataSet::LoadFile( const std::string &sFileName, std::vector<std::string> &OpNames,
                        std::vector<std::string> &ModelArgs, const std::string &Args )
{
  // First string is Filename
  Common::FileNameAtt Att( sFileName );
  if( Common::EqualIgnoreCase( Att.Type, Common::sFold ) )
  {
    // This is a correlator - load it
    Att.ParseOpNames( OpNames );
    const std::size_t i{ corr.size() };
    corr.emplace_back();
    corr[i].SetName( std::move( Att ) );
    corr[i].Read( "  " );
    // See whether this correlator is compatible with prior correlators
    if( i )
    {
      corr[0].IsCompatible( corr[i], &NSamples );
      if( MaxSamples > corr[i].NumSamples() )
        MaxSamples = corr[i].NumSamples();
    }
    else
      MaxSamples = corr[i].NumSamples();
    ModelArgs.push_back( Args );
  }
  else
  {
    // Break trailing arguments for this file into an array of strings
    std::vector<std::string> vThisArg{ Common::ArrayFromString( Args ) };
    // This is a pre-built model (i.e. the result of a previous fit)
    const std::size_t i{ constFile.size() };
    constFile.emplace_back();
    constFile[i].SetName( std::move( Att ) );
    constFile[i].Read( "  " );
    // Keep track of minimum number of replicas across all files
    if( NSamples == 0 )
      NSamples = constFile[i].NumSamples();
    else if( NSamples > constFile[i].NumSamples() )
      NSamples = constFile[i].NumSamples();
    // Keep track of minimum and maximum number of exponents across all files
    if( i )
    {
      if( MinExponents > constFile[i].NumExponents )
        MinExponents = constFile[i].NumExponents;
      if( MaxExponents < constFile[i].NumExponents )
        MinExponents = constFile[i].NumExponents;
    }
    else
    {
      MinExponents = constFile[i].NumExponents;
      MaxExponents = constFile[i].NumExponents;
    }
    // Now see which parameters we want to read from the model
    const vString &ColumnNames{ constFile[i].GetColumnNames() };
    const vString &OpNames{ constFile[i].OpNames };
    if( vThisArg.empty() )
    {
      // Load every constant in this file
      for( int j = 0; j < ColumnNames.size(); ++j )
      {
        // Add a map back to this specific parameter (making sure not already present)
        AddConstant( ColumnNames[j], i, j );
        // Now take best stab at whether this should be per exponent or individual
        std::string s{ ColumnNames[j] };
        int Exp;
        if( Common::ExtractTrailing( s, Exp )
           && ( Common::EqualIgnoreCase( E, s ) || Common::IndexIgnoreCase( OpNames, s ) != OpNames.size() ) )
        {
          ConstantNamesPerExp[s];
        }
        else
          ConstantNames[ColumnNames[j]];
      }
    }
    else
    {
      // Load only those constants specifically asked for
      for( const std::string &ThisArg : vThisArg )
      {
        std::vector<std::string> vPar{ Common::ArrayFromString( ThisArg, "=" ) };
        if( vPar.empty() || vPar[0].empty() || vPar.size() > 2 || ( vPar.size() > 1 && vPar[1].empty() ) )
          throw std::runtime_error( "Cannot interpret model parameter string \"" + Args + "\"" );
        const std::string &vLookFor{ vPar.back() };
        const std::string &vLoadAs{ vPar[0] };
        // Have we asked for a per-exponent constant (which includes energies)
        if( Common::EqualIgnoreCase( E, vLookFor ) || Common::IndexIgnoreCase( OpNames, vLookFor ) != OpNames.size() )
        {
          for( int e = 0; e < constFile[i].NumExponents; ++e )
            AddConstant( vLoadAs, i, constFile[i].GetColumnIndex( vLookFor, e ), e );
          ConstantNamesPerExp[vLoadAs];
        }
        else
        {
          AddConstant( vLoadAs, i, constFile[i].GetColumnIndex( vLookFor ) );
          ConstantNames[vLoadAs];
        }
      }
    }
  }
}

// Sort the operator names and renumber all the loaded correlators referring to them
void DataSet::SortOpNames( std::vector<std::string> &OpNames )
{
  int NumOps{ static_cast<int>( OpNames.size() ) };
  if( OpNames.size() > 1 )
  {
    // Sort the names
    UniqueNames OpSorted;
    for( int i = 0; i < NumOps; ++i )
      OpSorted.emplace( std::move( OpNames[i] ), i );
    // Extract the sorted names and indices (to renumber operators in correlator names)
    std::vector<std::string> SortedNames;
    std::vector<int> SortIndex( NumOps );
    SortedNames.reserve( NumOps );
    int idx{ 0 };
    for( UniqueNames::iterator it = OpSorted.begin(); it != OpSorted.end(); ++it )
    {
      SortedNames.emplace_back( it->first );
      SortIndex[it->second] = idx++;
    }
    // Renumber the operators and save the sorted operator names
    for( auto &f : corr )
      for( int i = 0; i < f.Name_.op.size(); ++i )
        f.Name_.op[i] = SortIndex[f.Name_.op[i]];
    OpNames = SortedNames;
  }
}
