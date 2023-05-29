// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeBoxReader.h"

#include "vtkDirectory.h"
#include "vtkDoubleArray.h"
#include "vtkGenericAttributeCollection.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"

#include "vtkDataRecords.h"
#include "vtkShoeAttribute.h"
#include "vtkShoeBox.h"
#include "vtkShoeBoxAlgorithm.h"
#include "vtkShoeBoxP.h"
#include "vtkShoeCellRecord.h"
#include "vtkShoeCellSpecies.h"
#include "vtkShoeEnums.h"
#include "vtkShoeOrderTuple.h"

//FIXME: Use NumberOfCellsByDimension to fill in appropriate blanks
//FIXME: Get rid of GeometryIndex
//FIXME: Binary data. It could *almost* work...
//FIXME: Cell arrays
//FIXME: More summary information in summary file (req'd for next 3 fixme items)
//FIXME: Time steps
//FIXME: Multiblock output?
//FIXME: Sharing attributes among meshes

vtkStandardNewMacro(vtkShoeBoxReader);
vtkCxxRevisionMacro(vtkShoeBoxReader,"$Revision: 6040 $");

vtkShoeBoxReader::vtkShoeBoxReader()
{
  this->FileName = 0;
  this->Title = 0;
  this->Description = 0;
  this->TimeStepRange[0] = this->TimeStepRange[1] = 0;
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  this->GeometryIndex = 0;
}

vtkShoeBoxReader::~vtkShoeBoxReader()
{
  this->SetFileName( 0 );
  this->SetTitle( 0 );
  this->SetDescription( 0 );
}

void vtkShoeBoxReader::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "FileName:         " << (this->FileName ? this->FileName : "(none)") << vtkstd::endl;
  os << indent << "Title:            " << (this->Title ? this->Title : "(none)") << vtkstd::endl;
  os << indent << "Description:      " << (this->Description ? this->Description : "(none)") << vtkstd::endl;
  os << indent << "NumberOfPoints:   " << this->NumberOfPoints << vtkstd::endl;
  os << indent << "NumberOfDofNodes: " << this->NumberOfDofNodes << vtkstd::endl;
  os << indent << "NumberOfCells:    " << this->NumberOfCells << vtkstd::endl;
  os << indent << "  dimension 0:    " << this->NumberOfCellsByDimension[0] << vtkstd::endl;
  os << indent << "  dimension 1:    " << this->NumberOfCellsByDimension[1] << vtkstd::endl;
  os << indent << "  dimension 2:    " << this->NumberOfCellsByDimension[2] << vtkstd::endl;
  os << indent << "  dimension 3:    " << this->NumberOfCellsByDimension[3] << vtkstd::endl;
  os << indent << "GeometryIndex:    " << this->GeometryIndex << vtkstd::endl;
  os << indent << "TimeStepRange:    (" << this->TimeStepRange[0] << ", " << this->TimeStepRange[1] << ")" << vtkstd::endl;
}

vtkIdType vtkShoeBoxReader::GetNumberOfCells( int dim )
{
  assert(dim >= 0 && dim <= 3);
  return this->NumberOfCellsByDimension[ dim ];
}

vtkIdType vtkShoeBoxReader::GetNumberOfTimeSteps() const
{
  return this->TimeValues.size();
}

double vtkShoeBoxReader::GetTimeValue( int t ) const
{
  if ( t < 0 || t >= this->GetNumberOfTimeSteps() )
    {
    return 0.;
    }
  return this->TimeValues[t];
}

int vtkShoeBoxReader::GetNumberOfPointArrays()
{
  return this->PointArrayNames.size();
}

const char* vtkShoeBoxReader::GetPointArrayName( int i )
{
  if ( i < 0 || i >= this->GetNumberOfPointArrays() )
    return 0;

  return this->PointArrayNames[i].c_str();
}

void vtkShoeBoxReader::SetPointArrayStatus( int i, int status )
{
  if ( i < 0 || i >= this->GetNumberOfPointArrays() )
    return;
  if ( status ^ this->PointArrayStatus[i] )
    {
    this->PointArrayStatus[i] = status;
    this->Modified();
    }
}

int vtkShoeBoxReader::GetPointArrayStatus( int i )
{
  if ( i < 0 || i >= this->GetNumberOfPointArrays() )
    return 0;
  return this->PointArrayStatus[ i ];
}

int vtkShoeBoxReader::GetNumberOfCellArrays()
{
  return this->CellArrayNames.size();
}

const char* vtkShoeBoxReader::GetCellArrayName( int i )
{
  if ( i < 0 || i >= this->GetNumberOfCellArrays() )
    return 0;

  return this->CellArrayNames[i].c_str();
}

void vtkShoeBoxReader::SetCellArrayStatus( int i, int status )
{
  if ( i < 0 || i >= this->GetNumberOfCellArrays() )
    return;
  if ( status ^ this->CellArrayStatus[i] )
    {
    this->CellArrayStatus[i] = status;
    this->Modified();
    }
}

int vtkShoeBoxReader::GetCellArrayStatus( int i )
{
  if ( i < 0 || i >= this->GetNumberOfCellArrays() )
    return 0;
  return this->CellArrayStatus[ i ];
}

int vtkShoeBoxReader::GetNumberOfComponentsInCellArray( int a )
{
  if ( a < 0 || a >= this->GetNumberOfCellArrays() )
    return 0;
  return this->CellArraySize[a];
}

int vtkShoeBoxReader::GetNumberOfComponentsInPointArray( int a )
{
  if ( a < 0 || a >= this->GetNumberOfPointArrays() )
    return 0;
  return this->PointArraySize[a];
}

int vtkShoeBoxReader::RequestInformation(
  vtkInformation* vtkNotUsed(i), vtkInformationVector** vtkNotUsed(ii), vtkInformationVector* vtkNotUsed(o) )
{
  return this->ReadSummary();
}

int vtkShoeBoxReader::RequestData(
  vtkInformation* vtkNotUsed(request), vtkInformationVector** vtkNotUsed(inputVector), vtkInformationVector* outputVector )
{
  // get the info objects
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  // get the ouptut
  this->Output = vtkShoeBox::SafeDownCast( outInfo->Get( vtkDataObject::DATA_OBJECT() ) );
  // read in the mesh
  this->ReadConnectivity();
  this->ReadPointLinks();
  this->ReadDOFLinks();
  this->ReadAttributes();
  this->ReadSpecies();
  this->ReadCells();

  return 1;
}


int vtkShoeBoxReader::ReadSummary()
{
  vtkstd::string fn( this->FileName );
  fn += "/Summary";
  ifstream summary;
  summary.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! summary.good() )
    {
    vtkWarningMacro("Unable to read summary information.");
    return 0;
    }

  vtkstd::string word;
  summary >> word;
  if ( word != "ASCII" )
    {
    vtkErrorMacro("Can't read binary files yet");
    return 0;
    }
  summary >> this->NumberOfPoints;
  summary >> this->NumberOfCells;
  summary >> this->NumberOfDofNodes;
  summary >> this->NumberOfCellsByDimension[0];
  summary >> this->NumberOfCellsByDimension[1];
  summary >> this->NumberOfCellsByDimension[2];
  summary >> this->NumberOfCellsByDimension[3];
  summary >> this->GeometryIndex;
  summary.close();

  vtkstd::string attrDirName( this->FileName );
  attrDirName += "/Attributes";
  vtkDirectory* attrDir = vtkDirectory::New();
  attrDir->Open( attrDirName.c_str() );
  int nf = attrDir->GetNumberOfFiles();
  for ( int i = 0; i < nf; ++i )
    {
    const char* aName = attrDir->GetFile( i );
    if ( aName && aName[0] != '.' )
      {
      fn = attrDirName + "/";
      fn += aName;
      fn += "/Info";
      summary.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
      if ( ! summary.good() )
        {
        vtkWarningMacro("Unable to open attribute info file \"" << fn.c_str() << "\"");
        continue;
        }

      int nc;
      int btype;
      vtkstd::vector<double> bds;
      summary >> nc;
      if ( nc <= 0 )
        {
        vtkWarningMacro("Attribute \"" << aName << "\" has an invalid number of components (" << nc << ")");
        continue;
        }

      summary >> word;
      if ( word == "None" )
        {
        btype = VTK_RANGE_NONE;
        }
      else if ( word == "Sloppy" )
        {
        btype = VTK_RANGE_SLOPPY;
        }
      else if ( word == "Proper" )
        {
        btype = VTK_RANGE_PROPER;
        }
      else if ( word == "Tight" )
        {
        btype = VTK_RANGE_TIGHT;
        }
      else
        {
        vtkWarningMacro("Attribute range not specified correctly.");
        continue;
        }
      if ( btype != VTK_RANGE_NONE )
        {
        double btmp;
        for ( int b = 0; b < 2 * nc; ++b )
          {
          summary >> btmp;
          bds.push_back( btmp );
          }
        }
      summary.close();
      this->PointArrayNames.push_back( aName );
      this->PointArrayStatus.push_back( 1 );
      this->PointArraySize.push_back( nc );
      this->PointArrayBoundsType.push_back( btype );
      this->PointArrayBounds.push_back( bds );
      }
    }

  return 1;
}

int vtkShoeBoxReader::ReadAttributes()
{
  for ( vtkIdType a = 0; a < this->GetNumberOfPointArrays(); ++a )
    {
    // Always read geometry attribute
    if ( this->GetPointArrayStatus( a ) || (this->PointArrayNames[a] == "Geometry") )
      {
      if ( ! this->ReadAttribute( a ) )
        {
        return 0;
        }
      }
    }
  return 1;
}

int vtkShoeBoxReader::ReadCells()
{
  vtkstd::string fn( this->FileName );
  fn += "/Cells";
  ifstream cellFile;
  cellFile.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! cellFile.good() )
    {
    vtkWarningMacro("Unable to read cell information.");
    return 0;
    }

  vtkstd::string word;
  vtkIdType unused;

  cellFile >> unused;
  vtkIdType total = unused + this->NumberOfCells;
  vtkIdType cell;
  for ( cell = 0; cell < total; ++cell )
    {
    this->Output->CellRecs->Cells.grab();
    }
  for ( cell = 0; cell < unused; ++cell )
    {
    vtkIdType tmp;
    cellFile >> tmp;
    this->Output->CellRecs->Cells.free( tmp );

    }
  for ( freelist<vtkShoeCellRecord,vtkIdType>::iterator it = this->Output->CellRecs->Cells.begin();
        it != this->Output->CellRecs->Cells.end(); ++it )
    {
    vtkIdType tmp;

    cellFile >> tmp;
    it->Species = vtkShoeCellSpecies::GetSpeciesById( tmp );
    cellFile >> it->Offset;
    cellFile >> it->NodePermutations;
    }

  return 1;
}

int vtkShoeBoxReader::ReadConnectivity()
{
  vtkstd::string fn( this->FileName );
  fn += "/Connectivity";
  ifstream connectivity;
  connectivity.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! connectivity.good() )
    {
    vtkWarningMacro("Unable to read connectivity information.");
    return 0;
    }

  return this->ReadDataRecordsIdType( connectivity, this->Output->GetConnectivity() );
}

int vtkShoeBoxReader::ReadDOFLinks()
{
  vtkstd::string fn( this->FileName );
  fn += "/DofLinks";
  ifstream dofLinks;
  dofLinks.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! dofLinks.good() )
    {
    vtkWarningMacro("Unable to read DOF links information.");
    return 0;
    }

  return this->ReadDataRecordsIdType( dofLinks, this->Output->GetDOFLinks() );
}

int vtkShoeBoxReader::ReadPointLinks()
{
  vtkstd::string fn( this->FileName );
  fn += "/PointLinks";
  ifstream pointLinks;
  pointLinks.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! pointLinks.good() )
    {
    vtkWarningMacro("Unable to read point links information.");
    return 0;
    }

  return this->ReadDataRecordsIdType( pointLinks, this->Output->GetPointLinks() );
}

int vtkShoeBoxReader::ReadSpecies()
{
  this->SpeciesSwizzle.clear();

  vtkstd::string fn( this->FileName );
  fn += "/Species";
  ifstream speciesFile;
  speciesFile.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! speciesFile.good() )
    {
    vtkWarningMacro("Unable to read species information.");
    return 0;
    }

  vtkstd::string word;
  vtkstd::vector<vtkPolyInterpolant> desiredInterpolants;
  vtkstd::vector<vtkPolyProductSpace> desiredProdSpaces;
  vtkstd::vector<vtkShoeOrderTuple> desiredOrders;
  vtkstd::vector<int> attributeIds;
  int nconstraints;
  int ns; // number of species
  vtkShoeCellGenus genus;
  speciesFile >> ns;
  for ( int s = 0; s < ns; ++s )
    {
    vtkIdType swizzle;
    vtkIdType na;
    int shape;

    speciesFile >> swizzle;
    speciesFile >> shape;
    speciesFile >> na;
    genus.Shape = shoe::CellShape( shape );

    desiredInterpolants.clear();
    desiredProdSpaces.clear();
    desiredOrders.clear();
    attributeIds.clear();
    nconstraints = 0;
    for ( int a = 0; a < na; ++a )
      {
      speciesFile >> word;
      vtkIdType aid; // id of attribute
      aid = this->Output->GetAttributes()->FindAttribute( word.c_str() );
      if ( aid < 0 )
        {
        vtkWarningMacro("Cannot create species because attribute \"" << word.c_str() << "\" isn't defined on the mesh.");
        return 0;
        }
      if ( this->GetPointArrayStatus(aid) )
        {
        nconstraints++;
        attributeIds.push_back( aid );

        speciesFile >> word;
        if ( word == "Lagrange" )
          desiredInterpolants.push_back( Lagrange );
        else if ( word == "Legendre" )
          desiredInterpolants.push_back( Legendre );
        else if ( word == "Hermite" )
          desiredInterpolants.push_back( Hermite );
        else if ( word == "BSpline" )
          desiredInterpolants.push_back( BSpline );
        else
          {
          vtkWarningMacro("Bad interpolant keyword \"" << word.c_str() << "\" reading species" << s << ".");
          return 0;
          }

        speciesFile >> word;
        if ( word == "Tensor" )
          desiredProdSpaces.push_back( Tensor );
        else if ( word == "MaxTotalDegree" )
          desiredProdSpaces.push_back( MaxTotalDegree );
        else if ( word == "TruncatedTotalDegree" )
          desiredProdSpaces.push_back( TruncatedTotalDegree );
        else
          {
          vtkWarningMacro("Bad product space keyword \"" << word.c_str() << "\" reading species" << s << ".");
          return 0;
          }

        vtkShoeOrderTuple o;
        speciesFile >> o.Order[0];
        speciesFile >> o.Order[1];
        speciesFile >> o.Order[2];
        desiredOrders.push_back( o );
        }
      }
    vtkShoeCellSpecies* species = vtkShoeCellSpecies::FindOrCreate( genus,
      &desiredInterpolants[0], &desiredProdSpaces[0], &desiredOrders[0], &attributeIds[0], nconstraints,
      this->Output->GetAttributes() );
    this->SpeciesSwizzle[ swizzle ] = species->GetId();
    }
  return 1;
}

int vtkShoeBoxReader::ReadAttribute( vtkIdType a )
{
  vtkShoeAttribute* attrib = vtkShoeAttribute::New();
  attrib->SetName( this->GetPointArrayName( a ) );
  attrib->SetNumberOfComponents( this->PointArraySize[a] );
  vtkDataRecords* dof = vtkDataRecords::New();
  vtkDoubleArray* pts = vtkDoubleArray::New();

  vtkstd::string fn( this->FileName );
  fn += "/Attributes/" + this->PointArrayNames[a] + "/Dof";
  ifstream dataFile;
  dataFile.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! dataFile.good() )
    {
    vtkWarningMacro("Unable to read DOF information.");
    dof->Delete();
    pts->Delete();
    attrib->Delete();
    return 0;
    }

  if ( ! this->ReadDataRecordsDouble( dataFile, dof ) )
    {
    vtkWarningMacro("Unable to read DOF data.");
    dof->Delete();
    pts->Delete();
    attrib->Delete();
    return 0;
    }
  dataFile.close();

  fn = this->FileName;
  fn += "/Attributes/" + this->PointArrayNames[a] + "/Points";
  ifstream ptsFile;
  ptsFile.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  if ( ! ptsFile.good() )
    {
    vtkWarningMacro("Unable to read point information.");
    dof->Delete();
    pts->Delete();
    attrib->Delete();
    return 0;
    }

  if ( ! this->ReadDoubleArray( ptsFile, pts ) )
    {
    vtkWarningMacro("Unable to read point data.");
    dof->Delete();
    pts->Delete();
    attrib->Delete();
    return 0;
    }
  ptsFile.close();

  attrib->SetPointData( pts );
  attrib->SetDOFData( dof );
  dof->FastDelete();
  pts->FastDelete();
  if ( this->PointArrayBoundsType[a] != VTK_RANGE_NONE )
    {
    attrib->SetRange( &this->PointArrayBounds[a][0], this->PointArrayBoundsType[a] );
    }

  if ( this->PointArrayNames[a] == "Geometry" )
    {
    this->Output->SetGeometry( attrib );
    }
  else
    {
    this->Output->GetAttributes()->InsertNextAttribute( attrib );
    }
  attrib->FastDelete();

  fn = this->FileName;
  fn += "/Attributes/" + this->PointArrayNames[a] + "/CriticalPoints";
  ifstream cptsFile;
  vtkDataRecords** critPts = attrib->GetCriticalPoints();
  cptsFile.open( fn.c_str(), this->Binary ? ifstream::in | ifstream::binary : ifstream::in );
  int c = 0;
  bool haveCP = false;
  while ( cptsFile.good() && c < attrib->GetNumberOfComponents() )
    {
    haveCP = true;
    if ( ! this->ReadDataRecordsDouble( cptsFile, critPts[c] ) )
      {
      haveCP = false;
      vtkWarningMacro("Unable to read critical points, though the file exists.");
      break;
      }
    ++c;
    }
  cptsFile.close();
  if ( haveCP )
    {
    attrib->SetCriticalPointsDirty( 0 );
    }

  return 1;
}

int vtkShoeBoxReader::ReadDataRecordsIdType( istream& in, vtkDataRecords* dr )
{
  vtkIdType nr;
  vtkIdType nc;
  in >> nr;
  in >> nc;
  vtkstd::vector<vtkIdType> tuple;
  tuple.reserve( nc );
  dr->SetNumberOfComponents( nc );
  dr->SetNumberOfRecords( nr );
  for ( vtkIdType r = 0; r < nr; ++r )
    {
    vtkIdType nt;
    in >> nt;
    for ( vtkIdType t = 0; t < nt; ++t )
      {
      for ( int i = 0; i < nc; ++i )
        {
        in >> tuple[i];
        }
      dr->InsertNextTuple( r, &tuple[0] );
      }
    }
  return 1;
}

int vtkShoeBoxReader::ReadDataRecordsDouble( istream& in, vtkDataRecords* dr )
{
  vtkIdType nr;
  vtkIdType nc;
  in >> nr;
  in >> nc;
  vtkstd::vector<double> tuple;
  tuple.reserve( nc );
  dr->SetNumberOfComponents( nc );
  dr->SetNumberOfRecords( nr );
  for ( vtkIdType r = 0; r < nr; ++r )
    {
    vtkIdType nt;
    in >> nt;
    for ( vtkIdType t = 0; t < nt; ++t )
      {
      for ( int i = 0; i < nc; ++i )
        {
        in >> tuple[i];
        }
      dr->InsertNextTuple( r, &tuple[0] );
      }
    }
  return 1;
}

int vtkShoeBoxReader::ReadDoubleArray( istream& in, vtkDoubleArray* da )
{
  vtkIdType nt;
  vtkIdType nc;
  in >> nt;
  in >> nc;
  vtkstd::vector<double> tuple;
  tuple.reserve( nc );
  da->SetNumberOfComponents( nc );
  da->SetNumberOfTuples( nt );
  for ( vtkIdType t = 0; t < nt; ++t )
    {
    for ( int i = 0; i < nc; ++i )
      {
      in >> tuple[i];
      }
    da->InsertNextTuple( &tuple[0] );
    }
  return 1;
}


