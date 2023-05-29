/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkObjectFactory.h>
#include <vtkCellType.h>
#include <vtkPoints.h>
#include <vtkDataSetAttributes.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>

#include <vtkShoeMeshPolygonizer.h>
#include <vtkShoeMesh.h>
#include <vtkShoeMeshIterator.h>
#include <vtkCellOps.h>
#include <vtkAdaptiveTessellator.h>
#include <vtkShoeMeshChordLengthSubdivision.h>
#include <vtkShoeMeshViewDependentSubdivision.h>
#include <vtkSubdivisionAlgorithm.h>



vtkCxxRevisionMacro(vtkShoeMeshPolygonizer, "$Revision: 4814 $");
vtkStandardNewMacro(vtkShoeMeshPolygonizer);

static double empty_normal[3] = { 0., 0., 0. };

void vtkShoeMeshPolygonizer::AddATetrahedron( const double* a, const double* b, const double* c, const double* d,
                                              vtkSubdivisionAlgorithm*, void* pd, const void* )
{
  vtkShoeMeshPolygonizer* self = (vtkShoeMeshPolygonizer*) pd;
  self->OutputTetrahedron( a, b, c, d );
}

void vtkShoeMeshPolygonizer::OutputTetrahedron( const double* a, const double* b, const double* c, const double* d )
{
  vtkIdType cellIds[4];

  cellIds[0] = this->OutputPoints->InsertNextPoint( a );
  cellIds[1] = this->OutputPoints->InsertNextPoint( b );
  cellIds[2] = this->OutputPoints->InsertNextPoint( c );
  cellIds[3] = this->OutputPoints->InsertNextPoint( d );

  this->OutputMesh->InsertNextCell( VTK_TETRA, 4, cellIds );

  const int* off = this->Subdivider->GetFieldOffsets();
  vtkDataArray** att = this->OutputAttributes;
  if ( this->GenerateNormals )
    {
    (*att)->InsertTuple( cellIds[0], empty_normal );
    (*att)->InsertTuple( cellIds[1], empty_normal );
    (*att)->InsertTuple( cellIds[2], empty_normal );
    (*att)->InsertTuple( cellIds[3], empty_normal );

    ++att;
    }

  for ( int at=0; at<this->Subdivider->GetNumberOfFields(); ++at, ++att, ++off )
    {
    (*att)->InsertTuple( cellIds[0], a + *off );
    (*att)->InsertTuple( cellIds[1], b + *off );
    (*att)->InsertTuple( cellIds[2], c + *off );
    (*att)->InsertTuple( cellIds[3], d + *off );
    }
}

void vtkShoeMeshPolygonizer::AddATriangle( const double* a, const double* b, const double* c,
                                           vtkSubdivisionAlgorithm*, void* pd, const void* )
{
  vtkShoeMeshPolygonizer* self = (vtkShoeMeshPolygonizer*) pd;
  self->OutputTriangle( a, b, c );
}

void vtkShoeMeshPolygonizer::OutputTriangle( const double* a, const double* b, const double* c )
{
  vtkIdType cellIds[3];

  cellIds[0] = this->OutputPoints->InsertNextPoint( a );
  cellIds[1] = this->OutputPoints->InsertNextPoint( b );
  cellIds[2] = this->OutputPoints->InsertNextPoint( c );

  this->OutputMesh->InsertNextCell( VTK_TRIANGLE, 3, cellIds );

  const int* off = this->Subdivider->GetFieldOffsets();
  vtkDataArray** att = this->OutputAttributes;
  if ( this->GenerateNormals )
    {
    double norm[3];

    this->CurrentCell.GetCellOps()->EvaluateNormalOnFace( norm, this->CurrentCell, &a[3], 6 );
    (*att)->InsertTuple( cellIds[0], norm );

    this->CurrentCell.GetCellOps()->EvaluateNormalOnFace( norm, this->CurrentCell, &b[3], 6 );
    (*att)->InsertTuple( cellIds[1], norm );

    this->CurrentCell.GetCellOps()->EvaluateNormalOnFace( norm, this->CurrentCell, &c[3], 6 );
    (*att)->InsertTuple( cellIds[2], norm );

    ++att;
    }

  for ( int at=0; at<this->Subdivider->GetNumberOfFields(); ++at, ++att, ++off )
    {
    (*att)->InsertTuple( cellIds[0], a + *off );
    (*att)->InsertTuple( cellIds[1], b + *off );
    (*att)->InsertTuple( cellIds[2], c + *off );
    }
}

void vtkShoeMeshPolygonizer::AddALine( const double* a, const double* b,
                                       vtkSubdivisionAlgorithm*, void* pd, const void* )
{
  vtkShoeMeshPolygonizer* self = (vtkShoeMeshPolygonizer*) pd;
  self->OutputLine( a, b );
}

void vtkShoeMeshPolygonizer::OutputLine( const double* a, const double* b )
{
  vtkIdType cellIds[2];

  cellIds[0] = this->OutputPoints->InsertNextPoint( a );
  cellIds[1] = this->OutputPoints->InsertNextPoint( b );

  this->OutputMesh->InsertNextCell( VTK_LINE, 2, cellIds );

  const int* off = this->Subdivider->GetFieldOffsets();
  vtkDataArray** att = this->OutputAttributes;
    if ( this->GenerateNormals )
    {
    (*att)->InsertTuple( cellIds[0], empty_normal );
    (*att)->InsertTuple( cellIds[1], empty_normal );

    ++att;
    }

  for ( int at=0; at<this->Subdivider->GetNumberOfFields(); ++at, ++att, ++off )
    {
    (*att)->InsertTuple( cellIds[0], a + *off );
    (*att)->InsertTuple( cellIds[1], b + *off );
    }
}

vtkShoeMeshPolygonizer::vtkShoeMeshPolygonizer()
  : Tessellator( 0 ), Subdivider( 0 ), GenerateNormals( true )
{
  this->SetTessellator( vtkAdaptiveTessellator::New() );
  this->SetSubdivider( vtkShoeMeshChordLengthSubdivision::New() );
}

vtkShoeMeshPolygonizer::~vtkShoeMeshPolygonizer()
{
  this->SetSubdivider( 0 );
  this->SetTessellator( 0 );
}

void vtkShoeMeshPolygonizer::PrintSelf(vtkstd::ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Tessellator: " << Tessellator << vtkstd::endl
     << indent << "Subdivider: " << Subdivider << " (" << Subdivider->GetClassName() << ")" << vtkstd::endl
     << indent << "GenerateNormals: " << (this->GenerateNormals ? "On" : "Off") << vtkstd::endl;
}

void vtkShoeMeshPolygonizer::SetTessellator( vtkAdaptiveTessellator* t )
{
  if ( this->Tessellator == t )
    return;

  if ( this->Tessellator )
    this->Tessellator->UnRegister( this );

  this->Tessellator = t;

  if ( this->Tessellator )
    {
    this->Tessellator->Register( this );
    this->Tessellator->SetSubdivisionAlgorithm( this->Subdivider );
    }

  this->Modified();
}

void vtkShoeMeshPolygonizer::SetSubdivider( vtkShoeMeshSubdivisionAlgorithm* s )
{
  if ( this->Subdivider == s )
    return;

  if ( this->Subdivider )
    this->Subdivider->UnRegister( this );

  this->Subdivider = s;

  if ( this->Subdivider )
    this->Subdivider->Register( this );

  if ( this->Tessellator )
    this->Tessellator->SetSubdivisionAlgorithm( this->Subdivider );

  this->Modified();
}

void vtkShoeMeshPolygonizer::SetGenerateNormals( bool n )
{
  if ( n == this->GenerateNormals )
    return;

  this->GenerateNormals = n;
  this->Modified();
}

unsigned long vtkShoeMeshPolygonizer::GetMTime()
{
  unsigned long mt = this->MTime;
  unsigned long tmp;

  if ( this->Tessellator )
    {
    tmp = this->Tessellator->GetMTime();
    if ( tmp > mt )
      mt = tmp;
    }

  if ( this->Subdivider )
    {
    tmp = this->Subdivider->GetMTime();
    if ( tmp > mt )
      mt = tmp;
    }

  return mt;
}

void vtkShoeMeshPolygonizer::SetupOutput()
{
  this->OutputMesh = this->GetOutput(); // avoid doing all the stupid checks on NumberOfOutputs for every triangle/line.
  this->OutputMesh->Reset();
  this->OutputMesh->Allocate(0,0);

  if ( ! (this->OutputPoints = OutputMesh->GetPoints()) )
    {
    this->OutputPoints = vtkPoints::New();
    this->OutputMesh->SetPoints( this->OutputPoints );
    this->OutputPoints->Delete();
    }

  int maxNumComponents = 0;

  // This returns the id numbers of arrays that are default scalars, vectors, normals, texture coords, and tensors.
  // These are the fields that will be interpolated and passed on to the output mesh, with the exception of normals,
  // which are specified by the GenerateNormals flag (which explains the indirection below).
  vtkFunctionData* fields = this->GetInput()->GetFunctionData();
  vtkDataSetAttributes* outarrays = this->OutputMesh->GetPointData();
  outarrays->Initialize(); // empty, turn off all attributes, and set CopyAllOn to true.

  //vtkFieldData::BasicIterator fit = outarrays->ComputeRequiredArrays( fields );
  this->OutputAttributes = new vtkDataArray* [ fields->GetNumberOfArrays() + 1 ]; // + 1: GenerateNormals is on and fields doesn't have Normals
  this->OutputAttributeIndices = new int [ fields->GetNumberOfArrays() + 1 ];

  // OK, we always add normals as the 0-th array so that there's less work to do inside the tight loop (OutputTriangle)
  int attrib = 0;
  if ( this->GenerateNormals )
    {
    this->OutputAttributes[ attrib ] = vtkFloatArray::New();
    this->OutputAttributes[ attrib ]->SetNumberOfComponents(3);
    this->OutputAttributes[ attrib ]->SetName( "Normals" );
    this->OutputAttributeIndices[ attrib ] = outarrays->AddArray( this->OutputAttributes[attrib] );
    this->OutputMesh->GetPointData()->SetActiveAttribute( this->OutputAttributeIndices[attrib], vtkDataSetAttributes::NORMALS );
    ++attrib;
    maxNumComponents = 3;
    }

  for ( int a = 0; a < fields->GetNumberOfArrays(); ++a )
    {
    if ( fields->IsArrayAnAttribute( a ) == vtkDataSetAttributes::NORMALS && this->GenerateNormals )
      continue;

    vtkDataArray* array = fields->GetArray( a );
    this->OutputAttributes[ attrib ] = vtkDataArray::CreateDataArray( array->GetDataType() );
    this->OutputAttributes[ attrib ]->SetNumberOfComponents( array->GetNumberOfComponents() );
    this->OutputAttributes[ attrib ]->SetName( array->GetName() );
    this->OutputAttributeIndices[ attrib ] = outarrays->AddArray( this->OutputAttributes[ attrib ] );
    this->OutputAttributes[ attrib ]->Delete(); // output mesh now owns the array
    int attribType;
    if ( (attribType = fields->IsArrayAnAttribute( a )) != -1 )
      outarrays->SetActiveAttribute( this->OutputAttributeIndices[ attrib ], attribType );

		this->Subdivider->PassField( a, array->GetNumberOfComponents(), this->Tessellator );
    ++attrib;
    }
}

void vtkShoeMeshPolygonizer::Teardown()
{
  this->OutputMesh = 0;
  this->OutputPoints = 0;
  if ( this->OutputAttributes )
    delete [] this->OutputAttributes;
  if ( this->OutputAttributeIndices )
    delete [] this->OutputAttributeIndices;
	this->Subdivider->ResetFieldList();
}

void vtkShoeMeshPolygonizer::Execute()
{
  vtkShoeMesh* mesh = this->GetInput();

  this->SetupOutput();

  this->Tessellator->SetEdgeCallback( AddALine );
  this->Tessellator->SetTriangleCallback( AddATriangle );
  this->Tessellator->SetTetrahedronCallback( AddATetrahedron );
  this->Tessellator->SetPrivateData( this );

  vtkShoeMeshIterator theEnd( mesh->End() );
  this->CurrentCell = vtkShoeMeshIterator( mesh->Begin(), MeshOrder );
  while ( this->CurrentCell != theEnd )
    {
    if ( this->CurrentCell.GetCellOps()->EmbeddingDimension <= 3 )
      {
      this->Subdivider->SetCell( this->CurrentCell );
      this->CurrentCell.GetCellOps()->Tessellate( this->Tessellator );
      }
    ++this->CurrentCell;
    }

  this->Teardown();
}

