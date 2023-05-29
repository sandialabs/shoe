/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include "vtkObjectFactory.h"
#include "vtkCellType.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkFieldData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkTetra.h"
#include "vtkGenericAttributeCollection.h"

#include "vtkShoeBoxContourFilter.h"

#include "vtkAdaptiveTessellator.h"
#include "vtkShoeAttribute.h"
#include "vtkShoeBox.h"
#include "vtkShoeBoxIsosurfaceSubdivision.h"
#include "vtkShoeBoxPartition.h"
#include "vtkShoeBoxPartitionIterator.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellIterator.h"

// Define VTK_DBG_SCF to print debug messages.
#undef VTK_DBG_SCF

vtkCxxRevisionMacro(vtkShoeBoxContourFilter, "$Revision: 6043 $");
vtkStandardNewMacro(vtkShoeBoxContourFilter);

vtkShoeBoxContourFilter::vtkShoeBoxContourFilter()
	: FieldId( 0 ), FieldComponent( 0 ), MaximumNumberOfSubdivisions( 0 ),
    Isovalue( 0. ), Tolerance( 1.e-4 ), Cell( 0 ), Attribute( 0 ),
    Partition( 0 ), Tessellator( 0 ), Subdivider( 0 ), ContourOut( 0 )
{
	// Override the default subdivision algorithm with our own
	this->SetSubdivider( vtkShoeBoxIsosurfaceSubdivision::New() );
  this->SetTessellator( vtkAdaptiveTessellator::New() );
}

vtkShoeBoxContourFilter::~vtkShoeBoxContourFilter()
{
	this->SetSubdivider( 0 );
	this->SetTessellator( 0 );
	this->SetPartition( 0 );
}

void vtkShoeBoxContourFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf( os, indent );
	os << indent << "FieldId:        " << this->FieldId << vtkstd::endl
	   << indent << "FieldComponent: " << this->FieldComponent << vtkstd::endl
     << indent << "MaxSubdivisions:" << this->MaximumNumberOfSubdivisions << vtkstd::endl
		 << indent << "Isovalue:       " << this->Isovalue << vtkstd::endl
		 << indent << "Tolerance:      " << this->Tolerance << vtkstd::endl
		 << indent << "Cell:           " << this->Cell << vtkstd::endl
		 << indent << "Attribute:      " << this->Attribute << vtkstd::endl
		 << indent << "Subdivider:     " << this->Subdivider << vtkstd::endl
		 << indent << "Tessellator:    " << this->Tessellator << vtkstd::endl
		 << indent << "Partition:      " << this->Partition << vtkstd::endl;
}

void vtkShoeBoxContourFilter::SetPartition( vtkShoeBoxPartition* partition )
{
  vtkSetObjectBodyMacro(Partition,vtkShoeBoxPartition,partition);
}

void vtkShoeBoxContourFilter::SetTessellator( vtkAdaptiveTessellator* tess )
{
  vtkSetObjectBodyMacro(Tessellator,vtkAdaptiveTessellator,tess);
	if ( this->Tessellator )
		this->Tessellator->SetSubdivisionAlgorithm( this->Subdivider );
}

void vtkShoeBoxContourFilter::SetSubdivider( vtkShoeBoxIsosurfaceSubdivision* v )
{
	if ( this->Subdivider == v )
		return;

	if ( this->Subdivider )
		this->Subdivider->UnRegister( this );

	this->Subdivider = v;
		if ( this->Subdivider )
		this->Subdivider->Register( this );

	if ( this->Tessellator )
		this->Tessellator->SetSubdivisionAlgorithm( this->Subdivider );

	this->Modified();
}

void vtkShoeBoxContourFilter::SetFieldId( int id )
{
	if ( id < 0 )
		{
		vtkWarningMacro( "Isocontouring requires a valid field id (>=0)" );
		return;
		}

	if ( id == this->FieldId )
		return;

	this->FieldId = id;
	this->Modified();
  this->SetPartition( 0 ); // we must construct a new partition
}

void vtkShoeBoxContourFilter::SetFieldComponent( int c )
{
	if ( c < 0 )
		{
		vtkWarningMacro( "Isocontouring requires a valid field component (>=0)" );
		return;
		}

	if ( c == this->FieldComponent )
		return;

	this->FieldComponent = c;
	this->Modified();
}

void vtkShoeBoxContourFilter::SetIsovalue( double value )
{
	if ( value == this->Isovalue )
		return;

	this->Isovalue=value;
	this->Modified();
}

int vtkShoeBoxContourFilter::RequestData(
  vtkInformation* vtkNotUsed(request),
  vtkInformationVector** inputVector,
  vtkInformationVector* outputVector )
{
  // get the info objects
  vtkInformation* inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  // get the input and ouptut
  vtkShoeBox* mesh = vtkShoeBox::SafeDownCast( inInfo->Get( vtkDataObject::DATA_OBJECT() ) );
  this->ContourOut = vtkUnstructuredGrid::SafeDownCast( outInfo->Get( vtkDataObject::DATA_OBJECT() ) );
  this->ContourOut->Reset();
  this->ContourOut->Allocate();
  this->ContourOut->SetPoints( vtkPoints::New() );

  this->Attribute = vtkShoeAttribute::SafeDownCast( mesh->GetAttributes()->GetAttribute( this->FieldId ) );
  if ( ! this->Attribute )
    {
    vtkErrorMacro("No valid attribute (" << this->FieldId << ") to isosurface.");
    }
  int nc = this->Attribute->GetNumberOfComponents();

  // Make sure the mesh is properly partitioned
  if ( ! this->Partition )
    {
    vtkGenericAttributeCollection* ac = vtkGenericAttributeCollection::New();
    ac->InsertNextAttribute( this->Attribute );
    //ac->InsertNextAttribute( mesh->GetGeometry() );
    vtkShoeBoxPartition* p = mesh->PartitionMesh( ac );
    p->Dump();
    this->SetPartition( p );
    p->Delete();
    ac->Delete();
    }

  // Set up the subdivider and tessellator
  this->Subdivider->SetAttribute( this->Attribute );
  this->Subdivider->SetComponent( 0 );
  this->Subdivider->SetTolerance( this->Tolerance );
  this->Subdivider->SetIsovalue( this->Isovalue );
  this->Tessellator->SetTetrahedronCallback( vtkShoeBoxContourFilter::AddTet );
  this->Tessellator->SetTriangleCallback( vtkShoeBoxContourFilter::AddTri );
  this->Tessellator->SetEdgeCallback( vtkShoeBoxContourFilter::AddLine );
  this->Tessellator->SetPrivateData( (void*) this );
  this->Tessellator->SetMaximumNumberOfSubdivisions( this->MaximumNumberOfSubdivisions );
  this->Tessellator->SetEmbeddingDimension( 2, 3 );
  // Should nc be 1? (so that we only pass the isocontoured component of this->Attribute)
  this->Subdivider->PassField( this->FieldId, nc, this->Tessellator );

  // Storage for geometry/parameter/attribute values at corners and edge midpoints
  vtkstd::vector<double> corners;
  vtkstd::vector<double> midpnts;
  int ptsz = nc + 6;
  corners.reserve( 4*ptsz );
  midpnts.reserve( 6*ptsz );
  double* v0 = &corners[0];
  double* v1 = &corners[ptsz];
  double* v2 = &corners[2*ptsz];
  double* v3 = &corners[3*ptsz];
  double* m0 = &midpnts[0];
  double* m1 = &midpnts[ptsz];
  double* m2 = &midpnts[2*ptsz];
  double* m3 = &midpnts[3*ptsz];
  double* m4 = &midpnts[4*ptsz];
  double* m5 = &midpnts[5*ptsz];
  vtkShoeCellIterator* scit = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator( 3 ) );
  vtkShoeBoxPartitionIterator* sbpit = vtkShoeBoxPartitionIterator::New();
  for ( scit->Begin(); ! scit->IsAtEnd(); scit->Next() )
    {
    this->Cell = vtkShoeCell::SafeDownCast( scit->GetCell() );
#ifdef VTK_DBG_SCF
    vtkstd::cout << "Cell " << this->Cell->GetId() << vtkstd::endl;
#endif // VTK_DBG_SCF
    this->Subdivider->SetCell( this->Cell );
    vtkIdType* conn = this->Cell->GetConnectivity();
    vtkIdType vdof = conn[ this->Cell->GetNumberOfPoints() + this->Cell->GetNumberOfDOFNodes() - 1 ];
    for ( sbpit->Begin( this->Partition, vdof ); ! sbpit->IsAtEnd(); sbpit->Next() )
      {
      // Evaluate attribute at each corner node.
      sbpit->GetTetrahedron( v0 + 3, v1 + 3, v2 + 3, v3 + 3 );
      this->Cell->InterpolateTuple( this->Attribute, v0 + 3, v0 + 6 );
      this->Cell->InterpolateTuple( this->Attribute, v1 + 3, v1 + 6 );
      this->Cell->InterpolateTuple( this->Attribute, v2 + 3, v2 + 6 );
      this->Cell->InterpolateTuple( this->Attribute, v3 + 3, v3 + 6 );
#ifdef VTK_DBG_SCF
      vtkstd::cout
        << "Tet " << sbpit->GetTetrahedronId() << " (" << (sbpit->GetCurrentEdgeFacet() / 8) << ", " << (sbpit->GetCurrentEdgeFacet() % 8)
        << "): " << v0[6] << ", " << v1[6] << ", " << v2[6] << ", " << v3[6]
        << vtkstd::endl;
#endif // VTK_DBG_SCF

      // Determine the case number from the values.
      // Each bit is an endpoint color (0 when the endpoint is below the isovalue, 1 when above it)
      int bitcode = 0;
      bitcode  = v0[6+this->FieldComponent] < this->Isovalue ? 0 : 1;
      bitcode |= v1[6+this->FieldComponent] < this->Isovalue ? 0 : 2;
      bitcode |= v2[6+this->FieldComponent] < this->Isovalue ? 0 : 4;
      bitcode |= v3[6+this->FieldComponent] < this->Isovalue ? 0 : 8;
      if ( !bitcode || bitcode == 15 )
        continue;

#ifdef VTK_DBG_SCF
      vtkstd::cout << "    bitcode " << bitcode << vtkstd::endl;
#endif // VTK_DBG_SCF
      // Find root along each relevant edge.
      // Exclusive-or guarantees that if the corresponding bits are both set or both unset in
      // bitcode, then both bits of result are on. Otherwise, only one will be on, indicating
      // the edge has endpoints that are colored differently.
      if ( ( ( bitcode ^ 3 ) & 3 ) != 0 )
        {
        this->EdgeContourIntersect( v0, v1, m0 ); // Edge 0-1
        }
      if ( ( ( bitcode ^ 6 ) & 6 ) != 0 )
        {
        this->EdgeContourIntersect( v1, v2, m1 ); // Edge 1-2
        }
      if ( ( ( bitcode ^ 5 ) & 5 ) != 0 )
        {
        this->EdgeContourIntersect( v2, v0, m2 ); // Edge 2-0
        }
      if ( ( ( bitcode ^ 9 ) & 9 ) != 0 )
        {
        this->EdgeContourIntersect( v0, v3, m3 ); // Edge 0-3
        }
      if ( ( ( bitcode ^ 10 ) & 10 ) != 0 )
        {
        this->EdgeContourIntersect( v1, v3, m4 ); // Edge 1-3
        }
      if ( ( ( bitcode ^ 12 ) & 12 ) != 0 )
        {
        this->EdgeContourIntersect( v2, v3, m5 ); // Edge 2-3
        }

      // Look up triangles.
      // Pass triangles to adaptive tessellator.
      switch ( bitcode )
        {
      case 1:
        this->Tessellator->AdaptivelySample2Facet( m0, m3, m2 );
        break;
      case 2:
        this->Tessellator->AdaptivelySample2Facet( m0, m1, m4 );
        break;
      case 3:
        this->Tessellator->AdaptivelySample2Facet( m1, m4, m2 );
        this->Tessellator->AdaptivelySample2Facet( m4, m3, m2 );
        break;
      case 4:
        this->Tessellator->AdaptivelySample2Facet( m1, m2, m5 );
        break;
      case 5:
        this->Tessellator->AdaptivelySample2Facet( m0, m3, m1 );
        this->Tessellator->AdaptivelySample2Facet( m3, m5, m1 );
        break;
      case 6:
        this->Tessellator->AdaptivelySample2Facet( m0, m2, m4 );
        this->Tessellator->AdaptivelySample2Facet( m2, m5, m4 );
        break;
      case 7:
        this->Tessellator->AdaptivelySample2Facet( m3, m5, m4 );
        break;
      case 8:
        this->Tessellator->AdaptivelySample2Facet( m3, m4, m5 );
        break;
      case 9:
        this->Tessellator->AdaptivelySample2Facet( m0, m4, m2 );
        this->Tessellator->AdaptivelySample2Facet( m2, m4, m5 );
        break;
      case 10:
        this->Tessellator->AdaptivelySample2Facet( m0, m1, m3 );
        this->Tessellator->AdaptivelySample2Facet( m3, m1, m5 );
        break;
      case 11:
        this->Tessellator->AdaptivelySample2Facet( m1, m5, m2 );
        break;
      case 12:
        this->Tessellator->AdaptivelySample2Facet( m1, m2, m4 );
        this->Tessellator->AdaptivelySample2Facet( m4, m2, m3 );
        break;
      case 13:
        this->Tessellator->AdaptivelySample2Facet( m0, m4, m1 );
        break;
      case 14:
        this->Tessellator->AdaptivelySample2Facet( m0, m2, m3 );
        break;
        }
      }
    }

  return 1;
}

int vtkShoeBoxContourFilter::FillInputPortInformation( int port, vtkInformation* info )
{
  info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkShoeBox" );
  return 1;
}

void vtkShoeBoxContourFilter::AddTet( const double* x, const double* y, const double* z, const double* w,
  vtkSubdivisionAlgorithm* sa, void* v, const void* cv )
{
  vtkShoeBoxContourFilter* cf = (vtkShoeBoxContourFilter*)v;
  cf->AddTet( x, y, z, w );
}

void vtkShoeBoxContourFilter::AddTri( const double* x, const double* y, const double* z, vtkSubdivisionAlgorithm* sa, void* v, const void* cv )
{
  vtkShoeBoxContourFilter* cf = (vtkShoeBoxContourFilter*)v;
#if 0
  fprintf( stdout, "%g %g %g\n%g %g %g\n%g %g %g\n\n",
    x[0], x[1], x[2], y[0], y[1], y[2], z[0], z[1], z[2] );
#endif // 0
  cf->AddTri( x, y, z );
}

void vtkShoeBoxContourFilter::AddLine( const double* x, const double* y, vtkSubdivisionAlgorithm* sa, void* v, const void* cv )
{
  vtkShoeBoxContourFilter* cf = (vtkShoeBoxContourFilter*)v;
  cf->AddLine( x, y );
}

void vtkShoeBoxContourFilter::AddTet( const double* x, const double* y, const double* z, const double* w )
{
  vtkIdType pts[4];
  pts[0] = this->ContourOut->GetPoints()->InsertNextPoint( x );
  pts[1] = this->ContourOut->GetPoints()->InsertNextPoint( y );
  pts[2] = this->ContourOut->GetPoints()->InsertNextPoint( z );
  pts[3] = this->ContourOut->GetPoints()->InsertNextPoint( w );
  this->ContourOut->InsertNextCell( VTK_TETRA, 4, pts );
}

void vtkShoeBoxContourFilter::AddTri( const double* x, const double* y, const double* z )
{
  // Don't insert degenerate triangles.
  bool degen[3] = { true, true, true };
  for ( int i = 0; i < 3; ++i )
    {
    if ( x[i] != y[i] ) { degen[0] = false; }
    if ( y[i] != z[i] ) { degen[1] = false; }
    if ( z[i] != x[i] ) { degen[2] = false; }
    }
  if ( degen[0] || degen[1] || degen[2] )
    {
    return;
    }
  vtkIdType pts[3];
  pts[0] = this->ContourOut->GetPoints()->InsertNextPoint( x );
  pts[1] = this->ContourOut->GetPoints()->InsertNextPoint( y );
  pts[2] = this->ContourOut->GetPoints()->InsertNextPoint( z );
  vtkIdType id = this->ContourOut->InsertNextCell( VTK_TRIANGLE, 3, pts );
  (void)id;
  //fprintf( stdout, "%lld added with %lld %lld %lld\n", id, pts[0], pts[1], pts[2] );
}

void vtkShoeBoxContourFilter::AddLine( const double* x, const double* y )
{
  vtkIdType pts[2];
  pts[0] = this->ContourOut->GetPoints()->InsertNextPoint( x );
  pts[1] = this->ContourOut->GetPoints()->InsertNextPoint( y );
  this->ContourOut->InsertNextCell( VTK_LINE, 2, pts );
}

int vtkShoeBoxContourFilter::EdgeContourIntersect( const double* p0, const double* p1, double* pint )
{
  double tleft = 0.;
  double tright = 1.;
  double fleft = p0[6 + this->FieldComponent] - this->Isovalue;
  double t = 0.5;
  double f;
  int count = 0;

  do {
    double t1 = 1. - t;
    pint[3] = t1 * p0[3] + t * p1[3];
    pint[4] = t1 * p0[4] + t * p1[4];
    pint[5] = t1 * p0[5] + t * p1[5];
    this->Cell->InterpolateTuple( this->Attribute, pint + 3, pint + 6 );
    f = pint[6 + this->FieldComponent] - this->Isovalue;
    if ( fleft*f > 0. )
      {
      tleft = t;
      fleft = f;
      }
    else
      {
      tright = t;
      }
    t = 0.5*(tleft + tright);
  } while ( count++ < 15 && fabs( f ) > this->Tolerance );
  if ( count >= 15 )
    {// If we failed to find a root, take the closest value
    if ( fabs( fleft ) < fabs( f ) )
      {
      pint[3] = (1. - tleft) * p0[3] + tleft * p1[3];
      pint[4] = (1. - tleft) * p0[4] + tleft * p1[4];
      pint[5] = (1. - tleft) * p0[5] + tleft * p1[5];
      this->Cell->InterpolateTuple( this->Attribute, pint + 3, pint + 6 );
      this->Cell->EvaluateLocation( -1, pint + 3, pint );
      return 0;
      }
    }
  this->Cell->EvaluateLocation( -1, pint + 3, pint );
  return 1;
}


