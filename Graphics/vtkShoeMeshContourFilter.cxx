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
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMath.h>
#include <vtkTetra.h>

#include <vtkShoeMeshContourFilter.h>

#include <vtkShoeMesh.h>
#include <vtkShoeMeshIterator.h>
#include <vtkCellOps.h>
#include <vtkAdaptiveTessellator.h>
#include <vtkIsosurfaceSubdivisionAlgorithm.h>

vtkCxxRevisionMacro(vtkShoeMeshContourFilter, "$Revision: 1205 $");
vtkStandardNewMacro(vtkShoeMeshContourFilter);

vtkShoeMeshContourFilter::vtkShoeMeshContourFilter()
	: FieldId( 0 ), FieldComponent( 0 ), Isovalue( 0. )
{
	// Override the default subdivision algorithm with our own
	this->SetSubdivider( vtkIsosurfaceSubdivisionAlgorithm::New() );
}

vtkShoeMeshContourFilter::~vtkShoeMeshContourFilter()
{
	this->SetSubdivider( 0 );
	this->SetTessellator( 0 );
}

void vtkShoeMeshContourFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf( os, indent );
	os << indent << "FieldId:        " << this->FieldId << endl
	   << indent << "FieldComponent: " << this->FieldComponent << endl
		 << indent << "Isovalue:       " << this->Isovalue << endl;
}

void vtkShoeMeshContourFilter::SetSubdivider( vtkSubdivisionAlgorithm* s )
{
	if ( this->Subdivider == s )
		return;

	vtkIsosurfaceSubdivisionAlgorithm *v = dynamic_cast<vtkIsosurfaceSubdivisionAlgorithm*>(s);
	if ( v )
		{
		v->SetFieldId( this->FieldId );
		v->SetFieldComponent( this->FieldComponent );
		v->SetIsovalue( this->Isovalue );
		}
	else
		{
		vtkWarningMacro( "Subdivision algorithm must be derived from vtkIsosurfaceSubdivisionAlgorithm" );
		return;
		}

	if ( this->Subdivider )
		this->Subdivider->UnRegister( this );

	this->Subdivider = v;
		if ( this->Subdivider )
		this->Subdivider->Register( this );

	if ( this->Tessellator )
		this->Tessellator->SetSubdivisionAlgorithm( this->Subdivider );

	this->Modified();
}

void vtkShoeMeshContourFilter::SetFieldId( int id )
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
}

void vtkShoeMeshContourFilter::SetFieldComponent( int c )
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

void vtkShoeMeshContourFilter::SetIsovalue( double value )
{
	if ( value == this->Isovalue )
		return;

	this->Isovalue=value;
	this->Modified();
}

void vtkShoeMeshContourFilter::Execute()
{
	// Eventally, this should look for/create seed cells and iterate over just
	// the seed cells. Right now, it just iterates over all the cells.
	
	// This is horribly wrong. We need to tessellate each cell,
	// subdivide where critical points exist, and then use the
	// adaptive tessellator on the output isosurface.
	this->vtkShoeMeshPolygonizer::Execute();
}

