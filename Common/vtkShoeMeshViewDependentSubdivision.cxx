/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkMatrix4x4.h>
#include <vtkObjectFactory.h>
#include <vtkWindow.h>
#include <vtkCamera.h>

#include <vtkShoeMesh.h>
#include <vtkCellOps.h>
#include <vtkShoeMeshViewDependentSubdivision.h>


vtkCxxRevisionMacro(vtkShoeMeshViewDependentSubdivision,"$Revision: 1169 $");
vtkStandardNewMacro(vtkShoeMeshViewDependentSubdivision);

void vtkShoeMeshViewDependentSubdivision::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
}

vtkShoeMeshViewDependentSubdivision::vtkShoeMeshViewDependentSubdivision()
{
  Camera = 0;
  Window = 0;
  Transform = 0;
  AllowableChordError = 1.; // in pixels
}

vtkShoeMeshViewDependentSubdivision::~vtkShoeMeshViewDependentSubdivision()
{
  if ( Camera )
    Camera->UnRegister( this );

  if ( Window )
    Window->UnRegister( this );
}

void vtkShoeMeshViewDependentSubdivision::SetCamera( vtkCamera* c )
{
  if ( c == this->Camera )
    return;

  if ( this->Camera )
    this->Camera->UnRegister( this );

  this->Camera = c;

  if ( this->Camera )
    this->Camera->Register( this );

  this->MemberSetup();
  this->Modified();
}

void vtkShoeMeshViewDependentSubdivision::SetWindow( vtkWindow* w )
{
  if ( w == this->Window )
    return;

  if ( this->Window )
    this->Window->UnRegister( this );

  this->Window = w;

  if ( this->Window )
    this->Window->Register( this );

  this->MemberSetup();
  this->Modified();
}

void vtkShoeMeshViewDependentSubdivision::SetAllowableChordError( double e )
{
  if ( e == this->AllowableChordError )
    return;

  if ( e > 0. )
    {
    this->AllowableChordError = e;
    }
  else
    {
    this->AllowableChordError = 1.;
    }

  this->Modified();
}

void vtkShoeMeshViewDependentSubdivision::MemberSetup()
{
  if ( ! this->Camera || ! this->Window )
    return;

  double nearz, farz;
  int* sz = this->Window->GetSize();
  int w=sz[0], h=sz[1];
  this->Camera->GetClippingRange( nearz, farz );
  double aspect = double(w)/double(h);
  this->Transform = this->Camera->GetCompositePerspectiveTransformMatrix( aspect, nearz, farz );
  this->Transform->Register(this);
  PixelSize[0] = 2./w;
  PixelSize[1] = 2./h;
}

bool vtkShoeMeshViewDependentSubdivision::EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start )
{
  double real_p1[4];
  double tmp_r = p1[3];

  this->CurrentCell.GetCellOps()->EvaluateGeometry( real_p1, this->CurrentCell, p1+3 );
  real_p1[3] = 1.0; // w=1
  p1[3] = 1.0; // w=1, smashing first parameter value (we'll fix it later...)
	bool subdivide = this->ViewDependentEval( p0, p1, real_p1, p2, field_start, this->Transform, this->PixelSize, this->AllowableChordError );
	if ( subdivide )
		{
    p1[3] = tmp_r; // don't forget we stomped the first parameter

    // Evaluate the field values at the new vertex-to-be.
    this->EvaluateFields( p1, field_start );
		return true;
		}
	else
		{
		// keep the linearly interpolated point in p1
		p1[3] = tmp_r; // don't forget we stomped the first parameter value
		return false;
		}
}

