/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkObjectFactory.h>

#include <vtkShoeMesh.h>
#include <vtkCellOps.h>
#include <vtkShoeMeshChordLengthSubdivision.h>


vtkCxxRevisionMacro(vtkShoeMeshChordLengthSubdivision,"$Revision: 8692 $");
vtkStandardNewMacro(vtkShoeMeshChordLengthSubdivision);

void vtkShoeMeshChordLengthSubdivision::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
}

vtkShoeMeshChordLengthSubdivision::vtkShoeMeshChordLengthSubdivision()
{
  this->AllowableChordError2 = .01; // in world coordinates
}

vtkShoeMeshChordLengthSubdivision::~vtkShoeMeshChordLengthSubdivision()
{
}

void vtkShoeMeshChordLengthSubdivision::SetAllowableChordError2( double e )
{
  if ( e == this->AllowableChordError2 )
    return;
  if ( e > 0. )
    {
    this->AllowableChordError2 = e;
    }
  else
    {
    this->AllowableChordError2 = 1.;
    }
  this->Modified();
}

bool vtkShoeMeshChordLengthSubdivision::EvaluateEdge(
  const double* vtkNotUsed(p0), double* p1, const double* vtkNotUsed(p2), int field_start )
{
  double real_p1[3];
  double dist=0.;
  double tmp;

  this->CurrentCell.GetCellOps()->EvaluateGeometry( real_p1, this->CurrentCell, p1+3 );
  for ( int i=0; i<3; ++i )
    {
    tmp = real_p1[i] - p1[i];
    dist += tmp*tmp;
    }

  bool rval = dist > this->AllowableChordError2;
  if ( rval )
    {
    for ( int j=0; j<3; ++j )
      p1[j] = real_p1[j];

		// Evaluate field values at the midpoint
		this->EvaluateFields( p1, field_start );
    }

  return rval;
}

