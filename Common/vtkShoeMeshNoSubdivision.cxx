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
#include <vtkShoeMeshNoSubdivision.h>


vtkCxxRevisionMacro(vtkShoeMeshNoSubdivision,"$Revision: 8692 $");
vtkStandardNewMacro(vtkShoeMeshNoSubdivision);

void vtkShoeMeshNoSubdivision::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
}

vtkShoeMeshNoSubdivision::vtkShoeMeshNoSubdivision()
{
}

vtkShoeMeshNoSubdivision::~vtkShoeMeshNoSubdivision()
{
}

bool vtkShoeMeshNoSubdivision::EvaluateEdge(
  const double* vtkNotUsed(p0), double* vtkNotUsed(p1), const double* vtkNotUsed(p2), int vtkNotUsed(field_start) )
{
	return false;
}

