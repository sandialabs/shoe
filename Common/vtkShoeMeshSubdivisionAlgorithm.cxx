/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkShoeMeshSubdivisionAlgorithm.h>
#include <vtkCellOps.h>


vtkCxxRevisionMacro(vtkShoeMeshSubdivisionAlgorithm,"$Revision: 1198 $");

void vtkShoeMeshSubdivisionAlgorithm::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
}

void vtkShoeMeshSubdivisionAlgorithm::SetCell( const vtkShoeMeshIterator& c )
{
  this->CurrentCell = c;
  this->Modified();
}

double* vtkShoeMeshSubdivisionAlgorithm::EvaluateFields( double* vertex, int field_start )
{
	const vtkCellOps* ops = this->CurrentCell.GetCellOps();
	for ( int f=0; f<this->GetNumberOfFields(); ++f )
		ops->EvaluateField( vertex + field_start + this->GetFieldOffsets()[f], this->CurrentCell, vertex + 3, this->GetFieldIds()[f] );
	return vertex;
}
