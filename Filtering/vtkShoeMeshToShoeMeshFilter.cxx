/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkShoeMesh.h>
#include <vtkShoeMeshToShoeMeshFilter.h>



vtkCxxRevisionMacro(vtkShoeMeshToShoeMeshFilter,"$Revision: 759 $");

vtkShoeMeshToShoeMeshFilter::vtkShoeMeshToShoeMeshFilter()
{
}

vtkShoeMeshToShoeMeshFilter::~vtkShoeMeshToShoeMeshFilter()
{
}

void vtkShoeMeshToShoeMeshFilter::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf( os, indent );
}

void vtkShoeMeshToShoeMeshFilter::SetInput( vtkShoeMesh* m )
{
	this->vtkProcessObject::SetNthInput( 0, m );
}

vtkShoeMesh* vtkShoeMeshToShoeMeshFilter::GetInput()
{
	return this->GetNumberOfInputs() ? dynamic_cast<vtkShoeMesh*>(this->Inputs[0]) : 0;
}

const vtkShoeMesh* vtkShoeMeshToShoeMeshFilter::GetInput() const
{
	return this->NumberOfInputs ? dynamic_cast<vtkShoeMesh*>(this->Inputs[0]) : 0;
}

