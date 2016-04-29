/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkShoeMesh.h>
#include <vtkShoeMeshToUnstructuredGridFilter.h>



vtkCxxRevisionMacro(vtkShoeMeshToUnstructuredGridFilter, "$Revision: 759 $");

void vtkShoeMeshToUnstructuredGridFilter::PrintSelf(vtkstd::ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf( os, indent );
}

void vtkShoeMeshToUnstructuredGridFilter::SetInput( vtkShoeMesh* input )
{
	this->vtkProcessObject::SetNthInput( 0, input );
}

vtkShoeMesh* vtkShoeMeshToUnstructuredGridFilter::GetInput()
{
	if ( this->NumberOfInputs < 1 )
		return 0;

	return (vtkShoeMesh*)(this->Inputs[0]);
}

const vtkShoeMesh* vtkShoeMeshToUnstructuredGridFilter::GetInput() const
{
	if ( this->NumberOfInputs < 1 )
		return 0;

	return (vtkShoeMesh*)(this->Inputs[0]);
}

void vtkShoeMeshToUnstructuredGridFilter::ComputeInputUpdateExtents( vtkDataObject* output )
{
	vtkDataObject* input = this->GetInput();

	if ( input ) {
		this->vtkUnstructuredGridSource::ComputeInputUpdateExtents( output );
		input->RequestExactExtentOn();
	}
}

