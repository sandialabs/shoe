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
#include <vtkShoeMeshSource.h>



vtkCxxRevisionMacro(vtkShoeMeshSource,"$Revision: 759 $");

vtkShoeMeshSource::vtkShoeMeshSource()
{
	vtkShoeMesh* out = vtkShoeMesh::New();
	this->SetNthOutput( 0, out );
	out->ReleaseData();
	out->Delete();
}

vtkShoeMeshSource::~vtkShoeMeshSource()
{
}

void vtkShoeMeshSource::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf( os, indent );
}

vtkShoeMesh* vtkShoeMeshSource::GetOutput()
{
	return this->GetNumberOfOutputs() ? dynamic_cast<vtkShoeMesh*>( this->Outputs[0] ) : 0;
}

const vtkShoeMesh* vtkShoeMeshSource::GetOutput() const
{
	return this->NumberOfOutputs ? dynamic_cast<vtkShoeMesh*>( this->Outputs[0] ) : 0;
}

vtkShoeMesh* vtkShoeMeshSource::GetOutput( int i )
{
	return this->GetNumberOfOutputs() > i ? dynamic_cast<vtkShoeMesh*>( this->Outputs[i] ) : 0;
}

const vtkShoeMesh* vtkShoeMeshSource::GetOutput( int i ) const
{
	return this->NumberOfOutputs > i ? dynamic_cast<vtkShoeMesh*>( this->Outputs[i] ) : 0;
}

void vtkShoeMeshSource::SetOutput( vtkShoeMesh* output_in )
{
	this->vtkSource::SetNthOutput( 0, output_in );
}

void vtkShoeMeshSource::ComputeInputUpdateExtents(vtkDataObject *data)
{
	int piece=0, numPieces=1, ghostLevel=0;
	vtkShoeMesh* output = dynamic_cast<vtkShoeMesh*>(data);

	if ( ! output )
	{
		vtkErrorMacro( "InputUpdateExtents of an object that's null or not a vtkShoeMesh!" );
		return;
	}
	output->GetUpdateExtent( piece, numPieces, ghostLevel );

	// make sure piece is valid
	if (piece < 0 || piece >= numPieces)
		return;

	// just vtkstd::copy the Update extent as default behavior.
	for ( int i = 0; i < this->GetNumberOfInputs(); ++i )
		if ( this->Inputs[i] )
			this->Inputs[i]->SetUpdateExtent(piece, numPieces, ghostLevel);
}

