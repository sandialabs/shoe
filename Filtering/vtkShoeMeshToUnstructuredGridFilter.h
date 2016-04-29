/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-09-30 19:04:45 -0700 (Tue, 30 Sep 2003) $
  Version:   $Revision: 759 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkShoeMeshToUnstructuredGridFilter_h
#define __vtkShoeMeshToUnstructuredGridFilter_h

#include <vtksnlConfigure.h>
#include <vtksnlFilteringWin32Header.h>
#include <vtkUnstructuredGridSource.h>

// .NAME vtkShoeMeshToUnstructuredGridFilter - abstract class that filters nonlinear FEM elements into polydata
// .SECTION Description
// This is an abstract filter class for converting nonlinear FEM elements
// into polygonal data for use with the gahoodius number of VTK algorithms
// for piecewise linear manifolds. Yes, I did say gahoodius.
//
// .SECTION See Also
// vtkShoeMeshPolygonizer vtkShoeMesh vtkUnstructuredGrid

class vtkShoeMesh;

class VTK_SNL_FILTERING_EXPORT vtkShoeMeshToUnstructuredGridFilter
	: public vtkUnstructuredGridSource
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshToUnstructuredGridFilter,vtkUnstructuredGridSource);
		void PrintSelf( ostream& os, vtkIndent indent );

		// Description:
		// Set and get the input data in a away that restricts it to
		// be a vtkShoeMesh.
		virtual void SetInput( vtkShoeMesh* input );
		vtkShoeMesh* GetInput();
		const vtkShoeMesh* GetInput() const;

		// Description:
		// Implementing this restricts downstream requests from
		// getting more data than they asked for.
		virtual void ComputeInputUpdateExtents( vtkDataObject* output );

	protected:
		vtkShoeMeshToUnstructuredGridFilter() { this->NumberOfRequiredInputs = 1; }
		~vtkShoeMeshToUnstructuredGridFilter() {}

	private:
		vtkShoeMeshToUnstructuredGridFilter( const vtkShoeMeshToUnstructuredGridFilter& ); // Not implemented.
		void operator = ( const vtkShoeMeshToUnstructuredGridFilter& ); // Not implemented.
};

#endif // __vtkShoeMeshToUnstructuredGridFilter_h

