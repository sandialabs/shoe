/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-12-31 21:15:12 -0800 (Wed, 31 Dec 2003) $
  Version:   $Revision: 1201 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkShoeMeshNoSubdivision_h
#define __vtkShoeMeshNoSubdivision_h
// .NAME vtkShoeMeshNoSubdivision - never subdivides an element during tessellation
//
// .SECTION Description
// A vtkAdaptiveTessellator subdivision algorithm that won't subdivide.
//
// .SECTION See Also
// vtkShoeMeshSubdivisionAlgorithm

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>

#include <vtkShoeMeshIterator.h>
#include <vtkShoeMeshSubdivisionAlgorithm.h>
#include <vtkShoeMeshIterator.h>

class vtkShoeMesh;

class VTK_SNL_COMMON_EXPORT vtkShoeMeshNoSubdivision
	: public vtkShoeMeshSubdivisionAlgorithm
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshNoSubdivision,vtkShoeMeshSubdivisionAlgorithm);
		static vtkShoeMeshNoSubdivision* New();
		virtual void PrintSelf( ostream& os, vtkIndent indent );

		virtual bool EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start );

	protected:
		vtkShoeMeshNoSubdivision();
		virtual ~vtkShoeMeshNoSubdivision();

	private:
		vtkShoeMeshNoSubdivision( const vtkShoeMeshNoSubdivision& ); // Not implemented.
		void operator = ( const vtkShoeMeshNoSubdivision& ); // Not implemented.

};

#endif // __vtkShoeMeshNoSubdivision_h

