/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-12-10 17:27:57 -0800 (Wed, 10 Dec 2003) $
  Version:   $Revision: 1130 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkShoeMeshViewDependentSubdivision_h
#define __vtkShoeMeshViewDependentSubdivision_h
// .NAME vtkShoeMeshViewDependentSubdivision - subdivides faces for pixel-accurate geometry
//
// .SECTION Description
// For each edge passed to this algorithm, we perform the following tests to
// see if subdivision can be avoided:
// (1) If both vertices fall out of the viewing frustum, we do not subdivide
// (2) If the projected chord error is less than 1 pixel width, we do not subdivide
// (3) In the future, we may test whether the surface normals at the endpoints are similar
// Otherwise, we subdivide.
//
// .SECTION See Also
// vtkShoeMeshSubdivisionAlgorithm

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>

#include <vtkShoeMeshIterator.h>
#include <vtkShoeMeshSubdivisionAlgorithm.h>
#include <vtkShoeMeshIterator.h>

class vtkCamera;
class vtkMatrix4x4;
class vtkWindow;

class vtkShoeMesh;

class VTK_SNL_COMMON_EXPORT vtkShoeMeshViewDependentSubdivision
	: public vtkShoeMeshSubdivisionAlgorithm
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshViewDependentSubdivision,vtkShoeMeshSubdivisionAlgorithm);
		static vtkShoeMeshViewDependentSubdivision* New();
		virtual void PrintSelf( ostream& os, vtkIndent indent );

		virtual void SetCamera( vtkCamera* );
		virtual void SetWindow( vtkWindow* );
		virtual void SetAllowableChordError( double inPixels );
		virtual bool EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start );

	protected:
		vtkShoeMeshViewDependentSubdivision();
		virtual ~vtkShoeMeshViewDependentSubdivision();

		// Description:
		// An internal routine called to precompute some constants
		// used in EvaluateEdge() after the Window or Camera has been
		// changed.
		virtual void MemberSetup();

		vtkCamera* Camera;
		vtkWindow* Window;

		// Description:
		// The <b>square</b> of the allowable chord error. By
		// storing the square, we reduce the need for a sqrt
		// call in the EvaluateEdge() member.
		double AllowableChordError;
		vtkMatrix4x4* Transform;
		double PixelSize[2]; // Pixel size in normalized device coordinate space ([-1,1]x[-1,1]x[near,far])

	private:
		vtkShoeMeshViewDependentSubdivision( const vtkShoeMeshViewDependentSubdivision& ); // Not implemented.
		void operator = ( const vtkShoeMeshViewDependentSubdivision& ); // Not implemented.

};

#endif // __vtkShoeMeshViewDependentSubdivision_h

