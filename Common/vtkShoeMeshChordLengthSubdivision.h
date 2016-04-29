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
#ifndef __vtkShoeMeshChordLengthSubdivision_h
#define __vtkShoeMeshChordLengthSubdivision_h
// .NAME vtkShoeMeshChordLengthSubdivision - subdivides edges based on chord-length error at their midpoints
//
// .SECTION Description
// For each edge passed to this algorithm, we evaluate the position at the
// parametric center of the edge and compare it to a linear interpolation of
// the endpoint vertices. If the two are further apart than the specified
// chord length error, we signal that subdivision should occur.
//
// .SECTION See Also
// vtkShoeMeshSubdivisionAlgorithm

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>

#include <vtkShoeMeshSubdivisionAlgorithm.h>

class VTK_SNL_COMMON_EXPORT vtkShoeMeshChordLengthSubdivision
	: public vtkShoeMeshSubdivisionAlgorithm
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshChordLengthSubdivision,vtkShoeMeshSubdivisionAlgorithm);
		static vtkShoeMeshChordLengthSubdivision* New();
		virtual void PrintSelf( ostream& os, vtkIndent indent );

		virtual void SetAllowableChordError2( double length_in );
		double GetAllowableChordError2() const;

		virtual bool EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start );

	protected:
		vtkShoeMeshChordLengthSubdivision();
		virtual ~vtkShoeMeshChordLengthSubdivision();

		// Description:
		// The <b>square</b> of the allowable chord error. By
		// storing the square, we reduce the need for a sqrt
		// call in the EvaluateEdge() member.
		double AllowableChordError2;

	private:
		vtkShoeMeshChordLengthSubdivision( const vtkShoeMeshChordLengthSubdivision& ); // Not implemented.
		void operator = ( const vtkShoeMeshChordLengthSubdivision& ); // Not implemented.

};

//BTX

inline double vtkShoeMeshChordLengthSubdivision::GetAllowableChordError2() const { return this->AllowableChordError2; }

//ETX

#endif // __vtkShoeMeshChordLengthSubdivision_h

