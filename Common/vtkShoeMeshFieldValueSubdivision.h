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
#ifndef __vtkShoeMeshFieldValueSubdivision_h
#define __vtkShoeMeshFieldValueSubdivision_h
// .NAME vtkShoeMeshFieldValueSubdivision - subdivides edges based on L2 error at their midpoints
//
// .SECTION Description
// For each edge passed to this algorithm, we evaluate the position at the
// parametric center of the edge and compare it to a linear interpolation of
// the endpoint vertices. If the two are further apart than the specified
// chord length error, we signal that subdivision should occur. Otherwise,
// we evaluate the scalar field values at the midpoint and compare those
// values to linear interpolation. If any of those values exceed the specified
// L2 error, then subdivision should occur. Only if both geometric and
// field value tests fail will the edge remain undivided.
//
// .SECTION See Also
// vtkShoeMeshSubdivisionAlgorithm

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>

#include <vtkShoeMeshSubdivisionAlgorithm.h>

class VTK_SNL_COMMON_EXPORT vtkShoeMeshFieldValueSubdivision
	: public vtkShoeMeshSubdivisionAlgorithm
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshFieldValueSubdivision,vtkShoeMeshSubdivisionAlgorithm);
		static vtkShoeMeshFieldValueSubdivision* New();
		virtual void PrintSelf( ostream& os, vtkIndent indent );

		// Description:
		// Set a limit on the square of the L2 norm of the difference between the linearly
		// interpolated field value and the calculated (nonlinear) field value.
		// Whew, that's a mouthful. Think of it this way; for each field \a f (where \a f is
		// the <i>source</i> mesh field id), we take the difference between the linearly
		// interpolated midpoint value and one that we calculate using the nonlinear shape
		// functions. We square that difference (if the difference is a vector, we take the dot
		// product of the difference with itself) and then compare the result to
		// \a max. If the result is larger than \a max, we require the edge to be subdivided.
		// 
		// There is a \a max value associated with each field, \a f.
		// By default, \a max is infinitely large and is not a criterion for subdivision.
		// When you call \p SetAllowableL2Error2( \a f, \p Inf ),
		// field \a f will no longer be considered as a criterion for subdivision.
		// Calling \p SetAllowableL2Error2 with any other value will cause field \a f to
		// be checked.
		//
		// This check is performed <i>in addition</i> to the geometry check (which uses
		// the AllowableChordError2 parameter). You must \p SetAllowableChordError2(0)
		// to eliminate the geometry check.
		//
		// .SECTION Subtle Issues
		// Note that when a subdivision is required, <b>all</b> field values at the midpoint
		// vertex will be replaced with their exact (nonlinear) values.
		// When no subdivision is required, <b>none</b> of the field values at the midpoint
		// vertex will be affected (they will all be linearly interpolated).
		// This behavior may not meet the principle of minimum astonishment, but is
		// required for continuity.
		//
		// Another subtle issue is that this algorithm allows you to specify <b>no</b>
		// criteria for edge subdivision; you can specify that neither field values
		// nor geometry should be allowed to cause subdivision. If you do, this algorithm
		// acts like a very expensive no-op.
		virtual void SetAllowableL2Error2( int f, double max );
		double GetAllowableL2Error2( int field ) const;

		// Description:
		// Set the maximum allowable distance between the linear interpolation and exact
		// evaluation of an edge's midpoint before subdivision is required.
		//
		// Setting this to \p Inf will ignore geometry during subdivision (only field values
		// will be calculated).
		//
		// See \p SetAllowableL2Error2 for a more detailed discussion.
		virtual void SetAllowableChordError2( double max );
		double GetAllowableChordError2() const;

		virtual bool EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start );

		// Description:
		// These methods are inherited from vtkSubdivisionAlgorithm so that we can keep the
		// list of AllowableL2Error2 values up to date with the attributes being passed.
		//virtual void PassFields( vtkDataSetAttributes* outarrays, vtkAdaptiveTessellator* t );
		virtual int PassField( int sourceId, int sourceSize, vtkAdaptiveTessellator* t );
		virtual void ResetFieldList();
		virtual bool DontPassField( int sourceId, vtkAdaptiveTessellator* t );

	protected:
		vtkShoeMeshFieldValueSubdivision();
		virtual ~vtkShoeMeshFieldValueSubdivision();

		// Description:
		// Reallocate the AllowableL2Error2 array, preserving any values
		// not past this->NumberOfFields.
		void ResizeL2Array( int old_length );

		// Description:
		// A bit vector of the fields that serve as criteria
		// for edge subdivision. (When a bit is on, the AllowableL2Error
		// for the field corresponding to that bit must be checked
		// against the calculated and interpolated midpoint values.
		// When a bit is off, the field may be skipped.)
		int FieldCriteria;

		// Description:
		// The <b>square</b> of the allowable chord error. By
		// storing the square, we reduce the need for a sqrt
		// call in the EvaluateEdge() member.
		double* AllowableL2Error2;

		// Description:
		// The <b>square</b> of the allowable chord error. By
		// storing the square, we reduce the need for a sqrt
		// call in the EvaluateEdge() member.
		double AllowableChordError2;

	private:
		vtkShoeMeshFieldValueSubdivision( const vtkShoeMeshFieldValueSubdivision& ); // Not implemented.
		void operator = ( const vtkShoeMeshFieldValueSubdivision& ); // Not implemented.

};

//BTX

inline double vtkShoeMeshFieldValueSubdivision::GetAllowableChordError2() const { return this->AllowableChordError2; }

//ETX

#endif // __vtkShoeMeshFieldValueSubdivision_h

