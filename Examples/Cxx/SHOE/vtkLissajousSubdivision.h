/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2004-07-24 21:56:19 -0700 (Sat, 24 Jul 2004) $
  Version:   $Revision: 2272 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkLissajousSubdivision_h
#define __vtkLissajousSubdivision_h
// .NAME vtkLissajousSubdivision - a subclass of
// vtkDataSetSubdivisionAlgorithm for testing the tessellator.
//
// .SECTION Description
// This is a subclass of vtkDataSetSubdivisionAlgorithm that is used for
// tessellating cells of a vtkDataSet using the implicit equation of a
// Lissajous curve as a size map specification. 
//
// .SECTION See Also
// vtkDataSetSubdivisionAlgorithm

#include "vtkDataSetSubdivisionAlgorithm.h"

class vtkImplicitFunction;

class vtkLissajousSubdivision : public vtkDataSetSubdivisionAlgorithm
{
  public:
    vtkTypeRevisionMacro(vtkLissajousSubdivision,vtkDataSetSubdivisionAlgorithm);
		static vtkLissajousSubdivision* New();
    virtual void PrintSelf( ostream& os, vtkIndent indent );

		virtual bool EvaluateEdge( const double* p0, double* midpt, const double* p1, int field_start );

    virtual void SetRefiner( vtkImplicitFunction* );
    vtkGetObjectMacro(Refiner,vtkImplicitFunction);

  protected:
    vtkLissajousSubdivision();
    virtual ~vtkLissajousSubdivision();

		vtkImplicitFunction* Refiner;

  private:
    vtkLissajousSubdivision( const vtkLissajousSubdivision& ); // Not implemented.
    void operator = ( const vtkLissajousSubdivision& ); // Not implemented.
};

#endif // __vtkLissajousSubdivision_h
