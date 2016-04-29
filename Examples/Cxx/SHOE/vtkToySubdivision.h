/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2004-06-17 00:56:22 -0700 (Thu, 17 Jun 2004) $
  Version:   $Revision: 2100 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkToySubdivision_h
#define __vtkToySubdivision_h
// .NAME vtkToySubdivision - a subclass of vtkDataSetSubdivisionAlgorithm for testing the tessellator.
//
// .SECTION Description
// This is a subclass of vtkDataSetSubdivisionAlgorithm that is used for
// tessellating cells of a vtkDataSet using an implicit function to
// specify a region of space that should be refined.
//
// .SECTION See Also
// vtkDataSetSubdivisionAlgorithm

#include "vtkDataSetSubdivisionAlgorithm.h"

class vtkImplicitFunction;

class vtkToySubdivision : public vtkDataSetSubdivisionAlgorithm
{
  public:
    vtkTypeRevisionMacro(vtkToySubdivision,vtkDataSetSubdivisionAlgorithm);
		static vtkToySubdivision* New();
    virtual void PrintSelf( ostream& os, vtkIndent indent );

		virtual bool EvaluateEdge( const double* p0, double* midpt, const double* p1, int field_start );

    virtual void SetRefiner( vtkImplicitFunction* );
    vtkGetObjectMacro(Refiner,vtkImplicitFunction);

  protected:
    vtkToySubdivision();
    virtual ~vtkToySubdivision();

		vtkImplicitFunction* Refiner;

  private:
    vtkToySubdivision( const vtkToySubdivision& ); // Not implemented.
    void operator = ( const vtkToySubdivision& ); // Not implemented.
};

#endif // __vtkToySubdivision_h
