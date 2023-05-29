/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2004-09-21 11:38:10 -0700 (Tue, 21 Sep 2004) $
  Version:   $Revision: 2526 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkFlowerSubdivision_h
#define __vtkFlowerSubdivision_h
// .NAME vtkFlowerSubdivision - a subclass of
// vtkDataSetSubdivisionAlgorithm for testing the tessellator.
//
// .SECTION Description
// This is a subclass of vtkDataSetSubdivisionAlgorithm that is used for
// tessellating cells of a vtkDataSet using the implicit equation of a
// 3-D polar flower surface as a size map specification. 
//
// .SECTION See Also
// vtkDataSetSubdivisionAlgorithm

#include "vtkDataSetSubdivisionAlgorithm.h"

class vtkImplicitFunction;

class vtkFlowerSubdivision : public vtkDataSetSubdivisionAlgorithm
{
  public:
    vtkTypeRevisionMacro(vtkFlowerSubdivision,vtkDataSetSubdivisionAlgorithm);
		static vtkFlowerSubdivision* New();
    virtual void PrintSelf( ostream& os, vtkIndent indent );

		virtual bool EvaluateEdge( const double* p0, double* midpt, const double* p1, int field_start );

    virtual void SetRefiner( vtkImplicitFunction* );
    vtkGetObjectMacro(Refiner,vtkImplicitFunction);

  protected:
    vtkFlowerSubdivision();
    virtual ~vtkFlowerSubdivision();

		vtkImplicitFunction* Refiner;

  private:
    vtkFlowerSubdivision( const vtkFlowerSubdivision& ); // Not implemented.
    void operator = ( const vtkFlowerSubdivision& ); // Not implemented.
};

#endif // __vtkFlowerSubdivision_h
