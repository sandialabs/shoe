/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-09-29 10:24:11 -0700 (Mon, 29 Sep 2003) $
  Version:   $Revision: 730 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef Element_Point_h
#define Element_Point_h

#include <vtkCellOps.h>

class vtkShoeElemPoint
{
	public:
		static void EvaluateShapeFunctionPoint( double* shape, vtkShoeMeshIterator&, const int[3], const double[3] );
};

#endif // Element_Point_h
