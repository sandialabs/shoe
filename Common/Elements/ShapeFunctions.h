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
#ifndef Elem_ShapeFunctions_h
#define Elem_ShapeFunctions_h

#include <math.h>

#include <vtksnlConfigure.h>

class vtkShoeElemShapeFunctions
{
	public:
		/// Function pointers to various polynomial shape functions.
		static double (*EvaluatePhi[7])(double *x);
		static double (*EvaluatePsi[7])(double *x);
		static double (*EvaluateLegendre[7])(double *x);
		static double (*EvaluatePhiDerivative[7])(double *x);
		static double (*EvaluatePsiDerivative[7])(double *x);
		static double (*EvaluateLegendreDerivative[7])(double *x);

		static double phiValue[6][6];
		static double dphiValue[6][6];
		static double dpsiValue[6][6];
};

#endif // Elem_ShapeFunctions_h
