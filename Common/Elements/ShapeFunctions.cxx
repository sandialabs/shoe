/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <math.h>

#include <vtksnlConfigure.h>
#include <Elements/ShapeFunctions.h>



#include <Elements/ShapeFunctionsImpl.h>

// Function pointers to the integrated legendre polynomial evaluation functions.
double (*vtkShoeElemShapeFunctions::EvaluatePhi[7])(double *x) =
{
	0, 0, phi2, phi3, phi4, phi5, phi6
} ;

double (*vtkShoeElemShapeFunctions::EvaluatePsi[7])(double *x) =
{
	0, 0, psi2, psi3, psi4, psi5, psi6
};

double (*vtkShoeElemShapeFunctions::EvaluateLegendre[7])(double *x) =
{
	legendre0, legendre1, legendre2, legendre3, legendre4, legendre5, legendre6
};

double (*vtkShoeElemShapeFunctions::EvaluatePhiDerivative[7])(double *x) =
{
	0, 0, dphi2, dphi3, dphi4, dphi5, dphi6
};

double (*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[7])(double *x) =
{
	0, 0, dpsi2, dpsi3, dpsi4, dpsi5, dpsi6
};

double (*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[7])(double *x) =
{
	dleg0, dleg1, dleg2, dleg3, dleg4, dleg5, dleg6
};

// Static storage used by interpolation routines.
// NOT THREADSAFE AT ALL! NOT REENTRANT!
double vtkShoeElemShapeFunctions::phiValue[6][6];
double vtkShoeElemShapeFunctions::dphiValue[6][6];
double vtkShoeElemShapeFunctions::dpsiValue[6][6];
