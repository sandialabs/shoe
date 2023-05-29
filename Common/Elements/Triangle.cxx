/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtksnlConfigure.h>
#include <vtkCellOps.h>
#include <Elements/Triangle.h>
#include <Elements/ShapeFunctions.h>

#ifdef FOUND_GINAC
#  include <ginac/ginac.h>
USING_NAMESPACE(GiNaC);
#endif // FOUND_GINAC

// Return total number of shape functions for a triangle
int vtkShoeElemTriangle::GetTriangleNumberOfShapeFunctions(const int order[3])
{
	int p = order[0];
	return 3 + 3*(p-1) + (p-1)*(p-2)/2;
}

// p space
void vtkShoeElemTriangle::EvaluateShapeFunctionsTriangleMaxTotalOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3])
{
	double l[3];
	l[0] = .5*(1-r[0]-r[1]/sqrt(3.));
	l[1] = .5*(1+r[0]-r[1]/sqrt(3.));
	l[2] = r[1]/sqrt(3.);
	int index=0;

	// Corner modes
	shape[index++] = l[0];
	shape[index++] = l[1];
	shape[index++] = l[2];
	// Edge modes
	//Edge 0-1
	double x[10];
	double y[10];
	x[0]=l[1]-l[0];
	for(int cnt1=1;cnt1<order[0];cnt1++)
	{
		x[cnt1] = x[cnt1-1]*(l[1]-l[0]);
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++] = l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
	}
	//Edge 1-2
	x[0]=l[2]-l[1];
	for(int cnt1=1;cnt1<order[0];cnt1++)
	{
		x[cnt1] = x[cnt1-1]*(l[2]-l[1]);
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++] = l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
	}
	//Edge 2-0
	x[0]=l[0]-l[2];
	for(int cnt1=1;cnt1<order[0];cnt1++)
	{
		x[cnt1] = x[cnt1-1]*(l[0]-l[2]);
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++] = l[0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
	}
	
	// Face Modes
	x[0]=l[1]-l[0];
	for(int xcnt=1;xcnt<order[0];xcnt++)
	{
		x[xcnt] = x[xcnt-1]*(l[1]-l[0]);
	}
	y[0]=2*l[2]-1;
	for(int xcnt=1;xcnt<order[0];xcnt++)
	{
		y[xcnt] = y[xcnt-1]*(2*l[2]-1);
	}
	for(int cnt=0;cnt<=order[0]-3;cnt++)
	{
		for(int cnt1=cnt;cnt1>=0;cnt1--)
		{
			shape[index++] = l[0]*l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y);
		}
	}
}

void vtkShoeElemTriangle::EvaluateShapeFunctionDerivativesTriangleMaxTotalOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3]) // p space
{
/*	int no = cell.GetCellOps()->GetNumberOfShapeFunctions(order);
	
	double l[3];
	l[0] = .5*(1-r[0]-r[1]/sqrt(3.));
	l[1] = .5*(1+r[0]-r[1]/sqrt(3.));
	l[2] = r[1]/sqrt(3.);

	double ld[3][2];
	ld[0][0] = -.5; ld[0][1] = -.5/sqrt(3.);
	ld[1][0] = .5; ld[1][1] = -.5/sqrt(3.);
	ld[2][0] = 0.; ld[2][1] = 1/sqrt(3.);
	
	int index=0;

	// Corner modes
	
	shape[index] = ld[0][0];
	shape[no+index++] = ld[0][1];
	shape[index] = ld[1][0];
	shape[no+index++] = ld[1][1];
	shape[index] = ld[2][0];
	shape[no+index++] = ld[2][1];


	// Edge modes
	
	//Edge 0-1
	double x[10];
	double y[10];
	x[0]=l[1]-l[0];
	for(int cnt1=1;cnt1<order[0];cnt1++)
	{
		x[cnt1] = x[cnt1-1]*(l[1]-l[0]);
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index] = ld[0][0]*l[1]*psi[cnt](x) + l[0]*ld[1][0]*psi[cnt](x) + l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[1][0]-ld[0][0]);
		shape[no+index++] = ld[0][1]*l[1]*psi[cnt](x) + l[0]*ld[1][1]*psi[cnt](x) + l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[1][1]-ld[0][1]);
	}
	//Edge 1-2
	x[0]=l[2]-l[1];
	for(int cnt1=1;cnt1<order[0];cnt1++)
	{
		x[cnt1] = x[cnt1-1]*(l[2]-l[1]);
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index] = ld[1][0]*l[2]*psi[cnt](x) + l[1]*ld[2][0]*psi[cnt](x) + l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[2][0]-ld[1][0]);
		shape[no+index++] = ld[1][1]*l[2]*psi[cnt](x) + l[1]*ld[2][1]*psi[cnt](x) + l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[2][1]-ld[1][1]);
	}
	//Edge 2-0
	x[0]=l[0]-l[2];
	for(int cnt1=1;cnt1<order[0];cnt1++)
	{
		x[cnt1] = x[cnt1-1]*(l[0]-l[2]);
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index] = ld[0][0]*l[2]*psi[cnt](x) + l[0]*ld[2][0]*psi[cnt](x) + l[0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[0][0]-ld[2][0]);
		shape[no+index++] = ld[0][1]*l[2]*psi[cnt](x) + l[0]*ld[2][1]*psi[cnt](x) + l[0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[0][1]-ld[2][1]);
	}
	
	// Face Modes
	x[0]=l[1]-l[0];
	for(int xcnt=1;xcnt<order[0];xcnt++)
	{
		x[xcnt] = x[xcnt-1]*(l[1]-l[0]);
	}
	y[0]=2*l[2]-1;
	for(int xcnt=1;xcnt<order[0];xcnt++)
	{
		y[xcnt] = y[xcnt-1]*(2*l[2]-1);
	}
	
	for(int cnt=0;cnt<=order[0]-3;cnt++)
	{
		for(int cnt1=cnt;cnt1>=0;cnt1--)
		{
//	shape[index++] = l[0]*l[1]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y);
			shape[index++] = ld[0][0]*l[1]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*ld[1][0]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*ld[2][0]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*l[2]*leg[cnt1](x)*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1](y);
			shape[index++] = ld[0][1]*l[1]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*ld[1][1]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*ld[2][1]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*l[2]*leg[cnt1](x)*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1](y);
			shape[index++] = ld[0][2]*l[1]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*ld[1][2]*l[2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*ld[2][2]*leg[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*leg[cnt-cnt1](y) + l[0]*l[1]*l[2]*leg[cnt1](x)*vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1](y);
		}
	}
*/
}

