/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vector>

#include <vtkMath.h>
#include <vtkDoubleArray.h>

#include <vtkCellOps.h>
#include <vtkShoeMeshIterator.h>
#include <vtkAdaptiveTessellator.h>
#include <vtkShoeMeshSubdivisionAlgorithm.h>
#include <Elements/ShapeFunctions.h>
#include <Elements/CriticalPoints.h>
#include <Elements/Curve.h>

#ifdef FOUND_GINAC
USING_NAMESPACE(GiNaC);
#endif // FOUND_GINAC

// It takes \a o + 1 points to specify an order \a o curve.
// Edges associate 2 shape functions with their endpoints (corner vertices).
// That leaves \a o - 1 shape functions to be stored in the edge node.
// Curves are 1D, so they only use the first entry of \a o.
int vtkShoeElemCurve::GetNumberOfEdgeModesPerNode( int id, const int o[3] )
{
	return o[0] - 1;
}

int vtkShoeElemCurve::GetNumberOfShapeFunctions(const int order[3])
{
	return order[0]+1;
}

// Permute (negate) polynomial amplitudes as necessary
// Both the Legendre MaxTotalOrder and the Lagrange Tensor product polynomials have the same rule! Cool-O!
void vtkShoeElemCurve::GetPermutedEdgeSigns( vtkstd::vector<bool>& signs, int, bool node_permutation, const int order[3] )
{
	int total_dof = order[0] - 1;
	signs.resize( total_dof );

	vtkstd::fill( signs.begin(), signs.end(), false );

	if ( node_permutation )
		for ( int i=1; i<total_dof; i+=2 )
			signs[i] = true;
}

void vtkShoeElemCurve::EvaluateShapeFunctionsMaxTotalOrderLegendre( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	static double rpow[15];
	//static double phi_eval[15];
	double R = r[0];

	rpow[0] = R;
	for( int cnt=1; cnt<order[0]; )
	{
		rpow[cnt] = rpow[cnt-1]*R;
		cnt++;
		shape[cnt] = vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
	}
	R /= 2.;
	shape[0] = .5-R;
	shape[1] = .5+R;
}

void vtkShoeElemCurve::EvaluateShapeFunctionDerivativesMaxTotalOrderLegendre( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	static double rpow[15];
	shape[0] = -0.5;
	shape[1] =  0.5;
	double R = r[0];

	rpow[0] = R;
	for ( int o=1; o<order[0]; )
	{
		rpow[o] = rpow[o-1]*R;
		++o;
		shape[o] = vtkShoeElemShapeFunctions::EvaluatePhiDerivative[o](rpow);
	}
}

void vtkShoeElemCurve::EvaluateShapeFunctionsTensorLagrange( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	double R = r[0];
	shape[0] = R*(R - 1.)/2;
	shape[1] = R*(R + 1.)/2;
	if ( order[0] > 1 )
		shape[2] = 1. - R*R;
}

void vtkShoeElemCurve::EvaluateShapeFunctionDerivativesTensorLagrange( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	double R = r[0];
	shape[0] = R - 0.5;
	shape[1] = R + 0.5;
	if ( order[0] > 1 )
	  {
		shape[2] = -2*R;
		}
}

void vtkShoeElemCurve::GetPolynomialMaxTotalOrderLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector< vtkstd::vector< vtkstd::vector<int> > >& powers,  vtkShoeMeshIterator& cell_in, const int field_num)
{
#ifdef FOUND_GINAC
	const int *orders = cell_in.GetCellFunctionOrder(field_num);
	int dofcounter=0;
	int maxorder=0;
	cell_in.CacheFunction(field_num);
	vtkDoubleArray *field = cell_in.GetCachedFunction(field_num);
	double *dof = field->GetPointer(0);			
	
	int components = field->GetNumberOfComponents();

	symbol x("x");
	ex xp = 1+x;
	ex xn = 1-x;

	for(int i=0;i<field->GetNumberOfTuples()*field->GetNumberOfComponents();i++)
	{
		vtkstd::cout<<dof[i]<<" ";

	}
	vtkstd::cout<<vtkstd::endl;

	for(int CompCnt=0;CompCnt<components;CompCnt++)
	{
		dofcounter=CompCnt;
		ex finalPoly;
		// The corner shape function
		finalPoly += (dof[dofcounter]*.5*xn).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.5*xp).expand();
		dofcounter += components;

		// The edge shape functions
		// edge y=-1
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand().evalf();
			dofcounter += components;

		}
		finalPoly = finalPoly.expand().evalf();

	//	if(finalPoly.nops()==0)
	//		return;
		int deg[2];
		//     generate the coefficient and power vector
		vtkstd::vector<double> Tcoeff;
		vtkstd::vector<vtkstd::vector<int> > Tpowers;
		for(int cnt=0;cnt<finalPoly.nops();cnt++)
		{
			ex temp;	
			// Get the monomial
			temp=finalPoly.op(cnt); 
			// Find the degrees of the variables
			vtkstd::vector<int> tempVec;
			tempVec.push_back(temp.degree(x));
			Tpowers.push_back(tempVec);
			// Find the coefficient
			numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.)).evalf());
			Tcoeff.push_back(tempnum.to_double());
		}
	//	cout<<Tpowers.size()<<"***"<<Tcoeff.size()<<endl;
		powers.push_back(Tpowers);
		coeff.push_back(Tcoeff);
	}
// Sort the powers
/*	int n = powers.size();
	for(int cnt=0;cnt<n-1;cnt++)
	{
		for(int cnt1=cnt+1;cnt1<powers.size();cnt1++)
		{
			if(powers[cnt1]<powers[cnt])
			{	
				int temp1;
				double temp2;
				temp1=powers[cnt];
				powers[cnt]=powers[cnt1];
				powers[cnt1]=temp1;
				temp2=coeff[cnt];
				coeff[cnt]=coeff[cnt1];
				coeff[cnt1]=temp2;
			}
		}
	}
*/	
#endif // FOUND_GINAC
}
