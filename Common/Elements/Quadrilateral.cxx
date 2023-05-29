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
#include <vtkPointData.h>
#include <vtkIdTypeArray.h>

#include <vtksnlConfigure.h>
#include <vtkCellOps.h>
#include <vtkShoeMeshIterator.h>
#include <vtkAdaptiveTessellator.h>
#include <vtkShoeMeshSubdivisionAlgorithm.h>
#include <Elements/CriticalPoints.h>
#include <Elements/ShapeFunctions.h>
#include <Elements/Quadrilateral.h>
#include <Elements/Hexahedron.h>
#include <Elements/Generic.h>

#ifdef FOUND_GINAC
USING_NAMESPACE(GiNaC);
#endif // FOUND_GINAC

using namespace shoe;

// Permute (negate) Lagrange polynomial amplitudes as necessary for quadrilateral edge with given ID
void vtkShoeElemQuadrilateral::GetPermutedEdgeSignsTensorLagrange( vtkstd::vector<bool>& signs, int edge_ID, bool node_permutation, const int order[3] )
{
	// FIXME: Assumes order is the same in both parameters AND function is either linear or quadratic
	if ( order[0] == 1 )
		return ; // it's linear

	signs.resize( order[0] - 1 );
	*signs.begin() = node_permutation;
}

void vtkShoeElemQuadrilateral::GetPermutedFaceIndicesTensorLagrange( vtkstd::vector<int>& indices, vtkstd::vector<bool>& signs, int face_ID, int node_permutation, const int order[3] )
{
	// FIXME: Assumes no face nodes
	indices.clear();
	signs.clear();
}

// Permute (negate) Legendre polynomial amplitudes as necessary for quadrilateral edge with given ID
void vtkShoeElemQuadrilateral::GetPermutedEdgeSignsLegendre( vtkstd::vector<bool>& signs, int edge_ID, bool node_permutation, const int order[3] )
{

	// Determine order for parameter for this edge
	int total_dof = (edge_ID % 2) ? order[1] : order[0];
	
	if ( (int) signs.size() != total_dof - 1 )
		signs.resize( total_dof - 1 );

	vtkstd::fill( signs.begin(), signs.end(), false );
	
	if (node_permutation)
		for (int i = 3; i <= total_dof; i += 2)
			signs[i-2] = true;
}

// Based on the permutation bits for the face, determine the indices and signs
// of quadrilateral face modes after amplitude swaps and sign changes.
void vtkShoeElemQuadrilateral::GetPermutedFaceIndicesPLegendre( vtkstd::vector<int>& indices, vtkstd::vector<bool>& signs, int face_ID, int node_permutation, const int order[3] )
{
	vtkShoeElemHexahedron::GetPermutedFaceIndicesPLeg( indices, signs, 0, node_permutation, order );
}

// Return total number of shape functions for quad with lagrange tensor product interpolation.
int vtkShoeElemQuadrilateral::GetNumberOfShapeFunctionsTensorLagrange(const int order[3])  // p space
{
	int no =  4 + 4*(order[0]-1);
	return no;
}

// Return total number of shape functions for quad with legendre interpolation and maximum total order space.
int vtkShoeElemQuadrilateral::GetNumberOfShapeFunctionsMaxTotalOrderLegendre(const int order[3])  // p space
{
	int no =  4 + 4*(order[0]-1);
	if(order[0]>=4)
		no += (order[0]-2)*(order[0]-3)/2;
	return no;
}

void vtkShoeElemQuadrilateral::EvaluateShapeFunctionsMaxTotalOrderLegendre(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3])
{
	static double rpow[15];
	static double spow[15];
	rpow[0] = r[0];
	spow[0] = r[1];
	for(int cnt=1;cnt<order[0];cnt++)
	{
		rpow[cnt] = rpow[cnt-1]*r[0];
		spow[cnt] = spow[cnt-1]*r[1];
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
		vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
	}
//	r=r[0];s=r[2];t=r[2];
	double terms[6];
	int index=0;
	terms[0] = 1+r[0];      //  1+r
	terms[1] = 1-r[0];	// 1-r
	terms[2] = 1+r[1];	// 1+s
	terms[3] = 1-r[1];	// 1-s
	// Corner Nodes
	shape[index++]=(.25*terms[1]*terms[3]);    // .125*(1-r)(1-s)
	shape[index++]=(.25*terms[0]*terms[3]);   //  .125*(1+r)(1-s)
	shape[index++]=(.25*terms[0]*terms[2]);   //  .125*(1+r)(1+s)
	shape[index++]=(.25*terms[1]*terms[2]);   //  .125*(1-r)(1+s)
	// Edge Modes
	// Edge 0 - 1  (1-2)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]);      //  sn*tn
	}
	// (dge 1 - 2  (2,3)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]);     // rp*tn
	}
	// Edge 2 - 3  (3,4)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]);    // sp*tn
	}
	// Edge 3 - 0  (4,1)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]);     // rn*tn
	}
	// Face Modes
	for(int cnt=4;cnt<=order[0];cnt++)
	{
		for(int cnt1=cnt-2;cnt1>=2;cnt1--)
		{
			//				cout<<cnt1<<" "<<cnt-cnt1<<endl;
			shape[index++]=(vtkShoeElemShapeFunctions::phiValue[0][cnt1]*vtkShoeElemShapeFunctions::phiValue[1][cnt-cnt1]);
		}
	}
}

void vtkShoeElemQuadrilateral::EvaluateShapeFunctionsTensorLegendre(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3])  // pq space
{
	static double rpow[15];
	static double spow[15];
	rpow[0] = r[0];
	spow[0] = r[1];
	for(int cnt=1;cnt<order[0];cnt++)
	{
		rpow[cnt] = rpow[cnt-1]*r[0];
	}
	for(int cnt=1;cnt<order[1];cnt++)
	{
		spow[cnt] = spow[cnt-1]*r[1];
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
	}
	for(int cnt=2;cnt<=order[1];cnt++)
	{
		vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
	}
	double terms[6];
	int index=0;
	terms[0] = 1+r[0];      //  1+r
	terms[1] = 1-r[0];	// 1-r
	terms[2] = 1+r[1];	// 1+s
	terms[3] = 1-r[1];	// 1-s
	// Corner Nodes
	shape[index++]=(.25*terms[1]*terms[3]);    // .125*(1-r)(1-s)
	shape[index++]=(.25*terms[0]*terms[3]);   //  .125*(1+r)(1-s)
	shape[index++]=(.25*terms[0]*terms[2]);   //  .125*(1+r)(1+s)
	shape[index++]=(.25*terms[1]*terms[2]);   //  .125*(1-r)(1+s)
	
	// Edge 0 - 1  (1-2)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]);      //  sn*tn
	}
	// Edge 1 - 2  (2,3)
	for(int cnt=2;cnt<=order[1];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]);     // rp*tn
	}
	// Edge 2 - 3  (3,4)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]);    // sp*tn
	}
	// Edge 3 - 0  (4,1)
	for(int cnt=2;cnt<=order[1];cnt++)
	{
		shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]);     // rn*tn
	}
	
	
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		for(int cnt1=2;cnt1<=order[1];cnt1++)
		{
			shape[index++]=(vtkShoeElemShapeFunctions::phiValue[0][cnt]*vtkShoeElemShapeFunctions::phiValue[1][cnt1]);
		}
	}
	
}

void vtkShoeElemQuadrilateral::EvaluateShapeFunctionsTruncatedOrderLegendre(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]) // ~pq space
{
	static double rpow[15];
	static double spow[15];
	rpow[0] = r[0];
	spow[0] = r[1];
	for(int cnt=1;cnt<order[0];cnt++)
	{
		rpow[cnt] = rpow[cnt-1]*r[0];
	}
	for(int cnt=1;cnt<order[1];cnt++)
	{
		spow[cnt] = spow[cnt-1]*r[1];
	}
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
	}
	for(int cnt=2;cnt<=order[1];cnt++)
	{
		vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
	}
//	r=r[0];s=r[2];t=r[2];
	double terms[6];
	int index=0;
	terms[0] = 1+r[0];      //  1+r
	terms[1] = 1-r[0];	// 1-r
	terms[2] = 1+r[1];	// 1+s
	terms[3] = 1-r[1];	// 1-s

	// Corner Nodes
	shape[index++]=(.25*terms[1]*terms[3]);    // .125*(1-r)(1-s)
	shape[index++]=(.25*terms[0]*terms[3]);   //  .125*(1+r)(1-s)
	shape[index++]=(.25*terms[0]*terms[2]);   //  .125*(1+r)(1+s)
	shape[index++]=(.25*terms[1]*terms[2]);   //  .125*(1-r)(1+s)

	// Edge Modes
		
	// Edge 0 - 1  (1-2)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]);      //  
	}
	// Edge 1 - 2  (2,3)
	for(int cnt=2;cnt<=order[1];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]);     // 
	}
	// Edge 2 - 3  (3,4)
	for(int cnt=2;cnt<=order[0];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]);    // 
	}
	// Edge 3 - 0  (4,1)
	for(int cnt=2;cnt<=order[1];cnt++)
	{
		shape[index++]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]);     // 
	}
	
	// Face Modes
	for(int cnt=4;cnt<=order[0];cnt++)
	{
		for(int cnt1=cnt-2;cnt1>=2;cnt1--)
		{
			if(cnt-cnt1>order[1])
				continue;
			shape[index++]=vtkShoeElemShapeFunctions::phiValue[0][cnt1]*vtkShoeElemShapeFunctions::phiValue[1][cnt-cnt1];
		}
	}
}

void vtkShoeElemQuadrilateral::EvaluateShapeFunctionDerivativesMaxTotalOrderLegendre( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	int no = cell.GetCellOps()->GetNumberOfShapeFunctions(order);
	static double rpow[15];
	static double spow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  for ( int cnt = 1; cnt < order[0]; ++cnt )
    {
    rpow[cnt] = rpow[cnt-1] * r[0];
    spow[cnt] = spow[cnt-1] * r[1];
    }
  for ( int cnt = 2; cnt <= order[0]; ++cnt)
    {
    vtkShoeElemShapeFunctions::phiValue[0][cnt] = vtkShoeElemShapeFunctions::EvaluatePhi[cnt]( rpow );
    vtkShoeElemShapeFunctions::phiValue[1][cnt] = vtkShoeElemShapeFunctions::EvaluatePhi[cnt]( spow );
    vtkShoeElemShapeFunctions::dphiValue[0][cnt] = vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt]( rpow );
    vtkShoeElemShapeFunctions::dphiValue[1][cnt] = vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt]( spow );
    }

	double terms[6];
	int index=0;
	terms[0] = 1 + r[0];	// 1+r
	terms[1] = 1 - r[0];	// 1-r
	terms[2] = 1 + r[1];	// 1+s
	terms[3] = 1 - r[1];	// 1-s
	
	// Corner Nodes
	shape[index] = (.25*-1*terms[3]);      //  .125*(1-r)(1-s)(1-t)
	shape[no+index++] = (.25*terms[1]*-1);
	
	shape[index]=(.25*terms[3]);           //  .125*(1+r)(1-s)(1-t)
	shape[no+index++]=(.25*terms[0]*-1);   //  .125*(1+r)(1-s)(1-t)
	
	shape[index]=(.25*terms[2]);           //  .125*(1+r)(1+s)(1-t)
	shape[no+index++]=(.25*terms[0]);      //  .125*(1+r)(1+s)(1-t)
	
	shape[index]=(.25*-1*terms[2]);        //  .125*(1-r)(1+s)(1-t)
	shape[no+index++]=(.25*terms[1]);      //  .125*(1-r)(1+s)(1-t)

  // Edge Modes
  // Edge 0 - 1  (1-2)
  for(int cnt=2;cnt<=order[0];cnt++)
    {
    shape[index]=(.5*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]);      //  sn*tn
    shape[no+index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1);      //  sn*tn
    }
  // (edge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[0];cnt++)
    {
    shape[index]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]);     // rp*tn
    shape[no+index++]=(.5*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]);     // rp*tn
    }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
    {
    shape[index]=(.5*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]);    // sp*tn
    shape[no+index++]=(.5*vtkShoeElemShapeFunctions::phiValue[0][cnt]);    // sp*tn
    }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[0];cnt++)
    {
    shape[index]=(.5*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1);     // rn*tn
    shape[no+index++]=(.5*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]);     // rn*tn
    }

	// Face Modes
	for ( int cnt = 4; cnt <= order[0]; ++cnt )
    {
    for ( int cnt1 = cnt - 2; cnt1 >= 2; --cnt1)
      {
      shape[index] = vtkShoeElemShapeFunctions::dphiValue[0][cnt1] * vtkShoeElemShapeFunctions::phiValue[1][cnt - cnt1];
      shape[no + index++] = vtkShoeElemShapeFunctions::phiValue[0][cnt1] * vtkShoeElemShapeFunctions::dphiValue[1][cnt - cnt1];
      }
    }
}

void vtkShoeElemQuadrilateral::EvaluateNormalOnFace( double *norm, vtkShoeMeshIterator& cell, const double r[3], const int FaceId )
{
  vtkShoeElemHexahedron::EvaluateNormalOnFace( norm, cell, r, FaceId );
}


void vtkShoeElemQuadrilateral::EvaluateShapeFunctionsTensorLagrange( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	// FIXME: Assumes order[0] = order[1] AND order is 1 or 2.

	double omr = 1. - r[0];
	double oms = 1. - r[1];
	double opr = 1. + r[0];
	double ops = 1. + r[1];

	shape[0] = omr*oms*(r[0] + r[1] + 1.)/-4;
	shape[1] = opr*oms*(r[0] - r[1] - 1.)/ 4;
	shape[2] = opr*ops*(r[0] + r[1] - 1.)/ 4;
	shape[3] = omr*ops*(r[0] - r[1] + 1.)/-4;

	if ( order[0] == 1 )
		return;

	double omr2 = 1. - r[0]*r[0];
	double oms2 = 1. - r[1]*r[1];

	shape[4] = omr2*oms/2;
	shape[5] = opr*oms2/2;
	shape[6] = omr2*ops/2;
	shape[7] = omr*oms2/2;
}

void vtkShoeElemQuadrilateral::EvaluateShapeFunctionDerivativesTensorLagrange( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] )
{
	// FIXME: Assumes order[0] = order[1] AND order is 1 or 2.

	int no = order[0] == 1 ? 4 : 8;

	double omr = 1. - r[0];
	double oms = 1. - r[1];
	double opr = 1. + r[0];
	double ops = 1. + r[1];
	double trms = 2*r[0] - r[1];
	double trps = 2*r[0] + r[1];
	double tsmr = 2*r[1] - r[0];
	double tspr = 2*r[1] + r[0];

	shape[0] = oms*trps/4; shape[no+0] = omr*tspr/4;
	shape[1] = oms*trms/4; shape[no+1] = opr*tsmr/4;
	shape[2] = ops*trps/4; shape[no+2] = opr*tspr/4;
	shape[3] = ops*trms/4; shape[no+3] = omr*tsmr/4;

	if ( order[0] == 1 )
		return;

	double omr2 = 1. - r[0]*r[0];
	double oms2 = 1. - r[1]*r[1];

	shape[4] = -r[0]*oms;  shape[no+4] = omr2/-2;
	shape[5] = oms2/ 2;    shape[no+5] = -r[1]*opr;
	shape[6] = -r[0]*ops;  shape[no+6] = omr2/ 2;
	shape[7] = oms2/-2;    shape[no+7] = -r[1]*omr;
}

// This function returns the polynomial function for quadrilateral for total maximum order space.
// Input :
// orders = orders in each direction
// dof = one array containing all the degrees of freedom
// output:
// coeff = coefficient vector
// powers = powers vector
void vtkShoeElemQuadrilateral::GetPolynomialMaxTotalOrderLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num)
{
#ifdef FOUND_GINAC
	const int *orders = cell_in.GetCellFunctionOrder(field_num);
	cell_in.CacheFunction(field_num);
	vtkDoubleArray *field = cell_in.GetCachedFunction(field_num);
	int components = field->GetNumberOfComponents();
	
 	double *dof = field->GetPointer(0);
	int dofcounter=0;
	int maxorder=0;
	symbol x("x");symbol y("y");
	ex xp = 1+x;
	ex xn = 1-x;
	ex yp = 1+y;
	ex yn = 1-y;

	for(int CompCnt=0;CompCnt<components;CompCnt++)
	{
		dofcounter=CompCnt;
		ex finalPoly;

		// The corner shape function
		finalPoly += (dof[dofcounter]*.25*xn*yn).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.25*xp*yn).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.25*xp*yp).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.25*xn*yp).expand();
		dofcounter += components;

		// The edge shape functions

		// edge y=-1
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
			dofcounter += components;

		}
		// edge x=1;
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*xp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
			dofcounter += components;
		}
		//	cout<<"edge y=1"<<endl;
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
			dofcounter += components;
		}
		//	cout<<"edge x=-1"<<endl;
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*xn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
			dofcounter += components;
		}
		//       cout<<"evaluating face modes "<<endl;
		// face modes exponent sequence is {2,2} {3,2} {2,3} {4,2}  {3,3} {2,4} .....
		for(int cnt=4;cnt<=orders[0];cnt++)
		{
			//		cout<<orders[4]<<endl;
			for(int cnt1=cnt-2;cnt1>=2;cnt1--)
			{
				//		cout<<cnt1<<" "<<cnt-cnt1<<" "<<dof[dofcounter]<<endl;
				finalPoly += (dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,x)*vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1,y)).expand();
				dofcounter += components;

			}
		}
		finalPoly = finalPoly.expand().evalf();
		vtkstd::vector<double> Tcoeff;
		vtkstd::vector<vtkstd::vector<int> > Tpowers;
		int deg[2];
		//     generate the coefficient and power vector
		for(int cnt=0;cnt<finalPoly.nops();cnt++)
		{
			vtkstd::vector<int> tempvec;
			ex temp;	
			// Get the monomial
			temp=finalPoly.op(cnt);  
			// Find the degrees of the variables
			tempvec.push_back(temp.degree(x));
			tempvec.push_back(temp.degree(y));
			// Find the coefficient
			numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.,y==1)).evalf());
			if(tempnum>-1.e-8 && tempnum<1.e-8)
			{
				tempvec.clear();
				continue;
			}
			Tcoeff.push_back(tempnum.to_double());
			Tpowers.push_back(tempvec);

		}
		coeff.push_back(Tcoeff);
		powers.push_back(Tpowers);
		

	}
#endif // FOUND_GINAC
}

// This function returns the polynomial function for a quadrilateral for tensor product space.
// Input : 
// orders = orders for each higher order node // The maximum order in each direction
// dof = one array containing all the degrees of freedom
// output:
// coeff = coefficient vector
// powers = powers vector

void vtkShoeElemQuadrilateral::GetPolynomialTensorLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector< vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num)
{
#ifdef FOUND_GINAC
	int dofcounter=0;
	int maxorder=0;
	const int *orders = cell_in.GetCellFunctionOrder(field_num);
	cell_in.CacheFunction(field_num);
	vtkDoubleArray *field = cell_in.GetCachedFunction(field_num);

	if (field->GetNumberOfComponents()>1)
	{
		vtkstd::cout<<"Vector fields are not supported by critical point functions yet "<<vtkstd::endl;
		return;
	}
 
	int components = field->GetNumberOfComponents();
 	double *dof = field->GetPointer(0);			
	symbol x("x");symbol y("y");
	ex xp = 1+x;
	ex xn = 1-x;
	ex yp = 1+y;
	ex yn = 1-y;

	for(int CompCnt=0;CompCnt<components;CompCnt++)
	{

		dofcounter=CompCnt;

		ex finalPoly;

		// The corner shape function
		finalPoly += (dof[dofcounter]*.25*xn*yn).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.25*xp*yn).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.25*xp*yp).expand();
		dofcounter += components;
		finalPoly += (dof[dofcounter]*.25*xn*yp).expand();
		dofcounter += components;

		// The edge shape functions

		// edge y=-1
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
			dofcounter += components;

		}
		// edge x=1;
		for(int cnt=2;cnt<=orders[1];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*xp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
			dofcounter += components;
		}
		//	cout<<"edge y=1"<<endl;
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
			dofcounter += components;
		}
		//	cout<<"edge x=-1"<<endl;
		for(int cnt=2;cnt<=orders[1];cnt++)
		{
			finalPoly += (0.5*dof[dofcounter]*xn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
			dofcounter += components;
		}
		// face modes exponent sequence is {2,2} {3,2} {2,3} {4,2}  {3,3} {2,4} .....
		for(int cnt=2;cnt<=orders[0];cnt++)
		{
			for(int cnt1=2;cnt1<=orders[1];cnt1++)
			{
				//			cout<<cnt<<" "<<cnt1<<endl;
				finalPoly += (dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt,x)*vtkShoeElemCriticalPoints::eval_phi(cnt1,y)).expand();
				dofcounter += components;

			}
		}
		finalPoly = finalPoly.expand().evalf();
		int deg[2];
		vtkstd::vector<double> Tcoeff;
		vtkstd::vector<vtkstd::vector<int> > Tpowers;


		//     generate the coefficient and power vector
		for(int cnt=0;cnt<finalPoly.nops();cnt++)
		{
			vtkstd::vector<int> tempvec;
			ex temp;	
			// Get the monomial
			temp=finalPoly.op(cnt);  
			//		cout<<temp<<endl;
			// Find the degrees of the variables
			tempvec.push_back(temp.degree(x));
			tempvec.push_back(temp.degree(y));
			Tpowers.push_back(tempvec);

			// Find the coefficient
			numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.,y==1.)).evalf());
			Tcoeff.push_back(tempnum.to_double());
		}	

		powers.push_back(Tpowers);
		coeff.push_back(Tcoeff);
		Tpowers.clear();
		Tcoeff.clear();
	}
#endif // FOUND_GINAC
}

void vtkShoeElemQuadrilateral::GetPolynomialTruncatedTotalOrderLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, int *orders, vtkstd::vector<double>& dof)
{
#ifdef FOUND_GINAC

#endif // FOUND_GINAC
}

void vtkShoeElemQuadrilateral::GetCriticalPoints(vtkShoeMeshIterator& iter, int field_no) 
{
	int FaceCorners[4][3] = {{-1,-1,0},
		        {1,-1,0},
	        	{1,1,0},
	          	{-1,1,0}};
	int EdgeCorners[2][3] = { {-1,0,0},
				  {1,0,0}};
	
        vtkstd::vector<vtkstd::vector<double> > coef;         //  Vector storing the coefficients.
	vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >powers;  // Vector storing the exponents of the field polynomial
	// Get the cellops
	const vtkCellOps *ops;        
	ops = iter.GetCellOps();
	

	// Get the degrees of freedom for the field
	vtkDoubleArray *data = iter.GetCachedFunction(field_no);
	int components = data->GetNumberOfComponents();

	shoe::CellShape eleShape = iter.GetCell()->Def->GetType().DomainShape;
	
	 //  Get the field polynomial
	ops->GetFieldPolynomial(coef,powers,iter,field_no);   	
	int max[3];
	
	// Find the critical points of the field polynomial
	
	double param_bounds[3] = {-1,1,0};  // the range is  -1 to 1.
	

	vtkstd::vector<vtkstd::vector<double> >::iterator ic = coef.begin();
	vtkstd::vector<vtkstd::vector<double> > cornerValues;    // Vector storing the corner values
	// Get the corner values.
	for(int cnt=0;cnt<ops->NumberOfPoints;cnt++)
	{
		vtkstd::vector<double> tempvec;
		for(int i=0;i<components;i++)
		{
			tempvec.push_back(data->GetValue(components*cnt+i));
		}
		cornerValues.push_back(tempvec);
	}
	
	// Vectors storing extrema and types for each of the components.
	vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > > Vextrema;
	vtkstd::vector<vtkstd::vector<int> > Vtypes;

	ops->FindFieldCriticalPoints(Vextrema,Vtypes,max,coef,powers,ops->EmbeddingDimension);
	ops->FindBoundedCriticalValues(Vextrema,Vtypes,ops->EmbeddingDimension,param_bounds,eleShape,cornerValues);
	
	vtkstd::vector<vtkstd::vector<double> > ListOfExtrema;

	// Change the corner point ids to the corner point parametric coordinates. Prepare a list
	// of extrema to pass to the DofNodeExtrema functions.
	for(int Ccnt=0;Ccnt<components;Ccnt++)
	{
		for(int i=0;i < (int) Vextrema[Ccnt].size();i++)
		{
	//		cout<<Vextrema[Ccnt][i][0]<<" "<<Vextrema[Ccnt][i][1]<<" "<<Vextrema[Ccnt][i][2]<<" "<<Vextrema[Ccnt][i][3]<<" "<<Vextrema[Ccnt][i][4]<<endl;
			if((((int)Vextrema[Ccnt][i][Vextrema[Ccnt][i].size()-1])&7)==4 || ((int)(Vextrema[Ccnt][i][Vextrema[Ccnt][i].size()-1])&7)==5 )
			{
				int id = (int)Vextrema[Ccnt][i][0];
				for(int j=0;j<3;j++)
				{
					Vextrema[Ccnt][i][j]=FaceCorners[id][j];
				}

			}

			ListOfExtrema.push_back(Vextrema[Ccnt][i]);
		}
	}
		vtkCellDefinition tempDef;
	
	vtkIdType conn[10];
	
	const vtkIdType* this_cell_conn = iter.GetMesh()->GetCellConnectivity( iter.GetCellId() );
	// Insert the extrema in the data structure
	if(ListOfExtrema.size()>0)
		iter.GetMesh()->GetFunctionData()->InsertDofNodeExtrema(field_no,this_cell_conn[8],ListOfExtrema);

	// clear everything.
	Vextrema.clear();
	Vtypes.clear();
	ListOfExtrema.clear();
	cornerValues.clear();
	Vextrema.clear();
	coef.clear();
	powers.clear();
	uint32_t perm_e;
	for(int Ecnt=0;Ecnt<4;Ecnt++)
	{
		if(iter.GetMesh()->DoFunctionDofNodeExtremaExist(field_no,this_cell_conn[8+Ecnt]))
			continue;
	
		vtkstd::vector<vtkstd::vector<double> > Ecoef;
		vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > > Epowers;
//              Get boundary edge.
		ops->GetBoundaryEdge(tempDef,conn,perm_e,iter,Ecnt);
//              Insert it in a temporary mesh.
		vtkShoeMesh *tempMesh = vtkShoeMesh::New();
		tempMesh->SetPoints(iter.GetMesh()->GetPoints());
		tempMesh->SetFunctionData(iter.GetMesh()->GetFunctionData());
		tempMesh->SetGeometryData(iter.GetMesh()->GetGeometryData());
		tempMesh->SetFieldData(iter.GetMesh()->GetFieldData());
		tempMesh->GetPointData()->ShallowCopy(iter.GetMesh()->GetPointData());
		vtkShoeMesh::CellDefIterator i = tempMesh->FindOrCreateCellDef(tempDef.GetType(),
						 tempDef.GetGeometricOrder(),tempDef.GetFieldOrder(field_no));
		tempMesh->InsertNextCell(i,conn,perm_e);
		vtkShoeMeshIterator Titer(tempMesh->Begin(),MeshOrder);
		Titer.GetCellOps()->GetFieldPolynomial(Ecoef,Epowers,Titer,field_no);
//		cout<<Ecoef.size()<<" "<<Epowers.size()<<endl;
		vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > >Eextrema;
		vtkstd::vector<vtkstd::vector<int> > Etypes;
//              Find the critical points
//
		for(int Ccnt=0;Ccnt<components;Ccnt++)
		{
			vtkstd::vector<vtkstd::vector<double> > ComponentExtrema;
			vtkstd::vector<int> ComponentTypes;


			ops->FindUnivariateCriticalPoints(ComponentExtrema,ComponentTypes,Ecoef[Ccnt],Epowers[Ccnt]);
			for(int i=0;i < (int) ComponentExtrema.size();i++)
			{
				for(int cnt1=0;cnt1<components;cnt1++)
				{
					double value=0.;
					for ( int j = 0; j < (int) Epowers[cnt1].size(); ++j )
					{
						value+=Ecoef[cnt1][j]*::pow(ComponentExtrema[i][0],Epowers[cnt1][j][0]);
					}

					ComponentExtrema[i].push_back(value);
				}
			}
//			tempvec.push_back(0);
//			tempvec.push_back(0);
			Eextrema.push_back(ComponentExtrema);
			Etypes.push_back(ComponentTypes);
			ComponentExtrema.clear();
			ComponentTypes.clear();
			
		}
		vtkstd::vector<vtkstd::vector<double> > EcornerValues;

		for(int cnt=0;cnt<2;cnt++)
		{
			vtkstd::vector<double> tempvec;
			for(int i=0;i<components;i++)
			{
				tempvec.push_back(data->GetValue(components*vtkShoeElemGeneric::edgeConnectivity[Ecnt][cnt]+i));
			}
			EcornerValues.push_back(tempvec);
		}
		//              Find bounded critical values.
		ops->FindBoundedCriticalValues(Eextrema,Etypes,ops->EmbeddingDimension-2,param_bounds,Curve,EcornerValues);

		for(int Ccnt=0;Ccnt<components;Ccnt++)
		{
		
			for ( int i = 0; i < (int) Eextrema[Ccnt].size(); ++i )
			{
				vtkstd::vector<double>::iterator tempVec = Eextrema[Ccnt][i].begin();
				tempVec++;
				Eextrema[Ccnt][i].insert(tempVec,0); 
				tempVec=Eextrema[Ccnt][i].begin();
				tempVec++;
				Eextrema[Ccnt][i].insert(tempVec,0);
				if(((int)Eextrema[Ccnt][i][Eextrema[Ccnt][i].size()-1]&7) == 4 || ((int)Eextrema[Ccnt][i][Eextrema[Ccnt][i].size()-1]&7)==5)
				{
					int id = (int)Eextrema[Ccnt][i][0];
					for(int j=0;j<3;j++)
						Eextrema[Ccnt][i][j]=EdgeCorners[id][j];
				}
				ListOfExtrema.push_back(Eextrema[Ccnt][i]);
			}
		}

		// Depermute the values .............................
		for ( int cnt = 0; cnt < (int) ListOfExtrema.size(); ++cnt )
		{
			if(iter.GetCell()->GetEdgePermutation( Ecnt ))
			{
				// cout<<"Permuting values "<<endl;
				ListOfExtrema[cnt][0] = -ListOfExtrema[cnt][0];
				
			}
		}
		iter.GetMesh()->GetFunctionData()->InsertDofNodeExtrema(field_no,this_cell_conn[8+Ecnt],ListOfExtrema);
		EcornerValues.clear();
		Eextrema.clear();
		Etypes.clear();
		Ecoef.clear();
		Epowers.clear();
		ListOfExtrema.clear();
		tempMesh->Delete();
	}
}

bool vtkShoeElemQuadrilateral::DoCriticalPointsExist(vtkShoeMeshIterator &iter, const vtkstd::vector<double>& coef, const vtkstd::vector<vtkstd::vector<int> >& powers)
{
	
	int dimension = iter.GetCellOps()->EmbeddingDimension;
	const int *order = iter.GetCellGeometricOrder();
	double sign=0;
	bool SignChanged=false;
	double para[2]= {-1.,-1.};
	double step = 1./order[0];
	for(int cnt=0;cnt<dimension;cnt++)
	{
		SignChanged=false;
		sign=0;
		for(para[0]=-1.;para[0]<=1.;para[0]+=step)
		{
			for(para[1]=-1.;para[1]<=1.;para[1]+=step)
			{
				double value=0.;
				for ( int cnt1 = 0; cnt1 < (int) coef.size(); ++cnt1 )
				{
					if(powers[cnt1][cnt]<1)
						continue;
					double term=1.;
					for(int l=0;l<2;l++)
					{
						if(l==cnt)
							term=term*powers[cnt1][l]*::pow(para[l],powers[cnt1][l]-1);
						else
							term=term*::pow(para[l],powers[cnt1][l]);
					}
					value += term*coef[cnt1];

				}
		//		cout<<"value "<<value<<endl;
				if(value==0)
					continue;
				if(sign==0)
				{
					sign=value;
					continue;
				}
				if(sign*value<0)
				{
					SignChanged=true;
					break;
				}
			}

			if(SignChanged==true)
				break;

		}
		if(SignChanged==false)
		{
			vtkstd::cout<<"No critical points "<<vtkstd::endl;
			return false;
		}
			
		SignChanged=false;

	}
	return true;
}

void vtkShoeElemQuadrilateral::GetBoundaryEdge( vtkCellDefinition& def_out, vtkIdType* connectivity_out, uint32_t& perm_out, vtkShoeMeshIterator& cell, int edge_number )
{
	// cell operations for current cell
	const vtkCellOps* cellOps = cell.GetCellOps();

	// get the current cell definition
	vtkShoeMesh::CellDefConstIterator this_def = cell.GetCell()->Def;

	// get geometric order of current cell and first field order as place holder
	const int* geom_order = this_def->GetGeometricOrder();
	const int* func_order = this_def->GetFieldOrder( 0 );
	int nfields = this_def->GetNumberOfFields();

	// adjust geometric orders
	int new_geom_order[3];
	new_geom_order[0] = geom_order[edge_number % 2];
	new_geom_order[1] = 0;
	new_geom_order[2] = 0;

	// make new cell definition
	def_out = vtkCellDefinition( Curve, this_def->GetInterpolant(), this_def->GetProductSpace(), 0, nfields, new_geom_order, func_order );

	// adjust function orders in new definition
	int new_func_order[3];
	for (int i = 0; i < nfields; i++)
	{
		func_order = this_def->GetFieldOrder( i );
		new_func_order[0] = func_order[edge_number% 2];
		new_func_order[1] = 0;
		new_func_order[2] = 0;
		def_out.SetFieldOrder( i, new_func_order );
	}

	// declare a cell spec for setting permutations
	vtkShoeMesh::CellSpec temp_spec;

	// permutation is the same as in the parent cell
	temp_spec.SetEdgePermutation( 0, cell.GetCell()->GetEdgePermutation( edge_number ) );
	perm_out = temp_spec.GetCellPermutation();

	// connectivity array of this cell
	const vtkIdType* this_cell_conn = cell.GetMesh()->GetCellConnectivity( cell.GetCellId() );
	
	// get the endpoints of the new cell
	const int* corners = cellOps->EdgeCornerNodes[ edge_number ];
	connectivity_out[0] = this_cell_conn[corners[0]];
	connectivity_out[1] = this_cell_conn[corners[1]];

	// add the edge node
	connectivity_out[2] = this_cell_conn[4+edge_number];
}

void vtkShoeElemQuadrilateral::GetPermutedCriticalPoints( vtkstd::vector<vtkstd::vector<double> >& extrema, vtkShoeMeshIterator& iter, const int field )
{
	if( ! iter.GetMesh()->DoFunctionExtremaExist(field) )
		{
		vtkGenericWarningMacro( << "Critical points for field " << field << " do not exist. Computing them now." ); 
		for ( vtkShoeMeshIterator it( iter.GetMesh()->Begin(),MeshOrder); it != iter.GetMesh()->End(); ++it )
			it.GetCellOps()->GetCellCriticalPoints( it, field );
		}

	// Get the extrema
	vtkFunctionData *extremaData = iter.GetMesh()->GetFunctionData()->GetExtrema();
	//Get the offset array
	vtkIdTypeArray *offsets = extremaData->GetOffsets(field);
	// Get the extrema values (coordinates, type and value)
	vtkDoubleArray *values = dynamic_cast<vtkDoubleArray *>( extremaData->GetValues(field) );
	int components = values->GetNumberOfComponents();
	double *exVals = new double[components];
	const vtkIdType* this_cell_conn = iter.GetMesh()->GetCellConnectivity( iter.GetCellId() );
	double tup[2];

	// ========= Volumetric critical points ============
	//
	// this code intentionally left blank
	//
	// ========= Face critical points ============
	//
	int node_permutation = iter.GetCell()->GetFacePermutation( 0 );
	int node_id = this_cell_conn[4+4]; // Get the node id of the higher order node
	offsets->GetTuple( node_id, tup ); // Get the tuple containing the offset and number of critical points

	// For each critical point in node node_id
	for ( int cnt=0; cnt < tup[1]; ++cnt )
		{
		values->GetTuple( int(tup[0]) + cnt, exVals );
		vtkShoeElemGeneric::DepermuteQuadParameters( exVals, node_permutation );

		vtkstd::vector<double> tempVec;
		for ( int i=0; i < 2; ++i )
			tempVec.push_back( exVals[i] );
		tempVec.push_back( 0. ); // last parametric coordinate is 0.

		for ( int cnt = 3; cnt < components; ++cnt )
			tempVec.push_back( exVals[cnt] );

		extrema.push_back(tempVec);
		tempVec.clear();
		}

	// ========= Edge critical points ============
	//
	bool perm_e;
	for ( int Ecnt=0; Ecnt < 4; Ecnt++ )
		{
		perm_e = iter.GetCell()->GetEdgePermutation( Ecnt );
		int node_id = this_cell_conn[ 4 + Ecnt ];
		offsets->GetTuple( node_id, tup );

		for ( int cnt=0; cnt < tup[1]; cnt++ )
			{
			vtkstd::vector<double> tempVec;

			if ( tup[0] == -1 )
				break;

			values->GetTuple( (int) tup[0]+cnt, exVals );
			if ( perm_e )
				exVals[0] = -exVals[0];
			
			if ( Ecnt % 2 )
				{ // odd  edges are vertical   (s parameter varies, r=+/-1)
				tempVec.push_back( double(1 - ((Ecnt>>1)<<1)) );
				tempVec.push_back( exVals[0] );
				}
			else
				{ // even edges are horizontal (r parameter varies, s=+/-1)
				tempVec.push_back( exVals[0] );
				tempVec.push_back( double(-1 + ((Ecnt>>1)<<1)) );
				}
			tempVec.push_back( 0. ); // t = 0. for all edges

			for ( int j=3; j < components; j++ )
				tempVec.push_back( exVals[j] );

			extrema.push_back( tempVec );
			tempVec.clear();
			}
		}

	delete [] exVals;
}

