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
#include <Elements/Generic.h>



// The coordinates which vary for each face node.
// For 5th node, the 0th and 1st (r,s) vary and t is constant. 
// // The 7th entry is a dummy entry for a quadrilateral cell.
const int vtkShoeElemGeneric::faceIds[][3] =
{
	{ 5,0,1 },
	{ 4,0,1 },
	{ 3,0,2 },
	{ 2,0,2 },
	{ 1,1,2 },
	{ 0,1,2 },
	{-1,0,1 }
};

// The edge connectivity
// First vertex,second vertex,the parametric coordinate which varies, the values of constant parametric cordinates
const int vtkShoeElemGeneric::edgeConnectivity[][5]=
{
	{ 0,1,0,-1,-1 },
	{ 1,2,1, 1,-1 },
	{ 3,2,1, 1,-1 },
	{ 0,3,1,-1,-1 },
	{ 0,4,2,-1,-1 },
	{ 1,5,2, 1,-1 },
	{ 2,6,2, 1, 1 },
	{ 3,7,2,-1, 1 },
	{ 4,5,0,-1, 1 },
	{ 5,6,1, 1, 1 },
	{ 7,6,0, 1, 1 },
	{ 4,7,0,-1, 1 }
};

const int vtkShoeElemGeneric::faceConnectivity[][8] =
{
	{ 0, 1, 2, 3, 0, 1, 2, 3 },
	{ 4, 5, 6, 7, 8, 9,10,11 },
	{ 0, 1, 5, 4, 0, 5, 8, 4 },
	{ 3, 2, 6, 7, 2, 6,10, 7 },
	{ 0, 3, 7, 4, 3, 7,11, 4 },
	{ 1, 2, 6, 5, 1, 6, 9, 5 }
};

const int vtkShoeElemGeneric::faceConstantPara[][2]=
{
	{ 2,-1 },
	{ 2, 1 },
	{ 1,-1 },
	{ 1, 1 },
	{ 0,-1 },
	{ 0, 1 }
};

static const double GaussianIntegrationAbcissas_n1[] = {  0.000000000000000 };
static const double GaussianIntegrationAbcissas_n2[] = { -0.577350269189626,  0.577350269189626 };
static const double GaussianIntegrationAbcissas_n3[] = { -0.774596669241483,  0.000000000000000,  0.774596669241483 };
static const double GaussianIntegrationAbcissas_n4[] = { -0.861136311594053, -0.339981043584856,  0.339981043584856,  0.861136311594053 };
static const double GaussianIntegrationAbcissas_n5[] = { -0.906179845938664, -0.538469310105683,  0.000000000000000,  0.538469310105683,  0.906179845938664 };
static const double GaussianIntegrationAbcissas_n6[] = { -0.932469514203152, -0.661209386466265, -0.238619186083197,  0.238619186083197,  0.661209386466265,  0.932469514203152 };
static const double GaussianIntegrationAbcissas_n7[] = { -0.949107912342759, -0.741531185599394, -0.405845151377397,  0.000000000000000,  0.405845151377397,  0.741531185599394, 0.949107912342759 };
static const double GaussianIntegrationAbcissas_n8[] = { -0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650,  0.183434642495650,  0.525532409916329,  0.796666477413627, 0.960289856497536 };

const double* vtkShoeElemGeneric::GaussianIntegrationAbcissas[] =
{
	0,
	GaussianIntegrationAbcissas_n1,
	GaussianIntegrationAbcissas_n2,
	GaussianIntegrationAbcissas_n3,
	GaussianIntegrationAbcissas_n4,
	GaussianIntegrationAbcissas_n5,
	GaussianIntegrationAbcissas_n6,
	GaussianIntegrationAbcissas_n7,
	GaussianIntegrationAbcissas_n8
};

static const double GaussianIntegrationWeights_n1[] = { 2.000000000000000 };
static const double GaussianIntegrationWeights_n2[] = { 1.000000000000000, 1.000000000000000 };
static const double GaussianIntegrationWeights_n3[] = { 0.555555555555556, 0.888888888888888, 0.555555555555556 };
static const double GaussianIntegrationWeights_n4[] = { 0.347854845137454, 0.652145154862456, 0.652145154862456, 0.347854845137454 };
static const double GaussianIntegrationWeights_n5[] = { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 };
static const double GaussianIntegrationWeights_n6[] = { 0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170 };
static const double GaussianIntegrationWeights_n7[] = { 0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870 };
static const double GaussianIntegrationWeights_n8[] = { 0.101228536290376, 0.222381034453374, 0.313706645877887, 0.362683783378362, 0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290376 };

const double* vtkShoeElemGeneric::GaussianIntegrationWeights[] =
{
	0,
	GaussianIntegrationWeights_n1,
	GaussianIntegrationWeights_n2,
	GaussianIntegrationWeights_n3,
	GaussianIntegrationWeights_n4,
	GaussianIntegrationWeights_n5,
	GaussianIntegrationWeights_n6,
	GaussianIntegrationWeights_n7,
	GaussianIntegrationWeights_n8
};

// Points have no modes... they only interpolate. This same function is used to
// return the number of face/volume modes on a 1D element and volume modes of 2D elements.
int vtkShoeElemGeneric::ZeroModesPerNode( int, const int[3] )
{
	return 0;
}

int vtkShoeElemGeneric::ZeroModes( const int[3] )
{
	return 0;
}

int vtkShoeElemGeneric::OneModePerNode( int, const int[3] )
{
	return 1;
}

int vtkShoeElemGeneric::OneMode( const int[3] )
{
	return 1;
}

// For Legendre polynomials in the product space \f${\cal S}^p(\Omega^{(h)}_{st})\f$,
// all edges should have the same order, so use order[0] to compute the number of nodes.
// This works for pretty much all of the element shapes since the edge modes are standard.
int vtkShoeElemGeneric::GetNumberOfEdgeModesPerNodePLeg( int ID, const int order[3] )
{
	return order[0] - 1;
}

double* vtkShoeElemGeneric::DepermuteQuadParameters( double* pcoord, int node_permutation )
{
	double temp;
	switch( node_permutation )
	  {
	case 0:// bottom left, right up: no transformations needed
		break;
	case 1: // bottom left, up right: swap
		temp      = pcoord[0];
		pcoord[0] = pcoord[1];
		pcoord[1] = temp;
		break;
	case 2: // top left, down right: swap and negate r
		temp = pcoord[0];
		pcoord[0] = -pcoord[1];
		pcoord[1] = temp;
		break;
	case 3: // top left, right down: negate s
		pcoord[1] = pcoord[1];
		break;
	case 4: // top right, left down: negate both
		pcoord[0] = -pcoord[0];
		pcoord[1] = -pcoord[1];
		break;
	case 5: // top right, down left: swap and negate both
		temp      =  pcoord[0];
		pcoord[0] = -pcoord[1];
		pcoord[1] = -temp;
		break;
	case 6: // bottom right, up left: swap and negate s
		temp      = pcoord[0];
		pcoord[0] = pcoord[1];
		pcoord[1] = -temp;
		break;
	case 7: // bottom right, left up: negate r
		pcoord[0] = -pcoord[0];
		break ;
		}
  return 0;
}

