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
#include <Elements/CriticalPoints.h>
#include <Elements/ShapeFunctions.h>
#include <Elements/Tetrahedron.h>

#ifdef FOUND_GINAC
USING_NAMESPACE(GiNaC);
#endif // FOUND_GINAC

// Return total number of shape functions for a tetrahedron
int vtkShoeElemTetrahedron::GetNumberOfShapeFunctionsPLeg(const int order[3])
{
  int p = order[0];
  return 4 + 6*(p-1) + 2*(p-1)*(p-2) + (p-1)*(p-2)*(p-3)/6;
}

int vtkShoeElemTetrahedron::GetNumberOfFaceModesPerNodePLeg( int ID, const int order[3] )
{
  return ((order[0] - 1) * (order[0] - 2) / 2);
}

int vtkShoeElemTetrahedron::GetNumberOfVolumeModesPerNodePLeg( const int order[3] )
{
  return ((order[0] - 1) * (order[0] - 2) * (order[0] - 3) / 6) ;
}

// return the permutation code for face_B of a tetrahedron based on its
// mating face's code and the offset of their starting vertices
int vtkShoeElemTetrahedron::FaceNodePermutation( int face_A, int face_B, int offset )
{
  int perm;
  switch (offset)
  {
    case 0:
      if (face_A == 3)
          perm = (face_B == 3) ? 2 : 0;
      else
        perm = (face_B == 3) ? 0 : 2;
      break;

    case 1:
      if (face_A == 3)
          perm = (face_B == 3) ? 4 : 6;
      else
        perm = (face_B == 3) ? 5 : 1;
      break;

    case 2:
      if (face_A == 3)
          perm = (face_B == 3) ? 1 : 5;
      else
        perm = (face_B == 3) ? 6 : 4;
      break;
  }

  return perm;
}

int vtkShoeElemTetrahedron::ConcatenateFacePermutations( int perm_A, int perm_B )
{
  int perm_codes[][3] =
  {
    {0, 1, 2},
    {1, 0, 2},
    {2, 1, 0}
  };

  int bit_A = perm_A & 3;
  int bit_B = perm_B & 3;
  int perm = perm_codes[bit_A][bit_B];

  bit_A = perm_A & 4;
  bit_B = perm_B & 4;
  perm |= bit_A ^ bit_B;

  if (((perm_A & 3) + (perm_B & 3)) == 3)
    if (perm & 4)
      perm &= 3;
    else
      perm |= 4;

  return perm;
}

// Permute (negate) Legendre polynomial amplitudes as necessary for tetrahedral edge with given ID.
// Assumes max total order product space, and order[0] contains the correct order.
void vtkShoeElemTetrahedron::GetPermutedEdgeSignsPLeg( vtkstd::vector<bool>& signs, int edge_ID, bool node_permutation, const int order[3] )
{
  if ( (int)signs.size() != order[0] - 1)
    signs.resize( order[0] - 1 );

  vtkstd::fill( signs.begin(), signs.end(), false );
  
  if (node_permutation)
    for (int i = 3; i <= order[0]; i += 2)
      signs[i-2] = true;
}

// Based on the permutation bits for the face, determine the indices and signs
// of tetrahedral face modes after amplitude swaps and sign changes.
void vtkShoeElemTetrahedron::GetPermutedFaceIndicesPLeg( vtkstd::vector<int>& indices, vtkstd::vector<bool>& signs, int face_ID, int node_permutation, const int order[3] )
{
  // compute number of face modes for this node
  // for this product space, assume orders are equal for all parameters
  // if (n_indices < 1) should declare an error
  int n_indices = (order[0] - 1) * (order[0] - 2) / 2;

  // make sure the vectors are big enough
  if ( (int) indices.size() != n_indices)
    indices.resize( n_indices );
  if ( (int) signs.size() != n_indices)
    signs.resize( n_indices );

  // initialize the signs vector -- assume no signs change
  vtkstd::fill( signs.begin(), signs.end(), false ) ;

  // determine if we need to negate odd powers
  int negate = node_permutation & 4;

  // counter for shape functions
  int m = 0;
  for (int p = 0; p <= order[0]-3; p++)
  {
    // This loop implements the indexing in Szabo and Babuska, with a little help from
    // intermediate variables k and p_local
    for (int j = 0, i = p; i >= 0; j++, i--)
    {
      // load them and set sign flags for odd values of i
      indices[m] = m;
      if (negate && (i % 2))
        signs[m] = true;
      m++;
    }
  }
}

// Evaluate the normal on a face of a tetrahedron. The normal is
// calculated by taking cross product of parametric
// tangents along two orthogonal vectors lying in that face and which
// span that face completely. These vectors vary for each face of the tetrahedron.
// Parametric tangent along a certain direction is given by the gradient dotted with that direction.
// This function produces the outward pointing normal of the face.
void vtkShoeElemTetrahedron::EvaluateNormalOnFace(double *norm,vtkShoeMeshIterator& cell, const double r[3],const int FaceId)
{
  // Mutually perpendicular directions for each face of the tetrahedron in parameter space.
  static const double direction[8][3] =
  { 
    {-0.5,      0.866,   0.      },
    {-0.28868, -0.16667, 0.94281 },
    {-0.5,     -0.86603, 0.      },
    { 0.28868, -0.1667,  0.94281 },
    { 1.,       0.,      0.      },
    { 0.,       0.333,   0.94281 },
    { 1.,       0.,      0.      },
    { 0.,      -1.,      0.      }
  };

  int dimension = cell.GetCellOps()->EmbeddingDimension;
  const int *order = cell.GetCellGeometricOrder();
  int no=cell.GetCellOps()->GetNumberOfShapeFunctions(order);
  vtkDataArray *geometry = cell.GetCachedGeometry();
  double *shape = new double[no*dimension];
  cell.GetCellOps()->EvaluateShapeFuncDerivatives(shape,cell,order,r);
  double n1[3] = {0.,0.,0.};
  double n2[3] = {0.,0.,0.};
  double n3[3];
  double geom[3];

  //     For every direction
  for(int j=0;j<dimension;j++)
  {
    // For every shape function
    for(int cnt=0;cnt<no;cnt++)
    {
      geometry->GetTuple(cnt,geom);
      //           For each derivative direction
      for(int i=0;i<dimension;i++)
      {
        n1[i] += geom[j]*shape[i*no+cnt];
      }
    }
    // n1 is now the gradient of the geometry function for j^th coordinate.
    // Dot it with the direction in which the derivative is required.

    n2[j] = vtkMath::Dot(n1,direction[2*FaceId]);
    n3[j] = vtkMath::Dot(n1,direction[2*FaceId+1]);


  }

  vtkMath::Cross(n2,n3,norm);
  vtkMath::Normalize(norm);
  delete [] shape;
}

static int triFaceIds[4][4] =
{
  {1,2,3},
  {2,0,3},
  {0,1,3},
  {0,1,2}
};

// Evaluate the shape functions for a tetrahedral element.
void vtkShoeElemTetrahedron::EvaluateShapeFunctionsMaxTotalOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3])
{
  double l[4];
  l[0] = .5*(1-r[0]-.57735027*r[1]-.408248290464*r[2]);
  l[1] = .5*(1+r[0]-.57735027*r[1]-.408248290464*r[2]);
  l[2] = .57735027*(r[1]-.35355339059327*r[2]);
  l[3] = .61237243569579*r[2];
  int index=0;
  double x[6];
  double y[6];
  double z[6];
  //vertex shape functions
  shape[index++] = l[0];
  shape[index++] = l[1];
  shape[index++] = l[2];
  shape[index++] = l[3];
  //Edge Modes 
  //Edge 0-1
  x[0]=l[1]-l[0];
  for(int cnt1=1;cnt1<order[0];cnt1++)
  {
    x[cnt1] = x[cnt1-1]*(l[1]-l[0]);
    
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++] = l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
  }
  // Edge 1-2
  x[0]=l[2]-l[1];
  for(int cnt1=1;cnt1<order[0];cnt1++)
  {
    x[cnt1] = x[cnt1-1]*(l[2]-l[1]);
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++] = l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    
  }
  // Edge 2-0
  x[0]=l[0]-l[2];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[0]-l[2]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++] = l[2]*l[0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
  }
  // edge 0-3
  x[0]=l[3]-l[0];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[3]-l[0]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++] = l[0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
  }
  //edge 1-3
  x[0]=l[3]-l[1];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[3]-l[1]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++] = l[1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
  }
  //edge 2-3
  x[0]=l[3]-l[2];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[3]-l[2]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++] = l[2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
  }

  // Face Modes - possibly permuted
  int shift[][3] = { {0, 1, 2}, {2, 1, 0}, {0, 2, 1} };
  for(int Fcnt=0;Fcnt<4;Fcnt++)
  {
    int perm = cell.GetCell()->GetFacePermutation(Fcnt) & 3;

    x[0]=l[triFaceIds[Fcnt][shift[perm][1]]]-l[triFaceIds[Fcnt][shift[perm][0]]];
    for(int xcnt=1;xcnt<order[0];xcnt++)
    {
      x[xcnt] = x[xcnt-1]*(l[triFaceIds[Fcnt][shift[perm][1]]]-l[triFaceIds[Fcnt][shift[perm][0]]]);
    }
    y[0]=2*l[triFaceIds[Fcnt][shift[perm][2]]]-1;
    for(int xcnt=1;xcnt<order[0];xcnt++)
    {
      y[xcnt] = y[xcnt-1]*(2*l[triFaceIds[Fcnt][shift[perm][2]]]-1);
    }
    for(int cnt=0;cnt<=order[0]-3;cnt++)
    {
      for(int cnt1=cnt;cnt1>=0;cnt1--)
      {
        shape[index++] = l[triFaceIds[Fcnt][shift[perm][0]]]*l[triFaceIds[Fcnt][shift[perm][1]]]*l[triFaceIds[Fcnt][shift[perm][2]]]
                          *vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y);
      }
    }
  }

  // Volume modes
  x[0]=l[1]-l[0];
  y[0]=2*l[2]-1;
  z[0]=2*l[3]-1;

  for(int cnt=1;cnt<order[0];cnt++)
  {
    x[cnt] = x[cnt-1]*x[0];
    y[cnt] = y[cnt-1]*y[0];
    z[cnt] = z[cnt-1]*z[0];
  }
  double factor = l[0]*l[1]*l[2]*l[3];
  for(int cnt=0;cnt<=order[0]-4;cnt++)
  {
    for(int cnt1=cnt;cnt1>=0;cnt1--)
    {
      for(int cnt2=cnt-cnt1;cnt2>=0;cnt2--)
      {
        shape[index++] = factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z);
      }
    }
  }
}

void vtkShoeElemTetrahedron::EvaluateShapeFunctionDerivativesMaxTotalOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3])
{
  double l[4];
  double ld[4][3];
  
  l[0] = .5*(1-r[0]-.57735027*r[1]-.408248290464*r[2]);
  l[1] = .5*(1+r[0]-.57735027*r[1]-.408248290464*r[2]);
  l[2] = .57735027*(r[1]-.35355339059327*r[2]);
  l[3] = .61237243569579*r[2];

  ld[0][0] = -.5; 
  ld[0][1] = -1/sqrt(3.)*.5; 
  ld[0][2] = -.5*1/sqrt(6.);
  ld[1][0] = .5;
  ld[1][1] = -1/sqrt(3.)*.5;
  ld[1][2] = -1/sqrt(6.)*.5;
  ld[2][0] = 0.;
  ld[2][1] = 1/sqrt(3.);
  ld[2][2] = -1/sqrt(3.)*1/sqrt(8.);
  ld[3][0] = 0.;
  ld[3][1] = 0.;
  ld[3][2] = sqrt(3./8.);
  int no = cell.GetCellOps()->GetNumberOfShapeFunctions(order);

  
  int index=0;
  double x[6];
  double y[6];
  double z[6];
  //vertex shape functions
  shape[index] = ld[0][0];
  shape[no+index]= ld[0][1];
  shape[2*no+index++] = ld[0][2];
  
  // shape[index++] = l[1];
  shape[index] = ld[1][0];
  shape[no+index] = ld[1][1];
  shape[2*no+index++] = ld[1][2];
  
  // shape[index++] = l[2];
  shape[index] = ld[2][0];
  shape[no+index] = ld[2][1];
  shape[2*no+index++] = ld[2][2];

    
  // shape[index++] = l[3];
  shape[index] = ld[3][0];
  shape[no+index] = ld[3][1];
  shape[2*no+index++] = ld[3][2];
  
  
  //Edge Modes 
  //Edge 0-1
  x[0]=l[1]-l[0];
  for(int cnt1=1;cnt1<order[0];cnt1++)
  {
    x[cnt1] = x[cnt1-1]*(l[1]-l[0]);
    
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
//    shape[index++] = l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    shape[index] =
      ld[0][0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[1][0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[1][0]-ld[0][0]);
    shape[no+index] =
      ld[0][1]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[1][1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[1][1]-ld[0][1]);
    shape[2*no+index++] =
      ld[0][2]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[1][2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[1]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[1][2]-ld[0][2]);
  }

  // Edge 1-2
  x[0]=l[2]-l[1];
  for(int cnt1=1;cnt1<order[0];cnt1++)
  {
    x[cnt1] = x[cnt1-1]*(l[2]-l[1]);
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
//    shape[index++] = l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    shape[index] =
      ld[1][0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*ld[2][0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[2][0]-ld[1][0]);
    shape[no+index] =
      ld[1][1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*ld[2][1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[2][1]-ld[1][1]);
    shape[2*no+index++] =
      ld[1][2]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*ld[2][2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[2][2]-ld[1][2]);
  }
  
  // Edge 2-0
  x[0]=l[0]-l[2];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[0]-l[2]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
//    shape[index++] = l[2]*l[0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    shape[index] =
      ld[0][0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[2][0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[0][0]-ld[2][0]);
    shape[no+index] =
      ld[0][1]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[2][1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[0][1]-ld[2][1]);
    shape[2*no+index++] =
      ld[0][2]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[2][2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[2]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[0][2]-ld[2][2]);
  }
  // edge 0-3
  x[0]=l[3]-l[0];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[3]-l[0]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
//    shape[index++] = l[3]*l[0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    shape[index] =
      ld[0][0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[3][0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][0]-ld[0][0]);
    shape[no+index] =
      ld[0][1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[3][1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][1]-ld[0][1]);
    shape[2*no+index++] =
      ld[0][2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*ld[3][2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][2]-ld[0][2]);
    
  }
  
  //edge 1-3
  x[0]=l[3]-l[1];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[3]-l[1]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
//    shape[index++] = l[1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    shape[index] =
      ld[1][0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*ld[3][0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][0]-ld[1][0]);
    shape[no+index] =
      ld[1][1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*ld[3][1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][1]-ld[1][1]);
    shape[2*no+index++] =
      ld[1][2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*ld[3][2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][2]-ld[1][2]);
  }
  
  //edge 2-3
  x[0]=l[3]-l[2];
  for(int cnt1=1;cnt1<order[0];cnt1++)
    x[cnt1] = x[cnt1-1]*(l[3]-l[2]);

  for(int cnt=2;cnt<=order[0];cnt++)
  {
//    shape[index++] = l[2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x);
    shape[index] =
      ld[2][0]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[2]*ld[3][0]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][0]-ld[2][0]);
    shape[no+index] =
      ld[2][1]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[2]*ld[3][1]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][1]-ld[2][1]);
    shape[2*no+index++] =
      ld[2][2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[2]*ld[3][2]*vtkShoeElemShapeFunctions::EvaluatePsi[cnt](x) +
      l[2]*l[3]*vtkShoeElemShapeFunctions::EvaluatePsiDerivative[cnt](x)*(ld[3][2]-ld[2][2]);
  }

  // Face Modes - possibly permuted
  int shift[][3] = { {0, 1, 2}, {2, 1, 0}, {0, 2, 1} };
  for(int Fcnt=0;Fcnt<4;Fcnt++)
  {
    int perm = cell.GetCell()->GetFacePermutation(Fcnt) & 3;

    x[0]=l[triFaceIds[Fcnt][shift[perm][1]]]-l[triFaceIds[Fcnt][shift[perm][0]]];
    for(int xcnt=1;xcnt<order[0];xcnt++)
    {
      x[xcnt] = x[xcnt-1]*(l[triFaceIds[Fcnt][shift[perm][1]]]-l[triFaceIds[Fcnt][shift[perm][0]]]);
    }
    y[0]=2*l[triFaceIds[Fcnt][shift[perm][2]]]-1;
    for(int xcnt=1;xcnt<order[0];xcnt++)
    {
      y[xcnt] = y[xcnt-1]*(2*l[triFaceIds[Fcnt][shift[perm][2]]]-1);
    }
    for(int cnt=0;cnt<=order[0]-3;cnt++)
    {
      for(int cnt1=cnt;cnt1>=0;cnt1--)
      {
        //        shape[index++] = l[triFaceIds[Fcnt][0]]*l[triFaceIds[Fcnt][1]]*l[triFaceIds[Fcnt][2]]*vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y);
        shape[index] =
          ld[triFaceIds[Fcnt][shift[perm][0]]][0]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) + 

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          ld[triFaceIds[Fcnt][shift[perm][1]]][0]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          ld[triFaceIds[Fcnt][shift[perm][2]]][0]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y)*
          (ld[triFaceIds[Fcnt][shift[perm][1]]][0]-ld[triFaceIds[Fcnt][shift[perm][0]]][0]) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1](y)*
          (2*ld[triFaceIds[Fcnt][shift[perm][2]]][0]);

        shape[no+index] =
          ld[triFaceIds[Fcnt][shift[perm][0]]][1]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          ld[triFaceIds[Fcnt][shift[perm][1]]][1]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          ld[triFaceIds[Fcnt][shift[perm][2]]][1]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y)*
          (ld[triFaceIds[Fcnt][shift[perm][1]]][1]-ld[triFaceIds[Fcnt][shift[perm][0]]][1]) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1](y)*
          (2*ld[triFaceIds[Fcnt][shift[perm][2]]][1]);

        shape[2*no+index++]=
          ld[triFaceIds[Fcnt][shift[perm][0]]][2]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) + 

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          ld[triFaceIds[Fcnt][shift[perm][1]]][2]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          ld[triFaceIds[Fcnt][shift[perm][2]]][2]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1](y)*
          (ld[triFaceIds[Fcnt][shift[perm][1]]][2]-ld[triFaceIds[Fcnt][shift[perm][0]]][2]) +

          l[triFaceIds[Fcnt][shift[perm][0]]]*
          l[triFaceIds[Fcnt][shift[perm][1]]]*
          l[triFaceIds[Fcnt][shift[perm][2]]]*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1](y)*
          (2*ld[triFaceIds[Fcnt][shift[perm][2]]][2]);
      }
    }
  }

  // Volume modes
  x[0]=l[1]-l[0];
  y[0]=2*l[2]-1;
  z[0]=2*l[3]-1;

  for(int cnt=1;cnt<order[0];cnt++)
  {
    x[cnt] = x[cnt-1]*x[0];
    y[cnt] = y[cnt-1]*y[0];
    z[cnt] = z[cnt-1]*z[0];
  }
  double factor = l[0]*l[1]*l[2]*l[3];
  double factord1 = ld[0][0]*l[1]*l[2]*l[3] + l[0]*ld[1][0]*l[2]*l[3]+ l[0]*l[1]*ld[2][0]*l[3]+l[0]*l[1]*l[2]*ld[3][0];
  double factord2 =  ld[0][1]*l[1]*l[2]*l[3] + l[0]*ld[1][1]*l[2]*l[3]+ l[0]*l[1]*ld[2][1]*l[3]+l[0]*l[1]*l[2]*ld[3][1];
  double factord3 =  ld[0][2]*l[1]*l[2]*l[3] + l[0]*ld[1][2]*l[2]*l[3]+ l[0]*l[1]*ld[2][2]*l[3]+l[0]*l[1]*l[2]*ld[3][2];

  for(int cnt=0;cnt<=order[0]-4;cnt++)
  {
    for(int cnt1=cnt;cnt1>=0;cnt1--)
    {
      for(int cnt2=cnt-cnt1;cnt2>=0;cnt2--)
      {
        shape[index] =
          factord1*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z)*
          (ld[1][0]-ld[0][0]) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z)*
          (2*ld[2][0]) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1-cnt2](z)*
          (2* ld[3][0]);
        shape[no+index] =
          factord2*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z)*
          (ld[1][1]-ld[0][1]) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z)*
          (2*ld[2][1]) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1-cnt2](z)*
          (2*ld[3][1]);
        shape[ 2*no+index++ ] =
          factord3*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z)*
          (ld[1][2]-ld[0][2]) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt-cnt1-cnt2](z)*
          (2*ld[2][2]) +

          factor*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt1](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendre[cnt2](x)*
          vtkShoeElemShapeFunctions::EvaluateLegendreDerivative[cnt-cnt1-cnt2](z)*
          (2*ld[3][2]);
      }
    }
  }
}

void vtkShoeElemTetrahedron::GetPolynomialMaxTotalOrder(vtkstd::vector<double>& coeff, vtkstd::vector<vtkstd::vector<int> >&powers, int *orders,vtkstd::vector<double>& dof)
{
#ifdef FOUND_GINAC
  int TrifaceIds[4][4] = { {1,2,3},
    {0,3,2},
    {0,1,3},
    {0,1,2} };


    symbol x("x");
    symbol y("y");
    symbol z("z");

    ex l[4];
    l[0] = .5*(1-x-.57735027*y-.408248290464*z);
    l[1] = .5*(1+x-.57735027*y-.408248290464*z);
    l[2] = .57735027*(y-.35355339059327*z);
    l[3] = .61237243569579*z;
    int index=0;
    //  double x[6];
    //  double y[6];
    //  double z[6];
    ex finalPoly;
    //vertex shape functions
    finalPoly = dof[index++]*l[0]; 
    finalPoly += dof[index++]*l[1]; 
    finalPoly += dof[index++]*l[2] ;
    finalPoly += dof[index++]*l[3];
    //Edge Modes 
    //Edge 1-2
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += dof[index++]*l[0]*l[1]*vtkShoeElemCriticalPoints::eval_psi(cnt,l[1]-l[0]);
    }
    // Edge 2-3
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
    finalPoly += dof[index++]*l[1]*l[2]*vtkShoeElemCriticalPoints::eval_psi(cnt,l[2]-l[1]);
    
  }
  // Edge 3-1
  for(int cnt=2;cnt<=orders[0];cnt++)
  {
    finalPoly += dof[index++]*l[2]*l[0]*vtkShoeElemCriticalPoints::eval_psi(cnt,l[2]-l[0]);
  }
  // edge 1-4
  for(int cnt=2;cnt<=orders[0];cnt++)
  {
    finalPoly += dof[index++]*l[3]*l[0]*vtkShoeElemCriticalPoints::eval_psi(cnt,l[3]-l[0]);
  }
  //edge 2-4

  for(int cnt=2;cnt<=orders[0];cnt++)
  {
    finalPoly += dof[index++]*l[1]*l[3]*vtkShoeElemCriticalPoints::eval_psi(cnt,l[3]-l[1]);
  }
  //edge 3-4
  for(int cnt=2;cnt<=orders[0];cnt++)
  {
        
    finalPoly += dof[index++]*l[2]*l[3]*vtkShoeElemCriticalPoints::eval_psi(cnt,l[2]-l[3]);
  }

  // Face Modes - 
  ex xx,yy,zz;
  for(int Fcnt=0;Fcnt<4;Fcnt++)
  {
    xx = l[TrifaceIds[Fcnt][1]]-l[TrifaceIds[Fcnt][0]];
    yy = 2*l[TrifaceIds[Fcnt][2]]-1;
    for(int cnt=0;cnt<=orders[0]-3;cnt++)
    {
      for(int cnt1=cnt;cnt1>=0;cnt1--)
      {
        finalPoly +=
          dof[index++]*l[TrifaceIds[Fcnt][0]]*l[TrifaceIds[Fcnt][1]]*l[TrifaceIds[Fcnt][2]]*
          vtkShoeElemCriticalPoints::legendre_poly(cnt1,xx)*
          vtkShoeElemCriticalPoints::legendre_poly(cnt-cnt1,yy);
      }
    }
  }
  // Volume modes
  xx=l[1]-l[0];
  yy=2*l[2]-1;
  zz=2*l[3]-1;
  ex factor = l[0]*l[1]*l[2]*l[3];
  for(int cnt=0;cnt<=orders[0]-4;cnt++)
  {
    for(int cnt1=cnt;cnt1>=0;cnt1--)
    {
      for(int cnt2=cnt-cnt1;cnt2>=0;cnt2--)
      {
        finalPoly +=
          dof[index++]*factor*
          vtkShoeElemCriticalPoints::legendre_poly(cnt1,xx)*
          vtkShoeElemCriticalPoints::legendre_poly(cnt2,x)*
          vtkShoeElemCriticalPoints::legendre_poly(cnt-cnt1-cnt2,z);

      }
    }
  }
  finalPoly = finalPoly.expand().evalf();
  // Generate the powers and coefficient vectors
  for(int cnt=0;cnt<finalPoly.nops();cnt++)
  {
    vtkstd::vector<int> tempvec;
    ex temp;  
    temp=finalPoly.op(cnt);
    tempvec.push_back(temp.degree(x));
    tempvec.push_back(temp.degree(y));
    tempvec.push_back(temp.degree(z));
    // cout<<tempvec[0]<<" "<<tempvec[1]<<" "<<tempvec[2]<<endl;
    powers.push_back(tempvec);
    numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.,y==1.,z==1.)).evalf());
    coeff.push_back(tempnum.to_double());
  } 
#endif // FOUND_GINAC
}

