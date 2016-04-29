/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-12-09 18:54:48 -0800 (Tue, 09 Dec 2003) $
  Version:   $Revision: 1119 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef Element_Tetrahedron_h
#define Element_Tetrahedron_h

#include <vtkCellOps.h>

class vtkShoeElemTetrahedron
{
	public:
		static int GetNumberOfFaceModesPerNodePLeg( int ID, const int order[3] );

		static int GetNumberOfVolumeModesPerNodePLeg( const int order[3] );

		static int GetNumberOfShapeFunctionsPLeg(const int order[3]);

		static int FaceNodePermutation( int face_A, int face_B, int offset );
		static int ConcatenateFacePermutations( int perm_A, int perm_B );
		static void GetPermutedEdgeSignsPLeg( vtkstd::vector <bool>& signs, int edge_ID, bool node_permutation, const int order[3] );
		static void GetPermutedFaceIndicesPLeg( vtkstd::vector <int>& indices, vtkstd::vector <bool>& signs, int face_ID, int node_permutation, const int order[3] );

		static void EvaluateNormalOnFace(double *norm,vtkShoeMeshIterator& cell, const double r[3],const int FaceId=6);

		static void EvaluateShapeFunctionsMaxTotalOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3]);
		static void EvaluateShapeFunctionDerivativesMaxTotalOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3]);
		static void GetPolynomialMaxTotalOrder(vtkstd::vector<double>& coeff, vtkstd::vector<vtkstd::vector<int> >&powers, int *orders,vtkstd::vector<double>& dof);
};

#endif // Element_Tetrahedron_h
