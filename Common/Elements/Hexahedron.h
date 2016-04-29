/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2005-06-06 18:49:35 -0700 (Mon, 06 Jun 2005) $
  Version:   $Revision: 4814 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef Element_Hexahedron_h
#define Element_Hexahedron_h

#include <vtkCellOps.h>

class vtkShoeElemHexahedron
{
	public:
		static int GetNumberOfEdgeModesPerNodePPQLeg( int ID, const int order[3] );

		static int GetNumberOfFaceModesPerNodePLeg( int ID, const int order[3] );
		static int GetNumberOfFaceModesPerNodePPQLeg( int ID, const int order[3] );

		static int GetNumberOfVolumeModesPerNodePLeg( const int order[3] );
		static int GetNumberOfVolumeModesPerNodePPQLeg( const int order[3] );

		static void GetBoundaryEdge( vtkCellDefinition& def_out, vtkIdType* connectivity_out, uint32_t& perm_out, vtkShoeMeshIterator& cell, int edge_number );
		static void GetBoundaryFace( vtkCellDefinition& def_out, vtkIdType* connectivity_out, uint32_t& perm_out, vtkShoeMeshIterator& cell, int face_number );

		static int FaceNodePermutation( int face_A, int face_B, int offset );
		static int ConcatenateFacePermutations( const int perm_A_in, const int perm_B_in );
		static void GetPermutedEdgeSignsLeg( vtkstd::vector <bool>& signs, int edge_ID, bool node_permutation, const int order[3] );
		static void GetPermutedFaceIndicesPLeg( vtkstd::vector <int>& indices, vtkstd::vector <bool>& signs, int face_ID, int node_permutation, const int order[3] );

		static int GetNumberOfShapeFunctionsTensorProduct(const int orders[3]); // ppq space
		static int GetNumberOfShapeFunctionsMaxTotalOrder(const int order[3]);  // p space
		static int GetNumberOfShapeFunctionsTruncatedOrder(const int order[3]);  // ~ space  // orders = {p,q}

		static void EvaluateNormalOnFace( double *norm, vtkShoeMeshIterator& cell, const double r[3], const int FaceId );

		static void GetPermutedCriticalPoints( vtkstd::vector<vtkstd::vector<double> >& extrema, vtkShoeMeshIterator& iter, const int field_no );
		static void GetCriticalPoints( vtkShoeMeshIterator& iter, int field_no );

		static void GetIsosurface( vtkShoeMeshIterator& iter, vtkAdaptiveTessellator* tess_in, const double Isovalue, const int field_no );

		static void EvaluateShapeFunctionsMaxTotalOrder(double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3]);
		static void EvaluateShapeFunctionsTensorProduct(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]);  // ppq space
		static void EvaluateShapeFunctionsTruncatedOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]); // ~pq space

		static void EvaluateShapeFunctionDerivativesMaxTotalOrder(double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3]);
		static void EvaluateShapeFunctionDerivativesTensorProduct(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]);  // ppq space
		static void EvaluateTruncatedOrderShapeFunctionsDerivatives(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]); // ~pq space

		static void GetFieldPolynomialMaxTotalOrder(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num);
		static void GetFieldPolynomialTruncatedOrder(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num);
		static void GetFieldPolynomialTensorProduct(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num);

		static bool DoCriticalPointsExist(vtkShoeMeshIterator &iter, const vtkstd::vector<double>& coef, const vtkstd::vector<vtkstd::vector<int> >& powers);
};

#endif // Element_Hexahedron_h
