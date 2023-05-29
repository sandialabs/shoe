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
#ifndef Element_Quadrilateral_h
#define Element_Quadrilateral_h

#include <vtkCellOps.h>

class vtkShoeElemQuadrilateral
{
	public:
		static void GetPermutedEdgeSignsTensorLagrange( vtkstd::vector <bool>& signs, int edge_ID, bool node_permutation, const int order[3] );
		static void GetPermutedEdgeSignsLegendre( vtkstd::vector <bool>& signs, int edge_ID, bool node_permutation, const int order[3] );

		static void GetPermutedFaceIndicesTensorLagrange( vtkstd::vector <int>& indices, vtkstd::vector <bool>& signs, int face_ID, int node_permutation, const int order[3] );
		static void GetPermutedFaceIndicesPLegendre( vtkstd::vector <int>& indices, vtkstd::vector <bool>& signs, int face_ID, int node_permutation, const int order[3] );

		static int  GetNumberOfShapeFunctionsTensorLagrange(const int order[3]);
		static int  GetNumberOfShapeFunctionsMaxTotalOrderLegendre(const int order[3]); // p space

		static void EvaluateShapeFunctionsTensorLagrange(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3]);
		static void EvaluateShapeFunctionDerivativesTensorLagrange(double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3]);

		static void EvaluateShapeFunctionsTensorLegendre(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]);  // pq space
		// FIXME: No derivatives for legendre tensor product

		static void EvaluateShapeFunctionsTruncatedOrderLegendre(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]); // ~pq space
		// FIXME: No derivatives for legendre truncated order product

		static void EvaluateShapeFunctionsMaxTotalOrderLegendre(double *shape,  vtkShoeMeshIterator& cell, const int order[3], const double r[3]);
		static void EvaluateShapeFunctionDerivativesMaxTotalOrderLegendre(double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3]);

    static void EvaluateNormalOnFace( double *norm, vtkShoeMeshIterator& cell, const double r[3], const int FaceId );

		static void GetPolynomialTensorLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector< vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num);
		static void GetPolynomialMaxTotalOrderLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num);
		static void GetPolynomialTruncatedTotalOrderLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, int *orders, vtkstd::vector<double>& dof);

		static void GetCriticalPoints( vtkShoeMeshIterator& iter, int field_no );
		static bool DoCriticalPointsExist( vtkShoeMeshIterator &iter, const vtkstd::vector<double>& coef, const vtkstd::vector<vtkstd::vector<int> >& powers );
		static void GetPermutedCriticalPoints( vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkShoeMeshIterator& iter_in, const int field_number_in );

		static void GetBoundaryEdge( vtkCellDefinition& def_out, vtkIdType* connectivity_out, uint32_t& perm_out, vtkShoeMeshIterator& cell, int edge_number );
};

#endif // Element_Quadrilateral_h
