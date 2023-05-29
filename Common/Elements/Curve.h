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
#ifndef Element_Curve_h
#define Element_Curve_h

#include <vtkCellOps.h>

class vtkShoeElemCurve
{
	public:
		static int  GetNumberOfEdgeModesPerNode( int id, const int o[3] );
		static int  GetNumberOfShapeFunctions(const int order[3]);
		static void GetPermutedEdgeSigns( vtkstd::vector <bool>& signs, int, bool node_permutation, const int order[3] );

		static void EvaluateShapeFunctionsMaxTotalOrderLegendre( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] );
		static void EvaluateShapeFunctionDerivativesMaxTotalOrderLegendre( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] );

		static void EvaluateShapeFunctionsTensorLagrange( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] );
		static void EvaluateShapeFunctionDerivativesTensorLagrange( double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3] );

		static void GetPolynomialMaxTotalOrderLegendre(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector< vtkstd::vector< vtkstd::vector<int> > >& powers,  vtkShoeMeshIterator& cell_in, const int field_num);
};

#endif // Element_Curve_h
