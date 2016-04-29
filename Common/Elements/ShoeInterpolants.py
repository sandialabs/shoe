#!/usr/bin/python
#
# Copyright 2012 Sandia Corporation.
# Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the
# U.S. Government. Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that this Notice and any
# statement of authorship are reproduced on all copies.
#

# ============================================================================
class vtkUnivariateLagrange:
    """
A class for constructing code that evaluates univariate Lagrange polynomial
shape functions and their derivatives
"""
    # ============
    def NumericShapeFunctions( OutputFile ):
        print >> OutputFile, """
  double v;
  for ( i=0; i<3; ++i )
    {
    v = order[i] * ( pcoords[i] + 1. ) / 2;
    for ( j=0; j<=order[i]; ++j )
      {
      l[i][j] = 1.;
      for ( k=0; k<=order[i]; ++k )
        {
        if ( j != k )
          { // FIXME: (j - k) could be a register decremented inside the k loop:
          l[i][j] *= ( v - k ) / ( j - k ); 
          }
        }
      }
    }
"""
    NumericShapeFunctions = staticmethod( NumericShapeFunctions )

    def SymbolicShapeFunctions( OutputFile ):
        print >> OutputFile, """
  vtkPolynomalExpanderOperation* V;
  for ( i=0; i<3; ++i )
    {
    V = NT( '*', NN( order[i] / 2. ), NT( '+', NS( 'r' + i, 1 ), NN( 1 ) ) );
    for ( j=0; j<=order[i]; ++j )
      {      
      l[i][j] = NN( 1. );
      for ( k=0; k<=order[i]; ++k )
        {
        if ( j != k )
          {
          l[i][j] = NT( '*', l[i][j], NT( '/', NT( '-', V, NN(k) ), NN( j - k ) ) );
          }
        }
      }
    }
"""
    SymbolicShapeFunctions = staticmethod( SymbolicShapeFunctions )

    # ============
    def NumericShapeFunctionDerivatives( OutputFile ):
        print >> OutputFile, 'std::cout << "UnivariateLagrange does not implement NumericShapeFunctionDerivatives.\\n";'
    NumericShapeFunctionDerivatives = staticmethod( NumericShapeFunctionDerivatives )

# ============================================================================
class vtkTensorProduct:
    """
A class that writes code to evaluate multivariate polynomials given  a
univariate polynomial through a product."""
    # ============
    def NumericShapeFunctions( OutputFile ):
        print >> OutputFile, """
  int sn=0;

  // Corners
  shape[sn++] = l[0][0]        * l[1][0]        * l[2][0];
  shape[sn++] = l[0][order[0]] * l[1][0]        * l[2][0];
  shape[sn++] = l[0][order[0]] * l[1][order[1]] * l[2][0];
  shape[sn++] = l[0][0]        * l[1][order[1]] * l[2][0];
  shape[sn++] = l[0][0]        * l[1][0]        * l[2][order[2]];
  shape[sn++] = l[0][order[0]] * l[1][0]        * l[2][order[2]];
  shape[sn++] = l[0][order[0]] * l[1][order[1]] * l[2][order[2]];
  shape[sn++] = l[0][0]        * l[1][order[1]] * l[2][order[2]];
  

  int sn1, sn2, sn3;
  sn1 = order[0] + order[1] - 2;
  sn2 = sn1 * 2;
  sn3 = sn + sn1 + sn2;
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[0]; ++i )
    {
    //cout << sn << ", " << sn1 << ", " << sn2 << ", " << sn3 << "\\n";
    shape[sn++ ] = l[0][i] * l[1][0] * l[2][0];               // Edge 0-1
    shape[sn1++] = l[0][i] * l[1][order[1]] * l[2][0];        // Edge 2-3
    shape[sn2++] = l[0][i] * l[1][0] * l[2][order[2]];        // Edge 4-5
    shape[sn3++] = l[0][i] * l[1][order[1]] * l[2][order[2]]; // Edge 6-7
    }
  for ( i=1; i<order[1]; ++i )
    {
    //cout << sn << ", " << sn1 << ", " << sn2 << ", " << sn3 << "\\n";
    shape[sn++ ] = l[0][order[0]] * l[1][i] * l[2][0];        // Edge 1-2
    shape[sn1++] = l[0][0] * l[1][i] * l[2][0];               // Edge 3-0
    shape[sn2++] = l[0][order[0]] * l[1][i] * l[2][order[2]]; // Edge 5-6
    shape[sn3++] = l[0][0] * l[1][i] * l[2][order[2]];        // Edge 7-4
    }
  sn = sn3;
  sn1 = order[2] - 1;
  sn2 = sn1 * 2;
  sn3 = sn + sn1 + sn2;
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[2]; ++i )
    {
    //cout << sn << ", " << sn1 << ", " << sn2 << ", " << sn3 << "\\n";
    shape[sn++ ] = l[0][0] * l[1][0] * l[2][i];                // Edge 0-4
    shape[sn1++] = l[0][order[0]] * l[1][0] * l[2][i];         // Edge 1-5
    shape[sn2++] = l[0][order[0]] * l[1][order[1]] * l[2][i];  // Edge 2-6
    shape[sn3++] = l[0][0] * l[1][order[1]] * l[2][i];         // Edge 3-7
    }

  sn = sn3;
  sn1 = (order[1] - 1)*(order[2] - 1);
  sn2 = sn1 * 2;
  sn3 = sn + sn2 + (order[2] - 1)*(order[0] - 1);
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[2]; ++i )
    {
    for ( j=1; j<order[1]; ++j )
      {
      //cout << sn << ", " << sn1 << "\\n";
      shape[sn++ ] = l[0][0] * l[1][j] * l[2][i];        // Face 0-4-7-3
      shape[sn1++] = l[0][order[0]] * l[1][j] * l[2][i]; // Face 1-2-6-5
      }
    for ( j=1; j<order[0]; ++j )
      {
      //cout << sn2 << ", " << sn3 << "\\n";
      shape[sn2++] = l[0][j] * l[1][0] * l[2][i];        // Face 0-1-5-4
      shape[sn3++] = l[0][j] * l[1][order[1]] * l[2][i]; // Face 2-3-7-6
      }
    }
  sn = sn3;
  sn1 = sn + (order[0] - 1)*(order[1] - 1);
  for ( i=1; i<order[1]; ++i )
    {
    for ( j=1; j<order[0]; ++j )
      {
      //cout << sn << ", " << sn1 << "\\n";
      shape[sn++ ] = l[0][j] * l[1][i] * l[2][0];        // Face 0-1-2-3
      shape[sn1++] = l[0][j] * l[1][i] * l[2][order[2]]; // Face 4-7-6-5
      }
    }
  sn = sn1;
  for ( i=1; i<order[2]; ++i )
    {
    for ( j=1; j<order[1]; ++j )
      {
      for ( k=1; k<order[0]; ++k )
        {
        //cout << sn << "\\n";
        shape[sn++] = l[0][k] * l[1][j] * l[2][i]; // Body
        }
      }
    }
"""
    NumericShapeFunctions = staticmethod( NumericShapeFunctions );

    # ============
    def NumericShapeFunctionDerivatives( OutputFile ):
        print >> OutputFile, \
              'std::cout << "TensorProduct does not implement NumericShapeFunctionDerivatives!\\n";'
    NumericShapeFunctionDerivatives = staticmethod( NumericShapeFunctionDerivatives );

# ============================================================================
class vtkShoeInterpolantGenerator:
    """
vtkShoeInterpolantGenerator is a class to output various SHOE interpolants
and derivatives with both FP and analytical evaluation routines
"""
    # ============
    def __init__( self, outputFileName, cellShape, interpolant, productSpace ):
        self.OutputFileName = outputFileName
        self.OutputFile = open( self.OutputFileName, 'w' )
        self.CellShape = cellShape
        self.Interpolant = interpolant
        self.ProductSpace = productSpace

    # ============
    def GenerateHeader( self ):
        print >> self.OutputFile, """/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "Elements/Shoe%s.h"
""" %( self.CellShape )

    # ============
    def GenerateInterpolants( self ):
        print >> self.OutputFile, """
void Shoe%s::%s%sShapeFunctions( const int* order, const double* pcoords, double* shape )
{
  // FIXME: Eventually needs to be varying length.
  double l[3][ShoeHexahedron::MaxDegree + 1];
  int i,j,k;
""" %( self.CellShape, self.Interpolant, self.ProductSpace )

        eval( 'vtkUnivariate%s.NumericShapeFunctions( self.OutputFile )' % self.Interpolant )
        eval( 'vtk%sProduct.NumericShapeFunctions( self.OutputFile )' % self.ProductSpace )
        print >> self.OutputFile, """
}
"""
        print >> self.OutputFile, """
void Shoe%s::Symbolic%s%sShapeFunctions( const int* order, const double* pcoords, double* shape )
{
  // FIXME: Eventually needs to be varying length.
  // vtkPolynomialExpanderOperation* l[3][ShoeHexahedron::MaxDegree + 1];
  int i,j,k;
""" %( self.CellShape, self.Interpolant, self.ProductSpace )

        #eval( 'vtkUnivariate%s.SymbolicShapeFunctions( self.OutputFile )' % self.Interpolant )
        #eval( 'vtk%sProduct.SymbolicShapeFunctions( self.OutputFile )' % self.ProductSpace )
        print >> self.OutputFile, """
}
"""

    # ============
    def GenerateDerivatives( self ):
        print >> self.OutputFile, \
              'void Shoe%s::%s%sShapeFunctionDerivatives( const int* order, const double* pcoords, double* sderiv )\n{' % ( self.CellShape, self.Interpolant, self.ProductSpace )
        eval( 'vtkUnivariate%s.NumericShapeFunctionDerivatives( self.OutputFile )' % self.Interpolant )
        eval( 'vtk%sProduct.NumericShapeFunctionDerivatives( self.OutputFile )' % self.ProductSpace )
        print >> self.OutputFile, '}'

    def GenerateTail( self ):
        """
Doing nothing
"""
        # Doing nothing

    # ============
    def Go( self ):
        self.GenerateHeader()
        self.GenerateInterpolants()
        self.GenerateDerivatives()
        self.GenerateTail()




gen = vtkShoeInterpolantGenerator( 'ShoeInterpolants.cxx', 'Hexahedron', 'Lagrange', 'Tensor' )
gen.Go()
