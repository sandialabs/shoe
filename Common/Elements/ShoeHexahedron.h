// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __ShoeHexahedron_h
#define __ShoeHexahedron_h

#include "vtksnlConfigure.h"
#include "vtkSystemIncludes.h"

class vtkShoeBox;
class vtkDataRecords;
class vtkPolynomialSystem;

class ShoeHexahedron
{
public:
  static int VerticesOnFace[6];
  static int FaceArray[6][4];
  static int EdgeArray[12][2];
  static int FacesOnEdge[12][2];
  static int EdgesOnFace[6][4];
  static double ParametricCenter[3];
  static double ParametricCorners[24];
  static double BoundaryEdgeTransforms[12][6];
  static double BoundaryFaceTransforms[6][12];

  static void LagrangeTensorPermuteEdge( double* CellCache, int NumberOfComponents, vtkDataRecords* dof, vtkIdType dconni, int eord, bool perm );
  static void LagrangeTensorPermuteFace( double* CellCache , int NumberOfComponents, vtkDataRecords* dof, vtkIdType dconni, int *ford, int fdir[2], uint32_t perm );
  static void LagrangeTensorShapeFunctions( const int* order, const double* pcoords, double* shape );
  static void LagrangeTensorShapeFunctionDerivatives( const int* order, const double* pcoords, double* sderiv );
  static void SymbolicLagrangeFieldGradient( const int* order, const double* phi_e, vtkPolynomialSystem* grad );
  static void EmbedEdgeCoords( int edge, double* tuple );
  static void EmbedEdgeCoordsInFace( int edge, int face, double* tuple );
  static void EmbedFaceCoords( int face, double* tuple );
  static void ProjectCoordsToEdge( int edge, double* tuple );
  static void ProjectCoordsToFace( int face, double* tuple );
  static void ProjectFaceCoordsToEdge( int face, int edge, double* tuple );
  static void TransformFaceCoordsToStorage( uint32_t perm, double* tuple );

  static int  IsParameterInDomain( double* param, int bdy );

  enum 
    {
    MaxDegree = 10
    };
  static double LagrangeScratch[3][MaxDegree + 1];
};

#endif // __ShoeHexahedron_h
