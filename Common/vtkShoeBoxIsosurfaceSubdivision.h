// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __vtkShoeBoxIsosurfaceSubdivision_h
#define __vtkShoeBoxIsosurfaceSubdivision_h

// .NAME vtkShoeBoxIsosurfaceSubdivision - Adaptively refine an isocontour
// .SECTION Description
// This subdivision algorithm evaluates a scalar field value (or a single component of a vector/tensor value) at
// the midpoint of each edge it is passed.
// When the value exceeds the given tolerance, it uses Newton iteration to converge to the zero of the field
// by following the gradient of the field.

#include "vtksnlCommonWin32Header.h"
#include "vtkSubdivisionAlgorithm.h"

#include <vtkstd/vector>

class vtkShoeBox;
class vtkShoeCell;
class vtkShoeAttribute;

class vtkShoeBoxIsosurfaceSubdivision : public vtkSubdivisionAlgorithm
{
public:
  static vtkShoeBoxIsosurfaceSubdivision* New();
  virtual void PrintSelf( ostream& os, vtkIndent indent );
  vtkTypeRevisionMacro(vtkShoeBoxIsosurfaceSubdivision,vtkSubdivisionAlgorithm);

  virtual bool EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start );

  // Description:
  // Evaluate each of the attributes at the updated edge midpoint.
  void EvaluateFields( double* vertex, int field_start );

  virtual void SetMesh( vtkShoeBox* m ) { this->Mesh = m; }
  virtual void SetCell( vtkShoeCell* c ) { this->Cell = c; }
  virtual void SetAttribute( vtkShoeAttribute* a );
  virtual void SetComponent( int c ) { this->Component = c; }
  virtual void SetIsovalue( double v ) { this->Isovalue = v; }
  virtual void SetTolerance( double t ) { this->Tolerance = t; }

protected:
  vtkShoeBoxIsosurfaceSubdivision();
  ~vtkShoeBoxIsosurfaceSubdivision();

  /** Determine the best feasible direction to move the midpoint.
   * Starts with the gradient and prevents movement along some
   * directions.
   */
  bool RestrictPathToCellInterior( const double e0[3], double r[3], const double e1[3], double dr[3], int obc ) const;
  /** Returns a bit vector indicating which cell boundaries \a lies on.
   */
  int BoundaryClassify( double r[3] ) const;

  vtkShoeCell* Cell;
  vtkShoeBox* Mesh;
  vtkShoeAttribute* Attribute;
  int Component;
  double Isovalue;
  double Tolerance;
  //BTX
  vtkstd::vector<double> MidPtScalar;
  vtkstd::vector<double> MidPtDerivs;
  vtkstd::vector<double> EndPt0Derivs;
  vtkstd::vector<double> EndPt1Derivs;
  //ETX
};

#endif // __vtkShoeBoxIsosurfaceSubdivision_h
