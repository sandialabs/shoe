// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeBoxIsosurfaceSubdivision.h"

#include "vtkObjectFactory.h"
#include "vtkGenericAttributeCollection.h"

#include "vtkShoeAttribute.h"
#include "vtkShoeBox.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellMetaData.h"

#undef VTK_DBG_ISOSUBDIV

vtkStandardNewMacro(vtkShoeBoxIsosurfaceSubdivision);
vtkCxxRevisionMacro(vtkShoeBoxIsosurfaceSubdivision,"$Revision: 6012 $");

void vtkShoeBoxIsosurfaceSubdivision::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Cell:      " << this->Cell << vtkstd::endl;
  os << indent << "Mesh:      " << this->Mesh << vtkstd::endl;
  os << indent << "Attribute: " << this->Attribute << vtkstd::endl;
  os << indent << "Isovalue:  " << this->Isovalue << vtkstd::endl;
}

vtkShoeBoxIsosurfaceSubdivision::vtkShoeBoxIsosurfaceSubdivision()
{
  this->Cell = 0;
  this->Mesh = 0;
  this->Attribute = 0;
  this->Component = 0;
  this->Isovalue = 0.;
  this->Tolerance = 1.e-5;
}

vtkShoeBoxIsosurfaceSubdivision::~vtkShoeBoxIsosurfaceSubdivision()
{
}

bool vtkShoeBoxIsosurfaceSubdivision::EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start )
{
  bool toSubdivideOrNotToSubdivideThatIsTheQuestion = false;
  int cnt = 0;
  int origBdyClass = this->BoundaryClassify( p1 + 3 );
  this->Cell->InterpolateTuple( this->Attribute, p1 + 3, &this->MidPtScalar[0] );
  while ( (cnt < 15) && fabs( this->MidPtScalar[ this->Component ] - this->Isovalue ) > this->Tolerance )
    {
    toSubdivideOrNotToSubdivideThatIsTheQuestion = true;
    this->Cell->ParametricDerivatives( p1 + 3, this->Attribute, &this->MidPtDerivs[0] );
    double dr[3]; // direction of line the midpoint will be moved along.
    // Fill dr with MidPtDerivs, clipping if p1 is on 1 or 2 cell boundary faces.
    // Then step along dr to expected isosurface location.
    if ( ! this->RestrictPathToCellInterior( p0 + 3, p1 + 3, p2 + 3, dr, origBdyClass ) )
      {
      vtkstd::cerr << "Unable to move from ("
        << p1[3] << ", " << p1[4] << ", " << p1[5] << "), grad ("
        << this->MidPtDerivs[this->Component*3 + 0] << ", "
        << this->MidPtDerivs[this->Component*3 + 1] << ", "
        << this->MidPtDerivs[this->Component*3 + 2] << ")" << vtkstd::endl;
      return false;
      }
    this->Cell->InterpolateTuple( this->Attribute, p1 + 3, &this->MidPtScalar[0] );
#ifdef VTK_DBG_ISOSUBDIV
    vtkstd::cout
      << "   " << "g*(" << p1[3] << ", " << p1[4] << ", " << p1[5] << ")=" << this->MidPtScalar[0]
      << vtkstd::endl;
#endif // VTK_DBG_ISOSUBDIV
    ++cnt;
    }
  if ( ! toSubdivideOrNotToSubdivideThatIsTheQuestion )
    {
    // OK, the midpoint is on the surface, but what about the gradients at each
    // endpoint? If they are off by much, the surface normals will look funny when
    // rendered...
    //double e0grad[9];
    //double e1grad[9];
    double v[3];
    double vmag;
    int m;
    //this->Cell->ParametricDerivatives( p0 + 3, this->Mesh->GetGeometry(), e0grad );
    this->Cell->ParametricDerivatives( const_cast<double*>(p0) + 3, this->Attribute, &this->EndPt0Derivs[0] );
    vmag = 0.;
    for ( m = 0; m < 3; ++m )
      {
      v[m] = this->EndPt0Derivs[3 * this->Component + m];
      vmag += v[m]*v[m];
      }
    if ( vmag > 1.e-6 )
      {
      //this->Cell->ParametricDerivatives( p2 + 3, this->Mesh->GetGeometry(), e1grad );
      this->Cell->ParametricDerivatives( const_cast<double*>(p2) + 3, this->Attribute, &this->EndPt1Derivs[0] );
      double v2mag = 0.;
      for ( m = 0; m < 3; ++m )
        {
        v[m] = this->EndPt1Derivs[3 * this->Component + m];
        v2mag += v[m]*v[m];
        }
      if ( v2mag > 1.e-6 )
        {
        v2mag *= vmag;
        vmag = 0;
        // OK, both normals are nonzero, but are they close to each other? The dot product will tell:
        for ( m = 0; m < 3; ++m )
          {
          vmag += this->EndPt0Derivs[3 * this->Component + m] * v[m];
          }
        if ( 1. - fabs( vmag / sqrt(v2mag) ) < 0.8 )
          {
          toSubdivideOrNotToSubdivideThatIsTheQuestion = true;
          }
        }
      else
        {
        toSubdivideOrNotToSubdivideThatIsTheQuestion = true;
        }
      }
    else
      { // zero gradient means an undefined normal. Always subdivide near these; they are degenerate.
      toSubdivideOrNotToSubdivideThatIsTheQuestion = true;
      }
    }
  if ( toSubdivideOrNotToSubdivideThatIsTheQuestion )
    {
    this->Cell->EvaluateLocation( -1, p1 + 3, p1 );
    this->EvaluateFields( p1, field_start );
    }
  return toSubdivideOrNotToSubdivideThatIsTheQuestion;
}

void vtkShoeBoxIsosurfaceSubdivision::EvaluateFields( double* vertex, int field_start )
{
  vtkGenericAttributeCollection* attribs = this->Cell->GetDataSet()->GetAttributes();
  for ( int f=0; f<this->GetNumberOfFields(); ++f )
    this->Cell->InterpolateTuple( attribs->GetAttribute(this->GetFieldIds()[f]), vertex + 3, vertex + field_start + this->GetFieldOffsets()[f] );
}

void vtkShoeBoxIsosurfaceSubdivision::SetAttribute( vtkShoeAttribute* a )
{
  vtkSetObjectBodyMacro(Attribute,vtkShoeAttribute,a);
  if ( a )
    {
    int nc = a->GetNumberOfComponents();
    this->MidPtScalar.reserve( nc );
    this->MidPtDerivs.reserve( 3*nc );
    this->EndPt0Derivs.reserve( 3*nc );
    this->EndPt1Derivs.reserve( 3*nc );
    if ( this->Component >= nc )
      {
      vtkWarningMacro("New attribute \"" << a->GetName() << "\" only has "
        << nc << " components, but Component was set to " << this->Component
        << ". Resetting Component to 0.");
      this->Component = 0;
      }
    }
}

// Fills dr with MidPtDerivs, but will clip if p1 is on 1 or 2 cell boundary faces.
bool vtkShoeBoxIsosurfaceSubdivision::RestrictPathToCellInterior(
  const double e0[3], double r[3], const double e1[3], double dr[3], int obc ) const
{
  double y;
  const vtkShoeCellMetaData* meta = this->Cell->GetMetaData();

  dr[0] = this->MidPtDerivs[3 * this->Component + 0];
  dr[1] = this->MidPtDerivs[3 * this->Component + 1];
  dr[2] = this->MidPtDerivs[3 * this->Component + 2];

  int to = obc;
  int f = 0;
  while ( obc )
    {
    if ( obc & 1 )
      {
      y = -( dr[0]*meta->BoundaryFaceTransforms[f][ 9] +
             dr[1]*meta->BoundaryFaceTransforms[f][10] +
             dr[2]*meta->BoundaryFaceTransforms[f][11] );
      for ( int i = 0; i < 3; ++i )
        dr[i] += y*meta->BoundaryFaceTransforms[f][9+i];
      }
    obc >>= 1;
    ++f;
    }
  if ( fabs( dr[0] ) < 1e-4 && fabs( dr[1] ) < 1e-4 && fabs( dr[2] ) < 1e-4 )
    {
    return 0; // no good directions
    }
  obc = to;

  int bc = this->BoundaryClassify( r );
  f = 0;
  while ( bc ^ obc )
    {
    if ( bc & 1 )
      {
      y = -( dr[0]*meta->BoundaryFaceTransforms[f][ 9] +
             dr[1]*meta->BoundaryFaceTransforms[f][10] +
             dr[2]*meta->BoundaryFaceTransforms[f][11] );
      // Only constrain the point from moving outside the boundary
      if ( y * this->MidPtScalar[this->Component] > 0 )
        {
        for ( int i = 0; i < 3; ++i )
          dr[i] += y*meta->BoundaryFaceTransforms[f][9+i];
        }
      }
    obc >>= 1;
    bc >>= 1;
    ++f;
    }
  if ( fabs( dr[0] ) < 1e-4 && fabs( dr[1] ) < 1e-4 && fabs( dr[2] ) < 1e-4 )
    {
    return 0; // no good directions
    }

  int z;
  double x[3], edir[3];
  y = 0.;
  for ( z = 0; z < 3; ++z )
    {
    edir[z] = e1[z] - e0[z];
    y += edir[z]*edir[z];
    }
  y = sqrt( y );
  double w = 0.;
  for ( z = 0; z < 3; ++z )
    {
    w += (edir[z] = edir[z] / y) * dr[z];
    }
  for ( z = 0; z < 3; ++z )
    {
    x[z] = dr[z] - w * edir[z];
    }
  if ( fabs( x[0] ) > 1e-4 || fabs( x[1] ) > 1e-4 || fabs( x[2] ) > 1e-4 )
    {
    // Try moving normal to the line as long as we can...
    dr[0] = x[0]; dr[1] = x[1]; dr[2] = x[2];
    }

  double graddirdotprod;
  graddirdotprod =
    this->MidPtDerivs[3 * this->Component + 0] * dr[0] +
    this->MidPtDerivs[3 * this->Component + 1] * dr[1] +
    this->MidPtDerivs[3 * this->Component + 2] * dr[2];
  double t = (this->MidPtScalar[this->Component] - this->Isovalue) / graddirdotprod;
#ifdef VTK_DBG_ISOSUBDIV
  vtkstd::cout
    << ": g(" << r[0] << ", " << r[1] << ", " << r[2] << ")=" << this->MidPtScalar[0]
    << ", gradgsq=" << graddirdotprod << ", t=" << t << vtkstd::endl;
#endif // VTK_DBG_ISOSUBDIV
  r[0] -= dr[0] * t;
  r[1] -= dr[1] * t;
  r[2] -= dr[2] * t;
  return true;
}

int vtkShoeBoxIsosurfaceSubdivision::BoundaryClassify( double r[3] ) const
{
  int bitvector = 0;
  int curbit = 1;
  const vtkShoeCellMetaData* meta = this->Cell->GetMetaData();
  int nf = meta->NumberOfBoundaries[2];
  for ( int f = 0; f < nf; ++f, curbit <<= 1 )
    {
    double d = 0.;
    for ( int c = 0; c < 3; ++c )
      {
      d += (r[c] - meta->BoundaryFaceTransforms[f][c])*meta->BoundaryFaceTransforms[f][c+9];
      }
    if ( d < 1.e-6 )
      {
      bitvector += curbit;
      }
    }
  return bitvector;
}
