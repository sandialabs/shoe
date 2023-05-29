/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include "vtkFlowerSubdivision.h"
#include "vtkAdaptiveTessellator.h"

#include <vtkObjectFactory.h>
#include <vtkImplicitFunction.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkDataSet.h>

vtkCxxRevisionMacro(vtkFlowerSubdivision,"$Revision: 3149 $");
vtkStandardNewMacro(vtkFlowerSubdivision);
vtkCxxSetObjectMacro(vtkFlowerSubdivision,Refiner,vtkImplicitFunction);

vtkFlowerSubdivision::vtkFlowerSubdivision()
  : Refiner( 0 )
{
}

vtkFlowerSubdivision::~vtkFlowerSubdivision()
{
	if ( this->Refiner )
		this->Refiner->UnRegister( this );
}

void vtkFlowerSubdivision::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
	os << indent << "Refiner: " << this->Refiner << endl;
}

bool vtkFlowerSubdivision::EvaluateEdge( const double* p0, double* midpt, const double* p1, int field_start )
{
	static double weights[27];
	static int dummySubId=-1;
	double realMidPt[ 3 ];
  int i;

  //cout << "( " << midpt[0] << ", " << midpt[1] << ", " << midpt[2] << " ) ";
    //<< " ) - ( " << p1[0] << ", " << p1[1] << ", " << p1[2]
    //<< " )  ";

  double x,y,z,x2,y2,z2;

  x = (midpt[0] -  97.) / 95.;
  y = (midpt[1] - 100.) / 56.;
  z = (midpt[2] -  60.) / 60.;
  x2 = x * x;
  y2 = y * y;
  z2 = z * z;

  double d = 0.;

  d = sqrt(x2 + y2 +z2); // first of all, compute r
  if ( d < 1.e-8 )
    {
    //cout << "d=" << d << endl;
    this->CurrentCellData->EvaluateLocation( dummySubId, midpt + 3, midpt, weights );
    this->EvaluateFields( midpt, weights, field_start );
    //cout << "1" << endl;
    return true;
    }
  // now that we know z/d isn't heading to infinity, use r to compute the real distance function:
  d = sin( 3.*atan2(y,x) )*sin( 4.*acos(z/d) ) - d;

  double l = 0.;
  for ( i=0; i<3; ++i )
    {
    double tmp = (p1[i] - p0[i]);
    l += tmp*tmp;
    }

  if ( l / d > 0.0001 )
    {
    this->CurrentCellData->EvaluateLocation( dummySubId, midpt + 3, midpt, weights );
    this->EvaluateFields( midpt, weights, field_start );
    //cout << "1" << endl;
    return true;
    }

  //cout << "0" << endl;
  //cout << endl;
  return false;
}
