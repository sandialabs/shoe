/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include "vtkToySubdivision.h"
#include "vtkAdaptiveTessellator.h"

#include <vtkObjectFactory.h>
#include <vtkImplicitFunction.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkDataSet.h>

vtkCxxRevisionMacro(vtkToySubdivision,"$Revision: 5466 $");
vtkStandardNewMacro(vtkToySubdivision);
vtkCxxSetObjectMacro(vtkToySubdivision,Refiner,vtkImplicitFunction);

vtkToySubdivision::vtkToySubdivision()
  : Refiner( 0 )
{
}

vtkToySubdivision::~vtkToySubdivision()
{
	if ( this->Refiner )
		this->Refiner->UnRegister( this );
}

void vtkToySubdivision::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
	os << indent << "Refiner: " << this->Refiner << endl;
}

bool vtkToySubdivision::EvaluateEdge( const double* p0, double* midpt, const double* p1, int field_start )
{
	static double weights[27];
	static int dummySubId=-1;
  int i;

  //cout << "( " << midpt[0] << ", " << midpt[1] << ", " << midpt[2] << " ) ";
    //<< " ) - ( " << p1[0] << ", " << p1[1] << ", " << p1[2]
    //<< " )  ";

  double d = this->Refiner->EvaluateFunction( midpt );
  double l = 0.;
  d *= d;
  if ( d < 1.e-8 )
    {
    //cout << "d=" << d << endl;
    this->CurrentCellData->EvaluateLocation( dummySubId, midpt + 3, midpt, weights );
    this->EvaluateFields( midpt, weights, field_start );
    //cout << "1" << endl;
    return true;
    }

  for ( i=0; i<3; ++i )
    {
    double tmp = (p1[i] - p0[i]);
    l += tmp*tmp;
    }

  if ( l / d > 0.5 )
    {
    //cout << "l/d**2 = " << (l/d) << ", l=" << sqrt(l) << ", d=" << sqrt(d) << endl;
    this->CurrentCellData->EvaluateLocation( dummySubId, midpt + 3, midpt, weights );
    this->EvaluateFields( midpt, weights, field_start );
    //cout << "1" << endl;
    return true;
    }

  //cout << "0" << endl;
  //cout << endl;
  return false;
}
