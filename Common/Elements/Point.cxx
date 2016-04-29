/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkCellOps.h>
#include <Elements/Point.h>



void vtkShoeElemPoint::EvaluateShapeFunctionPoint( double* shape, vtkShoeMeshIterator&, const int[3], const double[3] )
{
	shape[0] = 1.0;
}
