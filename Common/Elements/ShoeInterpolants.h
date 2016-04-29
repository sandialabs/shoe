// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __ShoeInterpolants_h
#define __ShoeInterpolants_h

#include "vtksnlConfigure.h"
#include "vtkShoeOrderTuple.h" // for LagrangePrepareForOrder

#include <vtkstd/vector>

class ShoeInterpolants
{
public:
  static void LagrangePrepareForOrder( const vtkShoeOrderTuple& );

  static vtkstd::vector< vtkstd::vector<int> > BinomialCoefficients;
  static vtkstd::vector< vtkstd::vector<double> > LagrangeNormalizationFactors;
  static vtkstd::vector< vtkstd::vector< vtkstd::vector< double > > > LagrangeInterpolants;
  static vtkstd::vector< vtkstd::vector< vtkstd::vector< double > > > LagrangeDerivatives;
};

#endif // __ShoeInterpolants_h

