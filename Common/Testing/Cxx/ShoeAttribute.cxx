// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeAttribute.h"
#include "vtkDataRecords.h"
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

static double ShoeAttributePts[] =
{
  0., 0., 0.,
  1., 0., 0.,
  0., 1., 0.,
  0., 0., 1.
};

static double ShoeAttributeDOF[] =
{
  0.5 , 0.  , 0.  ,
  0.5 , 0.5 , 0.  ,
  0.  , 0.5 , 0.  ,
  0.  , 0.  , 0.5 ,
  0.5 , 0.  , 0.5 ,
  0.  , 0.5 , 0.5 ,
  1/3., 1/3., 0.  ,
  1/3., 0.  , 1/3.,
  0.  , 1/3., 1/3.,
  0.5 , 0.5 , 0.5
};

int ShoeAttribute( int argc, char** argv )
{
  vtkShoeAttribute* att = vtkShoeAttribute::New();
  att->SetName("Velocity");
  att->SetNumberOfComponents(3);

  int i;
  vtkDoubleArray* pts = vtkDoubleArray::New();
  vtkDataRecords* dof = vtkDataRecords::New();

  pts->SetNumberOfComponents(3);
  pts->SetNumberOfTuples(4);

  for (i=0; i<4; i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof->SetNumberOfComponents(3);
  for (i=0; i<10; i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att->SetPointData( pts );
  att->SetDOFData( dof );

  return 0;
}
