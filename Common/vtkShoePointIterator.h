// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __vtkShoePointIterator_h
#define __vtkShoePointIterator_h

#include "vtkGenericPointIterator.h"

class vtkDataArray;
class vtkShoeCell;
class vtkShoeBox;

class vtkShoePointIterator : public vtkGenericPointIterator
{
public:
  vtkTypeRevisionMacro(vtkShoePointIterator,vtkGenericPointIterator);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkShoePointIterator* New();

  virtual void Begin();
  virtual int IsAtEnd();
  virtual void Next();
  virtual double* GetPosition();
  virtual void GetPosition(double*);
  virtual vtkIdType GetId();

  // Description:
  // These calls allow you to describe which set of points should be iterated.
  // In addition to the points of a dataset or cell, an arbitrary set of points may be traversed.
  virtual void SetDataSet( vtkShoeBox* );
  virtual void SetCell( vtkShoeCell* );
  virtual void SetSelection( vtkDataArray* Points, vtkIdType N, vtkIdType* Subset );

protected:
  vtkShoePointIterator();
  virtual ~vtkShoePointIterator();

  vtkShoeBox* DataSet;
  vtkDataArray* Points;
  vtkIdType Id;
  vtkIdType* List;
  vtkIdType ListSize;

private:
  vtkShoePointIterator( const vtkShoePointIterator& ); // Not implemented.
  void operator =( const vtkShoePointIterator& ); // Not implemented.
};

#endif // __vtkShoePointIterator_h
