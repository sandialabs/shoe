// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __vtkShoeCellIterator_h
#define __vtkShoeCellIterator_h

#include "vtkGenericCellIterator.h"
#include "vtkShoeEnums.h" // for traversal style

class vtkShoeBox;
class vtkShoeCell;
class vtkDataRecordsIterator;
class vtkShoeCellIteratorP;

class vtkShoeCellIterator : public vtkGenericCellIterator
{
public:
  vtkTypeRevisionMacro(vtkShoeCellIterator,vtkGenericCellIterator);
  void PrintSelf( ostream& os, vtkIndent indent );
  static vtkShoeCellIterator* New();

  virtual void Begin();
  virtual int IsAtEnd();
  virtual vtkGenericAdaptorCell* NewCell();
  virtual void GetCell( vtkGenericAdaptorCell* c );

  // Description:
  // Returns the ID of the cell currently referenced by the iterator (or -1 if
  // the iterator is not pointing to a valid cell, as can happen before initialization
  // or after a traversal has been completed).
  vtkIdType GetCellId();

  // Description:
  // Returns the cell currently referenced by the iterator.
  // You should <strong>not</strong> delete the cell when you're done; it
  // is managed by the iterator.
  virtual vtkGenericAdaptorCell* GetCell();

  // Description:
  // Advance to the next cell.
  virtual void Next();

  int GetTraversalStyle() const;

  int EndCheckMeshOrder();
  int EndCheckBoundary();
  int EndCheckExtBoundary();
  int EndCheckCellBoundary();

  void NextCellMeshOrder();
  void NextCellBoundary();
  void NextCellExtBoundary();
  void NextCellCellBoundary();

  virtual void SetDataSet( vtkShoeBox* );
  vtkShoeBox* GetDataSet() { return this->DataSet; }
  const vtkShoeBox* GetDataSet() const { return this->DataSet; }

protected:
  vtkShoeCellIterator();
  virtual ~vtkShoeCellIterator();

  //BTX
  friend class vtkShoeCell;
  friend class vtkShoeBox;
  //ETX

  void SetCellById( vtkIdType cellID );
  void SetBoundaryCellByDOF( vtkIdType DOF );
  void SetParentCell( vtkShoeCell* );

  vtkShoeBox* DataSet;
  vtkShoeCellIteratorP* Iter;
  vtkShoeCell* CurrentCell;
  vtkShoeCell* CurrentParent; // CurrentCell is facet CurrentFacet of CurrentParent when traversing boundaries
  int CurrentFacet; // for Boundary or ExtBoundary iteration, which face of CurrentParent are we on?
  int Dimension; // for Boundary or ExtBoundary iteration, which boundaries do we consider?

  MeshTraversalStyle TraversalStyle;
  int TraversalMask; // used with "custom" traversal style

private:
  vtkShoeCellIterator( const vtkShoeCellIterator& ); // Not implemented.
  void operator =( const vtkShoeCellIterator& ); // Not implemented.
};

//BTX
inline int vtkShoeCellIterator::GetTraversalStyle() const { return this->TraversalStyle; }
//ETX

#endif // __vtkShoeCellIterator_h
