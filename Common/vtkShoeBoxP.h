// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
// .NAME vtkShoeBoxP - private data used by vtkShoeBox to represent a mesh
// .SECTION Description
// Although each vtkShoeBox makes a vtkShoeBoxP object available to
// selected friend classes, you must never attempt to access this data
// structure directly from any code not contained in the vtksnlCommon
// library or Microsoft's incompetent compiler will become baffled.
#ifndef __vtkShoeBoxP_h
#define __vtkShoeBoxP_h

#include "vtksnlConfigure.h" // for uint32
#include "freelist" // for Cells
#include "vtkstd/map" // for Extant

class vtkShoeCellSpecies;
class vtkShoeCellRecord;
class vtkUnstructuredGrid;

class vtkShoeBoxP
{
public:
  vtkShoeBoxP();
  ~vtkShoeBoxP();

  // Description:
  // Insert a cell into the private mesh representation.
  // An integer "handle" for the cell is returned, but we
  // make no guarantee that all the cells in a mesh are
  // sequentially numbered.
  vtkIdType InsertNextCell( vtkShoeCellSpecies* species, vtkIdType connectivityEntry, uint32_t permutation );

  // Description:
  // A freelist of cells.
  // Because we use a freelist, cells may not be sequentially numbered.
  // But by giving up this constraint, deleting cells from the mesh
  // becomes an O(1) task, the same as inserting cells.
  freelist<vtkShoeCellRecord,vtkIdType> Cells;

  // Description:
  // A list of all the cell types currently in use by the mesh.
  vtkstd::map< int, vtkIdType > Extant; // pair(species ID,count) for each extant species

  // Description:
  // A triangulation of 2-D cells and any 2-D boundaries of 3-D cells.
  // One day, this mesh may also contain line segment approximations of 1-D cells, but not yet.
  vtkUnstructuredGrid* BdyTri;

  // Description:
  // This map takes a DOF node ID and returns the first triangle on
  // approximating the given DOF node plus the total number of triangles
  // in the approximation.
  vtkstd::map< vtkIdType, vtkstd::pair< vtkIdType, vtkIdType > > BdyTriAssoc;
};

#endif // __vtkShoeBoxP_h
