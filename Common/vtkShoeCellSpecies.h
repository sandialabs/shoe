// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __vtkShoeCellSpecies_h
#define __vtkShoeCellSpecies_h

#include "vtkObject.h"
#include "vtkShoeCellGenus.h"  // You won't mind including these...
#include "vtkShoeOrderTuple.h" // They're very lightweight.
#include "vtkShoeEnums.h"      // And yummy, too.

class vtkShoeCellSpeciesP;
class vtkGenericAttributeCollection;
class vtkShoeAttribute;
class vtkShoeCellMetaData;

class vtkShoeCellSpecies
{
public:
  // Description:
  // Find a species of cell that matches the requested genus and interpolant orders.
  // By creating a vector of vtkShoeOrderTuples, one for every attribute of a mesh, you
  // can add new a new cell type to an existing mesh with this call.
  static vtkShoeCellSpecies* FindOrCreate( vtkShoeCellGenus&,
    vtkPolyInterpolant* desiredInterpolants, vtkPolyProductSpace* desiredProdSpaces,
    vtkShoeOrderTuple*, int*, int, vtkGenericAttributeCollection* );
  static vtkShoeCellSpecies* GetSpeciesById( int id );

  // Description:
  // Set the order of the interpolant for this species for the requested attribute(s).
  void SetAttributeOrder( vtkGenericAttributeCollection*, vtkShoeOrderTuple*, int*, int );
  void SetAttributeOrder( vtkShoeAttribute*, vtkShoeOrderTuple& );
  int GetAttributeOrder( vtkShoeAttribute*, vtkShoeOrderTuple& );

  // Description:
  // Return the number of cells of this species
  vtkIdType GetNumberOfCells();

  // Description:
  // Increase the number of cells that reference this species.
  void Reference( vtkIdType count=1 );

  // Description:
  // Decrease the number of cells that reference this species.
  void Dereference( vtkIdType count=1 );

  // Description:
  // Get the unique ID associated with a species.
  int GetId();

  // Description:
  // Returns the number of connectivity entries (number of points + number of DOF nodes)
  // required to describe a cell of this type.
  int GetNumberOfConnectivityEntries();

  // Description:
  // Set the interpolant, product space, and order for
  // a whole collection of attributes.
  void SetInterpolantInfo( vtkGenericAttributeCollection*, int*, int,
    vtkPolyInterpolant*, vtkPolyProductSpace*, vtkShoeOrderTuple* );

  const vtkShoeCellMetaData* Meta;

protected:
  vtkShoeCellSpecies( vtkShoeCellGenus& );
  ~vtkShoeCellSpecies();

  vtkShoeCellSpeciesP* Data;
};

#endif // __vtkShoeCellSpecies_h
