/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2005-08-08 11:28:08 -0700 (Mon, 08 Aug 2005) $
  Version:   $Revision: 5405 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
// .NAME vtkShoeCellGenus - Represent the unique characteristics that define a class of cells
// .SECTION Description
// At the moment, a cell's "genus" is simply its shape -- i.e., the
// number and connectivity of "corner" vertices.
// For example, a cell with 8 corners arranged into a solid with 6 faces
// and 12 edges is a hexahedron.
// 
// Each genus of cell will have several species associated with it.
// Species are differentiated by the way in which they interpolate a
// given set of attributes.
// For instance, a hexahedron might have position and momentum attributes.
// One species of hexahedron might interpolate position and momentum with
// a tri-linear Lagrange basis while another might interpolant position
// with a tensor-product NURBs spline and momentum with a quadratic
// hierarchic Legenedre polynomial. Both species have the same genus.
#ifndef __vtkShoeCellGenus_h
#define __vtkShoeCellGenus_h

#include "vtkShoeEnums.h"

class vtkShoeCellGenus
{
public:

  shoe::CellShape     Shape;

	// Description:
	// Are the types equivalent?
	bool operator == ( const vtkShoeCellGenus& other ) const;
	bool operator != ( const vtkShoeCellGenus& other ) const;

  // Description:
  // Return the number of connectivity entries required to describe
  // a cell of this genus. This number depends only on the shape
  // and is equal to the number of points plus the number of DOF nodes.
  int GetNumberOfConnectivityEntries() const;

  // This is implemented in vtkShoeCellSpecies.cxx to avoid requiring
  // a C++ file for one class-static member.
  static const int NumberOfConnectivityEntriesByShape[];
};

//BTX
inline bool vtkShoeCellGenus::operator != ( const vtkShoeCellGenus& other ) const
{
	return (!( *this == other ));
}

inline bool vtkShoeCellGenus::operator == ( const vtkShoeCellGenus& other ) const
{
	return ( this->Shape == other.Shape );
}

inline int vtkShoeCellGenus::GetNumberOfConnectivityEntries() const
{
  return vtkShoeCellGenus::NumberOfConnectivityEntriesByShape[ this->Shape ];
}
//ETX

#endif // __vtkShoeCellGenus_h
