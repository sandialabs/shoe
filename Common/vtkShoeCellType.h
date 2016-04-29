/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-11-29 19:05:34 -0800 (Sat, 29 Nov 2003) $
  Version:   $Revision: 1007 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
//	.NAME vtkCellType - Represent the unique characteristics that define a class of cells
//	.SECTION Description
//	vtkCellType contains information about
//	<ul>
//	<li> the shape of the parameter-space domain over which the cell is defined,
//	<li> the family of polynomial interpolants used to map this domain into world space, and
//	<li> the product space of one-dimensional shape functions used to extend univariate
//		polynomials into a multivariate interpolant. (This is unused for 1D cells.)
//	</ul>
//	There is a 1-to-1 mapping from vtkCellType to vtkCellOps objects; in other words, the
//	cell type (and only the cell type) determines the set of functions that operate on
//	the class of cells the vtkCellType describes.
//
//	.SECTION See Also
//	vtkCellOps, vtkCellDefinition, vtkShoeMesh
#ifndef vtkCellType_h
#define vtkCellType_h

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>

#include <vtkCellEnums.h>

struct VTK_SNL_COMMON_EXPORT vtkCellType
{
	// Description:
	// The shape of the parameter-space domain of the
	// specified cells.
	shoe::CellShape DomainShape;

	// Description:
	// The family of polynomials used to interpolate
	// function values over the cell.
	PolyInterpolant Interpolant;

	// Description:
	// The rule by which univariate shape functions are
	// combined into multivariate shape functions.
	PolyProductSpace ProductSpace;

	// Description:
	// Are the types equivalent?
	bool operator == ( const vtkCellType& other ) const;
	bool operator != ( const vtkCellType& other ) const;
};

//BTX

inline bool vtkCellType::operator != ( const vtkCellType& other ) const
{
	return (!( *this == other ));
}
inline bool vtkCellType::operator == ( const vtkCellType& other ) const
{
	return ( (this->DomainShape == other.DomainShape) &&
	         (this->Interpolant == other.Interpolant) &&
	         (this->ProductSpace == other.ProductSpace ) );
}
//ETX

#endif // vtkCellType_h
