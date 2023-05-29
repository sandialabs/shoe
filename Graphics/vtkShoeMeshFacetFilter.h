/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-11-29 18:43:33 -0800 (Sat, 29 Nov 2003) $
  Version:   $Revision: 997 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
// .NAME vtkShoeMeshFacetFilter - Replace all $k$-d elements with their $k-1$-d bounding faces
// .SECTION Description
// Grab just the $k-1$-facets off either the boundary of the solid (in SetBoundaryElementsOnlyOn()
// mode) or all the faces of all the elements without repetition (in SetBoundaryElementsOnlyOff()
// mode). $k$ defaults to 3 but may be set using SetFacetDimension().

#ifndef __vtkShoeMeshFacetFilter_h
#define __vtkShoeMeshFacetFilter_h

#include <vtksnlConfigure.h>
#include <vtksnlGraphicsWin32Header.h>
#include <vtkShoeMeshToShoeMeshFilter.h>

class VTK_SNL_GRAPHICS_EXPORT vtkShoeMeshFacetFilter
	: public vtkShoeMeshToShoeMeshFilter
{
	public:
		static vtkShoeMeshFacetFilter* New();
		vtkTypeRevisionMacro(vtkShoeMeshFacetFilter,vtkShoeMeshToShoeMeshFilter);

		vtkGetMacro(AverageCellData, int);
		vtkSetMacro(AverageCellData, int);
		vtkBooleanMacro(AverageCellData, int);

		vtkGetMacro(BoundaryElementsOnly, int);
		vtkSetMacro(BoundaryElementsOnly, int);
		vtkBooleanMacro(BoundaryElementsOnly, int);

		vtkGetMacro(FacetDimension, int);
		vtkSetMacro(FacetDimension, int);

	protected:
		vtkShoeMeshFacetFilter();
		~vtkShoeMeshFacetFilter();

		virtual void Execute();

		virtual void ComputeFacesOfVolumes( vtkShoeMesh*, vtkShoeMesh* );
		virtual void ComputeEdgesOfFaces( vtkShoeMesh*, vtkShoeMesh* );
		virtual void ComputePointsOfEdges( vtkShoeMesh*, vtkShoeMesh* );

		int AverageCellData;
		int BoundaryElementsOnly;
		int FacetDimension;

	private:
		vtkShoeMeshFacetFilter( const vtkShoeMeshFacetFilter& ); // not implemented
		void operator = ( const vtkShoeMeshFacetFilter& ); // not implemented
};

#endif // __vtkShoeMeshFacetFilter_h
