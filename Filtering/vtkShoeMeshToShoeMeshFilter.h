/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-09-30 19:04:45 -0700 (Tue, 30 Sep 2003) $
  Version:   $Revision: 759 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef vtkShoeMeshToShoeMeshFilter_h
#define vtkShoeMeshToShoeMeshFilter_h

// .NAME vtkShoeMeshToShoeMeshFilter - Base class for filters that process vtkShoeMesh objects
// .SECTION Description
// This class should be inherited by all filters that create an altered vtkShoeMesh
// from some initial vtkShoeMesh.
//
// .SECTION See Also
// vtkShoeMeshFacetFilter, vtkShoeMeshPolygonizer, vtkShoeMeshWarp, vtkShoeMeshCullFilter

#include <vtksnlConfigure.h>
#include <vtksnlFilteringWin32Header.h>
#include <vtkShoeMeshSource.h>

class vtkShoeMesh;

class VTK_SNL_FILTERING_EXPORT vtkShoeMeshToShoeMeshFilter
	: public vtkShoeMeshSource
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshToShoeMeshFilter,vtkShoeMeshSource);
		void PrintSelf( ostream&, vtkIndent );

		// Description:
		// Set or get the filter input.
		void SetInput( vtkShoeMesh* );
		vtkShoeMesh* GetInput();
		const vtkShoeMesh* GetInput() const;

	protected:
		vtkShoeMeshToShoeMeshFilter();
		~vtkShoeMeshToShoeMeshFilter();

	private:
		vtkShoeMeshToShoeMeshFilter( const vtkShoeMeshToShoeMeshFilter& ); // Not implemented: verboten
		vtkShoeMeshToShoeMeshFilter& operator = ( const vtkShoeMeshToShoeMeshFilter& ); // Not implemented: verboten
};

#endif // vtkShoeMeshToShoeMeshFilter_h

