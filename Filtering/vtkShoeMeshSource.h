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
#ifndef vtkShoeMeshSource_h
#define vtkShoeMeshSource_h
// .NAME vtkShoeMeshSource - An abstract class for filters that generate vtkShoeMesh objects
// .SECTION Description
// Inherit this class if you generate vtkShoeMeshes.

#include <vtkSource.h>

#include <vtksnlConfigure.h>
#include <vtksnlFilteringWin32Header.h>

class vtkDataObject;

class vtkShoeMesh;

class VTK_SNL_FILTERING_EXPORT vtkShoeMeshSource
	: public vtkSource
{
	public:
		vtkTypeRevisionMacro(vtkShoeMeshSource,vtkSource);
		void PrintSelf( ostream&, vtkIndent );

		// Description:
		// Get or set the output of the source.
		// Without an argument, the 0-th output is referenced.
		vtkShoeMesh* GetOutput();
		const vtkShoeMesh* GetOutput() const;
		vtkShoeMesh* GetOutput( int );
		const vtkShoeMesh* GetOutput( int ) const;
		void SetOutput( vtkShoeMesh* );

	protected:
		vtkShoeMeshSource();
		~vtkShoeMeshSource();

		void ComputeInputUpdateExtents(vtkDataObject *output);

	private:
		vtkShoeMeshSource( const vtkShoeMeshSource& ); // Not implemented: verboten
		vtkShoeMeshSource& operator = ( const vtkShoeMeshSource& ); // Not implemented: verboten
};

#endif // vtkShoeMeshSource_h
