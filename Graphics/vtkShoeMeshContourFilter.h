/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2004-08-18 11:28:39 -0700 (Wed, 18 Aug 2004) $
  Version:   $Revision: 2326 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkShoeMeshContourFilter_h
#define __vtkShoeMeshContourFilter_h

#include <vtksnlConfigure.h>
#include <vtkShoeMesh.h>
#include <vtkShoeMeshPolygonizer.h>

// .NAME vtkShoeMeshContourFilter - generates isosurfaces in a SHOE mesh.
// .SECTION Description
// This class generates polygonal approximation of isosurfaces in a higher order mesh.
// This filter requires the critical points to be calculated beforehand. The output
// may be wrong if critical points are not identified.
//
// This filter uses the vtkCellOps::Tessellate() member and a
// vtkAdaptiveTessellator/vtkShoeMeshFieldValueSubdivision pair to generate tetrahedra
// that partition each cell into regions where the marching tetrahedra algorithm
// assumptions hold. A second vtkAdaptiveTessellator is then used to extract
// the isocontour from each tetrahedron.
//
// .SECTION See Also
// vtkShoeMeshPolygonizer vtkShoeMesh vtkUnstructuredGrid

class vtkAdaptiveTessellator;
class vtkIsosurfaceSubdivisionAlgorithm;

class VTK_SNL_GRAPHICS_EXPORT vtkShoeMeshContourFilter
	: public vtkShoeMeshPolygonizer
{
public:
	vtkTypeRevisionMacro(vtkShoeMeshContourFilter,vtkShoeMeshPolygonizer);
	void PrintSelf( ostream& os, vtkIndent indent );

	static vtkShoeMeshContourFilter* New();

	virtual void SetSubdivider( vtkIsosurfaceSubdivisionAlgorithm* );

	// Description:
	// Set/Get the field to isosurface.
	// The number is the ID of the vtkFunctionData array to be used for isosurfacing.
	virtual void SetFieldId( int field );
	int GetFieldId() const;

	// Description:
	// Set/Get the component to isosurface.
	// If a field has more than one component (such as a vector or tensor field), you
	// may specify which component of the field will serve as the scalar function to
	// isosurface.
	// By default, this value is 0.
	virtual void SetFieldComponent( int component );
	int GetFieldComponent() const;

	virtual void SetIsovalue(const double value);
	double GetIsovalue() const;

protected:
	vtkShoeMeshContourFilter();
	~vtkShoeMeshContourFilter();

	// Description:
	// Run the filter; produce a polygonal approximation to the isosurface.
	virtual void Execute();

	int FieldId;
	int FieldComponent;
	double Isovalue;

private:
	vtkShoeMeshContourFilter( const vtkShoeMeshContourFilter& ); // Not implemented.
	void operator = ( const vtkShoeMeshContourFilter& ); // Not implemented.
};

//BTX

inline int   vtkShoeMeshContourFilter::GetFieldId() const        { return this->FieldId; }
inline int   vtkShoeMeshContourFilter::GetFieldComponent() const { return this->FieldComponent; }
inline double vtkShoeMeshContourFilter::GetIsovalue() const      { return this->Isovalue; }

//ETX

#endif // __vtkShoeMeshContourFilter_h

