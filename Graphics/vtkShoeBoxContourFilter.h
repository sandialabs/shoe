/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2005-10-21 15:53:07 -0700 (Fri, 21 Oct 2005) $
  Version:   $Revision: 6004 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkShoeBoxContourFilter_h
#define __vtkShoeBoxContourFilter_h

#include "vtksnlConfigure.h"
#include "vtksnlGraphicsWin32Header.h"
#include "vtkShoeBox.h"
#include "vtkUnstructuredGridAlgorithm.h"

// .NAME vtkShoeBoxContourFilter - generates isosurfaces of a vtkShoeBox mesh.
// .SECTION Description
// This class generates polygonal approximation of isosurfaces in a higher order mesh.
// Given a scalar field \f$\Phi\f$ defined over the mesh, it constructs a \f$\Phi\f$-compatible
// tessellation of the domain (stored in a vtkShoeBoxPartition).
// This tessellation is a simplicial complex that is then used to enumerate the isosurface.
// Marching tetrahedra is applied to each tetrahedron in the complex and the resulting triangles
// are then adaptively sampled to generate a better approximation of the isosurface.
//
// .SECTION See Also
// vtkShoeBox vtkUnstructuredGrid

class vtkUnstructuredGrid;
class vtkAdaptiveTessellator;
class vtkShoeBoxIsosurfaceSubdivision;
class vtkShoeBoxPartition;
class vtkSubdivisionAlgorithm;

class VTK_SNL_GRAPHICS_EXPORT vtkShoeBoxContourFilter : public vtkUnstructuredGridAlgorithm
{
public:
	vtkTypeRevisionMacro(vtkShoeBoxContourFilter,vtkUnstructuredGridAlgorithm);
	void PrintSelf( ostream& os, vtkIndent indent );
	static vtkShoeBoxContourFilter* New();

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

  vtkSetMacro(MaximumNumberOfSubdivisions,int);
	int GetMaximumNumberOfSubdivisions() const;

	virtual void SetIsovalue(const double value);
	double GetIsovalue() const;

  virtual void SetPartition( vtkShoeBoxPartition* partition );
  vtkShoeBoxPartition* GetPartition() { return this->Partition; }
  const vtkShoeBoxPartition* GetPartition() const { return this->Partition; }

  virtual void SetTessellator( vtkAdaptiveTessellator* tess );
  vtkAdaptiveTessellator* GetTessellator() { return this->Tessellator; }
  const vtkAdaptiveTessellator* GetTessellator() const { return this->Tessellator; }

	virtual void SetSubdivider( vtkShoeBoxIsosurfaceSubdivision* );
  vtkShoeBoxIsosurfaceSubdivision* GetSubdivider() { return this->Subdivider; }
  const vtkShoeBoxIsosurfaceSubdivision* GetSubdivider() const { return this->Subdivider; }

protected:
	vtkShoeBoxContourFilter();
	~vtkShoeBoxContourFilter();

	// Description:
	// Run the filter; produce a polygonal approximation to the isosurface.
	virtual int RequestData( vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector* outputVector );

  // Description:
  // Tell the pipeline that we require a vtkShoeBox object as the input.
  virtual int FillInputPortInformation( int port, vtkInformation* info );

  // Description:
  // Callback routines used by the adaptive tessellator.
  // These static functions simply recast the void* argument as a pointer to a vtkShoeBoxContourFilter
  // and call the corresponding non-static method.
  static void AddTet( const double*, const double*, const double*, const double*, vtkSubdivisionAlgorithm*, void*, const void* );
  static void AddTri( const double*, const double*, const double*, vtkSubdivisionAlgorithm*, void*, const void* );
  static void AddLine( const double*, const double*, vtkSubdivisionAlgorithm*, void*, const void* );

  // Description:
  // Non-static callback routines used by the adaptive tessellator.
  void AddTet( const double* x, const double* y, const double* z, const double* w );
  void AddTri( const double* x, const double* y, const double* z );
  void AddLine( const double* x, const double* y );

  // Description:
  // Find a single edge-contour intersection (as opposed to all possible intersections) on a line segment in 3-D parameter-space.
  // The line segment is defined by \a p0 and \a p1.
  // The isocontour value stored in \a this->Isovalue is used.
  // The resulting intersection is stored in \a pint.
  // The return value is 1 for success and 0 for failure. Failure occurs when the algorithm cannot find a root
  // on the given line segment within 15 iterations.
  int EdgeContourIntersect( const double* p0, const double* p1, double* pint );

	int FieldId;
	int FieldComponent;
	int MaximumNumberOfSubdivisions;
	double Isovalue;
  double Tolerance; //!< How close (in parameter space) we must be to an isocontour before stopping. Currently hardwired to 1.e-4.
  vtkShoeCell* Cell; //!< Current cell being traversed. Only valid in RequestData() calls. Not reference counted.
  vtkShoeAttribute* Attribute; //!< Attribute corresponding to FieldId. Only valid in RequestData() calls. Not reference counted.
  vtkShoeBoxPartition* Partition; //!< A partition of the mesh compatible with the scalar field to be contoured.
  vtkAdaptiveTessellator* Tessellator; //!< The tessellator used to adaptively sample the isocontour (not the cell).
	vtkShoeBoxIsosurfaceSubdivision* Subdivider; //!< A subdivision algorithm that integrates edge midpoints to the isosurface.
  vtkUnstructuredGrid* ContourOut; //!< A mesh containing the output isosurface. Only valid during RequestData() calls.

private:
	vtkShoeBoxContourFilter( const vtkShoeBoxContourFilter& ); // Not implemented.
	void operator = ( const vtkShoeBoxContourFilter& ); // Not implemented.
};

//BTX

inline int   vtkShoeBoxContourFilter::GetFieldId() const        { return this->FieldId; }
inline int   vtkShoeBoxContourFilter::GetFieldComponent() const { return this->FieldComponent; }
inline int   vtkShoeBoxContourFilter::GetMaximumNumberOfSubdivisions() const { return this->MaximumNumberOfSubdivisions; }
inline double vtkShoeBoxContourFilter::GetIsovalue() const      { return this->Isovalue; }

//ETX

#endif // __vtkShoeBoxContourFilter_h

