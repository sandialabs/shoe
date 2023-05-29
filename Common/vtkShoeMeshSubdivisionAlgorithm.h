/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-12-31 21:15:12 -0800 (Wed, 31 Dec 2003) $
  Version:   $Revision: 1201 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef __vtkShoeMeshSubdivisionAlgorithm_h
#define __vtkShoeMeshSubdivisionAlgorithm_h
// .NAME vtkShoeMeshSubdivisionAlgorithm - a subclass of vtkSubdivisionAlgorithm for SHOE meshes
//
// .SECTION Description
// This is a subclass of vtkSubdivisionAlgorithm that is used for
// tessellating cells of a SHOE mesh.
//
// It provides functions for setting the current cell being tessellated and a
// convenience routine, \a EvaluateFields() to evaluate field values at a
// point. You should call \a EvaluateFields() from inside \a EvaluateEdge()
// whenever the result of \a EvaluateEdge() will be true. Otherwise, do
// not call \a EvaluateFields() as the midpoint is about to be discarded.
// (<i>Implementor's note</i>: This isn't true if UGLY_ASPECT_RATIO_HACK
// has been defined. But in that case, we don't want the exact field values;
// we need the linearly interpolated ones at the midpoint for continuity.)
//
// .SECTION See Also
// vtkShoeMeshChordLengthSubdivision vtkShoeMeshViewDependentSubdivision

#include <vtkSubdivisionAlgorithm.h>

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>
#include <vtkShoeMeshIterator.h>

class VTK_SNL_COMMON_EXPORT vtkShoeMeshSubdivisionAlgorithm
  : public vtkSubdivisionAlgorithm
{
  public:
    vtkTypeRevisionMacro(vtkShoeMeshSubdivisionAlgorithm,vtkSubdivisionAlgorithm);
    virtual void PrintSelf( ostream& os, vtkIndent indent );

    //BTX
    virtual void SetCell( const vtkShoeMeshIterator& c );
    vtkShoeMeshIterator& GetCell();
    const vtkShoeMeshIterator& GetCell() const;
    //ETX

		// Description:
		// Evaluate all of the fields that should be output with the
		// given \a vertex and store them just past the parametric coordinates
		// of \a vertex, at the offsets given by
		// \p vtkSubdivisionAlgorithm::GetFieldOffsets() plus \a field_start.
		// \a field_start contains the number of world-space coordinates (always 3)
		// plus the embedding dimension (the size of the parameter-space in which
		// the cell is embedded). It will range between 3 and 6, inclusive.
		//
		// You must have called SetCell() before calling this routine or there
		// will not be a mesh over which to evaluate the fields.
		//
		// You must have called \p vtkSubdivisionAlgorithm::PassDefaultFields()
		// or \p vtkSubdivisionAlgorithm::PassField() or there will be no fields
		// defined for the output vertex.
		//
		// This routine is public and returns its input argument so that it
		// may be used as an argument to
		// \p vtkAdaptiveTessellator::AdaptivelySamplekFacet():
		// @verbatim
		//   vtkAdaptiveTessellator* t = vtkAdaptiveTessellator::New();
		//   vtkSubdivisionAlgorithm* s;
		//   ...
		//   t->AdaptivelySample1Facet( s->EvaluateFields( p0 ), s->EvaluateFields( p1 ) );
		//   ...
		// @endverbatim
		// Although this will work, using \p EvaluateFields() in this manner
		// should be avoided. It's much more efficient to fetch the corner values
		// for each attribute and copy them into \a p0, \a p1, ... as opposed to
		// performing shape function evaluations.
		double* EvaluateFields( double* vertex, int field_start );

  protected:
    vtkShoeMeshSubdivisionAlgorithm() { }
    virtual ~vtkShoeMeshSubdivisionAlgorithm() { }

    vtkShoeMeshIterator CurrentCell;

  private:
    vtkShoeMeshSubdivisionAlgorithm( const vtkShoeMeshSubdivisionAlgorithm& ); // Not implemented.
    void operator = ( const vtkShoeMeshSubdivisionAlgorithm& ); // Not implemented.

};

//BTX

inline vtkShoeMeshIterator& vtkShoeMeshSubdivisionAlgorithm::GetCell() { return this->CurrentCell; }
inline const vtkShoeMeshIterator& vtkShoeMeshSubdivisionAlgorithm::GetCell() const { return this->CurrentCell; }

//ETX

#endif // __vtkShoeMeshSubdivisionAlgorithm_h
