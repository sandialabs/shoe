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
#ifndef __vtkShoeMeshPolygonizer_h
#define __vtkShoeMeshPolygonizer_h

// .NAME vtkShoeMeshPolygonizer - approximate nonlinear FEM elements with simplices
// .SECTION Description
// This class approximates nonlinear FEM elements with vtkUnstructuredGrid.
// Only 1-D and 2-D Cells are currently handled, but adding a volume appoximator
// is a simple exercise and existing subdivision algorithms should
// work without modification.
//
// This filter uses a vtkShoeMeshIterator to rifle through all the cells in an
// input vtkShoeMesh. It calls each cell's Tessellate operation (defined in
// a vtkCellOps structure). The tessellator uses the vtkAdaptiveTessellator
// and vtkShoeMeshSubdivisionAlgorithm classes to generate simplices that approximate the
// nonlinear mesh using some approximation metric (encoded in the particular
// vtkShoeMeshSubdivisionAlgorithm::EvaluateEdge implementation). The simplices are placed into
// the filter's output vtkUnstructuredGrid object by the callback routines
// AddATriangle and AddALine, which are registered with the tessellator.
//
// The output mesh will have geometry and any fields specified as attributes
// in the input mesh's function data (see vtkShoeMesh::GetFunctionData()).
// The attribute's copy flags are honored, except for normals. Normals may
// be generated by calling the SetGenerateNormals() member before updating
// the pipeline. This is because analytical normals exist for nonlinear
// elements (or on their boundary faces for 3-D elements).
//
// .SECTION Internals
//
// The filter's main member function is Execute(). This function first calls
// SetupOutput() which allocates arrays and some temporary variables for the
// primitive callbacks (OutputTriangle and OutputLine which are called by
// AddATriangle and AddALine, respectively). The CurrentCell iterator is
// then initialized and each cell's Tessellate() operation is called, which
// results in one or more calls to OutputTriangle and OutputLine to add
// elements to the OutputMesh. Finally, Teardown() is called to free the
// filter's working space.
//
// .SECTION See Also
// vtkShoeMeshToUnstructuredGridFilter vtkShoeMesh vtkUnstructuredGrid

#include <vtksnlConfigure.h>
#include <vtksnlGraphicsWin32Header.h>
#include <vtkShoeMeshIterator.h>
#include <vtkShoeMeshToUnstructuredGridFilter.h>

class vtkDataArray;
class vtkUnstructuredGrid;

class vtkAdaptiveTessellator;
class vtkSubdivisionAlgorithm;
class vtkShoeMeshSubdivisionAlgorithm;

class VTK_SNL_GRAPHICS_EXPORT vtkShoeMeshPolygonizer
  : public vtkShoeMeshToUnstructuredGridFilter
{
  public:
    vtkTypeRevisionMacro(vtkShoeMeshPolygonizer,vtkShoeMeshToUnstructuredGridFilter);
    void PrintSelf( ostream& os, vtkIndent indent );

    static vtkShoeMeshPolygonizer* New();

    virtual void SetTessellator( vtkAdaptiveTessellator* );
    vtkAdaptiveTessellator* GetTessellator();
    const vtkAdaptiveTessellator* GetTessellator() const;

    virtual void SetSubdivider( vtkShoeMeshSubdivisionAlgorithm* );
    vtkShoeMeshSubdivisionAlgorithm* GetSubdivider();
    const vtkShoeMeshSubdivisionAlgorithm* GetSubdivider() const;

    virtual void SetGenerateNormals( bool gen_normals_in );
    bool GetGenerateNormals() const;

    // Description:
    // This object references a tessellator and subdivider. When either is modified,
    // we consider the polygonizer to be modified (since it will produce different
    // output).
    virtual unsigned long GetMTime();

  protected:
    vtkShoeMeshPolygonizer();
    ~vtkShoeMeshPolygonizer();

    // Description:
    // Called by Execute to set up a multitude of member variables used by the per-primitive
    // output functions (OutputLine, OutputTriangle, and maybe one day... OutputTetrahedron).
    void SetupOutput();

    // Description:
    // Run the filter; produce a polygonal approximation to the grid.
    virtual void Execute();

    // Description:
    // Reset the temporary variables used during the filter's Execute() method.
    void Teardown();

    //BTX
    vtkAdaptiveTessellator* Tessellator;
    vtkShoeMeshSubdivisionAlgorithm* Subdivider;

    bool GenerateNormals;

    // Description:
    // These member variables are set by SetupOutput for use inside the
    // callback members OutputLine and OutputTriangle.
    vtkShoeMeshIterator CurrentCell;
    vtkUnstructuredGrid* OutputMesh;
    vtkPoints* OutputPoints;
    vtkDataArray** OutputAttributes;
    int* OutputAttributeIndices;

    static void AddALine( const double*, const double*, vtkSubdivisionAlgorithm*, void*, const void* );
    static void AddATriangle( const double*, const double*, const double*, vtkSubdivisionAlgorithm*, void*, const void* );
    static void AddATetrahedron( const double*, const double*, const double*, const double*, vtkSubdivisionAlgorithm*, void*, const void* );
    void OutputLine( const double*, const double* );
    void OutputTriangle( const double*, const double*, const double* );
    void OutputTetrahedron( const double*, const double*, const double*, const double* );
    //ETX

  private:
    vtkShoeMeshPolygonizer( const vtkShoeMeshPolygonizer& ); // Not implemented.
    void operator = ( const vtkShoeMeshPolygonizer& ); // Not implemented.
};

//BTX
inline vtkAdaptiveTessellator* vtkShoeMeshPolygonizer::GetTessellator() { return this->Tessellator; }
inline const vtkAdaptiveTessellator* vtkShoeMeshPolygonizer::GetTessellator() const { return this->Tessellator; }

inline vtkShoeMeshSubdivisionAlgorithm* vtkShoeMeshPolygonizer::GetSubdivider() { return this->Subdivider; }
inline const vtkShoeMeshSubdivisionAlgorithm* vtkShoeMeshPolygonizer::GetSubdivider() const { return this->Subdivider; }

inline bool vtkShoeMeshPolygonizer::GetGenerateNormals() const { return this->GenerateNormals; }
//ETX

#endif // __vtkShoeMeshPolygonizer_h
