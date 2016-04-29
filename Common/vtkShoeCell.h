// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
// .NAME vtkShoeCell - A class that effectively behaves like a mesh iterator by referencing a single cell in a mesh.
// .SECTION Description
//
// The new "generic" dataset API treats cells as iterators over a mesh, rather
// than as first-class datasets in their own right.
// Also, this one class (vtkShoeCell) represents <i>all</i> the possible
// cells in a mesh, where the old VTK API requires a different class
// for each cell shape (i.e., hexahedron, tetrahedron, triangle, ...).
//
// .SECTION Metadata
//
// This class uses function uses a helper class (vtkShoeCellMetaData) and
// function pointers to perform operations on cells with an efficiency that
// virtual functions and inheritance cannot.
//
// Inheritance is a bad model for cell types because many subroutines could be
// shared across cell types if static data was allowed to vary by cell type.
// For instance, let's say we have an algorithm that needs to access the
// list of connectivity entries that define a given face. If the algorithm
// only works for hexahedra, it's a simple matter... just look at
// ShoeHexahedron::FaceArray, which is a static class array.
// But what if the algorithm is general enough to work on any cell shape?
// We can't define a virtual static array in the parent class of ShoeHexahedron.
// If we insist on using inheritance, the only solution is to make a virtual
// function that returns the array.
// On its own, that would be OK, but now what if we have another member function
// that depends only on the interpolant family and product space of an element (for
// instance, a degree reduction operator). If we use multiple levels of
// inheritance, such as Cell->Hexahedron->TensorLagrangeHex, then our degree reduction
// routine would have to be copied once for each shape (Hex, Tet, Quad, Tri, ...).
// By keeping a "metadata" class, we can directly access "virtual static" arrays without
// and reuse the same "virtual" functions across all cell types.
// And it is very easy to write a script that creates metadata objects as needed.
//
// .SECTION See also
// vtkShoeCellMetaData vtkShoeBox

#ifndef __vtkShoeCell_h
#define __vtkShoeCell_h

#include "vtkGenericAdaptorCell.h"
#include "vtkShoeCellRecord.h" // this is not derived from vtkObject and contained by a ShoeCell
#include "vtksnlCommonWin32Header.h"

class vtkShoeCellMetaData;
class vtkShoeBox;
class vtkPolynomialSystem;
class vtkShoeAttribute;
class vtkDataRecords;

class VTK_SNL_COMMON_EXPORT vtkShoeCell : public vtkGenericAdaptorCell
{
public:
  vtkTypeRevisionMacro(vtkShoeCell,vtkGenericAdaptorCell);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkShoeCell* New();

  virtual int GetType();
  virtual int GetDimension();
  virtual int GetGeometryOrder();
  virtual int GetAttributeOrder(vtkGenericAttribute *a);
  //virtual int GetHighestOrderAttribute(vtkGenericAttributeCollection *ac);
  virtual int GetNumberOfPoints();
  virtual int GetNumberOfBoundaries(int dim=-1);
  virtual int GetNumberOfDOFNodes();

  virtual void Derivatives(int subId, double pcoords[3], vtkGenericAttribute *attribute, double *derivs);
  virtual void ParametricDerivatives( double pcoords[3], vtkShoeAttribute *attribute, double *derivs);
  virtual void GetSymbolicAttributeComponentGradient( vtkGenericAttribute* a, int c, vtkPolynomialSystem* s );
  virtual int GetParametricCenter(double pcoords[3]);
  virtual double GetParametricDistance(double pcoords[3]);
  virtual double* GetParametricCoords();
  virtual int GetNumberOfVerticesOnFace(int faceId);
  virtual int* GetFaceArray( int faceId );
  virtual int* GetEdgeArray( int edgeId );

  virtual vtkIdType GetId() { return this->CellID; }
  virtual int IsInDataSet();
  virtual int IsPrimary() { return 1; }

  virtual void GetPointIterator( vtkGenericPointIterator *it );
  virtual vtkGenericCellIterator* NewCellIterator();
  virtual void GetBoundaryIterator( vtkGenericCellIterator *boundaries, int dim=-1 );

  virtual int FindClosestBoundary( int subId, double pcoords[3], vtkGenericCellIterator* &boundary );

  virtual int CountNeighbors( vtkGenericAdaptorCell *boundary );
  virtual void CountEdgeNeighbors( int* sharing );
  virtual void GetNeighbors( vtkGenericAdaptorCell *boundary, vtkGenericCellIterator *neighbors );

  virtual int EvaluatePosition( double x[3], double *closestPoint, int &subId, double pcoords[3], double &dist2 );
  virtual void EvaluateLocation( int subId, double pcoords[3], double x[3] );
  virtual void InterpolateTuple( vtkGenericAttribute *a, double pcoords[3], double *val );
  virtual void InterpolateTuple( vtkGenericAttributeCollection *c, double pcoords[3], double *val );

#if 0
  virtual void Contour( vtkContourValues *values, vtkImplicitFunction *f, vtkGenericAttributeCollection *attributes,
                        vtkGenericCellTessellator *tess, vtkPointLocator *locator, vtkCellArray *verts,
                        vtkCellArray *lines, vtkCellArray *polys, vtkPointData *outPd, vtkCellData *outCd,
                        vtkPointData *internalPd, vtkPointData *secondaryPd, vtkCellData *secondaryCd );
  virtual void Clip( double value, vtkImplicitFunction *f, vtkGenericAttributeCollection *attributes,
                     vtkGenericCellTessellator *tess, int insideOut, vtkPointLocator *locator, 
                     vtkCellArray *connectivity, vtkPointData *outPd, vtkCellData *outCd, vtkPointData *internalPd,
                     vtkPointData *secondaryPd, vtkCellData *secondaryCd );
#endif
  virtual int IntersectWithLine( double p1[3], double p2[3], double tol, double &t,
                                 double x[3], double pcoords[3], int &subId );
  virtual void GetBounds( double bounds[6] );
#if 0
  virtual void Tessellate( vtkGenericAttributeCollection *attributes, vtkGenericCellTessellator *tess,
                           vtkPoints *points, vtkPointLocator *locator, vtkCellArray* cellArray,
                           vtkPointData *internalPd, vtkPointData *pd, vtkCellData* cd, vtkUnsignedCharArray *types);
#endif
  virtual int IsFaceOnBoundary( vtkIdType faceId );
  virtual int IsOnBoundary();
  virtual void GetPointIds( vtkIdType* id );
#if 0
  virtual void TriangulateFace( vtkGenericAttributeCollection *attributes, vtkGenericCellTessellator *tess, int index, 
                                vtkPoints *points, vtkPointLocator *locator, vtkCellArray *cellArray,
                                vtkPointData *internalPd, vtkPointData *pd, vtkCellData *cd );
#endif

  // Description:
  // Return the connectivity array for the cell
  vtkIdType* GetConnectivity();

  // Description:
  // Given a polynomial function defined over an element's parameter space (i.e.,
  // with variables in \f$\{ r, s, t \}\f$ ), compute the restriction of
  // that function to a boundary of the cell.
  // @param ps the system of equations defining the polynomial
  // @param bdy an offset into the cell's DOF connectivity describing the
  //            boundary of interest. This does not include points; 0 is
  //            the first edge of the cell and any value up to
  //            vtkShoeCell::GetNumberOfBoundaries(-1) is valid, although
  //            the highest allowed will simply return \a ps.
  vtkPolynomialSystem* RestrictToBoundary( vtkPolynomialSystem* ps, int bdy );

  // Description:
  // Compute critical points of an attribute interior to the cell, along with
  // any critical points on the restriction of the attribute to the boundaries
  // of the cell. The critical points are then associated with the DOF node
  // of the cell or the boundary, depending on whether a restriction of the
  // attribute was used. The resulting critical points are stored in a data
  // structure internal to vtkShoeAttribute. The results are used by
  // vtkShoeCell::TriangulateBoundary.
  virtual void ComputeDOFCriticalPoints( vtkShoeAttribute* att );

  const vtkShoeCellMetaData* GetMetaData() const { return this->Meta; }
  void CopyFrom( const vtkShoeCell* src );
  void SetDataSet( vtkShoeBox* );
  vtkShoeBox* GetDataSet() { return this->DataSet; }
  void InvalidateCache();

  // Description:
  // Return the permutation of the requested edge (face) of the current cell.
  // This should be a number between 0 and the number of edges (faces) minus 1.
  // The returned permutation bit(s) are the permutation of cell's coordinates
  // projected to the given edge (face) relative to the coordinate system used
  // for storage.
  int GetEdgePermutation( int e ) { return this->CellRecord.GetEdgePermutation( e ); }
  int GetFacePermutation( int f ) { return this->CellRecord.GetFacePermutation( f ); }

  // Description:
  // Please do not call this function... it exists only for testing purposes.
  void SetMetaData( const vtkShoeCellMetaData* m ) { this->Meta = m; }

  // Description:
  // Return a reference to the particulars of the cell's storage in the mesh.
  vtkShoeCellRecord& GetRecord() { return this->CellRecord; }
  const vtkShoeCellRecord& GetRecord() const { return this->CellRecord; }
  
protected:
  vtkShoeCell();
  virtual ~vtkShoeCell();

  //BTX
  friend class vtkShoeCellIterator;
  friend class vtkShoeBox;
  //ETX

  vtkShoeBox* DataSet;
  const vtkShoeCellMetaData* Meta;
  vtkIdType CellID;
  vtkShoeCellRecord CellRecord;

  //vtkDoubleArray* GeometryCache;
  vtkDoubleArray** AttributeCache;
  int* IsAttributeDirty;
  int AttributeCacheLength;

  // Description:
  // Working space for holding shape function (or derivatives) during InterpolateTuple, Derivatives, etc.
  // Warning! This does make the evaluate functions non-reentrant and thread UNsafe!
  double ShapeTmp[3000]; // FIXME: Eventually needs to be varying length.

  void UpdateCache( int );
  void GetPermutedAttributeValues( int at );

  // Description:
  // Place any roots inside the cell boundary into the given critical point storage, \a cpStorage.
  // @param perm is the permutation of the cell's coordinates relative to the storage.
  // @param conn is the connectivity array for the cell
  // @param connOffset is the offset into \a conn array specifying which DOF node \a cpStorage is associated with.
  //        Note that \a connOffset does <b>not</b> include the \a NumberOfPoints... that must be added before
  //        looking up a value in \a conn.
  // @param ps the polynomial system containing the critical points for the specified DOF node.
  void StoreRelevantRoots( vtkPolynomialSystem* ps, vtkIdType* conn, int connOffset, int perm, vtkDataRecords* cpStorage, vtkShoeAttribute* sa );

private:
  vtkShoeCell( const vtkShoeCell& ); // Not implemented.
  void operator =( const vtkShoeCell& ); // Not implemented.
};

#endif // __vtkShoeCell_h
