// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
// .NAME vtkShoeBox - Represent a higher order finite element mesh.
// .SECTION Description
// SHOE (Sandia Higher Order Elements) is a set of classes for representing
// piecewise polynomial functions over arbitrary regions of space.
// These polynomials are typically used by finite element analysis codes to
// approximate the solution to a partial differential equation.
//
// This class contains connectivity information that specifies how different
// finite elements (<em>cells</em> in VTK's terms) share coefficients.
// The actual coefficients are store in vtkShoeAttribute classes -- one per
// attribute of the mesh (e.g., pressure, temperature, velocity).
// The geometry of the mesh is treated just like any other attribute, but
// each vtkShoeBox keeps a separate pointer to it so that it can be referenced
// by its function as a geometric interpolant.
//
// .SECTION See Also
// vtkShoeAttribute vtkShoeCell vtkShoeCellIterator vtkShoePointIterator vtkShoeCellGenus vtkShoeCellSpecies
#ifndef __vtkShoeBox_h
#define __vtkShoeBox_h

#include <vtksnlConfigure.h> // for uint32_t
#include <vtkGenericDataSet.h>

#define VTK_SHOE_MESH 42

class vtkShoeBoxP;
class vtkDataRecords;
class vtkShoeAttribute;
class vtkGenericAttribute;
class vtkGenericAttributeCollection;
class vtkShoeCellSpecies;
class vtkShoeBoxPartition;
class vtkShoeCell;

class VTK_EXPORT vtkShoeBox : public vtkGenericDataSet
{
public:
  vtkTypeRevisionMacro(vtkShoeBox,vtkGenericDataSet);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkShoeBox *New();

  virtual vtkIdType GetNumberOfPoints();
  virtual vtkIdType GetNumberOfDOFNodes();
  virtual vtkIdType GetNumberOfCells( int dim );

  virtual int GetCellDimension() { return -1; } /* FIXME: Eventually, this might also return 0,1,2, or 3 */
  virtual void GetCellTypes( vtkCellTypes* types );

  virtual vtkGenericCellIterator* NewCellIterator( int dim=-1 );
  virtual vtkGenericCellIterator* NewBoundaryIterator( int dim=-1, int exteriorOnly=0 );
  virtual vtkGenericPointIterator* NewPointIterator();

  virtual int FindCell( double x[3], vtkGenericCellIterator*& cell, double tol2, int &subId, double pcoords[3] );
  virtual void FindPoint( double x[3], vtkGenericPointIterator *p );

  //virtual unsigned long int GetMTime();

  virtual void ComputeBounds();
  //virtual double* GetBounds();
  //virtual void GetBounds(double bounds[6]);

  //virtual double *GetCenter();
  //virtual void GetCenter(double center[3]);

  //virtual double GetLength();

  int GetDataObjectType() { return VTK_SHOE_MESH; }

  virtual unsigned long GetActualMemorySize();
  virtual vtkIdType GetEstimatedSize();

  // Description:
  // Returns whether the mesh is a "temporary" mesh, e.g., one that
  // is owned by an iterator for the purposes of enumerating cells
  // on the boundary of a cell.
  //
  // <b>Design Note:</b>
  // It may be enough to have one temporary mesh per actual mesh;
  // by building a smart iterator (that can look over an arbitrary
  // subset of a mesh), all the iterators could share a single
  // "scratch" mesh to hold boundary faces/edges as required.
  int IsTemporaryMesh() const { return this == this->TemporaryMesh; }

  // Description:
  // Return an integer identifying the given attribute.
  int AttributeId( vtkGenericAttribute* );

  // Description:
  // Set the attribute containing the geometric interpolant.
  // The ID of the attribute is returned.
  int SetGeometry( vtkShoeAttribute* );
  vtkShoeAttribute* GetGeometry() { return this->Geometry; }

  // Description:
  // Insert a cell into the mesh given its connectivity, the type of cell, and a permutation,
  // The ID of the cell is returned.
  // The permutation is a 4-byte bit-vector specifying the orientation of each DOF node's
  // coordinate system relative to the cell's coordinate system.
  // Edge DOFs require a single bit and are always stored in the 15 least significant bits (0-14).
  // Any face DOFs present require 2 bits and are always stored in bits 15-28.
  // If a volumetric DOF is present, it requires 3 bits and is stored in the most significant 3 bits (29-31).
  // If a bit is unset, then the coordinate to which it corresponds is unpermuted relative to the
  // cell's coordinate system. Otherwise, its orientation is reversed.
  //
  // You are responsible for insuring that the cell species has an
  // appropriate interpolant, product space, and order defined for each
  // attribute in the mesh before any interpolation is performed.
  vtkIdType InsertNextCell( int species, vtkIdType* connectivity, uint32_t permutation );

  // Description:
  // If required, rebuild the PointLinks and DOFLinks structures used for reverse lookup of cells.
  // This will run even if MaintainLinks is set to Off.
  virtual void UpdateLinks();

  int GetMaintainLinks() const;
  void SetMaintainLinks( int );
  vtkBooleanMacro(MaintainLinks,int);

  // Description:
  // Return the ID of the first cell in the mesh (which may not have an ID of 0).
  virtual vtkIdType GetFirstCellId();

  // Description:
  // Return the cell at a given ID. This will return NULL if no cell exists for the given ID.
  // You are responsible for deleting the cell.
  vtkShoeCell* GetCell( vtkIdType );

  // Description:
  // Return the ID of the first 2-D DOF node (face node) that is on the boundary of the mesh.
  virtual vtkIdType GetFirstBoundaryDOF();

  // Description:
  // Return the array holding the mesh connectivity.
  vtkDataRecords* GetConnectivity() { return this->Connectivity; }

  // Description:
  // Return the reverse lookup tables for identifying which cells refer to
  // a given point or DOF node.
  vtkDataRecords* GetPointLinks() { return this->PointLinks; }
  vtkDataRecords* GetDOFLinks() { return this->DOFLinks; }

  // Description:
  // Construct a \f$\kappa\f$-compatible tessellation.
  vtkShoeBoxPartition* PartitionMesh( vtkGenericAttributeCollection* kappa );

  virtual void Reset();

protected:
  vtkShoeBox();
  virtual ~vtkShoeBox();

  //BTX
  friend class vtkShoeCell;
  friend class vtkShoeCellIterator;
  friend class vtkShoeBoxReader;
  //ETX

  // Description:
  // This is used by iterators to prepare a mesh for boundary traversal.
  // It creates a new child mesh if this mesh is not already temporary (otherwise TemporaryMesh = this)
  // and does a shallow copy of the mesh's attributes.
  void CreateTemporaryMesh();

  // Description:
  // Loop over a cell's connectivity and add the cell ID to
  // each corresponding CornerLinks and DOFLinks entry.
  void AddLinksForCell( vtkIdType cellId );

  // Description:
  // Loop over a cell's connectivity and remove the cell ID
  // from each corresponding CornerLinks and DOFLinks entry.
  void DeleteLinksForCell( vtkIdType cellId );

  // Description:
  // Add (delete) a single ID to (from) a Link array
  // This is called by AddLinksForCell (DeleteLinksForCell).
  // DeleteLink will fail silently if the cell wasn't present in the Link's record.
  void AddLink( vtkDataRecords* links, vtkIdType id, vtkIdType cell );
  void DeleteLink( vtkDataRecords* links, vtkIdType id, vtkIdType cell );

  // Description:
  // Compute critical points for an attribute of the mesh.
  //
  // Note that this is an interesting problem, because multiple
  // meshes may reference the same attribute, but each mesh may
  // use a different subset of the DOF nodes in the mesh.
  // Thus, to be entirely certain that <em>every</em> DOF node's
  // critical points are computed, vtkShoeBox::ComputeCriticalPoints
  // must be called on <em>every</em> one of those meshes.
  // That is why this function is protected; it should only be
  // called by vtkShoeBox::PartitionMesh and
  // vtkShoeAttribute::GetRange.
  //
  // As critical points are computed for each DOF node of the attribute,
  // a bit array stored in the attribute is marked. Only when <em>all</em>
  // of the DOF nodes have been marked will the attribute's
  // CriticalPointsDirty flag be cleared.
  void ComputeCriticalPoints( vtkShoeAttribute* att );

  // Description:
  // Compute critical points for a collection of attributes defined over a mesh.
  // This is a convenience routine that calls
  // vtkShoeBox::ComputeCriticalPoints(vtkShoeAttribute*) for each attribute.
  void ComputeCriticalPoints( vtkGenericAttributeCollection* kappa );

  // Description:
  // Triangulate the boundaries of cells using critical point information.
  // This is called from vtkShoeBox::PartitionMesh() in order to create
  // a \f$\kappa\f$-compatible partition of the mesh.
  void TriangulateBoundaries( vtkGenericAttributeCollection* kappa, vtkShoeBoxPartition* bdyTri );

  // Description:
  // Tetrahedralize the interiors of cells using critical point information.
  // This is called from vtkShoeBox::PartitionMesh() in order to create
  // a \f$\kappa\f$-compatible partition of the mesh.
  void TetrahedralizeInteriors( vtkGenericAttributeCollection* kappa, vtkShoeBoxPartition* bdyTri, vtkShoeBoxPartition* tetInt );

  // Description:
  // Split tetrahedra in the partition where interior faces or edges have critical points.
  // This is called from vtkShoeBox::PartitionMesh() in order to create
  // a \f$\kappa\f$-compatible partition of the mesh.
  void CorrectPartitionTopology( vtkGenericAttributeCollection* kappa, vtkShoeBoxPartition* tetInt );

  vtkShoeBox* TemporaryMesh;

  int GeometryIndex;
  vtkShoeAttribute* Geometry;
  vtkDataRecords* Connectivity;
  vtkIdType NumberOfCellsOfDim[4]; // How many cells of a given dimension are present?

  // The cells are stored using templated classes,
  // so a private class is used.
  vtkShoeBoxP* CellRecs;

  vtkDataRecords* PointLinks;
  vtkDataRecords* DOFLinks;
  int LinkState;

private:
  vtkShoeBox(const vtkShoeBox&); // Not implemented.
  void operator=(const vtkShoeBox&); // Not implemented.

  static int ElementsRegistered;
  static void RegisterElements();
};

//BTX
inline int vtkShoeBox::GetMaintainLinks() const { return this->LinkState & 1 ; }
inline void vtkShoeBox::SetMaintainLinks( int b )
{
  if ( !(b ^ (this->LinkState & 1)) )
    {
    return;
    }
  this->LinkState = this->LinkState & (~1) | (b != 0);
  this->Modified();
}
//ETX
#endif // __vtkShoeBox_h
