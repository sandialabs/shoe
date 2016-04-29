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
// .NAME vtkShoeMeshIterator - Sandia Higher Order Elements (SHOE) mesh interator for mixed-order nonlinear finite elements
// .SECTION Description
// vtkShoeMeshIterator is an iterator class for traversing a ShoeMesh

#ifndef vtkShoeMeshIterator_h
#define vtkShoeMeshIterator_h

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>

#include <vtkCellEnums.h>
#include <vtkShoeMesh.h>

class vtkDoubleArray;
struct vtkCellOps;

class VTK_SNL_COMMON_EXPORT vtkShoeMeshIterator
{
	public:
		// Description:
		// The default constructor. It constructs an invalid iterator.
		vtkShoeMeshIterator();

		// Description:
		// The copy constructor.
		vtkShoeMeshIterator( const vtkShoeMeshIterator& iter_in );

		//BTX
		// Description:
		// This constructor is intended to be used by vtkShoeMesh.begin() and vtkShoeMesh.end()
		// for setting the bounds of traversal.
		vtkShoeMeshIterator( vtkShoeMesh* mesh_in, vtkShoeMesh::CellsType::iterator init_cell_in );

		// Description:
		// This constructor does not take a cell definition -- used for MeshOrder traversal
		vtkShoeMeshIterator( vtkShoeMesh::const_iterator it_in, MeshTraversalStyle traversal_in );

		// Description:
		// These constructors are used by functions for traversing the mesh, perhaps selectively.
		// The second differs from the first only by allowing a mask on the target cell definition.
		vtkShoeMeshIterator( vtkShoeMesh::const_iterator it_in, MeshTraversalStyle traversal_in,
												 vtkShoeMesh::CellDefsType::const_iterator celldef_in );
		vtkShoeMeshIterator( vtkShoeMesh::const_iterator it, MeshTraversalStyle traversal_in,
												 vtkShoeMesh::CellDefsType::const_iterator celldef_in, int mask_in );

		// Description:
		// Deletes vtkDoubleArrays used for caching function and geometry values
		virtual ~vtkShoeMeshIterator();
		//ETX

		// Description:
		// <i>Do <b>not</b> use this member function</i>.
		// vtkShoeMeshIterator is a lightweight class that is not derived from vtkObject;
		// it is meant to be instantiated directly.
		// This member function exists only for the Tcl wrapper.
		static vtkShoeMeshIterator* New();

		//BTX
		// Description:
		// Accessor functions for various members
		vtkShoeMesh::CellsType::const_iterator GetCell() const;
		vtkShoeMesh::CellsType::iterator& GetCell();
		const vtkShoeMesh* GetMesh() const;
		vtkShoeMesh* GetMesh();
		//ETX
		
		// Description:
		// Return geometric and field parametric orders for current cell
		const int* GetCellGeometricOrder() const;
		const int* GetCellFunctionOrder ( int field_in ) const;

		// Description:
		// The mesh iterator caches cell geometry and function data for use by the calling function.
		// These accessors request that data for the current cell be cached, check to see if it has
		// been cached, retrieve the data for use by the calling function, and set the data.
		void CacheGeometry();
		bool IsGeometryCached() const;
		vtkDoubleArray* GetCachedGeometry();
		void SetCachedGeometry( vtkDoubleArray* geometry_in );
		void CacheFunction( int field_in );
		bool IsFunctionCached( int field_in ) const;
		vtkDoubleArray* GetCachedFunction( int field_in );
		void SetCachedFunction( int field_in, vtkDoubleArray* function_in );

		//BTX
		// Description:
		// A faster way to get the cell ops for the current cell.
		// Use this instead of GetCell()->Def->GetCellOps().
		const vtkCellOps* GetCellOps() const;
		//ETX

		// Description:
		// Get the ID of the currently referenced cell.
		vtkIdType GetCellId() const;

		//BTX
		// Description:
		// The pre- and post-increment operators account for different traversal styles
		// and masks on the target cell definition.
		vtkShoeMeshIterator& operator ++();
		vtkShoeMeshIterator operator ++( int );

		// Description:
		// Logical equality operator. Two iterators are equal if
		// they point to the same cell and same mesh.
		bool operator ==( const vtkShoeMeshIterator& it_in );
		bool operator !=( const vtkShoeMeshIterator& it_in );

		// Description:
		// The assignment operator.
		vtkShoeMeshIterator& operator = ( const vtkShoeMeshIterator& );
		//ETX
	protected:
		//BTX
		// Description:
		// Allocate memory for function cache and set other values that are not pass as arguments
		void CommonConstructor();

		// Description:
		// Current cell in the mesh and back pointer to the mesh
		vtkShoeMesh::CellsType::iterator CurCell;
		vtkShoeMesh* Mesh;

		// Description:
		// Structure containing all operations that apply to current cell
		const vtkCellOps* CurCellOps;

		// Description:
		// The cached data
		vtkDoubleArray* GeomCache;
		vtkDoubleArray** FuncCache;
		int FuncCacheLen;

		// Description:
		// Indicator for type of cell traversal
		// See values in vtkCellEnums.h for details
		MeshTraversalStyle TraversalStyle;

		// Description:
		// The target cell definition for this traversal
		// The mask indictes which fields in the definition
		// we are interested in
		struct TraversalSpec
		{
			int Mask;
			vtkShoeMesh::CellDefsType::const_iterator CellDef;
		};
		TraversalSpec TraversalDef;
		//ETX
};

//BTX
inline const vtkShoeMesh* vtkShoeMeshIterator::GetMesh() const { return this->Mesh; }
inline vtkShoeMesh* vtkShoeMeshIterator::GetMesh() { return this->Mesh; }

inline vtkShoeMesh::CellsType::const_iterator vtkShoeMeshIterator::GetCell() const { return this->CurCell; }
inline vtkShoeMesh::CellsType::iterator& vtkShoeMeshIterator::GetCell() { return this->CurCell; }

inline const vtkCellOps* vtkShoeMeshIterator::GetCellOps() const { return this->CurCellOps; }

inline const int* vtkShoeMeshIterator::GetCellGeometricOrder() const { return this->CurCell->Def->GetGeometricOrder(); }
inline bool vtkShoeMeshIterator::IsGeometryCached() const { return this->GeomCache; }
inline vtkDoubleArray* vtkShoeMeshIterator::GetCachedGeometry()
	{ if ( ! IsGeometryCached() ) CacheGeometry(); return this->GeomCache; }

inline const int* vtkShoeMeshIterator::GetCellFunctionOrder ( int field_in ) const { return this->CurCell->Def->GetFieldOrder( field_in ); }
inline bool vtkShoeMeshIterator::IsFunctionCached( int field_in ) const { return this->FuncCache[field_in]; }
inline vtkDoubleArray* vtkShoeMeshIterator::GetCachedFunction( int field_in )
	{ if ( ! IsFunctionCached( field_in ) ) CacheFunction( field_in ); return this->FuncCache[field_in]; }

inline vtkIdType vtkShoeMeshIterator::GetCellId() const { return this->CurCell.where(); }
//ETX

#endif // vtkShoeMeshIterator_h

