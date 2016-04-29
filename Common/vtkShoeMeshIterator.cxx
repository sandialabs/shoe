/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <functional>
#include <algorithm>

#include <vtkDoubleArray.h>

#include <vtkShoeMeshIterator.h>
#include <vtkShoeMesh.h>
#include <vtkCellOps.h>



vtkShoeMeshIterator* vtkShoeMeshIterator::New()
{
	return new vtkShoeMeshIterator;
}

void vtkShoeMeshIterator::CommonConstructor()
{
	this->FuncCacheLen = this->Mesh->GetNumberOfNonlinearFunctions();

	this->FuncCache = new vtkDoubleArray*[this->FuncCacheLen];
	vtkstd::fill( this->FuncCache, this->FuncCache+this->FuncCacheLen, (vtkDoubleArray*)0 );
	GeomCache = 0;
}

vtkShoeMeshIterator::vtkShoeMeshIterator()
{
	this->Mesh = 0;
	this->CurCellOps = 0;
	this->TraversalStyle = MeshOrder;
	this->GeomCache = 0;
	this->FuncCache = 0;
	this->FuncCacheLen = 0;
}

vtkShoeMeshIterator::vtkShoeMeshIterator( const vtkShoeMeshIterator& iter_in )
{
	this->Mesh = iter_in.Mesh;
	this->CurCell = iter_in.CurCell;
	this->CurCellOps = iter_in.CurCellOps;
	this->TraversalStyle = iter_in.TraversalStyle;
	this->TraversalDef = iter_in.TraversalDef;
	this->CommonConstructor();

	// OK, the FuncCache has been allocated, copy over the caches
	this->FuncCacheLen = iter_in.FuncCacheLen;
	for ( int i=0; i< this->FuncCacheLen; i++ )
	{
		if ( iter_in.IsFunctionCached( i ) )
			this->SetCachedFunction( i, iter_in.FuncCache[i] );
	}
	if ( iter_in.IsGeometryCached() )
		this->SetCachedGeometry( iter_in.GeomCache );
}

vtkShoeMeshIterator::vtkShoeMeshIterator& vtkShoeMeshIterator::operator = ( const vtkShoeMeshIterator& iter_in )
{
	if ( this == &iter_in )
		return *this;

	this->Mesh = iter_in.Mesh;
	this->CurCell = iter_in.CurCell;
	this->CurCellOps = iter_in.CurCellOps;
	this->TraversalStyle = iter_in.TraversalStyle;
	this->TraversalDef = iter_in.TraversalDef;
	this->CommonConstructor();

	// OK, the FuncCache has been allocated, copy over the caches
	this->FuncCacheLen = iter_in.FuncCacheLen;
	for ( int i=0; i<this->FuncCacheLen; i++ )
	{
		if ( iter_in.IsFunctionCached( i ) )
			this->SetCachedFunction( i, iter_in.FuncCache[i] );
	}
	if ( iter_in.IsGeometryCached() )
		this->SetCachedGeometry( iter_in.GeomCache );

	return *this;
}

vtkShoeMeshIterator::vtkShoeMeshIterator( vtkShoeMesh* mesh, vtkShoeMesh::CellsType::iterator init_cell ) 
{
	this->Mesh = mesh;
	this->CurCell = init_cell;
	this->TraversalStyle = MeshOrder;
	this->CommonConstructor();
}

vtkShoeMeshIterator::vtkShoeMeshIterator( vtkShoeMesh::const_iterator it, MeshTraversalStyle traversal )
{
	this->Mesh = it.Mesh;
	this->CurCell = it.CurCell;
	this->TraversalStyle = traversal;
	this->TraversalDef.Mask = DefnMask;
	this->CurCellOps = this->CurCell->Def->GetCellOps();
	this->CommonConstructor();
}

vtkShoeMeshIterator::vtkShoeMeshIterator( vtkShoeMesh::const_iterator it, MeshTraversalStyle traversal,
														vtkShoeMesh::CellDefsType::const_iterator celldef )
{
	this->Mesh = it.Mesh;
	this->CurCell = it.CurCell;
	this->TraversalStyle = traversal;
	this->TraversalDef.CellDef = celldef;
	this->TraversalDef.Mask = DefnMask;
	this->CurCellOps = this->CurCell->Def->GetCellOps();
	this->CommonConstructor();
}

vtkShoeMeshIterator::vtkShoeMeshIterator( vtkShoeMesh::const_iterator it, MeshTraversalStyle traversal,
														vtkShoeMesh::CellDefsType::const_iterator celldef, int mask )
{
	this->Mesh = it.Mesh;
	this->CurCell = it.CurCell;
	this->TraversalStyle = traversal;
	this->TraversalDef.CellDef = celldef;
	this->TraversalDef.Mask = mask;
	this->CurCellOps = this->CurCell->Def->GetCellOps();
	this->CommonConstructor();
}

vtkShoeMeshIterator::~vtkShoeMeshIterator()
{
	if ( this->Mesh )
	{
		for (int i = 0; i < this->FuncCacheLen; i++)
			this->SetCachedFunction( i, 0 );
		delete [] this->FuncCache;
		this->SetCachedGeometry( 0 );
	}
}

void vtkShoeMeshIterator::SetCachedGeometry( vtkDoubleArray* geometry )
{ 
	if (this->GeomCache)
		this->GeomCache->UnRegister( 0 );
	
	this->GeomCache = geometry;

	if (this->GeomCache)
		this->GeomCache->Register( 0 );
}

void vtkShoeMeshIterator::SetCachedFunction( int field_in, vtkDoubleArray* function )
{ 
	if (this->FuncCache[field_in])
		this->FuncCache[field_in]->UnRegister( 0 );
	
	this->FuncCache[field_in] = function;

	if (this->FuncCache[field_in])
		this->FuncCache[field_in]->Register( 0 );
}

void vtkShoeMeshIterator::CacheFunction ( int field )
{
	if (!this->IsFunctionCached( field ))
			this->CurCellOps->GetPermutedFunctionValues( this, field );
}

void vtkShoeMeshIterator::CacheGeometry ()
{
	if (!this->IsGeometryCached())
			this->CurCellOps->GetPermutedGeometryValues( this );
}

struct CompareCellDefsMasked 
			: public vtkstd::binary_function <vtkShoeMesh::CellDefsType::const_iterator,vtkShoeMesh::CellDefsType::const_iterator,bool>
{
	MeshTraversalStyle Traversal;
	int Mask;
	int NumberOfFunctions;

	CompareCellDefsMasked( MeshTraversalStyle t, int m, int n ) { Traversal = t; Mask = m; NumberOfFunctions = n; }

	bool operator () ( const vtkShoeMesh::CellSpec cell1, vtkShoeMesh::CellDefsType::const_iterator def2 )
	{

		vtkShoeMesh::CellDefsType::const_iterator def1 = cell1.Def;
		switch (Traversal)
		{
			case OfSameType:
				return (def1->GetType() == def2->GetType());
				break;

			case OfSameDefn:
				return (def1 == def2);
				break;

			case CustomOrder:
				if ((Mask & int( ShapeMask)) == ShapeMask)
					return (def1->GetShape() == def2->GetShape());
				else if ((Mask & int( InterpMask )) == InterpMask)
					return (def1->GetInterpolant() == def2->GetInterpolant());
				else if ((Mask & int( ProdSpaceMask)) == ProdSpaceMask)
					return (def1->GetProductSpace() == def2->GetProductSpace());
				else if ((Mask & int( TypeMask )) == TypeMask)
					return def1->GetType() == def2->GetType();
				else if ((Mask & int( FunctionOrderMask)) == FunctionOrderMask)
				{
					for (int i = 0; i < NumberOfFunctions; i++)
						if (!def1->IsFieldOrder( i,  def2->GetFieldOrder( i )))
							return false;
						return true;
				}
				else if ((Mask & int( GeomOrderMask)) == GeomOrderMask)
					return (!def1->IsGeometricOrder( def2->GetGeometricOrder()));
				else if ((Mask & int( DefnMask )) == DefnMask)
					return def1 == def2;
				break;

			case MeshOrder:
				return true;
				break;
		}
    return false; // Things we don't know about can't be the same?
	}
};

vtkShoeMeshIterator& vtkShoeMeshIterator::operator ++()
{

	vtkShoeMesh::CellsType::iterator NextCell = ++(this->CurCell);
	CompareCellDefsMasked pred( this->TraversalStyle, this->TraversalDef.Mask, this->Mesh->GetNumberOfNonlinearFunctions() );

	//this->CurCell = find_if( NextCell, this->Mesh->EndCells(), bind2nd( pred, this->TraversalDef.CellDef) );
	for (; NextCell != this->Mesh->EndCells(); ++NextCell)
	{
		if (pred(*NextCell, this->TraversalDef.CellDef) )
			break;
	}

	if (NextCell != this->Mesh->EndCells())
	{
		this->CurCell = NextCell;
		this->CurCellOps = this->CurCell->Def->GetCellOps();
	}
	else
	{
		this->CurCell = NextCell;
		this->CurCellOps = 0;
	}
	this->SetCachedGeometry(0);
	for ( int i=0; i<this->Mesh->GetNumberOfNonlinearFunctions(); i++ )
		this->SetCachedFunction( i, 0 );

	return *this;
}

vtkShoeMeshIterator vtkShoeMeshIterator::operator ++( int )
{

	vtkShoeMeshIterator dummy( *this );
	++(*this);
	return dummy;
}

bool vtkShoeMeshIterator::operator ==( const vtkShoeMeshIterator& it)
{
	return (this->CurCell == it.CurCell && this->Mesh == it.Mesh);
}

bool vtkShoeMeshIterator::operator !=( const vtkShoeMeshIterator& it )
{
	return !(*this == it);
}

