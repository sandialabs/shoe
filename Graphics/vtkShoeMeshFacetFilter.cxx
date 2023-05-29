/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <iostream>
#include <iterator>

#include <vtkObjectFactory.h>
#include <vtkPointData.h>

#include <vtkShoeMesh.h>
#include <vtkShoeMeshIterator.h>
#include <vtkCellOps.h>
#include <vtkShoeMeshFacetFilter.h>



vtkStandardNewMacro(vtkShoeMeshFacetFilter);
vtkCxxRevisionMacro(vtkShoeMeshFacetFilter,"$Revision: 997 $");

vtkShoeMeshFacetFilter::vtkShoeMeshFacetFilter()
{
	AverageCellData = 0;
	BoundaryElementsOnly = 1;
	FacetDimension = 3;
}

vtkShoeMeshFacetFilter::~vtkShoeMeshFacetFilter()
{
}

void vtkShoeMeshFacetFilter::Execute()
{
	vtkDebugMacro("Execute");
	vtkDataSet* input = this->GetInput();
	vtkDataSet* output = this->GetOutput();

	vtkShoeMesh* src = vtkShoeMesh::SafeDownCast(input);
	vtkShoeMesh* dst = vtkShoeMesh::SafeDownCast(output);

	if ( !src || !dst )
		vtkErrorMacro( "This filter's input and ouput must both be vtkShoeMesh objects!" );

	// Clear out the destination
	dst->Reset();

	// Copy the points and point field data to the destination
	dst->SetPoints( src->GetPoints() );
	dst->SetFieldData( src->GetFieldData() );
	dst->GetPointData()->ShallowCopy( src->GetPointData() );
	dst->GetGeometryData()->ShallowCopy( src->GetGeometryData() );
	if ( src->GetFunctionData() )
		{
		if ( ! dst->GetFunctionData() )
			dst->SetFunctionData( vtkFunctionData::New() );
		dst->GetFunctionData()->ShallowCopy( src->GetFunctionData() );
		}
	else
		{
		dst->SetFunctionData( 0 );
		}

	switch ( FacetDimension )
		{
		case 0:
			vtkErrorMacro( "Semantically impossible to find bounding facets of a point set. Go away." );
			break;
		case 1:
			ComputePointsOfEdges( src, dst );
			break;
		case 2:
			ComputeEdgesOfFaces( src, dst );
			break;
		case 3:
			ComputeFacesOfVolumes( src, dst );
			break;
		default:
			vtkErrorMacro( "Invalid facet dimension." );
		}
}

void vtkShoeMeshFacetFilter::ComputePointsOfEdges( vtkShoeMesh* src, vtkShoeMesh* dst )
{
	vtkErrorMacro( "Bounding points are unimplemented. Why? Because they're useless! You're screwier than a wing nut!" );
}

void vtkShoeMeshFacetFilter::ComputeEdgesOfFaces( vtkShoeMesh* src, vtkShoeMesh* dst )
{
	vtkIdType bdyConn[27]; // No cell connectivity will be longer than this

	dst->GetDofNodeLinks().resize( src->GetNumberOfDofNodes() );
	for ( vtkShoeMeshIterator it( src->Begin(), MeshOrder ); it != src->End(); ++it )
		{
		const vtkCellOps* ops = it.GetCellOps();
		vtkIdType itId = it.GetCellId();
		vtkIdType id;
		vtkShoeMesh::CellsType::const_iterator cellSpec= it.GetCell();
		if ( ops->EmbeddingDimension < 2 )
			{
			id = dst->InsertNextCell( dst->FindOrCreateCellDef( cellSpec->Def ),
			                          it.GetMesh()->GetCellConnectivity( itId ),
			                          cellSpec->GetCellPermutation() );
			}
		else
			{
			vtkCellDefinition bdyDef;
			uint32_t bdyPerm;
			for ( int e=0; e<ops->NumberOfEdges; e++ )
				{
				// we can do this next step because src (and thus this iterator) and dst share DOF node ids
				ops->GetBoundaryEdge( bdyDef, bdyConn, bdyPerm, it, e );
				vtkIdType dof = ops->GetBoundaryEdgeDofNodeId( src, itId, e );
				if ( BoundaryElementsOnly && (src->GetDofNodeLinks()[ dof ].size() > 2) )
					continue;
				if ( dst->GetDofNodeLinks()[ dof].empty() )
					id = dst->InsertNextCell( dst->FindOrCreateCellDef( bdyDef ), bdyConn, bdyPerm );
				} // for each boundary edge of the cell
			} // check cell dimensionality < 2
		} // for each cell in src
}

void vtkShoeMeshFacetFilter::ComputeFacesOfVolumes( vtkShoeMesh* src, vtkShoeMesh* dst )
{
	vtkIdType bdyConn[27]; // No cell connectivity will be longer than this

	dst->GetDofNodeLinks().resize( src->GetNumberOfDofNodes() );
	for ( vtkShoeMeshIterator it( src->Begin(), MeshOrder ); it != src->End(); ++it )
		{
		const vtkCellOps* ops = it.GetCellOps();
		vtkIdType itId = it.GetCellId();
		vtkIdType id;
		vtkShoeMesh::CellsType::const_iterator cellSpec= it.GetCell();
		if ( ops->EmbeddingDimension < 3 )
			{
			id = dst->InsertNextCell( dst->FindOrCreateCellDef( cellSpec->Def ),
			                          it.GetMesh()->GetCellConnectivity( itId ),
			                          cellSpec->GetCellPermutation() );
			}
		else
			{
			vtkCellDefinition bdyDef;
			uint32_t bdyPerm;
			for ( int f=0; f<ops->NumberOfFaces; f++ )
				{
				// we can do this next step because src (and thus this iterator) and dst share DOF node ids
				ops->GetBoundaryFace( bdyDef, bdyConn, bdyPerm, it, f );
				vtkIdType dof = ops->GetBoundaryFaceDofNodeId( src, itId, f );
				if ( BoundaryElementsOnly && (src->GetDofNodeLinks()[ dof ].size() > 1) )
					continue;
				if ( dst->GetDofNodeLinks()[ dof].empty() )
					id = dst->InsertNextCell( dst->FindOrCreateCellDef( bdyDef ), bdyConn, bdyPerm );
				} // for each boundary face of the cell
			} // check cell dimensionality < 3
		} // for each cell in src
}

