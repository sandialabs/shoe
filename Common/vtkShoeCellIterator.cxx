// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeCellIterator.h"
#include "vtkShoeBoxP.h"
#include "vtkShoeBox.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellMetaData.h"
#include "vtkShoeCellRecord.h"
#include "vtkShoeCellSpecies.h"
#include "vtkDataRecordsIterator.h"

#include <vtkObjectFactory.h>

vtkStandardNewMacro(vtkShoeCellIterator);
vtkCxxRevisionMacro(vtkShoeCellIterator,"$Revision: 9531 $");

class vtkShoeCellIteratorP
{
public:
  vtkShoeCellIteratorP() { }
  ~vtkShoeCellIteratorP() { }

  freelist<vtkShoeCellRecord,vtkIdType>::iterator Cell;
};

void vtkShoeCellIterator::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << "DataSet: " << this->DataSet << vtkstd::endl
     << "Iter: " << this->Iter << vtkstd::endl
     << "CurrentCell: " << this->CurrentCell << vtkstd::endl
     << "CurrentParent: " << this->CurrentParent << vtkstd::endl
     << "CurrentFacet: " << this->CurrentFacet << vtkstd::endl
     << "Dimension: " << this->Dimension << vtkstd::endl
     << "TraversalStyle: " << this->TraversalStyle << vtkstd::endl
     << "TraversalMask: " << this->TraversalMask << vtkstd::endl;
}

vtkShoeCellIterator::vtkShoeCellIterator()
{
  this->DataSet = 0;
  this->CurrentCell = vtkShoeCell::New();
  this->Iter = new vtkShoeCellIteratorP;
}

vtkShoeCellIterator::~vtkShoeCellIterator()
{
  this->SetDataSet(0);
  this->CurrentCell->Delete();
  delete Iter;
}

void vtkShoeCellIterator::SetDataSet( vtkShoeBox* ds )
{
  if ( ds == this->DataSet )
    {
    return;
    }
  if ( this->DataSet )
    {
    this->DataSet->UnRegister( this );
    }
  this->DataSet = ds;
  this->CurrentCell->SetDataSet( ds );
  this->Modified();
  if ( this->DataSet )
    {
    this->Iter->Cell = this->DataSet->CellRecs->Cells.begin();
    this->DataSet->Register( this );
    }
}

void vtkShoeCellIterator::Begin()
{
  // FIXME. This should eventually handle mix-n-match iterators such as
  // "iterate over boundaries of the selected cells if their order is (2,3,3)"
  assert( this->DataSet );
  switch ( this->TraversalStyle )
    {
    case OfSameType:  //!< Make multiple passes through the mesh, covering only cells of a single type per pass
    case OfSameDefn:  //!< Traverse cells, visited only those with same CellDefinition (possible with a mask)
    case CustomOrder: //!< User defined traversal using MeshTraversalMask defined below
    case Selection: //!< Iterate over the selected cells
    case CellBoundary: //!< Traverse the boundary of a single cell.
      vtkWarningMacro("Can't traverse cells this way yet. Falling back to MeshOrder.");
    case MeshOrder:   //!< Traverse cells one by one as they occur in the mesh
      this->Iter->Cell = this->DataSet->CellRecs->Cells.begin();
      if ( this->Iter->Cell != this->DataSet->CellRecs->Cells.end() )
        this->SetCellById( this->Iter->Cell.where() );
      break;
    case Boundary:
      this->DataSet->UpdateLinks(); // need to find cells referring to DOF nodes used to iterate bdy
      this->DataSet->CreateTemporaryMesh(); // need to hold connectivity for bdy facets
      this->Iter->Cell = this->DataSet->CellRecs->Cells.begin();
      this->SetCellById( this->Iter->Cell.where() );
      break;
    case ExtBoundary:
      this->DataSet->UpdateLinks();
      this->DataSet->CreateTemporaryMesh();
      this->SetBoundaryCellByDOF( this->DataSet->GetFirstBoundaryDOF() );
      break;
    };
}

vtkGenericAdaptorCell* vtkShoeCellIterator::NewCell()
{
  return 0;
}

void vtkShoeCellIterator::GetCell( vtkGenericAdaptorCell* c )
{
  vtkShoeCell* sc = vtkShoeCell::SafeDownCast(c);
  assert(sc);
  sc->CopyFrom( this->CurrentCell );
}

vtkGenericAdaptorCell* vtkShoeCellIterator::GetCell()
{
  return this->CurrentCell;
}

vtkIdType vtkShoeCellIterator::GetCellId()
{
  return this->CurrentCell->CellID;
}

void vtkShoeCellIterator::Next()
{
  switch( this->TraversalStyle )
    {
  case MeshOrder:
    this->NextCellMeshOrder();
    break;
  case Boundary:
    this->NextCellBoundary();
    break;
  case ExtBoundary:
    this->NextCellExtBoundary();
    break;

  default:
    vtkWarningMacro("Can't traverse cells this way yet.");
    }
}

void vtkShoeCellIterator::NextCellMeshOrder()
{
  if ( this->Dimension < 0 )
    {
    this->Iter->Cell++;
    if ( this->Iter->Cell != this->DataSet->CellRecs->Cells.end() )
      {
      this->CurrentCell->CellID = this->Iter->Cell.where();
      this->CurrentCell->CellRecord = this->DataSet->CellRecs->Cells[ this->CurrentCell->CellID ];
      }
    }
  else
    {
    do
      {
      this->Iter->Cell++;
      if ( this->Iter->Cell != this->DataSet->CellRecs->Cells.end() )
        {
        this->CurrentCell->CellID = this->Iter->Cell.where();
        this->CurrentCell->CellRecord = this->DataSet->CellRecs->Cells[ this->CurrentCell->CellID ];
        }
      } while ( ! this->EndCheckMeshOrder() &&
                this->CurrentCell->CellRecord.Species->Meta->Dimension != this->Dimension );
    }
  this->CurrentCell->InvalidateCache();
}

void vtkShoeCellIterator::NextCellBoundary()
{
}

void vtkShoeCellIterator::NextCellExtBoundary()
{
}

void vtkShoeCellIterator::NextCellCellBoundary()
{
}

int vtkShoeCellIterator::IsAtEnd()
{
  switch( this->TraversalStyle )
    {
  case MeshOrder:
    return this->EndCheckMeshOrder();
  case Boundary:
    return this->EndCheckBoundary();
  case ExtBoundary:
    return this->EndCheckExtBoundary();

  default:
    vtkWarningMacro("Can't traverse cells this way yet.");
    }
  return 1;
}

int vtkShoeCellIterator::EndCheckMeshOrder()
{
  return this->Iter->Cell == this->DataSet->CellRecs->Cells.end();
}

int vtkShoeCellIterator::EndCheckBoundary()
{
  return 1;
}

int vtkShoeCellIterator::EndCheckExtBoundary()
{
  return 1;
}

int vtkShoeCellIterator::EndCheckCellBoundary()
{
  return 1;
}

void vtkShoeCellIterator::SetCellById( vtkIdType cellID )
{
  this->CurrentCell->CellID = cellID;
  this->CurrentCell->CellRecord = this->DataSet->CellRecs->Cells[ cellID ];
  this->CurrentCell->Meta = this->CurrentCell->CellRecord.Species->Meta;
  this->CurrentCell->InvalidateCache();
}

void vtkShoeCellIterator::SetBoundaryCellByDOF( vtkIdType DOF )
{
}

void vtkShoeCellIterator::SetParentCell( vtkShoeCell* c )
{
  if ( c == this->CurrentParent )
    {
    return;
    }
  if ( this->CurrentParent )
    {
    this->CurrentParent->UnRegister( this );
    }
  this->CurrentParent = c;
  if ( this->CurrentParent )
    {
    this->CurrentParent->Register( this );
    }
  this->Modified();
}
