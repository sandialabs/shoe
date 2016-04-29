// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoePointIterator.h"
#include "vtkShoeAttribute.h"
#include "vtkShoeBox.h"
#include "vtkShoeCell.h"

#include <vtkObjectFactory.h>
#include <vtkDataArray.h>

vtkStandardNewMacro(vtkShoePointIterator);
vtkCxxRevisionMacro(vtkShoePointIterator,"$Revision: 5295 $");

void vtkShoePointIterator::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Id: " << this->Id << vtkstd::endl;
  os << indent << "DataSet: " << this->DataSet << vtkstd::endl;
  os << indent << "Points: " << this->Points << vtkstd::endl;
  os << indent << "List: " << this->List << vtkstd::endl;
  os << indent << "ListSize: " << this->ListSize << vtkstd::endl;
}

vtkShoePointIterator::vtkShoePointIterator()
{
  this->Id = -1;
  this->DataSet = 0;
  this->Points = 0;
  this->List = 0;
  this->ListSize = -1;
}

vtkShoePointIterator::~vtkShoePointIterator()
{
  if ( this->List )
    {
    delete [] this->List;
    }
  if ( this->Points )
    {
    this->Points->UnRegister( this );
    }
  if ( this->DataSet )
    {
    this->DataSet->UnRegister( this );
    }
}

void vtkShoePointIterator::Begin()
{
  if ( ! this->Points )
    {
    vtkWarningMacro( "Cannot iterate over a NULL point set" );
    this->Id = -1;
    }

  this->Id = 0;
}

int vtkShoePointIterator::IsAtEnd()
{
  return this->Id < 0;
}

void vtkShoePointIterator::Next()
{
  this->Id++;
  if ( ((! this->ListSize) && (this->Id >= this->Points->GetNumberOfTuples()))
    || ( this->ListSize && ( this->Id >= this->ListSize )) )
    {
    this->Id = -1;
    }
}

double* vtkShoePointIterator::GetPosition()
{
  if ( ! this->ListSize )
    {
    return this->Points->GetTuple( this->Id );
    }
  return this->Points->GetTuple( this->List[ this->Id ] );
}

void vtkShoePointIterator::GetPosition( double* x )
{
  if ( ! this->ListSize )
    {
    this->Points->GetTuple( this->Id, x );
    }
  this->Points->GetTuple( this->List[ this->Id ], x );
}

vtkIdType vtkShoePointIterator::GetId()
{
  return this->ListSize ? this->List[ this->Id ] : this->Id;
}

void vtkShoePointIterator::SetDataSet( vtkShoeBox* d )
{
  if ( d == this->DataSet )
    {
    return;
    }
  if ( this->DataSet )
    {
    this->DataSet->UnRegister( this );
    this->Points->UnRegister( this );
    }
  if ( this->ListSize )
    {
    this->ListSize = 0;
    delete [] this->List;
    this->List = 0;
    }
  this->Id = -1;
  this->DataSet = d;
  this->Modified();
  if ( this->DataSet )
    {
    this->Points = this->DataSet->GetGeometry()->GetPointData();
    this->Points->Register( this );
    this->DataSet->Register( this );
    }
  else
    {
    this->Points = 0;
    }
}

void vtkShoePointIterator::SetCell( vtkShoeCell* c )
{
  if ( ! c )
    {
    this->SetDataSet( 0 );
    }
  else
    {
    this->SetDataSet( c->GetDataSet() );
    }

  if ( this->List )
    {
    delete [] this->List;
    this->List = 0;
    }
  this->ListSize = c ? c->GetNumberOfPoints() : 0 ;
  if ( this->ListSize )
    {
    this->List = new vtkIdType [ this->ListSize ];
    c->GetPointIds( this->List );
    }
  this->Id = -1;
  this->Modified();
}

void vtkShoePointIterator::SetSelection( vtkDataArray* Points, vtkIdType N, vtkIdType* Subset )
{
  if ( this->Points != Points )
    {
    this->SetDataSet(0);
    }
  this->Points = Points;
  this->Points->Register( this );
  if ( this->ListSize != N )
    {
    if ( this->List )
      {
      delete [] this->List;
      }
    this->List = new vtkIdType[ N ];
    this->ListSize = N;
    }
  for ( vtkIdType i=0; i<N; ++i )
    {
    this->List[i] = Subset[i];
    }
  this->Id = -1;
  this->Modified();
}

