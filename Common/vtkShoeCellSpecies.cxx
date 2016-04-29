// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeCellSpecies.h"
#include "vtkShoeCellMetaData.h"
#include "vtkShoeAttribute.h"
#include "vtkGenericAttributeCollection.h"
#include "vtkShoeOrderTuple.h"

#include <vtkstd/vector>
#include <vtkstd/map>

class vtkShoeCellSpeciesP
{
public:
  vtkShoeCellSpeciesP( shoe::CellShape, vtkShoeCellSpecies* );
  ~vtkShoeCellSpeciesP();

  vtkIdType NumberOfCells;
  vtkShoeCellGenus Genus;
  int ID;

  static vtkstd::map<int,vtkShoeCellSpecies*> AllSpecies;
  static int NextID;
};

const int vtkShoeCellGenus::NumberOfConnectivityEntriesByShape[] =
{
  1, // Point
  1, // Curve
  1, // Tube
  1, // Rod
  7, // Triangle
  9, // Quadrilateral
  7, // TriangleShell
  9, // QuadrilateralShell
  11, // Tetrahedron
  27, // Hexahedron
  21, // Wedge
  19 // Pyramid
};

vtkstd::map<int,vtkShoeCellSpecies*> vtkShoeCellSpeciesP::AllSpecies;
int vtkShoeCellSpeciesP::NextID = 0;

vtkShoeCellSpeciesP::vtkShoeCellSpeciesP(
  shoe::CellShape s,
  vtkShoeCellSpecies* parent )
{
  this->Genus.Shape = s;
  this->NumberOfCells = 0;
  this->ID = vtkShoeCellSpeciesP::NextID++;
  vtkShoeCellSpeciesP::AllSpecies.insert( vtkstd::pair<int,vtkShoeCellSpecies*>( this->ID, parent ) );
}

vtkShoeCellSpeciesP::~vtkShoeCellSpeciesP()
{
  vtkShoeCellSpeciesP::AllSpecies.erase( this->ID );
}

vtkShoeCellSpecies::vtkShoeCellSpecies( vtkShoeCellGenus& g )
{
  this->Meta = vtkShoeCellMetaData::ByShape[ g.Shape ];
  this->Data = new vtkShoeCellSpeciesP( g.Shape, this );
}

vtkShoeCellSpecies::~vtkShoeCellSpecies()
{
  delete this->Data;
}

vtkShoeCellSpecies* vtkShoeCellSpecies::GetSpeciesById( int id )
{
  vtkstd::map<int,vtkShoeCellSpecies*>::iterator it = vtkShoeCellSpeciesP::AllSpecies.find( id );
  if ( it == vtkShoeCellSpeciesP::AllSpecies.end() )
    {
    return 0;
    }
  return it->second;
}

int vtkShoeCellSpecies::GetId()
{
  return this->Data->ID;
}

int vtkShoeCellSpecies::GetNumberOfConnectivityEntries()
{
  return this->Data->Genus.GetNumberOfConnectivityEntries();
}

vtkShoeCellSpecies* vtkShoeCellSpecies::FindOrCreate(
  vtkShoeCellGenus& g,
  vtkPolyInterpolant* desiredInterpolants, vtkPolyProductSpace* desiredProdSpaces,
  vtkShoeOrderTuple* desiredOrders, int* desiredOrderIndices, int n,
  vtkGenericAttributeCollection* a )
{
  vtkstd::map<int,vtkShoeCellSpecies*>::iterator it = vtkShoeCellSpeciesP::AllSpecies.begin();
  for ( ; it != vtkShoeCellSpeciesP::AllSpecies.end(); ++it )
    {
    if ( it->second->Data->Genus != g )
      continue;
    int match = 1;
    int itID = it->second->Data->ID;
    for ( int i=0; i<n; ++i )
      {
      vtkShoeAttribute* at = vtkShoeAttribute::SafeDownCast(a->GetAttribute(desiredOrderIndices[i]));
      if ( (! at) ||
           ( at->GetInterpolant( itID ) != desiredInterpolants[i] ) ||
           ( at->GetProductSpace( itID ) != desiredProdSpaces[i] ) )
        {
        match = 0;
        break;
        }
      vtkShoeOrderTuple o;
      at->GetOrder( itID, o );
      if ( desiredOrders[i] != o )
        {
        match = 0;
        break;
        }
      }
    if ( match )
      {
      return it->second;
      }
    }
  // No match, create a new entry and add the order info to all the attributes
  vtkShoeCellSpecies* entry = new vtkShoeCellSpecies(g);
  entry->SetInterpolantInfo( a, desiredOrderIndices, n, desiredInterpolants, desiredProdSpaces, desiredOrders );

  return entry;
}

void vtkShoeCellSpecies::SetInterpolantInfo( vtkGenericAttributeCollection* a, int* indices, int n,
    vtkPolyInterpolant* interp, vtkPolyProductSpace* pspace, vtkShoeOrderTuple* order )
{
  int id = this->GetId();
  for ( int i=0; i<n; ++i )
    {
    vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast(a->GetAttribute( i ));
    if ( sa )
      {
      sa->SetInterpolant( id, interp[i] );
      sa->SetProductSpace( id, pspace[i] );
      sa->SetOrder( id, order[i] );
      }
    }
}

void vtkShoeCellSpecies::SetAttributeOrder( vtkGenericAttributeCollection* a, vtkShoeOrderTuple* desiredOrders, int* desiredOrderIndices, int n )
{
  for ( int i=0; i<n; ++i )
    {
    vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( a->GetAttribute( desiredOrderIndices[i] ) );
    if ( sa )
      {
      this->SetAttributeOrder( sa, desiredOrders[i] );
      }
    else
      {
      vtkGenericWarningMacro( "vtkShoeCellSpecies::SetAttributeOrder called with an attribute collection that contains a non-SHOE attribute, " << a );
      }
    }
}

void vtkShoeCellSpecies::SetAttributeOrder( vtkShoeAttribute* a, vtkShoeOrderTuple& o )
{
  a->SetOrder( this->Data->ID, o );
}

vtkIdType vtkShoeCellSpecies::GetNumberOfCells()
{
  return this->Data->NumberOfCells;
}

void vtkShoeCellSpecies::Reference( vtkIdType cnt )
{
  this->Data->NumberOfCells += cnt;
}

void vtkShoeCellSpecies::Dereference( vtkIdType cnt )
{
  this->Data->NumberOfCells -= cnt;
}

