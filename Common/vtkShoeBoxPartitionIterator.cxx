// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeBoxPartitionIterator.h"
#include "vtkShoeBoxPartition.h"
#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkShoeBoxPartitionIterator);
vtkCxxRevisionMacro(vtkShoeBoxPartitionIterator,"$Revision: 8684 $");

void vtkShoeBoxPartitionIterator::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Visited: " << this->Visited.size() << " bytes" << vtkstd::endl;
  os << indent << "Triangle: " << this->Triangle << vtkstd::endl;
  os << indent << "TriangleBegin: " << this->TriangleBegin << vtkstd::endl;
  os << indent << "TriangleEnd: " << this->TriangleEnd << vtkstd::endl;
  os << indent << "Partition: " << this->Partition << vtkstd::endl;
  os << indent << "Tetrahedron: " << this->Tetrahedron << vtkstd::endl;
}

vtkShoeBoxPartitionIterator::vtkShoeBoxPartitionIterator()
{
  this->Triangle = 0;
  this->TriangleBegin = 0;
  this->TriangleEnd = 0;
  this->TriangleSide = 0;
  this->Partition = 0;
  this->Tetrahedron = -1;
}

vtkShoeBoxPartitionIterator::~vtkShoeBoxPartitionIterator()
{
}

void vtkShoeBoxPartitionIterator::Begin( vtkShoeBoxPartition* partition, vtkIdType dof )
{
  this->Partition = partition;
  this->Visited.clear();
  vtkIdType N = this->Partition->GetNumberOfElements( dof );
  if ( N <= 0 )
    {
    this->Partition = 0;
    this->Triangle = 0;
    this->TriangleBegin = 0;
    this->TriangleEnd = -1;
    this->TriangleSide = 0;
    return;
    }
  N = ( 2 * N + 7 ) / 8;
  this->Visited.resize( N, 0 );
  this->Tetrahedron = -1; // a label for each tetrahedron in a DOF node
  this->Triangle = 0;
  this->TriangleBegin = this->Partition->GetFirstTriangle( dof );
  this->TriangleEnd = this->TriangleBegin + this->Partition->GetNumberOfElements( dof );
  this->Next(); // possibly advance and definitely mark this->Visited as needed
}

int vtkShoeBoxPartitionIterator::IsAtEnd()
{
  return this->Triangle + this->TriangleBegin >= this->TriangleEnd;
}

void vtkShoeBoxPartitionIterator::MarkInterior( vtkIdType* tris, int* sides )
{
  for ( int i = 0; i < 3; ++ i )
    {
    vtkIdType byte = ( 2 * (tris[i] - this->TriangleBegin) ) / 8;
    vtkIdType bit = ( sides[i] ? 2 : 1 ) << ( ( 2 * (tris[i] - this->TriangleBegin) ) % 8 );
    this->Visited[byte] |= bit;
    }
}

void vtkShoeBoxPartitionIterator::Next()
{
  vtkIdType tnum = this->Triangle + this->TriangleBegin;
  if ( tnum >= this->TriangleEnd )
    {
    return;
    }

  vtkIdType interiorTris[3];
  int       interiorSides[3];
  do {
    vtkIdType byte = ( 2 * this->Triangle ) / 8;
    int bt0 = 1 << ( ( 2 * this->Triangle ) % 8 );
    int bt1 = 2 << ( ( 2 * this->Triangle ) % 8 );
    if ( !(this->Visited[byte] & bt0) )
      {
      this->Visited[byte] |= bt0;
      this->TriangleSide = 0;
      if ( ! this->Partition->IsTriangleOnHull( tnum, this->TriangleSide ) )
        {
        //vtkstd::cerr << (tnum*8 + this->TriangleSide) << " (" << tnum << ", " << this->TriangleSide << ") is NOT on hull.\n";
        this->Partition->GetTetrahedron( tnum, this->TriangleSide, this->Coords, interiorTris, interiorSides );
        this->MarkInterior( interiorTris, interiorSides );
        ++this->Tetrahedron;
        return;
        }
#if 0
      else
        {
        vtkstd::cerr << (tnum*8 + this->TriangleSide) << " (" << tnum << ", " << this->TriangleSide << ") is on hull.\n";
        }
#endif // 0
      }
    if ( !(this->Visited[byte] & bt1) )
      {
      this->Visited[byte] |= bt1;
      this->TriangleSide = 1;
      if ( ! this->Partition->IsTriangleOnHull( tnum, this->TriangleSide ) )
        {
        //vtkstd::cerr << (tnum*8 + this->TriangleSide) << " (" << tnum << ", " << this->TriangleSide << ") is NOT on hull.\n";
        this->Partition->GetTetrahedron( tnum, this->TriangleSide, this->Coords, interiorTris, interiorSides );
        this->MarkInterior( interiorTris, interiorSides );
        ++this->Tetrahedron;
        return;
        }
#if 0
      else
        {
        vtkstd::cerr << (tnum*8 + this->TriangleSide) << " (" << tnum << ", " << this->TriangleSide << ") is on hull.\n";
        }
#endif // 0
      }
    ++tnum;
    ++this->Triangle;
  } while ( tnum < this->TriangleEnd );
}

void vtkShoeBoxPartitionIterator::GetTetrahedron( double* x, double* y, double* z, double* w )
{
  double* xs = this->Coords;
  double* ys = this->Coords + 3;
  double* zs = this->Coords + 6;
  double* ws = this->Coords + 9;
  for ( int i = 0; i < 3; ++i )
    {
    x[i] = xs[i];
    y[i] = ys[i];
    z[i] = zs[i];
    w[i] = ws[i];
    }
}

void vtkShoeBoxPartitionIterator::GetWorldTetrahedron( double* x, double* y, double* z, double* w )
{
#ifdef VTK_SBP_DUMP_FORMAT_U
  double* xs = this->Coords + 12;
  double* ys = this->Coords + 15;
  double* zs = this->Coords + 18;
  double* ws = this->Coords + 21;
  for ( int i = 0; i < 3; ++i )
    {
    x[i] = xs[i];
    y[i] = ys[i];
    z[i] = zs[i];
    w[i] = ws[i];
    }
#endif // VTK_SBP_DUMP_FORMAT_U
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartitionIterator::GetCurrentEdgeFacet() const
{
  return ( ( this->Triangle + this->TriangleBegin ) << 3 ) | ( this->TriangleSide ? 1 : 0 );
}
