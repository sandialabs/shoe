// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtksnlConfigure.h"
#include "vtkShoeCellRecord.h"

#include <assert.h>

int  vtkShoeCellRecord::GetFacePermutation( int face ) const
{
  assert( face < 6 && face >= 0 );
  return (this->NodePermutations >> (3*face)) & 0x7;
}

void vtkShoeCellRecord::SetFacePermutation( int face, int permutation )
{
  assert( face < 6 && face >= 0 );
  this->NodePermutations &= UINT32_MAX - (0x7 << (face*3));
  this->NodePermutations |= permutation << (face*3);
}

bool vtkShoeCellRecord::GetEdgePermutation( int edge ) const
{
  return (this->NodePermutations << edge) & (1 << 31);
}

void vtkShoeCellRecord::SetEdgePermutation( int edge, bool orientation )
{
  if ( orientation )
    this->NodePermutations |= 1 << (31 - edge);
  else
    this->NodePermutations &= UINT32_MAX - ( 1 << (31 - edge) );
}

void vtkShoeCellRecord::ParametricFaceCoordsToStorageOrder( vtkIdType N, double* rs, int perm )
{ 
  double tmp;
  switch ( perm )
    {
  case 1:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = -rs[1];
      rs[1] = tmp;
      rs += 2;
      --N;
      }
    break;
  case 2:
    while ( N )
      {
      rs[0] = -rs[0];
      rs[1] = -rs[1];
      rs += 2;
      --N;
      }
    break;
  case 3:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = rs[1];
      rs[1] = -tmp;
      rs += 2;
      --N;
      }
    break;
  case 4:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = rs[1];
      rs[1] = tmp;
      rs += 2;
      --N;
      }
    break;
  case 5:
    while ( N )
      {
      rs[1] = -rs[1];
      rs += 2;
      --N;
      }
    break;
  case 6:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = -rs[1];
      rs[1] = -tmp;
      rs += 2;
      --N;
      }
    break;
  case 7:
    while ( N )
      {
      rs[0] = -rs[0];
      rs += 2;
      --N;
      }
    break;
  default:
  case 0:
    // do nothing;
    break;
    }
}

void vtkShoeCellRecord::ParametricFaceCoordsFromStorageOrder( vtkIdType N, double* rs, int perm, int skip )
{
  double tmp;
  skip += 2;
  switch ( perm )
    {
  case 1:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = rs[1];
      rs[1] = -tmp;
      rs += skip;
      --N;
      }
    break;
  case 2:
    while ( N )
      {
      rs[0] = -rs[0];
      rs[1] = -rs[1];
      rs += skip;
      --N;
      }
    break;
  case 3:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = -rs[1];
      rs[1] = tmp;
      rs += skip;
      --N;
      }
    break;
  case 4:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = rs[1];
      rs[1] = tmp;
      rs += skip;
      --N;
      }
    break;
  case 5:
    while ( N )
      {
      rs[1] = -rs[1];
      rs += skip;
      --N;
      }
    break;
  case 6:
    while ( N )
      {
      tmp = rs[0];
      rs[0] = -rs[1];
      rs[1] = -tmp;
      rs += skip;
      --N;
      }
    break;
  case 7:
    while ( N )
      {
      rs[0] = -rs[0];
      rs += skip;
      --N;
      }
    break;
  default:
  case 0:
    // do nothing;
    break;
    }
}

void vtkShoeCellRecord::ParametricEdgeCoordsToStorageOrder( vtkIdType N, double* r, bool perm )
{
  if ( perm )
    {
    while ( N )
      {
      *r = -(*r);
      --N;
      }
    }
}

void vtkShoeCellRecord::ParametricEdgeCoordsFromStorageOrder( vtkIdType N, double* r, bool perm, int skip )
{
  skip += 1;
  if ( perm )
    {
    while ( N )
      {
      *r = -(*r);
      r += skip;
      --N;
      }
    }
}

