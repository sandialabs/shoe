/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include "vtksnlConfigure.h"
#include "vtkDataRecords.h"
#include "vtkPolynomialSystem.h"
#include "vtkShoeCellRecord.h"
#include "Elements/ShoeHexahedron.h"

int ShoeHexahedron::VerticesOnFace[6] = { 4, 4, 4, 4, 4, 4 };
int ShoeHexahedron::EdgeArray[12][2] = {{0, 1}, {1, 2}, {3, 2}, {0, 3}, {4, 5}, {5, 6}, {7, 6}, {4, 7}, {0, 4}, {1, 5}, {3, 7}, {2, 6}};
int ShoeHexahedron::FaceArray[6][4] = {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4}, {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};
int ShoeHexahedron::EdgesOnFace[6][4] = {{8, 7, 10, 3}, {1, 11, 5, 9}, {0, 9, 4, 8}, {10, 6, 11, 2}, {3, 2, 1, 0}, {4, 5, 6, 7}};
int ShoeHexahedron::FacesOnEdge[12][2] = {{2, 4}, {1, 4}, {3, 4}, {0, 4}, {2, 5}, {1, 5}, {3, 5}, {0, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3}};
double ShoeHexahedron::ParametricCenter[3] = { 0., 0., 0. };
double ShoeHexahedron::ParametricCorners[24] =
{
  -1., -1., -1.,
   1., -1., -1.,
   1.,  1., -1.,
  -1.,  1., -1.,
  -1., -1.,  1.,
   1., -1.,  1.,
   1.,  1.,  1.,
  -1.,  1.,  1.,
};
double ShoeHexahedron::BoundaryEdgeTransforms[12][6] = // FIXME: These should be derived from EdgeArray and ParametricCorners, not hand-written!
{
    {  0., -1., -1.,    1., 0., 0. }, // edge 0-1
    { -1.,  0., -1.,    0., 1., 0. }, // edge 1-2
    {  0.,  1., -1.,    1., 0., 0. }, // edge 3-2
    {  1.,  0., -1.,    0., 1., 0. }, // edge 0-3

    {  0., -1.,  1.,    1., 0., 0. }, // edge 4-5
    { -1.,  0.,  1.,    0., 1., 0. }, // edge 5-6
    {  0.,  1.,  1.,    1., 0., 0. }, // edge 7-6
    {  1.,  0.,  1.,    0., 1., 0. }, // edge 4-7

    { -1., -1.,  0.,    0., 0., 1. }, // edge 0-4
    {  1., -1.,  0.,    0., 0., 1. }, // edge 1-5
    { -1.,  1.,  0.,    0., 0., 1. }, // edge 3-7
    {  1.,  1.,  0.,    0., 0., 1. }  // edge 2-6
};
double ShoeHexahedron::BoundaryFaceTransforms[6][12] = // FIXME: These should be derived from FaceArray and ParametricCorners, not hand-written!
{
    { -1.,  0.,  0.,    0., 1., 0.,    0., 0., 1.,    1.,  0.,  0. }, // face 0-4-7-3 (-X)
    {  1.,  0.,  0.,    0., 1., 0.,    0., 0., 1.,   -1.,  0.,  0. }, // face 1-2-6-5 (+X)

    {  0., -1.,  0.,    1., 0., 0.,    0., 0., 1.,    0.,  1.,  0. }, // face 0-1-5-4 (-Y)
    {  0.,  1.,  0.,    1., 0., 0.,    0., 0., 1.,    0., -1.,  0. }, // face 3-7-6-2 (+Y)

    {  0.,  0., -1.,    1., 0., 0.,    0., 1., 0.,    0.,  0.,  1. }, // face 0-3-2-1 (-Z)
    {  0.,  0.,  1.,    1., 0., 0.,    0., 1., 0.,    0.,  0., -1. }, // face 4-5-6-7 (+Z)
};

void ShoeHexahedron::LagrangeTensorPermuteEdge( double* CellCache, int NumberOfComponents, vtkDataRecords* dof, vtkIdType dconni, int eord, bool perm )
{
  if ( perm )
    {
    int i,j;
    double* TupleData;
    TupleData = dof->GetRecordAsDouble( dconni );
    int nt = eord - 2;
    int sn = 0;
    
    for ( i = nt; i > -1; -- i )
      {
      for ( j = 0; j < NumberOfComponents; ++ j )
        {
        CellCache[sn++] = TupleData[i * NumberOfComponents + j];
        }
      }
    }
  else
    {
    dof->GetRecord( dconni, CellCache );
    }
}

inline int D4perm_0321( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;
  int n = np1 - 1;

  for ( i = 0; i < mp1; ++i )
    {
    for ( j = 0; j < np1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
        Barray[( i * np1 + j ) * noc + k] = Aarray[( ( n - j ) * mp1 + i ) * noc + k];
        }
      }
    }
  return np1;
}

inline int D4perm_0123( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;
  int m = mp1 - 1;

  for ( i = 0; i < mp1; ++i )
    {
    for ( j = 0; j < np1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
        Barray[( i * np1 + j ) * noc + k] = Aarray[( j * mp1 + m - i ) * noc + k];
        }
      }
    }
  return np1;
}

inline int D4perm_02_13( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;
  int m = mp1 - 1;
  int n = np1 - 1;

  for ( i = 0; i < np1; ++i )
    {
    for ( j = 0; j < mp1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
        Barray[( i * mp1 + j ) * noc + k] = Aarray[( ( n - i ) * mp1 + m - j ) * noc + k];
        }
      }
    }
  return mp1;
}

inline int D4perm_01_23( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;
  int m = mp1 - 1;

  for ( i = 0; i < np1; ++i )
    {
    for ( j = 0; j < mp1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
        Barray[( i * mp1 + j ) * noc + k] = Aarray[( i * mp1 + m - j ) * noc + k];
        }
      }
    }
  return mp1;
}

inline int D4perm_03_12( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;
  int n = np1 - 1;

  for ( i = 0; i < np1; ++i )
    {
    for ( j = 0; j < mp1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
        Barray[( i * mp1 + j ) * noc + k] = Aarray[( ( n - i ) * mp1 + j ) * noc + k];
        }
      }
    }
  return mp1;
}

inline int D4perm_02( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;

  for ( i = 0; i < mp1; ++i )
    {
    for ( j = 0; j < np1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
      Barray[( i * np1 + j ) * noc + k] = Aarray[( j * mp1 + i ) * noc + k];
        }
      }
    }
  return np1;
}

inline int D4perm_13( double *Aarray, double *Barray, int mp1, int np1, int noc )
{
  int i,j,k;
  int m = mp1 - 1;
  int n = np1 - 1;

  for ( i = 0; i < mp1; ++i )
    {
    for ( j = 0; j < np1; ++j )
      {
      for ( k = 0; k < noc; ++k )
        {
        Barray[( i * np1 + j ) * noc + k] = Aarray[( ( n - j ) * mp1 + m - i ) * noc + k];
        }
      }
    }
  return np1;
}

void ShoeHexahedron::LagrangeTensorPermuteFace( double* CellCache, int NumberOfComponents, vtkDataRecords* dof, vtkIdType dconni, int *ford, int fdir[2], uint32_t perm )
{
  if ( perm )
    {
    double* TupleData;
    TupleData = dof->GetRecordAsDouble( dconni );
    int mp1 = ford[fdir[0]] - 1;
    int np1 = ford[fdir[1]] - 1;
    int ncol;

    switch( perm )
      {
      case 1:
        ncol = D4perm_13( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      case 2:
        ncol = D4perm_0123( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      case 3:
        ncol = D4perm_03_12( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      case 4:
        ncol = D4perm_02_13( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      case 5:
        ncol = D4perm_02( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      case 6:
        ncol = D4perm_0321( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      case 7:
        ncol = D4perm_01_23( TupleData, CellCache, mp1, np1, NumberOfComponents );
        break;
      default:
        vtkGenericWarningMacro( "Undefined perm value: " << perm << "." );
        break;
      }

#if 0
    int i,j,k;

    cout << "perm " << perm << " on " << fdir[0] << " -- " << fdir[1] << " :" << endl;
    for ( i = 0; i < np1; ++i )
      {
      for ( j = 0; j < mp1; ++j )
        {
        for ( k = 0; k < NumberOfComponents; ++k )
          {
          cout << TupleData[i * mp1 * NumberOfComponents + j * NumberOfComponents + k] << " ";
          }
        cout << endl;
        }
      }
    cout << "permuted:" << endl;
    int a,b;
    if ( ncol == mp1)
      {
      a =  np1;
      b =  mp1;
      }
    else
      {
      a =  mp1;
      b =  np1;
      }

    for ( i = 0; i < a; ++i )
      {
      for ( j = 0; j < b; ++j )
        {
        for ( k = 0; k < NumberOfComponents; ++k )
          {
          cout << CellCache[i * ncol * NumberOfComponents + j * NumberOfComponents + k] << " ";
          }
        cout << endl;
        }
      }
#endif // 0

    }
  else
    {
    dof->GetRecord( dconni, CellCache );
    }
}

void ShoeHexahedron::EmbedEdgeCoords( int edge, double* tuple )
{
  if ( edge < 8 )
    {
    if ( edge % 2 )
      {
      tuple[1] = tuple[0];
      tuple[0] = ( ( edge / 2 ) % 2 ) ? -1. : 1.;
      }
    else
      {
      tuple[1] = ( edge % 4 ) ? 1. : -1.;
      }
    tuple[2] = edge / 4 ? 1. : -1.;
    }
  else
    {
    tuple[2] = tuple[0];
    tuple[0] = (edge == 8 || edge == 10) ? -1. : 1.;
    tuple[1] = edge < 10 ? -1. : 1.;
    }
}

void ShoeHexahedron::EmbedEdgeCoordsInFace( int edge, int face, double* tuple )
{
  ShoeHexahedron::EmbedEdgeCoords( edge, tuple );
  ShoeHexahedron::ProjectCoordsToFace( face, tuple );
#if 0
  // FIXME!! +/-X needs to be checked. Others are guaranteed not to handle the
  // case where an edge is perpendicular to a face properly.
  switch (face)
    {
    case 0: // -X
    case 1: // +X
      if ( edge < 8 )
	{
	  if ( edge % 2 )
	    {
	      if ( edge == 3 || edge == 7 )
		{ // edges 3,7 go in opposite direction of t coordinate
		  tuple[0] = -tuple[0];
		}
	    }
	  else
	    {
	      tuple[0] = edge % 4 ? +1. : -1.; // edges 0,4 have t=-1, edges 2,6 have t=+1
	    }
	  tuple[1] = edge < 5 ? -1. : +1.;
	}
      else
	{
	  tuple[1] = tuple[0];
	  tuple[0] = edge < 10 ? -1. : +1.;
	}
      break;
    case 2: // -Y
    case 3: // +Y
      if ( edge < 8 )
	{
	  tuple[1] = edge < 4 ? -1. : +1.;
	}
      else
	{
	  tuple[1] = tuple[0];
	  tuple[0] = edge % 2 ? +1. : -1.;
	}
      break;
    case 4: // -Z
    case 5: // +Z
      if ( edge % 2 )
	{
	  tuple[1] = tuple[0];
	  tuple[0] = (edge - 1) % 4 ? -1. : +1.;
	}
      else
	{
	  tuple[1] = edge % 4 ? +1. : -1.;
	}
      break;
    default:
      vtkGenericWarningMacro( "I am confused. And everybody knows that when you're confused, the best thing to do is pretend you're a deer staring at headlights." );
      break;
    }
#endif // 0
}

// This takes coordinates from a face DOF node and puts the
// values in the correct coordinate for the volume node. For
// instance, a DOF node might store a critical point on the
// 3rd face of a hexahedron as as ( 0.2, 0.5, x ), where x
// indicates x an unused cooordinate. You might then
// transform the coordinates from storage order into the local
// coordinate frame with
//    vtkShoeCellRecord::ParametricCoordsFromStorageOrder
// to obtain ( 0.5, -0.2, x ). The routine below performs the
// final step by mapping ( 0.5, -0.2, x ) to ( 0.5, -1.0, -0.2 ).
// It knows to do this because the 3rd face of a hex is always
// at s = -1 and the given "face" coordinates correspond to
// r and t. There is a similar routine for edges, but it is 
// much simpler.
void ShoeHexahedron::EmbedFaceCoords( int face, double* tuple )
{
  switch( face )
    {
    case 0:
      tuple[2] = tuple[1];
      tuple[1] = tuple[0];
      tuple[0] = -1.;
      break;
    case 1:
      tuple[2] = tuple[1];
      tuple[1] = tuple[0];
      tuple[0] =  1.;
      break;
    case 2:
      tuple[2] = tuple[1];
      tuple[1] = -1.;
      break;
    case 3:
      tuple[2] = tuple[1];
      tuple[1] =  1.;
      break;
    case 4:
      tuple[2] = -1.;
      break;
    case 5:
      tuple[2] =  1.;
      break;
    default:
      vtkGenericWarningMacro( "The is no face with index " << face << " in a ShoeHexahedron." );
      break;
    }
}

void ShoeHexahedron::ProjectCoordsToEdge( int edge, double* r )
{
  if ( edge < 8 )
    {
    if ( edge % 2 )
      {
      // it's an s-axis edge
      r[0] = r[1];
      }
    else
      {
      // do nothing: it's an r-axis edge
      }
    }
  else
    {
    // it's a t-axis edge
    r[0] = r[2];
    }
}

void ShoeHexahedron::ProjectCoordsToFace( int face, double* r )
{
  switch (face/2)
    {
    case 0:
      // s-t faces
      r[0] = r[1];
      r[1] = r[2];
      break;
    case 1:
      // r-t faces
      r[1] = r[2];
      break;
    case 2:
      // r-s faces: do nothing
      break;
    }
}

void ShoeHexahedron::ProjectFaceCoordsToEdge( int face, int edge, double* r )
{
  (void)face;
  if ( edge % 2 )
    {
    r[0] = r[1];
    }
  else
    {
    // do nothing, edge is on r axis
    }
}

void ShoeHexahedron::TransformFaceCoordsToStorage( uint32_t perm, double* tuple )
{
  switch ( perm )
    {
    case 0:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
      // FIXME: complete implementation below this line
    case 1:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    case 2:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    case 3:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    case 4:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    case 5:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    case 6:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    case 7:
      // local to storage is ()
      tuple[0] = tuple[2];
      tuple[1] = tuple[3];
      break;
    default:
      vtkGenericWarningMacro("Incorrect permutation: " << perm << ".");
    }
}

int ShoeHexahedron::IsParameterInDomain( double* param, int bdy )
{
  int result = 1;
  int dim;
  if ( bdy < 0 || bdy >= 18 )
    { // use the whole cell
    dim = 3;
    }
  else if ( bdy < 12 )
    { // it's an edge
    dim = 1;
    }
  else
    { // must be a face
    dim = 2;
    }
  double u;
  for ( int i = 0; i < dim; ++i, ++param )
    {
    u = fabs( *param );
    if ( u > 1. )
      {
      result = -1;
      }
    else if ( u > (1. - 1.e-8) ) // FIXME: why is 1.e-8 a good fudge factor?
      {
      result =  0; // on at least one boundary
      }
    }
  return result;
}
