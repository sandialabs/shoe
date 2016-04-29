// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __vtkShoeCellRecord_h
#define __vtkShoeCellRecord_h

#include "vtksnlConfigure.h" // for uint32_t
#include "vtkSystemIncludes.h" // for vtkIdType

class vtkShoeCellSpecies;

class vtkShoeCellRecord
{
public:
  // Description:
  // A reference to the definition of the cell:<p><ul>
  // <li> its shape,
  // <li> its interpolant,
  // <li> its product space,
  // <li> the interpolation order for each geometric and function field,
  // <li> its operations (functions performed on a cell that vary by cell type)
  // </ul>
  vtkShoeCellSpecies* Species;

  // Description:
  // An index into the Connectivity array that defines which corner and
  // higher order "nodes" make up the cell.
  vtkIdType Offset;

  // Description:
  // A bit array of permutations for each of the edge and face nodes referenced
  // by this cell in the vtkShoeMesh::Connectivity array.
  //
  // Because, for example, a higher order face "node" may be shared by two volume
  // cells that do not share the same coordinate system, we must store how each
  // of the face "nodes" is oriented with respect to each volume cell.
  //
  // For edges, the values associated with a single higher order edge "node" may be
  // stored in forwards or backwards order.
  //
  // For faces, there are up to 8 possible permutations:<p><ul>
  // <li> the starting location may be any of 4 corner nodes on a brick element's face and
  // <li> the values may be in row-first or column-first order.
  // </ul>
  // 
  // To use NodePermutations, first ask for a particular face or edge's permutation.
  // Then, for each degree of freedom at a "node", call GetFaceIndex or GetEdgeIndex
  // to retrieve an offset into the "node"'s list of values.
  uint32_t NodePermutations;

  int  GetFacePermutation(int) const;
  void SetFacePermutation(int,int);

  bool GetEdgePermutation(int) const;
  void SetEdgePermutation(int,bool);

  // Description:
  // Convert \f$(r,s)\f$ tuples specified in a cell's local coordinate frame to/from
  // \f$(r^{\prime},s^{\prime})\f$ specified in the storage coordinate frame.
  // This modifies values in the \a rs array.
  // @param N the number of tuples in the \a rs array.
  // @param skip is an additional number of doubles you would like to skip over as
  //        each tuple is processed. That is, instead of adding 2 to \a rs to get
  //        the next tuple to conver, add 2 + \a skip.
  static void ParametricFaceCoordsToStorageOrder( vtkIdType N, double* rs, int perm );
  static void ParametricFaceCoordsFromStorageOrder( vtkIdType N, double* rs, int perm, int skip=0 );
  void ParametricFaceCoordsToStorageOrder( int face, double* rs ) const;
  void ParametricFaceCoordsFromStorageOrder( int face, double* rs ) const;

  // Description:
  // Convert \f$r\f$ specified in a cell's local coordinate frame to/from
  // \f$r^{\prime}\f$ specified in the storage coordinate frame.
  // This modifies the value in the \a r array.
  // @param N the number of tuples in the \a r array.
  // @param skip is an additional number of doubles you would like to skip over as
  //        each tuple is processed. That is, instead of adding 2 to \a r to get
  //        the next tuple to convert, add 2 + \a skip.
  static void ParametricEdgeCoordsToStorageOrder( vtkIdType N, double* r, bool perm );
  static void ParametricEdgeCoordsFromStorageOrder( vtkIdType N, double* r, bool perm, int skip=0 );
  void ParametricEdgeCoordsToStorageOrder( int edge, double* r ) const;
  void ParametricEdgeCoordsFromStorageOrder( int edge, double* r ) const;
};

//BTX
inline void vtkShoeCellRecord::ParametricFaceCoordsToStorageOrder( int face, double* rs ) const
{ vtkShoeCellRecord::ParametricFaceCoordsToStorageOrder( 1, rs, this->GetFacePermutation( face ) ); }
inline void vtkShoeCellRecord::ParametricFaceCoordsFromStorageOrder( int face, double* rs ) const
{ vtkShoeCellRecord::ParametricFaceCoordsFromStorageOrder( 1, rs, this->GetFacePermutation( face ) ); }

inline void vtkShoeCellRecord::ParametricEdgeCoordsToStorageOrder( int edge, double* r ) const
{ vtkShoeCellRecord::ParametricEdgeCoordsToStorageOrder( 1, r, this->GetEdgePermutation( edge ) ); }
inline void vtkShoeCellRecord::ParametricEdgeCoordsFromStorageOrder( int edge, double* r ) const
{ vtkShoeCellRecord::ParametricEdgeCoordsFromStorageOrder( 1, r, this->GetEdgePermutation( edge ) ); }
//ETX

#endif // __vtkShoeCellRecord_h
