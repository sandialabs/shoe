// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#ifndef __vtkShoeOrderTuple_h
#define __vtkShoeOrderTuple_h

struct vtkShoeOrderTuple
{
  int Order[3];

  inline void Set( int* o ) { Order[0] = o[0]; Order[1] = o[1]; Order[2] = o[2]; }
  inline void Set( int ro, int so, int to ) { Order[0] = ro; Order[1] = so; Order[2] = to; }
  inline void Set( int ro, int so ) { Order[0] = ro; Order[1] = so; }
  inline void Set( int ro ) { Order[0] = ro; }

  inline bool operator == ( const vtkShoeOrderTuple& o ) const
    {
    return (o.Order[0] == this->Order[0]) && (o.Order[1] == this->Order[1]) && (o.Order[2] == this->Order[2]);
    }
  inline bool operator != ( const vtkShoeOrderTuple& o ) const
    {
    return ! ( o == *this );
    }
  inline bool operator < ( const vtkShoeOrderTuple& o ) const
    {
    if ( o.Order[0] > this->Order[0] )
      {
      return true;
      }
    else if ( o.Order[0] == this->Order[0] )
      {
      if ( (o.Order[1] > this->Order[1]) ||
        ( (o.Order[1] == this->Order[1]) && (o.Order[2] > this->Order[2]) ) )
        {
        return true;
        }
      }
    return false;
    }
};

#endif // __vtkShoeOrderTuple_h
