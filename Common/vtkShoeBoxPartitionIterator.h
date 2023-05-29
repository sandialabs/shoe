// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtksnlCommonWin32Header.h"
#include "vtkShoeBoxPartition.h"

#include <vtkstd/vector>

#include "vtkObject.h"

class VTK_SNL_COMMON_EXPORT vtkShoeBoxPartitionIterator : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkShoeBoxPartitionIterator,vtkObject);
  static vtkShoeBoxPartitionIterator* New();
  virtual void PrintSelf( ostream& os, vtkIndent indent );

  void Begin( vtkShoeBoxPartition* partition, vtkIdType dof );
  int IsAtEnd();
  void Next();

  void GetTetrahedron( double*, double*, double*, double* );
  void GetWorldTetrahedron( double*, double*, double*, double* );
  int GetTetrahedronId() const { return this->Tetrahedron; }

  //BTX
  vtkShoeBoxPartition::EleRef GetCurrentEdgeFacet() const;
  //ETX

protected:
  //BTX
  vtkstd::vector<int> Visited;
  //ETX
  vtkIdType Triangle;
  vtkIdType TriangleBegin;
  vtkIdType TriangleEnd;
  int Tetrahedron;
  int TriangleSide;
#define VTK_SBP_DUMP_FORMAT_U 
#ifdef VTK_SBP_DUMP_FORMAT_U
  double Coords[24];
#else // VTK_SBP_DUMP_FORMAT_U
  double Coords[12];
#endif // VTK_SBP_DUMP_FORMAT_U
  vtkShoeBoxPartition* Partition;

  vtkShoeBoxPartitionIterator();
  ~vtkShoeBoxPartitionIterator();

  void MarkInterior( vtkIdType* tris, int* sides );

private:
  void operator = ( const vtkShoeBoxPartitionIterator& ); // Not implemented.
  vtkShoeBoxPartitionIterator( const vtkShoeBoxPartitionIterator& ); // Not implemented.
};

