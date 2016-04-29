// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeEnums.h"
#include "vtkShoeAttribute.h"
#include "vtkGenericAdaptorCell.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellGenus.h"
#include "vtkShoeCellIterator.h"
#include "vtkShoeCellSpecies.h"
#include "vtkShoeCellMetaData.h"
#include "vtkShoeOrderTuple.h"
#include "vtkShoePointIterator.h"
#include "vtkShoeBox.h"
#include "vtkShoeBoxPartition.h"
#include "vtkDataRecords.h"
#include "vtkDataRecordsIterator.h"
#include "vtkRegressionTest.h"

#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkGenericAttributeCollection.h>
#include <vtkGenericContourFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGarbageCollector.h>

static double ShoeAttributePts[] =
{
-1., -1., -1. ,
 1., -1., -1. ,
 1.,  1., -1. ,
-1.,  1., -1. ,
-1., -1.,  1. ,
 1., -1.,  1. ,
 1.,  1.,  1. ,
-1.,  1.,  1. ,
 3., -1., -1. ,
 3.,  1., -1. ,
 3., -1.,  1. ,
 3.,  1.,  1. ,
};

static double ShoeAttributeDOF[] =
{
 0., -1., -1. , // 0
 1.,  0., -1. ,
 0.,  1., -1. ,
-1.,  0., -1. ,

 0., -1.,  1. , // 4
 1.,  0.,  1. ,
 0.,  1.,  1. ,
-1.,  0.,  1. ,

-1., -1.,  0. , // 8
 1., -1.,  0. ,
-1.,  1.,  0. ,
 1.,  1.,  0. ,

-1.,  0.,  0. , // 12
 1.,  0.,  0. ,
 0., -1.,  0. ,
 0.,  1.,  0. ,

 0.,  0., -1. ,
 0.,  0.,  1. ,

 0.,  0.,  0. , // 18

 2., -1., -1. , // 19
 3.,  0., -1. ,
 2.,  1., -1. ,

 2., -1.,  1. , // 22
 3.,  0.,  1. ,
 2.,  1.,  1. ,

 3., -1.,  0. , // 25
 3.,  1.,  0. ,

 3.,  0.,  0. , // 27

 2., -1.,  0. ,
 2.,  1.,  0. ,

 2.,  0., -1. , // 27
 2.,  0.,  1. ,

 2.,  0.,  0. , // 32
};

static double ShoeScalarPts[] =
{
  -1.,
   1.,
   2.,
   0.,

   0.,
   1.,
  -.5,
   .5,

  -1.,
   .9,
   .3,
   .5
};

static double ShoeScalarDOF[] =
{ // some nice random numbers
  -1.0511238,
   0.6419952,
   0.0203955,
   0.4717119,
   0.4631557,
   0.2172587,
  -1.0107735,
  -1.0625039,
  -0.5642809,
  -0.7541493,
  -1.8459450,
   0.7072287,
   0.0053725,
  -1.0461605,
   1.4775577,
   1.1940160,
   1.3181717,
   1.9951038,
   0.8072696,
   0.4845049,
  -1.8695624,
  -1.4175837,
  -1.5098062,
  -0.7453152,
   0.0809994,
  -1.9680787,
  -0.5888761,
   0.6309946,
  -0.9367993,
  -0.7582506,
   0.4309616,
   1.6155434,
   0.1517625
};

static vtkIdType ShoeTestCellConn[] =
{
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
  1,  8,  9,  2,  5, 10, 11,  6, 19, 20, 21,  1, 22, 23, 24,  5,  9, 25, 11, 26, 13, 27, 28, 29, 30, 31, 32
};

int ShoeBox( int argc, char** argv )
{
  vtkRegressionTest test( "ShoeBox" );
  //vtkGarbageCollector::SetGlobalDebugFlag(1);

  // CREATE A MESH ============================================================
  vtkShoeBox* mesh = vtkShoeBox::New();
  vtkShoeAttribute* att = vtkShoeAttribute::New();
  att->SetName( "Geometry" );
  att->SetNumberOfComponents(3);
  vtkShoeAttribute* sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents(1);

  int i;
  vtkIdType id;
  vtkDoubleArray* pts = vtkDoubleArray::New();
  vtkDataRecords* dof = vtkDataRecords::New();

  pts->SetNumberOfComponents(1);
  pts->SetNumberOfTuples(sizeof(ShoeScalarPts)/sizeof(ShoeScalarPts[0]));

  for (i=0; i<(int)(sizeof(ShoeScalarPts)/sizeof(ShoeScalarPts[0])); i++)
    {
    pts->SetTuple( i, ShoeScalarPts + i );
    }

  dof->SetNumberOfComponents(1);
  for (i=0; i<(int)(sizeof(ShoeScalarDOF)/sizeof(ShoeScalarDOF[0])); i++)
    {
    dof->InsertNextRecord( 1, ShoeScalarDOF + i );
    }

  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->FastDelete();
  dof->FastDelete();

  pts = vtkDoubleArray::New();
  dof = vtkDataRecords::New();

  pts->SetNumberOfComponents(3);
  pts->SetNumberOfTuples(sizeof(ShoeAttributePts)/sizeof(ShoeAttributePts[0])/3);

  for (i=0; i<(int)(sizeof(ShoeAttributePts)/sizeof(ShoeAttributePts[0])/3); i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof->SetNumberOfComponents(3);
  for (i=0; i<(int)(sizeof(ShoeAttributeDOF)/sizeof(ShoeAttributeDOF[0])/3); i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att->SetPointData( pts );
  att->SetDOFData( dof );
  pts->FastDelete();
  dof->FastDelete();

  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); // Need to call this to resize PointLinks, DOFLinks arrays
  att->FastDelete();
  sca->FastDelete();
  vtkShoeCellGenus g;
  g.Shape = shoe::Hexahedron;
  vtkPolyInterpolant in[2] = { Lagrange, Lagrange };
  vtkPolyProductSpace ps[2] = { Tensor, Tensor };
  vtkShoeOrderTuple o[2];
  o[0].Set( 2, 2, 2 );
  o[1].Set( 2, 2, 2 );
  int oi[2] = { 0, 1 };
  int n = 2;
  vtkShoeCellSpecies* sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, mesh->GetAttributes() );
  if ( ! sp )
    {
    test.StdOut() << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }
  vtkShoeOrderTuple od;
  att->GetOrder( sp->GetId(), od );
  test.StdOut() << "Interpolant " << att->GetInterpolant( sp->GetId() ) << vtkstd::endl
    << "Product space " << att->GetProductSpace( sp->GetId() ) << vtkstd::endl
    << "Order " << od.Order[0] << " " << od.Order[1] << " " << od.Order[2] << vtkstd::endl;
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn + 27, 0 /*permutation*/ );

  // TEST RECORD ITERATOR =====================================================
  test.StdOut() << "Test DataRecordsIterator" << vtkstd::endl;
  vtkDataRecordsIterator* drit = vtkDataRecordsIterator::New();
  drit->SetData( mesh->GetConnectivity() );
  for ( drit->Begin(); ! drit->IsAtEnd(); drit->Next() )
    {
    test.StdOut() << "  Record " << drit->GetRecord() << " (" << mesh->GetConnectivity()->GetNumberOfTuples( drit->GetRecord() ) << ")\n";
    }

  // TEST CELL ITERATOR =======================================================
  test.StdOut() << "Test ShoeCellIterator and ShoePointIterator (over cell)" << vtkstd::endl;
  vtkShoeCellIterator* it = vtkShoeCellIterator::SafeDownCast(mesh->NewCellIterator());
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }
  vtkGenericAdaptorCell* gc;
  vtkShoeCell* scp = 0;
  vtkShoeCell*& sc( scp );
  double pcoords[3] = { 0.3, 0.2, -.9 };
  //double pcoords[3] = { 1., 1., 0. };
  double xcoords[3];
  vtkShoePointIterator* pit = vtkShoePointIterator::SafeDownCast( mesh->NewPointIterator() );
  if ( ! pit )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoePointIterator." << vtkstd::endl;
    }
  it->Begin();
  gc = it->GetCell();
  scp = vtkShoeCell::SafeDownCast( gc );
  if ( sc && sc->GetDimension() == 3 )
    {
    test.StdOut() << "    Test hex edge and face embeddings for this genus of cell" << vtkstd::endl;
    for ( int e = 0; e < sc->GetNumberOfBoundaries(1); ++e )
      {
      xcoords[0] = 0.125; xcoords[1] = xcoords[2] = 0.;
      test.StdOut() << "        Edge " << e << "( " << xcoords[0] << " ) = < ";
      sc->GetMetaData()->EmbedEdgeCoords( e, xcoords );
      test.StdOut() << xcoords[0] << ", " << xcoords[1] << ", " << xcoords[2] << " > == < ";
      xcoords[0] = 0.125; xcoords[1] = xcoords[2] = 0.;
      int f = sc->GetMetaData()->FacesOnEdge[e][0];
      sc->GetMetaData()->EmbedEdgeCoordsInFace( e, f, xcoords );
      test.StdOut() << xcoords[0] << ", " << xcoords[1] << " > --> < ";
      sc->GetMetaData()->EmbedFaceCoords( f, xcoords );
      test.StdOut() << xcoords[0] << ", " << xcoords[1] << ", " << xcoords[2] << " >" << vtkstd::endl;
      }
    }
  for ( ; ! it->IsAtEnd(); it->Next() )
    {
    test.StdOut() << "  Cell " << it->GetCellId();
    gc = it->GetCell();
    scp = vtkShoeCell::SafeDownCast( gc );
    scp->ComputeDOFCriticalPoints( sca );
    gc->InterpolateTuple( att, pcoords, xcoords );
    test.StdOut() << " ( " << pcoords[0] << ", " << pcoords[1] << ", " << pcoords[2] << " ) = ( " << xcoords[0] << ", " << xcoords[1] << ", " << xcoords[2] << " ) x ( ";
    gc->InterpolateTuple( sca, pcoords, xcoords );
    test.StdOut() << xcoords[0] << " )\n";
    test.StdOut() << "    Cell points:";
    pit->SetCell( sc );
    for ( pit->Begin(); !pit->IsAtEnd(); pit->Next() )
      {
      test.StdOut() << " " << pit->GetId();
      }
    test.StdOut() << vtkstd::endl;
    test.StdOut() << "    Boundary faces:";
    for ( int f=0; f<gc->GetNumberOfBoundaries(2); ++f )
      {
      test.StdOut() << " " << gc->IsFaceOnBoundary(f);
      }
    test.StdOut() << vtkstd::endl;
    }
  // TEST MESH PARTITIONING ===================================================
  vtkGenericAttributeCollection* kappa = vtkGenericAttributeCollection::New();
  kappa->InsertNextAttribute( sca );
  double bogusCP[4] = { 0.1, 0.28, 0., 0. };
  sca->GetCriticalPoints()[0]->InsertNextTuple( 12, bogusCP );
  sca->SetCriticalPointsDirty( 0 );
  vtkShoeBoxPartition* ptn = mesh->PartitionMesh( kappa );
  ptn->Dump();
  // TEST POINT ITERATOR ======================================================
  test.StdOut() << "Test ShoePointIterator (over mesh)" << vtkstd::endl;
  pit->SetDataSet(0);
  pit->SetDataSet(mesh);
  for ( pit->Begin(); ! pit->IsAtEnd(); pit->Next() )
    {
    double* x = pit->GetPosition();
    test.StdOut() << "  Id " << pit->GetId() << " ( " << x[0] << ", " << x[1] << ", " << x[2] << " )\n";
    }
  // TEST CONTOURING ==========================================================
  vtkGenericContourFilter *contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetValue( 0, 0.1 );
  contour->Update();
  test.StdOut() << "Contour has " << contour->GetOutput()->GetNumberOfCells() << " triangles." << vtkstd::endl;

  vtkXMLPolyDataWriter* pdw = vtkXMLPolyDataWriter::New();
  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "ShoeBox.vtp" );
  pdw->Write();

  pdw->Delete();
  contour->SetInput( 0 );
  contour->Delete();
  pit->Delete();
  it->Delete();
  drit->Delete();
  mesh->Delete();

  vtkShoeCellMetaData::Shutdown(); // free metadata records
  return test.GetReturnCode( vtkRegressionTest::Passed );
}
