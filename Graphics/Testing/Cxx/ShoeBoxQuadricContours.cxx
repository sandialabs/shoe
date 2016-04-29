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
#include "vtkShoeOrderTuple.h"
#include "vtkShoePointIterator.h"
#include "vtkShoeBox.h"
#include "vtkDataRecords.h"
#include "vtkDataRecordsIterator.h"
#include "vtkRegressionTest.h"
#include "vtkPolynomialSystem.h"
#include "vtkShoeBoxContourFilter.h"

#include "vtkGenericDataSetTessellator.h"
#include "vtkSimpleCellTessellator.h"
#include "vtkGenericStreamingTessellator.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkUnstructuredGridWriter.h"

#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkGenericAttributeCollection.h>
#include <vtkGenericContourFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkUnstructuredGridWriter.h>

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
-1.,  1.,  0. , // Kitwareism
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

 3.,  0.,  0. ,
 2., -1.,  0. ,
 2.,  1.,  0. ,
 2.,  0., -1. , // 27
 2.,  0.,  1. ,

 2.,  0.,  0. , // 32
};

static double ShoeScalarPts[] =
{
   1.,
   1.,
   1.,
   1.,

   1.,
   1.,
   1.,
   1.,

   5.,
   5.,
   5.,
   5.
};

static double ShoeScalarDOF[] =
{ 
  2.,
  0.,
  2.,
  0.,
  2.,
  0.,
  2.,
  0.,
  0.,
  0.,
  0.,
  0.,
 -1.,
 -1.,
  1.,
  1.,
  1.,
  1.,
  0.,
  2., 
  4.,
  2.,
  2.,
  4.,
  2.,
  4.,
  4.,
  3.,
  1.,
  1.,
  1.,
  1.,
  0.
};

static double EScalarDOF[] =
{ 
  3.,
  3.,
  3.,
  3.,
  3.,
  3.,
  3.,
  3.,
  2.,
  2.,
  2.,
  2.,
  1.,
  1.,
  1.,
  1.,
  2.,
  2.,
  0.
};

static double H1ScalarDOF[] =
{ 
  0.,
  0.,
  0.,
  0.,
  0.,
  0.,
  0.,
  0.,
  2.,
  2.,
  2.,
  2.,
  1.,
  1.,
  1.,
  1.,
 -1.,
 -1.,
  0.
};

static double H2ScalarDOF[] =
{ 
  0.,
  0.,
  0.,
  0.,
  0.,
  0.,
  0.,
  0.,
 -2.,
 -2.,
 -2.,
 -2.,
 -1.,
 -1.,
 -1.,
 -1.,
  1.,
  1.,
  0.
};

static double EPScalarDOF[] =
{ 
  2.,
  2.,
  2.,
  2.,
  0.,
  0.,
  0.,
  0.,
  2.,
  2.,
  2.,
  2.,
  1.,
  1.,
  1.,
  1.,
  1.,
 -1.,
  0.
};

static double HPScalarDOF[] =
{ 
 -2.,
  0.,
 -2.,
  0.,
  0.,
  2.,
  0.,
  2.,
  0.,
  0.,
  0.,
  0.,
  1.,
  1.,
 -1.,
 -1.,
 -1.,
  1.,
  0.
};

static vtkIdType ShoeTestCellConn[] =
{
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
  1,  8,  9,  2,  5, 10, 11,  6, 19, 20, 21,  1, 22, 23, 24,  5,  9, 25, 11, 26, 13, 27, 28, 29, 30, 31, 32
};

int ShoeBoxQuadricContours( int vtkNotUsed(argc), char** vtkNotUsed(argv) )
{
  // DEFINE VARIABLES NEEDED FOR TESTING ================================
  vtkShoeCell* sc;
  vtkShoeCellIterator* it;
  vtkShoeBox* mesh = vtkShoeBox::New();
  vtkShoeAttribute* att = vtkShoeAttribute::New();
  vtkShoeAttribute* sca;
  vtkShoeCellGenus g;
  vtkGenericContourFilter *contour;

  int i;
  double range[2];
  vtkIdType id;
  vtkDoubleArray* pts;
  vtkDataRecords* dof;

  // CREATE A MESH ============================================================
  att->SetName( "Geometry" );
  att->SetNumberOfComponents( 3 );
  sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents( 1 );

  g.Shape = shoe::Hexahedron;
  vtkPolyInterpolant in[2] = { Lagrange, Lagrange };
  vtkPolyProductSpace ps[2] = { Tensor, Tensor };
  vtkShoeOrderTuple o[2];
  o[0].Set( 2, 2, 2 );
  o[1].Set( 2, 2, 2 );
  int oi[2] = { 0, 1 };
  int n = 2;

  // TEST 2-CELL CONTOURING ===================================================
  vtkRegressionTest test( "ShoeBoxQuadricContours" );
  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples(sizeof(ShoeScalarPts)/sizeof(ShoeScalarPts[0]));

  for (i=0; i<int(sizeof(ShoeScalarPts)/sizeof(ShoeScalarPts[0])); i++)
    {
    pts->SetTuple( i, ShoeScalarPts + i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents(1);
  for (i=0; i<int(sizeof(ShoeScalarDOF)/sizeof(ShoeScalarDOF[0])); i++)
    {
    dof->InsertNextRecord( 1, ShoeScalarDOF + i );
    }

  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  sca->SetRangeStyle( VTK_RANGE_TIGHT );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents(3);
  pts->SetNumberOfTuples(sizeof(ShoeAttributePts)/sizeof(ShoeAttributePts[0])/3);
  for (i=0; i<int(sizeof(ShoeAttributePts)/sizeof(ShoeAttributePts[0])/3); i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents(3);
  for (i=0; i<int(sizeof(ShoeAttributeDOF)/sizeof(ShoeAttributeDOF[0])/3); i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att->SetPointData( pts );
  att->SetDOFData( dof );
  att->SetRangeStyle( VTK_RANGE_TIGHT );
  pts->Delete();
  dof->Delete();

  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); // Need to call this to resize PointLinks, DOFLinks arrays
  att->Delete();
  sca->Delete();

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
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn + 27, 0 /*permutation*/ );

  contour = vtkGenericContourFilter::New();
  contour->SetNumberOfContours( 3 );
  contour->SetInput( mesh );
  contour->SetValue( 0, 0. );
  contour->SetValue( 1, 0.75 );
  contour->SetValue( 2, -0.5 );
  contour->Update();

  vtkXMLPolyDataWriter* pdw = vtkXMLPolyDataWriter::New();
  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "ShoeBoxQuadricContours.vtp" );
  pdw->Write();
  contour->Delete();

  // TEST BRUTE FORCE SAMPLING  ===================================
  vtkSimpleCellTessellator* simpleTess = vtkSimpleCellTessellator::SafeDownCast( mesh->GetTessellator() );
  if ( simpleTess )
    {
      simpleTess->SetMaxSubdivisionLevel( 2 );
      simpleTess->SetFixedSubdivisions( 2 );
      vtkstd::cout << "# Using the simple Tessellator\n";
    }
  vtkGenericStreamingTessellator* streamTess = vtkGenericStreamingTessellator::SafeDownCast( mesh->GetTessellator() );
  if ( streamTess )
    {
      streamTess->SetMaxSubdivisionLevel( 2 );
      streamTess->SetFixedSubdivisions( 2 );
      vtkstd::cout << "# Using the generic Tessellator\n";
    }

  vtkGenericDataSetTessellator* tessellator = vtkGenericDataSetTessellator::New();

  tessellator->SetInput( mesh );
  tessellator->MergingOff();
  tessellator->Update(); //So that we can call GetRange() on the scalars
  
  vtkXMLUnstructuredGridWriter* writer=vtkXMLUnstructuredGridWriter::New();
  writer->SetInput( mesh );
  writer->SetFileName("QuadricVolume0.vtu");
  writer->SetDataModeToAscii();
  writer->Write();
  writer->Delete();

  tessellator->Delete();

  mesh->Delete();

  // TEST ELLIPSOID  ================================================
  vtkRegressionTest testE( "EllipsoidContours" );

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples( 8 );
  for (i=0; i<8; i++)
    {
    pts->SetTuple1( i, 4. );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 1 );
  for (i=0; i<19 ; i++)
    {
    dof->InsertNextRecord( 1, EScalarDOF + i );
    }

  sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents( 1 );
  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 3 );
  pts->SetNumberOfTuples( 8 );

  for (i=0; i<8; i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 3 );
  for (i=0; i<19; i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att = vtkShoeAttribute::New();
  att->SetName( "Geometry" );
  att->SetNumberOfComponents( 3 );
  att->SetPointData( pts );
  att->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  mesh = vtkShoeBox::New();
  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); 
  att->Delete();
  sca->Delete();

  sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, mesh->GetAttributes() );
  if ( ! sp )
    {
    test.StdOut() << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }

  att->GetOrder( sp->GetId(), od );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );

  contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetNumberOfContours( 3 );
  contour->SetValue( 0, .25 );
  contour->SetValue( 1, 2. );
  contour->SetValue( 2, 3. );
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();

  testE.StdOut() << "Ellipsoid contours have " 
                << contour->GetOutput()->GetNumberOfCells() 
                << " triangles."
                << vtkstd::endl;

  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "Ellipsoid.vtp" );
  pdw->Write();
  contour->Delete();

  // TEST ELLIPSOID CRITICAL POINTS ===========================================
  it = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator() );
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }

  it->Begin();
  sc = vtkShoeCell::SafeDownCast( it->GetCell() );

  vtkstd::cout << "#################################\n"
               << "       Ellipsoid\n"
               << "       ---------\n";
  sc->ComputeDOFCriticalPoints( sca );

  sca->SetRangeStyle( 3 );
  sca->GetRange( 0, range );
  vtkstd::cout << "## Tight range : [ " 
               << range[0] << " , " 
               << range[1] << " ]\n";
  vtkstd::cout << "#################################\n";

  it->Delete();
  mesh->Delete();

  // TEST 1-SHEET HYPERBOLOID  ================================================
  vtkRegressionTest testH1( "Hyperboloid-1SheetContours" );

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples( 8 );
  for (i=0; i<8; i++)
    {
    pts->SetTuple1( i, 1. );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 1 );
  for (i=0; i<19 ; i++)
    {
    dof->InsertNextRecord( 1, H1ScalarDOF + i );
    }

  sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents( 1 );
  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 3 );
  pts->SetNumberOfTuples( 8 );

  for (i=0; i<8; i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 3 );
  for (i=0; i<19; i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att = vtkShoeAttribute::New();
  att->SetName( "Geometry" );
  att->SetNumberOfComponents( 3 );
  att->SetPointData( pts );
  att->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  mesh = vtkShoeBox::New();
  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); 
  att->Delete();
  sca->Delete();

  sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, mesh->GetAttributes() );
  if ( ! sp )
    {
    test.StdOut() << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }

  att->GetOrder( sp->GetId(), od );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );

  contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetNumberOfContours( 3 );
  contour->SetValue( 0, 0. );
  contour->SetValue( 1, .25 );
  contour->SetValue( 2, 1. );
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();

  testH1.StdOut() << "1-sheet hyperboloid contours have " 
                << contour->GetOutput()->GetNumberOfCells() 
                << " triangles."
                << vtkstd::endl;

  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "Hyperboloid-1Sheet.vtp" );
  pdw->Write();
  contour->Delete();

  // TEST 1-SHEET HYPERBOLOID CRITICAL POINTS==================================
  it = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator() );
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }

  it->Begin();
  sc = vtkShoeCell::SafeDownCast( it->GetCell() );

  vtkstd::cout << "#################################\n"
               << "       1-Sheet Hyperboloid\n"
               << "       -------------------\n";
  sc->ComputeDOFCriticalPoints( sca );

  sca->SetRangeStyle( 3 );
  sca->GetRange( 0, range );
  vtkstd::cout << "## Tight range : [ " 
               << range[0] << " , " 
               << range[1] << " ]\n";
  vtkstd::cout << "#################################\n";

  it->Delete();
  mesh->Delete();

  // TEST 2-SHEET HYPERBOLOID  ================================================
  vtkRegressionTest testH2( "Hyberboloid-2SheetContours" );

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples( 8 );
  for (i=0; i<8; i++)
    {
    pts->SetTuple1( i, -1. );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 1 );
  for (i=0; i<19 ; i++)
    {
    dof->InsertNextRecord( 1, H2ScalarDOF + i );
    }

  sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents( 1 );
  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 3 );
  pts->SetNumberOfTuples( 8 );

  for (i=0; i<8; i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 3 );
  for (i=0; i<19; i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att = vtkShoeAttribute::New();
  att->SetName( "Geometry" );
  att->SetNumberOfComponents( 3 );
  att->SetPointData( pts );
  att->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  mesh = vtkShoeBox::New();
  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); 
  att->Delete();
  sca->Delete();

  sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, mesh->GetAttributes() );
  if ( ! sp )
    {
    test.StdOut() << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }

  att->GetOrder( sp->GetId(), od );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );

  contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetNumberOfContours( 3 );
  contour->SetValue( 0, 0. );
  contour->SetValue( 1, .25 );
  contour->SetValue( 2, .9 );
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();

  testH2.StdOut() << "2-sheet hyperboloid contours have " 
                << contour->GetOutput()->GetNumberOfCells() 
                << " triangles."
                << vtkstd::endl;

  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "Hyperboloid-2Sheet.vtp" );
  pdw->Write();
  contour->Delete();

  // TEST 2-SHEET HYPERBOLOID CRITICAL POINTS==================================
  it = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator() );
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }

  it->Begin();
  sc = vtkShoeCell::SafeDownCast( it->GetCell() );

  vtkstd::cout << "#################################\n"
               << "       2-Sheet Hyperboloid\n"
               << "       -------------------\n";
  sc->ComputeDOFCriticalPoints( sca );

  sca->SetRangeStyle( 3 );
  sca->GetRange( 0, range );
  vtkstd::cout << "## Tight range : [ " 
               << range[0] << " , " 
               << range[1] << " ]\n";
  vtkstd::cout << "#################################\n";

  it->Delete();
  mesh->Delete();

  // TEST ELLIPTIC PARABOLOID  ================================================
  vtkRegressionTest testEP( "EllipticParaboloidContours" );

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples( 8 );
  for (i=0; i<4; i++)
    {
    pts->SetTuple1( i, 3. );
    pts->SetTuple1( i + 4, 1. );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 1 );
  for (i=0; i<19 ; i++)
    {
    dof->InsertNextRecord( 1, EPScalarDOF + i );
    }

  sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents( 1 );
  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 3 );
  pts->SetNumberOfTuples( 8 );

  for (i=0; i<8; i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 3 );
  for (i=0; i<19; i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att = vtkShoeAttribute::New();
  att->SetName( "Geometry" );
  att->SetNumberOfComponents( 3 );
  att->SetPointData( pts );
  att->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  mesh = vtkShoeBox::New();
  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); 
  att->Delete();
  sca->Delete();

  sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, mesh->GetAttributes() );
  if ( ! sp )
    {
    test.StdOut() << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }

  att->GetOrder( sp->GetId(), od );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );

  contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetNumberOfContours( 3 );
  contour->SetValue( 0, 0. );
  contour->SetValue( 1, .25 );
  contour->SetValue( 2, 1. );
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();

  testEP.StdOut() << "Elliptic paraboloid contours have " 
                << contour->GetOutput()->GetNumberOfCells() 
                << " triangles."
                << vtkstd::endl;

  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "EllipticParaboloid.vtp" );
  pdw->Write();
  contour->Delete();

  // TEST ELLIPTIC PARABOLOID CRITICAL POINTS==================================
  it = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator() );
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }

  it->Begin();
  sc = vtkShoeCell::SafeDownCast( it->GetCell() );

  vtkstd::cout << "#################################\n"
               << "       Elliptic Paraboloid\n"
               << "       -------------------\n";
  sc->ComputeDOFCriticalPoints( sca );

  sca->SetRangeStyle( 3 );
  sca->GetRange( 0, range );
  vtkstd::cout << "## Tight range : [ " 
               << range[0] << " , " 
               << range[1] << " ]\n";
  vtkstd::cout << "#################################\n";

  it->Delete();
  mesh->Delete();

  // TEST HYPERBOLIC PARABOLOID  ==============================================
  vtkRegressionTest testHP( "HyperbolicParaboloidContours" );

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples( 8 );
  for (i=0; i<4; i++)
    {
    pts->SetTuple1( i, -1. );
    pts->SetTuple1( i + 4, 1. );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 1 );
  for (i=0; i<19 ; i++)
    {
    dof->InsertNextRecord( 1, HPScalarDOF + i );
    }

  sca = vtkShoeAttribute::New();
  sca->SetName( "Scalar" );
  sca->SetNumberOfComponents( 1 );
  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 3 );
  pts->SetNumberOfTuples( 8 );

  for (i=0; i<8; i++)
    {
    pts->SetTuple( i, ShoeAttributePts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 3 );
  for (i=0; i<19; i++)
    {
    dof->InsertNextRecord( 1, ShoeAttributeDOF + 3*i );
    }

  att = vtkShoeAttribute::New();
  att->SetName( "Geometry" );
  att->SetNumberOfComponents( 3 );
  att->SetPointData( pts );
  att->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  mesh = vtkShoeBox::New();
  mesh->SetGeometry( att );
  mesh->GetAttributes()->InsertNextAttribute( sca );
  mesh->GetAttributes()->SetActiveAttribute( mesh->AttributeId( sca ) );
  mesh->UpdateLinks(); 
  att->Delete();
  sca->Delete();

  sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, mesh->GetAttributes() );
  if ( ! sp )
    {
    test.StdOut() << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }

  att->GetOrder( sp->GetId(), od );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );

  contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetNumberOfContours( 3 );
  contour->SetValue( 0, -1.5 );
  contour->SetValue( 1, 0. );
  contour->SetValue( 2, 1. );
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();

  testHP.StdOut() << "Hyperbolic paraboloid contours have " 
                << contour->GetOutput()->GetNumberOfCells() 
                << " triangles."
                << vtkstd::endl;

  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "HyperbolicParaboloid.vtp" );
  pdw->Write();
  contour->Delete();

  // TEST HYPERBOLIC PARABOLOID CRITICAL POINTS ================================
  it = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator() );
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }

  it->Begin();
  sc = vtkShoeCell::SafeDownCast( it->GetCell() );

  vtkstd::cout << "#################################\n"
               << "       Hyperbolic Paraboloid\n"
               << "       ---------------------\n";
  sc->ComputeDOFCriticalPoints( sca );

  sca->SetRangeStyle( 3 );
  sca->GetRange( 0, range );
  vtkstd::cout << "## Tight range : [ " 
               << range[0] << " , " 
               << range[1] << " ]\n";
  vtkstd::cout << "#################################\n";

  // TEST HYPERBOLIC PARABOLOID ISOCONTOURS ================================

  vtkstd::cout << "#################################\n"
               << "       Higer Order Isocontouring\n"
               << "       ---------------------\n";

  vtkShoeBoxContourFilter* iso = vtkShoeBoxContourFilter::New();
  iso->SetInput( mesh );
  iso->SetFieldId( mesh->GetAttributes()->FindAttribute( sca->GetName() ) );
  iso->SetFieldComponent( 0 );
  iso->SetMaximumNumberOfSubdivisions( 0 );

  vtkUnstructuredGridWriter* ugw = vtkUnstructuredGridWriter::New();
  ugw->SetInput( iso->GetOutput() );

  vtkstd::cout << "## -1.5-isocontour:\n";
  iso->SetIsovalue( -1.5 );
  ugw->SetFileName( "hyperbolic_paraboloid_m1_5.vtk" );
  ugw->Write();

  vtkstd::cout << "## 0-isocontour:\n";
  iso->SetIsovalue( 0. );
  ugw->SetFileName( "hyperbolic_paraboloid_0.vtk" );
  ugw->Write();

  vtkstd::cout << "## 1.1-isocontour:\n";
  iso->SetIsovalue( 1.1 );
  ugw->SetFileName( "hyperbolic_paraboloid_1_1.vtk" );
  ugw->Write();

  ugw->Delete();
  iso->Delete();
  it->Delete();
  mesh->Delete();

  vtkstd::cout << "#################################\n";

  // ==========================================================================

  pdw->Delete();

  return 0;
}
