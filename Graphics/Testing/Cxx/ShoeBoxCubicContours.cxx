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

static double CubicHexPts[] =
{
-1., -1., -1. ,
 1., -1., -1. ,
 1.,  1., -1. ,
-1.,  1., -1. ,
-1., -1.,  1. ,
 1., -1.,  1. ,
 1.,  1.,  1. ,
-1.,  1.,  1.
};

static double CubicHexDOF[] =
{
-0.333333, -1, -1,
 0.333333, -1, -1,
 1, -0.333333, -1,
 1, 0.333333, -1,
-0.333333, 1, -1,
 0.333333, 1, -1,
-1, -0.333333, -1,
-1, 0.333333, -1,
-0.333333, -1, 1,
 0.333333, -1, 1,
 1, -0.333333, 1,
 1, 0.333333, 1,
-0.333333, 1, 1,
 0.333333, 1, 1,
-1, -0.333333, 1,
-1, 0.333333, 1,
-1, -1, -0.333333,
-1, -1, 0.333333,
 1, -1, -0.333333,
 1, -1, 0.333333,
-1, 1, -0.333333, // Kitwareism ( 10 <-> 11 )
-1, 1, 0.333333,
 1, 1, -0.333333,
 1, 1, 0.333333,
-1, -0.333333, -0.333333,
-1, 0.333333, -0.333333,
-1, -0.333333, 0.333333,
-1, 0.333333, 0.333333,
 1, -0.333333, -0.333333,
 1, 0.333333, -0.333333,
 1, -0.333333, 0.333333,
 1, 0.333333, 0.333333,
-0.333333, -1, -0.333333,
 0.333333, -1, -0.333333,
-0.333333, -1, 0.333333,
 0.333333, -1, 0.333333,
-0.333333, 1, -0.333333,
 0.333333, 1, -0.333333,
-0.333333, 1, 0.333333,
 0.333333, 1, 0.333333,
-0.333333, -0.333333, -1,
 0.333333, -0.333333, -1,
-0.333333, 0.333333, -1,
 0.333333, 0.333333, -1,
-0.333333, -0.333333, 1,
 0.333333, -0.333333, 1,
-0.333333, 0.333333, 1,
 0.333333, 0.333333, 1,
-0.333333, -0.333333, -0.333333,
 0.333333, -0.333333, -0.333333,
-0.333333, 0.333333, -0.333333,
 0.333333, 0.333333, -0.333333,
-0.333333, -0.333333, 0.333333,
 0.333333, -0.333333, 0.333333,
-0.333333, 0.333333, 0.333333,
 0.333333, 0.333333, 0.333333
};

static double DDScalarDOF[] =
{ 
// Edge values:
-0.888889,
-0.888889,
-0.888889,
-0.888889,
-0.888889,
-0.888889,
-0.888889,
-0.888889,
1.11111,
1.11111,
1.11111,
1.11111,
1.11111,
1.11111,
1.11111,
1.11111,
1.85185,
1.92593,
1.85185,
1.92593,
1.85185,
1.92593,
1.85185,
1.92593,
// Face values:
0.962963,
0.962963,
1.03704,
1.03704,
0.962963,
0.962963,
1.03704,
1.03704,
0.962963,
0.962963,
1.03704,
1.03704,
0.962963,
0.962963,
1.03704,
1.03704,
-1.77778,
-1.77778,
-1.77778,
-1.77778,
0.222222,
0.222222,
0.222222,
0.222222,
// Body values:
0.0740741,
0.0740741,
0.0740741,
0.0740741,
0.148148,
0.148148,
0.148148,
0.148148
};

static vtkIdType ShoeTestCellConn[] =
{
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18
};

#if 0
static void PolynomialSystemPrintRoots( vtkPolynomialSystem* ps, vtkstd::ostream& os )
{
  os << ps->GetNumberOfRoots() << " total roots and ";
  int nr = ps->GetNumberOfRoots();
  os << nr << " real roots.\n";
  int i;
  int nv = ps->GetVariableOrder( 0 );
  double* root = new double[ nv ];
  for ( i=0; i<nr; ++i )
    {
    os << "    " << i << ": (";
    ps->GetRoot( i, root );
    for ( int v=0; v<nv; ++v )
      {
      os << " " << root[v];
      }
    os << " )\n";
    }
  delete [] root;
}
#endif // 0

int ShoeBoxCubicContours( int vtkNotUsed(argc), char** vtkNotUsed(argv) )
{
  vtkRegressionTest test( "ShoeBoxCubicContours" );

  // DEFINE VARIABLES NEEDED FOR TESTING ================================
  vtkShoeBox* mesh;
  vtkShoeAttribute* att;
  vtkShoeAttribute* sca;
  vtkShoeCellGenus g;
  vtkGenericContourFilter* contour;
  vtkXMLPolyDataWriter* pdw = vtkXMLPolyDataWriter::New();

  int i;
  double range[2];
  vtkIdType id;
  vtkDoubleArray* pts;
  vtkDataRecords* dof;

  // TEST DING-DONG SURFACE ==============================================
  vtkRegressionTest testDD( "DingDongSurface" );

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 1 );
  pts->SetNumberOfTuples( 8 );
  for (i=0; i<4; i++)
    {
    pts->SetTuple1( i, 0. );
    pts->SetTuple1( i + 4, 2. );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 1 );
  for (i=0; i<12 ; i++)
    {
    dof->InsertNextRecord( 2, DDScalarDOF + 2*i );
    }
  for (i=0; i<6 ; i++)
    {
    dof->InsertNextRecord( 4, DDScalarDOF + 4*i + 24 );
    }
  dof->InsertNextRecord( 8, DDScalarDOF + 48 );

  sca = vtkShoeAttribute::New();
  sca->SetPointData( pts );
  sca->SetDOFData( dof );
  pts->Delete();
  dof->Delete();

  pts = vtkDoubleArray::New();
  pts->SetNumberOfComponents( 3 );
  pts->SetNumberOfTuples( 8 );

  for (i=0; i<8; i++)
    {
    pts->SetTuple( i, CubicHexPts + 3*i );
    }

  dof = vtkDataRecords::New();
  dof->SetNumberOfComponents( 3 );
  for (i=0; i<12; i++)
    {
    dof->InsertNextRecord( 2, CubicHexDOF + 6*i );
    }
  for (i=0; i<6 ; i++)
    {
    dof->InsertNextRecord( 4, CubicHexDOF + 12*i + 72 );
    }
  dof->InsertNextRecord( 8, CubicHexDOF + 144 );

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
  sca->Delete();
  att->Delete();

  g.Shape = shoe::Hexahedron;
  vtkPolyInterpolant in[2] = { Lagrange, Lagrange };
  vtkPolyProductSpace ps[2] = { Tensor, Tensor };
  vtkShoeOrderTuple o[2];
  o[0].Set( 3, 3, 3 );
  o[1].Set( 3, 3, 3 );
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
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );

  // TEST DING-DONG VOLUME CRITICAL POINTS ==============================
#if 0
  vtkShoeCellIterator* it = vtkShoeCellIterator::SafeDownCast( mesh->NewCellIterator() );
  if ( ! it )
    {
    test.StdOut() << "Error: Iterator is not a vtkShoeCellIterator." << vtkstd::endl;
    }

  it->Begin();
  vtkShoeCell* sc = vtkShoeCell::SafeDownCast( it->GetCell() );

  vtkstd::cout << "#################################\n"
               << "       Ding Dong Function\n"
               << "       ------------------\n";
  sc->ComputeDOFCriticalPoints( sca );

  sca->SetRangeStyle( 3 );
  sca->GetRange( 0, range );
  vtkstd::cout << "## Tight range : [ " 
               << range[0] << " , " 
               << range[1] << " ]\n";
  vtkstd::cout << "#################################\n";

  it->Delete();
#endif

  // TEST DING-DONG SURFACE SAMPLING + CONTOURING ============================
  vtkSimpleCellTessellator* simpleTess = vtkSimpleCellTessellator::SafeDownCast( mesh->GetTessellator() );
  if ( simpleTess )
    {
      simpleTess->SetMaxSubdivisionLevel( 3 );
      simpleTess->SetFixedSubdivisions( 3 );
      vtkstd::cout << "# Using the simple Tessellator\n";
    }
  vtkGenericStreamingTessellator* streamTess = vtkGenericStreamingTessellator::SafeDownCast( mesh->GetTessellator() );
  if ( streamTess )
    {
      streamTess->SetMaxSubdivisionLevel( 3 );
      streamTess->SetFixedSubdivisions( 3 );
      vtkstd::cout << "# Using the generic Tessellator\n";
    }

  contour = vtkGenericContourFilter::New();
  contour->SetInput( mesh );
  contour->SetNumberOfContours( 1 );
  contour->SetValue( 0, 0. );
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();

  testDD.StdOut() << "Ding dong surface contour has " 
		<< contour->GetOutput()->GetNumberOfCells() 
		<< " triangles."
		<< vtkstd::endl;

  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "DingDongSurface.vtp" );
  pdw->Write();

  // TEST DING-DONG SURFACE HO ISOCONTOURING  ===================================
  vtkstd::cout << "#################################\n"
               << "       High Order Isocontouring of the Ding Dong Function\n"
               << "       ------------------\n";
  vtkShoeBoxContourFilter* iso = vtkShoeBoxContourFilter::New();
  iso->SetInput( mesh );
  iso->SetFieldId( mesh->GetAttributes()->FindAttribute( sca->GetName() ) );
  iso->SetFieldComponent( 0 );
  iso->SetMaximumNumberOfSubdivisions( 0 );
  iso->SetIsovalue( 0. );

  vtkUnstructuredGridWriter* ugw = vtkUnstructuredGridWriter::New();
  ugw->SetInput( iso->GetOutput() );
  ugw->SetFileName( "DingDongContour.vtk" );
  ugw->Write();
  iso->Delete();
  ugw->Delete();
  vtkstd::cout << "#################################\n";

  // TEST DING-DONG VOLUME RENDERING ===================================
  vtkGenericDataSetTessellator* tessellator = vtkGenericDataSetTessellator::New();

  tessellator->SetInput( mesh );
  tessellator->Update(); //So that we can call GetRange() on the scalars
  
  vtkXMLUnstructuredGridWriter* writer=vtkXMLUnstructuredGridWriter::New();
  writer->SetInput(tessellator->GetOutput());
  writer->SetFileName("DingDongVolume.vtu");
  writer->SetDataModeToAscii();
  writer->Write();
  writer->Delete();

  pdw->Delete();
  contour->Delete();

  tessellator->Delete();
  mesh->Delete();

  return 0;
}
