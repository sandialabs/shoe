// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeEnums.h"
#include "vtkShoeAttribute.h"
#include "vtkGenericStreamingTessellator.h"
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

#include <vtkCollection.h>
#include <vtkGenericCellTessellator.h>
#include <vtkGenericDataSetTessellator.h>
#include <vtkGenericSubdivisionErrorMetric.h>
#include <vtkAttributesErrorMetric.h>
#include <vtkGeometricErrorMetric.h>
#include <vtkGenericAdaptorCell.h>
#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkGenericAttributeCollection.h>
#include <vtkGenericContourFilter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSimpleCellTessellator.h>
#include <vtkTimerLog.h>

static double ShoeAttributePts[] =
{
-.7, -1., -1.,
 1., -1., -1. ,
 1.,  1., -1. ,
-.7,  1., -1. ,
-.7, -1.,  1. ,
 1., -1.,  1. ,
 1.,  1.,  1. ,
-.7,  1.,  1. ,
 3.5, -1., -1. ,
 3.5,  1., -1. ,
 3.5, -1.,  1. ,
 3.5,  1.,  1. ,
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

 -.5,  0.,  0. , // 12
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

static vtkIdType ShoeTestCellConn[] =
{
  0,  1,  2,  3,  4,  5,  6,  7,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
  1,  8,  9,  2,  5, 10, 11,  6, 19, 20, 21,  1, 22, 23, 24,  5,  9, 25, 11, 26, 13, 27, 28, 29, 30, 31, 32
};

int main( int argc, char** argv )
{
  // CREATE TIMER LOG ========================================================
  vtkTimerLog* tl = vtkTimerLog::New();

  // DEFINE VARIABLES NEEDED FOR TESTING ======================================
  vtkShoeBox* mesh = vtkShoeBox::New();
  vtkShoeAttribute* att = vtkShoeAttribute::New();
  vtkShoeAttribute* sca;
  vtkShoeCellGenus g;

  int i;
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

  // TEST SIMPLE TESSELLATOR ==================================================
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
      vtkstd::cout << "Unable to find or create a cell species that matches all the criteria:" << vtkstd::endl
      << "  Shape: " << g.Shape << vtkstd::endl
      << "  Number of attributes to match: " << n << vtkstd::endl
      << "Perhaps you've specified an interpolant not implemented?" << vtkstd::endl;
    }
  vtkShoeOrderTuple od;
  att->GetOrder( sp->GetId(), od );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn, 0 /*permutation*/ );
  id = mesh->InsertNextCell( sp->GetId(), ShoeTestCellConn + 27, 0 /*permutation*/ );

  // Set the error metric thresholds:
  // 1. for the geometric error metric
  vtkGeometricErrorMetric *geometricError=vtkGeometricErrorMetric::New();
  geometricError->SetRelativeGeometricTolerance( 0.01, mesh );

  mesh->GetTessellator()->GetErrorMetrics()->AddItem( geometricError );
  geometricError->Delete();

  // 2. for the attribute error metric
  vtkAttributesErrorMetric *attributesError=vtkAttributesErrorMetric::New();
  attributesError->SetAttributeTolerance( 0.01 ); // 0.11, 0.005

  mesh->GetTessellator()->GetErrorMetrics()->AddItem( attributesError );
  attributesError->Delete();

  vtkGenericDataSetTessellator* tessellator = vtkGenericDataSetTessellator::New();
  vtkXMLUnstructuredGridWriter* ugw = vtkXMLUnstructuredGridWriter::New();

  // TEST GENERIC STREAMING TESSELLATOR =====================================
  vtkstd::cout << mesh->GetTessellator()->GetClassName()
	       << " : ";
  vtkGenericStreamingTessellator* gst = vtkGenericStreamingTessellator::SafeDownCast( mesh->GetTessellator() );
  if ( gst ) 
    {
    gst->SetFixedSubdivisions( 4 );
    gst->SetMaximumNumberOfSubdivisions( 6 );
    }

  tessellator->SetInput( mesh );

  tl->StartTimer();
  tessellator->Update();
  tl->StopTimer();
  vtkstd::cout << tl->GetElapsedTime()
	       << vtkstd::endl;

  ugw->SetFileName( "GenericStreamingTessellatorBenchmark.vtu" );
  ugw->SetInput( tessellator->GetOutput() );
  ugw->Write();
  
  // TEST SIMPLE CELL TESSELLATOR ===========================================
  mesh->SetTessellator( vtkSimpleCellTessellator::New() );
  vtkstd::cout << mesh->GetTessellator()->GetClassName()
	       << " : ";
  vtkSimpleCellTessellator* sct = vtkSimpleCellTessellator::SafeDownCast( mesh->GetTessellator() );
  if ( sct ) sct->SetSubdivisionLevels( 4, 100 );

  tessellator->SetInput( mesh );

  tl->StartTimer();
  tessellator->Update();
  tl->StopTimer();
  vtkstd::cout << tl->GetElapsedTime()
	       << vtkstd::endl;

  ugw->SetFileName( "SimpleCellTessellatorBenchmark.vtu" );
  ugw->SetInput( tessellator->GetOutput() );
  ugw->Write();

  // CLEAN UP =================================================================
  mesh->Delete();
  ugw->Delete();
  tessellator->Delete();

  return 0;
}
