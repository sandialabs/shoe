/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#define DEBUG_HEX
#define TINY
#undef KITWARE_TESS
#define REFINE_FIRST

#include <vtkHexahedron.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPlane.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericContourFilter.h>
#include <vtkGenericAttributeCollection.h>
#include <vtkSimpleCellTessellator.h>
#include <vtkGenericDataSetTessellator.h>
#include <vtkGenericSubdivisionErrorMetric.h>
#include <vtkGeometricErrorMetric.h>
#include <vtkAttributesErrorMetric.h>
#include <vtkContourFilter.h>
#include <vtkCollection.h>

#include "vtkDataRecords.h"
#include "vtkDataRecordsIterator.h"
#include "vtkShoeEnums.h"
#include "vtkShoeAttribute.h"
#include "vtkShoeCellGenus.h"
#include "vtkShoeCellIterator.h"
#include "vtkShoeCellSpecies.h"
#include "vtkShoeOrderTuple.h"
#include "vtkShoeBox.h"
#include "vtkGenericStreamingTessellator.h"


vtkIdType NumberOfCornerNodes = 0;
vtkIdType NumberOfDOFNodes = 0;
vtkIdType NumberOfElements = 0;

vtkIdType* CornerIndex;
vtkDoubleArray* CornerNodes = vtkDoubleArray::New();
vtkDoubleArray* CornerValues = vtkDoubleArray::New();
vtkDoubleArray* CornerVectors = vtkDoubleArray::New();

vtkDataRecords* DOFNodes = vtkDataRecords::Create( VTK_DOUBLE );
vtkDataRecords* DOFValues = vtkDataRecords::Create( VTK_DOUBLE );
vtkDataRecords* DOFVectors = vtkDataRecords::Create( VTK_DOUBLE );

#ifdef DEBUG_HEX
vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
vtkPoints* ugpts = vtkPoints::New();
#endif // 0

vtkShoeBox* sb = 0;
vtkShoeAttribute* sga = 0;
vtkShoeAttribute* ssa = 0;
vtkShoeAttribute* sva = 0;
vtkShoeCellSpecies* sp = 0;

void ReadPoints( ifstream &fcor, ifstream &fdof )
{
  double p[3],f,v[3];
  vtkIdType i;
  int id;

  fcor >> NumberOfCornerNodes;
  cout << "Reading " << NumberOfCornerNodes << " corner nodes and attached values...";
  cout.flush();
  
  CornerNodes->SetNumberOfComponents( 3 );
  CornerNodes->SetNumberOfTuples( NumberOfCornerNodes );

  CornerValues->SetName( "VorticityMagnitude" );
  CornerValues->SetNumberOfComponents( 1 );
  CornerValues->SetNumberOfTuples( NumberOfCornerNodes );

  CornerVectors->SetName( "Velocity" );
  CornerVectors->SetNumberOfComponents( 3 );
  CornerVectors->SetNumberOfTuples( NumberOfCornerNodes );

  for ( i=0; i<NumberOfCornerNodes; i++ )
    {
    fcor >> p[0];
    fcor >> p[1];
    fcor >> p[2];
    fcor >> v[0];
    fcor >> v[1];
    fcor >> v[2];
    fcor >> f   ;
    fcor >> id  ;
    CornerNodes->SetTuple( i, p );
    CornerValues->SetTuple( i, &f );
    CornerVectors->SetTuple( i, v );
    }

  cout << " done." << vtkstd::endl;

  fdof >> NumberOfDOFNodes;
  cout << "Reading " << NumberOfDOFNodes << " DOF nodes and attached values...";
  cout.flush();
  
  DOFNodes->SetNumberOfComponents( 3 );
  DOFNodes->AllocateRecords( NumberOfDOFNodes );

  //  DOFValues->SetName( "VorticityMagnitude" );
  DOFValues->SetNumberOfComponents( 1 );
  DOFValues->AllocateRecords( NumberOfDOFNodes );

  //  DOFVectors->SetName( "Velocity" );
  DOFVectors->SetNumberOfComponents( 3 );
  DOFVectors->AllocateRecords( NumberOfDOFNodes );

  for ( i=0; i<NumberOfDOFNodes; i++ )
    {
    fdof >> p[0];
    fdof >> p[1];
    fdof >> p[2];
    fdof >> v[0];
    fdof >> v[1];
    fdof >> v[2];
    fdof >> f   ;
    fdof >> id  ;
    DOFNodes->InsertNextRecord( 1, p );
    DOFValues->InsertNextRecord( 1, &f );
    DOFVectors->InsertNextRecord( 1, v );
    }

  cout << " done." << vtkstd::endl;

  return;
}

void ReadHexahedraLagrangeQuadraticTensor( ifstream &fhex )
{
  vtkIdType conn[27],i,n;

  fhex >> NumberOfElements;
  cout << "Reading " << NumberOfElements << " quadratic Lagrange hexahedra...";
  cout.flush();

  int j;
  bool consecutive = true;

  for ( i=0; i<NumberOfElements; i++ )
    {
    fhex >> n ;

    if ( ! consecutive && i + 1 != n )
      {
      cout << vtkstd::endl << "WAR: non-consecutive indexing after element " 
	   << i << "." << vtkstd::endl;
      consecutive = false;
      }

    for ( j=0; j<27 ; j++ )
      {
      fhex >> conn[j];
      }
#ifdef DEBUG_HEX
    ug->InsertNextCell( VTK_HEXAHEDRON, 8, conn );
#endif // 0
    sb->InsertNextCell( sp->GetId(), conn, 0 /*permutation*/ );
    }


  cout << " done." << vtkstd::endl;

  return;
}

int main( int argc, char** argv )
{
  sb = vtkShoeBox::New();
  sga = vtkShoeAttribute::New();
  ssa = vtkShoeAttribute::New();
  sva = vtkShoeAttribute::New();

  ifstream fcor( "corners.dat" );
  ifstream fdof( "dof.dat" );
  ReadPoints( fcor, fdof );

#ifdef DEBUG_HEX  
  ugpts->SetData( CornerNodes );
  ug->SetPoints( ugpts );
  ug->GetPointData()->SetScalars( CornerValues );
  ug->GetPointData()->SetVectors( CornerVectors );
#endif // 0

  sga->SetName( "Geometry" );
  sga->SetNumberOfComponents( 3 );
  sga->SetPointData( CornerNodes );
  sga->SetDOFData( DOFNodes );

  ssa->SetName( "VorticityMagnitude" );
  ssa->SetNumberOfComponents( 1 );
  ssa->SetPointData( CornerValues );
  ssa->SetDOFData( DOFValues );

  sva->SetName( "Velocity" );
  sva->SetNumberOfComponents( 3 );
  sva->SetPointData( CornerVectors );
  sva->SetDOFData( DOFVectors );

  sb->SetGeometry( sga );
  sb->GetAttributes()->InsertNextAttribute( ssa );
  sb->GetAttributes()->SetActiveAttribute( sb->AttributeId( ssa ) );
  sb->UpdateLinks();

  vtkShoeCellGenus g;
  g.Shape = shoe::Hexahedron;
  vtkPolyInterpolant in[2] = { Lagrange, Lagrange };
  vtkPolyProductSpace ps[2] = { Tensor, Tensor };
  vtkShoeOrderTuple o[2];
  o[0].Set( 2, 2, 2 );
  o[1].Set( 2, 2, 2 );
  int oi[2] = { 0, 1 };
  int n = 2;
  sp = vtkShoeCellSpecies::FindOrCreate( g, in, ps, o, oi, n, sb->GetAttributes() );
   if ( ! sp )
    {
    cout << "Unable to find or create a cell species." << vtkstd::endl;
    }

  ifstream fhex( "hex27.dat" );
  ReadHexahedraLagrangeQuadraticTensor( fhex );  

#ifdef DEBUG_HEX
  vtkXMLUnstructuredGridWriter* wr = vtkXMLUnstructuredGridWriter::New();
  wr->SetInput( ug );
  wr->SetFileName( "linearized_hex27.vtu" );
  wr->Write();
#endif // 0

#ifdef KITWARE_TESS
  vtkSimpleCellTessellator* tess = vtkSimpleCellTessellator::New();
  tess->SetMaxSubdivisionLevel(10);
  tess->SetSubdivisionLevels(1,10);
#else
  vtkGenericStreamingTessellator* tess = vtkGenericStreamingTessellator::New();
  tess->SetMaximumNumberOfSubdivisions(4);
#endif
  tess->Initialize( sb );
  sb->SetTessellator( tess );

#ifdef REFINE_FIRST
  vtkGenericDataSetTessellator* dst = vtkGenericDataSetTessellator::New();
  dst->SetInput( sb );

  vtkGeometricErrorMetric* geometricError = vtkGeometricErrorMetric::New();
  geometricError->SetRelativeGeometricTolerance( 0.01, sb );
  sb->GetTessellator()->GetErrorMetrics()->AddItem( geometricError );
  geometricError->Delete();

  vtkAttributesErrorMetric* attributesError = vtkAttributesErrorMetric::New();
  attributesError->SetAttributeTolerance( 0.01 );
  sb->GetTessellator()->GetErrorMetrics()->AddItem( attributesError );
  attributesError->Delete();

  vtkContourFilter* contour = vtkContourFilter::New();
  contour->SetInput( dst->GetOutput() );

#else // REFINE_FIRST !
  vtkGenericContourFilter *contour = vtkGenericContourFilter::New();
  contour->SetInput( sb );
#endif // ! REFINE_FIRST

  contour->SetNumberOfContours( 3 );
#ifdef TINY
  contour->SetValue( 0, 0. );
  contour->SetValue( 1, .25 );
  contour->SetValue( 2, 1. );
#else
  contour->SetValue( 0, 150000 );
  contour->SetValue( 1, 272500 );
  contour->SetValue( 2, 617500 );
#endif // TINY
  contour->ComputeGradientsOn();
  contour->ComputeNormalsOn();
  contour->ComputeScalarsOn();
  contour->Update();
  vtkstd::cout << "Contour has " 
	       << contour->GetOutput()->GetNumberOfCells() 
	       << " triangles."
	       << vtkstd::endl;

  vtkXMLPolyDataWriter* pdw = vtkXMLPolyDataWriter::New();
  pdw->SetInput( contour->GetOutput() );
  pdw->SetFileName( "ShoeBox.vtp" );
  pdw->Write();
                                                                                               

  return 0;
}
