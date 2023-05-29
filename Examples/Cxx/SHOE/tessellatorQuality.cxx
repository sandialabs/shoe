/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "vtkToySubdivision.h"
#include "vtkINRIAMeshReader.h"
#include "vtkAdaptiveTessellator.h"
#include "vtkDataSetPolygonizer.h"
#include "vtkMeshQuality.h"

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPlane.h>
#include <vtkFieldData.h>

int main( int argc, char** argv )
{
  vtkINRIAMeshReader* mr = vtkINRIAMeshReader::New();
  vtkDataSetPolygonizer* dp = vtkDataSetPolygonizer::New();
  vtkToySubdivision* ts = vtkToySubdivision::New();
  vtkPlane* pl = vtkPlane::New();
  vtkUnstructuredGrid* ug;
  vtkPoints* pt;
  double av,bds[6];

  mr->SetFileName( argc > 1 ? argv[1] : VTKSNL_DATA_DIRECTORY VTKSNL_DIRSEP "samp_cuboid.mesh" );
        mr->SetReadCorners( 0 );
        mr->SetReadEdges( 0 );
        mr->SetReadTriangles( 0 );
        mr->SetReadTetrahedra( 1 );
  mr->Update();

  ug = mr->GetOutput();
  pt = ug->GetPoints();
  pt->GetBounds( bds );
  for ( int i=0; i<3; ++i )
    {
    bds[2*i] += bds[2*i+1];
    bds[2*i] /= 2.;
    }

  //pl->SetNormal( 1., 0., 0. );
  pl->SetNormal( sqrt(3.)/3., sqrt(3.)/3., sqrt(3.)/3. );
  pl->SetOrigin( bds[0], bds[2], bds[4] );
  ts->SetRefiner( pl );
  ts->SetMesh( ug );

  vtkMeshQuality* iq = vtkMeshQuality::New();
  iq->SetTetQualityMeasureToAspectRatio();
  iq->SetInput( ug );
  iq->Update();
  av = iq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 1 );
  cout << "Initial Mesh Triangle Quality: "
       << iq->GetOutput()->GetFieldData()->GetArray( "Mesh Triangle Quality" )->GetComponent( 0, 0 )
       << endl;
  cout << "Initial Mesh Quadrilateral Quality: "
       << iq->GetOutput()->GetFieldData()->GetArray( "Mesh Quadrilateral Quality" )->GetComponent( 0, 0 )
       << endl;
  cout << "Initial Best Mesh Tetrahedron Quality: "
       << iq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 0 )
       << endl;
  cout << "Initial Average Mesh Tetrahedron Quality: " << av
       << endl;
  cout << "Initial Worst Mesh Tetrahedron Quality: "
       << iq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 2 )
       << endl;
  cout << "Initial Mesh Tetrahedron Quality Standard Deviation: "
       << sqrt(iq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 3 ) - av * av)
       << endl;
  cout << "Initial Mesh Hexahedron Quality: "
       << iq->GetOutput()->GetFieldData()->GetArray( "Mesh Hexahedron Quality" )->GetComponent( 0, 0 )
       << endl;
  
  dp->SetInput( iq->GetOutput() );
  dp->SetSubdivider( ts );
  dp->GetTessellator()->SetMaximumNumberOfSubdivisions( argc > 3 ? atoi(argv[3]) : 3 );
  dp->Update();
  
  vtkMeshQuality* fq = vtkMeshQuality::New();
  fq->SetTetQualityMeasureToAspectRatio();
  fq->SetInput( dp->GetOutput() );
  fq->Update();
  av = fq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 1 );
  cout << "Final Mesh Triangle Quality: "
       << fq->GetOutput()->GetFieldData()->GetArray( "Mesh Triangle Quality" )->GetComponent( 0, 0 )
       << endl;
  cout << "Final Mesh Quadrilateral Quality: "
       << fq->GetOutput()->GetFieldData()->GetArray( "Mesh Quadrilateral Quality" )->GetComponent( 0, 0 )
       << endl;
  cout << "Final Best Mesh Tetrahedron Quality: "
       << fq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 0 )
       << endl;
  cout << "Final Average Mesh Tetrahedron Quality: " << av
       << endl;
  cout << "Final Worst Mesh Tetrahedron Quality: "
       << fq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 2 )
       << endl;
  cout << "Final Mesh Tetrahedron Quality Standard Deviation: "
       << sqrt(fq->GetOutput()->GetFieldData()->GetArray( "Mesh Tetrahedron Quality" )->GetComponent( 0, 3 ) - av * av)
       << endl;
  cout << "Final Mesh Hexahedron Quality: "
       << fq->GetOutput()->GetFieldData()->GetArray( "Mesh Hexahedron Quality" )->GetComponent( 0, 0 )
       << endl;
  
#ifdef DBG_CASE_COUNTS
  int c, s;
  for ( c = 0; c < 11; ++c )
    {
    cout << dp->GetTessellator()->GetCaseCount( c );
    for ( s = 0; s < 51; ++s )
      {
      cout << " " << dp->GetTessellator()->GetSubcaseCount( c, s );
      }
    cout << endl;
    }
#endif // DBG_CASE_COUNTS

  vtkXMLUnstructuredGridWriter* wr = vtkXMLUnstructuredGridWriter::New();
  //wr->SetDataModeToAscii();
  wr->SetInput( dp->GetOutput() );
  wr->SetFileName( argc > 2 ? argv[2] : "samp_inria_subdivided.vtu" );
  wr->Write();

  return 0;
}
