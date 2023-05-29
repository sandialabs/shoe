// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeBoxReader.h"

#include "vtkIndent.h"

#include "vtkRegressionTest.h"
#include "vtkShoeBox.h"

int ShoeBoxReader( int argc, char* argv[] )
{
  char* testFile = vtkRegressionTest::ExpandDataFileName(
        argc, argv, VTKSNL_DATA_DIRECTORY, "samp_shoe2.shoe" );
  vtkstd::cout << "Reading \"" << testFile << "\"" << vtkstd::endl;
  vtkShoeBoxReader* rdr = vtkShoeBoxReader::New();
  rdr->SetFileName( testFile );
  rdr->UpdateInformation();
  vtkstd::cout
    << rdr->GetNumberOfPoints() << " points, "
    << rdr->GetNumberOfDofNodes() << " DOF nodes, "
    << rdr->GetNumberOfCells() << " cells.\n";
  rdr->SetPointArrayStatus( "Scalar", 0 );
  for ( int a = 0; a < rdr->GetNumberOfPointArrays(); ++a )
    {
    vtkstd::cout << "Attribute \"" << rdr->GetPointArrayName( a ) << "\" [" << rdr->GetPointArrayStatus( a ) << "]\n";
    }

  rdr->Update();
  vtkIndent indent;
  vtkShoeBox* mesh = rdr->GetOutput();
  vtkstd::cout << "Resulting mesh has " << mesh->GetNumberOfCells(-1) << " cells." << vtkstd::endl;

  rdr->Delete();
  return 0;
}

