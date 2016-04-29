// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeCellMetaData.h"
#include "Elements/ShoeHexahedron.h"
#include "Elements/ShoeInterpolants.h"

#include "vtkSystemIncludes.h"
#include "vtkCellType.h"
#include "vtkDoubleArray.h"
#include "vtkObjectFactory.h"

vtkstd::set<vtkShoeCellMetaData*> vtkShoeCellMetaData::All;
const vtkShoeCellMetaData* vtkShoeCellMetaData::ByShape[shoe::NumberOfCellShapes+1];

vtkCxxRevisionMacro(vtkShoeCellMetaData,"$Revision: 10018 $");
vtkStandardNewMacro(vtkShoeCellMetaData);

void vtkShoeCellMetaData::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Type: " << this->Type << vtkstd::endl;
  os << indent << "Dimension:                " << this->Dimension << vtkstd::endl;
  os << indent << "NumberOfBoundaries:       "
    << this->NumberOfBoundaries[0] << ", "
    << this->NumberOfBoundaries[1] << ", "
    << this->NumberOfBoundaries[2] << vtkstd::endl;
  os << indent << "NumberOfDOFNodes:              " << this->NumberOfDOFNodes << vtkstd::endl;
  os << indent << "NumberOfVerticesOnFace:        " << this->NumberOfVerticesOnFace << vtkstd::endl;
  os << indent << "EdgeArray:                     " << this->EdgeArray << vtkstd::endl;
  os << indent << "FaceArray:                     " << this->FaceArray << vtkstd::endl;
  os << indent << "FacesOnEdge:                   " << this->FacesOnEdge << vtkstd::endl;
  os << indent << "EdgesOnFace:                   " << this->EdgesOnFace << vtkstd::endl;
  os << indent << "ParametricCenter:              " << this->ParametricCenter << vtkstd::endl;
  os << indent << "ParametricCorners:             " << this->ParametricCorners << vtkstd::endl;
  os << indent << "BoundaryFaceTransforms:        " << this->BoundaryFaceTransforms << vtkstd::endl;
  os << indent << "BoundaryEdgeTransforms:        " << this->BoundaryEdgeTransforms << vtkstd::endl;
  os << indent << "ShapeFunctions:                " << this->ShapeFunctions << vtkstd::endl;
  os << indent << "ShapeFunctionDerivatives:      " << this->ShapeFunctionDerivatives << vtkstd::endl;
  os << indent << "SymbolicFieldGradient:         " << this->SymbolicFieldGradient << vtkstd::endl;
  os << indent << "PrepareForOrder:               " << this->PrepareForOrder << vtkstd::endl;
  os << indent << "EmbedEdgeCoords:               " << this->EmbedEdgeCoords << vtkstd::endl;
  os << indent << "EmbedEdgeCoordsInFace:         " << this->EmbedEdgeCoordsInFace << vtkstd::endl;
  os << indent << "EmbedFaceCoords:               " << this->EmbedFaceCoords << vtkstd::endl;
  os << indent << "ProjectCoordsToEdge:           " << this->ProjectCoordsToEdge << vtkstd::endl;
  os << indent << "ProjectCoordsToFace:           " << this->ProjectCoordsToFace << vtkstd::endl;
  os << indent << "ProjectFaceCoordsToEdge:       " << this->ProjectFaceCoordsToEdge << vtkstd::endl;
  os << indent << "TransformFaceCoordsToStorage:  " << this->TransformFaceCoordsToStorage << vtkstd::endl;
  os << indent << "IsParameterInDomain:           " << this->IsParameterInDomain << vtkstd::endl;
}

vtkShoeCellMetaData::vtkShoeCellMetaData()
{
}

vtkShoeCellMetaData::~vtkShoeCellMetaData()
{
}

int vtkShoeCellMetaData::InitializeShoeCells()
{
  int i;
  // ===================================
  // Hexahedron - Lagrange - Tensor
  vtkShoeCellMetaData* m = vtkShoeCellMetaData::New();
  // FIXME!!! Hack for paper. m->Type should be a unique
  // integer for each vtkShoeCellMetaData instance.
  m->Type = VTK_HIGHER_ORDER_HEXAHEDRON;
  m->Dimension = 3;
  m->NumberOfBoundaries[0] = 8;
  m->NumberOfBoundaries[1] = 12;
  m->NumberOfBoundaries[2] = 6;
  m->NumberOfDOFNodes = 19;
  m->NumberOfVerticesOnFace = ShoeHexahedron::VerticesOnFace;
  m->ShapeFunctions = ShoeHexahedron::LagrangeTensorShapeFunctions;
  m->ShapeFunctionDerivatives = ShoeHexahedron::LagrangeTensorShapeFunctionDerivatives;
  m->SymbolicFieldGradient = ShoeHexahedron::SymbolicLagrangeFieldGradient;
  m->IsParameterInDomain = ShoeHexahedron::IsParameterInDomain;

  m->EdgeArray = new int*[ m->NumberOfBoundaries[1] ];
  m->FacesOnEdge = new int*[ m->NumberOfBoundaries[1] ];
  m->BoundaryEdgeTransforms = new double*[ m->NumberOfBoundaries[1] ];
  for (i=0; i<m->NumberOfBoundaries[1]; ++i)
    {
    m->EdgeArray[i] = ShoeHexahedron::EdgeArray[i];
    m->FacesOnEdge[i] = ShoeHexahedron::FacesOnEdge[i];
    m->BoundaryEdgeTransforms[i] = ShoeHexahedron::BoundaryEdgeTransforms[i];
    }

  m->FaceArray = new int*[ m->NumberOfBoundaries[2] ];
  m->EdgesOnFace = new int*[ m->NumberOfBoundaries[2] ];
  m->BoundaryFaceTransforms = new double*[ m->NumberOfBoundaries[2] ];
  for (i=0; i<m->NumberOfBoundaries[2]; ++i)
    {
    m->FaceArray[i] = ShoeHexahedron::FaceArray[i];
    m->EdgesOnFace[i] = ShoeHexahedron::EdgesOnFace[i];
    m->BoundaryFaceTransforms[i] = ShoeHexahedron::BoundaryFaceTransforms[i];
    }

  m->ParametricCenter = ShoeHexahedron::ParametricCenter;
  m->ParametricCorners = ShoeHexahedron::ParametricCorners;
  m->PrepareForOrder = ShoeInterpolants::LagrangePrepareForOrder;
  m->EmbedEdgeCoords = ShoeHexahedron::EmbedEdgeCoords;
  m->EmbedEdgeCoordsInFace = ShoeHexahedron::EmbedEdgeCoordsInFace;
  m->EmbedFaceCoords = ShoeHexahedron::EmbedFaceCoords;
  m->ProjectCoordsToEdge = ShoeHexahedron::ProjectCoordsToEdge;
  m->ProjectCoordsToFace = ShoeHexahedron::ProjectCoordsToFace;
  m->ProjectFaceCoordsToEdge = ShoeHexahedron::ProjectFaceCoordsToEdge;
  m->TransformFaceCoordsToStorage = ShoeHexahedron::TransformFaceCoordsToStorage;
  vtkShoeCellMetaData::All.insert( m );
  vtkShoeCellMetaData::ByShape[ shoe::Hexahedron ] = m;

  // ===================================
  // Point
  m = vtkShoeCellMetaData::New();
  m->Type = VTK_VERTEX; // FIXME: See note for m->Type above.
  m->Dimension = 0;
  m->NumberOfBoundaries[0] = 0;
  m->NumberOfBoundaries[1] = 0;
  m->NumberOfBoundaries[2] = 0;
  m->NumberOfDOFNodes = 0;
  m->NumberOfVerticesOnFace = 0;
  m->ShapeFunctions = 0; // FIXME: should return 1.0 for (0., 0., 0.) and 0. for everything else
  m->ShapeFunctionDerivatives = 0; // FIXME: should return -Inf in every direction
  m->SymbolicFieldGradient = 0;
  m->IsParameterInDomain = 0;
  m->EdgeArray = 0;
  m->FaceArray = 0;
  m->FacesOnEdge = 0;
  m->EdgesOnFace = 0;
  m->BoundaryFaceTransforms = 0;
  m->BoundaryEdgeTransforms = 0;
  m->ParametricCenter = ShoeHexahedron::ParametricCenter;
  m->ParametricCorners = 0;
  m->PrepareForOrder = vtkShoeCellMetaData::EmptyPrepareForOrder;
  m->EmbedEdgeCoords = 0;
  m->EmbedEdgeCoordsInFace = 0;
  m->EmbedFaceCoords = 0;
  m->ProjectCoordsToEdge = 0;
  m->ProjectCoordsToFace = 0;
  m->ProjectFaceCoordsToEdge = 0;
  m->TransformFaceCoordsToStorage = 0;
  vtkShoeCellMetaData::All.insert( m );
  vtkShoeCellMetaData::ByShape[ shoe::NumberOfCellShapes ] = m;

  return 1;
}

void vtkShoeCellMetaData::EmptyPrepareForOrder( const vtkShoeOrderTuple& o )
{
  // do nothing
  //cerr << "Prepare for " << o.Order[0] << ", " << o.Order[1] << ", " << o.Order[2] << vtkstd::endl; 
}

double vtkShoeCellMetaData::GetParametricDistance( double x[3] ) const
{
  cerr << "Oops. Not implemented.\n";
  return 0.0;
}

int vtkShoeCellMetaData::FindClosestBoundary( double pcoords[3], vtkShoeCellIterator*& boundary ) const
{
  cerr << "Oops. Not implemented.\n";
  return -1;
}

void vtkShoeCellMetaData::Derivatives( double*, vtkDoubleArray*, double* ) const
{
  cerr << "Oops. Not implemented.\n";
}

void vtkShoeCellMetaData::Shutdown()
{
  for (vtkstd::set<vtkShoeCellMetaData*>::iterator it = vtkShoeCellMetaData::All.begin(); it != vtkShoeCellMetaData::All.end(); ++it )
    {
    if ( (*it)->EdgeArray )   delete [] (*it)->EdgeArray;
    if ( (*it)->FaceArray )   delete [] (*it)->FaceArray;
    if ( (*it)->EdgesOnFace ) delete [] (*it)->EdgesOnFace;
    if ( (*it)->FacesOnEdge ) delete [] (*it)->FacesOnEdge;
    if ( (*it)->BoundaryEdgeTransforms) delete [] (*it)->BoundaryEdgeTransforms;
    if ( (*it)->BoundaryFaceTransforms) delete [] (*it)->BoundaryFaceTransforms;
    (*it)->Delete();
    }
}
