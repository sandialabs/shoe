// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeBox.h"
#include "vtkShoeBoxP.h"

#include "vtkObjectFactory.h"
#include "vtkGenericCellIterator.h"
#include "vtkGenericAttribute.h"
#include "vtkGenericAttributeCollection.h"
#include "vtkSimpleCellTessellator.h"
#include "vtkGenericStreamingTessellator.h"

#include "vtkShoeEnums.h"
#include "vtkShoeAttribute.h"
#include "vtkShoeBoxPartition.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellIterator.h"
#include "vtkShoeCellMetaData.h"
#include "vtkShoeCellSpecies.h"
#include "vtkShoeCellRecord.h"
#include "vtkShoePointIterator.h"
#include "vtkDataRecords.h"

#include <vtkDataArray.h>

#include <vtkstd/map>
#include "freelist"

// Uncomment this to debug the creation of mesh partitioning.
#undef VTK_DEBUG_SBPARTITION

// =SPM========================================================================
vtkShoeBoxP::vtkShoeBoxP()
{
}

vtkShoeBoxP::~vtkShoeBoxP()
{
}

vtkIdType vtkShoeBoxP::InsertNextCell( vtkShoeCellSpecies* sp, vtkIdType connectivityEntry, uint32_t permutation )
{
  assert(sp);
  vtkstd::map< int, vtkIdType >::iterator it = this->Extant.find( sp->GetId() );
  if ( it != this->Extant.end() )
    {
    this->Extant.insert( vtkstd::pair<int,vtkIdType>( sp->GetId(), 1 ) );
    // Check that interpolant/prodspace/order is defined for all attributes?
    // No. See docs in header.
    }
  else
    {
    it->second++;
    }

  vtkIdType id = this->Cells.grab();
  vtkShoeCellRecord* rec = &this->Cells[id];

  rec->Species = sp;
  rec->Offset = connectivityEntry;
  rec->NodePermutations = permutation;

  return id;
}

// =EPM========================================================================
vtkCxxRevisionMacro(vtkShoeBox,"$Revision: 9648 $");
vtkStandardNewMacro(vtkShoeBox);

int vtkShoeBox::ElementsRegistered = 0;

vtkShoeBox::vtkShoeBox()
{
  if ( ! vtkShoeBox::ElementsRegistered )
    {
    // Dynamically load our finite element definitions
    vtkShoeBox::RegisterElements();
    }
  this->TemporaryMesh = 0; // not temporary by default
  this->LinkState = 1; // By default, update links as cells are inserted

  this->Geometry = vtkShoeAttribute::New();
  this->Geometry->Register( this );
  this->Geometry->FastDelete();
  this->Geometry->SetName( "Geometry" );
  this->Geometry->SetNumberOfComponents( 3 );
  this->Geometry->SetStorageType( VTK_DOUBLE );
  // This way, geometry is just a special case of an attribute
  // and can be used by filters that traditionally would only
  // handle fields:
  this->GeometryIndex = this->Attributes->GetNumberOfAttributes();
  this->Attributes->InsertNextAttribute( this->Geometry );

  this->Connectivity = vtkDataRecords::New();
  this->Connectivity->Register( this );
  this->Connectivity->FastDelete();
  this->Connectivity->SetStorageType( VTK_ID_TYPE );
  this->Connectivity->SetNumberOfComponents( 1 );

  this->CellRecs = new vtkShoeBoxP;

  this->PointLinks = vtkDataRecords::New();
  this->PointLinks->Register( this );
  this->PointLinks->FastDelete();
  this->PointLinks->SetStorageType( VTK_ID_TYPE );
  this->PointLinks->SetNumberOfComponents( 1 );

  this->DOFLinks = vtkDataRecords::New();
  this->DOFLinks->Register( this );
  this->DOFLinks->FastDelete();
  this->DOFLinks->SetStorageType( VTK_ID_TYPE );
  this->DOFLinks->SetNumberOfComponents( 1 );

  //this->Tessellator = vtkSimpleCellTessellator::New();
  if ( this->Tessellator )
    {
    this->Tessellator->Delete();
    }
  this->Tessellator = vtkGenericStreamingTessellator::New();
  this->Tessellator->Register(this);
  this->Tessellator->FastDelete();
  this->Tessellator->Initialize( this );
}

vtkShoeBox::~vtkShoeBox()
{
  if ( this->TemporaryMesh && this->TemporaryMesh != this )
    {
    this->TemporaryMesh->Delete();
    }
  this->Geometry->UnRegister(this);
  this->PointLinks->UnRegister(this);
  this->DOFLinks->UnRegister(this);
  this->Connectivity->UnRegister(this);
  delete this->CellRecs;
}

void vtkShoeBox::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Geometry: " << this->Geometry << vtkstd::endl;
  os << indent << "CellRecs: " << this->CellRecs << vtkstd::endl;
  os << indent << "Connectivity: " << this->Connectivity << vtkstd::endl;
  os << indent << "PointLinks: " << this->PointLinks << vtkstd::endl;
  os << indent << "DOFLinks: " << this->DOFLinks << vtkstd::endl;
}

void vtkShoeBox::RegisterElements()
{
  //FIXME
  //This is a hack for now. Eventually, this will call a dynamic loader
  //that looks for libraries of cells...
  vtkShoeCellMetaData::InitializeShoeCells();
  vtkShoeBox::ElementsRegistered = 1;
}

vtkIdType vtkShoeBox::GetNumberOfPoints()
{
  return this->Geometry->GetNumberOfPoints();
}

vtkIdType vtkShoeBox::GetNumberOfDOFNodes()
{
  return this->Geometry->GetNumberOfDOFNodes();
}

vtkIdType vtkShoeBox::GetNumberOfCells( int dim )
{
  assert( dim < 4 );
  if ( dim < 0 )
    {
    return this->CellRecs->Cells.size();
    }
  return this->NumberOfCellsOfDim[dim];
}

void vtkShoeBox::GetCellTypes( vtkCellTypes* types )
{
  (void)types;
  vtkWarningMacro( "No GetCellTypes" );
}

vtkGenericCellIterator* vtkShoeBox::NewCellIterator( int dim )
{
  vtkShoeCellIterator* it = vtkShoeCellIterator::New();
  it->SetDataSet( this );
  it->TraversalStyle = MeshOrder;
  it->Dimension = dim;
  return it;
}

vtkGenericCellIterator* vtkShoeBox::NewBoundaryIterator( int dim, int exteriorOnly )
{
  (void)dim;
  (void)exteriorOnly;
  vtkWarningMacro( "No bdy iterator yet" );
  return 0;
}

vtkGenericPointIterator* vtkShoeBox::NewPointIterator()
{
  vtkShoePointIterator* it = vtkShoePointIterator::New();
  it->SetDataSet( this );
  return it;
}

int vtkShoeBox::FindCell( double x[3], vtkGenericCellIterator*& cell, double tol2, int &subId, double pcoords[3] )
{
  (void)x;
  (void)cell;
  (void)tol2;
  (void)subId;
  (void)pcoords;
  vtkWarningMacro( "No FindCell" );
  return -1;
}

void vtkShoeBox::FindPoint( double x[3], vtkGenericPointIterator *p )
{
  (void)x;
  (void)p;
  vtkWarningMacro( "No FindPoint" );
}

void vtkShoeBox::ComputeBounds()
{
  int c;
  int nc = this->Geometry->GetNumberOfComponents();
  assert( nc <= 3 );
  for ( c = 0; c < nc; ++ c )
    this->Geometry->GetRange( c, this->Bounds + 2*c );
  for ( c = nc; c < 3; ++ c )
    this->Bounds[2*c] = this->Bounds[2*c + 1] = 0.;
}

unsigned long vtkShoeBox::GetActualMemorySize()
{
  return 1024*1024;
}

vtkIdType vtkShoeBox::GetEstimatedSize()
{
  return 1024*1024;
}

int vtkShoeBox::AttributeId( vtkGenericAttribute* a )
{
  for ( int i=0; i<this->Attributes->GetNumberOfAttributes(); ++ i )
    {
    if ( this->Attributes->GetAttribute(i) == a )
      {
      return i;
      }
    }
  return -1;
}

int vtkShoeBox::SetGeometry( vtkShoeAttribute* a )
{
  if ( this->Geometry == a )
    {
    return -1;
    }

  if ( this->Geometry )
    {
    this->Attributes->RemoveAttribute( this->GeometryIndex );
    this->Geometry->UnRegister( this );
    }

  for ( int i=0; i<this->Attributes->GetNumberOfAttributes(); ++ i )
    {
    if ( vtkShoeAttribute::SafeDownCast(this->Attributes->GetAttribute(i)) == a )
      {
      this->GeometryIndex = i;
      this->Geometry = a;
      this->Geometry->Register( this );
      return this->GeometryIndex;
      }
    }

  // FIXME!!! InsertNextAttribute should return a number but it doesn't!
  this->GeometryIndex = this->Attributes->GetNumberOfAttributes();
  this->Attributes->InsertNextAttribute( a );
  this->Geometry = a;
  this->Geometry->Register( this );
  assert( this->Attributes->GetAttribute( this->GeometryIndex ) == this->Geometry );

  return this->GeometryIndex;
}

vtkIdType vtkShoeBox::InsertNextCell( int species, vtkIdType* connectivity, uint32_t permutation )
{
  vtkShoeCellSpecies* sp = vtkShoeCellSpecies::GetSpeciesById( species );
  if ( !sp )
    {
    vtkWarningMacro( "Could not find species with ID " << species );
    return -1;
    }

  vtkIdType connEntry = this->Connectivity->InsertNextRecord( sp->GetNumberOfConnectivityEntries(), connectivity );
  vtkIdType id = this->CellRecs->InsertNextCell( sp, connEntry, permutation );
  if ( this->GetMaintainLinks() )
    {
    this->AddLinksForCell( id );
    }
  return id;
}

void vtkShoeBox::UpdateLinks()
{
  this->DOFLinks->Reset();
  this->PointLinks->Reset();

  this->PointLinks->SetNumberOfRecords( this->Geometry->GetNumberOfPoints() );
  this->DOFLinks->SetNumberOfRecords( this->Geometry->GetNumberOfDOFNodes() );

  vtkShoeCellIterator* it = vtkShoeCellIterator::SafeDownCast( this->NewCellIterator() );
  for ( it->Begin(); ! it->IsAtEnd(); it->Next() )
    {
    this->AddLinksForCell( it->GetCellId() );
    }
  it->Delete();
}

void vtkShoeBox::AddLinksForCell( vtkIdType cellId )
{
  vtkShoeCellRecord* cr = &this->CellRecs->Cells[ cellId ];
  const vtkShoeCellMetaData* meta = cr->Species->Meta;
  vtkIdType* conn = this->Connectivity->GetRecordAsIdType( cr->Offset );

  for ( int c=0; c<meta->NumberOfBoundaries[0]; ++ c )
    {
    this->AddLink( this->PointLinks, *conn, cellId );
    ++ conn;
    }
  for ( int d=0; d<meta->NumberOfDOFNodes; ++ d )
    {
    this->AddLink( this->DOFLinks, *conn, cellId );
    ++ conn;
    }
}

void vtkShoeBox::DeleteLinksForCell( vtkIdType cellId )
{
  vtkShoeCellRecord* cr = &this->CellRecs->Cells[ cellId ];
  const vtkShoeCellMetaData* meta = cr->Species->Meta;
  vtkIdType* conn = this->Connectivity->GetRecordAsIdType( cr->Offset );

  for ( int c=0; c<meta->NumberOfBoundaries[0]; ++ c )
    {
    this->DeleteLink( this->PointLinks, *conn, cellId );
    ++ conn;
    }
  for ( int d=0; d<meta->NumberOfDOFNodes; ++ d )
    {
    this->DeleteLink( this->DOFLinks, *conn, cellId );
    ++ conn;
    }
}

void vtkShoeBox::AddLink( vtkDataRecords* links, vtkIdType id, vtkIdType cell )
{
  links->InsertNextTuple( id, &cell );
}

void vtkShoeBox::DeleteLink( vtkDataRecords* links, vtkIdType id, vtkIdType cell )
{
  vtkIdType* R = links->GetRecordAsIdType( id );
  vtkIdType N = links->GetNumberOfTuples( id );
  for ( vtkIdType i=0; i<N; ++ i )
    {
    if ( R[i] == cell )
      {
      links->RemoveTuple( id, i );
      break;
      }
    }
}

vtkShoeCell* vtkShoeBox::GetCell( vtkIdType id )
{
  if ( this->CellRecs->Cells.is_unused( id ) )
    {
    return 0;
    }
  vtkShoeCell* sc = vtkShoeCell::New();
  sc->SetDataSet( this );
  sc->CellID = id;
  sc->CellRecord = this->CellRecs->Cells[ id ];
  sc->Meta = sc->CellRecord.Species->Meta;
  sc->InvalidateCache();
  return sc;
}

void vtkShoeBox::ComputeCriticalPoints( vtkShoeAttribute* sa )
{
  // Not a smart loop. Try this instead:
  // Loop over DOF D not yet computed
  //   for each mesh M referenced by sa
  //     if M has a cell C that references D
  //       C->ComputeDOFCriticalPoints( sa )
  //       break;
  // Mark sa->CriticalPointsDirty = 0 (clean)
  if ( ! sa->GetCriticalPointsDirty() )
    {
    vtkstd::cout << "\n\n===============\n\nno recompute\n\n====================\n\n";
    return;
    }

  vtkGenericCellIterator* gcit = this->NewCellIterator();
  for ( gcit->Begin(); ! gcit->IsAtEnd(); gcit->Next() ) 
    {
    vtkShoeCell* sc = vtkShoeCell::SafeDownCast( gcit->GetCell() );
    sc->ComputeDOFCriticalPoints( sa );
    }
  sa->SetCriticalPointsDirty( 0 );
}

void vtkShoeBox::ComputeCriticalPoints( vtkGenericAttributeCollection* kappa )
{
  int na = kappa->GetNumberOfAttributes();
  for ( int i = 0; i < na; ++ i )
    {
    vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( i ) );
    if ( ! sa ) continue;
    this->ComputeCriticalPoints( sa );
    }
}

struct vtkCritPointSorter {
  bool operator() ( const vtkstd::vector<double>& a, const vtkstd::vector<double>&b )
    {
    return a[0] < b[0];
    }
};

void vtkShoeBox::TriangulateBoundaries( vtkGenericAttributeCollection* kappa, vtkShoeBoxPartition* bdyTri )
{
  vtkstd::cout << "#### TriangulateBoundaries:\n";

  int na = kappa->GetNumberOfAttributes();
  vtkShoeCellIterator* it = vtkShoeCellIterator::SafeDownCast(this->NewCellIterator());
  for ( it->Begin(); ! it->IsAtEnd(); it->Next() )
    {
    vtkShoeCell* sc = vtkShoeCell::SafeDownCast(it->GetCell());
    if ( !sc )
      continue;

    vtkIdType* conn = sc->GetConnectivity() + sc->GetNumberOfPoints();
    vtkIdType ne = sc->GetNumberOfBoundaries( 1 );
    vtkIdType nf = sc->GetNumberOfBoundaries( 2 );
    vtkShoeCellRecord& cellRec( sc->GetRecord() );
    int dim = sc->GetDimension();

    if ( (dim < 2)  || 
         (bdyTri->GetNumberOfElements( conn[ne + nf] ) != 0) )
      {
      // If this is an edge node, skip it.
      // If we've already processed this volume or face node, skip it.
      continue;
      }

    double rep;    // projected edge coordinates (cell-aligned)
    double rel;    // face-loop edge coordinates (aligned w/ edge loop bounding face)
    double res;    // storage-order edge coordinates
    double rfe[2]; // face-embedded version of rep
    //double rce[3]; // cell-embedded version of rep (or rfe, for that matter)
    vtkstd::vector< vtkstd::vector<double> > critPts;

    if ( dim == 2 )
      {
      for ( int e = 0; e < ne; ++ e )
        {
        // The start of the edge is always "critical".
        vtkstd::vector<double> edgeEndPoint;
        edgeEndPoint.push_back( -1.0 ); // rel
        bool c0 = sc->GetEdgePermutation( e );
        bool c1 = sc->GetMetaData()->FaceCoedgeReversed( 0 /*face*/, e );
        edgeEndPoint.push_back( (c0 ^ c1) ? 1.0 : -1.0 ); // res
        edgeEndPoint.push_back( c0 ? 1.0 : -1.0 ); // rep
        edgeEndPoint.push_back( edgeEndPoint.back() ); // rfe[0]
        edgeEndPoint.push_back( 0.0 ); // rfe[1]
        sc->GetMetaData()->EmbedEdgeCoords( e, &edgeEndPoint[3] );
        critPts.push_back( edgeEndPoint );

        for ( int a = 0; a < na; ++ a )
          {
          vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( a ) );
          if ( ! sa ) continue;

          vtkDataRecords* cp;
          int nc = sa->GetNumberOfComponents();
          for ( int c = 0; c < nc; ++ c )
            {
            cp = sa->GetCriticalPoints()[c];
            int cptsz = cp->GetNumberOfComponents(); // critical point "tuple size"
            vtkIdType ncp = cp->GetNumberOfTuples( conn[e] );
            double* param = cp->GetRecordAsDouble( conn[e] );
            vtkstd::vector<double> point;
            for ( int critpt = 0; critpt < ncp; ++ critpt )
              {
              res = param[critpt*cptsz];
              rel = res;
              if ( c0 ^ c1 )
                vtkShoeCellRecord::ParametricEdgeCoordsFromStorageOrder( 1, &rel, true, 0 );
              rep = res;
              cellRec.ParametricEdgeCoordsFromStorageOrder( e, &rep );
              rfe[0] = rep;
              sc->GetMetaData()->EmbedEdgeCoords( e, rfe ); // for 3-D cells, this should be EmbedEdgeCoordsInFace
              point.push_back( rel );
              point.push_back( res );
              point.push_back( rep );
              point.push_back( rfe[0] );
              point.push_back( rfe[1] );
              critPts.push_back( point );
              point.clear();
              }
            }
          }
        // Sort the critical points along the edge.
        vtkstd::sort( critPts.begin(), critPts.end(), vtkCritPointSorter() );

        // Insert the ordered points, ignoring duplicates
        double last = -2.;
        for ( vtkstd::vector< vtkstd::vector< double > >::iterator it = critPts.begin(); it != critPts.end(); ++ it )
          {
          if ( (*it)[0] == last )
            continue;
          bdyTri->AddPointToFaceLoop( sc, conn, e, &(*it)[1] );
          last = (*it)[0];
          }
        critPts.clear();
        }

      // Now insert face-interior critical points, using the first to create a triangle fan.
      bool haveStarredFace = false;
      for ( int a = 0; a < na; ++ a )
        {
        vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( a ) );
        if ( ! sa ) continue;

        vtkDataRecords* cp;
        int nc = sa->GetNumberOfComponents();
        for ( int c = 0; c < nc; ++ c )
          {
          cp = sa->GetCriticalPoints()[c];
          int cptsz = cp->GetNumberOfComponents(); // critical point "tuple size"
          vtkIdType ncp = cp->GetNumberOfTuples( conn[ne] );
          double* param = cp->GetRecordAsDouble( conn[ne] );
          vtkstd::vector<double> point;
          double rs[5]; // (r,s) projection of (r,s,t) to face
          for ( int icp = 0; icp < ncp; ++ icp )
            {
            if ( sc->GetMetaData()->IsParameterInDomain( param, -1 ) == -1 )
              {
              continue;
              }
            rs[0] = param[0];
            rs[1] = param[1];
            rs[2] = rs[0];
            rs[3] = rs[1];
            cellRec.ParametricFaceCoordsFromStorageOrder( 0 /* face */, &rs[2] );

            if ( !haveStarredFace )
              {
              haveStarredFace = true;
              bdyTri->TriangulateFaceLoop( sc, conn, ne, ne, rs, 0 /* IsBoundary */ );
              }
            else
              {
              bdyTri->InsertPointIntoTriangulation( sc, conn, ne, rs );
              }
            param += cptsz;
            }
          }
        }
      if ( ! haveStarredFace )
        {
        double pcenter[3];
        sc->GetParametricCenter( pcenter );
        pcenter[2] = pcenter[0];
        pcenter[3] = pcenter[1];
        cellRec.ParametricFaceCoordsToStorageOrder( 0, &pcenter[0] );
        // We don't have any interior critical points.
        // For now, assume cell center is visible from all points on boundary.
        bdyTri->TriangulateFaceLoop( sc, conn, ne, ne, pcenter, 0  /* IsBoundary */);
        }

      // At this point, we are guaranteed to have a triangulation.
      // Now correct its topology to satisy kappa-compatibility.
      bdyTri->CorrectTriTopology( sc, conn, ne, ne, kappa );
      }
    else if ( dim == 3 )
      {
      // all of the above, but for each face of the cell
      for ( int face = 0; face < nf; ++ face )
        {
        if ( bdyTri->GetNumberOfElements( conn[ne + face] ) != 0 )
          {
          // If this face has been triangulated by a cell that shares it,
          // skip this face.
          continue;
          }
        for ( int e = 0; e < sc->GetMetaData()->NumberOfVerticesOnFace[face]; ++ e )
          {
          int edge = sc->GetMetaData()->EdgesOnFace[face][e];
          // The start of the edge is always "critical".
          vtkstd::vector<double> edgeEndPoint;
          edgeEndPoint.push_back( -1.0 ); // rel
          bool c0 = sc->GetEdgePermutation( edge );
          bool c1 = sc->GetMetaData()->FaceCoedgeReversed( face, e );
          edgeEndPoint.push_back( (c0 ^ c1) ? 1.0 : -1.0 ); // res
          edgeEndPoint.push_back( c1 ? 1.0 : -1.0 ); // rep
          edgeEndPoint.push_back( c1 ? 1.0 : -1.0 ); // rfe[0]
          edgeEndPoint.push_back( 0.0 ); // rfe[1]
          sc->GetMetaData()->EmbedEdgeCoordsInFace( edge, face, &edgeEndPoint[3] );
          critPts.push_back( edgeEndPoint );
#ifdef VTK_DEBUG_SBPARTITION
          vtkstd::cout << "Face " << face << " Edge " << edge << " =("
            << sc->GetMetaData()->EdgeArray[edge][0] << "," << sc->GetMetaData()->EdgeArray[edge][1] << "), rfl "
            << edgeEndPoint[0] << " rst " << edgeEndPoint[1] << " rep " << edgeEndPoint[2]
            << " rfp " << edgeEndPoint[3] << "," << edgeEndPoint[4] << vtkstd::endl;
#endif // VTK_DEBUG_SBPARTITION

          for ( int a = 0; a < na; ++ a )
            {
            vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( a ) );
            if ( ! sa ) continue;

            vtkDataRecords* cp;
            int nc = sa->GetNumberOfComponents();
            for ( int c = 0; c < nc; ++ c )
              {
              cp = sa->GetCriticalPoints()[c];
              int cptsz = cp->GetNumberOfComponents(); // critical point "tuple size"
              vtkIdType ncp = cp->GetNumberOfTuples( conn[edge] );
              double* param = cp->GetRecordAsDouble( conn[edge] );
              vtkstd::vector<double> point;
              for ( int critpt = 0; critpt < ncp; ++ critpt )
                {
                res = param[critpt*cptsz];
                rel = res;
                if ( c0 ^ c1 )
                  vtkShoeCellRecord::ParametricEdgeCoordsFromStorageOrder( 1, &rel, true, 0 );
                rep = res;
                cellRec.ParametricEdgeCoordsFromStorageOrder( e, &rep );
                rfe[0] = rep;
                sc->GetMetaData()->EmbedEdgeCoordsInFace( edge, face, rfe );
                point.push_back( rel );
                point.push_back( res );
                point.push_back( rep );
                point.push_back( rfe[0] );
                point.push_back( rfe[1] );
#ifdef VTK_DEBUG_SBPARTITION
                vtkstd::cout << "Face " << face << " Edge " << edge << " =("
                  << sc->GetMetaData()->EdgeArray[edge][0] << "," << sc->GetMetaData()->EdgeArray[edge][1] << "), rfl "
                  << point[0] << " rst " << point[1] << " rep " << point[2]
                  << " rfp " << point[3] << "," << point[4] << vtkstd::endl;
#endif // VTK_DEBUG_SBPARTITION
                critPts.push_back( point );
                point.clear();
                }
              }
            }
          // Sort the critical points along the edge.
          vtkstd::sort( critPts.begin(), critPts.end(), vtkCritPointSorter() );

          // Insert the ordered points, ignoring duplicates
          double last = -2.;
          for ( vtkstd::vector< vtkstd::vector< double > >::iterator it = critPts.begin(); it != critPts.end(); ++ it )
            {
            if ( (*it)[0] == last )
              continue;

            bdyTri->AddPointToFaceLoop( sc, conn, edge, &(*it)[1] );
            last = (*it)[0];
            }
          critPts.clear();
          }

        // We've been all the way around the face now... time to create a triangle fan
        // Now insert face-interior critical points, using the first to create a triangle fan.
        bool haveStarredFace = false;
        for ( int a = 0; a < na; ++ a )
          {
          vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( a ) );
          if ( ! sa ) continue;

          vtkDataRecords* cp;
          int nc = sa->GetNumberOfComponents();
          for ( int c = 0; c < nc; ++ c )
            {
            cp = sa->GetCriticalPoints()[c];
            int cptsz = cp->GetNumberOfComponents(); // critical point "tuple size"
            vtkIdType ncp = cp->GetNumberOfTuples( conn[ne + face] );
            double* param = cp->GetRecordAsDouble( conn[ne + face] );
            vtkstd::vector<double> point;
            double rs[5]; // (r,s) projection of (r,s,t) to face
            for ( int icp = 0; icp < ncp; ++ icp )
              {
              if ( sc->GetMetaData()->IsParameterInDomain( param, ne + face ) == -1 )
                {
                // skip critical points outside the cell
                continue;
                }
              rs[0] = param[0];
              rs[1] = param[1];
              rs[2] = rs[0];
              rs[3] = rs[1];
              cellRec.ParametricFaceCoordsFromStorageOrder( face, &rs[2] );

              if ( !haveStarredFace )
                {
                haveStarredFace = true;
                bdyTri->TriangulateFaceLoop( sc, conn, ne, ne + face, rs, 1 /* IsBoundary */ );
                }
              else
                {
                bdyTri->InsertPointIntoTriangulation( sc, conn, ne + face, rs );
                }
              param += cptsz;
              }
            }
          }
        if ( ! haveStarredFace )
          {
          double pcenter[4];
          for ( int i = 0; i < 3; ++ i )
            pcenter[i] = sc->GetMetaData()->BoundaryFaceTransforms[face][i];
          sc->GetMetaData()->ProjectCoordsToFace( face, pcenter );
          pcenter[2] = pcenter[0];
          pcenter[3] = pcenter[1];
          cellRec.ParametricFaceCoordsToStorageOrder( face, &pcenter[0] );
          // We don't have any interior critical points.
          // For now, assume cell center is visible from all points on boundary.
          bdyTri->TriangulateFaceLoop( sc, conn, ne, ne + face, pcenter, 1  /* IsBoundary */);
          }

        // At this point, we are guaranteed to have a triangulation.
        // Now correct its topology to satisy kappa-compatibility.
        bdyTri->CorrectTriTopology( sc, conn, ne, ne + face, kappa );
        }
      }
    }
}

void vtkShoeBox::TetrahedralizeInteriors( vtkGenericAttributeCollection* kappa, vtkShoeBoxPartition* bdyTri, vtkShoeBoxPartition* tetInt )
{
  // All of the points in the triangulation will be used in the tetrahedralization.
  //tetInt->SetPoints( bdyTri->GetPoints() );

  vtkstd::cout << "#### TetrahedralizeInteriors:\n";

  vtkShoeCellIterator* it = vtkShoeCellIterator::SafeDownCast(this->NewCellIterator());
  for ( it->Begin(); ! it->IsAtEnd(); it->Next() )
    {
    vtkShoeCell* sc = vtkShoeCell::SafeDownCast(it->GetCell());
    if ( !sc || sc->GetDimension() < 3 )
      continue;

    vtkIdType* conn = sc->GetConnectivity() + sc->GetNumberOfPoints();
    vtkIdType ne = sc->GetNumberOfBoundaries( 1 );
    vtkIdType nf = sc->GetNumberOfBoundaries( 2 );
    vtkIdType vdof = conn[ne + nf];

    // Grab all of the interior critical points
    int na = kappa->GetNumberOfAttributes();
    bool haveStarredInterior = false;
    for ( int a = 0; a < na; ++ a )
      {
      vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( a ) );
      if ( ! sa ) continue;

      vtkDataRecords* cp;
      int nc = sa->GetNumberOfComponents();
      for ( int c = 0; c < nc; ++ c )
        {
        cp = sa->GetCriticalPoints()[c];
        int cptsz = cp->GetNumberOfComponents(); // critical point "tuple size"
        vtkIdType ncp = cp->GetNumberOfTuples( vdof );
        double* param = cp->GetRecordAsDouble( vdof );
        int icp = 0;

        if ( ( ! haveStarredInterior ) && ncp )
          {
          vtkstd::cout << "\n### Critical point " << icp << ": ( "
                       << param[0] << ", "
                       << param[1] << ", "
                       << param[2] << " )\n";
        
          haveStarredInterior = true;
          tetInt->StarInterior( bdyTri, sc, conn, ne, nf, param );
          param += cptsz;
          ++ icp;
          }

        for ( ; icp < ncp; ++ icp )
          {
          vtkstd::cout << "\n### Critical point " << icp << ": ( "
                       << param[0] << ", "
                       << param[1] << ", "
                       << param[2] << " )\n";

          tetInt->InsertPointIntoTetrahedralization( sc, conn, ne + nf, param );
          param += cptsz;
          }
        }
      }
    if ( ! haveStarredInterior )
      {
      // We don't have any interior critical points.
      // For now, assume cell center is visible from all points on boundary.
      double pcenter[3];
      sc->GetParametricCenter( pcenter );
      tetInt->StarInterior( bdyTri, sc, conn, ne, nf, pcenter );
      }

    // At this point, we are guaranteed to have a tetrahedralization.
    // Now correct its topology to satisy kappa-compatibility.
    tetInt->CorrectTetTopology( sc, conn, ne, nf, kappa );
    }
}

void vtkShoeBox::CorrectPartitionTopology( vtkGenericAttributeCollection* kappa, vtkShoeBoxPartition* tetInt )
{
  int na = kappa->GetNumberOfAttributes();
  for ( int i = 0; i < na; ++ i )
    {
    vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( i ) );
    if ( ! sa ) continue;
    //this->CorrectPartitionTopology( sa, tetInt );
    (void)tetInt; // FIXME
    (void)sa;
    }
}

vtkIdType vtkShoeBox::GetFirstCellId()
{
  // FIXME!
  return 0;
}

vtkIdType vtkShoeBox::GetFirstBoundaryDOF()
{
  // FIXME!
  return 0;
}

vtkShoeBoxPartition* vtkShoeBox::PartitionMesh( vtkGenericAttributeCollection* kappa )
{
  this->ComputeCriticalPoints( kappa );
  vtkShoeBoxPartition* btri = vtkShoeBoxPartition::New();
  vtkstd::cout << "\n\n========\n\n";
  this->TriangulateBoundaries( kappa, btri );
  btri->Dump();
  vtkstd::cout << "\n\n========\n\n";
  vtkShoeBoxPartition* btet = vtkShoeBoxPartition::New();
  this->TetrahedralizeInteriors( kappa, btri, btet );
  btri->Delete();
  // This is now done in TetrahedralizeInteriors()
  //this->CorrectPartitionTopology( kappa, btet );
  return btet;
}

void vtkShoeBox::Reset()
{
  this->Attributes->Reset();
  this->Geometry->Delete();
  this->Geometry = vtkShoeAttribute::New();
  this->Geometry->Register( this );
  this->Geometry->FastDelete();
  this->Geometry->SetName( "Geometry" );
  this->Geometry->SetNumberOfComponents( 3 );
  this->Geometry->SetStorageType( VTK_DOUBLE );
  // This way, geometry is just a special case of an attribute
  // and can be used by filters that traditionally would only
  // handle fields:
  this->GeometryIndex = this->Attributes->GetNumberOfAttributes();
  this->Attributes->InsertNextAttribute( this->Geometry );
  this->Connectivity->Reset();
  this->CellRecs->Cells.clear();
  for ( int i=0; i<4; ++ i )
    {
    this->NumberOfCellsOfDim[i] = 0;
    }
  this->PointLinks->Reset();
  this->DOFLinks->Reset();
  if ( this->TemporaryMesh )
    {
    this->TemporaryMesh->Delete();
    }
  this->TemporaryMesh = 0;
  this->LinkState = 1; // By default, update links as cells are inserted
}

void vtkShoeBox::CreateTemporaryMesh()
{
  if ( ! this->TemporaryMesh )
    {
    this->TemporaryMesh = vtkShoeBox::New();
    this->TemporaryMesh->TemporaryMesh = this->TemporaryMesh; // mark the mesh as temporary
    this->TemporaryMesh->Attributes->ShallowCopy( this->Attributes );
    this->TemporaryMesh->SetGeometry( this->Geometry );
    }
}
