// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkMath.h"
#include "vtkstd/list"

#include "vtkGenericAttribute.h"
#include "vtkGenericAttributeCollection.h"

#include "vtkShoeAttribute.h"
#include "vtkPolynomialSystem.h"
#include "vtkShoeBoxPartition.h"
#include "vtkShoeBoxPartitionIterator.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellMetaData.h"

#include <float.h>

#include "vtkObjectFactory.h"

// Choose a format for dumps:
// 'T' Text
/* #define VTK_SBP_DUMP_FORMAT 'T' */
// 'P' PLY (Stanford polygon)
/* #define VTK_SBP_DUMP_FORMAT 'P' */
// 'U' Unstructured grid (.vtu)
// #define VTK_SBP_DUMP_FORMAT 'U'
#define VTK_SBP_DUMP_FORMAT 'U'

#if VTK_SBP_DUMP_FORMAT == 'U'
#  include "vtkCellType.h"
#  include "vtkUnstructuredGrid.h"
#  include "vtkUnstructuredGridWriter.h"
#  include "vtkPoints.h"
#endif // VTK_SBP_DUMP_FORMAT == 'U'

// Define/undefine various debugging flags as needed
#define VTK_VERBOSE_SBP_FLOW
#define VTK_VERBOSE_CORRECT_EDGE
#define VTK_VERBOSE_CORRECT_TRI
#define VTK_VERBOSE_CORRECT_TET
#undef VTK_VERBOSE_MERGE
#undef VTK_DEBUG_INS_PT_TET
#undef VTK_DEBUG_BND_EDGE_CP
#define VTK_DEBUG_INT_EDGE_CP
#undef VTK_DEBUG_SOLVE_EDGE_CP
#undef VTK_DEBUG_FACE_CP
#undef VTK_DEBUG_SPLIT_INT_EDGE
#undef VTK_DEBUG_SPLIT_BND_EDGE
#undef VTK_DEBUG_INSERT_POINT_IN_TRI
#undef VTK_DEBUG_INSERT_POINT_IN_TET

//double vtkShoeBoxPartition::Epsilon = 2.2204460492503131e-016;
const double vtkShoeBoxPartition::Epsilon = DBL_EPSILON;

int vtkShoeBoxPartition::VO[] = { 0, 1, 1, 2, 2, 0 };
int vtkShoeBoxPartition::VE[] = { 2, 5, 4, 1, 0, 3 };

vtkCxxRevisionMacro(vtkShoeBoxPartition,"$Revision: 10051 $");
vtkStandardNewMacro(vtkShoeBoxPartition);

vtkShoeBoxPartition::vtkShoeBoxPartition()
{
  this->OwnPoints = 1;
  this->Points = new vtkShoeBoxPartition::PointSet;
  this->NextPointLabel = 0;
  this->MaximumSplitEdgeDepth = 5;
  this->MaximumSplitFaceDepth = 5;
  this->NullTriangleThreshold = .01;
  this->NullTriangleThreshold2 = this->NullTriangleThreshold * this->NullTriangleThreshold;
  this->NullTetThreshold = .01;
}

vtkShoeBoxPartition::~vtkShoeBoxPartition()
{
  if ( this->OwnPoints )
    {
    delete this->Points;
    }
}

void vtkShoeBoxPartition::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Points:     " << this->Points << vtkstd::endl;
  os << indent << "OwnPoints:  " << this->OwnPoints << vtkstd::endl;
  os << indent << "EdgePoints: " << &this->EdgePoints << vtkstd::endl;
  os << indent << "Triangles:  " << &this->Triangles << vtkstd::endl;
}

void vtkShoeBoxPartition::Dump()
{
#if VTK_SBP_DUMP_FORMAT == 'T'
  vtkIndent indent;
  this->PrintSelf( vtkstd::cout, indent );
  vtkstd::cout << "===== POINTS =====\n";
  for ( PointSet::iterator psit = this->Points->begin(); psit != this->Points->end(); ++ psit )
    {
    for ( vtkstd::set<vtkShoeBoxPartition::Point>::iterator dpit = psit->second.begin(); dpit != psit->second.end(); ++ dpit )
      {
      vtkstd::cout.width(5);
      vtkstd::cout << dpit->Label << ": "
        << "DOF " << dpit->DOF << ", Element " << dpit->Element
        << " (" << dpit->R[0] << ", " << dpit->R[1] << ", " << dpit->R[2] << ", " << dpit->R[3] << ") "
        << " (" << dpit->X[0] << ", " << dpit->X[1] << ", " << dpit->X[2] << ") "
        << vtkstd::endl;
      }
    }

  if ( ! this->Offsets.empty() )
    {
    vtkstd::cout << "===== DOF OFFSETS =====\n";
    for ( vtkstd::map<vtkIdType,DofLocator>::iterator it = this->Offsets.begin(); it != this->Offsets.end(); ++ it )
      {
      vtkstd::cout << it->first << ": " << it->second.Start << ", " << it->second.NumberOfElements << vtkstd::endl;
      }
    }

  if ( this->Triangles.size() )
    {
    vtkstd::cout << "===== TRIANGLES =====\n";
    vtkIdType tnum = 0;
    for ( vtkstd::vector<Triangle>::iterator trit = this->Triangles.begin(); trit != this->Triangles.end(); ++trit, ++tnum )
      {
      vtkstd::cout.width(4);
      vtkstd::cout << tnum << ": (";
      vtkstd::cout.width(3);
      vtkstd::cout
        << trit->Flags
        << ") " << trit->EndPoints[0]->Label << " - " << trit->EndPoints[1]->Label << " - " << trit->EndPoints[2]->Label
        << ", Neighbors "
        << trit->Neighbors[0] << "=(" << (trit->Neighbors[0]/8) << "," << (trit->Neighbors[0]%8) << "), "
        << trit->Neighbors[1] << "=(" << (trit->Neighbors[1]/8) << "," << (trit->Neighbors[1]%8) << "), "
        << trit->Neighbors[2] << "=(" << (trit->Neighbors[2]/8) << "," << (trit->Neighbors[2]%8) << "), "
        << trit->Neighbors[3] << "=(" << (trit->Neighbors[3]/8) << "," << (trit->Neighbors[3]%8) << "), "
        << trit->Neighbors[4] << "=(" << (trit->Neighbors[4]/8) << "," << (trit->Neighbors[4]%8) << "), "
        << trit->Neighbors[5] << "=(" << (trit->Neighbors[5]/8) << "," << (trit->Neighbors[5]%8) << ")"
        << vtkstd::endl;
      }
    }

#elif VTK_SBP_DUMP_FORMAT == 'P'

  static int plyserial = 8;
  char plyname[256];
  sprintf( plyname, "sbp%03d.ply", plyserial++ );
  vtkstd::ofstream outy;
  outy.open( plyname, ofstream::out );
  vtkstd::map<vtkIdType,PointRef> pointmap;
  for ( PointSet::iterator psit = this->Points->begin(); psit != this->Points->end(); ++ psit )
    {
    for ( vtkstd::set<vtkShoeBoxPartition::Point>::iterator dpit = psit->second.begin(); dpit != psit->second.end(); ++ dpit )
      {
      pointmap[ dpit->Label ] = dpit;
      }
    }
  outy << "ply\nformat ascii 1.0\ncomment vtkShoeBoxPartition\nelement vertex " << pointmap.size()
    << "\nproperty float x\nproperty float y\nproperty float z\nelement face " << this->Triangles.size()
    << "\nproperty list uchar int vertex_indices\nend_header\n";
  for ( vtkstd::map<vtkIdType,PointRef>::iterator it = pointmap.begin(); it != pointmap.end(); ++ it )
    {
    outy << it->second->X[0] << " " << it->second->X[1] << " " << it->second->X[2] << "\n";
    }
  for ( vtkstd::vector<Triangle>::iterator trit = this->Triangles.begin(); trit != this->Triangles.end(); ++trit )
    {
    outy
      << "3 " << trit->EndPoints[0]->Label << " " << trit->EndPoints[1]->Label << " " << trit->EndPoints[2]->Label
      << vtkstd::endl;
    }
  outy.close();

#elif VTK_SBP_DUMP_FORMAT == 'U'

  for ( vtkstd::map<vtkIdType,DofLocator>::iterator dlit = this->Offsets.begin(); dlit != this->Offsets.end(); ++ dlit )
    {
    vtkstd::cout << "DOF " << dlit->first << ": " << dlit->second.Start << ", " << dlit->second.NumberOfElements << vtkstd::endl;
    }
  vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
  vtkPoints* pts = vtkPoints::New();
  vtkShoeBoxPartitionIterator* sbpit = vtkShoeBoxPartitionIterator::New();
  for ( vtkShoeBoxPartition::PointSet::iterator psit = this->Points->begin(); psit != this->Points->end(); ++ psit )
    {
    vtkstd::cout << "Dumping DOF " << psit->first << vtkstd::endl;
    for ( sbpit->Begin( this, psit->first ); ! sbpit->IsAtEnd(); sbpit->Next() )
      {
      double x[3], y[3], z[3], w[3];
      sbpit->GetWorldTetrahedron( x, y, z, w );
      vtkIdType conn[4];
      conn[0] = pts->InsertNextPoint( x );
      conn[1] = pts->InsertNextPoint( y );
      conn[3] = pts->InsertNextPoint( z );
      conn[2] = pts->InsertNextPoint( w );
      ug->InsertNextCell( VTK_TETRA, 4, conn );
      }
    }
  ug->SetPoints( pts );
  pts->Delete();
  vtkUnstructuredGridWriter* ugw = vtkUnstructuredGridWriter::New();
  ugw->SetInput( ug );
  ug->Delete();
  static int serialnum = 8;
  char dumpname[256];
  sprintf( dumpname, "sbp%03d.vtk", serialnum++ );
  ugw->SetFileName( dumpname );
  ugw->Write();

  ugw->Delete();
  sbpit->Delete();

#else
#  error "No format for Dump() selected. Set VTK_SBP_DUMP_FORMAT to T, P, or U"
#endif // 0
}

#define VTK_TRI_DEMUX_FV(f,v,fv) \
  f = fv >> 3;                   \
  v = fv & 7

#define VTK_TRI_DDEMUX_FV(f,v,fv) \
  int f = fv >> 3;                \
  int v = fv & 7

#define VTK_TRI_MUX_FV(f,v) ( (f)<<3 | (v) )

#define VTK_TRI_FACE_FROM_FV(fv) this->Triangles[(fv) >> 3]

#define VTK_TRI_F_FROM_FV(fv) ( (fv) >> 3 )

#define VTK_TRI_SWAP_FACEPTR( fa, va, fb, vb )                              \
    tmp = this->Triangles[fa].Neighbors[va] ;                               \
    this->Triangles[fa].Neighbors[va] = this->Triangles[fb].Neighbors[vb] ; \
    this->Triangles[fb].Neighbors[vb] = tmp ;


int vtkShoeBoxPartition::UnSanity( int dof )
{
  vtkIdType tBegin = this->GetFirstTriangle( dof );
  vtkIdType tEnd = tBegin + this->GetNumberOfElements( dof );
  vtkIdType tCur = 0;
  int flag = 0;
  cout << "Sanity check: DOF " << dof << " contains triangles " << tBegin << " to " << tEnd << "\n";
  while ( tCur + tBegin < tEnd )
    {
    for ( int n = 0; n < 6; ++n )
      {
      EleRef fv = this->Triangles[tBegin + tCur].Neighbors[n];
      VTK_TRI_DDEMUX_FV(nfv,vfv,fv);
      (void)vfv;
      if ( nfv < tBegin || nfv >= tEnd )
        {
        cout << "    " << dof << " refers to " << nfv << " as fnext(" << (tBegin+tCur) << "," << n << ")\n";
        ++flag;
        }
      }
    ++tCur;
    }
  return flag;
}

int vtkShoeBoxPartition::MadITellYou()
{
  int problemCount = 0;
  vtkstd::map<vtkIdType, DofLocator>::iterator it;
  for ( it = this->Offsets.begin(); it != this->Offsets.end(); ++ it )
    {
    problemCount += this->UnSanity( it->first );
    }
  return problemCount;
}

bool vtkShoeBoxPartition::CheckTetOrientations( vtkIdType d )
{
  vtkstd::cout << "* Checking tetrahedral orientations: ";

  bool orientation = true;

  vtkShoeBoxPartitionIterator* sbpit = vtkShoeBoxPartitionIterator::New();
  for ( sbpit->Begin( this, d ); ! sbpit->IsAtEnd(); sbpit->Next() )
    {
    EleRef fv0 = sbpit->GetCurrentEdgeFacet();
    EleRef fv1E = this->TriENext( fv0 );
    EleRef fv1F = this->TriFNext( fv0 );

    PointRef pr0 = this->TriOrg( fv0 );
    PointRef pr1 = this->TriOrg( fv1E );
    PointRef pr2 = this->TriOrg( this->TriENext( fv1E ) );
    PointRef pr3 = this->TriDest( this->TriENext( fv1F ) );

    double AB[3], AC[3], AD[3];
    for ( int i = 0; i < 3 ; ++ i )
      {
      AB[i] = pr1->R[i] - pr0->R[i];
      AC[i] = pr2->R[i] - pr0->R[i];
      AD[i] = pr3->R[i] - pr0->R[i];
      }
    double ABCD = vtkMath::Determinant3x3( AB, AC, AD );
    if ( ABCD < 0 ) 
      {
      vtkstd::cout << "\n Inverted tetrahedron in the tessellation ( "
                   << fv0 << " ): " 
                   << pr0->Label << " "
                   << pr1->Label << " "
                   << pr2->Label << " "
                   << pr3->Label << " "
                   << "det= " << ABCD;
      orientation = false;
      }
    }

  if ( orientation ) vtkstd::cout << " ok.";
  vtkstd::cout << "\n";

  return orientation;
}

int vtkShoeBoxPartition::IsTriangleOnHull( vtkIdType t, int frontOrBack )
{
  return this->TriIsFaceOnHull( VTK_TRI_MUX_FV(t,frontOrBack) );
}

void vtkShoeBoxPartition::GetTetrahedron( vtkIdType t, int frontOrBack, double* coords, vtkIdType* tris, int* sides )
{
  int verts[4];
  EleRef fv = VTK_TRI_MUX_FV(t,frontOrBack);
  EleRef nf;
  PointRef p;
  for ( int i = 0; i < 3; ++ i )
    {
    // Enumerate other co-facets that form this tetrahedron
    nf = this->TriFNext(fv);
    VTK_TRI_DDEMUX_FV(fnf,vnf,nf);
    tris[i] = fnf;
    sides[i] = vnf % 2 ? 0 : 1; // We want TriSym(nf), since that's inward-pointing!

    // Add the coordinates of the 3 points on the input triangle
    p = this->TriOrg( fv );
    verts[i] = p->Label;
    for ( int c = 0; c < 3; ++ c )
      {
      coords[c] = p->R[c];
#if VTK_SBP_DUMP_FORMAT == 'U'
      coords[c + 12] = p->X[c];
#endif // VTK_SBP_DUMP_FORMAT == 'U'
      }
    coords += 3;

    // Step around the input triangle to the next edge
    fv = this->TriENext( fv );
    }

  // Fill in the remaining coordinate using the last neighbor face we visited
  fv = this->TriENext2( nf );
  p = this->TriOrg( fv );
  verts[3] = p->Label;
  for ( int c = 0; c < 3; ++ c )
    {
    coords[c] = p->R[c];
#if VTK_SBP_DUMP_FORMAT == 'U'
    coords[c+12] = p->X[c];
#endif // VTK_SBP_DUMP_FORMAT == 'U'
    }
#if 0
  vtkstd::cerr << "GetTetrahedron " << (VTK_TRI_MUX_FV(t,frontOrBack)) << "=(" << t << ", " << frontOrBack << ") is "
    << verts[0] << " - " << verts[1] << " - " << verts[2] << " - " << verts[3] << vtkstd::endl;
#endif // 0
}

vtkShoeBoxPartition::PointRef vtkShoeBoxPartition::TriOrg( EleRef fv ) const
{
  VTK_TRI_DDEMUX_FV(f,v,fv) ;
  return this->Triangles[f].EndPoints[vtkShoeBoxPartition::VO[v]];
}

vtkShoeBoxPartition::PointRef vtkShoeBoxPartition::TriDest( EleRef fv ) const
{
  int tmp = this->TriENext( fv ) ;
  VTK_TRI_DDEMUX_FV(f,v,tmp) ;
  return this->Triangles[f].EndPoints[vtkShoeBoxPartition::VO[v]];
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartition::TriSym( EleRef fv ) const
{
  return fv + ( fv % 2 ? -1 : 1 );
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartition::TriENext( EleRef fv ) const
{
  VTK_TRI_DDEMUX_FV(f,v,fv);
  return VTK_TRI_MUX_FV(f,vtkShoeBoxPartition::VE[v]);
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartition::TriENext2( EleRef fv ) const
{
  VTK_TRI_DDEMUX_FV(f,v,fv) ;
  return VTK_TRI_MUX_FV(f,vtkShoeBoxPartition::VE[vtkShoeBoxPartition::VE[v]]);
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartition::TriFNext( EleRef fv ) const
{
  VTK_TRI_DDEMUX_FV(f,v,fv) ;
  return this->Triangles[f].Neighbors[v];
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartition::TriTurn( EleRef fv ) const
{
  VTK_TRI_DDEMUX_FV(f,v,fv);
  v = vtkShoeBoxPartition::VE[v];
  fv = VTK_TRI_MUX_FV(f,v); // fv <-- enext(fv)
  fv = this->TriSym( fv );// fv <-- sym(fv)
  VTK_TRI_DEMUX_FV(f,v,fv);
  fv = this->Triangles[f].Neighbors[v]; // fv <-- fnext(fv)
  VTK_TRI_DEMUX_FV(f,v,fv);
  v = vtkShoeBoxPartition::VE[v]; // fv <-- enext(fv)
  return VTK_TRI_MUX_FV(f,v);
}

void vtkShoeBoxPartition::TriFSplice( EleRef fva, EleRef fvb )
{
  int tmp;
  int fa, va, fb, vb;
  VTK_TRI_DEMUX_FV(fa,va,fva);
  VTK_TRI_DEMUX_FV(fb,vb,fvb);
  VTK_TRI_SWAP_FACEPTR(fa,va,fb,vb);

  fva = this->TriSym( this->TriFNext( fva ) );
  VTK_TRI_DEMUX_FV(fa,va,fva);
  fvb = this->TriSym( this->TriFNext( fvb ) );
  VTK_TRI_DEMUX_FV(fb,vb,fvb);
  VTK_TRI_SWAP_FACEPTR(fa,va,fb,vb);
}

void vtkShoeBoxPartition::TriInvolutiveFMerge( EleRef fva, EleRef fvb )
{
  this->TriFSplice( fva, this->TriSym( this->TriFNext( this->TriSym( fvb ) ) ) );
}

void vtkShoeBoxPartition::TriIdempotentFMerge( EleRef fva, EleRef fvb )
{
#ifdef VTK_VERBOSE_MERGE
  vtkstd::cout << "Merging "
               << this->TriOrg( fva )->Label << "-"
               << this->TriDest( fva )->Label << "-"
               << this->TriDest( this->TriENext( fva ) )->Label << " to "
               << this->TriOrg( fvb )->Label << "-"
               << this->TriDest( fvb )->Label << "-"
               << this->TriDest( this->TriENext( fvb ) )->Label << "\n";
#endif // VTK_VERBOSE_MERGE

  if ( this->TriFNext( fva ) != fvb ) this->TriInvolutiveFMerge( fva, fvb );
#ifdef VTK_VERBOSE_MERGE
  else vtkstd::cout << "... already glued, doing nothing!\n";
#endif // VTK_VERBOSE_MERGE
}
        
void vtkShoeBoxPartition::TriFDel( EleRef fva )
{
  // int fa = VTK_TRI_F_FROM_FV(fva);
#if 0
  VTK_TRI_DDEMUX_FV(f,v,fva);
  int vi = org(a);
  if ( VTK_TRI_F_FROM_FV( m_vt[vi]->edfa) == f )
    updateVertEdfaAll( fnext( m_vt[vi]->edfa ) );
#endif // 0
  int fvb = this->TriENext( fva );
  int fvc = this->TriENext( fvb );
  this->TriInvolutiveFMerge( fva, fva );
  this->TriInvolutiveFMerge( fvb, fvb );
  this->TriInvolutiveFMerge( fvc, fvc );

#if 0
  switch ( this->Triangles[f]->getFaceColor() ) 
    {
    case vtkShoeBoxPartition::GREEN:
      m_flippable--;
      break;
    case vtkShoeBoxPartition::RED:
      m_unflippable--;
      break;
    default:
      break;
    }
  this->Triangles[f]->markDeleted();
#endif
}

vtkIdType vtkShoeBoxPartition::TriCreate( PointRef& a, PointRef& b, PointRef& c )
{
  vtkIdType tnum = this->Triangles.size();
  vtkShoeBoxPartition::Triangle tri;
  tri.EndPoints[0] = a;
  tri.EndPoints[1] = b;
  tri.EndPoints[2] = c;
  tri.Flags = 48;
#if 0
  vtkstd::cout << "Triangle:";
  vtkstd::cout << " (" << a->R[2] << ", " << a->R[3] << ") -"
               << " (" << b->R[2] << ", " << b->R[3] << ") -"
               << " (" << c->R[2] << ", " << c->R[3] << ")\n";
#endif // 0
  for ( int i = 0; i < 6; ++ i )
    tri.Neighbors[i] = VTK_TRI_MUX_FV( tnum, i );
  this->Triangles.push_back( tri );
  return tnum;
}

vtkShoeBoxPartition::EleRef vtkShoeBoxPartition::TriFindOrCreate( DofLocator& dl, 
                                                                  PointRef a, PointRef b, PointRef c )
{
  int i = dl.Start;
  for ( vtkstd::vector<Triangle>::iterator it = this->Triangles.begin() + i; it != this->Triangles.end(); ++ it, ++ i )
    {
    if ( it->EndPoints[0] == a )
      {
      if ( it->EndPoints[1] == b && it->EndPoints[2] == c )
        { // we found a-b-c
        return VTK_TRI_MUX_FV(i,0);
        }
      else if ( it->EndPoints[1] == c && it->EndPoints[2] == b )
        { // found a-c-b
        return VTK_TRI_MUX_FV(i,5);
        }
      }
    else if ( it->EndPoints[0] == b )
      {
      if ( it->EndPoints[1] == a && it->EndPoints[2] == c )
        { // found b-a-c
        return VTK_TRI_MUX_FV(i,1);
        }
      else if ( it->EndPoints[1] == c && it->EndPoints[2] == a )
        { // found b-c-a
        return VTK_TRI_MUX_FV(i,4);
        }
      }
    else if ( it->EndPoints[0] == c )
      {
      if ( it->EndPoints[1] == a && it->EndPoints[2] == b )
        { // found c-a-b
        return VTK_TRI_MUX_FV(i,2);
        }
      else if ( it->EndPoints[1] == b && it->EndPoints[2] == a )
        { // found c-b-a
        return VTK_TRI_MUX_FV(i,3);
        }
      }
    }
  // No triangles match... create a new one
  vtkIdType tnum = this->TriCreate( a, b, c );
  ++ dl.NumberOfElements;
  return VTK_TRI_MUX_FV(tnum,0);
}

void vtkShoeBoxPartition::TriMarkEdge( EleRef fv, bool mark )
{
  // use the origin of the even co-edge as the bit associated
  // with the edge:
  if ( fv % 2 )
    { 
    fv = this->TriSym( fv );
    }
  VTK_TRI_DDEMUX_FV(f,v,fv);
  if ( mark )
    {
    this->Triangles[f].Flags |= (1<<(v/2));
    }
  else
    {
    this->Triangles[f].Flags &= ~(1<<(v/2));
    }
}

void vtkShoeBoxPartition::TriMarkEdges( EleRef fv, bool mark )
{
  EleRef start = fv;
  int i = 0;
  do {
    this->TriMarkEdge( fv, mark );
    fv = this->TriFNext( fv );
    ++ i;
  } while ( fv != start && i < 50 );
  if ( i == 50 )
    {
    vtkErrorMacro("More than 50 faces on edge " << start);
    }
}

void vtkShoeBoxPartition::TriMarkFace( EleRef fv, bool mark )
{
  vtkIdType f = VTK_TRI_F_FROM_FV(fv);
  if ( mark )
    {
    this->Triangles[f].Flags |= 0x08;
    }
  else
    {
    this->Triangles[f].Flags &= ~(0x08);
    }
}

bool vtkShoeBoxPartition::TriGetEdgeMark( EleRef fv )
{
  // use the origin of the even co-edge as the bit associated
  // with the edge:
  if ( fv % 2 )
    { 
    fv = this->TriSym( fv );
    }
  VTK_TRI_DDEMUX_FV(f,v,fv);
  return (this->Triangles[f].Flags & (1<<(v/2))) != 0;
}

bool vtkShoeBoxPartition::TriGetFaceMark( EleRef fv )
{
  return (this->Triangles[VTK_TRI_F_FROM_FV(fv)].Flags & 0x08) != 0;
}

void vtkShoeBoxPartition::TriFaceIsOnHull( EleRef fv, bool onHull )
{
  VTK_TRI_DDEMUX_FV(f,v,fv);
  if ( onHull )
    {
    this->Triangles[f].Flags |= v % 2 ? 0x20 : 0x10;
    }
  else
    {
    this->Triangles[f].Flags &= ~(v % 2 ? 0x20 : 0x10);
    }
}

bool vtkShoeBoxPartition::TriIsFaceOnHull( EleRef fv )
{
  VTK_TRI_DDEMUX_FV(f,v,fv);
  int flags = this->Triangles[f].Flags;
  int vmask = v % 2 ? 0x20 : 0x10;
  int rint = flags & vmask;
  bool result = rint != 0;
  return result;
  //return (this->Triangles[f].Flags & (v % 2 ? 0x20 : 0x10) != 0);
}

void vtkShoeBoxPartition::SplitBoundaryEdge( vtkIdType d, EleRef eref,
                                             vtkShoeBoxPartition::PointRef p )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "# SplitBoundaryEdge\n";
#endif // VTK_VERBOSE_SBP_FLOW

  // First adjacent face: P0 - P1 - P2
  EleRef fv012 = eref;
  PointRef p1 = this->TriDest( fv012 );

#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  PointRef p0 = this->TriOrg( fv012 );
  vtkstd::cout << "# Inserting point " 
               << p->Label << " ( "
               << p->X[0] << ", "
               << p->X[1] << ", "
               << p->X[2] << " ) on edge "
               << p0->Label << "-" 
               << p1->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE

  EleRef fv120 = this->TriENext( fv012 );
  EleRef fv122next = this->TriFNext( fv120 );
  PointRef p2 = this->TriDest( fv120 );

#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << " on face " 
               << p0->Label << " "
               << p1->Label << " "
               << p2->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE

  // Replace face P1 - P2 - P0 with face P* - P2 - P0
  int f, v;
  VTK_TRI_DEMUX_FV( f, v, fv120 );
  this->Triangles[f].EndPoints[this->VO[v]] = p;
  
#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << "  created: "
               << this->TriOrg( fv012 )->Label << " "
               << this->TriDest( fv012 )->Label << " "
               << this->TriDest( this->TriENext( fv012 ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE
  
  // Create face P* - P1 - P2
  vtkIdType tnum = this->TriCreate( p, p1, p2 );
  DofLocator& dlout( this->Offsets[d] );
  ++ dlout.NumberOfElements;
  EleRef fvP12 = VTK_TRI_MUX_FV( tnum, 0 );
  
#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << "  created: "
               << this->TriOrg( fvP12 )->Label << " "
               << this->TriDest( fvP12 )->Label << " "
               << this->TriDest( this->TriENext( fvP12 ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE
  
  // Stick in P1 - P2 - P*
  EleRef fv12P = this->TriENext( fvP12 );
  this->TriIdempotentFMerge( fv12P, fv122next );
  this->TriIdempotentFMerge( fv122next, fv12P );
  //this->TriIdempotentFMerge( this->TriSym( this->TriFNext( this->TriSym( fv120 ) ) ), fv12P ); 
  EleRef fvP21 = this->TriSym ( this->TriENext( fv12P ) );
  this->TriIdempotentFMerge( fvP21, fv120 ); 
  this->TriIdempotentFMerge( fv120, fvP21 ); 

  // Mark edges that need it
  this->TriMarkEdges( fv120, 1 );
  this->TriMarkEdges( fv12P, this->TriGetEdgeMark( fv122next ) );
  
  // Second adjacent face: P0 - P1 - P3
  EleRef fv013 = this->TriFNext( fv012 );
  EleRef fv130 = this->TriENext( fv013 );
  EleRef fv133next = this->TriFNext( fv130 );
  PointRef p3 = this->TriDest( fv130 );

#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << " on face " 
               << p0->Label << " "
               << p1->Label << " "
               << p3->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE

  // Replace face P1 - P3 - P0 with face P* - P3 - P0
  VTK_TRI_DEMUX_FV( f, v, fv130 );
  this->Triangles[f].EndPoints[this->VO[v]] = p;
  this->TriMarkEdges( fv130, 1 );
  
#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << "  created: "
               << this->TriOrg( fv013 )->Label << " "
               << this->TriDest( fv013 )->Label << " "
               << this->TriDest( this->TriENext( fv013 ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE
  
  // Create face P1 - P* - P3
  tnum = this->TriCreate( p1, p, p3 );
  ++ dlout.NumberOfElements;
  EleRef fvP13 = VTK_TRI_MUX_FV( tnum, 1 );
  
#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << "  created: "
               << this->TriOrg( fvP13 )->Label << " "
               << this->TriDest( fvP13 )->Label << " "
               << this->TriDest( this->TriENext( fvP13 ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_BND_EDGE
  
  // Stick in P1 - P3 - P*
  EleRef fv13P = this->TriENext( fvP13 );
  this->TriIdempotentFMerge( fv13P, fv133next );
  this->TriIdempotentFMerge( fv133next, fv13P );
  //this->TriIdempotentFMerge( this->TriSym( this->TriFNext( this->TriSym( fv130 ) ) ), fv13P ); 
  EleRef fvP31 = this->TriSym ( this->TriENext( fv13P ) );
  this->TriIdempotentFMerge( fvP31, fv130 ); 
  this->TriIdempotentFMerge( fv130, fvP31 ); 

  // Mark edges that need it
  this->TriMarkEdges( fv130, 1 );
  this->TriMarkEdges( fv13P, this->TriGetEdgeMark( fv133next ) );

  // Finally, attach faces along edge P* - P1
  this->TriIdempotentFMerge( fvP12, fvP13 );
  this->TriIdempotentFMerge( fvP13, fvP12 );
  
#ifdef VTK_DEBUG_SPLIT_BND_EDGE
  vtkstd::cout << vtkstd::endl;
#endif // VTK_DEBUG_SPLIT_BND_EDGE
}

void vtkShoeBoxPartition::SplitInteriorEdge( vtkIdType d, EleRef eref, 
                                             vtkShoeBoxPartition::PointRef p )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "# SplitInteriorEdge\n";
#endif // VTK_VERBOSE_SBP_FLOW

  EleRef fv012 = eref;
  PointRef p1 = this->TriDest( fv012 );

#ifdef VTK_DEBUG_SPLIT_INT_EDGE
  PointRef p0 = this->TriOrg( fv012 );
  vtkstd::cout << "# Inserting point " 
               << p->Label << " ( "
               << p->X[0] << ", "
               << p->X[1] << ", "
               << p->X[2] << " ) on edge "
               << p0->Label << "-" 
               << p1->Label << "\n";
#endif // VTK_DEBUG_SPLIT_INT_EDGE

  DofLocator& dlout( this->Offsets[d] );
  EleRef fv012next, fv120;
  EleRef fv122next, fv22next1;
  EleRef fvP12prev, fvP12, fv12P;
  EleRef fvP2first0, fvP20;
  EleRef fvP21first, fvP21;
  EleRef fvP22prev, fvP22next, fvP22second, fv22nextP;
  PointRef p2, p2next;
  int f,v;
  vtkIdType tnum;
  bool mark, first = true;
  do
    {
    fv120 = this->TriENext( fv012 );
    fv122next = this->TriFNext( fv120 );
    p2 = this->TriDest( fv120 );
    fv012next = this->TriFNext( fv012 );

#ifdef VTK_DEBUG_SPLIT_INT_EDGE
    vtkstd::cout << " on face " 
                 << p0->Label << " "
                 << p1->Label << " "
                 << p2->Label << "\n";
#endif // VTK_DEBUG_SPLIT_INT_EDGE
    
    // Replace face P1 - P2 - P0 with face P* - P2 - P0
    VTK_TRI_DEMUX_FV( f, v, fv120 );
    this->Triangles[f].EndPoints[this->VO[v]] = p;
    this->TriMarkEdges( fv120, 1 );

#ifdef VTK_DEBUG_SPLIT_INT_EDGE
    vtkstd::cout << "  created: "
                 << this->TriOrg( fv012 )->Label << " "
                 << this->TriDest( fv012 )->Label << " "
                 << this->TriDest( this->TriENext( fv012 ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_INT_EDGE
    
    // Create face P* - P1 - P2
    tnum = this->TriCreate( p, p1, p2 );
    ++ dlout.NumberOfElements;
    fvP12 = VTK_TRI_MUX_FV( tnum, 0 );
    this->TriFaceIsOnHull( fvP12, 0 );

#ifdef VTK_DEBUG_SPLIT_INT_EDGE
    vtkstd::cout << "  created: "
                 << this->TriOrg( fvP12 )->Label << " "
                 << this->TriDest( fvP12 )->Label << " "
                 << this->TriDest( this->TriENext( fvP12 ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_INT_EDGE
    
    // Stick in P1 - P2 - P*
    fv12P = this->TriENext( fvP12 );
    this->TriIdempotentFMerge( fv12P, fv122next );
    this->TriIdempotentFMerge( this->TriSym( this->TriFNext( this->TriSym( fv120 ) ) ), fv12P ); 
    this->TriMarkEdges( fv12P, this->TriGetEdgeMark( fv122next ) );
    
    // Create face P* - P2 - P2next
    p2next = this->TriDest( this->TriENext( fv012next ) );
    tnum = this->TriCreate( p, p2, p2next );
    ++ dlout.NumberOfElements;
    fvP22next = VTK_TRI_MUX_FV( tnum, 0 );
    this->TriFaceIsOnHull( fvP22next, 0 );
    this->TriMarkEdges( fvP22next, 1 );

#ifdef VTK_DEBUG_SPLIT_INT_EDGE
    vtkstd::cout << "  created: "
                 << this->TriOrg( fvP22next )->Label << " "
                 << this->TriDest( fvP22next )->Label << " "
                 << this->TriDest( this->TriENext( fvP22next ) )->Label << "\n";
#endif // VTK_DEBUG_SPLIT_INT_EDGE
    
    fvP20 = this->TriENext( fv012 );
    this->TriMarkEdges( fvP20, 1 );
    fvP21 = this->TriSym( this->TriENext2( fvP12 ) );
    this->TriMarkEdges( fvP21, 1 );
    
    if ( first ) 
      {
      fvP2first0 = fvP20;
      fvP21first = fvP21;
      fvP22second = fvP22next;
      first = false;
      }
    else
      {
      // Stick in P* - P1 - P2
      this->TriIdempotentFMerge( fvP12prev, fvP12 );
      
      // Glue the P* - P2 mill
      this->TriIdempotentFMerge( fvP22prev, fvP20 );
      this->TriIdempotentFMerge( fvP20, fvP22next );
      this->TriIdempotentFMerge( fvP22next, fvP21 );
      }

    fvP12prev = fvP12;
    fvP22prev = this->TriSym( this->TriENext2( fvP22next ) );
    
    fv22nextP = this->TriENext( fvP22next );
    fv22next1 = this->TriENext( fv122next );
    this->TriIdempotentFMerge( fv22nextP, fv22next1 );
    
    mark = this->TriGetEdgeMark( this->TriENext( this->TriSym( this->TriFNext( this->TriENext( fv120 ) ) ) ) );
    this->TriMarkEdges( fv22nextP, mark );
    this->TriMarkEdges( fv22next1, mark );
    
    fv012 = fv012next;
    } while ( fv012 != eref );

  // Glue the P* - P2first mill
  this->TriIdempotentFMerge( fvP22prev, fvP2first0 );
  this->TriIdempotentFMerge( fvP2first0, fvP22second );
  this->TriIdempotentFMerge( fvP22second, fvP21first );

#ifdef VTK_DEBUG_SPLIT_INT_EDGE
  vtkstd::cout << vtkstd::endl;
#endif // VTK_DEBUG_SPLIT_INT_EDGE
}

void vtkShoeBoxPartition::SplitEdge( vtkShoeCell* sc, vtkIdType d, EleRef eref, double r[3], int face )
{
  vtkShoeBoxPartition::Point p;

  if ( face < 0 ) // this occurs iff the DOF is a volume DOF
    {
    for ( int i = 0; i < 3; ++ i ) p.R[i] = r[i];
    sc->EvaluateLocation( -1, p.R, p.X ); // This is p.X = Xi( p.R ), -1 is litter
    }
  else // this occurs iff the DOF is a face DOF
    {
    p.R[2] = r[0]; 
    p.R[3] = r[1];
    p.R[4] = r[2];
    sc->GetMetaData()->ProjectCoordsToFace( face, p.R + 2 );
    sc->GetMetaData()->TransformFaceCoordsToStorage( sc->GetFacePermutation( face ), p.R );
    sc->EvaluateLocation( -1, r, p.X ); // This is p.X = Xi( p.R ), -1 is litter
    }

  p.DOF = d;
  p.Label = this->NextPointLabel;
  p.Element = -1;
  vtkstd::pair<PointRef,bool> result = (*this->Points)[d].insert( p );
  if ( ! result.second ) return;
  ++ this->NextPointLabel;

  if ( face < 0 ) // this occurs iff the DOF is a volume DOF
    {
    this->SplitInteriorEdge( d, eref, result.first );
    }
  else // this occurs iff the DOF is a face DOF
    {
    this->SplitBoundaryEdge( d, eref, result.first );
    }
}

void vtkShoeBoxPartition::StarTriangle( vtkIdType d, vtkIdType i,
                                        vtkShoeBoxPartition::PointRef p )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "# StarTriangle\n";
#endif // VTK_VERBOSE_SBP_FLOW

  DofLocator& dlout( this->Offsets[d] );
  PointRef* pr = this->Triangles[i].EndPoints;

  // Get a hold of the neighbours
  EleRef fv01P = VTK_TRI_MUX_FV( i, 0 );
  EleRef fvN01 = this->TriFNext( fv01P );
  EleRef fvN12 = this->TriFNext( this->TriENext ( fv01P ) );
  EleRef fvN20 = this->TriFNext( this->TriENext2 ( fv01P ) );

  // Replace first vertex of first facet with the point to be inserted
  this->TriFDel( fv01P );
  pr[2] = p;
  
  // Create new facets
  vtkIdType tnum;
  dlout.NumberOfElements += 2;

  tnum = this->TriCreate( pr[1], pr[2], p );
  EleRef fv12P = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fv12P, 1 );
  this->TriFaceIsOnHull( this->TriSym( fv12P ), 1 );

  tnum = this->TriCreate( pr[2], pr[0], p );
  EleRef fv20P = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fv20P, 1 );
  this->TriFaceIsOnHull( this->TriSym( fv20P ), 1 );

  // Glue to the neighbours
  this->TriIdempotentFMerge( fvN01, fv01P );
  this->TriIdempotentFMerge( fvN12, fv12P );
  this->TriIdempotentFMerge( fvN20, fv20P );
  
  // Glue inner edges together and mark them for subsequent C.P. search
  EleRef fv0P1 = this->TriSym( this->TriENext2( fv01P ) );
  this->TriIdempotentFMerge( this->TriENext( fv20P ), fv0P1 );
  this->TriMarkEdges( fv0P1, 1 );

  EleRef fv1P2 = this->TriSym( this->TriENext2( fv12P ) );
  this->TriIdempotentFMerge( this->TriENext( fv01P ), fv1P2 );
  this->TriMarkEdges( fv1P2, 1 );

  EleRef fv2P0 = this->TriSym( this->TriENext2( fv20P ) );
  this->TriIdempotentFMerge( this->TriENext( fv12P ), fv2P0 );
  this->TriMarkEdges( fv2P0, 1 );
}

void vtkShoeBoxPartition::StarTetrahedron( vtkIdType d, EleRef fv012, 
                                           vtkShoeBoxPartition::PointRef p )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "# StarTetrahedron\n";
#endif // VTK_VERBOSE_SBP_FLOW

  EleRef fv120 = this->TriENext( fv012 );
  EleRef fv201 = this->TriENext( fv120 );
  EleRef fv013 = this->TriFNext( fv012 );
  EleRef fv103 = this->TriSym( fv013 );
  EleRef fv031 = this->TriENext( fv103 );
  EleRef fv310 = this->TriENext( fv031 );
  EleRef fv321 = this->TriSym( this->TriENext( this->TriFNext( fv120 ) ) );
    
  PointRef pr0 = this->TriOrg( fv012 );
  PointRef pr1 = this->TriOrg( fv120 );
  PointRef pr2 = this->TriOrg( fv201 );
  PointRef pr3 = this->TriDest( fv031 );
    
  vtkIdType tnum;
  EleRef fnum;
    
  // Create new facets, get a hold of them, attach them to the mother tet.
  tnum = this->TriCreate( pr0, pr1, p );
  fnum = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fnum, 0 );
  this->TriFaceIsOnHull( this->TriSym( fnum ), 0 );
  this->TriInvolutiveFMerge( fv012, fnum );
  EleRef fv1P0 = this->TriENext( fnum );
  this->TriFaceIsOnHull( fv1P0, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv1P0 ), 0 );
  EleRef fvP01 = this->TriENext( fv1P0 );
  this->TriFaceIsOnHull( fvP01, 0 );
  this->TriFaceIsOnHull( this->TriSym( fvP01 ), 0 );
    
  tnum = this->TriCreate( pr1, pr2, p );
  fnum = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fnum, 0 );
  this->TriFaceIsOnHull( this->TriSym( fnum ), 0 );
  this->TriInvolutiveFMerge( fv120, fnum );
  EleRef fv2P1 = this->TriENext( fnum );
  this->TriFaceIsOnHull( fv2P1, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv2P1 ), 0 );
  EleRef fv1P2 = this->TriSym( this->TriENext ( fv2P1 ) );
  this->TriFaceIsOnHull( fv1P2, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv1P2 ), 0 ); // FIXME
    
  tnum = this->TriCreate( pr2, pr0, p );
  fnum = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fnum, 0 );
  this->TriFaceIsOnHull( this->TriSym( fnum ), 0 );
  this->TriInvolutiveFMerge( fv201, fnum );
  EleRef fv2P0 = this->TriENext( this->TriSym( fnum ) );
  this->TriFaceIsOnHull( fv2P0, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv2P0 ), 0 ); // FIXME
  EleRef fvP02 = this->TriENext( fv2P0 );
  this->TriFaceIsOnHull( fvP02, 0 );
  this->TriFaceIsOnHull( this->TriSym( fvP02 ), 0 );
    
  tnum = this->TriCreate( pr0, pr3, p );
  fnum = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fnum, 0 );
  this->TriFaceIsOnHull( this->TriSym( fnum ), 0 );
  this->TriInvolutiveFMerge( fv031, fnum );
  EleRef fv3P0 = this->TriENext ( fnum );
  this->TriFaceIsOnHull( fv3P0, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv3P0 ), 0 );
  EleRef fvP03 = this->TriENext ( fv3P0 );
  this->TriFaceIsOnHull( fvP03, 0 );
  this->TriFaceIsOnHull( this->TriSym( fvP03 ), 0 );
    
  tnum = this->TriCreate( pr3, pr1, p );
  fnum = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fnum, 0 );
  this->TriFaceIsOnHull( this->TriSym( fnum ), 0 );
  this->TriInvolutiveFMerge( fv310, fnum );
  EleRef fv1P3 = this->TriENext ( fnum );
  this->TriFaceIsOnHull( fv1P3, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv1P3 ), 0 );
  EleRef fv3P1 = this->TriSym( this->TriENext( fv1P3 ) );
  this->TriFaceIsOnHull( fv3P1, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv3P1 ), 0 ); // FIXME

  tnum = this->TriCreate( pr3, pr2, p );
  fnum = VTK_TRI_MUX_FV( tnum, 0 );
  this->TriFaceIsOnHull( fnum, 0 );
  this->TriFaceIsOnHull( this->TriSym( fnum ), 0 );
  this->TriInvolutiveFMerge( fv321, fnum );
  EleRef fv2P3 = this->TriENext ( fnum );
  this->TriFaceIsOnHull( fv2P3, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv2P3 ), 0 );
  EleRef fv3P2 = this->TriSym( this->TriENext( fv2P3 ) );
  this->TriFaceIsOnHull( fv3P2, 0 );
  this->TriFaceIsOnHull( this->TriSym( fv3P2 ), 0 ); // FIXME

  this->Offsets[d].NumberOfElements += 6;
    
  // Glue inner edges together
  this->TriIdempotentFMerge( fvP03, fvP02 );
  this->TriIdempotentFMerge( fvP02, fvP01 );
#ifdef VTK_DEBUG_INS_PT_TET
  vtkstd::cout << this->TriOrg( fvP03 )->Label << " "
               << this->TriDest( fvP03 )->Label << " "
               << this->TriOrg( fvP02 )->Label << " "
               << this->TriDest( fvP02 )->Label << " "
               << this->TriOrg( fvP01 )->Label << " "
               << this->TriDest( fvP01 )->Label << "\n";
#endif // VTK_DEBUG_INS_PT_TET
    
  this->TriIdempotentFMerge( fv1P0, fv1P3 );
  this->TriIdempotentFMerge( fv1P3, fv1P2 );
#ifdef VTK_DEBUG_INS_PT_TET
  vtkstd::cout << this->TriOrg( fv1P0 )->Label << " "
               << this->TriDest( fv1P0 )->Label << " "
               << this->TriOrg( fv1P3 )->Label << " "
               << this->TriDest( fv1P3 )->Label << " "
               << this->TriOrg( fv1P2 )->Label << " "
               << this->TriDest( fv1P2 )->Label << "\n";
#endif // VTK_DEBUG_INS_PT_TET
    
  this->TriIdempotentFMerge( fv2P1, fv2P3 );
  this->TriIdempotentFMerge( fv2P3, fv2P0 );
#ifdef VTK_DEBUG_INS_PT_TET
  vtkstd::cout << this->TriOrg( fv2P1 )->Label << " "
               << this->TriDest( fv2P1 )->Label << " "
               << this->TriOrg( fv2P3 )->Label << " "
               << this->TriDest( fv2P3 )->Label << " "
               << this->TriOrg( fv2P0 )->Label << " "
               << this->TriDest( fv2P0 )->Label << "\n";
#endif // VTK_DEBUG_INS_PT_TET
    
  this->TriIdempotentFMerge( fv3P2, fv3P1 );
  this->TriIdempotentFMerge( fv3P1, fv3P0 );
#ifdef VTK_DEBUG_INS_PT_TET
  vtkstd::cout << this->TriOrg( fv3P2 )->Label << " "
               << this->TriDest( fv3P2 )->Label << " "
               << this->TriOrg( fv3P1 )->Label << " "
               << this->TriDest( fv3P1 )->Label << " "
               << this->TriOrg( fv3P0 )->Label << " "
               << this->TriDest( fv3P0 )->Label << "\n";
#endif // VTK_DEBUG_INS_PT_TET
    
  // Finally, mark inner edges and faces for subsequent C.P. search
  this->TriMarkEdges( fvP01, 1 );
  this->TriMarkEdges( fvP02, 1 );
  this->TriMarkEdges( fvP03, 1 );
    
  this->TriMarkEdges( fv1P0, 1 );
  this->TriMarkEdges( fv1P2, 1 );
  this->TriMarkEdges( fv1P3, 1 );
    
  this->TriMarkEdges( fv2P0, 1 );
  this->TriMarkEdges( fv2P1, 1 );
  this->TriMarkEdges( fv2P3, 1 );
    
  this->TriMarkEdges( fv3P0, 1 );
  this->TriMarkEdges( fv3P1, 1 );
  this->TriMarkEdges( fv3P2, 1 );
    
  this->TriMarkFace( fvP01, 1 );
  this->TriMarkFace( fv2P1, 1 );
  this->TriMarkFace( fv1P3, 1 );
  this->TriMarkFace( fv3P2, 1 );
  this->TriMarkFace( fvP02, 1 );
  this->TriMarkFace( fvP03, 1 );
}

void vtkShoeBoxPartition::InsertPointIntoTriangulation( vtkShoeCell* sc, vtkIdType* conn, 
                                                        int dof, double* r )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "## InsertPointIntoTriangulation:\n";
#endif // VTK_VERBOSE_SBP_FLOW

  vtkShoeBoxPartition::Point p;
  double px = p.R[0] = r[2];
  double py = p.R[1] = r[3];
  int ne = sc->GetMetaData()->NumberOfBoundaries[1];

  sc->GetMetaData()->EmbedFaceCoords( dof - ne, p.R ); // p.R is modified as a side effect
  sc->EvaluateLocation( -1, p.R, p.X ); // This is p.X = Xi( p.R ), -1 is litter
  p.R[0] = r[0]; // Restore storage order face coordinates
  p.R[1] = r[1];
  p.R[2] = r[2]; // Restore cell local face coordinates
  p.R[3] = r[3]; 

  vtkIdType d = conn[dof];
  p.DOF = d;
  p.Label = this->NextPointLabel;
  p.Element = -1;
  vtkstd::pair<PointRef,bool> result = (*this->Points)[d].insert( p );
  if ( ! result.second ) return;
  ++ this->NextPointLabel;

  bool triNotFound = true;
  DofLocator dl = this->Offsets[d];
  vtkIdType dlEnd = dl.Start + dl.NumberOfElements;
  for ( vtkIdType t = dl.Start; t < dlEnd; ++ t )
    {
    PointRef* pr = this->Triangles[t].EndPoints;

#ifdef VTK_DEBUG_INSERT_POINT_IN_TRI
    vtkstd::cout << "Triangle " << t  << " :\n"
                 << " " << pr[0]->Label << " ( " << pr[0]->R[2] << ", " << pr[0]->R[3] << " )"
                 << " [ " << pr[0]->X[0] << ", " << pr[0]->X[1] << ", " << pr[0]->X[2] << " ]\n"
                 << " " << pr[1]->Label << " ( " << pr[1]->R[2] << ", " << pr[1]->R[3] << " )"
                 << " [ " << pr[1]->X[0] << ", " << pr[1]->X[1] << ", " << pr[1]->X[2] << " ]\n"
                 << " " << pr[2]->Label << " ( " << pr[2]->R[2] << ", " << pr[2]->R[3] << " )"
                 << " [ " << pr[2]->X[0] << ", " << pr[2]->X[1] << ", " << pr[2]->X[2] << " ]\n";
#endif // VTK_DEBUG_INSERT_POINT_IN_TRI

#if 0 // use local parametric coordinates
    double AB[2], AC[2];
    AB[0] = pr[1]->R[2] - pr[0]->R[2];
    AB[1] = pr[1]->R[3] - pr[0]->R[3];
    AC[0] = pr[2]->R[2] - pr[0]->R[2];
    AC[1] = pr[2]->R[3] - pr[0]->R[3];
    double ABC = AB[0] * AC[1] - AB[1] * AC[0];
    if ( ! ABC )
      {
      vtkWarningMacro( "Triangle " << i << " is degenerate." );
      return;
      }

    double AM[2];
    AM[0] = px - pr[0]->R[2];
    AM[1] = py - pr[0]->R[3];
    double ABM = AB[0] * AM[1] - AB[1] * AM[0];
    if ( ABM * ABC < 0. ) continue;
    double AMC = AM[0] * AC[1] - AM[1] * AC[0];
    if ( AMC * ABC < 0. ) continue;

    double BC[2];
    BC[0] = pr[2]->R[2] - pr[1]->R[2];
    BC[1] = pr[2]->R[3] - pr[1]->R[3];
    double BCM = BC[0] * ( py - pr[1]->R[3] ) - BC[1] * ( px - pr[1]->R[2] );
    if ( BCM * ABC < 0. ) continue;

    // If we're here, then this means that the point is either in the interior
    // or on one of the edges of the triangle. 
    tetNotFound = false;
      
#ifdef VTK_DEBUG_INSERT_POINT_IN_TRI
    vtkstd::cout << "Point belongs to tri "
                 << pr[0]->Label << " "
                 << pr[1]->Label << " "
                 << pr[2]->Label << " "
                 << " ABC= " << ABC
                 << " ABM= " << ABM
                 << " AMC= " << AMC
                 << " BCM= " << BCM
                 << vtkstd::endl;
#endif // VTK_DEBUG_INSERT_POINT_IN_TRI

    // We must check if it makes sense to add one of the neighbours.
    
    double ABMtoABC = fabs ( ABM / ABC );
    int evilTriangle = ABMtoABC < this->NullTriangleThreshold ? 1 : 0;
    
    double BCMtoABC = fabs ( BCM / ABC );
    evilTriangle += BCMtoABC < this->NullTriangleThreshold ? 2 : 0;
    if ( evilTriangle == 3 ) continue;
    
    double AMCtoABC = 1. - ( ABMtoABC + BCMtoABC );
    evilTriangle += AMCtoABC < this->NullTriangleThreshold ? 4 : 0;
    if ( ( evilTriangle == 5 ) || ( evilTriangle == 6 ) ) continue;
    
#endif // use local parametric coordinates

#if 1 // use world coordinates
    double AB[3], AC[3], ABAC[3];
    for ( int ii = 0; ii < 3; ++ ii )
      {
      AB[ii] = pr[1]->X[ii] - pr[0]->X[ii];
      AC[ii] = pr[2]->X[ii] - pr[0]->X[ii];
      }
    vtkMath::Cross( AB, AC, ABAC );
    double ABC2 = vtkMath::Dot( ABAC, ABAC );
    if ( ! ABC2 )
      {
      vtkWarningMacro( "Triangle " << t << " is degenerate." );
      return;
      }

    double AM[3], BC[3], BM[3];
    for ( int ii = 0; ii < 3; ++ ii )
      {
      AM[ii] = p.X[ii] - pr[0]->X[ii];
      BC[ii] = pr[2]->X[ii] - pr[1]->X[ii];
      BM[ii] = p.X[ii] - pr[1]->X[ii];
      }
    double ABAM[3], AMAC[3], BCBM[3];
    vtkMath::Cross( AB, AM, ABAM );
    if ( vtkMath::Dot( ABAM, ABAC ) < 0. ) continue;
    vtkMath::Cross( AM, AC, AMAC );
    if ( vtkMath::Dot( AMAC, ABAC ) < 0. ) continue;
    vtkMath::Cross( BC, BM, BCBM );
    if ( vtkMath::Dot( BCBM, ABAC ) < 0. ) continue;

    double ABM2 = vtkMath::Dot( ABAM, ABAM );
    double AMC2 = vtkMath::Dot( AMAC, AMAC );
    double BCM2 = vtkMath::Dot( BCBM, BCBM );

    // If we're here, then this means that the point is either in the interior
    // or on one of the edges of the triangle. 
    triNotFound = false;
      
#ifdef VTK_DEBUG_INSERT_POINT_IN_TRI
    vtkstd::cout << "Point belongs to tri "
                 << pr[0]->Label << " "
                 << pr[1]->Label << " "
                 << pr[2]->Label << " "
                 << " ABC= " << sqrt( ABC2 )
                 << " ABM= " << sqrt( ABM2 )
                 << " AMC= " << sqrt( AMC2 )
                 << " BCM= " << sqrt( BCM2 )
                 << vtkstd::endl;
#endif // VTK_DEBUG_INSERT_POINT_IN_TRI

    // We must check if it makes sense to add one of the neighbours.
    
    double ABM2toABC2 = ABM2 / ABC2;
    int evilTriangle = ABM2toABC2 < this->NullTriangleThreshold2 ? 1 : 0;
    
    double BCM2toABC2 = BCM2 / ABC2;
    evilTriangle += BCM2toABC2 < this->NullTriangleThreshold2 ? 2 : 0;
    if ( evilTriangle == 3 ) continue;
    
    double AMC2toABC2 = AMC2 / ABC2;
    evilTriangle += AMC2toABC2 < this->NullTriangleThreshold2 ? 4 : 0;
    if ( ( evilTriangle == 5 ) || ( evilTriangle == 6 ) ) continue;
    
#endif // use world coordinates

    if ( evilTriangle )
      {
      // If we're here, then this means that the point is close to only
      // one edge of the triangle. We need to add the neighboring triangle
      // to the cavity.
      
      if ( evilTriangle == 1 ) evilTriangle = 0;

      EleRef fv = this->Triangles[t].Neighbors[evilTriangle];
      int f, v;
      VTK_TRI_DEMUX_FV(f,v,fv);

#ifdef VTK_DEBUG_INSERT_POINT_IN_TRI
      vtkstd::cout << "  -> Critical point close to edge " 
                   << this->TriOrg( fv )->Label
                   << " - "
                   << this->TriDest( fv )->Label
                   << ".\n";
#endif // VTK_DEBUG_INSERT_POINT_IN_TRI
      
      this->SplitBoundaryEdge( d, fv, result.first );
      break;
      }
    else
      {
      // If we're here, then this means that the point is interior to
      // a tetrahedron.
      this->StarTriangle( d, t, result.first );
      break;
      }
    }

  if ( triNotFound )
    {
    vtkErrorMacro("The point has not been found to belong to a triangle.");
    }

  this->Dump();
}

void vtkShoeBoxPartition::InsertPointIntoTetrahedralization( vtkShoeCell* sc, vtkIdType* conn, 
                                                             int dof, double* r )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "## InsertPointIntoTetrahedralization:\n";
#endif // VTK_VERBOSE_SBP_FLOW

  this->MadITellYou();

  Point p;
  for ( int i = 0; i < 3; ++ i) p.R[i] = r[i];
  sc->EvaluateLocation( -1, p.R, p.X ); // This is p.X = Xi( p.R ), -1 is litter

  vtkIdType d = conn[dof];
  p.DOF = d;
  p.Label = this->NextPointLabel;
  p.Element = -1;
  vtkstd::pair<PointRef,bool> result = (*this->Points)[d].insert( p );
  if ( ! result.second ) return;
  ++ this->NextPointLabel;

  bool tetNotFound = true;
  EleRef fv0;
  vtkShoeBoxPartitionIterator* sbpit = vtkShoeBoxPartitionIterator::New();
  for ( sbpit->Begin( this, d ); ! sbpit->IsAtEnd(); sbpit->Next() )
    {
    fv0 = sbpit->GetCurrentEdgeFacet();
    EleRef fv1E = this->TriENext( fv0 );
    EleRef fv1F = this->TriFNext( fv0 );

    PointRef pr[4]; 
    pr[0] = this->TriOrg( fv0 );
    pr[1] = this->TriOrg( fv1E );
    pr[2] = this->TriOrg( this->TriENext( fv1E ) );
    pr[3] = this->TriDest( this->TriENext( fv1F ) );

    double AB[3], AC[3], AD[3];
    for ( int i = 0; i < 3; ++ i)
      {
      AB[i] = pr[1]->R[i] - pr[0]->R[i];
      AC[i] = pr[2]->R[i] - pr[0]->R[i];
      AD[i] = pr[3]->R[i] - pr[0]->R[i];
      }
    double ABCD = vtkMath::Determinant3x3( AB, AC, AD );

    if ( ! ABCD ) 
      {
      vtkWarningMacro( "Tetrahedron " << pr[0]->Label << "-" << pr[1]->Label << "-" << pr[2]->Label << "-" << pr[3]->Label << " is degenerate." );
      return;
      }

    if ( ABCD < 0 ) 
      {
      vtkWarningMacro( "Tetrahedron " << pr[0]->Label << "-" << pr[1]->Label << "-" << pr[2]->Label << "-" << pr[3]->Label << " is inverted." );
      return;
      }

    double AM[3];
    for ( int i = 0; i < 3; ++ i) AM[i] = r[i] - pr[0]->R[i];
    double ABCM = vtkMath::Determinant3x3( AB, AC, AM );
    if ( ABCM * ABCD < 0. ) continue;
    double ADBM = vtkMath::Determinant3x3( AD, AB, AM );
    if ( ADBM * ABCD < 0. ) continue;
    double ACDM = vtkMath::Determinant3x3( AC, AD, AM );
    if ( ACDM * ABCD < 0. ) continue;

    double BC[3], BD[3], BM[3];
    for ( int i = 0; i < 3; ++ i)
      {
      BC[i] = pr[2]->R[i] - pr[1]->R[i];
      BD[i] = pr[3]->R[i] - pr[1]->R[i];
      BM[i] = r[i] - pr[1]->R[i];
      }
    double BDCM = vtkMath::Determinant3x3( BD, BC, BM );
    if ( BDCM * ABCD < 0. ) continue;

    // If we're here, then this means that the point is either in the interior
    // or on the boundary of the tetrahedron
    tetNotFound = false;
      
#ifdef VTK_DEBUG_INSERT_POINT_IN_TET
    vtkstd::cout << "Point belongs to tet "
                 << pr[0]->Label << " "
                 << pr[1]->Label << " "
                 << pr[2]->Label << " "
                 << pr[3]->Label << " "
                 << " ABCD= " << ABCD
                 << " ABCM= " << ABCM
                 << " ADBM= " << ADBM
                 << " ACDM= " << ACDM
                 << " BDCM= " << BDCM
                 << vtkstd::endl;
#endif // VTK_DEBUG_INSERT_POINT_IN_TET

    // We must check if it makes sense to add one ore more of the neighbours.

    double ABCMtoABCD = fabs ( ABCM / ABCD );
    int evilTet = ABCMtoABCD < this->NullTetThreshold ? 1 : 0;
    
    double ADBMtoABCD = fabs ( ADBM / ABCD );
    evilTet += ADBMtoABCD < this->NullTetThreshold ? 2 : 0;
    
    double ACDMtoABCD = fabs ( ACDM / ABCD );
    evilTet += ACDMtoABCD < this->NullTetThreshold ? 4 : 0;
    if ( evilTet == 7 ) continue;
    
    double BDCMtoABCD = 1. - ( ABCMtoABCD + ADBMtoABCD + ACDMtoABCD );
    evilTet += BDCMtoABCD < this->NullTetThreshold ? 8 : 0;
    if ( ( evilTet == 11 ) || ( evilTet == 13 ) || ( evilTet == 14 ) ) continue;

    if ( evilTet )
      {
      // If we're here, then this means that the point is close to one or
      // two faces of the tetrahedron. We need to add the neighboring 
      // tetrahedr(on/a) to the cavity.
      
      if ( ( evilTet == 1 ) || ( evilTet == 2 ) || 
           ( evilTet == 4 ) || ( evilTet == 8 ) )
        {

#ifdef VTK_DEBUG_INSERT_POINT_IN_TET
        vtkstd::cout << "  -> Critical point close to one face.\n";
#endif // VTK_DEBUG_INSERT_POINT_IN_TET
        
        // FIXME: complete implementation
        vtkErrorMacro("Case not implemented.\n");
        break;
        }
      else if ( ( evilTet == 3 ) || ( evilTet == 5 ) || 
                ( evilTet == 6 ) || ( evilTet == 9 ) ||
                ( evilTet == 10 ) || ( evilTet == 12 ) )
        {

        switch ( evilTet )
          {
          case 3:
            // Do nothing, fv0 is the edge to be split
            break;
          case 9:
            fv0 = fv1E;
            break;
          case 5:
            fv0 = this->TriENext( fv1E );
            break;
          case 6:
            fv0 = this->TriENext2( fv1F );
            break;
          case 10:
            fv0 = this->TriENext( fv1F );
            break;
          case 12:
            fv0 = this->TriENext( this->TriFNext( fv1E ) );
            break;
          default:
            vtkErrorMacro("Incorrect edge bitcode ("<< evilTet << ")."); 
            break;
          }

#ifdef VTK_DEBUG_INSERT_POINT_IN_TET
        vtkstd::cout << "  -> Critical point close to edge " 
                     << this->TriOrg( fv0 )->Label
                     << " - "
                     << this->TriDest( fv0 )->Label
                     << ".\n";
#endif // VTK_DEBUG_INSERT_POINT_IN_TET
          
        this->SplitInteriorEdge( d, fv0, result.first );
        break;
        }
      }
    else
      {
      // If we're here, then this means that the point is interior to
      // a tetrahedron.
      this->StarTetrahedron( d, fv0, result.first );
      break;
      }
    }

  if ( tetNotFound )
    {
    vtkErrorMacro("The point has not been found to belong to a tetrahedron.");
    }

  this->CheckTetOrientations( d );
  this->Dump();
}

void vtkShoeBoxPartition::AddPointToFaceLoop( vtkShoeCell* sc, vtkIdType* conn, 
                                              int dof, double* r )
{
  Point p;
  // Convert r into cell-local order for EvaluateLocation...
  p.R[0] = r[0];
  sc->GetRecord().ParametricEdgeCoordsFromStorageOrder( dof, p.R );
  sc->GetMetaData()->EmbedEdgeCoords( dof, p.R );
  sc->EvaluateLocation( -1, p.R, p.X );
  // ...but store as passed, not cell-local order.
  for ( int i = 0; i < 4; ++ i ) p.R[i] = r[i];
  vtkIdType d = conn[dof];
  p.DOF = d;
  p.Label = this->NextPointLabel;
  p.Element = -1;

  vtkstd::pair<vtkShoeBoxPartition::PointRef,bool> result;
  result = (*this->Points)[d].insert( p );
  this->EdgePoints.push_back( result.first );
  if ( result.second ) ++ this->NextPointLabel;
}

void vtkShoeBoxPartition::TriangulateFaceLoop( vtkShoeCell* sc, vtkIdType* conn, int ne, 
                                               int dof, double* r, int isBdy )
{
#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "\n### TriangulateFaceLoop (dof " << dof << ")\n";
#endif // VTK_VERBOSE_SBP_FLOW

  int i;
  Point p;
  // Convert r into cell-local order for EvaluateLocation...
  p.R[0] = r[2];
  p.R[1] = r[3];
  if ( isBdy )
    {
    sc->GetMetaData()->EmbedFaceCoords( dof - ne, p.R );
    }
  sc->EvaluateLocation( -1, p.R, p.X );
  // ...but store in DOF-node order, not cell-local order.
  p.R[0] = r[0];
  p.R[1] = r[1];
  p.R[2] = r[2];
  p.R[3] = r[3];
  vtkIdType d = conn[dof];
  p.DOF = d;
  p.Label = this->NextPointLabel;

  vtkstd::pair<vtkShoeBoxPartition::PointRef,bool> result;
  result = (*this->Points)[d].insert( p );
  if ( result.second ) ++ this->NextPointLabel;
  //this->EdgePoints->push_back( result.first );

  DofLocator dl;
  dl.Start = this->Triangles.size();
  dl.NumberOfElements = this->EdgePoints.size();

  EleRef tnum;
  if ( sc->GetFacePermutation( dof - ne ) % 2 )
    {
    // odd face permutations mean the face is reversed from storage order; triangles must be flipped
    (EleRef&)result.first->Element = -1;
    for ( i = 1; i < dl.NumberOfElements; ++ i )
      {
      tnum = this->TriCreate( result.first, this->EdgePoints[i-1], this->EdgePoints[i] );
      this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,0), 1 ); // Edge 0-1
      this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,4), 1 ); // Edge 2-0
      this->TriFaceIsOnHull( VTK_TRI_MUX_FV(tnum,0), 0 );
      EleRef last = result.first->Element;
      if ( last >= 0 )
        this->TriInvolutiveFMerge( last, VTK_TRI_MUX_FV(tnum,0) );
      (EleRef&)result.first->Element = tnum; // use common vertex to store previous results
      }
    tnum = this->TriCreate( result.first, this->EdgePoints[dl.NumberOfElements-1], this->EdgePoints[0] ); // close the loop
    this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,0), 1 ); // Edge 0-1
    this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,4), 1 ); // Edge 2-0
    this->TriFaceIsOnHull( VTK_TRI_MUX_FV(tnum,0), 0 );
    this->TriInvolutiveFMerge( result.first->Element, VTK_TRI_MUX_FV(tnum,0) );
    this->TriInvolutiveFMerge( VTK_TRI_MUX_FV(tnum,5), VTK_TRI_MUX_FV(dl.Start,0) );
    }
  else
    {
    // even face permutations mean the face is only rotated from storage order; no need to flip triangles.
    (EleRef&)result.first->Element = -1;
    for ( i = dl.NumberOfElements - 1; i > 0; --i )
      {
      tnum = this->TriCreate( result.first, this->EdgePoints[i], this->EdgePoints[i-1] );
      this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,0), 1 ); // Edge 0-1
      this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,4), 1 ); // Edge 2-0
      this->TriFaceIsOnHull( VTK_TRI_MUX_FV(tnum,0), 0 );
      EleRef last = result.first->Element;
      if ( last >= 0 )
        this->TriInvolutiveFMerge( last, VTK_TRI_MUX_FV(tnum,0) );
      (EleRef&)result.first->Element = VTK_TRI_MUX_FV(tnum,5); // use common vertex to store previous results
      }
    tnum = this->TriCreate( result.first, this->EdgePoints[0], this->EdgePoints[dl.NumberOfElements-1] ); // close the loop
    this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,0), 1 ); // Edge 0-1
    this->TriMarkEdges( VTK_TRI_MUX_FV(tnum,4), 1 ); // Edge 2-0
    this->TriFaceIsOnHull( VTK_TRI_MUX_FV(tnum,0), 0 );
    this->TriInvolutiveFMerge( result.first->Element, VTK_TRI_MUX_FV(tnum,0) );
    this->TriInvolutiveFMerge( VTK_TRI_MUX_FV(tnum,5), VTK_TRI_MUX_FV(dl.Start,0) );
    }
  // The face-point that is the center of the triangle fan should have a reference to
  // one element of the triangulation:
  (EleRef&)this->Triangles[tnum].EndPoints[0]->Element = VTK_TRI_MUX_FV(tnum,0);
  // Clear the EdgePoints for the next face to be processed.
  this->EdgePoints.clear();
  // Insert the DofOffset record into the map associating a DOF ID to triangles
  this->Offsets[conn[dof]] = dl;
}

void vtkShoeBoxPartition::CopyPointsToDof( vtkShoeBoxPartition::PointRef npts[3],
  vtkShoeCell* sc, vtkIdType* conn, int ne, int nf, int vdof,
  vtkShoeBoxPartition* bdyTri, vtkIdType facedof, int i )
{
  (void)nf; // FIXME: Do we really not need nf?
  Point p;
  PointRef oldpt;
  i += bdyTri->Offsets[facedof].Start;
  vtkShoeBoxPartition::Triangle& tri( bdyTri->Triangles[i] );

  for ( int j = 0; j < 3; ++j )
    {
    oldpt = tri.EndPoints[j];
    int k = 0;
    while ( k < vdof && conn[k] != oldpt->DOF ) ++k;
    if ( k == vdof )
      {
      vtkErrorMacro("Bad point in boundary triangulation, DOF " << oldpt->DOF << " isn't in connectivity array.");
      continue;
      }
    if ( k < ne )
      {
      p.R[0] = oldpt->R[0];
      sc->GetRecord().ParametricEdgeCoordsFromStorageOrder( k, p.R );
      sc->GetMetaData()->EmbedEdgeCoords( k, p.R );
      }
    else
      {
      p.R[0] = oldpt->R[0];
      p.R[1] = oldpt->R[1];
      k -= ne;
      sc->GetRecord().ParametricFaceCoordsFromStorageOrder( k, p.R );
      sc->GetMetaData()->EmbedFaceCoords( k, p.R );
      }
    for ( int x = 0; x < 3; ++x ) p.X[x] = oldpt->X[x];
    p.Label = this->NextPointLabel;
    p.DOF = conn[vdof];
    p.Element = -1;

    vtkstd::pair<PointRef,bool> result = (*this->Points)[conn[vdof]].insert( p );
    if ( result.second ) ++ this->NextPointLabel;
    npts[j] = result.first;
    }
}

void vtkShoeBoxPartition::StarInterior( vtkShoeBoxPartition* bdyTri,
  vtkShoeCell* sc, vtkIdType* conn, int ne, int nf, double* interiorPt )
{
  // This is going to be dog slow, but it will work.
  //typedef vtkstd::pair< PointRef, PointRef > Edge;
  //vtkstd::map< Edge, vtkstd::set<EleRef> > edgeMap;

#ifdef VTK_VERBOSE_SBP_FLOW
  vtkstd::cout << "### StarInterior\n";
#endif // VTK_VERBOSE_SBP_FLOW

  // Insert an entry in this->Offsets to record triangles associated with this DOF node.
  vtkIdType vdof = conn[ne + nf];
  vtkIdType d = conn[vdof];
  DofLocator tmp;
  tmp.Start = this->Triangles.size();
  tmp.NumberOfElements = 0;
  this->Offsets[d] = tmp;
  DofLocator& eleOff( this->Offsets[d] );

  // Insert the point interior to the cell
  vtkShoeBoxPartition::PointRef p[4];
  vtkShoeBoxPartition::Point ptmp;
  
  ptmp.DOF = d;
  ptmp.Label = this->NextPointLabel;
  ptmp.Element = -1;
  int m;
  for ( m = 0; m < 3; ++m )
    ptmp.R[m] = interiorPt[m];
  sc->EvaluateLocation( -1, ptmp.R, ptmp.X );
  vtkstd::pair<PointRef,bool> result = (*this->Points)[d].insert( ptmp );
  if ( result.second ) ++ this->NextPointLabel;
  else 
    {
    vtkErrorMacro("StarInterior: inserted point " << result.first->Label << " that already exists in DOF " << vdof << ".");
    return;
    }
  p[0] = result.first;

  vtkIdType b;
  for ( b = 0; b < nf; ++b )
    {
    vtkIdType i;
    vtkIdType faceDof = ne + b;
    vtkIdType nt = bdyTri->GetNumberOfElements( conn[faceDof] );
    int facePerm = sc->GetFacePermutation( b );
    for ( i = 0; i < nt; ++ i )
      {
      //vtkstd::cout << "About to insert Boundary " << b << " Tri " << i << vtkstd::endl;
      this->CopyPointsToDof( p + 1, sc, conn, ne, nf, vdof, bdyTri, conn[faceDof], i );

      double e0[3];
      double e1[3];
      double e2[3];
      int h;
      for ( h = 0; h < 3; ++ h )
        {
        e0[h] = p[1]->R[h] - p[0]->R[h];
        e1[h] = p[2]->R[h] - p[0]->R[h];
        e2[h] = p[3]->R[h] - p[0]->R[h];
        }

      //if ( facePerm % 2 )
      if ( vtkMath::Determinant3x3( e0, e1, e2) < 0 )
        {
        // Odd faces have reversed normals
        // Swap a pair of points on the triangle to get a positive arrangement of vertices
        PointRef tmpref = p[3];
        p[3] = p[2];
        p[2] = tmpref;
        }

      EleRef eref[6];
      // Do NOT reorder! These add the faces in the traditional order: 0-1-2, 0-3-1, 1-3-2, 2-3-0
      eref[2] = this->TriENext( this->TriSym( this->TriFindOrCreate( eleOff, p[0], p[1], p[2] ) ) ); // Edge 0-2
      eref[0] = this->TriENext( this->TriSym( this->TriFindOrCreate( eleOff, p[0], p[3], p[1] ) ) ); // Edge 0-1
      eref[1] = this->TriENext( this->TriSym( this->TriFindOrCreate( eleOff, p[1], p[3], p[2] ) ) ); // Edge 1-2
      eref[3] = this->TriSym( this->TriENext( this->TriFindOrCreate( eleOff, p[2], p[3], p[0] ) ) ); // Edge 0-3
      eref[4] = this->TriENext( eref[0] ); // Edge 1-3
      eref[5] = this->TriENext( eref[1] ); // Edge 2-3

      // Glue faces together
      this->TriIdempotentFMerge( this->TriSym( eref[0] ), this->TriENext2( eref[2] ) );
      this->TriIdempotentFMerge( this->TriSym( eref[1] ), this->TriENext ( eref[2] ) );
      this->TriIdempotentFMerge( this->TriSym( eref[2] ), this->TriENext2( eref[3] ) );
      this->TriIdempotentFMerge( this->TriSym( eref[3] ), this->TriENext ( eref[4] ) );
      this->TriIdempotentFMerge( this->TriSym( eref[5] ), this->TriENext ( eref[3] ) );
      this->TriIdempotentFMerge( this->TriSym( eref[4] ), this->TriENext ( eref[5] ) );

      // Mark edges/faces that need to be examined for critical points.
      // The three new edges/faces leading to the interior point must always be searched:
      this->TriMarkEdges( eref[0], 1 );
      this->TriMarkEdges( eref[2], 1 );
      this->TriMarkEdges( eref[3], 1 );
      this->TriMarkFace( eref[0], 1 );
      this->TriMarkFace( eref[2], 1 );
      this->TriMarkFace( eref[3], 1 );

      // The three edges/one face on the boundary triangulation *may* need to be searched:
      vtkIdType tnum = i + bdyTri->Offsets[conn[faceDof]].Start;
      EleRef btri = (facePerm % 2) ? VTK_TRI_MUX_FV(tnum,5) : VTK_TRI_MUX_FV(tnum,0);
      this->TriMarkEdges( eref[1], bdyTri->TriGetEdgeMark( btri ) );
      btri = bdyTri->TriENext( btri );
      this->TriMarkEdges( eref[5], bdyTri->TriGetEdgeMark( btri ) );
      btri = bdyTri->TriENext( btri );
      this->TriMarkEdges( eref[4], bdyTri->TriGetEdgeMark( btri ) );
      this->TriMarkFace( eref[1], bdyTri->TriGetFaceMark( btri ) );
      this->TriFaceIsOnHull( this->TriSym( eref[0] ), 0 );
      this->TriFaceIsOnHull( this->TriSym( eref[1] ), 0 );
      this->TriFaceIsOnHull( this->TriSym( eref[2] ), 0 );
      this->TriFaceIsOnHull( this->TriSym( eref[3] ), 0 );
      }
    }

  this->CheckTetOrientations( d );
}

void vtkShoeBoxPartition::ComputeEdgeCriticalPoints( vtkShoeCell* sc, 
                                                     vtkGenericAttributeCollection* kappa, 
                                                     EleRef eref,
                                                     double* p1, double* p2,
                                                     vtkstd::vector<vtkstd::pair<double[3],EleRef> >& edgeCPs )
{
  char vars[] = "rst";
  char bvar[] = "u";
  double range[] = {0., 1.};

  int na = kappa->GetNumberOfAttributes();
  for ( int k = 0; k < na; ++ k )
    {
    vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( k ) );
    if ( ! sa ) continue;
    
    vtkPolynomialSystem* ps = vtkPolynomialSystem::New();
    vtkPolynomialSystem* bps;
          
    int nc = sa->GetNumberOfComponents();
    for ( int c = 0; c < nc; ++ c )
      {
      double dir[3];
      for ( int i = 0; i < 3; ++ i ) dir[i] = p2[i] - p1[i];

      sc->GetSymbolicAttributeComponentGradient( sa, c, ps );
      bps = ps->RestrictToLine( vars, p1, dir );
      bps->SetVariableRange( bvar, range );

#ifdef VTK_DEBUG_SOLVE_EDGE_CP
      vtkstd::cout << " O( " << p1[0] << " " << p1[1] << " " << p1[2] << " ),";
      vtkstd::cout << " D( " << dir[0] << " " << dir[1] << " " << dir[2] << ")\n";
      bps->PrintSystem( vtkstd::cout );
#endif // VTK_DEBUG_SOLVE_EDGE_CP

      bps->SolveSystem();

      vtkstd::pair<double[3],EleRef> newCP;
      bool go = true;
      for ( vtkIdType n = 0; go && n < bps->GetNumberOfRoots(); ++ n )
        {
        double root;
        int mult;
        bps->GetRoot( n, &root, mult );

#ifdef VTK_DEBUG_SOLVE_EDGE_CP
        vtkstd::cout << " root " << n << ": " << root << " ( mult " << mult << " )\n";
#endif // VTK_DEBUG_SOLVE_EDGE_CP

        if ( mult > -1 )
          {
          // FIXME: optimize that so that all roots can be treated at once.
          go = false; // One root is enough for an edge
          for ( int i = 0; i < 3; ++ i ) newCP.first[i] = p1[i] + root * dir[i];
          newCP.second = eref;
          edgeCPs.push_back( newCP );
          }
        }
      }
    bps->Delete();
    ps->Delete();
    }
}

void vtkShoeBoxPartition::ComputeTriangulationEdgeCriticalPoints( vtkShoeCell* sc, 
                                                                  int face,
                                                                  vtkGenericAttributeCollection* kappa, 
                                                                  vtkIdType d,
                                                                  vtkstd::vector<vtkstd::pair<double[3],EleRef> >& bndEdgeCPs )
{
  DofLocator dl = this->Offsets[d];
  vtkIdType dlEnd = dl.Start + dl.NumberOfElements;
  for ( vtkIdType i = dl.Start; i < dlEnd; ++ i )
    {
#ifdef VTK_DEBUG_BND_EDGE_CP
    vtkstd::cout << "Triangle " << i << " ( "
                 << this->Triangles[i].EndPoints[0]->Label << " "
                 << this->Triangles[i].EndPoints[1]->Label << " "
                 << this->Triangles[i].EndPoints[2]->Label << " ):";
#endif // VTK_DEBUG_BND_EDGE_CP
    EleRef eref = VTK_TRI_MUX_FV(i,0);
    for ( int j = 0; j < 3; ++ j )
      {
      if ( this->TriGetEdgeMark( eref ) )
        {
#ifdef VTK_DEBUG_BND_EDGE_CP
        vtkstd::cout << " edge " << j << " (" << this->TriGetEdgeMark( eref ) << ") ";
#endif // VTK_DEBUG_BND_EDGE_CP

        double p1[3], p2[3];
        for ( int i = 0; i < 2; ++ i ) 
          {
          p1[i] = this->TriOrg( eref )->R[i + 2];
          p2[i] = this->TriDest( eref )->R[i + 2];
          }

        sc->GetMetaData()->EmbedFaceCoords( face, p1 ); // p1 is modified as a side effect
        sc->GetMetaData()->EmbedFaceCoords( face, p2 ); // p2 is modified as a side effect

        this->ComputeEdgeCriticalPoints( sc, kappa, eref, p1, p2, bndEdgeCPs);
        this->TriMarkEdges( eref, 0 );
        }
      eref = this->TriENext( eref );
      }

#ifdef VTK_DEBUG_BND_EDGE_CP
    vtkstd::cout << "\n";
#endif // VTK_DEBUG_BND_EDGE_CP

    }
}

void vtkShoeBoxPartition::ComputeTetrahedralizationEdgeCriticalPoints( vtkShoeCell* sc, 
                                                                       vtkGenericAttributeCollection* kappa, 
                                                                       vtkIdType d,
                                                                       vtkstd::vector<vtkstd::pair<double[3],EleRef> >& intEdgeCPs )
{
  DofLocator dl = this->Offsets[d];
  vtkIdType dlEnd = dl.Start + dl.NumberOfElements;
  for ( vtkIdType i = dl.Start; i < dlEnd; ++ i )
    {
#ifdef VTK_DEBUG_INT_EDGE_CP
    vtkstd::cout << "Triangle " << i << " ( "
                 << this->Triangles[i].EndPoints[0]->Label << " "
                 << this->Triangles[i].EndPoints[1]->Label << " "
                 << this->Triangles[i].EndPoints[2]->Label << " ):";
#endif // iVTK_DEBUG_INT_EDGE_CP
    EleRef eref = VTK_TRI_MUX_FV(i,0);
    for ( int j = 0; j < 3; ++ j )
      {
      if ( this->TriGetEdgeMark( eref ) )
        {

#ifdef VTK_DEBUG_INT_EDGE_CP
        vtkstd::cout << " edge " << j;
#endif // VTK_DEBUG_INT_EDGE_CP

        EleRef ring = eref;
        // Make sure this edge is not on the hull; if it is, this is wrong, as it should have already been treated
        // as an edge of a face triangulation.
        do
          {
          if ( this->TriIsFaceOnHull( ring ) ) 
            {
            
#ifdef VTK_DEBUG_INT_EDGE_CP
            vtkstd::cout << " is a boundary edge and thus should have already been treated.";
#endif // VTK_DEBUG_INT_EDGE_CP

            break;
            }
          ring = this->TriFNext( ring );
          } while ( ring != eref );

        double p1[3], p2[3];
        for ( int i = 0; i < 3; ++ i )
          {
            p1[i] = this->TriOrg( eref )->R[i];
            p2[i] = this->TriDest( eref )->R[i];
          }

        this->ComputeEdgeCriticalPoints( sc, kappa, eref, p1, p2, intEdgeCPs);
        this->TriMarkEdges( eref, 0 );
        }
      eref = this->TriENext( eref );
      }

#ifdef VTK_DEBUG_INT_EDGE_CP
    vtkstd::cout << "\n";
#endif // VTK_DEBUG_INT_EDGE_CP

    }
}

void vtkShoeBoxPartition::ComputeTetrahedralizationFaceCriticalPoints( vtkShoeCell* sc, 
                                                                       vtkGenericAttributeCollection* kappa, 
                                                                       vtkIdType d,
                                                                       vtkstd::vector<vtkstd::pair<double[3],EleRef> >& intFaceCPs )
{
  DofLocator dl = this->Offsets[d];
  vtkIdType dlEnd = dl.Start + dl.NumberOfElements;
  for ( vtkIdType i = dl.Start; i < dlEnd; ++ i )
    {
    PointRef p0 = this->Triangles[i].EndPoints[0];
    PointRef p1 = this->Triangles[i].EndPoints[1];
    PointRef p2 = this->Triangles[i].EndPoints[2];

#ifdef VTK_DEBUG_FACE_CP
    vtkstd::cout << "# Check face " << i << " ( "
                 << p0->Label << " "
                 << p1->Label << " "
                 << p2->Label << " ):\n";
#endif // VTK_DEBUG_FACE_CP
    EleRef eref = VTK_TRI_MUX_FV(i,0);
    if ( this->TriGetFaceMark( eref ) )
      {
      // FIXME: normally not needed as boundary faces should be compatible and marked accordingly
      bool boundary = false;
        
      int na = kappa->GetNumberOfAttributes();
      for ( int k = 0; k < na; ++ k )
        {
        vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( kappa->GetAttribute( k ) );
        if ( ! sa ) continue;
        
        vtkPolynomialSystem* ps = vtkPolynomialSystem::New();
        char vars[] = "rst";
        
        vtkPolynomialSystem* bps;
        char bvars[] = "uv";
        double range[4];
        range[0] = range[2] = 0.;
        range[1] = range[3] = 1.;
          
        int nc = sa->GetNumberOfComponents();
        for ( int l = 0; l < nc; ++ l )
          {
          sc->GetSymbolicAttributeComponentGradient( sa, l, ps );
          double ori[3];
          double dir1[3];
          double dir2[3];
          for ( int m = 0; m < 3; ++ m )
            {
            ori[m] = p0->R[m];
            dir1[m] = p1->R[m] - p0->R[m];
            dir2[m] = p2->R[m] - p0->R[m];
            }
          bps = ps->RestrictToPlane( vars, ori, dir1, dir2 );
          bps->SetVariableRange( bvars, range );
#ifdef VTK_DEBUG_FACE_CP
          vtkstd::cout << " O( " << ori[0] << " " << ori[1] << " " << ori[2] << " ),";
          vtkstd::cout << " D1( " << dir1[0] << " " << dir1[1] << " " << dir1[2] << "),";
          vtkstd::cout << " D2( " << dir2[0] << " " << dir2[1] << " " << dir2[2] << ")\n";
          bps->PrintSystem( vtkstd::cout );
#endif // VTK_DEBUG_FACE_CP
          bps->SolveSystem();
          
          vtkstd::pair<double[3],EleRef> newCP;
          bool go = true;
          for ( vtkIdType n = 0; go && n < bps->GetNumberOfRoots(); ++ n )
            {
            double root[2];
            int mult;
            bps->GetRoot( n, root, mult );
            // FIXME: allow for ( strict ) variable ranging in PolynomialSystem
            if ( ( root[0] > 0. ) && ( root[1] > 0. ) && ( root[0] + root[1] < 1. ) )
              {
#ifdef VTK_DEBUG_FACE_CP
              vtkstd::cout << "( " << root[0]
                           << ", " << root[1]
                           << " )\n";
#endif // VTK_DEBUG_FACE_CP
              go = false; // One root is enough for an edge
              for ( int m = 0; m < 3; ++ m ) newCP.first[m] = ori[m] + root[0] * dir1[m] + root[1] * dir2[m];
              newCP.second = eref;
              
                if ( boundary )
                  {
                  vtkErrorMacro(" Boundary face (" << ori[0] << " " << ori[1] << " " << ori[2] << " ) - (" 
                                << dir1[0] << " " << dir1[1] << " " << dir1[2] << ") - ("
                                << dir2[0] << " " << dir2[1] << " " << dir2[2] << ") has a critical point\n" );  
                  }
                else intFaceCPs.push_back( newCP );
              }
            }
          }
        bps->Delete();
        ps->Delete();
        }
      }
    }
}

bool vtkShoeBoxPartition::CorrectEdgeTopology( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, 
					       vtkIdType d, int face )
{
#ifdef VTK_VERBOSE_CORRECT_EDGE
  vtkstd::cout << "## CorrectEdgeTopology :\n";
#endif // VTK_VERBOSE_CORRECT_EDGE

  vtkstd::vector<vtkstd::pair<double[3],EleRef> > edgeCPs;
  for ( int i = 0; i < this->MaximumSplitEdgeDepth; ++ i )
    {
    if ( face < 0 ) // this occurs iff the DOF is a volume DOF
      {
      this->ComputeTetrahedralizationEdgeCriticalPoints( sc, kappa, d, edgeCPs );
      }
    else // this occurs iff the DOF is a face DOF
      {
      this->ComputeTriangulationEdgeCriticalPoints( sc, face, kappa, d, edgeCPs );
      }
    
#ifdef VTK_VERBOSE_CORRECT_EDGE
    vtkstd::cout << "# Edge depth level " << i << ":\n";
    vtkstd::cout << edgeCPs.size() << " edge critical points found.\n";
#endif // VTK_VERBOSE_CORRECT_EDGE
      
    if ( ! edgeCPs.size() )
      {
    
#ifdef VTK_VERBOSE_CORRECT_EDGE
      vtkstd::cout << "* Edge topology is kappa-compatible.\n";
#endif // VTK_VERBOSE_CORRECT_EDGE
      
      return false;
      }

    for ( vtkstd::vector<vtkstd::pair<double[3],EleRef> >::iterator eit = edgeCPs.begin(); eit != edgeCPs.end(); ++ eit )
      {

#ifdef VTK_VERBOSE_CORRECT_EDGE
      if ( face < -1 ) vtkstd::cout << " Interior edge ";
      else vtkstd::cout << " Boundary edge ";
      vtkstd::cout << this->TriOrg( eit->second )->Label << "-"
                   << this->TriDest( eit->second )->Label << " must be split.\n";
#endif // VTK_VERBOSE_CORRECT_EDGE

      this->SplitEdge( sc, d, eit->second, eit->first /* rst */, face );
      }
      
    edgeCPs.clear();
    }
  
  return true;
}

bool vtkShoeBoxPartition::CorrectFaceTopology( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, vtkIdType d )
{
#ifdef VTK_VERBOSE_CORRECT_FACE
  vtkstd::cout << "## CorrectFaceTopology :\n";
#endif // VTK_VERBOSE_CORRECT_FACE

  vtkstd::vector<vtkstd::pair<double[3],EleRef> > faceCPs;
  for ( int i = 0; i < this->MaximumSplitFaceDepth; ++ i )
    {
    this->ComputeTetrahedralizationFaceCriticalPoints( sc, kappa, d, faceCPs );
    
#ifdef VTK_VERBOSE_CORRECT_FACE
    vtkstd::cout << "# Face depth level " << i << ":\n";
    vtkstd::cout << faceCPs.size() << " face critical points found.\n";
#endif // VTK_VERBOSE_CORRECT_FACE
      
    if ( ! faceCPs.size() )
      {
    
#ifdef VTK_VERBOSE_CORRECT_EDGE
      vtkstd::cout << "* Face topology is kappa-compatible.\n";
#endif // VTK_VERBOSE_CORRECT_EDGE
      
      return false;
      }

    for ( vtkstd::vector<vtkstd::pair<double[3],EleRef> >::iterator eit = faceCPs.begin(); eit != faceCPs.end(); ++ eit )
      {

#ifdef VTK_VERBOSE_CORRECT_FACE
      vtkstd::cout << " Interior face "
                   << this->TriOrg( eit->second )->Label << "-"
                   << this->TriDest( eit->second )->Label << "-"
                   << this->TriDest( this->TriENext( eit->second ) )->Label << " must be split.\n";
#endif // VTK_VERBOSE_CORRECT_FACE

      //this->SplitFace( sc, d, eit->second, eit->first /* rst */ );
      }

    faceCPs.clear();
    }
  
  vtkstd::cout << "!! Face topology is NOT kappa-compatible.\n";
  return false;
}

bool vtkShoeBoxPartition::CorrectTriTopology( vtkShoeCell* sc, vtkIdType* conn, int ne, 
                                              int dof, vtkGenericAttributeCollection* kappa )
{
  vtkIdType d = conn[dof];

#ifdef VTK_VERBOSE_CORRECT_TRI
  vtkstd::cout << "\n### CorrectTriTopology (dof " << dof
               << ", " << this->Offsets[d].NumberOfElements << " interior triangles):\n";
#endif // VTK_VERBOSE_CORRECT_TRI

  if ( this->CorrectEdgeTopology( sc, kappa, d, dof - ne ) )
    {
    vtkWarningMacro( "Triangle edge topology not kappa-compatible at maximal allowed depth (" << this->MaximumSplitEdgeDepth << ")." );
    return true;
    }

#ifdef VTK_VERBOSE_CORRECT_TRI
  vtkstd::cout << "* Triangle edge topology is kappa-compatible.\n";
#endif // VTK_VERBOSE_CORRECT_TRI
  
  return false;
}

bool vtkShoeBoxPartition::CorrectTetTopology( vtkShoeCell* sc, vtkIdType* conn, int ne, 
                                              int nf, vtkGenericAttributeCollection* kappa )
{
  int dof = ne + nf;
  vtkIdType d = conn[dof];

#ifdef VTK_VERBOSE_CORRECT_TET
  vtkstd::cout << "\n### CorrectTetTopology (dof " << dof 
               << ", " << this->Offsets[d].NumberOfElements << " interior triangles):\n";
#endif // VTK_VERBOSE_CORRECT_TET

  // First, take care of edges
  if ( this->CorrectEdgeTopology( sc, kappa, d, -1 ) ) // -1 < 0 to indicate that this is a body DOF
    {
    vtkWarningMacro( "Tetrahedral edge topology not kappa-compatible at maximal allowed depth (" << this->MaximumSplitEdgeDepth << ")." );
    return true;
    }

#ifdef VTK_VERBOSE_CORRECT_TET
  vtkstd::cout << "* Tetrahedral edge topology is kappa-compatible.\n";
#endif // VTK_VERBOSE_CORRECT_TET

  // Last, take care of faces
  if ( this->CorrectFaceTopology( sc, kappa, d ) ) 
    {
    vtkWarningMacro( "Tetrahedral face topology not kappa-compatible at maximal allowed depth (" << this->MaximumSplitFaceDepth << ")." );
    return true;
    }
  
#ifdef VTK_VERBOSE_CORRECT_TET
  vtkstd::cout << "** Tetrahedralization is kappa-compatible!\n\n";
#endif // VTK_VERBOSE_CORRECT_TET

  this->Dump();

  return false;
}


