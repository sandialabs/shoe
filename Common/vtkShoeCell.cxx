// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeCell.h"

#include <assert.h>
#include <vtkstd/map>

#include "vtkDataRecords.h"
#include "vtkShoeBox.h"
#include "vtkShoeBoxP.h"
#include "vtkShoeAttribute.h"
#include "vtkPolynomialSystem.h"
#include "vtkShoeCellMetaData.h"
#include "vtkShoeCellSpecies.h"
#include "vtkShoeCellIterator.h"
#include "vtkShoePointIterator.h"

#include "vtkObjectFactory.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkBitArray.h"
#include "vtkGenericAttributeCollection.h"

#include "Elements/ShoeHexahedron.h"

#define VTK_PRINT_POLSYS

vtkCxxRevisionMacro(vtkShoeCell,"$Revision: 10017 $");
vtkStandardNewMacro(vtkShoeCell);

vtkShoeCell::vtkShoeCell()
{
  this->DataSet = 0;
  this->CellID = -1;
  this->Meta = vtkShoeCellMetaData::ByShape[ shoe::NumberOfCellShapes ]; // "empty" metadata
  this->AttributeCache = 0;
  this->AttributeCacheLength = 0;
  this->IsAttributeDirty = 0;
}

vtkShoeCell::~vtkShoeCell()
{
  if ( this->AttributeCache )
    {
    vtkDoubleArray** arr = this->AttributeCache ;
    while ( this->AttributeCacheLength > 0 )
      {
      if  ( arr[ --this->AttributeCacheLength ] )
        {
        arr[ this->AttributeCacheLength ]->UnRegister( this );
        }
      }
    delete [] this->AttributeCache;
    }
  if ( this->IsAttributeDirty )
    {
    delete [] this->IsAttributeDirty;
    }
  if ( this->DataSet )
    {
    this->DataSet->UnRegister( this );
    }
}

void vtkShoeCell::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Meta: " << this->Meta << vtkstd::endl;
  os << indent << "DataSet: " << this->DataSet << vtkstd::endl;
  os << indent << "AttributeCache: " << this->AttributeCache
    << " (" << this->AttributeCacheLength << ")" << vtkstd::endl;
  os << indent << "IsAttributeDirty: " << this->IsAttributeDirty
    << " (" << this->AttributeCacheLength << ")" << vtkstd::endl;
}

void vtkShoeCell::UpdateCache( int at )
{
  if ( ! this->DataSet )
    {
    return;
    }
  //assert( at < this->DataSet->GetAttributes()->GetNumberOfAttributes() );
  if ( this->AttributeCacheLength < this->DataSet->GetAttributes()->GetNumberOfAttributes() )
    {
    int N;
    if ( this->AttributeCacheLength )
      {
      vtkDoubleArray** aCache = this->AttributeCache;
      int* aDirty = this->IsAttributeDirty;
      N = this->DataSet->GetAttributes()->GetNumberOfAttributes();
      this->AttributeCache = new vtkDoubleArray* [ N ];
      this->IsAttributeDirty = new int [ N ];
      int i;
      for ( i=0; i<this->AttributeCacheLength; ++i )
        {
        this->AttributeCache[i] = aCache[i];
        this->IsAttributeDirty[i] = aDirty[i];
        }
      for ( i=this->AttributeCacheLength; i<N; ++i )
        {
        this->IsAttributeDirty[i] = 1;
        }
      if ( aCache )
        {
        delete [] aCache;
        }
      if ( aDirty )
        {
        delete [] aDirty;
        }
      }
    else
      { // no attributes
      N = this->DataSet->GetAttributes()->GetNumberOfAttributes();
      this->AttributeCache = new vtkDoubleArray* [ N ];
      this->IsAttributeDirty = new int [ N ];
      int i;
      for ( i=0; i<N; ++i )
        {
        this->AttributeCache[i] = 0;
        this->IsAttributeDirty[i] = 1;
        }
      }
      this->AttributeCacheLength = N;
    }
  if ( this->IsAttributeDirty[at] )
    {
    this->GetPermutedAttributeValues( at ); 
    this->IsAttributeDirty[at] = false;
    }
}

void vtkShoeCell::GetPermutedAttributeValues( int at )
{
  vtkShoeAttribute* aa = vtkShoeAttribute::SafeDownCast( this->DataSet->GetAttributes()->GetAttribute( at ) );
  vtkDoubleArray* ac = this->AttributeCache[at];
  if ( ! ac )
    {
    ac = vtkDoubleArray::New();
    ac->SetNumberOfComponents( aa->GetNumberOfComponents() );
    this->AttributeCache[at] = ac;
    }
  int sz = this->GetNumberOfPoints();
  vtkIdType* pconn = this->DataSet->Connectivity->GetRecordAsIdType( this->CellRecord.Offset );
  vtkIdType* dconn = pconn + sz;
  vtkDataArray* pts = aa->GetPointData();
  vtkDataRecords* dof = aa->GetDOFData();
  for ( int i=0; i<this->GetNumberOfDOFNodes(); ++i )
    {
    sz += dof->GetNumberOfTuples( dconn[i] );
    }
  ac->SetNumberOfTuples( sz );
  for ( int i=0; i<this->GetNumberOfPoints(); ++i )
    {
    ac->SetTuple( i, pts->GetTuple( pconn[i] ) );
    }
  vtkShoeOrderTuple ord;
  vtkPolyInterpolant interp;
  vtkPolyProductSpace pspace;
  aa->GetCellTypeInfo( this->CellRecord.Species->GetId(), ord, interp, pspace );

  int s = this->Meta->NumberOfBoundaries[0];
  int nc = ac->GetNumberOfComponents();
  int ne = this->Meta->NumberOfBoundaries[1];
  for ( int i=0; i<ne; ++i )
    {
    vtkIdType nt = dof->GetNumberOfTuples( i );
    // FIMXE: make this generic
    ShoeHexahedron::LagrangeTensorPermuteEdge( ac->GetPointer( s * nc ), nc, dof, dconn[i], ord.Order[i < 8 ? i % 2 : 2], this->CellRecord.GetEdgePermutation( i ) );
    s += nt;
    }
  for ( int i=0; i<this->Meta->NumberOfBoundaries[2]; ++i, ++ne )
    {
    vtkIdType nt = dof->GetNumberOfTuples( ne );
    int fdir[2];
    if ( i < 2 )
      {
      fdir[0] = 1;
      fdir[1] = 2;
      }
    else
      {
      fdir[0] = 0;
      fdir[1] = ( i > 3 ? 1 : 2 ); 
      }
    // FIMXE: make this generic
    ShoeHexahedron::LagrangeTensorPermuteFace( ac->GetPointer( s * nc ), nc, dof, dconn[ne], ord.Order, fdir, this->CellRecord.GetFacePermutation( i ) );
    s += nt;
    }
  if ( this->Meta->Dimension == 3 )
    {
    dof->GetRecord( dconn[this->Meta->NumberOfDOFNodes-1], ac->GetPointer( s * nc ) );
    }
}

int vtkShoeCell::GetType()
{
  return this->Meta->Type;
}

int vtkShoeCell::GetDimension()
{
  return this->Meta->Dimension;
}

int vtkShoeCell::GetGeometryOrder()
{
  vtkShoeOrderTuple ot;
  this->DataSet->Geometry->GetOrder( this->CellRecord.Species->GetId(), ot );
  return ot.Order[0] + ot.Order[1] + ot.Order[2];
}

int vtkShoeCell::GetAttributeOrder( vtkGenericAttribute* a )
{
  vtkShoeOrderTuple ot;
  vtkShoeAttribute* sa = static_cast<vtkShoeAttribute*>(a);
  sa->GetOrder( this->CellRecord.Species->GetId(), ot );
  return ot.Order[0] + ot.Order[1] + ot.Order[2];
}

//int vtkShoeCell::GetHighestOrderAttribute(vtkGenericAttributeCollection *ac)

int vtkShoeCell::GetNumberOfPoints()
{
  return this->Meta->NumberOfBoundaries[0];
}

int vtkShoeCell::GetNumberOfBoundaries( int dim )
{
  assert( dim >= 0 && dim < 3 );
  return this->Meta->NumberOfBoundaries[dim];
}

int vtkShoeCell::GetNumberOfDOFNodes()
{
  return this->Meta->NumberOfDOFNodes;
}

void vtkShoeCell::GetSymbolicAttributeComponentGradient( vtkGenericAttribute* a, int c, vtkPolynomialSystem* s )
{
  vtkShoeOrderTuple ot;
  vtkShoeAttribute* sa = static_cast<vtkShoeAttribute*>(a);
  sa->GetOrder( this->CellRecord.Species->GetId(), ot );

  int at = this->DataSet->AttributeId( sa );
  this->UpdateCache( at );

  vtkDoubleArray* ac = this->AttributeCache[at];
  int nc = ac->GetNumberOfComponents();
  if ( nc == 1 )
    {
    this->Meta->SymbolicFieldGradient( ot.Order, ac->GetPointer(0), s );
    }
  else
    {
    vtkIdType nt = ac->GetNumberOfTuples();
    double* componentCache = new double[ nt ];
    double* src = ac->GetPointer(0) + c;
    for ( vtkIdType i = 0; i < nt; ++i, src += nc )
      {
      componentCache[i] = *src;
      }
    this->Meta->SymbolicFieldGradient( ot.Order, componentCache, s );
    }
}

void vtkShoeCell::Derivatives( int subId, double pcoords[3], vtkGenericAttribute* attribute, double* derivs )
{
  int at = this->DataSet->AttributeId( attribute ); // FIXME!!!
  this->UpdateCache( at );
  this->Meta->Derivatives( pcoords, this->AttributeCache[ at ], derivs );
}

void vtkShoeCell::ParametricDerivatives( double pcoords[3], vtkShoeAttribute* attribute, double* derivs )
{
  // FIXME: Horribly inefficient hack. Uses forward differences and a fixed step size.
  // These derivatives are available analytically!
  double rPlusDr[3];
  vtkstd::vector<double> base;
  vtkstd::vector<double> fdr;
  vtkstd::vector<double> fds;
  vtkstd::vector<double> fdt;
  int nc = attribute->GetNumberOfComponents();
  base.reserve( nc );
  fdr.reserve( nc );
  fds.reserve( nc );
  fdt.reserve( nc );
  this->InterpolateTuple( attribute, pcoords, &base[0] );
  int i;
  for ( i = 0; i < 3; ++i ) rPlusDr[i] = pcoords[i];
  rPlusDr[0] += 1.e-5;
  this->InterpolateTuple( attribute, rPlusDr, &fdr[0] );
  rPlusDr[0] = pcoords[0];
  rPlusDr[1] += 1.e-5; 
  this->InterpolateTuple( attribute, rPlusDr, &fds[0] );
  rPlusDr[1] = pcoords[1];
  rPlusDr[2] += 1.e-5; 
  this->InterpolateTuple( attribute, rPlusDr, &fdt[0] );
  for ( int c = 0; c < nc; ++c )
    {
    *(derivs++) = (fdr[c] - base[c])/1.e-5;
    *(derivs++) = (fds[c] - base[c])/1.e-5;
    *(derivs++) = (fdt[c] - base[c])/1.e-5;
    }
}

int vtkShoeCell::GetParametricCenter( double pcoords[3] )
{
  for ( int i=0; i<3; ++i )
    pcoords[i] = this->Meta->ParametricCenter[i];
  return 0;
}

double vtkShoeCell::GetParametricDistance( double pcoords[3] )
{
  return this->Meta->GetParametricDistance( pcoords );
}

double* vtkShoeCell::GetParametricCoords()
{
  return this->Meta->ParametricCorners;
}

int vtkShoeCell::GetNumberOfVerticesOnFace( int faceId )
{
  return this->Meta->NumberOfVerticesOnFace[ faceId ];
}

int* vtkShoeCell::GetFaceArray( int faceId )
{
  return this->Meta->FaceArray[ faceId ];
}

int* vtkShoeCell::GetEdgeArray( int edgeId )
{
  return this->Meta->EdgeArray[ edgeId ];
}

int vtkShoeCell::IsInDataSet()
{
  return this->DataSet->IsTemporaryMesh();
}

void vtkShoeCell::GetPointIterator( vtkGenericPointIterator* it )
{
  vtkShoePointIterator* sit = vtkShoePointIterator::SafeDownCast( it );
  sit->SetDataSet( 0 );
  sit->SetCell( this );
}

vtkGenericCellIterator* vtkShoeCell::NewCellIterator()
{
  return this->DataSet->NewCellIterator();
}

void vtkShoeCell::GetBoundaryIterator( vtkGenericCellIterator* boundaries, int dim )
{
  assert( dim < 3 && dim >= -1 );
  vtkShoeCellIterator* bdyIt = vtkShoeCellIterator::SafeDownCast( boundaries );
  assert( bdyIt );

  bdyIt->SetDataSet( this->DataSet );
  bdyIt->Dimension = dim;
  bdyIt->SetParentCell( this );
  bdyIt->CurrentFacet = -1; // must call Begin to init.
}

int vtkShoeCell::FindClosestBoundary( int vtkNotUsed(subId), double pcoords[3], vtkGenericCellIterator*& boundary )
{
  return this->Meta->FindClosestBoundary( pcoords, (vtkShoeCellIterator*&)boundary );
}


int vtkShoeCell::CountNeighbors( vtkGenericAdaptorCell *boundaryRaw )
{
  vtkShoeCell* boundary = (vtkShoeCell*) boundaryRaw;
  // Last DOF node of boundary cell will have a list of all the neighbors we need.
  // Since reverse lookup table includes boundary cell itself, subtract one from count.
  vtkShoeCell* bdy = static_cast<vtkShoeCell*>(boundary);
  vtkIdType* bdyDof;
  boundary->DataSet->Connectivity->GetRecord( bdy->CellID, bdyDof );
  return this->DataSet->DOFLinks->GetNumberOfTuples( bdyDof[ boundary->Meta->NumberOfDOFNodes + boundary->Meta->NumberOfBoundaries[0] - 1 ] ) - 1;
}

void vtkShoeCell::CountEdgeNeighbors( int* sharing )
{
  int sz = this->GetNumberOfPoints();
  vtkIdType* dconn = this->DataSet->Connectivity->GetRecordAsIdType( this->CellRecord.Offset ) + sz;
  for ( int i=0; i<this->Meta->NumberOfBoundaries[1]; ++i )
    {
    sharing[i] = this->DataSet->GetDOFLinks()->GetNumberOfTuples( dconn[i] ) - 1;
    }
}

void vtkShoeCell::GetNeighbors( vtkGenericAdaptorCell *boundary, vtkGenericCellIterator *neighbors )
{
  // See CountNeighbors() for a hint on how to implement.
}

int vtkShoeCell::EvaluatePosition( double x[3], double *closestPoint, int &subId, double pcoords[3], double &dist2 )
{
  subId = 0;
  // err. *Is* x in this cell?
  return 1;
}

void vtkShoeCell::EvaluateLocation( int vtkNotUsed(subId), double pcoords[3], double x[3] )
{
  this->InterpolateTuple( this->DataSet->GetGeometry(), pcoords, x );
}

void vtkShoeCell::InterpolateTuple( vtkGenericAttribute* a, double pcoords[3], double* val )
{
  int at = this->DataSet->AttributeId( a );
  this->UpdateCache( at );
  vtkShoeOrderTuple ord;
  ((vtkShoeAttribute*)a)->GetOrder( this->CellRecord.Species->GetId(), ord );
  this->Meta->ShapeFunctions( ord.Order, pcoords, ShapeTmp );
  vtkDoubleArray* attrib = this->AttributeCache[at];
  double* cellVals = attrib->GetPointer(0);
  int c,s;
  for ( c=0; c<attrib->GetNumberOfComponents(); ++c )
    {
    val[c] = 0.;
    }
  for ( s=0; s<attrib->GetNumberOfTuples(); ++s )
    {
    //cout << ShapeTmp[s] << " * (";
    for ( c=0; c<attrib->GetNumberOfComponents(); ++c )
      {
      //cout << " " << *cellVals;
      val[c] += ShapeTmp[s]*(*cellVals);
      cellVals++;
      }
    //cout << " )\n";
    }
#if 0
  cout << "( " << pcoords[0] << ", " << pcoords[1] << ", " << pcoords[2] << " ): (";
  for ( c=0; c<attrib->GetNumberOfComponents(); ++c )
    {
    cout << " " << val[c];
    }
  cout << " )\n";
#endif // 0
 }

void vtkShoeCell::InterpolateTuple( vtkGenericAttributeCollection* c, double pcoords[3], double* val )
{
  vtkGenericAttribute* a;
  int na = c->GetNumberOfAttributes();
  for ( int i=0; i<na; ++i )
    {
    a = c->GetAttribute(i);
    this->InterpolateTuple( a, pcoords, val );
    val += a->GetNumberOfComponents();
    }
}

/*
void vtkShoeCell::Contour(
  vtkContourValues *values, vtkImplicitFunction *f, vtkGenericAttributeCollection *attributes,
  vtkGenericCellTessellator *tess, vtkPointLocator *locator, vtkCellArray *verts,
  vtkCellArray *lines, vtkCellArray *polys, vtkPointData *outPd, vtkCellData *outCd,
  vtkPointData *internalPd, vtkPointData *secondaryPd, vtkCellData *secondaryCd )
{
}

void vtkShoeCell::Clip(
  double value, vtkImplicitFunction *f, vtkGenericAttributeCollection *attributes,
  vtkGenericCellTessellator *tess, int insideOut, vtkPointLocator *locator, 
  vtkCellArray *connectivity, vtkPointData *outPd, vtkCellData *outCd, vtkPointData *internalPd,
  vtkPointData *secondaryPd, vtkCellData *secondaryCd )
{
}
*/

int vtkShoeCell::IntersectWithLine( double p1[3], double p2[3], double tol, double &t, double x[3], double pcoords[3], int &subId )
{
  vtkWarningMacro( "No IntersectWithLine" );
  return 1;
}

void vtkShoeCell::GetBounds( double bounds[6] )
{
  vtkWarningMacro( "No Cell GetBounds" );
}

/*
void vtkShoeCell::Tessellate( vtkGenericAttributeCollection *attributes, vtkGenericCellTessellator *tess,
                           vtkPoints *points, vtkPointLocator *locator, vtkCellArray* cellArray,
                           vtkPointData *internalPd, vtkPointData *pd, vtkCellData* cd, vtkUnsignedCharArray *types)
{
}
*/

int vtkShoeCell::IsFaceOnBoundary( vtkIdType faceId )
{
  if ( this->Meta->Dimension != 3 )
    { // volumetric cells cannot be on the boundary
    return 0;
    }

  vtkIdType* conn = this->DataSet->Connectivity->GetRecordAsIdType( this->CellID );
  conn += this->Meta->NumberOfBoundaries[0] + this->Meta->NumberOfDOFNodes - this->Meta->NumberOfBoundaries[2] + faceId - 1;
  return this->DataSet->GetDOFLinks()->GetNumberOfTuples( *conn ) == 1;
}

int vtkShoeCell::IsOnBoundary()
{
  if ( this->Meta->Dimension != 2 )
    {
    // our trick only works for faces.
    // How do we tell if an edge is on the boundary?
    // I guess either it's the only cell using its body-centered DOF
    // OR it is referenced by at least one face on the boundary
    // OR it is referenced by at least one volume with one face on the boundary...
    return 0;
    }

  vtkIdType* conn = this->DataSet->GetConnectivity()->GetRecordAsIdType( this->CellID );
  conn += this->Meta->NumberOfBoundaries[0] + this->Meta->NumberOfDOFNodes - 1;
  return this->DataSet->GetDOFLinks()->GetNumberOfTuples( *conn ) == 1;
}

void vtkShoeCell::GetPointIds( vtkIdType* id )
{
  vtkIdType* conn = this->DataSet->Connectivity->GetRecordAsIdType( this->CellID );
  for ( int i=0; i<this->Meta->NumberOfBoundaries[0]; ++i )
    {
    *id = *conn;
    ++id;
    ++conn;
    }
}

vtkIdType* vtkShoeCell::GetConnectivity()
{
  return this->DataSet->Connectivity->GetRecordAsIdType( this->CellRecord.Offset );
}

/*
void vtkShoeCell::TriangulateFace(
  vtkGenericAttributeCollection *attributes, vtkGenericCellTessellator *tess, int index, 
  vtkPoints *points, vtkPointLocator *locator, vtkCellArray *cellArray,
  vtkPointData *internalPd, vtkPointData *pd, vtkCellData *cd )
{
}
*/

void vtkShoeCell::SetDataSet( vtkShoeBox* d )
{
  if ( d == this->DataSet )
    {
    return;
    }
  if ( this->DataSet )
    {
    this->DataSet->UnRegister( this );
    }
  this->DataSet = d;
  if ( this->DataSet )
    {
    this->DataSet->Register( this );
    }
  this->Modified();
}

vtkPolynomialSystem* vtkShoeCell::RestrictToBoundary( vtkPolynomialSystem* ps, int bdy )
{
  vtkPolynomialSystem* bps;
  if ( bdy < 0 || bdy > this->Meta->NumberOfDOFNodes - 2 )
    {
    bps = 0;
    }
  else if ( bdy == this->Meta->NumberOfDOFNodes - 1 )
    { // Restricting a system to itself is itself
    bps = ps;
    }
  else if ( bdy < this->Meta->NumberOfBoundaries[1] )
    { // Get an edge from the volume or face
    double* bp = &(this->Meta->BoundaryEdgeTransforms[bdy][0]);
    double* dir = &(this->Meta->BoundaryEdgeTransforms[bdy][3]);
    char vars[] = "rst";
    bps = ps->RestrictToLine( vars, bp, dir );
    }
  else
    { // Must be a face
    // turn bdy from a connectivity index into a face index
    bdy -= this->Meta->NumberOfBoundaries[1];
    double* bp = &(this->Meta->BoundaryFaceTransforms[bdy][0]);
    double* dir1 = &(this->Meta->BoundaryFaceTransforms[bdy][3]);
    double* dir2 = &(this->Meta->BoundaryFaceTransforms[bdy][6]);
    // FIXME: In a perfect world, we would check that dir1 & dir2 are not
    // collinear or of 0 length -- if they were, we would move to a corner
    // node where these conditions hold.
    char vars[] = "rst";
    bps = ps->RestrictToPlane( vars, bp, dir1, dir2 );
    }
  return bps;
}

void vtkShoeCell::StoreRelevantRoots( vtkPolynomialSystem* ps, vtkIdType* conn, int connOffset, int perm, vtkDataRecords* cpStorage, vtkShoeAttribute* sa )
{
  int NumberOfComponents = cpStorage->GetNumberOfComponents();
  int NumberOfRoots = ps->GetNumberOfRoots();
  int NumberOfRealRoots; // real root counter
  int NumberOfValidRoots = 0; // valid root counter

  vtkIdType dof = conn[ this->Meta->NumberOfBoundaries[0] + connOffset ];
  char varOrder[3];
  double tuple[NumberOfComponents];
  // FIXME. We currently rely on alphabetical ordering rather than the contents of varOrder
  // as returned by GetVariableOrder:
  int dim = ps->GetVariableOrder( varOrder );

  for ( NumberOfRealRoots = 0; NumberOfRealRoots < NumberOfRoots; ++ NumberOfRealRoots )
    {
    ps->GetRoot( NumberOfRealRoots, tuple );
    if ( this->Meta->IsParameterInDomain( tuple, connOffset ) >= 0 )
      {
      if ( dim == 3 )
        {
        this->InterpolateTuple( sa, tuple, tuple + 3 );
        }
      else if ( dim == 2 )
        {
        double evaltuple[3];
        for ( int i = 0; i < 2; ++ i )
          {
          evaltuple[i] = tuple[i];
          }
        this->Meta->EmbedFaceCoords( connOffset - this->Meta->NumberOfBoundaries[1], evaltuple );
        this->InterpolateTuple( sa, evaltuple, tuple + 3 );
        
        vtkShoeCellRecord::ParametricFaceCoordsToStorageOrder( 1, tuple, perm );
        tuple[2] = 0.; // Pad unused coordinates
        }
      else if ( dim == 1 )
        {
        double evaltuple[3];
        evaltuple[0] = tuple[0];
        this->Meta->EmbedEdgeCoords( connOffset, evaltuple );
        this->InterpolateTuple( sa, evaltuple, tuple + 3 );

        vtkShoeCellRecord::ParametricEdgeCoordsToStorageOrder( 1, tuple, perm );
        tuple[1] = tuple[2] = 0.; // Pad unused coordinates
        }
       cpStorage->InsertNextTuple( dof, tuple );
       ++ NumberOfValidRoots;
      }
    }
}

void vtkShoeCell::ComputeDOFCriticalPoints( vtkShoeAttribute* att )
{
  // If the critical points are marked as valid (as a whole), don't
  // check individual bits in att->GetCriticalPointsComputed() array.
  if ( ! att->GetCriticalPointsDirty() )
    return;

  // Start at node with highest dimension (volume node, if this is a volumetric cell)
  // There is always exactly one DOF node of this dimension.
  // This node is immediately preceded by GetNumberOfDOFNodes(Dimension-1)...
  int connOffset = this->GetNumberOfDOFNodes() - 1;
  int npts = this->Meta->NumberOfBoundaries[0];
  vtkIdType* conn = this->DataSet->GetConnectivity()->GetRecordAsIdType( this->CellRecord.Offset );
  int perm;

  // First things first... see if the critical points have already been marked as computed
  // We need only check the highest-dimension node, since if it is done, all its boundaries are
  // done as well. Then, we can check each boundary as we go along.
  vtkDataRecords** cp = att->GetCriticalPoints();
  vtkBitArray* cpc = att->GetCriticalPointsComputed();
  if ( cpc->GetValue( conn[npts + connOffset] ) )
    {
    return;
    }

  vtkPolynomialSystem* ps = vtkPolynomialSystem::New();
  vtkPolynomialSystem* bps;
  int nc = att->GetNumberOfComponents();
  for ( int i = 0; i < nc; ++i )
    {
    this->GetSymbolicAttributeComponentGradient( att, i, ps );
#ifdef VTK_PRINT_POLSYS
    ps->PrintSystem( vtkstd::cout );
#endif // VTK_PRINT_POLSYS
    ps->SolveSystem();
    switch ( this->Meta->Dimension )
      {
      default:
      case 0:
        break; // do nothing
      case 1:
        perm = this->CellRecord.GetEdgePermutation( 0 );
        this->StoreRelevantRoots( ps, conn, connOffset, perm, cp[i], att );
        break;
      case 2:
          {
          perm = this->CellRecord.GetFacePermutation( 0 );
          this->StoreRelevantRoots( ps, conn, connOffset, perm, cp[i], att );
          int ne = this->GetNumberOfBoundaries(1);
          connOffset -= ne; // move to the first edge
          for ( int e=0; e<ne; ++e )
            {
            if ( cpc->GetValue( conn[connOffset+e] ) )
              continue; // we've already computed these critical points

            bps = this->RestrictToBoundary( ps, connOffset + e );
#ifdef VTK_PRINT_POLSYS
            bps->PrintSystem( vtkstd::cout );
#endif // VTK_PRINT_POLSYS
            bps->SolveSystem();
            perm = this->CellRecord.GetEdgePermutation( e );
            this->StoreRelevantRoots( bps, conn, connOffset + e, perm, cp[i], att );
            bps->Delete();
            if ( i == nc - 1 )
              {
              cpc->SetValue( conn[connOffset+e], 1 );
              }
            }
          }
        break;
      case 3:
          {
          // FIXME. The stuff below assumes there's a way to get from a face to an edge.
          // While this is easy for humans to do, it's not so simple to program -- at
          // least not in the general case. It is easy but tedious if done on a per-element
          // basis. But I don't want to do that.
#if 0
          this->StoreRelevantRoots( ps, conn, connOffset, 0, cp[i] );
          // Critical points on faces
          vtkPolynomialSystem** psFace = new vtkPolynomialSystem*[ this->GetNumberOfBoundaries(2) ];
          int nf = this->GetNumberOfBoundaries(1);
          connOffset -= nf; // move to the first face
          for ( int f=0; f<nf; ++f )
            {
            if ( cpc->GetValue( conn[connOffset+f] ) )
              continue; // we've already computed these critical points

            vtkPolynomialSystem* bps = this->RestrictToBoundary( ps, connOffset + f );
            bps->SolveSystem();
            perm = this->CellRecord.GetEdgePermutation( f );
            this->StoreRelevantRoots( bps, conn, connOffset + f, perm, cp[i] );
            if ( i == nc - 1 )
              {
              cpc->SetValue( conn[connOffset+f], 1 );
              }
            }

          // Critical points on edges
          int ne = this->GetNumberOfBoundaries(1);
          connOffset -= ne; // move to the first edge
          for ( int e=0; e<ne; ++e )
            {
            if ( cpc->GetValue( conn[connOffset+e] ) )
              continue; // we've already computed these critical points

            int faceWithEdge;
            for ( faceWithEdge = 0; ! psFace[this->Meta->FacesOnEdge[e][faceWithEdge]]; ++faceWithEdge )
              ; // guaranteed to terminate with a valid value because otherwise edge would have been handled
            vtkPolynomialSystem* bps = this->Meta->RestrictFaceToEdge( psFace[this->Meta->FacesOnEdge[e][faceWithEdge]], e, faceWithEdge );
            bps->SolveSystem();
            perm = this->CellRecord.GetEdgePermutation( e );
            this->StoreRelevantRoots( bps, conn, connOffset + e, perm, cp[i] );
            bps->Delete();
            if ( i == nc - 1 )
              {
              cpc->SetValue( conn[connOffset+e], 1 );
              }
            }
          for ( int j = 0; j < nf; ++j )
            {
            if ( psFace[j] )
              psFace[j]->Delete();
            }
          delete [] psFace;
#else // 0
          this->StoreRelevantRoots( ps, conn, connOffset, 0, cp[i], att );
          // Critical points on faces
          int nf = this->GetNumberOfBoundaries(2);
          connOffset -= nf; // move to the first face
          for ( int f=0; f<nf; ++f )
            {
            if ( cpc->GetValue( conn[npts + connOffset+f] ) )
              continue; // we've already computed these critical points

            vtkPolynomialSystem* bps = this->RestrictToBoundary( ps, connOffset + f );
#ifdef VTK_PRINT_POLSYS
            bps->PrintSystem( vtkstd::cout );
#endif // VTK_PRINT_POLSYS
            bps->SolveSystem();
            perm = this->CellRecord.GetEdgePermutation( f );
            this->StoreRelevantRoots( bps, conn, connOffset + f, perm, cp[i], att );
            bps->Delete();
            if ( i == nc - 1 )
              {
              cpc->SetValue( conn[npts + connOffset+f], 1 );
              }
            }

          // Critical points on edges
          int ne = this->GetNumberOfBoundaries(1);
          connOffset -= ne; // move to the first edge
          for ( int e=0; e<ne; ++e )
            {
            if ( cpc->GetValue( conn[npts + connOffset+e] ) )
              continue; // we've already computed these critical points

            vtkPolynomialSystem* bps = this->RestrictToBoundary( ps, connOffset + e );
#ifdef VTK_PRINT_POLSYS
            bps->PrintSystem( vtkstd::cout );
#endif // VTK_PRINT_POLSYS
            bps->SolveSystem();
            perm = this->CellRecord.GetEdgePermutation( e );
            this->StoreRelevantRoots( bps, conn, connOffset + e, perm, cp[i], att );
            bps->Delete();
            if ( i == nc - 1 )
              {
              cpc->SetValue( conn[npts + connOffset+e], 1 );
              }
            }
#endif // 0
          }
        break;
      }
    }
  ps->Delete();
}

void vtkShoeCell::CopyFrom( const vtkShoeCell* src )
{
  this->SetDataSet( src->DataSet );
  this->Meta = src->Meta; // not a vtkObject; not reference-counted
  this->CellID = src->CellID;

  // Copy the cached values if the cache exists and isn't dirty
  if ( this->AttributeCacheLength )
    {
    for ( int i=0; i<this->AttributeCacheLength; ++i )
      {
      if ( this->AttributeCache[i] )
        {
        this->AttributeCache[i]->Delete();
        }
      }
    if ( this->AttributeCacheLength != src->AttributeCacheLength )
      {
      delete [] this->AttributeCache;
      delete [] this->IsAttributeDirty;
      this->AttributeCacheLength = 0;
      this->AttributeCache = 0;
      this->IsAttributeDirty = 0;
      }
    }
  if ( src->AttributeCacheLength )
    {
    if ( ! this->AttributeCacheLength )
      {
      this->AttributeCacheLength = src->AttributeCacheLength;
      this->AttributeCache = new vtkDoubleArray* [this->AttributeCacheLength];
      this->IsAttributeDirty = new int [this->AttributeCacheLength];
      }
    for ( int c=0; c<src->AttributeCacheLength; ++c )
      {
      if ( src->AttributeCache[c] && !src->IsAttributeDirty[c] )
        {
        this->AttributeCache[c] = vtkDoubleArray::New();
        this->AttributeCache[c]->DeepCopy( src->AttributeCache[c] );
        this->IsAttributeDirty[c] = 0;
        }
      else
        {
        this->AttributeCache[c] = 0;
        this->IsAttributeDirty[c] = 1;
        }
      }
    }
}

void vtkShoeCell::InvalidateCache()
{
  for ( int i=0; i<this->AttributeCacheLength; ++i )
    {
    this->IsAttributeDirty[i] = 1;
    }
}

