/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <iostream>
#include <string>
#include <iterator>

#include <vtkDataArray.h>
#include <vtkIdTypeArray.h>
#include <vtkDoubleArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkPointData.h>
#include <vtkGenericCell.h>
#include <vtkObjectFactory.h>

#include <vtkCellOps.h>
#include <vtkShoeMesh.h>
#include <vtkFunctionData.h>
#include <vtkShoeMeshIterator.h>

#define DBG_LINKS
#define DBG_CELLS

using namespace shoe; // for CellShape

vtkCxxRevisionMacro( vtkShoeMesh, "$Revision: 8692 $" );
vtkStandardNewMacro( vtkShoeMesh );

vtkShoeMesh::vtkShoeMesh()
  :
  FunctionData( vtkFunctionData::New() ),
  GeometryData( vtkFunctionData::New() ),
  NumberOfDofNodes( 0 ),
  Name(0),
  Description(0)
{
  // Take official ownership of the memory
  this->FunctionData->Register( this );
  this->FunctionData->Delete();

  this->GeometryData->Register( this );
  this->GeometryData->Delete();

  // Geometry function data should always be a single (vector) array, which we allocate by default
  this->NumberOfDofNodes=0;
  this->GeometryData->AllocateArrays( 1 );
  vtkDoubleArray* dtmp = vtkDoubleArray::New();
  vtkIdTypeArray* itmp = vtkIdTypeArray::New();
  this->GeometryData->InsertValues( 0, dtmp );
  this->GeometryData->InsertOffsets( 0, itmp );
  dtmp->Delete();
  itmp->Delete();
  this->GeometryData->GetValues(0)->SetNumberOfComponents(3);
  this->GeometryData->GetOffsets(0)->SetNumberOfComponents(1);
  this->GeometryData->GetOffsets(0)->SetNumberOfTuples(this->GetNumberOfDofNodes());
}

vtkShoeMesh::~vtkShoeMesh()
{
  this->FunctionData->UnRegister( this );
  this->GeometryData->UnRegister( this );
}

void vtkShoeMesh::PrintSelf(vtkstd::ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Connectivity: " << this->Connectivity.size() << " entries" << vtkstd::endl;
#ifdef DBG_CELLS
  os << indent.GetNextIndent() ;
  for ( ConnectivityType::iterator it = this->Connectivity.begin(); it != this->Connectivity.end(); ++it )
    os << *it << " ";
  os << vtkstd::endl;
#endif // DBG_CELLS
  os << indent << "Cells: " << this->Cells.size() << " entries" << vtkstd::endl;
#ifdef DBG_CELLS
  os << indent.GetNextIndent() ;
  for ( CellsType::iterator it = this->Cells.begin(); it != this->Cells.end(); ++it )
    os << "( " << it->Offset << ", " << it->GetCellPermutation() << " ) ";
  os << vtkstd::endl;
#endif // DBG_CELLS
  os << indent << "FunctionData: " << this->FunctionData;
  if ( this->FunctionData )
    os << " ( " << this->FunctionData->GetNumberOfArrays() << " arrays )" << vtkstd::endl;
  else
    os << vtkstd::endl;
  this->FunctionData->PrintSelf(os,indent.GetNextIndent());

  os << indent << "GeometryData: " << this->GeometryData;
  if ( this->GeometryData )
    os << " ( " << this->GeometryData->GetNumberOfArrays() << " arrays )" << vtkstd::endl;
  else
    os << vtkstd::endl;
  this->GeometryData->PrintSelf(os,indent.GetNextIndent());

  os << indent << "CornerLinks:  " << this->CornerLinks.size() << " entries" << vtkstd::endl;
#ifdef DBG_LINKS
  int i;
  i=0;
  for ( LinkType::const_iterator it = CornerLinks.begin(); it != CornerLinks.end(); ++it, ++i ) {
    os << indent << "  " << i << " : ";
    copy( it->begin(), it->end(), vtkstd::ostream_iterator<vtkIdType>( os, " " ) );
    os << vtkstd::endl;
  }
#endif // DBG_LINKS
  os << indent << "DofNodeLinks: " << this->DofNodeLinks.size() << " entries" << vtkstd::endl;
#ifdef DBG_LINKS
  i=0;
  for ( LinkType::const_iterator it = DofNodeLinks.begin(); it != DofNodeLinks.end(); ++it, ++i ) {
    os << indent << "  " << i << " : ";
    copy( it->begin(), it->end(), vtkstd::ostream_iterator<vtkIdType>( os, " " ) );
    os << vtkstd::endl;
  }
#endif // DBG_LINKS
}

vtkShoeMesh::ConnectivityType::iterator vtkShoeMesh::GetCellConnectivityBegin( vtkIdType cell_id_in )
{
  return ConnectivityType::iterator( &this->Connectivity, this->Cells[cell_id_in].Offset );
}

vtkShoeMesh::ConnectivityType::iterator vtkShoeMesh::GetCellConnectivityEnd( vtkIdType cell_id_in )
{
  CellDefConstIterator cdi = this->Cells[cell_id_in].Def;
  const vtkCellOps* ops = cdi->GetCellOps();
  return ConnectivityType::iterator( &this->Connectivity, this->Cells[cell_id_in].Offset + ops->NumberOfPoints + ops->NumberOfDofNodes );
}

void vtkShoeMesh::SetNumberOfDofNodes( vtkIdType n )
{
  if ( n == this->NumberOfDofNodes )
    return;

  this->Modified();

  this->GeometryData->SetNumberOfDofNodes( n );
  if ( this->FunctionData )
  {
    this->FunctionData->SetNumberOfDofNodes( n );
  }

  this->NumberOfDofNodes = n;

  this->DofNodeLinks.resize( this->NumberOfDofNodes );
}

void vtkShoeMesh::SetFunctionData( vtkFunctionData* f )
{
  if ( f == this->FunctionData )
    return;

  if ( this->FunctionData )
    this->FunctionData->UnRegister( this );

  this->FunctionData = f;
  this->Modified();

  if ( this->FunctionData ) {
    this->FunctionData->Register( this );
    if ( this->NumberOfDofNodes == 0 )
      this->SetNumberOfDofNodes( this->FunctionData->GetNumberOfDofNodes() );
    int bad_field = this->FunctionData->VerifySize( this->NumberOfDofNodes );
    if ( bad_field >= 0 )
    {
      vtkErrorMacro( "Field " << bad_field << " has "
          << this->FunctionData->GetOffsets( bad_field )->GetNumberOfTuples() << " DOF nodes but this mesh has "
          << this->NumberOfDofNodes << " DOF nodes." );
    }

  }
  else
  {
    vtkErrorMacro( "vtkShoeMesh pretty much expects FunctionData to be valid at all times. Setting it to NULL is a Bad Idea." );
  }
}

void vtkShoeMesh::SetGeometryData( vtkFunctionData* f )
{
  if ( f == this->GeometryData )
    return;

  if ( this->GeometryData )
    this->GeometryData->UnRegister( this );

  this->GeometryData = f;
  this->Modified();

  if ( this->GeometryData->GetNumberOfArrays() != 1 )
  {
    vtkErrorMacro( "vtkShoeMesh pretty much expects GeometryData to have exactly one array." );
    if ( this->GeometryData->GetNumberOfArrays() == 0 )
    {
      this->GeometryData->AllocateOffsets( 1 );
      vtkDoubleArray* dtmp = vtkDoubleArray::New();
      vtkIdTypeArray* itmp = vtkIdTypeArray::New();
      this->GeometryData->InsertValues( 0, dtmp );
      this->GeometryData->InsertOffsets( 0, itmp );
      dtmp->FastDelete();
      itmp->FastDelete();
      this->GeometryData->GetValues(0)->SetNumberOfComponents(3);
      this->GeometryData->GetOffsets(0)->SetNumberOfComponents(1);
      this->GeometryData->GetOffsets(0)->SetNumberOfTuples(this->GetNumberOfDofNodes());
    }
  }
  else if ( this->GeometryData->GetValues(0)->GetNumberOfComponents() != 3 )
  {
    vtkErrorMacro( "vtkShoeMesh pretty much expects GeometryData to have three components per tuple." );
  }

  if ( this->GeometryData ) {
    this->GeometryData->Register( this );
    if ( this->NumberOfDofNodes == 0 )
      this->SetNumberOfDofNodes( this->GeometryData->GetNumberOfDofNodes() );
    int bad_field = this->GeometryData->VerifySize( this->NumberOfDofNodes );
    if ( bad_field >= 0 )
    {
      vtkErrorMacro( "Geometry field " << bad_field << " has "
          << this->GeometryData->GetOffsets( bad_field )->GetNumberOfTuples() << " DOF nodes but this mesh has "
          << this->NumberOfDofNodes << " DOF nodes." );
    }
  }
  else
  {
    vtkErrorMacro( "vtkShoeMesh pretty much expects GeometryData to be valid at all times. Setting it to NULL is a Bad Idea." );
  }
}

vtkIdType vtkShoeMesh::InsertNextDofNode( int num_modes_this_node_in )
{
  // Not much we can do except hope no one's dorked with our arrays 'cause there's not really
  // any way we can guarantee we get the same ID unless ya keep yer skifty mitts off the nums!
  vtkIdType idg = this->GeometryData->InsertNextDofNode( num_modes_this_node_in );
  vtkIdType idf = this->FunctionData->InsertNextDofNode( num_modes_this_node_in );
  assert( idg == idf || idf == -1 );
  this->NumberOfDofNodes++;

  if ( idg >= vtkIdType(this->DofNodeLinks.size()) )
    this->DofNodeLinks.resize( idg + 1 );
  return idg;
}

vtkShoeMesh::DofNodeType vtkShoeMesh::GetDofNodeType( vtkIdType dof_node_id_in ) const
{
  int nd = this->GetNumberOfDofNodes();
  if ( dof_node_id_in < 0 || dof_node_id_in >= nd )
    return Unused;

  CellSpec c( this->Cells[ DofNodeLinks[dof_node_id_in][0] ] );
  const vtkCellOps* ops = c.Def->GetCellOps();
  const vtkIdType* conn = &(this->Connectivity[ c.Offset ]) + ops->NumberOfPoints;
  for ( int e=0; e<ops->NumberOfEdges; ++e, ++conn )
    if ( *conn == dof_node_id_in )
      return EdgeNode;

  for ( int f=0; f<ops->NumberOfEdges; ++f, ++conn )
    if ( *conn == dof_node_id_in )
      return FaceNode;

  if ( *conn == dof_node_id_in )
      return VolumeNode;

  return Unused;
}

void vtkShoeMesh::GetDofNodeSizes( vtkstd::vector<int>& nodesize_out, int field_id_in ) const
{
  int nd = this->GetNumberOfDofNodes();
  for ( int i=0; i<nd; ++i )
  {
    if ( this->DofNodeLinks[i].empty() )
    {
      nodesize_out.push_back( 0 );
    }
    else
    {
      CellSpec c( this->Cells[ DofNodeLinks[i][0] ] );
      const vtkCellOps* ops = c.Def->GetCellOps();
      const vtkIdType* conn = &(this->Connectivity[ c.Offset ]) + ops->NumberOfPoints;
      int off;
      for ( off=0; (conn[off] != i) && (off < ops->NumberOfDofNodes); ++off )
        ;
      if ( field_id_in < 0 )
        nodesize_out.push_back( ops->GetNumberOfModesPerNode( off, c.Def->GetGeometricOrder() ) );
      else
        nodesize_out.push_back( ops->GetNumberOfModesPerNode( off, c.Def->GetFieldOrder( field_id_in ) ) );
    }
  }
}

bool vtkShoeMesh::DoGeometryExtremaExist()
{
  vtkFunctionData* extrema = this->GetGeometryData()->GetExtrema();
  if ( extrema == 0 )
    return false;

  vtkIdTypeArray* off = extrema->GetOffsets(0);
  if ( off == 0 )
    return false;

  return off->GetNumberOfTuples() != 0;
}

bool vtkShoeMesh::DoGeometryDofNodeExtremaExist( vtkIdType dof_node_id_in )
{
  if ( this->DoGeometryExtremaExist() )
  {
    vtkIdTypeArray* off = this->GetGeometryData()->GetExtrema()->Offsets[0];
    int nc = off->GetNumberOfComponents();
    if ( off->GetValue( nc*dof_node_id_in ) != -1 )
      return true;
    return false;
  }
  else
  {
    return false;
  }
}

void vtkShoeMesh::InsertGeometryExtrema( vtkIdType dof_node_id_in, const vtkstd::vector< vtkstd::vector< double > >& extrema_in )
{
  if ( extrema_in.empty() )
  {
    vtkErrorMacro( "You passed an empty list of extrema!" );
    return;
  }

  if ( extrema_in[0].empty() )
  {
    vtkErrorMacro( "You passed an empty tuple (no components) as an extremum!" );
    return;
  }

  this->GetGeometryData()->InsertDofNodeExtrema( 0, dof_node_id_in, extrema_in );
}

void vtkShoeMesh::GetGeometryExtrema( double*& extrema_out, vtkIdType& num_extrema_this_node_out, vtkIdType dof_node_id_in )
{
  vtkFunctionData* gd = this->GetGeometryData();
  vtkFunctionData* ge = gd->GetExtrema();
  vtkIdTypeArray* off;
  vtkDoubleArray* val;
  if ( !ge || ! (off = ge->GetOffsets(0)) )
  {
    num_extrema_this_node_out = 0;
    extrema_out = 0;
    return;
  }
  int nc = off->GetNumberOfComponents();
  vtkIdType start = off->GetValue( dof_node_id_in * nc );
  if ( start == -1 )
  {
    num_extrema_this_node_out = 0;
    extrema_out = 0;
    return;
  }

  if ( nc > 1 )
    num_extrema_this_node_out = off->GetValue( dof_node_id_in * nc + 1 );

  val = dynamic_cast<vtkDoubleArray*>(ge->GetValues(0));
  extrema_out =  val ? val->GetPointer( start * val->GetNumberOfComponents() ) : 0 ;
}

bool vtkShoeMesh::DoFunctionExtremaExist( int field_in )
{
  vtkFunctionData* extrema = this->GetFunctionData()->GetExtrema();
  if ( extrema == 0 )
    return false;

  vtkIdTypeArray* off = extrema->GetOffsets( field_in );
  if ( off == 0 )
    return false;

  return off->GetNumberOfTuples() != 0;
}

bool vtkShoeMesh::DoFunctionDofNodeExtremaExist( int field_in, vtkIdType dof_node_id_in )
{
  if ( this->DoFunctionExtremaExist( field_in ) )
  {
    vtkIdTypeArray* off = this->GetFunctionData()->GetExtrema()->GetOffsets( field_in );
    int nc = off->GetNumberOfComponents();
    if ( off->GetValue( nc*dof_node_id_in ) == -1 )
      return false;
    return true;
  }
  else
  {
    return false;
  }
}

void vtkShoeMesh::InsertFunctionExtrema( int field_in, vtkIdType dof_node_id_in, const vtkstd::vector< vtkstd::vector< double > >& extrema_in )
{
  if ( extrema_in.empty() )
  {
    vtkErrorMacro( "You passed an empty list of extrema!" );
    return;
  }

  if ( extrema_in[0].empty() )
  {
    vtkErrorMacro( "You passed an empty tuple (no components) as an extremum!" );
    return;
  }

  this->GetFunctionData()->InsertDofNodeExtrema( field_in, dof_node_id_in, extrema_in );
}

void vtkShoeMesh::GetFunctionExtrema( double*& extrema_out, vtkIdType& num_extrema_this_node_out, vtkIdType field_in, vtkIdType dof_node_id_in )
{
  vtkFunctionData* fd = this->GetFunctionData();
  vtkFunctionData* fe = fd->GetExtrema();
  vtkIdTypeArray* off;
  vtkDoubleArray* val;
  if ( !fe || ! (off = fe->GetOffsets( field_in )) )
  {
    num_extrema_this_node_out = 0;
    extrema_out = 0;
    return;
  }
  int nc = off->GetNumberOfComponents();
  vtkIdType start = off->GetValue( dof_node_id_in * nc );
  if ( start == -1 )
  {
    num_extrema_this_node_out = 0;
    extrema_out = 0;
    return;
  }

  if ( nc > 1 )
    num_extrema_this_node_out = off->GetValue( dof_node_id_in * nc + 1 );

  val = dynamic_cast<vtkDoubleArray*>(fe->GetValues( field_in ));
  extrema_out =  val ? val->GetPointer( start * val->GetNumberOfComponents() ) : 0 ;
}

void vtkShoeMesh::GetPermutedDofNodeGeometryExtrema(
  vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, vtkShoeMeshIterator& it )
{
#if 0
  if ( ! this->DoGeometryDofNodeExtremaExist( dof_node_id_in ) )
    it.GetCellOps()->GetCellCriticalPoints( it, -1 );

  this->GetPermutedDofNodeGeometryExtrema( extrema_out, dof_node_id_in, it.GetCellId() );
#else
  (void)extrema_out;
  (void)dof_node_id_in;
  (void)it;
#endif // 0
}


void vtkShoeMesh::GetPermutedDofNodeGeometryExtrema( 
  vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, vtkIdType cell_id )
{
  (void)extrema_out;
  (void)dof_node_id_in;
  (void)cell_id;
  vtkstd::cout<<"GetPermutedDofNodeGeometryExtrema  not implemented "<<vtkstd::endl;
  exit(1);
}

void vtkShoeMesh::GetPermutedDofNodeFunctionExtrema(vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, int field_in, vtkShoeMeshIterator& it)
{
  if ( ! this->DoFunctionDofNodeExtremaExist( field_in, dof_node_id_in ) )
    it.GetCellOps()->GetCellCriticalPoints( it, field_in );
    
  this->GetPermutedDofNodeFunctionExtrema( extrema_out, dof_node_id_in, field_in, it.GetCellId() );
}

void vtkShoeMesh::GetPermutedDofNodeFunctionExtrema(vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, int field_in, vtkIdType cell_id)
{
  
  double *extrema;
  vtkIdType num_extrema;
  // Get a pointer to the extrema array.
  this->GetFunctionExtrema(extrema,num_extrema,field_in,dof_node_id_in);

  
  CellSpec c( this->Cells[ cell_id ] );
  const vtkCellOps* ops = c.Def->GetCellOps();
  const vtkIdType* connectivity = &(this->Connectivity[ c.Offset ]);
  vtkFunctionData *ex = this->GetFunctionData()->GetExtrema();

  // Determine the tuple size of extrema entry.
  int tuple_size = ex->GetValues(field_in)->GetNumberOfComponents();

  // Find the edge id or face id of the dofnode.
  int local_id=0;

  // Determine the node type and find the dof node local edge id or face id.
  DofNodeType x = EdgeNode;
  for ( int cnt = ops->NumberOfPoints; cnt < ops->NumberOfPoints + ops->NumberOfDofNodes; ++cnt, ++local_id )
    if ( connectivity[cnt] == dof_node_id_in )
      break;

  // If the local id is greater than the number of edges then this node is a face node or a volume node
  if(local_id>=ops->NumberOfEdges)
    {
    if ( local_id == ops->NumberOfEdges + ops->NumberOfFaces )
      {
      x = VolumeNode;
      }
    else
      {
      x=FaceNode;
      local_id = local_id-ops->NumberOfEdges;
      }
    }
  else
    {
    x = EdgeNode;
    }
  
  double* exVals = new double[tuple_size];
  bool perm_e;
  int node_permutation;
  switch (x)
    {
  case EdgeNode:
    perm_e = c.GetEdgePermutation( local_id );
    for ( int cnt = 0; cnt < num_extrema * tuple_size; cnt += tuple_size )
      {
      vtkstd::vector<double> tempvec;
      if ( perm_e )
        tempvec.push_back( -extrema[cnt] );
      else
        tempvec.push_back( extrema[cnt] );
      for( int i = 1; i < tuple_size; ++i )
        tempvec.push_back( extrema[i + cnt] );
      extrema_out.push_back( tempvec );
      tempvec.clear();
      }
    break;

  case FaceNode:
    node_permutation = c.GetFacePermutation(local_id);
    for ( int cnt = 0; cnt < num_extrema * tuple_size; cnt += tuple_size )
      {
      vtkstd::vector<double> tempVec;
      for ( int i = 0; i < tuple_size; ++i )
        exVals[i] = extrema[cnt + i];

      double temp;
      switch (node_permutation) 
        {
        // starting point, order of traversal
      case 0: // bottom left, right up
              // no transformations needed
        break;

      case 1: // bottom left, up right
              // swap
        temp             = exVals[0];
        exVals[0] = exVals[1];
        exVals[1] = temp;
        break;

      case 2: // top left, down right
              // swap and negate r
        temp = exVals[0];
        exVals[0] = -exVals[1];  // r gets negated
        exVals[1] = temp;
        break;

      case 3: // top left, right down
              // negate s
        exVals[1] = exVals[1];
        break;

      case 4: // top right, left down
              // negate both
        exVals[0] = -exVals[0];
        exVals[1] = -exVals[1];
        break;

      case 5: // top right, down left
              // swap and negate both
        temp             = exVals[0];
        exVals[0] = -exVals[1];
        exVals[1] = -temp;
        break;

      case 6: // bottom right, up left
              // swap and negate s
              // To depermute negate s and then swap
        temp             = exVals[0];
        exVals[0] = exVals[1];
        exVals[1] = -temp;
        break;

      case 7: // bottom right, left up
              // negate r
        exVals[0] = -exVals[0];
        break ;
        } // switch (node_permutation)

      vtkstd::vector<double> tempvec;
      for ( int i = 0; i < tuple_size; ++i )
        tempvec.push_back( exVals[i] );
      extrema_out.push_back( tempvec );
      tempvec.clear();
      }
      break;

  case VolumeNode:
      for(int cnt=0;cnt<num_extrema*tuple_size;cnt+=tuple_size)
        {
        vtkstd::vector<double> tempvec;
        for ( int i = 0; i < tuple_size; ++i )
          {
          tempvec.push_back( extrema[i + cnt] );
          }
        extrema_out.push_back( tempvec );
        tempvec.clear();
        }
      break;

  case Unused:
      // avoid a compiler warning
      break;
    }
  delete [] exVals;
}

    
vtkShoeMesh::CellDefConstIterator vtkShoeMesh::FindCellDef( const vtkCellType& type_in, const int* geom_order_in, const int* func_order_in ) const
{
  vtkCellDefinition cd( type_in.DomainShape, type_in.Interpolant, type_in.ProductSpace, 0,
                        this->GetNumberOfNonlinearFunctions(), geom_order_in, func_order_in );
  return find( CellDefs.begin(), CellDefs.end(), cd );
}

vtkShoeMesh::CellDefIterator vtkShoeMesh::FindOrCreateCellDef( const vtkCellType& type_in, const int* geom_order_in, const int* func_order_in )
{
  vtkCellDefinition cd( type_in.DomainShape, type_in.Interpolant, type_in.ProductSpace, 0,
                        this->GetNumberOfNonlinearFunctions(), geom_order_in, func_order_in );
  CellDefIterator it = find( CellDefs.begin(), CellDefs.end(), cd );
  if ( it == CellDefs.end() ) {
    this->CellDefs.push_front( cd );
    it = this->CellDefs.begin();
  }
  return it;
}

vtkShoeMesh::CellDefIterator vtkShoeMesh::FindOrCreateCellDef( const vtkShoeMesh::CellDefConstIterator& celldef_in )
{
  CellDefIterator it = find( CellDefs.begin(), CellDefs.end(), *celldef_in );
  if ( it == CellDefs.end() ) {
    this->CellDefs.push_front( *celldef_in );
    it = this->CellDefs.begin();
  }
  return it;
}

vtkShoeMesh::CellDefIterator vtkShoeMesh::FindOrCreateCellDef( const vtkCellDefinition& celldef_in )
{
  CellDefIterator it = find( CellDefs.begin(), CellDefs.end(), celldef_in );
  if ( it == CellDefs.end() ) {
    this->CellDefs.push_front( celldef_in );
    it = this->CellDefs.begin();
  }
  return it;
}

vtkIdType vtkShoeMesh::InsertNextCell( const CellDefIterator& defn_in, const vtkIdType* connectivity_in, uint32_t node_permutations_in )
{
  const vtkCellOps* ops = defn_in->GetCellOps();
  int ncorner = ops->NumberOfPoints;
  int ndof = ops->NumberOfDofNodes;
  int nconn = ncorner + ndof;
  vtkIdType c = this->Cells.grab();
  this->Cells[c].SetCellPermutation( node_permutations_in );
  this->Cells[c].Offset = this->Connectivity.grab( nconn );
  this->Cells[c].Def = defn_in;
  defn_in->Reference( c );
  vtkstd::copy( connectivity_in, connectivity_in + nconn, &(this->Connectivity[ this->Cells[c].Offset ]) );
  for ( int cn=0; cn<ncorner; ++cn )
    this->CornerLinks[ connectivity_in[cn] ].push_back( c );
  for ( int dn=0; dn<ndof; ++dn )
    this->DofNodeLinks[ connectivity_in[ncorner + dn] ].push_back( c );
  return c;
}

int vtkShoeMesh::CreateNewField( int data_type_in, int num_components_in, int* order_in, const char* field_name_in, bool init_values_in )
{
  vtkDataArray*   corners = vtkDataArray::CreateDataArray( data_type_in );
  vtkDataArray*   values  = vtkDataArray::CreateDataArray( data_type_in );
  vtkIdTypeArray* offsets = vtkIdTypeArray::New();

  corners->SetName( field_name_in );
  corners->SetNumberOfComponents( num_components_in );
  corners->SetNumberOfTuples( this->GetNumberOfPoints() );

  offsets->SetName( field_name_in );
  offsets->SetNumberOfComponents( 1 );
  offsets->SetNumberOfTuples( this->GetNumberOfDofNodes() );

  values->SetName( field_name_in );
  values->SetNumberOfComponents( num_components_in );

  int field = this->GetPointData()->AddArray( corners );
  vtkFunctionData* fd = this->GetFunctionData();
  fd->InsertValues( field, values );
  fd->InsertOffsets( field, offsets );
  offsets->FastDelete();
  corners->FastDelete();
  values->FastDelete();

  for ( CellDefIterator di = this->BeginCellDefs(); di != this->EndCellDefs(); ++di )
    di->SetFieldOrder( field, order_in );

  vtkIdType valueOff = 0;
  for ( int dn=0; dn<this->GetNumberOfDofNodes(); ++dn )
  {
    if ( this->DofNodeLinks[ dn ].empty() )
    {
      offsets->SetValue( dn, -1 );
    }
    else
    {
      CellSpec c( this->Cells[ DofNodeLinks[ dn ][0] ] );
      const vtkCellOps* ops = c.Def->GetCellOps();
      const vtkIdType* conn = &(this->Connectivity[ c.Offset ]) + ops->NumberOfPoints;
      int off;
      for ( off=0; (conn[off] != dn) && (off < ops->NumberOfDofNodes); ++off )
        ;
      offsets->SetValue( dn, valueOff );
      valueOff += ops->GetNumberOfModesPerNode( off, order_in );
    }
  }
  fd->NextEntry->SetValue( field, valueOff );
  values->SetNumberOfTuples( valueOff );

  if ( init_values_in )
  {
    for ( vtkIdType t=0; t<num_components_in; ++t )
    {
      values->FillComponent( t, 0. );
      corners->FillComponent( t, 0. );
    }
  }

  return field;
}

void vtkShoeMesh::Reset()
{
  this->Cells.clear();
  this->Connectivity.clear();
  this->DofNodeLinks.clear();
  this->CornerLinks.clear();
  this->CellDefs.clear();
  if ( this->GetGeometryData() )
    this->GetGeometryData()->Reset();
  if ( this->GetFunctionData() )
    this->GetFunctionData()->Reset();
  this->NumberOfDofNodes = 0;
  this->SetName( 0 );
  this->SetDescription( 0 );
}

void vtkShoeMesh::SetUpdateExtent( int vtkNotUsed(piece_in), int vtkNotUsed(num_pieces_in), int vtkNotUsed(ghostLevel) )
{
}

int vtkShoeMesh::GetNumberOfNonlinearFunctions() const
{
  return this->GetFunctionData()->GetNumberOfArrays();
}


void vtkShoeMesh::SqueezeCells()
{
  while ( this->Cells.unused_size() )
  {
    // we know that id is not at the end of the list, so grab_and_assign never changes max_used_entry.
    vtkIdType old_id;
    vtkIdType new_id = this->Cells.grab_and_assign( this->Cells[ (old_id = this->Cells.max_used_entry()) ] );
    const vtkCellOps* ops = this->Cells[old_id].Def->GetCellOps();
    vtkIdType* conn = this->GetCellConnectivity( old_id );

    // update all the CornerLinks and DofNodeLinks entries for the cell
    for ( int p=0; p<ops->NumberOfPoints; p++ )
    {
      LinkEntryType::iterator it = find( CornerLinks[ *conn ].begin(), CornerLinks[ *conn ].end(), old_id ) ;
      if ( it == CornerLinks[ *conn ].end() )
      {
        vtkErrorMacro( "Missing CornerLinks entry for cell " << old_id << " at point " << *conn );
      }
      else
      {
        *it = new_id;
      }
      conn++;
    }

    for ( int h=0; h<ops->NumberOfDofNodes; h++ )
    {
      if ( *conn < 0 )
        continue ; // there are no coefficients associated with this DOF node
      LinkEntryType::iterator it = find( DofNodeLinks[ *conn ].begin(), DofNodeLinks[ *conn ].end(), old_id ) ;
      if ( it == DofNodeLinks[ *conn ].end() )
      {
        vtkErrorMacro( "Missing DofNodeLinks entry for cell " << old_id << " at point " << *conn );
      }
      else
      {
        *it = new_id;
      }
      conn++;
    }

    // free the old cell
    this->Cells.free( old_id );
  }
}

void vtkShoeMesh::SqueezeConnectivity()
{
  vtkWarningMacro( "vtkShoeMesh::SqueezeConnectivity not yet implemented" );
}

vtkCell* vtkShoeMesh::GetCell( vtkIdType id )
{
  static vtkGenericCell* TempCell = vtkGenericCell::New();

  this->GetCell( id, TempCell );
  return TempCell;
}

void vtkShoeMesh::GetCell( vtkIdType id, vtkGenericCell* cell )
{
  if ( (id < 0) || (id > this->Cells.size() ) || Cells.is_unused( id ) )
  {
    cell->SetCellTypeToEmptyCell();
  }
  else
  {
    CellShape s = Cells[id].Def->GetShape();
    vtkIdType* conn = &(this->Connectivity[ this->Cells[id].Offset ]);
    int npts = 0;

    switch ( s )
    {
      case Point:
        cell->SetCellTypeToVertex();
        npts = 1;
        break;
      case Curve:
      case Rod:
      case Tube:
        cell->SetCellTypeToLine();
        npts = 2;
        break;
      case Triangle:
      case TriangleShell:
        cell->SetCellTypeToTriangle();
        npts = 3;
        break;
      case Quadrilateral:
      case QuadrilateralShell:
        cell->SetCellTypeToQuad();
        npts = 4;
        break;
      case Tetrahedron:
        cell->SetCellTypeToTetra();
        npts = 4;
        break;
      case Hexahedron:
        cell->SetCellTypeToHexahedron();
        npts = 8;
        break;
      case Wedge:
        cell->SetCellTypeToWedge();
        npts = 6;
        break;
      case Pyramid:
        cell->SetCellTypeToPyramid();
        npts = 5;
        break;
    }

    cell->PointIds->SetNumberOfIds( npts );
    cell->Points->SetNumberOfPoints( npts );

    double x[3];
    for ( int i=0; i<npts; i++ )
    {
      cell->PointIds->SetId( i, conn[i] );
      this->Points->GetPoint( conn[i], x );
      cell->Points->SetPoint(i, x);
    }
  }

  if ( cell->RequiresInitialization() )
    cell->Initialize();
}

int vtkShoeMesh::GetCellType( vtkIdType id )
{
  if ( (id < 0) || (id > this->Cells.size() ) || Cells.is_unused( id ) )
    return VTK_EMPTY_CELL;

  CellShape s = Cells[id].Def->GetShape();
  switch ( s )
  {
    case Point:
      return VTK_VERTEX;
    case Curve:
    case Rod:
    case Tube:
      return VTK_LINE;
    case Triangle:
    case TriangleShell:
      return VTK_TRIANGLE;
    case Quadrilateral:
    case QuadrilateralShell:
      return VTK_QUAD;
    case Tetrahedron:
      return VTK_TETRA;
    case Hexahedron:
      return VTK_HEXAHEDRON;
    case Wedge:
      return VTK_WEDGE;
    case Pyramid:
      return VTK_PYRAMID;
  }

  return VTK_EMPTY_CELL;
}

void vtkShoeMesh::GetCellPoints( vtkIdType id, vtkIdList* points )
{
  points->Reset();
  if ( (id < 0) || (id > this->Cells.size() ) || Cells.is_unused( id ) )
  {
    return;
  }

  const vtkCellOps* ops = Cells[id].Def->GetCellOps();
  register int npts = ops->NumberOfPoints;
  points->SetNumberOfIds( npts );

  vtkIdType* conn = &(this->Connectivity[ this->Cells[id].Offset ]);
  for ( int i=0; i<npts; i++ )
    points->SetId( i, conn[i] );
}

void vtkShoeMesh::GetPointCells( vtkIdType id, vtkIdList* cells)
{
  cells->Reset();
  if ( (id < 0) || (id > this->GetPoints()->GetNumberOfPoints()) )
    return;

  for ( LinkEntryType::const_iterator it=this->CornerLinks[ id ].begin(); it != this->CornerLinks[ id ].end(); ++it )
    cells->InsertNextId( *it );
}

int vtkShoeMesh::GetMaxCellSize()
{
  // a nonlinear hexahedron has the largest number of connectivity entries
  return 27;
}

void vtkShoeMesh::Squeeze()
{
  if ( this->FunctionData )
    this->FunctionData->Squeeze();

  if ( this->GeometryData )
    this->GeometryData->Squeeze();

  this->SqueezeCells();
  this->SqueezeConnectivity();
}

void vtkShoeMesh::SetPoints( vtkPoints* points_in )
{
  this->Superclass::SetPoints( points_in );

  if ( points_in )
    this->CornerLinks.resize( points_in->GetNumberOfPoints() );
}

vtkShoeMesh::CellSpec::CellSpec()
{
  Offset = -1;
}

int vtkShoeMesh::CellSpec::GetFacePermutation( int face ) const
{
  assert( face < 6 && face >= 0 );
  return (this->NodePermutations >> (3*face)) & 0x7;
}

void vtkShoeMesh::CellSpec::SetFacePermutation( int face, int permutation )
{
  assert( face < 6 && face >= 0 );
  this->NodePermutations &= UINT32_MAX - (0x7 << (face*3));
  this->NodePermutations |= permutation << (face*3);
}

bool vtkShoeMesh::CellSpec::GetEdgePermutation( int edge ) const
{
  return (NodePermutations << edge) & (1 << 31);
}

void vtkShoeMesh::CellSpec::SetEdgePermutation( int edge, bool orientation )
{
  if ( orientation )
    NodePermutations |= 1 << (31 - edge);
  else
    NodePermutations &= UINT32_MAX - ( 1 << (31 - edge) );
}

vtkShoeMesh::const_iterator vtkShoeMesh::Begin()
{
  return vtkShoeMeshIterator( this, this->Cells.begin() );
}

vtkShoeMesh::const_iterator vtkShoeMesh::End()
{
  return vtkShoeMeshIterator( this, this->Cells.end() );
}

vtkShoeMesh::CellsType::iterator vtkShoeMesh::BeginCells()
{
  return this->Cells.begin();
}

vtkShoeMesh::CellsType::iterator vtkShoeMesh::EndCells()
{
  return this->Cells.end();
}

