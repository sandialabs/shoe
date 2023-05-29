/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vector>

#include <vtkMath.h>
#include <vtkDataArray.h>
#include <vtkIdTypeArray.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDelaunay3D.h>

#include <vtkCellOps.h>
#include <vtkShoeMeshIterator.h>
#include <vtkAdaptiveTessellator.h>
#include <vtkShoeMeshSubdivisionAlgorithm.h>
#include <Elements/VelhoFunctions.h>
#include <Elements/Generic.h>
#include <Elements/CriticalPoints.h>
#include <Elements/ShapeFunctions.h>
#include <Elements/Hexahedron.h>


#ifdef FOUND_GINAC
USING_NAMESPACE(GiNaC);
#endif // FOUND_GINAC
using namespace shoe;

// For Legendre polynomials in the product space \f${\cal S}^{p,p,q}(\Omega^{(h)}_{st})\f$,
// we can mask off the edges to determine which parameter determines the number of nodes.
// These edge IDs, with binary equivalents, can be masked with 5.
// 
//           0        0000
//           2        0010
//           8        1000
//          10        1010
//
// These edge IDs, with binary equivalents, can then be masked with 4.
// 
//           1        0001
//           3        0011
//           9        1001
//          11        1011
//
// Finally, these edge IDs, with binary equivalents, are left.  Note that
// this order must be strictly followed.
// 
//           4        0101
//           5        0110
//           6        0110
//           7        0111
//
int vtkShoeElemHexahedron::GetNumberOfEdgeModesPerNodePPQLeg( int ID, const int order[3] )
{
  return ((ID & 5) == 0) ? order[0] - 1 : ((ID & 4) == 0) ? order[1] - 1 : order[2] - 1;
}

// For Legendre polynomials in the product space \f${\cal S}^p(\Omega^{(h)}_{st})\f$,
// all faces should have the same order, so use order[0] to compute the number of nodes.
int vtkShoeElemHexahedron::GetNumberOfFaceModesPerNodePLeg( int ID, const int order[3] )
{
  return ((order[0] - 2) * (order[0] - 3) / 2);
}

// For Legendre polynomials in the product space \f${\cal S}^{p,p,q}(\Omega^{(h)}_{st})\f$,
// we can mask off the faces to determine which parameters determine the number of nodes.
// These face IDs, with binary equivalents, can be masked with 4,
// and have shape functions that are functions of s and t.
// 
//           4        100
//           5        110
//
// The shape functions for these faces are functions of s and t.
// These face IDs, with binary equivalents, can then be masked with 2,
// and have shape functions that are functions of r and t.
// 
//           2        010
//           3        011
//
// Finally, these face IDs, with binary equivalents, are left.
// These faces have shape functions that are functions of r and s.
// 
//           0        000
//           1        001
//
// Note that this order must be strictly followed.
int vtkShoeElemHexahedron::GetNumberOfFaceModesPerNodePPQLeg( int ID, const int order[3] )
{
  int order_r;
  int order_s;

  if (ID & 4)
  {
    order_r = order[1];
    order_s = order[2];
  }
  else if (ID & 2) 
  {
    order_r = order[0];
    order_s = order[2];
  }
  else
  {
    order_r = order[0];
    order_s = order[1];
  }

  if (order_r == order_s)
    return ((order_r - 2) * (order_r - 3) / 2);
  else
    return ((order_r - 1) * (order_s - 1)) ;
}

// For Legendre polynomials in the product space \f${\cal S}^p(\Omega^{(h)}_{st})\f$,
// all parameters should have the same order, so use order[0] to compute the number of nodes.
int vtkShoeElemHexahedron::GetNumberOfVolumeModesPerNodePLeg( const int order[3] )
{
  return ((order[0] - 3) * (order[0] - 4) * (order[0] - 5) / 6) ;
}

// For Legendre polynomials in the product space \f${\cal S}^{p,p,q}(\Omega^{(h)}_{st})\f$,
// two parameter orders are equal, so we just have to identify the one that is different
// to compute the number of nodes.
int vtkShoeElemHexahedron::GetNumberOfVolumeModesPerNodePPQLeg( const int order[3] )
{
  int orderq = (order[0] == order[1]) ? order[2] : (order[0] == order[2]) ? order[1] : order[0];
  int orderp = (order[0] == orderq) ? order[1] : order[0];
  return ((orderp - 3) * (orderp - 2) * (orderq - 1) / 2) ;
}

void vtkShoeElemHexahedron::GetBoundaryEdge( vtkCellDefinition& def_out, vtkIdType* connectivity_out, uint32_t& perm_out, vtkShoeMeshIterator& cell, int edge_number )
{
  int order_swap[] = { 0, 1, 0, 1, 2, 2, 2, 2, 0, 1, 0, 1 };

  // cell operations for current cell
  const vtkCellOps* cellOps = cell.GetCellOps();

  // get the current cell definition
  vtkShoeMesh::CellDefConstIterator this_def = cell.GetCell()->Def;

  // get geometric order of current cell and first field order as place holder
  const int* geom_order = this_def->GetGeometricOrder();
  const int* func_order = this_def->GetFieldOrder( 0 );
  int nfields = this_def->GetNumberOfFields();

  // adjust geometric orders
  int new_geom_order[3];
  new_geom_order[0] = geom_order[order_swap[edge_number]];
  new_geom_order[1] = 0;
  new_geom_order[2] = 0;

  // make new cell definition
  def_out = vtkCellDefinition( Curve, this_def->GetInterpolant(), this_def->GetProductSpace(), 0, nfields, new_geom_order, func_order );

  // adjust function orders in new definition
  int new_func_order[3];
  for (int i = 0; i < nfields; i++)
  {
    func_order = this_def->GetFieldOrder( i );
    new_func_order[0] = func_order[order_swap[edge_number]];
    new_func_order[1] = 0;
    new_func_order[2] = 0;
    def_out.SetFieldOrder( i, new_func_order );
  }

  // declare a cell spec for setting permutations
  vtkShoeMesh::CellSpec temp_spec;

  // permutation is the same as in the parent cell
  temp_spec.SetEdgePermutation( 0, cell.GetCell()->GetEdgePermutation( edge_number ) );
  perm_out = temp_spec.GetCellPermutation();

  // connectivity array of this cell
  const vtkIdType* this_cell_conn = cell.GetMesh()->GetCellConnectivity( cell.GetCellId() );
  
  // get the endpoints of the new cell
  const int* corners = cellOps->EdgeCornerNodes[ edge_number ];
  connectivity_out[0] = this_cell_conn[corners[0]];
  connectivity_out[1] = this_cell_conn[corners[1]];

  // add the edge node
  connectivity_out[2] = this_cell_conn[8+edge_number];
}

// return a quadrilateral cell definition corresponding to a hexadron face with given ID.
void vtkShoeElemHexahedron::GetBoundaryFace( vtkCellDefinition& def_out, vtkIdType* connectivity_out, uint32_t& perm_out, vtkShoeMeshIterator& cell, int face_number )
{
  static const int order_swap[][2] =
  {
    {0, 1},
    {1, 2},
    {2, 0},
    {0, 2},
    {1, 2},
    {2, 1}
  };

  static const bool odd_permutation[] = { false, true, true, false, false, true };

  // cell operations for current cell
  const vtkCellOps* cellOps = cell.GetCellOps();

  // get the current cell definition
  vtkShoeMesh::CellDefConstIterator this_def = cell.GetCell()->Def;

  // get geometric order of current cell and first field order as place holder
  const int* geom_order = this_def->GetGeometricOrder();
  const int* func_order = this_def->GetFieldOrder( 0 );
  int nfields = this_def->GetNumberOfFields();

  // adjust geometric orders
  int new_geom_order[3];
  new_geom_order[0] = geom_order[order_swap[face_number][0]];
  new_geom_order[1] = geom_order[order_swap[face_number][1]];
  new_geom_order[2] = 0;

  // make new cell definition
  def_out = vtkCellDefinition( Quadrilateral, this_def->GetInterpolant(), this_def->GetProductSpace(), 0, nfields, new_geom_order, func_order );

  // adjust function orders in new definition
  int new_func_order[3];
  for (int i = 0; i < nfields; i++)
  {
    func_order = this_def->GetFieldOrder( i );
    new_func_order[0] = func_order[order_swap[face_number][0]];
    new_func_order[1] = func_order[order_swap[face_number][1]];
    new_func_order[2] = 0;
    def_out.SetFieldOrder( i, new_func_order );
  }

  // declare a cell spec for setting permutations
  vtkShoeMesh::CellSpec temp_spec;

  // set face permutation
  int this_perm = (face_number % 2) ? 4 : 0;
  if (odd_permutation[face_number])
    this_perm |= 1;
  this_perm = cellOps->ConcatenateFacePermutations( cell.GetCell()->GetFacePermutation( face_number ), this_perm );
  temp_spec.SetFacePermutation( 0, this_perm );

  // set edge permutations
  const int* face_edges = cellOps->EdgesOfFaces[ face_number ];
  for (int i = 0; i < 4; i++)
    if (face_number % 2)
      temp_spec.SetEdgePermutation( i, !cell.GetCell()->GetEdgePermutation( face_edges[i] ) );
    else
      temp_spec.SetEdgePermutation( i, cell.GetCell()->GetEdgePermutation( face_edges[i] ) );

  // return the permutation
  perm_out = temp_spec.GetCellPermutation();

  // connectivity array of this cell
  const vtkIdType* this_cell_conn = cell.GetMesh()->GetCellConnectivity( cell.GetCellId() );
  
  // get the corners of the new cell
  const int* corner_nodes = cellOps->FaceCornerNodes[ face_number ];
  for (int i = 0; i < 4; i++)
    connectivity_out[i] = this_cell_conn[corner_nodes[i]];

  // add in the edges
  const int* edges = cellOps->EdgesOfFaces[ face_number ];
  for (int i = 0; i < 4; i++)
    connectivity_out[4+i] = this_cell_conn[edges[i]+8];

  // add the face
  connectivity_out[4+4] = this_cell_conn[8+12+face_number];
}

// return the permutation code for face_B of a hexahedron based on its
// mating face's code and the offset of their starting vertices
int vtkShoeElemHexahedron::FaceNodePermutation( int face_A, int face_B, int offset )
{
  // do the right thing, based on the offset
  int perm;
  switch (offset)
  {
    case 0:
      {
        static const bool odd[] = {false, true, true, false, false, true};

        if ((face_A + face_B) % 2)
          perm = 4;
        else
          perm = 0;

        if (odd[face_A] == odd[face_B])
          perm += 1;
      }
      break;

    case 1:
      {
        static const int seq1[] = {7, 2, 6, 3};
        static const int seq2[] = {2, 7, 3, 6};

        if (face_A < 2 || face_A > 3)
          perm = seq1[(face_A * 6 + face_B) % 4];
        else
          perm = seq2[(face_A * 6 + face_B) % 4];
      }
      break;

    case 2:
      {
        static const bool odd[] = {true, false, false, true, true, false};

        if ((face_A + face_B) % 2)
          perm = 0;
        else
          perm = 4;

        if (odd[face_A] == odd[face_B])
          perm += 1;
      }
      break;

    case 3:
      {
        static const int seq1[] = {3, 6, 2, 7};
        static const int seq2[] = {6, 3, 7 ,2};

        if (face_A < 2 || face_A > 3)
          perm = seq1[(face_A * 6 + face_B) % 4];
        else
          perm = seq2[(face_A * 6 + face_B) % 4];
      }
      break;
  }

  return perm;
}

// concatenate two face permutaions
int vtkShoeElemHexahedron::ConcatenateFacePermutations( const int perm_A_in, const int perm_B_in )
{
  return (perm_A_in % 2) ? (perm_A_in + 8 - perm_B_in) % 8 : (perm_A_in + perm_B_in) % 8;
}

// Permute (negate) Legendre polynomial amplitudes as necessary for hexahedral edge with given ID
void vtkShoeElemHexahedron::GetPermutedEdgeSignsLeg( vtkstd::vector<bool>& signs, int edge_ID, bool node_permutation, const int order[3] )
{

  // Determine order for parameter for this edge
  int total_dof = ((edge_ID & 5) == 0) ? order[0] : ((edge_ID & 4) == 0) ? order[1] : order[2];
  
  if ( (int) signs.size() != total_dof - 1)
    signs.resize( total_dof - 1 );

  vtkstd::fill( signs.begin(), signs.end(), false );
  
  if (node_permutation)
    for (int i = 3; i <= total_dof; i += 2)
      signs[i-2] = true;
}

// Based on the permutation bits for the face, determine the indices and signs
// of hexahedral face modes after amplitude swaps and sign changes.
void vtkShoeElemHexahedron::GetPermutedFaceIndicesPLeg( vtkstd::vector<int>& indices, vtkstd::vector<bool>& signs, int face_ID, int node_permutation, const int order[3] )
{
  // for this product space, assume orders are equal for all parameters
  int total_dof = order[0];

  // compute number of face modes for this node
  // if (n_indices < 1) should declare an error
  int n_indices = (total_dof - 2) * (total_dof - 3) / 2;

  // make sure the vectors are big enough
  if ( (int) indices.size() != n_indices)
    indices.resize( n_indices );
  if ( (int) signs.size() != n_indices)
    signs.resize( n_indices );

  // initialize the signs vector -- assume no signs change
  vtkstd::fill( signs.begin(), signs.end(), false );

  // counter for shape functions
  int m = 0;
 
  // this loop steps through lower orders up to total_dof
  // keeps us in the right place in the vectors
  for (int p_local = 3; p_local < total_dof; p_local++)
  {
    // number of higher order shape functions for current value of p
    int k = (p_local - 2) * (p_local - 3) / 2;

    // This loop implements the indexing in Szabo and Babuska, with a little help from
    // intermediate variables k and p_local
    for (int j = m - k + 2, i = p_local - j + 1; ((i + j) == p_local+1 && i >= 2); j++, i--)
    {
      // index for swapping if needed
      int n = i + j - 3;

      // switching on the permutation, handling each case separately
      switch (node_permutation) 
        {
           //       starting point, order of traversal
        case 0: //     bottom left, right up
          // no transformations needed
          indices[m] = m;
          break;
          
        case 1: //    bottom left, up right
          // swap
          indices[m] = n * n - m - 1;
          break;

        case 2: //    top left, down right
          // swap and negate r
          indices[m] = n * n - m - 1;
          if (i % 2)
            signs[m] = true;
          break;

        case 3: //    top left, right down
          // negate s
          indices[m] = m;
          if (j % 2)
            signs[m] = true;
          break;

        case 4: //     top right, left down
          // negate both
          indices[m] = m ;
          if ((i + j) % 2)
            signs[m] = true;
          break;
          
        case 5: //    top right, down left
          // swap and negate both
          indices[m] = n * n - m - 1;
          if ((i + j) % 2)
            signs[m] = true;
          break;

        case 6: //    bottom right, up left
          // swap and negate s
          indices[m] = n * n - m - 1;
          if (j % 2)
            signs[m] = true;
          break;

        case 7: //   bottom right, left up
          // negate r
          indices[m] = m;
          if (i % 2)
            signs[m] = true;
          break ;
        }  // switch
      m++ ;
    } // for loop i and j
  } // for loop on p_local
}

// Return total number of shape functions for hex using tensor product space and legendre shape functions
int vtkShoeElemHexahedron::GetNumberOfShapeFunctionsTensorProduct(const int orders[3]) // ppq space
{
    return 8 + 4*(orders[0]-1) + 4*(orders[1]-1) + 4*(orders[2]-1) + 4*(orders[0]-1)*(orders[2]-1) + (orders[0]-2)*(orders[0]-3) + (orders[0]-3)*(orders[0]-2)*(orders[2]-1)/2;
}

// Return total number of shape functions for hex with legendre interpolation and maximum total order space.
int vtkShoeElemHexahedron::GetNumberOfShapeFunctionsMaxTotalOrder(const int order[3])  // p space
{
  int no =  8 + 12*(order[0]-1);
  if(order[0]>=4)
    no += 3*(order[0]-2)*(order[0]-3);
  if(order[0]>=6)
    no += (order[0]-3)*(order[0]-4)*(order[0]-5)/6;
  return no;
}


// Return total number of shape functions for hex with truncated total order space (~ space)
int vtkShoeElemHexahedron::GetNumberOfShapeFunctionsTruncatedOrder(const int order[3])  // ~ space  // orders = {p,q}
{
  int p=order[0];
  int q=order[1];
  //Number of vertex modes and edge modes
  int no = 8 + 8*(order[0]-1) + 4*(order[1]-1);
  //Number of Face modes
  if(order[0]>=4 && order[2]>=2)
  {
    no += 4*(q*(q-1)/2+(p-q-2)*(q-1)) + (p-2)*(p-3);

  }
  //Number of Volume modes
  for(int cnt=6;cnt<=order[0];cnt++)
  {
    for(int cnt1 = cnt-4;cnt1>=2;cnt1--)
    {
      for(int cnt2 = cnt-cnt1-2;cnt2>=2;cnt2--)
      {
        if(cnt-cnt1-cnt2>order[1])
          continue;
        no++;
      }
    }
  }
  return no;
}

// Evaluate normal to a face of a hexahedron. The normal is calculated by
// taking cross product of the parametric tangents in the two parametric
// directions which span the face. 
// FaceId==6 by default if the function is called for a quadrilaterial element.
void vtkShoeElemHexahedron::EvaluateNormalOnFace( double *norm, vtkShoeMeshIterator& cell, const double r[3], const int FaceId )
{
  int dimension = cell.GetCellOps()->EmbeddingDimension;
  const int *order = cell.GetCellGeometricOrder();
  int no=cell.GetCellOps()->GetNumberOfShapeFunctions(order);
  vtkDataArray *geometry = cell.GetCachedGeometry();
//  cout<<dimension<<" "<<no<<endl;
  double *shape = new double[no*dimension];

  cell.GetCellOps()->EvaluateShapeFuncDerivatives(shape,cell,order,r);
  double n1[3] = {0.,0.,0.};
  double n2[3] = {0.,0.,0.};
  double geom[3];
  for(int cnt=0;cnt<no;cnt++)
  {
    geometry->GetTuple(cnt,geom);
    for(int i=0;i<3;i++)
    {
      n1[i] += geom[i]*shape[vtkShoeElemGeneric::faceIds[FaceId][1]*no+cnt];
      n2[i] += geom[i]*shape[vtkShoeElemGeneric::faceIds[FaceId][2]*no+cnt];
    }
  }
  if(FaceId==1||FaceId==2||FaceId==5 || FaceId==6)
  {
    vtkMath::Cross(n1,n2,norm);
  }
  else
  {
    vtkMath::Cross(n2,n1,norm);
  }

  vtkMath::Normalize(norm);
  delete [] shape;
}

// Function to permute the extrema cordinates stored at the node to match with coordinate system of a cell.
void vtkShoeElemHexahedron::GetPermutedCriticalPoints( vtkstd::vector<vtkstd::vector<double> >& extrema, vtkShoeMeshIterator& iter, int field_no )
{
  if(!iter.GetMesh()->DoFunctionExtremaExist(field_no))
    {
    vtkGenericWarningMacro( << "Critical points for field " << field_no << " do not exist. Computing them now." );
    for ( vtkShoeMeshIterator it( iter.GetMesh()->Begin(),MeshOrder); it != iter.GetMesh()->End(); ++it )
      it.GetCellOps()->GetCellCriticalPoints( it, field_no );
    }

  // Get the extrema 
  vtkFunctionData *extremaData = iter.GetMesh()->GetFunctionData()->GetExtrema();
  //Get the offset array
  vtkIdTypeArray *offsets = extremaData->GetOffsets(field_no);
  // Get the extrema values (coordinates, type and value)
  vtkDoubleArray *values = dynamic_cast<vtkDoubleArray *>( extremaData->GetValues(field_no) );
  int components = values->GetNumberOfComponents();
  //  Add the volume critical points.
  double tup[2];
  const vtkIdType* this_cell_conn = iter.GetMesh()->GetCellConnectivity( iter.GetCellId() );
  
  // Get the offset tuple (offset,number of critical points)
  offsets->GetTuple( this_cell_conn[26], tup );
  double *exVals = new double[components];
  
  for ( int cnt=0; cnt < tup[1]; cnt++ )
    {
    values->GetTuple( (int)tup[0]+cnt, exVals );
    vtkstd::vector<double> tempVec;
    for( int i=0; i < components; i++ )
      {
      tempVec.push_back( exVals[i] );
      }
    extrema.push_back( tempVec );
    tempVec.clear();
    }
  
  // Add face critical points.
  for ( int Fcnt=0; Fcnt < 6; Fcnt++ )
    {
    int node_permutation = iter.GetCell()->GetFacePermutation( Fcnt );
    // Get the node id of the higher order node
    int node_id = this_cell_conn[8+12+Fcnt];
    // Get the tuple containing the offset and number of critical points
    offsets->GetTuple( node_id, tup );
    //For each critical point for that node
    for ( int cnt=0; cnt < tup[1]; cnt++ )
      {
      values->GetTuple( int(tup[0])+cnt, exVals );
      vtkShoeElemGeneric::DepermuteQuadParameters( exVals, node_permutation );

      vtkstd::vector<double> tempVec;
      int tempVal=0;

      for ( int i = 0; i < 3; ++i )
        {
        if ( i == vtkShoeElemGeneric::faceConstantPara[Fcnt][0] )
          tempVec.push_back( vtkShoeElemGeneric::faceConstantPara[Fcnt][1] );
        else
          tempVec.push_back( exVals[tempVal++] );
        }

      for ( int cnt = 3; cnt < components; ++cnt )
        tempVec.push_back(exVals[cnt]);

      extrema.push_back(tempVec);
      tempVec.clear();
      }
    }


  // Add edge points 
  bool perm_e;
  for ( int Ecnt=0; Ecnt < 12; Ecnt++ )
    {
    perm_e = iter.GetCell()->GetEdgePermutation( Ecnt );
    int node_id = this_cell_conn[8+Ecnt];
    offsets->GetTuple( node_id, tup );
    for ( int cnt=0; cnt < tup[1]; cnt++ )
      {
      vtkstd::vector<double> tempVec;
      //cout<<tup[0]<<" "<<cnt<<" "<<(int)tup[0]/5+cnt<<endl;
      if ( tup[0] == -1 )
        break;
      values->GetTuple( (int) tup[0]+cnt, exVals );
      if ( perm_e )
        {
        exVals[0] = -exVals[0];
        }
      
      int temp=0;
      for ( int j=0; j < 3; j++ )
        {
        if ( j != vtkShoeElemGeneric::edgeConnectivity[Ecnt][2] )
          {
          tempVec.push_back( vtkShoeElemGeneric::edgeConnectivity[Ecnt][3+temp++] );
          }
        else
          {
          tempVec.push_back( exVals[0] );
          }
        }
      for ( int j=3; j < components; j++ )
        {
        tempVec.push_back( exVals[j] );
        }
      extrema.push_back( tempVec );
      tempVec.clear();
      }
    }

  /*
  for(int cnt=0;cnt<extrema.size();cnt++)
    {
    vtkstd::cout<<"extrema "<<extrema[cnt][0]<<" "<<extrema[cnt][1]<<" "<<extrema[cnt][2]<<vtkstd::endl;
    }
  */

  delete [] exVals;
}

void vtkShoeElemHexahedron::GetCriticalPoints( vtkShoeMeshIterator& iter, int field_no )
{
  int FaceCorners[4][3] = {
    {-1,-1, 0},
    { 1,-1, 0},
    { 1, 1, 0},
    {-1, 1, 0}
  };
  int EdgeCorners[2][3] = {
    {-1, 0, 0},
    { 1, 0, 0}
  };
  int VolumeCorners[8][3] = {
    {-1,-1,-1},
    { 1,-1,-1},
    { 1, 1,-1},
    {-1, 1,-1},
    {-1,-1, 1},
    { 1,-1, 1},
    { 1, 1, 1},
    {-1, 1, 1}
  };

  vtkstd::vector<vtkstd::vector<double> > coef; //  Vector storing the coefficients.
  vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >powers; // Vector storing the exponents of the field polynomial
  // Get the cellops
  const vtkCellOps *ops;
  ops = iter.GetCellOps();

  // Get the degrees of freedom for the field
  vtkDoubleArray *data = iter.GetCachedFunction(field_no);
  int components = data->GetNumberOfComponents();

  shoe::CellShape eleShape = iter.GetCell()->Def->GetType().DomainShape;

   //  Get the field polynomial
  ops->GetFieldPolynomial(coef,powers,iter,field_no);     
  int max[3];
  
  // Find the critical points of the field polynomial
  double param_bounds[3] = {-1,1,0};  // the range is  -1 to 1.
  
  vtkstd::vector<vtkstd::vector<double> >::iterator ic = coef.begin();
  vtkstd::vector<vtkstd::vector<double> > cornerValues;    // Vector storing the corner values
  // Get the corner values.
  for(int cnt=0;cnt<ops->NumberOfPoints;cnt++)
    {
    vtkstd::vector<double> tempvec;
    for(int i=0;i<components;i++)
      {
      tempvec.push_back(data->GetValue(components*cnt+i));
      }
    cornerValues.push_back(tempvec);
    }
  
  // Vectors storing extrema and types for each of the components.
  vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > > Vextrema;
  vtkstd::vector<vtkstd::vector<int> > Vtypes;

  ops->FindFieldCriticalPoints( Vextrema, Vtypes, max, coef, powers, ops->EmbeddingDimension );
  ops->FindBoundedCriticalValues( Vextrema, Vtypes, ops->EmbeddingDimension, param_bounds, eleShape, cornerValues );
  
  vtkstd::vector<vtkstd::vector<double> > ListOfExtrema;

  // Change the corner point ids to the corner point parametric coordinates. Prepare a list
  // of extrema to pass to the DofNodeExtrema functions.
  for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
    {
    for ( int i = 0; i < (int) Vextrema[Ccnt].size(); ++i )
      {
      vtkstd::cout
        << Vextrema[Ccnt][i][0] << " " << Vextrema[Ccnt][i][1] << " "
        << Vextrema[Ccnt][i][2] << " " << Vextrema[Ccnt][i][3] << " "
        << Vextrema[Ccnt][i][4] << vtkstd::endl;
      if ( (((int)Vextrema[Ccnt][i][Vextrema[Ccnt][i].size()-1])&7)==4 ||
           (((int)Vextrema[Ccnt][i][Vextrema[Ccnt][i].size()-1])&7)==5 )
        {
        int id = (int)Vextrema[Ccnt][i][0];
        for(int j=0;j<3;j++)
          {
          Vextrema[Ccnt][i][j]=VolumeCorners[id][j];
          }
        }
      ListOfExtrema.push_back(Vextrema[Ccnt][i]);
      }
    }
  vtkCellDefinition tempDef;
  
  vtkIdType conn[10];
  
  uint32_t perm;
  const vtkIdType* this_cell_conn = iter.GetMesh()->GetCellConnectivity( iter.GetCellId() );
  // Insert the extrema in the data structure
  if ( ListOfExtrema.size() > 0 )
    iter.GetMesh()->GetFunctionData()->InsertDofNodeExtrema(field_no,this_cell_conn[26],ListOfExtrema);

  // clear everything.
  Vextrema.clear();
  Vtypes.clear();
  ListOfExtrema.clear();
  cornerValues.clear();
  Vextrema.clear();
  coef.clear();
  powers.clear();
  
  // For each face of the cell.
  for(int Fcnt=0;Fcnt<ops->NumberOfFaces;Fcnt++)
    { 
    vtkstd::cout<<"Face "<<Fcnt<<vtkstd::endl;
    // Skip the face if the field extrema exist.
    if ( iter.GetMesh()->DoFunctionDofNodeExtremaExist(field_no,this_cell_conn[8+12+Fcnt]) )
    {
      continue;
    }
    vtkstd::vector<vtkstd::vector<double> > Fcoef;          // Vector of coefficients of the polynomial
    vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > > Fpower;  // Vector of powers of the polynomial
    vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > >Fextrema;  // Vector of critical point coordinated
    vtkstd::vector<vtkstd::vector<int> > Ftypes; //           Vector giving the types of the critical points
    // Get the boundary face
    ops->GetBoundaryFace(tempDef,conn,perm,iter,Fcnt);
    // Insert it in temparary mesh.
    vtkShoeMesh *tempMesh = vtkShoeMesh::New();
    tempMesh->SetPoints( iter.GetMesh()->GetPoints() );
    tempMesh->SetFunctionData( iter.GetMesh()->GetFunctionData() );
    tempMesh->SetGeometryData( iter.GetMesh()->GetGeometryData() );
    tempMesh->SetFieldData( iter.GetMesh()->GetFieldData() );
    tempMesh->GetPointData()->ShallowCopy( iter.GetMesh()->GetPointData() );
    vtkShoeMesh::CellDefIterator i = tempMesh->FindOrCreateCellDef( tempDef.GetType(), tempDef.GetGeometricOrder(),tempDef.GetFieldOrder(field_no) );
    tempMesh->InsertNextCell( i, conn, perm );
    vtkShoeMeshIterator Titer( tempMesh->Begin(), MeshOrder );
    Titer.GetCellOps()->GetFieldPolynomial( Fcoef, Fpower, Titer, field_no );
    vtkstd::vector<vtkstd::vector<double> > FcornerValues;
    
    // Get The corner values
    for(int cnt=0;cnt<4;cnt++)
      {
      vtkstd::vector<double> tempvec;
      for(int i=0;i<components;i++)
        {
        tempvec.push_back(data->GetValue(components*vtkShoeElemGeneric::faceConnectivity[Fcnt][cnt]+i));
        }
      FcornerValues.push_back(tempvec);
      }

    ops->FindFieldCriticalPoints(Fextrema,Ftypes,max,Fcoef,Fpower,ops->EmbeddingDimension-1);
    ops->FindBoundedCriticalValues(Fextrema,Ftypes,ops->EmbeddingDimension-1,
              param_bounds,Quadrilateral,FcornerValues);
    
    //cout<<"Number of face extrema "<<Fextrema[0].size()<<" "<<Fextrema[1].size()<<" "<<Fextrema[2].size()<<endl;
    for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
      {
      for ( int i = 0; i < (int) Fextrema[Ccnt].size(); ++i )
        {
        vtkstd::vector<double>::iterator tempvec = Fextrema[Ccnt][i].begin();
        tempvec+=2;
        Fextrema[Ccnt][i].insert(tempvec,0);
        if ( (((int)Fextrema[Ccnt][i][Fextrema[Ccnt][i].size()-1])&7)==4 || ((int)(Fextrema[Ccnt][i][Fextrema[Ccnt][i].size()-1])&7)==5 )
          {
          int id = (int)Fextrema[Ccnt][i][0];
          for(int j=0;j<3;j++)
            {
            //cout<<"id "<<id<<" j "<<j<<" "<<FaceCorners[id][j]<<endl;
            Fextrema[Ccnt][i][j]=FaceCorners[id][j];
            }

          }


        double para[3];
        int temp1=0;
        //double grad[3];
        for(int cnt=0;cnt<3;cnt++)
          {
          if(cnt==vtkShoeElemGeneric::faceConstantPara[Fcnt][0])
            {
            para[cnt]=vtkShoeElemGeneric::faceConstantPara[Fcnt][1];
            continue;
            }
          para[cnt]=Fextrema[Ccnt][i][temp1++];
          }
        /*
        iter.GetCellOps()->EvaluateIsosurfaceGradient(grad,iter,0,para);
        vtkstd::cout<<"Face Critical Point "<<para[0]<<" "<<para[1]<<" "<<para[2]<<" "<<Fextrema[Ccnt][i][3]<<" "<<Fextrema[Ccnt][i][4]<<vtkstd::endl;
        vtkstd::cout<<"Gradient "<<grad[0]<<" "<<grad[1]<<" "<<grad[2]<<vtkstd::endl;
        for(int cnt=0;cnt<3;cnt++)
          {
          para[cnt] += grad[cnt]*.1;
          }
        iter.GetCellOps()->EvaluateIsosurfaceGradient(grad,iter,0.,para);
        vtkstd::cout<<"Gradient after step"<<grad[0]<<" "<<grad[1]<<" "<<grad[2]<<vtkstd::endl;
        */
        ListOfExtrema.push_back(Fextrema[Ccnt][i]);
        }
      }

    // Have to depermute the values ........................................................
    for ( int cnt = 0; cnt < (int) ListOfExtrema.size(); ++cnt )
      {
      int node_permutation = Titer.GetCell()->GetFacePermutation(0);
      double temp;
      switch (node_permutation) 
        {
        //       starting point, order of traversal
      case 0: //     bottom left, right up
        // no transformations needed
        break;
      case 1: //    bottom left, up right
        // swap
        temp = ListOfExtrema[cnt][0];
        ListOfExtrema[cnt][0] = ListOfExtrema[cnt][1];
        ListOfExtrema[cnt][1] = temp;
        break;
      case 2: //    top left, down right
        // swap and negate r // So negate r and swap to depermute
        temp = -ListOfExtrema[cnt][0];
        ListOfExtrema[cnt][0] = ListOfExtrema[cnt][1];
        ListOfExtrema[cnt][1] = temp;
        break;
      case 3: //    top left, right down
        // negate s
        ListOfExtrema[cnt][1] = - ListOfExtrema[cnt][1];
        break;
      case 4: //     top right, left down
        // negate both
        ListOfExtrema[cnt][0] = -ListOfExtrema[cnt][0];
        ListOfExtrema[cnt][1] = -ListOfExtrema[cnt][1];
        break;
      case 5: //    top right, down left
        // swap and negate both
        temp = ListOfExtrema[cnt][0];
        ListOfExtrema[cnt][0] = -ListOfExtrema[cnt][1];
        ListOfExtrema[cnt][1] = -temp;
        break;
      case 6: //    bottom right, up left
        // swap and negate s  // To depermute negate s and then swap
        temp = ListOfExtrema[cnt][0];
        ListOfExtrema[cnt][0] = -ListOfExtrema[cnt][1];
        ListOfExtrema[cnt][1] = temp;
        break;
      case 7: //   bottom right, left up
        // negate r
        ListOfExtrema[cnt][0] = - ListOfExtrema[cnt][0];
        break ;
        }  // switch
      }

    // Insert the extrema into the data structure.
    if(ListOfExtrema.size()>0)
      iter.GetMesh()->GetFunctionData()->InsertDofNodeExtrema(field_no,this_cell_conn[8+12+Fcnt],ListOfExtrema);

    //  Clear everything
    Fcoef.clear();
    ListOfExtrema.clear();
    Fpower.clear();
    Fextrema.clear();
    Ftypes.clear();
    tempMesh->Delete();
    FcornerValues.clear();
  }
  
  uint32_t perm_e;
  for ( int Ecnt=0; Ecnt < 12; Ecnt++ )
    {
    vtkstd::cout<<"Edge Number "<<Ecnt<<vtkstd::endl;
    if(iter.GetMesh()->DoFunctionDofNodeExtremaExist(field_no,this_cell_conn[8+Ecnt]))
      {
      //continue;
      }
  
    vtkstd::vector<vtkstd::vector<double> > Ecoef;
    vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > > Epowers;

    // Get boundary edge.
    ops->GetBoundaryEdge( tempDef,conn,perm_e,iter,Ecnt );
    // Insert it in a temporary mesh.
    vtkShoeMesh *tempMesh = vtkShoeMesh::New();
    tempMesh->SetPoints( iter.GetMesh()->GetPoints() );
    tempMesh->SetFunctionData( iter.GetMesh()->GetFunctionData() );
    tempMesh->SetGeometryData( iter.GetMesh()->GetGeometryData() );
    tempMesh->SetFieldData( iter.GetMesh()->GetFieldData() );
    tempMesh->GetPointData()->ShallowCopy( iter.GetMesh()->GetPointData() );
    vtkShoeMesh::CellDefIterator i = tempMesh->FindOrCreateCellDef( tempDef.GetType(), tempDef.GetGeometricOrder(),tempDef.GetFieldOrder(field_no) );
    tempMesh->InsertNextCell( i, conn, perm_e );
    vtkShoeMeshIterator Titer( tempMesh->Begin(),MeshOrder );
    Titer.GetCellOps()->GetFieldPolynomial( Ecoef, Epowers, Titer, field_no );

    //cout<<Ecoef.size()<<" "<<Epowers.size()<<endl;
    vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > > Eextrema;
    vtkstd::vector<vtkstd::vector<int> > Etypes;

    // Find the critical points
    for(int Ccnt=0;Ccnt<components;Ccnt++)
    {
      vtkstd::vector<vtkstd::vector<double> > ComponentExtrema;
      vtkstd::vector<int> ComponentTypes;


      ops->FindUnivariateCriticalPoints( ComponentExtrema, ComponentTypes, Ecoef[Ccnt], Epowers[Ccnt] );
      for ( int i = 0; i < (int) ComponentExtrema.size(); ++i )
      {
        for ( int cnt1 = 0; cnt1 < components; ++cnt1 )
        {
          double value=0.;
          for ( int j = 0; j < (int) Epowers[cnt1].size(); ++j )
          {
            value+=Ecoef[cnt1][j]*::pow( ComponentExtrema[i][0], Epowers[cnt1][j][0] );
          }
          ComponentExtrema[i].push_back( value );
        }
      }
      //tempvec.push_back(0);
      //tempvec.push_back(0);
      Eextrema.push_back(ComponentExtrema);
      Etypes.push_back(ComponentTypes);
      ComponentExtrema.clear();
      ComponentTypes.clear();
    }

    vtkstd::vector<vtkstd::vector<double> > EcornerValues;

    for ( int cnt=0; cnt < 2; cnt++ )
      {
      vtkstd::vector<double> tempvec;
      for ( int i=0; i < components; i++ )
        {
        tempvec.push_back(data->GetValue(components*vtkShoeElemGeneric::edgeConnectivity[Ecnt][cnt]+i));
        }
      EcornerValues.push_back(tempvec);
      }

    // Find bounded critical values.
    ops->FindBoundedCriticalValues(Eextrema,Etypes,ops->EmbeddingDimension-2,param_bounds,Curve,EcornerValues);

    for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
      {
      for ( int i = 0; i < (int) Eextrema[Ccnt].size(); ++i )
        {
        vtkstd::vector<double>::iterator tempVec = Eextrema[Ccnt][i].begin();
        tempVec++;
        Eextrema[Ccnt][i].insert(tempVec,0); 
        tempVec=Eextrema[Ccnt][i].begin();
        tempVec++;
        Eextrema[Ccnt][i].insert(tempVec,0);
        if ( ((int)Eextrema[Ccnt][i][Eextrema[Ccnt][i].size()-1]&7) == 4 || ((int)Eextrema[Ccnt][i][Eextrema[Ccnt][i].size()-1]&7)==5 )
          {
          int id = (int)Eextrema[Ccnt][i][0];
          for (int j=0; j < 3; j++ )
            Eextrema[Ccnt][i][j]=EdgeCorners[id][j];
          }
        ListOfExtrema.push_back(Eextrema[Ccnt][i]);
        }
      }

    // Depermute the values .............................
    for ( int cnt = 0; cnt < (int) ListOfExtrema.size(); ++cnt )
      {
      if ( iter.GetCell()->GetEdgePermutation( Ecnt ) )
        {
        //cout<<"Permuting values "<<endl;
        ListOfExtrema[cnt][0] = -ListOfExtrema[cnt][0];
        }
      }

    iter.GetMesh()->GetFunctionData()->InsertDofNodeExtrema(field_no,this_cell_conn[8+Ecnt],ListOfExtrema);
    EcornerValues.clear();
    Eextrema.clear();
    Etypes.clear();
    Ecoef.clear();
    Epowers.clear();
    ListOfExtrema.clear();
    tempMesh->Delete();
    }
}

bool vtkShoeElemHexahedron::DoCriticalPointsExist(vtkShoeMeshIterator& iter,const vtkstd::vector<double>& coef, const vtkstd::vector<vtkstd::vector<int> >& powers)
{
  int order = 0;
  int dimension=iter.GetCellOps()->EmbeddingDimension;
  /*
  for(int i=0;i<dimension;i++)
    order=maxExponent[i];    FIXME :  Uncomment this when maxExponents are used. Currently they are all zero.
  */  
  order = 8;
  double sign=0;
  bool SignChanged=false;
  double para[3]=  {-1.,-1.,-1.};
  double step = 1./order;
  vtkstd::cout << "Step " << step << vtkstd::endl;

  // For every partial derivative
  for(int cnt=0;cnt<dimension;cnt++)
  {
    sign=0;
    SignChanged=false;
    for(para[0]=-1.;para[0]<=1.;para[0]+=step)
    {
      for(para[1]=-1.;para[1]<=1.;para[1]+=step)
      {
        for(para[2]=-1.;para[2]<=1;para[2]+=step)
        {
        //  cout<<"Para "<<para[0]<<" "<<para[1]<<" "<<para[2]<<endl;
          double value=0.;
          for ( int cnt1 = 0; cnt1 < (int) coef.size(); ++cnt1 )
          {
            if(powers[cnt1][cnt]<1)
              continue;
            double term=1.;
            for(int l=0;l<3;l++)
            {
              if(l==cnt)
                term=term*powers[cnt1][l]*::pow(para[l],powers[cnt1][l]-1);
              else
                term=term*::pow(para[l],powers[cnt1][l]);
            }
//            cout<<"term and coeff "<<term<<" "<<coef[cnt1]<<endl;
            
            value += term*coef[cnt1];

          }
//          cout<<"value "<<value<<endl;
          if(value==0)
            continue;
          if(sign==0)
          {
            sign=value;
            continue;
          }
          if(sign*value<0)
          {
            SignChanged=true;
            break;
          }
        }

        if(SignChanged==true)
          break;

      }
      if(SignChanged==true)
        break;
    }
    if(SignChanged==false)
      return false;
    SignChanged=false;

  }
  return true;
}


void vtkShoeElemHexahedron::GetIsosurface(vtkShoeMeshIterator& iter,vtkAdaptiveTessellator* tess_in, const double Isovalue, const int field_no)
{ 
  double VolumeCorners[8][3] =
  {
    {-1,-1,-1},
    { 1,-1,-1},
    { 1, 1,-1},
    {-1, 1,-1},
    {-1,-1, 1},
    { 1,-1, 1},
    { 1, 1, 1},
    {-1, 1, 1}
  };

  // Polydata for tessalation.
  vtkPolyData *cPoints = vtkPolyData::New();
  vtkPoints *cp = vtkPoints::New();
  cPoints->Allocate(1,1);
  cp->Allocate(1,1);
  iter.CacheGeometry();
  iter.CacheFunction(field_no);
  //vtkDoubleArray* geometry = iter.GetCachedGeometry();
  vtkDoubleArray* field = iter.GetCachedFunction(field_no);
  vtkDoubleArray* values = vtkDoubleArray::New();
  // Add corner nodes
  for(int i=0;i<iter.GetCellOps()->NumberOfPoints;i++)
  {
    double vec[3];
    cp->InsertNextPoint(VolumeCorners[i]);
    field->GetTuple(i,vec);
    values->InsertNextTuple(vec);
  }

  cPoints->SetPoints(cp);
  // Get Extrema
  vtkstd::vector<vtkstd::vector<double> > extrema;
  iter.GetCellOps()->GetPermutedCriticalPoints(extrema,iter,field_no);
  // Add extrema points
  for ( int cnt = 0; cnt < (int) extrema.size(); ++cnt )
  {
    if(extrema[cnt][4]==4)   // Do not add extrema if it is a corner node
      continue;
    double pt[3];
    pt[0]=extrema[cnt][0];pt[1]=extrema[cnt][1];pt[2]=extrema[cnt][2];
    double pt1[3];
    iter.GetCellOps()->EvaluateGeometry(pt1,iter,pt);
    int id=cp->InsertNextPoint(pt);
    (void)id;
    iter.GetCellOps()->EvaluateField(pt1,iter,pt,field_no);
    values->InsertNextTuple(pt1);
  }

  vtkUnstructuredGrid *tetraGrid;
  vtkDelaunay3D *filter = vtkDelaunay3D::New();
  filter->SetInput(cPoints);
  tetraGrid=filter->GetOutput();
  filter->Update();
  
  vtkstd::cout << "Number of Tetrahedra = " << tetraGrid->GetNumberOfCells() << vtkstd::endl;
  // int ii=0;
  for(int ii=0;ii<tetraGrid->GetNumberOfCells();ii++)
  {
    vtkCell *cell = tetraGrid->GetCell(ii);
    if(cell->GetCellType()==VTK_TETRA)
    {
      vtkShoeElemVelhoFunctions::Velho( (vtkTetra *)(cell), iter, Isovalue, field_no,tess_in);
    }
  } 

  extrema.clear();
  filter->Delete();
  cPoints->Delete();
}

void vtkShoeElemHexahedron::EvaluateShapeFunctionsMaxTotalOrder(double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3])
{
  // cout<<"evaluating "<<endl;
  //
  // Calculate required phi functions 
  // cout<<"Order "<<order[0]<<endl;
  static double rpow[15];
  static double spow[15];
  static double tpow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  tpow[0] = r[2];
  for(int cnt=1;cnt<order[0];cnt++)
  {
    rpow[cnt] = rpow[cnt-1]*r[0];
    spow[cnt] = spow[cnt-1]*r[1];
    tpow[cnt] = tpow[cnt-1]*r[2];
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
    vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
    vtkShoeElemShapeFunctions::phiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](tpow);
  }

  double terms[6];
  int index=0;
  terms[0] = 1+r[0];  // 1+r
  terms[1] = 1-r[0];  // 1-r
  terms[2] = 1+r[1];  // 1+s
  terms[3] = 1-r[1];  // 1-s
  terms[4] = 1+r[2];  // 1+t
  terms[5] = 1-r[2];  // 1-t

  // Corner Nodes
  shape[index++]=(.125*terms[1]*terms[3]*terms[5]);   //  .125*(1-r)(1-s)(1-t)
  shape[index++]=(.125*terms[0]*terms[3]*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[index++]=(.125*terms[0]*terms[2]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[index++]=(.125*terms[1]*terms[2]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[index++]=(.125*terms[1]*terms[3]*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[index++]=(.125*terms[0]*terms[3]*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[index++]=(.125*terms[0]*terms[2]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[index++]=(.125*terms[1]*terms[2]*terms[4]);   //  .125*(1-r)(1+s)(1+t)

  // Edge Modes
  // Edge 0 - 1  (1-2)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*terms[5]);      //  sn*tn
  }
  // (dge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*terms[5]);     // rp*tn
  }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*terms[5]);    // sp*tn
  }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*terms[5]);     // rn*tn
  }
  // edge 0 - 4 (1,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*terms[3]);     // rn*sn
  }
  // edge 1-5 (2,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*terms[3]);    // rp *sn
  }
  // edge 2-6 (3,7)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*terms[2]);  
  }
  // edge 3-7 (4,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*terms[2]);     
  }
  // edge 4 - 5 (5,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*terms[4]);   
  }
  // edge 5-6 (4,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*terms[4]);    
  }
  // edge 6-7 (7,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*terms[4]);
  }
  // edge 7-4 (8,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*terms[4]);
  }

  // Face Modes
  // cout<<"Number of shape functions calculated "<<index<<endl;
  for(int Fcnt=0;Fcnt<6;Fcnt++)
  {
    for(int cnt=4;cnt<=order[0];cnt++)
    {
      for(int cnt1=cnt-2;cnt1>=2;cnt1--)
      {
        // cout<<cnt1<<" "<<cnt-cnt1<<endl;
        shape[index++]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      }
    }
  }
  //  cout<<"Number of shape functions calculated "<<index<<endl;

  // Volume modes
  for(int cnt=6;cnt<=order[0];cnt++)
  {
    for(int cnt1 = cnt-4;cnt1>=2;cnt1--)
    {
      for(int cnt2 = cnt-cnt1-2;cnt2>=2;cnt2--)
      {
        // cout<<cnt1<<" "<<cnt2<<" "<<cnt-cnt1-cnt2<<endl;
        shape[index++]=(vtkShoeElemShapeFunctions::phiValue[0][cnt1]*vtkShoeElemShapeFunctions::phiValue[1][cnt2]*vtkShoeElemShapeFunctions::phiValue[2][cnt-cnt1-cnt2]);
      }
    }
  }
/*  for(int i=0;i<index;i++)
    if(shape[i]<1.e-12 && shape[i]>-1.e-12)
      shape[i] = 0.;
*/
  //  cout<<"Number of shape functions calculated "<<index<<endl;
}

void vtkShoeElemHexahedron::EvaluateShapeFunctionsTensorProduct(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3])  // ppq space
{
  static double rpow[15];
  static double spow[15];
  static double tpow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  tpow[0] = r[2];
  for(int cnt=1;cnt<order[0];cnt++)
  {
    rpow[cnt] = rpow[cnt-1]*r[0];
    spow[cnt] = spow[cnt-1]*r[1];
  }
  for(int cnt=1;cnt<order[2];cnt++)
  {
    tpow[cnt] = tpow[cnt-1]*r[2];
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
    vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
  }
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](tpow);
  }

  double terms[6];
  int index=0;
  terms[0] = 1+r[0];  // 1+r
  terms[1] = 1-r[0];  // 1-r
  terms[2] = 1+r[1];  // 1+s
  terms[3] = 1-r[1];  // 1-s
  terms[4] = 1+r[2];  // 1+t
  terms[5] = 1-r[2];  // 1-t

  // Corner Nodes
  shape[index++]=(.125*terms[1]*terms[3]*terms[5]);   //  .125*(1-r)(1-s)(1-t)
  shape[index++]=(.125*terms[0]*terms[3]*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[index++]=(.125*terms[0]*terms[2]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[index++]=(.125*terms[1]*terms[2]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[index++]=(.125*terms[1]*terms[3]*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[index++]=(.125*terms[0]*terms[3]*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[index++]=(.125*terms[0]*terms[2]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[index++]=(.125*terms[1]*terms[2]*terms[4]);   //  .125*(1-r)(1+s)(1+t)

  // Edge Modes
  // Edge 0 - 1  (1-2)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*terms[5]);      //  sn*tn
  }
  // (dge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*terms[5]);     // rp*tn
  }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*terms[5]);    // sp*tn
  }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*terms[5]);     // rn*tn
  }
  // edge 0 - 4 (1,5)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*terms[3]);     // rn*sn
  }
  // edge 1-5 (2,6)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*terms[3]);    // rp *sn
  }
  // edge 2-6 (3,7)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*terms[2]);  
  }
  // edge 3-7 (4,8)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*terms[2]);     
  }
  // edge 4 - 5 (5,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*terms[4]);   
  }
  // edge 5-6 (4,5)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*terms[4]);    
  }
  // edge 6-7 (7,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*terms[4]);
  }
  // edge 7-4 (8,5)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*terms[4]);
  }

  // Face Modes
  // For the z = constant faces, the shape function are similar to max. total order space
  for(int Fcnt=0;Fcnt<2;Fcnt++)
  {
    for(int cnt=4;cnt<=order[0];cnt++)
    {
      for(int cnt1=cnt-2;cnt1>=2;cnt1--)
      {
        // cout<<cnt1<<" "<<cnt-cnt1<<endl;
        shape[index++] =
          .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
          vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      }
    }
  }

  // For the remaining faces, the shape functions are tensor products.
  for(int Fcnt=2;Fcnt<6;Fcnt++)
  {
    for(int cnt=2;cnt<=order[vtkShoeElemGeneric::faceIds[Fcnt][1]];cnt++)
    {
      for(int cnt1=2;cnt1<=order[vtkShoeElemGeneric::faceIds[Fcnt][2]];cnt1++)
      {
        shape[index++] =
          .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
          vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
          vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      }
    }
  }

  // Volume modes i,j = 2,3,4....p-2 . k = 2,3,...,p, i+j = 4,5,6,..,p
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      for(int cnt2=2;cnt2<=order[2];cnt2++)
      {
        shape[index++] =
          vtkShoeElemShapeFunctions::phiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[1][cnt-cnt1]*
          vtkShoeElemShapeFunctions::phiValue[2][cnt2];
      }
    }
  }
  //cout<<index<<endl;
}

void vtkShoeElemHexahedron::EvaluateShapeFunctionsTruncatedOrder(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]) // ~pq space
{
  static double rpow[15];
  static double spow[15];
  static double tpow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  tpow[0] = r[2];
  for(int cnt=1;cnt<order[0];cnt++)
  {
    rpow[cnt] = rpow[cnt-1]*r[0];
    spow[cnt] = spow[cnt-1]*r[1];
    tpow[cnt] = tpow[cnt-1]*r[2];
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
    vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
    vtkShoeElemShapeFunctions::phiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](tpow);
  }

  double terms[6];
  int index=0;
  terms[0] = 1+r[0];  // 1+r
  terms[1] = 1-r[0];  // 1-r
  terms[2] = 1+r[1];  // 1+s
  terms[3] = 1-r[1];  // 1-s
  terms[4] = 1+r[2];  // 1+t
  terms[5] = 1-r[2];  // 1-t

  // Corner Nodes
  shape[index++]=(.125*terms[1]*terms[3]*terms[5]);    // .125*(1-r)(1-s)(1-t)
  shape[index++]=(.125*terms[0]*terms[3]*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[index++]=(.125*terms[0]*terms[2]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[index++]=(.125*terms[1]*terms[2]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[index++]=(.125*terms[1]*terms[3]*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[index++]=(.125*terms[0]*terms[3]*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[index++]=(.125*terms[0]*terms[2]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[index++]=(.125*terms[1]*terms[2]*terms[4]);   //  .125*(1-r)(1+s)(1+t)

  // Edge Modes
  // Edge 0 - 1  (1-2)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*terms[5]);      //  sn*tn
  }
  // (dge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*terms[5]);     // rp*tn
  }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*terms[5]);    // sp*tn
  }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*terms[5]);     // rn*tn
  }
  // edge 0 - 4 (1,5)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*terms[3]);     // rn*sn
  }
  // edge 1-5 (2,6)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*terms[3]);    // rp *sn
  }
  // edge 2-6 (3,7)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*terms[2]);  
  }
  // edge 3-7 (4,8)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*terms[2]);     
  }
  // edge 4 - 5 (5,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*terms[4]);   
  }
  // edge 5-6 (4,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*terms[4]);    
  }
  // edge 6-7 (7,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*terms[4]);
  }
  // edge 7-4 (8,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*terms[4]);
  }

  // Face Modes
  for(int Fcnt=0;Fcnt<6;Fcnt++)
  {
    for(int cnt=4;cnt<=order[0];cnt++)
    {
      for(int cnt1=cnt-2;cnt1>=2;cnt1--)
      {
        // cout<<cnt1<<" "<<cnt-cnt1<<endl;
        if(vtkShoeElemGeneric::faceIds[Fcnt][2]==2 && cnt-cnt1>order[1])
          continue;

        shape[index++] =
          .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
          vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      }
    }
  }

  // Volume modes
  for(int cnt=6;cnt<=order[0];cnt++)
  {

    for(int cnt1 = cnt-4;cnt1>=2;cnt1--)
    {
      for(int cnt2 = cnt-cnt1-2;cnt2>=2;cnt2--)
      {
        // cout<<cnt1<<" "<<cnt2<<" "<<cnt-cnt1-cnt2<<endl;
        if(cnt-cnt1-cnt2>order[1])
          continue;

        shape[index++] =
          vtkShoeElemShapeFunctions::phiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[1][cnt2]*
          vtkShoeElemShapeFunctions::phiValue[2][cnt-cnt1-cnt2];
      }
    }
  }
}

void vtkShoeElemHexahedron::EvaluateShapeFunctionDerivativesMaxTotalOrder(double* shape, vtkShoeMeshIterator& cell, const int order[3], const double r[3])
{
//  Calculate required phi functions 


  int no = cell.GetCellOps()->GetNumberOfShapeFunctions(order);
  vtkstd::fill(shape,shape+3*no,0.);
  static double rpow[15];
  static double spow[15];
  static double tpow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  tpow[0] = r[2];
  for(int cnt=1;cnt<order[0];cnt++)
  {
    rpow[cnt] = rpow[cnt-1]*r[0];
    spow[cnt] = spow[cnt-1]*r[1];
    tpow[cnt] = tpow[cnt-1]*r[2];
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
    vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
    vtkShoeElemShapeFunctions::phiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](tpow);
    vtkShoeElemShapeFunctions::dphiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt](rpow);
    vtkShoeElemShapeFunctions::dphiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt](spow);
    vtkShoeElemShapeFunctions::dphiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt](tpow);
  }
  double terms[6];
  int index=0;
  terms[0] = 1+r[0];      //  1+r
  terms[1] = 1-r[0];  // 1-r
  terms[2] = 1+r[1];  // 1+s
  terms[3] = 1-r[1];  // 1-s
  terms[4] = 1+r[2];  // 1+t
  terms[5] = 1-r[2];  // 1-t
  
  // Corner Node .
  
  shape[index] = (.125*-1*terms[3]*terms[5]);    // .125*(1-r)(1-s)(1-t)
  shape[no+index] = (.125*terms[1]*-1*terms[5]);
  shape[2*no+index++] = (.125*terms[1]*terms[3]*-1);
  
  shape[index]=(.125*terms[3]*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[no+index]=(.125*terms[0]*-1*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[2*no+index++]=(.125*terms[0]*terms[3]*-1);   //  .125*(1+r)(1-s)(1-t)
  
  shape[index]=(.125*terms[2]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[no+index]=(.125*terms[0]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[2*no+index++]=(.125*terms[0]*terms[2]*-1);   //  .125*(1+r)(1+s)(1-t)
  
  shape[index]=(.125*-1*terms[2]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[no+index]=(.125*terms[1]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[2*no+index++]=(.125*terms[1]*terms[2]*-1);   //  .125*(1-r)(1+s)(1-t)
  
  shape[index]=(.125*-1*terms[3]*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[no+index]=(.125*terms[1]*-1*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[2*no+index++]=(.125*terms[1]*terms[3]);   //  .125*(1-r)(1-s)(1+t)
  
  shape[index]=(.125*terms[3]*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[no+index]=(.125*terms[0]*-1*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[2*no+index++]=(.125*terms[0]*terms[3]);   //  .125*(1+r)(1-s)(1+t)
  
  shape[index]=(.125*terms[2]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[no+index]=(.125*terms[0]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[2*no+index++]=(.125*terms[0]*terms[2]);   //  .125*(1+r)(1+s)(1+t)
  
  shape[index]=(.125*-1*terms[2]*terms[4]);   //  .125*(1-r)(1+s)(1+t)
  shape[no+index]=(.125*terms[1]*terms[4]);   //  .125*(1-r)(1+s)(1+t)
  shape[2*no+index++]=(.125*terms[1]*terms[2]);   //  .125*(1-r)(1+s)(1+t)

  // Edge Modes
  // Edge 0 - 1  (1-2)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]*terms[5]);      //  sn*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1*terms[5]);      //  sn*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*-1);      //  sn*tn
  }
  // (edge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[5]);     // rp*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]*terms[5]);     // rp*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*-1);     // rp*tn
  }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]*terms[5]);    // sp*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[5]);    // sp*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*-1);    // sp*tn
  }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1*terms[5]);     // rn*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]*terms[5]);     // rn*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*-1);     // rn*tn
  }
  // edge 0 - 4 (1,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*-1*terms[3]);     // rn*sn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*-1);     // rn*sn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[1]*terms[3]);     // rn*sn
  }
  // edge 1-5 (2,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[3]);    // rp *sn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*-1);    // rp *sn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[0]*terms[3]);    // rp *sn
  }
  // edge 2-6 (3,7)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[2]);   // rp*sp
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]);  
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[0]*terms[2]);  
  }
  // edge 3-7 (4,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*-1*terms[2]);     //rn*sp
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]);     
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[1]*terms[2]);     
  }
  // edge 4 - 5 (5,6)
  for(int cnt=2;cnt<=order[0];cnt++)      // sn*tp
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]*terms[4]);   
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1*terms[4]);   
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]);   
  }
  // edge 5-6 (4,5)              
  for(int cnt=2;cnt<=order[0];cnt++)   // rp*tp
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[4]);    
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]*terms[4]);    
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]);    
  }
  // edge 6-7 (7,8)
  for(int cnt=2;cnt<=order[0];cnt++)    //sp*tp
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]*terms[4]);
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[4]);
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]);
  }
  // edge 7-4 (8,5)
  for(int cnt=2;cnt<=order[0];cnt++) 
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1*terms[4]);
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]*terms[4]);
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]);
  }

  // Face Modes
  
  //  Face 1 : 
  int Fcnt=0;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      // cout<<cnt1<<" "<<cnt-cnt1<<endl;
      shape[index]=
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index]=
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++]=
        -.5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  //     Face : 2
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      // cout<<cnt1<<" "<<cnt-cnt1<<endl;
      shape[index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        .5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  //  Face : 3
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      shape[index] = 
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        -.5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  //  Face : 4
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      shape[index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        .5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  //  Face : 5
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      shape[index] =
        -.5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  //  Face : 6
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      shape[index] =
        .5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  // Volume modes

  for(int cnt=6;cnt<=order[0];cnt++)
  {
    for(int cnt1 = cnt-4;cnt1>=2;cnt1--)
    {
      for(int cnt2 = cnt-cnt1-2;cnt2>=2;cnt2--)
      {
        shape[index] =
          vtkShoeElemShapeFunctions::dphiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[1][cnt2]*
          vtkShoeElemShapeFunctions::phiValue[2][cnt-cnt1-cnt2];
        shape[no+index] =
          vtkShoeElemShapeFunctions::phiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::dphiValue[1][cnt2]*
          vtkShoeElemShapeFunctions::phiValue[2][cnt-cnt1-cnt2];
        shape[2*no+index++] =
          vtkShoeElemShapeFunctions::phiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[1][cnt2]*
          vtkShoeElemShapeFunctions::dphiValue[2][cnt-cnt1-cnt2];
      }
    }
  }

  /*
   // Occasional trouble with shape functions not being close enough to zero forced this hack.
   // But it's commented out because that was causing other problems.
  for(int cnt=0;cnt<3*no;cnt++)
    if(shape[cnt]<1.e-8 && shape[cnt]>-1.e-8)
      shape[cnt]=0.;
  */
  //  cout<<index<<endl;
}

void vtkShoeElemHexahedron::EvaluateShapeFunctionDerivativesTensorProduct(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3])  // ppq space
{
  int no = cell.GetCellOps()->GetNumberOfShapeFunctions(order);
  static double rpow[15];
  static double spow[15];
  static double tpow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  tpow[0] = r[2];
  for(int cnt=1;cnt<order[0];cnt++)
  {
    rpow[cnt] = rpow[cnt-1]*r[0];
    spow[cnt] = spow[cnt-1]*r[1];
  }
  for(int cnt=1;cnt<order[2];cnt++)
  {
    tpow[cnt] = tpow[cnt-1]*r[2];
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
    vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
    vtkShoeElemShapeFunctions::dphiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt](rpow);
    vtkShoeElemShapeFunctions::dphiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt](spow);
  }
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](tpow);
    vtkShoeElemShapeFunctions::dphiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhiDerivative[cnt](tpow);
  }

//  r=r[0];s=r[2];t=r[2];
  double terms[6];
  int index=0;
  terms[0] = 1+r[0];      //  1+r
  terms[1] = 1-r[0];  // 1-r
  terms[2] = 1+r[1];  // 1+s
  terms[3] = 1-r[1];  // 1-s
  terms[4] = 1+r[2];  // 1+t
  terms[5] = 1-r[2];  // 1-t
  // Corner Nodes
  // 
  shape[index] = (.125*-1*terms[3]*terms[5]);    // .125*(1-r)(1-s)(1-t)
  shape[no+index] = (.125*terms[1]*-1*terms[5]);
  shape[2*no+index++] = (.125*terms[1]*terms[3]*-1);
  
  shape[index]=(.125*terms[3]*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[no+index]=(.125*terms[0]*-1*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[2*no+index++]=(.125*terms[0]*terms[3]*-1);   //  .125*(1+r)(1-s)(1-t)
  
  shape[index]=(.125*terms[2]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[no+index]=(.125*terms[0]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[2*no+index++]=(.125*terms[0]*terms[2]*-1);   //  .125*(1+r)(1+s)(1-t)
  
  shape[index]=(.125*-1*terms[2]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[no+index]=(.125*terms[1]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[2*no+index++]=(.125*terms[1]*terms[2]*-1);   //  .125*(1-r)(1+s)(1-t)
  
  shape[index]=(.125*-1*terms[3]*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[no+index]=(.125*terms[1]*-1*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[2*no+index++]=(.125*terms[1]*terms[3]);   //  .125*(1-r)(1-s)(1+t)
  
  shape[index]=(.125*terms[3]*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[no+index]=(.125*terms[0]*-1*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[2*no+index++]=(.125*terms[0]*terms[3]);   //  .125*(1+r)(1-s)(1+t)
  
  shape[index]=(.125*terms[2]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[no+index]=(.125*terms[0]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[2*no+index++]=(.125*terms[0]*terms[2]);   //  .125*(1+r)(1+s)(1+t)
  
  shape[index]=(.125*-1*terms[2]*terms[4]);   //  .125*(1-r)(1+s)(1+t)
  shape[no+index]=(.125*terms[1]*terms[4]);   //  .125*(1-r)(1+s)(1+t)
  shape[2*no+index++]=(.125*terms[1]*terms[2]);   //  .125*(1-r)(1+s)(1+t)

  // Edge Modes
  // Edge 0 - 1  (1-2)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]*terms[5]);      //  sn*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1*terms[5]);      //  sn*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*-1);      //  sn*tn
  }
  // (edge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[5]);     // rp*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]*terms[5]);     // rp*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*-1);     // rp*tn
  }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]*terms[5]);    // sp*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[5]);    // sp*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*-1);    // sp*tn
  }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1*terms[5]);     // rn*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]*terms[5]);     // rn*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*-1);     // rn*tn
  }
  // edge 0 - 4 (1,5)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*-1*terms[3]);     // rn*sn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*-1);     // rn*sn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[1]*terms[3]);     // rn*sn
  }
  // edge 1-5 (2,6)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[3]);    // rp *sn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*-1);    // rp *sn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[0]*terms[3]);    // rp *sn
  }
  // edge 2-6 (3,7)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[2]);  
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]);  
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[0]*terms[2]);  
  }
  // edge 3-7 (4,8)
  for(int cnt=2;cnt<=order[2];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*-1*terms[2]);     
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]);     
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[1]*terms[2]);     
  }
  // edge 4 - 5 (5,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]*terms[4]);   
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1*terms[4]);   
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]);   
  }
  // edge 5-6 (4,5)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[4]);    
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]*terms[4]);    
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]);    
  }
  
  // edge 6-7 (7,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]*terms[4]);
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[4]);
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]);
  }
  // edge 7-4 (8,5)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1*terms[4]);
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]*terms[4]);
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]);
  }

  //  Face 1 : 
  int Fcnt=0;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      shape[index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        -.5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  //     Face : 2
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      shape[index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
      shape[2*no+index++] =
        .5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1];
    }
  }

  // For the remaining faces, the shape functions are tensor products.
  // Face 3:
  Fcnt++;
  for(int cnt=2;cnt<=order[vtkShoeElemGeneric::faceIds[Fcnt][1]];cnt++)
  {
    for(int cnt1=2;cnt1<=order[vtkShoeElemGeneric::faceIds[Fcnt][2]];cnt1++)
    {
      shape[index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[no+index] =
        -.5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      
    }
  }
  //  Face 4:
  Fcnt++;
  for(int cnt=2;cnt<=order[vtkShoeElemGeneric::faceIds[Fcnt][1]];cnt++)
  {
    for(int cnt1=2;cnt1<=order[vtkShoeElemGeneric::faceIds[Fcnt][2]];cnt1++)
    {
      shape[index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[no+index] =
        .5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      
    }
  }
  // Face 5:
  Fcnt++;
  for(int cnt=2;cnt<=order[vtkShoeElemGeneric::faceIds[Fcnt][1]];cnt++)
  {
    for(int cnt1=2;cnt1<=order[vtkShoeElemGeneric::faceIds[Fcnt][2]];cnt1++)
    {
      shape[index] =
        -.5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      
    }
  }
  // Face 6:
  Fcnt++;
  for(int cnt=2;cnt<=order[vtkShoeElemGeneric::faceIds[Fcnt][1]];cnt++)
  {
    for(int cnt1=2;cnt1<=order[vtkShoeElemGeneric::faceIds[Fcnt][2]];cnt1++)
    {
      shape[index] =
        .5*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[no+index] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
      shape[2*no+index++] =
        .5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*
        vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt]*
        vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt1];
    }
  }

  // Volume modes i,j = 2,3,4....p-2 . k = 2,3,...,p, i+j = 4,5,6,..,p
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      for(int cnt2=2;cnt2<=order[2];cnt2++)
      {
        shape[index] =
          vtkShoeElemShapeFunctions::dphiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[1][cnt-cnt1]*
          vtkShoeElemShapeFunctions::phiValue[2][cnt2];
        shape[no+index] =
          vtkShoeElemShapeFunctions::phiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::dphiValue[1][cnt-cnt1]*
          vtkShoeElemShapeFunctions::phiValue[2][cnt2];
        shape[2*no+index++] =
          vtkShoeElemShapeFunctions::phiValue[0][cnt1]*
          vtkShoeElemShapeFunctions::phiValue[1][cnt-cnt1]*
          vtkShoeElemShapeFunctions::dphiValue[2][cnt2];
      }
    }
  }
}




void vtkShoeElemHexahedron::EvaluateTruncatedOrderShapeFunctionsDerivatives(double *shape,  vtkShoeMeshIterator& cell, const int order[3], double r[3]) // ~pq space
{
  static double rpow[15];
  int no = cell.GetCellOps()->GetNumberOfShapeFunctions(order);
  static double spow[15];
  static double tpow[15];
  rpow[0] = r[0];
  spow[0] = r[1];
  tpow[0] = r[2];
  for(int cnt=1;cnt<order[0];cnt++)
  {
    rpow[cnt] = rpow[cnt-1]*r[0];
    spow[cnt] = spow[cnt-1]*r[1];
    tpow[cnt] = tpow[cnt-1]*r[2];
  }
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    vtkShoeElemShapeFunctions::phiValue[0][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](rpow);
    vtkShoeElemShapeFunctions::phiValue[1][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](spow);
    vtkShoeElemShapeFunctions::phiValue[2][cnt]=vtkShoeElemShapeFunctions::EvaluatePhi[cnt](tpow);
  }
//  r=r[0];s=r[2];t=r[2];
  double terms[6];
  int index=0;
  terms[0] = 1+r[0];      //  1+r
  terms[1] = 1-r[0];  // 1-r
  terms[2] = 1+r[1];  // 1+s
  terms[3] = 1-r[1];  // 1-s
  terms[4] = 1+r[2];  // 1+t
  terms[5] = 1-r[2];  // 1-t
  // Corner Nodes
  shape[index] = (.125*-1*terms[3]*terms[5]);    // .125*(1-r)(1-s)(1-t)
  shape[no+index] = (.125*terms[1]*-1*terms[5]);
  shape[2*no+index++] = (.125*terms[1]*terms[3]*-1);
  
  shape[index]=(.125*terms[3]*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[no+index]=(.125*terms[0]*-1*terms[5]);   //  .125*(1+r)(1-s)(1-t)
  shape[2*no+index++]=(.125*terms[0]*terms[3]*-1);   //  .125*(1+r)(1-s)(1-t)
  
  shape[index]=(.125*terms[2]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[no+index]=(.125*terms[0]*terms[5]);   //  .125*(1+r)(1+s)(1-t)
  shape[2*no+index++]=(.125*terms[0]*terms[2]*-1);   //  .125*(1+r)(1+s)(1-t)
  
  shape[index]=(.125*-1*terms[2]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[no+index]=(.125*terms[1]*terms[5]);   //  .125*(1-r)(1+s)(1-t)
  shape[2*no+index++]=(.125*terms[1]*terms[2]*-1);   //  .125*(1-r)(1+s)(1-t)
  
  shape[index]=(.125*-1*terms[3]*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[no+index]=(.125*terms[1]*-1*terms[4]);   //  .125*(1-r)(1-s)(1+t)
  shape[2*no+index++]=(.125*terms[1]*terms[3]);   //  .125*(1-r)(1-s)(1+t)
  
  shape[index]=(.125*terms[3]*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[no+index]=(.125*terms[0]*-1*terms[4]);   //  .125*(1+r)(1-s)(1+t)
  shape[2*no+index++]=(.125*terms[0]*terms[3]);   //  .125*(1+r)(1-s)(1+t)
  
  shape[index]=(.125*terms[2]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[no+index]=(.125*terms[0]*terms[4]);   //  .125*(1+r)(1+s)(1+t)
  shape[2*no+index++]=(.125*terms[0]*terms[2]);   //  .125*(1+r)(1+s)(1+t)
  
  shape[index]=(.125*-1*terms[2]*terms[4]);   //  .125*(1-r)(1+s)(1+t)
  shape[no+index]=(.125*terms[1]*terms[4]);   //  .125*(1-r)(1+s)(1+t)
  shape[2*no+index++]=(.125*terms[1]*terms[2]);   //  .125*(1-r)(1+s)(1+t)

  // Edge Modes
  
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]*terms[5]);      //  sn*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1*terms[5]);      //  sn*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]*-1);      //  sn*tn
  }
  // (edge 1 - 2  (2,3)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[5]);     // rp*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]*terms[5]);     // rp*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]*-1);     // rp*tn
  }
  // Edge 2 - 3  (3,4)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]*terms[5]);    // sp*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[5]);    // sp*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]*-1);    // sp*tn
  }
  // Edge 3 - 0  (4,1)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1*terms[5]);     // rn*tn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]*terms[5]);     // rn*tn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]*-1);     // rn*tn
  }
  // edge 0 - 4 (1,5)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*-1*terms[3]);     // rn*sn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]*-1);     // rn*sn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[1]*terms[3]);     // rn*sn
  }
  // edge 1-5 (2,6)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[3]);    // rp *sn
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]*-1);    // rp *sn
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[0]*terms[3]);    // rp *sn
  }
  // edge 2-6 (3,7)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[2]);  
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[0]);  
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[0]*terms[2]);  
  }
  // edge 3-7 (4,8)
  for(int cnt=2;cnt<=order[1];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*-1*terms[2]);     
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[2][cnt]*terms[1]);     
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::dphiValue[2][cnt]*terms[1]*terms[2]);     
  }
  // edge 4 - 5 (5,6)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[3]*terms[4]);   
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*-1*terms[4]);   
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[3]);   
  }
  // edge 5-6 (4,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[4]);    
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[0]*terms[4]);    
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[0]);    
  }
  // edge 6-7 (7,8)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::dphiValue[0][cnt]*terms[2]*terms[4]);
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[4]);
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[0][cnt]*terms[2]);
  }
  // edge 7-4 (8,5)
  for(int cnt=2;cnt<=order[0];cnt++)
  {
    shape[index]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*-1*terms[4]);
    shape[no+index]=(.25*vtkShoeElemShapeFunctions::dphiValue[1][cnt]*terms[1]*terms[4]);
    shape[2*no+index++]=(.25*vtkShoeElemShapeFunctions::phiValue[1][cnt]*terms[1]);
  }
  // Face Modes
//  Face 1 : 
  int Fcnt=0;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
//        cout<<cnt1<<" "<<cnt-cnt1<<endl;
      shape[index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[no+index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[2*no+index++]=(.5*-1*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
    }
  }
//     Face : 2
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
//        cout<<cnt1<<" "<<cnt-cnt1<<endl;
      shape[index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[no+index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[2*no+index++]=(.5*1*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
    }
  }
//  Face : 3
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      if(cnt-cnt1>order[1])
        continue;
      shape[index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[no+index]=(.5*-1*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[2*no+index++]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
    }
  }
//  Face : 4
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      if(cnt-cnt1>order[1])
        continue;
      shape[index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[no+index]=(.5*1*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[2*no+index++]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
    }
  }
//  Face : 5
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      if(cnt-cnt1>order[1])
        continue;
      shape[index]=(.5*-1*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[no+index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[2*no+index++]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
    }
  }
//  Face : 6
  Fcnt++;
  for(int cnt=4;cnt<=order[0];cnt++)
  {
    for(int cnt1=cnt-2;cnt1>=2;cnt1--)
    {
      if(cnt-cnt1>order[1])
        continue;
      shape[index]=(.5*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[no+index]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      shape[2*no+index++]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::dphiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
    }
  }

/*
  for(int Fcnt=0;Fcnt<6;Fcnt++)
  {
    for(int cnt=4;cnt<=order[0];cnt++)
    {
      for(int cnt1=cnt-2;cnt1>=2;cnt1--)
      {
//        cout<<cnt1<<" "<<cnt-cnt1<<endl;
        if(vtkShoeElemGeneric::faceIds[Fcnt][2]==2 && cnt-cnt1>order[1])
          continue;
        shape[index++]=(.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][1]][cnt1]*vtkShoeElemShapeFunctions::phiValue[vtkShoeElemGeneric::faceIds[Fcnt][2]][cnt-cnt1]);
      }
    }
  }
*/
  // Volume modes

  for(int cnt=6;cnt<=order[0];cnt++)
  {

    for(int cnt1 = cnt-4;cnt1>=2;cnt1--)
    {
      for(int cnt2 = cnt-cnt1-2;cnt2>=2;cnt2--)
      {
//        cout<<cnt1<<" "<<cnt2<<" "<<cnt-cnt1-cnt2<<endl;
        if(cnt-cnt1-cnt2>order[1])
          continue;
        shape[index]=(vtkShoeElemShapeFunctions::phiValue[0][cnt1]*vtkShoeElemShapeFunctions::dphiValue[1][cnt2]*vtkShoeElemShapeFunctions::phiValue[2][cnt-cnt1-cnt2]);
        shape[no+index]=(vtkShoeElemShapeFunctions::phiValue[0][cnt1]*vtkShoeElemShapeFunctions::dphiValue[1][cnt2]*vtkShoeElemShapeFunctions::phiValue[2][cnt-cnt1-cnt2]);
        shape[2*no+index++]=(vtkShoeElemShapeFunctions::phiValue[0][cnt1]*vtkShoeElemShapeFunctions::phiValue[1][cnt2]*vtkShoeElemShapeFunctions::dphiValue[2][cnt-cnt1-cnt2]);
      }
    }
  }

}

// Static array used in the following function.
// First number = the index of the term (1-x),(1+x),(1-y).. to be multiplied for that face.
// Next two numbers = Index of the variables to be passed to phi. 
// Get polynomial function for hexahedron constant max order space.
// orders[0] has the Constant Total order.
// Input - dof - degrees of freedom arranged in the proper order for all the nodes.
// Output - coeff array and powers array to be passed to GetPolynomialCriticalPoints.
void vtkShoeElemHexahedron::GetFieldPolynomialMaxTotalOrder(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num)
{
#ifdef FOUND_GINAC
  int dofcounter=0;
  int maxorder=0;
  const int* orders = cell_in.GetCellFunctionOrder(field_num);
  cell_in.CacheFunction(field_num);
  vtkDoubleArray* field = cell_in.GetCachedFunction(field_num);


  double* dof = field->GetPointer(0);     
  
  int components = field->GetNumberOfComponents();
  symbol x("x");symbol y("y");symbol z("z");
  symbol sym[3];
  sym[0]=x;sym[1]=y;sym[2]=z;
  ex xp = 1+x;
  ex xn = 1-x;
  ex yp = 1+y;
  ex yn = 1-y;
  ex zn = 1-z;
  ex zp = 1+z;
  ex terms[6];
  terms[0]=xp;terms[1]=xn;terms[2]=yp;terms[3]=yn;terms[4]=zp;terms[5]=zn;
/*  for(int cnt=0;cnt<field->GetNumberOfTuples()*field->GetNumberOfComponents();cnt++)
  {
    vtkstd::cout<<dof[cnt]<<" ";
  }
  vtkstd::cout<<vtkstd::endl;
*/  
  for(int CompCnt=0;CompCnt<components;CompCnt++)
  {

    dofcounter=CompCnt;
    ex finalPoly;

    // The corner shape functions 
    finalPoly += (dof[dofcounter]*.125*xn*yn*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yn*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yp*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yp*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yn*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yn*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yp*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yp*zp).expand();
    dofcounter += components;


    // The edge shape functions


    // edge y=-1,z=-1   edge between 1-2
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yn*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;

    }

    // edge x=1,z=-1;    2-3
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }

    //  cout<<"edge y=1, z=-1"<<endl;   3-4
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yp*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }

    //  cout<<"edge x=-1,z=-1"<<endl;  4-1
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }

    //  edge x=-1,y=-1       1-5
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=1,y=-1        2-6
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=1,y=1        3-7
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=-1,y=1        4-8
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }

    //  edge y=-1,z=1     5-6
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yn*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }


    //  edge x=1,z=1;    6-7
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }
    //  edge y=1,z=1    7-8
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yp*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }
    //  edge x=-1,z=1   8-5
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }


    // face modes the sequence of exponents is {2,2} {3,2} {2,3} {4,2} {3,3} {2,4}..... for each face.
    for(int Fcnt=0;Fcnt<6;Fcnt++)
    {
      for(int cnt=4;cnt<=orders[0];cnt++)
      {
        for(int cnt1=cnt-2;cnt1>=2;cnt1--)
        {
          finalPoly += (.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][1]])*
                                              vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][2]])).expand();
          dofcounter += components;
        }
      }
    }


    //      Volume modes 
    for(int cnt=6;cnt<=orders[0];cnt++)
    {
      for(int cnt1=cnt-4;cnt1>=2;cnt1--)
      {
        for(int cnt2=cnt-cnt1-2;cnt2>=2;cnt2--)
        {
          //        cout<<cnt1<<" "<<cnt2<<" "<<cnt-cnt1-cnt2<<endl;
          
          finalPoly+=(dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,x)*vtkShoeElemCriticalPoints::eval_phi(cnt2,y)*vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1-cnt2,z)).expand();

          dofcounter += components;

        }

      }

    }
    
    finalPoly = finalPoly.expand().evalf();
    vtkstd::cout<<"Final Field Polynomial "<<vtkstd::endl<<finalPoly<<vtkstd::endl;
    int deg[2];
    double tempval=0.;

    //    Generate the powers and coefficient vectors
    vtkstd::vector<double> Tcoeff;
    vtkstd::vector<vtkstd::vector<int> > Tpowers;
    bool coef_too_small = false;

    for(int cnt=0;cnt<finalPoly.nops();cnt++)
    {
      vtkstd::vector<int> tempvec;
      ex temp;  
      temp=finalPoly.op(cnt);
      tempvec.push_back(temp.degree(x));
      tempvec.push_back(temp.degree(y));
      tempvec.push_back(temp.degree(z));
      numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.,y==1.,z==1.)).evalf());
/*      if(tempnum<1.e-8 && tempnum>-1.e-8)
      {
        tempvec.clear();
        continue;
      }*/
      Tcoeff.push_back(tempnum.to_double());
      Tpowers.push_back(tempvec);
//      cout<<"Coef "<<Tcoeff[cnt]<<" "<<Tpowers[cnt][0]<<" "<<Tpowers[cnt][1]<<" "<<Tpowers[cnt][2]<<endl;
      tempvec.clear();
    }
  
  
    powers.push_back(Tpowers);
    coeff.push_back(Tcoeff);
    Tpowers.clear();
    Tcoeff.clear();
  }
//  cout<<"Tempval = "<<tempval<<endl;
#endif // FOUND_GINAC
}

// Get polynomial function for hexahedron constant max order space.
// Input - dof - degrees of freedom arranged in the proper order for all the nodes.
// orders - maximum orders in r,s,t.
// Output - coeff array and powers array to be passed to GetPolynomialCriticalPoints.

void vtkShoeElemHexahedron::GetFieldPolynomialTensorProduct(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num)
{
#ifdef FOUND_GINAC
  int dofcounter=0;
  int maxorder=0;
  const int *orders = cell_in.GetCellFunctionOrder(field_num);
  cell_in.CacheFunction(field_num);
  vtkDoubleArray *field = cell_in.GetCachedFunction(field_num);
//  cout<<"Orders "<<orders[0]<<" "<<orders[1]<<" "<<orders[2]<<endl;
  if (field->GetNumberOfComponents()>1)
  {
    vtkstd::cout<<"Vector fields are not supported by critical point functions yet "<<vtkstd::endl;
    return;
  }

  double *dof = field->GetPointer(0);     
  int components = field->GetNumberOfComponents();
  symbol x("x");symbol y("y");symbol z("z");
  symbol sym[3];
  sym[0]=x;sym[1]=y;sym[2]=z;
  ex xp = 1+x;
  ex xn = 1-x;
  ex yp = 1+y;
  ex yn = 1-y;
  ex zn = 1-z;
  ex zp = 1+z;
  ex terms[6];
  terms[0]=xp;terms[1]=xn;terms[2]=yp;terms[3]=yn;terms[4]=zp;terms[5]=zn;
  
  for(int CompCnt=0;CompCnt<components;CompCnt++)
  {

    dofcounter=CompCnt;

    ex finalPoly;

    // The corner shape functions 
    finalPoly += (dof[dofcounter]*.125*xn*yn*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yn*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yp*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yp*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yn*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yn*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yp*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yp*zp).expand();
    dofcounter += components;


    // The edge shape functions
    // edge y=-1,z=-1   edge between 1-2
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yn*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;

    }
    // edge x=1,z=-1;    2-3
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }

    //  cout<<"edge y=1, z=-1"<<endl;   3-4
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yp*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }

    //  cout<<"edge x=-1,z=-1"<<endl;  4-1
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }

    //  edge x=-1,y=-1       1-5
    for(int cnt=2;cnt<=orders[2];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=1,y=-1        2-6
    for(int cnt=2;cnt<=orders[2];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=1,y=1        3-7
    for(int cnt=2;cnt<=orders[2];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=-1,y=1        4-8
    for(int cnt=2;cnt<=orders[2];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }

    //  edge y=-1,z=1     5-6
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yn*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }
    //  edge x=1,z=1;    6-7
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }
    //  edge y=1,z=1    7-8
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yp*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }
    //  edge x=-1,z=1   8-5
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }


    for(int Fcnt=0;Fcnt<2;Fcnt++)
    {
      for(int cnt=4;cnt<=orders[0];cnt++)
      {
        for(int cnt1=cnt-2;cnt1>=2;cnt1--)
        {
          //        cout<<cnt1<<" "<<cnt-cnt1<<endl;
          finalPoly += (.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][1]])*vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][2]])).expand();
          dofcounter += components;
        }
      }
    }


    for(int Fcnt=2;Fcnt<6;Fcnt++)
    {
      for(int cnt=2;cnt<=orders[vtkShoeElemGeneric::faceIds[Fcnt][1]];cnt++)
      {
        for(int cnt1=2;cnt1<=orders[vtkShoeElemGeneric::faceIds[Fcnt][2]];cnt1++)
        {
          finalPoly += (.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt,sym[vtkShoeElemGeneric::faceIds[Fcnt][1]])*vtkShoeElemCriticalPoints::eval_phi(cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][2]])).expand();
          dofcounter += components;
        }
      }
    }


    //      Volume modes 

    for(int cnt=4;cnt<=orders[0];cnt++)
    {
      for(int cnt1=cnt-2;cnt1>=2;cnt1--)
      {
        for(int cnt2=2;cnt2<=orders[2];cnt2++)
        {
          //        cout<<cnt1<<" "<<cnt2<<" "<<cnt-cnt1-cnt2<<endl;
          finalPoly+=(dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,x)*vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1,y)*vtkShoeElemCriticalPoints::eval_phi(cnt2,z)).expand();
          dofcounter += components;

        }

      }

    }
    finalPoly = finalPoly.expand().evalf();     

    //  cout<<finalPoly<<endl;
    int deg[2];
    vtkstd::vector<double> Tcoeff;
    vtkstd::vector<vtkstd::vector<int> > Tpowers;
    //    Generate the powers and coefficient vectors
    for(int cnt=0;cnt<finalPoly.nops();cnt++)
    {
      vtkstd::vector<int> tempvec;
      ex temp;  
      temp=finalPoly.op(cnt);
      tempvec.push_back(temp.degree(x));
      tempvec.push_back(temp.degree(y));
      tempvec.push_back(temp.degree(z));
      //    cout<<tempvec[0]<<" "<<tempvec[1]<<" "<<tempvec[2]<<endl;

      Tpowers.push_back(tempvec);
      numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.,y==1.,z==1.)).evalf());
      Tcoeff.push_back(tempnum.to_double());
    } 
    powers.push_back(Tpowers);
    coeff.push_back(Tcoeff);
    
    

  }
#endif // FOUND_GINAC
}

//This function gives polynomial for truncated order space for hex element.
//It is an intersection of Const Max Order for p and tensor product of (p,p,q) where q<p
//Thus we neglect the shape functions which are greater in degree q in z.
//orders array has {the maximum total order, the limiting order for z}
//
void vtkShoeElemHexahedron::GetFieldPolynomialTruncatedOrder(vtkstd::vector<vtkstd::vector<double> >& coeff, vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& powers, vtkShoeMeshIterator& cell_in, const int field_num)
{
#ifdef FOUND_GINAC
  const int *orders = cell_in.GetCellFunctionOrder(field_num);
  cell_in.CacheFunction(field_num);
  vtkDoubleArray *field = cell_in.GetCachedFunction(field_num);
  if (field->GetNumberOfComponents()>1)
  {
    vtkstd::cout<<"Vector fields are not supported by critical point functions yet "<<vtkstd::endl;
    return;
  }
  double *dof =  field->GetPointer(0);      
  
  int components = field->GetNumberOfComponents();
  int dofcounter=0;
  int maxorder=0;

  symbol x("x");symbol y("y");symbol z("z");
  symbol sym[3];
  sym[0]=x;sym[1]=y;sym[2]=z;
  ex xp = 1+x;
  ex xn = 1-x;
  ex yp = 1+y;
  ex yn = 1-y;
  ex zn = 1-z;
  ex zp = 1+z;
  ex terms[6];
  terms[0]=xp;terms[1]=xn;terms[2]=yp;terms[3]=yn;terms[4]=zp;terms[5]=zn;
  for(int CompCnt=0;CompCnt<components;CompCnt++)
  {
    dofcounter=CompCnt;


    ex finalPoly;

    // The corner shape functions 
    finalPoly += (dof[dofcounter]*.125*xn*yn*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yn*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yp*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yp*zn).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yn*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yn*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xp*yp*zp).expand();
    dofcounter += components;
    finalPoly += (dof[dofcounter]*.125*xn*yp*zp).expand();
    dofcounter += components;


    // The edge shape functions


    // edge y=-1,z=-1   edge between 1-2
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yn*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;

    }

    // edge x=1,z=-1;    2-3
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }

    //  cout<<"edge y=1, z=-1"<<endl;   3-4
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yp*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }

    //  cout<<"edge x=-1,z=-1"<<endl;  4-1
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*zn*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }

    //  edge x=-1,y=-1       1-5
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=1,y=-1        2-6
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*yn*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=1,y=1        3-7
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }
    //  edge x=-1,y=1        4-8
    for(int cnt=2;cnt<=orders[1];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*yp*vtkShoeElemCriticalPoints::eval_phi(cnt,z)).expand();
      dofcounter += components;
    }

    //  edge y=-1,z=1     5-6
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yn*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }


    //  edge x=1,z=1;    6-7
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xp*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }
    //  edge y=1,z=1    7-8
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*yp*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,x)).expand();
      dofcounter += components;
    }
    //  edge x=-1,z=1   8-5
    for(int cnt=2;cnt<=orders[0];cnt++)
    {
      finalPoly += (.25*dof[dofcounter]*xn*zp*vtkShoeElemCriticalPoints::eval_phi(cnt,y)).expand();
      dofcounter += components;
    }


    // face modes the sequence of exponents is {2,2} {3,2} {2,3} {4,2} {3,3} {2,4}..... for each face.
    for(int Fcnt=0;Fcnt<6;Fcnt++)
    {
      for(int cnt=4;cnt<=orders[0];cnt++)
      {
        for(int cnt1=cnt-2;cnt1>=2;cnt1--)
        {
          if(vtkShoeElemGeneric::faceIds[Fcnt][2]==2 && cnt-cnt1>orders[1])
            continue;
          //    cout<<cnt1<<" "<<cnt-cnt1<<endl;
          finalPoly += (.5*terms[vtkShoeElemGeneric::faceIds[Fcnt][0]]*dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][1]])*vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1,sym[vtkShoeElemGeneric::faceIds[Fcnt][2]])).expand();
          dofcounter += components;
        }
      }
    }

    //      Volume modes 
    for(int cnt=6;cnt<=orders[0];cnt++)
    {
      for(int cnt1=cnt-4;cnt1>=2;cnt1--)
      {
        for(int cnt2=cnt-cnt1-2;cnt2>=2;cnt2--)
        {
          if(cnt-cnt1-cnt2>orders[1])
            continue;
          //  cout<<cnt1<<" "<<cnt2<<" "<<cnt-cnt1-cnt2<<endl;
          finalPoly+=(dof[dofcounter]*vtkShoeElemCriticalPoints::eval_phi(cnt1,x)*vtkShoeElemCriticalPoints::eval_phi(cnt2,y)*vtkShoeElemCriticalPoints::eval_phi(cnt-cnt1-cnt2,z)).expand();
          dofcounter += components;
        }
      }

    }
    finalPoly = finalPoly.expand().evalf();     

    //  cout<<finalPoly<<endl;
    int deg[2];

    //    Generate the powers and coefficient vectors
    vtkstd::vector<double> Tcoeff;
    vtkstd::vector<vtkstd::vector<int> > Tpowers;

    for(int cnt=0;cnt<finalPoly.nops();cnt++)
    {
      vtkstd::vector<int> tempvec;
      ex temp;  
      temp=finalPoly.op(cnt);
      tempvec.push_back(temp.degree(x));
      tempvec.push_back(temp.degree(y));
      tempvec.push_back(temp.degree(z));
      //    cout<<tempvec[0]<<" "<<tempvec[1]<<" "<<tempvec[2]<<endl;

      Tpowers.push_back(tempvec);
      numeric tempnum = ex_to<numeric>(temp.subs(lst(x==1.,y==1.,z==1.)).evalf());
      Tcoeff.push_back(tempnum.to_double());
    } 
    powers.push_back(Tpowers);
    coeff.push_back(Tcoeff);
    Tpowers.clear();
    Tcoeff.clear();
  }
#endif // FOUND_GINAC
}

