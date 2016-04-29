// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
// .NAME vtkShoeCellMetaData - record for data specific to a genus of cells
// .SECTION Description
// One instance of this class is created for each genus (shape) of cell.
#ifndef __vtkShoeCellMetaData_h
#define __vtkShoeCellMetaData_h

#include <vtkstd/set>
#include "vtksnlCommonWin32Header.h" // for Windows filth
#include "vtkShoeEnums.h" // for NumberOfCellShapes
#include "vtkShoeOrderTuple.h" // for PrepareForOrder arg
#include "vtkObject.h"

class vtkShoeCellIterator;
class vtkDoubleArray;
class vtkPolynomialSystem;
class vtkDataRecords;

class VTK_SNL_COMMON_EXPORT vtkShoeCellMetaData : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkShoeCellMetaData,vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkShoeCellMetaData* New();

  // Description:
  // Cell type as a unique integer
  int Type;

  // Description:
  // 0 - 3 (embedding dimension)
  int Dimension;

  // Description:
  // Number of 0-, 1-, or 2-D boundaries of this cell
  int NumberOfBoundaries[3];

  // Description:
  // sum of NumberOfBoundaries[1 - 2] + 1 (for cell DOF node)
  int NumberOfDOFNodes;

  // Description:
  // Number of topological corner points on a given 2-D boundary face (length NumberOfBoundaries[2])
  int* NumberOfVerticesOnFace;

  // Description:
  // Array of offsets into cell connectivity describing corner vertices bounding a given edge
  int** EdgeArray;

  // Description:
  // Array of offsets into cell connectivity describing corner vertices bounding a given face
  int** FaceArray;

  // Description:
  // FacesOnEdge[ee][ff] is the ff-th face that is bounded by edge ee
  int** FacesOnEdge;

  // Description:
  // EdgesOnFace[ff][ee] is the ee-th edge that bounds face ff
  int** EdgesOnFace;


  // Description:
  // Parametric coordinates of the center of the cell
  double* ParametricCenter;

  // Description:
  // An array of the parametric coordinates of the topological corner points (length 3*NumberOfBoundaries[0])
  double* ParametricCorners;

  // Description:
  // An array used by vtkShoeCell::RestrictToBoundary to transform (r,s,t)
  // coordinates into a proper (u,v) parameterization of the face.
  // There are 9 numbers stored for each face:
  // <ul>
  // <li> the first 3 are the parametric center of the face
  // <li> the next 3 are the "u" direction in the face's local coordinates
  // <li> the last 3 are the "v" direction in the face's local coordinates
  // </ul>
  // An example: for face 0 of a hexahedron (defined by the plane r = -1),
  // we store ( -1, 0, 0,  0, 1, 0,  0, 0, 1 ).
  double** BoundaryFaceTransforms;

  // Description:
  // An array used by vtkShoeCell::RestrictToBoundary to transform (r,s,t)
  // coordinates into a proper (u) parameterization of the edge.
  // There are 6 numbers stored for each edge:
  // <ul>
  // <li> the first 3 are the parametric center of the edge
  // <li> the next 3 are the "u" direction in the edge's local coordinates
  // </ul>
  // An example: for edge 0 of a hexahedron (defined by the planes s = t = -1),
  // we store ( 0, -1, -1,  1, 0, 0 ).
  double** BoundaryEdgeTransforms;

  double GetParametricDistance( double x[3] ) const;
  int FindClosestBoundary( double pcoords[3], vtkShoeCellIterator*& boundary ) const;
  void Derivatives( double* pcoords, vtkDoubleArray* cellData, double* derivs ) const;

  /** Is coedge \a e on face \a f reversed with respect to the underlying edge of the cell?
    * Since cells have all their edges oriented so that the projection onto the cell's coordinate
    * axes is positive, half the edges of a face loop will be reversed (for the quadrilateral and
    * hexahedral elements; life is different for simplices).
    * Note that \a e < \a NumberOfVerticesOnFace[f] == \a NumberOfEdgesOnFace[f].
    * Computed by computing \a this->EdgeArray[this->EdgesOnFace[f][e]][0] == \a this->FaceArray[f][e];
    * if they are equal, the coedge is positively oriented and 0 is returned; otherwise 1 is returned.
    */
  int FaceCoedgeReversed( int f, int e ) const;

  //BTX
  /** Return the value of each shape function at the given parametric coordinates (\a pcoords).
    */
  void (*ShapeFunctions)( const int* order, const double* pcoords, double* shape );
  /** Return the derivatives of each shape function at the given parametric coordinates (\a pcoords).
    */
  void (*ShapeFunctionDerivatives)( const int* order, const double* pcoords, double* sderivs );
  /** Construct a polynomial system of equations representing the gradient of a scalar field.
    *
    * For 3-D elements, this will be a set of 3 equations, 2 for 2-D cells, and 1 for 1-D cells.
    *
    * Note that this <strong>only</strong> accepts <em>scalar</em> (i.e., a single real component) fields,
    * not vector or tensor fields. If you have a vector or tensor field, you are expected to call
    * this function once for each component.
    */
  void (*SymbolicFieldGradient)( const int* order, const double* permuted_field_vals, vtkPolynomialSystem* grad );

  /** Prepare to handle a given order of interpolant.
    * This is a pointer to a function that will be called whenever a new species of cell
    * is created.
    * It will be called once for each specified attribute, and the attribute's
    * order will be passed in the \a order parameter.
    * This is used to allocate and evaluate constants used in polynomials of
    * the given order.
    *
    * As an example, see how the Lagrange interpolation routines use this to
    * compute the coefficients that convert a univariate Lagrange basis into
    * the power basis.
    */
  void (*PrepareForOrder)( const vtkShoeOrderTuple& );

  /** Embed edge parametric coordinates into element parametric coordinates.
   * @param edge - The edge of the element along which a coordinate is given.
   * @param coords - Initially, coords[0] contains the coordinate along the edge. 
   *        Upon completion, all three parametric coordinates will be set to
   *        valid values (coord[2] will be 0.0 for 2-dimensional elements).
   */
  void (*EmbedEdgeCoords)( int edge, double* coords ); 

  /** Embed edge parametric coordinates into the parametric coordinates of one of an element's faces.
   * This is defined only for 3-dimensional elements; 2-dimensional elements should use EmbedEdgeCoords
   * instead.
   * @param edge - The edge of the element along which a coordinate is given.
   * @param face - The face of the element into which the edge coordinate should be embedded.
   *               Note that \a edge must be a bounding edge of \a face.
   * @param coords - Initially, coords[0] contains the coordinate along the edge. 
   *        Upon completion, the first two parametric coordinates will be set to
   *        valid values (coord[2] is unused).
   *        Note that the result is undefined if you ask for coordinates along an edge
   *        to be embedded into a face the edge does not bound.
   */
  void (*EmbedEdgeCoordsInFace)( int edge, int face, double* coords ); 

  /** Embed face parametric coordinates into element parametric coordinates
   * @param coords - Initially, coords[0] contains the coordinate along the edge. 
   *        Upon completion, all three parametric coordinates will be set to
   *        valid values (coord[2] will be 0.0 for 2-dimensional elements).
   */
  void (*EmbedFaceCoords)( int face, double* coords );

  /** Project parameteric coordinates in a cell's local coordinate space (3D) to an edge of the cell (1D).
    * This routine uses the cell's local coordinate system projected to the edge of interest.
    * Output edge coordinates are again in cell-local coordinates, not face-loop(coedge) or storage coordinates.
    * @param tuple - the coordinates; the first entry of the array is overwritten with the output coordinate.
    */
  void (*ProjectCoordsToEdge)( int edge, double* tuple );

  /** Project parameteric coordinates in a cell's local coordinate space (3D) to a face of the cell (2D).
    * This routine uses the cell's local coordinate system projected to the face of interest.
    * @param tuple - the coordinates; the first 2 entries of the array are overwritten with the output coordinates.
    */
  void (*ProjectCoordsToFace)( int face, double* tuple );

  /** Project parameteric coordinates in a cell's face's local coordinate space (2D) to an edge of the cell face (1D).
    * This routine uses the cell's local coordinate system projected to the face of interest, not storage coordinates.
    * Output edge coordinates are again in cell-local coordinates, not face-loop(coedge) or storage coordinates.
    * @param tuple - the face coordinates; the first entry of the array is overwritten with the output coordinate.
    */
  void (*ProjectFaceCoordsToEdge)( int face, int edge, double* tuple );

  /** Transform parameteric coordinates in a face's local coordinate space (2D) to storage coordinates (2D).
    * This routine uses the face's permutation to transform its local coordinate system into storage coordinates.
    * @param tuple - the coordinates; the first 2 entries of the array are overwritten with the output coordinates.
    */
  void (*TransformFaceCoordsToStorage)( uint32_t perm, double* tuple );

  /** Decide whether a given parameter-space tuple is inside, outside, or on the boundary of an element or its boundary.
    *
    * Returns -1, 0, or 1 depending on whether the tuple is outside, on the boundary of, or inside the region of interest, respectively.
    * The region of interest is specified by \a bdy.
    * If \a bdy < 0 or \a bdy >= NumberOfDOFNodes, then the region of interest is the entire finite element.
    * Otherwise, the region of interest is the edge or face specified by \a bdy.
    * For example, a hexahedron will take \a bdy = 12 to be the first face (vertices 0-4-7-3, with plane equation r = -1).
    * Since the dimension of region of interest is the number of coordinates required to specify a point in that region, only
    * that many coordinates of \a param will be examined.
    * Staying with our example of the hexahedron, only the first 2 coordinates will be examined for \a bdy 12 and they
    * correspond to (s,t).
    */
  int (*IsParameterInDomain)( double* param, int bdy );
  //ETX

  static int InitializeShoeCells();
  static const vtkShoeCellMetaData* ByShape[shoe::NumberOfCellShapes+1];

  // Description:
  // Call this static method to release all of the metadata records representing
  // the different element types.
  static void Shutdown();

  // Description:
  // Function pointer used by some instances that need no implementation.
  static void EmptyPrepareForOrder( const vtkShoeOrderTuple& );
  static void EmptyStoreRelevantRoots( vtkPolynomialSystem*, vtkIdType*, int, int, vtkDataRecords* );

protected:
  vtkShoeCellMetaData();
  ~vtkShoeCellMetaData();

private:
  vtkShoeCellMetaData( const vtkShoeCellMetaData& ); // Not implemented.
  void operator = ( const vtkShoeCellMetaData& ); // Not implemented.

  //BTX
  static vtkstd::set<vtkShoeCellMetaData*> All;
  //ETX
};

//BTX
inline int vtkShoeCellMetaData::FaceCoedgeReversed( int f, int e ) const
{
  return this->EdgeArray[this->EdgesOnFace[f][e]][0] == this->FaceArray[f][e] ? 0 : 1;
}
//ETX
#endif // __vtkShoeCellMetaData_h
