/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2005-04-19 09:22:52 -0700 (Tue, 19 Apr 2005) $
  Version:   $Revision: 4272 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef shoe_enums_h
#define shoe_enums_h

/*!\brief Specify how a mesh iterator should traverse cells in the mesh.
 */
enum MeshTraversalStyle
{
	OfSameType,  //!< Make multiple passes through the mesh, covering only cells of a single type per pass
	OfSameDefn,  //!< Traverse cells, visited only those with same CellDefinition (possible with a mask)
	CustomOrder, //!< User defined traversal using MeshTraversalMask defined below
	MeshOrder,   //!< Traverse cells one by one as they occur in the mesh
  Selection,   //!< Traverse cells in a specified list
  Boundary,    //!< Traverse the n-dimensional boundary of each cell of dimension n+1 or higher in the mesh
  ExtBoundary, //!< Traverse the n-dimensional boundary of each cell of dimension n+1 or higher that resides on the exterior of the mesh
  CellBoundary //!< Traverse the n-dimensional boundary of a single cell (of dimension n+1 or higher)
};

/*!\brief Specify special mesh traversal criteria
 */
enum MeshTraversalMask
{
	ShapeMask = 0x1,         //!< Visit all cells with same shape
	InterpMask = 0x2,        //!< Visit all cells using same interpolation functions
	ProdSpaceMask = 0x4,     //!< Visit all cells using same product space
	TypeMask = 0x7,          //!< Visit all cells that are of same type (all of the above)
	FunctionOrderMask = 0x8, //!< Visit all cells that have all functions of the same respective orders
	GeomOrderMask = 0x10,    //!< Visit all cells with same order for geometry shape functions
	DefnMask = 0x1f          //!< Visit all cells that have the same definition (all of the above)
};

// CellShape must be in a private namespace on MacOS X since Carbon
// defines struct Point and Qt typedefs it to Point. Please don't
// check in lameness that renames Point to POINT (yelling, how rude!)
// or Node (not true!) or vtkSnlPointEnum (simply boorish). Use a
// namespace for once in your life.
namespace shoe {

/*!\brief Specify the parametric domain of a cell
 */
enum CellShape
{
	Point,               //!< A 0-D, infinitesimal point
	Curve,               //!< A 1-D, infinitesimally thin cell
	Rod,                 //!< A solid 1-D cell with finite thickness
	Tube,                //!< A hollow 1-D cell with a shell of finite thickness
	Triangle,            //!< A 2-D, infinitesimally thin triangular patch
	Quadrilateral,       //!< A 2-D, infinitesimally thin quadrilateral patch
	TriangleShell,       //!< A solid 2-D triangular patch with a finite thickness
	QuadrilateralShell,  //!< A solid 2-D quadrilateral patch with a finite thickness
	Tetrahedron,         //!< A solid 3-D tetrahedral patch
	Hexahedron,          //!< A solid 3-D hexahedral patch
	Wedge,               //!< A solid 3-D wedge-shaped patch
	Pyramid,             //!< A solid 3-D pyramidal patch
  NumberOfCellShapes   // THIS MUST BE THE LAST ENTRY
};

} // shoe namespace

/*!\brief Families of orthogonal polynomials and their interpolation schemes.
 * All polynomial interpolants operate in the range [-1,1].
 */
enum vtkPolyInterpolant
{
	Lagrange, //!< Lagrange interpolants
	Legendre, //!< Legendre polynomials
  Hermite,  //!< Hermitian polynomials
	BSpline   //!< B-Spline interpolants
};

/*!\brief The shape of the polynomial product space.
 *
 * This enumerant specifies how a univariate polynomial family specified
 * by PolyInterpolant is to be combined into a multivariate polynomial
 * basis. Call the univariate polynomial family \f$\phi_q(u)\f$. Then basis
 * polynomials are formed (in 3-D) by taking products of these polynomials
 * in each parameter-space coordinate: \f$\phi_{q_r}(r)\phi_{q_s}(s)\phi_{q_t}(t)\f$,
 * where \f$q_r\f$, \f$q_s\f$, and \f$q_t\in Z^+\f$.
 * Each product space specifies which combinations of \f$(q_r,q_s,q_t)\f$
 * are chosen for inclusion in the multivariate polynomial.
 */
enum vtkPolyProductSpace
{
	/*! A tensor product of univariate polynomials
	 *  Given integers \f$k,\ell,m > 0,\f$ choose \f$q_r\in{0,1,\ldots,k}, q_s\in{0,1,\ldots,\ell}, q_t\in{0,1,\ldots,m}\f$
	 */
	Tensor,
	/*! Products of polynomials whose total degree is bounded.
	 * Given \f$k > 0\f$, choose \f$q_r,q_s,q_t\in{0,1,\ldots,k}\f$ st \f$q_r+q_s+q_t\leq k\f$.
	 * The side modes \f$(k,1,1), (1,k,1), (1,1,k)\f$ are added.
	 */
	MaxTotalDegree,
	/*! A subset of \a MaxTotalDegree.
	 * Given integers \f$k,\ell,m > 0,\f$ assume without loss of generality that \f$\ell=max(k,\ell,m)\f$.
	 * Then no shape function in \f$r\f$ or \f$t\f$ is used if its degree is larger than \f$k\f$ or \f$m\f$, respectively.
	 * This can be used to create a cell of nominal order \f$k\f$ while truncating higher order terms on boundaries in the \f$r\f$ 
	 * and \f$t\f$ directions where boundary constraints prevent the full interpolant.
	 */
	TruncatedTotalDegree
};

#endif // shoe_enums_h
