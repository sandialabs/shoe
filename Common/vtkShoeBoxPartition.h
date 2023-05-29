// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
// .NAME vtkShoeBoxPartition - A partition of a higher order mesh used for processing
// .SECTION Description
// This is an internal class used to store a partition of a mesh into regions
// that have no interior extrema (over a given set of attributes).

#ifndef __vtkShoeBoxPartition_h
#define __vtkShoeBoxPartition_h

#include "vtksnlCommonWin32Header.h" // for microsoft crud
#include "vtkObject.h"

#include <vtkstd/vector> // for std::vector
#include <vtkstd/list> // for std::list
#include <vtkstd/map> // for std::map
#include <vtkstd/set> // for std::set

class vtkShoeCell;
class vtkGenericAttributeCollection;

class VTK_SNL_COMMON_EXPORT vtkShoeBoxPartition : public vtkObject
{
public:
  vtkTypeRevisionMacro(vtkShoeBoxPartition,vtkObject);
  static vtkShoeBoxPartition* New();
  virtual void PrintSelf( ostream& os, vtkIndent indent );
  void Dump(); // for use in gdb

  // Description:
  // The tolerance on how close two distinct coordinates may be to each other
  // and still be considered identical.
  static const double Epsilon;

  void AddPointToFaceLoop( vtkShoeCell* containingCell, vtkIdType* conn, int dof, double* r );
  void TriangulateFaceLoop( vtkShoeCell* containingCell, vtkIdType* conn, int ne, int dof, double* r, int isBdy );

  // Description:
  // Generate the initial tessellation of a DOF, using either a C.P. (if there is one) or the element center.
  void StarInterior( vtkShoeBoxPartition* bdyTri, vtkShoeCell* sc, vtkIdType* conn, int ne, int nf, double* interiorPt );

  // Description:
  // Insert a point into the partition -- after the initial triangulation/tetrahedralization
  // has been performed.
  void InsertPointIntoTriangulation( vtkShoeCell* cell, vtkIdType* conn, int dof, double* r );
  void InsertPointIntoTetrahedralization( vtkShoeCell* cell, vtkIdType* conn, int dof, double* r );

  // Description:
  // Correct the edge topology of a tetrahedralization (iff face < 0) or a triangulation (otherwise).
  // Return true if the correction failed, false otherwise.
  bool CorrectEdgeTopology( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, vtkIdType d, int face );

  // Description:
  // Correct the face topology of a tetrahedralization.
  // Return true if the correction failed, false otherwise.
  bool CorrectFaceTopology( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, vtkIdType d );

  // Description:
  // Correct the topology of a triangulation.
  // Return true if the correction failed, false otherwise.
  bool CorrectTriTopology( vtkShoeCell* sc, vtkIdType* conn, int ne, int dof, vtkGenericAttributeCollection* kappa );

  // Description:
  // Correct the topology of a tetrahedralization.
  // Return true if the correction failed, false otherwise.
  bool CorrectTetTopology( vtkShoeCell* sc, vtkIdType* conn, int ne, int nf, vtkGenericAttributeCollection* kappa );

  // Description:
  // Return the first triangle in the partition of the given DOF node's domain.
  vtkIdType GetFirstTriangle( int dof ) const;

  // Description:
  // Return the number of triangles/tetrahedra in the partition of a DOF.
  vtkIdType GetNumberOfElements( int dof ) const;

  // Description:
  // Returns 1 if the co-facet indicated by \a frontOrBack of triangle \a t is on the hull, 0 otherwise.
  // The direction of the triangle does matter (an inwards-pointing triangle is not
  // considered to be on the hull even if its outwards face is on the hull).
  int IsTriangleOnHull( vtkIdType t, int frontOrBack );

  // Description:
  // Return the parameter-space coordinates of the tetrahedron described by
  // triangle \a t and \a frontOrBack, along with the other 3 triangular
  // co-facets (triangle number and "sidedness") that belong to the tetrahedron.
  void GetTetrahedron( vtkIdType t, int frontOrBack, double* coords, vtkIdType* tris, int* sides );

  //BTX
  typedef vtkTypeInt32 EleRef;

  class Point {
  public:
    double X[3]; // world coordinates
    double R[4]; // parametric coordinates
    // What's going to be stored in R:
    // for points on cell edges, R contains r_storage, r_projected,  (r,s)_current_cell_embedded_in_face
    // for points on cell faces, R contains (r,s)_storage, (r,s)_current_cell_projected_to_face
    // for points inside a cell, R contains (r,s,t)_storage
    vtkIdType DOF; // DOF node on which point lives (-1 for corners)
    vtkIdType Label; // A unique ID for the point, but not one used to index it here.
    EleRef Element; // Offset to an element that references this point 

    int Compare( const vtkShoeBoxPartition::Point& other ) const
      {
#if 0
      if ( DOF < other.DOF )
        return -1;
      else if ( DOF == other.DOF )
        {
        double d;
        d = X[0] - other.X[0];
        if ( d < -vtkShoeBoxPartition::Epsilon )
          return -1;
        else if ( fabs(d) < vtkShoeBoxPartition::Epsilon )
          {
          d = X[1] - other.X[1];
          if ( d < -vtkShoeBoxPartition::Epsilon )
            return -1;
          else if ( fabs(d) < vtkShoeBoxPartition::Epsilon )
            {
            d = X[2] - other.X[2];
            if ( d < -vtkShoeBoxPartition::Epsilon )
              return -1;
            else if ( fabs(d) < vtkShoeBoxPartition::Epsilon )
              return 0;
            }
          }
        }
      return 1;
#else
        double d;
        d = X[0] - other.X[0];
        if ( d < -vtkShoeBoxPartition::Epsilon )
          return -1;
        else if ( fabs(d) < vtkShoeBoxPartition::Epsilon )
          {
          d = X[1] - other.X[1];
          if ( d < -vtkShoeBoxPartition::Epsilon )
            return -1;
          else if ( fabs(d) < vtkShoeBoxPartition::Epsilon )
            {
            d = X[2] - other.X[2];
            if ( d < -vtkShoeBoxPartition::Epsilon )
              return -1;
            else if ( fabs(d) < vtkShoeBoxPartition::Epsilon )
              return 0;
            }
          }
        return 1;
#endif // 0
      }
    vtkShoeBoxPartition::Point& operator = ( const vtkShoeBoxPartition::Point& src )
      {
      if ( this == &src )
        return *this;
      for ( int i = 0; i < 3; ++i )
        {
        this->X[i] = src.X[i];
        this->R[i] = src.R[i];
        }
      this->R[3] = src.R[3];
      this->DOF = src.DOF;
      this->Label = src.Label;
      this->Element = src.Element;
      return *this;
      }
    bool operator < ( const vtkShoeBoxPartition::Point& other ) const
      { return this->Compare( other ) < 0; }
    bool operator > ( const vtkShoeBoxPartition::Point& other ) const
      { return this->Compare( other ) > 0; }
    bool operator == ( const vtkShoeBoxPartition::Point& other ) const
      { return this->Compare( other ) == 0; }
  };

  typedef vtkstd::map<vtkIdType,vtkstd::set<vtkShoeBoxPartition::Point> > PointSet;
  typedef vtkstd::set<vtkShoeBoxPartition::Point>::iterator PointRef;

  struct Triangle {
    PointRef EndPoints[3];
    EleRef Neighbors[6];
    int Flags;
  };
  //ETX

  // Description:
  // Set and get the value of \a MaximumSplitEdgeDepth.
  // By default, \a MaximumSplitEdgeDepth is 5.
  vtkGetMacro(MaximumSplitEdgeDepth, int);
  vtkSetMacro(MaximumSplitEdgeDepth, int);

  // Description:
  // Set and get the value of \a MaximumSplitFaceDepth.
  // By default, \a MaximumSplitFaceDepth is 5.
  vtkGetMacro(MaximumSplitFaceDepth, int);
  vtkSetMacro(MaximumSplitFaceDepth, int);

  // Description:
  // Set and get the value of \a NullTriangleThreshold.
  // By default, \a NullTriangleThreshold is 0.01.
  vtkGetMacro(NullTriangleThreshold, double);
  vtkSetMacro(NullTriangleThreshold, double);

  // Description:
  // Set and get the value of \a NullTriangleThreshold2.
  // By default, \a NullTriangleThreshold2 is 0.0001.
  vtkGetMacro(NullTriangleThreshold2, double);
  vtkSetMacro(NullTriangleThreshold2, double);

  // Description:
  // Set and get the value of \a NullTetThreshold.
  // By default, \a NullTetThreshold is 0.01.
  vtkGetMacro(NullTetThreshold, double);
  vtkSetMacro(NullTetThreshold, double);

protected:
  vtkShoeBoxPartition();
  virtual ~vtkShoeBoxPartition();

  //BTX
  typedef struct DofLocator {
    vtkIdType Start;
    vtkIdType NumberOfElements;
    DofLocator() { }
    DofLocator( const vtkIdType& s, const vtkIdType l ) : Start( s ), NumberOfElements( l ) { }
    DofLocator( const DofLocator& src ) : Start( src.Start ), NumberOfElements( src.NumberOfElements ) { }
  };

  typedef vtkstd::map<vtkIdType, DofLocator> DofOffsets;

  PointRef TriOrg( EleRef fv ) const;
  PointRef TriDest( EleRef fv ) const;
  EleRef   TriSym( EleRef fv ) const;
  EleRef   TriENext( EleRef fv ) const;
  EleRef   TriENext2( EleRef fv ) const;
  EleRef   TriFNext( EleRef fv ) const;
  EleRef   TriTurn( EleRef fv ) const;
  void     TriFSplice( EleRef fva, EleRef fvb );
  void     TriInvolutiveFMerge( EleRef fva, EleRef fvb );
  void     TriIdempotentFMerge( EleRef fva, EleRef fvb );
  void     TriFDel( EleRef fva );
  vtkIdType TriCreate( PointRef& a, PointRef& b, PointRef& c);
  EleRef   TriFindOrCreate( DofLocator& eleRef, PointRef a, PointRef b, PointRef c );
  void     TriMarkEdge( EleRef fv, bool mark );
  void     TriMarkEdges( EleRef fv, bool mark );
  void     TriMarkFace( EleRef fv, bool mark );
  bool     TriGetEdgeMark( EleRef fv );
  bool     TriGetFaceMark( EleRef fv );
  void     TriFaceIsOnHull( EleRef fv, bool onHull );
  bool     TriIsFaceOnHull( EleRef fv );

  void CopyPointsToDof( PointRef npts[3],
    vtkShoeCell* sc, vtkIdType* conn, int ne, int nf, int vdof,
    vtkShoeBoxPartition* bdyTri, vtkIdType facedof, int i );

  // Description:
  // Insert a point interior to a triangle in the tessellation of a face DOF.
  void StarTriangle( vtkIdType d, vtkIdType i,vtkShoeBoxPartition::PointRef p );

  // Description:
  // Insert a point interior to a tetrahedron in the tessellation of a body DOF.
  void StarTetrahedron( vtkIdType d, EleRef eref, vtkShoeBoxPartition::PointRef p );

  // Description:
  // Compute the critical points of the restrictions of the kappa fields of a cell to an edge.
  void ComputeEdgeCriticalPoints( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, EleRef eref, double* p1, double* p2, vtkstd::vector<vtkstd::pair<double[3],EleRef> >& edgeCPs );

  // Description:
  // Compute the critical points of the restrictions of the kappa fields of a cell to the edges of a face dof triangulation.
  void ComputeTriangulationEdgeCriticalPoints( vtkShoeCell* sc, int face, vtkGenericAttributeCollection* kappa, vtkIdType d, vtkstd::vector<vtkstd::pair<double[3],EleRef> >& intEdgeCPs );

  // Description:
  // Compute the critical points of the restrictions of the kappa fields of a cell to the interior edges of a body dof tetrahedralization, and stores them in a vector.
  void ComputeTetrahedralizationEdgeCriticalPoints( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, vtkIdType d, vtkstd::vector<vtkstd::pair<double[3],EleRef> >& intEdgeCPs );

  // Description:
  // Compute the critical points of the restrictions of the kappa fields of a cell to the interior faces of a body dof tetrahedralization, and stores them in a vector.
  void ComputeTetrahedralizationFaceCriticalPoints( vtkShoeCell* sc, vtkGenericAttributeCollection* kappa, vtkIdType d, vtkstd::vector<vtkstd::pair<double[3],EleRef> >& intFaceCPs );

  // Description:
  // Split an edge interior to the triangulation of a face DOF using a given point.
  void SplitBoundaryEdge( vtkIdType d, EleRef eref, vtkShoeBoxPartition::PointRef p );

  // Description:
  // Split an edge interior to the tetrahedralization of a body DOF using a given point.
  void SplitInteriorEdge( vtkIdType d, EleRef eref, vtkShoeBoxPartition::PointRef p );

  // Description:
  // Split an edge interior to the triangulation of a body (iff face < 0) or face (otherwise) DOF 
  // and create the point that splits the edge.
  void SplitEdge( vtkShoeCell* sc, vtkIdType d, EleRef eref, double r[3], int face );

  int UnSanity( int dof );
  int MadITellYou();

  // Description:
  // Check the orientation of the tetrahedralization attached to a DOF.
  // Returns true if all tetrahedra are positively oriented, false otherwise.
  bool CheckTetOrientations( vtkIdType d );

  // Description:
  // Edge-facet references.
  static int VO[];
  static int VE[];

  // Description:
  // Points in the partition.
  vtkShoeBoxPartition::PointSet* Points;
  int OwnPoints;
  vtkIdType NextPointLabel;

  // Description:
  // Some temporary storage used to hold points on the edge loop
  // bounding a face before we triangulate it. Points are inserted
  // with AddPointOnFaceLoop() and then consumed by calling
  // TriangulateFaceLoop().
  vtkstd::vector<vtkShoeBoxPartition::PointRef> EdgePoints;

  // Description:
  // Elements in the partition.
  vtkstd::vector<Triangle> Triangles;

  // Description:
  // Offset (into Elements) and size (number of elements) for a given DOF ID.
  DofOffsets Offsets;

  // Description:
  // Maximum number of depth levels when iteratively correcting the edge topology w.r.t. the given kappa attributes.
  int MaximumSplitEdgeDepth;

  // Description:
  // Maximum number of depth levels when iteratively correcting the face topology w.r.t. the given kappa attributes.
  int MaximumSplitFaceDepth;

  // Description:
  // Threshold for the relative area for a subtriangle to its parent below which the former is declared to be null.
  double NullTriangleThreshold;

  // Description:
  // Threshold for the squared relative area for a subtriangle to its parent below which the former is declared to be null.
  double NullTriangleThreshold2;

  // Description:
  // Threshold for the relative volume for a subtetrahedron to its parent below which the former is declared to be null.
  double NullTetThreshold;
  //ETX


private:
  void operator = ( const vtkShoeBoxPartition& ); // Not implemented.
  vtkShoeBoxPartition( const vtkShoeBoxPartition& ); // Not implemented.
};

//BTX
inline vtkIdType vtkShoeBoxPartition::GetFirstTriangle( int dof ) const
{
  DofOffsets::const_iterator it;
  it = this->Offsets.find( dof );
  if ( it == this->Offsets.end() )
    return 0;
  return it->second.Start;
}

inline vtkIdType vtkShoeBoxPartition::GetNumberOfElements( int dof ) const
{
  DofOffsets::const_iterator it;
  it = this->Offsets.find( dof );
  if ( it == this->Offsets.end() )
    return 0;
  return it->second.NumberOfElements;
}
//ETX

#endif // __vtkShoeBoxPartition_h
