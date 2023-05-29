// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
// .NAME vtkShoeAttribute - Storage for a field to be interpolated over a vtkShoeBox mesh.
// .SECTION Description
// This class is a container for the coefficients of a polynomial interpolant
// used to represent a field over a mesh.
// Usually, the polynomial is an approximant to the solution of a finite element solver.
// In SHOE, this representation takes the form of values taken on at topological "corner"
// nodes -- which are always exactly interpolated -- and coefficients that are associated
// with edges, faces, or volumetric regions of the finite element.
// The values of these coefficients may or may not be exactly interpolated.
// In any event, all of the coefficients are sometimes called <em>degrees of
// freedom</em> (DOFs for short).
// DOFs frequently represent the magnitude of a spectral function, so each coefficient
// may be associated with a "mode shape" or "DOF mode."
// Because there may be multiple DOF modes for a single edge, face, or volume, these
// modes are grouped into "DOF nodes" -- the collection of coefficients whose mode shapes
// act upon a given edge, face, or volume.
//
// There are some sticky implementation details for attributes.
// First of all, the coefficients that are stored for edges, faces, and volumes are not
// all that is required to uniquely define an interpolant;
// the set of polynomials that the coefficients weight is also required.
// However, these polynomials vary by element species (see vtkShoeCellSpecies) and so
// the mesh connectivity and cell species must be known in order to reconstruct the
// interpolating polynomial.
// That means the range of values that an attribute may take on is not defined without
// a mesh that completely specifies the polynomial interpolants.
// For this reason, each attribute maintains a list of meshes that define polynomials
// using the attribute's coefficients.
// You are responsible for ensuring that all the meshes that define polynomials using
// an attribute's coefficients define identical polynomials from the same coefficients;
// otherwise, the concept of an attribute's range is not well-defined.
// 
// In spite of the problem outlined above, we allow multiple meshes to refer to the
// same attribute because many visualization operations would otherwise require exact
// copies of an attribute to work.
// Consider a thresholding operation that simply subsets a mesh;
// a cell is inserted into the subset when the values of an attribute fall within a given interval.
// Clearly, no coefficients of any of the thresholded mesh are any different from the
// initial mesh (its superset).
// So we allow both meshes to reference the same coefficients.
// Note that the thresholded mesh may not use all of the DOFs stored in its attributes.
// This is another sticky situation, because the range of values taken on by the attribute
// may have been calculated from the initial mesh, and this range will be larger than the
// thresholded mesh by definition.
// You should be aware that the range of an attribute is the convex hull of the range
// of all the polynomials defined over all the meshes that refer to that attribute,
// <strong>assuming</strong> that each mesh interprets the same coefficients in the same
// way. (Otherwise, all bets are off.)
//
// .SECTION See also
// vtkDataRecords
#ifndef __vtkShoeAttribute_h
#define __vtkShoeAttribute_h

#include "vtkGenericAttribute.h"
#include "vtkShoeEnums.h" // for vtkPolyInterpolant, vtkPolyProductSpace
#include "vtkShoeOrderTuple.h" // very lightweight wrapper for int[3]

#define VTK_RANGE_NONE 0
#define VTK_RANGE_SLOPPY 1
#define VTK_RANGE_PROPER 2
#define VTK_RANGE_TIGHT 3

class vtkShoeAttributeP;
class vtkDataArray;
class vtkBitArray;
class vtkDataRecords;
class vtkShoeBox;

class vtkShoeAttribute : public vtkGenericAttribute
{
public:
  vtkTypeRevisionMacro(vtkShoeAttribute,vtkGenericAttribute);
  virtual void PrintSelf(ostream& os, vtkIndent indent);
  static vtkShoeAttribute* New();

  // ================= Inherited methods =====================================
  virtual const char* GetName();
  virtual int GetNumberOfComponents();
  virtual int GetCentering();
  virtual int GetType();
  virtual int GetComponentType();
  virtual vtkIdType GetSize();
  virtual unsigned long GetActualMemorySize();
  virtual double* GetRange( int component );
  virtual void GetRange( int component, double range[2] );
  virtual double GetMaxNorm();
  virtual double* GetTuple( vtkGenericAdaptorCell* c );
  virtual void GetTuple( vtkGenericAdaptorCell* c, double* tuple );
  virtual double* GetTuple( vtkGenericCellIterator* c );
  virtual void GetTuple( vtkGenericCellIterator* c, double* tuple );
  virtual double* GetTuple( vtkGenericPointIterator* p );
  virtual void GetTuple( vtkGenericPointIterator* p, double* tuple );
  virtual void GetComponent( int i, vtkGenericCellIterator* c, double* values );
  virtual double GetComponent( int i, vtkGenericPointIterator* p );
  virtual void DeepCopy( vtkGenericAttribute* other );
  virtual void ShallowCopy( vtkGenericAttribute* other );

  // ================= New methods ===========================================
  virtual vtkIdType GetNumberOfPoints();
  virtual vtkIdType GetNumberOfDOFNodes();
  virtual void SetPointData( vtkDataArray* );
  virtual void SetDOFData( vtkDataRecords* );
  virtual void SetName( const char* );
  virtual void SetNumberOfComponents( int );
  virtual void SetNumberOfPoints( vtkIdType );
  virtual void SetNumberOfDOFNodes( vtkIdType );
  virtual void SetType( int );
  virtual void SetStorageType( int );
  // Description:
  // Set the range style and the range of all the components.
  // This method is only intended for use by readers.
  // This does <strong>not</strong> set the desired range type, only the style of the computed range.
  virtual void SetRange( double* range, int style );
  //BTX
  virtual void SetOrder( int cellSpeciesId, vtkShoeOrderTuple& order );
  int GetOrder( int cellSpeciesId, vtkShoeOrderTuple& orderOut );
  //ETX
  virtual void SetInterpolant( int cellSpeciesId, int interpolant );
  int GetInterpolant( int cellSpeciesId );
  virtual void SetProductSpace( int cellSpeciesId, int productSpace );
  int GetProductSpace( int cellSpeciesId );
  //BTX
  int GetCellTypeInfo( int cellSpeciesId, vtkShoeOrderTuple& orderOut, vtkPolyInterpolant& interpOut, vtkPolyProductSpace& psOut );
  void SetCriticalPointsDirty( int i );
  int GetCriticalPointsDirty() const;
  //ETX
  vtkDataArray* GetPointData();
  vtkDataRecords* GetDOFData();

  // Description:
  // Return an array of vtkDataRecords pointers, one for each component of the attribute.
  // These vtkDataRecords objects contain the critical points associated with each DOF node.
  // For an attribute with \f$N\f$ components, each vtkDataRecords object has \f$3+N\f$
  // components: three parametric coordinates (in storage order) locating the critical point
  // and \f$N\f$ values specifying the value of the attribute at the critical point.
  vtkDataRecords** GetCriticalPoints();

  // Description:
  // Return a bit array with one bit per DOF node.
  // The bit corresponding to a particular DOF is set when the critical points for that DOF
  // node have been computed.
  vtkBitArray* GetCriticalPointsComputed();

  // Description:
  // A debug subroutine that prints the critical points currently stored with an attribute.
  void PrintCriticalPoints( ostream& os );

  // Description:
  // Determine the way in which the bounds of an attribute are computed.
  // By default, \a RangeStyle is VTK_RANGE_SLOPPY.
  // The following enumerations can be used: <ul>
  // <li> VTK_RANGE_SLOPPY - Compute the range quickly, even if it means that
  //      some attribute values are not contained by the bounds.
  // <li> VTK_RANGE_PROPER - Compute the range quickly, but insure that the
  //      bounds are conservative. This can return very loose bounds.
  // <li> VTK_RANGE_TIGHT - Compute the exact range of the attribute. This
  //      can take longer to compute, but will be the exact range the
  //      attribute takes on over <strong>all</strong> meshes that reference
  //      this attribute.
  // </ul>
  vtkSetClampMacro(RangeStyle,int,1,3);
  vtkGetMacro(RangeStyle,int);
  void SetRangeStyleSloppy() { this->SetRangeStyle( VTK_RANGE_SLOPPY ); }
  void SetRangeStyleProper() { this->SetRangeStyle( VTK_RANGE_PROPER ); }
  void SetRangeStyleTight()  { this->SetRangeStyle( VTK_RANGE_TIGHT );  }

  virtual void InsertMeshReference( vtkShoeBox* );
  virtual void RemoveMeshReference( vtkShoeBox* );

protected:
  vtkShoeAttribute();
  virtual ~vtkShoeAttribute();

  //BTX
  friend class vtkShoeAttributeP;
  //ETX

  vtkShoeAttributeP* Data;
  int RangeStyle;

  // Description:
  // Mark the attribute to indicate that any critical points in storage are
  // invalid and should be recomputed.
  void MarkCriticalPointsDirty();

private:
  vtkShoeAttribute( const vtkShoeAttribute& ); // Not implemented.
  void operator=( const vtkShoeAttribute& ); // Not implemented.
};

#endif // __vtkShoeAttribute_h
