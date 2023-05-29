/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <math.h>

#include <iostream>
#include <algorithm>
#include <iterator>

#ifdef __APPLE__
// Yeah, this sucks, but std::isinf() is broken.
// Thanks, Mr. Jobs, I guess I didn't want IEEE floats after all.
#  ifndef isinf
#    define  isinf( x )         ( ( sizeof ( x ) == sizeof(double) ) ?           \
	                              __isinfd ( x ) :                                 \
																 ( sizeof ( x ) == sizeof( float) ) ?            \
																 __isinff ( x ) :                                \
																 __isinf  ( x ) )
#  endif // isinf
#endif // __APPLE__

#include <vtkObjectFactory.h>
#include <vtkMath.h>

#include <vtkShoeMesh.h>
#include <vtkCellOps.h>
#include <vtkAdaptiveTessellator.h>
#include <vtkShoeMeshFieldValueSubdivision.h>

using vtkstd::ostream_iterator;
using vtkstd::copy;

vtkCxxRevisionMacro(vtkShoeMeshFieldValueSubdivision,"$Revision: 8644 $");
vtkStandardNewMacro(vtkShoeMeshFieldValueSubdivision);

void vtkShoeMeshFieldValueSubdivision::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
	os << indent << "FieldCriteria:        " << this->FieldCriteria << endl;
	os << indent << "AllowableL2Error2:    [ ";
	copy( this->AllowableL2Error2, this->AllowableL2Error2 + this->NumberOfFields, ostream_iterator<double>( os, " " ) );
	os << " ]" << endl;
	os << indent << "AllowableChordError2: " << this->AllowableChordError2 << endl;
}

vtkShoeMeshFieldValueSubdivision::vtkShoeMeshFieldValueSubdivision()
{
  this->FieldCriteria = 0;

  if ( this->NumberOfFields )
    this->AllowableL2Error2 = new double[ this->NumberOfFields ];
  else
    this->AllowableL2Error2 = 0;

  this->AllowableChordError2 = .01; // in world coordinates
}

vtkShoeMeshFieldValueSubdivision::~vtkShoeMeshFieldValueSubdivision()
{
  if ( this->AllowableL2Error2 )
    delete [] this->AllowableL2Error2;
}

#if 0
void vtkShoeMeshFieldValueSubdivision::PassFields( vtkDataSetAttributes* outarrays, vtkAdaptiveTessellator* t )
{
  this->Superclass::PassFields( outarrays, t );
}
#endif // 0

void vtkShoeMeshFieldValueSubdivision::ResizeL2Array( int old_length )
{
  if ( this->NumberOfFields == 0 )
    {
    if ( this->AllowableL2Error2 )
      delete [] this->AllowableL2Error2;
    this->AllowableL2Error2 = 0;
    }
  else
    {
      double* allocated = new double [ this->NumberOfFields ];
      int len = old_length > this->NumberOfFields ? this->NumberOfFields : old_length;
      for ( int i=0; i<len; ++i )
        allocated[i] = this->AllowableL2Error2[i];
      delete [] this->AllowableL2Error2;
      this->AllowableL2Error2 = allocated;
    }
}

int vtkShoeMeshFieldValueSubdivision::PassField( int sourceId, int sourceSize, vtkAdaptiveTessellator* t )
{
  int nf = this->NumberOfFields;
  int val = this->Superclass::PassField( sourceId, sourceSize, t );
  if ( nf != this->NumberOfFields )
    this->ResizeL2Array( nf );

  if ( val >= 0 )
    this->AllowableL2Error2[ val ] = vtkMath::Inf(); // To infinity.. and beyond!

  return val;
}

void vtkShoeMeshFieldValueSubdivision::ResetFieldList()
{
  this->Superclass::ResetFieldList();

  if ( this->AllowableL2Error2 )
    {
    delete [] this->AllowableL2Error2;
    this->AllowableL2Error2 = 0;
    }

  this->FieldCriteria = 0;
}

bool vtkShoeMeshFieldValueSubdivision::DontPassField( int sourceId, vtkAdaptiveTessellator* t )
{
  int id = this->GetOutputField( sourceId );
  if ( id >= 0 )
    {
    // Update the FieldCriteria
    int missing_bit = 1 << id; // the bit to be excised
    int low_mask = missing_bit - 1 < 0 ? 0 : missing_bit -1;
    int high_mask = (~low_mask) & (~1);
    int low_part = this->FieldCriteria & low_mask;
    int high_part = this->FieldCriteria & high_mask;
    high_part >>= 1;
    this->FieldCriteria = high_part | low_part;

    // Update AllowableL2Error2
    for ( int j=id; j<this->NumberOfFields-1; ++j )
      this->AllowableL2Error2[j] = this->AllowableL2Error2[j+1];
    }

  // Cut things down to size
  int nf = this->NumberOfFields;
  bool val = this->Superclass::DontPassField( sourceId, t );
  if ( val )
    this->ResizeL2Array( nf );

  return val;
}

void vtkShoeMeshFieldValueSubdivision::SetAllowableChordError2( double e )
{
  if ( e == this->AllowableChordError2 )
    return;

  // NOTE: This is different than vtkShoeMeshViewDependentSubdivision!!!
  // We allow e == 0, so that you can subdivide purely on field values.
  if ( e >= 0. )
    {
    this->AllowableChordError2 = e;
    }
  else
    {
      vtkWarningMacro( "Negative values (" << e << ") for AllowableChordError2 are forbidden." );
      return;
    }
  this->Modified();
}

void vtkShoeMeshFieldValueSubdivision::SetAllowableL2Error2( int sourceFieldId, double e )
{
  int id = this->GetOutputField( sourceFieldId );
  if ( id < 0 )
    {
      vtkWarningMacro( "Field " << sourceFieldId << " is not currently in the output mesh." );
      return;
    }

  if ( e == this->AllowableL2Error2[id] )
    return;

  // NOTE: This is different than vtkShoeMeshViewDependentSubdivision!!!
  // We allow e == 0, so that you can subdivide purely on field values.
  if ( e >= 0. )
    {
    this->AllowableL2Error2[id] = e;

    if ( isinf( e ) )
        this->FieldCriteria &= ~( 1<<id );
    else
        this->FieldCriteria |=  ( 1<<id );
    }
  else
    {
      vtkWarningMacro( "Negative values (" << e << ") for AllowableL2Error2 are forbidden." );
      return;
    }
  this->Modified();
}

double vtkShoeMeshFieldValueSubdivision::GetAllowableL2Error2( int sourceFieldId ) const
{
	int id = this->GetOutputField( sourceFieldId );
	if ( id < 0 )
		{

		// Can't call warning macro since it uses non-const methods. Aaaagh!
		//vtkWarningMacro( "Field " << sourceFieldId << " is not currently in the output mesh." );
		return -1.0;
		}

	return this->AllowableL2Error2[ id ];
}

bool vtkShoeMeshFieldValueSubdivision::EvaluateEdge( const double* p0, double* p1, const double* p2, int field_start )
{
  double real_p1[3];

  double dist=0.;
  double tmp;
  const vtkCellOps* ops = this->CurrentCell.GetCellOps();

  ops->EvaluateGeometry( real_p1, this->CurrentCell, p1+3 );
  for ( int i=0; i<3; ++i )
    {
    tmp = real_p1[i] - p1[i];
    dist += tmp*tmp;
    }

  bool rval = dist > this->AllowableChordError2;
  if ( rval )
    {
    for ( int j=0; j<3; ++j )
      p1[j] = real_p1[j];

		// Evaluate field values at the midpoint
		this->EvaluateFields( p1, field_start );
    }

  // If the geometry test requires a subdivision, we have already
  // evaluated the fields and there's no need to do any tests
  if ( (!rval) && this->FieldCriteria )
    {
    // Otherwise, we need to evaluate the fields but save the
    // linearly interpolated values, so we need some temporary
    // storage.
    double real_pf[6+vtkAdaptiveTessellator::MaxFieldSize];
    copy( p1, p1 + field_start, real_pf );
    this->EvaluateFields( real_pf, field_start );

		rval = this->FixedFieldErrorEval( p0, p1, real_pf, p2, field_start, this->FieldCriteria, this->AllowableL2Error2 );
		if ( rval )
			{
			copy( real_pf+field_start, real_pf+field_start+this->FieldOffsets[this->NumberOfFields], p1+field_start );
			}
		}

	return rval;
}

