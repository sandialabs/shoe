// Copyright 2012 Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the
// U.S. Government. Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that this Notice and any
// statement of authorship are reproduced on all copies.
#include "vtkShoeAttribute.h"

#include "vtkShoeBox.h"
#include "vtkShoeCell.h"
#include "vtkShoeCellSpecies.h"
#include "vtkShoeCellMetaData.h"
#include "vtkShoeCellIterator.h"
#include "vtkShoePointIterator.h"
#include "vtkDataRecords.h"
#include "vtkDataRecordsIterator.h"

#include "vtkObjectFactory.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkBitArray.h"
#include "vtkDataSetAttributes.h"

#include <vtkstd/vector>
#include <vtkstd/string>
#include <vtkstd/set>
#include <vtkstd/map>

#include <assert.h>

vtkCxxRevisionMacro(vtkShoeAttribute,"$Revision:");
vtkStandardNewMacro(vtkShoeAttribute);

enum vtkRangeStyle {
  vtkNoRange,     //!< Range has not been computed
  vtkSloppyRange, //!< Range may not be conservative... function may take on values outside the range.
  vtkProperRange, //!< Range serves as proper bounds, but may not be tight.
  vtkTightRange   //!< Range is both proper and tight.
};

typedef struct _PerCellTypeRec
{
  vtkShoeOrderTuple Order;
  vtkPolyInterpolant Interpolant;
  vtkPolyProductSpace ProductSpace;
} PerCellTypeRec;

// ============================================================================

class vtkShoeAttributeP
{
public:
  vtkShoeAttributeP( vtkShoeAttribute* parent );
  ~vtkShoeAttributeP();

  void SetNumberOfPoints( vtkIdType );
  void SetNumberOfDOFNodes( vtkIdType );

  void SetNumberOfComponents( int );
  int GetNumberOfComponents();
  int GetComponentType();
  void SetStorageType( vtkShoeAttribute* parent, int );

  double* GetRange( vtkRangeStyle rangeStyle, int component );
  void ResetRange();
  double GetMaxNorm();

  vtkIdType GetSize();
  unsigned long GetActualMemorySize();

  void GetTuple( vtkShoeCell*, double* );
  void GetTuple( vtkShoeCellIterator*, double* );
  void GetTuple( vtkShoePointIterator*, double* );
  void GetComponent( int, vtkShoeCellIterator*, double* );
  double GetComponent( int, vtkShoePointIterator* );
  int SetOrder( int ID, const vtkShoeOrderTuple& );
  int GetOrder( int ID, vtkShoeOrderTuple& );
  int SetInterpolant( int ID, vtkPolyInterpolant );
  vtkPolyInterpolant GetInterpolant( int ID );
  int SetProductSpace( int ID, vtkPolyProductSpace );
  vtkPolyProductSpace GetProductSpace( int ID );
  int GetCellTypeInfo( int cellSpeciesId, vtkShoeOrderTuple& orderOut, vtkPolyInterpolant& interpOut, vtkPolyProductSpace& psOut );
  vtkIdType GetNumberOfPoints();
  vtkIdType GetNumberOfDOFNodes();
  void UpdateRange( vtkRangeStyle rangeStyle, int component );

  void DeepCopy( vtkShoeAttributeP* other );
  void ShallowCopy( vtkShoeAttributeP* other );

  // basic bookkeeping information
  vtkstd::string Name;
  int Type; // one of vtkDataSetAttributes::AttributeTypes or -1
  int Centering;
  vtkShoeAttribute* Parent;
  vtkstd::set<vtkShoeBox*> References; // meshes referencing this attribute

  // coefficient storage
  vtkDataArray* PointData; // Corner values
  vtkDataRecords* DOFData; // DOF node values

  // polynomial map specification: given a cell type (int), return how it is interpolated over this attribute
  vtkstd::map<int,PerCellTypeRec> MapSpec;

  // critical point storage
  vtkstd::vector<vtkDataRecords*> CriticalPoints; // Pre-computed critical points (if they exist) per component
  vtkBitArray* CriticalPointsComputed; // Bit flag denoting whether a given DOF's critical points have been computed.
  int CriticalPointsDirty; // Do we need to update the critical points? 0 = no, 1 = only for new values, 2 = recompute

  // range data for each attribute component
  vtkstd::vector<vtkstd::pair<double,double> > Range; // Range of each component... may be imprecise
  vtkRangeStyle RangeStyle; // May take on any value of vtkRangeStyle
  vtkstd::vector<double> TupleCache;

  // used to assign a unique name when none is specified
  static unsigned UniqueId;
};

// =BPM========================================================================

unsigned vtkShoeAttributeP::UniqueId = 0;

vtkShoeAttributeP::vtkShoeAttributeP( vtkShoeAttribute* parent )
{
  char DefaultName[48];
  sprintf( DefaultName, "UnnamedAttribute%02u", vtkShoeAttributeP::UniqueId );
  this->Name = vtkstd::string( DefaultName );

  //this->Centering = vtkBoundaryCentered;
  this->Centering = vtkPointCentered;
  this->Parent = parent;
  this->PointData = vtkDoubleArray::New();
  this->DOFData = vtkDataRecords::Create( VTK_DOUBLE );
  this->CriticalPoints.push_back( vtkDataRecords::Create( VTK_DOUBLE ) );
  this->CriticalPoints[0]->SetNumberOfComponents(4);
  this->CriticalPointsComputed = vtkBitArray::New();
  this->CriticalPointsComputed->SetNumberOfComponents(1);
  this->CriticalPointsDirty = 2;
  this->SetNumberOfComponents( 1 );
  this->Type = vtkDataSetAttributes::SCALARS;
  this->ResetRange();
}

vtkShoeAttributeP::~vtkShoeAttributeP()
{
  this->CriticalPointsComputed->Delete();
  for ( int i = 0; i < this->GetNumberOfComponents(); ++i )
    {
    this->CriticalPoints[i]->Delete();
    }
  this->DOFData->Delete();
  // Warning: GetNumberOfComponents() above uses PointData->GetNumberOfComponents.
  this->PointData->Delete();
}

vtkIdType vtkShoeAttributeP::GetNumberOfPoints()
{
  return this->PointData->GetNumberOfTuples();
}

vtkIdType vtkShoeAttributeP::GetNumberOfDOFNodes()
{
  return this->DOFData->GetNumberOfRecords();
}
void vtkShoeAttributeP::SetNumberOfPoints( vtkIdType N )
{
  this->PointData->SetNumberOfTuples( N );
}

void vtkShoeAttributeP::SetNumberOfDOFNodes( vtkIdType N )
{
  this->DOFData->SetNumberOfRecords( N );
  for ( vtkstd::vector<vtkDataRecords*>::iterator cpit = this->CriticalPoints.begin(); cpit != this->CriticalPoints.end(); ++cpit )
    (*cpit)->SetNumberOfRecords( N );
  this->CriticalPointsComputed->SetNumberOfTuples( N );
}

int vtkShoeAttributeP::GetNumberOfComponents()
{
  return this->PointData->GetNumberOfComponents();
}

void vtkShoeAttributeP::SetNumberOfComponents( int N )
{
  this->PointData->SetNumberOfComponents( N );
  this->DOFData->SetNumberOfComponents( N );
  // CriticalPoints have:
  // 1 component specifying the type of extremum (index of Hessian)
  // 3 component specifying the parametric cooordinates of the critical point
  // N components for the value at the critical point
  int cpsz = int(this->CriticalPoints.size());
  if ( N > cpsz )
    {
    int k;
    for ( k = 0; k < cpsz; ++k )
      {
      this->CriticalPoints[k]->Delete();
      vtkDataRecords* cpComponent = vtkDataRecords::Create( VTK_DOUBLE );
      cpComponent->SetNumberOfComponents( N + 3 );
      this->CriticalPoints[k] = cpComponent;
      }
    for ( k = cpsz; k < N; ++k )
      {
      vtkDataRecords* cpComponent = vtkDataRecords::Create( VTK_DOUBLE );
      cpComponent->SetNumberOfComponents( N + 3 );
      this->CriticalPoints.push_back( cpComponent );
      }
    }
  else if ( N < cpsz )
    {
    for ( int k = 0; k < cpsz; ++k )
      {
      this->CriticalPoints[k]->Delete();
      }
    this->CriticalPoints.erase( this->CriticalPoints.begin() + N, this->CriticalPoints.end() );
    for ( int k = 0; k < N; ++k )
      {
      vtkDataRecords* cpComponent = vtkDataRecords::Create( VTK_DOUBLE );
      cpComponent->SetNumberOfComponents( N + 3 );
      this->CriticalPoints[k] = cpComponent;
      }
    }
  this->Parent->MarkCriticalPointsDirty();
  this->ResetRange();
  switch (N)
    {
  case 1:
    this->Type = vtkDataSetAttributes::SCALARS;
    break;
  // NOTE: Do not assume N = 2 => TCOORDS
  case 3:
    this->Type = vtkDataSetAttributes::VECTORS;
    break;
  case 6:
    this->Type = vtkDataSetAttributes::TENSORS;
    break;
  default:
    this->Type = -1;
    break;
    }
}

int vtkShoeAttributeP::GetComponentType()
{
  return this->PointData->GetDataType();
}

void vtkShoeAttributeP::SetStorageType( vtkShoeAttribute* parent, int t )
{
  this->PointData->UnRegister( parent );
  this->DOFData->UnRegister( parent );

  this->PointData = vtkDataArray::CreateDataArray( t );
  this->DOFData = vtkDataRecords::Create( t );

  this->PointData->Register( parent );
  this->DOFData->Register( parent );

  this->PointData->FastDelete();
  this->DOFData->FastDelete();

  this->SetNumberOfComponents( 1 );
  // the above calls this->MarkCriticalPointsDirty();
}

vtkIdType vtkShoeAttributeP::GetSize()
{
  return this->PointData->GetNumberOfTuples() + this->DOFData->GetNumberOfTuples();
}

unsigned long vtkShoeAttributeP::GetActualMemorySize()
{
  unsigned long sz = 
    this->PointData->GetActualMemorySize() + this->DOFData->GetActualMemorySize() + \
    (this->Name.size() + this->TupleCache.size()*sizeof(double) + sizeof(this) + 4 /*parent's pointer to this*/)/1024;
  for ( vtkstd::vector<vtkDataRecords*>::iterator cpit = this->CriticalPoints.begin(); cpit != this->CriticalPoints.end(); ++cpit )
    {
    sz += (*cpit)->GetActualMemorySize();
    }
  sz += this->CriticalPointsComputed->GetActualMemorySize();
  return sz;
}

double* vtkShoeAttributeP::GetRange( vtkRangeStyle rangeStyle, int component )
{
  if ( component == -1 )
    {
    component = this->GetNumberOfComponents();
    }
  
  this->UpdateRange( rangeStyle, component );

  // FIXME: This relies on pair<double,double> being packed tightly (i.e., no padding)
  return &this->Range[component].first;
}

void vtkShoeAttributeP::ResetRange()
{
  // last pair is space for euclidean norm returned when component -1 is requested
  for ( int i=0; i<=this->GetNumberOfComponents(); ++i)
    this->Range.push_back( vtkstd::pair<double,double>(VTK_DOUBLE_MIN,VTK_DOUBLE_MAX) );

  this->RangeStyle = vtkNoRange;
}

double vtkShoeAttributeP::GetMaxNorm()
{
  // FIXME: Eventually, this could be a tighter bound, but is
  // extremely difficult to compute since the max norm is really
  // the max of a derived field (for which we would have to
  // compute new extrema).
  double n = 0.;
  for ( int i=0; i<this->GetNumberOfComponents(); ++i)
    {
    double v;
    v = this->Range[i].first > this->Range[i].second ? this->Range[i].first : this->Range[i].second;
    n += v*v;
    }

  return sqrt(n);
}

void vtkShoeAttributeP::GetTuple( vtkShoeCell*, double* )
{
}

void vtkShoeAttributeP::GetTuple( vtkShoeCellIterator*, double* )
{
}

void vtkShoeAttributeP::GetTuple( vtkShoePointIterator*, double* )
{
}

void vtkShoeAttributeP::GetComponent( int, vtkShoeCellIterator*, double* )
{
}

double vtkShoeAttributeP::GetComponent( int, vtkShoePointIterator* )
{
  return 0.;
}

int vtkShoeAttributeP::SetOrder( int ID, const vtkShoeOrderTuple& o )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( ID );
  if ( result != this->MapSpec.end() )
    { // the cell species already has an order... is it different
    if ( result->second.Order == o )
      {
      return 0;
      }
    }

  vtkShoeCellSpecies* sp = vtkShoeCellSpecies::GetSpeciesById( ID );
  assert(sp && sp->Meta);
  sp->Meta->PrepareForOrder( o );

  result->second.Order = o;
  return 1;
}

int vtkShoeAttributeP::GetOrder( int ID, vtkShoeOrderTuple& orderOut )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( ID );
  if ( result == this->MapSpec.end() )
    { // Could set orderOut to (-1,-1,-1) here to be mean - err... thorough.
    return 0;
    }
  orderOut = result->second.Order;
  return 1;
}

int vtkShoeAttributeP::SetInterpolant( int ID, vtkPolyInterpolant i )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( ID );
  if ( result != this->MapSpec.end() )
    { // the cell species already has an interpolant... is it different
    if ( result->second.Interpolant == i )
      {
      return 0;
      }
    }
  this->MapSpec[ ID ].Interpolant = i;
  return 1;
}

vtkPolyInterpolant vtkShoeAttributeP::GetInterpolant( int ID )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( ID );
  if ( result == this->MapSpec.end() )
    { // Could set orderOut to (-1,-1,-1) here to be mean - err... thorough.
    return (vtkPolyInterpolant)-1;
    }
  return result->second.Interpolant;
}

int vtkShoeAttributeP::SetProductSpace( int ID, vtkPolyProductSpace i )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( ID );
  if ( result != this->MapSpec.end() )
    { // the cell species already has a product space... is it different
    if ( result->second.ProductSpace == i )
      {
      return 0;
      }
    }
  this->MapSpec[ ID ].ProductSpace = i;
  return 1;
}

vtkPolyProductSpace vtkShoeAttributeP::GetProductSpace( int ID )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( ID );
  if ( result == this->MapSpec.end() )
    { // Could set orderOut to (-1,-1,-1) here to be mean - err... thorough.
    return (vtkPolyProductSpace)-1;
    }
  return result->second.ProductSpace;
}

int vtkShoeAttributeP::GetCellTypeInfo( int cellSpeciesId, vtkShoeOrderTuple& orderOut, vtkPolyInterpolant& interpOut, vtkPolyProductSpace& psOut )
{
  vtkstd::map<int,PerCellTypeRec>::iterator result = this->MapSpec.find( cellSpeciesId );
  if ( result == this->MapSpec.end() )
    { // Could set orderOut to (-1,-1,-1) here to be mean - err... thorough.
    return 0;
    }
  orderOut = result->second.Order;
  interpOut = result->second.Interpolant;
  psOut = result->second.ProductSpace;
  return 1;
}

void vtkShoeAttributeP::UpdateRange( vtkRangeStyle rangeStyle, int vtkNotUsed(component) )
{
  if ( this->RangeStyle >= rangeStyle )
    { // we have tighter bounds than requested or at least meet the request
    return;
    }

  if ( ! this->PointData || this->PointData->GetNumberOfTuples() < 1 )
    { // no data
    return;
    }

  // FIXME
  // - ignore component for now and recompute everything
  // - assume values taken on at corners are "good enough"
  int c;
  int nc = this->PointData->GetNumberOfComponents();
  vtkstd::vector<double> tuple;
  tuple.reserve( nc );

  // Save ourselves some work by initializing the range
  // to the first tuple (means we can use "else if"
  // below instead of testing every point against both
  // min and max)
  this->PointData->GetTuple( 0, &tuple[0] );
  for ( c = 0; c < nc; ++ c )
    {
    this->Range[c].first = this->Range[c].second = tuple[c];
    }

  // Sift through the corners
  for ( vtkIdType i = 1; i < this->PointData->GetNumberOfTuples(); ++ i )
    {
    this->PointData->GetTuple( i, &tuple[0] );
    for ( c = 0; c < nc; ++ c )
      {
      if ( tuple[c] < this->Range[c].first )
        {
        this->Range[c].first = tuple[c];
        }
      else if ( tuple[c] > this->Range[c].second )
        {
        this->Range[c].second = tuple[c];
        }
      }
    }

  if ( rangeStyle == vtkTightRange )
    {
    int nc = this->GetNumberOfComponents();
    vtkDataRecords* cp;
    vtkDataRecordsIterator* cpit = vtkDataRecordsIterator::New();
    cpit->SetData( this->DOFData );
    for ( cpit->Begin(); ! cpit->IsAtEnd(); cpit->Next() )
      {
      vtkIdType r = cpit->GetRecord(); // this is a DOF ID
      if ( ! this->CriticalPointsComputed->GetValue( r ) )
        {
      // find a mesh whose DOFLookup contains an entry E for r
      // Take cell C(E) and call ComputeDOFCriticalPoints
      // meshes referencing this attribute
        for( vtkstd::set<vtkShoeBox*>::iterator mit =  this->References.begin(); mit != this->References.end(); ++ mit )
          {
          if ( (*mit)->GetDOFLinks()->GetNumberOfTuples( r ) )
            {
            vtkIdType* cells = (*mit)->GetDOFLinks()->GetRecordAsIdType( r );
            vtkShoeCell* sc = (*mit)->GetCell( cells[0] );
            sc->ComputeDOFCriticalPoints( this->Parent );
            sc->Delete();
            break;
            }
          }
        }

      for ( int c = 0; c < nc; ++ c )
        {
        cp = this->CriticalPoints[c];
        int nt = cp->GetNumberOfTuples( r );
        double* tuples = cp->GetRecordAsDouble( r );
        for ( int t = 0; t < nt; ++ t )
          {
          if ( tuples[3 + c] < this->Range[c].first )
            {
            this->Range[c].first = tuples[3 + c];
            }
          else if ( tuples[3 + c] > this->Range[c].second )
            {
            this->Range[c].second = tuples[3 + c];
            }
          tuples += nc + 3;
          }
        }
      }
    cpit->Delete();
    }

  this->RangeStyle = rangeStyle;
}

void vtkShoeAttributeP::DeepCopy( vtkShoeAttributeP* other )
{
  this->Name = other->Name;
  this->Type = other->Type;
  this->Centering = other->Centering;
  // Don't copy Parent
  // Don't copy References
  if ( this->PointData != other->PointData )
    {
    if ( this->PointData )
      this->PointData->UnRegister( this->Parent );

    if ( other->PointData )
      {
      this->PointData = vtkDataArray::CreateDataArray( other->PointData->GetDataType() );
      this->PointData->DeepCopy( other->PointData );
      }
    else
      {
      this->PointData = 0;
      }
    }

  if ( this->DOFData != other->DOFData )
    {
    if ( other->DOFData )
      {
      if ( ! this->DOFData )
        {
        this->DOFData = vtkDataRecords::New();
        }
      this->DOFData->DeepCopy( other->DOFData );
      }
    else
      {
      this->DOFData = 0;
      }
    }

  this->MapSpec = other->MapSpec;

  // Copy CriticalPoints
  vtkstd::vector<vtkDataRecords*>::iterator cit;
  for ( cit = this->CriticalPoints.begin(); cit != this->CriticalPoints.end(); ++cit )
    {
    if ( *cit )
      {
      (*cit)->UnRegister( this->Parent );
      }
    }
  this->CriticalPoints.clear();
  for ( cit = other->CriticalPoints.begin(); cit != other->CriticalPoints.end(); ++cit )
    {
    vtkDataRecords* cprec = *cit;
    if ( cprec )
      {
      vtkDataRecords* arr = vtkDataRecords::New();
      arr->DeepCopy( cprec );
      cprec = arr;
      }
    this->CriticalPoints.push_back( cprec );
    }
  this->CriticalPointsComputed->DeepCopy( other->CriticalPointsComputed );
  this->CriticalPointsDirty = other->CriticalPointsDirty;
  this->Range = other->Range;
  this->RangeStyle = other->RangeStyle;
  this->TupleCache = other->TupleCache;
}

void vtkShoeAttributeP::ShallowCopy( vtkShoeAttributeP* other )
{
  this->Name = other->Name;
  this->Type = other->Type;
  this->Centering = other->Centering;
  // Don't copy Parent
  // Don't copy References
  if ( this->PointData )
    this->PointData->UnRegister( this->Parent );
  this->PointData = other->PointData;
  if ( this->PointData )
    this->PointData->Register( this->Parent );

  if ( this->DOFData )
    this->DOFData->UnRegister( this->Parent );
  this->DOFData = other->DOFData;
  if ( this->DOFData )
    this->DOFData->Register( this->Parent );

  this->MapSpec = other->MapSpec;

  // Copy CriticalPoints
  vtkstd::vector<vtkDataRecords*>::iterator cit;
  for ( cit = this->CriticalPoints.begin(); cit != this->CriticalPoints.end(); ++cit )
    {
    if ( *cit )
      {
      (*cit)->UnRegister( this->Parent );
      }
    }
  this->CriticalPoints.clear();
  for ( cit = other->CriticalPoints.begin(); cit != other->CriticalPoints.end(); ++cit )
    {
    vtkDataRecords* cprec = *cit;
    this->CriticalPoints.push_back( cprec );
    if ( cprec )
      {
      cprec->Register( this->Parent );
      }
    }
  this->CriticalPointsComputed->DeepCopy( other->CriticalPointsComputed );
  this->Range = other->Range;
  this->RangeStyle = other->RangeStyle;
  this->TupleCache = other->TupleCache;
}

// =EPM========================================================================

vtkShoeAttribute::vtkShoeAttribute()
{
  this->Data = 0;
  this->Data = new vtkShoeAttributeP( this );
  this->RangeStyle = VTK_RANGE_SLOPPY;
}

vtkShoeAttribute::~vtkShoeAttribute()
{
  delete this->Data;
}

void vtkShoeAttribute::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os, indent );
  os << indent << "Data: " << this->Data << vtkstd::endl;
}

const char* vtkShoeAttribute::GetName()
{
  return this->Data->Name.c_str();
}

int vtkShoeAttribute::GetNumberOfComponents()
{
  return this->Data->GetNumberOfComponents();
}

void vtkShoeAttribute::SetNumberOfComponents( int N )
{
  if ( N <= 0 )
    {
    vtkErrorMacro( "Cannot have " << N << " components" );
    return;
    }

  if ( this->Data->GetNumberOfComponents() == N )
    {
    return;
    }

  this->Modified();
  this->Data->SetNumberOfComponents( N );
}

int vtkShoeAttribute::GetCentering()
{
  return this->Data->Centering;
}

int vtkShoeAttribute::GetType()
{
  return this->Data->Type;
}

void vtkShoeAttribute::SetType( int T )
{
  if ( T == this->Data->Type )
    {
    return;
    }

  switch (T)
    {
  case vtkDataSetAttributes::VECTORS:
    if ( this->Data->GetNumberOfComponents() == 3 )
      {
      this->Data->Type = T;
      this->Modified();
      }
    else
      {
      vtkErrorMacro( "VECTORS must have 3 components" );
      }
    break;

  case vtkDataSetAttributes::NORMALS:
    if ( this->Data->GetNumberOfComponents() == 3 )
      {
      this->Data->Type = T;
      this->Modified();
      }
    else
      {
      vtkErrorMacro( "NORMALS must have 3 components" );
      }
    break;

  case vtkDataSetAttributes::TCOORDS:
    if ( this->Data->GetNumberOfComponents() == 2 || this->Data->GetNumberOfComponents() == 3 )
      {
      this->Data->Type = T;
      this->Modified();
      }
    else
      {
      vtkErrorMacro( "TCOORDS must have 2 or 3 components" );
      }
    break;

  default:
    // This can't be SCALARS or TENSORS because Type == T is guaranteed for them.
    // The only question here is whether to allow users to set the type to -1 or not.
    // Certainly, it is an odd situation when there are 2 components, since those
    // may be marked as TCOORDS but not reverted to their original -1.
    vtkErrorMacro( "Invalid attribute type (" << T << ")" );
    break;
    }
}

int vtkShoeAttribute::GetComponentType()
{
  return this->Data->GetComponentType();
}

void vtkShoeAttribute::SetStorageType( int t )
{
  if ( t == this->Data->GetComponentType() )
    {
    return;
    }

  if ( t < 0 || t > VTK_ID_TYPE )
    {
    vtkErrorMacro( "Invalid storage type (" << t << ")" );
    return;
    }

  this->Data->SetStorageType( this, t );
}

void vtkShoeAttribute::SetRange( double* range, int style )
{
  this->Data->ResetRange();
  this->Data->RangeStyle = vtkRangeStyle( style );
  for ( int c = 0; c < this->GetNumberOfComponents(); ++c )
    {
    this->Data->Range[c].first = *(range++);
    this->Data->Range[c].second = *(range++);
    }
}

vtkIdType vtkShoeAttribute::GetSize()
{
  return this->Data->GetSize();
}

unsigned long vtkShoeAttribute::GetActualMemorySize()
{
  return this->Data->GetActualMemorySize();
}

double* vtkShoeAttribute::GetRange( int component )
{
  if ( component < -1 || component >= this->Data->GetNumberOfComponents() )
    {
    vtkErrorMacro( "Invalid component (" << component << ")" );
    return 0;
    }

  vtkRangeStyle rangeStyle;
  bool criticalPointsNeeded = false;
  switch( this->RangeStyle )
    {
    case 1:
      rangeStyle = vtkSloppyRange;
      break;
    case 2:
      rangeStyle = vtkProperRange;
      break;
    case 3:
      rangeStyle = vtkTightRange;
      criticalPointsNeeded = true;
      break;
    default:
      rangeStyle = vtkNoRange;
    }

  return this->Data->GetRange( rangeStyle, component );
}

void vtkShoeAttribute::GetRange( int component, double range[2] )
{
  if ( component < -1 || component >= this->Data->GetNumberOfComponents() )
    {
    vtkErrorMacro( "Invalid component (" << component << ")" );
    return;
    }

  vtkRangeStyle rangeStyle;
  bool criticalPointsNeeded = false;
  switch( this->RangeStyle )
    {
    case 1:
      rangeStyle = vtkSloppyRange;
      break;
    case 2:
      rangeStyle = vtkProperRange;
      break;
    case 3:
      rangeStyle = vtkTightRange;
      criticalPointsNeeded = true;
      break;
    default:
      rangeStyle = vtkNoRange;
    }

  double* R = this->Data->GetRange( rangeStyle, component );
  range[0] = R[0];
  range[1] = R[1];
}

double vtkShoeAttribute::GetMaxNorm()
{
  // FIXME: Eventually, this could be a tighter bound, but is
  // extremely difficult to compute since the max norm is really
  // the max of a derived field (for which we would have to
  // compute new extrema).
  return this->Data->GetMaxNorm();
}

double* vtkShoeAttribute::GetTuple( vtkGenericAdaptorCell* c )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoeCell* C = dynamic_cast<vtkShoeCell*>(c);
  this->Data->GetTuple(C,&this->Data->TupleCache[0]);
  return &this->Data->TupleCache[0];
}

void vtkShoeAttribute::GetTuple( vtkGenericAdaptorCell* c, double* tuple )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoeCell* C = dynamic_cast<vtkShoeCell*>(c);
  this->Data->GetTuple(C,tuple);
}

double* vtkShoeAttribute::GetTuple( vtkGenericCellIterator* c )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoeCellIterator* C = dynamic_cast<vtkShoeCellIterator*>(c);
  this->Data->GetTuple(C,&this->Data->TupleCache[0]);
  return &this->Data->TupleCache[0];
}

void vtkShoeAttribute::GetTuple( vtkGenericCellIterator* c, double* tuple )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoeCellIterator* C = dynamic_cast<vtkShoeCellIterator*>(c);
  this->Data->GetTuple(C,tuple);
}

double* vtkShoeAttribute::GetTuple( vtkGenericPointIterator* p )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoePointIterator* P = dynamic_cast<vtkShoePointIterator*>(p);
  this->Data->GetTuple(P,&this->Data->TupleCache[0]);
  return &this->Data->TupleCache[0];
}

void vtkShoeAttribute::GetTuple( vtkGenericPointIterator* p, double* tuple )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoePointIterator* P = dynamic_cast<vtkShoePointIterator*>(p);
  this->Data->GetTuple(P,tuple);
}

void vtkShoeAttribute::GetComponent( int i, vtkGenericCellIterator* c, double* values )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoeCellIterator* C = dynamic_cast<vtkShoeCellIterator*>(c);
  this->Data->GetComponent( i, C, values );
}

double vtkShoeAttribute::GetComponent( int i, vtkGenericPointIterator* p )
{
  // FIXME: Should be static_cast<>() for speed.
  vtkShoePointIterator* P = dynamic_cast<vtkShoePointIterator*>(p);
  return this->Data->GetComponent( i, P );
}

void vtkShoeAttribute::DeepCopy( vtkGenericAttribute* other )
{
  if ( ! other )
    {
    return;
    }

  vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( other );
  if ( ! sa )
    {
    vtkErrorMacro( "Cannot copy a non-SHOE attribute to a SHOE attribute" );
    return;
    }

  if ( sa == this )
    {
    return;
    }

  this->Data->DeepCopy( sa->Data );
  this->RangeStyle = sa->RangeStyle;
}

void vtkShoeAttribute::ShallowCopy( vtkGenericAttribute* other )
{
  if ( ! other )
    {
    return;
    }

  vtkShoeAttribute* sa = vtkShoeAttribute::SafeDownCast( other );
  if ( ! sa )
    {
    vtkErrorMacro( "Cannot copy a non-SHOE attribute to a SHOE attribute" );
    return;
    }

  if ( sa == this )
    {
    return;
    }

  this->Data->ShallowCopy( sa->Data );
  this->RangeStyle = sa->RangeStyle;
}

void vtkShoeAttribute::SetPointData( vtkDataArray* a )
{
  if ( this->Data->PointData == a || a == 0 )
    return;

  this->Data->PointData->Delete();
  this->Data->PointData = a;
  this->Data->PointData->Register( this );
  this->Data->ResetRange();
  this->Modified();
}

void vtkShoeAttribute::SetDOFData( vtkDataRecords* r )
{
  if ( this->Data->DOFData == r || r == 0 )
    return;

  this->Data->DOFData->Delete();
  this->Data->DOFData = r;
  this->Data->DOFData->Register( this );
  this->Data->ResetRange();
  vtkIdType nr = r->GetNumberOfRecords();
  int nc = r->GetNumberOfComponents();
  int csz = this->Data->CriticalPoints.size();
  this->Data->CriticalPointsComputed->SetNumberOfTuples( nr );
  this->MarkCriticalPointsDirty();
  if ( csz < nc )
    {
    for ( int i = csz; i < nc; ++i )
      {
      this->Data->CriticalPoints.push_back( vtkDataRecords::Create( VTK_DOUBLE ) );
      this->Data->CriticalPoints[i]->SetNumberOfComponents(4);
      }
    }
  else if ( csz > nc )
    {
    for ( int i = csz; i < nc; ++i )
      {
      this->Data->CriticalPoints[i]->Delete();
      }
    this->Data->CriticalPoints.erase(
      this->Data->CriticalPoints.begin() + csz, this->Data->CriticalPoints.end() );
    }
  for ( int j = 0; j < nc; ++j )
    {
    this->Data->CriticalPoints[j]->SetNumberOfRecords(nr);
    }
  this->Modified();
}

void vtkShoeAttribute::SetName( const char* name )
{
  if ( !strcmp( name, this->Data->Name.c_str() ) )
    {
    return;
    }

  this->Modified();
  this->Data->Name = vtkstd::string( name );
}

void vtkShoeAttribute::SetNumberOfPoints( vtkIdType N )
{
  if ( this->Data->PointData->GetNumberOfTuples() == N )
    {
    return;
    }

  this->Data->PointData->SetNumberOfTuples( N );
}

void vtkShoeAttribute::SetNumberOfDOFNodes( vtkIdType N )
{
  if ( N < 0 )
    {
    vtkErrorMacro("Cannot request a negative number (" << N << ") of DOF nodes");
    }

  if ( this->Data->DOFData->GetNumberOfRecords() == N )
    {
    return;
    }

  this->Data->SetNumberOfDOFNodes( N );
}

void vtkShoeAttribute::SetOrder( int cellSpeciesId, vtkShoeOrderTuple& order )
{
  if ( this->Data->SetOrder( cellSpeciesId, order ) )
    {
    this->Modified();
    }
}

int vtkShoeAttribute::GetOrder( int cellSpeciesId, vtkShoeOrderTuple& orderOut )
{
  return this->Data->GetOrder( cellSpeciesId, orderOut );
}

void vtkShoeAttribute::SetInterpolant( int cellSpeciesId, int interpolant )
{
  if ( this->Data->SetInterpolant( cellSpeciesId, (vtkPolyInterpolant) interpolant ) )
    {
    this->Modified();
    }
}

int vtkShoeAttribute::GetInterpolant( int cellSpeciesId )
{
  return (int) this->Data->GetInterpolant( cellSpeciesId );
}

void vtkShoeAttribute::SetProductSpace( int cellSpeciesId, int productSpace )
{
  if ( this->Data->SetProductSpace( cellSpeciesId, (vtkPolyProductSpace) productSpace ) )
    {
    this->Modified();
    }
}

int vtkShoeAttribute::GetProductSpace( int cellSpeciesId )
{
  return (int) this->Data->GetProductSpace( cellSpeciesId );
}

int vtkShoeAttribute::GetCellTypeInfo( int cellSpeciesId, vtkShoeOrderTuple& orderOut, vtkPolyInterpolant& interpOut, vtkPolyProductSpace& psOut )
{
  return this->Data->GetCellTypeInfo( cellSpeciesId, orderOut, interpOut, psOut );
}

void vtkShoeAttribute::SetCriticalPointsDirty( int i )
{
  this->Data->CriticalPointsDirty = i;
}

int vtkShoeAttribute::GetCriticalPointsDirty() const
{
  return this->Data->CriticalPointsDirty;
}

vtkDataArray* vtkShoeAttribute::GetPointData()
{
  return this->Data->PointData;
}

vtkDataRecords* vtkShoeAttribute::GetDOFData()
{
  return this->Data->DOFData;
}

vtkIdType vtkShoeAttribute::GetNumberOfPoints()
{
  return this->Data->GetNumberOfPoints();
}

vtkIdType vtkShoeAttribute::GetNumberOfDOFNodes()
{
  return this->Data->GetNumberOfDOFNodes();
}

vtkDataRecords** vtkShoeAttribute::GetCriticalPoints()
{
  return &this->Data->CriticalPoints[0];
}

vtkBitArray* vtkShoeAttribute::GetCriticalPointsComputed()
{
  return this->Data->CriticalPointsComputed;
}

void vtkShoeAttribute::PrintCriticalPoints( ostream& os )
{
  int nc = this->GetNumberOfComponents();
  vtkDataRecordsIterator* dit = vtkDataRecordsIterator::New();
  dit->SetData( this->Data->DOFData );
  for ( dit->Begin(); ! dit->IsAtEnd(); dit->Next() )
    {
    vtkIdType r = dit->GetRecord();
    os << "DOF " << r << ": " << vtkstd::endl;
    for ( int c = 0; c < nc; ++c )
      {
      int nt = this->Data->CriticalPoints[c]->GetNumberOfTuples( r );
      double* tuples = this->Data->CriticalPoints[c]->GetRecordAsDouble( r );
      for ( int t = 0; t < nt; ++t )
        {
        os << "  " << c << ", " << t << ": " << tuples[0] << ", " << tuples[1] << ", " << tuples[2] << vtkstd::endl;
        tuples += 4;
        }
      }
    }
  dit->Delete();
}

void vtkShoeAttribute::InsertMeshReference( vtkShoeBox* sb )
{
  this->Data->References.insert( sb );
}
  
void vtkShoeAttribute::RemoveMeshReference( vtkShoeBox* sb )
{
  this->Data->References.erase( sb );
}

void vtkShoeAttribute::MarkCriticalPointsDirty()
{
  if ( this->Data )
    {
    this->Data->CriticalPointsDirty = 2;
    // Quickly set all the bits to 0
    vtkIdType mid = this->Data->CriticalPointsComputed->GetMaxId();
    if ( mid >= 0 )
      {
      memset( this->Data->CriticalPointsComputed->GetPointer(0), 0, mid / 8 + 1 );
      }
    }
}


