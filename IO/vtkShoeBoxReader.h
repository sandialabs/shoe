/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#ifndef __vtkShoeBoxReader_h
#define __vtkShoeBoxReader_h

#include "vtksnlIOWin32Header.h"

#include <vtkShoeBoxAlgorithm.h>

#include <vtkstd/vector> // for CellArrayNames, PointArrayNames
#include <vtkstd/string> // for CellArrayNames, PointArrayNames
#include <vtkstd/map> // for Species list

class vtkUnstructuredAlgorithm;
class vtkDataRecords;
class vtkDoubleArray;
class vtkShoeCellSpecies;

class VTK_SNL_IO_EXPORT vtkShoeBoxReader : public vtkShoeBoxAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkShoeBoxReader,vtkShoeBoxAlgorithm);
  virtual void PrintSelf(ostream &os, vtkIndent indent);
  static vtkShoeBoxReader *New();

  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  vtkGetStringMacro(Title);
  vtkGetStringMacro(Description);
  vtkGetMacro(NumberOfPoints,vtkIdType);
  vtkGetMacro(NumberOfDofNodes,vtkIdType);
  vtkGetMacro(NumberOfCells,vtkIdType);
  vtkIdType GetNumberOfCells( int dim );

  vtkIdType GetNumberOfTimeSteps() const;
  vtkGetVector2Macro(TimeStepRange,int);
  vtkSetMacro(TimeStep,int);
  vtkGetMacro(TimeStep,int);
  double GetTimeValue( int ) const;

  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int);
  virtual void SetPointArrayStatus( int arr, int status );
  virtual void SetPointArrayStatus( const char* arrName, int status );
  int GetPointArrayStatus( int arr );
  int GetPointArrayStatus( const char* arrName );

  int GetNumberOfCellArrays();
  const char* GetCellArrayName(int);
  virtual void SetCellArrayStatus( int arr, int status );
  virtual void SetCellArrayStatus( const char* arrName, int status );
  int GetCellArrayStatus( int arr );
  int GetCellArrayStatus( const char* arrName );

  int GetNumberOfComponentsInCellArray( int a );
  int GetNumberOfComponentsInCellArray( const char* arrName );

  int GetNumberOfComponentsInPointArray( int a );
  int GetNumberOfComponentsInPointArray( const char* arrName );

  vtkGetMacro(Binary,int);

protected:
  char* FileName;
  char* Title;
  char* Description;
  int TimeStepRange[2];
  int TimeStep;
  //BTX
  vtkstd::vector<double> TimeValues;
  //ETX
  vtkIdType NumberOfPoints;
  vtkIdType NumberOfCells;
  vtkIdType NumberOfDofNodes;
  vtkIdType NumberOfCellsByDimension[4];
  vtkIdType GeometryIndex;
  vtkShoeBox* Output; // Only valid during RequestData()
  //BTX
  vtkstd::vector<vtkstd::string> CellArrayNames;
  vtkstd::vector<vtkstd::string> PointArrayNames;
  vtkstd::vector<vtkstd::vector<double> > CellArrayBounds;
  vtkstd::vector<vtkstd::vector<double> > PointArrayBounds;
  vtkstd::vector<int> CellArrayBoundsType;
  vtkstd::vector<int> PointArrayBoundsType;
  vtkstd::vector<int> CellArrayStatus;
  vtkstd::vector<int> PointArrayStatus;
  vtkstd::vector<int> CellArraySize;
  vtkstd::vector<int> PointArraySize;
  vtkstd::map<vtkIdType,vtkIdType> SpeciesSwizzle;
  //ETX
  int Binary;

  vtkShoeBoxReader();
  virtual ~vtkShoeBoxReader();

  virtual int RequestInformation( vtkInformation*, vtkInformationVector**, vtkInformationVector* );
  virtual int RequestData( vtkInformation*, vtkInformationVector**, vtkInformationVector* );

  vtkSetStringMacro(Title);
  vtkSetStringMacro(Description);

  // Description:
  // Read summary information (without reading the actual mesh).
  // This is called by vtkShoeBoxReader::RequestInformation().
  virtual int ReadSummary();

  // Description:
  // Read bits of the mesh.
  // These routines are all called by vtkShoeBoxReader::RequestData().
  virtual int ReadAttributes();
  virtual int ReadCells();
  virtual int ReadConnectivity();
  virtual int ReadDOFLinks();
  virtual int ReadPointLinks();
  virtual int ReadSpecies();

  // Description:
  // Read in an attribute.
  // This routine is called by vtkShoeBoxReader::ReadAttributes() and assumes that
  // the array status is set to 1 (i.e., the user has requested that this attribute be loaded).
  virtual int ReadAttribute( vtkIdType a );

  // Description:
  // Utility routines to read arrays used to hold mesh data.
  virtual int ReadDataRecordsIdType( istream&, vtkDataRecords* );
  virtual int ReadDataRecordsDouble( istream&, vtkDataRecords* );
  virtual int ReadDoubleArray( istream&, vtkDoubleArray* );

private:
  vtkShoeBoxReader( const vtkShoeBoxReader& ); // Not implemented.
  void operator = ( const vtkShoeBoxReader& ); // Not implemented.
};

//BTX
inline void vtkShoeBoxReader::SetPointArrayStatus( const char* arrName, int status )
{
  for ( int a=0; a<this->GetNumberOfPointArrays(); ++a )
    {
    if ( strcmp( arrName, this->GetPointArrayName(a) ) == 0 )
      {
      this->SetPointArrayStatus( a, status );
      return;
      }
    }
  vtkWarningMacro( "Point array \"" << arrName << "\" does not exist" );
}

inline int vtkShoeBoxReader::GetPointArrayStatus( const char* arrName )
{
  for ( int a=0; a<this->GetNumberOfPointArrays(); ++a )
    {
    if ( strcmp( arrName, this->GetPointArrayName(a) ) == 0 )
      {
      return this->GetPointArrayStatus( a );
      }
    }
  //vtkWarningMacro( "Point array \"" << arrName << "\" does not exist" );
  return 0;
}

inline void vtkShoeBoxReader::SetCellArrayStatus( const char* arrName, int status )
{
  for ( int a=0; a<this->GetNumberOfCellArrays(); ++a )
    {
    if ( strcmp( arrName, this->GetCellArrayName(a) ) == 0 )
      {
      this->SetCellArrayStatus( a, status );
      return;
      }
    }
  vtkWarningMacro( "Cell array \"" << arrName << "\" does not exist" );
}

inline int vtkShoeBoxReader::GetCellArrayStatus( const char* arrName )
{
  for ( int a=0; a<this->GetNumberOfCellArrays(); ++a )
    {
    if ( strcmp( arrName, this->GetCellArrayName(a) ) == 0 )
      {
      return this->GetCellArrayStatus( a );
      }
    }
  //vtkWarningMacro( "Cell array \"" << arrName << "\" does not exist" );
  return 0;
}

inline int vtkShoeBoxReader::GetNumberOfComponentsInCellArray( const char* arrName )
{
  for ( int a=0; a<this->GetNumberOfCellArrays(); ++a )
    {
    if ( strcmp( arrName, this->GetCellArrayName(a) ) == 0 )
      {
     return this->GetNumberOfComponentsInCellArray( a );
      }
    }
  //vtkWarningMacro( "Cell array \"" << arrName << "\" does not exist" );
  return 0;
}

inline int vtkShoeBoxReader::GetNumberOfComponentsInPointArray( const char* arrName )
{
  for ( int a=0; a<this->GetNumberOfPointArrays(); ++a )
    {
    if ( strcmp( arrName, this->GetPointArrayName(a) ) == 0 )
      {
     return this->GetNumberOfComponentsInPointArray( a );
      }
    }
  //vtkWarningMacro( "Point array \"" << arrName << "\" does not exist" );
  return 0;
}
//ETX

#endif // __vtkShoeBoxReader_h
