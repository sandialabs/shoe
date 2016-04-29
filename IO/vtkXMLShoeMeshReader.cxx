/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <iostream>

#include <vtkSystemIncludes.h>
#include <vtkObjectFactory.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLDataParser.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkDataArray.h>
#include <vtkExecutive.h>

#include <vtkXMLShoeMeshReader.h>
#include <vtkShoeMesh.h>
#include <vtkFunctionData.h>

using namespace shoe;

vtkCxxRevisionMacro(vtkXMLShoeMeshReader,"$Revision:");
vtkStandardNewMacro(vtkXMLShoeMeshReader);

vtkXMLShoeMeshReader::vtkXMLShoeMeshReader()
{
}

vtkXMLShoeMeshReader::~vtkXMLShoeMeshReader()
{
}

vtkIdType vtkXMLShoeMeshReader::GetNumberOfCellsInPiece( int piece_in )
{
	return this->NumberOfCells[ piece_in ];
}

void vtkXMLShoeMeshReader::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf( os, indent );
}

void vtkXMLShoeMeshReader::SetOutput( vtkShoeMesh* mesh_in )
{
	this->GetExecutive()->SetOutputData( 0, mesh_in );
}

vtkShoeMesh* vtkXMLShoeMeshReader::GetOutput()
{
	return static_cast<vtkShoeMesh*>( this->GetOutputDataObject(0) );
}

const char* vtkXMLShoeMeshReader::GetDataSetName()
{
	return "ShoeMesh";
}

void vtkXMLShoeMeshReader::GetOutputUpdateExtent( int& piece_out, int& num_pieces_out, int& ghost_level_out )
{
	vtkShoeMesh* mesh = this->GetOutput();

	piece_out = mesh->GetUpdatePiece();
	num_pieces_out = mesh->GetUpdateNumberOfPieces();
	ghost_level_out = mesh->GetUpdateGhostLevel();
}

vtkIdType vtkXMLShoeMeshReader::GetNumberOfDofNodes() const
{
	return TotalNumberOfDofNodes;
}

void vtkXMLShoeMeshReader::SetupOutputTotals()
{
	this->Superclass::SetupOutputTotals();

	// Find the total size of the output.
	vtkIdType i;
	this->TotalNumberOfCells = 0;
	this->TotalNumberOfDofNodes = 0;
	for( i=this->StartPiece; i<this->EndPiece; ++i )
	{
		this->TotalNumberOfCells += this->NumberOfCells[i];
		this->TotalNumberOfDofNodes += this->NumberOfDofNodes[i];
	}
}

void vtkXMLShoeMeshReader::SetupPieces( int num_pieces_in )
{
	this->Superclass::SetupPieces( num_pieces_in );

	this->NumberOfFunctionFields.resize( num_pieces_in );
	this->NumberOfPointFields.resize( num_pieces_in );
	this->NumberOfDofNodes.resize( num_pieces_in );
	this->FieldDofElements.resize( num_pieces_in );
	this->GeometryDofElements.resize( num_pieces_in );
	this->CellElements.resize( num_pieces_in );
	this->NumberOfCells.resize( num_pieces_in );
	for ( int i=0; i<num_pieces_in; i++ )
	{
		this->FieldDofElements[i] = 0;
		this->GeometryDofElements[i] = 0;
		this->CellElements[i] = 0;
	}
}

void vtkXMLShoeMeshReader::DestroyPieces()
{
	this->NumberOfCells.clear();
	this->NumberOfDofNodes.clear();
	this->FieldDofElements.clear();
	this->GeometryDofElements.clear();
	this->CellElements.clear();
	this->NumberOfPointFields.clear();
	this->NumberOfFunctionFields.clear();

	this->Superclass::DestroyPieces();
}

int vtkXMLShoeMeshReader::FillOutputPortInformation( int port, vtkInformation* info )
{
  (void)port;
  info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkShoeMesh");
  return 1;
}

void vtkXMLShoeMeshReader::SetupOutputData()
{
	this->Superclass::SetupOutputData();

	vtkShoeMesh* mesh = this->GetOutput();
	mesh->SetNumberOfDofNodes( this->GetNumberOfDofNodes() );
}

int vtkXMLShoeMeshReader::ReadPiece( vtkXMLDataElement* ePiece )
{
	if ( !this->Superclass::ReadPiece( ePiece ) )
		return 0;

	if( !ePiece->GetScalarAttribute("NumberOfCells", this->NumberOfCells[ this->Piece ]) )
	{
		vtkErrorMacro( "Piece " << this->Piece << " is missing its NumberOfCells attribute." );
		this->NumberOfCells[ this->Piece ] = 0;
		return 0;
	}

	if( !ePiece->GetScalarAttribute("NumberOfDofNodes", this->NumberOfDofNodes[ this->Piece ]) )
	{
		vtkErrorMacro( "Piece " << this->Piece << " is missing its NumberOfDofNodes attribute." );
		this->NumberOfDofNodes[ this->Piece ] = 0;
		return 0;
	}

	this->CellElements[ this->Piece ] = 0;
	for ( int i=0; i<ePiece->GetNumberOfNestedElements(); ++i )
	{
		vtkXMLDataElement* ele = ePiece->GetNestedElement( i );
		if( (strcmp(ele->GetName(), "Cells") == 0 ) && (ele->GetNumberOfNestedElements() > 0) )
		{
			this->CellElements[ this->Piece ] = ele;
		}
		else if ( (strcmp( ele->GetName(), "PointData" ) == 0) )
		{
			this->NumberOfPointFields[ this->Piece ] = 0;
			for ( int j=0; j<ele->GetNumberOfNestedElements(); ++j )
			{
				if ( strcmp( ele->GetNestedElement(j)->GetName(), "DataArray" ) == 0 )
					this->NumberOfPointFields[ this->Piece ]++;
			}
		}
		else if ( (strcmp( ele->GetName(), "FunctionData" ) == 0) )
		{
			const char* func_type = ele->GetAttribute( "Type" );
			if ( ! func_type )
				continue;

			if ( strcmp( func_type, "Geometry" ) == 0 )
			{
				this->GeometryDofElements[ this->Piece ] = ele;
			}
			else if ( strcmp( func_type, "Fields" ) == 0 )
			{
				this->FieldDofElements[ this->Piece ] = ele;

				this->NumberOfFunctionFields[ this->Piece ] = 0;
				for ( int j=0; j<ele->GetNumberOfNestedElements(); ++j )
				{
					if ( strcmp( ele->GetNestedElement(j)->GetName(), "VaryingDataArray" ) == 0 )
					{
						this->NumberOfFunctionFields[ this->Piece ]++;
					}
				}
			}
		}
	}

	if( !this->CellElements[ this->Piece ] )
	{
		vtkErrorMacro("A piece is missing its Cells element.");
		return 0;
	}

	return 1;
}

void vtkXMLShoeMeshReader::SetupNextPiece()
{
	this->Superclass::SetupNextPiece();
	// Change starting cell? What good does this do? It doesn't look like
	// vtkUnstructuredGrid does anything but add the next piece to the
	// existing mesh, consuming more memory, ...
}

bool IsFreeRange( vtkXMLDataElement* ele )
{
	return (strcmp( ele->GetName(), "FreeRange" ) == 0);
}

bool IsFreeList( vtkXMLDataElement* ele )
{
	return (strcmp( ele->GetName(), "FreeList" ) == 0);
}

bool IsList( vtkXMLDataElement* ele )
{
	return (strcmp( ele->GetName(), "List" ) == 0);
}

int vtkXMLShoeMeshReader::ReadPieceData()
{
	// FIXME: Figure out how frequently to update progress bar

	// Let vtkXMLUnstructuredDataReader load the points
	if ( !this->Superclass::ReadPieceData() )
		return 0;

	vtkShoeMesh* mesh = this->GetOutput();

	// Find the different arrays that define the cells.
	vtkXMLDataElement* xmlCells = this->CellElements[ this->Piece ];
	vtkXMLDataElement* xmlConn=0;
	vtkXMLDataElement* xmlCellSpec=0;
	vtkXMLDataElement* xmlCellDefn=0;
	vtkXMLDataElement* xmlCornerLinks=0;
	vtkXMLDataElement* xmlDofNodeLinks=0;

	for( int i=0; i<xmlCells->GetNumberOfNestedElements(); ++i )
	{
		vtkXMLDataElement* ele = xmlCells->GetNestedElement(i);
		if ( (strcmp(ele->GetAttribute("Name"), "Connectivity") == 0)
		     && (ele->GetNumberOfNestedElements() > 1))
			xmlConn = ele;
		else if ( (strcmp(ele->GetAttribute("Name"), "CellSpecs") == 0)
		          && (ele->GetNumberOfNestedElements() > 1))
			xmlCellSpec = ele;
		else if ( strcmp(ele->GetAttribute("Name"), "CellDefinitions") == 0 )
			xmlCellDefn = ele;
		else if ( (strcmp(ele->GetAttribute("Name"), "CornerLinks") == 0)
		          && (ele->GetNumberOfNestedElements() > 1))
			xmlCornerLinks = ele;
		else if ( (strcmp(ele->GetAttribute("Name"), "DofNodeLinks") == 0)
		          && (ele->GetNumberOfNestedElements() > 1))
			xmlDofNodeLinks = ele;
	}

	if ( (! xmlConn) || (! IsFreeRange( xmlConn )) )
	{
		vtkErrorMacro( "The current Piece is missing (or has improper) Connectivity information" );
		return 0;
	}

	if ( (! xmlCellSpec) || (! IsFreeList( xmlCellSpec )) )
	{
		vtkErrorMacro( "The current Piece is missing (or has improper) CellSpecs information" );
		return 0;
	}

	if ( (! xmlCellDefn) || (! IsList( xmlCellDefn )) )
	{
		vtkErrorMacro( "The current Piece is missing (or has improper) CellDefinitions information" );
		return 0;
	}

	if ( ! xmlCornerLinks )
	{
		vtkErrorMacro( "The current Piece is missing (or has improper) CornerLinks information" );
		return 0;
	}

	if ( ! xmlDofNodeLinks )
	{
		vtkErrorMacro( "The current Piece is missing (or has improper) DofNodeLinks information" );
		return 0;
	}

	// Allocate spots for all the FunctionData arrays
	mesh->GetFunctionData()->AllocateArrays( this->NumberOfPointFields[ this->Piece ] );

	// Store iterators to the vtkCellDefinitions for use when loading the vtkCellSpecs.
	vtkstd::vector< vtkShoeMesh::CellDefsType::const_iterator > cdIterators;
	if ( ! this->ReadCellDefinitions( mesh->CellDefs, cdIterators, xmlCellDefn ) )
	{
		vtkErrorMacro( "Could not read CellDefs list" );
		return 0;
	}

	if ( ! this->ReadCellSpecs( mesh->Cells, cdIterators, xmlCellSpec, this->GetNumberOfCells() ) )
	{
		vtkErrorMacro( "Could not read Cells array" );
		return 0;
	}

	if ( ! this->ReadFreerangeData( mesh->Connectivity, xmlConn ) )
	{
		vtkErrorMacro( "Could not read Connectivity array" );
		return 0;
	}

	if ( ! this->ReadVaryingDataArray( mesh->GetCornerLinks(), xmlCornerLinks, this->NumberOfPoints[ this->Piece ] ) )
	{
		vtkErrorMacro( "Could not read CornerLinks" );
		return 0;
	}

	if ( ! this->ReadVaryingDataArray( mesh->GetDofNodeLinks(), xmlDofNodeLinks, this->NumberOfDofNodes[ this->Piece ] ) )
	{
		vtkErrorMacro( "Could not read DofNodeLinks" );
		return 0;
	}

	if ( ! this->ReadFunctionData( mesh->GetGeometryData(), this->GeometryDofElements[ this->Piece ], 1 ) )
	{
		vtkErrorMacro( "Could not read FunctionData for fields" );
		return 0;
	}

	if ( this->FieldDofElements[ this->Piece] && ! this->ReadFunctionData( mesh->GetFunctionData(), this->FieldDofElements[ this->Piece], this->NumberOfFunctionFields[ this->Piece ] ) )
	{
		vtkErrorMacro( "Could not read FunctionData for fields" );
		return 0;
	}

	return 1;
}

vtkIdType vtkXMLShoeMeshReader::ReadVaryingDataArray( vtkFunctionData* data_out, int n,
                                                      vtkXMLDataElement* xml_node_in, vtkIdType num_offsets_in )
{
	vtkXMLDataElement* data = 0;
	vtkXMLDataElement* offs = 0;
	vtkIdType count = 0;
	vtkXMLDataElement* ele;

	int id;
	if ( xml_node_in->GetScalarAttribute( "Id", id ) )
		n = id; // Override the field number if the XML file specifies an array position

	for ( int i=0; i<xml_node_in->GetNumberOfNestedElements(); i++ )
	{
		ele = xml_node_in->GetNestedElement( i );

		if ( (strcmp(ele->GetName(), "Values") == 0))
			data = ele;
		else if ( (strcmp(ele->GetName(), "Offsets") == 0))
			offs = ele;
	}

	if ( ! data || ! offs ) {
		vtkErrorMacro( "Varying-length data array is missing Values or Offsets elements" );
		return 0;
	}

	int data_type;
	if ( ! xml_node_in->GetWordTypeAttribute( "Type", data_type ) )
	{
		vtkErrorMacro( "No type specified for data array" );
		return 0;
	}
	vtkDataArray*   arr = vtkDataArray::CreateDataArray( data_type );
	vtkIdTypeArray* off = vtkIdTypeArray::New();
	arr->SetName( xml_node_in->GetAttribute( "Name" ) );
	off->SetName( "Offsets" );
	//cout << "Reading " << ( data_out->IsExtrema() ? "extrema" : "field" ) << " named " << arr->GetName();

	int num_components;
	if ( xml_node_in->GetScalarAttribute( "NumberOfComponents", num_components ) )
		arr->SetNumberOfComponents( num_components );
	else
		arr->SetNumberOfComponents( 1 );

	int num_off_components;
	if ( xml_node_in->GetScalarAttribute( "NumberOfOffsetComponents", num_off_components ) )
	{
		off->SetNumberOfComponents( num_off_components );
	}
	else
	{
		if ( data_out->IsExtrema() ) {
			off->SetNumberOfComponents( 2 );
			num_off_components = 2;
		}
		else
		{
			off->SetNumberOfComponents( 1 );
			num_off_components = 1;
		}
	}

	//cout << " with " << num_offsets_in << " offset tuples (offsets have " << num_off_components << " components)" << endl;

	data_out->InsertValues( n, arr );
	data_out->InsertOffsets( n, off );
	arr->Delete();
	off->Delete();

	vtkIdType data_size;
	if ( ! data->GetScalarAttribute( "Size", data_size ) )
	{
		vtkErrorMacro( "Length of VaryingDataArray Values entry was not specified" );
		return 0;
	}
	data_out->GetOffsets( n )->SetNumberOfTuples( num_offsets_in );
	data_out->GetValues( n )->SetNumberOfTuples( data_size / num_components );

	count += this->ReadArrayValues( data, 0, data_out->GetValues( n ), 0, data_size );
	count += this->ReadArrayValues( offs, 0, data_out->GetOffsets( n ), 0, num_offsets_in*num_off_components );

	return count;
}

// FIXME!!!
// This is the worst way to do this, but it was easier to program than rewriting a lot of the XML parser.
// It should be rewritten the right way at some point or the reader will never be fast or compact.
vtkIdType vtkXMLShoeMeshReader::ReadVaryingDataArray( vtkShoeMesh::LinkType& data_out,
                                                      vtkXMLDataElement* xml_node_in, vtkIdType num_offsets_in )
{
	vtkXMLDataElement* data = 0;
	vtkXMLDataElement* offs = 0;
	vtkIdType count = 0;
	vtkXMLDataElement* ele;

	for ( int i=0; i<xml_node_in->GetNumberOfNestedElements(); i++ )
	{
		ele = xml_node_in->GetNestedElement( i );

		if ( (strcmp(ele->GetName(), "Values") == 0))
			data = ele;
		else if ( (strcmp(ele->GetName(), "Offsets") == 0))
			offs = ele;
	}

	if ( ! data || ! offs ) {
		vtkErrorMacro( "Varying-length data array is missing Values or Offsets elements" );
		return 0;
	}

	int data_type;
	if ( xml_node_in->GetWordTypeAttribute( "Type", data_type ) )
	{
#if 0
		if ( data_type != VTK_ID_TYPE )
		{
			vtkErrorMacro( "Type specified for data array (" << xml_node_in->GetAttribute("Type") << ") is not this machine's vtkIdType!" );
			return 0;
		}
#endif
	}

	//data_out->InsertValues( n, arr );

	vtkIdType data_size;
	if ( ! data->GetScalarAttribute( "Size", data_size ) )
	{
		vtkErrorMacro( "Length of VaryingDataArray Values entry was not specified" );
		return 0;
	}

	vtkIdTypeArray* datv = vtkIdTypeArray::New();
	vtkIdTypeArray* offv = vtkIdTypeArray::New();

	datv->SetNumberOfComponents(1);
	offv->SetNumberOfComponents(1);

	datv->SetNumberOfTuples( data_size );
	offv->SetNumberOfTuples( num_offsets_in + 1 );

	data_out.resize( num_offsets_in );

	count += this->ReadArrayValues( data, 0, datv, 0, data_size );
	count += this->ReadArrayValues( offs, 0, offv, 0, num_offsets_in );

	offv->InsertValue( num_offsets_in, data_size );

	// OK, I'm sure there's an easier, less mind-bending way to do this... but...
	// there is a special case to handle.
	// Consider the case when we are reading values for the first vector
	// in the array of vectors. If the _second_ vector in the array is empty, there
	// will be a -1 entry in the offsets. So we must increment i until offv->GetValue(i) >= 0
	// in order to get the number of entries in the first vector.
	// But we must remember the number of times we've incremented i so we can skip over
	// that many vectors in the array.
	int i=1;
	int c=0;
	int skip=0;
	// FIXME: Assumes offv->GetValue(0) == 0
	for ( vtkShoeMesh::LinkType::iterator it = data_out.begin(); it != data_out.end(); ++it, ++i )
	{
		while ( i <= num_offsets_in && offv->GetValue(i) == -1 )
		{
			++i;
			++skip;
		}
		if ( i > num_offsets_in )
			break;
		for ( ; c<offv->GetValue(i); c++ )
			it->push_back( datv->GetValue( c ) );
		while ( skip && it != data_out.end() )
		{
			++it;
			--skip;
		}
		if ( it == data_out.end() )
			break;
	}

	datv->Delete();
	offv->Delete();
	
	return count;
}

int vtkXMLShoeMeshReader::ReadArrayForCells( vtkXMLDataElement* da, vtkDataArray* outArray )
{
	vtkIdType startCell = 0;
	vtkIdType numCells = this->NumberOfCells[this->Piece];
	vtkIdType components = outArray->GetNumberOfComponents();
	return this->ReadArrayValues( da, startCell*components, outArray, 0, numCells*components);
}

int vtkXMLShoeMeshReader::ReadFunctionData( vtkFunctionData* data_out, vtkXMLDataElement* data_in, int num_fields_in )
{
	// FIXME: We currently allocate enough slots, but don't guarantee that the arrays
	// are ordered to match the vtkPointData.
	// For real XML files, point data and function data may not have the same field
	// orderings, so we should use array names instead of the \p field counter below.
	int field=0;
	int count=0;
	vtkXMLDataElement* extr = 0;

	for ( int e=0; e < data_in->GetNumberOfNestedElements(); e++ )
	{
		vtkXMLDataElement* ele = data_in->GetNestedElement( e );
		if ( (strcmp( ele->GetName(), "VaryingDataArray" ) == 0) && (field < num_fields_in) )
			count += this->ReadVaryingDataArray( data_out, field++, ele, this->GetNumberOfDofNodes() );
		else if ( strcmp( ele->GetName(), "Extrema" ) == 0 && !extr ) // only take the first extrema
			extr = ele;
	}

	const char* activeAttrib;
	if ( (activeAttrib = data_in->GetAttribute( "Scalars" )) != 0 )
		data_out->SetActiveScalars( activeAttrib );
	if ( (activeAttrib = data_in->GetAttribute( "Vectors" )) != 0 )
		data_out->SetActiveVectors( activeAttrib );
	if ( (activeAttrib = data_in->GetAttribute( "Tensors" )) != 0 )
		data_out->SetActiveTensors( activeAttrib );
	if ( (activeAttrib = data_in->GetAttribute( "TextureCoords" )) != 0 )
		data_out->SetActiveTCoords( activeAttrib );

	// OK, if there are extrema, try to read that
	if ( extr )
	{
		vtkFunctionData* extrema = data_out->GetExtrema();
		if ( ! extrema )
		{
			extrema = vtkFunctionData::New();
			extrema->MarkAsExtrema();
			extrema->AllocateOffsets( data_out->GetNumberOfArrays() );
			data_out->SetExtrema( extrema );
			extrema->FastDelete();
		}
		count += this->ReadFunctionData( extrema, extr, num_fields_in );
	}

	return count;
}

vtkIdType vtkXMLShoeMeshReader::ReadFreerangeData( vtkShoeMesh::ConnectivityType& data_out, vtkXMLDataElement* xml_node_in )
{
	vtkIdType count = 0;

	vtkXMLDataElement* xmlValues=0;
	vtkXMLDataElement* xmlHoles=0;
	vtkIdType val_size;
	vtkIdType capacity;
	vtkIdType holes_runlength;

	if ( !xml_node_in->GetScalarAttribute( "Size", val_size ) )
	{
		vtkErrorMacro( "FreeRange entry has a bad or missing Size attribute" );
		return 0;
	}

	if ( !xml_node_in->GetScalarAttribute( "Capacity", capacity ) )
	{
		vtkErrorMacro( "FreeRange entry has a bad or missing Capacity attribute" );
		return 0;
	}

	for ( int i=0; i<xml_node_in->GetNumberOfNestedElements(); i++ )
	{
		vtkXMLDataElement* ele = xml_node_in->GetNestedElement( i );
		if ( strcmp( ele->GetName(), "Values" ) == 0 )
			xmlValues = ele;
		else if ( strcmp( ele->GetName(), "Holes" ) == 0 )
			xmlHoles = ele;
	}

	if ( ! xmlValues )
	{
		vtkErrorMacro( "FreeRange has a bad or missing Values entry" );
		return 0;
	}

	if ( ! xmlHoles )
	{
		vtkErrorMacro( "FreeRange has a bad or missing Holes entry" );
		return 0;
	}

	if ( !xmlHoles->GetScalarAttribute( "RunLength", holes_runlength ) )
	{
		vtkErrorMacro( "FreeRange Hole entry has a bad or missing RunLength Attribute" );
		return 0;
	}
	
	data_out.clear();
	data_out.resize( capacity );

	count += this->ReadIdData( xmlValues, &(data_out[data_out.grab(val_size)]), val_size );

	vtkIdTypeArray* holes = vtkIdTypeArray::New();

	holes->SetNumberOfComponents( 1 );
	holes->SetNumberOfTuples( holes_runlength );
	count += this->ReadArrayValues( xmlHoles, 0, holes, 0, holes_runlength );

	int hsz = 1;
	for ( int h=0; h<holes_runlength; h++ )
	{
		int nh = holes->GetValue( h++ );
		for ( ; nh ; --nh )
			data_out.free( holes->GetValue( h++ ), hsz );
		hsz++;
	}

	holes->Delete();

	return count;
}

vtkIdType vtkXMLShoeMeshReader::ReadIdData( vtkXMLDataElement* da, vtkIdType* ptr, vtkIdType numValues )
{
  if ( !ptr )
    return 0;
  vtkIdType num = numValues;
  int result;
  if ( da->GetAttribute( "offset" ))
    {
    int offset = 0;
    da->GetScalarAttribute( "offset", offset );
    result = (this->XMLParser->ReadAppendedData( offset, (void*) ptr, 0, numValues, VTK_ID_TYPE ) == num );
    }
  else
    {
    int isAscii = 1;
    const char* fmt = da->GetAttribute( "format" );
    if ( fmt && ( strcmp( fmt, "binary" ) == 0 ) )
      isAscii = 0;
    result = (this->XMLParser->ReadInlineData( da, isAscii, (void*) ptr, 0, numValues, VTK_ID_TYPE ) == num );
    }
  return numValues;
}

vtkIdType vtkXMLShoeMeshReader::ReadCellDefinitions( vtkShoeMesh::CellDefsType& data_out, vtkstd::vector<vtkShoeMesh::CellDefsType::const_iterator>& cell_defs_out, vtkXMLDataElement* xml_node_in )
{
	vtkIdType count = 0;

	vtkIdType num_entries;

	if ( !xml_node_in->GetScalarAttribute( "NumberOfEntries", num_entries ) )
	{
		vtkErrorMacro( "FreeRange entry has a bad or missing NumberOfEntries attribute" );
		return 0;
	}

	data_out.clear();

	vtkIdTypeArray* ia = vtkIdTypeArray::New();

	int num_int_per_entry = 4 /*type & number*/ + 3*(1 /*geom*/ + this->NumberOfFunctionFields[ this->Piece ] /*funcs*/);
	ia->SetNumberOfComponents( num_int_per_entry );
	ia->SetNumberOfTuples( num_entries );
	count += this->ReadArrayValues( xml_node_in, 0, ia, 0, num_entries*num_int_per_entry );

	vtkIdType* iap = ia->GetPointer(0);
  int osize = 3*(1 + this->NumberOfFunctionFields[ this->Piece ]);
  int* orders = new int[ osize ];
	for ( int i=0; i<num_entries; i++, iap+=num_int_per_entry )
    {
    for ( int j=0; j<osize; ++j )
      {
      orders[j] = iap[4+j];
      }
		data_out.push_back( vtkCellDefinition( CellShape(iap[0]), PolyInterpolant(iap[1]), PolyProductSpace(iap[2]), iap[3],
		                    this->NumberOfFunctionFields[ this->Piece ], orders, orders + 3 ) );
    }

  delete [] orders;
	ia->Delete();

	for ( vtkShoeMesh::CellDefsType::const_iterator it = data_out.begin(); it != data_out.end(); ++it )
		cell_defs_out.push_back( it );

	return count;
}

vtkIdType vtkXMLShoeMeshReader::ReadCellSpecs( vtkShoeMesh::CellsType& data_out, vtkstd::vector<vtkShoeMesh::CellDefsType::const_iterator>& cell_defs_in, vtkXMLDataElement* xml_node_in, vtkIdType num_cells_in )
{
	vtkIdType count = 0;

	vtkXMLDataElement* xmlValues=0;
	vtkXMLDataElement* xmlHoles=0;
	vtkIdType val_size;
	vtkIdType capacity;
	vtkIdType num_entries;
	vtkIdType next_entry;
	vtkIdType holes_length;

	if ( !xml_node_in->GetScalarAttribute( "NumberOfEntries", num_entries ) )
	{
		vtkErrorMacro( "FreeRange entry has a bad or missing NumberOfEntries attribute" );
		return 0;
	}

	if ( !xml_node_in->GetScalarAttribute( "NextEntry", next_entry ) )
	{
		vtkErrorMacro( "FreeRange entry has a bad or missing NextEntry attribute" );
		return 0;
	}

	if ( !xml_node_in->GetScalarAttribute( "Capacity", capacity ) )
	{
		vtkErrorMacro( "FreeRange entry has a bad or missing Capacity attribute" );
		return 0;
	}

	for ( int i=0; i<xml_node_in->GetNumberOfNestedElements(); i++ )
	{
		vtkXMLDataElement* ele = xml_node_in->GetNestedElement( i );
		if ( strcmp( ele->GetName(), "Values" ) == 0 )
			xmlValues = ele;
		else if ( strcmp( ele->GetName(), "Holes" ) == 0 )
			xmlHoles = ele;
	}

	if ( ! xmlValues )
	{
		vtkErrorMacro( "FreeRange has a bad or missing Values entry" );
		return 0;
	}

	if ( ! xmlHoles )
	{
		vtkErrorMacro( "FreeRange has a bad or missing Holes entry" );
		return 0;
	}

	if ( !xmlValues->GetScalarAttribute( "Size", val_size ) )
	{
		vtkErrorMacro( "FreeRange Values entry has a bad or missing Size Attribute" );
		return 0;
	}

	if ( !xmlHoles->GetScalarAttribute( "Size", holes_length ) )
	{
		vtkErrorMacro( "FreeRange Hole entry has a bad or missing Size Attribute" );
		return 0;
	}

	vtkIdTypeArray* cell = vtkIdTypeArray::New();
	cell->SetNumberOfComponents(3);
	cell->SetNumberOfTuples( val_size/3 );

	count += this->ReadArrayValues( xmlValues, 0, cell, 0, val_size );

	data_out.clear();
	data_out.resize( capacity );

	vtkIdType* cv = cell->GetPointer(0);
	for ( int c=0; c<(val_size/3); c++, cv += 3 )
	{
		c = data_out.grab(); // this had better not change c -- we called data_out.clear()!
		data_out[c].Def = cell_defs_in[ cv[0] ]; // FIXME: This is an O(cv[0]), operation, not O(1)
		data_out[c].Offset = cv[1];
		data_out[c].SetCellPermutation( cv[2] );
	}

	cell->Delete();

	vtkIdTypeArray* holes = vtkIdTypeArray::New();

	holes->SetNumberOfComponents( 1 );
	holes->SetNumberOfTuples( holes_length );
	count += this->ReadArrayValues( xmlHoles, 0, holes, 0, holes_length );

	for ( int h=0; h<holes_length; h++ )
		data_out.free( holes->GetValue(h) );

	holes->Delete();

	return count;
}


