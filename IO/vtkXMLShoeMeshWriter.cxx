/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtkObjectFactory.h>
#include <vtkDataArray.h>
#include <vtkIdTypeArray.h>
#include <vtkOutputStream.h>

#include <vtkCellDefinition.h>
#include <vtkExecutive.h>
#include <vtkXMLShoeMeshWriter.h>
#include <vtkInformation.h>



vtkCxxRevisionMacro(vtkXMLShoeMeshWriter, "$Revision: 7653 $");
vtkStandardNewMacro(vtkXMLShoeMeshWriter);

vtkXMLShoeMeshWriter::vtkXMLShoeMeshWriter()
{
}

vtkXMLShoeMeshWriter::~vtkXMLShoeMeshWriter()
{
}

void vtkXMLShoeMeshWriter::PrintSelf( vtkstd::ostream& os, vtkIndent indent )
{
	this->Superclass::PrintSelf( os, indent );
	os << indent << "CellDefSwizzle: " << this->CellDefSwizzle.size() << " entries" << vtkstd::endl;
}

void vtkXMLShoeMeshWriter::SetInput( vtkShoeMesh* mesh_in )
{
	this->Superclass::SetInput( 0, mesh_in );
}

vtkShoeMesh* vtkXMLShoeMeshWriter::GetInput()
{
  if ( this->GetNumberOfInputConnections(0) < 1 )
    return 0;

  return vtkShoeMesh::SafeDownCast(this->GetExecutive()->GetInputData( 0, 0 ));
}

void vtkXMLShoeMeshWriter::SetStream( vtkstd::ostream* stream_in )
{
	this->Stream =  stream_in ;
}

void vtkXMLShoeMeshWriter::SetInputUpdateExtent( int piece_in, int num_pieces_in, int ghost_level_in )
{
	this->GetInput()->SetUpdateExtent( piece_in, num_pieces_in, ghost_level_in );
}

#if 0
void vtkXMLShoeMeshWriter::WriteInlineMode( vtkIndent indent )
{
}
#endif // 0

void vtkXMLShoeMeshWriter::WriteInlinePieceAttributes()
{
	this->Superclass::WriteInlinePieceAttributes();
	vtkShoeMesh* input = this->GetInput();
	this->WriteScalarAttribute( "NumberOfCells", input->GetNumberOfCells() );
	this->WriteScalarAttribute( "NumberOfDofNodes", input->GetNumberOfDofNodes() );
}

void vtkXMLShoeMeshWriter:: WriteInlinePiece(vtkIndent indent)
{
	this->Superclass::WriteInlinePiece(indent);

	vtkShoeMesh* input = this->GetInput();
	this->WriteFunctionDataInline( indent, input->GetGeometryData(), "Geometry" );
	this->WriteFunctionDataInline( indent, input->GetFunctionData(), "Fields" );

	vtkstd::ostream& os = *(this->Stream);
	os << indent << "<Cells>" << vtkstd::endl;
	vtkIndent cellIndent = indent.GetNextIndent();

	this->WriteConnectivityInline( cellIndent );
	this->WriteCellDefinitionsInline( cellIndent );
	this->WriteCellSpecsInline( cellIndent );
	this->WriteVaryingDataArrayInline( cellIndent, "CornerLinks", input->GetCornerLinks() );
	this->WriteVaryingDataArrayInline( cellIndent, "DofNodeLinks", input->GetDofNodeLinks() );

	os << indent << "</Cells>" << vtkstd::endl;
}



void vtkXMLShoeMeshWriter::WriteAppendedPieceAttributes(int index)
{
}

void vtkXMLShoeMeshWriter::WriteAppendedPiece(int index, vtkIndent indent)
{
}

void vtkXMLShoeMeshWriter::WriteAppendedPieceData(int index)
{
}


vtkIdType vtkXMLShoeMeshWriter::GetNumberOfInputCells()
{
	vtkShoeMesh* input = this->GetInput();
	if ( ! input )
	{
		vtkWarningMacro( "GetNumberOfInputCells called without an input dataset (or input of wrong type)" );
		return 0;
	}
	return input->GetNumberOfCells();
}

void vtkXMLShoeMeshWriter::WriteFunctionDataInline( vtkIndent indent, vtkFunctionData* data_in, const char* type_in )
{
	vtkShoeMesh* input = this->GetInput();
	vtkstd::ostream& os = *(this->Stream);
	bool is_extrema = data_in->IsExtrema();
	if ( is_extrema )
	{
		os << indent << "<Extrema>" << vtkstd::endl;
	}
	else
	{
		os << indent << "<FunctionData";
		os << " Type=\"" << type_in << "\"";
		// Here's where you might call this->WriteInlineFunctionDataAttributes( data_in ); so people
		// could subclass vtkFunctionData...
		os << ">" << vtkstd::endl;
	}
	vtkIndent funcIndent = indent.GetNextIndent();

	for ( int i=0; i<input->GetNumberOfNonlinearFunctions(); ++i )
		if ( data_in->GetOffsets( i ) && (data_in->GetOffsets(i)->GetSize() > 0) )
			this->WriteVaryingDataArrayInline( funcIndent, data_in, i );

	if ( (! is_extrema) && data_in->GetExtrema() )
		this->WriteFunctionDataInline( funcIndent, data_in->GetExtrema(), "Extrema" );

	if ( is_extrema )
		os << indent << "</Extrema>" << vtkstd::endl;
	else
		os << indent << "</FunctionData>" << vtkstd::endl;
}


void vtkXMLShoeMeshWriter::WriteVaryingDataArrayInline( vtkIndent indent, vtkFunctionData* data, int n )
{
	//vtkShoeMesh* input = this->GetInput();
	vtkstd::ostream& os = *(this->Stream);
	vtkIdTypeArray* off = data->GetOffsets(n);
	vtkDataArray* val = data->GetValues(n);
	if ( ! off || ! val )
		return;
	os << indent << "<VaryingDataArray";
	if ( val->GetName() )
		this->WriteStringAttribute( "Name", val->GetName() );
	this->WriteScalarAttribute( "Id", n );
	this->WriteWordTypeAttribute( "Type", val->GetDataType() );
	this->WriteScalarAttribute( "NumberOfComponents", val->GetNumberOfComponents() );
	this->WriteDataModeAttribute( "format" );
	os << ">" << vtkstd::endl;

	vtkIndent arrIndent = indent.GetNextIndent();

	vtkIdType n_entries = val->GetNumberOfTuples()*val->GetNumberOfComponents();
	os << arrIndent << "<Values";
	this->WriteScalarAttribute( "Size", n_entries );
	os << ">" << vtkstd::endl;
	this->WriteInlineData( val, arrIndent.GetNextIndent() );
	os << arrIndent << "</Values>" << vtkstd::endl;

	os << arrIndent << "<Offsets>" << vtkstd::endl;
	this->WriteInlineData( off, arrIndent.GetNextIndent() );
	os << arrIndent << "</Offsets>" << vtkstd::endl;

	os << indent << "</VaryingDataArray>" << vtkstd::endl;
}

void vtkXMLShoeMeshWriter::WriteVaryingDataArrayInline( vtkIndent indent, const char* name, vtkShoeMesh::LinkType& data )
{
	vtkstd::ostream& os = *(this->Stream);

	os << indent << "<VaryingDataArray";
	this->WriteStringAttribute( "Name", name );
	this->WriteScalarAttribute( "NumberOfComponents", 1 );
	switch ( this->GetIdType() )
	{
		case vtkXMLWriter::Int32:
			this->WriteStringAttribute( "Type", "Int32" );
			break;
		case vtkXMLWriter::Int64:
			this->WriteStringAttribute( "Type", "Int64" );
			break;
		default:
			os << "/>" << vtkstd::endl;
			vtkErrorMacro( "Undefined vtkIdType storage type" );
			return;
	}
	this->WriteDataModeAttribute( "format" );
	os << ">" << vtkstd::endl;

	vtkIndent arrIndent = indent.GetNextIndent();

	vtkIdType nval = 0;
	for ( vtkShoeMesh::LinkType::const_iterator lit = data.begin(); lit != data.end(); ++lit )
		nval += lit->size();
	os << arrIndent << "<Values";
	this->WriteScalarAttribute( "Size", nval );
	os << ">" << vtkstd::endl << arrIndent.GetNextIndent();
	for ( vtkShoeMesh::LinkType::const_iterator lit = data.begin(); lit != data.end(); ++lit )
		copy( lit->begin(), lit->end(), vtkstd::ostream_iterator<vtkIdType>( os, " " ) );
	os << vtkstd::endl << arrIndent << "</Values>" << vtkstd::endl;

	nval = 0;
	os << arrIndent << "<Offsets>" << vtkstd::endl << arrIndent.GetNextIndent();
	for ( vtkShoeMesh::LinkType::const_iterator lit = data.begin(); lit != data.end(); ++lit )
	{
		if ( lit->size() )
			os << nval << " ";
		else
			os << "-1 ";
		nval += lit->size();
	}
	os << vtkstd::endl << arrIndent << "</Offsets>" << vtkstd::endl;

	os << indent << "</VaryingDataArray>" << vtkstd::endl;
}

void vtkXMLShoeMeshWriter::WriteConnectivityInline( vtkIndent indent )
{
	vtkShoeMesh* input = this->GetInput();
	vtkstd::ostream& os = *(this->Stream);

	os << indent << "<FreeRange";
	this->WriteStringAttribute( "Name", "Connectivity" );
	switch ( this->GetIdType() ) {
		case vtkXMLWriter::Int32:
			this->WriteStringAttribute( "Type", "Int32" );
			break;
		case vtkXMLWriter::Int64:
			this->WriteStringAttribute( "Type", "Int64" );
			break;
		default:
			vtkErrorMacro( "Unknown vtkIdType specified for output" );
			this->WriteStringAttribute( "Type", "***ERROR***" );
			os << "/>" << vtkstd::endl;
			return;
	}
	this->WriteScalarAttribute( "EmptyEntry", input->GetConnectivity().empty_entry_value() );
	this->WriteScalarAttribute( "Size", input->GetConnectivity().max_id() + 1 );
	this->WriteScalarAttribute( "Capacity", input->GetConnectivity().capacity() );
	this->WriteDataModeAttribute( "format" );
	os << ">" << vtkstd::endl;

	vtkShoeMesh::ConnectivityType& conn( input->GetConnectivity() );

	vtkIndent arrIndent = indent.GetNextIndent();
	os << arrIndent << "<Values>" << vtkstd::endl;
	os << arrIndent.GetNextIndent();
	for ( int i = 0; i <= conn.max_id(); ++i )
		os << conn[i] << " ";
	os << vtkstd::endl << arrIndent << "</Values>" << vtkstd::endl;

	os << arrIndent << "<Holes";
	this->WriteScalarAttribute( "MaxDeadSize", conn.max_chunked_grab_size() );
	vtkIdType dead_total = 0;
	vtkShoeMesh::ConnectivityType::dead_list_t::iterator it;
	for ( it = conn.dead_begin(); it != conn.dead_end(); ++it )
		dead_total += it->size()+1;
	this->WriteScalarAttribute( "RunLength", dead_total );
	os << ">" << vtkstd::endl << arrIndent.GetNextIndent() << "0 ";
	for ( it = conn.dead_begin(); it != conn.dead_end(); ++it )
	{
		os << it->size() << " ";
		for ( vtkShoeMesh::ConnectivityType::dead_list_entry_t::iterator hit = it->begin(); hit != it->end(); ++hit )
			os << *hit << " ";
	}
	os << vtkstd::endl << arrIndent << "</Holes>" << vtkstd::endl;

	os << indent << "</FreeRange>" << vtkstd::endl;
}

// This is a hack so that the CellDefSwizzle map can work... We simply compare
// two iterators by comparing the addresses of the CellDefinition objects they
// reference. Who's yer daddy!
bool operator < ( const vtkShoeMesh::CellDefConstIterator& a, const vtkShoeMesh::CellDefConstIterator& b )
{
	return &(*a) < &(*b) ;
}

void vtkXMLShoeMeshWriter::WriteCellDefinitionsInline( vtkIndent indent )
{
	vtkShoeMesh* input = this->GetInput();
	vtkstd::ostream& os = *(this->Stream);

	this->CellDefSwizzle.clear();
	int ndefs = 0;
	for ( vtkShoeMesh::CellDefConstIterator it = input->BeginCellDefs(); it != input->EndCellDefs(); ++it )
		++ndefs;
	os << indent << "<List";
	this->WriteStringAttribute( "Name", "CellDefinitions" );
	this->WriteStringAttribute( "Type", "IdType*13" ); // FIXME: HACK for now.
	this->WriteScalarAttribute( "NumberOfEntries", ndefs );
	this->WriteDataModeAttribute( "format" );
	os << ">" << vtkstd::endl;

	int cd = 0;
	vtkIndent itemIndent = indent.GetNextIndent();	
	for ( vtkShoeMesh::CellDefConstIterator it = input->BeginCellDefs(); it != input->EndCellDefs(); ++it )
	{
		this->CellDefSwizzle.insert( vtkstd::pair<vtkShoeMesh::CellDefConstIterator,int>( it, cd++ ) );
		os << itemIndent << it->GetShape() << " " << it->GetInterpolant() << " " << it->GetProductSpace()
		   << "   " << it->GetNumberOfRefs() << "   ";
		for ( int i=0; i<3; ++i )
			os << it->GetGeometricOrder()[i] << " ";
		os << "  ";
		for ( int f=0; f<it->GetNumberOfFields(); ++f )
		{
			for ( int j=0; j<3; ++j )
				os << it->GetFieldOrder(f)[j] << " ";
			os << "  ";
		}
		os << vtkstd::endl;
	}

	os << indent << "</List>" << vtkstd::endl;
}

void vtkXMLShoeMeshWriter::WriteCellSpecsInline( vtkIndent indent )
{
	vtkShoeMesh* input = this->GetInput();
	vtkstd::ostream& os = *(this->Stream);

	os << indent << "<FreeList";
	this->WriteStringAttribute( "Name", "CellSpecs" );
	this->WriteStringAttribute( "Type", "Int32,Int32,Int32" ); // FIXME: HACK for now
	this->WriteScalarAttribute( "NumberOfEntries", input->Cells.size() );
	this->WriteScalarAttribute( "NextEntry", input->Cells.max_used_entry() + 1 );
	this->WriteScalarAttribute( "Capacity", input->Cells.capacity() );
	this->WriteDataModeAttribute( "format" );
	os << ">" << vtkstd::endl;

	vtkIndent arrIndent = indent.GetNextIndent();
	os << arrIndent << "<Values";
	this->WriteScalarAttribute( "Size", (input->Cells.max_used_entry() + 1 + input->Cells.unused_size())*3 );
	os << ">" << vtkstd::endl;
	os << arrIndent.GetNextIndent();
	vtkIdType cur = 0;
	for ( vtkShoeMesh::CellsType::iterator it = input->BeginCells(); it != input->EndCells(); ++it )
	{
		for ( vtkIdType tmp=cur; tmp<it.where(); ++tmp )
			os << "0 0 0   "; // Write out entries for holes
		cur = it.where()+1;
		vtkstd::map<vtkShoeMesh::CellDefConstIterator,int>::iterator swiz = this->CellDefSwizzle.find( it->Def );
		if ( swiz == this->CellDefSwizzle.end() )
		{
			os << " ***ERROR*** </Values></FreeList>" << vtkstd::endl;
			vtkErrorMacro( "Bad cell definition swizzle. File will be unreadable." );
			return;
		}
		os << swiz->second << " " << it->Offset << " " << it->GetCellPermutation() << "   ";
	}
	os << vtkstd::endl << arrIndent << "</Values>" << vtkstd::endl;

	os << arrIndent << "<Holes";
	this->WriteScalarAttribute( "Size", input->Cells.unused_size() );
	os << ">" << vtkstd::endl;
	os << arrIndent.GetNextIndent();
	copy( input->Cells.unused_begin(), input->Cells.unused_end(), vtkstd::ostream_iterator<vtkIdType>( os, " " ) );
	os << vtkstd::endl << arrIndent << "</Holes>" << vtkstd::endl;

	os << indent << "</FreeList>" << vtkstd::endl;
}

int vtkXMLShoeMeshWriter::FillInputPortInformation( int port, vtkInformation* info )
{
  info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkShoeMesh" );
  info->Set( vtkAlgorithm::INPUT_IS_REPEATABLE(), 1 );
  info->Set( vtkAlgorithm::INPUT_IS_OPTIONAL(), port > 0 ? 1 : 0 );

  return 1;
}

