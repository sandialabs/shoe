/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2007-01-08 11:06:20 -0800 (Mon, 08 Jan 2007) $
  Version:   $Revision: 8642 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
// .NAME vtkXMLShoeMeshReader - Read Sandia Higher Order Element Mesh files.
// .SECTION Description
// vtkXMLShoeMeshReader reads an XML file containing a single vtkShoeMesh
// object. 
//
// The standard extension for this reader's file format is "shoe".
// This reader will also eventually be used to read a single piece
// for a parallel version of the reader.
//
// Since FunctionData is derived from DataSetAttributes, we look for
// XML attributes specifying the default scalars, vectors, tensors, and
// so on. The XML attribute tags we look for are named: Scalars, Vectors,
// Tensors, TextureCoords. We do <i>not</i> look for normals because
// analytical normals exist.

#ifndef vtkXMLShoeMeshReader_h
#define vtkXMLShoeMeshReader_h

#include <vtkXMLUnstructuredDataReader.h>

#include <vtksnlIOWin32Header.h>
#include <vtkShoeMesh.h>

class vtkFunctionData;
class vtkCellDefinition;

class VTK_SNL_IO_EXPORT vtkXMLShoeMeshReader
	: public vtkXMLUnstructuredDataReader
{
	public:
		vtkTypeRevisionMacro(vtkXMLShoeMeshReader,vtkXMLUnstructuredDataReader);
		void PrintSelf( ostream& os, vtkIndent indent );
		static vtkXMLShoeMeshReader* New();

		// Description:
		// Get or set the reader's output mesh.
		virtual void SetOutput( vtkShoeMesh* mesh_in );
		virtual vtkShoeMesh* GetOutput();
		virtual vtkIdType GetNumberOfCellsInPiece( int piece_in );

		virtual vtkIdType GetNumberOfDofNodes() const;

	protected:
		vtkXMLShoeMeshReader();
		virtual ~vtkXMLShoeMeshReader();

		const char* GetDataSetName();
		void GetOutputUpdateExtent( int& piece_out, int& num_pieces_out, int& ghost_level_out );
		void SetupOutputTotals();
		void SetupPieces( int num_pieces_in );
		void DestroyPieces();

    virtual int FillOutputPortInformation( int port, vtkInformation* info );
		void SetupOutputData();
		int ReadPiece( vtkXMLDataElement* ePiece );
		void SetupNextPiece();
		int ReadPieceData();

		//BTX
		// Description:
		// Set the \a n_in -th array of \a data_out to the values specified by \a xml_node_in.
		// The total number of entries in the list of offsets must be \a num_offsets_in.
		vtkIdType ReadVaryingDataArray( vtkFunctionData* data_out, int n_in, vtkXMLDataElement* xml_node_in, vtkIdType num_offsets_in );
		vtkIdType ReadVaryingDataArray( vtkShoeMesh::LinkType& data_out, vtkXMLDataElement* xml_node_in, vtkIdType num_offsets_in );

		// Description:
		// Read in data for the storage type freerange<vtkIdType> (currently only used for connectivity).
		// The size of the storage is included in the XML node and is not directly related to the
		// number of cells, points, or dof nodes present.
		vtkIdType ReadFreerangeData( vtkShoeMesh::ConnectivityType& data_out, vtkXMLDataElement* xml_node_in );

		// Description:
		// Read in cell definitions
		vtkIdType ReadCellDefinitions( vtkShoeMesh::CellDefsType& data_out, vtkstd::vector<vtkShoeMesh::CellDefsType::const_iterator>& cell_defs_out, vtkXMLDataElement* xml_node_in );
		
		// Description:
		// Read in cell definitions
		vtkIdType ReadCellSpecs( vtkShoeMesh::CellsType& data_out, vtkstd::vector<vtkShoeMesh::CellDefsType::const_iterator>& cell_defs_in, vtkXMLDataElement* xml_node_in, vtkIdType num_cells_in );
		
		virtual int ReadArrayForCells( vtkXMLDataElement* da, vtkDataArray* outArray );
		virtual int ReadFunctionData( vtkFunctionData* data_out, vtkXMLDataElement* data_in, int num_fields_in );

    // Description:
    // Low-level routine to read vtkIdType data from an XML file without a vtkArrayIterator.
    // This is called by ReadFreerangeData to read the connectivity.
    vtkIdType ReadIdData( vtkXMLDataElement* da, vtkIdType* ptr, vtkIdType numValues );

		// Description:
		// The total number of function data dof nodes in the mesh.
		vtkIdType TotalNumberOfDofNodes;
		// Description:
		// The number of dof nodes per piece.
		vtkstd::vector<vtkIdType> NumberOfDofNodes;

		// Description:
		// Pointers to the <FunctionData> xml tags containing field and geometric
		// dof coefficients, respectively.
		vtkstd::vector<vtkXMLDataElement*> FieldDofElements;
		vtkstd::vector<vtkXMLDataElement*> GeometryDofElements;

		vtkstd::vector<vtkIdType> NumberOfCells;
		vtkstd::vector<vtkXMLDataElement*> CellElements;

		vtkstd::vector<int> NumberOfFunctionFields;
		vtkstd::vector<int> NumberOfPointFields;
		//ETX
};

#endif // vtkXMLShoeMeshReader_h

