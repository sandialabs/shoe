/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2005-02-23 13:37:58 -0800 (Wed, 23 Feb 2005) $
  Version:   $Revision: 3807 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
// .NAME vtkXMLShoeMeshWriter - Write Sandia Higher Order Element Mesh files.
// .SECTION Description
// vtkXMLShoeMeshWriter reads an XML file containing a single vtkShoeMesh
// object. 
// The standard extension for this reader's file format is "shoe".
// This reader will also eventually be used to read a single piece
// for a parallel version of the reader.
//
// .SECTION See Also
// vtkXMLPShoeMeshWriter

#ifndef vtkXMLShoeMeshWriter_h
#define vtkXMLShoeMeshWriter_h

#include <map>

#include <vtkXMLUnstructuredDataWriter.h>

#include <vtksnlIOWin32Header.h>
#include <vtkShoeMesh.h>

class vtkOutputStream;

class vtkFunctionData;
class vtkCellDefinition;

class VTK_SNL_IO_EXPORT vtkXMLShoeMeshWriter
	: public vtkXMLUnstructuredDataWriter
{
	public:
		vtkTypeRevisionMacro(vtkXMLShoeMeshWriter,vtkXMLUnstructuredDataWriter);
		void PrintSelf( ostream& os, vtkIndent indent );
		static vtkXMLShoeMeshWriter* New();

		// Description:
		// The default file extension for this type of file.
		// Required by vtkXMLWriter.
		// (Returns ".shoe")
		virtual const char* GetDefaultFileExtension();

		// Description:
		// Get or set the reader's output mesh.
		virtual void SetInput( vtkShoeMesh* mesh_in );
		virtual vtkShoeMesh* GetInput();
		
		virtual void SetStream( ostream* stream_in );

	protected:
		//BTX

		vtkXMLShoeMeshWriter();
		virtual ~vtkXMLShoeMeshWriter();

		// Description:
		// Returns a string used as the XML tag name for the dataset.
		// This class returns "ShoeMesh".
		const char* GetDataSetName();

		// Description:
		// Set information about what part of the dataset is to be written.
		void SetInputUpdateExtent( int piece_in, int num_pieces_in, int ghost_level_in );

		// Description:
		// Versioning for the ShoeMesh file type.
		// Could be used by the reader to switch between multiple formats, but I don't see
		// support in the reader classes yet.
		virtual int GetDataSetMajorVersion();
		virtual int GetDataSetMinorVersion();

		// Description:
		// Write the file with inline or appended data.
		// Only one of these two methods will be called per write.
		//virtual void WriteInlineMode( vtkIndent indent );
		//virtual void WriteAppendedMode( vtkIndent indent );

		// Description:
		// Helper routines for the \p WriteInlineMode and \p WriteAppendedMode methods.
		virtual void WriteInlinePieceAttributes();
		virtual void WriteInlinePiece( vtkIndent indent );

		virtual void WriteAppendedPieceAttributes( int index );
		virtual void WriteAppendedPiece( int index, vtkIndent indent );
		virtual void WriteAppendedPieceData( int index );

		// Description:
		// Return the number of cells in the mesh.
		// A pure virtual method inherited from vtkXMLUnstructuredDataWriter.
		virtual vtkIdType GetNumberOfInputCells();

		// Description:
		// Write all the fields in a vtkFunctionData object.
		virtual void WriteFunctionDataInline( vtkIndent indent_in, vtkFunctionData* data_in, const char* type_in );

		// Description:
		// Set the \a n_in -th array of \a data_out to the values specified by \a xml_node_in.
		// The total number of entries in the list of offsets must be \a num_offsets_in.
		void WriteVaryingDataArrayInline( vtkIndent indent_in, vtkFunctionData* data_in, int n_in );
		void WriteVaryingDataArrayInline( vtkIndent indent_in, const char* name_in, vtkShoeMesh::LinkType& data_out );

		// Description:
		// Write in data for the storage type freerange<vtkIdType> (currently only used for connectivity).
		// The size of the storage is included in the XML node and is not directly related to the
		// number of cells, points, or dof nodes present.
		void WriteConnectivityInline( vtkIndent indent_in );

		// Description:
		// Write cell definitions
		void WriteCellDefinitionsInline( vtkIndent indent_in );
		
		// Description:
		// Write cell specifications
		void WriteCellSpecsInline( vtkIndent indent_in );

    // Description:
    // Force all inputs to be of type vtkShoeMesh. Duh.
    virtual int FillInputPortInformation( int port, vtkInformation* info );

		vtkstd::map<vtkShoeMesh::CellDefConstIterator,int> CellDefSwizzle;

		//ETX
};

//BTX

inline const char* vtkXMLShoeMeshWriter::GetDefaultFileExtension() { return ".shoe"; }
inline const char* vtkXMLShoeMeshWriter::GetDataSetName() { return "ShoeMesh"; }
inline int vtkXMLShoeMeshWriter::GetDataSetMajorVersion() { return 0; }
inline int vtkXMLShoeMeshWriter::GetDataSetMinorVersion() { return 1; }

//ETX

#endif // vtkXMLShoeMeshWriter_h

