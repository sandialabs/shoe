/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2004-11-08 18:26:06 -0800 (Mon, 08 Nov 2004) $
  Version:   $Revision: 2875 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef vtkShoeMesh_h
#define vtkShoeMesh_h

// .NAME vtkShoeMesh - Sandia Higher Order Elements (SHOE) mesh representation for mixed-order nonlinear finite elements
// .SECTION Description
// vtkShoeMesh is a concrete implementation of vtkPointSet tailored to represent
// nonlinear finite elements with mixed-order interpolation and fast cell
// insertion/deletion for hp-adaptive refinement.
// <p>
//
// The vtkShoeMesh class introduces a new genre of field data to the data model.
// All objects that inherit from vtkDataSet have both point data (values defined
// at points) and cell data (values defined over the entire cell). Most filters
// assume that point data may be linearly interpolated over a cell. However,
// there is an increasing need for nonlinear maps ???
//
// Because fast deletion is a goal for this class, it does not store cells in a
// contiguous array like vtkUnstructuredGrid; instead, a freelist is used.
// A freelist is essentially an array of values with some entries marked as
// unused. The list of unused entries is stored in a second array.
// Because we allow blank entries in the middle of our list, cell deletion only
// requires the addition of an entry to the unused pool -- no movement of cells
// beyond the deleted cell is required. This also means that all cell ids are preserved
// when a deletion occurs.
// <p>
//
// This technique for storing cells has other ramifications. The most notable is
// that existing VTK filters assume that cells are numbered from 0 to GetNumberOfCells(),
// which is no longer the case. Before running any vtkPointSetToXXX or vtkDataSetToXXX
// filter on a vtkShoeMesh, you should call vtkShoeMesh::SqueezeCells(). This is not a fast
// operation, but it will remove blank entries from the list of Cells.
// <p>
//
// Finally, you should be aware that vtkShoeMesh does not support nonlinear cell queries
// through the vtkCell class. Filters that use the GetCell methods of this class will
// receive linear versions of the cells with no higher order information. If you need
// nonlinear cell information, use a vtkShoeMeshIterator. The functionality is similar
// but made available through a more efficient model that does not involve copying
// connectivity information for each cell (which can reach up to 27 entries per cell
// for a hexahedron compared to 8 entries for the linear case).
// <p>
//
// .SECTION See Also
// vtkShoeMeshIterator, freerange, freelist, vtkCell, vtkUnstructuredGrid

#ifndef	__STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS // explicitly request MIN/MAX macros as per C99 standard
#endif
#include <stdint.h>

#include <list>
#include <vector>
#include <stack>

#include <freelist>
#include <freerange>

#include <vtkPointSet.h>
#include <vtkIndent.h>

#include <vtksnlConfigure.h>
#include <vtksnlCommonWin32Header.h>
#include <vtkCellDefinition.h>
#include <vtkFunctionData.h>

class vtkDataSetAttributes;
class vtkIdTypeArray;
class vtkCell;

class vtkFunctionData;
class vtkShoeMeshIterator;

class VTK_SNL_COMMON_EXPORT vtkShoeMesh
	: public vtkPointSet
{
	public:
		vtkTypeRevisionMacro(vtkShoeMesh,vtkPointSet);
		void PrintSelf(ostream& os, vtkIndent indent);

		static vtkShoeMesh *New();

		//BTX
		/** A degree of freedom (DOF) node may contain edge, face, or volumetric mode coefficients.
		 */
		enum DofNodeType
		{
			Unused=-1,    //!< No cells refer to this node, so it is impossible to know its type
			EdgeNode=1,   //!< The coefficients at this node are for edge modes
			FaceNode=2,   //!< The coefficients at this node are for face modes
			VolumeNode=3  //!< The coefficients at this node are for volumetric modes
		};

		typedef vtkstd::vector<vtkIdType> LinkEntryType;
		typedef vtkstd::vector< LinkEntryType > LinkType;
		typedef freerange<vtkIdType,vtkIdType,-3> ConnectivityType;
		typedef vtkstd::list<vtkCellDefinition> CellDefsType;
		typedef CellDefsType::iterator CellDefIterator;
		typedef CellDefsType::const_iterator CellDefConstIterator;
		typedef vtkShoeMeshIterator const_iterator;
		/** Represent the specification of a cell
		 * A cell is specified by <ul>
		 * <li> its defining type and interpolation orders;
		 * <li> the permutations of its higher order node entries; and 
		 * <li> its offset into the connectivity array.
		 * </ul>
		 */
		class CellSpec {
			public:
				// Description:
				// Construct a default cell specification. A default cell spec is
				// equivalent to a unused cell.
				CellSpec();

				// Description:
				// A reference to the definition of the cell:<p><ul>
				// <li> its shape,
				// <li> its interpolant,
				// <li> its product space,
				// <li> the interpolation order for each geometric and function field,
				// <li> its operations (functions performed on a cell that vary by cell type)
				// </ul>
				CellDefConstIterator Def;

				// Description:
				// An index into the Connectivity array that defines which corner and
				// higher order "nodes" make up the cell.
				vtkIdType Offset;
				
				int  GetFacePermutation(int) const;
				void SetFacePermutation(int,int);

				bool GetEdgePermutation(int) const;
				void SetEdgePermutation(int,bool);

				void SetCellPermutation(uint32_t);
				uint32_t GetCellPermutation() const;
			protected:
				// Description:
				// A bit array of permutations for each of the edge and face nodes referenced
				// by this cell in the vtkShoeMesh::Connectivity array.
				//
				// Because, for example, a higher order face "node" may be shared by two volume
				// cells that do not share the same coordinate system, we must store how each
				// of the face "nodes" is oriented with respect to each volume cell.
				//
				// For edges, the values associated with a single higher order edge "node" may be
				// stored in forwards or backwards order.
				//
				// For faces, there are up to 8 possible permutations:<p><ul>
				// <li> the starting location may be any of 4 corner nodes on a brick element's face and
				// <li> the values may be in row-first or column-first order.
				// </ul>
				// 
				// To use NodePermutations, first ask for a particular face or edge's permutation.
				// Then, for each degree of freedom at a "node", call GetFaceIndex or GetEdgeIndex
				// to retrieve an offset into the "node"'s list of values.
				uint32_t NodePermutations;
		};

		typedef freelist<CellSpec,vtkIdType> CellsType;
		//ETX

		// Description:
		// Get/set the field data associated with higher order "nodes".
		int GetNumberOfNonlinearFunctions() const;
		vtkFunctionData* GetFunctionData();
		const vtkFunctionData* GetFunctionData() const;
		void SetFunctionData( vtkFunctionData* );

		// Description:
		// Get/set the geometric degrees of freedom associated with higher "nodes".
		vtkFunctionData* GetGeometryData();
		const vtkFunctionData* GetGeometryData() const;
		void SetGeometryData( vtkFunctionData* );

		// Description:
		// Get the number of higher order "nodes".
		// 
		// All of the vtkShoeMesh::FunctionData and vtkShoeMesh::GeometryData
		// Offsets arrays should have this many entries.
		vtkIdType GetNumberOfDofNodes() const;

		//BTX
		// Description:
		// Return the underlying links representation for higher order nodes.
		LinkType& GetDofNodeLinks();
		LinkType  GetDofNodeLinks() const;
		
		// Description:
		// Return the underlying links representation for corner nodes.
		LinkType& GetCornerLinks();
		LinkType  GetCornerLinks() const;

		// Description:
		// Return a pointer to the connectivity entries for a cell.
		vtkIdType* GetCellConnectivity( vtkIdType );
		const vtkIdType* GetCellConnectivity( vtkIdType ) const;
		ConnectivityType::iterator GetCellConnectivityBegin( vtkIdType );
		ConnectivityType::iterator GetCellConnectivityEnd( vtkIdType );
		ConnectivityType& GetConnectivity();
		const ConnectivityType& GetConnectivity() const;

		// Description:
		// Return iterator pointing to beginning/end of CellSpecs in this mesh.
		// Use of these functions is discouraged.  Please use Begin() and End().
		CellsType::iterator BeginCells();
		CellsType::iterator EndCells();
		CellsType::iterator GetCellSpec( vtkIdType cell_ID_in );
		CellsType::const_iterator GetCellSpec( vtkIdType cell_ID_in ) const;
		
		// Description:
		// Return iterator pointing to beginning/end of vtkCellDefinitions in this mesh.
		// Use of these functions is discouraged.  Please use vtkShoeMesh::FindCellDef.
		CellDefIterator BeginCellDefs();
		CellDefConstIterator BeginCellDefs() const;
		CellDefIterator EndCellDefs();
		CellDefConstIterator EndCellDefs() const;
		
		// Description:
		// Return iterator pointing to beginning/end of Cells in this mesh.
		const_iterator Begin();
		const_iterator End();
		//ETX
	
		// Description:
		// Free any unused memory at the tail end of the connectivity or the list of cells.
		virtual void SqueezeCells();
		virtual void SqueezeConnectivity();

		// Description:
		// Insert a DOF node at the tail of both the geometry and function data.
		vtkIdType InsertNextDofNode( int num_modes_this_node_in );

		//BTX
		// Description:
		// Return the type of modes stored at the given DOF node.
		// This function uses the DofNodeLinks lookup table to find a cell that refers to \a dof_node_id_in
		// and uses the position of the DOF node in that cell's connectivity array to discern the type of
		// node.
		// If there is no DofNodeLink for the given DOF node, then Unused is returned.
		// If the DOF node ID is out of bounds, then Unused is returned.
		DofNodeType GetDofNodeType( vtkIdType dof_node_id_in ) const;

		// Description:
		// Find the number of entries (modes) associated with each DOF node used by any cell in the mesh.
		// The number of modes per node depends on the order of the interpolant, which can vary from field
		// to field. You must specify a field id (or -1 for the geometry) for which the corresponding mode
		// sizes will be returned.
		//
		// For instance, consider a mesh containing 2 edge DOF nodes: the first node having 6 modes and
		// the second node having 3 modes for field 2. Calling \p GetDofNodeSizes( \a nodesize_out, 2 )
		// would <i>append</i> [ 6, 3 ] to \a nodesize_out.
		//
		// The number of modes per node is found by asking one of the cells listed in DofNodeLinks for
		// each node, so if a DOF node has become orphaned (i.e., DofNodeLinks is empty for that node),
		// a 0 will be entered.
		//
		// .SECTION Warning
		// This does <b>not</b> empty \a nodesize_out before appending the results. You are
		// responsible for calling \a nodesize_out.clear() before calling this member function.
		// (Actually, this behavior is a sneaky hack that will let you prepend a leading 0 to the
		// list of sizes... you can then add the sum of all previous entries to a given entry
		// and get a mapping between DOF node IDs and a global numbering of all mode coefficients.
		// This is used by the Analysis/SolidMechanics code to number modes in the stiffness matrix. Heh.)
		virtual void GetDofNodeSizes( vtkstd::vector<int>& nodesize_out, int field_id_in ) const;
		//ETX

		/**
		 * Accessing extrema of higher order functions
		 */
		//@{

		// Description:
		// Have extrema been calculated for the geometry?
		// This is a convenience routine that checks to see whether \a GetGeometryData()->GetExtrema()->GetValues(0) is non-null.
		bool DoGeometryExtremaExist();

		// Description:
		// Have extrema been calculated for the given DOF node?
		bool DoGeometryDofNodeExtremaExist( vtkIdType dof_node_id_in );

		//BTX
		// Description:
		// Insert the geometric extrema related to \a dof_node_id_in into the mesh.
		void InsertGeometryExtrema( vtkIdType dof_node_id_in, const vtkstd::vector< vtkstd::vector< double > >& extrema_in );
		//ETX

		// Description:
		// Get the geometric extrema over the DOF node \a dof_node_id_in .
		// The \a extrema_out array contains a record for each extremum.
		// Each record consists of a parametric coordinate (in the DOF node's coordinate system for permutation 0),
		// followed by the field value(s) at that parametric coordinate, followed by the type of extremum
		// (its index computed from the Hessian or 3 if the extremum is on the boundary of the k-cell associated
		// with the DOF node).
		void GetGeometryExtrema( double*& extrema_out, vtkIdType& num_extrema_this_node_out, vtkIdType dof_node_id_in );

		//BTX
		// Description:
		// Return the extrema values for a given DofNode and field.
		void GetPermutedDofNodeFunctionExtrema(vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, int field_in, vtkShoeMeshIterator& it_in);
		void GetPermutedDofNodeFunctionExtrema(vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, int field_in, vtkIdType cell_id);

		// Description:
		// Return the geometric extrema values for a given DofNode.
		void GetPermutedDofNodeGeometryExtrema(vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, vtkShoeMeshIterator& it_in);
		void GetPermutedDofNodeGeometryExtrema(vtkstd::vector<vtkstd::vector<double> >& extrema_out, vtkIdType dof_node_id_in, vtkIdType cell_id);
		//ETX

		// Description:
		// Have extrema been calculated for the specified field?
		// This is a convenience routine that checks to see whether \a GetFunctionData()->GetExtrema()->GetValues( \p field_in \a ) is non-null.
		bool DoFunctionExtremaExist( int field_in );

		// Description:
		// Have extrema of \a field_in been calculated for the given DOF node?
		bool DoFunctionDofNodeExtremaExist( int field_in, vtkIdType dof_node_id_in );

		//BTX
		// Description:
		// Insert the field \a field_in extrema related to \a dof_node_id_in into the mesh.
		void InsertFunctionExtrema( int field_in, vtkIdType dof_node_id_in, const vtkstd::vector< vtkstd::vector< double > >& extrema_in );
		//ETX

		// Description:
		// Get the extrema of a field over the DOF node \a dof_node_id_in .
		// The \a extrema_out array contains a record for each extremum.
		// Each record consists of a parametric coordinate (in the DOF node's coordinate system for permutation 0),
		// followed by the field value(s) at that parametric coordinate, followed by the type of extremum
		// (its index computed from the Hessian or 3 if the extremum is on the boundary of the k-cell associated
		// with the DOF node).
		void GetFunctionExtrema( double*& extrema_out, vtkIdType& num_extrema_this_node_out, vtkIdType field_in, vtkIdType dof_node_id_in );

		//@}

		//BTX
		// Description:
		// Find a cell definition.
		// FindCellDef will return EndCellDefs() if there is no match.
		CellDefConstIterator FindCellDef( const vtkCellType& type_in, const int* geom_order_in, const int* func_order_in ) const;

		// Description:
		// Find a cell definition.
		// FindOrCreateCellDef will create a cell definition and return it to you if there isn't one that
		// matches in the mesh.
		//
		// The variant that accepts an iterator to a vtkCellDefinition as input assumes that the iterator
		// refers to a vtkCellDefinition associated with some other mesh.
		CellDefIterator FindOrCreateCellDef( const vtkCellType& type_in, const int* geom_order_in, const int* func_order_in );
		CellDefIterator FindOrCreateCellDef( const CellDefConstIterator& celldef_in );
		CellDefIterator FindOrCreateCellDef( const vtkCellDefinition& celldef_in );

		// Description:
		// Insert a new cell.
		// If you don't have a CellDefIterator, use FindOrCreateCellDef to grab one.
		// The cell ID of the new cell is returned.
		vtkIdType InsertNextCell( const CellDefIterator& defn_in, const vtkIdType* connectivity_in, uint32_t node_permutations_in );
		//ETX

		// Description:
		// Insert a new field into the mesh. This will create arrays for both the corner nodes and DOF nodes.
		// Returns the ID of the newly created field (which you can use with \p GetValues(), \p GetOffsets(), or
		// \p GetPointData()->GetArray() to retrieve the new arrays).
		//
		// The DOF node offsets array will be sized correctly (to match the number of DOF nodes) and filled
		// with valid offsets into the values array.
		// The corner node array and the DOF node
		// values array will be sized correctly (memory allocated) but uninitialized.
		// If \a init_values_in is true, then the DOF node values array will be initialized to 0.
		// \a init_values_in defaults to false.
		//
		// The order of the newly created field will be set to \a order_in for all cell definitions. \a order_in
		// is also used to determine the number of modes per node that should be allocated in the DOF node values
		// array.
		virtual int CreateNewField( int data_type_in, int num_components_in, int* order_in, const char* field_name_in, bool init_values_in=false );

		// Description:
		// Each mesh can (<i>should</i>) have a short description associated with it.
		virtual void SetName( const char* );
		virtual const char* GetName() const ;

		// Description:
		// Each mesh can (<i>should</i>) have a short description associated with it.
		virtual void SetDescription( const char* );
		virtual const char* GetDescription() const ;

		// Description:
		// Reset a mesh to its initial state.
		// This will not free memory in connectivity and links lists, but will make the mesh appear completely empty.
		virtual void Reset();

		//@{
		/**
		 * Pure virtual vtkDataObject Methods
		 */
		virtual void SetUpdateExtent( int piece_in, int num_pieces_in, int ghost_level_in );
		//@}

		virtual void GetUpdateExtent( int &piece_in, int &num_pieces_in, int &ghost_level_in );

		// Description:
		// Call superclass method to avoid hiding
		// Since this data type does not use 3D extents, this set method
		// is useless but necessary since vtkDataSetToDataSetFilter does not
		// know what type of data it is working on.
		void SetUpdateExtent( int x1, int x2, int y1, int y2, int z1, int z2 );
		void SetUpdateExtent( int ext[6] );

		//@{
		/**
		 * Pure virtual vtkDataSet methods
		 */
		virtual vtkIdType GetNumberOfCells();
		virtual vtkCell* GetCell(vtkIdType);
		virtual void GetCell(vtkIdType, vtkGenericCell*);
		virtual int GetCellType(vtkIdType);
		virtual void GetCellPoints(vtkIdType, vtkIdList*);
		virtual void GetPointCells(vtkIdType, vtkIdList*);
		virtual int GetMaxCellSize();
		//@}

		//@{
		/**
		 * Virtual vtkDataSet methods we override.
		 */
		virtual void Squeeze();
		//@}
		
		//@{
		/**
		 * Virtual vtkPointSet methods we override.
		 */
		virtual void SetPoints( vtkPoints* points_in );
		//@}

	protected:
		//BTX
		vtkShoeMesh();
		virtual ~vtkShoeMesh();

		// Description:
		// This member stores extra degrees of freedom for higher order "nodes" for each field.
		vtkFunctionData* FunctionData;

		// Description:
		// This member stores extra geometric degrees of freedom associated with higher order "nodes".
		vtkFunctionData* GeometryData;

		// Description:
		// This member stores connectivity information for all the cells in the mesh.
		// Each entry is either:<p><ul>
		// <li> a corner node index (and thus an offset into a Points->Data or PointData->Data array), or
		// <li> a higher order node index (and thus an offset into a GeometryData->Offsets or FunctionData->Offsets array).
		// </ul>
		// The cell type (specified in the Cells array) determines the number of each type of entry associated
		// with a particular cell.
		ConnectivityType Connectivity;

		// Description:
		// Each cell has an entry in this array specifying<p><ul>
		// <li> the cell type (through a reference to a CellDefinition),
		// <li> its offset into the Connectivity array, and
		// <li> permutations of any higher order nodal values.
		// </ul>
		CellsType Cells;

		// Description:
		// This lists all the CellDefinitions present in the mesh.
		// A CellDefinition defines a cell's<p><ul>
		// <li> shape,
		// <li> interpolant,
		// <li> product space,
		// <li> interpolation order for each geometric and function field,
		// <li> operations (functions performed on a cell that vary by cell type)
		// </ul>
		CellDefsType CellDefs;

		// Description:
		// We need to keep the number of higher order "nodes" somewhere because
		// it occurs in too many places for there to be a clear, canonical reference.
		// Any time the FunctionData or GeometryData is set, we should verify that
		// every Offsets array is of the same length and that length is NumberOfDofNodes.
		vtkIdType NumberOfDofNodes;

		// Description:
		// A reverse lookup table that lists all the cells associated with a given corner node ID.
		LinkType CornerLinks;

		// Description:
		// A reverse lookup table that lists all the cells associated with a given higher order node ID.
		LinkType DofNodeLinks;

		// Description:
		// The name of the mesh
		char* Name;

		// Description:
		// A short description of the mesh
		char* Description;

		// Description:
		// Set the number of higher order nodes.
		// Each vtkFunctionData::Offsets array must have a length equal to \p n.
		virtual void SetNumberOfDofNodes( vtkIdType n );

		friend class vtkXMLShoeMeshReader;
		friend class vtkXMLShoeMeshWriter;

	private:
		vtkShoeMesh(const vtkShoeMesh&);  // Not implemented.
		void operator=(const vtkShoeMesh&);  // Not implemented.
		//ETX
};

//BTX
inline vtkIdType vtkShoeMesh::GetNumberOfCells() { return this->Cells.size(); }
inline vtkIdType vtkShoeMesh::GetNumberOfDofNodes() const { return this->NumberOfDofNodes; }
inline vtkShoeMesh::LinkType& vtkShoeMesh::GetDofNodeLinks() { return this->DofNodeLinks; }
inline vtkShoeMesh::LinkType  vtkShoeMesh::GetDofNodeLinks() const { return this->DofNodeLinks; }
inline vtkShoeMesh::LinkType& vtkShoeMesh::GetCornerLinks() { return this->CornerLinks; }
inline vtkShoeMesh::LinkType  vtkShoeMesh::GetCornerLinks() const { return this->CornerLinks; }
inline uint32_t vtkShoeMesh::CellSpec::GetCellPermutation() const { return this->NodePermutations; }
inline void vtkShoeMesh::CellSpec::SetCellPermutation(uint32_t p) { this->NodePermutations = p; }
inline vtkFunctionData* vtkShoeMesh::GetGeometryData() { return this->GeometryData; }
inline const vtkFunctionData* vtkShoeMesh::GetGeometryData() const { return this->GeometryData; }
inline vtkFunctionData* vtkShoeMesh::GetFunctionData() { return this->FunctionData; }
inline const vtkFunctionData* vtkShoeMesh::GetFunctionData() const { return this->FunctionData; }
inline vtkShoeMesh::CellDefIterator vtkShoeMesh::BeginCellDefs() { return this->CellDefs.begin(); }
inline vtkShoeMesh::CellDefConstIterator vtkShoeMesh::BeginCellDefs() const { return this->CellDefs.begin(); }
inline vtkShoeMesh::CellDefIterator vtkShoeMesh::EndCellDefs() { return this->CellDefs.end(); }
inline vtkShoeMesh::CellDefConstIterator vtkShoeMesh::EndCellDefs() const { return this->CellDefs.end(); }

inline vtkIdType* vtkShoeMesh::GetCellConnectivity( vtkIdType id )
{ return ( (id < 0) || (id > this->Cells.size()) ) ? 0 : &(this->Connectivity[ this->Cells[id].Offset ]) ; }
inline const vtkIdType* vtkShoeMesh::GetCellConnectivity( vtkIdType id ) const
{ return ( (id < 0) || (id > this->Cells.size()) ) ? 0 : &(this->Connectivity[ this->Cells[id].Offset ]) ; }
inline vtkShoeMesh::ConnectivityType& vtkShoeMesh::GetConnectivity() { return this->Connectivity; }
inline const vtkShoeMesh::ConnectivityType& vtkShoeMesh::GetConnectivity() const { return this->Connectivity; }

inline vtkShoeMesh::CellsType::iterator vtkShoeMesh::GetCellSpec( vtkIdType cell_ID_in ) { return CellsType::iterator( &this->Cells, cell_ID_in ); }
inline vtkShoeMesh::CellsType::const_iterator vtkShoeMesh::GetCellSpec( vtkIdType cell_ID_in ) const
{ return CellsType::const_iterator( &this->Cells, cell_ID_in ); }
inline void vtkShoeMesh::SetDescription( const char* d ) { if ( d == this->Description ) return; if ( this->Description ) delete [] this->Description; if ( d ) { this->Description = new char[ strlen(d) + 1 ]; strcpy( this->Description, d ); } }
inline const char* vtkShoeMesh::GetDescription() const { return this->Description; }
inline void vtkShoeMesh::SetName( const char* d ) { if ( d == this->Name ) return; if ( this->Name ) delete [] this->Name; if ( d ) { this->Name = new char[ strlen(d) + 1 ]; strcpy( this->Name, d ); } }
inline const char* vtkShoeMesh::GetName() const { return this->Name; }

inline void vtkShoeMesh::GetUpdateExtent( int& piece_in, int& num_pieces_in, int& ghost_level_in )
{
	piece_in = this->GetUpdatePiece();
	num_pieces_in = this->GetUpdateNumberOfPieces();
	ghost_level_in = this->GetUpdateGhostLevel();
}
inline void vtkShoeMesh::SetUpdateExtent( int x1, int x2, int y1, int y2, int z1, int z2 )
{ this->Superclass::SetUpdateExtent( x1, x2, y1, y2, z1, z2 ); };
inline void vtkShoeMesh::SetUpdateExtent( int ext[6] ) { this->Superclass::SetUpdateExtent( ext ); };
//ETX

#endif // vtkShoeMesh_h
