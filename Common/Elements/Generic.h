/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2007-01-11 15:03:33 -0800 (Thu, 11 Jan 2007) $
  Version:   $Revision: 8688 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef Elements_Generic_h
#define Elements_Generic_h

/**
 * Routines used across multiple cell shapes.
 */
class vtkShoeElemGeneric
{
	public:
		static const int faceIds[][3];
		static const int edgeConnectivity[][5];
		static const int faceConnectivity[][8];
		static const int faceConstantPara[][2];

		/** Integration points for the Legendre shape functions.
		 * \a GaussianIntegrationAbcissas and \a GaussianIntegrationWeights contain the
		 * constants needed to perform Gaussian quadrature.
		 *
		 * There are \a n entries in the \a n -th array. There are no entries in
		 * \a GaussianIntegrationAbcissas[0] or \a GaussianIntegrationWeights[0].
		 */
		//@{
		static const double* GaussianIntegrationAbcissas[];
		static const double* GaussianIntegrationWeights[];
		//@}

		static int ZeroModesPerNode( int, const int[3] );
		static int ZeroModes( const int[3] );
		static int OneModePerNode( int, const int[3] );
		static int OneMode( const int[3] );
		static int GetNumberOfEdgeModesPerNodePLeg( int ID, const int order[4] );

		// Description:
		// Convert a set of parametric coordinates from their DOF node storage order
		// into the current coordinates.
		// @param parametric_coords is a pointer to parametric coordinates to be
		//   converted. Only the first 2 entries are accessed for a quadrilateral
		//   element.
		// @param node_permutation is the permutation of the storage <i>relative</i>
		//   to the current coordinates.
		static double* DepermuteQuadParameters( double* parametric_coords, int node_permutation );
};

extern "C" int vtksnlRegisterElements();

#endif // Elements_Generic_h
