/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-12-31 21:15:12 -0800 (Wed, 31 Dec 2003) $
  Version:   $Revision: 1201 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef VELHO_FUNCTIONS_H
#define VELHO_FUNCTIONS_H

#include <vtkObjectFactory.h>
#include <vtkCellType.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMath.h>
#include <vtkTetra.h>

#include <vtkShoeMesh.h>
#include <vtkShoeMeshIterator.h>
#include <vtkCellOps.h>
#include <vtkAdaptiveTessellator.h>

class vtkShoeElemVelhoFunctions
{
	public:
		// Description:
		// Evaluate the implicit function
		static double ifunction(double coord[3],vtkShoeMeshIterator& cell, double Isovalue,int field_no);

		// Description:
		// Find the intersection of an edge with the isosurface. Uses 
		static void rootfind(double p11[3], double p22[3], double p3[3], double Isovalue, vtkShoeMeshIterator& cell,int field_no);

		// Description:
		// Find intrsection of surface and a tetrahedron. Gives a triangle as an output.
		static void tri(int v0,int v1,int v2,int v3,double tetraPts[4][3],vtkShoeMeshIterator& cell,double Isovalue,int field_no,vtkAdaptiveTessellator* tess_in);

		// Description:
		// Find intersection of surface with a tetrahedron. Gives quadrilateral as an output.
		static void quad(int v0,int v1,int v2,int v3,double tetraPts[4][3],vtkShoeMeshIterator& cell,double Isovalue,int field_no,vtkAdaptiveTessellator* tess_in);

		// Description:
		// Determine whether the isosurface intersects a tetrahedron forming a triangle or quadrilateral.
		static void Velho(vtkTetra *tetra,vtkShoeMeshIterator& cell,double Isovalue,int field_no,vtkAdaptiveTessellator* tess_in);
};

#endif

