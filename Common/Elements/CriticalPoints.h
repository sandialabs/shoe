/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-12-12 01:01:22 -0800 (Fri, 12 Dec 2003) $
  Version:   $Revision: 1145 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef Elements_CriticalPoints_h
#define Elements_CriticalPoints_h

#include <vtksnlConfigure.h>

#ifdef FOUND_GINAC
#include <ginac/ginac.h>

class vtkShoeElemCriticalPoints
{
	public:
		static GiNaC::ex eval_phi(const GiNaC::ex& j,const GiNaC::ex& xi);
		static GiNaC::ex eval_psi(int j, const GiNaC::ex& xi);
		static GiNaC::ex legendre_poly( const GiNaC::ex& n, const GiNaC::ex& x );
};
#endif // FOUND_GINAC
#endif // Elements_CriticalPoints_h
