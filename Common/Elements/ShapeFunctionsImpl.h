/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$
  Language:  C++
  Date:      $Date: 2003-09-29 10:24:11 -0700 (Mon, 29 Sep 2003) $
  Version:   $Revision: 730 $

  Copyright 2012 Sandia Corporation.
  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
  license for use of this work by or on behalf of the
  U.S. Government. Redistribution and use in source and binary forms, with
  or without modification, are permitted provided that this Notice and any
  statement of authorship are reproduced on all copies.

=========================================================================*/
#ifndef ShapeFunctionsImpl_h
#define ShapeFunctionsImpl_h

// Legendre polynomials
static double legendre0(double* x)
{
	return 1.0000000000000000e+00;
}
static double legendre1(double* x)
{
	return x[0];
}
static double legendre2(double* x)
{
	return 1.5000000000000000e+00*(x[1])-5.0000000000000000e-01;
}
static double legendre3(double* x)
{
	return 2.5000000000000000e+00*(x[2])-1.5000000000000000e+00*x[0];
}
static double legendre4(double* x)
{
	return -3.7500000000000000e+00*(x[1])+4.3750000000000000e+00*(x[3])+3.7500000000000000e-01;
}
static double legendre5(double* x)
{
	return -8.7500000000000000e+00*(x[2])+7.8750000000000000e+00*(x[4])+1.8750000000000000e+00*x[0];
}
static double legendre6(double* x)
{
	return 6.5625000000000000e+00*(x[1])-1.9687500000000000e+01*((x[3]))+1.4437500000000000e+01*(x[5])-3.1250000000000000e-01;
}

// Legendre polynomial derivatives
static double dleg0(double *x)
{
	return 0.0;
}
static double dleg1(double *x)
{
	return 1.0000000000000000e+00;
}
static double dleg2(double *x)
{
	return 3.0000000000000000e+00*x[0];
}
static double dleg3(double *x)
{
	return 7.5000000000000000e+00*(x[1])-1.5000000000000000e+00;
}
static double dleg4(double *x)
{
	return 1.7500000000000000e+01*(x[2])-7.5000000000000000e+00*x[0];
}
static double dleg5(double *x)
{
	return -2.6250000000000000e+01*(x[1])+3.9375000000000000e+01*((x[3]))+1.8750000000000000e+00;
}
static double dleg6(double *x)
{
	return -7.8750000000000000e+01*(x[2])+8.6625000000000000e+01*(x[4])+1.3125000000000000e+01*x[0];
}

// Phi
static double phi2(double *x)
{
	return 6.1237243569579447e-01*(x[1]-1.);
}
static double phi3(double *x)
{
	return 7.9056941504209488e-01*(x[2]-x[0]);
}
static double phi4(double *x)
{
	return -1.4031215200402281e+00*(x[1])+1.1692679333668567e+00*(x[3])+2.3385358667337133e-01;
}
static double phi5(double *x)
{
	return -2.6516504294495533e+00*(x[2])+1.8561553006146871e+00*(x[4])+7.9549512883486595e-01*x[0];
}
static double phi6(double *x)
{
	return 2.1986323874172324e+00*(x[1])-5.1301422373068757e+00*(x[3])+3.0780853423841257e+00*(x[5])-1.4657549249448218e-01;
}

// Phi derivatives
static double dphi2(double *x)
{
	return 1.2247448713915889e+00*x[0];
}
static double dphi3(double *x)
{
	return 2.3717082451262845e+00*(x[1])-7.9056941504209488e-01;
}
static double dphi4(double *x)
{
	return 4.6770717334674270e+00*(x[2])-2.8062430400804561e+00*x[0];
}
static double dphi5(double *x)
{
	return -7.9549512883486599e+00*(x[1])+9.2807765030734366e+00*((x[3]))+7.9549512883486595e-01;
}
static double dphi6(double *x)
{
	return -2.0520568949227503e+01*(x[2])+1.8468512054304753e+01*(x[4])+4.3972647748344649e+00*x[0];
}

// Psi
static double psi2(double *x)
{
	return -2.4494897427831779e+00;
}
static double psi3(double *x)
{
	return -3.1622776601683795e+00*x[0];
}
static double psi4(double *x)
{
	return -4.6770717334674270e+00*(x[1])+9.3541434669348533e-01;
}
static double psi5(double *x)
{
	return -7.4246212024587486e+00*(x[2])+3.1819805153394638e+00*x[0];
}
static double psi6(double *x)
{
return 8.2082275796910018e+00*(x[1])-1.2312341369536503e+01*((x[3]))-5.8630196997792872e-01;
}

// Psi derivatives
static double dpsi2(double *x)
{
return 0.0;
}
static double dpsi3(double *x)
{
return -3.1622776601683795e+00;
}
static double dpsi4(double *x)
{
return -9.3541434669348540e+00*x[0];
}
static double dpsi5(double *x)
{
return -2.2273863607376246e+01*(x[1])+3.1819805153394638e+00;
}
static double dpsi6(double *x)
{
return -4.9249365478146011e+01*(x[2])+1.6416455159382004e+01*x[0];
}

#endif // ShapeFunctionsImpl_h
