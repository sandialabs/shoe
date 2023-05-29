/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */

#include "Elements/ShoeInterpolants.h"
#include "Elements/ShoeHexahedron.h"
#include "vtkPolynomialSystem.h"

#include "vtkMath.h"

#include <string.h>

#undef VTK_DEBUG_LAGRANGE_UNI
#undef VTK_DEBUG_LAGRANGE_MULTI
#undef VTK_DEBUG_PSS_TEST

vtkstd::vector< vtkstd::vector<int> > ShoeInterpolants::BinomialCoefficients;
vtkstd::vector< vtkstd::vector<double> > ShoeInterpolants::LagrangeNormalizationFactors;
vtkstd::vector< vtkstd::vector< vtkstd::vector< double > > > ShoeInterpolants::LagrangeInterpolants;
vtkstd::vector< vtkstd::vector< vtkstd::vector< double > > > ShoeInterpolants::LagrangeDerivatives;

static void ComputeBinomialCoefficients( int order, vtkstd::vector<int>& N )
{
  if ( order < 0 )
    return;

  int fac = 1;
  N.push_back( fac );
 
  int sz = order / 2 + 1;
  for ( int o = 1; o < sz; ++ o )
    {
    fac *= order --;
    fac /= o;
    N.push_back( fac );
    }
}

static void ComputeLagrangeNormalization( int order, vtkstd::vector<double>& N )
{
  if ( order < 1 )
    return;

  int sz = order / 2 + 1;
  int sign = order % 2 ? -1 : 1;

  int NpN = order;
  for ( int i = 1; i < order; ++ i )
    {
    NpN *= order;
    }

  double fac1 = NpN / double( 1 << order );
  int fac2 = vtkMath::Factorial( order );
  for ( int o = 0; o < sz; ++ o )
    {
    N.push_back( sign * fac1 / fac2  );
    sign *= -1;
    fac2 *= o + 1;
    fac2 /= order - o;
    }
}

static void ComputeLagrangeCoefficients( int order, vtkstd::vector< vtkstd::vector< double > >& LC, vtkstd::vector< vtkstd::vector< double > >& LD )
{
  if ( order < 1 )
    return;

  int sz = order / 2 + 1;

  double xi[order + 1];
  xi[0] = -1;
  double dxi = 2. / order;
  for ( int i = 0; i < order; ++ i )
    {
    xi[i + 1] = xi[i] + dxi;
    }

  vtkstd::vector<double> LCi;
  vtkstd::vector<double> LDi;
  int orderm1 = order - 1;
  for ( int i = 0; i < orderm1; ++ i )
    {
    LCi.push_back( 0. );
    LDi.push_back( 0. );
    }
  LCi.push_back( 0. );
  LCi.push_back( 1. );
  LDi.push_back( order );

  for ( int o = 0; o < sz; ++ o )
    {
    LC.push_back(LCi);
    LD.push_back(LDi);

    for ( int d = 1; d < ( 1 << order ); ++ d )
      {
      int v = 0;
      int tmp = d;
      int j = 0;
      double fac = 1.;

      while ( tmp )
        {
        if ( tmp & 1 )
          {
          ++v;
          fac *= ( j < o ? xi[j] : xi[j + 1] );
          }
        tmp >>= 1;
        ++ j;
        }
      fac *= ( v % 2 ? -1 : 1 );
      LC[o][order - v] += fac;
      if ( order - v )
        {
        LD[o][order - v - 1] += ( order - v ) * fac;
        }
      }
    }
  
#ifdef VTK_DEBUG_LAGRANGE_UNI
  for ( int o = 0; o < sz; ++ o ) {
    vtkstd::cout << " L_" << o << "^" << order << "(X) = ";
    for ( int i = order; i > 0; -- i )
      vtkstd::cout << LC[o][i] << " X^" << i << " + ";
    vtkstd::cout << LC[o][0] << vtkstd::endl;

     vtkstd::cout << " L_" << o << "^" << order << "'(X) = ";
     for ( int i = order - 1; i > 0; -- i )
       vtkstd::cout << LD[o][i] << " X^" << i << " + ";
     vtkstd::cout << LD[o][0] << vtkstd::endl;
  }
#endif // VTK_DEBUG_LAGRANGE_UNI
}

void ShoeInterpolants::LagrangePrepareForOrder( const vtkShoeOrderTuple& o )
{
  for ( int i = 0; i < 3; ++ i )
    {
    if ( o.Order[i] < int(LagrangeNormalizationFactors.size()) )
      continue;
    for ( int j = LagrangeNormalizationFactors.size(); j <= o.Order[i]; ++ j )
      {
      vtkstd::vector<int> BC;
      BinomialCoefficients.push_back( BC );
      ComputeBinomialCoefficients( j, BinomialCoefficients[j] );
      
      vtkstd::vector<double> LN;
      LagrangeNormalizationFactors.push_back( LN );
      ComputeLagrangeNormalization( j, LagrangeNormalizationFactors[j] );

      vtkstd::vector< vtkstd::vector< double > > LC;
      LagrangeInterpolants.push_back( LC );
      vtkstd::vector< vtkstd::vector< double > > LD;
      LagrangeDerivatives.push_back( LD );
      ComputeLagrangeCoefficients( j, LagrangeInterpolants[j], LagrangeDerivatives[j] );
      }
    }
  
#ifdef VTK_DEBUG_PSS_TEST
  vtkPolynomialSystem* grad = vtkPolynomialSystem::New();
  int order[] = { 2, 2, 2 };
  ShoeHexahedron::SymbolicLagrangeFieldGradient( order, testH2, grad );
  grad->PrintSystem( vtkstd::cout );
  grad->SolveSystem();
  vtkstd::cout << grad->GetNumberOfRoots() << " total roots and ";
  grad->PruneComplexRoots( 1.0e-8 );
  int nr = grad->GetNumberOfRoots();
  vtkstd::cout << nr << " real roots.\n";
  double root[3];
  int i;
  for ( i=0; i<nr; ++i )
    {
    grad->GetRealRoot( i, root );
    vtkstd::cout << "    " << i << ": ( " << root[0] << ", " << root[1] << ", " << root[2] << " )\n";
    }
#endif // VTK_DEBUG_PSS_TEST
}

void ShoeHexahedron::LagrangeTensorShapeFunctions( const int* order, const double* pcoords, double* shape )
{
  // FIXME: Eventually needs to be varying length.
  double l[3][ShoeHexahedron::MaxDegree + 1];
  int i,j,k;

  double v;
  for ( i=0; i<3; ++i )
    {
    v = order[i] * ( pcoords[i] + 1. ) / 2;
    for ( j=0; j<=order[i]; ++j )
      {
      l[i][j] = 1.;
      for ( k=0; k<=order[i]; ++k )
        {
        if ( j != k )
          { // FIXME: (j - k) could be a register decremented inside the k loop:
            // and even better: normalization could be pre-computed and stored
            // somehow for each order that is actually used.
          l[i][j] *= ( v - k ) / ( j - k ); 
          }
        }
      }
    }

  int sn=0;

  // Corners
  shape[sn++] = l[0][0]        * l[1][0]        * l[2][0];
  shape[sn++] = l[0][order[0]] * l[1][0]        * l[2][0];
  shape[sn++] = l[0][order[0]] * l[1][order[1]] * l[2][0];
  shape[sn++] = l[0][0]        * l[1][order[1]] * l[2][0];
  shape[sn++] = l[0][0]        * l[1][0]        * l[2][order[2]];
  shape[sn++] = l[0][order[0]] * l[1][0]        * l[2][order[2]];
  shape[sn++] = l[0][order[0]] * l[1][order[1]] * l[2][order[2]];
  shape[sn++] = l[0][0]        * l[1][order[1]] * l[2][order[2]];

  int sn1, sn2, sn3;
  sn1 = order[0] + order[1] - 2;
  sn2 = sn1 * 2;
  sn3 = sn + sn1 + sn2;
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[0]; ++i )
    {
    //cout << sn << ", " << sn1 << ", " << sn2 << ", " << sn3 << "\n";
    shape[sn++ ] = l[0][i] * l[1][0] * l[2][0];               // Edge 0-1
    shape[sn1++] = l[0][i] * l[1][order[1]] * l[2][0];        // Edge 2-3
    shape[sn2++] = l[0][i] * l[1][0] * l[2][order[2]];        // Edge 4-5
    shape[sn3++] = l[0][i] * l[1][order[1]] * l[2][order[2]]; // Edge 6-7
    }
  for ( i=1; i<order[1]; ++i )
    {
    //cout << sn << ", " << sn1 << ", " << sn2 << ", " << sn3 << "\n";
    shape[sn++ ] = l[0][order[0]] * l[1][i] * l[2][0];        // Edge 1-2
    shape[sn1++] = l[0][0] * l[1][i] * l[2][0];               // Edge 3-0
    shape[sn2++] = l[0][order[0]] * l[1][i] * l[2][order[2]]; // Edge 5-6
    shape[sn3++] = l[0][0] * l[1][i] * l[2][order[2]];        // Edge 7-4
    }
  sn = sn3;
  sn1 = order[2] - 1;
  sn2 = sn1 * 2;
  sn3 = sn + sn1 + sn2;
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[2]; ++i )
    {
    //cout << sn << ", " << sn1 << ", " << sn2 << ", " << sn3 << "\n";
    shape[sn++ ] = l[0][0] * l[1][0] * l[2][i];                // Edge 0-4
    shape[sn1++] = l[0][order[0]] * l[1][0] * l[2][i];         // Edge 1-5
    // Kitware insists on swapping edges 10 and 11 as follows:
    shape[sn3++] = l[0][order[0]] * l[1][order[1]] * l[2][i];  // Edge 2-6
    shape[sn2++] = l[0][0] * l[1][order[1]] * l[2][i];         // Edge 3-7
    }

  sn = sn3;
  sn1 = (order[1] - 1)*(order[2] - 1);
  sn2 = sn1 * 2;
  sn3 = sn + sn2 + (order[2] - 1)*(order[0] - 1);
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[2]; ++i )
    {
    for ( j=1; j<order[1]; ++j )
      {
      //cout << sn << ", " << sn1 << "\n";
      shape[sn++ ] = l[0][0] * l[1][j] * l[2][i];        // Face 0-4-7-3
      shape[sn1++] = l[0][order[0]] * l[1][j] * l[2][i]; // Face 1-2-6-5
      }
    for ( j=1; j<order[0]; ++j )
      {
      //cout << sn2 << ", " << sn3 << "\n";
      shape[sn2++] = l[0][j] * l[1][0] * l[2][i];        // Face 0-1-5-4
      shape[sn3++] = l[0][j] * l[1][order[1]] * l[2][i]; // Face 2-3-7-6
      }
    }
  sn = sn3;
  sn1 = sn + (order[0] - 1)*(order[1] - 1);
  for ( i=1; i<order[1]; ++i )
    {
    for ( j=1; j<order[0]; ++j )
      {
      //cout << sn << ", " << sn1 << "\n";
      shape[sn++ ] = l[0][j] * l[1][i] * l[2][0];        // Face 0-1-2-3
      shape[sn1++] = l[0][j] * l[1][i] * l[2][order[2]]; // Face 4-7-6-5
      }
    }
  sn = sn1;
  for ( k=1; k<order[2]; ++k )
    {
    for ( j=1; j<order[1]; ++j )
      {
      for ( i=1; i<order[0]; ++i )
        {
        //cout << sn << "\n";
        shape[sn++] = l[0][i] * l[1][j] * l[2][k]; // Body
        }
      }
    }
}

#define VTK_LAGRANGE_DERIV(cmt,dder,aa,ab,ac,ba,bb,bc,ca,cb,cc) \
        VTK_DEBUG_LAGRANGE_PRINT_RST(cmt); \
        sgnm = sgnm ? ba : 0; \
        for ( int m = 0; m < aa; ++ m ) \
          { \
          sgnm *= -1; \
          sgnn = sgnn ? bb : 0; \
          for ( int n = 0; n < ab; ++ n ) \
            { \
            sgnn *= -1; \
            sgnp = sgnp ? bc : 0; \
            for ( int p = 0; p < ac; ++ p ) \
              { \
              sgnp *= -1; \
              VTK_DEBUG_LAGRANGE_PRINT_PWR(ca,cb,cc); \
              coef =  norm *  \
                      (sgnm ? sgnm : 1) * ca[order[0]][ip][m] * \
                      (sgnn ? sgnn : 1) * cb[order[1]][jp][n] * \
                      (sgnp ? sgnp : 1) * cc[order[2]][kp][p]; \
              dder[ p + op1[2] * ( n + op1[1] * m ) ] += coef; \
              } \
            } \
          }

#define VTK_LAGRANGE_DERIVS \
        VTK_LAGRANGE_DERIV("Dr:",rder,order[0],op1[1],op1[2],sgnM,-sgnN,-sgnP, \
          ShoeInterpolants::LagrangeDerivatives,ShoeInterpolants::LagrangeInterpolants,ShoeInterpolants::LagrangeInterpolants);      \
                                                                               \
        VTK_LAGRANGE_DERIV("Ds:",sder,op1[0],order[1],op1[2],-sgnM,sgnN,-sgnP, \
          ShoeInterpolants::LagrangeInterpolants,ShoeInterpolants::LagrangeDerivatives,ShoeInterpolants::LagrangeInterpolants);      \
                                                                               \
        VTK_LAGRANGE_DERIV("Dt:",tder,op1[0],op1[1],order[2],-sgnM,-sgnN,sgnP, \
          ShoeInterpolants::LagrangeInterpolants,ShoeInterpolants::LagrangeInterpolants,ShoeInterpolants::LagrangeDerivatives)

#ifdef VTK_DEBUG_LAGRANGE_MULTI
#  define VTK_DEBUG_LAGRANGE_PRINT_PWR(a,b,c) \
              vtkstd::cout << "  " << norm << " x "                           \
                           << "( " << ( sgnm ? sgnm : 1 ) << " ) "            \
                           << a [order[0]][ip][m]  << " x " \
                           << "( " << ( sgnn ? sgnn : 1 ) << " ) "            \
                           << b [order[1]][jp][n] << " x "  \
                           << "( " << ( sgnp ? sgnp : 1 ) << " ) "            \
                           << c [order[2]][kp][p]           \
                           << " X^" << m                                      \
                           << " Y^" << n                                      \
                           << " Z^" << p                                      \
                           << vtkstd::endl
#  define VTK_DEBUG_LAGRANGE_PRINT_LOOPVAR \
        vtkstd::cout << "At " << i << " " << j << " " << k                            \
                     << " ( proper: " << ip << " " << jp << " " << kp << " )"         \
                     << " ( signs: " << sgnM << " " << sgnN << " " << sgnP << " ):\n";
#  define VTK_DEBUG_LAGRANGE_PRINT_RST(x) vtkstd::cout << x "\n"
#else // VTK_DEBUG_LAGRANGE_MULTI
#  define VTK_DEBUG_LAGRANGE_PRINT_PWR(a,b,c)
#  define VTK_DEBUG_LAGRANGE_PRINT_LOOPVAR
#  define VTK_DEBUG_LAGRANGE_PRINT_RST(x)
#endif // VTK_DEBUG_LAGRANGE_MULTI

void ShoeHexahedron::LagrangeTensorShapeFunctionDerivatives( const int* order, const double* pcoords, double* sderiv )
{
  (void)order;
  (void)pcoords;
  (void)sderiv;
  std::cout << "UnivariateLagrange does not implement NumericShapeFunctionDerivatives.\n";
  std::cout << "TensorProduct does not implement NumericShapeFunctionDerivatives!\n";
}

void ShoeHexahedron::SymbolicLagrangeFieldGradient( const int* order, const double* phi_e, vtkPolynomialSystem* grad )
{
  // FIXME: Eventually needs to be varying length.
  // vtkPolynomialExpanderOperation* l[3][ShoeHexahedron::MaxDegree + 1];

  int i, j, k, sz[3], op1[3];
  for ( i = 0; i < 3; ++i )
    {
    sz[i] = order[i] / 2;
    op1[i] = order[i] + 1;
    }
  
  double prenorm,norm,coef;
  int llen = op1[0] * op1[1] * op1[2];
  double* rder = new double[ 3*llen ];
  double* sder = rder + llen;
  double* tder = sder + llen;
  memset( rder, 0, 3 * llen * sizeof(double) );

#ifdef VTK_DEBUG_LAGRANGE_MULTI
  vtkstd::cout << "#############################\n";
#endif // VTK_DEBUG_LAGRANGE_MULTI

  int ip,jp,kp,sgnm,sgnn,sgnp,sgnM,sgnN,sgnP;
  int sn = 0;

  // Corners
  ip = jp = kp = 0;

  sgnm = 0; sgnn = 0; sgnp = 0;
  sgnM = 1; 
  sgnN = 1; 
  sgnP = 1;
  norm = phi_e[sn++] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnm = 1;
  sgnM = order[0] % 2 ? -1 : 1; 
  norm = phi_e[sn++] * sgnM *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnn = 1;
  sgnN = order[1] % 2 ? -1 : 1;
  norm = phi_e[sn++] * sgnM * sgnN *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnm = 0;
  sgnM = 1; 
  norm = phi_e[sn++] * sgnN *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnn = 0; sgnp = 1;
  sgnN = 1; 
  sgnP = order[2] % 2 ? -1 : 1;
  norm = phi_e[sn++] * sgnP *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnm = 1;
  sgnM = order[0] % 2 ? -1 : 1; 
  norm = phi_e[sn++] * sgnM * sgnP *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnn = 1;
  sgnN = order[1] % 2 ? -1 : 1;
  norm = phi_e[sn++] * sgnM * sgnN * sgnP *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  sgnm = 0;
  sgnM = 1; 
  norm = phi_e[sn++] * sgnN * sgnP *
    ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
    ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
  VTK_LAGRANGE_DERIVS;

  int sn1, sn2, sn3;
  sn1 = order[0] + order[1] - 2;
  sn2 = sn1 * 2;
  sn3 = sn + sn1 + sn2;
  sn1 += sn;
  sn2 += sn;
  for ( i=1; i<order[0]; ++i )
    {
    if ( i  > sz[0] )
      {
      ip = order[0] - i;
      sgnM = order[0] % 2 ? -1 : 1;
      sgnm = 1;
      }
    else
      {
      ip = i;
      sgnM = 1;
      sgnm = 0;
      }

    // Edge 0-1
    sgnn = 0; sgnp = 0;
    sgnN = 1; 
    sgnP = 1;
    norm = phi_e[sn++] * sgnM *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;

    // Edge 2-3
    sgnn = 1;
    sgnN = order[1] % 2 ? -1 : 1; 
    norm = phi_e[sn1++] * sgnM * sgnN *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;

    // Edge 6-7
    sgnp = 1;
    sgnP = order[2] % 2 ? -1 : 1;
    norm = phi_e[sn3++] * sgnM * sgnN * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;

    // Edge 4-5
    sgnn = 0;
    sgnN = 1; 
    norm = phi_e[sn2++] * sgnM * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;
    } // i loop for r-direction edges
  ip = 0;

  for ( j=1; j<order[1]; ++j )
    {
    if ( j  > sz[1] )
      {
      jp = order[1] - j;
      sgnN = order[1] % 2 ? -1 : 1;
      sgnn = 1;
      }
    else
      {
      jp = j;
      sgnN = 1;
      sgnn = 0;
      }

    // Edge 1-2
    sgnm = 1; sgnp = 0;
    sgnM = order[0] % 2 ? -1 : 1; 
    sgnP = 1;
    norm = phi_e[sn++] * sgnM * sgnN *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;

    // Edge 3-0
    sgnm = 0;
    sgnM = 1;
    norm = phi_e[sn1++] * sgnN *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;

    // Edge 7-4
    sgnp = 1;
    sgnP = order[2] % 2 ? -1 : 1;
    norm = phi_e[sn3++] * sgnN * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;

    // Edge 5-6
    sgnm = 1;
    sgnM = order[0] % 2 ? -1 : 1; 
    norm = phi_e[sn2++] * sgnM * sgnN * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
    VTK_LAGRANGE_DERIVS;
    } // j loop for s-direction edges
  jp = 0;

  sn = sn3;
  sn1 = order[2] - 1;
  sn2 = sn1 * 2;
  sn3 = sn + sn1 + sn2;
  sn1 += sn;
  sn2 += sn;
  for ( k=1; k<order[2]; ++k )
    {
    if ( k  > sz[2] )
      {
      kp = order[2] - k;
      sgnP = order[2] % 2 ? -1 : 1;
      sgnp = 1;
      }
    else
      {
      kp = k;
      sgnP = 1;
      sgnp = 0;
      }

    // Edge 0-4
    sgnm = 0; sgnn = 0;
    sgnM = 1; 
    sgnN = 1;
    norm = phi_e[sn++] * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
    VTK_LAGRANGE_DERIVS;

    // Edge 1-5
    sgnm = 1;
    sgnM = order[0] % 2 ? -1 : 1; 
    norm = phi_e[sn1++] * sgnM * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
    VTK_LAGRANGE_DERIVS;

    // Kitware insists on swapping edges 10 and 11 as follows:
    // Edge 2-6
    sgnn = 1;
    sgnN = order[1] % 2 ? -1 : 1;
    norm = phi_e[sn3++] * sgnM * sgnN * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
    VTK_LAGRANGE_DERIVS;

    // Edge 3-7
    sgnm = 0;
    sgnM = 1; 
    norm = phi_e[sn2++] * sgnN * sgnP *
      ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
      ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
    VTK_LAGRANGE_DERIVS;
    } // k loop for t-direction edges
  kp = 0;

  sn = sn3;
  sn1 = (order[1] - 1) * (order[2] - 1);
  sn2 = sn1 * 2;
  sn3 = sn + sn2 + (order[2] - 1)*(order[0] - 1);
  sn1 += sn;
  sn2 += sn;
  for ( k=1; k<order[2]; ++k )
    {
    if ( k  > sz[2] )
      {
      kp = order[2] - k;
      sgnP = order[2] % 2 ? -1 : 1 ;
      sgnp = 1;
      }
    else
      {
      kp = k;
      sgnP = 1;
      sgnp = 0;
      }

    for ( j=1; j<order[1]; ++j )
      {
      if ( j  > sz[1] )
        {
        jp = order[1] - j;
        sgnN = order[1] % 2 ? -1 : 1;
        sgnn = 1;
        }
      else
        {
        jp = j;
        sgnN = 1;
        sgnn = 0;
        }

      // Face 0-4-7-3
      sgnM = 1;
      sgnm = 0;
      norm = phi_e[sn++] * sgnN * sgnP *
        ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
        ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
        ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
      VTK_LAGRANGE_DERIVS;
      
      // Face 1-2-6-5
      sgnM = order[0] % 2 ? -1 : 1;
      sgnm = 1;
      norm = phi_e[sn1++] * sgnM * sgnN * sgnP *
        ShoeInterpolants::LagrangeNormalizationFactors[order[0]][0] *
        ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
        ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
      VTK_LAGRANGE_DERIVS;
      }
    jp = 0;

    for ( i=1; i<order[0]; ++i )
      {
      if ( i  > sz[0] )
        {
        ip = order[0] - i;
        sgnM = order[0] % 2 ? -1 : 1;
        sgnm = 1;
        }
      else
        {
        ip = i;
        sgnM = 1;
        sgnm = 0;
        }

      // Face 0-1-5-4
      sgnN = 1;
      sgnn = 0;
      norm = phi_e[sn2++] * sgnM * sgnP *
          ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
      VTK_LAGRANGE_DERIVS;

      // Face 2-3-7-6
      sgnN = order[1] % 2 ? -1 : 1;
      sgnn = 1;
      norm = phi_e[sn3++] * sgnM * sgnN * sgnP *
          ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[1]][0] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp];
      VTK_LAGRANGE_DERIVS;
      }
    ip = 0;
    }
  kp = 0;

  sn = sn3;
  sn1 = sn + (order[0] - 1) * (order[1] - 1);
  for ( j=1; j<order[1]; ++j )
    {
    if ( j  > sz[1] )
      {
      jp = order[1] - j;
      sgnN = order[1] % 2 ? -1 : 1;
      sgnn = 1;
      }
    else
      {
      jp = j;
      sgnN = 1;
      sgnn = 0;
      }
    for ( i=1; i<order[0]; ++i )
      {
      if ( i  > sz[0] )
        {
        ip = order[0] - i;
        sgnM = order[0] % 2 ? -1 : 1;
        sgnm = 1;
        }
      else
        {
        ip = i;
        sgnM = 1;
        sgnm = 0;
        }

      // Face 0-1-2-3
      sgnP = 1;
      sgnp = 0;
      norm = phi_e[sn++] * sgnM * sgnN *
          ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
      VTK_LAGRANGE_DERIVS;

      // Face 4-7-6-5
      sgnP = order[2] % 2 ? -1 : 1;
      sgnp = 1;
      norm = phi_e[sn1++] * sgnM * sgnN * sgnP *
          ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[2]][0];
      VTK_LAGRANGE_DERIVS;
      }
    }
  ip = jp = 0;

  // Body modes
  sn = sn1;

  for ( k = 1; k < order[2]; ++k )
    {

    if ( k  > sz[2] )
      {
      kp = order[2] - k;
      sgnP = order[2] % 2 ? -1 : 1 ;
      sgnp = 1;
      }
    else
      {
      kp = k;
      sgnP = 1;
      sgnp = 0;
      }

    for ( j = 1; j < order[1]; ++j )
      {
      if ( j  > sz[1] )
        {
        jp = order[1] - j;
        sgnN = order[1] % 2 ? -1 : 1;
        sgnn = 1;
        }
      else
        {
        jp = j;
        sgnN = 1;
        sgnn = 0;
        }

      prenorm = sgnP * ShoeInterpolants::LagrangeNormalizationFactors[order[2]][kp] *
                sgnN * ShoeInterpolants::LagrangeNormalizationFactors[order[1]][jp];

      for ( i = 1; i < order[0]; ++i )
        {
        if ( i  > sz[0] )
          {
          ip = order[0] - i;
          sgnM = order[0] % 2 ? -1 : 1;
          sgnm = 1;
          }
        else
          {
          ip = i;
          sgnM = 1;
          sgnm = 0;
          }

        VTK_DEBUG_LAGRANGE_PRINT_LOOPVAR;

        norm = prenorm * sgnM * phi_e[sn++] *
          ShoeInterpolants::LagrangeNormalizationFactors[order[0]][ip];

        VTK_LAGRANGE_DERIVS;

        } // i loop
      } // j loop
    } // k loop


  // Now that we have coefficients, create the system of equations.
  int eqn[3];  
  eqn[0] = grad->InsertNextPolynomial();
  eqn[1] = grad->InsertNextPolynomial();
  eqn[2] = grad->InsertNextPolynomial();

  for ( i = 0; i < llen; ++i )
    {
      int m, n, p, q;
      p = i % op1[2];
      q = i / op1[2];
      n = q % op1[1];
      m = q / op1[1];

      if ( rder[i] )
        grad->InsertMonomial( eqn[0], rder[i], 3, 'r', m, 's', n, 't', p );
      if ( sder[i] )
        grad->InsertMonomial( eqn[1], sder[i], 3, 'r', m, 's', n, 't', p );
      if ( tder[i] )
        grad->InsertMonomial( eqn[2], tder[i], 3, 'r', m, 's', n, 't', p );
    }

  delete [] rder;
}


