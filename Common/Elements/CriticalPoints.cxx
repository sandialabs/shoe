/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include <vtksnlConfigure.h>

#include <stdio.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <iostream>
#include <vector>

#include <vtkMath.h>
#include <vtkDoubleArray.h>

#ifdef FOUND_GINAC
#include <ginac/ginac.h>
#endif

#include <gsl/gsl_poly.h>

#include <vtkCellEnums.h>
#include <vtkCellOps.h>
#include <vtkShoeMeshIterator.h>
#include <Elements/CriticalPoints.h>
#include "pss.h"

using namespace shoe;

#ifdef FOUND_GINAC
USING_NAMESPACE(GiNaC);
#endif // FOUND_GINAC

// pss function for solving a system of polynomial equations.
extern "C" int solveSystem(pss_polytope *p, pss_complex *coefs, pss_complex **roots);
static int iPrint=0;

// Function which calls pss polynomial solver routine.
void pssSolver(vtkstd::vector<vtkstd::vector<double> >& sysroots, int neqn,int nvar,vtkstd::vector<vtkstd::vector<int> >& exponents, vtkstd::vector<double>& coefficients, vtkstd::vector<int> &number_of_terms)
{
  vtkstd::cout<<"Number of coefficients "<<coefficients.size()<<" "<<exponents.size()<<vtkstd::endl;
  int cnt=1;
  if ( iPrint == 1 )
  {
    vtkstd::cout << "pssSolver; Coeff Input" << vtkstd::endl;
    for ( int cnt = 0; cnt < (int) coefficients.size(); ++cnt )
      vtkstd::cout << coefficients[cnt] << vtkstd::endl;

    vtkstd::cout << "pssSolver: number of terms " << vtkstd::endl;
    for ( int cnt = 0; cnt < (int) number_of_terms.size(); ++cnt )
      vtkstd::cout << number_of_terms[cnt] << vtkstd::endl;
  }


  // the polytope structure.
  pss_polytope p;
  p.nvar= nvar;   // Number of variables
  p.neq = neqn;    // Number of equations
  p.firstrow = (int *)malloc(sizeof(int)*(p.neq+1));  
  // each monomial term is a row. If all the rows of all the polynomials
  // are arranged in a list then the array firstrow gives the offset into that list for the starting row
  // of each polynomial. One extra entry is there which is equal to total number of rows.
  
  //array giving number of rows in each polynomial
  p.numrows = (int *)malloc(sizeof(int)*(p.neq));
  cnt=0;
  for(vtkstd::vector<int>::iterator i = number_of_terms.begin();i!=number_of_terms.end();i++)
  {
    p.numrows[cnt++] = *i;
  }
  
  // Type  = 1.
  p.type = (int *)malloc(sizeof(int)*(p.neq));
  for(cnt=0;cnt<p.neq;cnt++)
    p.type[cnt]=1;
  
  // setting the firstrow offsets
  p.firstrow[0]=0;  // first entry is allways zero
  for(int cnt=1;cnt<=p.neq;cnt++)
  {
      p.firstrow[cnt] = p.firstrow[cnt-1]+p.numrows[cnt-1];
      // add the number of rows for the previous polynomial to the previous entry
  }
  
  // the degree exponents. Size =  total number of monomials * number of variables
  p.A = (double *)malloc(sizeof(double)*p.firstrow[p.neq]*p.nvar);
  int ind = 0;

  for( vtkstd::vector<vtkstd::vector<int> >::iterator i = exponents.begin(); i!=exponents.end(); i++ )
  {
    vtkstd::vector<int>::iterator iint = (*i).begin();
    for(;iint!=(*i).end();iint++)
    {
      p.A[ind++] = *iint;
    }
  }
  pss_complex coefs[p.firstrow[p.neq]+p.nvar];
        cnt=0;
  if ( (int) coefficients.size() != p.firstrow[p.neq] )
  {
    vtkstd::cout << "Number of coefficients not equal to number of terms" << vtkstd::endl;
    exit(1);
  }

  
  for ( vtkstd::vector<double>::iterator i=coefficients.begin(); i!=coefficients.end(); i++ )
    coefs[cnt++].re=*i;

  for ( int cnt=p.firstrow[p.neq]; cnt<p.firstrow[p.neq]+p.nvar; cnt++ )
    coefs[cnt].re=1.;

  for ( cnt=0; cnt<p.firstrow[p.neq]+p.nvar; cnt++ )
    coefs[cnt].im=0.;
  pss_complex *roots;
  pss_complex **roots1=&roots;
  int numroots;
  
  // call the pss function.
  numroots=solveSystem(&p, coefs, roots1);
  vtkstd::cout<<"numroots "<<numroots<<vtkstd::endl;

  //cout<<"Number of roots "<<numroots<<endl<<endl<<endl;
  int imaginaryCount=0;
  for ( int cnt=0; cnt<numroots*p.nvar; cnt++ )
  {
    vtkstd::vector<double> tempVec;
    imaginaryCount=0;
    for ( int i=0; i<p.nvar; i++,cnt++ )
    {
      //cout<<(*roots1)[cnt].re<<" i"<<(*roots1)[cnt].im<<" ";
      if ( ((*roots1)[cnt].im<1.e-08) && ((*roots1)[cnt].im>-1.e-08) )
      {
        tempVec.push_back((*roots1)[cnt].re);
        imaginaryCount++;
      }
    }
    if ( imaginaryCount == p.nvar )
    {
      //cout<<tempVec[0]<<" "<<tempVec[1]<<" "<<tempVec[2]<<endl;
      sysroots.push_back(tempVec);
    }

    cnt--;
    // cout<<endl;
  }
  free(p.A);
  free(p.type);
  free(p.numrows);
  free(p.firstrow);
}

#ifdef FOUND_GINAC
ex vtkShoeElemCriticalPoints::legendre_poly( const ex& n, const ex& x )
{
  if ( n < 0 )
    throw "Legendre Polynomials only defined for non-negative n";

  // Use Bonnet's recursion formula
  if ( n == 0 )
    return 1;
  else if ( n == 1 )
    return x;

  return ( (2*n-1)*x*legendre_poly(n-1,x) - (n-1)*legendre_poly(n-2,x) ).expand().collect(x) / n ;
}

ex vtkShoeElemCriticalPoints::eval_phi(const ex& j,const ex& xi)
{
  if ( j < 2 )
    throw "Phi function indices start at 2";

  return 1/GiNaC::sqrt(2*(2*j-1))*(legendre_poly(j,xi) - legendre_poly(j-2,xi)).expand().collect(xi);
}

ex vtkShoeElemCriticalPoints::eval_psi(int j, const ex& xi)
{
  if ( j < 2 )
    throw "Psi function indices start at 2";
  switch ( j )
  {
    case 2:
      return -2.4494897427831779e+00;
      break;
    case 3:
      return -3.1622776601683795e+00*xi;
      break;
    case 4:
      return -4.6770717334674270e+00*(xi*xi)+9.3541434669348533e-01;
      break;
    case 5:
      return -7.4246212024587486e+00*(xi*xi*xi)+3.1819805153394638e+00*xi;
      break;
    case 6:
      return 8.2082275796910018e+00*(xi*xi)-1.2312341369536503e+01*((xi*xi*xi*xi))-5.8630196997792872e-01;
      break;
  }
}
#endif // FOUND_GINAC
  

// Find critical points for a field. A field may be a vector field or a scalar field. This function calls 
// FindPolynomialCriticalPoints for every component of the field.

void vtkCellOps::FindFieldCriticalPoints(
  vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > >& fieldRoots,
  vtkstd::vector<vtkstd::vector<int> >& fieldTypes,
  int maxExponents[3],
  const vtkstd::vector<vtkstd::vector<double> >& coef,
  const vtkstd::vector<vtkstd::vector<vtkstd::vector<int> > >& exp,
  int dimension) const
{

  int components = coef.size();

  for (int Ccnt=0; Ccnt<components; Ccnt++ )   // Ccnt=ComponentCount.
  {
    vtkstd::vector<vtkstd::vector<double> > sysroots;
    vtkstd::vector<int> types;
    // Temporary arrays storing extrema for each component
    this->FindPolynomialCriticalPoints(sysroots,types,maxExponents,coef[Ccnt],exp[Ccnt],dimension);
    // Add the type and field value at the end of each root entry.
    for ( int i = 0; i < (int) sysroots.size(); ++i )
    {
      // For each field component
      for ( int j = 0; j < components; ++j )
      {
        double value=0.;
        // for each monomial
        for ( int cnt1 = 0; cnt1 < (int) coef[j].size(); ++cnt1 )
        {
          double term=1.;
          // Cacluate the term r^a*s^b..
          for ( int cnt2 = 0; cnt2 < dimension; ++cnt2 )
          {
            term = term*::pow(sysroots[i][cnt2],exp[j][cnt1][cnt2]);
          }
          // add to the value after multiplying with coefficients     
          value+=term*coef[j][cnt1];
        }
        //cout<<"value "<<value<<endl;
        sysroots[i].push_back(value);
        
        //cout<<sysroots[i][0]<<" "<<sysroots[i][1]<<" "<<sysroots[i][2]<<" "<<sysroots[i][3]<<endl;
      }
      types[i]=types[i]+(8*Ccnt); //   Store the component in the type
    }

    // Sort the values in ascending order.
    if ( sysroots.size() > 0 )
    {
      for ( int i = 0; i < (int) sysroots.size() - 1; ++i )
      {
        for ( int j = i + 1; j < (int) sysroots.size(); ++j )
        {
          if ( sysroots[i][dimension + Ccnt] > sysroots[j][dimension + Ccnt] )
          {
            double temp;
            for ( int k = 0; k < (int) sysroots[j].size(); ++k )
            {
              temp = sysroots[i][k];
              sysroots[i][k] = sysroots[j][k];
              sysroots[j][k] = temp;
            }
            temp = types[i];
            types[i] = types[j];
            types[j] = (int) temp;
          }
        }
      }
    }

#if 0
    for(int i=0;i<sysroots.size();i++)
      {
      vtkstd::cout
        << sysroots[i][0] << " " << sysroots[i][1] << " "
        << sysroots[i][2] << " " << sysroots[i][3] << " "
        << types[i] << vtkstd::endl;
      }
    vtkstd::cout << vtkstd::endl;
#endif // 0

    fieldRoots.push_back( sysroots );
    fieldTypes.push_back( types );
    sysroots.clear();
    types.clear();
  }
}



// Find critical points of a polynomial
// Input : 
// maxExponents = maximum exponent in each direction. /// ????????
// coefficients = the array of coef
// exp = the tuples of exponents
// dimension = 1D/2D/3D
// Output :
// vals = tuples storing (r,s,t,val) giving the critical points of the polynomials
// types = type of each critical point min/max/saddle  0 -> Min , 1->Saddle, 2->Max, 3->Singular Hessian, 4->Corner Value.
// 
void vtkCellOps::FindPolynomialCriticalPoints(vtkstd::vector<vtkstd::vector<double> >& sysroots, vtkstd::vector<int>& types, int maxExponents[3], const vtkstd::vector<double>& coef, const vtkstd::vector<vtkstd::vector<int> >& exp, int dimension) const
{
  if(dimension<=1)
  {
    vtkstd::cout<<"FindPolynomialCriticalPoints: Dimension should be greater than 1 "<<vtkstd::endl;
  }
  if(coef.size()!=exp.size())
  {
    vtkstd::cout<<"FindPolynomialCriticalPoints: The size of expoenent and coeff arrays are not equal"<<vtkstd::endl;
    exit(1);
  }
  
  if(coef.size()==0)
    return;
  
  //construct the system of derrivatives of the functions
  
  //Array of exponents of the system
  vtkstd::vector<vtkstd::vector<int> >exponents;
  
  //Coefficiants of the systems
  
  vtkstd::vector<double> coefficients;

  // Number of terms in each system
  
  vtkstd::vector<int> number_of_terms;

  // Critical point coordinates from the solver
  
  // Create the system coefficients and exponents.
  
  
  int terms_counter=0;

  // For each direction
  for ( int cnt = 0; cnt < dimension; ++cnt )
  {
    terms_counter = 0;    
    // For each coefficient
    for ( int cnt2 = 0; cnt2 < (int) coef.size(); ++cnt2 )
    {
      
      if ( exp[cnt2][cnt] == 0 )
        continue;
      
      vtkstd::vector<int> temp;
      // For each variable
      for ( int cnt3 = 0; cnt3 < dimension; ++cnt3 )
      {
        if ( cnt3 != cnt )
          temp.push_back( exp[cnt2][cnt3] );
        else
          temp.push_back( exp[cnt2][cnt] - 1 );
      }
      exponents.push_back( temp );
      coefficients.push_back( coef[cnt2] * exp[cnt2][cnt] );
      terms_counter++;
    }
    if ( terms_counter <= 1 )
    {
      vtkstd::cout << "There  is a degenerate equation ***********" << vtkstd::endl;
      return;
    }
 
    number_of_terms.push_back( terms_counter );
  }
      
  // pass the system of derrivatives to psssolver
  
  pssSolver(sysroots,dimension,dimension,exponents,coefficients,number_of_terms);
  vtkstd::cout<<"number of roots "<<sysroots.size()<<vtkstd::endl;
  
  for ( int rootCnt = 0; rootCnt < (int) sysroots.size(); ++rootCnt )
  {
    vtkstd::cout
      << "size of the tuple " << sysroots[rootCnt].size() << vtkstd::endl
      << sysroots[rootCnt][0] << " " <<sysroots[rootCnt][1] << " " << sysroots[rootCnt][2] << vtkstd::endl;

    for ( int cnt = 0; cnt < dimension; ++cnt )
    {
      double value=0.;
      // For each coefficient
      for ( int cnt2 = 0; cnt2 < (int) coef.size(); ++cnt2 )
      {
        if ( exp[cnt2][cnt] == 0 )
          continue;
        
        double term=1.;

        // For each variable
        for ( int cnt3 = 0; cnt3 < dimension; ++cnt3 )
        {
          if ( cnt3 != cnt )
            term *= ::pow( sysroots[rootCnt][cnt3], exp[cnt2][cnt3] );
          else
            term *= ::pow( sysroots[rootCnt][cnt3], exp[cnt2][cnt] - 1 );
        }
        value += term * coef[cnt2] * exp[cnt2][cnt];
      }
    }

  }

  // calculate the value and add it to the coordinates
  double value;

  // For each root
  for ( int cnt = 0; cnt < (int) sysroots.size(); ++cnt )
  {
    value = 0.;
    types.push_back( 1 );
  }
  //cout<<"calculated values "<<endl;
    

  // Calculate the hessian and determine types
  double hess[3][3];
  // For each root
  for ( int cnt4 = 0; cnt4 < (int) sysroots.size(); ++cnt4 )
    {
    // For each variable take the partial derrivative
    for ( int cnt = 0; cnt < dimension; ++cnt )
      {
      // For each variable take partial derrivatuive
      for ( int cnt1 = 0; cnt1 < dimension; ++cnt1 )
        {
        double value = 0.;
        //For each monomial term
        for ( int cnt2 = 0; cnt2 < (int) coef.size(); ++cnt2 )
          {
          // If differntiated twice by same variable, if the exponent of that variable in the
          // monomial term is less than 2, neglect that term, it will be zero
          if ( cnt == cnt1 )
            {
            if ( exp[cnt2][cnt] < 2 )
              continue;
            }  // If differntiated by different variables, if exponent of any of them is less than 1, neglect the term
          else 
            {
            if ( exp[cnt2][cnt] == 0 || exp[cnt2][cnt1] == 0 )
              continue;
            }

          double term=1.;
          // calculate the term product

          for ( int cnt3 = 0; cnt3 < dimension; ++cnt3 )
            {
            // If the term is not diferntiated, raise by full power
            if ( cnt3 != cnt && cnt3 != cnt1 )
              {
              term *= ::pow( sysroots[cnt4][cnt3], exp[cnt2][cnt3] );
              }
            else 
              {
              // If the term is diffentiated only once by this variable, reduce power by one
              if ( cnt != cnt1 )
                term *= ::pow( sysroots[cnt4][cnt3], exp[cnt2][cnt3] - 1 );
              else // if twice differentiated, reduce power by 2
                term *= ::pow( sysroots[cnt4][cnt3], exp[cnt2][cnt3] - 2 );
              }
            }
          // If the term is twice differentiated in the same variable, multiply by the exp*(exp-1) 
          if ( cnt == cnt1 )
            value += coef[cnt2] * exp[cnt2][cnt] * (exp[cnt2][cnt] - 1) * term;
          else
            value += coef[cnt2] * exp[cnt2][cnt] * exp[cnt2][cnt1] * term;
          // If the term is differentiated by different variables, multiply by the exponents of those two variables.
          }
        //cout<<value<<endl;
        hess[cnt][cnt1] = value;
        }
      }
    //  Calculate the eigenvalues of the hessian matrix if all the values are +ve definite  then it is minima
    //  if -ve definite then it is maxima. If the hessian is 0 then it is degenerate. Otherwise it is saddle.
    double determinant;
    if ( iPrint == 1 )
      {
      vtkstd::cout << "hessian matrix " << vtkstd::endl;
      for ( int temp = 0; temp < dimension; ++temp )
        {
        for ( int temp1 = 0; temp1 < dimension; ++temp1 )
          {
          vtkstd::cout << hess[temp][temp1] << " ";
          }
        vtkstd::cout << vtkstd::endl;
        }
      vtkstd::cout << vtkstd::endl;
      }

    if ( dimension == 2 )
      {
      determinant=hess[0][0]*hess[1][1] - hess[0][1]*hess[1][0];
      }
    else if ( dimension == 3 )
      {
      determinant = vtkMath::Determinant3x3( hess );
      }
    else 
      {
      vtkstd::cout << "FindPolynomialCriticalPoints : dimension upto 3  are supported " << vtkstd::endl;
      exit( 1 );
      }

    if ( determinant < 1.e-6 && determinant > -1.e-6 )
      {
      vtkstd::cout << "singular hessian " << vtkstd::endl;
      //sysroots[cnt4].push_back(3);
      types[cnt4] = 3;
      continue;
      }
    double eigenvalues[3]; 
    double** eigenvectors = (double**) malloc( dimension * sizeof(double*) );
    double** hessian = (double**) malloc( dimension * sizeof(double*) );
    for ( int temp = 0; temp < dimension; ++temp )
      {
      eigenvectors[temp] = (double*) malloc( dimension * sizeof(double) );
      hessian[temp] = hess[temp]; // This is required to call jacobiN. It takes double** as argument and complains for hess which is
      }
    
    // Find the eigenvalues and eigenvectors.
    vtkMath::JacobiN( hessian, dimension, eigenvalues, (double**) eigenvectors );
    int flag=0;
    double oldvalue = eigenvalues[0];
    for ( int temp = 0; temp < dimension; ++temp )
      {
      // Check for zero eigenvalue. If any value is zero, set the type to saddle point.
      if ( eigenvalues[temp] < .000001 && eigenvalues[temp] > -.000001 )
        {
        //cout<<"saddle point "<<endl;
        //sysroots[cnt4].push_back(1);
        types[cnt4] = 1;
        flag = 1;
        break;
        }
      // check for sign change. If sign changes the point is saddle point.
      if ( eigenvalues[temp] * oldvalue < 0 )
        {
        //cout<<"saddle point "<<endl;
        types[cnt4] = 1;
        //sysroots[cnt4].push_back(1);
        flag = 1;
        break;
        }
      oldvalue = eigenvalues[temp];
      }

    if( flag == 1 )
      {
      continue;
      }

    // There is no sign change as well as no zero value. Thus the point is maxima or minima.
    if ( eigenvalues[0] > 0 )
      {
      types[cnt4] = 0;
      }
    else if ( eigenvalues[0] < 0 )
      {
      types[cnt4] = 2;
      }
    //sysroots[cnt4].push_back(types[cnt4]);

    for ( int i = 0; i < dimension; ++i )
      {
      free( eigenvectors[i] );
      }
    free( eigenvectors );
    free( hessian );
    }
}

// Reject the critical values which fall out of the bounds.
// input : extrema - vector consisting of all critical points and values for each component of a field.
// types - type of the extrema points for each component of a field.
// dimension = dimension of the proble.
// param_bounds - array consisting of parameter bounds.
// Domainshape - hex/tet or prism 
// CornerValues - field values at the corner nodes of the domain.
// Output :
// extrema = extrema vales lying inside the bounds.
// types - the types of extrema  0 -> Minimum, 1 -> Saddle, 2->Maximum 4->Corner Value 
void vtkCellOps::FindBoundedCriticalValues(
  vtkstd::vector<vtkstd::vector<vtkstd::vector<double> > >& extrema,
  vtkstd::vector<vtkstd::vector<int> >& types,
  int dimension,
  const double* param_bounds,
  shoe::CellShape const DomainShape,
  const vtkstd::vector<vtkstd::vector<double> >& CornerValues) const
{
  if ( dimension > 3 )
    {
    vtkstd::cout<<"dimensions greater than 3 are not supported by FindBoundedCriticalValues"<<vtkstd::endl;
    return;
    }
  int components = extrema.size();
  // Remove entries which are out of bounds.

  // For every component of field
  for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
    {
    int cnt1 = 0, flag = 0;
    vtkstd::vector<vtkstd::vector<double> >::iterator i = extrema[Ccnt].begin();
    vtkstd::vector<int>::iterator itype = types[Ccnt].begin();
    // For every extrema 
    for ( ; i != extrema[Ccnt].end(); )
      {
      flag=0;
      if ( DomainShape==Hexahedron || DomainShape==Quadrilateral || DomainShape == Curve )
        {
        // for each variable
        for ( int cnt = 0; cnt < dimension; ++cnt )
          {
          if ( (*i)[cnt] < param_bounds[0] || (*i)[cnt] > param_bounds[1] )     
            {
            // Remove that extrema
            extrema[Ccnt].erase( i ); // Erase the extrema
            types[Ccnt].erase( itype ); // Erase the type entry
            i = extrema[Ccnt].begin() + cnt1; // Reset the iterator to the next element
            itype = types[Ccnt].begin() + cnt1;
            flag = 1;
            break;
            }
          }
        }
      else if ( DomainShape == Tetrahedron )  
        {
        double l1 = 0.5 * ( 1 - (*i)[0] - (1 / sqrt( 3. )) * (*i)[1] - (1 / sqrt( 6. )) * (*i)[2] );
        double l2 = 0.5 * ( 1 + (*i)[0] - (1 / sqrt( 3. )) * (*i)[1] - (1 / sqrt( 6. )) * (*i)[2] );
        double l3 = ( sqrt( 3. ) / 3.) * ( (*i)[1] - ( 1 / sqrt( 8. )) * (*i)[2] );
        vtkstd::cout << "the bary coordinates " << l1<< " " << l2 << " " << l3 << vtkstd::endl;
        vtkstd::cout << (*i)[0] << " " << (*i)[1] << " "  << (*i)[2] << vtkstd::endl;
        if ( l1 < 0 || l2 < 0 || l3 < 0 || (*i)[2] < 0 )
          {
          extrema[Ccnt].erase( i );
          types[Ccnt].erase( itype );
          i = extrema[Ccnt].begin() + cnt1;
          itype = types[Ccnt].begin() + cnt1;
          flag = 1;
          }
        }
      else if ( DomainShape == Triangle )
        {
        double l1 = .5 * ( 1 - (*i)[0] - (1 / sqrt( 3. )) * (*i)[1] );
        double l2 = .5 * ( 1 + (*i)[0] - (1 / sqrt( 3. )) * (*i)[1] );
        vtkstd::cout << "Triangle " << vtkstd::endl;
        if ( l1 < 0 || l2 < 0 || (*i)[1] < 0 )
          {
          extrema[Ccnt].erase( i );
          types[Ccnt].erase( itype );
          i = extrema[Ccnt].begin() + cnt1;
          itype = types[Ccnt].begin() + cnt1;
          flag = 1;
          }
        }
      else
        {
        vtkstd::cout << "This element shape is not yet supported by FindBoundedCriticalValues" << vtkstd::endl;
        exit(1);
        }

      if ( flag == 0 ) // if the extrema was not removed, increament the counter and the iterator
        {
        ++cnt1;
        ++i;
        ++itype;
        }
      }
    }

  for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
    {
    if ( extrema[Ccnt].size() != types[Ccnt].size() )
      {
      vtkstd::cout << "FindBoundedCriticalValues: Number of entries in extrema array and types array are not equal" << vtkstd::endl;
      exit(1);
      }
    }

  // Find min and max
  // Compare with corner values..if corner values are extremum or minimum, add them 
  int cornerMin=0,cornerMax=0;
  // Find the minima and maxima of the corner values
  for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
    {
    for ( int cnt = 0; cnt < (int) CornerValues.size(); ++cnt )
      {
      if ( CornerValues[cnt][Ccnt] < CornerValues[cornerMin][Ccnt] )
        {
        cornerMin=cnt;
        }
      else if ( CornerValues[cnt][Ccnt] > CornerValues[cornerMax][Ccnt] )
        {
        cornerMax=cnt;
        }
      }
    vtkstd::vector<double> tempCorner;
    vtkstd::vector<double> tempCorner1;
    if ( extrema[Ccnt].size() == 0 )
      {
      for ( int cnt = 0; cnt < dimension; ++cnt )
        tempCorner.push_back( cornerMin );

      for ( int i = 0; i < components; ++i )
        tempCorner.push_back( CornerValues[cornerMin][i] );

      extrema[Ccnt].push_back( tempCorner );
      types[Ccnt].push_back( 8 * Ccnt + 4 );
      for ( int cnt = 0; cnt < dimension; ++cnt )
        tempCorner1.push_back( cornerMax );
      
      for ( int i = 0; i < components; ++i )
        tempCorner1.push_back( CornerValues[cornerMax][i] );
      
      extrema[Ccnt].push_back( tempCorner1 );
      types[Ccnt].push_back( 8 * Ccnt + 5 );
      continue;
      }
    else
      {
      if ( CornerValues[cornerMin][Ccnt] <= extrema[Ccnt][0][dimension + Ccnt] )
        {
        for ( int cnt = 0; cnt < dimension; ++cnt )
          tempCorner.push_back( cornerMin );

        for ( int i=0;i<components;i++)
          tempCorner.push_back( CornerValues[cornerMin][i] );

        extrema[Ccnt].insert( extrema[Ccnt].begin(), tempCorner );
        types[Ccnt].insert( types[Ccnt].begin(), 8 * Ccnt + 4 );
        }

      if ( CornerValues[cornerMax][Ccnt] >= extrema[Ccnt][extrema[Ccnt].size() - 1][dimension + Ccnt] )
        {
        for ( int cnt = 0; cnt < dimension; ++cnt )
          tempCorner1.push_back( cornerMax );

        for( int i = 0; i < components; ++i )
          tempCorner1.push_back( CornerValues[cornerMax][i] );

        extrema[Ccnt].push_back( tempCorner1 );
        types[Ccnt].push_back( 8 * Ccnt + 5 );
        }
      }
    tempCorner.clear();
    tempCorner1.clear();
    }

  // Move the maxim to the 2nd place.
  for ( int Ccnt = 0; Ccnt < components; ++Ccnt )
    {
    if ( extrema[Ccnt].size() > 2 )
      {
      double temp;
      for ( int i = 0; i < (int) extrema[Ccnt][0].size(); ++i )
        {
        temp = extrema[Ccnt][1][i];
        extrema[Ccnt][1][i] = extrema[Ccnt][extrema[Ccnt].size() - 1][i];
        extrema[Ccnt][extrema[Ccnt].size() - 1][i] = temp;
        }
      temp = types[Ccnt][1];
      types[Ccnt][1] = types[Ccnt][types[Ccnt].size() - 1];
      types[Ccnt][types[Ccnt].size() - 1] = int(temp);
      }

    for ( int cnt = 0; cnt < (int) extrema[Ccnt].size(); ++cnt )
      {
      extrema[Ccnt][cnt].push_back( types[Ccnt][cnt] ); // Add the type at the end of extrema tuple.
      }
    }
}

// This function finds the critical points of a univariate polynomial using the Gnu Scientific Library (gsl).
void vtkCellOps::FindUnivariateCriticalPoints(
  vtkstd::vector<vtkstd::vector<double> >& sysroots,
  vtkstd::vector<int>& types,
  vtkstd::vector<double>& coef,
  vtkstd::vector< vtkstd::vector<int> >& powers) const
{
  int temp=0;
  // Sort the powers..lowest to highest power..
  int n = powers.size();
  if ( n == 0 )
    {
    return;
    }

  for ( int cnt = 0; cnt < n - 1; ++cnt )
    {
    for ( int cnt1 = cnt + 1; cnt1 < (int) powers.size(); ++cnt1 )
      {
      if ( powers[cnt1][0] < powers[cnt][0] )
        { 
        int temp1;
        double temp2;
        temp1 = powers[cnt][0];
        powers[cnt][0] = powers[cnt1][0];
        powers[cnt1][0] = temp1;
        temp2 = coef[cnt];
        coef[cnt] = coef[cnt1];
        coef[cnt1] = temp2;
        }
      }
    }

  if ( powers[powers.size() - 1][0] < 2 )
    {
    // The polynomial is linear or a constant. It does not have any critical points.
    return;
    }

  vtkstd::vector<double> coef_deriv;
  vtkstd::vector<int> powers_deriv;
  // Build the derivative of the polynomial
  for ( int cnt = 0; cnt < (int) powers.size(); ++cnt )
    {
    if ( powers[cnt][0] > 0 )
      {
      coef_deriv.push_back( powers[cnt][0] * coef[cnt] );
      powers_deriv.push_back( powers[cnt][0] - 1 );
      }
    }
  // Generate the array to pass to the solver.
  // The array needs to be n+1 long where n is the degree of the univariate polynomial.
  // The terms which are absent have zero coefficient.
  temp = 0;
  double* c = new double[2 * powers_deriv[powers_deriv.size() - 1]];
  for ( int i = 0; i <= powers_deriv[powers_deriv.size() - 1]; ++i )
    {
    // Put zero for the coefficient term if the term with this power does not exist.
    if ( i != powers_deriv[temp] )
      {
      c[i] = 0;
      continue;
      }
    c[i] = coef_deriv[temp];
    vtkstd::cout << "coef " << c[i] << vtkstd::endl;
    ++temp;
    }
  
  double* z = new double[2 * powers_deriv[powers_deriv.size() - 1]];
  gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc( powers_deriv[powers_deriv.size() - 1] + 1 );
  gsl_poly_complex_solve( c, powers_deriv[powers_deriv.size() - 1] + 1, w, z );
  gsl_poly_complex_workspace_free( w );

  // Store the roots in output array
  for ( int i = 0; i < powers_deriv[powers_deriv.size() - 1]; ++i )
    {
    if ( z[2 * i + 1] > -0.000001 && z[2 *i + 1] < 0.000001) // If imaginary part is zero
      {
      vtkstd::vector<double> tempvec;
      tempvec.push_back( z[2 * i] );
#if 0
      double value = 0.;
#endif // 0
      double dderiv = 0.;

      // Calculate the function value and the double derivative value.
      for ( int j = 0; j < (int) powers_deriv.size(); ++j )
        {
        if ( powers_deriv[j] < 1 )
          continue;
        
        dderiv += coef_deriv[j] * ::pow( z[2 * j], powers_deriv[j] - 1 );
        }
#if 0
      for ( int j = 0; j < (int) powers.size(); ++j )
        {
        value + =coef[j] * ::pow( z[2 * j], powers[j][0]);
        }
      tempvec.push_back( value );
#endif // 0

      if ( dderiv < 0 )
        {
        types.push_back( 2 );
        //tempvec.push_back(2);
        }
      else if ( dderiv > 0 )
        {
        types.push_back( 0 );
        //tempvec.push_back(0);
        }
      else 
        {
        types.push_back( 1 );
        //tempvec.push_back(1);
        }
      sysroots.push_back( tempvec );
      tempvec.clear();
      }
    }

  delete [] c;
  delete [] z;
  coef_deriv.clear();
  powers_deriv.clear();
}

