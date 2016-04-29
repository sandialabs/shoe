/*
 * Copyright 2012 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
#include "Elements/VelhoFunctions.h"
#include "Elements/Generic.h"

double vtkShoeElemVelhoFunctions::ifunction(double coord[3],vtkShoeMeshIterator& cell, double Isovalue,int field_no)
{
  double x[3] = {0.,0.,0.};
  cell.GetCellOps()->EvaluateField(x,cell,coord,field_no);
//  cout<<x[0]<<endl;
  return x[0]-Isovalue;
  
}
void vtkShoeElemVelhoFunctions::rootfind(double p11[3], double p22[3], double p3[3], double Isovalue, vtkShoeMeshIterator& cell,int field_no)
{
  
  double midpt[3],p1[3],p2[3];
  double product;
  for(int cnt=0;cnt<3;cnt++)
  {
    p1[cnt]=p11[cnt];
    p2[cnt]=p22[cnt];
  }

  while(1)
  {
    for(int cnt=0;cnt<3;cnt++)
      midpt[cnt]=(p1[cnt]+p2[cnt])/2.;
    if(vtkMath::Distance2BetweenPoints(midpt,p1)<.00001)
      break;
    
    product = ifunction(midpt,cell,Isovalue,field_no)*ifunction(p1,cell,Isovalue,field_no);
//    cout<<product<<endl;
    if(product==0)
      break;
    if(product<0)
    {
      p2[0]=midpt[0];
      p2[1]=midpt[1];
      p2[2]=midpt[2];
    }
    else
    {
//      cout<<"changing p1"<<endl;
      p1[0]=midpt[0];
      p1[1]=midpt[1];
      p1[2]=midpt[2];
    }
  }
  p3[0]=midpt[0];
  p3[1]=midpt[1];
  p3[2]=midpt[2];
}

void vtkShoeElemVelhoFunctions::tri(int v0,int v1,int v2,int v3,double tetraPts[4][3],vtkShoeMeshIterator& cell,double Isovalue,int field_no,vtkAdaptiveTessellator* tess_in)
{
  double pt[3][6];
  rootfind(tetraPts[v0],tetraPts[v1],&pt[0][3],Isovalue,cell,field_no);
  rootfind(tetraPts[v0],tetraPts[v2],&pt[1][3],Isovalue,cell,field_no);
  rootfind(tetraPts[v0],tetraPts[v3],&pt[2][3],Isovalue,cell,field_no);

  for(int cnt=0;cnt<3;cnt++)
  {
    cell.GetCellOps()->EvaluateGeometry(&pt[cnt][0],cell,&pt[cnt][3]);
  }

  tess_in->AdaptivelySample2Facet(pt[0],pt[1],pt[2]);
}

void vtkShoeElemVelhoFunctions::quad(int v0,int v1,int v2,int v3,double tetraPts[4][3],vtkShoeMeshIterator& cell,double Isovalue,int field_no,vtkAdaptiveTessellator* tess_in)
{
  double pt[4][6];
  int intPtIds[3];
  rootfind(tetraPts[v0],tetraPts[v2],&pt[0][3],Isovalue,cell,field_no);
  rootfind(tetraPts[v0],tetraPts[v3],&pt[1][3],Isovalue,cell,field_no);
  rootfind(tetraPts[v1],tetraPts[v3],&pt[2][3],Isovalue,cell,field_no);
  for(int cnt=0;cnt<3;cnt++)
  {
    cell.GetCellOps()->EvaluateGeometry(&pt[cnt][0],cell,&pt[cnt][3]);
  }
  tess_in->AdaptivelySample2Facet(pt[0],pt[1],pt[2]);
  rootfind(tetraPts[v1],tetraPts[v2],&pt[3][3],Isovalue,cell,field_no);
  cell.GetCellOps()->EvaluateGeometry(&pt[3][0],cell,&pt[3][3]);
  int intPtIds1[3];
  intPtIds1[0]=intPtIds[0];
  intPtIds1[2] = intPtIds[2];

  tess_in->AdaptivelySample2Facet(pt[0],pt[3],pt[2]);
  

}




void vtkShoeElemVelhoFunctions::Velho(vtkTetra *tetra,vtkShoeMeshIterator& cell,double Isovalue,int field_no,vtkAdaptiveTessellator* tess_in)
{
  double tetraPts[4][3];
  double scalPts[4];
  for(int cnt=0;cnt<4;cnt++)
  {
    tetra->GetPoints()->GetPoint(cnt,tetraPts[cnt]);
  }
  
  for(int cnt=0;cnt<4;cnt++)
  {
    cell.GetCellOps()->EvaluateField(&(scalPts[cnt]),cell,tetraPts[cnt],field_no);
    scalPts[cnt]+=-Isovalue;
  }

  if(scalPts[0]<0)
  {
      if(scalPts[1]<0)
      {
         if(scalPts[2]<0)
         {
      if(scalPts[3]<0)
      { 
        return;
      }
      else
        tri(3,0,1,2,tetraPts,cell,Isovalue,field_no,tess_in);
         }
         else
         {
      if(scalPts[3]<0)
        tri(2,3,1,0,tetraPts,cell,Isovalue,field_no,tess_in);
      else
        quad(2,3,0,1,tetraPts,cell,Isovalue,field_no,tess_in);
         }
      }
      else
      {
         if(scalPts[2]<0)
         {
      if(scalPts[3]<0) 
        tri(1,3,2,0,tetraPts,cell,Isovalue,field_no,tess_in);
      else
        quad(1,3,2,0,tetraPts,cell,Isovalue,field_no,tess_in);
         }
         else
         {
      if(scalPts[3]<0)
        quad(1,2,0,3,tetraPts,cell,Isovalue,field_no,tess_in);
      else
        tri(0,3,2,1,tetraPts,cell,Isovalue,field_no,tess_in);
         }
      }
  }
  else
  {
      if(scalPts[1]<0)
      {
         if(scalPts[2]<0)
         {
      if(scalPts[3]<0)
        tri(0,3,2,1,tetraPts,cell,Isovalue,field_no,tess_in);
      else
        quad(0,3,1,2,tetraPts,cell,Isovalue,field_no,tess_in);
         }
         else
         {
      if(scalPts[3]<0)
        quad(0,2,3,1,tetraPts,cell,Isovalue,field_no,tess_in);
      else
        tri(1,3,0,2,tetraPts,cell,Isovalue,field_no,tess_in);
         }
      }
      else
      {
         if(scalPts[2]<0)
         {
      if(scalPts[3]<0) 
        quad(0,1,2,3,tetraPts,cell,Isovalue,field_no,tess_in);
      else
        tri(2,3,1,0,tetraPts,cell,Isovalue,field_no,tess_in);
         }
         else
         {
      if(scalPts[3]<0)
        tri(3,0,1,2,tetraPts,cell,Isovalue,field_no,tess_in);
      else
      {
        return;
      }
         }
      }
  }

}
