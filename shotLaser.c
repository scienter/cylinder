#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <complex.h>
#include "mpi.h"

void shotLaser_Gaussian(Domain *D,LaserList *L);

void shotLaser(Domain *D,LaserList *L)
{
    if(L->mode==Gaussian)
      shotLaser_Gaussian(D,L);
    else ;
}

void shotLaser_Gaussian(Domain *D,LaserList *L)
{
   double rU,rD,longitudinal,t0,flat,x1,x2,elliptic;
   double zR,w0,w,phi,omega,kx,pphi,amp,focus,a0;
   double x,r,r2,w2,retard,positionX,dz,dr,dt,dtBydr;
   double upPrR,dnPrR,PrR,upSrI,dnSrI,SrI;
   double ***field1,***field2;	
   double ***field3,***field4;	
   int istart,iend,jstart,jend;
   int rank,i,j,minXSub,minYSub;
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   minXSub=D->minXSub; 
   minYSub=D->minYSub;
   dz=D->dz; dr=D->dr; dt=D->dt;

   rU=L->rU*L->lambda/D->lambda;  //D->lambda;
   rD=L->rD*L->lambda/D->lambda; //*D->lambda;
   flat=L->flat*L->lambda/D->lambda;
   retard=L->retard*L->lambda/D->lambda;

   a0=L->amplitude;
   focus=L->focus;
   zR=L->rayleighLength;	//normalized
   w0=L->beamWaist;  		//normalized
   positionX=L->loadPointX*D->dz;
   x1=positionX-flat*0.5-retard;
   x2=positionX+flat*0.5-retard;
   kx=2.0*M_PI/L->lambda*D->lambda;
   dtBydr=dt/dr;

   switch(D->fieldType)  {
   case Yee:
     field1=D->ErR; field2=D->EpI;
     field3=D->BrI; field4=D->BpR;
     break;
   case Split:
     field1=D->PrR; field2=D->SrI;
     field3=D->PrR; field4=D->SrI;
     break;
   }

   switch (L->add)  {
   case OFF :
     if(D->fieldType==Split) {
       for(i=0; i<iend+3; i++)  {
         x=(minXSub+i)*dz+retard;
         w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
         w2=w*w;
         if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
         else if(x>x1 && x<=x2) longitudinal=a0;
         else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
         phi=atan((x-focus)/zR);
         for(j=jstart; j<jend+3; j++)         {
           r=(j-jstart+minYSub+0.5)*dr;
           r2=r*r;
           pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
           amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
           field1[1][i][j]=2*amp;            
           field2[1][i][j]=2*amp;
         }
       }
/*
       for(i=1; i<iend+2; i++)  {
         for(j=jstart+1; j<jend+2; j++)         {
           r=j-jstart;
           upPrR=(field1[1][i-1][j]  +field1[1][i][j])*0.5; 
           dnPrR=(field1[1][i-1][j-1]+field1[1][i][j-1])*0.5; 
           upSrI=(field2[1][i-1][j]  +field2[1][i][j])*0.5; 
           dnSrI=(field2[1][i-1][j-1]+field2[1][i][j-1])*0.5; 
           PrR=0.5*(upPrR+dnPrR);
           SrI=0.5*(upSrI+dnSrI);
           field3[1][i][j]=(0.5/r*dtBydr*(PrR-SrI)+0.5*dtBydr*(upPrR-dnPrR));
           field4[1][i][j]=(0.5/r*dtBydr*(PrR-SrI)-0.5*dtBydr*(upSrI-dnSrI));
           field3[1][i][j]=(0.5/r*dtBydr*(PrR-SrI)+0.5*dtBydr*(upPrR-dnPrR));
           field4[1][i][j]=(0.5/r*dtBydr*(PrR-SrI)-0.5*dtBydr*(upSrI-dnSrI));
         }
       }
*/
     } else if(D->fieldType==Yee) {
       for(i=0; i<iend+3; i++)  {
         x=(minXSub+i)*dz+retard;
         w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
         w2=w*w;
         if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
         else if(x>x1 && x<=x2) longitudinal=a0;
         else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
         phi=atan((x-focus)/zR);
         for(j=jstart; j<jend+3; j++)         {
           r=(j-jstart+minYSub+0.5)*dr;
           r2=r*r;
           pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
           amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
           field1[1][i][j]=amp;            
           field2[1][i][j]=amp;
         }
       }
       for(i=1; i<iend+2; i++)
         for(j=jstart; j<jend+3; j++)         {
           field3[1][i][j]=-(field1[1][i][j]+field1[1][i+1][j])*0.5;
           field4[1][i][j]= (field2[1][i][j]+field2[1][i+1][j])*0.5;
         }
     }
     break;
   case ON :
     for(i=0; i<iend+3; i++)  {
       x=(minXSub+i)*dz+retard;
       w=w0*sqrt(1.0+(x-focus)*(x-focus)/zR/zR);
       w2=w*w;
       if(x<=x1)      longitudinal=exp(-(x-positionX)*(x-positionX)/rU/rU)*a0;
       else if(x>x1 && x<=x2) longitudinal=a0;
       else if(x>x2)  longitudinal=exp(-(x-positionX)*(x-positionX)/rD/rD)*a0;
       phi=atan((x-focus)/zR);
       for(j=jstart; j<jend+3; j++)         {
         r=(j-jstart+minYSub+0.5)*dr;
         r2=r*r;
         pphi=(x-focus)/zR*r2/w2-0.5*phi+kx*x;           
         amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
         field1[1][i][j]+=2*amp;            
         field2[1][i][j]+=2*amp;            
       }
     }
     break;
   default :
     printf("In shotLaser.c, what laser mode?\n");
     break;
   }		//End of switch
}

