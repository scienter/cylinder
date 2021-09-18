#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void loadLaser2D_Yee_Pukhov(Domain *D,LaserList *L,double t);
void loadLaser_Split(Domain *D,LaserList *L,double t);
void calGaussianPulse(Domain *D,LaserList *L,double ***real,double ***img,double longitudinal,double t);
void calGaussianPulse_Add(Domain *D,LaserList *L,double ***real,double ***img,double ampli,double t);
void calFilePulse(Domain *D,LaserList *L,double ***ErR,double ***ErI,double ***EpR,double ***EpI,double longitudinal,double t);

void loadLaser(Domain *D,LaserList *L,double t)
{
  if(D->boostOn==OFF)
  {
    switch(D->fieldType)  {
    case Yee :
    case NoCherenkov :
      loadLaser2D_Yee_Pukhov(D,L,t);
      break;
    case Split :
      loadLaser_Split(D,L,t);
      break;
    default :
      printf("In loadLaser, what is field_type? and what is dimension?\n");
    }
  }
}

void loadLaser_Split(Domain *D,LaserList *L,double t)
{
   double rU,rD,longitudinal,t0,flat,phi0;
   double zR,w0,w,phi,omega,kx,pphi,amp,dt;
   double x,y,z,r2,w2,sinPh,cosPh;
   int istart,iend,jstart,jend;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   istart=D->istart; iend=D->iend;
   jstart=D->jstart; jend=D->jend;
   dt=D->dt/D->dtRatio;

   rU=L->rU/D->divisionLambda*dt;
   rD=L->rD/D->divisionLambda*dt;
   flat=L->flat/D->divisionLambda*L->lambda/D->lambda*dt;

   t0=2*rU+L->retard;
   zR=L->rayleighLength;
   w0=L->beamWaist;  
   omega=2*pi*L->omega/D->omega;
   kx=2*pi*D->lambda/L->lambda;

   if(t<2*rU)
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   else if(t>=2*rU && t<2*rU+flat) 
      longitudinal=L->amplitude*1.0;
   else if(t>=2*rU+flat && t<2*rU+flat+2*rD) 
      longitudinal=L->amplitude*exp(-(t-t0)*(t-t0)/rD/rD);
   else if(t>=2*rU+2*rD+flat) 
      longitudinal=0.0;

   positionX=L->loadPointX+istart-D->minXSub;
   if(positionX>=D->minXSub && positionX<D->maxXSub)
     laserOK=1;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;

   if(laserOK==1)
   {
     if(L->polarity==2) {
       x=-L->focus;   
       w=w0*sqrt(1.0+x*x/zR/zR);
       phi=atan(x/zR);
       phi0=0.0;
       cosPh=cos(phi0); sinPh=sin(phi0);
       w2=w*w;
       for(j=jstart; j<jend; j++) {
         y=(j-jstart+D->minYSub+0.5)*D->dr;
         r2=y*y;
         pphi=x/zR*r2/w2-phi+kx*x-omega*t;
         amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
         D->PrR[1][positionX][j]=2.0*amp*cosPh;
//         D->PrI[1][positionX][j]=0.0;
//         D->PlR[1][positionX][j]=0.0;
//         D->PlI[1][positionX][j]=0.0;
//         D->SrR[1][positionX][j]=0.0;
//         D->SlR[1][positionX][j]=0.0;
         D->SrI[1][positionX][j]=2.0*amp*cosPh;
//         D->SlI[1][positionX][j]=0.0;
       }

     } else if(L->polarity==3) {
/*
       k=0;
       w2=w*w;
       for(j=0; j<jend+3; j++) {
         y=(j-jstart+D->minYSub-jC)*D->dy;
         r2=y*y;
         pphi=x/zR*r2/w2-0.5*phi+kx*x-omega*t;
         amp=longitudinal*sqrt(w0/w)*exp(-r2/w2)*sin(pphi);
         D->Ez[positionX][j][k]=amp;            
       }
*/
     }  
   } else ;    //end of field is OK

}

void loadLaser2D_Yee_Pukhov(Domain *D,LaserList *L,double t)
{
   double rU,rD,amp=0.0,t0,flat,retard,phi0;
   double zR,w0,w,phi,omega,kx,pphi,dt,gdd,alpha,beta,gamma,tt,A;
   double x,y,z,r2,w2,x0,x1,x2,sig,sig1,sig2;
   int positionX,rank,i,j,k,jC,kC,laserOK=0;    
   int myrank, nTasks;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   dt=D->dt/D->dtRatio;

   rU=L->rU/D->divisionLambda*dt;
   rD=L->rD/D->divisionLambda*dt;
   retard=L->retard/D->divisionLambda*dt;
   flat=L->flat/D->divisionLambda*L->lambda/D->lambda*dt;

   t0=2*rU+retard;
   if(t>=retard && t<t0) {
      amp=L->amplitude*exp(-(t-t0)*(t-t0)/rU/rU);
   }
   else if(t>=t0 && t<t0+flat) {
      amp=L->amplitude*1.0;
   }
   else if(t>=t0+flat && t<t0+flat+2*rD) {
      amp=L->amplitude*exp(-(t-t0-flat)*(t-t0-flat)/rD/rD);
   }
   else if(t>=t0+2*rD+flat) {
      amp=0.0;
   }


   positionX=L->loadPointX+D->istart-D->minXSub;
   if(positionX>=D->minXSub && positionX<D->maxXSub) laserOK=1; else ;
//   if(positionX>D->minXSub && positionX<=D->maxXSub &&
//      jC-D->minYSub>=jstart && jC-D->minYSub<jend &&
//      kC-D->minZSub>=kstart && kC-D->minZSub<kend)
//     laserOK=1;

   // frequency chirping
   gdd=L->gdd;
   if(gdd==0) tt=t;
   else {
     A=1.0+gdd-2.0*gdd/(rU+rD)*rU;
     beta=A-gdd;
     alpha=gdd/(rU+rD);
     gamma=2.0*rU*(1.0-beta-alpha*rU);
     tt=0.5*alpha*t*t+beta*t+gamma;
   }
   if(laserOK==1) {
     if(L->add==OFF) {
       if(L->polarity==2) calGaussianPulse(D,L,D->ErR,D->EpI,amp,tt);
       else ;
     } else if(L->add==ON) {
       if(L->polarity==2) calGaussianPulse_Add(D,L,D->ErR,D->EpI,amp,tt);
       else ;
     } else if(L->add==File) {
       if(L->polarity==2) calFilePulse(D,L,D->ErR,D->ErI,D->EpR,D->EpI,amp,tt);
       else ;
     }
   } else ;    //end of field is OK

}

void calFilePulse(Domain *D,LaserList *L,double ***ErR,double ***ErI,double ***EpR,double ***EpI,double longitudinal,double t)
{
   int positionX,i,j,k,istart,iend,jstart,jend;
   double omega,tmp,tmp1;

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   omega=2*pi*L->omega/D->omega;
   positionX=L->loadPointX+istart-D->minXSub;

   for(j=0; j<jend+3; j++) {     
     tmp=longitudinal*D->laserI[j]*cos(-omega*t+D->laserPhase[j]);
     ErR[1][positionX][j]=tmp;
     EpI[1][positionX][j]=tmp;
//     tmp1=longitudinal*(D->laserI[j]+D->laserI[j-1])*0.5*cos(-omega*t+0.5*(D->laserPhase[j]+D->laserPhase[j-1]));
   }
}

void calGaussianPulse(Domain *D,LaserList *L,double ***real,double ***img,double longitudinal,double t)
{
   double zR,w0,w,phi,omega,pphi,phi0,amp;
   double x,y,r2,w2,x0,x1,x2,sig,sig1,sig2;
   int positionX,i,j,k,istart,iend,jstart,jend;    

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   zR=L->rayleighLength;
   w0=L->beamWaist;  
   omega=2*pi*L->omega/D->omega;

   positionX=L->loadPointX+istart-D->minXSub;
   
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   phi0=0.0;
   w2=w*w;
   for(j=jstart; j<jend; j++) {
     y=(j-jstart+D->minYSub+0.5)*D->dr;
     r2=y*y;
     pphi=x/zR*r2/w2-phi-omega*t;
     amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
     real[1][positionX][j]=amp;      
   }
   for(j=jstart; j<jend; j++) {
     y=(j-jstart+D->minYSub)*D->dr;
     r2=y*y;
     pphi=x/zR*r2/w2-phi-omega*t;
     amp=longitudinal*w0/w*exp(-r2/w2)*sin(pphi);
     img[1][positionX][j]=amp;            
   }
}

void calGaussianPulse_Add(Domain *D,LaserList *L,double ***real,double ***img,double ampli,double t)
{
   double zR,w0,w,phi,omega,pphi,phi0,amp;
   double x,y,r2,w2,x0,x1,x2,sig,sig1,sig2;
   int positionX,i,j,k,istart,iend,jstart,jend;    

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   zR=L->rayleighLength;
   w0=L->beamWaist;  
   omega=2*pi*L->omega/D->omega;

   positionX=L->loadPointX+istart-D->minXSub;
   
   x=-L->focus;   
   w=w0*sqrt(1.0+x*x/zR/zR);
   phi=atan(x/zR);
   phi0=0.0;
   w2=w*w;
   for(j=jstart; j<jend+3; j++) {
     y=(j-jstart+D->minYSub)*D->dr;
     r2=y*y;
     pphi=x/zR*r2/w2-phi-omega*t;
     amp=ampli*w0/w*exp(-r2/w2)*sin(pphi);
     real[1][positionX][j]+=amp*cos(phi0);            
//if(j==jstart) printf("Er,longitudinal=%g, amp=%g\n",longitudinal,amp); 
   }
   for(j=jstart+1; j<jend+3; j++) {
     y=(j-jstart+D->minYSub-0.5)*D->dr;
     r2=y*y;
     pphi=x/zR*r2/w2-phi-omega*t;
     amp=ampli*w0/w*exp(-r2/w2)*sin(pphi);
     img[1][positionX][j]+=amp*cos(phi0);            
   }
}

