#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void absorb_U(Domain *D)
{
  double tmp,y,LdU,rr,rd,upL;
  int j,jstart,jend;

  jstart=D->jstart;    jend=D->jend;

  upL=(double)(D->ny-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;
  rd=D->pmld;

  for(j=jstart; j<jend; j++) {
    y=j-jstart;
    if(y>=upL && LdU>0)  {
      tmp=(y-upL)/LdU;
      D->upr[j]=1.0-rr*rr*tmp*tmp;
      D->upd[j]=1.0-rd*rd*tmp*tmp;
    } else {
      D->upr[j]=1.0;
      D->upd[j]=1.0;
    }
  }

}

void absorb_R(Domain *D,double *rtr,double *rtd,double x,double rL,double LdR,double rr,double rd)
{
  double tmp;

  if(x>=rL)  {
    tmp=(x-rL)/LdR;
    *rtr=1.0-rr*rr*tmp*tmp;
    *rtd=1.0-rd*rd*tmp*tmp;
  } else  {  
    *rtr=1.0;
    *rtd=1.0;
  }
}
/*
void PML_Esolve2D_Left(Domain *D)
{
  int i,j,m,istart,iend,jstart,jend,left,numMode;
  double sigZ,coef1,coef2,dtBydr,dtBydz,dt,[mlr,r;
 [
  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;
  
  left=D->pmlCellLeft;  pmlr=D->pmlr;
  dtBydz=D->dt/D->dz; dtBydr=D->dt/D->dr; dt=D->dt;

  for(m=0; m<numMode; m++)
    for(i=0; i<left-1; i++)  
      for(j=jstart+1; j<jend; j++)
      {
        r=(double)(j-jstart);
        sigZ=(pmlr*i)*(pmlr*i);
        coef2=1.0/(1.0+dt*sigZ); coef1=(1.0-dt*sigZ)*coef2;
        D->lpml[m][i][j].ErR=coef1*D->lpml[m][i][j].ErR+coef2*(-dtBydz*(D->lpml[m][i][j].BpzR+D->lpml[m][i][j].BprR-D->lpml[m][i+1][j].BpzR-D->lpml[m][i+1][j].BprR)-dtBydr/r*m*D->lpml[m][i][j].BzI);
        D->lpml[m][i][j].ErI=coef1*D->lpml[m][i][j].ErI+coef2*(-dtBydz*(D->lpml[m][i][j].BpzI+D->lpml[m][i][j].BprI-D->lpml[m][i+1][j].BpzI-D->lpml[m][i+1][j].BprI)+dtBydr/r*m*D->lpml[m][i][j].BzR);
  
        D->lpml[m][i][j].EprR+=-dtBydr*(D->lpml[m][i][j].BzR-D->lpml[m][i][j-1].BzR);
        D->lpml[m][i][j].EprI+=-dtBydr*(D->lpml[m][i][j].BzI-D->lpml[m][i][j-1].BzI);
        D->lpml[m][i][j].EpzR=coef1*D->lpml[m][i][j].EpzR+dtBydz*coef2*(D->lpml[m][i][j].BrR-D->lpml[m][i+1][j].BrR);
        D->lpml[m][i][j].EpzI=coef1*D->lpml[m][i][j].EpzI+dtBydz*coef2*(D->lpml[m][i][j].BrI-D->lpml[m][i+1][j].BrI);

        D->lpml[m][i][j].EzR+=dtBydr*(1.0/(r-0.5)*m*D->lpml[m][i][j].BrI+1.0/(r-0.5)*(r*(D->lpml[m][i][j].BpzR+D->lpml[m][i][j].BprR)-(r-1.0)*(D->lpml[m][i][j-1].BpzR+D->lpml[m][i][j-1].BprR)));
        D->lpml[m][i][j].EzI+=dtBydr*(-1.0/(r-0.5)*m*D->lpml[m][i][j].BrR+1.0/(r-0.5)*(r*(D->lpml[m][i][j].BpzI+D->lpml[m][i][j].BprI)-(r-1.0)*(D->lpml[m][i][j-1].BpzI+D->lpml[m][i][j-1].BprI)));
      }

  //treatment of center
  m=1;
  j=jstart;
  for(i=0; i<left-1; i++)  
  {
    sigZ=(pmlr*i)*(pmlr*i);
    coef2=1.0/(1.0+dt*sigZ); coef1=(1.0-dt*sigZ)*coef2;
    D->lpml[m][i][j].ErR=coef1*D->lpml[m][i][j].ErR+coef2*(-dtBydz*(D->lpml[m][i][j].BpzR+D->lpml[m][i][j].BprR-D->lpml[m][i+1][j].BpzR-D->lpml[m][i+1][j].BprR)-dtBydr*m*D->lpml[m][i][j+1].BzI);
    D->lpml[m][i][j].ErI=coef1*D->lpml[m][i][j].ErI+coef2*(-dtBydz*(D->lpml[m][i][j].BpzI+D->lpml[m][i][j].BprI-D->lpml[m][i+1][j].BpzI-D->lpml[m][i+1][j].BprI)+dtBydr*m*D->lpml[m][i][j+1].BzR);
  }
}

void copyPML_Left_E(Domain *D,int left)
{
  int m,i,j,istart,iend,jstart,jend,numMode;

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;

  for(m=0; m<numMode; m++)  {
    for(i=0; i<2; i++) 
      for(j=jstart; j<jend; j++) {
        D->ErR[m][istart-i-1][j]=D->lpml[m][i][j].ErR;
        D->ErI[m][istart-i-1][j]=D->lpml[m][i][j].ErI;
        D->EpR[m][istart-i-1][j]=D->lpml[m][i][j].EpzR+D->lpml[m][i][j].EprR;
        D->EpI[m][istart-i-1][j]=D->lpml[m][i][j].EpzI+D->lpml[m][i][j].EprI;
        D->EzR[m][istart-i-1][j]=D->lpml[m][i][j].EzR;
        D->EzI[m][istart-i-1][j]=D->lpml[m][i][j].EzI;
      }
  }
}

void PML_Bsolve2D_Left(Domain *D)
{
  int i,j,m,istart,iend,jstart,jend,left,numMode;
  double sigZ,coef1,coef2,dtBydr,dtBydz,dt,pmlr,r;
  
  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;
  
  left=D->pmlCellLeft;  pmlr=D->pmlr;
  dtBydz=D->dt/D->dz; dtBydr=D->dt/D->dr; dt=D->dt;

  for(m=0; m<numMode; m++)
    for(i=0; i<left; i++)  
      for(j=jstart+1; j<jend; j++)
      {
        r=(double)(j-jstart);
        D->lpml[m][i][j].BzR+=dtBydr*(-1.0/r*m*D->lpml[m][i][j].ErI-1.0/r*((r+0.5)*(D->lpml[m][i][j+1].EpzR+D->lpml[m][i][j+1].EprR)-(r-0.5)*(D->lpml[m][i][j].EpzR+D->lpml[m][i][j].EprR)));
        D->lpml[m][i][j].BzI+=dtBydr*(1.0/r*m*D->lpml[m][i][j].ErR-1.0/r*((r+0.5)*(D->lpml[m][i][j+1].EpzI+D->lpml[m][i][j+1].EprI)-(r-0.5)*(D->lpml[m][i][j].EpzI+D->lpml[m][i][j].EprI)));

        D->lpml[m][i][j].BprR+=dtBydr*(D->lpml[m][i][j+1].EzR-D->lpml[m][i][j].EzR);
        D->lpml[m][i][j].BprI+=dtBydr*(D->lpml[m][i][j+1].EzI-D->lpml[m][i][j].EzI);
      }

  for(m=0; m<numMode; m++)
    for(i=0; i<left-1; i++)  
      for(j=jstart+1; j<jend; j++)
      {
        r=(double)(j-jstart);
        sigZ=(pmlr*(i+0.5))*(pmlr*(i+0.5));
        coef2=1.0/(1.0+dt*sigZ); coef1=(1.0-dt*sigZ)*coef2;
        D->lpml[m][i][j].BrR=coef1*D->lpml[m][i][j].BrR+coef2*(dtBydz*(D->lpml[m][i][j].EpzR+D->lpml[m][i][j].EprR-D->lpml[m][i+1][j].EpzR-D->lpml[m][i+1][j].EprR)+dtBydr/(r-0.5)*m*D->lpml[m][i][j].EzI);
        D->lpml[m][i][j].BrI=coef1*D->lpml[m][i][j].BrI+coef2*(dtBydz*(D->lpml[m][i][j].EpzI+D->lpml[m][i][j].EprI-D->lpml[m][i+1][j].EpzI-D->lpml[m][i+1][j].EprI)-dtBydr/(r-0.5)*m*D->lpml[m][i][j].EzR);
  
        D->lpml[m][i][j].BpzR=coef1*D->lpml[m][i][j].BpzR-dtBydz*coef2*(D->lpml[m][i][j].ErR-D->lpml[m][i+1][j].ErR);
        D->lpml[m][i][j].BpzI=coef1*D->lpml[m][i][j].BpzI-dtBydz*coef2*(D->lpml[m][i][j].ErI-D->lpml[m][i+1][j].ErI);
      }

  i=0;
  for(m=0; m<numMode; m++)
    for(j=jstart+1; j<jend; j++)
      {
        r=(double)(j-jstart);
        sigZ=(pmlr*(i+0.5))*(pmlr*(i+0.5));
        coef2=1.0/(1.0+dt*sigZ); coef1=(1.0-dt*sigZ)*coef2;
        D->lpml[m][i][j].BrR=coef1*D->lpml[m][i][j].BrR+coef2*(dtBydz*(D->EpR[m][istart][j]-D->lpml[m][i][j].EpzR-D->lpml[m][i][j].EprR)+dtBydr/(r-0.5)*m*D->lpml[m][i][j].EzI);
        D->lpml[m][i][j].BrI=coef1*D->lpml[m][i][j].BrI+coef2*(dtBydz*(D->EpI[m][istart][j]-D->lpml[m][i][j].EpzI-D->lpml[m][i][j].EprI)-dtBydr/(r-0.5)*m*D->lpml[m][i][j].EzR);
  
        D->lpml[m][i][j].BpzR=coef1*D->lpml[m][i][j].BpzR-dtBydz*coef2*(D->ErR[m][istart][j]-D->lpml[m][i][j].ErR);
        D->lpml[m][i][j].BpzI=coef1*D->lpml[m][i][j].BpzI-dtBydz*coef2*(D->ErI[m][istart][j]-D->lpml[m][i][j].ErI);
      }

  //treatment of center
  j=jstart;
  for(i=0; i<left; i++)  
    D->lpml[0][i][j].BzR+=-4.0*dtBydr*(D->lpml[0][i][j+1].BpzR+D->lpml[0][i][j+1].BprR);

  i=0;
  D->lpml[1][i][j].BprR+=dtBydr*(4.0*D->lpml[1][i][j+1].EzR-D->lpml[1][i][j+2].EzR);
  D->lpml[1][i][j].BprI+=dtBydr*(4.0*D->lpml[1][i][j+1].EzI-D->lpml[1][i][j+2].EzI);
  sigZ=(pmlr*(i+0.5))*(pmlr*(i+0.5));
  coef2=1.0/(1.0+dt*sigZ); coef1=(1.0-dt*sigZ)*coef2;
  D->lpml[1][i][j].BpzR=coef1*D->lpml[1][i][j].BpzR-dtBydz*coef2*(D->ErR[1][istart][j]-D->lpml[1][i][j].ErR);
  D->lpml[1][i][j].BpzI=coef1*D->lpml[1][i][j].BpzI-dtBydz*coef2*(D->ErI[1][istart][j]-D->lpml[1][i][j].ErI);

  for(i=1; i<left-1; i++) {
    sigZ=(pmlr*(i+0.5))*(pmlr*(i+0.5));
    coef2=1.0/(1.0+dt*sigZ); coef1=(1.0-dt*sigZ)*coef2;
    D->lpml[1][i][j].BpzR=coef1*D->lpml[1][i][j].BpzR-dtBydz*coef2*(D->lpml[1][i][j].ErR-D->lpml[1][i+1][j].ErR);
    D->lpml[1][i][j].BpzI=coef1*D->lpml[1][i][j].BpzI-dtBydz*coef2*(D->lpml[1][i][j].ErI-D->lpml[1][i+1][j].ErI);
  }

}

//lala
void copyPML_Left_B(Domain *D,int left)
{
  int m,i,j,istart,iend,jstart,jend,numMode;

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;

  for(m=0; m<numMode; m++)  {
    for(i=0; i<2; i++) 
      for(j=jstart; j<jend; j++) {
        D->BzR[m][istart-i-1][j]=D->lpml[m][i][j].BzR;
        D->BzI[m][istart-i-1][j]=D->lpml[m][i][j].BzI;
        D->BrR[m][istart-i-1][j]=D->lpml[m][i][j].BrR;
        D->BrI[m][istart-i-1][j]=D->lpml[m][i][j].BrI;
        D->BpR[m][istart-i-1][j]=D->lpml[m][i][j].BpzR+D->lpml[m][i][j].BprR;
        D->BpI[m][istart-i-1][j]=D->lpml[m][i][j].BpzI+D->lpml[m][i][j].BprI;
      }
  }
}
*/

/*
void PML_Esolve2D_Up(Domain *D)
{
  int i,j,m,istart,iend,jstart,jend,up,numMode;
  double sigY,coef1,coef2,dtBydr,dtBydz,dt,pmlr,r;
  
  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;
  
  up=D->pmlCellUp;  pmlr=D->pmlr;
  dtBydz=D->dt/D->dz; dtBydr=D->dt/D->dr; dt=D->dt;

  //Here pml.Exy=Ex, pml.Eyz=Eyi, pml.Bxy=Bx, pml.Byz=By
  j=0;
  r=(double)(jend-1-jstart);
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      sigY=pmlr*j*j;
      coef2=1.0/(1+dt*sigY); coef1=(1.0-dt*sigY)*coef2;

      D->upml[m][i][j].EzR=coef1*D->upml[m][i][j].EzR+dtBydr*coef2*(1.0/(r-0.5)*m*D->BrI[m][i][j]+(0.5/(r-0.5)+1.0)*(D->upml[m][i][j].BpzR+D->upml[m][i][j].BprR)+(0.5/(r-0.5)-1.0)*D->BpR[m][i][jend-2]);
      D->upml[m][i][j].ErR+=dtBydz*(-(D->upml[m][i][j].BpzR+D->upml[m][i][j].BprR-D->upml[m][i-1][j].BpzR-D->upml[m][i-1][j].BprR)-1.0/r*m*D->upml[m][i][j].BzI);
      D->upml[m][i][j].EpzR+=dtBydz*(D->BrR[m][i][jend-1]-D->BrR[m][i-1][jend-1]);
      D->upml[m][i][j].EprR=coef1*D->upml[m][i][j].EprR+dtBydr*coef2*(-(D->upml[m][i][j].BzR-D->BzR[m][i][jend-2]));
  
      D->upml[m][i][j].EzI=coef1*D->upml[m][i][j].EzI+dtBydr*coef2*(-1.0/(r-0.5)*m*D->BrR[m][i][j]+(0.5/(r-0.5)+1.0)*(D->upml[m][i][j].BpzI+D->upml[m][i][j].BprI)+(0.5/(r-0.5)-1.0)*D->BzI[m][i][jend-2]);
      D->upml[m][i][j].ErI+=dtBydz*(-(D->upml[m][i][j].BpzI+D->upml[m][i][j].BprI-D->upml[m][i-1][j].BpzI-D->upml[m][i-1][j].BprI)+1.0/r*m*D->upml[m][i][j].BzR);
      D->upml[m][i][j].EpzI+=dtBydz*(D->BrI[m][i][jend-1]-D->BrI[m][i-1][jend-1]);
      D->upml[m][i][j].EprI=coef1*D->upml[m][i][j].EprI+dtBydr*coef2*(-(D->upml[m][i][j].BzI-D->BzI[m][i][jend-2]));
    }

  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)  
      for(j=1; j<up; j++)
      {
        r=(double)(jend-1-jstart+j);
        sigY=pmlr*j*j;
        coef2=1.0/(1+dt*sigY); coef1=(1.0-dt*sigY)*coef2;

        D->upml[m][i][j].EzR=coef1*D->upml[m][i][j].EzR+dtBydr*coef2*(1.0/(r-0.5)*m*D->upml[m][i][j].BrI+(0.5/(r-0.5)+1.0)*(D->upml[m][i][j].BpzR+D->upml[m][i][j].BprR)+(0.5/(r-0.5)-1.0)*(D->upml[m][i][j-1].BpzR+D->upml[m][i][j-1].BprR));
        D->upml[m][i][j].ErR+=dtBydz*(-(D->upml[m][i][j].BpzR+D->upml[m][i][j].BprR-D->upml[m][i-1][j].BpzR-D->upml[m][i-1][j].BprR)-1.0/r*m*D->upml[m][i][j].BzI);
        D->upml[m][i][j].EpzR+=dtBydz*(D->upml[m][i][j].BrR-D->upml[m][i-1][j].BrR);
        D->upml[m][i][j].EprR=coef1*D->upml[m][i][j].EprR+dtBydr*coef2*(-(D->upml[m][i][j].BzR-D->upml[m][i][j-1].BzR));
 
        D->upml[m][i][j].EzI=coef1*D->upml[m][i][j].EzI+dtBydr*coef2*(-1.0/(r-0.5)*m*D->upml[m][i][j].BrR+(0.5/(r-0.5)+1.0)*(D->upml[m][i][j].BpzI+D->upml[m][i][j].BprI)+(0.5/(r-0.5)-1.0)*(D->upml[m][i][j-1].BpzI+D->upml[m][i][j-1].BprI));
        D->upml[m][i][j].ErI+=dtBydz*(-(D->upml[m][i][j].BpzI+D->upml[m][i][j].BprI-D->upml[m][i-1][j].BpzI-D->upml[m][i-1][j].BprI)+1.0/r*m*D->upml[m][i][j].BzR);
        D->upml[m][i][j].EpzI+=dtBydz*(D->upml[m][i][j].BrI-D->upml[m][i-1][j].BrI);
        D->upml[m][i][j].EprI=coef1*D->upml[m][i][j].EprI+dtBydr*coef2*(-(D->upml[m][i][j].BzI-D->upml[m][i][j-1].BzI));
      }
}

void copyPML_Up_E(Domain *D)
{
  int m,i,j,istart,iend,jstart,jend,numMode;

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;

 for(m=0; m<numMode; m++)  
   for(i=istart; i<iend; i++) 
     for(j=0; j<3; j++) {
      D->EzR[m][i][jend-1+j]=D->upml[m][i][j].EzR;
      D->ErR[m][i][jend-1+j]=D->upml[m][i][j].ErR;
      D->EpR[m][i][jend-1+j]=D->upml[m][i][j].EpzR+D->upml[m][i][j].EprR;
      D->EzI[m][i][jend-1+j]=D->upml[m][i][j].EzI;
      D->ErI[m][i][jend-1+j]=D->upml[m][i][j].ErI;
      D->EpI[m][i][jend-1+j]=D->upml[m][i][j].EpzI+D->upml[m][i][j].EprI;
    }
}


void PML_Bsolve2D_Up(Domain *D)
{
  int i,j,m,istart,iend,jstart,jend,up,numMode;
  double sigY,coef1,coef2,dtBydr,dtBydz,dt,pmlr,r;
  double oldBzR,oldBzI,oldBrR,oldBrI,oldBpR,oldBpI;
  
  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;
  
  up=D->pmlCellUp;  pmlr=D->pmlr;
  dtBydz=D->dt/D->dz; dtBydr=D->dt/D->dr; dt=D->dt;

  //Here pml.Exy=Ex, pml.Eyz=Eyi, pml.Bxy=Bx, pml.Byz=By
  j=0;
  r=(double)(jend-1-jstart);
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      sigY=pmlr*(j+0.5)*(j+0.5);
      coef2=1.0/(1+dt*sigY); coef1=(1.0-dt*sigY)*coef2;

      oldBzR=D->BzR[m][i][jend-1];
      oldBzI=D->BzI[m][i][jend-1];
      oldBrR=D->BrR[m][i][jend-1];
      oldBrI=D->BrI[m][i][jend-1];
      oldBpR=D->BpR[m][i][jend-1];
      oldBpI=D->BpI[m][i][jend-1];
      D->upml[m][i][j].BzR=coef1*D->upml[m][i][j].BzR+dtBydr*coef2*(1.0/r*m*D->ErI[m][i][j]-(0.5/r+1.0)*(D->upml[m][i+1][j].EpzR+D->upml[m][i+1][j].EprR)-(0.5/r-1.0)*D->EpR[m][i][jend-1]);
      D->upml[m][i][j].BrR+=dtBydz*((D->EpR[m][i+1][jend-1]-D->EpR[m][i][jend-1])+1.0/(r-0.5)*m*D->EzI[m][i][jend-1]);
      D->upml[m][i][j].BpzR+=-dtBydz*(D->upml[m][i+1][j].ErR-D->upml[m][i][j].ErR);
      D->upml[m][i][j].BprR=coef1*D->upml[m][i][j].BprR+dtBydr*coef2*(D->upml[m][i][j+1].EzR-D->EzR[m][i][jend-1]);
  
      D->upml[m][i][j].BzI=coef1*D->upml[m][i][j].BzI+dtBydr*coef2*(-1.0/r*m*D->ErR[m][i][j]-(0.5/r+1.0)*(D->upml[m][i+1][j].EpzI+D->upml[m][i+1][j].EprI)-(0.5/r-1.0)*D->EpI[m][i][jend-1]);
      D->upml[m][i][j].BrI+=dtBydz*((D->EpI[m][i+1][jend-1]-D->EpI[m][i][jend-1])-1.0/(r-0.5)*m*D->EzR[m][i][jend-1]);
      D->upml[m][i][j].BpzI+=-dtBydz*(D->upml[m][i+1][j].ErI-D->upml[m][i][j].ErI);
      D->upml[m][i][j].BprI=coef1*D->upml[m][i][j].BprI+dtBydr*coef2*(D->upml[m][i][j+1].EzI-D->EzI[m][i][jend-1]);
  
      D->BzNowR[m][i][jend-1]=0.5*(oldBzR+D->upml[m][i][j].BzR);
      D->BzNowI[m][i][jend-1]=0.5*(oldBzI+D->upml[m][i][j].BzI);
      D->BrNowR[m][i][jend-1]=0.5*(oldBrR+D->upml[m][i][j].BrR);
      D->BrNowI[m][i][jend-1]=0.5*(oldBrI+D->upml[m][i][j].BrI);
      D->BpNowR[m][i][jend-1]=0.5*(oldBpR+D->upml[m][i][j].BpzR+D->upml[m][i][j].BprR);
      D->BpNowI[m][i][jend-1]=0.5*(oldBpI+D->upml[m][i][j].BpzI+D->upml[m][i][j].BprI);
    }

  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)  
      for(j=1; j<up; j++)
      {
        r=(double)(jend-1-jstart+j);
        sigY=pmlr*(j+0.5)*(j+0.5);
        coef2=1.0/(1+dt*sigY); coef1=(1.0-dt*sigY)*coef2;

        D->upml[m][i][j].BzR=coef1*D->upml[m][i][j].BzR+dtBydr*coef2*(-1.0/r*m*D->upml[m][i][j].ErI-(0.5/r+1.0)*(D->upml[m][i][j+1].EpzR+D->upml[m][i][j+1].EprR)-(0.5/r-1.0)*(D->upml[m][i][j].EpzR+D->upml[m][i][j].EprR));
        D->upml[m][i][j].BrR+=dtBydz*((D->upml[m][i+1][j].EpzR+D->upml[m][i+1][j].EprR-D->upml[m][i][j].EpzR-D->upml[m][i][j].EprR)+1.0/(r-0.5)*m*D->upml[m][i][j].EzI);
        D->upml[m][i][j].BpzR+=-dtBydz*(D->upml[m][i+1][j].ErR-D->upml[m][i][j].ErR);
        D->upml[m][i][j].BprR=coef1*D->upml[m][i][j].BprR+dtBydr*coef2*(D->upml[m][i][j+1].EzR-D->upml[m][i][j].EzR);
 
        D->upml[m][i][j].BzI=coef1*D->upml[m][i][j].BzI+dtBydr*coef2*(1.0/r*m*D->upml[m][i][j].ErR-(0.5/r+1.0)*(D->upml[m][i][j+1].EpzI+D->upml[m][i][j+1].EprI)-(0.5/r-1.0)*(D->upml[m][i][j].EpzI+D->upml[m][i][j].EprI));
        D->upml[m][i][j].BrI+=dtBydz*((D->upml[m][i+1][j].EpzI+D->upml[m][i+1][j].EprI-D->upml[m][i][j].EpzI-D->upml[m][i][j].EprI)-1.0/(r-0.5)*m*D->upml[m][i][j].EzR);
        D->upml[m][i][j].BpzI+=-dtBydz*(D->upml[m][i+1][j].ErI-D->upml[m][i][j].ErI);
        D->upml[m][i][j].BprI=coef1*D->upml[m][i][j].BprI+dtBydr*coef2*(D->upml[m][i][j+1].EzI-D->upml[m][i][j].EzI);
      }
}

void copyPML_B(Domain *D)
{
  int m,i,j,istart,iend,jstart,jend,numMode;

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode;

 for(m=0; m<numMode; m++)  
   for(i=istart; i<iend; i++) 
     for(j=0; j<3; j++) {
      D->BzR[m][i][jend-1+j]=D->upml[m][i][j].BzR;
      D->BrR[m][i][jend-1+j]=D->upml[m][i][j].BrR;
      D->BpR[m][i][jend-1+j]=D->upml[m][i][j].BpzR+D->upml[m][i][j].BprR;
      D->BzI[m][i][jend-1+j]=D->upml[m][i][j].BzI;
      D->BrI[m][i][jend-1+j]=D->upml[m][i][j].BrI;
      D->BpI[m][i][jend-1+j]=D->upml[m][i][j].BpzI+D->upml[m][i][j].BprI;
    }
}
*/
/*
void PML_Esolve2D_Right(Domain *D)
{
  int i,j,k,istart,iend,jstart,jend,rt;
  double sigX,sigY,coef1,coef2,dtBydx,dtBydy,dt,pmlr;
  
  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  rt=D->pmlCellRight;  pmlr=D->pmlr;
  dtBydx=D->dt/D->dx; dtBydy=D->dt/D->dy; dt=D->dt;

  //Here pml.Exy=Ex, pml.Eyz=Eyi, pml.Bxy=Bx, pml.Byz=By
  k=0; i=0;
  for(j=jstart; j<jend; j++)  {
    sigX=pmlr*i*i;
    coef2=1.0/(1+dt*sigX); coef1=(1.0-dt*sigX)*coef2;
    D->rpml[i][j][k].Exy+=dtBydy*(D->rpml[i][j][k].Bzx+D->rpml[i][j][k].Bzy-D->rpml[i][j-1][k].Bzx-D->rpml[i][j-1][k].Bzy);
    D->rpml[i][j][k].Eyz=coef1*D->rpml[i][j][k].Eyz-dtBydx*coef2*(D->rpml[i][j][k].Bzx+D->rpml[i][j][k].Bzy-D->Bz[iend-2][j][k]);
    D->rpml[i][j][k].Ezx=coef1*D->rpml[i][j][k].Ezx+dtBydx*coef2*(D->rpml[i][j][k].Byz-D->By[iend-2][j][k]);
    D->rpml[i][j][k].Ezy+=-dtBydy*(D->Bx[iend-1][j][k]-D->Bx[iend-1][j-1][k]);
  }

  for(i=1; i<rt; i++)  
    for(j=jstart; j<jend; j++)  {
      sigX=pmlr*i*i;
      coef2=1.0/(1+dt*sigX); coef1=(1.0-dt*sigX)*coef2;
      D->rpml[i][j][k].Exy+=dtBydy*(D->rpml[i][j][k].Bzx+D->rpml[i][j][k].Bzy-D->rpml[i][j-1][k].Bzx-D->rpml[i][j-1][k].Bzy);
      D->rpml[i][j][k].Eyz=coef1*D->rpml[i][j][k].Eyz-dtBydx*coef2*(D->rpml[i][j][k].Bzx+D->rpml[i][j][k].Bzy-D->rpml[i-1][j][k].Bzx-D->rpml[i-1][j][k].Bzy);
      D->rpml[i][j][k].Ezx=coef1*D->rpml[i][j][k].Ezx+dtBydx*coef2*(D->rpml[i][j][k].Byz-D->rpml[i-1][j][k].Byz);
      D->rpml[i][j][k].Ezy+=-dtBydy*(D->rpml[i][j][k].Bxy-D->rpml[i][j-1][k].Bxy);
    }
  if(D->L>1) {
    MPI_TransferPML_Yminus(D,D->rpml,D->pmlCellRight,1,3);
    MPI_TransferPML_Yplus(D,D->rpml,D->pmlCellRight,1,3);
  } else    ;
  //update calculation domain
  for(i=0; i<3; i++) 
    for(j=jstart; j<jend; j++)  {
      D->Ex[iend-1+i][j][k]=D->rpml[i][j][k].Exy;
      D->Ey[iend-1+i][j][k]=D->rpml[i][j][k].Eyz;
      D->Ez[iend-1+i][j][k]=D->rpml[i][j][k].Ezx+D->rpml[i][j][k].Ezy;
    }
}

void PML_Bsolve2D_Right(Domain *D)
{
  int i,j,k,istart,iend,jstart,jend,rt;
  double sigX,sigY,coef1,coef2,dtBydx,dtBydy,dt,pmlr;
  
  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  rt=D->pmlCellRight;  pmlr=D->pmlr;
  dtBydx=D->dt/D->dx; dtBydy=D->dt/D->dy; dt=D->dt;

  //Here pml.Exy=Ex, pml.Eyz=Eyi, pml.Bxy=Bx, pml.Byz=By
  for(i=0; i<rt-1; i++)  
    for(j=jstart; j<jend; j++)  {
      sigX=pmlr*i*i;
      coef2=1.0/(1+dt*sigX); coef1=(1.0-dt*sigX)*coef2;
      D->rpml[i][j][k].Bxy+=-dtBydy*(D->rpml[i][j+1][k].Ezx+D->rpml[i][j+1][k].Ezy-D->rpml[i][j][k].Ezx-D->rpml[i][j][k].Ezy);
      D->rpml[i][j][k].Byz=coef1*D->rpml[i][j][k].Byz+dtBydx*coef2*(D->rpml[i+1][j][k].Ezx+D->rpml[i+1][j][k].Ezy-D->rpml[i][j][k].Ezx-D->rpml[i][j][k].Ezy);
      D->rpml[i][j][k].Bzx=coef1*D->rpml[i][j][k].Bzx-dtBydx*coef2*(D->rpml[i+1][j][k].Eyz-D->rpml[i][j][k].Eyz);
      D->rpml[i][j][k].Bzy+=dtBydy*(D->rpml[i][j+1][k].Exy-D->rpml[i][j][k].Exy);
    }
  if(D->L>1) {
    MPI_TransferPML_Yminus(D,D->rpml,D->pmlCellRight,1,3);
    MPI_TransferPML_Yplus(D,D->rpml,D->pmlCellRight,1,3);
  } else    ;
  //update calculation domain
  for(i=0; i<3; i++) 
    for(j=jstart; j<jend; j++)  {
      D->Bx[iend-1+i][j][k]=D->upml[i][j][k].Bxy;
      D->By[iend-1+i][j][k]=D->upml[i][j][k].Byz;
      D->Bz[iend-1+i][j][k]=D->upml[i][j][k].Bzx+D->upml[i][j][k].Bzy;
    }
}
*/

/*
void MPI_TransferPML_Xplus(Domain *D,PML ***f1,int ny,int share)
{
    int i,j,m,num;
    int istart,iend,jstart,jend,numMode;
    int myrank, nTasks,rank,start;
    double *data;

    MPI_Status status;

    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=myrank/D->M;

    num=(share-1)*numMode*ny*16;  //6 is pml component
    data = (double *)malloc(num*sizeof(double ));
    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++) 
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)      { 
          data[start+0]=f1[m][iend-i][j].EzR; 
          data[start+1]=f1[m][iend-i][j].ErR; 
          data[start+2]=f1[m][iend-i][j].EpzR; 
          data[start+3]=f1[m][iend-i][j].EprR; 
          data[start+4]=f1[m][iend-i][j].EzI; 
          data[start+5]=f1[m][iend-i][j].ErI; 
          data[start+6]=f1[m][iend-i][j].EpzI; 
          data[start+7]=f1[m][iend-i][j].EprI; 
          data[start+8]=f1[m][iend-i][j].BzR; 
          data[start+9]=f1[m][iend-i][j].BrR; 
          data[start+10]=f1[m][iend-i][j].BpzR; 
          data[start+11]=f1[m][iend-i][j].BprR; 
          data[start+12]=f1[m][iend-i][j].BzI; 
          data[start+13]=f1[m][iend-i][j].BrI; 
          data[start+14]=f1[m][iend-i][j].BpzI; 
          data[start+15]=f1[m][iend-i][j].BprI; 
          start+=16;
        }

    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);
      start=0;
      for(m=0; m<numMode; m++) 
        for(i=1; i<share; i++)
          for(j=0; j<ny; j++)    {
            f1[m][istart-i][j].EzR=data[start+0];
            f1[m][istart-i][j].ErR=data[start+1];
            f1[m][istart-i][j].EpzR=data[start+2];
            f1[m][istart-i][j].EprR=data[start+3];
            f1[m][istart-i][j].EzI=data[start+4];
            f1[m][istart-i][j].ErI=data[start+5];
            f1[m][istart-i][j].EpzI=data[start+6];
            f1[m][istart-i][j].EprI=data[start+7];
            f1[m][istart-i][j].BzR=data[start+8];
            f1[m][istart-i][j].BrR=data[start+9];
            f1[m][istart-i][j].BpzR=data[start+10];
            f1[m][istart-i][j].BprR=data[start+11];
            f1[m][istart-i][j].BzI=data[start+12];
            f1[m][istart-i][j].BrI=data[start+13];
            f1[m][istart-i][j].BpzI=data[start+14];
            f1[m][istart-i][j].BprI=data[start+15];
            start+=16;
          }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++) 
      for(i=1; i<share; i++)
        for(j=0; j<ny; j++)      { 
          data[start+0]=f1[m][iend-i][j].EzR; 
          data[start+1]=f1[m][iend-i][j].ErR; 
          data[start+2]=f1[m][iend-i][j].EpzR; 
          data[start+3]=f1[m][iend-i][j].EprR; 
          data[start+4]=f1[m][iend-i][j].EzI; 
          data[start+5]=f1[m][iend-i][j].ErI; 
          data[start+6]=f1[m][iend-i][j].EpzI; 
          data[start+7]=f1[m][iend-i][j].EprI; 
          data[start+8]=f1[m][iend-i][j].BzR; 
          data[start+9]=f1[m][iend-i][j].BrR; 
          data[start+10]=f1[m][iend-i][j].BpzR; 
          data[start+11]=f1[m][iend-i][j].BprR; 
          data[start+12]=f1[m][iend-i][j].BzI; 
          data[start+13]=f1[m][iend-i][j].BrI; 
          data[start+14]=f1[m][iend-i][j].BpzI; 
          data[start+15]=f1[m][iend-i][j].BprI; 
          start+=16;
        }

    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
      start=0;
      for(m=0; m<numMode; m++) 
        for(i=1; i<share; i++)
          for(j=0; j<ny; j++)    {
            f1[m][istart-i][j].EzR=data[start+0];
            f1[m][istart-i][j].ErR=data[start+1];
            f1[m][istart-i][j].EpzR=data[start+2];
            f1[m][istart-i][j].EprR=data[start+3];
            f1[m][istart-i][j].EzI=data[start+4];
            f1[m][istart-i][j].ErI=data[start+5];
            f1[m][istart-i][j].EpzI=data[start+6];
            f1[m][istart-i][j].EprI=data[start+7];
            f1[m][istart-i][j].BzR=data[start+8];
            f1[m][istart-i][j].BrR=data[start+9];
            f1[m][istart-i][j].BpzR=data[start+10];
            f1[m][istart-i][j].BprR=data[start+11];
            f1[m][istart-i][j].BzI=data[start+12];
            f1[m][istart-i][j].BrI=data[start+13];
            f1[m][istart-i][j].BpzI=data[start+14];
            f1[m][istart-i][j].BprI=data[start+15];
            start+=16;
          }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(data);
}

void MPI_TransferPML_Xminus(Domain *D,PML ***f1,int ny,int share)
{
    int i,j,m,num,start,end;
    int istart,iend,jstart,jend,numMode;
    int myrank, nTasks,rank;
    double *data;

    MPI_Status status;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    rank=myrank/D->M;

    num=share*numMode*ny*16;  //6 is pml component
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++) 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++) {
          data[start+0]=f1[m][i+istart][j].EzR;
          data[start+1]=f1[m][i+istart][j].ErR;
          data[start+2]=f1[m][i+istart][j].EpzR;
          data[start+3]=f1[m][i+istart][j].EprR;
          data[start+4]=f1[m][i+istart][j].EzI;
          data[start+5]=f1[m][i+istart][j].ErI;
          data[start+6]=f1[m][i+istart][j].EpzI;
          data[start+7]=f1[m][i+istart][j].EprI;
          data[start+8]=f1[m][i+istart][j].BzR;
          data[start+9]=f1[m][i+istart][j].BrR;
          data[start+10]=f1[m][i+istart][j].BpzR;
          data[start+11]=f1[m][i+istart][j].BprR;
          data[start+12]=f1[m][i+istart][j].BzI;
          data[start+13]=f1[m][i+istart][j].BrI;
          data[start+14]=f1[m][i+istart][j].BpzI;
          data[start+15]=f1[m][i+istart][j].BprI;
          start+=16;
        }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
      for(m=0; m<numMode; m++) 
        for(i=0; i<share; i++)
          for(j=0; j<ny; j++) {
            f1[m][iend+i][j].EzR=data[start+0];
            f1[m][iend+i][j].ErR=data[start+1];
            f1[m][iend+i][j].EpzR=data[start+2];
            f1[m][iend+i][j].EprR=data[start+3];
            f1[m][iend+i][j].EzI=data[start+4];
            f1[m][iend+i][j].ErI=data[start+5];
            f1[m][iend+i][j].EpzI=data[start+6];
            f1[m][iend+i][j].EprI=data[start+7];
            f1[m][iend+i][j].BzR=data[start+8];
            f1[m][iend+i][j].BrR=data[start+9];
            f1[m][iend+i][j].BpzR=data[start+10];
            f1[m][iend+i][j].BprR=data[start+11];
            f1[m][iend+i][j].BzI=data[start+12];
            f1[m][iend+i][j].BrI=data[start+13];
            f1[m][iend+i][j].BpzI=data[start+14];
            f1[m][iend+i][j].BprI=data[start+15];
            start+=16;
          }
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++) 
      for(i=0; i<share; i++)
        for(j=0; j<ny; j++) {
          data[start+0]=f1[m][i+istart][j].EzR;
          data[start+1]=f1[m][i+istart][j].ErR;
          data[start+2]=f1[m][i+istart][j].EpzR;
          data[start+3]=f1[m][i+istart][j].EprR;
          data[start+4]=f1[m][i+istart][j].EzI;
          data[start+5]=f1[m][i+istart][j].ErI;
          data[start+6]=f1[m][i+istart][j].EpzI;
          data[start+7]=f1[m][i+istart][j].EprI;
          data[start+8]=f1[m][i+istart][j].BzR;
          data[start+9]=f1[m][i+istart][j].BrR;
          data[start+10]=f1[m][i+istart][j].BpzR;
          data[start+11]=f1[m][i+istart][j].BprR;
          data[start+12]=f1[m][i+istart][j].BzI;
          data[start+13]=f1[m][i+istart][j].BrI;
          data[start+14]=f1[m][i+istart][j].BpzI;
          data[start+15]=f1[m][i+istart][j].BprI;
          start+=16;
        }

    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);
      start=0;
      for(m=0; m<numMode; m++) 
        for(i=0; i<share; i++)
          for(j=0; j<ny; j++) {
            f1[m][iend+i][j].EzR=data[start+0];
            f1[m][iend+i][j].ErR=data[start+1];
            f1[m][iend+i][j].EpzR=data[start+2];
            f1[m][iend+i][j].EprR=data[start+3];
            f1[m][iend+i][j].EzI=data[start+4];
            f1[m][iend+i][j].ErI=data[start+5];
            f1[m][iend+i][j].EpzI=data[start+6];
            f1[m][iend+i][j].EprI=data[start+7];
            f1[m][iend+i][j].BzR=data[start+8];
            f1[m][iend+i][j].BrR=data[start+9];
            f1[m][iend+i][j].BpzR=data[start+10];
            f1[m][iend+i][j].BprR=data[start+11];
            f1[m][iend+i][j].BzI=data[start+12];
            f1[m][iend+i][j].BrI=data[start+13];
            f1[m][iend+i][j].BpzI=data[start+14];
            f1[m][iend+i][j].BprI=data[start+15];
            start+=16;
          }
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(data);
}
*/

/*
void MPI_TransferPML_Yminus(Domain *D,PML ***f1,int nx,int nz,int share)
{
    int i,j,k,numberData,start,end,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank;
    double *data;

    MPI_Status status;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    rank=(myrank%(D->M*D->N))%D->M;

    num=share*nx*nz*12;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)
        for(k=0; k<nz; k++) { 
          data[start]=f1[i][j+jstart][k].Exy;
          data[start+1]=f1[i][j+jstart][k].Exz;
          data[start+2]=f1[i][j+jstart][k].Eyz;
          data[start+3]=f1[i][j+jstart][k].Eyx;
          data[start+4]=f1[i][j+jstart][k].Ezx;
          data[start+5]=f1[i][j+jstart][k].Ezy;
          data[start+6]=f1[i][j+jstart][k].Bxy;
          data[start+7]=f1[i][j+jstart][k].Bxz;
          data[start+8]=f1[i][j+jstart][k].Byz;
          data[start+9]=f1[i][j+jstart][k].Byx;
          data[start+10]=f1[i][j+jstart][k].Bzx;
          data[start+11]=f1[i][j+jstart][k].Bzy;
          start+=12;
        }

    if(rank%2==0 && rank!=D->M-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)
          for(k=0; k<nz; k++) {
            f1[i][jend+j][k].Exy=data[start+0];
            f1[i][jend+j][k].Exz=data[start+1];
            f1[i][jend+j][k].Eyz=data[start+2];
            f1[i][jend+j][k].Eyx=data[start+3];
            f1[i][jend+j][k].Ezx=data[start+4];
            f1[i][jend+j][k].Ezy=data[start+5];
            f1[i][jend+j][k].Bxy=data[start+6];
            f1[i][jend+j][k].Bxz=data[start+7];
            f1[i][jend+j][k].Byz=data[start+8];
            f1[i][jend+j][k].Byx=data[start+9];
            f1[i][jend+j][k].Bzx=data[start+10];
            f1[i][jend+j][k].Bzy=data[start+11];
            start+=12;
          }
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevYrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(i=0; i<nx; i++)
      for(j=0; j<share; j++)
        for(k=0; k<nz; k++) { 
          data[start]=f1[i][j+jstart][k].Exy;
          data[start+1]=f1[i][j+jstart][k].Exz;
          data[start+2]=f1[i][j+jstart][k].Eyz;
          data[start+3]=f1[i][j+jstart][k].Eyx;
          data[start+4]=f1[i][j+jstart][k].Ezx;
          data[start+5]=f1[i][j+jstart][k].Ezy;
          data[start+6]=f1[i][j+jstart][k].Bxy;
          data[start+7]=f1[i][j+jstart][k].Bxz;
          data[start+8]=f1[i][j+jstart][k].Byz;
          data[start+9]=f1[i][j+jstart][k].Byx;
          data[start+10]=f1[i][j+jstart][k].Bzx;
          data[start+11]=f1[i][j+jstart][k].Bzy;
          start+=12;
        }

    if(rank%2==1 && rank!=D->M-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextYrank,D->nextYrank, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<nx; i++)
        for(j=0; j<share; j++)
          for(k=0; k<nz; k++) {
            f1[i][jend+j][k].Exy=data[start+0];
            f1[i][jend+j][k].Exz=data[start+1];
            f1[i][jend+j][k].Eyz=data[start+2];
            f1[i][jend+j][k].Eyx=data[start+3];
            f1[i][jend+j][k].Ezx=data[start+4];
            f1[i][jend+j][k].Ezy=data[start+5];
            f1[i][jend+j][k].Bxy=data[start+6];
            f1[i][jend+j][k].Bxz=data[start+7];
            f1[i][jend+j][k].Byz=data[start+8];
            f1[i][jend+j][k].Byx=data[start+9];
            f1[i][jend+j][k].Bzx=data[start+10];
            f1[i][jend+j][k].Bzy=data[start+11];
            start+=12;
          }
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevYrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_TransferPML_Yplus(Domain *D,PML ***f1,int nx,int nz,int share)
{
    int i,j,k,num;
    int istart,iend,jstart,jend,kstart,kend;
    int myrank, nTasks,rank,start;
    double *data;

    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    rank=(myrank%(D->M*D->N))%D->M;
    num=(share-1)*nx*nz*12;
    data = (double *)malloc(num*sizeof(double ));

   //Transferring even ~ odd cores 
    start=0;
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)
        for(k=0; k<nz; k++) {
          data[start+0]=f1[i][jend-j][k].Exy;
          data[start+1]=f1[i][jend-j][k].Exz;
          data[start+2]=f1[i][jend-j][k].Eyz;
          data[start+3]=f1[i][jend-j][k].Eyx;
          data[start+4]=f1[i][jend-j][k].Ezx;
          data[start+5]=f1[i][jend-j][k].Ezy;
          data[start+6]=f1[i][jend-j][k].Bxy;
          data[start+7]=f1[i][jend-j][k].Bxz;
          data[start+8]=f1[i][jend-j][k].Byz;
          data[start+9]=f1[i][jend-j][k].Byx;
          data[start+10]=f1[i][jend-j][k].Bzx;
          data[start+11]=f1[i][jend-j][k].Bzy;
          start+=12;
        }
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevYrank,D->prevYrank, MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)
          for(k=0; k<nz; k++) {
            f1[i][jstart-j][k].Exy=data[start+0];
            f1[i][jstart-j][k].Exz=data[start+1];
            f1[i][jstart-j][k].Eyz=data[start+2];
            f1[i][jstart-j][k].Eyx=data[start+3];
            f1[i][jstart-j][k].Ezx=data[start+4];
            f1[i][jstart-j][k].Ezy=data[start+5];
            f1[i][jstart-j][k].Bxy=data[start+6];
            f1[i][jstart-j][k].Bxz=data[start+7];
            f1[i][jstart-j][k].Byz=data[start+8];
            f1[i][jstart-j][k].Byx=data[start+9];
            f1[i][jstart-j][k].Bzx=data[start+10];
            f1[i][jstart-j][k].Bzy=data[start+11];
            start+=12;
          }
    }
    else if(rank%2==0 && rank!=D->M-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextYrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(i=0; i<nx; i++)
      for(j=1; j<share; j++)
        for(k=0; k<nz; k++) {
          data[start+0]=f1[i][jend-j][k].Exy;
          data[start+1]=f1[i][jend-j][k].Exz;
          data[start+2]=f1[i][jend-j][k].Eyz;
          data[start+3]=f1[i][jend-j][k].Eyx;
          data[start+4]=f1[i][jend-j][k].Ezx;
          data[start+5]=f1[i][jend-j][k].Ezy;
          data[start+6]=f1[i][jend-j][k].Bxy;
          data[start+7]=f1[i][jend-j][k].Bxz;
          data[start+8]=f1[i][jend-j][k].Byz;
          data[start+9]=f1[i][jend-j][k].Byx;
          data[start+10]=f1[i][jend-j][k].Bzx;
          data[start+11]=f1[i][jend-j][k].Bzy;
          start+=12;
        }
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevYrank,D->prevYrank,MPI_COMM_WORLD,&status);
      start=0;
      for(i=0; i<nx; i++)
        for(j=1; j<share; j++)
          for(k=0; k<nz; k++) {
            f1[i][jstart-j][k].Exy=data[start+0];
            f1[i][jstart-j][k].Exz=data[start+1];
            f1[i][jstart-j][k].Eyz=data[start+2];
            f1[i][jstart-j][k].Eyx=data[start+3];
            f1[i][jstart-j][k].Ezx=data[start+4];
            f1[i][jstart-j][k].Ezy=data[start+5];
            f1[i][jstart-j][k].Bxy=data[start+6];
            f1[i][jstart-j][k].Bxz=data[start+7];
            f1[i][jstart-j][k].Byz=data[start+8];
            f1[i][jstart-j][k].Byx=data[start+9];
            f1[i][jstart-j][k].Bzx=data[start+10];
            f1[i][jstart-j][k].Bzy=data[start+11];
            start+=12;
          }
    }
    else if(rank%2==1 && rank!=D->M-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextYrank,myrank,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}
*/


/*
double dampingU(double x,double limitR,double pmlCell,double r)
{
  double result,tmp;

  if(x<=limitR)  result=1.0;
  else      {
    tmp=r*(x-limitR)/pmlCell;
    result=1.0-tmp*tmp;
  }

  return result;
}

double dampingD(double x,double limitR,double pmlCell,double r)
{
  double result,tmp;

  if(x>=limitR)  result=1.0;
  else      {
    tmp=r*(x-limitR)/pmlCell;
    result=1.0-tmp*tmp;
  }

  return result;
}

*/
