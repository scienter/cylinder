#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "mesh.h"
#include "constants.h"
#include "math.h"

void absorb_R(Domain *D,double *rtr,double *rtd,double x,double rL,double LdR,double rr,double rd);
void absorb_L(Domain *D,double *lftr,double *lftd,double x,double lL,double LdL,double rr,double rd);
void absorb_U(Domain *D,double *upr,double *upd,double y,double upL,double LdU,double rr,double rd);
void Bsolve2D_Yee(Domain *D);
void Esolve2D_Yee(Domain *D,double dF,int iteration);
void Bsolve2D_NoCherenkov(Domain *D,int iteration);
void solve_NDFX_C(Domain *D);
void solve_NDFX(Domain *D);

void fieldSolve1(Domain D,double t,int iteration,double dF)
{
  int rankX,rankY;
  float limit;
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D.fieldType) {

  case Yee :
    Bsolve2D_Yee(&D);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    } else	;

/*
    if(D.pmlLeft==ON && iteration>D.pmlStart) {
      if(myrank==0) PML_Bsolve2D_Left(&D); else ;
      MPI_Barrier(MPI_COMM_WORLD);
      if(myrank==0) copyPML_Left_B(&D,D.pmlCellLeft); else ;
    } else ;
*/

    break ;

  case NoCherenkov :
    Bsolve2D_NoCherenkov(&D,iteration);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
//    MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        limit=((L->rU+L->rD)*2.1+L->retard)*1.0/D.divisionLambda;
        if(iteration<=limit) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    } else ;
    break ;

  case NDFX :
     ;
    break ;
  }
}

void fieldSolve2(Domain D,double t,int iteration,double dF)
{
  int rankX,rankY;
  float limit;
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D.fieldType) {
  case Yee :
    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        limit=((L->rU+L->rD)*2.1+L->retard)*1.0/D.divisionLambda;
        if(iteration<=limit) loadLaser(&D,L,t); else ;
        L=L->next;
      }
    } else ;

    Esolve2D_Yee(&D,dF,iteration);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
      MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    } else	;

    break ;

  case NoCherenkov :
    Esolve2D_Yee(&D,dF,iteration);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
      MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    } else	;
    break ;

  case NDFX :
    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
        L=L->next;
      }
    } else ;

    solve_NDFX_C(&D);
    if(D.L>1)  {
      MPI_Transfer4F_NDFX_Xminus(&D,D.PlCR,D.PlCI,D.SlCR,D.SlCI,D.nySub+5,3);
      MPI_Transfer4F_NDFX_Xplus(&D,D.PlCR,D.PlCI,D.SlCR,D.SlCI,D.nySub+5,3);
      MPI_Transfer8F_Xplus(&D,D.PrCR,D.PrCI,D.SrCR,D.SrCI,D.EzCR,D.EzCI,D.BzCR,D.BzCI,D.nySub+5,3);
      MPI_Transfer8F_Xminus(&D,D.PrCR,D.PrCI,D.SrCR,D.SrCI,D.EzCR,D.EzCI,D.BzCR,D.BzCI,D.nySub+5,3);
    } else	;

    solve_NDFX(&D);
    if(D.L>1)  {
      MPI_Transfer4F_NDFX_Xminus(&D,D.PlR,D.PlI,D.SlR,D.SlI,D.nySub+5,3);
      MPI_Transfer4F_NDFX_Xplus(&D,D.PlR,D.PlI,D.SlR,D.SlI,D.nySub+5,3);
      MPI_Transfer8F_Xminus(&D,D.PrR,D.PrI,D.SrR,D.SrI,D.EzR,D.BzR,D.EzI,D.BzI,D.nySub+5,3);
      MPI_Transfer8F_Xplus(&D,D.PrR,D.PrI,D.SrR,D.SrI,D.EzR,D.BzR,D.EzI,D.BzI,D.nySub+5,3);
    } else	;
    break ;
  }
}

void Bsolve2D_NoCherenkov(Domain *D,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double x,minZSub,dtBydr,dtBydz,r,dr,dz,dt,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double deltaZ,betaPZ,betaRZ,alphaZ,alphaP,alphaR;
  double upr,upd,lftr,lftd,upL,leftL,LdU,LdL,rr,rd,tmp,tmpr,tmpd;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  minZSub=D->minXSub;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt;

  leftL=(double)(D->minXDomain+D->pmlCellLeft);
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;    rd=D->pmld;
  lftr=lftd=upr=upd=1.0;

  tmp=sin(0.5*pi*dtBydz);
  deltaZ=0.25*(1.0-tmp*tmp/dtBydz/dtBydz);
  betaPZ=0.25;
  betaRZ=0.25;
  alphaZ=1.0-3.0*deltaZ;
  alphaP=1.0-2.0*betaPZ;
  alphaR=1.0-2.0*betaRZ;

  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      x=(i-istart)+minZSub;
      if(D->pmlLeft==ON && iteration>D->pmlStart)
          absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
      else ;
      for(j=jstart+1; j<jend; j++)  
      {
        r=(double)(j-jstart);
        if(D->pmlUp==ON)
          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
        else    ;
        oldBzR=D->BzR[m][i][j];
        oldBrR=D->BrR[m][i][j];
        oldBpR=D->BpR[m][i][j];
        oldBzI=D->BzI[m][i][j];
        oldBrI=D->BrI[m][i][j];
        oldBpI=D->BpI[m][i][j];

        tmpr=upr*lftr;
        tmpd=upd*lftd;
        tmp=tmpr*dtBydr/r*(-alphaR*((r+0.5)*D->EpR[m][i][j+1]-(r-0.5)*D->EpR[m][i][j])-betaRZ*((r+0.5)*D->EpR[m][i+1][j+1]-(r-0.5)*D->EpR[m][i+1][j]+(r+0.5)*D->EpR[m][i-1][j+1]-(r-0.5)*D->EpR[m][i-1][j])-m*(alphaP*D->ErI[m][i][j]+betaPZ*(D->ErI[m][i+1][j]+D->ErI[m][i-1][j])));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);
        tmp=tmpr*dtBydr/r*(-alphaR*((r+0.5)*D->EpI[m][i][j+1]-(r-0.5)*D->EpI[m][i][j])-betaRZ*((r+0.5)*D->EpI[m][i+1][j+1]-(r-0.5)*D->EpI[m][i+1][j]+(r+0.5)*D->EpI[m][i-1][j+1]-(r-0.5)*D->EpI[m][i-1][j])+m*(alphaP*D->ErR[m][i][j]+betaPZ*(D->ErR[m][i+1][j]+D->ErR[m][i-1][j])));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(alphaZ*(D->EpR[m][i+1][j]-D->EpR[m][i][j])+deltaZ*(D->EpR[m][i+2][j]-D->EpR[m][i-1][j]))+m*dtBydr/(r-0.5)*(alphaP*D->EzI[m][i][j]+betaPZ*(D->EzI[m][i+1][j]+D->EzI[m][i-1][j])));
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(alphaZ*(D->EpI[m][i+1][j]-D->EpI[m][i][j])+deltaZ*(D->EpI[m][i+2][j]-D->EpI[m][i-1][j]))-m*dtBydr/(r-0.5)*(alphaP*D->EzR[m][i][j]+betaPZ*(D->EzR[m][i+1][j]+D->EzR[m][i-1][j])));
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(alphaR*(D->EzR[m][i][j+1]-D->EzR[m][i][j])+betaRZ*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j]+D->EzR[m][i-1][j+1]-D->EzR[m][i-1][j]))-dtBydz*(alphaZ*(D->ErR[m][i+1][j]-D->ErR[m][i][j])+deltaZ*(D->ErR[m][i+2][j]-D->ErR[m][i-1][j])));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpr*(dtBydr*(alphaR*(D->EzI[m][i][j+1]-D->EzI[m][i][j])+betaRZ*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j]+D->EzI[m][i-1][j+1]-D->EzI[m][i-1][j]))-dtBydz*(alphaZ*(D->ErI[m][i+1][j]-D->ErI[m][i][j])+deltaZ*(D->ErI[m][i+2][j]-D->ErI[m][i-1][j])));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  lftr=lftd=upr=upd=1.0;
  j=jstart;
    for(i=istart; i<iend; i++) {
      x=(i-istart)+minZSub;
      if(D->pmlLeft==ON && iteration>D->pmlStart)
          absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
      else ;

      oldBzR=D->BzR[0][i][j];
      oldBpR=D->BpR[1][i][j];
      oldBpI=D->BpI[1][i][j];

      tmpr=upr*lftr;
      tmpd=upd*lftd;

      D->BzR[0][i][j]+=-4.0*dtBydr*D->EpR[0][i][j+1];

      tmp=tmpr*(dt*(4.0*D->EzR[1][i][j+1]-D->EzR[1][i][j+2])-dtBydz*(alphaZ*(D->ErR[1][i+1][j]-D->ErR[1][i][j])+deltaZ*(D->ErR[1][i+2][j]-D->ErR[1][i-1][j])));
      D->BpR[1][i][j]=tmpd*(oldBpR+tmp);
      tmp=tmpr*(dt*(4.0*D->EzI[1][i][j+1]-D->EzI[1][i][j+2])-dtBydz*(alphaZ*(D->ErI[1][i+1][j]-D->ErI[1][i][j])+deltaZ*(D->ErI[1][i+2][j]-D->ErI[1][i-1][j])));
      D->BpI[1][i][j]=tmpd*(oldBpI+tmp);
      D->BzNowR[0][i][j]=0.5*(D->BzR[0][i][j]+oldBzR);
      D->BpNowR[1][i][j]=0.5*(D->BpR[1][i][j]+oldBpR);
      D->BpNowI[1][i][j]=0.5*(D->BpI[1][i][j]+oldBpI);
    }
}

//lala
void Bsolve2D_Yee(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double upr,upd,upL,LdU,rr,rd,tmp,tmpr,tmpd,coef1,coef2,coef3,coef4;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;
  upr=upd=1.0;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
      for(j=jstart+1; j<jend; j++)  
      {
        r=(double)(j-jstart);
        coef1=m/r*dtBydr;
        coef2=m/(r-0.5)*dtBydr;
        coef3=(2.0*r)/(2.0*r-1.0)*dtBydr;
        coef4=(2.0-2.0*r)/(2.0*r-1.0)*dtBydr;
//        if(D->pmlUp==ON)
//          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
//        else    ;
        oldBzR=D->BzR[m][i][j];
        oldBrR=D->BrR[m][i][j];
        oldBpR=D->BpR[m][i][j];
        oldBzI=D->BzI[m][i][j];
        oldBrI=D->BrI[m][i][j];
        oldBpI=D->BpI[m][i][j];

        D->BzR[m][i][j]+=-dtBydr*m/r*D->ErI[m][i][j]-dtBydr/r*((r+0.5)*D->EpR[m][i][j]-(r-0.5)*D->EpR[m][i][j-1]);
        D->BzI[m][i][j]+=dtBydr*m/r*D->ErR[m][i][j]-dtBydr/r*((r+0.5)*D->EpI[m][i][j]-(r-0.5)*D->EpI[m][i][j-1]);

        D->BrR[m][i][j]+=dtBydz*(D->EpR[m][i][j]-D->EpR[m][i-1][j])+dtBydr*m/(r+0.5)*D->EzI[m][i][j];
        D->BrI[m][i][j]+=dtBydz*(D->EpI[m][i][j]-D->EpI[m][i-1][j])-dtBydr*m/(r+0.5)*D->EzR[m][i][j];

        D->BpR[m][i][j]+=dtBydr*(D->EzR[m][i][j]-D->EzR[m][i][j-1])-dtBydz*(D->ErR[m][i][j]-D->ErR[m][i-1][j]);
        D->BpI[m][i][j]+=dtBydr*(D->EzI[m][i][j]-D->EzI[m][i][j-1])-dtBydz*(D->ErI[m][i][j]-D->ErI[m][i-1][j]);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  j=jstart; m=1;
    for(i=istart; i<iend; i++) {
      oldBzR=D->BzR[0][i][j];
      oldBpR=D->BpR[1][i][j];
      oldBpI=D->BpI[1][i][j];

      D->BzR[0][i][j]+=-4.0*dtBydr*D->EpR[0][i][j+1];
      D->BpR[1][i][j]+=2.0*dtBydr*D->EzR[1][i][j]-dtBydz*(D->ErR[1][i][j]-D->ErR[1][i-1][j]);
      D->BpI[1][i][j]+=2.0*dtBydr*D->EzI[1][i][j]-dtBydz*(D->ErI[1][i][j]-D->ErI[1][i-1][j]);
      D->BzNowR[0][i][j]=0.5*(D->BzR[0][i][j]+oldBzR);
      D->BpNowR[1][i][j]=0.5*(D->BpR[1][i][j]+oldBpR);
      D->BpNowI[1][i][j]=0.5*(D->BpI[1][i][j]+oldBpI);
    }
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      oldBrR=D->BrR[m][i][j];
      oldBrI=D->BrI[m][i][j];
      D->BrR[m][i][j]+=dtBydz*(D->EpR[m][i][j]-D->EpR[m][i-1][j])+2.0*dtBydr*m*D->EzI[m][i][j];
      D->BrI[m][i][j]+=dtBydz*(D->EpI[m][i][j]-D->EpI[m][i-1][j])-2.0*dtBydr*m*D->EzR[m][i][j];
      D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
      D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
    }
}


void Esolve2D_Yee(Domain *D,double dF,int iteration)
{
  int i,j,m,numMode,istart,iend,jstart,jend,a;  
  double dtBydr,dtBydz,r,dt,coef1,coef2,coef3,coef4;
  double EzR,ErR,EpR,EzI,ErI,EpI,x,minZSub;
  double upr,upd,lftr,lftd,upL,leftL,LdU,LdL,rr,rd,tmp,tmpr,tmpd,alpha[2];
  double oldEzR,oldEzI,oldErR,oldErI,oldEpR,oldEpI;
  double beforeEpR,beforeEpI,beforeEzR,beforeEzI,beforeErR,beforeErI;
  double nowEpR,nowEpI,nowEzR,nowEzI,nowErR,nowErI;
  double beforeBpR,beforeBpI,beforeBzR,beforeBzI,beforeBrR,beforeBrI;
  double nowBpR,nowBpI,nowBzR,nowBzI,nowBrR,nowBrI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;
  minZSub=D->minXSub;

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  leftL=(double)(D->minXDomain+D->pmlCellLeft);
  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  LdL=D->pmlCellLeft;
  rr=D->pmlr;    rd=D->pmld;
  lftr=lftd=upr=upd=1.0;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz; dt=D->dt;
  for(m=0; m<numMode; m++) 
    for(i=istart; i<iend; i++) 
    {
//      x=(i-istart)+minZSub;
//      if(D->pmlLeft==ON && iteration>D->pmlStart)  
//        absorb_L(D,&lftr,&lftd,x,leftL,LdL,rr,rd);
//      else   ;
      for(j=jstart+1; j<jend; j++) 
      {
        r=(double)(j-jstart);
        coef1=m/(r-0.5)*dtBydr;
        coef2=m/r*dtBydr;
        coef3=(2.0*r+1)/(2.0*r)*dtBydr;
        coef4=(1-2.0*r)/(2.0*r)*dtBydr;
 
//        if(D->pmlUp==ON)
//          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
//        else    ;
        oldEzR=D->EzR[m][i][j];
        oldEzI=D->EzI[m][i][j];
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];
        oldEpR=D->EpR[m][i][j];
        oldEpI=D->EpI[m][i][j];

//        tmpr=1.0; //upr*lftr;
//        tmpd=1.0; //upd*lftd;
        D->EzR[m][i][j]+=dtBydr*m/(r+0.5)*D->BrI[m][i][j]+dtBydr/(r+0.5)*((r+1)*D->BpR[m][i][j+1]-r*D->BpR[m][i][j])-2.0*pi*dt*D->JzR[m][i][j];
        D->EzI[m][i][j]+=-dtBydr*m/(r+0.5)*D->BrR[m][i][j]+dtBydr/(r+0.5)*((r+1)*D->BpI[m][i][j+1]-r*D->BpI[m][i][j])-2.0*pi*dt*D->JzI[m][i][j];
       
        D->ErR[m][i][j]+=-dtBydz*(D->BpR[m][i+1][j]-D->BpR[m][i][j])-dtBydr*m/r*D->BzI[m][i][j]-2.0*pi*dt*D->JrR[m][i][j];
        D->ErI[m][i][j]+=-dtBydz*(D->BpI[m][i+1][j]-D->BpI[m][i][j])+dtBydr*m/r*D->BzR[m][i][j]-2.0*pi*dt*D->JrI[m][i][j];
        
        D->EpR[m][i][j]+=-dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])+dtBydz*(D->BrR[m][i+1][j]-D->BrR[m][i][j])-2.0*pi*dt*D->JpR[m][i][j];
        D->EpI[m][i][j]+=-dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])+dtBydz*(D->BrI[m][i+1][j]-D->BrI[m][i][j])-2.0*pi*dt*D->JpI[m][i][j];
      }
    }

  j=jstart; m=1;
    for(i=istart; i<iend; i++) {
//      D->ErR[m][i][j]=dtBydz*(D->BpR[m][i+1][j]-D->BpR[m][i][j])-dtBydr*(2*D->BzI[m][i][j+1]-0.5*D->BzI[m][i][j+2]);
//      D->ErI[m][i][j]=dtBydz*(D->BpI[m][i+1][j]-D->BpI[m][i][j])+dtBydr*(2*D->BzR[m][i][j+1]-0.5*D->BzR[m][i][j+2]);
      D->ErR[m][i][j]=dtBydz*(D->BpR[m][i+1][j]-D->BpR[m][i][j])-dtBydr*D->BzI[m][i][j+1];
if(i==istart+2) printf("ErR=%g, dif BpR=%g, BzI=%g\n",D->ErR[m][i][j],D->BpR[m][i+1][j]-D->BpR[m][i][j],D->BzI[m][i][j+1]);
      D->ErI[m][i][j]=dtBydz*(D->BpI[m][i+1][j]-D->BpI[m][i][j])+dtBydr*D->BzR[m][i][j+1];
    }
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      D->EzR[m][i][j]+=2.0*dtBydr*m*D->BrI[m][i][j]+2.0*dtBydr*(D->BpR[m][i][j+1]-D->BpR[m][i][j])-2.0*pi*dt*D->JzR[m][i][j];
      D->EzI[m][i][j]+=-2.0*dtBydr*m*D->BrR[m][i][j]+2.0*dtBydr*(D->BpI[m][i][j+1]-D->BpI[m][i][j])-2.0*pi*dt*D->JzI[m][i][j];
      D->EpR[m][i][j]+=-dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])+dtBydz*(D->BrR[m][i+1][j]-D->BrR[m][i][j])-2.0*pi*dt*D->JpR[m][i][j];
      D->EpI[m][i][j]+=-dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])+dtBydz*(D->BrI[m][i+1][j]-D->BrI[m][i][j])-2.0*pi*dt*D->JpI[m][i][j];
    }
}

void solve_NDFX_C(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,coefEB1,coefEB2,coefEB3,coefRL;
  double nowPrCR,nowSrCR,prevPrCR,prevSrCR;
  double nowPrCI,nowSrCI,prevPrCI,prevSrCI;
  double PrR0,PrR1,PrR2,PrR3,PlR0,PlR1,PlR2,PlR3;
  double PrI0,PrI1,PrI2,PrI3,PlI0,PlI1,PlI2,PlI3;
  double SrR0,SrR1,SrR2,SrR3,SlR0,SlR1,SlR2,SlR3;
  double SrI0,SrI1,SrI2,SrI3,SlI0,SlI1,SlI2,SlI3;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt;

  for(m=0; m<numMode; m++) 
    for(j=jstart+1; j<jend; j++) {
      nowPrCR=D->PrCR[m][istart-1][j];
      nowSrCR=D->SrCR[m][istart-1][j];
      nowPrCI=D->PrCI[m][istart-1][j];
      nowSrCI=D->SrCI[m][istart-1][j];
      coefEB1=m*0.125/((double)(j-jstart-0.5))*dtBydr;
      coefEB2=0.5*(j-jstart)/(2.0*(j-jstart)-1.0)*dtBydr;
      coefEB3=0.5*(1.0-j+jstart)/(2.0*(j-jstart)-1.0)*dtBydr;
      coefRL=0.5*m/((double)(j-jstart))*dtBydr;
      for(i=istart; i<iend; i++) {
        PrR0=D->PrR[m][i-1][j-1]; PrR1=D->PrR[m][i][j-1]; PrR2=D->PrR[m][i-1][j]; PrR3=D->PrR[m][i][j];
        PlR0=D->PlR[m][i-1][j-1]; PlR1=D->PlR[m][i][j-1]; PlR2=D->PlR[m][i-1][j]; PlR3=D->PlR[m][i][j];
        PrI0=D->PrI[m][i-1][j-1]; PrI1=D->PrI[m][i][j-1]; PrI2=D->PrI[m][i-1][j]; PrI3=D->PrI[m][i][j];
        PlI0=D->PlI[m][i-1][j-1]; PlI1=D->PlI[m][i][j-1]; PlI2=D->PlI[m][i-1][j]; PlI3=D->PlI[m][i][j];
        SrR0=D->SrR[m][i-1][j-1]; SrR1=D->SrR[m][i][j-1]; SrR2=D->SrR[m][i-1][j]; SrR3=D->SrR[m][i][j];
        SlR0=D->SlR[m][i-1][j-1]; SlR1=D->SlR[m][i][j-1]; SlR2=D->SlR[m][i-1][j]; SlR3=D->SlR[m][i][j];
        SrI0=D->SrI[m][i-1][j-1]; SrI1=D->SrI[m][i][j-1]; SrI2=D->SrI[m][i-1][j]; SrI3=D->SrI[m][i][j];
        SlI0=D->SlI[m][i-1][j-1]; SlI1=D->SlI[m][i][j-1]; SlI2=D->SlI[m][i-1][j]; SlI3=D->SlI[m][i][j];

        D->EzCR[m][i][j]+=coefEB1*(SlI0+SlI1+SlI2+SlI3-SrI0-SrI1-SrI2-SrI3)+coefEB2*(PrR3+PrR2-PlR3-PlR2)+coefEB3*(PrR1+PrR0-PlR1-PlR0)-0.5*pi*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j]+D->JzCR[m][i-1][j]+D->JzCR[m][i][j]);
        D->EzCI[m][i][j]+=-1.0*coefEB1*(SlR0+SlR1+SlR2+SlR3-SrR0-SrR1-SrR2-SrR3)+coefEB2*(PrI3+PrI2-PlI3-PlI2)+coefEB3*(PrI1+PrI0-PlI1-PlI0)-0.5*pi*dt*(D->JzI[m][i-1][j]+D->JzI[m][i][j]+D->JzCI[m][i-1][j]+D->JzCI[m][i][j]);
        D->BzCR[m][i][j]+=-1.0*coefEB1*(PrI0+PrI1+PrI2+PrI3+PlI0+PlI1+PlI2+PlI3)-coefEB2*(SrR3+SrR2+SlR3+SlR2)-coefEB3*(SrR1+SrR0+SlR1+SlR0);
        D->BzCI[m][i][j]+=coefEB1*(PrR0+PrR1+PrR2+PrR3+PlR0+PlR1+PlR2+PlR3)-coefEB2*(SrI3+SrI2+SlI3+SlI2)-coefEB3*(SrI1+SrI0+SlI1+SlI0);

        prevPrCR=nowPrCR;
        nowPrCR=D->PrCR[m][i][j];
        D->PrCR[m][i][j]=prevPrCR-coefRL*(D->BzI[m][i][j]+D->BzI[m][i][j+1])+dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])-pi*dt*(D->JrR[m][i][j]+D->JrCR[m][i][j]);
        prevPrCI=nowPrCI;
        nowPrCI=D->PrCI[m][i][j];
        D->PrCI[m][i][j]=prevPrCI+coefRL*(D->BzR[m][i][j]+D->BzR[m][i][j+1])+dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])-pi*dt*(D->JrI[m][i][j]+D->JrCI[m][i][j]);
        prevSrCR=nowSrCR;
        nowSrCR=D->SrCR[m][i][j];
        D->SrCR[m][i][j]=prevSrCR-coefRL*(D->EzI[m][i][j]+D->EzI[m][i][j+1])-dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])-0.5*pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1]+D->JpCR[m][i][j]+D->JpCR[m][i][j+1]);
        prevSrCI=nowSrCI;
        nowSrCI=D->SrCI[m][i][j];
        D->SrCI[m][i][j]=prevSrCI+coefRL*(D->EzR[m][i][j]+D->EzR[m][i][j+1])-dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])-0.5*pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1]+D->JpCI[m][i][j]+D->JpCI[m][i][j+1]);

        D->PlCR[m][i-1][j]=D->PlCR[m][i][j]-coefRL*(D->BzI[m][i][j]+D->BzI[m][i][j+1])-dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])-pi*dt*(D->JrR[m][i][j]+D->JrCR[m][i][j]);
        D->PlCI[m][i-1][j]=D->PlCI[m][i][j]+coefRL*(D->BzR[m][i][j]+D->BzR[m][i][j+1])-dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])-pi*dt*(D->JrI[m][i][j]+D->JrCI[m][i][j]);
        D->SlCR[m][i-1][j]=D->SlCR[m][i][j]+coefRL*(D->EzI[m][i][j]+D->EzI[m][i][j+1])-dtBydr*(D->BzR[m][i][j+1]-D->BzR[m][i][j])-0.5*pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1]+D->JpCR[m][i][j]+D->JpCR[m][i][j+1]);
        D->SlCI[m][i-1][j]=D->SlCI[m][i][j]-coefRL*(D->EzR[m][i][j]+D->EzR[m][i][j+1])-dtBydr*(D->BzI[m][i][j+1]-D->BzI[m][i][j])-0.5*pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1]+D->JpCI[m][i][j]+D->JpCI[m][i][j+1]);
      }       //End of i
    }         //End of j
  
  j=jstart; m=1;
  nowPrCR=D->PrCR[m][istart-1][j];
  nowSrCR=D->SrCR[m][istart-1][j];
  nowPrCI=D->PrCI[m][istart-1][j];
  nowSrCI=D->SrCI[m][istart-1][j];
  for(i=istart; i<iend; i++) {
    prevPrCR=nowPrCR;
    nowPrCR=D->PrCR[m][i][j];
    D->PrCR[m][i][j]=prevPrCR-2.0*dtBydr*D->BzI[m][i][j+1]+2.0*dtBydr*D->EzR[m][i][j+1]-pi*dt*(D->JrR[m][i][j]+D->JrCR[m][i][j]);
    prevPrCI=nowPrCI;
    nowPrCI=D->PrCI[m][i][j];
    D->PrCI[m][i][j]=prevPrCI+2.0*dtBydr*D->BzR[m][i][j+1]+2.0*dtBydr*D->EzI[m][i][j+1]-pi*dt*(D->JrI[m][i][j]+D->JrCI[m][i][j]);
    prevSrCR=nowSrCR;
    nowSrCR=D->SrCR[m][i][j];
    D->SrCR[m][i][j]=prevSrCR-2.0*dtBydr*D->EzI[m][i][j+1]-2.0*dtBydr*D->BzR[m][i][j+1]-pi*dt*(1.5*D->JpR[m][i][j+1]-0.5*D->JpR[m][i][j+2]+1.5*D->JpCR[m][i][j+1]-0.5*D->JpCR[m][i][j+2]);
    prevSrCI=nowSrCI;
    nowSrCI=D->SrCI[m][i][j];
    D->SrCI[m][i][j]=prevSrCI+2.0*dtBydr*D->EzR[m][i][j+1]-2.0*dtBydr*D->BzI[m][i][j+1]-0.5*pi*dt*(D->JpI[m][i][j+1]+D->JpCI[m][i][j+1]);

    D->PlCR[m][i-1][j]=D->PlCR[m][i][j]-2.0*dtBydr*D->BzI[m][i][j+1]-2.0*dtBydr*D->EzR[m][i][j+1]-pi*dt*(D->JrR[m][i][j]+D->JrCR[m][i][j]);
    D->PlCI[m][i-1][j]=D->PlCI[m][i][j]+2.0*dtBydr*D->BzR[m][i][j+1]-2.0*dtBydr*D->EzI[m][i][j+1]-pi*dt*(D->JrI[m][i][j]+D->JrCI[m][i][j]);
    D->SlCR[m][i-1][j]=D->SlCR[m][i][j]+2.0*dtBydr*D->EzI[m][i][j+1]-2.0*dtBydr*D->BzR[m][i][j+1]-0.5*pi*dt*(D->JpR[m][i][j+1]+D->JpCR[m][i][j+1]);
    D->SlCI[m][i-1][j]=D->SlCI[m][i][j]-2.0*dtBydr*D->EzR[m][i][j+1]-2.0*dtBydr*D->BzI[m][i][j+1]-0.5*pi*dt*(D->JpI[m][i][j+1]+D->JpCI[m][i][j+1]);
  }

}

//lala
void solve_NDFX(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,coefEB1,coefEB2,coefEB3,coefRL,coefF;
  double nowPrR,nowSrR,prevPrR,prevSrR;
  double nowPrI,nowSrI,prevPrI,prevSrI;
  double PrR0,PrR1,PrR2,PrR3,PlR0,PlR1,PlR2,PlR3;
  double PrI0,PrI1,PrI2,PrI3,PlI0,PlI1,PlI2,PlI3;
  double SrR0,SrR1,SrR2,SrR3,SlR0,SlR1,SlR2,SlR3;
  double SrI0,SrI1,SrI2,SrI3,SlI0,SlI1,SlI2,SlI3;
  double FR0,FR1,FR2,FR3,FI0,FI1,FI2,FI3,dF;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz; dt=D->dt; dF=D->dF;

  for(m=0; m<numMode; m++) 
    for(j=jstart+1; j<jend; j++) {
      r=(double)(j-jstart);
      nowPrR=D->PrR[m][istart-1][j];
      nowSrR=D->SrR[m][istart-1][j];
      nowPrI=D->PrI[m][istart-1][j];
      nowSrI=D->SrI[m][istart-1][j];
      coefEB1=m*0.125/(r-0.5)*dtBydr;
      coefEB2=0.5*(j-jstart)/(2.0*r-1.0)*dtBydr;
      coefEB3=0.5*(1.0-r)/(2.0*r-1.0)*dtBydr;
      coefRL=0.5*m/r*dtBydr;
      coefF=dF*m/(r-0.5)*dtBydr;
      for(i=istart; i<iend; i++) {
        PrR0=D->PrCR[m][i-1][j-1]; PrR1=D->PrCR[m][i][j-1]; PrR2=D->PrCR[m][i-1][j]; PrR3=D->PrCR[m][i][j];
        PlR0=D->PlCR[m][i-1][j-1]; PlR1=D->PlCR[m][i][j-1]; PlR2=D->PlCR[m][i-1][j]; PlR3=D->PlCR[m][i][j];
        PrI0=D->PrCI[m][i-1][j-1]; PrI1=D->PrCI[m][i][j-1]; PrI2=D->PrCI[m][i-1][j]; PrI3=D->PrCI[m][i][j];
        PlI0=D->PlCI[m][i-1][j-1]; PlI1=D->PlCI[m][i][j-1]; PlI2=D->PlCI[m][i-1][j]; PlI3=D->PlCI[m][i][j];
        SrR0=D->SrCR[m][i-1][j-1]; SrR1=D->SrCR[m][i][j-1]; SrR2=D->SrCR[m][i-1][j]; SrR3=D->SrCR[m][i][j];
        SlR0=D->SlCR[m][i-1][j-1]; SlR1=D->SlCR[m][i][j-1]; SlR2=D->SlCR[m][i-1][j]; SlR3=D->SlCR[m][i][j];
        SrI0=D->SrCI[m][i-1][j-1]; SrI1=D->SrCI[m][i][j-1]; SrI2=D->SrCI[m][i-1][j]; SrI3=D->SrCI[m][i][j];
        SlI0=D->SlCI[m][i-1][j-1]; SlI1=D->SlCI[m][i][j-1]; SlI2=D->SlCI[m][i-1][j]; SlI3=D->SlCI[m][i][j];
        FR0=D->FR[m][i][j]; FR1=D->FR[m][i+1][j]; FR2=D->FR[m][i][j+1]; FR3=D->FR[m][i+1][j+1];
        FI0=D->FI[m][i][j]; FI1=D->FI[m][i+1][j]; FI2=D->FI[m][i][j+1]; FI3=D->FI[m][i+1][j+1];

        D->EzR[m][i][j]+=coefEB1*(SlI0+SlI1+SlI2+SlI3-SrI0-SrI1-SrI2-SrI3)+coefEB2*(PrR3+PrR2-PlR3-PlR2)+coefEB3*(PrR1+PrR0-PlR1-PlR0)-pi*dt*(D->JzR[m][i-1][j]+D->JzR[m][i][j])+0.5*dF*dtBydz*(FR1-D->FR[m][i-1][j]);
        D->EzI[m][i][j]+=-coefEB1*(SlR0+SlR1+SlR2+SlR3-SrR0-SrR1-SrR2-SrR3)+coefEB2*(PrI3+PrI2-PlI3-PlI2)+coefEB3*(PrI1+PrI0-PlI1-PlI0)-pi*dt*(D->JzI[m][i-1][j]+D->JzI[m][i][j])+0.5*dF*dtBydz*(FI1-D->FI[m][i-1][j]);
        D->BzR[m][i][j]+=-coefEB1*(PrI0+PrI1+PrI2+PrI3+PlI0+PlI1+PlI2+PlI3)-coefEB2*(SrR3+SrR2+SlR3+SlR2)-coefEB3*(SrR1+SrR0+SlR1+SlR0);
        D->BzI[m][i][j]+=coefEB1*(PrR0+PrR1+PrR2+PrR3+PlR0+PlR1+PlR2+PlR3)-coefEB2*(SrI3+SrI2+SlI3+SlI2)-coefEB3*(SrI1+SrI0+SlI1+SlI0);

        prevPrR=nowPrR;
        nowPrR=D->PrR[m][i][j];
        D->PrR[m][i][j]=prevPrR-coefRL*(D->BzCI[m][i][j]+D->BzCI[m][i][j+1])+dtBydr*(D->EzCR[m][i][j+1]-D->EzCR[m][i][j])-2.0*pi*dt*D->JrR[m][i][j]+0.5*dF*dtBydr*(FR3+FR2-FR1-FR0);
        prevPrI=nowPrI;
        nowPrI=D->PrI[m][i][j];
        D->PrI[m][i][j]=prevPrI+coefRL*(D->BzCR[m][i][j]+D->BzCR[m][i][j+1])+dtBydr*(D->EzCI[m][i][j+1]-D->EzCI[m][i][j])-2.0*pi*dt*D->JrI[m][i][j]+0.5*dF*dtBydr*(FI3+FI2-FI1-FI0);
        prevSrR=nowSrR;
        nowSrR=D->SrR[m][i][j];
        D->SrR[m][i][j]=prevSrR-coefRL*(D->EzCI[m][i][j]+D->EzCI[m][i][j+1])-dtBydr*(D->BzCR[m][i][j+1]-D->BzCR[m][i][j])-pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1])-coefF*0.25*(FI0+FI1+FI2+FI3);
        prevSrI=nowSrI;
        nowSrI=D->SrI[m][i][j];
        D->SrI[m][i][j]=prevSrI+coefRL*(D->EzCR[m][i][j]+D->EzCR[m][i][j+1])-dtBydr*(D->BzCI[m][i][j+1]-D->BzCI[m][i][j])-pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1])+coefF*0.25*(FR0+FR1+FR2+FR3);

        D->PlR[m][i-1][j]=D->PlR[m][i][j]-coefRL*(D->BzCI[m][i][j]+D->BzCI[m][i][j+1])-dtBydr*(D->EzCR[m][i][j+1]-D->EzCR[m][i][j])-2.0*pi*dt*D->JrR[m][i][j]+0.5*dF*dtBydr*(FR3+FR2-FR1-FR0);
if(isnan(D->PlR[m][i-1][j])) printf("PlR=%g,i-1=%d,j=%d\n",D->PlR[m][i-1][j],i-1,j);
        D->PlI[m][i-1][j]=D->PlI[m][i][j]+coefRL*(D->BzCR[m][i][j]+D->BzCR[m][i][j+1])-dtBydr*(D->EzCI[m][i][j+1]-D->EzCI[m][i][j])-2.0*pi*dt*D->JrI[m][i][j]+0.5*dF*dtBydr*(FI3+FI2-FI1-FI0);
        D->SlR[m][i-1][j]=D->SlR[m][i][j]+coefRL*(D->EzCI[m][i][j]+D->EzCI[m][i][j+1])-dtBydr*(D->BzCR[m][i][j+1]-D->BzCR[m][i][j])-pi*dt*(D->JpR[m][i][j]+D->JpR[m][i][j+1])-coefF*0.25*(FI0+FI1+FI2+FI3);
        D->SlI[m][i-1][j]=D->SlI[m][i][j]-coefRL*(D->EzCR[m][i][j]+D->EzCR[m][i][j+1])-dtBydr*(D->BzCI[m][i][j+1]-D->BzCI[m][i][j])-pi*dt*(D->JpI[m][i][j]+D->JpI[m][i][j+1])+coefF*0.25*(FR0+FR1+FR2+FR3);
      }       //End of i
    }         //End of j
  
  j=jstart; m=1;
  nowPrR=D->PrR[m][istart-1][j];
  nowSrR=D->SrR[m][istart-1][j];
  nowPrI=D->PrI[m][istart-1][j];
  nowSrI=D->SrI[m][istart-1][j];
  for(i=istart; i<iend; i++) {
    coefF=dF*m*dtBydr;
    FR0=D->FR[m][i][j]; FR1=D->FR[m][i+1][j]; FR2=D->FR[m][i][j+1]; FR3=D->FR[m][i+1][j+1];
    FI0=D->FI[m][i][j]; FI1=D->FI[m][i+1][j]; FI2=D->FI[m][i][j+1]; FI3=D->FI[m][i+1][j+1];

    prevPrR=nowPrR;
    nowPrR=D->PrR[m][i][j];
    D->PrR[m][i][j]=prevPrR-2.0*dtBydr*D->BzCI[m][i][j+1]+2.0*dtBydr*D->EzCR[m][i][j+1]-2.0*pi*dt*D->JrR[m][i][j]+dF*dtBydr*(FR3+FR2);
    prevPrI=nowPrI;
    nowPrI=D->PrI[m][i][j];
    D->PrI[m][i][j]=prevPrI+2.0*dtBydr*D->BzCR[m][i][j+1]+2.0*dtBydr*D->EzCI[m][i][j+1]-2.0*pi*dt*D->JrI[m][i][j]+dF*dtBydr*(FI3+FI2);
    prevSrR=nowSrR;
    nowSrR=D->SrR[m][i][j];
    D->SrR[m][i][j]=prevSrR-2.0*dtBydr*D->EzCI[m][i][j+1]-2.0*dtBydr*D->BzCR[m][i][j+1]-coefF*0.5*(FI2+FI3)-pi*dt*(D->JpR[m][i][j+1]);
    prevSrI=nowSrI;
    nowSrI=D->SrI[m][i][j];
    D->SrI[m][i][j]=prevSrI+2.0*dtBydr*D->EzCR[m][i][j+1]-2.0*dtBydr*D->BzCI[m][i][j+1]+coefF*0.5*(FR2+FR3)-pi*dt*(D->JpI[m][i][j+1]);

    D->PlR[m][i-1][j]=D->PlR[m][i][j]-2.0*dtBydr*D->BzCI[m][i][j+1]-2.0*dtBydr*D->EzCR[m][i][j+1]-2.0*pi*dt*D->JrR[m][i][j]+dF*dtBydr*(FR3+FR2);
    D->PlI[m][i-1][j]=D->PlI[m][i][j]+2.0*dtBydr*D->BzCR[m][i][j+1]-2.0*dtBydr*D->EzCI[m][i][j+1]-2.0*pi*dt*D->JrI[m][i][j]+dF*dtBydr*(FI3+FI2);
    D->SlR[m][i-1][j]=D->SlR[m][i][j]+2.0*dtBydr*D->EzCI[m][i][j+1]-2.0*dtBydr*D->BzCR[m][i][j+1]-coefF*0.5*(FI2+FI3)-pi*dt*(D->JpR[m][i][j+1]);
    D->SlI[m][i-1][j]=D->SlI[m][i][j]-2.0*dtBydr*D->EzCR[m][i][j+1]-2.0*dtBydr*D->BzCI[m][i][j+1]+coefF*0.5*(FR2+FR3)-pi*dt*(D->JpI[m][i][j+1]);
  }
}




void filter(Domain *D)
{
  int a,i,j,m,istart,iend,jstart,jend,numMode;
  double alpha[2];
  double nowEpR,nowEzR,nowErR,nowBpR,nowBzR,nowBrR;
  double nowEpI,nowEzI,nowErI,nowBpI,nowBzI,nowBrI;
  double beforeEpR,beforeEzR,beforeErR,beforeBpR,beforeBzR,beforeBrR;
  double beforeEpI,beforeEzI,beforeErI,beforeBpI,beforeBzI,beforeBrI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  alpha[0]=0.5;
  alpha[1]=1.5;
  j=jstart;
  for(a=0; a<2; a++) 
    for(m=0; m<numMode; m++)
    {
      nowEpR=D->EpR[m][0][j];  nowEpI=D->EpI[m][0][j];
      nowEzR=D->EzR[m][0][j];  nowEzI=D->EzI[m][0][j];
      nowErR=D->ErR[m][0][j];  nowErI=D->ErI[m][0][j];
      nowBpR=D->BpR[m][0][j];  nowBpI=D->BpI[m][0][j];
      nowBzR=D->BzR[m][0][j];  nowBzI=D->BzI[m][0][j];
      nowBrR=D->BrR[m][0][j];  nowBrI=D->BrI[m][0][j];
      for(i=istart-1; i<iend+2; i++) {
        beforeEpR=nowEpR;      beforeEpI=nowEpI;
        beforeEzR=nowEzR;      beforeEzI=nowEzI;
        beforeErR=nowErR;      beforeErI=nowErI;
        beforeBpR=nowBpR;      beforeBpI=nowBpI;
        beforeBzR=nowBzR;      beforeBzI=nowBzI;
        beforeBrR=nowBrR;      beforeBrI=nowBrI;
        nowEpR=D->EpR[m][i][j];    nowEpI=D->EpI[m][i][j];
        nowEzR=D->EzR[m][i][j];    nowEzI=D->EzI[m][i][j];
        nowErR=D->ErR[m][i][j];    nowErI=D->ErI[m][i][j];

        D->EpR[m][i][j+1]=alpha[a]*nowEpR+(1.0-alpha[a])*(beforeEpR+D->EpR[m][i+1][j+1])*0.5;
        D->EpI[m][i][j+1]=alpha[a]*nowEpI+(1.0-alpha[a])*(beforeEpI+D->EpI[m][i+1][j+1])*0.5;
        D->EzR[m][i][j+1]=alpha[a]*nowEzR+(1.0-alpha[a])*(beforeEzR+D->EzR[m][i+1][j+1])*0.5;
        D->EzI[m][i][j+1]=alpha[a]*nowEzI+(1.0-alpha[a])*(beforeEzI+D->EzI[m][i+1][j+1])*0.5;
        D->ErR[m][i][j]=alpha[a]*nowErR+(1.0-alpha[a])*(beforeErR+D->ErR[m][i+1][j])*0.5;
        D->ErI[m][i][j]=alpha[a]*nowErI+(1.0-alpha[a])*(beforeErI+D->ErI[m][i+1][j])*0.5;
        D->BpR[m][i][j]=alpha[a]*nowBpR+(1.0-alpha[a])*(beforeBpR+D->BpR[m][i+1][j])*0.5;
        D->BpI[m][i][j]=alpha[a]*nowBpI+(1.0-alpha[a])*(beforeBpI+D->BpI[m][i+1][j])*0.5;
        D->BzR[m][i][j]=alpha[a]*nowBzR+(1.0-alpha[a])*(beforeBzR+D->BzR[m][i+1][j])*0.5;
        D->BzI[m][i][j]=alpha[a]*nowBzI+(1.0-alpha[a])*(beforeBzI+D->BzI[m][i+1][j])*0.5;
        D->BrR[m][i][j+1]=alpha[a]*nowBrR+(1.0-alpha[a])*(beforeBrR+D->BrR[m][i+1][j+1])*0.5;
        D->BrI[m][i][j+1]=alpha[a]*nowBrI+(1.0-alpha[a])*(beforeBrI+D->BrI[m][i+1][j+1])*0.5;
      }
    }
}
