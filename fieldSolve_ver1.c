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
void Esolve2D_Yee(Domain *D,double dF);
void Bsolve2D_NoCherenkov(Domain *D);
void solve_NDFX_C(Domain *D);
void solve_NDFX(Domain *D);
void PML_Esolve2D_Up(Domain *D);
void PML_Bsolve2D_Up(Domain *D);
void MPI_TransferPML_Xplus(Domain *D,PML ***f1,int ny,int share);
void MPI_TransferPML_Xminus(Domain *D,PML ***f1,int ny,int share);
void copyPML_E(Domain *D);
void copyPML_B(Domain *D);


void fieldSolve(Domain D,double t,int iteration,double dF)
{
  int rankX,rankY;
  LaserList *L;
  int myrank, nTasks,rank,rankM,rankN;
  void MPI_Transfer6F_Xminus();
  void MPI_Transfer6F_Xplus();
  void MPI_Transfer8F_Xminus();
  void MPI_Transfer8F_Xplus();
  void MPI_Transfer12F_Xminus();
  void MPI_Transfer12F_Xplus();

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  switch(D.fieldType) {

  case Yee :
    Esolve2D_Yee(&D,dF);
/*
    if(D.pmlUp==ON) {
      rankY=myrank%D.M;
      if(rankY==D.M-1) PML_Esolve2D_Up(&D); else ;
      MPI_Barrier(MPI_COMM_WORLD);
      if(D.L>1) {
        MPI_TransferPML_Xminus(&D,D.upml,D.pmlCellUp,3);
        MPI_TransferPML_Xplus(&D,D.upml,D.pmlCellUp,3);
      } else    ;
      if(rankY==D.M-1) copyPML_E(&D); else ;
    } else ;
*/
    MPI_Barrier(MPI_COMM_WORLD);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
      MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
//    MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
        L=L->next;
      }
    } else ;

    Bsolve2D_Yee(&D);
/*
    if(D.pmlUp==ON) {
      rankY=myrank%D.M;
      if(rankY==D.M-1) PML_Bsolve2D_Up(&D); else ;
      MPI_Barrier(MPI_COMM_WORLD);
      if(D.L>1) {
        MPI_TransferPML_Xminus(&D,D.upml,D.pmlCellUp,3);
        MPI_TransferPML_Xplus(&D,D.upml,D.pmlCellUp,3);
      } else    ;
      if(rankY==D.M-1) copyPML_B(&D); else ;
    } else ;
*/
    MPI_Barrier(MPI_COMM_WORLD);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
//    MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
    } else	;
    break ;

  case NoCherenkov :
    Esolve2D_Yee(&D,dF);
    if(D.L>1)  {
      MPI_Transfer6F_Xminus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
      MPI_Transfer6F_Xplus(&D,D.EzR,D.ErR,D.EpR,D.EzI,D.ErI,D.EpI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
//    MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
        L=L->next;
      }
    } else ;

    Bsolve2D_NoCherenkov(&D);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.BzR,D.BrR,D.BpR,D.BzI,D.BrI,D.BpI,D.BzNowR,D.BrNowR,D.BpNowR,D.BzNowI,D.BrNowI,D.BpNowI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
//    MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
    } else	;
    break ;

  case NDFX :
    solve_NDFX_C(&D);
    if(D.L>1)  {
      MPI_Transfer8F_Xminus(&D,D.PrCR,D.PlCR,D.SrCR,D.SlCR,D.PrCI,D.PlCI,D.SrCI,D.SlCI,D.nySub+5,3);
      MPI_Transfer4F_Xplus(&D,D.PrCR,D.SrCR,D.PrCI,D.SrCI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer3F_Yminus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
//    MPI_Transfer3F_Yplus(&D,D.Ex,D.Ey,D.Ez,D.nxSub+5,1,3);
    } else	;

    //load laser
    if(D.boostOn==OFF)   {
      L=D.laserList;
      while(L->next)  {
        loadLaser(&D,L,t); 
        L=L->next;
      }
    } else ;


    solve_NDFX(&D);
    if(D.L>1)  {
      MPI_Transfer12F_Xminus(&D,D.PrR,D.PlR,D.SrR,D.SlR,D.PrI,D.PlI,D.SrI,D.SlI,D.EzR,D.EzI,D.BzR,D.BzI,D.nySub+5,3);
      MPI_Transfer12F_Xplus(&D,D.PrR,D.PlR,D.SrR,D.SlR,D.PrI,D.PlI,D.SrI,D.SlI,D.EzR,D.EzI,D.BzR,D.BzI,D.nySub+5,3);
    } else	;
    if(D.M>1)  {
//    MPI_Transfer6F_Yminus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
//    MPI_Transfer6F_Yplus(&D,D.Bx,D.By,D.Bz,D.BxNow,D.ByNow,D.BzNow,D.nxSub+5,1,3);
    } else	;
    break ;
  }
}

//lala
void Bsolve2D_NoCherenkov(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double deltaZ,betaPZ,betaRZ,alphaZ,alphaP,alphaR;
  double upr,upd,upL,LdU,rr,rd,tmp,tmpr,tmpd;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr; dz=D->dz;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;

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

        tmpr=upr;
        tmpd=upd;
        tmp=tmpr*(-0.5*dtBydr/r*(D->EpR[m][i][j+1]+D->EpR[m][i][j])-dtBydr*(alphaR*(D->EpR[m][i][j+1]-D->EpR[m][i][j])+betaRZ*(D->EpR[m][i+1][j+1]-D->EpR[m][i+1][j]+D->EpR[m][i-1][j+1]-D->EpR[m][i-1][j]))-m*dtBydr/r*(alphaP*D->ErI[m][i][j]+betaPZ*(D->ErI[m][i+1][j]+D->ErI[m][i-1][j])));
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);
        tmp=tmpr*(-0.5*dtBydr/r*(D->EpI[m][i][j+1]+D->EpI[m][i][j])-dtBydr*(alphaR*(D->EpI[m][i][j+1]-D->EpI[m][i][j])+betaRZ*(D->EpI[m][i+1][j+1]-D->EpI[m][i+1][j]+D->EpI[m][i-1][j+1]-D->EpI[m][i-1][j]))+m*dtBydr/r*(alphaP*D->ErR[m][i][j]+betaPZ*(D->ErR[m][i+1][j]+D->ErR[m][i-1][j])));
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(alphaZ*(D->EpR[m][i+1][j]-D->EpR[m][i][j])+deltaZ*(D->EpR[m][i+2][j]-D->EpR[m][i-1][j]))+m*dtBydr/(r-0.5)*(alphaP*D->EzI[m][i][j]+betaPZ*(D->EzI[m][i+1][j]+D->EzI[m][i-1][j])));
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(alphaZ*(D->EpI[m][i+1][j]-D->EpI[m][i][j])+deltaZ*(D->EpI[m][i+2][j]-D->EpI[m][i-1][j]))-m*dtBydr/(r-0.5)*(alphaP*D->EzR[m][i][j]+betaPZ*(D->EzR[m][i+1][j]+D->EzR[m][i-1][j])));
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(alphaR*(D->EzR[m][i][j+1]-D->EzR[m][i][j])+betaRZ*(D->EzR[m][i+1][j+1]-D->EzR[m][i+1][j]+D->EzR[m][i-1][j+1]-D->EzR[m][i-1][j]))-dtBydz*(alphaZ*(D->ErR[m][i+1][j]-D->ErR[m][i][j])+deltaZ*(D->ErR[m][i+2][j]-D->ErR[m][i-1][j])));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpd*(dtBydr*(alphaR*(D->EzI[m][i][j+1]-D->EzI[m][i][j])+betaRZ*(D->EzI[m][i+1][j+1]-D->EzI[m][i+1][j]+D->EzI[m][i-1][j+1]-D->EzI[m][i-1][j]))-dtBydz*(alphaZ*(D->ErI[m][i+1][j]-D->ErI[m][i][j])+deltaZ*(D->ErI[m][i+2][j]-D->ErI[m][i-1][j])));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  j=jstart;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      oldBzR=D->BzR[m][i][j];
      oldBrR=D->BrR[m][i][j];
      oldBpR=D->BpR[m][i][j];
      oldBzI=D->BzI[m][i][j];
      oldBrI=D->BrI[m][i][j];
      oldBpI=D->BpI[m][i][j];
      if(m==0) {
        D->BzR[m][i][j]+=-4.0*dtBydr*D->EpR[m][i][j+1];
        D->BpR[m][i][j]=0.0;
        D->BrR[m][i][j]=-D->BrR[m][i][j];
      } else if (m==1) {
        D->BzR[m][i][j]=0.0;
        D->BzI[m][i][j]=0.0;
        D->BpR[m][i][j]+=-dtBydz*(alphaZ*(D->ErR[m][i+1][j]-D->ErR[m][i][j])+deltaZ*(D->ErR[m][i+2][j]-D->ErR[m][i-1][j]));
        D->BpI[m][i][j]+=-dtBydz*(alphaZ*(D->ErI[m][i+1][j]-D->ErI[m][i][j])+deltaZ*(D->ErI[m][i+2][j]-D->ErI[m][i-1][j]));
        D->BrR[m][i][j]=D->BrR[m][i][j+1];
        D->BrI[m][i][j]=D->BrI[m][i][j+1];
      } else ;
      D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
      D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
      D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
      D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
      D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
      D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
    }
}

void Bsolve2D_Yee(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,BzR,BrR,BpR,BzI,BrI,BpI;
  double oldBzR,oldBrR,oldBpR,oldBzI,oldBrI,oldBpI;
  double upr,upd,upL,LdU,rr,rd,tmp,tmpr,tmpd;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  numMode=D->numMode;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;

  dtBydr=D->dt/D->dr;
  dtBydz=D->dt/D->dz;
  dr=D->dr;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++)
    {
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

        tmpr=upr;
        tmpd=upd;
        tmp=tmpr*(-0.5*dtBydr/r*(D->EpR[m][i][j+1]+D->EpR[m][i][j])-dtBydr*(D->EpR[m][i][j+1]-D->EpR[m][i][j])-m*dtBydr/r*D->ErI[m][i][j]);
        D->BzR[m][i][j]=tmpd*(oldBzR+tmp);
        tmp=tmpr*(-0.5*dtBydr/r*(D->EpI[m][i][j+1]+D->EpI[m][i][j])-dtBydr*(D->EpI[m][i][j+1]-D->EpI[m][i][j])+m*dtBydr/r*D->ErR[m][i][j]);
        D->BzI[m][i][j]=tmpd*(oldBzI+tmp);

        tmp=tmpr*(dtBydz*(D->EpR[m][i+1][j]-D->EpR[m][i][j])+m*dtBydr/(r-0.5)*D->EzI[m][i][j]);
        D->BrR[m][i][j]=tmpd*(oldBrR+tmp);
        tmp=tmpr*(dtBydz*(D->EpI[m][i+1][j]-D->EpI[m][i][j])-m*dtBydr/(r-0.5)*D->EzR[m][i][j]);
        D->BrI[m][i][j]=tmpd*(oldBrI+tmp);

        tmp=tmpr*(dtBydr*(D->EzR[m][i][j+1]-D->EzR[m][i][j])-dtBydz*(D->ErR[m][i+1][j]-D->ErR[m][i][j]));
        D->BpR[m][i][j]=tmpd*(oldBpR+tmp);
        tmp=tmpd*(dtBydr*(D->EzI[m][i][j+1]-D->EzI[m][i][j])-dtBydz*(D->ErI[m][i+1][j]-D->ErI[m][i][j]));
        D->BpI[m][i][j]=tmpd*(oldBpI+tmp);

        D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
        D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
        D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
        D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
        D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
        D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
      }
    }

  j=jstart;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      oldBzR=D->BzR[m][i][j];
      oldBrR=D->BrR[m][i][j];
      oldBpR=D->BpR[m][i][j];
      oldBzI=D->BzI[m][i][j];
      oldBrI=D->BrI[m][i][j];
      oldBpI=D->BpI[m][i][j];
      if(m==0) {
        D->BzR[m][i][j]+=-4.0*dtBydr*D->EpR[m][i][j+1];
        D->BpR[m][i][j]=0.0;
        D->BrR[m][i][j]=-D->BrR[m][i][j];
      } else if (m==1) {
        D->BzR[m][i][j]=0.0;
        D->BzI[m][i][j]=0.0;
        D->BpR[m][i][j]+=-dtBydz*(D->ErR[m][i+1][j]-D->ErR[m][i][j]);
        D->BpI[m][i][j]+=-dtBydz*(D->ErI[m][i+1][j]-D->ErI[m][i][j]);
        D->BrR[m][i][j]=D->BrR[m][i][j+1];
        D->BrI[m][i][j]=D->BrI[m][i][j+1];
      } else ;
      D->BzNowR[m][i][j]=0.5*(D->BzR[m][i][j]+oldBzR);
      D->BrNowR[m][i][j]=0.5*(D->BrR[m][i][j]+oldBrR);
      D->BpNowR[m][i][j]=0.5*(D->BpR[m][i][j]+oldBpR);
      D->BzNowI[m][i][j]=0.5*(D->BzI[m][i][j]+oldBzI);
      D->BrNowI[m][i][j]=0.5*(D->BrI[m][i][j]+oldBrI);
      D->BpNowI[m][i][j]=0.5*(D->BpI[m][i][j]+oldBpI);
    }
}

//lala
void Esolve2D_Yee(Domain *D,double dF)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dt,dr;
  double EzR,ErR,EpR,EzI,ErI,EpI;
  double upr,upd,upL,LdU,rr,rd,tmp,tmpr,tmpd;
  double oldEzR,oldEzI,oldErR,oldErI,oldEpR,oldEpI;

  istart=D->istart;    iend=D->iend;
  jstart=D->jstart;    jend=D->jend;
  dt=D->dt;	dr=D->dr;
  numMode=D->numMode;

  upL=(double)(D->ny+D->minYDomain-D->pmlCellUp);
  LdU=D->pmlCellUp;
  rr=D->pmlr;    rd=D->pmld;

  dtBydr=dt/dr;
  dtBydz=dt/D->dz;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) 
    {
      for(j=jstart+1; j<jend; j++) 
      {
        r=(double)(j-jstart);
        if(D->pmlUp==ON)
          absorb_U(D,&upr,&upd,r,upL,LdU,rr,rd);
        else    ;
        oldEzR=D->EzR[m][i][j];
        oldEzI=D->EzI[m][i][j];
        oldErR=D->ErR[m][i][j];
        oldErI=D->ErI[m][i][j];
        oldEpR=D->EpR[m][i][j];
        oldEpI=D->EpI[m][i][j];

        tmpr=upr;
        tmpd=upd;
        tmp=tmpr*(0.5*dtBydr/(r-0.5)*(D->BpR[m][i][j]+D->BpR[m][i][j-1])+dtBydr*(D->BpR[m][i][j]-D->BpR[m][i][j-1])+m*dtBydr/(r-0.5)*D->BrI[m][i][j]-2.0*pi*dt*D->JzR[m][i][j]+dF*dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]));
        D->EzR[m][i][j]=tmpd*(oldEzR+tmp);
        tmp=tmpr*(0.5*dtBydr/(r-0.5)*(D->BpI[m][i][j]+D->BpI[m][i][j-1])+dtBydr*(D->BpI[m][i][j]-D->BpI[m][i][j-1])-m*dtBydr/(r-0.5)*D->BrR[m][i][j]-2.0*pi*dt*D->JzI[m][i][j]+dF*dtBydz*(D->FI[m][i+1][j]-D->FI[m][i][j]));
        D->EzI[m][i][j]=tmpd*(oldEzI+tmp);
       
        tmp=tmpr*(-dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])-m*dtBydr/r*D->BzI[m][i][j]-2.0*pi*dt*D->JrR[m][i][j]+dF*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]));
        D->ErR[m][i][j]=tmpd*(oldErR+tmp);
        tmp=tmpr*(-dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])+m*dtBydr/r*D->BzR[m][i][j]-2.0*pi*dt*D->JrI[m][i][j]+dF*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]));
        D->ErI[m][i][j]=tmpd*(oldErI+tmp);
        
        tmp=tmpr*(-dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])+dtBydz*(D->BrR[m][i][j]-D->BrR[m][i-1][j])-2.0*pi*dt*D->JpR[m][i][j]-dF*dtBydr*m/(r-0.5)*D->FI[m][i][j]);
        D->EpR[m][i][j]=tmpd*(oldEpR+tmp);
        tmp=tmpr*(-dtBydr*(D->BzI[m][i][j]-D->BzI[m][i][j-1])+dtBydz*(D->BrI[m][i][j]-D->BrI[m][i-1][j])-2.0*pi*dt*D->JpI[m][i][j]+dF*dtBydr*m/(r-0.5)*D->FR[m][i][j]);
        D->EpI[m][i][j]=tmpd*(oldEpI+tmp);
      }
    }

  j=jstart;
  r=0.0;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      if(m==0) {
        D->EzR[m][i][j]=D->EzR[m][i][j];
        D->EpR[m][i][j]=-D->EpR[m][i][j];
        D->ErR[m][i][j]-=2.0*pi*dt*D->JrR[m][i][j];
      } else if (m==1) {
        D->EzR[m][i][j]=-D->EzR[m][i][j];
        D->EzI[m][i][j]=-D->EzI[m][i][j];
        D->EpR[m][i][j]=D->EpR[m][i][j];
        D->EpI[m][i][j]=D->EpI[m][i][j];
        D->ErR[m][i][j]+=-dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])-2.0*pi*dt*D->JrR[m][i][j];
        D->ErI[m][i][j]+=-dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])-2.0*pi*dt*D->JrI[m][i][j];
      } else ;
    }

/*
  j=jstart;
  r=0.0;
  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      D->EzR[m][i][j]-=2.0*pi*dt*D->JzR[m][i][j]; //dtBydz*(D->FR[m][i+1][j]-D->FR[m][i][j]); //2.0*D->EzR[m][i][j+1]-D->EzR[m][i][j+2];
      D->EzI[m][i][j]-=2.0*pi*dt*D->JzI[m][i][j]; //dtBydz*(D->FI[m][i+1][j]-D->FI[m][i][j]); //2.0*D->EzI[m][i][j+1]-D->EzI[m][i][j+2];
      D->EpR[m][i][j]-=2.0*pi*dt*D->JpR[m][i][j]; //dtBydr/(r-0.5)*D->FI[m][i][j]; //2.0*D->EpR[m][i][j+1]-D->EpR[m][i][j+2];
      D->EpI[m][i][j]-=2.0*pi*dt*D->JpI[m][i][j]; //2.0*D->EpI[m][i][j+1]-D->EpI[m][i][j+2];
      if(m==1) {
        D->ErR[m][i][j]+=-dtBydz*(D->BpR[m][i][j]-D->BpR[m][i-1][j])-dtBydr*D->BzI[m][i][j+1]-2.0*pi*dt*D->JrR[m][i][j]+df*dtBydr*(D->FR[m][i][j+1]-D->FR[m][i][j]);
        D->ErI[m][i][j]+=-dtBydz*(D->BpI[m][i][j]-D->BpI[m][i-1][j])+dtBydr*D->BzR[m][i][j+1]-2.0*pi*dt*D->JrI[m][i][j]+df*dtBydr*(D->FI[m][i][j+1]-D->FI[m][i][j]);
      } else {
        D->ErR[m][i][j]-=2.0*pi*dt*D->JrR[m][i][j];
        D->ErI[m][i][j]-=2.0*pi*dt*D->JrI[m][i][j];
      }
    }
*/
}

void solve_NDFX_C(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,coef1,coef2,ms,ps;
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
      coef1=m*0.125/((double)(j-jstart))*dtBydr;
      ps=0.25*(0.5/((double)(j-jstart))+1.0)*dtBydr;
      ms=0.25*(0.5/((double)(j-jstart))-1.0)*dtBydr;
      coef2=m*0.5/((double)(j-jstart-0.5))*dtBydr;
      for(i=istart; i<iend; i++) {
        PrR0=D->PrR[m][i][j]; PrR1=D->PrR[m][i][j+1]; PrR2=D->PrR[m][i+1][j]; PrR3=D->PrR[m][i+1][j+1];
        PlR0=D->PlR[m][i][j]; PlR1=D->PlR[m][i][j+1]; PlR2=D->PlR[m][i+1][j]; PlR3=D->PlR[m][i+1][j+1];
        SrR0=D->SrR[m][i][j]; SrR1=D->SrR[m][i][j+1]; SrR2=D->SrR[m][i+1][j]; SrR3=D->SrR[m][i+1][j+1];
        SlR0=D->SlR[m][i][j]; SlR1=D->SlR[m][i][j+1]; SlR2=D->SlR[m][i+1][j]; SlR3=D->SlR[m][i+1][j+1];
        PrI0=D->PrI[m][i][j]; PrI1=D->PrR[m][i][j+1]; PrI2=D->PrI[m][i+1][j]; PrI3=D->PrI[m][i+1][j+1];
        PlI0=D->PlI[m][i][j]; PlI1=D->PlR[m][i][j+1]; PlI2=D->PlI[m][i+1][j]; PlI3=D->PlI[m][i+1][j+1];
        SrI0=D->SrI[m][i][j]; SrI1=D->SrR[m][i][j+1]; SrI2=D->SrI[m][i+1][j]; SrI3=D->SrI[m][i+1][j+1];
        SlI0=D->SlI[m][i][j]; SlI1=D->SlR[m][i][j+1]; SlI2=D->SlI[m][i+1][j]; SlI3=D->SlI[m][i+1][j+1];

        D->EzCR[m][i][j]+=-coef1*(SrI0+SrI1+SrI2+SrI3-SlI0-SlI1-SlI2-SlI3)+ps*(PrR1+PrR3-PlR0-PlR2)+ms*(PrR0+PrR2-PlR1-PlR3)-pi*dt*(D->JzR[m][i][j]+D->JzR[m][i+1][j]);
        D->BzCR[m][i][j]+=-coef1*(PrI0+PrI1+PrI2+PrI3+PlI0+PlI1+PlI2+PlI3)-ps*(SrR3+SrR1+SlR3+SlR1)-ms*(SrR2+SrR0+SlR2+SlR0);
        D->EzCI[m][i][j]+= coef1*(SrR0+SrR1+SrR2+SrR3-SlR0-SlR1-SlR2-SlR3)+ps*(PrI1+PrI3-PlI0-PlI2)+ms*(PrI0+PrI2-PlI1-PlI3)-pi*dt*(D->JzI[m][i][j]+D->JzI[m][i+1][j]);
        D->BzCI[m][i][j]+= coef1*(PrR0+PrR1+PrR2+PrR3+PlR0+PlR1+PlR2+PlR3)-ps*(SrI3+SrI1+SlI3+SlI1)-ms*(SrI2+SrI0+SlI2+SlI0);

        prevPrCR=nowPrCR;
        nowPrCR=D->PrCR[m][i][j];
        D->PrCR[m][i][j]=prevPrCR          -coef2*(D->BzI[m][i][j]+D->BzI[m][i][j-1])+dtBydr*(D->EzR[m][i][j]-D->EzR[m][i][j-1])-2.0*pi*dt*D->JrR[m][i][j];
        D->PlCR[m][i][j]=D->PlCR[m][i+1][j]-coef2*(D->BzI[m][i][j]+D->BzI[m][i][j-1])-dtBydr*(D->EzR[m][i][j]-D->EzR[m][i][j-1])-2.0*pi*dt*D->JrR[m][i][j];
        prevPrCI=nowPrCI;
        nowPrCI=D->PrCI[m][i][j];
        D->PrCI[m][i][j]=prevPrCI          +coef2*(D->BzR[m][i][j]+D->BzR[m][i][j-1])+dtBydr*(D->EzI[m][i][j]-D->EzI[m][i][j-1])-2.0*pi*dt*D->JrI[m][i][j];
        D->PlCI[m][i][j]=D->PlCI[m][i+1][j]+coef2*(D->BzR[m][i][j]+D->BzR[m][i][j-1])-dtBydr*(D->EzI[m][i][j]-D->EzI[m][i][j-1])-2.0*pi*dt*D->JrI[m][i][j];

        prevSrCR=nowSrCR;
        nowSrCR=D->SrCR[m][i][j];
        D->SrCR[m][i][j]=prevSrCR          -coef2*(D->EzI[m][i][j]+D->EzI[m][i][j-1])-dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])-pi*dt*(D->JpR[m][i][j]+D->JpR[m][i+1][j]);
        D->SlCR[m][i][j]=D->SlCR[m][i+1][j]+coef2*(D->EzI[m][i][j]+D->EzI[m][i][j-1])-dtBydr*(D->BzR[m][i][j]-D->BzR[m][i][j-1])-pi*dt*(D->JpR[m][i][j]+D->JpR[m][i+1][j]);
        prevSrCI=nowSrCI;
        nowSrCI=D->SrCI[m][i][j];
        D->SrCI[m][i][j]=prevSrCI          +coef2*(D->EzR[m][i][j]+D->EzR[m][i][j-1])-dtBydr*(D->BzI[m][i][j]-D->BzI[m][i][j-1])-pi*dt*(D->JpI[m][i][j]+D->JpI[m][i+1][j]);
        D->SlCI[m][i][j]=D->SlCI[m][i+1][j]-coef2*(D->EzR[m][i][j]+D->EzR[m][i][j-1])-dtBydr*(D->BzI[m][i][j]-D->BzI[m][i][j-1])-pi*dt*(D->JpI[m][i][j]+D->JpI[m][i+1][j]);
      }       //End of i
    }         //End of j

  j=jstart;
  for(m=0; m<numMode; m++) {
    for(i=istart; i<iend; i++) {
      if(m==0) {
        D->BzCR[m][i][j]+=-dtBydr*(D->SrR[m][i][j+1]+D->SlR[m][i][j+1]+D->SrR[m][i+1][j+1]+D->SlR[m][i+1][j+1]);
        D->EzCR[m][i][j]+= dtBydr*(D->PrR[m][i][j+1]+D->PrR[m][i+1][j+1]-D->PlR[m][i][j+1]-D->PlR[m][i+1][j+1])-2.0*pi*dt*D->JzR[m][i][j];
      } else {
        D->PrCR[m][i][j]=-D->PrCR[m][i][j+1];
        D->PlCR[m][i][j]=-D->PlCR[m][i][j+1];
        D->SrCR[m][i][j]=-D->SrCR[m][i][j+1];
        D->SlCR[m][i][j]=-D->SlCR[m][i][j+1];
        D->PrCI[m][i][j]=-D->PrCI[m][i][j+1];
        D->PlCI[m][i][j]=-D->PlCI[m][i][j+1];
        D->SrCI[m][i][j]=-D->SrCI[m][i][j+1];
        D->SlCI[m][i][j]=-D->SlCI[m][i][j+1];
      }
    }
  }
}

void solve_NDFX(Domain *D)
{
  int i,j,m,numMode,istart,iend,jstart,jend;  
  double dtBydr,dtBydz,r,dr,dz,dt,coef1,coef2,ms,ps;
  double nowPrR,nowSrR,prevPrR,prevSrR;
  double nowPrI,nowSrI,prevPrI,prevSrI;
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
      nowPrR=D->PrR[m][istart-1][j];
      nowSrR=D->SrR[m][istart-1][j];
      nowPrI=D->PrI[m][istart-1][j];
      nowSrI=D->SrI[m][istart-1][j];
      coef1=m*0.125/((double)(j-jstart))*dtBydr;
      ps=0.25*(0.5/((double)(j-jstart))+1.0)*dtBydr;
      ms=0.25*(0.5/((double)(j-jstart))-1.0)*dtBydr;
      coef2=m*0.5/((double)(j-jstart-0.5))*dtBydr;
      for(i=istart; i<iend; i++) {
        PrR0=D->PrCR[m][i][j]; PrR1=D->PrCR[m][i][j+1]; PrR2=D->PrCR[m][i+1][j]; PrR3=D->PrCR[m][i+1][j+1];
        PlR0=D->PlCR[m][i][j]; PlR1=D->PlCR[m][i][j+1]; PlR2=D->PlCR[m][i+1][j]; PlR3=D->PlCR[m][i+1][j+1];
        SrR0=D->SrCR[m][i][j]; SrR1=D->SrCR[m][i][j+1]; SrR2=D->SrCR[m][i+1][j]; SrR3=D->SrCR[m][i+1][j+1];
        SlR0=D->SlCR[m][i][j]; SlR1=D->SlCR[m][i][j+1]; SlR2=D->SlCR[m][i+1][j]; SlR3=D->SlCR[m][i+1][j+1];
        PrI0=D->PrCI[m][i][j]; PrI1=D->PrCR[m][i][j+1]; PrI2=D->PrCI[m][i+1][j]; PrI3=D->PrCI[m][i+1][j+1];
        PlI0=D->PlCI[m][i][j]; PlI1=D->PlCR[m][i][j+1]; PlI2=D->PlCI[m][i+1][j]; PlI3=D->PlCI[m][i+1][j+1];
        SrI0=D->SrCI[m][i][j]; SrI1=D->SrCR[m][i][j+1]; SrI2=D->SrCI[m][i+1][j]; SrI3=D->SrCI[m][i+1][j+1];
        SlI0=D->SlCI[m][i][j]; SlI1=D->SlCR[m][i][j+1]; SlI2=D->SlCI[m][i+1][j]; SlI3=D->SlCI[m][i+1][j+1];

        D->EzR[m][i][j]+=-coef1*(SrI0+SrI1+SrI2+SrI3-SlI0-SlI1-SlI2-SlI3)+ps*(PrR1+PrR3-PlR0-PlR2)+ms*(PrR0+PrR2-PlR1-PlR3)-pi*dt*(D->JzR[m][i][j]+D->JzR[m][i+1][j]);
        D->BzR[m][i][j]+=-coef1*(PrI0+PrI1+PrI2+PrI3+PlI0+PlI1+PlI2+PlI3)-ps*(SrR3+SrR1+SlR3+SlR1)-ms*(SrR2+SrR0+SlR2+SlR0);
        D->EzI[m][i][j]+= coef1*(SrR0+SrR1+SrR2+SrR3-SlR0-SlR1-SlR2-SlR3)+ps*(PrI1+PrI3-PlI0-PlI2)+ms*(PrI0+PrI2-PlI1-PlI3)-pi*dt*(D->JzI[m][i][j]+D->JzI[m][i+1][j]);
        D->BzI[m][i][j]+= coef1*(PrR0+PrR1+PrR2+PrR3+PlR0+PlR1+PlR2+PlR3)-ps*(SrI3+SrI1+SlI3+SlI1)-ms*(SrI2+SrI0+SlI2+SlI0);

        prevPrR=nowPrR;
        nowPrR=D->PrR[m][i][j];
        D->PrR[m][i][j]=prevPrR          -coef2*(D->BzCI[m][i][j]+D->BzCI[m][i][j-1])+dtBydr*(D->EzCR[m][i][j]-D->EzCR[m][i][j-1])-2.0*pi*dt*D->JrR[m][i][j];
        D->PlR[m][i][j]=D->PlR[m][i+1][j]-coef2*(D->BzCI[m][i][j]+D->BzCI[m][i][j-1])-dtBydr*(D->EzCR[m][i][j]-D->EzCR[m][i][j-1])-2.0*pi*dt*D->JrR[m][i][j];
        prevSrR=nowSrR;
        nowSrR=D->SrR[m][i][j];
        D->SrR[m][i][j]=prevSrR          -coef2*(D->EzCI[m][i][j]+D->EzCI[m][i][j-1])-dtBydr*(D->BzCR[m][i][j]-D->BzCR[m][i][j-1])-pi*dt*(D->JpR[m][i][j]+D->JpR[m][i+1][j]);
        D->SlR[m][i][j]=D->SlR[m][i+1][j]+coef2*(D->EzCI[m][i][j]+D->EzCI[m][i][j-1])-dtBydr*(D->BzCR[m][i][j]-D->BzCR[m][i][j-1])-pi*dt*(D->JpR[m][i][j]+D->JpR[m][i+1][j]);
        prevPrI=nowPrI;
        nowPrI=D->PrI[m][i][j];
        D->PrI[m][i][j]=prevPrI          +coef2*(D->BzCR[m][i][j]+D->BzCR[m][i][j-1])+dtBydr*(D->EzCI[m][i][j]-D->EzCI[m][i][j-1])-2.0*pi*dt*D->JrI[m][i][j];
        D->PlI[m][i][j]=D->PlI[m][i+1][j]+coef2*(D->BzCR[m][i][j]+D->BzCR[m][i][j-1])-dtBydr*(D->EzCI[m][i][j]-D->EzCI[m][i][j-1])-2.0*pi*dt*D->JrI[m][i][j];
        prevSrI=nowSrI;
        nowSrI=D->SrI[m][i][j];
        D->SrI[m][i][j]=prevSrI          +coef2*(D->EzCR[m][i][j]+D->EzCR[m][i][j-1])-dtBydr*(D->BzCI[m][i][j]-D->BzCI[m][i][j-1])-pi*dt*(D->JpI[m][i][j]+D->JpI[m][i+1][j]);
        D->SlI[m][i][j]=D->SlI[m][i+1][j]-coef2*(D->EzCR[m][i][j]+D->EzCR[m][i][j-1])-dtBydr*(D->BzCI[m][i][j]-D->BzCI[m][i][j-1])-pi*dt*(D->JpI[m][i][j]+D->JpI[m][i+1][j]);
      }       //End of i
    }         //End of j

  j=jstart;
  for(m=0; m<numMode; m++) {
    for(i=istart; i<iend; i++) {
      if(m==0) {
        D->BzR[m][i][j]+=-dtBydr*(D->SrCR[m][i][j+1]+D->SlCR[m][i][j+1]+D->SrCR[m][i+1][j+1]+D->SlCR[m][i+1][j+1]);
        D->EzR[m][i][j]+=dtBydr*(D->PrCR[m][i][j+1]+D->PrCR[m][i+1][j+1]-D->PlCR[m][i][j+1]-D->PlCR[m][i+1][j+1])-2.0*pi*dt*D->JzR[m][i][j];
      } else {
        D->PrR[m][i][j]=-D->PrR[m][i][j+1];
        D->PlR[m][i][j]=-D->PlR[m][i][j+1];
        D->SrR[m][i][j]=-D->SrR[m][i][j+1];
        D->SlR[m][i][j]=-D->SlR[m][i][j+1];
        D->PrI[m][i][j]=-D->PrI[m][i][j+1];
        D->PlI[m][i][j]=-D->PlI[m][i][j+1];
        D->SrI[m][i][j]=-D->SrI[m][i][j+1];
        D->SlI[m][i][j]=-D->SlI[m][i][j+1];
     } 
    }
  }

}

