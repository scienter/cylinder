#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void MPI_TransferDen_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share);
void MPI_TransferDen_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share);
void solveF_Yee(Domain *D);
void solveF_Split(Domain *D);


void solveF(Domain D)
{
  switch (D.fieldType) {
  case Yee :
    solveF_Yee(&D);
    break;
  case Split :
    solveF_Split(&D);
    break;
  }
}

void solveF_Split(Domain *D)
{
  int i,j,m,s,numMode,istart,iend,jstart,jend,minRSub,n,iter;
  double invDr,invDz,r,dr,***val;
  double upPrR,upPrI,upPlR,upPlI,upSrR,upSrI,upSlR,upSlI;
  double dnPrR,dnPrI,dnPlR,dnPlI,dnSrR,dnSrI,dnSlR,dnSlI;
  char name[100];
  FILE *out;
  LoadList *LL;

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode; minRSub=D->minYSub;
  iter=D->filterIter;
  invDr=1.0/D->dr; invDz=1.0/D->dz;
  dr=D->dr;

  //initializing density
  for(m=0; m<numMode; m++)
    for(i=0; i<iend+3; i++)
      for(j=0; j<jend+3; j++) {
        D->RhoNoPairR[m][i][j]=0.0;
        D->RhoNoPairI[m][i][j]=0.0;
      }

//  if(myrank==0) istart=istart-2; else ;
//  if(myrank==D->L-1) iend=iend-3; else ;

  LL=D->loadList; s=0;
  while(LL->next)      {
    solveCharge(D,LL,D->RhoNoPairR,D->RhoNoPairI,istart,iend,jstart,jend,s,1.0);
    LL=LL->next; s++;
  }
  if(D->L>1)  {
    MPI_TransferDen_Xplus(D,D->RhoNoPairR,D->RhoNoPairI,D->nySub+5,3);
    MPI_TransferDen_Xminus(D,D->RhoNoPairR,D->RhoNoPairI,D->nySub+5,3);
  }  else ;

//  filter(D,D->RhoNoPairR,D->RhoNoPairI);

  m=0;
  for(i=istart; i<iend; i++)
    for(j=jstart+1; j<jend; j++) {
      r=(j-jstart)*dr;
      upPrR=0.5*(D->PrR[m][i][j]  +D->PrR[m][i-1][j]);
      dnPrR=0.5*(D->PrR[m][i][j-1]+D->PrR[m][i-1][j-1]);
      upPlR=0.5*(D->PlR[m][i][j]  +D->PlR[m][i-1][j]);
      dnPlR=0.5*(D->PlR[m][i][j-1]+D->PlR[m][i-1][j-1]);
      upSrR=0.5*(D->SrR[m][i][j]  +D->SrR[m][i-1][j]);
      dnSrR=0.5*(D->SrR[m][i][j-1]+D->SrR[m][i-1][j-1]);
      upSlR=0.5*(D->SlR[m][i][j]  +D->SlR[m][i-1][j]);
      dnSlR=0.5*(D->SlR[m][i][j-1]+D->SlR[m][i-1][j-1]);

      D->FR[m][i][j]=
         0.25/r*(upPrR+dnPrR+upPlR+dnPlR)
        +0.5*invDr*(upPrR-dnPrR+upPlR-dnPlR)
        +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
        -2.0*M_PI*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
    }

  for(m=1; m<numMode; m++)
    for(i=istart; i<iend; i++)
      for(j=jstart+1; j<jend; j++) {
        r=(j-jstart)*dr;
        upPrR=0.5*(D->PrR[m][i][j]  +D->PrR[m][i-1][j]);
        dnPrR=0.5*(D->PrR[m][i][j-1]+D->PrR[m][i-1][j-1]);
        upPlR=0.5*(D->PlR[m][i][j]  +D->PlR[m][i-1][j]);
        dnPlR=0.5*(D->PlR[m][i][j-1]+D->PlR[m][i-1][j-1]);
        upSrR=0.5*(D->SrR[m][i][j]  +D->SrR[m][i-1][j]);
        dnSrR=0.5*(D->SrR[m][i][j-1]+D->SrR[m][i-1][j-1]);
        upSlR=0.5*(D->SlR[m][i][j]  +D->SlR[m][i-1][j]);
        dnSlR=0.5*(D->SlR[m][i][j-1]+D->SlR[m][i-1][j-1]);
        upPrI=0.5*(D->PrI[m][i][j]  +D->PrI[m][i-1][j]);
        dnPrI=0.5*(D->PrI[m][i][j-1]+D->PrI[m][i-1][j-1]);
        upPlI=0.5*(D->PlI[m][i][j]  +D->PlI[m][i-1][j]);
        dnPlI=0.5*(D->PlI[m][i][j-1]+D->PlI[m][i-1][j-1]);
        upSrI=0.5*(D->SrI[m][i][j]  +D->SrI[m][i-1][j]);
        dnSrI=0.5*(D->SrI[m][i][j-1]+D->SrI[m][i-1][j-1]);
        upSlI=0.5*(D->SlI[m][i][j]  +D->SlI[m][i-1][j]);
        dnSlI=0.5*(D->SlI[m][i][j-1]+D->SlI[m][i-1][j-1]);

        D->FR[m][i][j]=
           0.25/r*(upPrR+dnPrR+upPlR+dnPlR)
          +0.5*invDr*(upPrR-dnPrR+upPlR-dnPlR)
          +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
          -0.25/r*m*(upSrI+dnSrI+upSlI+dnSlI)
          -2*M_PI*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
        D->FI[m][i][j]=
           0.25/r*(upPrI+dnPrI+upPlI+dnPlI)
          +0.5*invDr*(upPrI-dnPrI+upPlI-dnPlI)
          +invDz*(D->EzI[m][i][j]-D->EzI[m][i-1][j])
          +0.25/r*m*(upSrR+dnSrR+upSlR+dnSlR)
          -2*M_PI*(D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]);
      }

  //axis
  m=0; j=jstart;
  for(i=istart; i<iend; i++) {
      upPrR=0.5*(D->PrR[m][i][j]  +D->PrR[m][i-1][j]);
      upPlR=0.5*(D->PlR[m][i][j]  +D->PlR[m][i-1][j]);
      upSrR=0.5*(D->SrR[m][i][j]  +D->SrR[m][i-1][j]);
      upSlR=0.5*(D->SlR[m][i][j]  +D->SlR[m][i-1][j]);

      D->FR[m][i][j]=
//        2.0*invDr*(upPrR+upPlR)
        0+invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
        -2.0*M_PI*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
    }
  m=1; j=jstart;
    for(i=istart; i<iend; i++) {
        upPrR=0.5*(D->PrR[m][i][j]  +D->PrR[m][i-1][j]);
        upPlR=0.5*(D->PlR[m][i][j]  +D->PlR[m][i-1][j]);
        upSrR=0.5*(D->SrR[m][i][j]  +D->SrR[m][i-1][j]);
        upSlR=0.5*(D->SlR[m][i][j]  +D->SlR[m][i-1][j]);
        upPrI=0.5*(D->PrI[m][i][j]  +D->PrI[m][i-1][j]);
        upPlI=0.5*(D->PlI[m][i][j]  +D->PlI[m][i-1][j]);
        upSrI=0.5*(D->SrI[m][i][j]  +D->SrI[m][i-1][j]);
        upSlI=0.5*(D->SlI[m][i][j]  +D->SlI[m][i-1][j]);

        D->FR[m][i][j]=
//          2.0*invDr*(0.25*(D->PrR[m][i][j+1]+D->PrR[m][i-1][j+1]+D->PlR[m][i][j+1]+D->PlR[m][i-1][j+1])-0.5*(upPrR+upPlR))
          0+invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
//          -invDr*m*(0.25*(D->SrI[m][i][j+1]+D->SrI[m][i-1][j+1]+D->SlI[m][j][j+1]+D->SlI[m][i-1][j+1])-0.5*(upSrI+upSlI))
          -2*M_PI*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
        D->FI[m][i][j]=
//          2.0*invDr*(0.25*(D->PrI[m][i][j+1]+D->PrI[m][i-1][j+1]+D->PlI[m][i][j+1]+D->PlI[m][i-1][j+1])-0.5*(upPrI+upPlI))
          0+invDz*(D->EzI[m][i][j]-D->EzI[m][i-1][j])
//          +invDr*m*(0.25*(D->SrR[m][i][j+1]+D->SrR[m][i-1][j+1]+D->SlR[m][j][j+1]+D->SlR[m][i-1][j+1])-0.5*(upSrR+upSlR))
          -2*M_PI*(D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]);
      }
  for(m=2; m<numMode; m++)
    for(i=istart; i<iend; i++) {
        upPrR=0.5*(D->PrR[m][i][j]  +D->PrR[m][i-1][j]);
        upPlR=0.5*(D->PlR[m][i][j]  +D->PlR[m][i-1][j]);
        upSrR=0.5*(D->SrR[m][i][j]  +D->SrR[m][i-1][j]);
        upSlR=0.5*(D->SlR[m][i][j]  +D->SlR[m][i-1][j]);
        upPrI=0.5*(D->PrI[m][i][j]  +D->PrI[m][i-1][j]);
        upPlI=0.5*(D->PlI[m][i][j]  +D->PlI[m][i-1][j]);
        upSrI=0.5*(D->SrI[m][i][j]  +D->SrI[m][i-1][j]);
        upSlI=0.5*(D->SlI[m][i][j]  +D->SlI[m][i-1][j]);

        D->FR[m][i][j]=
//          2.0*invDr*(upPrR+upPlR)
          0+invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
//          -invDr*m*(upSrI+upSlI)
          -2*M_PI*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
        D->FI[m][i][j]=
//          2.0*invDr*(upPrI+upPlI)
          0+invDz*(D->EzI[m][i][j]-D->EzI[m][i-1][j])
//          +invDr*m*(upSrR+upSlR)
          -2*M_PI*(D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]);
      }

  if(D->L>1)  {
    MPI_Transfer2F_Xplus(D,D->FR,D->FI,D->nySub+5,3);
    MPI_Transfer2F_Xminus(D,D->FR,D->FI,D->nySub+5,3);
  }  else ;
}


void solveF_Yee(Domain *D)
{
  int i,j,m,s,numMode,istart,iend,jstart,jend,minRSub,n,iter;
  double invDr,invDz,r,***val;
  char name[100];
  FILE *out;
  LoadList *LL;

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart; iend=D->iend;
  jstart=D->jstart; jend=D->jend;
  numMode=D->numMode; minRSub=D->minYSub;
  iter=D->filterIter;
  invDr=1.0/D->dr; invDz=1.0/D->dz;

  //initializing density
  for(m=0; m<numMode; m++)
    for(i=0; i<iend+3; i++)
      for(j=0; j<jend+3; j++) {
        D->RhoNoPairR[m][i][j]=0.0;
        D->RhoNoPairI[m][i][j]=0.0;
      }

//  if(myrank==0) istart=istart-2; else ;
//  if(myrank==D->L-1) iend=iend-3; else ;

  LL=D->loadList; s=0;
  while(LL->next)      {
    solveCharge(D,LL,D->RhoNoPairR,D->RhoNoPairI,istart,iend,jstart,jend,s,1.0);
    LL=LL->next; s++;
  }
  if(D->L>1)  {
    MPI_TransferDen_Xplus(D,D->RhoNoPairR,D->RhoNoPairI,D->nySub+5,3);
    MPI_TransferDen_Xminus(D,D->RhoNoPairR,D->RhoNoPairI,D->nySub+5,3);
  }  else ;
//  filter(D,D->RhoNoPairR,D->RhoNoPairI);

  m=0;
  for(i=istart; i<iend; i++)
    for(j=jstart+1; j<jend; j++) {
      r=(double)(j-jstart);
      D->FR[m][i][j]=
         0.5*invDr/r*(D->ErR[m][i][j]+D->ErR[m][i][j-1])
        +invDr*(D->ErR[m][i][j]-D->ErR[m][i][j-1])
        +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
        -2.0*pi*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
    }

  for(m=1; m<numMode; m++)
    for(i=istart; i<iend; i++)
      for(j=jstart+1; j<jend; j++) {
        r=(double)(j-jstart);
        D->FR[m][i][j]=
           0.5*invDr/r*(D->ErR[m][i][j]+D->ErR[m][i][j-1])
          +invDr*(D->ErR[m][i][j]-D->ErR[m][i][j-1])
          +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
          -invDr/r*m*D->EpI[m][i][j];
          -2.0*pi*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
        D->FI[m][i][j]=
           0.5*invDr/r*(D->ErI[m][i][j]+D->ErI[m][i][j-1])
          +invDr*(D->ErI[m][i][j]-D->ErI[m][i][j-1])
          +invDz*(D->EzI[m][i][j]-D->EzI[m][i-1][j])
          +invDr/r*m*D->EpR[m][i][j];
          -2.0*pi*(D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]);
      }

  //axis
  m=0; j=jstart;
  for(i=istart; i<iend; i++)  {
      D->FR[m][i][j]=
         4.0*invDr*D->ErR[m][i][j]
        +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
        -2.0*pi*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
    }
  for(m=1; m<numMode; m++)
    for(i=istart; i<iend; i++) {
      D->FR[m][i][j]=0.0;
      D->FI[m][i][j]=0.0;
    }
/*
  m=1;
    for(i=istart; i<iend; i++) {
        D->FR[m][i][j]=2.0*invDr*D->ErR[m][i][j]
          +invDr*(D->ErR[m][i][j]-D->ErR[m][i][j-1])
          +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
          -2.0*invDr*m*D->EpI[m][i][j]
          -2.0*pi*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
        D->FI[m][i][j]=2.0*invDr*D->ErI[m][i][j]
          +invDr*(D->ErI[m][i][j]-D->ErI[m][i][j-1])
          +invDz*(D->EzI[m][i][j]-D->EzI[m][i-1][j])
          +2.0*invDr*m*D->EpR[m][i][j]
          -2.0*pi*(D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]);
      }
  for(m=2; m<numMode; m++)
    for(i=istart; i<iend; i++) {
        D->FR[m][i][j]=
           4.0*invDr*(D->ErR[m][i][j])
          +invDz*(D->EzR[m][i][j]-D->EzR[m][i-1][j])
          -2.0*invDr*m*D->EpI[m][i][j]
          -2.0*pi*(D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]);
        D->FI[m][i][j]=
           4.0*invDr*(D->ErI[m][i][j])
          +invDz*(D->EzI[m][i][j]-D->EzI[m][i-1][j])
          +2.0*invDr*m*D->EpR[m][i][j]
          -2.0*pi*(D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]);
      }
*/
  if(D->L>1)  {
    MPI_Transfer2F_Xplus(D,D->FR,D->FI,D->nySub+5,3);
    MPI_Transfer2F_Xminus(D,D->FR,D->FI,D->nySub+5,3);
  }  else ;
   
}

void solveCharge(Domain *D,LoadList *LL,double ***rhoR,double ***rhoI,int istart,int iend,int jstart,int jend,int s,double coef)
{
  int i,j,n,m,index,ii,jj,i1,j1,numMode,minRSub,cnt;
  double z,x,y,r,phi,Wz[2],Wr[2],weight,factor,tmp,value,invR,alpha;
  double coss[D->numMode],sins[D->numMode],rho0;
  Particle **particle;
  particle=D->particle;
  ptclList *p;
  FILE *out;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  numMode=D->numMode;
  minRSub=D->minYSub;

  rho0=LL->density/LL->criticalDensity;

  alpha=2.0;
  for(i=istart; i<iend; i++)
    for(j=jstart+1; j<jend; j++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          weight=coef*p->weight*rho0*p->charge;
          z=p->z; x=p->x; y=p->y;
          r=sqrt(x*x+y*y);   invR=1.0/r;
          index=j-jstart;

          Wr[0]=((index+1)*(index+1)-r*r)/(2.0*index+1.0);
          Wr[1]=1.0-Wr[0];
          Wz[1]=z-(int)(z);              Wz[0]=1.0-Wz[1];

          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }

          factor=weight/(2.0*index+1.0);
          for(ii=0; ii<2; ii++)
            for(jj=0; jj<2; jj++) {
//              factor=weight/(2.0*(j+jj-jstart));
              tmp=Wr[jj]*Wz[ii]*factor;
              rhoR[0][i+ii][j+jj]+=tmp;
              for(m=1; m<numMode; m++) {
                rhoR[m][i+ii][j+jj]+=tmp*coss[m]*alpha;
                rhoI[m][i+ii][j+jj]-=tmp*sins[m]*alpha;
              }
            }

          p=p->next;
        }
      }		//End of for(i,j)

  double center=1.5;
  for(i=istart; i<iend; i++)
    for(j=jstart; j<jstart+1; j++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          weight=coef*p->weight*rho0*p->charge;
          z=p->z; x=p->x; y=p->y;
          r=sqrt(x*x+y*y);   invR=1.0/r;
          index=j-jstart;

          Wr[0]=((index+1)*(index+1)-r*r)/(2.0*index+1.0);
          Wr[1]=1.0-Wr[0];
          Wz[1]=z-(int)(z);              Wz[0]=1.0-Wz[1];

          coss[1]=x*invR; sins[1]=y*invR;
          for(m=2; m<numMode; m++) {
            coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
            sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
          }

          factor=weight/(2.0*index+1.0);
          m=0;
          for(ii=0; ii<2; ii++) {
            jj=0; 
            tmp=Wr[jj]*Wz[ii]*factor*center;
            rhoR[0][i+ii][j+jj]+=tmp;
            jj=1;
            tmp=Wr[jj]*Wz[ii]*factor;
            rhoR[0][i+ii][j+jj]+=tmp;
          }
          for(m=1; m<numMode; m++) 
            for(ii=0; ii<2; ii++) {
//              for(jj=0; jj<1; jj++) {
//                rhoR[m][i+ii][j+jj]+=tmp*coss[m]*alpha*center;
//                rhoI[m][i+ii][j+jj]-=tmp*sins[m]*alpha*center;
//              }
              for(jj=1; jj<2; jj++) {
                rhoR[m][i+ii][j+jj]+=tmp*coss[m]*alpha;
                rhoI[m][i+ii][j+jj]-=tmp*sins[m]*alpha;
              }
            }

          p=p->next;
        }
      }		//End of for(i,j)
		
}

