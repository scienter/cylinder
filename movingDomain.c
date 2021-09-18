#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"

void movingDomain_Cylind(Domain *D,int iteration);
void movingPML_Left(Domain *D,PML ***lpml,int left);
void MPI_Transfer2F_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share);
void MPI_Transfer2F_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share);



void movingDomain(Domain *D,int iteration)
{
  int shiftDuration,myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  void MPI_Transfer6F_Xminus();
  void MPI_Transfer6F_Xplus();

  if(D->L>1)  {
//    MPI_Transfer6F_Xplus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->nySub+5,3);
    if(D->fieldType!=Split)
      MPI_Transfer6F_Xminus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->nySub+5,3);
    else 
      MPI_Transfer12F_Xminus(D,D->JzR,D->JrR,D->JpR,D->JzI,D->JrI,D->JpI,D->JzCR,D->JrCR,D->JpCR,D->JzCI,D->JrCI,D->JpCI,D->nySub+5,3);

    MPI_Transfer2F_Xminus(D,D->RhoPairR,D->RhoPairI,D->nySub+5,3);
//    MPI_Transfer2F_Xplus(D,D->RhoPairR,D->RhoPairI,D->nySub+5,3);
    MPI_Transfer2F_Xminus(D,D->RhoNoPairR,D->RhoNoPairI,D->nySub+5,3);
//    MPI_Transfer2F_Xplus(D,D->RhoNoPairR,D->RhoNoPairI,D->nySub+5,3);
  }  else     ;

//  if(D->pmlLeft==ON && myrank==0) movingPML_Left(D,D->lpml,D->pmlCellLeft); else ;
  movingDomain_Cylind(D,iteration);
}

void movingPML_Left(Domain *D,PML ***lpml,int left)
{
  int i,j,m,numMode,istart,jend;
  istart=D->istart;
  jend=D->jend;
  numMode=D->numMode;

  for(m=0; m<numMode; m++) {
    for(i=left-1; i>0; i--)
      for(j=0; j<jend+3; j++)  {
        lpml[m][i][j].EzR=0;//lpml[m][i-1][j].EzR;
        lpml[m][i][j].EzI=0;//lpml[m][i-1][j].EzI;
        lpml[m][i][j].ErR=0;//lpml[m][i-1][j].ErR;
        lpml[m][i][j].ErI=0;//lpml[m][i-1][j].ErI;
        lpml[m][i][j].EpzR=0;//lpml[m][i-1][j].EpzR;
        lpml[m][i][j].EpzI=0;//lpml[m][i-1][j].EpzI;
        lpml[m][i][j].EprR=0;//lpml[m][i-1][j].EprR;
        lpml[m][i][j].EprI=0;//lpml[m][i-1][j].EprI;
        lpml[m][i][j].BzR=0;//lpml[m][i-1][j].BzR;
        lpml[m][i][j].BzI=0;//lpml[m][i-1][j].BzI;
        lpml[m][i][j].BrR=0;//lpml[m][i-1][j].BrR;
        lpml[m][i][j].BrI=0;//lpml[m][i-1][j].BrI;
        lpml[m][i][j].BpzR=0;//lpml[m][i-1][j].BpzR;
        lpml[m][i][j].BpzI=0;//lpml[m][i-1][j].BpzI;
        lpml[m][i][j].BprR=0;//lpml[m][i-1][j].BprR;
        lpml[m][i][j].BprI=0;//lpml[m][i-1][j].BprI;
      } 
    i=0; 
      for(j=0; j<jend+3; j++)  {
        lpml[m][i][j].EzR=0;//D->EzR[m][istart][j];
        lpml[m][i][j].EzI=0;//D->EzI[m][istart][j];
        lpml[m][i][j].EpzR=0;//D->EpR[m][istart][j]*0.5;
        lpml[m][i][j].EprR=0;//D->EpR[m][istart][j]*0.5;
        lpml[m][i][j].EpzI=0;//D->EpI[m][istart][j]*0.5;
        lpml[m][i][j].EprI=0;//D->EpI[m][istart][j]*0.5;
        lpml[m][i][j].ErR=0;//D->ErR[m][istart][j];
        lpml[m][i][j].ErI=0;//D->ErI[m][istart][j];
        lpml[m][i][j].BzR=0;//D->BzR[m][istart][j];
        lpml[m][i][j].BzI=0;//D->BzI[m][istart][j];
        lpml[m][i][j].BpzR=0;//D->BpR[m][istart][j]*0.5;
        lpml[m][i][j].BprR=0;//D->BpR[m][istart][j]*0.5;
        lpml[m][i][j].BpzI=0;//D->BpI[m][istart][j]*0.5;
        lpml[m][i][j].BprI=0;//D->BpI[m][istart][j]*0.5;
        lpml[m][i][j].BrR=0;//D->BrR[m][istart][j];
        lpml[m][i][j].BrI=0;//D->BrI[m][istart][j];
      }   
  }

}  

void movingDomain_Cylind(Domain *D,int iteration)
{
    int i,j,m,s,istart,iend,jstart,jend,shiftDuration,numMode;
    ptclList *p;      
    Particle **particle;
    particle=D->particle;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    switch (D->fieldType) {
    case Yee :
    case NoCherenkov :
      for(m=0; m<numMode; m++)
        for(i=0; i<iend+2; i++)
          for(j=0; j<jend+3; j++)  {
            D->EzR[m][i][j]=D->EzR[m][i+1][j];
            D->ErR[m][i][j]=D->ErR[m][i+1][j];
            D->EpR[m][i][j]=D->EpR[m][i+1][j];
            D->BzR[m][i][j]=D->BzR[m][i+1][j];
            D->BrR[m][i][j]=D->BrR[m][i+1][j];
            D->BpR[m][i][j]=D->BpR[m][i+1][j];
            D->EzI[m][i][j]=D->EzI[m][i+1][j];
            D->ErI[m][i][j]=D->ErI[m][i+1][j];
            D->EpI[m][i][j]=D->EpI[m][i+1][j];
            D->BzI[m][i][j]=D->BzI[m][i+1][j];
            D->BrI[m][i][j]=D->BrI[m][i+1][j];
            D->BpI[m][i][j]=D->BpI[m][i+1][j];
            D->JzR[m][i][j]=D->JzR[m][i+1][j];
            D->JrR[m][i][j]=D->JrR[m][i+1][j];
            D->JpR[m][i][j]=D->JpR[m][i+1][j];
            D->JzI[m][i][j]=D->JzI[m][i+1][j];
            D->JrI[m][i][j]=D->JrI[m][i+1][j];
            D->JpI[m][i][j]=D->JpI[m][i+1][j];
            D->RhoNoPairR[m][i][j]=D->RhoNoPairR[m][i+1][j];
            D->RhoNoPairI[m][i][j]=D->RhoNoPairI[m][i+1][j];
            D->RhoPairR[m][i][j]=D->RhoPairR[m][i+1][j];
            D->RhoPairI[m][i][j]=D->RhoPairI[m][i+1][j];
            D->FR[m][i][j]=D->FR[m][i+1][j];
            D->FI[m][i][j]=D->FI[m][i+1][j];
            D->BzNowR[m][i][j]=D->BzNowR[m][i+1][j];
            D->BrNowR[m][i][j]=D->BrNowR[m][i+1][j];
            D->BpNowR[m][i][j]=D->BpNowR[m][i+1][j];
            D->BzNowI[m][i][j]=D->BzNowI[m][i+1][j];
            D->BrNowI[m][i][j]=D->BrNowI[m][i+1][j];
            D->BpNowI[m][i][j]=D->BpNowI[m][i+1][j];
          }   
//      if(D->pmlUp==ON) {
//        for(m=0; m<numMode; m++)
//          for(i=istart-1; i<iend+1; i++) {
//            //upml
//            for(j=0; j<D->pmlCellUp; j++)
//              D->upml[m][i][j]=D->upml[m][i+1][j];
//          }
//      } else ;
      break;
    case Split :
      for(m=0; m<numMode; m++)
        for(i=0; i<iend+2; i++)
          for(j=0; j<jend+3; j++)  {
            D->EzR[m][i][j]=D->EzR[m][i+1][j];
            D->EzI[m][i][j]=D->EzI[m][i+1][j];
            D->BzR[m][i][j]=D->BzR[m][i+1][j];
            D->BzI[m][i][j]=D->BzI[m][i+1][j];
            D->PrR[m][i][j]=D->PrR[m][i+1][j];
            D->PrI[m][i][j]=D->PrI[m][i+1][j];
            D->PlR[m][i][j]=D->PlR[m][i+1][j];
            D->PlI[m][i][j]=D->PlI[m][i+1][j];
            D->SrR[m][i][j]=D->SrR[m][i+1][j];
            D->SrI[m][i][j]=D->SrI[m][i+1][j];
            D->SlR[m][i][j]=D->SlR[m][i+1][j];
            D->SlI[m][i][j]=D->SlI[m][i+1][j];
            D->EzCR[m][i][j]=D->EzCR[m][i+1][j];
            D->EzCI[m][i][j]=D->EzCI[m][i+1][j];
            D->BzCR[m][i][j]=D->BzCR[m][i+1][j];
            D->BzCI[m][i][j]=D->BzCI[m][i+1][j];
            D->PrCR[m][i][j]=D->PrCR[m][i+1][j];
            D->PrCI[m][i][j]=D->PrCI[m][i+1][j];
            D->PlCR[m][i][j]=D->PlCR[m][i+1][j];
            D->PlCI[m][i][j]=D->PlCI[m][i+1][j];
            D->SrCR[m][i][j]=D->SrCR[m][i+1][j];
            D->SrCI[m][i][j]=D->SrCI[m][i+1][j];
            D->SlCR[m][i][j]=D->SlCR[m][i+1][j];
            D->SlCI[m][i][j]=D->SlCI[m][i+1][j];
            D->JzR[m][i][j]=D->JzR[m][i+1][j];
            D->JzI[m][i][j]=D->JzI[m][i+1][j];
            D->JrR[m][i][j]=D->JrR[m][i+1][j];
            D->JrI[m][i][j]=D->JrI[m][i+1][j];
            D->JpR[m][i][j]=D->JpR[m][i+1][j];
            D->JpI[m][i][j]=D->JpI[m][i+1][j];
            D->JzCR[m][i][j]=D->JzCR[m][i+1][j];
            D->JzCI[m][i][j]=D->JzCI[m][i+1][j];
            D->JrCR[m][i][j]=D->JrCR[m][i+1][j];
            D->JrCI[m][i][j]=D->JrCI[m][i+1][j];
            D->JpCR[m][i][j]=D->JpCR[m][i+1][j];
            D->JpCI[m][i][j]=D->JpCI[m][i+1][j];
            D->RhoNoPairR[m][i][j]=D->RhoNoPairR[m][i+1][j];
            D->RhoNoPairI[m][i][j]=D->RhoNoPairI[m][i+1][j];
            D->RhoPairR[m][i][j]=D->RhoPairR[m][i+1][j];
            D->RhoPairI[m][i][j]=D->RhoPairI[m][i+1][j];
            D->FR[m][i][j]=D->FR[m][i+1][j];
            D->FI[m][i][j]=D->FI[m][i+1][j];
          }   
//      if(D->pmlUp==ON) {
//        for(m=0; m<numMode; m++)
//          for(i=istart-1; i<iend+1; i++) {
//            //upml
//            for(j=0; j<D->pmlCellUp; j++)
//              D->upml[m][i][j]=D->upml[m][i+1][j];
//          }
//      } else ;
      break;
    }

    for(i=istart-1; i<iend; i++) {
      for(j=jstart; j<jend; j++)
        for(s=0; s<D->nSpecies; s++)      {
          p=particle[i][j].head[s]->pt;
          while(p)  {
            p->z-=1.0;  
            p->oldZ-=1.0;
            p=p->next;
          }
        }	
    }
    D->minXSub+=1;
    D->maxXSub+=1;
    D->minXDomain+=1;
    D->maxXDomain+=1;

}  
