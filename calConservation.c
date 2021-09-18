#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void calConservation_0th(Domain *D);
void compensate_Jz(Domain *D);

void calConservation(Domain D,int iteration)
{
  int nSpecies;

  nSpecies=D.nSpecies;

  if(D.L>1)  {
    MPI_TransferDen_Xplus(&D,D.CnR,D.CnI,D.nySub+5,3);
    MPI_TransferDen_Xminus(&D,D.CnR,D.CnI,D.nySub+5,3);
  }  else	;
  if(D.L>1)  {
    MPI_Transfer2F_Xplus(&D,D.CnR,D.CnI,D.nySub+5,3);
    MPI_Transfer2F_Xminus(&D,D.CnR,D.CnI,D.nySub+5,3);
  }  else	;

  calConservation_0th(&D);
  if(D.L>1)  {
    MPI_Transfer2F_Xplus(&D,D.CnR,D.CnI,D.nySub+5,3);
    MPI_Transfer2F_Xminus(&D,D.CnR,D.CnI,D.nySub+5,3);
  }  else	;
  compensate_Jz(&D);

}

void calConservation_0th(Domain *D)
{
    int i,j;
    int istart,iend,jstart,jend,minRSub;
    double R,invDt,invDz,invDr;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    minRSub=D->minYSub;

    invDt=1.0/D->dt; invDz=1.0/D->dz; invDr=1.0/D->dr;

    int myrank,nTasks;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    if(myrank==0) { istart=D->istart+2; } else ;
//    if(myrank==D->L-1) { iend=D->iend-3; } else ;

    for(i=istart; i<iend; i++)
      for(j=jstart+2; j<jend; j++)
      {
        R=j-jstart;
        D->CnR[0][i][j]+=((R+0.5)*D->JrR[0][i][j]-(R-0.5)*D->JrR[0][i][j-1])/R*invDr+(D->JzR[0][i][j]-D->JzR[0][i-1][j])*invDz;
      }

    // for Axis
    for(i=istart; i<iend; i++)
      for(j=jstart; j<jstart+2; j++)
      {
        if(j==jstart) {
          D->CnR[0][i][j]+=4.0*D->JrR[0][i][j]*invDr+(D->JzR[0][i][j]-D->JzR[0][i-1][j])*invDz;
        } else {
          R=j-jstart;
          D->CnR[0][i][j]+=((R+0.5)*D->JrR[0][i][j]-(R-0.5)*D->JrR[0][i][j-1])/R*invDr+(D->JzR[0][i][j]-D->JzR[0][i-1][j])*invDz;
        }
      }      //End of for(i,j)

}

void compensate_Jz(Domain *D)
{
    int i,j;
    int istart,iend,jstart,jend,minRSub;
    double dz,newJ,nowJ;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    minRSub=D->minYSub;

    dz=1.0/D->dz; 

    int myrank,nTasks;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//    if(myrank==0) { istart=D->istart+2; } else ;
//    if(myrank==D->L-1) { iend=D->iend-3; } else ;

    for(j=jstart; j<jend; j++) {
      nowJ=D->JzR[0][iend-1][j];
      for(i=iend-2; i>=istart; i--)
      {
        newJ=nowJ;
        nowJ=D->JzR[0][i][j];       
        D->JzR[0][i][j]=D->JzR[0][i+1][j]-newJ+nowJ-dz*D->CnR[0][i+1][j];
      }
    }
}

