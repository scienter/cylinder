#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void countParticle(Domain *D,int iteration,int *reIstart,int *reIend,int *reNxSub,int *reMinXSub,int *reMaxXSub);
void saveDumpParticleHDF(Domain *D,int iteration,char *fileName);
void saveDumpFieldHDF(Domain *D,int iteration,char *fileName);
void restoreRedistHDF(Domain *D,int iteration,char *partFile,char *fieldFile);

void particle_redist(Domain *D,int iteration,External *Ext)
{
  char partFile[100],fieldFile[100];
  int reIstart,reIend,reNxSub,reMinXSub,reMaxXSub;
  int myrank,i,j;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(iteration%D->redistStep==0) {
    sprintf(partFile,"redistPart%d.h5",iteration);
    sprintf(fieldFile,"redistField%d.h5",iteration);
    saveDumpFieldHDF(D,iteration,fieldFile);
    saveDumpParticleHDF(D,iteration,partFile);
    if(myrank==0) 
      printf("iteration=%d, %s and %s is made.\n",iteration,fieldFile,partFile);
    else ;

    countParticle(D,iteration,&reIstart,&reIend,&reNxSub,&reMinXSub,&reMaxXSub);

    tmpCleanMemory(D);
    D->istart=reIstart;
    D->iend=reIend;
    D->nxSub=reNxSub;
    D->minXSub=reMinXSub;
    D->maxXSub=reMaxXSub;

//printf("after rebound,myrank=%d, minXSub=%d,nxSub=%d,minXSub=%d\n",myrank,reMinXSub,D->nxSub,D->minXSub);
    reBoundary(D,Ext);


    restoreRedistHDF(D,iteration,partFile,fieldFile);
    D->minXSub+=D->minXDomain;
    D->maxXSub+=D->minXDomain;
    if(myrank==0) { ;
      remove(partFile);
      remove(fieldFile);
      printf("redistributed particles.\n");
    } else ;
    MPI_Barrier(MPI_COMM_WORLD);

  }
}

void countParticle(Domain *D,int iteration,int *reIstart,int *reIend,int *reNxSub,int *reMinXSub,int *reMaxXSub)
{
  int i,j,n,s,index,nx,minXSub,minXDomain,targetCore,fromCore,rankX;
  int total,subCnt,remainX,divisionFlag=OFF;
  int *numPart,*share,*divCore;
  int myrank,nTasks,istart,iend,jstart,jend,nSpecies,cnt;  
  double x; 
  ptclList *p;
  char name[100];
  FILE *out; 
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
//  MPI_Comm_rank(MPI_COMM_WORLD, &nTasks);
  nTasks=D->L;

  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  nSpecies=D->nSpecies;
  nx=D->nx;   minXSub=D->minXSub;   minXDomain=D->minXDomain;

  numPart=(int *)malloc(nx*sizeof(int ));
  share=(int *)malloc(nx*sizeof(int ));
  divCore=(int *)malloc((nTasks+1)*sizeof(int ));
  for(i=0; i<nx; i++) { numPart[i]=0; share[i]=0; }
  for(i=0; i<=nTasks; i++) { divCore[i]=0; }
  

  for(i=istart; i<iend; i++) {
    index=i-istart+minXSub-minXDomain;
    for(j=jstart; j<jend; j++) 
      for(s=0; s<nSpecies; s++) {
        cnt=0;
        p=D->particle[i][j].head[s]->pt;
        while(p) {
          cnt++;
          p=p->next;
        }
        numPart[index]+=cnt;
      }
  }

  for(i=0; i<nx; i++) share[i]=numPart[i];

  targetCore=myrank%D->M;
  rankX=myrank/D->M;
  for(n=1; n<D->L; n++)  {
    fromCore=targetCore+n*D->M;
    if(myrank==fromCore)
      MPI_Send(share,nx,MPI_INT,targetCore,myrank,MPI_COMM_WORLD);
    else	;
  } 

  if(myrank==targetCore) {
    for(n=1; n<D->L; n++)  {
      fromCore=targetCore+n*D->M;
      MPI_Recv(share,nx,MPI_INT,fromCore,fromCore,MPI_COMM_WORLD,&status);
      for(i=0; i<nx; i++) numPart[i]+=share[i];    
    }
  } else	;  //End of if(myrank=targetCore)
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(numPart,nx,MPI_INT,0,MPI_COMM_WORLD);

  total=0;
  for(i=0; i<nx; i++) total+=numPart[i];    
  subCnt=total/D->L;

  cnt=0;
  index=0;
  divCore[index]=0;
  for(i=0; i<nx; i++) {
    cnt+=numPart[i];
    if(cnt>=subCnt) {
      index=index+1;
      divCore[index]=i;
//      i-=1;
      cnt-=subCnt;
    } else ;
  }
/*
if(myrank==D->L-1) {
for(i=0; i<=D->L; i++) {
  printf("divCore[%d]=%d, total=%d\n",i,divCore[i],total);
 }
}
*/
  if(subCnt>1) divisionFlag=ON; else divisionFlag=OFF;

  if(divisionFlag==ON)  {
    divCore[D->L]=nx;
    *reIstart=2;
    *reNxSub=divCore[rankX+1]-divCore[rankX];
    *reIend=*reNxSub+2;
    *reMinXSub=divCore[rankX];  
    *reMaxXSub=divCore[rankX+1];      
  }
  else {
    divCore[D->L]=0;
    *reIstart=D->istart;
    *reNxSub=D->nxSub;
    *reIend=D->iend;
    *reMinXSub=D->minXSub;  
    *reMaxXSub=D->maxXSub;  
  }

  free(share);
  free(numPart);
  free(divCore);
}

