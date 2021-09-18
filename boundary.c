#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>

PML ***PMLmemoryAsign(int nx, int ny, int nz);
double ***memoryAsign(int nx, int ny, int nz);
double ***memoryAsignJ(int nx, int ny, int nz);

void boundary(Domain *D,External *Ext)
{
     FILE *in;
     char name[100];
     double x,y,yy,tmp,dR,minR,*intensity,*phase,*dataX;
     double intensity1,intensity2,intensity3;
     double phase1,phase2,phase3,aa,bb,cc;
     int i,j,s,rank,rankX,rankY,trackStart,cnt,index;
     int remainX,remainY,subX,subY,tmpX,tmpY;
     int nx,ny,nxSub,nySub,numdataUp,numdataBt,numberData;
     int startj,nxSub1D,nySub2D;
     int minX,maxX,minY,maxY;
     int myrank, nTasks,a,add;
     LaserList *L;
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     if(D->nx%D->resolX!=0 || D->ny%D->resolY!=0)  {
       printf("check nx, ny, nz. The remain of resol value is not zero.\n");
       exit(0);
     }  else ;
     nx=D->nx/D->resolX;
     ny=D->ny/D->resolY;

     D->nxSub=nx/D->L;
     subX=D->nxSub;
     remainX=nx%D->L;
     minX=maxX=0;

     D->nySub=ny/D->M;
     subY=D->nySub;
     remainY=ny%D->M;
     minY=maxY=0;

     minX=maxX=0;
     for(rankX=0; rankX<D->L; rankX++)
     {
       if(rankX<remainX)   tmpX=subX+1;
       else                tmpX=subX;
       minX=maxX;
       maxX=minX+tmpX;

       minY=maxY=0;
       for(rankY=0; rankY<D->M; rankY++)
       {
         if(rankY<remainY)   tmpY=subY+1;
         else                tmpY=subY;
         minY=maxY;
         maxY=minY+tmpY;

         rank=rankY+rankX*D->M;
         if(myrank==rank)
         {
            D->minXSub=minX*D->resolX;
            D->maxXSub=maxX*D->resolX;
            D->nxSub=tmpX*D->resolX;
            D->minYSub=minY*D->resolY;
            D->maxYSub=maxY*D->resolY;
            D->nySub=tmpY*D->resolY;
         }
       }
     }
     D->minXSub+=D->minXDomain;    
     D->minYSub+=D->minYDomain;    
     D->maxXSub+=D->minXDomain;    
     D->maxYSub+=D->minYDomain;    

     for(rank=0; rank<nTasks; rank++)
     {
       if(myrank==rank)
         printf("rank=%d, minXSub=%d,maxSub=%d,minYSub=%d,maxYSub=%d\n",myrank,D->minXSub,D->maxXSub,D->minYSub,D->maxYSub);
       else; 
     }

     //defining next and prev domain
     rankX=myrank/D->M;
     rankY=myrank%D->M;

     D->nextXrank=rankY+(rankX+1)*D->M;
     D->prevXrank=rankY+(rankX-1)*D->M;
     if(rankX==D->L-1) 
       D->nextXrank=rankY;
     else if(rankX==0) 
       D->prevXrank=rankY+(D->L-1)*D->M;

     D->nextYrank=(rankY+1)+rankX*D->M;
     D->prevYrank=(rankY-1)+rankX*D->M;
     if(rankY==D->M-1) 
       D->nextYrank=rankX*D->M;
     else if(rankY==0) 
       D->prevYrank=(D->M-1)+rankX*D->M;

     D->istart=2;     D->iend=D->nxSub+2;
     D->jstart=2;     D->jend=D->nySub+2;

     // Field setting
     nxSub1D=D->nxSub+5;
     nySub2D=D->nySub+5;

     //Density memory
     D->RhoNoPairR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->RhoNoPairI=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->RhoPairR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->RhoPairI=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->FR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->FI=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->CnR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->CnI=memoryAsign(D->numMode,nxSub1D,nySub2D);

     switch(D->fieldType)  {
     case Yee :
     case NoCherenkov :
       D->EzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->ErR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzNowR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrNowR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpNowR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->ErI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzNowI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrNowI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpNowI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       break;
     case Split :
       D->EzCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       break;
     default :
       printf("what field type?\n");
     }

     // PML
     D->upr=(double *)malloc(nySub2D*sizeof(double ));
     D->upd=(double *)malloc(nySub2D*sizeof(double ));
     D->rtr=(double *)malloc(nxSub1D*sizeof(double ));
     D->rtd=(double *)malloc(nxSub1D*sizeof(double ));
     D->ltr=(double *)malloc(nxSub1D*sizeof(double ));
     D->ltd=(double *)malloc(nxSub1D*sizeof(double ));

     // Particle setting
     nxSub1D=D->nxSub+3;
     nySub2D=D->nySub+3;
     startj=D->jstart-1;

     D->particle = (Particle **)malloc((nxSub1D)*sizeof(Particle *));
     for(i=0; i<nxSub1D; i++) 
       D->particle[i] = (Particle *)malloc((nySub2D)*sizeof(Particle ));

     // setting up particle's pointer
     for(i=0; i<D->iend+1; i++)	//i starts at 0 because of boost frame
       for(j=D->jstart-1; j<D->jend+1; j++)   {
         D->particle[i][j].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)           {
           D->particle[i][j].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j].head[s]->pt = NULL;
         }
       }

     D->laserI = (double *)malloc((D->nySub+5)*sizeof(double ));
     D->laserPhase = (double *)malloc((D->nySub+5)*sizeof(double ));
     for(i=0; i<D->nySub+5; i++) {
       D->laserI[i]=0.0;
       D->laserPhase[i]=0.0;
     }

     add=OFF;
     L=D->laserList;
     while(L->next) {
       if(L->add==File) add=ON; else ;
       L=L->next;
     }

     if(myrank==0) {
       sprintf(name,"laserIn");
       if(add==ON) {
         in=fopen(name,"r");
         cnt=0;
         while(fscanf(in,"%lf %lf %lf %lf %lf %lf %lf",&x,&tmp,&tmp,&tmp,&tmp,&tmp,&tmp)!=EOF) cnt++;
         fclose(in);
         dataX = (double *)malloc(cnt*sizeof(double ));
         intensity = (double *)malloc(cnt*sizeof(double ));
         phase = (double *)malloc(cnt*sizeof(double ));
         in=fopen("laserIn","r");
         for(i=0; i<cnt; i++) {
           fscanf(in,"%lf %lf %lf %lf %lf %lf %lf",&x,&tmp,&tmp,&tmp,&intensity[i],&tmp,&phase[i]);
           dataX[i]=x/D->lambda/D->dr;
         }
         fclose(in);

         dR=dataX[1]-dataX[0];     
         minR=dataX[0];
         for(j=D->jstart; j<D->jend; j++) {
           y=j-D->jstart;
           index=(int)((y-minR)/dR);
           yy=(y-minR)/dR-index;
           if(index<cnt-3) {
             intensity1=intensity[index];
             intensity2=intensity[index+1];
             intensity3=intensity[index+2];
             phase1=phase[index];
             phase2=phase[index+1];
             phase3=phase[index+2];
             aa=(intensity3+intensity1)*0.5-intensity2;
             bb=2.0*intensity2-1.5*intensity1-0.5*intensity3;
             cc=intensity1;
             D->laserI[j]=aa*yy*yy+bb*yy+cc;
             aa=(phase3+phase1)*0.5-phase2;
             bb=2.0*phase2-1.5*phase1-0.5*phase3;
             cc=phase1;
             D->laserPhase[j]=aa*yy*yy+bb*yy+cc;
           } else ;  
         }

         free(dataX); free(intensity); free(phase);
       } else ;
     }	//End of myrank==0
     MPI_Bcast(D->laserI,D->nySub+5,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(D->laserPhase,D->nySub+5,MPI_INT,0,MPI_COMM_WORLD);

/*
    //Making track Particle memory
    trackStart=D->dumpStep/D->trackSaveStep;
    if(trackStart==0)  D->trackStart=D->dumpStep;
    else 
      D->trackStart=(trackStart+1)*D->trackSaveStep;

    numberData=(D->maxStep-D->trackStart)/D->trackSaveStep+1;
    D->track = (Track **)malloc(D->idNums*sizeof(Track *));
    for(i=0; i<D->idNums; i++)
      D->track[i] = (Track *)malloc(numberData*sizeof(Track ));
    for(i=0; i<D->idNums; i++)
      for(j=0; j<numberData; j++)
      {
        D->track[i][j].x=0.0;
        D->track[i][j].y=0.0;
        D->track[i][j].z=0.0;
        D->track[i][j].px=0.0;
        D->track[i][j].py=0.0;
        D->track[i][j].pz=0.0;
        D->track[i][j].step=0;
        D->track[i][j].id=0;
        D->track[i][j].core=0;
      }
*/
/*
     D->probe = (Probe **)malloc(D->probeNum*sizeof(Probe *));
     for(i=0; i<D->probeNum; i++)
       D->probe[i] = (Probe *)malloc((D->maxStep+1)*sizeof(Probe ));
     for(i=0; i<D->probeNum; i++)
       for(j=0; j<=D->maxStep; j++)
       {
         D->probe[i][j].Ex=0.0;
         D->probe[i][j].Ey=0.0;
         D->probe[i][j].Ez=0.0;
         D->probe[i][j].Bx=0.0;
         D->probe[i][j].By=0.0;
         D->probe[i][j].Bz=0.0;
       }
*/

}

void reBoundary(Domain *D,External *Ext)
{
     FILE *out;
     int i,j,s,rank,rankX,rankY,trackStart;
     int startj,nxSub1D,nySub2D;
     int myrank, nTasks;
     MPI_Status status;

     MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
     MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

     // Field setting
     nxSub1D=D->nxSub+5;
     nySub2D=D->nySub+5;

     //Density memory
     D->RhoNoPairR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->RhoNoPairI=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->RhoPairR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->RhoPairI=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->FR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->FI=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->CnR=memoryAsign(D->numMode,nxSub1D,nySub2D);
     D->CnI=memoryAsign(D->numMode,nxSub1D,nySub2D);

     switch(D->fieldType)  {
     case Yee :
     case NoCherenkov :
       D->EzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->ErR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzNowR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrNowR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpNowR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->ErI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzNowI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BrNowI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BpNowI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       break;
     case Split :
       D->EzCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpCR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpCI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpR=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->EzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->BzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->PlI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->SlI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JzI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JrI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       D->JpI=memoryAsign(D->numMode,nxSub1D,nySub2D);
       break;
     default :
       printf("what field type?\n");
     }

     // PML
     D->upr=(double *)malloc(nySub2D*sizeof(double ));
     D->upd=(double *)malloc(nySub2D*sizeof(double ));
     D->rtr=(double *)malloc(nxSub1D*sizeof(double ));
     D->rtd=(double *)malloc(nxSub1D*sizeof(double ));
     D->ltr=(double *)malloc(nxSub1D*sizeof(double ));
     D->ltd=(double *)malloc(nxSub1D*sizeof(double ));


     // Particle setting
     nxSub1D=D->nxSub+3;
     nySub2D=D->nySub+3;
     startj=D->jstart-1;

     D->particle = (Particle **)malloc((nxSub1D)*sizeof(Particle *));
     for(i=0; i<nxSub1D; i++) 
       D->particle[i] = (Particle *)malloc((nySub2D)*sizeof(Particle ));

     // setting up particle's pointer
     for(i=0; i<D->iend+1; i++)	//i starts at 0 because of boost frame
       for(j=D->jstart-1; j<D->jend+1; j++)   {
         D->particle[i][j].head = (ptclHead **)malloc(D->nSpecies*sizeof(ptclHead *));
         for(s=0; s<D->nSpecies; s++)           {
           D->particle[i][j].head[s] = (ptclHead *)malloc(sizeof(ptclHead));
           D->particle[i][j].head[s]->pt = NULL;
         }
       }
}

double ***memoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   double ***field;

   field = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)   {
     field[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       field[i][j] = (double *)malloc((nz)*sizeof(double ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         field[i][j][k]=0.0;

   return field;
}

PML ***PMLmemoryAsign(int nx, int ny, int nz)
{
   int i,j,k;
   PML ***field;

   field = (PML ***)malloc((nx)*sizeof(PML **));
   for(i=0; i<nx; i++)   {
     field[i] = (PML **)malloc((ny)*sizeof(PML *));
     for(j=0; j<ny; j++)
       field[i][j] = (PML *)malloc((nz)*sizeof(PML ));
   }

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++) {
         field[i][j][k].EzR=0.0; field[i][j][k].ErR=0.0;
         field[i][j][k].EpzR=0.0; field[i][j][k].EprR=0.0;
         field[i][j][k].BzR=0.0; field[i][j][k].BrR=0.0;
         field[i][j][k].BpzR=0.0; field[i][j][k].BprR=0.0;
         field[i][j][k].EzI=0.0; field[i][j][k].ErI=0.0;
         field[i][j][k].EpzI=0.0; field[i][j][k].EprI=0.0;
         field[i][j][k].BzI=0.0; field[i][j][k].BrI=0.0;
         field[i][j][k].BpzI=0.0; field[i][j][k].BprI=0.0;
       }
   return field;
}


double ***memoryAsignJ(int nx, int ny, int nz)
{
   int i,j,k;
   double ***current;

   current = (double ***)malloc((nx)*sizeof(double **));
   for(i=0; i<nx; i++)   {
     current[i] = (double **)malloc((ny)*sizeof(double *));
     for(j=0; j<ny; j++)
       current[i][j] = (double *)malloc((nz)*sizeof(double ));
   }
   
   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         current[i][j][k]=0.0;

   return current;
}
