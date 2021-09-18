#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include <mpi.h>
#include "hdf5.h"
#include "hdf5_hl.h"


void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,int *offSet);
void restoreDumpHDF(Domain *D,int iteration,char *partFile,char *fieldFile);

void restoreDump(Domain D,int iteration)
{
  char partFile[100],fieldFile[100];

  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  if(D.saveDumpMode==HDF) {
    sprintf(fieldFile,"dumpField%d.h5",iteration);
    sprintf(partFile,"dumpParticle%d.h5",iteration);
    restoreDumpHDF(&D,iteration,partFile,fieldFile);
  } else ;

}


void restoreRedistHDF(Domain *D,int iteration,char *partFile,char *fieldFile)
{
    int i,j,s,n,istart,iend,jstart,jend,nx,ny,numMode;
    int indexI,indexJ;
    int nxSub,nySub,nSpecies,totalCnt,start,core,dataCnt;
    int minXDomain,minYDomain,minZDomain,maxXDomain,maxYDomain;
    int L,M,flag,shareCnt,fromCore,targetCore;
    int startIndexX,startIndexY,unitX,unitY;
    int rankX,rankY,tmp,shift,rank,cntSub,remain,sub;
    int *recv,*minXSub,*maxXSub,*minYSub,*maxYSub;
    int *sharePNum,*recvDataCnt,*dataCore,*coreCnt,*dataIndex,*dataCores,*nxCore,*share;
    double *data;
    double **sendData,**recvData,*shareData;
    double x,y,z,r;
    int offset[3];
    char name[100],dataName[100];
    void restoreIntMeta();
    ptclList *p;
    Particle **particle;
    particle=D->particle;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    nxSub=D->nxSub;      nySub=D->nySub;
    L=D->L;  M=D->M;
    numMode=D->numMode;

    rankX=myrank/D->M;
    rankY=myrank%D->M;

//    sprintf(name,"dumpField%d.h5",iteration);
    if(myrank==0)  {
      restoreIntMeta(fieldFile,"/minXDomain",&D->minXDomain,1);
    }    else	;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&D->minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);

    nx=D->nx+5; nxSub+=5; istart=0; iend+=3; 
    ny=D->ny+5; nySub+=5; jstart=0; jend+=3;
    offset[0]=D->minXSub;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;
//printf("myrank=%d,offset[0]=%d,offset[1]=%d,nxSub=%d,nySub=%d\n",myrank,offset[0],offset[1],nxSub,nySub);

    restoreFieldComp(D->EzR,fieldFile,"/EzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->EzI,fieldFile,"/EzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->BzR,fieldFile,"/BzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->BzI,fieldFile,"/BzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JzR,fieldFile,"/JzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JzI,fieldFile,"/JzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JrR,fieldFile,"/JrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JrI,fieldFile,"/JrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JpR,fieldFile,"/JpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JpI,fieldFile,"/JpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->RhoPairR,fieldFile,"/RhoPairR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->RhoPairI,fieldFile,"/RhoPairI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->FR,fieldFile,"/FR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->FI,fieldFile,"/FI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    if(D->fieldType!=Split) {
      restoreFieldComp(D->ErR,fieldFile,"/ErR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->ErI,fieldFile,"/ErI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EpR,fieldFile,"/EpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EpI,fieldFile,"/EpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BrR,fieldFile,"/BrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BrI,fieldFile,"/BrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BpR,fieldFile,"/BpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BpI,fieldFile,"/BpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    } 
    else {
      restoreFieldComp(D->PrR,fieldFile,"/PrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PrI,fieldFile,"/PrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlR,fieldFile,"/PlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlI,fieldFile,"/PlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrR,fieldFile,"/SrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrI,fieldFile,"/SrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlR,fieldFile,"/SlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlI,fieldFile,"/SlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EzCR,fieldFile,"/EzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EzCI,fieldFile,"/EzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BzCR,fieldFile,"/BzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BzCI,fieldFile,"/BzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PrCR,fieldFile,"/PrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PrCI,fieldFile,"/PrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlCR,fieldFile,"/PlCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlCI,fieldFile,"/PlCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrCR,fieldFile,"/SrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrCI,fieldFile,"/SrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlCR,fieldFile,"/SlCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlCI,fieldFile,"/SlCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JzCR,fieldFile,"/JzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JzCI,fieldFile,"/JzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JpCR,fieldFile,"/JpCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JpCI,fieldFile,"/JpCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JrCR,fieldFile,"/JrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JrCI,fieldFile,"/JrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    }

    //particle resotre
    nx=D->nx;
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    recv=(int *)malloc(nTasks*sizeof(int ));
    minXSub=(int *)malloc(nTasks*sizeof(int ));
    maxXSub=(int *)malloc(nTasks*sizeof(int ));
    minYSub=(int *)malloc(nTasks*sizeof(int ));
    maxYSub=(int *)malloc(nTasks*sizeof(int ));
    sharePNum=(int *)malloc(nTasks*sizeof(int ));
    recvDataCnt=(int *)malloc(nTasks*sizeof(int ));
    recvData=(double **)malloc(nTasks*sizeof(double *));
    sendData=(double **)malloc(nTasks*sizeof(double *));
    coreCnt=(int *)malloc(nTasks*sizeof(int ));

    nxCore=(int *)malloc(nx*sizeof(int ));
    share=(int *)malloc(nx*sizeof(int ));
    for(i=0; i<nx; i++) { nxCore[i]=0; share[i]=0; }
    for(i=istart; i<iend; i++) {
      indexI=i-istart+D->minXSub;
      nxCore[indexI]=myrank;
    }
    for(i=0; i<nx; i++) share[i]=nxCore[i];

    targetCore=myrank%D->M;
    rankX=myrank/D->M;
    for(n=1; n<D->L; n++)  {
      fromCore=targetCore+n*D->M;
      if(myrank==fromCore)
        MPI_Send(share,nx,MPI_INT,targetCore,myrank,MPI_COMM_WORLD);
      else ;
    }
    if(myrank==targetCore) {
      for(n=1; n<D->L; n++)  {
        fromCore=targetCore+n*D->M;
        MPI_Recv(share,nx,MPI_INT,fromCore,fromCore,MPI_COMM_WORLD,&mpi_status);
        for(i=0; i<nx; i++) nxCore[i]+=share[i];
      }
    } else ;  //End of if(myrank=targetCore)
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(nxCore,nx,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Gather(&D->minXSub,1,MPI_INT,minXSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxXSub,1,MPI_INT,maxXSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->minYSub,1,MPI_INT,minYSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxYSub,1,MPI_INT,maxYSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);

    //restore particle
//    sprintf(name,"dumpParticle%d.h5",iteration);
    if(myrank==0)  {
      restoreIntMeta(partFile,"/nSpecies",&nSpecies,1);
    }    else	;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

    int cntList[nSpecies];

    for(s=0; s<nSpecies; s++)  {
      sprintf(dataName,"%dtotalCnt",s);
      if(myrank==0)  
        restoreIntMeta(partFile,dataName,&cntList[s],1);
      else	;
      MPI_Bcast(&cntList[s],1,MPI_INT,0,MPI_COMM_WORLD);
    }

    //set HDF5 parameter
    hid_t file_id,dset_id,plist_id,attr_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
    herr_t ierr;

    minXDomain=0;
    maxXDomain=D->nx;
    minYDomain=D->minYDomain;
    maxYDomain=D->minYDomain+D->ny;
  
    for(s=0; s<nSpecies; s++)
    {
      totalCnt=cntList[s];
      if(totalCnt>0)
      {
        //open file
        plist_id=H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
        file_id=H5Fopen(partFile,H5F_ACC_RDWR,plist_id);
        H5Pclose(plist_id);

        //set dataset
        sprintf(dataName,"%d",s);
        dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

        sub=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        for(rank=0; rank<nTasks; rank++) {
          if(rank<remain)  tmp=sub+1;
          else	         tmp=sub;
          if(myrank==rank)
            cntSub=tmp;
        }
        MPI_Gather(&cntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        for(i=0; i<nTasks; i++)  {
          sharePNum[i]=0;
          recvDataCnt[i]=0;
          coreCnt[i]=0;
        }

        dataCnt=10;
        data = (double *)malloc(cntSub*dataCnt*sizeof(double ));

        //file space
        dimsf[0]=totalCnt;
        dimsf[1]=dataCnt;
        filespace=H5Screate_simple(2,dimsf,NULL);

        //memory space
        dimsf[0]=cntSub;
        dimsf[1]=dataCnt;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;
        stride[1]=1;
        count[0]=1;
        count[1]=1;

        //hyperslab in file space
        block[0]=cntSub;
        block[1]=dataCnt;
        offSet[0]=start;
        offSet[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);

        //hyperslab in memory space
        offSet[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);

        //read data
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);

        H5Dclose(dset_id);
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Fclose(file_id);

        //for testing the core of each particle
        dataCore = (int *)malloc(cntSub*sizeof(int ));       
        for(i=0; i<cntSub; i++)
          dataCore[i]=-1;

        for(i=0; i<cntSub; i++)  {          
          z=data[i*dataCnt+0]-D->minXDomain;         
          x=data[i*dataCnt+1];          
          y=data[i*dataCnt+2];      
          r=sqrt(x*x+y*y)+D->minYSub;    
          indexI=(int)z;
//          startIndexX=((int)z)/unitX;
//          startIndexY=((int)(r-D->minYDomain))/unitY;
//          rank=startIndexY+startIndexX*M;
          rank=nxCore[indexI];
          flag=0;

//          while(flag==0 && rank<nTasks)  {
            if(minXSub[rank]<=z && z<maxXSub[rank] &&
               minYSub[rank]<=r && r<maxYSub[rank]) { 
              flag=1;           
              if(rank==myrank)  {
                indexI=(int)(z-D->minXSub)+D->istart;
                indexJ=(int)(r-D->minYSub)+D->jstart;

                p = (ptclList *)malloc(sizeof(ptclList));
                p->next = particle[indexI][indexJ].head[s]->pt;
                particle[indexI][indexJ].head[s]->pt=p;
            
                p->z=(z-D->minXSub)-(int)(z-D->minXSub);
                p->x=x;
                p->y=y;
                p->pz=data[i*dataCnt+3];
                p->px=data[i*dataCnt+4];
                p->py=data[i*dataCnt+5];
                p->index=data[i*dataCnt+6];
                p->core=data[i*dataCnt+7];
                p->weight=data[i*dataCnt+8];
                p->charge=data[i*dataCnt+9];
              }
              dataCore[i]=rank;
              sharePNum[rank]+=1;
            } else if(z<minXDomain || z>=maxXDomain 
                   || r<minYDomain || r>=maxYDomain) { 
              dataCore[i]=-1;
              flag=1;
            } else ;
        }

        //set count '0' at myrank due to no reason for sharing
        for(i=0; i<nTasks; i++)
          if(myrank==i)  sharePNum[i]=0;

        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    
            MPI_Send(&sharePNum[i],1,MPI_INT,i,myrank,MPI_COMM_WORLD);         
        }
        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    {
            MPI_Recv(&recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&mpi_status);
          }  else	;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //set memory for send and recving data
        for(i=0; i<nTasks; i++)   {          
          sendData[i]=(double *)malloc(sharePNum[i]*dataCnt*sizeof(double ));
          recvData[i]=(double *)malloc(recvDataCnt[i]*dataCnt*sizeof(double ));
        }        

        for(i=0; i<cntSub; i++)  
        {          
          core=dataCore[i];
          if(myrank!=core && core>=0 && core<nTasks)  {
            n=coreCnt[core];

            sendData[core][n*dataCnt+0]=data[i*dataCnt+0]-D->minXDomain;
            sendData[core][n*dataCnt+1]=data[i*dataCnt+1];
            sendData[core][n*dataCnt+2]=data[i*dataCnt+2];
            sendData[core][n*dataCnt+3]=data[i*dataCnt+3];
            sendData[core][n*dataCnt+4]=data[i*dataCnt+4];
            sendData[core][n*dataCnt+5]=data[i*dataCnt+5];
            sendData[core][n*dataCnt+6]=data[i*dataCnt+6];
            sendData[core][n*dataCnt+7]=data[i*dataCnt+7];
            sendData[core][n*dataCnt+8]=data[i*dataCnt+8];
            sendData[core][n*dataCnt+9]=data[i*dataCnt+9];
            coreCnt[core]+=1;
          }	else	;
        }

        for(i=0; i<nTasks; i++)
        {
          if(myrank==i)  {
            for(j=0; j<nTasks; j++)
              if(i!=j)
                MPI_Send(sendData[j],sharePNum[j]*dataCnt,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);   
          } 
          else  {
            MPI_Recv(recvData[i],recvDataCnt[i]*dataCnt,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&mpi_status);
            for(j=0; j<recvDataCnt[i]; j++)  {
              z=recvData[i][j*dataCnt+0];
              x=recvData[i][j*dataCnt+1];
              y=recvData[i][j*dataCnt+2];
              r=sqrt(x*x+y*y);
              indexI=(int)(z-D->minXSub)+D->istart;
              indexJ=(int)(r)+D->jstart;
              p = (ptclList *)malloc(sizeof(ptclList));
              p->next = particle[indexI][indexJ].head[s]->pt;
              particle[indexI][indexJ].head[s]->pt=p;
            
              p->z=(z-D->minXSub)-(int)(z-D->minXSub);
              p->x=x;
              p->y=y;
              p->pz=recvData[i][j*dataCnt+3];
              p->px=recvData[i][j*dataCnt+4];
              p->py=recvData[i][j*dataCnt+5];
              p->index=recvData[i][j*dataCnt+6];
              p->core=recvData[i][j*dataCnt+7];
              p->weight=recvData[i][j*dataCnt+8];
              p->charge=recvData[i][j*dataCnt+9];
            }
          } 
          MPI_Barrier(MPI_COMM_WORLD);
        }

        free(data);
        free(dataCore);
        for(i=0; i<nTasks; i++)   {          
          free(recvData[i]);
          free(sendData[i]);
        }        

      }	else ;	//End of totalCnt>0
    } 		//End of nSpecies 

    free(recv);
    free(minXSub);
    free(maxXSub);
    free(minYSub);
    free(maxYSub);
    free(sharePNum);
    free(recvDataCnt);
    free(coreCnt);
    free(recvData);
    free(sendData);
    free(nxCore);
    free(share);

}

void restoreDumpHDF(Domain *D,int iteration,char *partFile,char *fieldFile)
{
    int i,j,s,n,istart,iend,jstart,jend,nx,ny,numMode;
    int indexI,indexJ;
    int nxSub,nySub,nSpecies,totalCnt,start,core,dataCnt;
    int minXDomain,minYDomain,minZDomain,maxXDomain,maxYDomain;
    int L,M,flag,shareCnt,fromCore,targetCore;
    int startIndexX,startIndexY,unitX,unitY;
    int rankX,rankY,tmp,shift,rank,cntSub,remain,sub;
    int *recv,*minXSub,*maxXSub,*minYSub,*maxYSub;
    int *sharePNum,*recvDataCnt,*dataCore,*coreCnt,*dataIndex,*dataCores,*nxCore,*share;
    double *data;
    double **sendData,**recvData,*shareData;
    double x,y,z,r;
    int offset[3];
    char name[100],dataName[100];
    void restoreIntMeta();
    ptclList *p;
    Particle **particle;
    particle=D->particle;

    int myrank, nTasks;    
    MPI_Status mpi_status;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    nxSub=D->nxSub;      nySub=D->nySub;
    L=D->L;  M=D->M;
    numMode=D->numMode;

    rankX=myrank/D->M;
    rankY=myrank%D->M;

//    sprintf(name,"dumpField%d.h5",iteration);
    if(myrank==0)  {
      restoreIntMeta(fieldFile,"/minXDomain",&D->minXDomain,1);
    }    else	;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&D->minXDomain,1,MPI_INT,0,MPI_COMM_WORLD);

    nx=D->nx+5; nxSub+=5; istart=0; iend+=3; 
    ny=D->ny+5; nySub+=5; jstart=0; jend+=3;
    offset[0]=D->minXSub;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;
      
    restoreFieldComp(D->EzR,fieldFile,"/EzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->EzI,fieldFile,"/EzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JzR,fieldFile,"/JzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JzI,fieldFile,"/JzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JrR,fieldFile,"/JrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JrI,fieldFile,"/JrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JpR,fieldFile,"/JpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->JpI,fieldFile,"/JpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->RhoPairR,fieldFile,"/RhoPairR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    restoreFieldComp(D->RhoPairI,fieldFile,"/RhoPairI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    if(D->fieldType!=Split) {
      restoreFieldComp(D->ErR,fieldFile,"/ErR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->ErI,fieldFile,"/ErI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EpR,fieldFile,"/EpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EpI,fieldFile,"/EpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BzR,fieldFile,"/BzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BzI,fieldFile,"/BzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BrR,fieldFile,"/BrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BrI,fieldFile,"/BrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BpR,fieldFile,"/BpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->BpI,fieldFile,"/BpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    } else {
      restoreFieldComp(D->PrR,fieldFile,"/PrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PrI,fieldFile,"/PrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlR,fieldFile,"/PlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlI,fieldFile,"/PlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrR,fieldFile,"/SrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrI,fieldFile,"/SrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlR,fieldFile,"/SlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlI,fieldFile,"/SlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EzCR,fieldFile,"/EzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->EzCI,fieldFile,"/EzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PrCR,fieldFile,"/PrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PrCI,fieldFile,"/PrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlCR,fieldFile,"/PlCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->PlCI,fieldFile,"/PlCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrCR,fieldFile,"/SrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SrCI,fieldFile,"/SrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlCR,fieldFile,"/SlCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->SlCI,fieldFile,"/SlCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JzCR,fieldFile,"/JzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JzCI,fieldFile,"/JzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JpCR,fieldFile,"/JpCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JpCI,fieldFile,"/JpCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JrCR,fieldFile,"/JrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      restoreFieldComp(D->JrCI,fieldFile,"/JrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    }

    //particle resotre
    nx=D->nx;
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;

    recv=(int *)malloc(nTasks*sizeof(int ));
    minXSub=(int *)malloc(nTasks*sizeof(int ));
    maxXSub=(int *)malloc(nTasks*sizeof(int ));
    minYSub=(int *)malloc(nTasks*sizeof(int ));
    maxYSub=(int *)malloc(nTasks*sizeof(int ));
    sharePNum=(int *)malloc(nTasks*sizeof(int ));
    recvDataCnt=(int *)malloc(nTasks*sizeof(int ));
    recvData=(double **)malloc(nTasks*sizeof(double *));
    sendData=(double **)malloc(nTasks*sizeof(double *));
    coreCnt=(int *)malloc(nTasks*sizeof(int ));

    nxCore=(int *)malloc(nx*sizeof(int ));
    share=(int *)malloc(nx*sizeof(int ));
    for(i=0; i<nx; i++) { nxCore[i]=0; share[i]=0; }
    for(i=istart; i<iend; i++) {
      indexI=i-istart+D->minXSub;
      nxCore[indexI]=myrank;
    }
    for(i=0; i<nx; i++) share[i]=nxCore[i];

    targetCore=myrank%D->M;
    rankX=myrank/D->M;
    for(n=1; n<D->L; n++)  {
      fromCore=targetCore+n*D->M;
      if(myrank==fromCore)
        MPI_Send(share,nx,MPI_INT,targetCore,myrank,MPI_COMM_WORLD);
      else ;
    }
    if(myrank==targetCore) {
      for(n=1; n<D->L; n++)  {
        fromCore=targetCore+n*D->M;
        MPI_Recv(share,nx,MPI_INT,fromCore,fromCore,MPI_COMM_WORLD,&mpi_status);
        for(i=0; i<nx; i++) nxCore[i]+=share[i];
      }
    } else ;  //End of if(myrank=targetCore)
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(nxCore,nx,MPI_INT,0,MPI_COMM_WORLD);

    MPI_Gather(&D->minXSub,1,MPI_INT,minXSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxXSub,1,MPI_INT,maxXSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxXSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->minYSub,1,MPI_INT,minYSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(minYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Gather(&D->maxYSub,1,MPI_INT,maxYSub,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxYSub,nTasks,MPI_INT,0,MPI_COMM_WORLD);

    //restore particle
//    sprintf(name,"dumpParticle%d.h5",iteration);
    if(myrank==0)  {
      restoreIntMeta(partFile,"/nSpecies",&nSpecies,1);
    }    else	;
    MPI_Bcast(&nSpecies,1,MPI_INT,0,MPI_COMM_WORLD);

    int cntList[nSpecies];

    for(s=0; s<nSpecies; s++)  {
      sprintf(dataName,"%dtotalCnt",s);
      if(myrank==0)  
        restoreIntMeta(partFile,dataName,&cntList[s],1);
      else	;
      MPI_Bcast(&cntList[s],1,MPI_INT,0,MPI_COMM_WORLD);
    }

    //set HDF5 parameter
    hid_t file_id,dset_id,plist_id,attr_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offSet[2],block[2],stride[2];
    herr_t ierr;

    minXDomain=0;
    maxXDomain=D->nx;
    minYDomain=D->minYDomain;
    maxYDomain=D->minYDomain+D->ny;
  
    for(s=0; s<nSpecies; s++)
    {
      totalCnt=cntList[s];
      if(totalCnt>0)
      {
        //open file
        plist_id=H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
        file_id=H5Fopen(partFile,H5F_ACC_RDWR,plist_id);
        H5Pclose(plist_id);

        //set dataset
        sprintf(dataName,"%d",s);
        dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);

        sub=totalCnt/nTasks;
        remain=totalCnt%nTasks;
        for(rank=0; rank<nTasks; rank++) {
          if(rank<remain)  tmp=sub+1;
          else	         tmp=sub;
          if(myrank==rank)
            cntSub=tmp;
        }
        MPI_Gather(&cntSub,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);
        start=0;
        for(i=0; i<myrank; i++)
          start+=recv[i];

        for(i=0; i<nTasks; i++)  {
          sharePNum[i]=0;
          recvDataCnt[i]=0;
          coreCnt[i]=0;
        }

        dataCnt=10;
        data = (double *)malloc(cntSub*dataCnt*sizeof(double ));

        //file space
        dimsf[0]=totalCnt;
        dimsf[1]=dataCnt;
        filespace=H5Screate_simple(2,dimsf,NULL);

        //memory space
        dimsf[0]=cntSub;
        dimsf[1]=dataCnt;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;
        stride[1]=1;
        count[0]=1;
        count[1]=1;

        //hyperslab in file space
        block[0]=cntSub;
        block[1]=dataCnt;
        offSet[0]=start;
        offSet[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offSet,stride,count,block);

        //hyperslab in memory space
        offSet[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offSet,stride,count,block);

        //read data
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        ierr=H5Dread(dset_id,H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);

        H5Dclose(dset_id);
        H5Pclose(plist_id);
        H5Sclose(memspace);
        H5Sclose(filespace);
        H5Fclose(file_id);

        //for testing the core of each particle
        dataCore = (int *)malloc(cntSub*sizeof(int ));       
        for(i=0; i<cntSub; i++)
          dataCore[i]=-1;

        for(i=0; i<cntSub; i++)  {          
          z=data[i*dataCnt+0]-D->minXDomain;         
          x=data[i*dataCnt+1];          
          y=data[i*dataCnt+2];          
          r=sqrt(x*x+y*y);
          indexI=(int)z;
//          startIndexX=((int)z)/unitX;
//          startIndexY=((int)(r-D->minYDomain))/unitY;
//          rank=startIndexY+startIndexX*M;
          rank=nxCore[indexI];
          flag=0;

//          while(flag==0 && rank<nTasks)  {
            if(minXSub[rank]<=z && z<maxXSub[rank] &&
               minYSub[rank]<=r && r<maxYSub[rank]) { 
              flag=1;           
              if(rank==myrank)  {
                indexI=(int)(z-D->minXSub)+D->istart;
                indexJ=(int)(r)+D->jstart;

                p = (ptclList *)malloc(sizeof(ptclList));
                p->next = particle[indexI][indexJ].head[s]->pt;
                particle[indexI][indexJ].head[s]->pt=p;
            
                p->z=(z-D->minXSub)-(int)(z-D->minXSub);
                p->x=x;
                p->y=y;
                p->pz=data[i*dataCnt+3];
                p->px=data[i*dataCnt+4];
                p->py=data[i*dataCnt+5];
                p->index=data[i*dataCnt+6];
                p->core=data[i*dataCnt+7];
                p->weight=data[i*dataCnt+8];
                p->charge=data[i*dataCnt+9];
              }
              dataCore[i]=rank;
              sharePNum[rank]+=1;
            } else if(z<minXDomain || z>=maxXDomain 
                   || r<minYDomain || r>=maxYDomain) { 
              dataCore[i]=-1;
              flag=1;
            } else ;
        }

        //set count '0' at myrank due to no reason for sharing
        for(i=0; i<nTasks; i++)
          if(myrank==i)  sharePNum[i]=0;

        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    
            MPI_Send(&sharePNum[i],1,MPI_INT,i,myrank,MPI_COMM_WORLD);         
        }
        for(i=0; i<nTasks; i++)   {          
          if(myrank!=i)    {
            MPI_Recv(&recvDataCnt[i],1,MPI_INT,i,i,MPI_COMM_WORLD,&mpi_status);
          }  else	;
        }
        MPI_Barrier(MPI_COMM_WORLD);

        //set memory for send and recving data
        for(i=0; i<nTasks; i++)   {          
          sendData[i]=(double *)malloc(sharePNum[i]*dataCnt*sizeof(double ));
          recvData[i]=(double *)malloc(recvDataCnt[i]*dataCnt*sizeof(double ));
        }        

        for(i=0; i<cntSub; i++)  
        {          
          core=dataCore[i];
          if(myrank!=core && core>=0 && core<nTasks)  {
            n=coreCnt[core];

            sendData[core][n*dataCnt+0]=data[i*dataCnt+0]-D->minXDomain;
            sendData[core][n*dataCnt+1]=data[i*dataCnt+1];
            sendData[core][n*dataCnt+2]=data[i*dataCnt+2];
            sendData[core][n*dataCnt+3]=data[i*dataCnt+3];
            sendData[core][n*dataCnt+4]=data[i*dataCnt+4];
            sendData[core][n*dataCnt+5]=data[i*dataCnt+5];
            sendData[core][n*dataCnt+6]=data[i*dataCnt+6];
            sendData[core][n*dataCnt+7]=data[i*dataCnt+7];
            sendData[core][n*dataCnt+8]=data[i*dataCnt+8];
            sendData[core][n*dataCnt+9]=data[i*dataCnt+9];
            coreCnt[core]+=1;
          }	else	;
        }

        for(i=0; i<nTasks; i++)
        {
          if(myrank==i)  {
            for(j=0; j<nTasks; j++)
              if(i!=j)
                MPI_Send(sendData[j],sharePNum[j]*dataCnt,MPI_DOUBLE,j,myrank,MPI_COMM_WORLD);   
          } 
          else  {
            MPI_Recv(recvData[i],recvDataCnt[i]*dataCnt,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&mpi_status);
            for(j=0; j<recvDataCnt[i]; j++)  {
              z=recvData[i][j*dataCnt+0];
              x=recvData[i][j*dataCnt+1];
              y=recvData[i][j*dataCnt+2];
              r=sqrt(x*x+y*y);
              indexI=(int)(z-D->minXSub)+D->istart;
              indexJ=(int)(r)+D->jstart;
              p = (ptclList *)malloc(sizeof(ptclList));
              p->next = particle[indexI][indexJ].head[s]->pt;
              particle[indexI][indexJ].head[s]->pt=p;
            
              p->z=(z-D->minXSub)-(int)(z-D->minXSub);
              p->x=x;
              p->y=y;
              p->pz=recvData[i][j*dataCnt+3];
              p->px=recvData[i][j*dataCnt+4];
              p->py=recvData[i][j*dataCnt+5];
              p->index=recvData[i][j*dataCnt+6];
              p->core=recvData[i][j*dataCnt+7];
              p->weight=recvData[i][j*dataCnt+8];
              p->charge=recvData[i][j*dataCnt+9];
            }
          } 
          MPI_Barrier(MPI_COMM_WORLD);
        }

        free(data);
        free(dataCore);
        for(i=0; i<nTasks; i++)   {          
          free(recvData[i]);
          free(sendData[i]);
        }        

      }	else ;	//End of totalCnt>0
    } 		//End of nSpecies 

    free(recv);
    free(minXSub);
    free(maxXSub);
    free(minYSub);
    free(maxYSub);
    free(sharePNum);
    free(recvDataCnt);
    free(coreCnt);
    free(recvData);
    free(sendData);
    free(nxCore);
    free(share);

}



void restoreFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,int *offSet)
{
  int i,j,m,start;
  double *field;
  char name[100];
  int myrank, nTasks;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,plist_id;
  herr_t status;
  hid_t subfilespace,filespace,memspace;
  hsize_t dimsf[3],count[3],offset[3];

  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//  H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//  MPI_Barrier(MPI_COMM_WORLD);

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
  H5Pclose(plist_id);
  dimsf[0]=ny; dimsf[1]=nx; dimsf[2]=numMode;     
  filespace=H5Screate_simple(3,dimsf,NULL);

  count[0]=nySub; count[1]=nxSub; count[2]=numMode;
  offset[0]=offSet[1]; offset[1]=offSet[0]; offset[2]=offSet[2];
  memspace=H5Screate_simple(3,count,NULL);

  field = (double *)malloc(nxSub*nySub*numMode*sizeof(double ));

  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  subfilespace=H5Dget_space(dset_id);
  H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
  status = H5Dread(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
  start=0;
  for(j=jstart; j<jend; j++)
    for(i=istart; i<iend; i++) {
      for(m=0; m<numMode; m++)
        data[m][i][j]=field[start+m];
      start+=numMode;
    }
  H5Pclose(plist_id);
  H5Sclose(subfilespace);
  H5Dclose(dset_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fclose(file_id);
  free(field);
}

void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dopen2(file_id,dataName,H5P_DEFAULT);
  status=H5Dread(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}

