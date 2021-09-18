#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);
void saveParticleComp_Double(double *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveParticleComp_Int(int *data,char *fileName,char *dataName,int totalCnt,int cnt,int offSet);
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
double calHighInter(double y1,double y2,double y3,double resolX);
double calLowInter(double y1,double y2,double y3,int resolX);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,int *offSet);
void saveDumpParticleHDF(Domain *D,int iteration,char *fileName);
void saveDumpFieldHDF(Domain *D,int iteration,char *fileName);


void saveDump(Domain D,int iteration)
{
  char particleName[100],fieldName[100];
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  if(D.saveDumpMode==HDF)  {
    sprintf(particleName,"dumpParticle%d.h5",iteration);
    sprintf(fieldName,"dumpField%d.h5",iteration);
    saveDumpFieldHDF(&D,iteration,fieldName);
    saveDumpParticleHDF(&D,iteration,particleName);
    if(myrank==0)  {
      printf("%s\n",particleName);
      printf("%s\n",fieldName);
    }   else	;
  }

}

void saveDumpParticleHDF(Domain *D,int iteration,char *fileName)
{
    int i,j,s,istart,iend,jstart,jend,dataCnt=10;
    int cnt,totalCnt,index,start,cntList[D->nSpecies];
    int minXSub,minYSub,nxSub,nySub;
    char dataName[100];
    double *data;
    int *recv,*offSetRank;
    Particle **particle;
    particle=D->particle;
    ptclList *p;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,plist_id,dset_id,attr_id,as_id;
    hid_t filespace,memspace;
    hsize_t dimsf[2],count[2],offset[2],block[2],stride[2],a_dims,metaDim[1];
    herr_t ierr;

    nxSub=D->nxSub;      nySub=D->nySub;
    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    minXSub=D->minXSub;  minYSub=D->minYSub;

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    ierr=H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
    //create file
//    sprintf(name,"dumpParticle%d.h5",iteration);
    file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
    ierr=H5Pclose(plist_id);

    recv = (int *)malloc(nTasks*sizeof(int ));

    for(s=0; s<D->nSpecies; s++)
    {
      cnt=0;
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++) {
          p=particle[i][j].head[s]->pt;
          while(p)   {
            cnt++;
            p=p->next;
          }
        }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Gather(&cnt,1,MPI_INT,recv,1,MPI_INT,0,MPI_COMM_WORLD);
      MPI_Bcast(recv,nTasks,MPI_INT,0,MPI_COMM_WORLD);

      start=0;
      for(i=0; i<myrank; i++)        start+=recv[i];
      totalCnt=0;
      for(i=0; i<nTasks; i++)        totalCnt+=recv[i];
      cntList[s]=totalCnt;

      //file space
      dimsf[0]=totalCnt;
      dimsf[1]=dataCnt;
      filespace=H5Screate_simple(2,dimsf,NULL);
      sprintf(dataName,"%d",s);
      dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

      if(totalCnt>0)
      {
        data = (double *)malloc(cnt*dataCnt*sizeof(double ));

        index=0;
        for(i=istart; i<iend; i++)
          for(j=jstart; j<jend; j++) {
            p=particle[i][j].head[s]->pt;
            while(p)    {
              data[index*dataCnt+0]=p->z+i-istart+minXSub;
              data[index*dataCnt+1]=p->x;
              data[index*dataCnt+2]=p->y;
              data[index*dataCnt+3]=p->pz;
              data[index*dataCnt+4]=p->px;
              data[index*dataCnt+5]=p->py;
              data[index*dataCnt+6]=p->index;
              data[index*dataCnt+7]=p->core;
              data[index*dataCnt+8]=p->weight;
              data[index*dataCnt+9]=p->charge;
              index++;
              p=p->next;
            }
          }

        //memory space
        dimsf[0]=cnt;
        dimsf[1]=dataCnt;
        memspace=H5Screate_simple(2,dimsf,NULL);

        stride[0]=1;       stride[1]=1;
        count[0]=1;        count[1]=1;

        //hyperslab in file space
        block[0]=cnt;        block[1]=dataCnt;
        offset[0]=start;     offset[1]=0;
        H5Sselect_hyperslab(filespace,H5S_SELECT_SET,offset,stride,count,block);

        //hyperslab in memory space
        offset[0]=0;
        H5Sselect_hyperslab(memspace,H5S_SELECT_SET,offset,stride,count,block);
      
        plist_id=H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_INDEPENDENT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,filespace,plist_id,data);
        H5Pclose(plist_id);
        H5Sclose(memspace);
      
        free(data);
      }	else	; 	//End of totalCnt>0
      MPI_Barrier(MPI_COMM_WORLD);
/*
      //write meta
      a_dims=1;
      as_id=H5Screate_simple(1,&a_dims,NULL);
      sprintf(dataName,"%dtotalCnt",s);
      attr_id=H5Acreate2(dset_id,dataName,H5T_NATIVE_INT,as_id,H5P_DEFAULT,H5P_DEFAULT);
      H5Awrite(attr_id,H5T_NATIVE_INT,&totalCnt);
      H5Aclose(attr_id);
      H5Sclose(as_id);
*/
      H5Dclose(dset_id);

      H5Sclose(filespace);
    }	//End of nSpecies
    free(recv);
    H5Fclose(file_id);

    if(myrank==0)  {
      sprintf(dataName,"nSpecies");
      saveIntMeta(fileName,dataName,&D->nSpecies,1);
      for(s=0; s<D->nSpecies; s++)  {
        sprintf(dataName,"%dtotalCnt",s);
        saveIntMeta(fileName,dataName,&cntList[s],1);
      }
      saveIntMeta(fileName,"/nx",&D->nx,1);
      saveIntMeta(fileName,"/ny",&D->ny,1);
      saveIntMeta(fileName,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(fileName,"/minYDomain",&D->minYDomain,1);
    }  else	;
}

void saveDumpFieldHDF(Domain *D,int iteration,char *fileName)
{
    int i,j,m,s,istart,iend,jstart,jend,nx,ny;
    int biasX,biasY,offSetY,rankX,rankY;
    int nxSub,nySub,numMode;
    int offset[3];
    char name[100],name2[100];

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;
    void saveFieldComp();

    nxSub=D->nxSub;    nySub=D->nySub;
    istart=D->istart;  iend=D->iend;
    jstart=D->jstart;  jend=D->jend;
    numMode=D->numMode;

//    sprintf(name,"dumpField%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(fileName,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/D->M;
    rankY=myrank%D->M;

    if(myrank==0)  {
      saveIntMeta(fileName,"/minXDomain",&D->minXDomain,1);
      saveIntMeta(fileName,"/minYDomain",&D->minYDomain,1);
      saveIntMeta(fileName,"/nSpecies",&D->nSpecies,1);
      saveIntMeta(fileName,"/nx",&D->nx,1);
      saveIntMeta(fileName,"/ny",&D->ny,1);
    } else	;
    MPI_Barrier(MPI_COMM_WORLD);
    
    nx=D->nx+5;
    calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D->L);
    ny=D->ny+5;
    calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D->M);

    offset[0]=(D->minXSub-D->minXDomain)+biasX;
    offset[1]=(D->minYSub-D->minYDomain)+biasY;
    offset[2]=0;
 
    saveFieldComp(D->EzR,fileName,"/EzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->EzI,fileName,"/EzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->BzR,fileName,"/BzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->BzI,fileName,"/BzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JzR,fileName,"/JzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JzI,fileName,"/JzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JrR,fileName,"/JrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JrI,fileName,"/JrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JpR,fileName,"/JpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JpI,fileName,"/JpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->RhoPairR,fileName,"/RhoPairR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->RhoPairI,fileName,"/RhoPairI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->FR,fileName,"/FR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->FI,fileName,"/FI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);

    if(D->fieldType!=Split) {
      saveFieldComp(D->ErR,fileName,"/ErR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->ErI,fileName,"/ErI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->EpR,fileName,"/EpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->EpI,fileName,"/EpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->BrR,fileName,"/BrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->BrI,fileName,"/BrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->BpR,fileName,"/BpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->BpI,fileName,"/BpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    } else {
      saveFieldComp(D->PrR,fileName,"/PrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PrI,fileName,"/PrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PlR,fileName,"/PlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PlI,fileName,"/PlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SrR,fileName,"/SrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SrI,fileName,"/SrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SlR,fileName,"/SlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SlI,fileName,"/SlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->EzCR,fileName,"/EzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->EzCI,fileName,"/EzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->BzCR,fileName,"/BzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->BzCI,fileName,"/BzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PrCR,fileName,"/PrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PrCI,fileName,"/PrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PlCR,fileName,"/PlCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->PlCI,fileName,"/PlCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SrCR,fileName,"/SrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SrCI,fileName,"/SrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SlCR,fileName,"/SlCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->SlCI,fileName,"/SlCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->JzCR,fileName,"/JzCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->JzCI,fileName,"/JzCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->JrCR,fileName,"/JrCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->JrCI,fileName,"/JrCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->JpCR,fileName,"/JpCR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
      saveFieldComp(D->JpCI,fileName,"/JpCI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    }
}


void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L)
{
  if(L==1) {
    *istart=0; *iend+=3; *biasX=0; *nxSub=nx;
  } else  {
    if(rankX==0)  {
      *istart=0; *nxSub+=2; *biasX=0;
    }  else if(rankX==L-1)  {
      *iend+=3; *nxSub+=3; *biasX=2;
    } else
      *biasX=2;
  }
}

void deleteFieldData(Domain *D,double ***field,int *nxSub,int *nySub,int *nzSub,int rankX,int rankY)
{
   int i,j,k;

   if(rankX==0)  *nxSub-=2;
   else if(rankX==D->L-1)  *nxSub-=3;
   if(rankY==0)  *nySub-=2;
   else if(rankY==D->M-1)  *nySub-=3;

   for(i=0; i<*nxSub+5; i++)
   {
     for(j=0; j<*nySub+5; j++)
       free(field[i][j]);
     free(field[i]);
   }
   free(field);
}


void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,int *offSet)
{
    int ii,i,j,m,start;
    double *field;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[3],count[3],offset[3];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=ny; dimsf[1]=nx; dimsf[2]=numMode;
    filespace=H5Screate_simple(3,dimsf,NULL);

    count[0]=nySub;      count[1]=nxSub;      count[2]=numMode;
    offset[0]=offSet[1]; offset[1]=offSet[0]; offset[2]=offSet[2];
    memspace=H5Screate_simple(3,count,NULL);

    field = (double *)malloc(nxSub*nySub*numMode*sizeof(double ));

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

    start=0;
    for(j=jstart; j<jend; j++) 
      for(i=istart; i<iend; i++)  {
        for(m=0; m<numMode; m++) 
          field[start+m]=data[m][i][j];
        start+=numMode;     
      }

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace,subfilespace,plist_id,field);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
    free(field);
}

void MPI_saveIntArray(int *data,char *fileName,char *dataName,int offSet)
{
    int i;
    FILE *out;
    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,dset_id,plist_id,tic_id;
    herr_t status;
    hid_t total_file_space,subfilespace,filespace,memspace,ticspace;
    hsize_t dimsf[1],count[1],offset[1];

    plist_id=H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id,MPI_COMM_WORLD,MPI_INFO_NULL);
//    H5Pset_fclose_degree(plist_id,H5F_CLOSE_SEMI);
//    MPI_Barrier(MPI_COMM_WORLD);

    file_id=H5Fopen(fileName,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    dimsf[0]=nTasks;
    filespace=H5Screate_simple(1,dimsf,NULL);

    count[0]=1;
    offset[0]=offSet;
    memspace=H5Screate_simple(1,count,NULL);

    dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    subfilespace=H5Dget_space(dset_id);
    H5Sselect_hyperslab(subfilespace,H5S_SELECT_SET,offset,NULL,count,NULL);

    plist_id=H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id,H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_INT,memspace,subfilespace,plist_id,data);
    H5Pclose(plist_id);
    H5Sclose(subfilespace);
    H5Dclose(dset_id);

    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Fclose(file_id);
}

void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_INT,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
      
void saveDoubleMeta(char *fileName,char *dataName,double *data,int dataCnt)
{
  hid_t file_id,dset_id,filespace;
  hsize_t metaDim[1];
  herr_t status;

  metaDim[0]=dataCnt;

  file_id=H5Fopen(fileName,H5F_ACC_RDWR,H5P_DEFAULT);
  filespace=H5Screate_simple(1,metaDim,NULL);
  dset_id=H5Dcreate2(file_id,dataName,H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  status=H5Dwrite(dset_id,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dset_id);
  H5Sclose(filespace);
  H5Fclose(file_id);
}
      

