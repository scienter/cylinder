#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveCoordHDF(Domain *D,char *fileName);
void Efield_xdmf(char *fileName,int nx,int ny);
void Bfield_xdmf(char *fileName,int nx,int ny);
void calCenter(Domain *D,double ***data,double ***field,int shI,int shJ,int istart,int iend,int jstart,int jend,int numMode);

void saveJHDF(Domain *D,int iteration)
{
    int i,j,m,s,istart,iend,jstart,jend,nx,ny,nz,nxSub,nySub,numMode;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    double ***data,sin2,cos2,cosP,sinP,phi=0;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;    nySub=D->nySub;
    istart=D->istart;  iend=D->iend;
    jstart=D->jstart;  jend=D->jend;
    nx=D->nx;          ny=D->ny; nz=1;
    numMode=D->numMode;
    double coss[numMode],sins[numMode];

    sprintf(name,"J%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
      saveIntMeta(name,"/nx",&nx,1);
      saveIntMeta(name,"/ny",&ny,1);
      saveIntMeta(name,"/nz",&nz,1);
      saveIntMeta(name,"/numMode",&numMode,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;


    data=(double ***)malloc(1*sizeof(double **));
    for(m=0; m<1; m++) {
      data[m]=(double **)malloc((nxSub+5)*sizeof(double *));
      for(i=0; i<nxSub+5; i++) 
        data[m][i]=(double *)malloc((nySub+5)*sizeof(double ));
    }

//    cosP=cos(phi*pi/180.0); sinP=sin(phi*pi/180.0);
//    coss[0]=1.0; sins[0]=0.0;
//    for(m=1; m<numMode; m++) {
//      cos2=coss[m-1]; sin2=sins[m-1];
//      coss[m]=cosP*cos2-sinP*sin2;
//      sins[m]=sinP*cos2+cosP*sin2;
//    }

//    calField(data,D->EzR,D->EzI,nxSub+5,nySub+5,istart,iend,jstart,jend,numMode,coss,sins);
    saveFieldComp(D->JzR,name,"/JzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JzI,name,"/JzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JrR,name,"/JrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JrI,name,"/JrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JpR,name,"/JpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->JpI,name,"/JpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);

    for(m=0; m<1; m++) {
      for(i=0; i<nxSub+5; i++) free(data[m][i]);
      free(data[m]);
    } free(data);

    if(myrank==0)   {
      saveCoordHDF(D,name);
      sprintf(fileName,"J%d",iteration);
//      Efield_xdmf(fileName,nx,ny);
      printf("%s\n",name); 
    }	 else  ;
}

void saveSplitFieldHDF(Domain *D,int iteration)
{
    int i,j,m,s,istart,iend,jstart,jend,nx,ny,nz,nxSub,nySub,numMode;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    double ***data,sin2,cos2,cosP,sinP,phi=0;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;    nySub=D->nySub;
    istart=D->istart;  iend=D->iend;
    jstart=D->jstart;  jend=D->jend;
    nx=D->nx;          ny=D->ny; nz=1;
    numMode=D->numMode;
    double coss[numMode],sins[numMode];

//lala
    sprintf(name,"fieldSplit%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
      saveIntMeta(name,"/nx",&nx,1);
      saveIntMeta(name,"/ny",&ny,1);
      saveIntMeta(name,"/nz",&nz,1);
      saveIntMeta(name,"/numMode",&numMode,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;

    saveFieldComp(D->PrR,name,"/PrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->PlR,name,"/PlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->SrR,name,"/SrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->SlR,name,"/SlR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->PrI,name,"/PrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->PlI,name,"/PlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->SrI,name,"/SrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->SlI,name,"/SlI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->EzR,name,"/EzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->EzI,name,"/EzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->BzR,name,"/BzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->BzI,name,"/BzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);

    saveFieldComp(D->FR,name,"/FR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->FI,name,"/FI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    for(m=0; m<numMode; m++) {
//      for(i=0; i<nxSub+5; i++) free(data[m][i]);
//      free(data[m]);
//    } free(data);
    if(myrank==0)   {
      saveCoordHDF(D,name);
      sprintf(fileName,"fieldSplit%d",iteration);
//      Efield_xdmf(fileName,nx,ny);
      printf("%s\n",name); 
    }	 else  ;
}

void saveEFieldHDF(Domain *D,int iteration)
{
    int i,j,m,s,istart,iend,jstart,jend,nx,ny,nz,nxSub,nySub,numMode;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    double ***data,sin2,cos2,cosP,sinP,phi=0;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;    nySub=D->nySub;
    istart=D->istart;  iend=D->iend;
    jstart=D->jstart;  jend=D->jend;
    nx=D->nx;          ny=D->ny; nz=1;
    numMode=D->numMode;
    double coss[numMode],sins[numMode];

    sprintf(name,"field%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
      saveIntMeta(name,"/nx",&nx,1);
      saveIntMeta(name,"/ny",&ny,1);
      saveIntMeta(name,"/nz",&nz,1);
      saveIntMeta(name,"/numMode",&numMode,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;


    data=(double ***)malloc(numMode*sizeof(double **));
    for(m=0; m<numMode; m++) {
      data[m]=(double **)malloc((nxSub+5)*sizeof(double *));
      for(i=0; i<nxSub+5; i++) 
        data[m][i]=(double *)malloc((nySub+5)*sizeof(double ));
    }

    cosP=cos(phi*M_PI/180.0); sinP=sin(phi*M_PI/180.0);
    coss[0]=1.0; sins[0]=0.0;
    for(m=1; m<numMode; m++) {
      cos2=coss[m-1]; sin2=sins[m-1];
      coss[m]=cosP*cos2-sinP*sin2;
      sins[m]=sinP*cos2+cosP*sin2;
    }

//    calCenter(D,data,D->EzR,0,-1,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->EzR,name,"/EzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->EzI,0,-1,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->EzI,name,"/EzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->ErR,-1,0,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->ErR,name,"/ErR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->ErI,-1,0,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->ErI,name,"/ErI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->EpR,-1,-1,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->EpR,name,"/EpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->EpI,-1,-1,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->EpI,name,"/EpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);

//    calCenter(D,data,D->BzNowR,-1,0,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->BzNowR,name,"/BzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->BzNowI,-1,0,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->BzNowI,name,"/BzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->BrNowR,0,-1,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->BrNowR,name,"/BrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//    calCenter(D,data,D->BrNowI,0,-1,istart,iend,jstart,jend,numMode);
    saveFieldComp(D->BrNowI,name,"/BrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->BpNowR,name,"/BpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->BpNowI,name,"/BpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);

    saveFieldComp(D->FR,name,"/FR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    saveFieldComp(D->FI,name,"/FI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
    for(m=0; m<numMode; m++) {
      for(i=0; i<nxSub+5; i++) free(data[m][i]);
      free(data[m]);
    } free(data);

    if(myrank==0)   {
      saveCoordHDF(D,name);
      sprintf(fileName,"field%d",iteration);
//      Efield_xdmf(fileName,nx,ny);
      printf("%s\n",name); 
    }	 else  ;
}
/*
void saveBFieldHDF(Domain *D,int iteration)
{
    int i,j,m,s,istart,iend,jstart,jend,nx,ny,nz,nxSub,nySub,numMode;
    int offset[3];
    char name[100],dataName[100],fileName[100];
    double ***data,sin2,cos2,cosP,sinP,phi=0;
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;    nySub=D->nySub;
    istart=D->istart;  iend=D->iend;
    jstart=D->jstart;  jend=D->jend;
    nx=D->nx;          ny=D->ny; nz=1;
    numMode=D->numMode;
    double coss[numMode],sins[numMode];

    sprintf(name,"fieldB%d.h5",iteration);
    if(myrank==0)     {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
      saveIntMeta(name,"/nx",&nx,1);
      saveIntMeta(name,"/ny",&ny,1);
      saveIntMeta(name,"/nz",&nz,1);
      saveIntMeta(name,"/numMode",&numMode,1);
    } else      ;
    MPI_Barrier(MPI_COMM_WORLD);

    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;


    data=(double ***)malloc(numMode*sizeof(double **));
    for(m=0; m<numMode; m++) {
      data[m]=(double **)malloc((nxSub+5)*sizeof(double *));
      for(i=0; i<nxSub+5; i++) 
        data[m][i]=(double *)malloc((nySub+5)*sizeof(double ));
    }

    calCenter(D,data,D->BzNowR,0,-1);
    saveFieldComp(data,name,"/BzR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,1,offset);
    calCenter(D,data,D->BzNowI,0,-1);
    saveFieldComp(data,name,"/BzI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,1,offset);
    calCenter(D,data,D->BrNowR,-1,0);
    saveFieldComp(data,name,"/BrR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,1,offset);
    calCenter(D,data,D->BrNowI,-1,0);
    saveFieldComp(data,name,"/BrI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,1,offset);
    calCenter(D,data,D->BpNowR,-1,-1);
    saveFieldComp(data,name,"/BpR",nx,ny,nxSub,nySub,istart,iend,jstart,jend,1,offset);
    calCenter(D,data,D->BpNowI,-1,-1);
    saveFieldComp(data,name,"/BpI",nx,ny,nxSub,nySub,istart,iend,jstart,jend,1,offset);

    for(m=0; m<1; m++) {
      for(i=0; i<nxSub+5; i++) free(data[m][i]);
      free(data[m]);
    } free(data);

    if(myrank==0)   {
      saveCoordHDF(D,name);
      sprintf(fileName,"fieldB%d",iteration);
      Efield_xdmf(fileName,nx,ny);
      printf("%s\n",name); 
    }	 else  ;
}
*/
//lala
void calCenter(Domain *D,double ***data,double ***field,int shI,int shJ,int istart,int iend,int jstart,int jend,int numMode)
{
  int i,j,m;

  for(m=0; m<numMode; m++)
    for(i=istart; i<iend; i++) 
      for(j=jstart; j<jend; j++) { 
        data[m][i][j]=0.25*(field[m][i+shI][j+shJ]+field[m][i+shI+1][j+shJ]+field[m][i+shI][j+shJ-1]+field[m][i+shI+1][j+shJ-1]);
      }
  
}

void Efield_xdmf(char *fileName,int nx,int ny)
{
    FILE *xmf = 0;
    char name[100];
    const char *Names[] = {"/Ez","/Er","/Ep"};
    int i;
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"%s.xmf",fileName);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
    fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
    fprintf(xmf, "        %s.h5:/X\n",fileName);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
    fprintf(xmf, "        %s.h5:/Y\n",fileName);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    for(i=0; i<3; i++)
    {
      fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx,ny);
      fprintf(xmf, "        %s.h5:/%s\n",fileName,Names[i]);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
    }
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

void Bfield_xdmf(char *fileName,int nx,int ny)
{
    FILE *xmf = 0;
    char name[100];
    const char *Names[] = {"/Bz","/Br","/Bp"};
    int i;
    //
     // Open the file and write the XML description of the mesh..
     //
    sprintf(name,"%s.xmf",fileName);
    xmf = fopen(name,"w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");

    fprintf(xmf, "     <Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"%d %d\"/>\n", ny,nx);
    fprintf(xmf, "     <Geometry GeometryType=\"VXVY\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
    fprintf(xmf, "        %s.h5:/X\n",fileName);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ny);
    fprintf(xmf, "        %s.h5:/Y\n",fileName);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    for(i=0; i<3; i++)
    {
      fprintf(xmf, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Node\">\n",Names[i]);
      fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx,ny);
      fprintf(xmf, "        %s.h5:/%s\n",fileName,Names[i]);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
    }
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

void calField(double ***data,double ***fieldR,double ***fieldI,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,double *coss,double *sins)
{
  int i,j,m;
  
  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++) {
      data[0][i][j]=fieldR[0][i][j];
      for(m=1; m<numMode; m++)
        data[0][i][j]+=fieldR[m][i][j]*coss[m]-fieldI[m][i][j]*sins[m];
    }
}
