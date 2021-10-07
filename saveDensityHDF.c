#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void deleteField(double ***field,int nx,int ny,int nz);
double ***memoryAsign(int nx, int ny, int nz);

void MPI_TransferDen_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share);
void MPI_TransferDen_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nxSub,int nySub,int istart,int iend,int jstart,int jend,int numMode,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void saveCoordHDF(Domain *D,char *fileName);
void density_xdmf(char *fileName,int nx,int ny,int s);
void setZero(double ***den,int nx, int ny, int nz);
void deleteField(double ***field,int nx,int ny,int nz);
void deleteParticle(Domain *D,int i,int j,int s);
void solveCharge(Domain *D,LoadList *LL,double ***rhoR,double ***rhoI,int istart,int iend,int jstart,int jend,int s,double coef);



void saveDensityHDF(Domain *D,int iteration)
{
    int i,j,m,s,istart,iend,jstart,jend,nx,ny,nz,numMode;
    int biasX,biasY,offSetY,rankX,rankY,minRSub;
    int nxSub,nySub,iter;
    int offset[3];
    double ***dataR,***dataI,***val;
    double charge,coef;
    char name[100],dataName[100],fileName[100];
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D->nxSub;    nySub=D->nySub;
    istart=D->istart;  iend=D->iend;
    jstart=D->jstart;  jend=D->jend;
    numMode=D->numMode;
	 iter=D->filterIter;

    rankX=myrank/D->M;
    rankY=myrank%D->M;
    dataR=memoryAsign(numMode,D->nxSub+5,D->nySub+5);
    dataI=memoryAsign(numMode,D->nxSub+5,D->nySub+5);


    nx=D->nx;  ny=D->ny; nz=1;
    offset[0]=D->minXSub-D->minXDomain;
    offset[1]=D->minYSub-D->minYDomain;
    offset[2]=0;

    s=0;
    LL=D->loadList;
    while(LL->next)  
    {
      for(m=0; m<numMode; m++)
        for(i=0; i<iend+3; i++)
          for(j=0; j<jend+3; j++) {
            dataR[m][i][j]=0.0;
            dataI[m][i][j]=0.0;
          }

//      rho0[s]=LL->criticalDensity*LL->density;
      coef=LL->criticalDensity;

      sprintf(name,"%ddensity%d.h5",s,iteration);

      if(myrank==0)        {
        file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5Fclose(file_id);
        saveIntMeta(name,"/nx",&nx,1);
        saveIntMeta(name,"/ny",&ny,1);
        saveIntMeta(name,"/nz",&nz,1);
        saveIntMeta(name,"/numMode",&numMode,1);
      } else      ;
      MPI_Barrier(MPI_COMM_WORLD);

      solveCharge(D,LL,dataR,dataI,istart,iend,jstart,jend,s,coef);
      if(D->L>1)  {
        MPI_TransferDen_Xplus(D,dataR,dataI,nySub+5,3);
        MPI_TransferDen_Xminus(D,dataR,dataI,nySub+5,3);
      }  else ;

//      filter(D,dataR,dataI);
/*
      val=(double ***)malloc(1*sizeof(double ** ));
      for(m=0; m<1; m++) {
        val[m]=(double **)malloc((D->nxSub+5)*sizeof(double * ));
        for(i=0; i<D->nxSub+5; i++)
          val[m][i]=(double *)malloc((D->nySub+5)*sizeof(double  ));
      }
      for(i=0; i<D->nxSub+5; i++)
        for(j=0; j<D->nySub+5; j++)
          val[0][i][j]=0.0;
    
      filter_current(D,val,dataR,iter);
      filter_current(D,val,dataI,iter);

      for(i=0; i<D->nxSub+5; i++) free(val[0][i]);
      free(val[0]); free(val);
*/




//      solveCharge(D,dataR,dataI,rho0[s],s);

//      for(m=0; m<numMode; m++)      
//       for(i=istart; i<iend; i++)      
//          for(j=jstart; j<jend; j++)    {  
//            dataR[m][i][j]=D->RhoNoPairR[m][i][j]+D->RhoPairR[m][i][j]; 
//            dataI[m][i][j]=D->RhoNoPairI[m][i][j]+D->RhoPairI[m][i][j]; 
//          }
      sprintf(dataName,"R");
      saveFieldComp(dataR,name,dataName,nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//      saveFieldComp(D->RhoIonR,name,dataName,nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);
//      MPI_Barrier(MPI_COMM_WORLD);
//      for(m=0; m<numMode; m++)      
//        for(i=istart; i<iend; i++)      
//          for(j=jstart; j<jend; j++)      
//            data[m][i][j]=D->RhoElcI[m][i][j]*rho0[s]; 
      sprintf(dataName,"I");
      saveFieldComp(dataI,name,dataName,nx,ny,nxSub,nySub,istart,iend,jstart,jend,numMode,offset);

      if(myrank==0)   {
        saveCoordHDF(D,name);
        sprintf(fileName,"%ddensity%d",s,iteration);
//        density_xdmf(fileName,nx,ny,s);
        printf("%s\n",name);
      }  else  ;

      LL=LL->next;
      s++;
    }		//End of for(s)
    deleteField(dataR,numMode,D->nxSub+5,D->nySub+5);
    deleteField(dataI,numMode,D->nxSub+5,D->nySub+5);
}


void saveCoordHDF(Domain *D,char *fileName)
{
  int ii,i,nx,ny;
  char name[100];
  double *xtic,*ytic;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

  hid_t file_id,dset_id,tic_id;
  herr_t status;
  hid_t filespace;
  hsize_t dimy[1],dimx[1],dimz[1];
  const char *coorName[] = {"/Y","/X"};

  nx=D->nx;  ny=D->ny;

  if(myrank==0)
  {
    sprintf(name,"%s",fileName);
    file_id=H5Fopen(name,H5F_ACC_RDWR,H5P_DEFAULT);

    dimx[0]=nx;      dimy[0]=ny;
    xtic=(double *)malloc(nx*sizeof(double));
    for(i=0;i<nx;i++) xtic[i]=(i+0.5+D->minXDomain)*D->lambda*D->dz;
    ytic=(double *)malloc(ny*sizeof(double));
    for(i=0;i<ny;i++) ytic[i]=(i+0.5+D->minYDomain)*D->lambda*D->dr;
    for(ii=0; ii<2; ii++)
    {
      if(ii==0) filespace=H5Screate_simple(1,dimy,NULL);
      else if(ii==1) filespace=H5Screate_simple(1,dimx,NULL);
      dset_id=H5Dcreate2(file_id,coorName[ii],H5T_NATIVE_DOUBLE,filespace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ii==0 ? ytic : xtic);
      H5Dclose(dset_id);
      H5Sclose(filespace);
    }
    free(xtic);
    free(ytic);
    H5Fclose(file_id);
  }
  else ;
}

// The number of cells in the X, Y dimensions
void density_xdmf(char *fileName,int nx,int ny,int s)
{
    FILE *xmf = 0;
    char name[100];
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

      fprintf(xmf, "     <Attribute Name=\"/0%d\" AttributeType=\"Scalar\" Center=\"Node\">\n",s);
      fprintf(xmf, "       <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",nx);
      fprintf(xmf, "        %s.h5:/%d\n",fileName,s);
      fprintf(xmf, "       </DataItem>\n");
      fprintf(xmf, "     </Attribute>\n");
      fprintf(xmf, "   </Grid>\n");
      fprintf(xmf, " </Domain>\n");
      fprintf(xmf, "</Xdmf>\n");
      fclose(xmf);
}
                                                          
void setZero(double ***den,int nx, int ny, int nz)
{
   int i,j,k;

   for(i=0; i<nx; i++)
     for(j=0; j<ny; j++)
       for(k=0; k<nz; k++)
         den[i][j][k]=0.0;
}

void deleteParticle(Domain *D,int i,int j,int s)
{
  ptclList *p,*tmp;

  p=D->particle[i][j].head[s]->pt;
  while(p)  {
    tmp=p->next;
    D->particle[i][j].head[s]->pt=tmp;
    p->next=NULL;
    free(p);
    p=D->particle[i][j].head[s]->pt;
  }
}

