#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "constants.h"
#include "mpi.h"
#include "hdf5.h"
#include "hdf5_hl.h"

void MPI_TransferJ_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share);
void MPI_TransferJ_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,int ny,int nz,int share);
void MPI_TransferJ_Yplus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share);
void MPI_TransferJ_Yminus(Domain *D,double ***f1,double ***f2,double ***f3,int nx,int nz,int share);
void MPI_TransferRho_Xplus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Xminus(Domain *D,double ***f1,int ny,int nz,int share);
void MPI_TransferRho_Yminus(Domain *D,double ***f1,int nx,int nz,int share);
void MPI_TransferRho_Yplus(Domain *D,double ***f1,int nx,int nz,int share);

void calParameter(int nx,int *istart,int *iend,int *nxSub,int rankX,int *biasX,int L);
void saveFieldComp(double ***data,char *fileName,char *dataName,int nx,int ny,int nz,int nxSub,int nySub,int nzSub,int istart,int iend,int jstart,int jend,int kstart,int kend,int *offSet);
void saveIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void solvePMap2D(Domain *D,double ***Px,double ***Py,double ***Pz,double ***den,int s);
//void solvePMap3D(Domain *D,int s);
double ***memoryAsign(int nx, int ny, int nz);
void deleteField(double ***field,int nx,int ny,int nz);



void saveDumpPMapResolHDF(Domain D,int iteration)
{
    int i,j,k,s,istart,iend,jstart,jend,kstart,kend,nx,ny,nz;
    int biasX,biasY,biasZ,offSetY,rankX,rankY,rankZ;
    int nxSub,nySub,nzSub;
    int offset[3];
    double ***Px,***Py,***Pz,***den;
    double x,y,z,px,py,pz;
    FILE *out;
    char name[100],dataName[100];
    LoadList *LL;

    int myrank, nTasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

    hid_t file_id,group_id,dset_id,filespace;
    herr_t status;

    nxSub=D.nxSub;
    nySub=D.nySub;
    nzSub=D.nzSub;
    istart=D.istart;
    iend=nxSub+2;
    jstart=D.jstart;
    jend=nySub+2;
    kstart=D.kstart;
    kend=nzSub+2;

    sprintf(name,"dumpPMapResol%d.h5",iteration);
    if(myrank==0)
    {
      file_id=H5Fcreate(name,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      H5Fclose(file_id);
    }
    else ;

    rankX=myrank/(D.M*D.N);
    rankZ=(myrank%(D.M*D.N))/D.M;
    rankY=(myrank%(D.M*D.N))%D.M;

    if(myrank==0)  {
      saveIntMeta(name,"/minXDomain",&D.minXDomain,1);
      saveIntMeta(name,"/minYDomain",&D.minYDomain,1);
      saveIntMeta(name,"/minZDomain",&D.minZDomain,1);
      saveIntMeta(name,"/nSpecies",&D.nSpecies,1);
      saveIntMeta(name,"/nx",&D.nx,1);
      saveIntMeta(name,"/ny",&D.ny,1);
      saveIntMeta(name,"/nz",&D.nz,1);
    } else      ;

    switch(D.dimension) {
    //2D
    case 2:
      nx=D.nx+5;
      calParameter(nx,&istart,&iend,&nxSub,rankX,&biasX,D.L);
      ny=D.ny+5;
      calParameter(ny,&jstart,&jend,&nySub,rankY,&biasY,D.M);

      offset[0]=(D.minXSub-D.minXDomain)+biasX;
      offset[1]=(D.minYSub-D.minYDomain)+biasY;
      offset[2]=0;

      Px=memoryAsign(D.nxSub+5,D.nySub+5,1);
      Py=memoryAsign(D.nxSub+5,D.nySub+5,1);
      Pz=memoryAsign(D.nxSub+5,D.nySub+5,1);
      den=memoryAsign(D.nxSub+5,D.nySub+5,1);
      for(s=0; s<D.nSpecies; s++)
      {
        solvePMap2D(&D,Px,Py,Pz,den,s);
        if(D.L>1)  {
          MPI_TransferJ_Xplus(&D,Px,Py,Pz,D.nySub+5,1,3);
          MPI_TransferJ_Xminus(&D,Px,Py,Pz,D.nySub+5,1,3);
          MPI_TransferRho_Xplus(&D,den,D.nySub+5,1,3);
          MPI_TransferRho_Xminus(&D,den,D.nySub+5,1,3);
        }  else     ;
        if(D.M>1)  {
          MPI_TransferJ_Yplus(&D,Px,Py,Pz,D.nxSub+5,1,3);
          MPI_TransferJ_Yminus(&D,Px,Py,Pz,D.nxSub+5,1,3);
          MPI_TransferRho_Yplus(&D,den,D.nxSub+5,1,3);
          MPI_TransferRho_Yminus(&D,den,D.nxSub+5,1,3);
        }  else     ;

        for(i=D.istart; i<D.iend; i++)
          for(j=D.jstart; j<D.jend; j++)  {   
            if(den[i][j][0]==0)  {            
              Px[i][j][0]=den[i][j][0];
              Py[i][j][0]=den[i][j][0];
              Pz[i][j][0]=den[i][j][0];
            }
            else  {            
              Px[i][j][0]/=den[i][j][0];
              Py[i][j][0]/=den[i][j][0];
              Pz[i][j][0]/=den[i][j][0];
            }
          }
        sprintf(dataName,"%dPx",s);
        saveFieldComp(Px,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%dPy",s);
        saveFieldComp(Py,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
        sprintf(dataName,"%dPz",s);
        saveFieldComp(Pz,name,dataName,nx,ny,1,nxSub,nySub,1,istart,iend,jstart,jend,0,1,offset);
      }

      if(myrank==0)   printf("%s\n",name);  else  ;

    sprintf(name,"pMap%d_%d",iteration,myrank);
    out = fopen(name,"w");
    k=0;
    for(i=D.istart; i<D.iend; i++)
    {
      for(j=D.jstart; j<D.jend; j++) 
      {
        x=(i-2+D.minXSub)*D.dx*D.lambda;
        y=(j-2+D.minYSub)*D.dy*D.lambda;
        px=Px[i][j][k];
        py=Py[i][j][k];
        pz=Pz[i][j][k];
        fprintf(out,"%g %g %g %g %g\n",x,y,px,py,pz);
      }
      fprintf(out,"\n");
    }
    fclose(out);
//    printf("%s is made\n");
      
      
      deleteField(Px,D.nxSub+5,D.nySub+5,1);
      deleteField(Py,D.nxSub+5,D.nySub+5,1);
      deleteField(Pz,D.nxSub+5,D.nySub+5,1);
      deleteField(den,D.nxSub+5,D.nySub+5,1);
      break;

    }   //End of switch (dimension)

}


void solvePMap2D(Domain *D,double ***Px,double ***Py,double ***Pz,double ***den,int s)
{
  int i,j,k,ii,jj,istart,iend,jstart,jend;
  double Wx[4],Wy[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;

  k=0;
  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)  {
      Px[i][j][k]=0.0;
      Py[i][j][k]=0.0;
      Pz[i][j][k]=0.0;
      den[i][j][k]=0.0;
    }

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
    {
      p=particle[i][j][k].head[s]->pt;
      while(p)
      {
        x=p->x;
        x1=1+x;
        x2=x;
        x3=1-x;
        x4=2-x;
        Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
        Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
        Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
        Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
        y=p->y;
        y1=1+y;
        y2=y;
        y3=1-y;
        y4=2-y;
        Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
        Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
        Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
        Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
        for(ii=0; ii<4; ii++)
          for(jj=0; jj<4; jj++)  {
            Px[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj]*p->p1;
            Py[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj]*p->p2;
            Pz[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj]*p->p3;
            den[i-1+ii][j-1+jj][k]+=Wx[ii]*Wy[jj];
          }
        p=p->next;
      }
    }
}   

void solvePMap3D(Domain *D,int s,double coef)
{
  int i,j,k,ii,jj,kk,istart,iend,jstart,jend,kstart,kend;
  double Wx[4],Wy[4],Wz[4];
  double x,x1,x2,x3,x4,y,y1,y2,y3,y4,z,z1,z2,z3,z4;
  Particle ***particle;
  particle=D->particle;
  ptclList *p;

  istart=D->istart;
  iend=D->iend;
  jstart=D->jstart;
  jend=D->jend;
  kstart=D->kstart;
  kend=D->kend;

  for(i=0; i<iend+3; i++)
    for(j=0; j<jend+3; j++)
      for(k=0; k<kend+3; k++)
        D->Rho[i][j][k]=0.0;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      for(k=kstart; k<kend; k++)
      {
        p=particle[i][j][k].head[s]->pt;
        while(p)
        {
          x=p->x;
          x1=1+x;
          x2=x;
          x3=1-x;
          x4=2-x;
          Wx[0]=(2-x1)*(2-x1)*(2-x1)/6.0;
          Wx[1]=(4-6*x2*x2+3*x2*x2*x2)/6.0;
          Wx[2]=(4-6*x3*x3+3*x3*x3*x3)/6.0;
          Wx[3]=(2-x4)*(2-x4)*(2-x4)/6.0;
          y=p->y;
          y1=1+y;
          y2=y;
          y3=1-y;
          y4=2-y;
          Wy[0]=(2-y1)*(2-y1)*(2-y1)/6.0;
          Wy[1]=(4-6*y2*y2+3*y2*y2*y2)/6.0;
          Wy[2]=(4-6*y3*y3+3*y3*y3*y3)/6.0;
          Wy[3]=(2-y4)*(2-y4)*(2-y4)/6.0;
          z=p->z;
          z1=1+z;
          z2=z;
          z3=1-z;
          z4=2-z;
          Wz[0]=(2-z1)*(2-z1)*(2-z1)/6.0;
          Wz[1]=(4-6*z2*z2+3*z2*z2*z2)/6.0;
          Wz[2]=(4-6*z3*z3+3*z3*z3*z3)/6.0;
          Wz[3]=(2-z4)*(2-z4)*(2-z4)/6.0;
          for(ii=0; ii<4; ii++)
            for(jj=0; jj<4; jj++)
              for(kk=0; kk<4; kk++)
                D->Rho[i-1+ii][j-1+jj][k-1+kk]+=Wx[ii]*Wy[jj]*Wz[kk]*coef;
          p=p->next;
        }
      }
}  

