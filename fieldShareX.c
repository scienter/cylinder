#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"


void MPI_Transfer1F_Xminus(Domain *D,double ***f1,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*1*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer1F_Xplus(Domain *D,double ***f1,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*1*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer2F_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*2*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer2F_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer4F_Split_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*4*ny*(share+1);
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=-1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=-1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=-1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=-1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer4F_Split_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*4*ny*(share-2);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=2; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=2; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=2; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=2; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}


void MPI_Transfer4F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*4*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*6*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer6F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*6*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer8F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*8*ny*3;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer8F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*8*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer12F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share)
{
    int i,j,m,num,start,end,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank; 
    double *data;

    MPI_Status status;         

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*12*ny*share;
    data = (double *)malloc(num*sizeof(double ));
 
    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][i+istart][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==1)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][i+istart][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][i+istart][j]; start+=ny;
      }
        
    if(rank%2==1 && rank!=D->L-1)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank, MPI_COMM_WORLD,&status);  
      start=0; 
      for(m=0; m<numMode; m++)
        for(i=0; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][iend+i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][iend+i][j]=data[start+j]; start+=ny;
        }  
    }
    else if(rank%2==0 && rank!=0)
       MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_Transfer12F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share)
{
    int i,j,m,num,numMode;
    int istart,iend,jstart,jend;
    int myrank, nTasks,rank,start; 
    double *data;

    MPI_Status status;         
   
    istart=D->istart;   iend=D->iend;
    jstart=D->jstart;   jend=D->jend;
    numMode=D->numMode;
  
    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);            

    rank=myrank/D->M;

    num=numMode*12*ny*(share-1);
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][iend-i][j]; start+=ny;
      }
      
    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank, MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==0 && rank!=D->L-1)
       MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0; 
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++)  {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f7[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f8[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f9[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f10[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f11[m][iend-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f12[m][iend-i][j]; start+=ny;
      }
        
    if(rank%2==0 && rank!=0)
    {
      MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);  
      start=0;
      for(m=0; m<numMode; m++)
        for(i=1; i<share; i++) {
          for(j=0; j<ny; j++) f1[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f2[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f3[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f4[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f5[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f6[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f7[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f8[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f9[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f10[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f11[m][istart-i][j]=data[start+j]; start+=ny;
          for(j=0; j<ny; j++) f12[m][istart-i][j]=data[start+j]; start+=ny;
        }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank,MPI_COMM_WORLD);             
    MPI_Barrier(MPI_COMM_WORLD);
    free(data);
}

void MPI_TransferJ_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*6*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend+i][j]; start+=ny;
      }

    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][iend+i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(data);
}
                                                
void MPI_TransferJ_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*6*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][istart-i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1)
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f3[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f4[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f5[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f6[m][istart-i][j]; start+=ny;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f3[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f4[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f5[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f6[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=0)
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   
    free(data);
}
/*
void MPI_TransferRho_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
      }

    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(data);
}
                                                
void MPI_TransferRho_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1)
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=0)
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   
    free(data);
}
*/
void MPI_TransferDen_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*3;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
      }

    if(rank%2==1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++)
      for(i=0; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][iend+i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][iend+i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=0)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->prevXrank,D->prevXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=0; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][istart+i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][istart+i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1 && rank!=D->L-1)
      MPI_Send(data,num,MPI_DOUBLE,D->nextXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    free(data);
}
                                                
void MPI_TransferDen_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share)
{
    int i,j,m,start,numMode,num;
    int istart,iend,jstart,jend;
    int myrank, nTasks, rank;
    double *data;

    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    rank=myrank/D->M;
    num=numMode*2*ny*2;
    data = (double *)malloc(num*sizeof(double ));

    //Transferring even ~ odd cores 
    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
      }

    if(rank%2==0 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==1)
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    //Transferring odd ~ even cores             
    start=0;
    for(m=0; m<numMode; m++)
      for(i=1; i<share; i++) {
        for(j=0; j<ny; j++) data[start+j]=f1[m][istart-i][j]; start+=ny;
        for(j=0; j<ny; j++) data[start+j]=f2[m][istart-i][j]; start+=ny;
      }

    if(rank%2==1 && rank!=D->L-1)
    {
       MPI_Recv(data,num,MPI_DOUBLE,D->nextXrank,D->nextXrank,MPI_COMM_WORLD,&status);
       start=0;
       for(m=0; m<numMode; m++)
         for(i=1; i<share; i++) {
           for(j=0; j<ny; j++) f1[m][iend-i][j]+=data[j+start]; start+=ny;
           for(j=0; j<ny; j++) f2[m][iend-i][j]+=data[j+start]; start+=ny;
         }
    }
    else if(rank%2==0 && rank!=0)
      MPI_Send(data,num,MPI_DOUBLE,D->prevXrank,myrank, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
   
    free(data);
}
                                                
