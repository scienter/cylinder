#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void saveParticleHDF(Domain *D,int iteration,int s,double minPx,double density);
void saveDensityHDF(Domain *D,int iteration);
void saveEFieldHDF(Domain *D,int iteration);
void saveBFieldHDF(Domain *D,int iteration);
void saveSplitFieldHDF(Domain *D,int iteration);
void saveJHDF(Domain *D,int iteration);


void saveFile(Domain D,int iteration)
{
  int myrank, nTasks,s;
  double minPx,density;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  void saveField();
  void saveParticle();
  void saveDensity();
  void saveCurrent();
  LoadList *LL;


  //save current
  if(D.currentSave==ON)   {
    if(D.saveCurrentMode==TXT)  
      ;//      saveCurrent(&D,iteration);
    else if(D.saveCurrentMode==HDF)   {
      saveJHDF(&D,iteration);
    }  else	;

    if(myrank==0)
      printf("J%d is made.\n",iteration); 
    else	;
  }	else	;


  //save field
  if(D.fieldSave==ON)   {
    if(D.saveFieldMode==TXT)  {
      printf("before saveField\n");
      saveField(&D,iteration);
    }
    else if(D.saveFieldMode==HDF)   {
      switch (D.fieldType) {
      case Yee :
      case NoCherenkov :
        saveEFieldHDF(&D,iteration);
//        saveBFieldHDF(&D,iteration);
        break;
      case Split :
        saveSplitFieldHDF(&D,iteration);
        break;
      }
    }  else	;

    if(myrank==0)
      printf("field%d is made.\n",iteration); 
    else	;
  }	else	;

  //save particle
  if(D.particleSave==ON)
  {
    LL=D.loadList;
    s=0;
    while(LL->next)  {
      minPx=LL->givenMinPx;
      density=LL->density;
      if(D.saveParticleMode==TXT) { 
        saveParticle(&D,iteration);
      }
      else if(D.saveParticleMode==HDF) 
        saveParticleHDF(&D,iteration,s,minPx,density);
      
      else	;

      LL=LL->next;
      s++;
    }

    if(myrank==0)
      printf("particle%d is made.\n",iteration);  
    else	;
  }	else	;

  //save density
  if(D.densitySave==ON)
  {
    if(D.saveDensityMode==TXT)  ;//  saveDensity(&D,iteration);
    else if(D.saveDensityMode==HDF)
      saveDensityHDF(&D,iteration);
    else ;
    if(myrank==0)
      printf("density%d is made.\n",iteration); 
  }  else	;
 
  //save current
  if(D.currentSave==ON)
  {
    if(D.saveCurrentMode==TXT)
//      saveCurrent(&D,iteration);
//    else if(D.saveCurrentMode==HDF)
//      saveCurrentHDF(&D,iteration);
    if(myrank==0)
      printf("current%d is made.\n",iteration); 
  }  else	;
}

/*
void saveProbe(Domain *D,int iteration)
{
    int i,j,n;
    char name[100];
    double t,Ex,Ey,Ez,Bx,By,Bz,Pr,Pl,Sr,Sl,x,y;
    double omega,frequency,dt;
    FILE *out;
    int myrank, nprocs;    
    Probe **probe;
    probe=D->probe;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    omega=2*pi*velocityC/D->lambda;
    frequency=omega/2.0/pi;   
    dt=1.0/frequency*D->dt;

    for(n=0; n<D->probeNum; n++)
    {
      if(D->probeX[n]>=D->minXSub && D->probeX[n]<D->maxXSub && 
         D->probeY[n]>=D->minYSub && D->probeY[n]<D->maxYSub)
      {
        x=D->probeX[n]*D->dx*D->lambda;
        y=D->probeY[n]*D->dy*D->lambda;
        sprintf(name,"probeRaman%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Pr=probe[n][i].Pr;
          Pl=probe[n][i].Pl;
          Sr=probe[n][i].Sr;
          Sl=probe[n][i].Sl;
          fprintf(out,"%g %g %g %g %g %g %g\n",t,Pr,Pl,Sr,Sl,x,y);
        }
        fclose(out);

        sprintf(name,"probe%d_%d",iteration,n);
        out = fopen(name,"w");
        for(i=0; i<=iteration; i++)
        {
          t=i*dt;
          Ex=probe[n][i].E1;
          Bx=probe[n][i].B1;
          Ey=probe[n][i].Pr+probe[n][i].Pl;
          Ez=probe[n][i].Sr+probe[n][i].Sl;
          By=probe[n][i].Sl-probe[n][i].Sr;
          Bz=probe[n][i].Pr-probe[n][i].Pl;
          fprintf(out,"%g %g %g %g %g %g %g %g %g\n",t,Ex,Ey,Ez,Bx,By,Bz,x,y);
        }             
        fclose(out);
      }
    }
}
*/


/*
void saveDensity(Domain *D,int iteration)
{
  int i,j,s,index,m,l,istart,iend,jstart,jend,angle[2];
  int ii,jj,minZSub,minRSub,numMode;
  double x,y,lambda,dz,dr,cosP[4],sinP[4],rho,tmp1,tmp2;
  double R,inv,factor,z,r,phi,Wz[2],Wr[2],weight,charge;
  char name[100];
  FILE *out;
  Particle **particle;
  particle=D->particle;
  ptclList *p;
  int myrank, nTasks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;  iend=D->iend;
  jstart=D->jstart;  jend=D->jend;
  minZSub=D->minXSub; minRSub=D->minYSub;
  lambda=D->lambda;  dr=D->dr; dz=D->dz;
  numMode=D->numMode;

  for(i=istart; i<iend; i++)
    for(j=jstart; j<jend; j++)
      D->Rho[0][i][j]=0.0;

  s=0;
  charge=1.0;
    for(i=istart; i<iend; i++)
      for(j=jstart+1; j<jend; j++)
      {
        p=particle[i][j].head[s]->pt;
        while(p) {
          weight=p->weight;
          z=p->z;  r=p->r;  phi=p->phi;
          Wz[1]=z;                            Wz[0]=1.0-Wz[1];
          R=j+r-jstart+minRSub-0.5;
          index=j-jstart+minRSub;
          inv=(double)(2.0*index);
          Wr[1]=(R*R-(index-0.5)*(index-0.5))/inv;  Wr[0]=1.0-Wr[1];
          factor=weight*charge/(4.0*index);
          for(ii=0; ii<2; ii++)
            for(jj=0; jj<2; jj++)
              D->Rho[0][i+ii][j+jj]+=Wr[jj]*Wz[ii]*factor;
          p=p->next;
        }
      }           
    sprintf(name,"density%d_%d",iteration,myrank);
    out = fopen(name,"w");    
    for(i=istart; i<iend; i++)
    {
      x=(i-istart+minZSub)*dz*lambda;
      for(j=jstart+1; j<jend; j++)
      {
        y=(j-jstart+minRSub-0.5)*dr*lambda;

        fprintf(out,"%g %g %g\n",x,y,D->Rho[0][i][j]);    
      }           
      fprintf(out,"\n");    
    }
    fclose(out);
  
}
*/

/*
void boostSaveField(Domain *D,int labSaveStep)
{
    int i,j,istart,iend,jstart,jend,show;
    char name[100];
    double x,y,e1,pr,pl,b1,sr,sl;
    double factor,dx;
    FILE *out;
    int myrank, nprocs;    
    Boost **boost;
    boost=D->boost;

    istart=D->istart;
    iend=D->iend;
    jstart=D->jstart;
    jend=D->jend;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    sprintf(name,"bField%d_%d",labSaveStep,myrank);
    out = fopen(name,"w");

    for(i=istart; i<iend; i++)
    {
      for(j=jstart; j<jend; j++)
      {
        x=boost[i][j].x;
        y=boost[i][j].y;
        e1=boost[i][j].E1;
        pr=boost[i][j].Pr;
        pl=boost[i][j].Pl;
        b1=boost[i][j].B1;
        sr=boost[i][j].Sr;
        sl=boost[i][j].Sl;
        if(x>0) {
          show=1;
          fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,e1,pr,pl,b1,sr,sl);
        }
        else show=0;
      }  
      if(show==1)           
        fprintf(out,"\n");
    }
    fclose(out);
    
    if(myrank==0)
      printf("bField%d is saved.\n",labSaveStep);
}
*/


void saveCenterField(Domain *D,int iteration)
{
   int myrank;    
   MPI_Status status;
   void calculationCenter();
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   switch (D->fieldType) {
     case Yee :
       calculationCenter(D,D->EzR,D->ErR,D->EpR,D->BzR,D->BrR,D->BpR,D->RhoNoPairR,"",iteration);
     break;
     case Split :
       calculationCenter(D,D->EzR,D->PrR,D->PlR,D->BzR,D->SrR,D->SlR,D->RhoNoPairR,"",iteration);
     break;
   }
}


void calculationCenter(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,char *fieldType,int iteration)
{
  int i,j,m,n,rankX,targetCore,fromCore,index,nx,minXSub,minXDomain,numMode,pick;
  char name[100];
  double x,Ex,Ey,Ez,Bx,By,Bz,den,nc,factor,*field,*share;
  double ff1,ff2,ff3,ff4,ff5,ff6,ff7,angle;
  double cosP[D->numMode],sinP[D->numMode];
  FILE *out;
  int myrank,number,istart,iend,jstart,jend;    
  MPI_Status status;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  numMode=D->numMode;
  nx=D->nx;   minXSub=D->minXSub;   minXDomain=D->minXDomain;
  nc=D->loadList->criticalDensity;

  number=7*nx*numMode;
  field=(double *)malloc(number*sizeof(double ));
  share=(double *)malloc(number*sizeof(double ));
  for(i=0; i<number; i++)  {
     share[i]=0.0;     field[i]=0.0;
  }
  pick=(int)(D->centerPick/D->dr);

  angle=D->centerAngle*pi/180.0;
  cosP[0]=1.0; sinP[0]=0.0;
  for(m=1; m<numMode; m++) { 
    cosP[m]=cos(m*angle); sinP[m]=sin(m*angle);
  }

  if(D->minYSub<=0 && D->maxYSub>=0)
  {
    targetCore=myrank%D->M;
    rankX=myrank/D->M;

    j=jstart+pick;
    for(i=istart; i<iend; i++)  {
      index=i-istart+minXSub-minXDomain;
      ff1=0.0; ff2=0.0; ff3=0.0;  ff4=0.0;
      ff5=0.0; ff6=0.0; ff7=0.0;
      for(m=0; m<numMode; m++)  {
        ff1+=f1[m][i][j]*cosP[m];
        ff1-=f1[m][i][j]*sinP[m];
        ff2+=f2[m][i][j]*cosP[m];
        ff2-=f2[m][i][j]*sinP[m];
        ff3+=f3[m][i][j]*cosP[m];
        ff3-=f3[m][i][j]*sinP[m];
        ff4+=f4[m][i][j]*cosP[m];
        ff4-=f4[m][i][j]*sinP[m];
        ff5+=f5[m][i][j]*cosP[m];
        ff5-=f5[m][i][j]*sinP[m];
        ff6+=f6[m][i][j]*cosP[m];
        ff6-=f6[m][i][j]*sinP[m];
        ff7+=f7[m][i][j]*cosP[m];
        ff7-=f7[m][i][j]*sinP[m];
      }
      share[nx*0+index]=ff1;
      share[nx*1+index]=ff2;
      share[nx*2+index]=ff3;
      share[nx*3+index]=ff4;
      share[nx*4+index]=ff5;
      share[nx*5+index]=ff6;
      share[nx*6+index]=ff7;
      field[nx*0+index]=share[nx*0+index];    
      field[nx*1+index]=share[nx*1+index];    
      field[nx*2+index]=share[nx*2+index];    
      field[nx*3+index]=share[nx*3+index];    
      field[nx*4+index]=share[nx*4+index];    
      field[nx*5+index]=share[nx*5+index];    
      field[nx*6+index]=share[nx*6+index];    
    }

    for(n=1; n<D->L; n++)  {
      fromCore=targetCore+n*D->M;
      if(myrank==fromCore)
        MPI_Send(share,number,MPI_DOUBLE,targetCore,myrank,MPI_COMM_WORLD);
      else	;
    } 

    if(myrank==targetCore) {
      for(n=1; n<D->L; n++)  {
        fromCore=targetCore+n*D->M;
        MPI_Recv(share,number,MPI_DOUBLE,fromCore,fromCore,MPI_COMM_WORLD,&status);
        for(i=0; i<number; i++)
          field[i]+=share[i];    
      }
    } else	;  //End of if(myrank=targetCore)

    if(myrank==targetCore)       {
      sprintf(name,"cen%d_%d",iteration,myrank);
      out = fopen(name,"w");
      for(i=0; i<nx; i++)  {
        x=(i+minXDomain)*D->dz*D->lambda;
        Ex=field[nx*0+i]; Ey=field[nx*1+i]; Ez=field[nx*2+i];
        Bx=field[nx*3+i]; By=field[nx*4+i]; Bz=field[nx*5+i]; den=field[nx*6+i];
        fprintf(out,"%.7g %g %g %g %g %g %g %g\n",x,Ex,Ey,Ez,Bx,By,Bz,den*nc);
      }
      fclose(out);
    }	else	;

  }   else	;  //End of minYSub<0<maxYSub
  free(share);
  free(field);
}
/*
void saveCurrent(Domain *D,int iteration)
{
    int i,j,k,istart,iend,jstart,jend,kstart,kend;
    char name[100];
    double x,y,z,Jx,Jy,Jz;
    FILE *out;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    kstart=D->kstart;    kend=D->kend;

    sprintf(name,"current%d_%d",iteration,myrank);
    out = fopen(name,"w");

    switch(D->dimension) {

    case 1:
      j=k=0;
      for(i=istart; i<iend; i++)      {
        x=(i-2+D->minXSub)*D->dx*D->lambda;
        Jx=D->Jx[i][j][k]; Jy=D->Jy[i][j][k]; Jz=D->Jz[i][j][k];
        fprintf(out,"%g %g %g %g\n",x,Jx,Jy,Jz);
      }
      fclose(out);
      break;

    case 2:
      k=0;
      for(i=istart; i<iend; i++)    {
        for(j=jstart; j<jend; j++)        {
          x=(i-2+D->minXSub)*D->dx*D->lambda;
          y=(j-2+D->minYSub)*D->dy*D->lambda;
          Jx=D->Jx[i][j][k]; Jy=D->Jy[i][j][k]; Jz=D->Jz[i][j][k];
          fprintf(out,"%g %g %g %g %g\n",x,y,Jx,Jy,Jz);
        }
        fprintf(out,"\n");                 
      }
      fclose(out);
      break;
    default :
      printf("In saveCurrent, what dimension?\n");
    }
}
*/

void saveField(Domain *D,int iteration)
{
    int i,j,m,istart,iend,jstart,jend,numMode;
    char name[100];
    double x,y,Ez,Er,Ep,Bz,Br,Bp,phi,factor,cos2,sin2;
    double cosP,sinP,coss[D->numMode],sins[D->numMode];
    FILE *out,*out1;
    int myrank, nprocs;    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    sprintf(name,"field%d_%d",iteration,myrank);
    out = fopen(name,"w");
    factor=D->gamma*(1+D->beta);

    phi=0.0;
    cosP=cos(phi*pi/180.0);
    sinP=sin(phi*pi/180.0);
    coss[0]=1.0; sins[0]=0.0;
    for(m=1; m<numMode; m++) {
      cos2=coss[m-1]; sin2=sins[m-1];
      coss[m]=cosP*cos2-sinP*sin2;
      sins[m]=sinP*cos2+cosP*sin2;
    }

    for(i=istart; i<iend; i++)  {    
      for(j=jstart; j<jend; j++)  {
        x=(i-2+D->minXSub)*D->dz*D->lambda;
        y=(j-2+D->minYSub)*D->dr*D->lambda;

        Ez=D->EzR[0][i][j];
        Er=D->ErR[0][i][j];
        Ep=D->EpR[0][i][j];
        Bz=D->BzR[0][i][j];
        Br=D->BrR[0][i][j];
        Bp=D->BpR[0][i][j];

        for(m=0; m<numMode; m++)  {
          Ez+=D->EzR[m][i][j]*coss[m]-D->EzI[m][i][j]*sins[m];
          Er+=D->ErR[m][i][j]*coss[m]-D->ErI[m][i][j]*sins[m];
          Ep+=D->EpR[m][i][j]*coss[m]-D->EpI[m][i][j]*sins[m];
          Bz+=D->BzR[m][i][j]*coss[m]-D->BzI[m][i][j]*sins[m];
          Br+=D->BrR[m][i][j]*coss[m]-D->BrI[m][i][j]*sins[m];
          Bp+=D->BpR[m][i][j]*coss[m]-D->BpI[m][i][j]*sins[m];
        }
        fprintf(out,"%g %g %g %g %g %g %g %g\n",x,y,Ez,Er,Ep,Bz,Br,Bp);
      }
      fprintf(out,"\n");
    }
    fclose(out);
}


void saveParticle(Domain *D,int iteration)
{
  int i,j,istart,iend,jstart,jend,s,core,index;
  int minZSub,minRSub;
  char name[100];
  double x,y,z,r,pz,px,py,gamma,mc,weight,cosP,sinP,phi;
  double minPz[D->nSpecies],dr,dz,lambda;
  Particle **particle;
  particle=D->particle;
  ptclList *p;
  LoadList *LL;
  FILE *out;
  int myrank, nprocs;    
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  istart=D->istart;  iend=D->iend;
  jstart=D->jstart;  jend=D->jend;
  dz=D->dz; dr=D->dr; lambda=D->lambda;
  minZSub=D->minXSub;
  minRSub=D->minYSub;

  LL=D->loadList;
  s=0;
  while(LL->next)
  {
    minPz[s]=LL->givenMinPx;
    LL=LL->next;
    s++;
  }

    for(s=0; s<D->nSpecies; s++)
    {
      sprintf(name,"%dParticle%d_%d",s,iteration,myrank);
      out = fopen(name,"w");    
      for(i=istart; i<iend; i++)
        for(j=jstart; j<jend; j++)
        {
          p=particle[i][j].head[s]->pt;
          while(p)
          {
            z=(p->z+i-istart+minZSub)*dz*lambda; 
            x=p->x;
            y=p->y;
            pz=p->pz;
            px=p->px;    
            py=p->py;
//            pz=p->Ez;
//            pr=p->Er;    
//            pphi=p->Ep;
            index=p->index;
            core=p->core;
            weight=p->weight;
//            if(pz>=minPz[s])
              fprintf(out,"%g %g %g %g %g %g %d %d %g\n",z,x,y,pz,px,py,index,core,weight);               
//              fprintf(out,"%g %g %g %g %g %g %g %g %g\n",z,x,y,p->Ez,p->Er,p->Ep,p->Bz,p->Br,p->Bp);               
            p=p->next;
          }	//End of while(p)
        }	//End of for(i,j)
      fclose(out);
    }				//End of for(s)
}

