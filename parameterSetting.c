#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>
#include <time.h>

int whatONOFF(char *str);
double randomV();
int FindParameters (char *block, int rank, char *options, char *input, char *ret);
int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input);
int findLaserParameters(int rank, LaserList *L,Domain *D,char *input);
int whatSaveMode(char *str);
int whatFieldType(char *str);
int whatPlasmaType(char *str);
int whatSpecies(char *str);
double whatMass(int species);
int whatCharge(int species);
int whatFunctionMode(char *str);
int whatDefineMode(char *str);
int whatCurrentCons(char *str);
int whatLoadingMethod(char *str);
int whatLaserMode(char *str);

void parameterSetting(Domain *D,External *Ext, char *input)
{
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   FILE *in=NULL;
   double minX,maxX,minY,maxY,minZ,maxZ;
   double x,y,z,px,py,pz,gamma,vg;
   double positionX,factor,pMinX,pMaxX,pPosition;
   double normalB,normalE,Ex,Ey,Ez,Bx,By,Bz,droverdz,dzoverdx;
   char str[100],name[100],fileName[100];
   int rank,minT,maxT,tmpInt,fail=0,cnt;
   int i,j,k,n,numProbeX,numProbeY,numProbeZ,probeType,id,core,species;
   double lambda,tmpDouble,probeDx,probeDy,probeDz,maxProbeX,minProbeX,maxProbeY,minProbeY,maxProbeZ,minProbeZ,tmp;
   LoadList *LL, *New;
   LaserList *L, *LNew;


   //initially
   if(FindParameters("Domain",1,"L",input,str)) D->L=atoi(str);
   else  {
      printf("in [Domain], L=?  (Sorry. Fix as L=1)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"M",input,str)) D->M=atoi(str);
   else  {
      printf("in [Domain], M=?  (y directionally dividing number)\n");
      fail=1;
   }
   if(D->M*D->L!=nTasks)  {
     printf("L=%d, M=%d, check the values of L and M.\n",D->L,D->M);
     fail=1;
   }

   //low pass filter
   if(FindParameters("Domain",1,"filter",input,str)) D->filter=whatONOFF(str); else D->filter=OFF;
   if(FindParameters("Domain",1,"filter_iteration",input,str)) D->filterIter=atoi(str); else D->filterIter=1e5;
   if(D->filterIter==0) D->filterIter=1e5; else ;

   //Field Type
   if(FindParameters("Domain",1,"field_type",input,str)) D->fieldType=whatFieldType(str);
   else { printf("in [Domain], field_type=?  (Split/Yee/NoCherenkov)\n"); fail=1; }
   if(FindParameters("Domain",1,"center_compensation",input,str)) D->cenComp=atof(str);
   else D->cenComp=0.95;
   if(FindParameters("Domain",1,"mode_number",input,str)) 
     D->numMode=atoi(str);
   else  {
      printf("in [Domain], mode_number=?  (m=?)\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"dF",input,str)) D->dF=atof(str);
   else  D->dF=0.0;
   //Current Type
   if(FindParameters("Domain",1,"current_order",input,str)) D->currentType=atoi(str);
   else  
      D->currentType=1;
   if(FindParameters("Domain",1,"interpolation_order",input,str)) D->interpolationType=atoi(str);
   else 
      D->interpolationType=1;
   if(FindParameters("Domain",1,"current_conservation",input,str)) D->currentCons=whatCurrentCons(str);
   else D->currentCons=Lifschitz;

   //Boost frame
   if(FindParameters("Domain",1,"boost_gamma",input,str)) D->gamma=atof(str);
   else D->gamma=1;
   if(FindParameters("Domain",1,"boost_ion",input,str)) D->boostIon=whatONOFF(str);
   else D->boostIon=ON;
   if(D->gamma>1)   D->boostOn=ON;
   else             D->boostOn=OFF;
   D->beta=sqrt(1-1.0/D->gamma/D->gamma);

   //Domain parameter setting
   if(FindParameters("Domain",1,"max_time",input,str)) D->maxTime=atoi(str);
   else  D->maxTime=525600;
   if(FindParameters("Domain",1,"max_step",input,str)) D->maxStep=atoi(str);
   else  {
      printf("In [Domain], maxStep=? \n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_step",input,str)) D->saveStep=atoi(str);
   else  {
      printf("In [Domain], save_step=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"save_start",input,str)) D->saveStart=atoi(str);
   else  {
      printf("In [Domain], save_start=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"center_save_step",input,str)) D->centerStep=atoi(str);
   else  {
      printf("In [Domain], center_save_step=?\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"center_pick",input,str)) D->centerPick=atof(str);
   else  D->centerPick=0e-6; 
   if(FindParameters("Domain",1,"center_angle",input,str)) D->centerAngle=atof(str);
   else  D->centerAngle=90; 
   //save options
   if(FindParameters("Save",1,"field_format",input,str)) 
     D->saveFieldMode=whatSaveMode(str);
   else  
     D->saveFieldMode=TXT;
   if(FindParameters("Save",1,"current_format",input,str)) 
     D->saveCurrentMode=whatSaveMode(str);
   else  
     D->saveCurrentMode=TXT;
   if(FindParameters("Save",1,"particle_format",input,str)) 
     D->saveParticleMode=whatSaveMode(str);
   else  
     D->saveParticleMode=TXT;
   if(FindParameters("Save",1,"density_format",input,str)) 
     D->saveDensityMode=whatSaveMode(str);
   else  
     D->saveDensityMode=TXT;
   if(FindParameters("Save",1,"dump_format",input,str)) 
     D->saveDumpMode=whatSaveMode(str);
   else  
     D->saveDumpMode=TXT;
   if(FindParameters("Save",1,"dump_save",input,str)) 
     D->dumpSave=whatONOFF(str);
   else  
     D->dumpSave=OFF;
   if(FindParameters("Save",1,"dump_start",input,str)) 
     D->dumpStart=atoi(str);
   else  
     D->dumpStart=D->saveStart;
   if(FindParameters("Save",1,"dump_save_step",input,str))
     D->dumpSaveStep=atoi(str);
   else D->dumpSaveStep=D->saveStep;
   D->dumpStep=0;
   if(FindParameters("Save",1,"field_save",input,str)) 
     D->fieldSave=whatONOFF(str);
   else  
     D->fieldSave=ON;
   if(FindParameters("Save",1,"current_save",input,str)) 
     D->currentSave=whatONOFF(str);
   else  
     D->currentSave=ON;
   if(FindParameters("Save",1,"particle_save",input,str)) 
     D->particleSave=whatONOFF(str);
   else  
     D->particleSave=ON;
   if(FindParameters("Save",1,"density_save",input,str)) 
     D->densitySave=whatONOFF(str);
   else  
     D->densitySave=ON;

   //particle redistributing option
   if(FindParameters("Domain",1,"particle_redist",input,str)) 
      D->redist=whatONOFF(str);
   else  D->redist=OFF;
   if(FindParameters("Domain",1,"particle_redist_step",input,str)) 
     D->redistStep=atoi(str);
   else  
     D->redistStep=100000000;

   //resolution Option
   if(FindParameters("Save",1,"resolution_change",input,str)) 
      D->resolChange=whatONOFF(str);
   else  D->resolChange=OFF;
   if(FindParameters("Save",1,"resolution_high",input,str)) 
      D->resolHigh=whatONOFF(str);
   else  D->resolHigh=OFF;
   if(FindParameters("Save",1,"resolution_low",input,str)) 
      D->resolLow=whatONOFF(str);
   else  D->resolLow=OFF;
   if(FindParameters("Save",1,"resolution_change_step",input,str)) 
      D->resolStep=atoi(str);
   else  D->resolStep=D->maxTime;
   if(FindParameters("Save",1,"resolution_rate_X",input,str)) 
      D->resolX=atoi(str);
   else  D->resolX=1;
   if(FindParameters("Save",1,"resolution_rate_Y",input,str)) 
      D->resolY=atoi(str);
   else  D->resolY=1;
   if(FindParameters("Domain",1,"minX",input,str)) 
   {
      minX=atof(str);
      minX*=D->gamma*(1+D->beta);
   }
   else  {
      printf("In [Domain], minX=? [m].\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"maxX",input,str)) 
   {
      maxX=atof(str);
      maxX*=D->gamma*(1+D->beta);
   }   else  {
      printf("In [Domain], maxX=? [m].\n");
      fail=1;
   }

   if(FindParameters("Domain",1,"minY",input,str)) 
     minY=atof(str);
   else  {
     printf("In [Domain], minY=? [m].\n");
     fail=1;
   }
   if(FindParameters("Domain",1,"maxY",input,str)) 
     maxY=atof(str);
   else  {
     printf("In [Domain], maxY=? [m].\n");
     fail=1;
   }
   if(FindParameters("Domain",1,"dr_over_dz",input,str)) 
     droverdz=atof(str);
   else  {
     printf("In [Domain], dr_over_dz=?  [dr/dz]\n");
     fail=1;
   }
   if(FindParameters("Domain",1,"lambda",input,str)) 
   {
      D->lambda=atof(str);
      D->lambda*=D->gamma*(1+D->beta);
   }
   else  {
      printf("In [Domain], lambda=? [m]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"division_lambda",input,str)) 
      D->divisionLambda=atof(str);
   else  {
      printf("In [Domain], divisionLambda=? [1/number of devided wavelength]\n");
      fail=1;
   }
	D->divisionLambda=1.0/D->divisionLambda;
   if(FindParameters("Domain",1,"dt_ratio",input,str)) 
      D->dtRatio=atof(str);
   else  {
      printf("In [Domain], dt_ratio=? [<1.0]\n");
      fail=1;
   }
   if(FindParameters("Domain",1,"moving_domain",input,str)) D->moving=whatONOFF(str);
   else  {
      printf("In [Domain], moving_domain=? [ON/OFF].\n");
      fail=1;
   }   
   if(FindParameters("Domain",1,"move_now",input,str)) D->moveIt=whatONOFF(str);
   else D->moveIt=OFF;
   if(FindParameters("Domain",1,"moving_velocity",input,str)) 
     D->movingV=atof(str);
   else D->movingV=D->dtRatio;

/*
   //External field parameter setting
   if(FindParameters("External",1,"Ex",input,str)) Ex=atof(str);
   else  {
      printf("In [External], Ex=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Ey",input,str)) Ey=atof(str);
   else  {
      printf("In [External], Ey=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Ez",input,str)) Ez=atof(str);
   else  {
      printf("In [External], Ez=? [V/m]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Bx",input,str)) Bx=atof(str);
   else  {
      printf("In [External], Bx=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"By",input,str)) By=atof(str);
   else  {
      printf("In [External], By=? [Tesla]\n");
      fail=1;
   }
   if(FindParameters("External",1,"Bz",input,str)) Bz=atof(str);
   else  {
      printf("In [External], Bz=? [Tesla]\n");
      fail=1;
   }
*/
   //pml
   if(FindParameters("PML",1,"pml",input,str)) D->pml=whatONOFF(str);
   else D->pml=OFF;
   if(FindParameters("PML",1,"up_pml_cells",input,str)) D->pmlCellUp=atoi(str);
   else  D->pmlCellUp=1;
   if(FindParameters("PML",1,"pml_start",input,str))    D->pmlStart=atoi(str);
   else  { printf("In [PML], pml_start=? [what step]\n"); fail=1; }
   if(FindParameters("PML",1,"right_pml_cells",input,str))  D->pmlCellRight=atoi(str);
   else  D->pmlCellRight=1;
   if(FindParameters("PML",1,"left_pml_cells",input,str))   D->pmlCellLeft=atoi(str);
   else  D->pmlCellLeft=1;
   if(FindParameters("PML",1,"pml_r",input,str))    D->pmlr=atof(str);
   else  D->pmlr=0.04;
   if(FindParameters("PML",1,"pml_d",input,str))    D->pmld=atof(str);
   else  D->pmld=0.02;

   //field ionization
   if(FindParameters("Domain",1,"field_ionization",input,str)) D->fieldIonization=whatONOFF(str);
   else  D->fieldIonization=OFF;


   //additional Domain parameters  
   D->centerPick/=D->lambda;
   D->dz=D->divisionLambda;
   D->nx=((int)((maxX-minX)/D->lambda/D->dz));
   tmpDouble=D->dz/(D->gamma*(1+D->beta));
   D->dt=D->dz*D->dtRatio;
   if(D->fieldType==Split) {
     D->dtRatio=1.0;
     D->movingV=1.0;
     D->dt=D->dz;
     vg=1.0;
   } else if(D->fieldType==NoCherenkov || D->fieldType==Yee)
     vg=1.0+0.5*(1.0-D->movingV)*2.0*pi*D->dz;
   else ;
	D->shiftStart=0;
//   if(D->movingV==1.0) D->shiftDuration=D->maxStep+1;
//   else                D->shiftDuration=(int)(1.0/(vg-D->movingV));


   D->dr=D->dz*droverdz;
   D->ny=((int)((maxY-minY)/D->lambda/D->dr));

//   if(D->fieldType==Yee)
//     D->dt=D->dtRatio/sqrt(1.0/D->dr/D->dr+1.0/D->dz/D->dz);
//   else	;

   D->minXDomain=D->minYDomain=0;
   D->maxXDomain=D->nx;

   D->omega=2*M_PI*velocityC/D->lambda;
   if(D->boostOn==ON)   {
     D->minXSub=-D->nx;
     D->minYDomain=(int)(minY/D->lambda/D->dr);
   }
   else   {
     D->minXSub=0;
     D->minYDomain=(int)(minY/D->lambda/D->dr);
   }

   if(myrank==0)
   {
     printf("dz=%g,dr=%g,dt=%g,gamma=%g,divisionLambda=%g\n",D->dz,D->dr,D->dt,D->gamma,D->divisionLambda);
     printf("current_conservation=%d, 1:Lifschitz,  2:Davidson\n",D->currentCons);
   }
   else ;
   MPI_Barrier(MPI_COMM_WORLD);

   //ID track
   if(FindParameters("Domain",1,"tracking",input,str)) 
     D->tracking=whatONOFF(str);
   else  
     D->tracking=OFF;
   if(D->tracking==ON)
   {
     if(FindParameters("Domain",1,"track_save_step",input,str)) 
       D->trackSaveStep=atoi(str);
     else  
       D->trackSaveStep=1;
     if(myrank==0)   
     {
       in=fopen("trackFile","r");
       if(in!=NULL)  
       {
         cnt=0;
         while(fscanf(in,"%lf %lf %lf %lf %lf %d %d %d"
                 ,&x,&y,&px,&py,&pz,&id,&core,&species)!=EOF)
           cnt++;
         D->idNums=cnt;
         fclose(in);
       }  else   {
         printf("'trackFile' is missing.\n");
         fail=1;
       }
     } else	;
     MPI_Bcast(&(D->idNums),1,MPI_INT,0,MPI_COMM_WORLD);
   }
   else D->idNums=0;

   if(D->idNums>0)
   {
     D->trackID=(int *)malloc(D->idNums*sizeof(int));
     D->trackCore=(int *)malloc(D->idNums*sizeof(int));
     D->trackS=(int *)malloc(D->idNums*sizeof(int));
    
     if(myrank==0)
     {
       in=fopen("trackFile","r");
       for(i=0; i<D->idNums; i++)
         fscanf(in,"%lf %lf %lf %lf %lf %d %d %d"
                     ,&x,&y,&px,&py,&pz
                     ,&(D->trackID[i]),&(D->trackCore[i]),&(D->trackS[i]));
       fclose(in);
     }	
     else	;
     MPI_Barrier(MPI_COMM_WORLD);
     MPI_Bcast(D->trackID,D->idNums,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(D->trackCore,D->idNums,MPI_INT,0,MPI_COMM_WORLD);
     MPI_Bcast(D->trackS,D->idNums,MPI_INT,0,MPI_COMM_WORLD);
   }     
/*
   //Probe parameter
   if(FindParameters("Probe",1,"probeType",input,str)) probeType=atoi(str);
   else probeType=0;
   D->probeNum=0;

   if(probeType==0)
   {
     if(FindParameters("Probe",1,"probeNum",input,str)) D->probeNum=atoi(str);
     else  {
       printf("in [Probe], probeNum=?  [ea]\n");
       fail=1;
     }
     if(D->probeNum>0)
     {
       D->probeX=(int *)malloc(D->probeNum*sizeof(int));
       D->probeY=(int *)malloc(D->probeNum*sizeof(int));
       D->probeZ=(int *)malloc(D->probeNum*sizeof(int));
       for(i=0; i<D->probeNum; i++)
       {
         sprintf(name,"probeX%d",i);
         if(FindParameters("Probe",1,name,input,str))   
           D->probeX[i]=((int)(atof(str)/D->lambda/D->dx));      
         else  {
           printf("in [Probe], probeX%d=?\n",i);
           fail=1;
         }
         sprintf(name,"probeY%d",i);
         if(FindParameters("Probe",1,name,input,str))      
           D->probeY[i]=((int)(atof(str)/D->lambda/D->dy));      
         else  {
           printf("in [Probe], probeY%d=?\n",i);
           fail=1;
         }
         sprintf(name,"probeZ%d",i);
         if(FindParameters("Probe",1,name,input,str))      
           D->probeZ[i]=((int)(atof(str)/D->lambda/D->dz));      
         else  {
           printf("in [Probe], probeZ%d=?\n",i);
           fail=1;
         }
       }
     }  
   }
   else if(probeType==1)
   {
     if(FindParameters("Probe",1,"minProbeX",input,str)) minProbeX=atof(str);
     else  {
       printf("in [Probe], minProbeX=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeX",input,str)) maxProbeX=atof(str);
     else  {
       printf("in [Probe], maxProbeX=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeX",input,str)) numProbeX=atoi(str);
     else 
       numProbeX=1;
     if(FindParameters("Probe",1,"minProbeY",input,str)) minProbeY=atof(str);
     else  {
       printf("in [Probe], minProbeY=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeY",input,str)) maxProbeY=atof(str);
     else  {
       printf("in [Probe], maxProbeY=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeY",input,str)) numProbeY=atoi(str);
     else 
       numProbeY=1;
     if(FindParameters("Probe",1,"minProbeZ",input,str)) minProbeZ=atof(str);
     else  {
       printf("in [Probe], minProbeZ=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"maxProbeZ",input,str)) maxProbeZ=atof(str);
     else  {
       printf("in [Probe], maxProbeZ=? [m]\n");
       fail=1;
     }
     if(FindParameters("Probe",1,"numProbeZ",input,str)) numProbeZ=atoi(str);
     else 
       numProbeZ=1;
     if(numProbeX==0 || numProbeY==0 || numProbeZ==0)  {
       printf("in [Probe], it must be that numProbeX or numProbeY or numProbeZ > 0 !!\n");
       fail=1;
     }
           
     probeDx=(maxProbeX-minProbeX)/((double)numProbeX);
     probeDy=(maxProbeY-minProbeY)/((double)numProbeY);
     probeDz=(maxProbeZ-minProbeZ)/((double)numProbeZ);
     D->probeNum=numProbeX*numProbeY*numProbeZ;
     D->probeX=(int *)malloc(D->probeNum*sizeof(int));
     D->probeY=(int *)malloc(D->probeNum*sizeof(int));
     D->probeZ=(int *)malloc(D->probeNum*sizeof(int));

     n=0;
     for(i=0; i<numProbeX; i++)
       for(j=0; j<numProbeY; j++)
         for(k=0; k<numProbeZ; k++)
         {       
           tmpDouble=minProbeX+i*probeDx;
           D->probeX[n]=((int)(tmpDouble/D->lambda/D->dx));      
           tmpDouble=minProbeY+j*probeDy;
           D->probeY[n]=((int)(tmpDouble/D->lambda/D->dy));     
           tmpDouble=minProbeZ+k*probeDz;
           D->probeZ[n]=((int)(tmpDouble/D->lambda/D->dz));     
           n++;
         } 
   }			//End of else if(probeType=2)
*/    
   //additional Boost parameters
   factor=D->gamma*(1+D->beta);
   D->minT=(int)(D->maxStep/factor/factor); 	//boost frame iteration
   D->maxT=(int)((D->maxStep+D->beta*D->nx)/(1+D->beta)-factor*D->gamma*D->minT*D->beta);	//boost frame iteration
   if(myrank==0)
     printf("maxT=%d, nx=%d, ny=%d\n",D->maxT,D->nx,D->ny);
   else	;

   //additional external field parameters
   normalB=eMass*D->omega/(-eCharge);
   normalE=normalB*velocityC;
   Ext->E1=Ex/normalE;
   Ext->E2=Ey/normalE;
   Ext->E3=Ez/normalE;
   Ext->B1=Bx/normalB;
   Ext->B2=By/normalB;
   Ext->B3=Bz/normalB;

   //Laser parameter setting
   D->laserList = (LaserList *)malloc(sizeof(LaserList));
   D->laserList->next = NULL;
   L = D->laserList;
   rank = 1;
   while(findLaserParameters(rank,L,D,input)) 
   {
      LNew = (LaserList *)malloc(sizeof(LaserList));
      LNew->next = NULL;
      L->next=LNew;
      L=L->next;
      rank ++;
   }
   D->nLaser = rank-1;

/*
      D->laserl = cnt;
        D->laserX = (double *)malloc(cnt*sizeof(double ));
        D->laserI = (double *)malloc(cnt*sizeof(double ));
        D->laserPhase = (double *)malloc(cnt*sizeof(double ));
        in=fopen("laserIn","r");
        for(i=0; i<cnt; i++) {
          fscanf(in,"%lf %lf %lf %lf %lf %lf",&x,&tmp,&tmp,&tmp,&D->laserI[i],&D->laserPhase[i]);
          D->laserX[i]=x/D->lambda/D->dr;
        }
        fclose(in);
      } else ;
*/
   //Plasma parameter setting
   D->loadList = (LoadList *)malloc(sizeof(LoadList));
   D->loadList->next = NULL;
	D->index=1;
   LL = D->loadList;
   rank = 1;
   while(findLoadParameters(rank, LL, D,input)) 
   {
      New = (LoadList *)malloc(sizeof(LoadList));
      New->next = NULL;
      LL->next=New;
      LL=LL->next;
      rank ++;
   }
   D->nSpecies = rank-1;

   if(fail==1)
      exit(0);
   else	;

}

int findLaserParameters(int rank, LaserList *L,Domain *D,char *input)
{
   double positionX,positionY,positionZ;
   char name[100], str[100];
   int fail=0,polarity;

   if(FindParameters("Laser",rank,"loading_method",input,str)) L->loadMethod=whatLoadingMethod(str);
   else  { L->loadMethod=0; fail=1; }

   if(L->loadMethod)
   {
     if(FindParameters("Laser",rank,"polarity",input,str)) L->polarity=atoi(str);
     else  L->polarity=0;
     if(FindParameters("Laser",rank,"wavelength",input,str)) {
        L->lambda=atof(str);
        L->lambda*=D->gamma*(1.0+D->beta);
     } else  L->lambda=D->lambda;
     if(FindParameters("Laser",rank,"a0",input,str)) L->amplitude=atof(str);
     else  { printf("in [Laser], a0=??\n");  fail=1;  }
     if(FindParameters("Laser",rank,"focus",input,str)) L->focus=atof(str);
     else  { printf("in [Laser], focus=?  [m]\n"); fail=1;  }
     if(FindParameters("Laser",rank,"loadPositionX",input,str)) positionX=atof(str);
     else  positionX=0;
     if(FindParameters("Laser",rank,"add",input,str)) L->add=whatONOFF(str);
     else  L->add=OFF;

     if(FindParameters("Laser",rank,"mode",input,str)) L->mode=whatLaserMode(str);
     else  L->mode=Gaussian;
     if(L->mode==Gaussian) {
       if(FindParameters("Laser",rank,"rU",input,str)) L->rU=atof(str);
       else  { printf("in [Laser], rU=? [# of basic wavelength]\n"); fail=1; }
       if(FindParameters("Laser",rank,"rD",input,str)) L->rD=atof(str);
       else  { printf("in [Laser], rD=? [# of basic wavelength]\n"); fail=1; }
       if(FindParameters("Laser",rank,"retard",input,str)) L->retard=atof(str);
       else  { printf("in [Laser], retard=? [# of basic wavelength]\n"); fail=1; }
       if(FindParameters("Laser",rank,"flat",input,str)) L->flat=atof(str);
       else  L->flat=0.0;
       if(FindParameters("Laser",rank,"direction",input,str)) L->direction=atoi(str);
       else  L->direction=1;
       if(FindParameters("Laser",rank,"gdd",input,str)) L->gdd=atof(str);
       else  L->gdd=0.0;
       if(FindParameters("Laser",rank,"beamWaist",input,str)) L->beamWaist=atof(str);
       else  { printf("in [Laser], beamWaist=?  [m]\n"); fail=1;  }
     }  else ;

     //additional laser parameters
     L->omega=2*M_PI*velocityC/L->lambda;
     L->loadPointX=((int)(positionX/D->lambda/D->dz));   
     L->rayleighLength=M_PI/(L->lambda/D->gamma/(1.0+D->beta))*L->beamWaist*L->beamWaist/D->lambda;
     L->beamWaist=L->beamWaist/D->lambda;
     L->focus=L->focus/D->lambda;

     if(L->loadMethod==Shot) D->shiftStart=D->nx; else ;

     if(fail==1) exit(0); else ;
   }
   return L->loadMethod;
}


int findLoadParameters(int rank, LoadList *LL,Domain *D,char *input)
{
   LoadList *New;
   char name[100], str[100];
   int i,n,cnt,fail=0,tmpInt;
   double pointPosition,wp,pDt;
   double tmp,max,min;
	double gamma0,emitR,gammaR,sigR;
   double *shareDouble;
   int myrank, nTasks;
   MPI_Status status;

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);     
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);     

   if(FindParameters("Plasma",rank,"type",input,name)) {
     LL->type = whatPlasmaType(name);
     if(D->boostOn==ON) LL->type = BoostFrame; else ;
   } else LL->type=0;

   if(LL->type>0)
   {

/*    
      //testing optimal dx(divisionLambda) size
      wp=sqrt(LL->density*eCharge*eCharge/eMass/eps0);
      pDt=2*pi/wp;   
      pDt=pDt/(2*pi/D->omega)/20.0;
      if(D->dt>pDt)
      {
         printf("dt must be less then %g!\n",pDt);
         printf("So, divisionLambda>%g!\n",1/pDt);
         fail=1;
      }
*/

      if(FindParameters("Plasma",rank,"species",input,name)) LL->species = whatSpecies(name);
      else  LL->species = 0;
      if(FindParameters("Plasma",rank,"pair",input,name)) LL->pair=whatONOFF(name);
      else LL->pair=OFF;
      if(FindParameters("Plasma",rank,"numberRZ",input,str)) LL->numberRZ=atoi(str);
      else  { printf("in [Plasma], numberRZ=? \n"); fail=1; }
      if(FindParameters("Plasma",rank,"numberPhi",input,str)) LL->numberPhi=atoi(str);
      else  { printf("in [Plasma], numberPhi=? \n"); fail=1; }
      if(FindParameters("Plasma",rank,"temperature_r",input,str)) LL->temperatureR=atof(str);
      else   LL->temperatureR=0.0;	
      if(FindParameters("Plasma",rank,"temperature_z",input,str)) LL->temperatureZ=atof(str);
      else   LL->temperatureZ=0.0;	
      if(FindParameters("Plasma",rank,"given_min_pz",input,str)) LL->givenMinPx = atof(str);
      else  LL->givenMinPx = -1e9;
      LL->mass=whatMass(LL->species);
      LL->charge=whatCharge(LL->species);
      LL->criticalDensity=eps0*eMass*D->omega*D->omega/eCharge/eCharge;
//      LL->superP=LL->density*D->lambda*D->dx*D->lambda*D->dy*D->lambda*D->dz/LL->numberInCell;

      //setting ionization
      ionizationSetup(LL,LL->species);
      if(FindParameters("Plasma",rank,"given_min_a0",input,str)) 
         LL->givenMinA0 = atof(str);
      else  LL->givenMinA0 = 0.0;

      srand(1*(myrank+1));
      switch (LL->type)  {
      case Polygon :
        if(FindParameters("Plasma",rank,"density",input,str)) {
          LL->density=atof(str);
          LL->density*=D->gamma;
        } else  { printf("in [Plasma], density=? [m-3]\n"); fail=1; }
        if(FindParameters("Plasma",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
        else  {
          printf("in [Plasma], Xnodes=?\n");
          printf("Each nodes indicates the point of plasma density changing.\n");
          fail=1;
        }
        if(LL->xnodes>0)
        {
          LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
          LL->xn = (double *)malloc(LL->xnodes*sizeof(double));   
          for(i=0; i<LL->xnodes; i++)
          {
            sprintf(name,"X%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xpoint[i] = atof(str)/D->gamma/D->lambda/D->dz;
            else 
            { printf("X%d should be defined.\n",i);  fail=1; }

            sprintf(name,"Xn%d",i);
            if(FindParameters("Plasma",rank,name,input,str)) 
              LL->xn[i] = atof(str);
            else 
            { printf("Xn%d should be defined.\n",i);  fail=1; } 
          }
          LL->delZ=LL->xpoint[LL->xnodes-1]-LL->xpoint[0];
          LL->z0=(LL->xpoint[LL->xnodes-1]+LL->xpoint[0])*0.5;
        }
          if(FindParameters("Plasma",rank,"Ynodes",input,str)) LL->ynodes=atoi(str);
          else  {
            printf("in [Plasma], Ynodes=?\n");
            printf("Each nodes indicates the point of plasma density changing.\n");
            fail=1;
          }
          if(LL->ynodes>0)
          {
            LL->ypoint = (double *)malloc(LL->ynodes*sizeof(double));
            LL->yn = (double *)malloc(LL->ynodes*sizeof(double));   
            for(i=0; i<LL->ynodes; i++)
            {
              sprintf(name,"Y%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) {
                LL->ypoint[i] = atof(str)/D->lambda/D->dr;
              }
              else 
              { printf("Y%d should be defined.\n",i);  fail=1; }
 
              sprintf(name,"Yn%d",i);
              if(FindParameters("Plasma",rank,name,input,str)) 
                LL->yn[i] = atof(str);
              else 
              { printf("Yn%d should be defined.\n",i);  fail=1; } 
            }
          }
        if(FindParameters("Plasma",rank,"centerX",input,str))  
          LL->centerX=atof(str)/D->lambda/D->dz;
        else   LL->centerX=0.0;	
        LL->centerY=0.0;	
        if(FindParameters("Plasma",rank,"gauss_coef_X",input,str))  
          LL->gaussCoefX=atof(str)/D->lambda/D->dz;
        else   LL->gaussCoefX=1.0;
        if(FindParameters("Plasma",rank,"poly_coef_X",input,str))  
          LL->polyCoefX=atof(str)/D->lambda/D->dz;
        else   LL->polyCoefX=0.0;	
        if(FindParameters("Plasma",rank,"function_mode_X",input,str))  
          LL->modeX=whatFunctionMode(str);
        else   LL->modeX=0;	
        if(FindParameters("Plasma",rank,"function_mode_YZ",input,str))  
          LL->modeYZ=whatFunctionMode(str);
        else   LL->modeYZ=0;	
          if(FindParameters("Plasma",rank,"centerY",input,str))  
            LL->centerY=atof(str)/D->lambda/D->dr;
          else   LL->centerY=0.0;	
          if(FindParameters("Plasma",rank,"gauss_coef_YZ",input,str))  
            LL->gaussCoefYZ=atof(str)/D->lambda/D->dr;
          else   LL->gaussCoefYZ=1.0;	
          if(FindParameters("Plasma",rank,"poly_coef_YZ",input,str))  
            LL->polyCoefYZ=atof(str)*D->lambda*D->dr*D->lambda*D->dr;
          else   LL->polyCoefYZ=0.0;	
        break;
      case Beam :
        if(FindParameters("Plasma",rank,"loading_step",input,str)) LL->loadingStep=atoi(str);
        else LL->loadingStep=0;
        if(FindParameters("Plasma",rank,"energy",input,str)) LL->energy=atof(str);
        else { printf("In [Plasma], energy=? [MeV] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"energy_spread",input,str)) LL->spread=atof(str);
        else { printf("In [Plasma], energy_spread=? [<1.0] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"energy_chirp",input,str)) LL->eChirp=atof(str);
        else LL->eChirp=0.0;
        if(FindParameters("Plasma",rank,"norm_emittance_r",input,str)) LL->emitR=atof(str)*1e-6;
        else { printf("In [Plasma], norm_emittance_r=? [mm mrad] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"beta_r",input,str)) LL->betaR=atof(str);
        else { printf("In [Plasma], beta_r=? [m] \n"); fail=1;}
        if(FindParameters("Plasma",rank,"alpha_r",input,str)) LL->alphaR=atof(str);
        else { printf("In [Plasma], alpha_r=?  \n"); fail=1;}
        if(FindParameters("Plasma",rank,"focal_length",input,str)) LL->focalL=atof(str);
		  else LL->focalL=0.0;
        if(FindParameters("Plasma",rank,"peak_current",input,str)) LL->peakCurr=atof(str);
        else { printf("In [Plasma], peak_current=? [A] \n"); fail=1;}

        gamma0=LL->energy/mc2;
        emitR=LL->emitR/gamma0;
        gammaR=(1+LL->alphaR*LL->alphaR)/LL->betaR;
        sigR=sqrt(emitR/gammaR);
        LL->density=LL->peakCurr*0.5/M_PI/sigR/sigR/velocityC/eCharge;
        if(myrank==0) {
          printf("density=%g,sigR=%g,gammaR=%g,emitR=%g,alphaR=%g,loading_step=%d\n",LL->density,sigR,gammaR,emitR,LL->alphaR,LL->loadingStep);
        } else ;

        if(FindParameters("Plasma",rank,"gaussian_ONOFF",input,str)) LL->gaussMode=whatONOFF(str);
	else LL->gaussMode=OFF;
	if(LL->gaussMode==ON) {
          if(FindParameters("Plasma",rank,"position_z",input,str)) LL->posZ=atof(str);
          else LL->posZ=(D->maxX-D->minX)*0.5*D->lambda;
	  if(FindParameters("Plasma",rank,"gaussian_power",input,str)) LL->gaussPower=atof(str);
	  else LL->gaussPower=2;
	  if(FindParameters("Plasma",rank,"gaussian_sigma_z",input,str)) LL->sigZ=atof(str);
	  else LL->sigZ=1.0;
        } else {
          if(FindParameters("Plasma",rank,"Xnodes",input,str)) LL->xnodes=atoi(str);
          else  LL->xnodes=0;
          if(LL->xnodes>0) {
            LL->xpoint = (double *)malloc(LL->xnodes*sizeof(double));
            LL->xn = (double *)malloc(LL->xnodes*sizeof(double));
            for(i=0; i<LL->xnodes; i++)  {
              sprintf(name,"X%d",i);
              if(FindParameters("Plasma",rank,name,input,str))
                LL->xpoint[i] = atof(str)/D->gamma/D->lambda/D->dz;
              else { printf("X%d should be defined.\n",i);  fail=1; }

              sprintf(name,"Xn%d",i);
              if(FindParameters("Plasma",rank,name,input,str))LL->xn[i] = atof(str);
              else { printf("Xn%d should be defined.\n",i);  fail=1; }
            }
          } else ;
        }
        if(FindParameters("Plasma",rank,"num_energy_interval",input,str)) LL->numGam=atoi(str);                              else  LL->numGam=1;

		  break;




/*
      case Defined :

//        srand(time(NULL)*(myrank+1));

        if(FindParameters("Plasma",rank,"define_mode",input,str))  
          LL->defineMode=whatDefineMode(str);
        else   LL->defineMode=byNumber;	
        if(FindParameters("Plasma",rank,"number_defined",input,str))  
          LL->numDefined=atoi(str);
        else   LL->numDefined=0;	
        if(LL->defineMode==byDensity)
        {
          if(FindParameters("Plasma",rank,"minX",input,str))  
            LL->minX=atof(str);
          else   {   printf("minX=?  [m]\n");  exit(0);   }	
          if(FindParameters("Plasma",rank,"maxX",input,str))  
            LL->maxX=atof(str);
          else   {   printf("maxX=?  [m]\n");  exit(0);   }	
          if(D->dimension>1)
          {
            if(FindParameters("Plasma",rank,"minY",input,str))  
              LL->minY=atof(str);
            else   {   printf("minY=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxY",input,str))  
              LL->maxY=atof(str);
            else   {   printf("maxY=?  [m]\n");  exit(0);   }	
          }
          if(D->dimension>2)
          {
            if(FindParameters("Plasma",rank,"minZ",input,str))  
              LL->minZ=atof(str);
            else   {   printf("minZ=?  [m]\n");  exit(0);   }	
            if(FindParameters("Plasma",rank,"maxZ",input,str))  
              LL->maxZ=atof(str);
            else   {   printf("maxZ=?  [m]\n");  exit(0);   }	
          }
        }
        else 	;

        if(FindParameters("Plasma",rank,"xlength_particle",input,str))  
          LL->xLengDef=atof(str)/D->lambda;
        else   LL->xLengDef=0.0;
        if(LL->numDefined>0)	
        {
          LL->xPosition=(double *)malloc(LL->numDefined*sizeof(double));
          shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
          for(i=0; i<LL->numDefined; i++)
          {
            if(LL->defineMode==byNumber)
            {
              sprintf(name,"xPosition%d",i);
              if(FindParameters("Plasma",rank,name,input,str))  
                LL->xPosition[i]=atof(str)/D->lambda;
              else
              { printf("xPosition%d should be defined.\n",i); fail=1;}
            }
            else if(LL->defineMode==byDensity)
            {
              if(myrank==0)
                shareDouble[i]=(LL->minX+randomV()*(LL->maxX-LL->minX))/D->lambda;
              else 	;
            }
            else 	;
          }
          if(LL->defineMode==byDensity)
          {
            MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            for(i=0; i<LL->numDefined; i++)
              LL->xPosition[i]=shareDouble[i];
          }
          else	;
          free(shareDouble);
        }
        if(D->dimension>1)
        {
          if(FindParameters("Plasma",rank,"ylength_particle",input,str))  
            LL->yLengDef=atof(str)/D->lambda;
          else   LL->yLengDef=0.0;	
          if(LL->numDefined>0)
          {
            LL->yPosition=(double *)malloc(LL->numDefined*sizeof(double));
            shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
            for(i=0; i<LL->numDefined; i++)
            {
              if(LL->defineMode==byNumber)
              {
                sprintf(name,"yPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str))  
                  LL->yPosition[i]=atof(str)/D->lambda;
                else
                { printf("yPosition%d should be defined.\n",i); fail=1;}
              }
              else if(LL->defineMode==byDensity)
              {
                if(myrank==0)
                  shareDouble[i]=(LL->minY+randomV()*(LL->maxY-LL->minY))/D->lambda;
                else 	;
              }
              else 	;
            }
            if(LL->defineMode==byDensity)
            {
              MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
              for(i=0; i<LL->numDefined; i++)
                LL->yPosition[i]=shareDouble[i];
            }
            else	;
            free(shareDouble);
          }
        }		//End of demension>1
        if(D->dimension>2)
        {
          if(FindParameters("Plasma",rank,"zlength_particle",input,str))  
            LL->zLengDef=atof(str)/D->lambda;
          else   LL->zLengDef=0.0;	
          if(LL->numDefined>0)
          {
            LL->zPosition=(double *)malloc(LL->numDefined*sizeof(double));
            shareDouble=(double *)malloc(LL->numDefined*sizeof(double));
            for(i=0; i<LL->numDefined; i++)
            {
              if(LL->defineMode==byNumber)
              {
                sprintf(name,"zPosition%d",i);
                if(FindParameters("Plasma",rank,name,input,str))  
                  LL->zPosition[i]=atof(str)/D->lambda;
                else
                { printf("zPosition%d should be defined.\n",i); fail=1;}
              }
              else if(LL->defineMode==byDensity)
              {
                if(myrank==0)
                  shareDouble[i]=(LL->minZ+randomV()*(LL->maxZ-LL->minZ))/D->lambda;
                else 	;
              }
              else 	;
            }
            if(LL->defineMode==byDensity)
            {
              MPI_Bcast(shareDouble,LL->numDefined,MPI_DOUBLE,0,MPI_COMM_WORLD);
              MPI_Barrier(MPI_COMM_WORLD);
              for(i=0; i<LL->numDefined; i++)
                LL->zPosition[i]=shareDouble[i];
            }
            else	;
            free(shareDouble);  
          }
        }
        if(D->dimension==2)
        {
          LL->numDefPtcls=(int)(LL->xLengDef*LL->yLengDef/D->dx/D->dy*LL->numberInCell);
          cnt=LL->numDefPtcls;
          if(LL->numDefined>0)
          {
            LL->define=(double **)malloc(LL->numDefined*sizeof(double *));
            for(i=0; i<LL->numDefined; i++)
              LL->define[i]=(double *)malloc((cnt*2)*sizeof(double ));
 
            max=0;
            min=(double)D->maxStep;
            for(i=0; i<LL->numDefined; i++)
              for(n=0; n<cnt; n++)
              {
                tmp=(double)(randomV());
                LL->define[i][n]=LL->xPosition[i]-0.5*LL->xLengDef+tmp*LL->xLengDef;
                if(LL->define[i][n]>max) max=LL->define[i][n];
                if(LL->define[i][n]<min) min=LL->define[i][n];
                tmp=(double)(randomV());
                LL->define[i][n+cnt]=LL->yPosition[i]-0.5*LL->yLengDef+tmp*LL->yLengDef;
              }
          }
        }  
        else if(D->dimension==3)
        {
          LL->numDefPtcls=(int)(LL->xLengDef*LL->yLengDef*LL->zLengDef/D->dx/D->dy/D->dz*LL->numberInCell);
          cnt=LL->numDefPtcls;
          if(LL->numDefined>0)
          {
            LL->define=(double **)malloc(LL->numDefined*sizeof(double *));
            for(i=0; i<LL->numDefined; i++)
              LL->define[i]=(double *)malloc((LL->numDefPtcls*3)*sizeof(double ));
            max=0;
            min=(double)D->maxStep;
            for(i=0; i<LL->numDefined; i++)
              for(n=0; n<cnt; n++)
              {
                tmp=(double)(randomV());
                LL->define[i][n]=LL->xPosition[i]-0.5*LL->xLengDef+tmp*LL->xLengDef;
                if(LL->define[i][n]>max) max=LL->define[i][n];
                if(LL->define[i][n]<min) min=LL->define[i][n];
                tmp=(double)(randomV());
                LL->define[i][n+cnt]=LL->yPosition[i]-0.5*LL->yLengDef+tmp*LL->yLengDef;
                tmp=(double)(randomV());
                LL->define[i][n+2*cnt]=LL->zPosition[i]-0.5*LL->zLengDef+tmp*LL->zLengDef;
              }
          }
        } 	//End fo dimension=3 : defined Plasma 
//        tmpInt=(int)(max*D->divisionLambda);
//        tmpInt=tmpInt/D->shiftDuration;
//        LL->maxLoadTime=(int)((max+tmpInt)*D->divisionLambda);
//        LL->minLoadTime=(int)((min-tmpInt)*D->divisionLambda);
        if(FindParameters("Plasma",rank,"min_load_step",input,str))  
          LL->minLoadTime=atoi(str);
        else   LL->minLoadTime=0;
        if(FindParameters("Plasma",rank,"max_load_step",input,str))  
          LL->maxLoadTime=atoi(str);
        else   LL->maxLoadTime=0;
        break;
*/
      }
   
   }	//end of if(species)

   if(fail==1) exit(0); else ;

   return LL->type;
}



int whatDefineMode(char *str)
{
   if(strstr(str,"by_number")) 		return byNumber;
   else if(strstr(str,"by_density"))   	return byDensity;
   else 				return byNumber;
}

int whatONOFF(char *str)
{
   if(strstr(str,"ON")) 		return ON;
   else if(strstr(str,"OFF"))   	return OFF;
   else if(strstr(str,"File"))   return File;
   else 				return OFF;
}

int whatSaveMode(char *str)
{
   if(strstr(str,"TXT")) 		return TXT;
   else if(strstr(str,"HDF"))   	return HDF;
   else 				return TXT;
}

int whatFieldType(char *str)
{
   if(strstr(str,"Split")) 		return Split;
   else if(strstr(str,"Yee"))   	return Yee;
   else if(strstr(str,"NoCherenkov"))   return NoCherenkov;
   else if(strstr(str,"Split"))  	return Split;
   else 				return 0;
}


int whatPlasmaType(char *str)
{
   if(strstr(str,"Polygon"))    		return Polygon;
   else if(strstr(str,"Beam"))   	return Beam;
   else if(strstr(str,"Defined"))  	return Defined;
   else  {
     printf("No Plasma type!\n"); 
     exit(0);
   }
   return 0;
}


int whatFunctionMode(char *str)
{
   if(strstr(str,"Constant")) 		return Constant;
   else if(strstr(str,"Gaussian"))   	return Gaussian;
   else if(strstr(str,"Polynomial"))   	return Polynomial;
   else return 0;
}


int whatCurrentCons(char *str)
{
   if(strstr(str,"Lifschitz")) 		   return Lifschitz;
   else if(strstr(str,"Davidson"))   	return Davidson;
   else return Lifschitz;
}

double randomV()
{
   double r;
   int intRand, randRange=1000, rangeDev;

   intRand = rand() % randRange;
   r = ((double)intRand)/randRange;

   return r;
}

int whatLoadingMethod(char *str)
{
  if(strstr(str,"Boundary"))       return Boundary;
  else if(strstr(str,"Shot"))      return Shot;
  else                             return 0;
}

int whatLaserMode(char *str)
{
  if(strstr(str,"Gaussian"))       return Gaussian;
  else           return 0;
}
