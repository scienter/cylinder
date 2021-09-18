#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include "mpi.h"

void fieldIonization1D(Domain *D);
void fieldIonization2D(Domain *D);
void fieldIonization3D(Domain *D);
double randomValue(double beta);
void RungeKutta(double *W,double *prob,int iter,double dt,int start,int end,int flag);

void fieldIonization(Domain *D,int iteration)
{
  int n,i,j,m,s,Z,i1,j1,flag,ss;
  int initZ,cnt,ionLevel,index,sCnt,*sList;
  int istart,iend,jstart,jend,numMode,nSpecies,minRSub;
  double Ea,Edc,wa,Efield,c2,n_eff,l_eff,ionE,phase,dt;
  double **ionEnergy,**W,**prob;
  double z,r,phi,unitEBydt,value,weight,z1,r1,testProb;
  Particle **particle;
  particle=D->particle;
  LoadList *LL;

  ptclList *p,*New;
  
  istart=D->istart;   iend=D->iend;
  jstart=D->jstart;   jend=D->jend;
  minRSub=D->minYSub; numMode=D->numMode;
  dt=D->dt; nSpecies=D->nSpecies;

  Ea=eCharge*5.1e11/eMass/D->omega/velocityC;
  wa=4.13e16/D->omega;

  unitEBydt=13.59844/0.511*1.0e-6/dt;

  double rho[nSpecies];
  double minA0[nSpecies];
  int levels[nSpecies],species[nSpecies],ionFinal[nSpecies];
  ionEnergy=(double **)malloc(nSpecies*sizeof(double *));
  W=(double **)malloc(nSpecies*sizeof(double *));
  prob=(double **)malloc(nSpecies*sizeof(double *));

  LL=D->loadList;
  s=0; sCnt=0;
  while(LL->next)
  {
    rho[s]=LL->density/LL->criticalDensity;
    levels[s]=LL->levels;
    species[s]=LL->species;
    ionFinal[s]=LL->ionFinal;
    minA0[s]=LL->givenMinA0;
    //ionEnergy Z number : 0,1,2,...,Z-1
    ionEnergy[s]=(double *)malloc(LL->levels*sizeof(double ));
    W[s]=(double *)malloc(LL->levels*sizeof(double ));
    prob[s]=(double *)malloc((LL->levels+1)*sizeof(double ));
    prob[s][LL->levels]=1.0;

    for(i=0; i<LL->levels; i++) 
      ionEnergy[s][i]=LL->ionEnergy[i];

    if(LL->species!=Electron && LL->ionFinal==OFF) sCnt++; else ;
    LL=LL->next;
    s++;
  }

  sList=(int *)malloc(sCnt*sizeof(int ));
  LL=D->loadList;
  s=0; index=0;
  while(LL->next)
  {
    if(LL->species!=Electron && LL->ionFinal==OFF) {
      sList[index]=s; index++;
    }  else ;
    LL=LL->next;
    s++;
  }

/*
  //initialize J
  if(D->fieldType==Split)  {
    for(m=0; m<numMode; m++)  
      for(i=0; i<iend+3; i++) 
        for(j=0; j<jend+3; j++)  {  
          D->JzCR[m][i][j]=D->JzR[m][i][j];
          D->JzCI[m][i][j]=D->JzI[m][i][j];
          D->JpCR[m][i][j]=D->JpR[m][i][j];
          D->JpCI[m][i][j]=D->JpI[m][i][j];
          D->JrCR[m][i][j]=D->JrR[m][i][j];
          D->JrCI[m][i][j]=D->JrI[m][i][j];
          D->JzR[m][i][j]=0.0;
          D->JzI[m][i][j]=0.0;
          D->JrR[m][i][j]=0.0;
          D->JrI[m][i][j]=0.0;
          D->JpR[m][i][j]=0.0;
          D->JpI[m][i][j]=0.0;
        }
  } else {
    for(m=0; m<numMode; m++)  
      for(i=0; i<iend+3; i++) 
        for(j=0; j<jend+3; j++)  {  
          D->JzR[m][i][j]=0.0;
          D->JrR[m][i][j]=0.0;
          D->JpR[m][i][j]=0.0;
          D->JzI[m][i][j]=0.0;
          D->JrI[m][i][j]=0.0;
          D->JpI[m][i][j]=0.0;
        }
  }
*/

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
      {
        for(ss=0; ss<sCnt; ss++) {
          s=sList[ss];

          p=particle[i][j].head[s]->pt;

          while(p)  {
//            Edc=sqrt(p->Ez*p->Ez+p->Ex*p->Ex+p->Ey*p->Ey);
            Edc=sqrt(p->Ex*p->Ex+p->Ey*p->Ey);
            initZ=p->charge;
            weight=p->weight;
 
            if(initZ<levels[s] && Edc>minA0[s])  {
              testProb=randomValue(1.0);
  
              for(Z=initZ; Z<levels[s]; Z++) {
                Efield=Ea/Edc;
                ionE=ionEnergy[s][Z];
                n_eff=(Z+1)*sqrt(1.0/ionE);
                c2=pow(5.43656/n_eff,2.0*n_eff)/(n_eff);

                phase=2.0*Efield*pow(ionE,1.5);
                W[s][Z]=wa*c2*(0.5*ionE)*pow(phase,2.0*n_eff-1)*exp(-1.0*phase/3.0);
              }

              RungeKutta(&W[s][0],&prob[s][0],10,dt,initZ,levels[s],flag);
              flag=0;

              Z=initZ; ionLevel=-1;
              while(Z<levels[s]) {
                if(prob[s][Z]<testProb && testProb<prob[s][Z+1]) {
                  ionLevel=Z;
                  break;
                }
                Z++;
              }

              if(ionLevel>=initZ) { 
                value=ionEnergy[s][ionLevel]*unitEBydt*rho[s]/Edc/Edc*weight;
                //creating electrons
                p->charge=ionLevel+1;
                cnt=p->charge-initZ;
                for(n=0; n<cnt; n++) {
                  New = (ptclList *)malloc(sizeof(ptclList));
                  New->next = particle[i][j].head[s-1]->pt;
                  particle[i][j].head[s-1]->pt = New;
//                  z=randomValue(1.0); r=randomValue(1.0);
                  New->z = p->z; New->oldZ= p->z;
                  New->x = p->x; New->oldX= p->x;
                  New->y = p->y; New->oldY= p->y;
                  New->Ez = p->Ez; New->Ex=p->Ex; New->Ey=p->Ey;
                  New->Bz = p->Bz; New->Bx=p->Bx; New->By=p->By;
                  New->weight=p->weight;

                  New->pz=p->pz; New->px=p->px; New->py=p->py;
                  New->charge=-1;
                  New->index=p->index;
                  New->core=p->core;
/*
                  x=p->x;
                    x1=(x+0.5)-(int)(x+0.5);
                    i1=i+(int)(x+0.5)-1;
                    D->Jx[i1][j][k]+=(1.0-x1)*p->E1*value;
                    D->Jx[i1+1][j][k]+=x1*p->E1*value;
                    D->Jy[i][j][k]+=(1.0-x)*p->E2*value;
                    D->Jy[i+1][j][k]+=x*p->E2*value;
                    D->Jz[i][j][k]+=(1.0-x)*p->E3*value;
                    D->Jz[i+1][j][k]+=x*p->E3*value;
*/

                }
              } else ;	//End of electron generation 
            }  	//End of if(initZ<levels[s])

            p=p->next;

          }  //End of while(p)

      }      //End of for(s)
    }        //End of for(i,j,k)

  free(sList);
  for(n=0; n<nSpecies; n++) {
    free(ionEnergy[n]); free(W[n]); free(prob[n]);
  } free(ionEnergy); free(W); free(prob);
}



void ionizationSetup(LoadList *LL,int species)
{

  switch (species)  {
    case Electron :
      LL->levels=0;
    case HPlus0 :
    case HPlus1 :
      LL->levels=1;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.0;
      if(species==HPlus1) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case HePlus0 :
    case HePlus1 :
    case HePlus2 :
      LL->levels=2;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.808;
      LL->ionEnergy[1]=4.002;
      if(species==HePlus2) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case CPlus0 :
    case CPlus1 :
    case CPlus2 :
    case CPlus3 :
    case CPlus4 :
    case CPlus5 :
    case CPlus6 :
      LL->levels=6;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=0.828;
      LL->ionEnergy[1]=1.793;
      LL->ionEnergy[2]=3.522;
      LL->ionEnergy[3]=4.743;
      LL->ionEnergy[4]=28.833;
      LL->ionEnergy[5]=36.033;
      if(species==CPlus6) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case NPlus0 :
    case NPlus1 :
    case NPlus2 :
    case NPlus3 :
    case NPlus4 :
    case NPlus5 :
    case NPlus6 :
    case NPlus7 :
      LL->levels=7;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.069;
      LL->ionEnergy[1]=2.177;
      LL->ionEnergy[2]=3.489;
      LL->ionEnergy[3]=5.697;
      LL->ionEnergy[4]=7.199;
      LL->ionEnergy[5]=40.598;
      LL->ionEnergy[6]=49.053;
      if(species==NPlus7) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
    case NePlus0 :
    case NePlus1 :
    case NePlus2 :
    case NePlus3 :
    case NePlus4 :
    case NePlus5 :
    case NePlus6 :
    case NePlus7 :
    case NePlus8 :
    case NePlus9 :
    case NePlus10 :
      LL->levels=10;
      LL->ionEnergy=(double *)malloc(LL->levels*sizeof(double) );
      LL->ionEnergy[0]=1.586;
      LL->ionEnergy[1]=3.012;
      LL->ionEnergy[2]=4.666;
      LL->ionEnergy[3]=7.142;
      LL->ionEnergy[4]=9.281;
      LL->ionEnergy[5]=11.614;
      LL->ionEnergy[6]=15.243;
      LL->ionEnergy[7]=17.583;
      LL->ionEnergy[8]=87.939;
      LL->ionEnergy[9]=100.173;
      if(species==NPlus7) LL->ionFinal=ON; else LL->ionFinal=OFF;
      break;
  }

}

double calW(double ionEnergy,double charge,double omega)
{
  double c2,n;
  double Ea=5.1e11, wa=4.13e16, euler=2.71828;
  
  Ea=eCharge*Ea/eMass/omega/velocityC;
  n=charge*sqrt(1.0/ionEnergy);

  c2=0.5/pi/n*pow(2*euler/n,2*n);
}
