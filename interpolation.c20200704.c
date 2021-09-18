#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
#include "math.h"

void interpolation_NDFX_1st(Domain *D,External *Ext,int iteration);
void interpolation_Yee_1st(Domain *D,External *Ext);

void interpolation(Domain *D,External *Ext,int iteration)
{
  switch((D->fieldType-1)*2+(D->interpolationType-1)) {  
  case ((Yee-1)*2+(FIRST-1)) :
  case ((NoCherenkov-1)*2+(FIRST-1)) :
    interpolation_Yee_1st(D,Ext);	//2D
    break;
  case ((NDFX-1)*2+(FIRST-1)) :
    interpolation_NDFX_1st(D,Ext,iteration);	//2D
    break;
  default :
    printf("In interpolation, what interpolationType(%d)?\n",D->interpolationType);
  }
}

void interpolation_NDFX_1st(Domain *D,External *Ext,int iteration)
{
   int ii,jj,i,j,m,i1,j1,istart,iend,jstart,jend,s,rank,numMode;
   double Bz,Br,Bp,Bx,By,Ez,Er,Ep,Ex,Ey,Pr,Pl,Sr,Sl;
   double coss[D->numMode],sins[D->numMode];
   double coss2[D->numMode],sins2[D->numMode];
   double extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,R,r,invR;
   double wz1,wz2,wr1,wr2,WZ1,WZ2,WR1,WR2,y1,y2;
   double Pr0,Pr1,Pr2,Pr3,Pl0,Pl1,Pl2,Pl3;
   double Sr0,Sr1,Sr2,Sr3,Sl0,Sl1,Sl2,Sl3;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   numMode=D->numMode;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;
   extE3=Ext->E3;   extB1=Ext->B1;
   extB2=Ext->B2;   extB3=Ext->B3;

   for(i=istart; i<iend; i++)
     for(j=jstart+1; j<jend; j++) 
     { 
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j].head[s]->pt;
         while(p) {
           z=p->z;  x=p->x; y=p->y;
           R=sqrt(x*x+y*y); invR=1.0/R;
           r=R-(j-jstart);

           wz2=z;             wz1=1.0-wz2;
           wr2=r;             wr1=1.0-wr2;
           i1=((int)(i+z+0.5));
           j1=((int)(j+r+0.5));

           WZ2=z+0.5-((int)(z+0.5));  WZ1=1.0-WZ2;
           WR2=r+0.5-((int)(r+0.5));  WR1=1.0-WR2;

           coss[0]=1.0; sins[0]=0.0;
           coss[1]=x*invR; sins[1]=y*invR;
           for(m=2; m<numMode; m++) {
             coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
             sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
           }

           Ez=wz1*wr1*D->EzR[0][i][j]
             +wz2*wr1*D->EzR[0][i+1][j]
             +wz1*wr2*D->EzR[0][i][j+1]
             +wz2*wr2*D->EzR[0][i+1][j+1];
           Bz=wz1*wr1*D->BzR[0][i][j]
             +wz2*wr1*D->BzR[0][i+1][j]
             +wz1*wr2*D->BzR[0][i][j+1]
             +wz2*wr2*D->BzR[0][i+1][j+1];
//if(iteration>=89) printf("before Pr[0], i1=%d, j1=%d,z=%g,r=%g, istart=%d,jend=%d\n",i1,j1,z,r,iend,jend);
           Pr=WZ1*WR1*D->PrR[0][i1-1][j1-1]
             +WZ2*WR1*D->PrR[0][i1][j1-1]
             +WZ1*WR2*D->PrR[0][i1-1][j1]
             +WZ2*WR2*D->PrR[0][i1][j1];
           Pl=WZ1*WR1*D->PlR[0][i1-1][j1-1]
             +WZ2*WR1*D->PlR[0][i1][j1-1]
             +WZ1*WR2*D->PlR[0][i1-1][j1]
             +WZ2*WR2*D->PlR[0][i1][j1];
           Sr=WZ1*WR1*D->SrR[0][i1-1][j1-1]
             +WZ2*WR1*D->SrR[0][i1][j1-1]
             +WZ1*WR2*D->SrR[0][i1-1][j1]
             +WZ2*WR2*D->SrR[0][i1][j1];
           Sl=WZ1*WR1*D->SlR[0][i1-1][j1-1]
             +WZ2*WR1*D->SlR[0][i1][j1-1]
             +WZ1*WR2*D->SlR[0][i1-1][j1]
             +WZ2*WR2*D->SlR[0][i1][j1];

           for(m=1; m<numMode; m++) {
             Ez+=((wz1*wr1*D->EzR[m][i][j]
                  +wz2*wr1*D->EzR[m][i+1][j]
                  +wz1*wr2*D->EzR[m][i][j+1]
                  +wz2*wr2*D->EzR[m][i+1][j+1])*coss[m]
                 -(wz1*wr1*D->EzI[m][i][j]
                  +wz2*wr1*D->EzI[m][i+1][j]
                  +wz1*wr2*D->EzI[m][i][j+1]
                  +wz2*wr2*D->EzI[m][i+1][j+1])*sins[m]);
             Bz+=((wz1*wr1*D->BzR[m][i][j]
                  +wz2*wr1*D->BzR[m][i+1][j]
                  +wz1*wr2*D->BzR[m][i][j+1]
                  +wz2*wr2*D->BzR[m][i+1][j+1])*coss[m]
                 -(wz1*wr1*D->BzI[m][i][j]
                  +wz2*wr1*D->BzI[m][i+1][j]
                  +wz1*wr2*D->BzI[m][i][j+1]
                  +wz2*wr2*D->BzI[m][i+1][j+1])*sins[m]);
             Pr+=((WZ1*WR1*D->PrR[m][i1-1][j1-1]
                  +WZ2*WR1*D->PrR[m][i1][j1-1]
                  +WZ1*WR2*D->PrR[m][i1-1][j1]
                  +WZ2*WR2*D->PrR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->PrI[m][i1-1][j1-1]
                  +WZ2*WR1*D->PrI[m][i1][j1-1]
                  +WZ1*WR2*D->PrI[m][i1-1][j1]
                  +WZ2*WR2*D->PrI[m][i1][j1])*sins[m]);
             Pl+=((WZ1*WR1*D->PlR[m][i1-1][j1-1]
                  +WZ2*WR1*D->PlR[m][i1][j1-1]
                  +WZ1*WR2*D->PlR[m][i1-1][j1]
                  +WZ2*WR2*D->PlR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->PlI[m][i1-1][j1-1]
                  +WZ2*WR1*D->PlI[m][i1][j1-1]
                  +WZ1*WR2*D->PlI[m][i1-1][j1]
                  +WZ2*WR2*D->PlI[m][i1][j1])*sins[m]);
             Sr+=((WZ1*WR1*D->SrR[m][i1-1][j1-1]
                  +WZ2*WR1*D->SrR[m][i1][j1-1]
                  +WZ1*WR2*D->SrR[m][i1-1][j1]
                  +WZ2*WR2*D->SrR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->SrI[m][i1-1][j1-1]
                  +WZ2*WR1*D->SrI[m][i1][j1-1]
                  +WZ1*WR2*D->SrI[m][i1-1][j1]
                  +WZ2*WR2*D->SrI[m][i1][j1])*sins[m]);
             Sl+=((WZ1*WR1*D->SlR[m][i1-1][j1-1]
                  +WZ2*WR1*D->SlR[m][i1][j1-1]
                  +WZ1*WR2*D->SlR[m][i1-1][j1]
                  +WZ2*WR2*D->SlR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->SlI[m][i1-1][j1-1]
                  +WZ2*WR1*D->SlI[m][i1][j1-1]
                  +WZ1*WR2*D->SlI[m][i1-1][j1]
                  +WZ2*WR2*D->SlI[m][i1][j1])*sins[m]);
           }

           Er=(Pr+Pl)*0.5; Bp=(Pr-Pl)*0.5;
           Ep=(Sl+Sr)*0.5; Br=(Sl-Sr)*0.5;
           Ex=Er*coss[1]-Ep*sins[1];
           Ey=Er*sins[1]+Ep*coss[1];
           Bx=Br*coss[1]-Bp*sins[1];
           By=Br*sins[1]+Bp*coss[1];

           p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
           p->Bz=Bz; p->Bx=Bx; p->By=By;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)

   j=jstart;
   for(i=istart; i<iend; i++)
     { 
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j].head[s]->pt;
         while(p) {
           z=p->z;  x=p->x; y=p->y;
           r=sqrt(x*x+y*y); invR=1.0/r;

           coss[0]=1.0; sins[0]=0.0;
           coss[1]=x*invR; sins[1]=y*invR;
           for(m=2; m<numMode; m++) {
             coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
             sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
           }

           wz2=z;             wz1=1.0-wz2;
           WZ2=z+0.5-((int)(z+0.5));  WZ1=1.0-WZ2;
           i1=((int)(i+z+0.5));
           j1=((int)(j+r+0.5));

           wr2=r;                 wr1=1.0-wr2;
           
           //Pr,Pl,Sr,Sl
           if(j1==jstart) {
             coss2[0]=1.0;    sins2[0]=0.0;
             coss2[1]=-x*invR; sins2[1]=-y*invR;
             for(m=2; m<numMode; m++) {
               coss2[m]=coss2[m-1]*coss2[1]-sins2[m-1]*sins2[1];
               sins2[m]=sins2[m-1]*coss2[1]+coss2[m-1]*sins2[1];
             }

             Pr=Pl=Sr=Sl=0.0;
             WR2=r/0.5;     WR1=1.0-WR2;
				 m=0;
             Pr+=WZ1*WR2*D->PrR[m][i1-1][j1]
				    +WZ2*WR2*D->PrR[m][i1][j1];
             Pl+=WZ1*WR2*D->PlR[m][i1-1][j1]
				    +WZ2*WR2*D->PlR[m][i1][j1];
             Sr+=WZ1*WR2*D->SrR[m][i1-1][j1]
				    +WZ2*WR2*D->SrR[m][i1][j1];
             Sl+=WZ1*WR2*D->SlR[m][i1-1][j1]
				    +WZ2*WR2*D->SlR[m][i1][j1];
             m=1;
             Pr+=WZ1*D->PrR[m][i1-1][j1]*coss[m]-WZ1*D->PrI[m][i1-1][j1]*sins[m]
                +WZ2*D->PrR[m][i1][j1]*coss[m]  -WZ2*D->PrI[m][i1][j1]*sins[m];
             Pl+=WZ1*D->PlR[m][i1-1][j1]*coss[m]-WZ1*D->PlI[m][i1-1][j1]*sins[m]
                +WZ2*D->PlR[m][i1][j1]*coss[m]  -WZ2*D->PlI[m][i1][j1]*sins[m];
             Sr+=WZ1*D->SrR[m][i1-1][j1]*coss[m]-WZ1*D->SrI[m][i1-1][j1]*sins[m]
                +WZ2*D->SrR[m][i1][j1]*coss[m]  -WZ2*D->SrI[m][i1][j1]*sins[m];
             Sl+=WZ1*D->SlR[m][i1-1][j1]*coss[m]-WZ1*D->SlI[m][i1-1][j1]*sins[m]
                +WZ2*D->SlR[m][i1][j1]*coss[m]  -WZ2*D->SlI[m][i1][j1]*sins[m];
             for(m=2; m<numMode; m++) {
               Pr+=WZ1*WR2*D->PrR[m][i1-1][j1]*coss[m]-WZ1*WR2*D->PrI[m][i1-1][j1]*sins[m]
                  +WZ2*WR2*D->PrR[m][i1][j1]*coss[m]  -WZ2*WR2*D->PrI[m][i1][j1]*sins[m];
               Pl+=WZ1*WR2*D->PlR[m][i1-1][j1]*coss[m]-WZ1*WR2*D->PlI[m][i1-1][j1]*sins[m]
                  +WZ2*WR2*D->PlR[m][i1][j1]*coss[m]  -WZ2*WR2*D->PlI[m][i1][j1]*sins[m];
               Sr+=WZ1*WR2*D->SrR[m][i1-1][j1]*coss[m]-WZ1*WR2*D->SrI[m][i1-1][j1]*sins[m]
                  +WZ2*WR2*D->SrR[m][i1][j1]*coss[m]  -WZ2*WR2*D->SrI[m][i1][j1]*sins[m];
               Sl+=WZ1*WR2*D->SlR[m][i1-1][j1]*coss[m]-WZ1*WR2*D->SlI[m][i1-1][j1]*sins[m]
                  +WZ2*WR2*D->SlR[m][i1][j1]*coss[m]  -WZ2*WR2*D->SlI[m][i1][j1]*sins[m];
				 }
/*
             Pr0=D->PrR[0][i1-1][j1]; Pr1=D->PrR[0][i1][j1];
             Pr2=D->PrR[0][i1-1][j1]; Pr3=D->PrR[0][i1][j1];
             Pl0=D->PlR[0][i1-1][j1]; Pl1=D->PlR[0][i1][j1];
             Pl2=D->PlR[0][i1-1][j1]; Pl3=D->PlR[0][i1][j1];
             Sr0=D->SrR[0][i1-1][j1]; Sr1=D->SrR[0][i1][j1];
             Sr2=D->SrR[0][i1-1][j1]; Sr3=D->SrR[0][i1][j1];
             Sl0=D->SlR[0][i1-1][j1]; Sl1=D->SlR[0][i1][j1];
             Sl2=D->SlR[0][i1-1][j1]; Sl3=D->SlR[0][i1][j1];
             for(m=1; m<numMode; m++) {
               Pr0+=D->PrR[m][i1-1][j1]*coss[m] -D->PrI[m][i1-1][j1]*sins[m];
               Pr1+=D->PrR[m][i1][j1]*coss[m]   -D->PrI[m][i1][j1]*sins[m];
               Pr2+=D->PrR[m][i1-1][j1]*coss2[m]-D->PrI[m][i1-1][j1]*sins2[m];
               Pr3+=D->PrR[m][i1][j1]*coss2[m]  -D->PrI[m][i1][j1]*sins2[m];

               Pl0+=D->PlR[m][i1-1][j1]*coss[m] -D->PlI[m][i1-1][j1]*sins[m];
               Pl1+=D->PlR[m][i1][j1]*coss[m]   -D->PlI[m][i1][j1]*sins[m];
               Pl2+=D->PlR[m][i1-1][j1]*coss2[m]-D->PlI[m][i1-1][j1]*sins2[m];
               Pl3+=D->PlR[m][i1][j1]*coss2[m]  -D->PlI[m][i1][j1]*sins2[m];

               Sr0+=D->SrR[m][i1-1][j1]*coss[m] -D->SrI[m][i1-1][j1]*sins[m];
               Sr1+=D->SrR[m][i1][j1]*coss[m]   -D->SrI[m][i1][j1]*sins[m];
               Sr2+=D->SrR[m][i1-1][j1]*coss2[m]-D->SrI[m][i1-1][j1]*sins2[m];
               Sr3+=D->SrR[m][i1][j1]*coss2[m]  -D->SrI[m][i1][j1]*sins2[m];

               Sl0+=D->SlR[m][i1-1][j1]*coss[m] -D->SlI[m][i1-1][j1]*sins[m];
               Sl1+=D->SlR[m][i1][j1]*coss[m]   -D->SlI[m][i1][j1]*sins[m];
               Sl2+=D->SlR[m][i1-1][j1]*coss2[m]-D->SlI[m][i1-1][j1]*sins2[m];
               Sl3+=D->SlR[m][i1][j1]*coss2[m]  -D->SlI[m][i1][j1]*sins2[m];
             }

             WR2=0.5+r;     WR1=1.0-WR2;
             Pr=WZ1*WR2*Pr0+WZ2*WR2*Pr1+WZ1*WR1*Pr2+WZ2*WR1*Pr3;
             Pl=WZ1*WR2*Pl0+WZ2*WR2*Pl1+WZ1*WR1*Pl2+WZ2*WR1*Pl3;
             Sr=WZ1*WR2*Sr0+WZ2*WR2*Sr1+WZ1*WR1*Sr2+WZ2*WR1*Sr3;
             Sl=WZ1*WR2*Sl0+WZ2*WR2*Sl1+WZ1*WR1*Sl2+WZ2*WR1*Sl3;
*/				 
           } else {
             WR2=r+0.5-((int)(r+0.5));  WR1=1.0-WR2;

             Pr=WZ1*WR1*D->PrR[0][i1-1][j1-1]
               +WZ2*WR1*D->PrR[0][i1][j1-1]
               +WZ1*WR2*D->PrR[0][i1-1][j1]
               +WZ2*WR2*D->PrR[0][i1][j1];
             Pl=WZ1*WR1*D->PlR[0][i1-1][j1-1]
               +WZ2*WR1*D->PlR[0][i1][j1-1]
               +WZ1*WR2*D->PlR[0][i1-1][j1]
               +WZ2*WR2*D->PlR[0][i1][j1];
             Sr=WZ1*WR1*D->SrR[0][i1-1][j1-1]
               +WZ2*WR1*D->SrR[0][i1][j1-1]
               +WZ1*WR2*D->SrR[0][i1-1][j1]
               +WZ2*WR2*D->SrR[0][i1][j1];
             Sl=WZ1*WR1*D->SlR[0][i1-1][j1-1]
               +WZ2*WR1*D->SlR[0][i1][j1-1]
               +WZ1*WR2*D->SlR[0][i1-1][j1]
               +WZ2*WR2*D->SlR[0][i1][j1];
             for(m=1; m<numMode; m++) {
               Pr+=((WZ1*WR1*D->PrR[m][i1-1][j1-1]
                  +WZ2*WR1*D->PrR[m][i1][j1-1]
                  +WZ1*WR2*D->PrR[m][i1-1][j1]
                  +WZ2*WR2*D->PrR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->PrI[m][i1-1][j1-1]
                  +WZ2*WR1*D->PrI[m][i1][j1-1]
                  +WZ1*WR2*D->PrI[m][i1-1][j1]
                  +WZ2*WR2*D->PrI[m][i1][j1])*sins[m]);
               Pl+=((WZ1*WR1*D->PlR[m][i1-1][j1-1]
                  +WZ2*WR1*D->PlR[m][i1][j1-1]
                  +WZ1*WR2*D->PlR[m][i1-1][j1]
                  +WZ2*WR2*D->PlR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->PlI[m][i1-1][j1-1]
                  +WZ2*WR1*D->PlI[m][i1][j1-1]
                  +WZ1*WR2*D->PlI[m][i1-1][j1]
                  +WZ2*WR2*D->PlI[m][i1][j1])*sins[m]);
               Sr+=((WZ1*WR1*D->SrR[m][i1-1][j1-1]
                  +WZ2*WR1*D->SrR[m][i1][j1-1]
                  +WZ1*WR2*D->SrR[m][i1-1][j1]
                  +WZ2*WR2*D->SrR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->SrI[m][i1-1][j1-1]
                  +WZ2*WR1*D->SrI[m][i1][j1-1]
                  +WZ1*WR2*D->SrI[m][i1-1][j1]
                  +WZ2*WR2*D->SrI[m][i1][j1])*sins[m]);
               Sl+=((WZ1*WR1*D->SlR[m][i1-1][j1-1]
                  +WZ2*WR1*D->SlR[m][i1][j1-1]
                  +WZ1*WR2*D->SlR[m][i1-1][j1]
                  +WZ2*WR2*D->SlR[m][i1][j1])*coss[m]
                 -(WZ1*WR1*D->SlI[m][i1-1][j1-1]
                  +WZ2*WR1*D->SlI[m][i1][j1-1]
                  +WZ1*WR2*D->SlI[m][i1-1][j1]
                  +WZ2*WR2*D->SlI[m][i1][j1])*sins[m]);
             }
           }


           //Ex, Bx
           Ez=wz1*wr1*D->EzR[0][i][j]
             +wz2*wr1*D->EzR[0][i+1][j]
             +wz1*wr2*D->EzR[0][i][j+1]
             +wz2*wr2*D->EzR[0][i+1][j+1];
           Bz=wz1*wr1*D->BzR[0][i][j]
             +wz2*wr1*D->BzR[0][i+1][j]
             +wz1*wr2*D->BzR[0][i][j+1]
             +wz2*wr2*D->BzR[0][i+1][j+1];
           for(m=1; m<numMode; m++) {
             Ez+=((wz1*wr1*D->EzR[m][i][j]
                  +wz2*wr1*D->EzR[m][i+1][j]
                  +wz1*wr2*D->EzR[m][i][j+1]
                  +wz2*wr2*D->EzR[m][i+1][j+1])*coss[m]
                 -(wz1*wr1*D->EzI[m][i][j]
                  +wz2*wr1*D->EzI[m][i+1][j]
                  +wz1*wr2*D->EzI[m][i][j+1]
                  +wz2*wr2*D->EzI[m][i+1][j+1])*sins[m]);
             Bz+=((wz1*wr1*D->BzR[m][i][j]
                  +wz2*wr1*D->BzR[m][i+1][j]
                  +wz1*wr2*D->BzR[m][i][j+1]
                  +wz2*wr2*D->BzR[m][i+1][j+1])*coss[m]
                 -(wz1*wr1*D->BzI[m][i][j]
                  +wz2*wr1*D->BzI[m][i+1][j]
                  +wz1*wr2*D->BzI[m][i][j+1]
                  +wz2*wr2*D->BzI[m][i+1][j+1])*sins[m]);
           }

           Er=(Pr+Pl)*0.5; Bp=(Pr-Pl)*0.5;
           Ep=(Sl+Sr)*0.5; Br=(Sl-Sr)*0.5;
           Ex=Er*coss[1]-Ep*sins[1];
           Ey=Er*sins[1]+Ep*coss[1];
           Bx=Br*coss[1]-Bp*sins[1];
           By=Br*sins[1]+Bp*coss[1];

           p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
           p->Bz=Bz; p->Bx=Bx; p->By=By;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

void interpolation_Yee_1st(Domain *D,External *Ext)
{
   int i,j,m,i1,j1,ii,jj,istart,iend,jstart,jend,s,rank,numMode;
   double Bz,Br,Bp,Bx,By,Ez,Er,Ep,Ex,Ey,yR1,yR2,yI1,yI2;
   double Bz0,Bz1,Bz2,Bz3;
   double Bp0,Bp1,Bp2,Bp3;
   double Er0,Er1,Er2,Er3,diag;
   double coss[D->numMode],sins[D->numMode];
   double cos0[D->numMode],sin0[D->numMode];
   double cos1[D->numMode],sin1[D->numMode];
   double cos2[D->numMode],sin2[D->numMode];
   double cos3[D->numMode],sin3[D->numMode];
   double extE1,extE2,extE3,extB1,extB2,extB3,x,y,z,R,r,invR;
   double wz[2],wr[2],WZ[2],WR[2];
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;
   numMode=D->numMode;

   Particle **particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;
   extE3=Ext->E3;   extB1=Ext->B1;
   extB2=Ext->B2;   extB3=Ext->B3;


   for(i=istart; i<iend; i++)
     for(j=jstart+1; j<jend; j++) 
     { 
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j].head[s]->pt;
         while(p) {
           z=p->z;  x=p->x; y=p->y;
           R=sqrt(x*x+y*y); invR=1.0/R;
           r=R-(j-jstart);

           wz[1]=z;           wz[0]=1.0-wz[1];
           wr[1]=r;           wr[0]=1.0-wr[1];
           i1=((int)(i+z+0.5));
           j1=((int)(j+r+0.5));
           WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];
           WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];

           coss[0]=1.0; sins[0]=0.0;
           coss[1]=x*invR; sins[1]=y*invR;
           for(m=2; m<numMode; m++) {
             coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
             sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
           }

           Bz=Br=Bp=Ez=Er=Ep=0.0;
           for(ii=0; ii<2; ii++)
             for(jj=0; jj<2; jj++) {
               Bz+=wz[ii]*WR[jj]*D->BzNowR[0][i+ii][j1-1+jj];
               Br+=WZ[ii]*wr[jj]*D->BrNowR[0][i1-1+ii][j+jj];
               Bp+=WZ[ii]*WR[jj]*D->BpNowR[0][i1-1+ii][j1-1+jj];
               Ez+=WZ[ii]*wr[jj]*D->EzR[0][i1-1+ii][j+jj];
               Er+=wz[ii]*WR[jj]*D->ErR[0][i+ii][j1-1+jj];
               Ep+=wz[ii]*wr[jj]*D->EpR[0][i+ii][j+jj];
             }

           for(m=1; m<numMode; m++) {
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 Bz+=wz[ii]*WR[jj]*D->BzNowR[m][i+ii][j1-1+jj]*coss[m];
                 Bz-=wz[ii]*WR[jj]*D->BzNowI[m][i+ii][j1-1+jj]*sins[m];
                 Br+=WZ[ii]*wr[jj]*D->BrNowR[m][i1-1+ii][j+jj]*coss[m];
                 Br-=WZ[ii]*wr[jj]*D->BrNowI[m][i1-1+ii][j+jj]*sins[m];
                 Bp+=WZ[ii]*WR[jj]*D->BpNowR[m][i1-1+ii][j1-1+jj]*coss[m];
                 Bp-=WZ[ii]*WR[jj]*D->BpNowI[m][i1-1+ii][j1-1+jj]*sins[m];
                 Ez+=WZ[ii]*wr[jj]*D->EzR[m][i1-1+ii][j+jj]*coss[m];
                 Ez-=WZ[ii]*wr[jj]*D->EzI[m][i1-1+ii][j+jj]*sins[m];
                 Er+=wz[ii]*WR[jj]*D->ErR[m][i+ii][j1-1+jj]*coss[m];
                 Er-=wz[ii]*WR[jj]*D->ErI[m][i+ii][j1-1+jj]*sins[m];
                 Ep+=wz[ii]*wr[jj]*D->EpR[m][i+ii][j+jj]*coss[m];
                 Ep-=wz[ii]*wr[jj]*D->EpI[m][i+ii][j+jj]*sins[m];
               }
           }
           Ex=Er*coss[1]-Ep*sins[1];
           Ey=Er*sins[1]+Ep*coss[1];
           Bx=Br*coss[1]-Bp*sins[1];
           By=Br*sins[1]+Bp*coss[1];

           p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
           p->Bz=Bz; p->Bx=Bx; p->By=By;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)

   j=jstart;
   diag=sqrt(2.0)*0.5;
   cos0[0]=1.0; sin0[0]=0.0;
   cos1[0]=1.0; sin1[0]=0.0;
   cos2[0]=1.0; sin2[0]=0.0;
   cos3[0]=1.0; sin3[0]=0.0;
   cos0[1]=-diag; sin0[1]=-diag;
   cos1[1]=diag;  sin1[1]=-diag;
   cos2[1]=-diag; sin2[1]=diag;
   cos3[1]=diag;  sin3[1]=diag;
   for(m=2; m<numMode; m++) {
     cos0[m]=cos0[m-1]*cos0[1]-sin0[m-1]*sin0[1];
     sin0[m]=sin0[m-1]*cos0[1]+cos0[m-1]*sin0[1];
     cos1[m]=cos1[m-1]*cos1[1]-sin1[m-1]*sin1[1];
     sin1[m]=sin1[m-1]*cos1[1]+cos1[m-1]*sin1[1];
     cos2[m]=cos2[m-1]*cos2[1]-sin2[m-1]*sin2[1];
     sin2[m]=sin2[m-1]*cos2[1]+cos2[m-1]*sin2[1];
     cos3[m]=cos3[m-1]*cos3[1]-sin3[m-1]*sin3[1];
     sin3[m]=sin3[m-1]*cos3[1]+cos3[m-1]*sin3[1];
   }
   for(i=istart; i<iend; i++)
     { 
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j].head[s]->pt;
         while(p) {
           z=p->z;  x=p->x; y=p->y;
           R=sqrt(x*x+y*y); invR=1.0/R;
           r=R-(j-jstart);

           coss[0]=1.0; sins[0]=0.0;
           coss[1]=x*invR; sins[1]=y*invR;
           for(m=2; m<numMode; m++) {
             coss[m]=coss[m-1]*coss[1]-sins[m-1]*sins[1];
             sins[m]=sins[m-1]*coss[1]+coss[m-1]*sins[1];
           }

           i1=((int)(i+z+0.5));
           j1=((int)(j+r+0.5));
           wz[1]=z;                     wz[0]=1.0-wz[1];
           WZ[1]=z+0.5-((int)(z+0.5));  WZ[0]=1.0-WZ[1];

           Bz=Br=Bp=Ez=Er=Ep=0.0;
           //Bz, Bp, Er
           if(j1==jstart) {
             WR[1]=R*2.0;                 WR[0]=1.0-WR[1];
	     Bz0=Bz1=Bz2=Bz3=0.0;
	     Bp0=Bp1=Bp2=Bp3=0.0;
	     Er0=Er1=Er2=Er3=0.0;
             m=0;
             for(ii=0; ii<2; ii++) {
               Bz0+=wz[ii]*D->BzNowR[m][i+ii][j1];
               Bp0+=WZ[ii]*D->BpNowR[m][i1-1+ii][j1];
               Er0+=wz[ii]*D->ErR[m][i+ii][j1];
             }
	     Bz1=Bz2=Bz3=Bz0;
	     Bp1=Bp2=Bp3=Bp0;
	     Er1=Er2=Er3=Er0;
             for(m=1; m<numMode; m++) 
               for(ii=0; ii<2; ii++) {
                 Bz0+=wz[ii]*D->BzNowR[m][i+ii][j1]*cos0[m]
                     -wz[ii]*D->BzNowR[m][i+ii][j1]*sin0[m];
                 Bz1+=wz[ii]*D->BzNowR[m][i+ii][j1]*cos1[m]
                     -wz[ii]*D->BzNowR[m][i+ii][j1]*sin1[m];
                 Bz2+=wz[ii]*D->BzNowR[m][i+ii][j1]*cos2[m]
                     -wz[ii]*D->BzNowR[m][i+ii][j1]*sin2[m];
                 Bz3+=wz[ii]*D->BzNowR[m][i+ii][j1]*cos3[m]
                     -wz[ii]*D->BzNowR[m][i+ii][j1]*sin3[m];
                 Bp0+=WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*cos0[m]
                     -WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*sin0[m];
                 Bp1+=WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*cos1[m]
                     -WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*sin1[m];
                 Bp2+=WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*cos2[m]
                     -WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*sin2[m];
                 Bp3+=WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*cos3[m]
                     -WZ[ii]*D->BpNowR[m][i1-1+ii][j1]*sin3[m];
                 Er0+=wz[ii]*D->ErR[m][i+ii][j1]*cos0[m]
                     -wz[ii]*D->ErR[m][i+ii][j1]*sin0[m];
                 Er1+=wz[ii]*D->ErR[m][i+ii][j1]*cos1[m]
                     -wz[ii]*D->ErR[m][i+ii][j1]*sin1[m];
                 Er2+=wz[ii]*D->ErR[m][i+ii][j1]*cos2[m]
                     -wz[ii]*D->ErR[m][i+ii][j1]*sin2[m];
                 Er3+=wz[ii]*D->ErR[m][i+ii][j1]*cos3[m]
                     -wz[ii]*D->ErR[m][i+ii][j1]*sin3[m];
               }
             Bz=(x+0.5)*(y+0.5)*Bz3+(0.5-x)*(y+0.5)*Bz2+(x+0.5)*(0.5-y)*Bz1*(0.5-x)*(0.5-y)*Bz0;
             Bp=(x+0.5)*(y+0.5)*Bp3+(0.5-x)*(y+0.5)*Bp2+(x+0.5)*(0.5-y)*Bp1*(0.5-x)*(0.5-y)*Bp0;
             Er=(x+0.5)*(y+0.5)*Er3+(0.5-x)*(y+0.5)*Er2+(x+0.5)*(0.5-y)*Er1*(0.5-x)*(0.5-y)*Er0;
           } else {
             WR[1]=r+0.5-((int)(r+0.5));  WR[0]=1.0-WR[1];

             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 Bz+=wz[ii]*WR[jj]*D->BzNowR[0][i+ii][j1-1+jj];
                 Bp+=WZ[ii]*WR[jj]*D->BpNowR[0][i1-1+ii][j1-1+jj];
                 Er+=wz[ii]*WR[jj]*D->ErR[0][i+ii][j1-1+jj];
               }
             for(m=1; m<numMode; m++) 
               for(ii=0; ii<2; ii++)
                 for(jj=0; jj<2; jj++) {
                   Bz+=wz[ii]*WR[jj]*D->BzNowR[m][i+ii][j1-1+jj]*coss[m];
                   Bz-=wz[ii]*WR[jj]*D->BzNowI[m][i+ii][j1-1+jj]*sins[m];
                   Bp+=WZ[ii]*WR[jj]*D->BpNowR[m][i1-1+ii][j1-1+jj]*coss[m];
                   Bp-=WZ[ii]*WR[jj]*D->BpNowI[m][i1-1+ii][j1-1+jj]*sins[m];
                   Er+=wz[ii]*WR[jj]*D->ErR[m][i+ii][j1-1+jj]*coss[m];
                   Er-=wz[ii]*WR[jj]*D->ErI[m][i+ii][j1-1+jj]*sins[m];
                 }             
           }

           //Ez, Ep, Br
           wr[1]=r;                  wr[0]=1.0-wr[1];

           for(ii=0; ii<2; ii++)
             for(jj=0; jj<2; jj++) {
               Br+=WZ[ii]*wr[jj]*D->BrNowR[0][i1-1+ii][j+jj];
               Ez+=WZ[ii]*wr[jj]*D->EzR[0][i1-1+ii][j+jj];
               Ep+=wz[ii]*wr[jj]*D->EpR[0][i+ii][j+jj];
             }

           for(m=1; m<numMode; m++) 
             for(ii=0; ii<2; ii++)
               for(jj=0; jj<2; jj++) {
                 Br+=WZ[ii]*wr[jj]*D->BrNowR[m][i1-1+ii][j+jj]*coss[m]
                    -WZ[ii]*wr[jj]*D->BrNowI[m][i1-1+ii][j+jj]*sins[m];
                 Ez+=WZ[ii]*wr[jj]*D->EzR[m][i1-1+ii][j+jj]*coss[m]
                    -WZ[ii]*wr[jj]*D->EzI[m][i1-1+ii][j+jj]*sins[m];
                 Ep+=wz[ii]*wr[jj]*D->EpR[m][i+ii][j+jj]*coss[m]
                    -wz[ii]*wr[jj]*D->EpI[m][i+ii][j+jj]*sins[m];
               }
           
           Ex=Er*coss[1]-Ep*sins[1];
           Ey=Er*sins[1]+Ep*coss[1];
           Bx=Br*coss[1]-Bp*sins[1];
           By=Br*sins[1]+Bp*coss[1];

           p->Ez=Ez; p->Ex=Ex; p->Ey=Ey;
           p->Bz=Bz; p->Bx=Bx; p->By=By;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)

}

/*
void interpolation2D_Yee_Pukhov_2nd(Domain *D,External *Ext)
{
   int s,i,j,k,ii,jj,i,j+1,k1,istart,iend,jstart,jend,kstart,kend;
   double E1,E2,E3,B1,B2,B3,extE1,extE2,extE3,extB1,extB2,extB3;
   double x,y,z,xx,yy,zz,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;
   
   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           //edge
           ii=(int)(i+x+0.5);           jj=(int)(j+y+0.5);
           xx=1.0+x-((int)(0.5+x));     yy=1.0+y-((int)(0.5+y));
           //side
           i=i;           j+1=j;
           x1=0.5+x;       y1=0.5+y;

           E1=calBi+1D(D->Ex,i,jj,x1,yy);
           E2=calBi+1D(D->Ey,ii,j+1,xx,y1);
           E3=calBi+1D(D->Ez,ii,jj,xx,yy);
           B1=calBi+1D(D->BxNow,ii,j+1,xx,y1);
           B2=calBi+1D(D->ByNow,i,jj,x1,yy);
           B3=calBi+1D(D->BzNow,i,j+1,x1,y1);

           p->E1=E1+extE1; p->E2=E2+extE2; p->E3=E3+extE3;
           p->B1=B1+extB1; p->B2=B2+extB2; p->B3=B3+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

void interpolation1D_Split_1st(Domain *D,External *Ext)
{
   int i,j,k,i,istart,iend,s,cnt;
   double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3,x,x1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart; iend=D->iend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;
   
   j=k=0;
   for(i=istart; i<iend; i++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x; 
           i=(int)(i+x+0.5);
           x1=x+0.5-((int)(x+0.5));

           B1=0;
           E1=(1-x1)*D->Ex[i-1][j][k] + x1*D->Ex[i][j][k];
           Pr=(1-x1)*D->Pr[i-1][j][k] + x1*D->Pr[i][j][k];
           Pl=(1-x1)*D->Pl[i-1][j][k] + x1*D->Pl[i][j][k];
           Sr=(1-x1)*D->Sr[i-1][j][k] + x1*D->Sr[i][j][k];
           Sl=(1-x1)*D->Sl[i-1][j][k] + x1*D->Sl[i][j][k];

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

void interpolation2D_Split_1st(Domain *D,External *Ext)
{
   int i,j,k,i,j+1,k1,istart,iend,jstart,jend,s,cnt;
   double E1,Pr,Pl,B1,Sr,Sl,x,y,z,x1,y1,z1;
   double extE1,extE2,extE3,extB1,extB2,extB3;
   ptclList *p;
   int myrank, nprocs;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;

   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         cnt=0;
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           i=(int)(i+x+0.5);           j+1=(int)(j+y+0.5);
           x1=x+0.5-((int)(x+0.5));     y1=y+0.5-((int)(y+0.5));

           B1=(1-x1)*(1-y1)*D->Bx[i-1][j][k1]
             +    x1*(1-y1)*D->Bx[i][j][k1]
             +(1-x1)*    y1*D->Bx[i-1][j+1][k1]
             +    x1*    y1*D->Bx[i][j+1][k1];
           E1=(1-x1)*(1-y)*D->Ex[i-1][j][k]
             +    x1*(1-y)*D->Ex[i][j][k]
             +(1-x1)*    y*D->Ex[i-1][j+1][k]
             +    x1*    y*D->Ex[i][j+1][k];
           Pr=(1-x1)*(1-y1)*D->Pr[i-1][j][k]
             +    x1*(1-y1)*D->Pr[i][j][k]
             +(1-x1)*    y1*D->Pr[i-1][j+1][k]
             +    x1*    y1*D->Pr[i][j+1][k];
           Pl=(1-x1)*(1-y1)*D->Pl[i-1][j][k]
             +    x1*(1-y1)*D->Pl[i][j][k]
             +(1-x1)*    y1*D->Pl[i-1][j+1][k]
             +    x1*    y1*D->Pl[i][j+1][k];
           Sr=(1-x1)*(1-y)*D->Sr[i-1][j][k1]
             +    x1*(1-y)*D->Sr[i][j][k1]
             +(1-x1)*    y*D->Sr[i-1][j+1][k1]
             +    x1*    y*D->Sr[i][j+1][k1];
           Sl=(1-x1)*(1-y)*D->Sl[i-1][j][k1]
             +    x1*(1-y)*D->Sl[i][j][k1]
             +(1-x1)*    y*D->Sl[i-1][j+1][k1]
             +    x1*    y*D->Sl[i][j+1][k1];

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
           cnt++;
         }
       }                //for(s)        
     }             //for(i,j)
}

void interpolation2D_Split_2nd(Domain *D,External *Ext)
{
   int s,i,j,k,ii,jj,i,j+1,k1,istart,iend,jstart,jend;
   double E1,Pr,Pl,B1,Sr,Sl,extE1,extE2,extE3,extB1,extB2,extB3;
   double x,y,z,xx,yy,zz,x1,y1,z1;
   ptclList *p;
   int myrank, nprocs;    

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;   iend=D->iend;
   jstart=D->jstart;   jend=D->jend;

   Particle ***particle;
   particle=D->particle;

   extE1=Ext->E1;   extE2=Ext->E2;   extE3=Ext->E3;
   extB1=Ext->B1;   extB2=Ext->B2;   extB3=Ext->B3;
   
   k=k1=0;
   for(i=istart; i<iend; i++)
     for(j=jstart; j<jend; j++)
     {
       for(s=0; s<D->nSpecies; s++)
       {
         p=particle[i][j][k].head[s]->pt;
         while(p)
         {
           x=p->x;  y=p->y; // z=p->z;
           //edge
           ii=(int)(i+x+0.5);           jj=(int)(j+y+0.5);
           xx=1.0+x-((int)(0.5+x));     yy=1.0+y-((int)(0.5+y));
           //side
           i=i;           j+1=j;
           x1=0.5+x;       y1=0.5+y;

           E1=calBi+1D(D->Ex,ii,jj,xx,yy);
           Pr=calBi+1D(D->Pr,i,j+1,x1,y1);
           Pl=calBi+1D(D->Pl,i,j+1,x1,y1);
           B1=calBi+1D(D->Bx,ii,jj,xx,yy);
           Sr=calBi+1D(D->Sr,i,j+1,x1,y1);
           Sl=calBi+1D(D->Sl,i,j+1,x1,y1);

           p->E1=E1+extE1; p->E2=Pr+Pl+extE2; p->E3=Sr+Sl+extE3;
           p->B1=B1+extB1; p->B2=Sl-Sr+extB2; p->B3=Pr-Pl+extB3;

           p=p->next;
         }
       }		//for(s)        
     }		   //for(i,j)
}

*/

