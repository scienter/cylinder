#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>
#include <mpi.h>

double maxwellianVelocity(double temterature); 

void loadPlasma_crystal(Domain *D,LoadList *LL,int s)
{
   int i,j,k,l,istart,iend,intNum,cnt,np,nc,leftIndex,rightIndex;
   double space,positionX,x,n0,n1,nc1,xL,weight;
   double wp,pDt,v1,v2,v3,gamma,mass,refX;
   Particle ***particle;
   particle=D->particle;
   ptclList *New,*p;   
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   weight=1.0/LL->numberInCell;

   j=k=0;
   for(l=0; l<LL->xnodes-1; l++)
   {
     if(LL->xn[l+1]-LL->xn[l]>0 || LL->xn[l+1]-LL->xn[l]<0)
     {
       n0=LL->xn[l];
       n1=LL->xn[l+1];
       nc=LL->numberInCell;
       xL=LL->xpoint[l+1]-LL->xpoint[l];

       np=((int)((n1+n0)*xL*0.5*nc));
       cnt=0;
       x=(sqrt(n0*n0*(np-1)*(np-1)+(n1+n0)*(n1-n0)*cnt*(np-1))-n0*(np-1))/(n1-n0)/(np-1)*xL+LL->xpoint[l]-D->minXSub;
       i=(int)(x);
       refX=0;

       p = particle[i][j][k].head[s]->pt;
       while(p)  {
         if(refX<=p->x)
           refX=p->x;
         p=p->next;
       }
       space=1.0/LL->numberInCell*LL->xn[l];
       refX+=space;
       refX-=(int)refX;     
       while(cnt<np)
       {
          x=(sqrt(n0*n0*(np-1)*(np-1)+(n1+n0)*(n1-n0)*cnt*(np-1))-n0*(np-1))/(n1-n0)/(np-1)*xL+LL->xpoint[l]-D->minXSub+refX;
          i=(int)(x);
          if(i>=istart && i<iend)
          {
             positionX=x-i;
             New = (ptclList *)malloc(sizeof(ptclList)); 
             New->next = particle[i][j][k].head[s]->pt;
             particle[i][j][k].head[s]->pt = New;

             New->x = positionX; New->y=0.0; New->z=0.0;
             New->oldX=i+positionX;
             New->E1=New->E2=New->E3=New->B1=New->B2=New->B3=0.0;
             v1=maxwellianVelocity(LL->temperature)/velocityC;
             v2=maxwellianVelocity(LL->temperature)/velocityC;
             v3=maxwellianVelocity(LL->temperature)/velocityC;
             New->p1=-D->gamma*D->beta+v1;
             New->p2=v2;
             New->p3=v3;
             LL->index++;
             New->index=LL->index;            
             New->core=myrank;            
             New->weight=weight;            
          }
          cnt++;
       }	//end of while(cnt)
     }	//end if (l)

     else
     {
        nc=LL->numberInCell*LL->xn[l];
        space=1.0/((double)nc);
        leftIndex=(int)(LL->xpoint[l]-D->minXSub);
        rightIndex=(int)(LL->xpoint[l+1]-D->minXSub);
        if(rightIndex>iend-1)
           rightIndex=iend; 
        np=(int)((rightIndex-leftIndex)*nc);
        cnt=0;
        refX=0;
        x=leftIndex+space*(cnt+0.5);

        i=(int)x-1;
        if(i>=istart && i<iend)
        {
          p = particle[i][j][k].head[s]->pt;
          while(p)  {
            if(refX<=p->x)
              refX=p->x;
             p=p->next;
          }
          refX+=space*1.5;
          refX=fabs(1-refX);
          if(refX==space) { 
             cnt=1;
          }
        }

        while(cnt<np)
        {
          x=leftIndex+space*(cnt+0.5)+refX;
          i=((int)(x));
          if(i>=istart && i<iend)
          {
             positionX=x-i;
              
             New = (ptclList *)malloc(sizeof(ptclList)); 
             New->next = particle[i][j][k].head[s]->pt;
             particle[i][j][k].head[s]->pt = New;

             New->x = positionX; New->y=0.0; New->z=0.0;
             New->oldX=i+positionX;
             New->E1=New->E2=New->E3=New->B1=New->B2=New->B3=0.0;
             v1=maxwellianVelocity(LL->temperature)/velocityC;
             v2=maxwellianVelocity(LL->temperature)/velocityC;
             v3=maxwellianVelocity(LL->temperature)/velocityC;
             New->p1=-D->gamma*D->beta+v1;
             New->p2=v2;
             New->p3=v3;
             LL->index++;
             New->index=LL->index;            
             New->core=myrank;
             New->weight=weight;            
          }
           cnt++;
        }

     }	//end of else

   }

}




void loadMovingPlasma_crystal(Domain *D,LoadList *LL,int s)
{
   int i,j,k,l,intNum,cnt,np,istart,iend;
   double space,position,positionX,x,ne,minX;
   double leftIndex,rightIndex,nc;
   double v1,v2,v3,gamma,mass,weight;
   Particle ***particle;
   particle=D->particle;
   ptclList *New,*p;
   int myrank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

   istart=D->istart;
   iend=D->iend;
   weight=1.0/LL->numberInCell;

   i=iend-1;
   j=k=0;
//   while(LL->next)
//   {  
     for(l=0; l<LL->xnodes-1; l++)
     {
       position=i+D->minXSub;
 
       if(position>=LL->xpoint[l] && position<LL->xpoint[l+1])
       {
         nc=((LL->xn[l+1]-LL->xn[l])/(LL->xpoint[l+1]-LL->xpoint[l])
            *(position-LL->xpoint[l])+LL->xn[l]);
         nc*=LL->numberInCell;	//it is the double number of superparticles.
         space=1.0/nc;
         p=particle[iend-1][j][k].head[s]->pt;
         x=-2;
         while(p)
         {
           if(p->x>x)  x=p->x;
           p=p->next;
         }
   
         while(x<1)
         {
           x+=space;
           positionX=x;

           New = (ptclList *)malloc(sizeof(ptclList)); 
           New->next = particle[i][j][k].head[s]->pt;
           particle[i][j][k].head[s]->pt = New;

           New->x = positionX; New->y=0.0; New->z=0.0;
           New->oldX=i+positionX;
           New->E1=New->E2=New->E3=New->B1=New->B2=New->B3=0.0;
           v1=maxwellianVelocity(LL->temperature)/velocityC;
           v2=maxwellianVelocity(LL->temperature)/velocityC;
           v3=maxwellianVelocity(LL->temperature)/velocityC;
//           New->p1=-D->gamma*D->beta+v1;
           New->p1=v1;
           New->p2=v2;
           New->p3=v3;
           LL->index++;
           New->index=LL->index;            
           New->core=myrank;
           New->weight=weight;
         }
       }
     }
//     LL=LL->next;
//     s++;
//   }		//End of while(LL)  
}

