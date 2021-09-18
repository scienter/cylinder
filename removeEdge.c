#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "mpi.h"
#include "constants.h"

void removeEdge(Domain *D)
{
    int m,i,j,istart,iend,jstart,jend,s,cnt,numMode;
    double r,phi;
    Particle **particle;    
    particle=D->particle;     
    ptclList *p,*tmp;

    int myrank,nTasks;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &nTasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    numMode=D->numMode;

    //remove left boundary
    for(i=0; i<istart; i++)
      for(j=jstart-1; j<jend+1; j++) 
        for(s=0; s<D->nSpecies; s++)
        {
          p=particle[i][j].head[s]->pt;
          while(p!=NULL)
          {
            particle[i][j].head[s]->pt=p->next; 
            p->next=NULL;
            free(p);
            p=particle[i][j].head[s]->pt;
          }
        }
    if(myrank==0) {
      for(m=0; m<numMode; m++)
          for(j=0; j<jend+3; j++) {
            D->RhoNoPairR[m][istart][j]=D->RhoNoPairR[m][istart+1][j];
            D->RhoNoPairI[m][istart][j]=D->RhoNoPairI[m][istart+1][j];
            D->RhoPairR[m][istart][j]=D->RhoPairR[m][istart+1][j];
            D->RhoPairI[m][istart][j]=D->RhoPairI[m][istart+1][j];
          }
    } else ;

    //remove right boundary
    i=iend;
    for(j=jstart-1; j<=jend; j++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          particle[i][j].head[s]->pt=p->next; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }

    //remove plusY boundary
    j=jend;
    for(i=0; i<=iend; i++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][j].head[s]->pt;
        while(p)
        {	
          particle[i][j].head[s]->pt=p->next; 
          p->next=NULL;
          free(p);
          p=particle[i][j].head[s]->pt;
        }
      }
/*
    //remove minusY boundary
    for(i=0; i<=iend; i++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][jstart].head[s]->pt; 
        while(p)
        {
          r=p->r; phi=p->phi; 
          if(r<0.5)  {
            p->r=1.0-r;
            p->phi=phi+pi;
//            p->pphi*=-1.0;
//            p->pr*=-1.0;
          } else ;
          p=p->next;
        }
      }
    for(i=0; i<=iend; i++)
      for(s=0; s<D->nSpecies; s++)
      { 
        p=particle[i][jstart-1].head[s]->pt; 
        while(p)
        {
printf("Oh,no! particles are beneth axis\n");

//          r=p->r; phi=p->phi; 
//          if(r<0.5)  {
//            p->r=1.0-r;
//            p->phi=phi+pi;
//          } else ;

          p=p->next;

        }
      }
*/
}

