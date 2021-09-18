#include <stdio.h>
#include <stdlib.h>
#include "mesh.h"
#include "constants.h"
#include <math.h>


void particlePush(Domain *D,int iteration)
{
    int i,j,istart,iend,jstart,jend,l,m,s,minRSub;
    double x,y,z,R,shiftZ,shiftX,shiftY,dt,dz,dr,sqrT,coef,coef1;
    double dtOverdz,dtOverdr,gamma,invGamma,numeric=1e-6;
    double pMinus[3],T[3],S[3],operate[3][3],pPlus[3];
    Particle **particle;
    particle=D->particle;
    LoadList *LL;
    ptclList *p, *New, *tmp, *prev;

    istart=D->istart;    iend=D->iend;
    jstart=D->jstart;    jend=D->jend;
    minRSub=D->minYSub;

    for(i=0; i<3; i++)   {
       pMinus[i]=0.0;
       T[i]=0.0;
       S[i]=0.0;
       pPlus[i]=0.0;
    }

    dt=D->dt; dz=D->dz; dr=D->dr;
    dtOverdz=dt/dz; dtOverdr=dt/dr;

    double mass[D->nSpecies];
    LL=D->loadList;
    s=0;

    while(LL->next)
    {
       mass[s]=LL->mass;
//       charge[s]=LL->charge;
       s++;
       LL=LL->next;
    }

    for(i=istart; i<iend; i++)
      for(j=jstart; j<jend; j++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            coef1=M_PI/mass[s]*dt;
            p=particle[i][j].head[s]->pt;     
            while(p)
            {    
              coef=coef1*p->charge;
              //Calculate vector P- 
              pMinus[0]=p->pz+coef*(p->Ez);
              pMinus[1]=p->px+coef*(p->Ex);
              pMinus[2]=p->py+coef*(p->Ey);    
 
             //Calculate vector T 
             invGamma=1.0/sqrt(1.0+pMinus[0]*pMinus[0]+pMinus[1]*pMinus[1]+pMinus[2]*pMinus[2]);           
             T[0]=coef*invGamma*(p->Bz);
             T[1]=coef*invGamma*(p->Bx);   
             T[2]=coef*invGamma*(p->By);

             //Calculate vector S
             sqrT=1.0+T[0]*T[0]+T[1]*T[1]+T[2]*T[2];
             for(l=0; l<3; l++)  
                S[l]=2.0*T[l]/sqrT;
  
             //Calculate operator A from P+=A.P-
             operate[0][0]=1.0-S[2]*T[2]-S[1]*T[1];      
             operate[0][1]=S[1]*T[0]+S[2];    
             operate[0][2]=S[2]*T[0]-S[1];     
             operate[1][0]=S[0]*T[1]-S[2];        
             operate[1][1]=1.0-S[0]*T[0]-S[2]*T[2];          
             operate[1][2]=S[2]*T[1]+S[0];         
             operate[2][0]=S[0]*T[2]+S[1];    
             operate[2][1]=S[1]*T[2]-S[0];    
             operate[2][2]=1.0-S[0]*T[0]-S[1]*T[1]; 
             //Calculate vector P+
             for(l=0; l<3; l++)  {
                pPlus[l]=0.0;
                for(m=0; m<3; m++)   
                   pPlus[l]+=operate[l][m]*pMinus[m];
                }
             //Updated momentum              
             p->pz=pPlus[0]+coef*(p->Ez); 
             p->px=pPlus[1]+coef*(p->Ex); 
             p->py=pPlus[2]+coef*(p->Ey); 
    
             //Translation
             gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
             shiftZ=p->pz/gamma*dtOverdz;
             shiftX=p->px/gamma*dtOverdr;
             shiftY=p->py/gamma*dtOverdr;

if(shiftZ>=1.0 || shiftX>=1.0 || shiftY>=1.0) { 
  printf("shiftZ=%g, shiftX=%g, shiftY=%g\n",shiftZ,shiftX,shiftY);
  exit(0);
} else ;
             p->oldZ=i+p->z; p->z+=shiftZ;
             p->oldX=p->x; p->x+=shiftX;
             p->oldY=p->y; p->y+=shiftY;
             
             p=p->next;
          }		//End of while(p)
        }		//End of for(s)
      }      	//End of for(i,j)
/*
    j=jstart;
    for(i=istart; i<iend; i++)
        {
          for(s=0; s<D->nSpecies; s++)
          {
            coef1=M_PI/mass[s]*dt;
            p=particle[i][j].head[s]->pt;     
            while(p)
            {    
              coef=coef1*p->charge;
              //Calculate vector P- 
              pMinus[0]=p->pz+coef*(p->Ez);
              pMinus[1]=p->px+coef*(p->Ex);
              pMinus[2]=p->py+coef*(p->Ey);    
 
             //Calculate vector T 
             invGamma=1.0/sqrt(1.0+pMinus[0]*pMinus[0]+pMinus[1]*pMinus[1]+pMinus[2]*pMinus[2]);           
             T[0]=coef*invGamma*(p->Bz);
             T[1]=coef*invGamma*(p->Bx);   
             T[2]=coef*invGamma*(p->By);

             //Calculate vector S
             sqrT=1.0+T[0]*T[0]+T[1]*T[1]+T[2]*T[2];
             for(l=0; l<3; l++)  
                S[l]=2.0*T[l]/sqrT;
  
             //Calculate operator A from P+=A.P-
             operate[0][0]=1.0-S[2]*T[2]-S[1]*T[1];      
             operate[0][1]=S[1]*T[0]+S[2];    
             operate[0][2]=S[2]*T[0]-S[1];     
             operate[1][0]=S[0]*T[1]-S[2];        
             operate[1][1]=1.0-S[0]*T[0]-S[2]*T[2];          
             operate[1][2]=S[2]*T[1]+S[0];         
             operate[2][0]=S[0]*T[2]+S[1];    
             operate[2][1]=S[1]*T[2]-S[0];    
             operate[2][2]=1.0-S[0]*T[0]-S[1]*T[1]; 
             //Calculate vector P+
             for(l=0; l<3; l++)  {
                pPlus[l]=0.0;
                for(m=0; m<3; m++)   
                   pPlus[l]+=operate[l][m]*pMinus[m];
                }
             //Updated momentum              
             p->pz=pPlus[0]+coef*(p->Ez); 
             p->px=pPlus[1]+coef*(p->Ex); 
             p->py=pPlus[2]+coef*(p->Ey); 
    
             //Translation
             gamma=sqrt(1.0+p->pz*p->pz+p->px*p->px+p->py*p->py);
             shiftZ=p->pz/gamma*dtOverdz;
             shiftX=p->px/gamma*dtOverdr;
             shiftY=p->py/gamma*dtOverdr;

//if(shiftZ>=1.0 || shiftX>=1.0 || shiftY>=1.0) 
//  printf("shiftZ=%g, shiftX=%g, shiftY=%g\n",shiftZ,shiftX,shiftY);
             p->oldZ=i+p->z; p->z+=shiftZ;
             p->oldX=p->x; p->x+=shiftX;
             p->oldY=p->y; p->y+=shiftY;

//             R=sqrt(p->x*p->x+p->y*p->y);
//             if(R==0.0) { p->x+=numeric; p->y+=numeric; } else ;
//             if(R==0.0) { printf("x=%g, y=%g, i=%d, j=%d\n",p->x,p->y,i,j); } else ;
             
             p=p->next;
          }		//End of while(p)
        }		//End of for(s)
      }      	//End of for(i,j)
*/		
}

             
