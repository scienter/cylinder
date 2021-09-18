#include "particle.h"
#include "laser.h"
#include "plasma.h"

#define FIRST 	1
#define SECOND 	2
#define THIRD 	3

#define ON	1
#define OFF	0
#define File 2

#define TXT	0
#define HDF	1

#define Split	0
#define Yee	1
#define NoCherenkov	2

#define UP		1
#define DOWN		2
#define FRONT		3
#define BACK		4
#define UPFRONT		5
#define UPBACK		6
#define DOWNFRONT	7
#define DOWNBACK	8

#define Lifschitz	1
#define Davidson	2

typedef struct _Domain 
{
   int dimension;
   int filter,filterIter;

   int fieldType;
   int currentType,currentCons;
   int interpolationType;
   double cenComp;

   int maxTime;
   int saveFieldMode;
   int saveParticleMode;
   int saveDensityMode;
   int saveCurrentMode;
   int saveDumpMode;
   int maxStep;
   int saveStep;
   int centerStep;
   int saveStart;
   int dumpStart;
   int dumpSave;
   int dumpSaveStep;
   int dumpStep;
   int fieldSave;
   int ramanSave;
   int particleSave;
   int densitySave;
   int currentSave;
   double centerPick,centerAngle;

   int nx;             //Total domain
   int ny;             //Total domain
   int nxSub;          //Each core has sub domain        
   int nySub;          //Each core has sub domain        
   int istart;
   int iend;
   int jstart;
   int jend;
   //Each core has start mesh point in total domain
	double minX,maxX;
   int minXSub,maxXSub,minYSub,maxYSub;
   int minXDomain,minYDomain,maxXDomain;
   int numberInCell;
   int moving,shiftStart,moveIt;         //Moving domain option. 1:on
   double movingV;
//   int shiftDuration;

   double lambda;
   double omega;
   double divisionLambda;
   double dt;
   double dtRatio;
   double dz;
   double dr;
   int resolChange;
   int resolHigh;
   int resolLow;
   int resolStep;
   int resolX;
   int resolY;
    

   //MPI parameter
   int L;
   int M;
   int nextXrank;
   int prevXrank;
   int nextYrank;
   int prevYrank;
   
   //sharing mesh
   double *XplusJ,*XminusJ,*YplusJ,*YminusJ;
   double *minusDenY,*plusDenY,*minusDenZ,*plusDenZ;
   int numPlusXJ,numMinusXJ,numPlusYJ,numMinusYJ;
   int numPlusDenY,numMinusDenY;

   int numMode;
   double dF;
   double ***FR,***FI,***CnR,***CnI;
   //Yee
   double ***RhoNoPairR,***RhoNoPairI,***RhoPairR,***RhoPairI;
   double ***EzR,***ErR,***EpR,***BzR,***BrR,***BpR;    
   double ***EzI,***ErI,***EpI,***BzI,***BrI,***BpI;    
   double ***JzR,***JrR,***JpR,***JzI,***JrI,***JpI;    
   double ***BzNowR,***BrNowR,***BpNowR,***BzNowI,***BrNowI,***BpNowI;
   //NDFX
   double ***PrCR,***PlCR,***SrCR,***SlCR,***EzCR,***BzCR;
   double ***PrCI,***PlCI,***SrCI,***SlCI,***EzCI,***BzCI;
   double ***PrR,***PlR,***SrR,***SlR;
   double ***PrI,***PlI,***SrI,***SlI;
   double ***JzCR,***JrCR,***JpCR,***JzCI,***JrCI,***JpCI;
   
   struct _Particle **particle;    
   struct _Boost **boost;    

   //Plasma load
   struct _LoadList *loadList;
   int nSpecies;
	double index;

   //Laser load
   struct _LaserList *laserList;
   int nLaser;
   double *laserI,*laserPhase;

   //Boost
   int boostOn;
   int boostIon;
   double gamma;
   double beta;
   int minT;	//boost frame's step
   int maxT;	//boost frame's step
   int boostSaveStep;	//lab frame's step
   
   //Probe
   int probeNum;
   int *probeX;
   int *probeY;
   int *probeZ;
   struct _Probe **probe;
   
   //ID Track
   int tracking;
   int trackSaveStep;
   int trackStart;
   int idNums;
   int *trackID;
   int *trackCore;
   int *trackS;
   struct _Track **track;

   //PML
   int pml;
   int pmlStart;
   int pmlCellRight,pmlCellLeft;   
   int pmlCellUp;   
//   struct _PML ***upml,***lpml; 
   double pmlr,pmld;
   double *upr,*upd,*rtr,*rtd,*ltr,*ltd;
   int period;

   //Field ionization
   int fieldIonization;

   //Particle redistributing
   int redist;
   int redistStep;

}  Domain; 

typedef struct _Boost
{
   double x;
   double y;
   double E1;
   double B1;
   double Pr;
   double Pl;
   double Sr;
   double Sl;   
}  Boost;

typedef struct _PML 
{
   double EzR,ErR,EpzR,EprR;
   double BzR,BrR,BpzR,BprR;
   double EzI,ErI,EpzI,EprI;
   double BzI,BrI,BpzI,BprI;
}  PML;

typedef struct _Particle 
{
   // Particle List Header
   ptclHead **head;            
}  Particle;

typedef struct _External 
{
   double E1;
   double E2;
   double E3;
   double B1;
   double B2;
   double B3;
}  External;

typedef struct _Probe
{
   double E1;
   double Pr;
   double Pl;
   double B1;
   double Sr;
   double Sl;
}  Probe;

typedef struct _Track
{
   double x;
   double y;
   double z;
   double px;
   double py;
   double pz;
   int step;
   int id;
   int core;
   double wp;
   double kp;
}  Track;

void cleanMemory(Domain *D);
void tmpCleanMemory(Domain *D);
void saveTracking(Domain *D);
void removeEdge(Domain *D);
void particleShareZ(Domain *D);
void particleShareY(Domain D);
void particleShareX(Domain D);
void rearrangeParticles(Domain *D);
void movingDomain(Domain *D,int iteration);
void updateCurrent(Domain D,int iteration);
void particlePush(Domain *D,int iteration);
void interpolation(Domain *D,External *Ext,int iteration);
void fieldSolve1(Domain D,double t,int iteration,double dF);
void fieldSolve2(Domain D,double t,int iteration,double dF);
void loadLaser(Domain *D,LaserList *L,double t);
void saveDump(Domain D,int iteration);
void saveBDump(Domain D,int iteration);
void saveEDump(Domain D,int iteration);
void saveJDump(Domain D,int iteration);
void saveDumpParticleResolHDF(Domain *D,int iteration);
void saveDumpDensityResolHDF(Domain D,int iteration);
void saveP_GridHDF(Domain D,int iteration);
void saveFile(Domain D,int iteration);
void firstFieldShare(Domain D);
void secondFieldShare(Domain D);
void trackID(Domain *D,int iteration,int istart,int iend,int jstart,int jend,int kstart,int kend);
void loadPlasma(Domain *D,LoadList *LL,int s,int iteration,int istart,int iend,int jstart,int jend);
void loadBeam(Domain *D,LoadList *LL,int s,int iteration);
void restoreDump(Domain D,int iteration);
void boundary(Domain *D,External *Ext);
void reBoundary(Domain *D,External *Ext);
void parameterSetting(Domain *D,External *Ext, char *input);
int FindParameters (char *block, int rank, char *options, char *input, char *ret);
void saveFieldHDF(Domain D,int iteration);
void saveCenterDensity(Domain *D,int iteration);
void saveCenterField(Domain *D,int iteration);
double ***memoryAsign(int nx, int ny, int nz);
void restoreIntMeta(char *fileName,char *dataName,int *data,int dataCnt);
void solveF(Domain D);
void solveCharge(Domain *D,LoadList *LL,double ***rhoR,double ***rhoI,int istart,int iend,int jstart,int jend,int s,double coef);
void movingPairCharge(Domain *D);
void filter(Domain *D,double ***dataR,double ***dataI);
void filterFieldC(Domain *D);
void ionizationSetup(LoadList *LL,int species);
void fieldIonization(Domain *D,int iteration);
double randomValue(double beta);
void particle_redist(Domain *D,int iteration,External *Ext);
void MPI_Transfer2F_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share);
void MPI_Transfer2F_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share);
void MPI_Transfer4F_NDFX_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share);
void MPI_Transfer4F_NDFX_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share);
void MPI_Transfer4F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share);
void MPI_Transfer4F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,int ny,int share);
void MPI_Transfer6F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share);
void MPI_Transfer6F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,int ny,int share);
void MPI_Transfer8F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share);
void MPI_Transfer8F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,int ny,int share);
void MPI_Transfer12F_Xminus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share);
void MPI_Transfer12F_Xplus(Domain *D,double ***f1,double ***f2,double ***f3,double ***f4,double ***f5,double ***f6,double ***f7,double ***f8,double ***f9,double ***f10,double ***f11,double ***f12,int ny,int share);
void MPI_TransferDen_Xplus(Domain *D,double ***f1,double ***f2,int ny,int share);
void MPI_TransferDen_Xminus(Domain *D,double ***f1,double ***f2,int ny,int share);
void calConservation(Domain D,int iteration);
void filter_current(Domain *D,double ***val,double ***J,int iter);
void MPI_filter_Xminus(Domain *D,double ***f1,int ny,int share,int m);
void MPI_filter_Xplus(Domain *D,double ***f1,int ny,int share,int m);
void saveDensityProfile(Domain *D);
void checkEnergyCons(Domain *D,int iteration);
