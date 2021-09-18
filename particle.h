

typedef struct _ptclHead  {
    struct _ptclList *pt;
}   ptclHead;

typedef struct _ptclList  {
    double z; 
    double oldZ;   
    double x; 
    double oldX;   
    double y; 
    double oldY;   
    double pz;    //momentum  
    double px;
    double py;
    double Ez;    
    double Ex;    
    double Ey;    
    double Bz;    
    double Bx;    
    double By;    
    double weight;    
    double charge;
    int index; 
    int core; 
    struct _ptclList *next;
} ptclList;

