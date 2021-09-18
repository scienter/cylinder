#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "plasma.h"
#include "constants.h"



int whatSpecies(char *str)
{
   if(strstr(str,"Electron")) 		return Electron;
   else if(strstr(str,"HPlus0"))   	return HPlus0;
   else if(strstr(str,"HPlus1"))   	return HPlus1;
   else if(strstr(str,"HePlus0"))   	return HePlus0;
   else if(strstr(str,"HePlus1"))   	return HePlus1;
   else if(strstr(str,"HePlus2"))   	return HePlus1;
   else if(strstr(str,"CPlus0"))   	return CPlus0;
   else if(strstr(str,"CPlus1"))   	return CPlus1;
   else if(strstr(str,"CPlus2"))   	return CPlus2;
   else if(strstr(str,"CPlus3"))   	return CPlus3;
   else if(strstr(str,"CPlus4"))   	return CPlus4;
   else if(strstr(str,"CPlus5"))   	return CPlus5;
   else if(strstr(str,"CPlus6"))   	return CPlus6;
   else if(strstr(str,"NPlus0"))    return NPlus0;
   else if(strstr(str,"NPlus1"))    return NPlus1;
   else if(strstr(str,"NPlus2"))    return NPlus2;
   else if(strstr(str,"NPlus3"))    return NPlus3;
   else if(strstr(str,"NPlus4"))    return NPlus4;
   else if(strstr(str,"NPlus5"))    return NPlus5;
   else if(strstr(str,"NPlus6"))    return NPlus6;
   else if(strstr(str,"NPlus7"))    return NPlus7;
   else if(strstr(str,"NePlus0"))    return NePlus0;
   else if(strstr(str,"NePlus1"))    return NePlus1;
   else if(strstr(str,"NePlus2"))    return NePlus2;
   else if(strstr(str,"NePlus3"))    return NePlus3;
   else if(strstr(str,"NePlus4"))    return NePlus4;
   else if(strstr(str,"NePlus5"))    return NePlus5;
   else if(strstr(str,"NePlus6"))    return NePlus6;
   else if(strstr(str,"NePlus7"))    return NePlus7;
   else if(strstr(str,"NePlus8"))    return NePlus8;
   else if(strstr(str,"NePlus9"))    return NePlus9;
   else if(strstr(str,"NePlus10"))    return NePlus10;
   else if(strstr(str,"AlPlus4"))   	return AlPlus4;
   else return 0;
}

double whatMass(int species)
{
   if(species == Electron) 		return 1;
   else if(species == HPlus0)  		return 1.00794/eMassU;
   else if(species == HPlus1)  		return (1.00794-1*eMassU)/eMassU;
   else if(species == HePlus0)  	return (4.00260-0*eMassU)/eMassU;
   else if(species == HePlus1)  	return (4.00260-1*eMassU)/eMassU;
   else if(species == HePlus2)  	return (4.00260-2*eMassU)/eMassU;
   else if(species == CPlus0)           return (12.0111-0*eMassU)/eMassU;
   else if(species == CPlus1)           return (12.0111-1*eMassU)/eMassU;
   else if(species == CPlus2)           return (12.0111-2*eMassU)/eMassU;
   else if(species == CPlus3)           return (12.0111-3*eMassU)/eMassU;
   else if(species == CPlus4)           return (12.0111-4*eMassU)/eMassU;
   else if(species == CPlus5)           return (12.0111-5*eMassU)/eMassU;
   else if(species == CPlus6)           return (12.0111-6*eMassU)/eMassU;
   else if(species == NPlus0)           return (14.0064-0*eMassU)/eMassU;
   else if(species == NPlus1)           return (14.0064-1*eMassU)/eMassU;
   else if(species == NPlus2)           return (14.0064-2*eMassU)/eMassU;
   else if(species == NPlus3)           return (14.0064-3*eMassU)/eMassU;
   else if(species == NPlus4)           return (14.0064-4*eMassU)/eMassU;
   else if(species == NPlus5)           return (14.0064-5*eMassU)/eMassU;
   else if(species == NPlus6)           return (14.0064-6*eMassU)/eMassU;
   else if(species == NPlus7)           return (14.0064-7*eMassU)/eMassU;
   else if(species == NePlus0)          return (20.1767-0*eMassU)/eMassU;
   else if(species == NePlus1)          return (20.1767-1*eMassU)/eMassU;
   else if(species == NePlus2)          return (20.1767-2*eMassU)/eMassU;
   else if(species == NePlus3)          return (20.1767-3*eMassU)/eMassU;
   else if(species == NePlus4)          return (20.1767-4*eMassU)/eMassU;
   else if(species == NePlus5)          return (20.1767-5*eMassU)/eMassU;
   else if(species == NePlus6)          return (20.1767-6*eMassU)/eMassU;
   else if(species == NePlus7)          return (20.1767-7*eMassU)/eMassU;
   else if(species == NePlus8)          return (20.1767-8*eMassU)/eMassU;
   else if(species == NePlus9)          return (20.1767-9*eMassU)/eMassU;
   else if(species == NePlus10)         return (20.1767-10*eMassU)/eMassU;
   else if(species == AlPlus4)          return (26.9815-4*eMassU)/eMassU;
   else {  printf("Species' mass not defined\n");  exit(0);  }
}

int whatCharge(int species)
{
   int fail;

   if(species == Electron) 		return -1;
   else if(species == HPlus0)  		return 0;
   else if(species == HPlus1)  		return 1;
   else if(species == HePlus0)  	return 0;
   else if(species == HePlus1)  	return 1;
   else if(species == HePlus2)  	return 2;
   else if(species == CPlus0)           return 0;
   else if(species == CPlus1)           return 1;
   else if(species == CPlus2)           return 2;
   else if(species == CPlus3)           return 3;
   else if(species == CPlus4)           return 4;
   else if(species == CPlus5)           return 5;
   else if(species == CPlus6)           return 6;
   else if(species == NPlus0)           return 0;
   else if(species == NPlus1)           return 1;
   else if(species == NPlus2)           return 2;
   else if(species == NPlus3)           return 3;
   else if(species == NPlus4)           return 4;
   else if(species == NPlus5)           return 5;
   else if(species == NPlus6)           return 6;
   else if(species == NPlus7)           return 7;
   else if(species == NePlus0)           return 0;
   else if(species == NePlus1)           return 1;
   else if(species == NePlus2)           return 2;
   else if(species == NePlus3)           return 3;
   else if(species == NePlus4)           return 4;
   else if(species == NePlus5)           return 5;
   else if(species == NePlus6)           return 6;
   else if(species == NePlus7)           return 7;
   else if(species == NePlus8)           return 8;
   else if(species == NePlus9)           return 9;
   else if(species == NePlus10)          return 10;
   else if(species == AlPlus4)           return 4;
   else {  printf("Species' charge not defined\n");  exit(0);  }
}

