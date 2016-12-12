#ifndef GRIDCONSTRUCTION
#define GRIDCONSTRUCTION

#include "vectorXYZ.hpp"
 
int MakeRegularGrid(vector <vectorXYZ> *point, vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax,int *dimx, int *dimy, int *dimz);

#endif
