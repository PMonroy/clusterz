#ifndef VELOCITY
#define VELOCITY

#include "eqdate.hpp"
#include "vectorXYZ.hpp"



/* FUNCTIONS */
int SetupVflow(const string & vflowParamsFileName);
int LoadVLatLonGrid(eqdate rdate);
int LoadVelocities(eqdate startdate, int ntau);


void FreeMemoryVelocities(int ntau);

int GetVelocity(double t,vectorXYZ point, vectorXYZ *vint);
int GetVzPartialDeriv(double t,vectorXYZ point, double *dvzdz);

#endif
