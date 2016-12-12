#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>

using namespace std;

#include "vflow.hpp" // Functions to read velocities 
#include "constants.hpp"
#include "rparameters.hpp"
#include "mt19937ar.hpp"

double vsink;
double Dh;
double Dz;

//Functions
double gasdev(double (*nrand)(void ));

int SetupLagrangianEngine(const string & LagEngParamsFileName) {
  vsink=getDoubleParam(LagEngParamsFileName, "vsink");
  Dh=getDoubleParam(LagEngParamsFileName, "Dh");
  Dh=Dh*secondsday;
  Dz=getDoubleParam(LagEngParamsFileName, "Dz");
  Dz=Dz*secondsday;

#ifdef DEBUG
  cout << "Lagragian  parameters";
  cout << " vsink = "<< vsink <<endl;
  cout << " Dh = " << Dh <<endl;
  cout << " Dz = " << Dz <<endl;
  cout << " [Complete]" <<endl;
  cout << endl;
#endif
  
  return 0;
}

int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* )) {
  vectorXYZ point2,point3,point4;
  vectorXYZ v1,v2,v3,v4;
  
  double tstep2,tstep6;
  double h;//scale factor for spherical coordinates
  double t;

  /* Time increments */
  tstep2 = intstep*0.5;
  tstep6 = intstep/6.0;
  
  /* Calculate V1: */
  if(velocity(t0,*point, &v1))
    return 1;
  h = rearth * cos(rads*(point->y));
  v1.x = degrees*(v1.x / h );
  v1.y = degrees*(v1.y / rearth);


  /* Calculate V2: */
  t = t0 + tstep2;
  point2 = *point + (tstep2 * v1);

  if(velocity(t,point2, &v2))
    return 1;

  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity


  /* Calculate V3: */
  point3 = *point + (tstep2 * v2);

  if(velocity(t,point3, &v3))
    return 1;

  h = rearth * cos(rads*(point3.y));
  v3.x = degrees*(v3.x / h);
  v3.y = degrees*(v3.y / rearth);

  
  /* Calculate V4: */
  t = t0 + intstep;
  point4 = *point + (intstep * v3);
  
  if(velocity(t,point4, &v4))
    return 1;

  h = rearth * cos(rads*(point4.y));
  v4.x = degrees*(v4.x / h);
  v4.y = degrees*(v4.y / rearth);

  /* Calculate Final point */  
  *point += (tstep6 * (v1 + v4 + 2.0*v2 + 2.0*v3));

  return 0;
}
int Heun2(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* )) {
  vectorXYZ point1,point2, stddev;
  vectorXYZ v1,v2;
  vectorXYZ g1,g2;
  vectorXYZ u, uh, aux;
  vectorXYZ k,l;
  double h;
  static int rset=0;
  static unsigned long seed=329571112321;
  
  stddev.x = sqrt(2.0*Dh);
  stddev.y = stddev.x;
  stddev.z = sqrt(2.0*Dz);

  if(rset == 0){
      rset = 1;
      init_genrand(seed);
  }

  u.x = gasdev(genrand_real2); 
  u.y = gasdev(genrand_real2); 
  u.z = gasdev(genrand_real2); 
  
  /* Calculate V1: */
  point1 = *point;
  if(velocity(t0, point1, &v1))
    return 1;
  h = rearth * cos(rads*(point1.y));
  v1.x = degrees*(v1.x / h ); 
  v1.y = degrees*(v1.y / rearth);

  /* Calculate G1*/
  g1 = stddev;
  g1.x = degrees*(g1.x / h );
  g1.y = degrees*(g1.y / rearth);
 

  /* Calculate V2: */
  uh=sqrt(abs(intstep))*u;
  aux = intstep * v1 + uh * g1;
  point2 = point1 + aux;

  if(velocity(t0+intstep,point2, &v2))
    return 1;
  h = rearth * cos(rads*(point2.y));
  v2.x = degrees*(v2.x / h ); // rads velocity
  v2.y = degrees*(v2.y / rearth); // rads velocity

  /* Calculate G2*/
  g2 = stddev;
  g2.x = degrees*(g2.x / h ); // rads velocity
  g2.y = degrees*(g2.y / rearth); // rads velocity

  /* Compute the new position */
  *point += 0.5*(aux + intstep*v2+ uh*g2);

  return 0;
}

double gasdev(double (*nrand)(void ))
{

  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;
  
  if (iset == 0) 
    {
      do 
	{
	  v1=2.0*nrand()-1.0;
	  v2=2.0*nrand()-1.0;
	  rsq=v1*v1+v2*v2;
	} 
      while (rsq >= 1.0 || rsq == 0.0);

      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      return v2*fac;
    } 
  else 
    {
      iset=0;
      return gset;
    }
}
int GetVflowplusVsink(double t,vectorXYZ point, vectorXYZ *vint) {
  if(GetVelocity( t, point, vint))
    return 1;
  
  vint->z = vint->z - vsink;
  return 0;
}

/* TEST RANDOM NUMBER GENERATOR

int mttest(void)
{
  int i;
  unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;
  init_by_array(init, length);
  printf("1000 outputs of genrand_int32()\n");
  for (i=0; i<1000; i++) {
    printf("%10lu ", genrand_int32());
    if (i%5==4) printf("\n");
  }
  printf("\n1000 outputs of genrand_real2()\n");
  for (i=0; i<1000; i++) {
    printf("%10.8f ", genrand_real2());
    if (i%5==4) printf("\n");
  }

  return 0;
}

*/
