#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

#include "vectorXYZ.hpp"

int MakeRegularGrid(vector<vectorXYZ> *point, vectorXYZ domainmin, vectorXYZ intergrid, vectorXYZ domainmax, int *dimx, int *dimy, int *dimz){

  double x,y,z;
  int i,j,k;

  int nx, ny, nz, numpoints;

  // INTEGER SPATIAL DIMENSIONS
  nz = floor((domainmax.z-domainmin.z)/intergrid.z);
  ny = floor((domainmax.y-domainmin.y)/intergrid.y);
  nx = floor((domainmax.x-domainmin.x)/intergrid.x);
  
  numpoints=nx*ny*nz;
  
  if(numpoints<1){
    cout <<"Error: Number of grid points is negative"<<endl;
    return 1;
  }

  (*point).reserve(numpoints);
  
  for(k=0; k<nz; k++){
    z=domainmin.z+intergrid.z*double(k);
    for(j=0; j<ny; j++){
      y=domainmin.y+intergrid.y*double(j);
      for(i=0; i<nx; i++){
	x=domainmin.x + intergrid.x*double(i);
	(*point).push_back(vectorXYZ(x,y,z));
      }
    }
  }
  
  *dimx=nx;
  *dimy=ny;
  *dimz=nz;  
  
  return 0;
}
