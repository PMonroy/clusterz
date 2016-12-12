#include <iostream>
#include <iomanip>
#include <netcdfcpp.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

#include "vflow.hpp"
#include "constants.hpp"
#include "eqdate.hpp"
#include "rparameters.hpp"
#include "VTK.hpp"

static const int NC_ERR = 2;

//GLOBAL VARIABLES
unsigned int DimLon;
unsigned int DimLonU;
unsigned int DimLat;
unsigned int DimLatV;
unsigned int DimDepth;
unsigned int DimTime;

vector <double> lon;
vector <double> lat;

vectorXYZ**** vflow;
double**** depth;

//char velocitydir[] = "/scratch/pmonroy/";
string vdir;

// Degree resolution
double degree_resolution;


/* Structures: Esta estructura tiene que ser interna a velocity.cpp */
struct  vectorIJK {
    int i;
    int j;
    int k[2][2];
};
//FUNCTIONS PROTOTYPES
vectorIJK GetIndices(vectorXYZ point);

//INLINE FUNCTIONS 
inline int LocateIndex(double x, double start, double step){
  return floor((x-start)/step);
}
inline int LocateIndex(double x, const vector <double> &xx){
  /* Given a vector xx[0...n-1], and given a value x, returns a value i such that x is between xx[i]
   * and xx[i+1]. xx must be monotonic, either increasing or decreasing. -1 or n is returned
   * to indicate that x is out of range.
   */

  if (x < xx.front()) 
    return -1;
  else if(x >= xx.back()) 
    return xx.size();
  else { 		
    unsigned long jl,jm,ju;
    jl=0;
    ju=xx.size()+1;
    int ascnd = (xx.back()>xx.front());
    while((ju-jl)>1) {
      jm =(ju+jl)>>1;
      if((x>=xx.at(jm-1))==ascnd) jl=jm;
      else ju=jm;
    }
    return int(jl-1);
  }
}

//FUNCTIONS
int SetupVflow(const string & vflowParamsFileName) {
  
  vdir=getStringParam(vflowParamsFileName, "vdir");
  degree_resolution=getDoubleParam(vflowParamsFileName, "degree_resolution");
  
  return 0;
}
int LoadVLatLonGrid(eqdate rdate) {

  NcError err(NcError::verbose_nonfatal);

  // Open Netcdf file correspond to date rdate
  char ncfile[256];
  sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",vdir.c_str(), rdate.GetYear(),rdate.GetMonth());
  
  NcFile dataFile(ncfile, NcFile::ReadOnly);
  // Check to see if the file was opened.
  if(!dataFile.is_valid()){
    cout << "Error opening file:"<< ncfile <<endl; 
   return NC_ERR;
  }

  NcDim *NcDimLon;
  NcDim *NcDimLonU;
  NcDim *NcDimLat;
  NcDim *NcDimLatV;
  NcDim *NcDimDepth;
  NcDim *NcDimTime;

  if(!(NcDimLon = dataFile.get_dim("xi_rho"))
     || !(NcDimLonU = dataFile.get_dim("xi_u"))
     || !(NcDimLat = dataFile.get_dim("eta_rho"))
     || !(NcDimLatV = dataFile.get_dim("eta_v"))
     || !(NcDimDepth = dataFile.get_dim("s_rho"))
     || !(NcDimTime = dataFile.get_dim("time"))){
    cout << "Error getting NcDimensions"<<endl;
    return NC_ERR;
  }

  DimLon = NcDimLon->size();
  DimLonU = NcDimLonU->size();
  DimLat = NcDimLat->size();
  DimLatV = NcDimLatV->size();
  DimDepth = NcDimDepth->size();
  DimTime = NcDimTime->size();

  lon.resize(DimLon);
  lat.resize(DimLat);

  // Get latitude 
  NcVar *NcVarLon;
  if(!(NcVarLon = dataFile.get_var("lon_rho"))
     || !NcVarLon->set_cur(0,0)
     || !NcVarLon->get(&lon[0], 1, DimLon)){
    cout << "Error reading longitude variable"<<endl;
    return NC_ERR;
  }

  // Get longitude 
  NcVar *NcVarLat;
  if(!(NcVarLat = dataFile.get_var("lat_rho"))
     || !NcVarLat->set_cur(0,0) 
     || !NcVarLat->get(&lat[0], DimLat, 1)){
    cout << "Error reading latitude variable"<<endl;
    return NC_ERR;
  }

  return 0;
}
int LoadVelocities(eqdate startdate, int ntau) {

  unsigned int t,i,j,k;
  int ft_center,fk_center,fj_center,q_center;
  int fk_vmon,fj_vmon,q_vmon;
  int fj_umon,q_umon;

  char ncfile[256];

  NcError err(NcError::verbose_nonfatal);

  // Buffer variables

  vector <double> Hbuffer(DimTime*DimLon*DimLat*DimDepth);
  vector <double> Ubuffer(DimTime*DimLonU*DimLat*DimDepth);
  vector <double> Vbuffer(DimTime*DimLon*DimLatV*DimDepth);
  vector <double> Wbuffer(DimTime*DimLon*DimLat*DimDepth);

  NcVar *NcVarDepth, *NcVarU, *NcVarV, *NcVarW;
  unsigned int time, final_time;
  eqdate date;
  unsigned int ndays, sumndays=0;

  time = startdate.GetDay();
  date.SetDay(time);
  final_time = startdate.GetDay() + ntau;


  /* Allocate memory */
  vflow = new vectorXYZ ***[ntau];
  depth = new double ***[ntau];

  for(t=0; t<(unsigned int) ntau; t++){
    vflow[t] = new vectorXYZ **[DimLon];
    depth[t] = new double **[DimLon];
    for(i=0; i<DimLon; i++){
      vflow[t][i] = new vectorXYZ *[DimLat];
      depth[t][i] = new double *[DimLat];
      for(j=0; j<DimLat; j++){
	vflow[t][i][j] = new vectorXYZ [DimDepth];
	depth[t][i][j] = new double [DimDepth];
      }
    }
  }



  sumndays=0;
  while(time < final_time)
    {

      ndays = DimTime - (date.GetMday()-1);
      if((time + ndays) > final_time)
	ndays = final_time - time;
      
      sprintf(ncfile, "%sextract_roms_avg_Y%1uM%1u.nc.1",vdir.c_str(), date.GetYear(), date.GetMonth());
      cout << "reading netcdf file: " <<ncfile <<" start day="<< date.GetMday() <<" number of days=" << ndays<<endl;

      NcFile dataFile(ncfile, NcFile::ReadOnly);

      // Check to see if the file was opened.
      if(!dataFile.is_valid())
	return NC_ERR;
 
      //Read depth
      if (!(NcVarDepth = dataFile.get_var("depth"))
	  || !NcVarDepth->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarDepth->get(&Hbuffer[0], ndays, DimDepth, DimLat, DimLon))
	return NC_ERR;

      //Read u
      if (!(NcVarU = dataFile.get_var("u"))
	  || !NcVarU->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarU->get(&Ubuffer[0], ndays, DimDepth, DimLat, DimLonU))
	return NC_ERR;

      //Read v
      if (!(NcVarV = dataFile.get_var("v"))
	  || !NcVarV->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarV->get(&Vbuffer[0], ndays, DimDepth, DimLatV, DimLon))
	return NC_ERR;

      //Read w
      if (!(NcVarW = dataFile.get_var("w"))
	  || !NcVarW->set_cur(date.GetMday()-1, 0, 0, 0)
	  || !NcVarW->get(&Wbuffer[0], ndays, DimDepth, DimLat, DimLon))
	return NC_ERR;

      dataFile.close();


      /* copy values and averages it */

      for(t=0; t<ndays; t++) {
	ft_center=t*DimDepth;
	
	for(k=0; k<DimDepth; k++){
	  fk_center=DimLat*(ft_center+k);
	  fk_vmon=DimLatV*(ft_center+k);
	    
	  for(j=0; j<DimLat; j++){
	    fj_center=DimLon*(fk_center+j);
	    fj_vmon=DimLon*(fk_vmon+j);
	    fj_umon=DimLonU*(fk_center+j);
	    for(i=0; i<DimLon; i++){
	      q_center=fj_center+i;
	      q_umon=fj_umon+i;
	      q_vmon=fj_vmon+i;
	      depth[t+sumndays][i][j][k]=Hbuffer[q_center];
	      vflow[t+sumndays][i][j][k].z=Wbuffer[q_center];
	      vflow[t+sumndays][i][j][k].x=(Ubuffer[q_umon-(i==(DimLon-1))]+Ubuffer[q_umon-(i!=0)])/2.0;
	      vflow[t+sumndays][i][j][k].y=(Vbuffer[q_vmon-(j==(DimLat-1))*DimLon]+Vbuffer[q_vmon-(j!=0)*DimLon])/2.0;

	      vflow[t+sumndays][i][j][k]=secondsday*vflow[t+sumndays][i][j][k];
	      
	    }
	  }
	}
      }

      time += ndays;
      date.SetDay(time);
      sumndays += ndays;
    }


  MakeVTKStructuredGrid_Vfield(&lon[0], &lat[0], depth[0], vflow[0], DimLon, DimLat, DimDepth, "velocity.vtk");

  return 0;
}

void FreeMemoryVelocities(int ntau) {
  int t,i,j;
  
  for(t=0; t<ntau; t++){
    for(i = 0; i < int(DimLon); i++){
      for(j = 0; j < int(DimLat); j++){
	delete[] depth[t][i][j];
	delete[] vflow[t][i][j];
      }
      delete[] depth[t][i];
      delete[] vflow[t][i]; 
    }
    delete[] depth[t];
    delete[] vflow[t]; 
  }
  delete[] depth;
  delete[] vflow;   
}


int GetIndices(unsigned long time, vectorXYZ point, vectorIJK *index) {

  /* point.x -> longitude
   * point.y -> latitude
   * point.z -> depth
   */

  int i,j,ki,kj;

  /* Locate index longitude*/
  index->i = LocateIndex(point.x, lon[0], degree_resolution);

  if(index->i < 0 || index->i >= int(DimLon-1))
      return 1;

  /* Locate index latitude*/
  index->j = LocateIndex(point.y, lat);

  if(index->j < 0 || index->j >= int(DimLat-1))
      return 1;
  
  /* Locate index depth */
  for(i=index->i,ki=0; i<index->i+2; i++,ki++)
    {
      for(j=index->j,kj=0; j<index->j+2; j++,kj++)
	{
	  //vector<double> dpt (depth[time][i][j], depth[time][i][j] + sizeof(depth[time][i][j])/sizeof(depth[time][i][j][0]));
	  vector<double> dpt;// (depth[time][i][j], depth[time][i][j] + DimDepth);
	  dpt.reserve (DimDepth); // Reserve memory not to allocate it 10 times...
	  for (unsigned int k=0; k<DimDepth; k++){
            dpt.push_back(depth[time][i][j][k]);
	  }
	  index->k[ki][kj] = LocateIndex(point.z, dpt);
	  if(index->k[ki][kj] < 0 || index->k[ki][kj] >= int(DimDepth-1))
	    return 1;
	  
	}
    }

  return 0;
}
int GetVelocity(double t,vectorXYZ point, vectorXYZ *vint) {
  vectorXYZ vgrid[16];
  vectorXYZ vcomp[16];
  vectorIJK index[2];
  unsigned long time;

  double alpha, beta;

  int h,i,j,k;
  int deltatime,deltai,deltaj,deltak;
  unsigned int q;

  /* Calculates the vectors vb[15] and points ptmb[15] */
  time = (unsigned long) t;

  if(GetIndices(time, point, &index[0])==1)
    return 1;

  if(GetIndices(time+1, point, &index[1])==1)  //increment 1 day,depend on vflow 
    return 1;

  /* Vectors and points with only one int index*/
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++)
    {
      for(deltai=0; deltai<2; deltai++)
	{
	  for(deltaj = 0; deltaj<2; deltaj++)
	    {
	      for(deltak=0; deltak<2; deltak++)
		{
		  i = index[deltatime].i + deltai;
		  j = index[deltatime].j + deltaj;
		  k = index[deltatime].k[deltai][deltaj] + deltak;
		  h = time+deltatime;
		  
		  vgrid[q].x = rads*lon[i];
		  vgrid[q].y = log(fabs((1.0/cos(rads*lat[j]))+tan(rads*lat[j])));
		  vgrid[q].z = depth[h][i][j][k];
		  
		  vcomp[q] = vflow[h][i][j][k];
	
		  q++; 
		}
	    }
	}
    }

  point.x = rads*point.x;
  point.y = log(fabs((1.0/cos(rads*point.y))+tan(rads*point.y)));

  /* COAST CHECKING: I think is only need to check land mask 
  if( (point.z < bathymetry[index[0].i  ][index[0].j  ]) &&
      (point.z < bathymetry[index[0].i+1][index[0].j  ]) &&
      (point.z < bathymetry[index[0].i  ][index[0].j+1]) &&
      (point.z < bathymetry[index[0].i+1][index[0].j+1]))
      return 1;// The particle has reached the coast */

  /* Depth Interpolation: */
  for(q=0; q<16; q+=2)
    {
      alpha = (vgrid[q+1].z - point.z)/(vgrid[q+1].z - vgrid[q].z);
      beta = 1 - alpha;
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q]; 
    }

  /* Lat Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(q = 0; q < 8; q+=2)
    {
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Phi Interpolation */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1.0 - alpha;
  for(q = 0; q < 4; q+=2)
    {
      vcomp[q>>1] = alpha * vcomp[q] + beta * vcomp[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Time Interpolation: */ 
  alpha = ((double) (time + 1)) - t;  
  beta = 1.0 - alpha;   
  vcomp[0] = alpha * vcomp[0] + beta * vcomp[1];
  /* Interpolated V*/
  *vint = vcomp[0];

  return 0;
}

int GetVzPartialDeriv(double t,vectorXYZ point, double *dvzdz){

  unsigned int time;
  
  vectorXYZ delta;
  vectorIJK index[2];

  vectorXYZ vgrid[16];
  double partialderiv[16];
  double alpha,beta;

  unsigned int l,i,j,k;
  unsigned int k0;
  unsigned int k1;
  unsigned int deltai,deltaj,deltak,deltatime;
  unsigned int q;

  double h1,h2;

  /* Calculates the vectors vb[15] and points ptmb[15] */
  time = (unsigned long) t;

  if(GetIndices(time, point, &index[0])==1)
    return 1;

  if(GetIndices(time+1, point, &index[1])==1)  //increment 1 day,depend on vflow 
    return 1;

  /* Vectors and points with only one int index*/
  q=0;  
  for(deltatime=0; deltatime<2; deltatime++)
    {
      for(deltai=0; deltai<2; deltai++)
	{
	  for(deltaj = 0; deltaj<2; deltaj++)
	    {
	      for(deltak=0; deltak<2; deltak++)
		{
		  i = index[deltatime].i + deltai;
		  j = index[deltatime].j + deltaj;
		  k = index[deltatime].k[deltai][deltaj] + deltak;
		  l = time+deltatime;
		  
		  k0 = k-1;
		  k1 = k+1;

		  if(k0<0 || k0> int(DimDepth-1))
		    return 1;

		  if(k1<0 || k1>int(DimDepth-1))
		    return 1;

		  vgrid[q].x = rads*lon[i];
		  vgrid[q].y = log(fabs((1.0/cos(rads*lat[j]))+tan(rads*lat[j])));
		  vgrid[q].z = depth[l][i][j][k];
		  
		  h1=depth[l][i][j][k]-depth[l][i][j][k0];
		  h2=depth[l][i][j][k1]-depth[l][i][j][k];
		  
		  partialderiv[q] = (vflow[l][i][j][k1].z*(h1/h2)-vflow[l][i][j][k0].z*(h2/h1))/(h1+h2);
		  partialderiv[q] += vflow[l][i][j][k].z*((1.0/h1)-(1.0/h2));
	
		  q++; 
		}
	    }
	}
    }

  point.x = rads*point.x;
  point.y = log(fabs((1.0/cos(rads*point.y))+tan(rads*point.y)));

  /* Depth Interpolation: */
  for(q=0; q<16; q+=2)
    {
      alpha = (vgrid[q+1].z - point.z)/(vgrid[q+1].z - vgrid[q].z);
      beta = 1 - alpha;
      partialderiv[q>>1] = alpha * partialderiv[q] + beta * partialderiv[q+1];
      vgrid[q>>1] = vgrid[q]; 
    }

  /* Lat Interpolation: */
  alpha = (vgrid[1].y - point.y)/(vgrid[1].y - vgrid[0].y);
  beta = 1.0 - alpha;
  for(q = 0; q < 8; q+=2)
    {
      partialderiv[q>>1] = alpha * partialderiv[q] + beta * partialderiv[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Phi Interpolation */
  alpha = (vgrid[1].x - point.x)/(vgrid[1].x - vgrid[0].x);
  beta = 1.0 - alpha;
  for(q = 0; q < 4; q+=2)
    {
      partialderiv[q>>1] = alpha * partialderiv[q] + beta * partialderiv[q+1];
      vgrid[q>>1] = vgrid[q];
    }

  /* Time Interpolation: */

  alpha = ((double) (time + 1)) - t;
  beta = 1 - alpha;

  partialderiv[0] = alpha*partialderiv[0]+beta*partialderiv[1];

  *dvzdz = partialderiv[0];
  return 0;
}
