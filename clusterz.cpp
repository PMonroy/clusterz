#include <cstdio>
#include <cstdlib>
#include <iostream> // cout, endl...  
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <iomanip> // setfill, setw
#include <string>
#include <omp.h>

using namespace std;

#include "rparameters.hpp"
#include "gridconstruction.hpp"
#include "lagrangian_engine.hpp"
#include "constants.hpp"
#include "vflow.hpp"
#include "eqdate.hpp"
#include "vectorXYZ.hpp"

string numprintf(int ndigits, int ndecimals, double number);
struct clusterzParameters {

  const vectorXYZ domainmin;
  const vectorXYZ domainmax;
  const vectorXYZ intergrid;
  const double density;
  const eqdate seeddate;
  const double intstep;
  const double vsink;
  const double maxdepth;
  
  // Define a constructor that will load stuff from a configuration file.
  clusterzParameters(const string & clusterzParamsFileName)
  :domainmin(getVectorXYZParam(clusterzParamsFileName, "domainmin"))
  ,domainmax(getVectorXYZParam(clusterzParamsFileName, "domainmax"))
  ,intergrid(getVectorXYZParam(clusterzParamsFileName, "intergrid"))
  ,density(getDoubleParam(clusterzParamsFileName, "density"))  
  ,seeddate(getEqDateParam(clusterzParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(clusterzParamsFileName, "intstep")) 
  ,vsink(getDoubleParam(clusterzParamsFileName, "vsink"))
  ,maxdepth(getDoubleParam(clusterzParamsFileName, "maxdepth"))
{}
};

int main(int argc, char **argv){

  /*********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string fnameparams;// File name that stores the parameters
  int namefileflag;
  int vequation;
  string me=argv[0];
  
  if(GetcmdlineParameters(argc, argv, &fnameparams, &namefileflag, &vequation)) {//Get cmd line parameters
    cout << me << ": Error getting parameters from command line" <<endl;
    return 1;
  }
#ifdef DEBUG
  cout << "Clusterz Parameters from command line: "<<endl;
  cout <<"ini file ="<<fnameparams<<": "<<endl; 
  cout << "V. equation = "<< vequation <<endl;
#endif
  
  const clusterzParameters clusterzParams(fnameparams);

#ifdef DEBUG
  cout << "Clusterz Parameters from file ";
  cout <<"\""<<fnameparams<<": "<<endl; 
  cout << " domainmin = "<< clusterzParams.domainmin<<endl;
  cout << " intergrid = " <<clusterzParams.intergrid<<endl;
  cout << " domainmax = "<< clusterzParams.domainmax<<endl;
  cout << " density = "<< clusterzParams.density<<endl;
  cout << " seeddate = "<< clusterzParams.seeddate.GetMday()<<"-";
  cout<<clusterzParams.seeddate.GetMonth()<<"-";
  cout<<clusterzParams.seeddate.GetYear()<<endl;
  cout << " intstep = "<<clusterzParams.intstep<<endl;
  cout << " vsink = "<<clusterzParams.vsink<<endl;
  cout << " maxdepth = "<<clusterzParams.maxdepth<<endl ;  
  cout << " [Complete]" <<endl;
  cout << endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> grid;
  int GridDimx, GridDimy,GridDimz;
 
  if(MakeRegularGrid(&grid, clusterzParams.domainmin, 
		     clusterzParams.intergrid, 
		     clusterzParams.domainmax,
		     &GridDimx,
		     &GridDimy,
		     &GridDimz)){//Grid construction 
    cout << "Regular Grid: [Fail]" << endl;
    return 1;
  }
  unsigned int numgridpoints=grid.size();

#ifdef DEBUG
  cout << "Grid construction: "<<endl; 
  cout << " num. nodes = "<< grid.size() <<endl;
  cout << "Dim x = "<< GridDimx << endl;
  cout << "Dim y = "<< GridDimy << endl;
  cout << "Dim z = "<< GridDimz << endl;
  cout << "[Complete]" << endl;
  cout << endl;
#endif


  /**********************************************
   * SETUP TIME PARAMETERS
   **********************************************/
  int tau=2*int((clusterzParams.domainmax.z-clusterzParams.maxdepth)/clusterzParams.vsink);
  
  eqdate inidate;

  double tstart;
  double tend;
  double h;

  int ntime = abs(tau);
  int ascnd = tau > 0;
  unsigned int seedtime = clusterzParams.seeddate.GetDay();

  if(ascnd){
    unsigned int initime = seedtime;
    inidate.SetDay(initime);
    tstart = 0.0;
    tend = (double) ntime;
    h=clusterzParams.intstep;
  }
  else{
    unsigned int initime = seedtime-ntime;
    inidate.SetDay(initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*clusterzParams.intstep;
  }

#ifdef DEBUG
  cout << "Setup time parameters: "<<endl;
  cout << " inidate = "<< inidate.GetMday() <<"-";
  cout<< inidate.GetMonth() <<"-";
  cout<< inidate.GetYear() <<endl;
  cout << " tstart = "<< tstart <<endl;
  cout << " tend = " << tend <<endl;
  cout << " tau = " << tau <<endl;
  cout << " h = "<< h <<endl;
#endif

  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/

#ifdef DEBUG
  cout<<"Loading velocity grid:"<<endl; //Verbose: loading vel. grid
#endif

  SetupVflow(fnameparams);

  if((LoadVLatLonGrid(clusterzParams.seeddate))!=0){//Load velocity grid
      cout<<"[Fail]"<<endl;
      return 1;
  }

#ifdef DEBUG
  cout<<"[Complete]" << endl;//Success loading vel. grid
  cout<<"Loading velocities from model:" << endl;//Verbose: loading velocities from model 
#endif

  if(LoadVelocities(inidate,ntime+2)!=0){// Load velocity field
    cout << "[Fail]"<< endl;
    return 1;
  }

#ifdef DEBUG//Succes loading velocities
  cout << "[Complete]"<<endl;
#endif


  /**********************************************
   * LANGRANGIAN ENGINE
   **********************************************/


  SetupLagrangianEngine(fnameparams);

  // Tracer grid
  vector<vectorXYZ> tracer,itracer;
  vectorXYZ rshift;
  vectorXYZ tracerintergrid;
  
  tracerintergrid.x=clusterzParams.intergrid.x/sqrt(clusterzParams.density);
  tracerintergrid.y=clusterzParams.intergrid.y/sqrt(clusterzParams.density);
  tracerintergrid.z=clusterzParams.intergrid.z;

  vectorXYZ tracerdomainmin;
  tracerdomainmin.x=clusterzParams.domainmin.x;
  tracerdomainmin.y=clusterzParams.domainmin.y;
  tracerdomainmin.z=clusterzParams.domainmin.z;

  vectorXYZ tracerdomainmax;
  tracerdomainmax.x=clusterzParams.domainmax.x+clusterzParams.intergrid.x;
  tracerdomainmax.y=clusterzParams.domainmax.y+clusterzParams.intergrid.y;
  tracerdomainmax.z=clusterzParams.domainmax.z;

  int TracerDimx, TracerDimy, TracerDimz;
 
  
  if(MakeRegularGrid(&itracer, tracerdomainmin, 
		     tracerintergrid, 
		     tracerdomainmax,
		     &TracerDimx,
		     &TracerDimy,
		     &TracerDimz)){//Tracer initial grid construction 
    cout << "Regular Tracer grid: [Fail]" << endl;
    return 1;
  }
  
  rshift.x=tracerintergrid.x/2.0;
  rshift.y=tracerintergrid.y/2.0;
  rshift.z=tracerintergrid.z/2.0;

  unsigned int numtracers=itracer.size();
  unsigned int q;

  for(q=0;q<numtracers; q++)
    itracer[q]=itracer[q]+rshift;
  
  tracer=itracer;


  
  vector<int> qflaggrid(numgridpoints,0);
  vector<double> dilationgrid(numgridpoints,0);
  vector<double> atimegrid(numgridpoints,abs(tau));  
  vector<int> nbox(numgridpoints,0);
  vector<double> efactor(numgridpoints,0.0);    

    
  vector<int> qflagtracer(numtracers,0);
  vector<double> dilationtracer(numtracers,0);
  vector<double> atimetracer(numtracers,abs(tau));  

  
  unsigned int step;
  double t;

  double dvbuffer;
  
  int depthdescnd=(clusterzParams.domainmax.z-clusterzParams.maxdepth)>0;

  cout << "depthdescnd=" << depthdescnd<<endl;
  omp_set_num_threads(7);
#pragma omp parallel for default(shared) private(t) // Parallelizing the code for computing trajectories

  for (q=0; q<numtracers; q++) {

    vector<double> dvzdz;

    dvzdz.reserve(abs(tau/clusterzParams.intstep));

    for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++) {
      //new point integration
      if(vequation==0){
	qflagtracer[q]=RK4(t, h, &tracer[q], GetVflowplusVsink);
      }
      else if(vequation==1){
	qflagtracer[q]=Heun2(t, h, &tracer[q], GetVflowplusVsink);
      }
	
      if(qflagtracer[q]==1){
	atimetracer[q]=abs(t-tstart);
	break;
      }
      //calc deriv
      qflagtracer[q]=GetVzPartialDeriv(t, tracer[q], &dvbuffer);
      if(qflagtracer[q]==1){
	atimetracer[q]=abs(t-tstart);
	break;
      }
      dvzdz.push_back((2*ascnd-1)*dvbuffer);

      // Check maxdepth
      if(depthdescnd==1?(tracer[q].z<=clusterzParams.maxdepth):(tracer[q].z>=clusterzParams.maxdepth)){
	atimetracer[q]=abs(t-tstart);
	break;
      }
    }


    if(qflagtracer[q]==0){

      efactor[q]=cos(rads*itracer[q].y)/cos(rads*tracer[q].y);

      for(unsigned int p=1; p<(dvzdz.size()-1); p++)
	dilationtracer[q]+=dvzdz[p];

      dilationtracer[q]=dvzdz.front() + 2*dilationtracer[q] + dvzdz.back();
      dilationtracer[q]=(clusterzParams.intstep/2.0)*dilationtracer[q];
    }  
  }
  
 // find index grid of initial position and final and add one to nini and nend in the corresponding grid box

  vectorXYZ rini;
  
  int iini,jini;
  int qini;

  for (q=0; q<numtracers; q++){
    rini=(itracer[q]-clusterzParams.domainmin)/clusterzParams.intergrid;

    // locate initial cell
    iini=floor(rini.x);
    jini=floor(rini.y);

    if(iini>=GridDimx || iini<0) continue;
    if(jini>=GridDimy || jini<0) continue;

    qini=iini+(jini*GridDimx);

    if(qflagtracer[q]==0){ 
      nbox[qini]++;
      dilationgrid[qini]+=dilationtracer[q];
      atimegrid[qini]+=atimetracer[q];
    }
    else{
      qflaggrid[qini]+=qflagtracer[q];
    }
  }

  for(q=0; q<numgridpoints; q++) {
    if(nbox[q]>0){
	dilationgrid[q]=dilationgrid[q]/nbox[q];
	atimegrid[q]=atimegrid[q]/nbox[q];
      }
  }
  
  /**********************************************
   * FREE RESOURCES
   **********************************************/

#ifdef DEBUG
  cout<<"freeing velocity :";
#endif

  FreeMemoryVelocities(ntime+2);

#ifdef DEBUG
  cout<<"[OUYEAH]"<<endl;
#endif

 
  /********************************************************************************
   * WRITING RESULT IN VTK FILE
   ********************************************************************************/
 // save clusterz values in a file

  string rawname;

  if(namefileflag==1){
    size_t lastdot = fnameparams.find_last_of(".");
    if(lastdot == string::npos){
      rawname=fnameparams;
    }
    else{
      rawname=fnameparams.substr(0,lastdot);
    }
  }
  else{
    rawname = "_lon"    + numprintf(4,1,clusterzParams.domainmin.x) 
      + numprintf(4,1,clusterzParams.domainmax.x)
      + numprintf(4,3,clusterzParams.intergrid.x)
      + "_lat"       + numprintf(4,1,clusterzParams.domainmin.y)
      + numprintf(4,1,clusterzParams.domainmax.y)
      + numprintf(4,3,clusterzParams.intergrid.y)
      + "_dpt"       + numprintf(4,0,clusterzParams.domainmin.z)
      + numprintf(4,0,clusterzParams.domainmax.z)
      + numprintf(4,0,clusterzParams.intergrid.z)
      + "_ts"        + numprintf(3,2,clusterzParams.intstep) 
      + "_mdpt"         + numprintf(4,0,clusterzParams.maxdepth)
      + "_den" + numprintf(4,0,clusterzParams.density)
      + "_v"         + numprintf(3,0,clusterzParams.vsink)
      + "_veq"         + numprintf(1,0,vequation)
      + "_d"         + numprintf(2,0,clusterzParams.seeddate.GetMday())
      + numprintf(2,0,clusterzParams.seeddate.GetMonth())
      + numprintf(2,0,clusterzParams.seeddate.GetYear());
  }

  string namevtkfilegrid = "clusterz"+ rawname + ".vtk";
  string namefiletracer = "clus_tracer"+ rawname + ".dat";
  string namefilecluster = "clusterz"+ rawname + ".dat";

#ifdef DEBUG
  cout << "Save cluster field in vtk file: " << namevtkfilegrid <<endl;
#endif

  ofstream vtkfile(namevtkfilegrid.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"grid cluster 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< GridDimx+1 <<" "<< GridDimy+1 <<" "<< GridDimz+1 <<endl;
  vtkfile<<"ORIGIN "<<grid[0].x<<" "<<grid[0].y<<" "<< grid[0].z <<endl;
  vtkfile<<"SPACING "<<clusterzParams.intergrid.x<<" "<<clusterzParams.intergrid.y<<" "<< clusterzParams.intergrid.z <<endl;
  vtkfile<<"CELL_DATA "<<(GridDimx)*(GridDimy)*(GridDimz)<<endl;
  vtkfile<<"SCALARS cluster float 1"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<dilationgrid.size(); q++) {
      vtkfile<<dilationgrid[q]<<endl;
  }
  vtkfile<<"SCALARS efactor float 1"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<efactor.size(); q++) {
      vtkfile<<efactor[q]<<endl;
  }
  vtkfile<<"SCALARS qflag int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<qflaggrid.size(); q++) {
    vtkfile<<qflaggrid[q]<<endl;
  }
  vtkfile.close();

#ifdef DEBUG
  cout << "Save data tracer in dat file: " << namefiletracer <<endl;
#endif

  ofstream filetracer(namefiletracer.c_str());
  for(q=0; q<numgridpoints; q++) {
    filetracer<<itracer[q]<<" ";
    filetracer<<tracer[q]<<" ";
    filetracer<<atimetracer[q]<<" ";
    filetracer<<qflagtracer[q]<<endl;
  }
  filetracer.close();

#ifdef DEBUG
  cout << "Save data cluster in dat file: " << namefilecluster <<endl;
#endif

  ofstream filecluster(namefilecluster.c_str());
  for(q=0; q<numgridpoints; q++) {
    filecluster<<grid[q]<<" ";
    filecluster<<qflaggrid[q]<<" ";
    filecluster<<dilationgrid[q]<<endl;
  }
  filecluster.close();
  return 0;     
}




string numprintf(int ndigits, int ndecimals, double number){
  ostringstream stream;// it needs to include <sstream>
  double factor=pow(10,ndecimals); // it needs to include <cmath>
  stream << fixed; //it needs to include <iostream>
  stream << setfill('0') << setw(ndigits);
  stream << setprecision(0) << internal<< (number*factor);
  return stream.str();
}
