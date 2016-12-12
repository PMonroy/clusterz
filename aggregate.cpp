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
struct clusterParameters {

  const vectorXYZ domainmin;
  const vectorXYZ domainmax;
  const vectorXYZ intergrid;
  const double density;
  const eqdate seeddate;
  const double intstep;
  const double vsink;
  const double maxdepth;
  
  // Define a constructor that will load stuff from a configuration file.
  clusterParameters(const string & clusterParamsFileName)
  :domainmin(getVectorXYZParam(clusterParamsFileName, "domainmin"))
  ,domainmax(getVectorXYZParam(clusterParamsFileName, "domainmax"))
  ,intergrid(getVectorXYZParam(clusterParamsFileName, "intergrid"))
  ,density(getDoubleParam(clusterParamsFileName, "density"))  
  ,seeddate(getEqDateParam(clusterParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(clusterParamsFileName, "intstep")) 
  ,vsink(getDoubleParam(clusterParamsFileName, "vsink"))
  ,maxdepth(getDoubleParam(clusterParamsFileName, "maxdepth"))
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
    cout << "Error getting parameters from file" <<endl;
    return 1;
  }
#ifdef DEBUG
  cout << "Clusterz Parameters from command line: "<<endl;
  cout <<"ini file ="<<fnameparams<<": "<<endl; 
  cout << "V. equation = "<< vequation <<endl;
#endif
  
  const clusterParameters clusterParams(fnameparams);

#ifdef DEBUG
  cout << "Cluster Parameters from file ";
  cout <<"\""<<fnameparams<<": "<<endl; 
  cout << " domainmin = "<< clusterParams.domainmin<<endl;
  cout << " intergrid = " <<clusterParams.intergrid<<endl;
  cout << " domainmax = "<< clusterParams.domainmax<<endl;
  cout << " density = "<< clusterParams.density<<endl;
  cout << " seeddate = "<< clusterParams.seeddate.GetMday()<<"-";
  cout<<clusterParams.seeddate.GetMonth()<<"-";
  cout<<clusterParams.seeddate.GetYear()<<endl;
  cout << " intstep = "<<clusterParams.intstep<<endl;
  cout << " vsink = "<<clusterParams.vsink<<endl;
  cout << " maxdepth = "<<clusterParams.maxdepth<<endl ;  
  cout << " [Complete]" <<endl;
  cout << endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> grid;
  int GridDimx, GridDimy,GridDimz;
 
  if(MakeRegularGrid(&grid, clusterParams.domainmin, 
		     clusterParams.intergrid, 
		     clusterParams.domainmax,
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
  int tau=2*int((clusterParams.domainmax.z-clusterParams.maxdepth)/clusterParams.vsink);
  
  eqdate inidate;

  double tstart;
  double tend;
  double h;

  int ntime = abs(tau);
  int ascnd = tau > 0;
  unsigned int seedtime = clusterParams.seeddate.GetDay();

  if(ascnd){
    unsigned int initime = seedtime;
    inidate.SetDay(initime);
    tstart = 0.0;
    tend = (double) ntime;
    h=clusterParams.intstep;
  }
  else{
    unsigned int initime = seedtime-ntime;
    inidate.SetDay(initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*clusterParams.intstep;
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

  if((LoadVLatLonGrid(clusterParams.seeddate))!=0){//Load velocity grid
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

  //Tracer Variables
  vector<vectorXYZ> itracer;
  vector<vectorXYZ> ftracer;


  vectorXYZ tracerintergrid;
  tracerintergrid.x=clusterParams.intergrid.x/sqrt(clusterParams.density);
  tracerintergrid.y=clusterParams.intergrid.y/sqrt(clusterParams.density);
  tracerintergrid.z=clusterParams.intergrid.z;

  vectorXYZ tracerdomainmin;
  tracerdomainmin.x=clusterParams.domainmin.x;
  tracerdomainmin.y=clusterParams.domainmin.y;
  tracerdomainmin.z=clusterParams.domainmin.z;

  vectorXYZ tracerdomainmax;
  tracerdomainmax.x=clusterParams.domainmax.x+clusterParams.intergrid.x;
  tracerdomainmax.y=clusterParams.domainmax.y+clusterParams.intergrid.y;
  tracerdomainmax.z=clusterParams.domainmax.z;

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
  vectorXYZ rshift;
  rshift.x=tracerintergrid.x/2.0;
  rshift.y=tracerintergrid.y/2.0;
  rshift.z=tracerintergrid.z/2.0;
  unsigned int numtracer=itracer.size();
  unsigned int q;
  for(q=0;q<numtracer; q++)
    itracer[q]=itracer[q]+rshift;
  
  ftracer=itracer;
  
  cout << "number of tracers =" << numtracer<< endl;
  vector<int> qflag(numtracer,0);
  vector<int> timetraj(numtracer,abs(tau));

  //Cell Variables
  vector<double> qflagini(numgridpoints,0);
  vector<double> qflagend(numgridpoints,0);
  vector<double> nini(numgridpoints,0);
  vector<double> nend(numgridpoints,0);
  vector<double> sumtend(numgridpoints,0);
  vector<double> sumt2end(numgridpoints,0);
  vector<double> cluster(numgridpoints,0);
  unsigned int step;
  double t;
  

  
  int depthdescnd=(clusterParams.domainmax.z-clusterParams.maxdepth)>0;
  cout << "depthdescnd=" << depthdescnd<<endl;
  

  
  omp_set_num_threads(10);
#pragma omp parallel for default(shared) private(t) // Parallelizing the code for computing trajectories

  for (q=0; q<numtracer; q++) {
    for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++) {
      //new point integration
      if(vequation==0){
	qflag[q]=RK4(t, h, &ftracer[q], GetVflowplusVsink);
      }
       else if(vequation==1){
	 qflag[q]=Heun2(t, h, &ftracer[q], GetVflowplusVsink);
       }
      
      if(qflag[q]==1){
	timetraj[q]=abs(t-tstart);
	break;
      }
      // Check maxdepth
      if(depthdescnd==1?(ftracer[q].z<=clusterParams.maxdepth):(ftracer[q].z>=clusterParams.maxdepth)){
	timetraj[q]=abs(t-tstart);
	break;
      }    
    }
  }
  
  // find index grid of initial position and final and add one to nini and nend in the corresponding grid box 
  
  vectorXYZ rini,rend;
  int iini,iend;
  int jini,jend;    
  int qini,qend;
  for (q=0; q<numtracer; q++){
    rini=(itracer[q]-clusterParams.domainmin)/clusterParams.intergrid;
    rend=(ftracer[q]-clusterParams.domainmin)/clusterParams.intergrid;
    // locate initial cell
    iini=floor(rini.x);
    jini=floor(rini.y);
    if(iini>=GridDimx || iini<0) continue;
    if(jini>=GridDimy || jini<0) continue;
    qini=iini+(jini*GridDimx);

    // locate arrival cell
    iend=floor(rend.x);
    jend=floor(rend.y);

    if(iend>=GridDimx || iend<0){
      qflagini[qini]+=1;
      continue;
    }
    if(jend>=GridDimy || jend<0){
      qflagini[qini]+=1;
      continue;
    } 
    qend=iend+(jend*GridDimx);

    if(qflag[q]==0){ 
      nini[qini]++;
      nend[qend]++;
      sumtend[qend]+=timetraj[q];
      sumt2end[qend]+=timetraj[q]*timetraj[q];      
    }
    else{
      qflagini[qini]+=qflag[q];
      qflagend[qend]+=qflag[q];
    }
  }
  
  for(q=0; q<numgridpoints; q++) {
    if(nini[q]>0)
      cluster[q]=nend[q]/nini[q];
    else
      cluster[q]=0.0;
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
 // save cluster values in a file

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
    rawname = "_lon"    + numprintf(4,1,clusterParams.domainmin.x) 
      + numprintf(4,1,clusterParams.domainmax.x)
      + numprintf(4,3,clusterParams.intergrid.x)
      + "_lat"       + numprintf(4,1,clusterParams.domainmin.y)
      + numprintf(4,1,clusterParams.domainmax.y)
      + numprintf(4,3,clusterParams.intergrid.y)
      + "_dpt"       + numprintf(4,0,clusterParams.domainmin.z)
      + numprintf(4,0,clusterParams.domainmax.z)
      + numprintf(4,0,clusterParams.intergrid.z)
      + "_ts"        + numprintf(3,2,clusterParams.intstep) 
      + "_mdpt"         + numprintf(4,0,clusterParams.maxdepth)
      + "_den" + numprintf(4,0,clusterParams.density)
      + "_v"         + numprintf(3,0,clusterParams.vsink)
      + "_veq"         + numprintf(1,0,vequation)
      + "_d"         + numprintf(2,0,clusterParams.seeddate.GetMday())    
      + numprintf(2,0,clusterParams.seeddate.GetMonth())
      + numprintf(2,0,clusterParams.seeddate.GetYear());
  }

  string namevtkfilegrid = "aggregation"+ rawname + ".vtk";
  string namefiletracer = "agg_tracer"+ rawname + ".dat";
  string namefileaggregation = "aggregation"+ rawname + ".dat";


#ifdef DEBUG
  cout << "Save cluster field in vtk file: " << namevtkfilegrid <<endl;
#endif

  ofstream vtkfile(namevtkfilegrid.c_str());
  vtkfile<<"# vtk DataFile Version 3.0"<<endl;
  vtkfile<<"grid aggregation 2D"<<endl; 
  vtkfile<<"ASCII"<<endl;
  vtkfile<<"DATASET STRUCTURED_POINTS"<< endl;
  vtkfile<<"DIMENSIONS "<< GridDimx+1 <<" "<< GridDimy+1 <<" "<< GridDimz+1 <<endl;
  vtkfile<<"ORIGIN "<<grid[0].x<<" "<<grid[0].y<<" "<< grid[0].z <<endl;
  vtkfile<<"SPACING "<<clusterParams.intergrid.x<<" "<<clusterParams.intergrid.y<<" "<< clusterParams.intergrid.z <<endl;
  vtkfile<<"CELL_DATA "<<(GridDimx)*(GridDimy)*(GridDimz)<<endl;
  vtkfile<<"SCALARS cluster float 1"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<cluster.size(); q++) {
      vtkfile<<cluster[q]<<endl;
  }
  vtkfile<<"SCALARS nini int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<nini.size(); q++) {
      vtkfile<<nini[q]<<endl;
  }
  vtkfile<<"SCALARS nend int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<nend.size(); q++) {
      vtkfile<<nend[q]<<endl;
  }
  vtkfile<<"SCALARS qflagini int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<qflagini.size(); q++) {
      vtkfile<<qflagini[q]<<endl;
  }
  vtkfile<<"SCALARS qflagend int"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<qflagend.size(); q++) {
      vtkfile<<qflagend[q]<<endl;
  }
  vtkfile<<"SCALARS sumt float 1"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++) {
      vtkfile<<sumtend[q]<<endl;
  }
  vtkfile<<"SCALARS sumt2 float 1"<<endl;
  vtkfile<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<numgridpoints; q++) {
      vtkfile<<sumt2end[q]<<endl;
  }
  vtkfile.close();

#ifdef TRACER
  cout << "Save data tracers in dat file: " << namefiletracer <<endl;

  ofstream filetracer(namefiletracer.c_str());
  for(q=0; q<numtracer; q++) {
    filetracer<<itracer[q]<<" ";
    filetracer<<ftracer[q]<<" ";
    filetracer<<timetraj[q]<<" ";
    filetracer<<qflag[q]<<endl;
  }
  filetracer.close();
#endif

#ifdef DEBUG
  cout << "Save data aggregation in dat file: " << namefileaggregation <<endl;
#endif

  ofstream fileaggregation(namefileaggregation.c_str());
  for(q=0; q<numgridpoints; q++) {
    fileaggregation<<grid[q]<<" ";
    fileaggregation<<nini[q]<<" ";
    fileaggregation<<nend[q]<<" ";
    fileaggregation<<qflagini[q]<<" ";
    fileaggregation<<qflagend[q]<<" ";
    fileaggregation<<sumtend[q]<<" ";
    fileaggregation<<sumt2end[q]<<" ";
    fileaggregation<<cluster[q]<<endl;
  }
  fileaggregation.close();

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
