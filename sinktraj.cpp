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
struct sinktrajParameters {

  const vectorXYZ domainmin;
  const vectorXYZ domainmax;
  const vectorXYZ intergrid;
  const eqdate seeddate;
  const double intstep;
  const double vsink;
  const double maxdepth;

  // Define a constructor that will load stuff from a configuration file.
  sinktrajParameters(const string & sinktrajParamsFileName)
  :domainmin(getVectorXYZParam(sinktrajParamsFileName, "domainmin"))
  ,domainmax(getVectorXYZParam(sinktrajParamsFileName, "domainmax"))
  ,intergrid(getVectorXYZParam(sinktrajParamsFileName, "intergrid")) 
  ,seeddate(getEqDateParam(sinktrajParamsFileName, "seeddate")) 
  ,intstep(getDoubleParam(sinktrajParamsFileName, "intstep")) 
  ,vsink(getDoubleParam(sinktrajParamsFileName, "vsink"))
  ,maxdepth(getDoubleParam(sinktrajParamsFileName, "maxdepth"))
{}
};

int main(int argc, char **argv){

  /*********************************************************************************
   * READ PARAMETERS FROM FILE
   ********************************************************************************/

  string fnameparams;// File name that stores the parameters
  int namefileflag;
  int equation;

  if(GetcmdlineParameters(argc, argv, &fnameparams, &namefileflag, &equation)) {//Get cmd line parameters
    cout << "Error getting parameters from file" <<endl;
    return 1;
  }
  const sinktrajParameters sinktrajParams(fnameparams);


#ifdef DEBUG
  cout << "Sinktraj Parameters from file ";
  cout <<"\""<<fnameparams<<": "<<endl; 
  cout << " domainmin = "<< sinktrajParams.domainmin<<endl;
  cout << " intergrid = " <<sinktrajParams.intergrid<<endl;
  cout << " domainmax = "<< sinktrajParams.domainmax<<endl;
  cout << " seeddate = "<< sinktrajParams.seeddate.GetMday()<<"-";
  cout<<sinktrajParams.seeddate.GetMonth()<<"-";
  cout<<sinktrajParams.seeddate.GetYear()<<endl;
  cout << " intstep = "<<sinktrajParams.intstep<<endl;
  cout << " [Complete]" <<endl;
  cout << endl;
#endif

  /********************************************************************************
   * GRID CONSTRUCTION
   ********************************************************************************/

  vector<vectorXYZ> grid;
  int GridDimx, GridDimy,GridDimz;
 
  if(MakeRegularGrid(&grid, sinktrajParams.domainmin, 
		     sinktrajParams.intergrid, 
		     sinktrajParams.domainmax,
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

  int tau=2*int((sinktrajParams.domainmax.z-sinktrajParams.maxdepth)/sinktrajParams.vsink);

  cout << tau <<endl;
  eqdate inidate;

  double tstart;
  double tend;
  double h;

  int ntime = abs(tau);
  int ascnd = tau > 0;
  unsigned int seedtime = sinktrajParams.seeddate.GetDay();

  if(ascnd){
    unsigned int initime = seedtime;
    inidate.SetDay(initime);
    tstart = 0.0;
    tend = (double) ntime;
    h=sinktrajParams.intstep;
  }
  else{
    unsigned int initime = seedtime-ntime;
    inidate.SetDay(initime); 
    tend = 0.0;
    tstart = (double) ntime;
    h=-1.0*sinktrajParams.intstep;
  }

#ifdef DEBUG
  cout << "Setup time parameters: "<<endl;
  cout << " inidate = "<< inidate.GetMday() <<"-";
  cout<< inidate.GetMonth() <<"-";
  cout<< inidate.GetYear() <<endl;
  cout << " tstart = "<< tstart <<endl;
  cout << " tend = " << tend <<endl;
  cout << " h = "<< h <<endl;
#endif

  /**********************************************
   * LOAD VELOCITY FIELD
   **********************************************/

#ifdef DEBUG
  cout<<"Loading velocity grid:"<<endl; //Verbose: loading vel. grid
#endif

  SetupVflow(fnameparams);

  if((LoadVLatLonGrid(sinktrajParams.seeddate))!=0){//Load velocity grid
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

  vector<vectorXYZ> tracer = grid;
  vector<int> qflag(numgridpoints,0);

  unsigned int step;
  double t;
  unsigned int q;

  int depthdescnd=(sinktrajParams.domainmax.z-sinktrajParams.maxdepth)>0;

  omp_set_num_threads(7);
#pragma omp parallel for default(shared) private(t) // Parallelizing the code for computing trajectories

  for (q=0; q<numgridpoints; q++) {

    for(t=tstart, step=0; ascnd==1?(t<tend):(t>=tend); t+=h,step++) {
      //new point integration
      qflag[q]=RK4(t, h, &tracer[q], GetVflowplusVsink);
      if(qflag[q]==1) break;
      if(depthdescnd==1?(tracer[q].z<sinktrajParams.maxdepth):(tracer[q].z>=sinktrajParams.maxdepth)){
	qflag[q]=2;      
	break;
      }
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

  /**********************************************
   * DEBUGGING
   **********************************************/

#ifdef DEBUG
  ofstream fileotracers;

  fileotracers.open("final.pos");
  for (q=0; q<numgridpoints; q++) {
    if(qflag[q]==0)
      fileotracers << tracer[q] <<endl;
  }

  fileotracers.close();

#endif
  /********************************************************************************
   * WRITING RESULT IN VTK FILE
   ********************************************************************************/
 // save sinktraj values in a file

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
    rawname = "_lon"    + numprintf(4,1,sinktrajParams.domainmin.x) 
      + numprintf(4,1,sinktrajParams.domainmax.x)
      + numprintf(4,3,sinktrajParams.intergrid.x)
      + "_lat"       + numprintf(4,1,sinktrajParams.domainmin.y)
      + numprintf(4,1,sinktrajParams.domainmax.y)
      + numprintf(4,3,sinktrajParams.intergrid.y)
      + "_dpt"       + numprintf(4,0,sinktrajParams.domainmin.z)
      + numprintf(4,0,sinktrajParams.domainmax.z)
      + numprintf(4,0,sinktrajParams.intergrid.z)
      + "_ts"        + numprintf(3,2,sinktrajParams.intstep) 
      + "_d"         + numprintf(2,0,sinktrajParams.seeddate.GetMday())
      + numprintf(2,0,sinktrajParams.seeddate.GetMonth())
      + numprintf(2,0,sinktrajParams.seeddate.GetYear());
  }

  string namevtkfileitrac = "itrac"+ rawname + ".vtk";
  string namevtkfileftrac = "ftrac"+ rawname + ".vtk";

#ifdef DEBUG
  cout << "Save initial positions in vtk file: " << namevtkfileitrac <<endl;
#endif

  vectorXYZ scalefactor(1.0,1.0,0.001);

  ofstream vtkfileitrac(namevtkfileitrac.c_str());
  vtkfileitrac<<"# vtk DataFile Version 3.0"<<endl;
  vtkfileitrac<<"Initial positions"<<endl; 
  vtkfileitrac<<"ASCII"<<endl;
  vtkfileitrac<<"DATASET POLYDATA"<< endl;
  vtkfileitrac<<"POINTS "<< grid.size() << " float"<<endl;
  for(q=0; q<grid.size(); q++) {
      vtkfileitrac<<grid[q]*scalefactor<<endl;
  }
  vtkfileitrac<<"VERTICES "<< grid.size() <<" " << 2*grid.size() << endl;
  for(q=0; q<grid.size(); q++) {
    vtkfileitrac<<"1 "<<q<<endl;
  }
  vtkfileitrac<<"POINT_DATA "<< tracer.size() << endl;
  vtkfileitrac<<"SCALARS qflaq int"<<endl;
  vtkfileitrac<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<qflag.size(); q++) {
    vtkfileitrac<<qflag[q]<<endl;
  }

  vtkfileitrac.close();

#ifdef DEBUG
  cout << "Save final positions in vtk file: " << namevtkfileftrac <<endl;
#endif

  ofstream vtkfileftrac(namevtkfileftrac.c_str());
  vtkfileftrac<<"# vtk DataFile Version 3.0"<<endl;
  vtkfileftrac<<"Initial positions"<<endl; 
  vtkfileftrac<<"ASCII"<<endl;
  vtkfileftrac<<"DATASET POLYDATA"<< endl;
  vtkfileftrac<<"POINTS "<< tracer.size() << " float"<<endl;
  for(q=0; q<tracer.size(); q++) {
    vtkfileftrac<<tracer[q]*scalefactor<<endl;
  }
  vtkfileftrac<<"VERTICES "<< tracer.size() <<" " << 2*tracer.size() << endl;
  for(q=0; q<tracer.size(); q++) {
    vtkfileftrac<<"1 "<<q<<endl;
  }
  vtkfileftrac<<"POINT_DATA "<< tracer.size() << endl;
  vtkfileftrac<<"SCALARS qflag int"<<endl;
  vtkfileftrac<<"LOOKUP_TABLE default"<<endl;
  for(q=0; q<qflag.size(); q++) {
    vtkfileftrac<<qflag[q]<<endl;
  }
  vtkfileftrac.close();

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
