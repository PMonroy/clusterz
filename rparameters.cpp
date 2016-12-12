#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>
#include <unistd.h>

using namespace std;
#include "vectorXYZ.hpp"
#include "eqdate.hpp"

string getStringParam(const string & paramsFileName, const string & key) {
  ifstream ifile(paramsFileName.c_str());
  if(!ifile.is_open()){
    cout<<"string Skipping unreadable file " << paramsFileName.c_str() <<" "<<endl; 
    exit(EXIT_FAILURE);
  }
  string sline;
  stringstream ssline;
  string nameVar;
  string Var;

  while(!ifile.eof()){
    getline(ifile,sline,'\n');
    if (sline[0]=='#') {
      continue;  /* Lines starts with # are comments and are ignored*/
    }
    ssline<<sline;
    ssline>>nameVar;
    if(nameVar==key){
      ssline>>Var;
      ifile.close();
      return Var;
      break;
    }
    else{
      ssline.clear();//clear any bits set
      ssline.str(string());
    }
  }
  ifile.close();
  cout << "Parameter "<<key<<" not set in file "<< paramsFileName<< endl;
  exit(EXIT_FAILURE);
}
int getIntParam(const string & paramsFileName, const string & key) {
  ifstream ifile(paramsFileName.c_str());
  if(!ifile.is_open()){
    cout<<" Skipping unreadable file " << paramsFileName.c_str() <<" "<<endl; 
    exit(EXIT_FAILURE);
  }

  string sline;
  stringstream ssline;
  string nameVar;
  int Var;

  while(!ifile.eof()){
    getline(ifile,sline,'\n');
    if (sline[0]=='#') {
      continue;  /* Lines starts with # are comments and are ignored*/
    }
    ssline<<sline;
    ssline>>nameVar;
    if(nameVar==key){
      ssline>>Var;
      ifile.close();
      return Var;
      break;
    }
    else{
      ssline.clear();//clear any bits set
      ssline.str(string());
    }
  }
  ifile.close();
  cout << "Parameter "<<key<<" not set in file "<< paramsFileName<< endl;
  exit(EXIT_FAILURE);
}
double getDoubleParam(const string & paramsFileName, const string & key) {
  ifstream ifile(paramsFileName.c_str());
  if(!ifile.is_open()){
    cout<<"Skipping unreadable file " << paramsFileName.c_str() <<" "<<endl; 
    exit(EXIT_FAILURE);
  }
  string sline;
  stringstream ssline;
  string nameVar;
  double Var;

  while(!ifile.eof()){
    getline(ifile,sline,'\n');
    if (sline[0]=='#') {
      continue;  /* Lines starts with # are comments and are ignored*/
    }
    ssline<<sline;
    ssline>>nameVar;
    if(nameVar==key){
      ssline>>Var;
      ifile.close();
      return Var;
      break;
    }
    else{
      ssline.clear();//clear any bits set
      ssline.str(string());
    }
  }
  ifile.close();
  cout << "Parameter "<<key<<" not set in file "<< paramsFileName<< endl;
  exit(EXIT_FAILURE);
}
vectorXYZ getVectorXYZParam(const string & paramsFileName, const string & key) {
  ifstream ifile(paramsFileName.c_str());
  if(!ifile.is_open()){
    cout<< key <<" Vector Skipping unreadable file " << paramsFileName.c_str() <<" "<<endl; 
    exit(EXIT_FAILURE);
  }
  string sline;
  stringstream ssline;
  string nameVar;
  vectorXYZ Var;

  while(!ifile.eof()){
    getline(ifile,sline,'\n');
    if (sline[0]=='#') {
      continue;  /* Lines starts with # are comments and are ignored*/
    }
    ssline<<sline;
    ssline>>nameVar;
    if(nameVar==key){
      ssline>>Var;
      ifile.close();
      return Var;
      break;
    }
    else{
      ssline.clear();//clear any bits set
      ssline.str(string());
    }
  }
  ifile.close();
  cout << "Parameter "<<key<<" not set in file "<< paramsFileName<< endl;
  exit(EXIT_FAILURE);
}
struct tm getDateParam(const string & paramsFileName, const string & key) {
  ifstream ifile(paramsFileName.c_str());
  if(!ifile.is_open()){
    cout<<"Skipping unreadable file " << paramsFileName.c_str() <<" "<<endl; 
    exit(EXIT_FAILURE);
  }
  string sline;
  stringstream ssline;
  string nameVar;
  struct tm Var;

  while(!ifile.eof()){
    getline(ifile,sline,'\n');
    if (sline[0]=='#') {
      continue;  /* Lines starts with # are comments and are ignored*/
    }
    ssline<<sline;
    ssline>>nameVar;
    if(nameVar==key){
      ssline>>Var.tm_mday>>Var.tm_mon>>Var.tm_year;
      Var.tm_mon = Var.tm_mon - 1;
      Var.tm_hour= 11;
      Var.tm_min = 0;
      Var.tm_sec = 0;
      ifile.close();
      return Var;
      break;
    }
    else{
      ssline.clear();//clear any bits set
      ssline.str(string());
    }
  }
  ifile.close();
  cout << "Parameter "<<key<<" not set in file "<< paramsFileName<< endl;  
  exit(EXIT_FAILURE);
}
eqdate getEqDateParam(const string & paramsFileName, const string & key) {
  ifstream ifile(paramsFileName.c_str());
  if(!ifile.is_open()){
    cout<<"Skipping unreadable file " << paramsFileName.c_str() <<" "<<endl; 
    exit(EXIT_FAILURE);
  }
  string sline;
  stringstream ssline;
  string nameVar;
  eqdate Var;

  unsigned int mday,month,year;

  while(!ifile.eof()){
    getline(ifile,sline,'\n');
    if (sline[0]=='#') {
      continue;  /* Lines starts with # are comments and are ignored*/
    }
    ssline<<sline;
    ssline>>nameVar;
    if(nameVar==key){
      ssline>>mday>>month>>year;
      Var.SetDate(year,month,mday);
      ifile.close();
      return Var;
      break;
    }
    else{
      ssline.clear();//clear any bits set
      ssline.str(string());
    }
  }
  ifile.close();
  cout << "Parameter "<<key<<" not set in file "<< paramsFileName<< endl;  
  exit(EXIT_FAILURE);
}


int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams, int *namefileflag, int *equation){

  int opt;
  int fnameparamsflag=0; // File parameters flag option
  int equationflag=0; // Equation type flag option
  *namefileflag=0;

  while((opt = getopt(narg, cmdarg, "f:ne:h")) != -1){
    switch(opt){
    case 'f':
      fnameparamsflag++;
      if(optarg) 
	*fnameparams = optarg;
      else 
	return 1;
      break;
    case 'n':
      *namefileflag=1;
      break;
    case 'e':
      equationflag++;
      if(optarg) 
	*equation = atoi(optarg);
      else 
	return 1;
      break;
    case 'h':
      cout<<"Usage: "<< cmdarg[0] <<" [OPTIONS]"<<endl;
      cout<<" -f [file]      Where [file] is the input file parmeters (mandatory)" <<endl;
      cout<<" -e [int]       Where [file] is the velocity equation (mandatory): 0=not random(RK4), 1=random(Heun)" <<endl;
      cout<<" -n             Output file names will write with input file without extension (optional)"<<endl;
      cout<<" -h             Print this help and exit (optional)"<<endl;
      return 1;

    case '?':
      cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
      return 1;
    }
  }
  
  /* Check mandatory or repeated options*/

  if(fnameparamsflag==0){//Check if file option is set  
    cout << "Error: Option -f <file> is mandatory" <<endl;
    cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
    return 1;
  }
  else if (fnameparamsflag>1){//Check if file option is repeated
    cout << "Warning: Option -f is repeated. Parameters file now is the last one ("<<*fnameparams<<")"<<endl;
  }

  if(equationflag==0){//Check if equation option is set  
    cout << "Error: Option -e [int] is mandatory" <<endl;
    cout << "Try "<<cmdarg[0]<<" -h for more information"<<endl;
    return 1;
  }
  else if (equationflag>1){//Check if equation option is repeated
    cout << "Warning: Option -e is repeated. Equation type now is the last one ("<<*equation<<")"<<endl;
  }

  return 0;
}
