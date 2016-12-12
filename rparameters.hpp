#ifndef PARAMETERS
#define PARAMETERS

#include <ctime>
#include <string>

using namespace std;
#include "vectorXYZ.hpp"
#include "eqdate.hpp"

int GetcmdlineParameters(int narg,char ** cmdarg, string *fnameparams,int *namefileflag, int *equation);

int getIntParam(const string & paramsFileName, const string & key);
string getStringParam(const string & paramsFileName, const string & key);
double getDoubleParam(const string & paramsFileName, const string & key);
vectorXYZ getVectorXYZParam(const string & paramsFileName, const string & key);
struct tm getDateParam(const string & paramsFileName, const string & key);
eqdate getEqDateParam(const string & paramsFileName, const string & key);

#endif
