


int SetupLagrangianEngine(const string & LagEngParamsFileName);
int GetVflowplusVsink(double t,vectorXYZ point, vectorXYZ *vint);
int RK4(double t0, double intstep, vectorXYZ *point, int (*velocity)(double t,vectorXYZ point, vectorXYZ *vint));
int Heun2(double t0, double intstep, vectorXYZ *point, int (*velocity)(double ,vectorXYZ , vectorXYZ* ));
/* TEST RANDOM NUMBER GENERATOR
int mttest(void);
*/
