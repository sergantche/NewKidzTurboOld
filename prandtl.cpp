#include <stdio.h>
#include "optimization.h"

using namespace alglib;


//Task parameters
struct TaskParams {
	double nuLeft; //Right point coordinate
	double nuRight; //Left point coordinate
	int NCells; //Number of cells in grid
	double* nu; //Grid points coordinates
	double dNu; //Grid step
};

//Right hand side
double R(double x) {
	return 0;
};

//Target function
void  FRes(const real_1d_array &fi, real_1d_array &F, void *ptr)
{	
    TaskParams* taskParams = (TaskParams*)ptr;
	for (int i = 0; i<=taskParams->NCells; i++) {
		//Boundary conditions
		if (i == 0) {
			F[0] = fi[0];
			continue;
		}

		//Compute derivatives
		//Second
		double d2;
		if(i == taskParams->NCells-1) {
			d2 = (fi[i-1] - fi[i] + 2 * taskParams->dNu) / (taskParams->dNu * taskParams->dNu);
		} else d2 = (fi[i-1] - 2*fi[i] + fi[i+1]) / (taskParams->dNu * taskParams->dNu);

		//Third
		double d3 = 0;
		if (i == 1) {
			d3 = fi[3] - fi[1] - 2*fi[2];
			d3 /=2 * taskParams->dNu * taskParams->dNu * taskParams->dNu;
			F[1] = fi[1];
			continue;
		} else if (i == taskParams->NCells-2) {
			double f0 = (fi[i-2] + fi[i-1]) / 2.0;
			double f1 = (fi[i-1] + fi[i]) / 2.0;
			double f2 = (fi[i] + fi[i+1]) / 2.0;
			double f3 = f2 + 2*taskParams->dNu;
			d3 = -f0 + 3*f1 - 3*f2 + f3;
			d3 /= taskParams->dNu * taskParams->dNu * taskParams->dNu;
		} else if (i == taskParams->NCells-1) {
			double f0 = (fi[i-2] + fi[i-1]) / 2.0;
			double f1 = (fi[i-1] + fi[i]) / 2.0;
			double f2 = f1 + 2*taskParams->dNu;
			double f3 = f2 + 2*taskParams->dNu;
			d3 = -f0 + 3*f1 - 3*f2 + f3;
			d3 /= taskParams->dNu * taskParams->dNu * taskParams->dNu;
		} else {
			double f0 = (fi[i-2] + fi[i-1]) / 2.0;
			double f1 = (fi[i-1] + fi[i]) / 2.0;
			double f2 = (fi[i] + fi[i+1]) / 2.0;
			double f3 = (fi[i+1] + fi[i+2]) / 2.0;
			d3 = -f0 + 3*f1 - 3*f2 + f3;
			d3 /= taskParams->dNu * taskParams->dNu * taskParams->dNu;
		};

		//Equation to solve
		//F[i] = d3 - fi[i]*d2 - R(taskParams->nu[i]);	//with right-hand side
		F[i] = d3 - fi[i]*d2;
	};
	return;
}

//Main program ))
int main_(int argc, char *argv[]) {
	//Task parameters
	TaskParams taskParams;
	taskParams.NCells = 1000; //Number of cells
	taskParams.nuLeft = 0.0;
	taskParams.nuRight = 5.0;
	taskParams.dNu = (taskParams.nuRight - taskParams.nuLeft) / (taskParams.NCells-1);

	//Set computational grid
	taskParams.nu = new double[taskParams.NCells];
	double deltaNu = taskParams.dNu;
	for (int i = 0; i<taskParams.NCells; i++) {
		taskParams.nu[i] = taskParams.nuLeft + deltaNu * i;
	};	

	//Set initial guess
	real_1d_array fi;
	fi.setlength(taskParams.NCells); 
	for (int i = 0; i<taskParams.NCells; i++) {
		fi[i] = (2.0/3.0)* taskParams.nu[i] * taskParams.nu[i];
	};	

	//Set optimization parameters
    double epsg = 0.000001;
    double epsf = 0;
    double epsx = 0;
    ae_int_t maxits = 10000;
    minlmstate state;
    minlmreport rep;
	int size = taskParams.NCells;
	//int size = taskParams.NCells+1;

	//Run optimization
    minlmcreatev(size, fi, 0.0001, state);
    minlmsetcond(state, epsg, epsf, epsx, maxits);
    minlmoptimize(state, FRes, NULL, (void *)&taskParams);
    minlmresults(state, fi, rep);

	//Output result
    printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
    //printf("%s\n", fi.tostring(20).c_str()); // 

	FILE* fout = fopen("res.dat", "w");
	for (int i = 0; i<taskParams.NCells; i++) {
		fprintf(fout, "%lf %lf\n", taskParams.nu[i], fi[i]); // 
	};
	fclose(fout);
	std::getchar();
    return 0;
}