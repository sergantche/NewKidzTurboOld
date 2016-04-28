#ifndef TURBO_TEMPORARY_TESTS
#define TURBO_TEMPORARY_TESTS

#include "model.h"

#define TEMPORARY_TESTS_LY_TR 1.0		////for RunVortexTrTest2D test


				// TROSHKIN VORTEX TEST
struct Rectangle_Area_TR {
	//coordinates of rectangle nodes
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
	//specific heat ratio
	double gamma;
	//initial pressure
	double P0;
	//initial density
	double RO0;
	//initial temperature
	double T0;
};

//values in dummy cells for 1D case
Vector BoundaryDistribution1D(Vector r, double time) {
	double time_scale = 0.0005;
	double init_vel = 500;
	return Vector(init_vel*exp((-1)*time/time_scale), 0, 0);
};

//Troshkin task 1D (impulse collision of the gas)
void RunVortexTrTest1D(int N_x, double delta_time_to_save, double Max_time)
{
	//set task parameters
	Rectangle_Area_TR info;

	//Fe Pe discrete perturbation
	info.gamma = 1.4;
	info.P0 = 1.0e5;	//initial pressure
	info.RO0 = 1.16;	//initial density
	info.T0 = 300.0;	//initial temperature
	info.x_max = 1.0;
	info.x_min = 0;
	//Shock tube problem setting
	Vector direction = Vector(1,0,0);
	Grid grid = GenGrid1D(N_x, info.x_min, info.x_max, direction);
	Model<Roe3DSolverPerfectGas> model;
	
	//Set fluid properties	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.DisableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.01);
	model.SetSchemeOrder(2);

	//Bind computational grid
	model.BindGrid(grid);
	if(model.SchemeOrder > 1) model.EnableLimiter();
	else model.DisableLimiter();

	//Set initial conditions
	ConservativeVariables initValues = model.PrimitiveToConservativeVariables(Vector(0, 0 , 0), info.P0, info.T0, model.medium);
	model.SetInitialConditions(initValues);

	//free boundary condition
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition Symmetry(model);
	Model<Roe3DSolverPerfectGas>::DummyTimeSpaceDistribution DummyDistribution(model);
	DummyDistribution.setParams(info.P0, info.T0, BoundaryDistribution1D);

	model.SetBoundaryCondition("left", DummyDistribution);
	model.SetBoundaryCondition("right", NaturalBC);
	//Load solution
	std::string outputSolutionFile = "vort1D";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	int n_delta = (int)(model.totalTime/delta_time_to_save) + 1;
	double TimeToSave = n_delta*delta_time_to_save;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
		};
		if (model.totalTime > Max_time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_time_to_save;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

//values in dummy cells for 2D case
Vector BoundaryDistribution2D(Vector r, double time) {
	double time_scale = 0.0005;
	double init_main_vel = 100;
	double init_pert_vel = -0.1*init_main_vel*cos(2*PI*r.y/TEMPORARY_TESTS_LY_TR);
	double u = init_main_vel + init_pert_vel;

	return Vector(u*exp((-1)*time/time_scale), 0, 0);
};

//Troshkin task 2D (impulse collision of the gas cause by velocity distribution in dummy cells)
void RunVortexTrTest2D(int N_x, int N_y, double delta_time_to_save, double Max_time)
{
	//set task parameters
	Rectangle_Area_TR info;

	info.gamma = 1.4;
	info.P0 = 1.0e5;	//initial pressure
	info.RO0 = 1.16;	//initial density
	info.T0 = 300.0;	//initial temperature
	info.x_max = 1.0;
	info.x_min = 0;
	info.y_min = -0.5*TEMPORARY_TESTS_LY_TR;
	info.y_max = -info.y_min;
	
	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);
	
	//Set fluid properties	
	model.SetGamma(1.4);
	model.SetCv(1006.43 / 1.4);
	model.SetMolecularWeight(28.966);
	model.DisableViscous();

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.01);
	model.SetSchemeOrder(2);

	//Bind computational grid
	model.BindGrid(grid);
	if(model.SchemeOrder > 1) model.EnableLimiter();
	else model.DisableLimiter();

	//Set initial conditions
	ConservativeVariables initValues = model.PrimitiveToConservativeVariables(Vector(0, 0 , 0), info.P0, info.T0, model.medium);
	model.SetInitialConditions(initValues);

	//free boundary condition
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition Symmetry(model);
	Model<Roe3DSolverPerfectGas>::DummyTimeSpaceDistribution DummyDistribution(model);
	DummyDistribution.setParams(info.P0, info.T0, BoundaryDistribution2D);

	model.SetBoundaryCondition("left", DummyDistribution);
	model.SetBoundaryCondition("right", NaturalBC);
	model.SetBoundaryCondition("top", Symmetry);
	model.SetBoundaryCondition("bottom", Symmetry);
	//Load solution
	std::string outputSolutionFile = "vort2D";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	int n_delta = (int)(model.totalTime/delta_time_to_save) + 1;
	double TimeToSave = n_delta*delta_time_to_save;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
		};
		if (model.totalTime > Max_time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_time_to_save;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

// Demchenko's article test

//2D initializing of init conditions for Test 2 (Demchenko)
ConservativeVariables FePlCollision_Demchenko_Init2D(Vector CellPos, void* params) {
	Rectangle_Area_RT info = *(Rectangle_Area_RT*)params;
	
	//define density pressure and velocity in the cell
	double ro, p;
	Vector v;

	//common border
	double l_y = info.y_max;		//a half of y size of the domain
	double x_zero_level = 0;		//position of free border without disturbances
	double bound_pos = 0;				//initial border distribution around zero level
	//distort bottom half of the domain
	if (CellPos.y < 0) {
		bound_pos = 1 - cos(CellPos.y*PI/l_y);
		bound_pos *= 0.5*l_y;
	};
	bound_pos += x_zero_level;	//position of distorted free porder
	
	//density left and right states
	if (CellPos.x <= bound_pos) ro = info.ro_bot;		
	else ro = info.ro_top;

	//pressure
	p = info.P0;

	//velocity
	v.y = 0;
	v.z = 0;
	if(CellPos.x <= bound_pos) v.x = 0.5*info.A;
	else v.x = -0.5*info.A;

	ConservativeVariables U;
	U.ro = ro;
	U.rou = ro * v.x;
	U.rov = ro * v.y;
	U.row = ro * v.z;
	U.roE = p/(info.gamma - 1.0) + ro * v.mod() * v.mod() / 2.0;
	
	return U;	
};

//2D collision test 2 from Demchenko
void RunCollisionDemchenkoTest2D(int N_x, int N_y, double delta_time_to_save, double Max_time)
{
	//set task parameters
	Rectangle_Area_RT info;

	//Fe Pb discrete perturbation
	info.A = 5.0e2;		//velocity of right plate
	info.gamma = 5.0/3.0;		
	info.P0 = 4.0e8;		
	info.ro_bot = 7.9e3;		//left state - Steel
	info.ro_top = 11.34e3;	//	right state - Lead
	//info.x_max = 3.9e-3;			//right plate width
	//info.x_min = -4.68e-3;			//left plate width
	info.x_max = 3.0e-3;
	info.x_min = - info.x_max;
	info.y_max = 1.0e-3;			//radius of perturbation
	info.y_min = - info.y_max;

	Model<Roe3DSolverPerfectGas> model;
	Grid grid = GenGrid2D(N_x, N_y, info.x_min, info.x_max, info.y_min, info.y_max, 1.0, 1.0);

	//Set fluid properties
	model.SetGamma(info.gamma);
	model.SetCv(1006.43 / info.gamma);
	model.SetMolecularWeight(28.966);
	model.SetThermalConductivity(0.0);
	model.DisableViscous();

	Vector g(0, -info.g, 0);
	model.SetUniformAcceleration(g);

	//Set computational settings
	model.SetCFLNumber(0.35);
	model.SetHartenEps(0.1);
	model.SchemeOrder = 2;
	model.EnableLimiter();
	//model.DisableLimiter();

	//Bind computational grid
	model.BindGrid(grid);	

	//Set initial conditions
	model.SetInitialConditions(FePlCollision_Demchenko_Init2D, &info);
		
	//Boundary conditions
	//Symmetry boundary
	Model<Roe3DSolverPerfectGas>::SymmetryBoundaryCondition SymmetryBC(model);
	//No slip boundary
	Model<Roe3DSolverPerfectGas>::NaturalCondition NaturalBC(model);

	//Set boundary conditions
	model.SetBoundaryCondition("left", NaturalBC);
	model.SetBoundaryCondition("right", NaturalBC);
	model.SetBoundaryCondition("bottom", SymmetryBC);
	model.SetBoundaryCondition("top", SymmetryBC);

	//Save initial solution
	model.ComputeGradients();
	model.SaveToTechPlot("init.dat");

	//Load solution
	std::string outputSolutionFile = "collision2D";
	//model.LoadSolution("solution_to_load.txt");

	//Run simulation
	int n_delta = (int)(model.totalTime/delta_time_to_save) + 1;
	double TimeToSave = n_delta*delta_time_to_save;
	
	bool isSave = true;	
	for (int i = 0; i < 1000000; i++) {
		model.Step();	
		if (i % 10 == 0) {
			std::cout<<"Iteration = "<<i<<"\n";
			std::cout<<"TimeStep = "<<model.stepInfo.TimeStep<<"\n";
			for (int k = 0; k<5; k++) std::cout<<"Residual["<<k<<"] = "<<model.stepInfo.Residual[k]<<"\n";
			std::cout<<"TotalTime = "<<model.totalTime<<"\n";
		};
		if ((i % 100 == 0) && (isSave)) {
			model.SaveSolution(outputSolutionFile+".txt");
			model.SaveToTechPlot(outputSolutionFile+".dat");
		};
		if (model.totalTime > Max_time) break;
		//save solutions during computation
		if(model.totalTime > TimeToSave) {
			std::stringstream ss;
			ss << n_delta;
			std::string timer = ss.str();
			model.SaveToTechPlot(outputSolutionFile+timer+".dat");
			TimeToSave += delta_time_to_save;
			n_delta++;
		};
	};

	//Save result to techplot
	if (isSave) {
		model.SaveSolution(outputSolutionFile+".txt");
		model.SaveToTechPlot(outputSolutionFile+".dat");
	};

	return;
};

#endif