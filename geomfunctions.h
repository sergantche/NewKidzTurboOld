#ifndef TURBO_GEOMFUNCTIONS
#define TURBO_GEOMFUNCTIONS

//Functions operating geometrical objects

#include <math.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "basetypes.h"
#include "cgnslib.h"
#include "optimization.h"
#include <Eigen\Dense>


//Compute cgns element measure ( volume, square, length ) 
double ComputeElementMeasure() {
	return 0;
};

//Compute gradient using least squares
Vector ComputeGradientByPoints(Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values, int nDims) {
	Vector grad;

	//Input	
	int m = points.size(); //Number of equations
	int nrhs = 1; //Number of right hand side
	Eigen::MatrixXd A(m, nDims);
	Eigen::VectorXd rhs(nrhs*m);

	//Compose matrices 				
	for (int i = 0; i<points.size(); i++) {
		Vector dr = point - points[i];
		A(i, 0) = dr.x;
		if (nDims > 1) A(i, 1) = dr.y;
		if (nDims > 2) A(i, 2) = dr.z;
	};

	for (int i = 0; i<points.size(); i++) {
		double dU = value - values[i];
		rhs(i) = dU;
	};

	//Solve problem	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
	Eigen::VectorXd x = svd.solve(rhs);

	//Check rank deficiency case
	if (svd.rank() < nDims) throw std::exception("Could not solve for gradient due to rank deficiency");

	//Output
	grad = Vector(x(0), 0, 0);
	if (nDims > 1) grad.y = x(1);
	if (nDims > 2) grad.z = x(2);

	return grad;
};

#endif