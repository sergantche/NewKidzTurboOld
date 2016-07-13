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

// Compute gradient using least squares
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

// !Rough algorithm for space reconstruction
Vector ComputeLinearGradient(Vector point, double value, const std::vector<Vector>& points, const std::vector<double>& values, int nDims) {
	
	auto eps_min{ std::numeric_limits<double>::min() };
	Vector res;

	// Find the most "right" and the most "left" points in all direction and compute components of the gradient
	auto pos_max{ 0.0 };
	auto pos_min{ 0.0 };
	auto val_r{ 0.0 };
	auto val_l{ 0.0 };
	for (auto i = 0; i < points.size(); i++) {
		auto delta_x = points[i].x - point.x;
		if (delta_x > pos_max + eps_min) {
			pos_max = delta_x;
			val_r = values[i];
		};
		if (delta_x < pos_min - eps_min) {
			pos_min = delta_x;
			val_l = values[i];
		};
	};
	if(pos_max - pos_min > 5.0 * eps_min)	res.x = (val_r - val_l) / (pos_max - pos_min);

	// Repeat for 2D and 3D case (TO DO)
	
	
	return res;
};

#endif