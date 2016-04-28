#ifndef TURBO_RIEMANNSOLVERS_PERFECT_GAS_HLLC_3D
#define TURBO_RIEMANNSOLVERS_PERFECT_GAS_HLLC_3D

#include "datatypes.h"
#include "basetypes.h"
#include "grid.h"

class HLLC3DSolverPerfectGas {
	//Required info
	double gamma;
public:	
	HLLC3DSolverPerfectGas(){};

	void SetGamma(double _g) {
		gamma = _g;
	};

	//Numerical flux
	std::vector<double> F(ConservativeVariables U, Vector n)
	{		
		std::vector<double> res(5,0);		
		double ro = U.ro;
		double vx = U.rou/ro;
		double vy = U.rov/ro;
		double vz = U.row/ro;
		double roE = U.roE;	//ro*e
		double p = (gamma - 1.0)*(roE - ro*(vx*vx + vy*vy + vz*vz)/2.0);
		double vn = vx*n.x + vy*n.y + vz*n.z;

		res[0] = ro*vn;
		res[1] = ro*vn*vx + n.x*p;
		res[2] = ro*vn*vy + n.y*p;
		res[3] = ro*vn*vz + n.z*p;
		res[4] = vn*(roE + p);
		return res;
	};

	//Solve riemann problem by HLLC method
	std::vector<double> ComputeFlux(const ConservativeVariables& UL, const ConservativeVariables& UR, const Face& f) {
		std::vector<double> res(5,0);
		
		//left variables
		double rol = UL.ro;
		double ro_el = UL.roE - 0.5*(UL.rou*UL.rou + UL.rov*UL.rov + UL.row*UL.row)/rol;
		double pl = (gamma - 1.0)*ro_el;
		double hl = (UL.roE + pl)/rol;	//enthalpy

		//right variables
		double ror = UR.ro;
		double ro_er = UR.roE - 0.5*(UR.rou*UR.rou + UR.rov*UR.rov + UR.row*UR.row)/ror;
		double pr = (gamma - 1.0)*ro_er;
		double hr = (UR.roE + pr)/ror;	//enthalpy

		//compute normal and tangential components of velosities
		Vector velocity_l, velocity_r;
		velocity_l.x = UL.rou/UL.ro;
		velocity_l.y = UL.rov/UL.ro;
		velocity_l.z = UL.row/UL.ro;
		velocity_r.x = UR.rou/UR.ro;
		velocity_r.y = UR.rov/UR.ro;
		velocity_r.z = UR.row/UR.ro;
		double uln = velocity_l*f.FaceNormal;
		double urn = velocity_r*f.FaceNormal;
		double ult = (velocity_l - uln*f.FaceNormal).mod();
		double urt = (velocity_r - urn*f.FaceNormal).mod();

		//step 1 - compute minimum and maximum wave velocity S_min S_max
		double ql = sqrt(rol)/(sqrt(rol) + sqrt(ror));
		double qr = sqrt(ror)/(sqrt(rol) + sqrt(ror));
		double u_roe = ql*uln + qr*urn;	//Roe everage normal velocity
		double h_roe = ql*hl + qr*hr;	//average enthalpy
		double a_roe = sqrt((gamma - 1.0)*(h_roe - 0.5*u_roe*u_roe));	//average speed of sound
		double SL = u_roe - a_roe;
		double SR = u_roe + a_roe;

		//step 2 - compute wave velocity in star region S_st
		double Sst = pr - pl + pl*uln*(SL - uln) - pr*urn*(SR - urn);
		Sst /= rol*(SL - uln) - ror*(SR - urn);
		
		//check where is S* directed to
		if(Sst>=0) {
			//in new basis (y_new directed as velocity_l's tangential component)
			ConservativeVariables ULnew = UL;
			ULnew.rou = rol*uln;
			ULnew.rov = rol*ult;
			ULnew.row = 0;
			std::vector<double> FL = F(ULnew, f.FaceNormal);	//flux in left area in new coordinate systems f
			//trivial case
			if(SL>=0) {
				res = FL;
			} else {
				std::vector<double> UstL;
				double coef_l = rol*(SL - uln)/(SL - Sst);
				UstL.push_back(coef_l);
				UstL.push_back(coef_l*Sst);
				UstL.push_back(coef_l*ult);
				UstL.push_back(0);
				UstL.push_back(coef_l*(ULnew.roE/ULnew.ro + (Sst - uln)*(Sst + pl/(rol*(SL - uln)))));
				for(int i=0; i<res.size(); i++)
					res[i] = FL[i] + SL*(UstL[i] - ULnew[i]);
			};

			//recompute flux in original coordinate system
			Vector j_new = velocity_l - uln*f.FaceNormal;	//tangetial direction of left velocity vector
			if(j_new.mod()>0) j_new *= (1.0/j_new.mod());
			Vector vel_flux = res[1]*f.FaceNormal + res[2]*j_new;
			res[1] = vel_flux.x;
			res[2] = vel_flux.y;
			res[3] = vel_flux.z;
		};

		//second case
		if(Sst<0) {
			//in new basis (y_new directed as velocity_l's tangential component)
			ConservativeVariables URnew = UR;
			URnew.rou = ror*urn;
			URnew.rov = ror*urt;
			URnew.row = 0;
			std::vector<double> FR = F(URnew, f.FaceNormal);	//flux in right area in new coordinate systems f
			//trivial case
			if(SR<=0) {
				res = FR;
			} else {
				std::vector<double> UstR;
				double coef_r = ror*(SR - urn)/(SR - Sst);
				UstR.push_back(coef_r);
				UstR.push_back(coef_r*Sst);
				UstR.push_back(coef_r*urt);
				UstR.push_back(0);
				UstR.push_back(coef_r*(URnew.roE/URnew.ro + (Sst - urn)*(Sst + pr/(ror*(SR - urn)))));
				for(int i=0; i<res.size(); i++)
					res[i] = FR[i] + SR*(UstR[i] - URnew[i]);
			};

			//recompute flux in original coordinate system
			Vector j_new = velocity_r - urn*f.FaceNormal;	//tangetial direction of left velocity vector
			if(j_new.mod()>0) j_new *= (1.0/j_new.mod());
			Vector vel_flux = res[1]*f.FaceNormal + res[2]*j_new;		//MISTAKE
			res[1] = vel_flux.x;
			res[2] = vel_flux.y;
			res[3] = vel_flux.z;
		};
			
		MaxEigenvalue = max(abs(SL), abs(SR));
		return res;
	};

	//Public properties
	double MaxEigenvalue;
};

#endif