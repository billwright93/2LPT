//COMPILE + RUN COMMANDS//

//g++ -I/opt/apps/libs/gsl/2.1/gcc-4.4.7/include -I/users/wrightw/work/cpplibs/headers -L/opt/apps/libs/gsl/2.1/gcc-4.4.7/lib -lgsl -lstdc++ -lgslcblas kernelb32.cpp /users/wrightw/work/cpplibs/cpps/*.cpp
// ./a.out


#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include "spline.h"

#include "stdafx.h"
#include <stdio.h>
#include "interpolation.h"



// CONSTANTS //
const double pipi = 3.1415926535897;
const	double H0 = (3.09e+19);//70.;//
double omega0 = 0.25;//0.99;//


// MG SWITCHES //
int MG = 1; //0=GR, 1=f(R), 2=DGP


// SAMPLING RANGES FOR a, k(and k1), x //
// Number of time steps you wish to output between a_ini and a_fin
const int Na = 10;
double a_ini = 0.1;
double a_fin = 1.1;

// Number of k-modes you wish to output between kmin and kmax
const int Nk = 4; //10
double kmin = 0.01;//0.01;
double kmax = 10;//0.2;

// Number of x=cos(theta) modes you wish to output between xmin and xmax
const int Nx = 1;//4; //10
double xmin = -1.+1.e-10;
double xmax = 1.-1.e-10;



// DECLARING VARIABLES AND ARRAYS //
double a_arr[Na];
double k_arr[Nk];
double D1_arr[Na*Nk];

alglib::real_1d_array aa_arr;
alglib::real_1d_array kk_arr;
alglib::real_1d_array D1D1_arr;

//D1/D1_dot(k)
double D1;
double dD1;

//k.L2(k, k1, theta)
double kL2[Nk*Nx];
double dkL2[Nk*Nx];



// FUNCTIONS //

//Easy power definitions
static double pow2(double base){
	return pow(base, 2.);}

static double pow3(double base){
	return pow(base, 3.);}

static double pow4(double base){
	return pow(base, 4.);}


//D1 2d spline interpolation function
double calc_D1_2d(double a, double k, alglib::real_1d_array& a_arr, alglib::real_1d_array& k_arr, alglib::real_1d_array& D1_arr){
	// k=x, a=y
	double v;
	alglib::spline2dinterpolant s;
	alglib::spline2dbuildbicubicv(k_arr, k_arr.length(), a_arr, a_arr.length(), D1_arr, 1, s); //1=dimension of D1 values since scalar field
	v = alglib::spline2dcalc(s, k, a);
	return v;
}


// Normalized Hubble parameter H, a*H*dH/da, and Omega_M(a) for the background (LCDM)

//H(a)
static double HA(double a, double OmegaM0){
	double OmegaL0 = 1.- OmegaM0;
	return  H0*sqrt((OmegaM0/pow3(a))+OmegaL0);}

// dH(a)/da
static double dHA(double a, double OmegaM0){
	double OmegaL0 = 1.- OmegaM0;
	return -1.5*(OmegaM0/pow4(a))*(1./((OmegaM0/pow3(a))+OmegaL0))*HA(a, OmegaM0);}

//
static double Omega_M(double a, double OmegaM0){
	//return OmegaM0*pow2(H0)/(pow2(HA(a, OmegaM0))*pow3(a));}
	double OmegaL0 = 1. - OmegaM0;
	return OmegaM0/(OmegaM0+(OmegaL0*pow3(a)));}


  //MODIFIED GRAVITY MODEL FUNCTIONS:
  //p1,p2,p3 are theory parameters
  //k1,k2,k3 are the mode magnitudes
  //u1 angle between k2 and k3; u2 between k1 and k2; u3 between k1 and k3


  // default is 3 model params. Add in parameters p4,p5,... as required.
  // Note on theories with more than 3 parameters:
  // 1) Add in extra parameters to param_type3 structure below
  // 2) Add in extra parameters to params of ode1 system of ODEs
  // 3) Add in extra parameters to initn function and in params passed by function to solver (this is right after initial conditions are passed)

//1-mu^2
static double ker1(double u1){
	return 1.-u1*u1;
}

//beta function for nDGP
static double beta(double a, double OmegaM0, double omegarc){
	return 1.+(HA(a, OmegaM0)/H0)/sqrt(omegarc)*(1.+(-1.5*OmegaM0/pow(a,3))/(3.*(HA(a,OmegaM0)*HA(a,OmegaM0)/pow2(H0))));
}

//mu
static double mu(double a, double k0, double OmegaM0, double p1, double p2, double p3, int mg){
	double h0 = 1./2997.9;
	if(mg == 0){
		return 1.; // GR
	}
	else if(mg == 1){
		return 1. + pow2(k0/a)/(3.*(pow2(k0/a)+pow3(OmegaM0/pow3(a)-4.*(OmegaM0-1.))/(2.*p1/pow2(h0)*pow2(4-3.*OmegaM0)))); //f(R) Hu- Sawicki
	}
	else if(mg == 2){
		return 1.+1./(3.*beta(a,OmegaM0,p1)); //nDGP
	}
	else{
		std::cout << "WARNING: MG SWITCH NOT DEFINED CORRECTLY, REVERTING TO GR" << '\n';
		return 1;
	}
}

	//gamma 2
static double gamma2(double a, double OmegaM0, double k0, double k1, double k2, double u1, double p1, double p2, double p3, int mg){
  double h0 = 1./2997.9;
	if(mg == 0){
		return 0. ; // GR
	}
	else if(mg == 1){
		return -(9.*pow2(k0/a)*pow2(OmegaM0/pow3(a))*pow(OmegaM0-4.*pow3(a)*(-1+OmegaM0),5))/
	  			    (48.*pow(a,15)*pow2(p1/pow2(h0))*pow2(HA(a,OmegaM0))*pow4(3.*OmegaM0-4.)
	  			   *(pow2(k0/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
	  			   *(pow2(k1/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
	  			   *(pow2(k2/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))); //f(R) Hu- Sawicki
	}
	else if(mg == 2){
		return -1.*p2/(HA(a,OmegaM0)*HA(a,OmegaM0)*24.*pow(beta(a,OmegaM0,p1),3)*p1)*pow(OmegaM0/(a*a*a),2)*ker1(u1); //nDGP
	}
	else{
		std::cout << "WARNING: MG SWITCH NOT DEFINED CORRECTLY, REVERTING TO GR" << '\n';
		return 0;
	}
}


//Jacobian of the system required when calling the system evolver, the below is neither needed, nor correct for solving
int jac (double a, const double D[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}




	//// 1ST ORDER ////

  struct param_type3 {
    double kk;
  	double OmegaM00;
  	double par1;
  	double par2;
  	double par3;
		int mmgg;
  };

  int ode1(double a, const double D[], double dDda[], void *params)
  {
  	param_type3 p = *(param_type3 *)(params);
  	double k = p.kk;
  	double OmegaM0 = p.OmegaM00;
  	double p1 = p.par1;
  	double p2 = p.par2;
  	double p3 = p.par3;
		int mg = p.mmgg;

		//Equations for 1st order growth factor D1/dD1(k1)
		dDda[0] = D[1];
		dDda[1] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3, mg)*D[0];

  	return GSL_SUCCESS;
  }


  // EDIT : Add in new gravity parameters to function's input as required (par1 is the nDGP parameter omega_rc as default)
  // EDIT : if more than 1 gravity parameter is needed, one must also edit initn function in SPT.h file.

  int solve1(double A, double k, double OmegaM0, double par1, double par2, double par3, int mg)
  {
  //#pragma omp parallel for schedule(dynamic)

		// Initial scale factor for solving system of equations
		double a = a_ini;

		// Einstein de Sitter initial condtions for 1st order growth factor and derivative
		double D[2] = { 1., 1./a };//{a, 1.}; //

		/*Parameters passed to system of equations */
	  // EDIT : If more than one gravity parameter is used, add them after p1
		struct param_type3 my_params1 = {k, OmegaM0, par1, par2, par3, mg};

		gsl_odeiv2_system sys = {ode1, jac, 2, &my_params1};

		//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
		gsl_odeiv2_driver * d =
		gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										1e-6, 1e-6, 0.0);

		int status1 = gsl_odeiv2_driver_apply (d, &a, A, D);

  	/*Allocation of array values */
		//D1/D1'(k)
		D1  = D[0];
		dD1 = D[1];

  	gsl_odeiv2_driver_free(d);

  	return 0;
  }





	///// 2ND ORDER ////

	struct param_type4 {
    double kk;
		double kk1;
		double kk2;
		double kk1dotkk2;
  	double OmegaM00;
  	double par1;
  	double par2;
  	double par3;
		int mmgg;

  };

	int ode2(double a, const double D[], double dDda[], void *params)
  {
  	param_type4 p = *(param_type4 *)(params);
		double k      = p.kk;
		double k1     = p.kk1;
  	double k2     = p.kk2;
  	double k1dotk2  = p.kk1dotkk2;
  	double OmegaM0 = p.OmegaM00;
  	double p1     = p.par1;
  	double p2     = p.par2;
  	double p3     = p.par3;
		int mg = p.mmgg;


		//Use 2d interpoaltion to find D1 for any (a, k) combination
		double D1k1 = calc_D1_2d(a, k1, aa_arr, kk_arr, D1D1_arr);
		double D1k2 = calc_D1_2d(a, k2, aa_arr, kk_arr, D1D1_arr);


		//Equations for 2nd order kernels: k.L2(k, k1, theta)
		dDda[0] = D[1];
		dDda[1] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3, mg)*D[0] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3, mg)*D1k1*D1k2*(1-pow2(k1dotk2)) + ( 2./(pow4(a)*pow2(HA(a, OmegaM0))) )*gamma2(a, OmegaM0, k, k1, k2, 1., p1, p2, p3, mg); //u1 not needed for f(R), need to FIX for DGP

  	return GSL_SUCCESS;
  }


  int solve2(double A, double k, double OmegaM0, double par1, double par2, double par3, int mg)
  {
  //#pragma omp parallel for schedule(dynamic)
		for(int k1_num = 0; k1_num < Nk; k1_num++){

			//double k1 = kmin + k1_num*(kmax-kmin)/Nk; //Linear sampling
			double k1 = kmin * exp(k1_num*log(kmax/kmin)/(Nk*1.-1.)); //Exponential sampling

			for(int x_num = 0; x_num < Nx; x_num++){

				//x=cos(theta_k,k1)
				double x = 0;//xmin + x_num*(xmax-xmin)/Nx;

				double k2 = sqrt(pow2(k)+pow2(k1)-(2.*k*k1*x));

				double alpha = pipi - acos((pow2(k1) + pow2(k2) - pow2(k))/(2*k1*k2));

				double k1dotk2 = cos(alpha);


				// Initial scale factor for solving system of equations
				double a = a_ini;

				// Einstein de Sitter initial condtions for 2nd order
				double D[2] = {(3./7.)*(1-pow2(k1dotk2)), (6./7.)*(1-pow2(k1dotk2))/a}; //{(3./7.)*(1-pow2(k1dotk2))*pow2(a), (6./7.)*(1-pow2(k1dotk2))*a}; //{ (3./7.)*(1.-pow2(k1dotk2)), (6./7.)*(1.-pow2(k1dotk2))/a };

  			/*Parameters passed to system of equations */
  			// EDIT : If more than one gravity parameter is used, add them after p1
  			struct param_type4 my_params2 = {k, k1, k2, k1dotk2, OmegaM0, par1, par2, par3, mg};


  			gsl_odeiv2_system sys = {ode2, jac, 2, &my_params2};

  			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
  			gsl_odeiv2_driver * d =
  			gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
  										  1e-6, 1e-6, 0.0);

  			int status1 = gsl_odeiv2_driver_apply (d, &a, A, D);

  			/*Allocation of array values */

				//(k.L2)/d(k.L2)
				kL2[k1_num*Nx + x_num]  = D[0];
				dkL2[k1_num*Nx + x_num] = D[1];

  			gsl_odeiv2_driver_free(d);
  		}
  	}

  	return 0;
  }





/// OUTPUT SECTION ////

int main(int argc, char* argv[]) {

	//output file name
	//Could create if statements here so that different MG gives different file names
	const char* output = "kernelb3.dat";

	/* Open output file */
	FILE* fp = fopen(output, "w");

	//Declare variables to store values
	double p1, p2, p3, p4, p5;

	// 1ST ORDER //

	for(int a_num = 0 ; a_num <Na; a_num++){

		// Example a value
		double a_val = a_ini + a_num*(a_fin-a_ini)/Na; //Linear sampling
		a_arr[a_num] = a_val;

		for(int k_num = 0 ; k_num < Nk;  k_num ++){

			// Example k value
			//double k = kmin + k_num*(kmax-kmin)/Nk; //Linear sampling
			double k = kmin * exp(k_num*log(kmax/kmin)/(Nk*1.-1.)); //Exponential sampling
			k_arr[k_num] = k;

			// Format : solve1(scale factor, k, omega_(total matter), theory parameter1, theory parameter2 ...)
			solve1(a_val, k, omega0, 0.00001, 1., 1., MG);
			D1_arr[a_num*Nk + k_num] = D1;

			fprintf(fp,"%e %e %e \n", a_val, k, D1);

			//std::cout << "a, k, D1: " << a_val << " " << k << " " << D1 << '\n';

			}

	}

	aa_arr.setcontent(sizeof(a_arr)/sizeof(*a_arr), a_arr);
	kk_arr.setcontent(sizeof(k_arr)/sizeof(*k_arr), k_arr);
	D1D1_arr.setcontent(sizeof(D1_arr)/sizeof(*D1_arr), D1_arr);

	std::cout << "1ST ORDER COMPLETED" << '\n';

	// 2ND ORDER //

	for(int a_num = 0 ; a_num <Na; a_num++){

		//Example a value
		double a_val = a_ini + a_num*(a_fin-a_ini)/Na; //Linear sampling
		std::cout << "a:" << a_num << "  " << a_val << '\n';

		for(int k_num = 0 ; k_num < Nk;  k_num ++){

				// Example k value
				//double k = kmin + k_num*(kmax-kmin)/Nk; //Linear sampling
				double k = kmin * exp(k_num*log(kmax/kmin)/(Nk*1.-1.)); //Exponential sampling


				solve2( a_val, k, omega0, 0.00001, 1., 1., MG);

				for(int k1_num = 0; k1_num < Nk; k1_num++){

					//Example k1 value
					//double k1 = kmin + k1_num*(kmax-kmin)/Nk;  //Linear sampling
					double k1 = kmin * exp(k1_num*log(kmax/kmin)/(Nk*1.-1.)); //Exponential sampling

					for(int x_num = 0; x_num < Nx; x_num++){

						//Example x=cos(theta_k,k1) value
						double x = 0.5;//xmin + x_num*(xmax-xmin)/Nx;  //Linear sampling

						double k2 = sqrt(pow2(k)+pow2(k1)-2.*k*k1*x);

						double alpha = pipi - acos((pow2(k1) + pow2(k2) - pow2(k))/(2*k1*k2));

						double k1dotk2 = cos(alpha);


						//Store 1st and 2nd order values
						p1 = D1_arr[a_num*Nk + k_num];
						p2 = kL2[k1_num*Nx+x_num];
						p3 = kL2[k1_num*Nx+x_num]/D1_arr[a_num*Nk + k_num];
						p4 = -1.*kL2[k1_num*Nx+x_num]/(1.-pow2(k1dotk2)); //D2 = -(k.L2)/(1 - k1dotk2^2)
						p5 = -7.*p4/(3.*pow2(D1_arr[a_num*Nk + k_num])); //Expect D2 = -3/7 D1^2 for EdS

						/*
						if( fabs(p5-1.) > 1e-10 ){
							std::cout << "Warning bad D2/D1:" << fabs(p5-1.) << '\n';
						}*/

						//std::cout << "a:" << a_val <<  " k:" << k << " k1:" << k1 << " x:" << x << " k1dotk2:" << k1dotk2 << " D1:" <<  p1 << " kL2:" << p2 << " D2/D1:" << p5 << '\n';
						std::cout << "a:" << a_val <<  " k:" << k << " k1:" << k1 << " k2:" << k2 << " D1:" <<  p1 << " D2:" << p4 << " D2/D1:" << p5 << '\n';

					  // a, k, k1, k2, D2
					  fprintf(fp,"%e %e %e %e %e \n", a_val, k, k1, k2, p4);

			 	  	}
			  	}
			}
	}

	std::cout << "2ND ORDER COMPLETED" << '\n';

	fclose(fp);


	return 0;
}
