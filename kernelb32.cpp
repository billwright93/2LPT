// g++ -I/opt/apps/libs/gsl/2.1/gcc-4.4.7/include I/users/wrightw/work/cpplibs -L/opt/apps/libs/gsl/2.1/gcc-4.4.7/lib -lgsl -lstd++ -lgslcblas kernelb32.cpp
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




const double pipi = 3.1415926535897;
const	double H0 = 70.;///(3.09e+19);
double omega0 = 0.99;//0.99;//

//F2 calc counter
int count = 0;

// Number of time steps you wish to output between a_ini and a_fin
const int Na = 10;
double a_ini = 0.01;
double a_fin = 1.1;

// Number of k-modes you wish to output between kmin and kmax
const int Nk = 4; //10
double kmin = 0.01;//0.01;
double kmax = 0.2;//0.2;

// Number of x=cos(theta) modes you wish to output between xmin and xmax
const int Nx = 4; //10
double xmin = -1.+1.e-10;
double xmax = 1.-1.e-10;



//D1 array to store values for 2nd order solver
//std::vector<double> k_arr(Nk);
//std::vector<double> a_arr(Na*Nk);
//std::vector<double> D1_arr(Na*Nk);
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


////////  NUMERICAL KERNEL CONSTRUCTION //////////

//Easy power definitions
static double pow2(double base){
	return pow(base, 2.);}

static double pow3(double base){
	return pow(base, 3.);}

static double pow4(double base){
	return pow(base, 4.);}

/*
double round_to_digits(double value, int digits){
	if(value == 0.0){ // otherwise it will return 'nan' due to the log10() of zero
	  return 0.0;
	}
	double factor = pow(10.0, digits - ceil(log10(fabs(value))));
	return round(value * factor) / factor;
}*/

/*
//D1 spline interpolation function
double calc_D1(double k, std::vector<double>& k_arr, std::vector<double>& D1_arr){
	tk::spline s;
	s.set_points(k_arr, D1_arr);
	return s(k);
}*/

double calc_D1_2d(double a, double k, alglib::real_1d_array& a_arr, alglib::real_1d_array& k_arr, alglib::real_1d_array& D1_arr, alglib::spline2dinterpolant s){ //double a_arr[], double k_arr[], double D1_arr[]){
	// k=x, a=y
	double v;
	//double length_k_arr = sizeof(k_arr)/sizeof(*k_arr);
	//double length_a_arr = sizeof(a_arr)/sizeof(*a_arr);
	//alglib::spline2dinterpolant s;
	alglib::spline2dbuildbicubicv(k_arr, k_arr.length(), a_arr, a_arr.length(), D1_arr, 1, s);//length_k_arr, a_arr, length_a_arr, D1_arr, 1, s); //1=dimension of D1 as scalar field
	v = alglib::spline2dcalc(s, k, a);
	return v;
}

/*double calc_D1_2d(double a, double k, alglib::real_1d_array& a_arr, alglib::real_1d_array& k_arr, alglib::real_1d_array& D1_arr){ //double a_arr[], double k_arr[], double D1_arr[]){
	// k=x, a=y
	using namespace alglib;
	double v;
	//double length_k_arr = sizeof(k_arr)/sizeof(*k_arr);
	//double length_a_arr = sizeof(a_arr)/sizeof(*a_arr);
  spline2dbuildbicubicv(k_arr, k_arr.length(), a_arr, a_arr.length(), D1_arr, 1, s);//length_k_arr, a_arr, length_a_arr, D1_arr, 1, s); //1=dimension of D1 as scalar field
	v = spline2dcalc(s, k, a);
	return v;
}*/


// Normalized Hubble parameter H and a*H*dH/da for the background (LCDM)

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





                                          ////
                                      /////////////
                              //////////////////////////////
                        /////////////////////////////////////////
                  ///////////////////////////////////////////////////////
              /////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////

  //MODIFIED GRAVITY MODEL FUNCTIONS:
  //p1,p2,p3 are theory parameters
  //k1,k2,k3 are the mode magnitudes
  //u1 angle between k2 and k3; u2 between k1 and k2; u3 between k1 and k3


  // default is 3 model params. Add in parameters p4,p5,... as required.
  // Note on theories with more than 3 parameters:
  // 1) Add in extra parameters to param_type3 structure below
  // 2) Add in extra parameters to params of ode1 system of ODEs
  // 3) Add in extra parameters to initn function and in params passed by function to solver (this is right after initial conditions are passed)

  static double mu(double a, double k0, double OmegaM0, double p1, double p2, double p3 ){
  	double h0 = 1./2997.9;
  //	return 1.; // GR
  	return 1. + pow2(k0/a)/(3.*(pow2(k0/a)+pow3(OmegaM0/pow3(a)-4.*(OmegaM0-1.))/(2.*p1/pow2(h0)*pow2(4-3.*OmegaM0)))); //f(R) Hu- Sawicki
  //	return 1.+1./(3.*beta(a,OmegaM0,p1)); //nDGP
  }


  static double gamma2(double a, double OmegaM0, double k0, double k1, double k2, double u1, double p1, double p2, double p3 ){
  	double h0 = 1./2997.9;
  	 return 0. ; // GR

  	 return -(9.*pow2(k0/a)*pow2(OmegaM0/pow3(a))*pow(OmegaM0-4.*pow3(a)*(-1+OmegaM0),5))/
  			    (48.*pow(a,15)*pow2(p1/pow2(h0))*pow2(HA(a,OmegaM0))*pow4(3.*OmegaM0-4.)
  			   *(pow2(k0/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
  			   *(pow2(k1/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
  			   *(pow2(k2/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))); //f(R) Hu- Sawicki


  //  	return -1.*p2/(HA(a,OmegaM0)*HA(a,OmegaM0)*24.*pow(beta(a,OmegaM0,p1),3)*p1)*pow(OmegaM0/(a*a*a),2)*ker1(u1); //nDGP
  }

  ///////// NON-SEPARABLE //////////////

	//Jacobian of the system required when calling the system evolver, the below is neither needed, nor correct for solving
	int jac (double a, const double D[], double *dfdy, double dfdt[], void *params)
	{
		return GSL_SUCCESS;
	}

	//1st order//

  /* Parameters passed to system of Euler and continuity equations*/
  // k (magnitude) and x (angular) values for the system of equations
  // k1.k2=k1k2x2 , k1.k3 = k1k3x3, k2.k3=k2k3x1
  // EDIT : Add in new gravity parameters to this list as required (par1 is the nDGP parameter omega_rc as default)
  struct param_type3 {
    double kk;
  	double OmegaM00;
  	double par1;
  	double par2;
  	double par3;
  };

  int ode1(double a, const double D[], double dDda[], void *params)
  {
  	param_type3 p = *(param_type3 *)(params);
  	double k = p.kk;
  	double OmegaM0 = p.OmegaM00;
  	double p1 = p.par1;
  	double p2 = p.par2;
  	double p3 = p.par3;

		double OmegaL0 = 1. - OmegaM0;

		//1st order: D1/dD1(k1)
		dDda[0] = D[1];
		dDda[1] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3)*D[0]; //+ 1.5*(OmegaM0/(pow2(a)*(OmegaM0 + OmegaL0*pow3(a))))*mu(a, k, OmegaM0, p1, p2, p3)*D[0];

  	return GSL_SUCCESS;
  }


  //Loops to construct numerical kernels
  // Kmax and Kmin indicates the maximum and minimum values of K
  // k loop dictates the angular parameter x =k1.k/(k^2*r)
  // j loop dictates the magnitude of k1 (k1=rk)

  // Default initialization at only scale factor A
  // k is your k vector dependence (See mu(a,k))
  // YMAX = QMAX/kmin and YMIN=QMIN/kmax where kmin and kmax are the k-limits and QMIN and QMAX are k1 limits (integrated over)
  // QMAX => 5, QMIN ~ 0.0005  FOR PS INTEGRALS
  // QMAX => 20, QMIN ~ 0.001 FOR XI INTEGRALS
  // set k = 1 and YMIN = QMIN/kmax, YMAX = QMAX/kmin for no k -dependence (initialize once)
  // set kmin = kmax = k for k dependence (initialize for each k)

  // EDIT : Add in new gravity parameters to function's input as required (par1 is the nDGP parameter omega_rc as default)
  // EDIT : if more than 1 gravity parameter is needed, one must also edit initn function in SPT.h file.

  int solve1(double A, double k, double OmegaM0, double par1, double par2, double par3)
  {
  //#pragma omp parallel for schedule(dynamic)

			// Initial scale factor for solving system of equations
			//Sets 1-3
				double a = a_ini; //0.0001;

				// Einstein de Sitter initial condtions for 1st order growth factor and derivative
				double D[2] = {a, 1.};//{ 1., 1./a };//{ 1., 1./a };// //{a, a*HA(a, OmegaM0)};//, -(3./7.)*a, -(3./7.)*a*HA(a, OmegaM0)}; // initial conditions

  		/*Parameters passed to system of equations */
  // EDIT : If more than one gravity parameter is used, add them after p1
  		struct param_type3 my_params1 = {k, OmegaM0, par1, par2, par3};

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


	//2nd order

	struct param_type4 {
    double kk;
		double kk1;
		double kk2;
		double kk1dotkk2;
  	double OmegaM00;
  	double par1;
  	double par2;
  	double par3;
		int aa_num;
		int kk1_num;
		alglib::spline2dinterpolant ss;
		//int kk2_num;

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
		int a_num     = p.aa_num;
		int k1_num 		= p.kk1_num;
		alglib::spline2dinterpolant s = p.ss;
		//int k2_num		= p.kk2_num;


		double OmegaL0 = 1. - OmegaM0;

		/*std::vector<double> D1a_arr(D1_arr.begin() + (a_num)*Nk , D1_arr.begin() + (a_num+1)*Nk  );

		int k_size = k_arr.size();//sizeof(k_arr)/sizeof(*k_arr);
		int D1_size = D1a_arr.size();//sizeof(D1a_arr)/sizeof(*D1a_arr);

		std::cout << "a_num: "  << a_num << " a: " << a << '\n';

		for(int i=0; i<Nk; i++){
			double test_k = kmin + i*(kmax-kmin)/Nk;
			std::cout << "D1_old[" << a << ", "  << test_k << "]= " << D1_arr[(a_num*Nk)+i] << '\n';
		}

		for(int j=0; j<D1_size; j++){
			std::cout << "D1_vec[" << a << ", "  << k_arr[j] << "]= " << D1a_arr[j] << '\n';
		}*/

		double D1k1 = calc_D1_2d(a, k1, aa_arr, kk_arr, D1D1_arr, s);
		double D1k2 = calc_D1_2d(a, k2, aa_arr, kk_arr, D1D1_arr, s);

		/*int sw = 0;
		if(a_num==2){
			//std::cout << a_num << "  " << a << '\n';
			for(int i=0; i<k_size; i++){
				//std::cout << "k: " << k_arr[i] << '\n';
				if(round_to_digits(k2, 4)==round_to_digits(k_arr[i], 4)){
					//std::cout << k2 << "k2==k_arr[i]" << k_arr[i] << '\n';
					sw++;
				}
			}
			for(int j=0; j<D1_size; j++){
				//std::cout << "D1: " << D1a_arr[j] << '\n';
			}
			if(sw == 0){
				//std::cout << a_num << "  " << a << '\n';
				for(int i=0; i<k_size; i++){
					//std::cout << "k: " << k_arr[i] << '\n';
				}
				for(int j=0; j<D1_size; j++){
					//std::cout << "D1: " << D1a_arr[j] << '\n';
				}
				//std::cout << "D1(" << k2 << "): " << D1k2 << '\n';
			}

		}

		//std::cout << D1_arr[a_num*Nk + k1_num] << "  " << D1_arr[a_num*Nk + k2_num] << '\n';
		//std::cout << "D1(k1):" << D1k1 << "  D1(k2):" << D1k2 << '\n';*/


		/*double D11_old = round_to_digits(D1_arr[a_num*Nk + k1_num], 4);//round(D1_arr[a_num*Nk + k1_num]);//
		//double D12_old = round_to_digits(D1_arr[a_num*Nk + k2_num], 10);//round(D1_arr[a_num*Nk + k2_num]);//
		double D11_new = round_to_digits(D1k1, 4);//round(D1k1);//
		double D12_new = round_to_digits(D1k2, 4);//round(D1k2);//

		//std::cout.precision(99);

		if(D11_old != D11_new){
			std::cout << D11_old << " old1 != new1 " << D11_new << '\n';
		}*/

		/*if(abs(D12_old) != abs(D12_new)){
			std::cout << D12_old << " old2 != new2 " << D12_new << '\n';
		}*/


		//2nd order kernels: k.L2(k, k1, theta)
		dDda[0] = D[1];
		//dDda[1] = - ( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + (3./2.)*(OmegaM0/(pow2(a)*(OmegaM0 + OmegaL0*pow3(a))))*mu(a, k, OmegaM0, p1, p2, p3)*D[0]; + (3./2.)*(OmegaM0/(pow2(a)*(OmegaM0 + OmegaL0*pow3(a))))*mu(a, k, OmegaM0, p1, p2, p3)*D1_arr[a_num*Nk + k1_num]*D1_arr[a_num*Nk + k2_num]*(1-pow2(k1dotk2)) + ( 2./(pow4(a)*pow2(HA(a, OmegaM0))) )*gamma2(a, OmegaM0, k, k1, k2, 1., p1, p2, p3); //u1 not needed for f(R), need to FIX for DGP
		dDda[1] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3)*D[0] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3)*D1k1*D1k2*(1-pow2(k1dotk2)) + ( 2./(pow4(a)*pow2(HA(a, OmegaM0))) )*gamma2(a, OmegaM0, k, k1, k2, 1., p1, p2, p3); //u1 not needed for f(R), need to FIX for DGP

  	return GSL_SUCCESS;
  }


  int solve2(int A_num, double A, double k, double OmegaM0, double par1, double par2, double par3, alglib::spline2dinterpolant s)//solve2(int A_num, double A, double k, double OmegaM0, double par1, double par2, double par3)
  {
  //#pragma omp parallel for schedule(dynamic)
		for(int k1_num = 0; k1_num < Nk; k1_num++){

			double k1 = kmin + k1_num*(kmax-kmin)/Nk;

			for(int x_num = 0; x_num < Nx; x_num++){

				//x=cos(theta_k,k1)
				double x = xmin + x_num*(xmax-xmin)/Nx;

				double k2 = sqrt(pow2(k)+pow2(k1)-(2.*k*k1*x));

				double alpha = pipi - acos((pow2(k1) + pow2(k2) - pow2(k))/(2*k1*k2));

				double k1dotk2 = cos(alpha);//(pi - acos( (pow2(k1)+pow2(k2)-pow2(k))/(2.*k1*k2) ) );

				//std::cout << "D1["<< A << ", " << k1 << "]=" << D1_arr[A_num*Nk + k1_num] << '\n';
				//std::cout << "D1["<< A << ", " << k2 << "]=" << D1_arr[A_num*Nk + k2_num] << '\n';

				// Initial scale factor for solving system of equations
				double a = a_ini;

				// Einstein de Sitter initial condtions for 2nd order
				double D[2] = {(3./7.)*(1-pow2(k1dotk2))*pow2(a), (6./7.)*(1-pow2(k1dotk2))*a};//{ (3./7.)*(1.-pow2(k1dotk2)), (6./7.)*(1.-pow2(k1dotk2))/a }; // initial conditions

  			/*Parameters passed to system of equations */
  			// EDIT : If more than one gravity parameter is used, add them after p1
  			struct param_type4 my_params2 = {k, k1, k2, k1dotk2, OmegaM0, par1, par2, par3, A_num, k1_num, s};//, k2_num};//


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

				//std::cout << "Done: " << k << "/" << k1<< "/" << k2 << '\n';

  			gsl_odeiv2_driver_free(d);
  		}
  	}

  	return 0;
  }

/// OUTPUT SECTION ////

int main(int argc, char* argv[]) {

	//output file name
	const char* output = "kernelb32.dat";

	/* Open output file */
	FILE* fp = fopen(output, "w");

	alglib::spline2dinterpolant s;

	double p1, p2, p3, p4, p5;

	for(int a_num = 0 ; a_num <Na; a_num++){

		// Example a value
		double a_val = a_ini + a_num*(a_fin-a_ini)/Na; //Linear sampling
		std::cout << "a:" << a_num << "  " << a_val << '\n';

		a_arr[a_num] = a_val;

		for(int k_num = 0 ; k_num < Nk;  k_num ++){

			// Example k value
			double k = kmin + k_num*(kmax-kmin)/Nk;
			k_arr[k_num] = k;

			//std::cout << "k:" << k << "\n";
			//std::cout << "mu(k):" << mu(0.1, k, 0.24, 0.0001, 1., 1.) << "\n";

			// Format : solve1(scale factor, ymin,ymax, k., omega_(total matter), theory parameter1, theory parameter2 ...)
			// Parameterized magnitude integration q = ky
			// set k = 1. for scale independant and take the function out of the loop
			solve1(a_val, k, omega0, 0.0001, 1., 1.);//solve1(1., ymax, ymin , k, 0.24, 0.0001, 1., 1.);
			D1_arr[a_num*Nk + k_num] = D1;

			if(a_num == (Na-1)){
				fprintf(fp,"%e %e \n", k, D1);}

			//std::cout << "a, k, D1: " << a_val << " " << k << " " << D1 << '\n';

			}

		aa_arr.setcontent(sizeof(a_arr)/sizeof(*a_arr), a_arr);
		kk_arr.setcontent(sizeof(k_arr)/sizeof(*k_arr), k_arr);
		D1D1_arr.setcontent(sizeof(D1_arr)/sizeof(*D1_arr), D1_arr);

		for(int k_num = 0 ; k_num < Nk;  k_num ++){

				// Example k value
				double k = kmin + k_num*(kmax-kmin)/Nk;
				//std::cout << "k:" << k << "\n";
				//std::cout << "mu(k):" << mu(0.1, k, 0.24, 0.0001, 1., 1.) << "\n";

				//std::cout << "D1(" << a_val << ")=" << D1_arr[a_num*Nk + k_num] << '\n';

				solve2(a_num, a_val, k, omega0, 0.0001, 1., 1., s);

				for(int k1_num = 0; k1_num < Nk; k1_num++){

					double k1 = kmin + k1_num*(kmax-kmin)/Nk;

					for(int x_num = 0; x_num < Nx; x_num++){

						//x=cos(theta_k,k1)
						double x = xmin + x_num*(xmax-xmin)/Nx;

						double k2 = sqrt(pow2(k)+pow2(k1)-2.*k*k1*x);

						double alpha = pipi - acos((pow2(k1) + pow2(k2) - pow2(k))/(2*k1*k2));

						double k1dotk2 = cos(alpha);//(pi - acos( (pow2(k1)+pow2(k2)-pow2(k))/(2.*k1*k2) ) );

						//solve2(a_num, a_val, k, 0.24, 0.0001, 1., 1.); //solve2(a_num, a_val, k, 0.24, 0.0001, 1., 1.)

						//F2/G2(p,k)
					  //p3 = F2A_nk[j*n2+m];
						p1 = D1_arr[a_num*Nk + k_num];
						p2 = kL2[k1_num*Nx+x_num];
						p3 = kL2[k1_num*Nx+x_num]/D1_arr[a_num*Nk + k_num];
						p4 = -1.*kL2[k1_num*Nx+x_num]/(1.-pow2(k1dotk2)); //D2 = -(k.L2)/(1 - k1dotk2^2)
						p5 = -7.*p4/(3.*pow2(D1_arr[a_num*Nk + k_num])); //Expect D2 = -3/7 D1^2 for EdS

						std::cout << "a:" << a_val <<  " k:" << k << " k1:" << k1 << " x:" << x << " k1dotk2:" << k1dotk2 << " D1:" <<  p1 << " kL2:" << p2 << " D2/D1:" << p5 << '\n';
						//std::cout << "a, k, k1, x, k1dotk2, D1, kL2, ratio: " << a_val << " " << k << " " << k1 << " " << x << " " << p1 << " " << p2 << " " << p3 << '\n';
					  //std::cout << "a, k, k1, k2, k1dotk2, D1, kL2: " << a_val << " " << k << " " << k1 << " " << k2 << " " << k1dotk2 << " " << p1 << " " << p2 << '\n';

					  // a, k, k1, x, D1(a, k), kL2(a, k, k1, x)
					  //fprintf(fp,"%e %e %e %e %e %e \n", a_val, k, k1, x, p1, p2);
			 	  	}
			  	}
			}
	}

	/*close output file*/
	fclose(fp);
	return 0;
}
