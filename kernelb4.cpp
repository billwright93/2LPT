//COMPILE + RUN COMMANDS//

// g++ -I/opt/apps/libs/gsl/2.1/gcc-4.4.7/include -L/opt/apps/libs/gsl/2.1/gcc-4.4.7/lib -lgsl -lstdc++ -lgslcblas kernelb4.cpp
// ./a.out



#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <stdlib.h>
#include <iostream>



// CONSTANTS //
const double pipi = 3.1415926535897;
//const	double H0 = (3.09e+19);//70.;//
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
const int Nx = 4; //10
double xmin = -1.+1.e-10;
double xmax = 1.-1.e-10;



// DECLARING VARIABLES AND ARRAYS //

//D1/D1_dot(k)
double D1;
double dD1;

//D1(k)
double D1k[Nk*Nx];

//D1(k1)
double D1k1[Nk*Nx];

//D1(k2)
double D1k2[Nk*Nx];

//k.L2(k, k1, theta)
double kL2[Nk*Nx];
//double dkL2[Nk*Nx];



// FUNCTIONS //

//Easy power definitions
static double pow2(double base){
	return pow(base, 2.);}

static double pow3(double base){
	return pow(base, 3.);}

static double pow4(double base){
	return pow(base, 4.);}


// Normalized Hubble parameter H, a*H*dH/da, and Omega_M(a) for the background (LCDM)

//H(a)
static double HA(double a, double OmegaM0){
	double OmegaL0 = 1.- OmegaM0;
	//return  H0*sqrt((OmegaM0/pow3(a))+OmegaL0);}
	return  sqrt((OmegaM0/pow3(a))+OmegaL0);}

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
	return 1.+(HA(a, OmegaM0))/sqrt(omegarc)*(1.+(-1.5*OmegaM0/pow(a,3))/(3.*pow2(HA(a,OmegaM0))));
}

//mu
static double mu(double a, double k, double OmegaM0, double p1, double p2, double p3, int mg){
	double h0 = 1./2997.9;
	if(mg == 0){
		return 1.; // GR
	}
	else if(mg == 1){
		return 1. + pow2(k/a)/(3.*(pow2(k/a)+pow3(OmegaM0/pow3(a)-4.*(OmegaM0-1.))/(2.*p1/pow2(h0)*pow2(4-3.*OmegaM0)))); //f(R) Hu- Sawicki
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
static double gamma2(double a, double OmegaM0, double k, double k1, double k2, double p1, double p2, double p3, int mg){
  double h0 = 1./2997.9;
	if(mg == 0){
		return 0. ; // GR
	}
	else if(mg == 1){
		return -(9.*pow2(k/a)*pow2(OmegaM0/pow3(a))*pow(OmegaM0-4.*pow3(a)*(-1+OmegaM0),5))/
	  			    (48.*pow(a,15)*pow2(p1/pow2(h0))*pow2(HA(a,OmegaM0))*pow4(3.*OmegaM0-4.)
	  			   *(pow2(k/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
	  			   *(pow2(k1/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
	  			   *(pow2(k2/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))); //f(R) Hu- Sawicki
	}
	else if(mg == 2){
		double u1 = pow2(k1)+pow2(k2)-pow2(k)/(2.*k1*k2);
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



	///// 1ST + 2ND ORDER ////

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


		//Equations for 1st order growth factor D1/dD1(k)
		dDda[0] = D[1];
		dDda[1] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3, mg)*D[0];

		//Equations for 1st order growth factor D1/dD1(k1)
		dDda[2] = D[3];
		dDda[3] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[3] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k1, OmegaM0, p1, p2, p3, mg)*D[2];

		//Equations for 1st order growth factor D1/dD1(k2)
		dDda[4] = D[5];
		dDda[5] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[5] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k2, OmegaM0, p1, p2, p3, mg)*D[4];

		//Equations for 2nd order kernels: k.L2(k, k1, theta)
		dDda[6] = D[7];
		dDda[7] = -1.*( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[7] + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3, mg)*D[6]
		         + 1.5*(Omega_M(a, OmegaM0)/pow2(a))*mu(a, k, OmegaM0, p1, p2, p3, mg)*D[2]*D[4]*(1-pow2(k1dotk2)) + ( 2.*D[2]*D[4]/(pow4(a)*pow2(HA(a, OmegaM0))) )*gamma2(a, OmegaM0, k, k1, k2, p1, p2, p3, mg);

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
				double x = xmin + x_num*(xmax-xmin)/Nx; //0.5;

				double k2 = sqrt(pow2(k)+pow2(k1)-(2.*k*k1*x));

				double alpha = pipi - acos((pow2(k1) + pow2(k2) - pow2(k))/(2*k1*k2));

				double k1dotk2 = cos(alpha);


				// Initial scale factor for solving system of equations
				double a = a_ini;

				// Einstein de Sitter initial condtions for 2nd order
				double D[8] = { 1., 1./a, 1., 1./a, 1., 1./a, (3./7.)*(1-pow2(k1dotk2)), (6./7.)*(1-pow2(k1dotk2))/a };//  1., 1./a}; //{(3./7.)*(1-pow2(k1dotk2))*pow2(a), (6./7.)*(1-pow2(k1dotk2))*a}; //{ (3./7.)*(1.-pow2(k1dotk2)), (6./7.)*(1.-pow2(k1dotk2))/a };

  			/*Parameters passed to system of equations */
  			// EDIT : If more than one gravity parameter is used, add them after p1
  			struct param_type4 my_params2 = {k, k1, k2, k1dotk2, OmegaM0, par1, par2, par3, mg};


  			gsl_odeiv2_system sys = {ode2, jac, 8, &my_params2};

  			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
  			gsl_odeiv2_driver * d =
  			gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
  										  1e-6, 1e-6, 0.0);

  			int status1 = gsl_odeiv2_driver_apply (d, &a, A, D);

  			/*Allocation of array values */

				//D1(k)
				D1k[k1_num*Nx + x_num] = D[0];

				//D1(k1)
				D1k1[k1_num*Nx + x_num] = D[2];

				//D1(k2)
				D1k2[k1_num*Nx + x_num] = D[4];

				//(k.L2)/d(k.L2)
				kL2[k1_num*Nx + x_num]  = D[6];
				//dkL2[k1_num*Nx + x_num] = D[7];

  			gsl_odeiv2_driver_free(d);
  		}
  	}

  	return 0;
  }





/// OUTPUT SECTION ////

int main(int argc, char* argv[]) {

	//output file name
	//Could create if statements here so that different MG gives different file names
	const char* output = "kernelb4.dat";

	/* Open output file */
	FILE* fp = fopen(output, "w");

	//Declare variables to store values
	double p1, p2, p3, p4, p5;

	// 1ST ORDER //
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
						double x = xmin + x_num*(xmax-xmin)/Nx;  //Linear sampling //0.5;

						double k2 = sqrt(pow2(k)+pow2(k1)-2.*k*k1*x);

						double alpha = pipi - acos((pow2(k1) + pow2(k2) - pow2(k))/(2*k1*k2));

						double k1dotk2 = cos(alpha);


						//Store 1st and 2nd order values
						p1 = D1k[k1_num*Nx+x_num]; //D1(a, k)
						p2 = kL2[k1_num*Nx+x_num]; //[k.L2](a, k, k1, x)
						p3 = p2/p1; //[k.L2](a, k, k1, x)/D1(a, k)
						p4 = -1.*p2/(1.-pow2(k1dotk2)); //D2 = -(k.L2)/(1 - k1dotk2^2)
						p5 = -7.*p4/(3.*pow2(p1)); //Expect D2 = -3/7 D1^2 for EdS

						/*
						if( fabs(p5-1.) > 1e-10 ){
							std::cout << "Warning bad D2/D1:" << fabs(p5-1.) << '\n';
						}*/

						//std::cout << "a:" << a_val <<  " k:" << k << " k1:" << k1 << " x:" << x << " k1dotk2:" << k1dotk2 << " D1:" <<  p1 << " kL2:" << p2 << " D2/D1:" << p5 << '\n';
						std::cout << "a:" << a_val <<  " k:" << k << " k1:" << k1 << " k2:" << k2 << " x:" << x << " D1:" <<  p1 << " D2:" << p4 << " D2/D1:" << p5 << '\n';

					  // a, k, k1, k2, D1, D2
					  fprintf(fp,"%e %e %e %e %e %e \n", a_val, k, k1, k2, p1, p4);

			 	  	}
			  	}
			}
	}

	//std::cout << "2ND ORDER COMPLETED" << '\n';

	fclose(fp);


	return 0;
}
