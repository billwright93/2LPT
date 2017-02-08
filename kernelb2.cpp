#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <stdlib.h>

#include <iostream>

const double pi = 3.1415926535897;
const	double H0 = 70.;///(3.09e+19);

//F2 calc counter
int count = 0;

// Number of k-modes you wish to output between kmin and kmax
const int Nk = 2; //10
double kmin = 0.01;//0.01;
double kmax = 0.2;//0.2;

// Number of time steps you wish to output between a_ini and a_fin
const int Na = 1000;
double a_ini = 0.01;
double a_fin = 1.1;

//D1 array to store values for 2nd order solver

double D1_arr[Na*Nk];


//D1/D1_dot(k)
double D1;
double dD1;

//k.L2(k, k1, theta)
double kL2[Nk*Nk];
double dkL2[Nk*Nk];


////////  NUMERICAL KERNEL CONSTRUCTION //////////

//Easy power definitions
static double pow2(double base){
	return pow(base, 2.);}

static double pow3(double base){
	return pow(base, 3.);}

static double pow4(double base){
	return pow(base, 4.);}

// Normalized Hubble parameter H and a*H*dH/da for the background (LCDM)

//H(a)
static double HA(double a, double OmegaM0){
	double OmegaL0 = 1.-OmegaM0;
	return  H0*sqrt(OmegaM0/pow(a,3)+OmegaL0);}

// dH(a)/da
static double dHA(double a, double OmegaM0){
	double OmegaL0 = 1.-OmegaM0;
	return -(3./2.)*(OmegaM0/pow4(a))*(1./(OmegaM0/pow(a,3)+OmegaL0))*HA(a, OmegaM0);}

//
static double Omega_M(double a, double OmegaM0){
	return OmegaM0*pow2(H0)/(pow2(HA(a, OmegaM0))*pow3(a));}





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
  	return 1.; // GR
  //	return 1. + pow2(k0/a)/(3.*(pow2(k0/a)+pow3(OmegaM0/pow3(a)-4.*(OmegaM0-1.))/(2.*p1/pow2(h0)*pow2(4-3.*OmegaM0)))); //f(R) Hu- Sawicki
  //	return 1.+1./(3.*beta(a,OmegaM0,p1)); //nDGP
  }


  static double gamma2(double a, double OmegaM0, double k0, double k1, double k2, double u1, double p1, double p2, double p3 ){
  	double h0 = 1./2997.9;
  	 return 0. ; // GR

  /*	 return -(9.*pow2(k0/a)*pow2(OmegaM0/pow3(a))*pow(OmegaM0-4.*pow3(a)*(-1+OmegaM0),5))/
  			    (48.*pow(a,15)*pow2(p1/pow2(h0))*pow2(HA(a,OmegaM0))*pow4(3.*OmegaM0-4.)
  			   *(pow2(k0/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
  			   *(pow2(k1/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))
  			   *(pow2(k2/a)+pow3(OmegaM0-4.*pow3(a)*(OmegaM0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*OmegaM0-4.)))); //f(R) Hu- Sawicki
	*/

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
		dDda[1] = - ( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + (3./2.)*(OmegaM0/(pow2(a)*(OmegaM0 + OmegaL0*pow3(a))))*mu(a, k, OmegaM0, p1, p2, p3)*D[0];

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
				double D[2] = { 1., 1./a };//{ 1., 1./a };// //{a, a*HA(a, OmegaM0)};//, -(3./7.)*a, -(3./7.)*a*HA(a, OmegaM0)}; // initial conditions

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
		int kk1_num;
		int kk2_num;
		int aa_num;
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
		int k1_num 		= p.kk1_num;
		int k2_num		= p.kk2_num;
		int a_num     = p.aa_num;

		double OmegaL0 = 1. - OmegaM0;

		//2nd order kernels: k.L2(k, k1, theta)
		dDda[0] = D[1];
		dDda[1] = - ( (3./a) + (dHA(a, OmegaM0)/HA(a, OmegaM0)) )*D[1] + (3./2.)*(OmegaM0/(pow2(a)*(OmegaM0 + OmegaL0*pow3(a))))*mu(a, k, OmegaM0, p1, p2, p3)*D[0]; + (3./2.)*(OmegaM0/(pow2(a)*(OmegaM0 + OmegaL0*pow3(a))))*mu(a, k, OmegaM0, p1, p2, p3)*D1_arr[a_num*Nk + k1_num]*D1_arr[a_num*Nk + k2_num]*(1-pow2(k1dotk2)) + ( 2./(pow4(a)*pow2(HA(a, OmegaM0))) )*gamma2(a, OmegaM0, k, k1, k2, 1., p1, p2, p3); //u1 not needed for f(R), need to FIX for DGP

  	return GSL_SUCCESS;
  }

  //Loops to construct numerical kernels
  // Kmax and Kmin indicates the maximum and minimum values of K
  // k loop dictates the angular parameter x =k1.k/(k^2*r)
  // j loop dictates the magnitude of k1 (k1=rk)

  // Default initialization at only scale factor Bill
  // k is your k vector dependence (See mu(a,k))
  // YMAX = QMAX/kmin and YMIN=QMIN/kmax where kmin and kmax are the k-limits and QMIN and QMAX are k1 limits (integrated over)
  // QMAX => 5, QMIN ~ 0.0005  FOR PS INTEGRALS
  // QMAX => 20, QMIN ~ 0.001 FOR XI INTEGRALS
  // set k = 1 and YMIN = QMIN/kmax, YMAX = QMAX/kmin for no k -dependence (initialize once)
  // set kmin = kmax = k for k dependence (initialize for each k)

  // EDIT : Add in new gravity parameters to function's input as required (par1 is the nDGP parameter omega_rc as default)
  // EDIT : if more than 1 gravity parameter is needed, one must also edit initn function in SPT.h file.

  int solve2(int A_num, double A, double k, double OmegaM0, double par1, double par2, double par3)//solve2(int A_num, double A, double k, double OmegaM0, double par1, double par2, double par3)
  {
  //#pragma omp parallel for schedule(dynamic)
		for(int k1_num = 0; k1_num < Nk; k1_num++){

			for(int k2_num = 0; k2_num < Nk; k2_num++){

				double k1 = kmin + k1_num*(kmax-kmin)/Nk;
				double k2 = kmin + k2_num*(kmax-kmin)/Nk;

				//k1dotk2/k1k2
				double k1dotk2 = cos(pi - acos( (pow2(k1)+pow2(k2)-pow2(k))/(2.*k1*k2) ) ); //pig C++ language here - FIX ?

				//std::cout << "D1["<< A << ", " << k1 << "]=" << D1_arr[A_num*Nk + k1_num] << '\n';
				//std::cout << "D1["<< A << ", " << k2 << "]=" << D1_arr[A_num*Nk + k2_num] << '\n';

				// Initial scale factor for solving system of equations
				double a = a_ini;

				// Einstein de Sitter initial condtions for 2nd order
				double D[2] = { -(3./7.)*(1.-pow2(k1dotk2)), -(6./7.)*(1.-pow2(k1dotk2))/a }; // initial conditions

  			/*Parameters passed to system of equations */
  			// EDIT : If more than one gravity parameter is used, add them after p1
  			struct param_type4 my_params2 = {k, k1, k2, k1dotk2, OmegaM0, par1, par2, par3, k1_num, k2_num, A_num};//


  			gsl_odeiv2_system sys = {ode2, jac, 2, &my_params2};


  			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
  			gsl_odeiv2_driver * d =
  			gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
  										  1e-6, 1e-6, 0.0);

  			int status1 = gsl_odeiv2_driver_apply (d, &a, A, D);

  			/*Allocation of array values */

				//(k.L2)/d(k.L2)
				kL2[k1_num*Nk + k2_num]  = D[0];
				dkL2[k1_num*Nk + k2_num] = D[1];

				//std::cout << "Done: " << k << "/" << k1<< "/" << k2 << '\n';

  			gsl_odeiv2_driver_free(d);
  		}
  	}

  	return 0;
  }

/// OUTPUT SECTION ////

int main(int argc, char* argv[]) {

	//output file name
	const char* output = "kernelb2.dat";

	/* Open output file */
	FILE* fp = fopen(output, "w");

	double p1, p2, p3;

	for(int a_num = 0 ; a_num <Na; a_num++){

		// Example a value
		double a_val = a_ini + a_num*(a_fin-a_ini)/Na; //Linear sampling
		std::cout << "a:" << a_num << "  " << a_val << '\n';

		for(int k_num = 0 ; k_num < Nk;  k_num ++){

			// Example k value
			double k = kmin + k_num*(kmax-kmin)/Nk;
			//std::cout << "k:" << k << "\n";
			//std::cout << "mu(k):" << mu(0.1, k, 0.24, 0.0001, 1., 1.) << "\n";

			// Format : solve1(scale factor, ymin,ymax, k., omega_(total matter), theory parameter1, theory parameter2 ...)
			// Parameterized magnitude integration q = ky
			// set k = 1. for scale independant and take the function out of the loop
			solve1(a_val, k, 0.24, 0.0001, 1., 1.);//solve1(1., ymax, ymin , k, 0.24, 0.0001, 1., 1.);
			D1_arr[a_num*Nk + k_num] = D1;

			//std::cout << "a, k, D1: " << a_val << " " << k << " " << D1 << '\n';

			}

		for(int k_num = 0 ; k_num < Nk;  k_num ++){

				// Example k value
				double k = kmin + k_num*(kmax-kmin)/Nk;
				//std::cout << "k:" << k << "\n";
				//std::cout << "mu(k):" << mu(0.1, k, 0.24, 0.0001, 1., 1.) << "\n";

				//std::cout << "D1(" << a_val << ")=" << D1_arr[a_num*Nk + k_num] << '\n';

				solve2(a_num, a_val, k, 0.24, 0.0001, 1., 1.);

				for(int k1_num = 0; k1_num < Nk; k1_num++){

					double k1 = kmin + k1_num*(kmax-kmin)/Nk;

					for(int k2_num = 0; k2_num < Nk; k2_num++){

						double k2 = kmin + k2_num*(kmax-kmin)/Nk;

						double k1dotk2 = cos(pi - acos( (pow2(k1)+pow2(k2)-pow2(k))/(2.*k1*k2) ) );

						//solve2(a_num, a_val, k, 0.24, 0.0001, 1., 1.); //solve2(a_num, a_val, k, 0.24, 0.0001, 1., 1.)

						//F2/G2(p,k)
					  //p3 = F2A_nk[j*n2+m];
						p1 = D1_arr[a_num*Nk + k_num];
						p2 = kL2[k1_num*Nk+k2_num];
						p3 = kL2[k1_num*Nk+k2_num]/D1_arr[a_num*Nk + k_num];

						std::cout << "a, k, k1, k2, k1dotk2, D1, kL2, ratio: " << a_val << " " << k << " " << k1 << " " << k2 << " " << p1 << " " << p2 << " " << p3 << '\n';
					  //std::cout << "a, k, k1, k2, k1dotk2, D1, kL2: " << a_val << " " << k << " " << k1 << " " << k2 << " " << k1dotk2 << " " << p1 << " " << p2 << '\n';

					  // k, y, x=cos(theta), F2(k, y, x)
					  // printf("%e %e %e %e \n", k, p1, p2, p3);
					  fprintf(fp,"%e %e %e %e %e %e \n", a_val, k, k1, k2, p1, p2);//fprintf(fp,"%e %e %e %e \n", k, p1, p2, p3);
			 	  	}
			  	}
			}
	}

	/*close output file*/
	fclose(fp);
	return 0;
}
