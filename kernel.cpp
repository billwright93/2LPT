#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <cfloat>
#include <cmath>
#include <stdlib.h>



//Kernel array sizes

//angle (x) sampling
const int n1=200;

// magnitude (r) sampling
const int n2=200;

//F1[k],G1[k]

double F1_nk;
double G1_nk;

//F1[k-p],G1[k-p]
double F1kmp_nk[(n1+1)*n2];
double G1kmp_nk[(n1+1)*n2];

//F1[p], G1[p]
double F1p_nk[(n1+1)*n2];
double G1p_nk[(n1+1)*n2];


//2nd order used in Bispectrum term : F2/G2[p,k]
double F2A_nk[(n1+1)*n2];
double G2A_nk[(n1+1)*n2];

//2nd order used in Bispectrum term : F2/G2[-p,k]
double F2B_nk[(n1+1)*n2];
double G2B_nk[(n1+1)*n2];

//2nd order used in Bispectrum term : F2/G2[-k,k-p]
double F2C_nk[(n1+1)*n2];
double G2C_nk[(n1+1)*n2];

//2nd order used in P22 : F2/G2[k-p, p]
double F2_nk[(n1+1)*n2];
double G2_nk[(n1+1)*n2];




////////  NUMERICAL KERNEL CONSTRUCTION //////////

// angular limits
//XMAX is set to less than 1 to avoid 2nd order kernel singularity at r->1, x = 1.
//XMIN is set to greater than -1 because of numerical singularity
const double XMAX = 0.999999;
const double XMIN = -1+1e-10;



// Normalized Hubble parameter H and a*H*dH/da for the background (LCDM)

static double HA(double a, double omega0){
	double omegaL= 1.-omega0;
	return  sqrt(omega0/pow(a,3)+omegaL);}

static double HA1(double a, double omega0){
	return -3*omega0/(2.*pow(a,3));}


  //1-mu^2
  static double ker1(double u1){
  	return 1.-u1*u1;
  }

  //alpha(k1,k2)
  static double alpha(double k1, double k2, double u1){
  	return 1.+k2*u1/k1;
  }

  //Symmetrized alpha : [alpha(k1,k2)+alpha(k2,k1)]/2
  static double alphas(double k1, double k2, double u1){
  	return (alpha(k1,k2,u1)+alpha(k2,k1,u1))/2.;
  }

  //beta(k1,k2)
  static double beta1(double k1, double k2, double u1){
  	return u1*(k1*k1+k2*k2+2*k1*k2*u1)/(2*k1*k2);
  }


  //beta function for nDGP
  static double beta(double a, double omega0, double omegarc){
  	return 1.+HA(a,omega0)/sqrt(omegarc)*(1.+HA1(a,omega0)/(3.*HA(a,omega0)*HA(a,omega0)));}


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
  // 2) Add in extra parameters to params of funcn1 system of ODEs
  // 3) Add in extra parameters to initn function and in params passed by function to solver (this is right after initial conditions are passed)

  static double mu(double a, double k0, double omega0, double p1, double p2, double p3 ){
  	double h0 = 1./2997.9;
  //	return 1.; // GR
  //	return 1. + pow2(k0/a)/(3.*(pow2(k0/a)+pow3(omega0/pow3(a)-4.*(omega0-1.))/(2.*p1/pow2(h0)*pow2(4-3.*omega0)))); //f(R) Hu- Sawicki
  	return 1.+1./(3.*beta(a,omega0,p1)); //nDGP
  }


  static double gamma2(double a, double omega0, double k0, double k1, double k2, double u1, double p1, double p2, double p3 ){
  	double h0 = 1./2997.9;
  //	 return 0. ; // GR

  /*	 return -(9.*pow2(k0/a)*pow2(omega0/pow3(a))*pow(omega0-4.*pow3(a)*(-1+omega0),5))/
  			    (48.*pow(a,15)*pow2(p1/pow2(h0))*pow2(HA(a,omega0))*pow4(3.*omega0-4.)
  			   *(pow2(k0/a)+pow3(omega0-4.*pow3(a)*(omega0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*omega0-4.)))
  			   *(pow2(k1/a)+pow3(omega0-4.*pow3(a)*(omega0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*omega0-4.)))
  			   *(pow2(k2/a)+pow3(omega0-4.*pow3(a)*(omega0-1.))/(2.*pow(a,9)*p1/pow2(h0)*pow2(3.*omega0-4.)))); //f(R) Hu- Sawicki
  */
       	return -1.*p2/(HA(a,omega0)*HA(a,omega0)*24.*pow(beta(a,omega0,p1),3)*p1)*pow(omega0/(a*a*a),2)*ker1(u1); //nDGP
  }



  ///////// NON-SEPARABLE //////////////


	//Jacobian of the system required when calling the system evolver, the below is neither needed, nor correct for solving
	int jac (double a, const double G[], double *dfdy, double dfdt[], void *params)
	{
		return GSL_SUCCESS;
	}


  /* Euler and Continuity equations for numerical kernels */

  //HA2 = -dH/dt/H^2
  static double HA2(double a, double omega0){
  	return 3.*omega0/(2.*HA(a,omega0)*HA(a,omega0)*pow(a,3))
  	;}

  /* Parameters passed to system of Euler and continuity equations*/
  // k (magnitude) and x (angular) values for the system of equations
  // k1.k2=k1k2x2 , k1.k3 = k1k3x3, k2.k3=k2k3x1
  // EDIT : Add in new gravity parameters to this list as required (par1 is the nDGP parameter omega_rc as default)
  struct param_type3 {
    double kk1;
    double xx1;
  	double kk2;
    double xx2;
  	double kk3;
    double xx3;
  	double omega00;
  	double par1;
  	double par2;
  	double par3;
  };


  // F1(k1,k2)=G[0], G1(k1,k2)=G[1], F3(k1,k2,k3)=G[10], G3[k1,k2,k3]=G[11]
  // Other G[i] are used to construct symmetrized kernels used in G3,F3 equations as well as RSD A-term

  int funcn1(double a, const double G[], double F[], void *params)
  {
  	param_type3 p = *(param_type3 *)(params);
  	double k1 = p.kk1;
  	double x1= p.xx1;
  	double k2 = p.kk2;
  	double x2= p.xx2;
  	double k3 = p.kk3;
  	double x3= p.xx3;
  	double omega0 = p.omega00;
  	double p1 = p.par1;
  	double p2 = p.par2;
  	double p3 = p.par3;

		/* 1st order */
		//1. F1/G1(k)
		F[0] = -G[1]/a;
		F[1] =1./a*(-(2.-HA2(a,omega0))*G[1]-HA2(a,omega0)*G[0]*mu(a,k1,omega0,p1,p2,p3));

		//2. F1/G1(k-p)
		F[2] = -G[3]/a;
		F[3] =1./a*(-(2.-HA2(a,omega0))*G[3]-HA2(a,omega0)*G[2]*mu(a,sqrt(k2*k2+k1*k1-2*k2*k1*x2),omega0,p1,p2,p3));

		//3. F1/G1(p)
		F[4] = -G[5]/a;
		F[5] =1./a*(-(2.-HA2(a,omega0))*G[5]-HA2(a,omega0)*G[4]*mu(a,k2,omega0,p1,p2,p3));

		/* 2nd order */

		//4. F2/G2(p,k-p) (P22)
		F[6] =1/a*(-(alpha(k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))*G[5]*G[2]+alpha(sqrt(k2*k2+k1*k1-2*k2*k1*x2),k2,(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))*G[3]*G[4])/2.-G[7]);
		F[7] =1/a*(-(2.-HA2(a,omega0))*G[7]-HA2(a,omega0)*G[6]*mu(a,k1,omega0,p1,p2,p3) - gamma2(a, omega0, k1, k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2),p1,p2,p3)*G[4]*G[2] - beta1(k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))*G[5]*G[3]);


		//5. F2/G2(-k,k-p)
		F[8] =1/a*(-(alpha(k1,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))*G[1]*G[2]+alpha(sqrt(k2*k2+k1*k1-2*k2*k1*x2),k1,(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))*G[3]*G[0])/2.-G[9]);
		F[9] =1/a*(-(2.-HA2(a,omega0))*G[9]-HA2(a,omega0)*G[8]*mu(a,k2,omega0,p1,p2,p3) - gamma2(a, omega0, k1, k1, sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2),p1,p2,p3)*G[2]*G[0] - beta1(k1,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))*G[3]*G[1]);



		//6. F2/G2(p,-p)
		F[10] =1./a*(-alpha(k2,k3,x1)*G[5]*G[4]-G[11]) ;
		F[11] =1./a*(-(2.-HA2(a,omega0))*G[11]-HA2(a,omega0)*G[10]*mu(a,0,omega0,p1,p2,p3) - gamma2(a, omega0, k1, k2,k3,x1,p1,p2,p3)*G[4]*G[4]-beta1(k2,k3,x1)*G[5]*G[5]);



		//7. F2/G2(p,k)
		F[12] =1/a*(-(alpha(k2,k1,x2)*G[5]*G[0]+alpha(k1,k2,x2)*G[1]*G[4])/2.-G[13]) ;
		F[13] =1/a*(-(2.-HA2(a,omega0))*G[13]-HA2(a,omega0)*G[12]*mu(a,sqrt(k2*k2+k1*k1+2*k2*k1*x2),omega0,p1,p2,p3) - gamma2(a, omega0, k1, k2,k1,x2,p1,p2,p3)*G[4]*G[0]-beta1(k2,k1,x2)*G[5]*G[1]);



		//8. F2/G2(-p,k)=F2/G2(p,-k)
		F[14] =1/a*(-(alpha(k3,k1,x3)*G[5]*G[0]+alpha(k1,k3,x3)*G[1]*G[4])/2.-G[15]) ;
		F[15] =1/a*(-(2.-HA2(a,omega0))*G[15]-HA2(a,omega0)*G[14]*mu(a,sqrt(k2*k2+k1*k1-2*k2*k1*x2),omega0,p1,p2,p3) -gamma2(a, omega0, k1, k3,k1,x3,p1,p2,p3)*G[4]*G[0]-beta1(k3,k1,x3)*G[5]*G[1]);


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

  int initn(double A, double YMIN, double YMAX, double k, double omega0, double par1, double par2, double par3)
  {
  //#pragma omp parallel for schedule(dynamic)
  	for ( int i = 0; i <= n1; ++i)
  	{
  		for (int j = 0; j< n2; ++j)
  		{
  				double x1,x2,x3,k1,k2,k3;

  				/* k - magnitude sampling */
  				    double R = (YMIN) * exp(j*log(YMAX/YMIN)/(n2*1.-1.)); //exponential sampling
  				//	double R = (j*1.)/(n2*1.-1)*YMAX + YMIN; // linear sampling
  				//	double R = pow2(j*1./(n2*1.-1))*YMAX+YMIN; //quadratic sampling

  				//k,k1,-k1

  				/*Parameter selection for kernels */

  				k1 = k;
  		    k2 = k*R;
  				k3 = k2;

  				x1 = XMIN;
  		    	if (i>=n1/2.+1.) {
  					x2 = XMAX - pow((i-n1*1./2.-1.)/(n1/2.-1.),2)*(XMAX-2./n1); // quadratic sampling
  				//	x2 = XMAX*exp((k-n1/2.-1)*log(1./(n1/2.))/(n1/2-1)); //exponential sampling
  				}
  				else if (i<=n1/2.-1.) {
  					x2 = XMIN + pow(i/(n1/2.-1.),2)*(XMAX-2./n1);
  				//	x2 = XMIN*exp(k*log(1./(n1/2))/(n1/2-1));// exponential sampling
  				}
  				else {
  					x2=0.;
  				}
  				//x2 = XMIN + (XMAX-XMIN)*k/n1; //linear sampling
  				x3 = -x2;

					// Initial scale factor for solving system of equations
					//Sets 1-3
						double a = 0.0001;

						// Einstein de Sitter initial condtions for 2nd order kernels

						//4. F2/G2(p,k-p) (P22)
						double Edsf2=a*a*(10./14.*alphas(k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))+ 2./7.*beta1(k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2)));
						double Edsg2=a*a*(6./14.*alphas(k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))+ 4./7.*beta1(k2,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k1*x2-k2)/sqrt(k2*k2+k1*k1-2*k2*k1*x2)));

						//5. F2/G2(-k,k-p)
						double EdsF2C=a*a*(10./14.*alphas(k1,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2))+ 2./7.*beta1(k1,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2)));
						double EdsG2C=a*a*(6./14.*alphas(k1,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2)) + 4./7.*beta1(k1,sqrt(k2*k2+k1*k1-2*k2*k1*x2),(k2*x2-k1)/sqrt(k2*k2+k1*k1-2*k2*k1*x2)));

						//6. F2/G2(p,-p)
						double EdsF2A= a*a*(10./14.*alphas(k2,k3,x1)+ 2./7.*beta1(k2,k3,x1));
						double EdsG2A= a*a*(6./14.*alphas(k2,k3,x1)+ 4./7.*beta1(k2,k3,x1));
						//7. F2/G2(p,k)
						double EdsF2B= a*a*(10./14.*alphas(k2,k1,x2)+ 2./7.*beta1(k2,k1,x2));
						double EdsG2B= a*a*(6./14.*alphas(k2,k1,x2) + 4./7.*beta1(k2,k1,x2));
						//8. F2/G2(-p,k)=F2/G2(p,-k)
						double Edsf2d= a*a*(10./14.*alphas(k3,k1,x3)+ 2./7.*beta1(k3,k1,x3));
						double Edsg2d= a*a*(6./14.*alphas(k3,k1,x3)+ 4./7.*beta1(k3,k1,x3));

						double G[16] = { a,-a,a,-a,a,-a,  Edsf2 ,Edsg2,  EdsF2C, EdsG2C, EdsF2A, EdsG2A , EdsF2B,  EdsG2B, Edsf2d, Edsg2d}; // initial conditions
  // Non-Eds ICs
  		//	  double G[16] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  			/*Parameters passed to system of equations */
  // EDIT : If more than one gravity parameter is used, add them after p1
  				struct param_type3 my_params1 = {k1,x1,k2,x2,k3,x3,omega0, par1, par2, par3};

  				gsl_odeiv2_system sys = {funcn1, jac, 16, &my_params1};

  			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
  				gsl_odeiv2_driver * d =
  				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
  										   1e-6, 1e-6, 0.0);

  				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


  				/*Allocation of array values */

  			//F1(k;a), G1(k;a)
  			F1_nk = G[0];
  			G1_nk = G[1];

  			// F1(k-p;a), G1(k-p;a)
  			F1kmp_nk[i*n2 + j] = G[2];
  			G1kmp_nk[i*n2 + j] = G[3];

  			//F1(p;a), G1(p;a)
  			F1p_nk[i*n2 + j] = G[4];
  			G1p_nk[i*n2 + j] = G[5];

  			/*2nd order*/

  			//F2/G2(p,k-p) (P22)
  			F2_nk[i*n2 + j] =  G[6];
  			G2_nk[i*n2 + j] =  G[7];

  			//F2/G2(-k,k-p)
  			F2C_nk[i*n2 + j] =  G[8];
  			G2C_nk[i*n2 + j] =  G[9];

  			//F2/G2(p,k)
  			F2A_nk[i*n2 + j] =  G[12];
  			G2A_nk[i*n2 + j] =  G[13];

  			//F2/G2(p,-k)
  			F2B_nk[i*n2 + j] =  G[14];
  			G2B_nk[i*n2 + j] =  G[15];

  			gsl_odeiv2_driver_free(d);
  		}
  	}

  	return 0;
  }



/// OUTPUT SECTION ////

int main(int argc, char* argv[]) {
//output file name
	const char* output = "kernelb.dat";

/* Open output file */
	FILE* fp = fopen(output, "w");


	// Number of k-modes you wish to output between kmin and kmax
		int Nk = 1;
		double kmin = 0.01;
		double kmax = 0.2;

		double qmin = 0.001;
		double qmax = 5.;

double p1,p2,p3;

for(int i = 0 ; i < Nk;  i ++) {

// Example k value
double k = kmin + i*(kmax-kmin)/Nk;

double ymin = qmin/k;
double ymax = qmax/k;
// Format : initn(scale factor, ymin,ymax, k., omega_(total matter), theory parameter1, theory parameter2 ...)
// Parameterized magnitude integration q = ky
// set k = 1. for scale independant and take the function out of the loop
 initn(1., ymax, ymin ,k, 0.24, 0.0001, 1., 1.);

for(int j = 0; j<1; j++){
// Value of y where k2 = yk
p1 = (ymin) * exp(j*log(ymax/ymin)/(n2*1.-1.)); //exponential sampling

for(int m = 0; m<1; m++){

//Value of k1.k2/k1k2 = cos(theta)
      	if (i>=n1/2.+1.) {
  					p2 = XMAX - pow((m-n1*1./2.-1.)/(n1/2.-1.),2)*(XMAX-2./n1); // with quadratic sampling
  				}
  				else if (i<=n1/2.-1.) {
  					p2 = XMIN + pow(m/(n1/2.-1.),2)*(XMAX-2./n1);
  				}
  				else {
  					p2=0.;
  				}

//F2/G2(p,k)
p3 = F2A_nk[j*n2+m];
// printf("%e %e %e %e \n", k, p1, p2, p3);
 fprintf(fp,"%e %e %e %e \n", k, p1, p2, p3);
}
}
}

/*close output file*/
	fclose(fp);
	return 0;
}
