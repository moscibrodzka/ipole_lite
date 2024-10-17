/**********************************************************/
/*** all you need to make a polarized radiative transfer***/
/******** used in ipole to evolve complex tensor N ********/
/******* along with standard evolution for I scalar *******/
/**********************************************************/
/**** written by Monika Moscibrodzka on 09 July 2014 ******/
/************ @ Eindhoven Airport *************************/
/************  last update: 9 May 2017   ******************/
/*************** co-author: C.F. Gammie *******************/
/**********************************************************/

#include "decs.h"
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>

/* the following definitions are used only locally */
#define S2     (1.41421356237310)	//sqrt(2)
#define S3     (1.73205080756888)	//sqrt(3)

/* declarations of local functions */
/* thermal plasma emissivity, absorptivity and Faraday conversion and rotation */
double g(double Xe);
double h(double Xe);
double Je(double Xe);
double jffunc(double Xe);
double I_I(double x);
double I_Q(double x);
double I_V(double x);
double besselk_asym(int n, double x);


/**** optional THERMAL rotativities one could use ****/
/* original Dexter 2016 */
/*
 *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
       (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
*/
/* Dexter 2016 + my correction */
/*
 *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
       (gsl_sf_bessel_Kn(0,1./Thetae) - Je(Xe)) / gsl_sf_bessel_Kn(2,1./Thetae) * cos(theta);
*/
/* Hung and Scherbakov fit */
/*
*rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
      gsl_sf_bessel_Kn(0, Thetaer) / (gsl_sf_bessel_Kn(2, Thetaer)+SMALL) * g(Xe) * cos(theta);
*/


/*invariant plasma emissivities/abs/rot in tetrad frame */
void jar_calc(double X[NDIM], double Kcon[NDIM],
	      double *jI, double *jQ, double *jU, double *jV,
	      double *aI, double *aQ, double *aU, double *aV,
	      double *rQ, double *rU, double *rV)
{
   
    double nu, Thetae, Ne, B, theta, nusq;
    double x, Xe, omega0, nuc;
    double Bnuinv;
    double Ucov[NDIM],Ucon[NDIM],Bcon[NDIM],Bcov[NDIM];
    double Thetaer, wp2;
    double Gcov[NDIM][NDIM];
    
    gcov_func(X, Gcov);
    get_fluid_params(X, Gcov, &Ne, &Thetae, &B, Ucon, Ucov, Bcon, Bcov);
    theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b  */
    nu = get_fluid_nu(Kcon, Ucov, X);	/* freqcgs1;  freq in Hz */
    
    if (Ne <= SMALL || nu==1. || nu > 1e24) {  

	 *jI = 0.0; *jQ = 0.0;
	 *jU = 0.0; *jV = 0.0;
	 
	 *aI = 0.0; *aQ = 0.0;
	 *aU = 0.0; *aV = 0.0;
	 
	 *rQ = 0.0; *rU = 0.0;
	 *rV = 0.0;

	 return;

    } else if (theta <= 0. || theta >= M_PI) {	/* no emission/absorption along field only Faraday rotation */

	*jI = 0.0; *jQ = 0.0;
	*jU = 0.0; *jV = 0.0;

	*aI = 0.0; *aQ = 0.0;
	*aU = 0.0; *aV = 0.0;

	*rQ = 0.0; *rU = 0.0;

	nu = get_fluid_nu(Kcon, Ucov, X);	/* freqcgs1;  freq in Hz */
	wp2 = 4. * M_PI * Ne * EE * EE / ME;
	omega0 = EE * B / ME / CL;
	/* Faraday rotativities for thermal plasma */

	  
	  Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));
	  Thetaer = 1. / Thetae;
	  //this is used for EHT-library only
	  if (Thetae > 3.0) {
	    // High temperature: use approximations to bessel
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	  } else if (0.2 < Thetae && Thetae <= 3.0) {
	    // Mid temperature: use real bessel functions (TODO fit?)
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	      (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	    //(gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
	  } else if (Thetae <= 0.2) {
	    // Use the constant low-temperature limit
	    *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
	  }
	
	  *rV *= nu;
	
	return;
	
    } else {

      //nu = get_fluid_nu(Kcon, Ucov, X);	/* freqcgs1;  freq in Hz */
      nusq = nu * nu;
      //      B = get_model_b(X);	/* field in G                */
      //Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
      omega0 = EE * B / ME / CL;
      wp2 = 4. * M_PI * Ne * EE * EE / ME;
      /* Faraday rotativities for thermal plasma */
      //Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));
	
        Xe = Thetae * sqrt(S2 * sin(theta) * (1.e3 * omega0 / 2. / M_PI / nu));
	Thetaer = 1. / Thetae;
	/* Here I use approximate bessel functions to match rhoqv with grtrans */
	*rQ = 2. * M_PI * nu / 2. / CL * wp2 * omega0 * omega0 / pow(2 * M_PI * nu, 4) *
	  jffunc(Xe) * (besselk_asym(1, Thetaer) / besselk_asym(2, Thetaer) +
			6. * Thetae) * sin(theta) * sin(theta);
	*rU = 0.0;
	//this is used for EHT-library only
	if (Thetae > 3.0) {
	  // High temperature: use approximations to bessel
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	    (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	} else if (0.2 < Thetae && Thetae <= 3.0) {
	  // Mid temperature: use real bessel functions (TODO fit?)
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) *
	    (besselk_asym(0, Thetaer) - Je(Xe)) / besselk_asym(2, Thetaer) * cos(theta);
	  //	  (gsl_sf_bessel_Kn(0, Thetaer) - Je(Xe)) / gsl_sf_bessel_Kn(2, Thetaer) * cos(theta);
	} else if (Thetae <= 0.2) {
	  // Use the constant low-temperature limit
	  *rV = 2.0 * M_PI * nu / CL * wp2 * omega0 / pow(2. * M_PI * nu, 3) * cos(theta);
	}
	
	/* invariant rotativities */
	*rQ *= nu;
	*rV *= nu;

      
      //// end of Faraday effects ////

	/* synchrotron emissivity */ 
	nuc = 3.0 * EE * B * sin(theta) / 4.0 / M_PI / ME / CL * Thetae * Thetae + 1.0;
	x = nu / nuc;
	//user is responsible for emissivities, there is no guarantee that this will work for your problem
	*jI =  jnu_synch(nu,Ne,Thetae,B,theta);
	//*jI = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_I(x); // [g/s^2/cm = ergs/s/cm^3]
	if( *jI > 0.0 ){
	  *jQ = Ne * EE * EE * nu / 2. / S3 / CL / Thetae / Thetae * I_Q(x); // here mistake in grtrans
	  *jU = 0.0;	                                        	 // convention; depends on tetrad
	  *jV = 2. * Ne * EE * EE * nu / tan(theta) / 3. / S3 / CL / Thetae / Thetae / Thetae * I_V(x);
	}else{
	  *jQ=0.0;
	  *jU=0.0;
	  *jV=0.0;
	}
	/* invariant emissivity */
	*jI /= nusq;
	*jQ /= nusq;
	*jU /= nusq;
	*jV /= nusq;

	/*invariant synchrotron absorptivity */
	Bnuinv = Bnu_inv(nu, Thetae) + 1e-100;   /* Planck function */
	*aI = *jI / Bnuinv;
	*aQ = *jQ / Bnuinv;
	*aU = *jU / Bnuinv;
	*aV = *jV / Bnuinv;
      
      return;
      
    }
    


}
/*emissivity functions and functions used for Faraday conversion and rotation*/
/*from J. Dexter PhD thesis (checked with Leung harmony program, and Huang & Shcherbakov 2011*/
double g(double Xe)
{
    return 1. - 0.11 * log(1 + 0.035 * Xe);
}

double h(double Xe)
{
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
	cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
	0.011 * exp(-Xe / 47.2);
}

double Je(double Xe)
{
    return 0.43793091 * log(1. + 0.00185777 * pow(Xe, 1.50316886));
}

double jffunc(double Xe)
{
    double extraterm;
    extraterm =
	(0.011 * exp(-Xe / 47.2) -
	 pow(2., -1. / 3.) / pow(3.,
				 23. / 6.) * M_PI * 1e4 * pow(Xe + 1e-16,
							      -8. / 3.)) *
	(0.5 + 0.5 * tanh((log(Xe) - log(120.)) / 0.1));
    return 2.011 * exp(-pow(Xe, 1.035) / 4.7) -
	cos(Xe * 0.5) * exp(-pow(Xe, 1.2) / 2.73) -
	0.011 * exp(-Xe / 47.2) + extraterm;
}

double I_I(double x)
{
    return 2.5651 * (1 + 1.92 * pow(x, -1. / 3.) +
		     0.9977 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
								    1. /
								    3.));
}

double I_Q(double x)
{
    return 2.5651 * (1 + 0.93193 * pow(x, -1. / 3.) +
		     0.499873 * pow(x, -2. / 3.)) * exp(-1.8899 * pow(x,
								      1. /
								      3.));
}

double I_V(double x)
{
    return (1.81348 / x + 3.42319 * pow(x, -2. / 3.) +
	    0.0292545 * pow(x, -0.5) + 2.03773 * pow(x,
						     -1. / 3.)) *
	exp(-1.8899 * pow(x, 1. / 3.));
}

double besselk_asym(int n, double x)
{

    if (n == 0)
	return -log(x / 2.) - 0.5772;

    if (n == 1)
	return 1. / x;

    if (n == 2)
	return 2. / x / x;

    fprintf(stderr,"this cannot happen\n");
    exit(1);
}


/* end of emissivity functions */

#undef S2
#undef S3

//optional
int radiating_region(double X[4])
{

    int i,j,k;
    double del[NDIM],Ne,B,sigma_m;
    //cut RT based on plasma density
    Ne=get_model_ne(X);
    if(Ne > SMALL*10 && X[1]<log(100.)) return(1);
    else return(0);
}


double jnu_synch(double nu, double Ne, double Thetae, double B, double theta)
{
    double K2,nuc,nus,x,f,j,sth ;
    
    if (Thetae < THETAE_MIN){
      return 0.;
    }
    if (Thetae > 1e2){
      K2=2. * Thetae * Thetae;
    }else{
      //K2 = besselk_asym(2, 1./Thetae);
      //to be accurate one should use below
      K2 = gsl_sf_bessel_Kn(2,1./Thetae);
    }

    nuc = EE*B/(2.*M_PI*ME*CL);
    sth = sin(theta);
    nus = (2./9.)*nuc*Thetae*Thetae*sth;
    if (nu > 1.e12 * nus) return (0.);
    x = nu / nus;
    double xxx,xp1;
    xp1 = pow(x, 1. / 3.);
    xxx = sqrt(x) + pow(2.,11./12.) * sqrt(xp1);
    f = xxx * xxx;
    j = (M_SQRT2 * M_PI * EE * EE * Ne * nus / (3. * CL * K2)) * f *
      exp(-xp1);
    
    return(j) ;

}

/* get the invariant emissivity and opacity at a given position
   for a given wavevector */
void get_jkinv(double X[NDIM], double Kcon[NDIM], double *jnuinv,
	       double *knuinv, double col_v[3])
{
    double nu, theta, B, Thetae, Ne, Bnuinv;
    double Ucov[NDIM], Bcov[NDIM];
    double Ucon[NDIM], Bcon[NDIM];
    double Kcov[NDIM], gcov[NDIM][NDIM];

    /* get fluid parameters */
    Ne = get_model_ne(X);	/* check to see if we're outside fluid model */
    B = get_model_b(X);		/* field in G */

    if (Ne <= 0.){
      *jnuinv = 0.;
      *knuinv = 0.;
      return;
    }

    /* get covariant four-velocity of fluid for use in get_bk_angle and get_fluid_nu */
    get_model_ucov(X, Ucov);
    get_model_ucon(X, Ucon);
    
    gcov_func(X, gcov);
    lower(Kcon, gcov, Kcov);

    theta = get_bk_angle(X, Kcon, Ucov);	/* angle between k & b */
    if (theta <= 0. || theta >= M_PI) {	/* no emission along field */
	*jnuinv = 0.;
	*knuinv = 0.;
	return;
    }

    nu = get_fluid_nu(Kcon, Ucov,X);	 /* freq in Hz */
    Thetae = get_model_thetae(X);	/* temp in e rest-mass units */
    
    Bnuinv = Bnu_inv(nu, Thetae);
    *jnuinv = jnu_inv(nu, Thetae, Ne, B, theta);
    
    if (Bnuinv < SMALL)
	*knuinv = SMALL;
    else{
	*knuinv = *jnuinv / Bnuinv;
    }
    

    return;

}

