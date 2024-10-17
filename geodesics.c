
#include "decs.h"

/*

this is the main photon orbit integrator 

*/
#define FAST_CPY(in,out) {out[0] = in[0]; out[1] = in[1]; out[2] = in[2]; out[3] = in[3];}
void push_photon(double X[NDIM], double Kcon[NDIM], double dl,double Xhalf[NDIM],double Kconhalf[NDIM])
{
  	double lconn[NDIM][NDIM][NDIM];
	double dKcon[NDIM];
	double Xh[NDIM], Kconh[NDIM];
	int i, j, k;

	
	/* 2nd order: scheme
	   take half-step and then evaluate derivatives 
	   at half-step. */

	/** half-step **/
	get_connection(X, lconn) ;
	
	/* advance K */
	
	for (k = 0; k < 4; k++)
		dKcon[k] = 0.;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				dKcon[k] -=
				    0.5 * dl * lconn[k][i][j] * Kcon[i] *
				    Kcon[j];
	/////////////////////////////////////////////////
	//taking advantage of symmetrices in connection instead of above, some hope for speedup here
	
	/* for (k = 0; k < 4; k++) { */
	/*   dKcon[k] = */
	/*     -2. * (Kcon[0] * */
	/* 	   (lconn[k][0][1] * Kcon[1] + */
	/* 	    lconn[k][0][2] * Kcon[2] + */
	/* 	    lconn[k][0][3] * Kcon[3]) */
	/* 	   + */
	/* 	   Kcon[1] * (lconn[k][1][2] * Kcon[2] + */
	/* 		       lconn[k][1][3] * Kcon[3]) */
	/* 	   + lconn[k][2][3] * Kcon[2] * Kcon[3] */
	/* 	   ); */
	  
	/*   dKcon[k] -= */
	/*     (lconn[k][0][0] * Kcon[0] * Kcon[0] + */
	/*      lconn[k][1][1] * Kcon[1] * Kcon[1] + */
	/*      lconn[k][2][2] * Kcon[2] * Kcon[2] + */
	/*      lconn[k][3][3] * Kcon[3] * Kcon[3] */
	/*      ); */
	/*   dKcon[k] *= 0.5*dl; */
	/* } */
	
	////////////////////////////////////////////
	
	for (k = 0; k < 4; k++)
	  Kconh[k] = Kcon[k] + dKcon[k];

	
	/* advance X */
	for (i = 0; i < 4; i++)
		Xh[i] = X[i] + 0.5 * dl * Kcon[i] ;

	/********************/
	//FAST_CPY(Xh,Xhalf);
        //FAST_CPY(Kconh,Kconhalf);
	for (i = 0; i < 4; i++) Xhalf[i]=Xh[i];
	for (i = 0; i < 4; i++) Kconhalf[i]=Kconh[i];
	/********************/

	/** full step **/
	get_connection(Xh, lconn) ;
	
	/* advance K */
	
	for (k = 0; k < 4; k++)
		dKcon[k] = 0.;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			for (k = 0; k < 4; k++)
				dKcon[k] -=
				    dl * lconn[k][i][j] * Kconh[i] *
				  Kconh[j];
	/////////////////////////////////////////////////
	//taking advantage of symmetrices in connection instead of above, some hope for speedup here
	
	/* for (k = 0; k < 4; k++) { */
	/*   dKcon[k] = */
	/*     -2. * (Kconh[0] * */
	/* 	   (lconn[k][0][1] * Kconh[1] + */
	/* 	    lconn[k][0][2] * Kconh[2] + */
	/* 	    lconn[k][0][3] * Kconh[3]) */
	/* 	   + */
	/* 	   Kconh[1] * (lconn[k][1][2] * Kconh[2] + */
	/* 		       lconn[k][1][3] * Kconh[3]) */
	/* 	   + lconn[k][2][3] * Kconh[2] * Kconh[3] */
	/* 	   ); */
	  
	/*   dKcon[k] -= */
	/*     (lconn[k][0][0] * Kconh[0] * Kconh[0] + */
	/*      lconn[k][1][1] * Kconh[1] * Kconh[1] + */
	/*      lconn[k][2][2] * Kconh[2] * Kconh[2] + */
	/*      lconn[k][3][3] * Kconh[3] * Kconh[3] */
	/*      ); */
	/*   dKcon[k] *= dl; */
	/* } */
	
	////////////////////////////////////////////////
	
	for (k = 0; k < 4; k++)
		Kcon[k] += dKcon[k];

	/* advance X */
	for (k = 0; k < 4; k++)
		X[k] += dl * Kconh[k];

	/*
	double Kcov[NDIM];
	double gcov[NDIM][NDIM];
	gcov_func(X,gcov);
	lower(Kcon, gcov, Kcov);
	double E0 = -Kcov[0];
	double L0 = Kcov[3];
	fprintf(stdout,"along g: r=%g E0=%g L0=%g \n",exp(X[1]),E0,L0);
	*/
	/* done! */
}


