/* 
model-dependent routines for integrating geodesics, 
including:
stop_backward_integration
stepsize
*/
#include "decs.h"

/* condition for stopping the backwards-in-lambda integration of the photon geodesic, default values */
#define LRMAX (log(50.))
#define LRMIN (log(1.05*Rh)) 
int stop_backward_integration(double X[NDIM], double Kcon[NDIM], double Xcam[NDIM])
{
    if ((X[1] > LRMAX && Kcon[1] < 0.) ||       /* out far */
        X[1] < LRMIN            /* in deep */
	){ 
      return (1);
    }else{
      return (0);             /* neither out far nor in deep */
    }
}
#undef LRMIN
#undef LRMAX

double stepsize(double X[NDIM], double Kcon[NDIM])
{
        double dl, dlx1, dlx2, dlx3;
        double idlx1,idlx2,idlx3 ;
	//stepsize paramter, default value
	double eps=0.03;
	
        dlx1 = eps / (fabs(Kcon[1]) + SMALL*SMALL) ;
	dlx2 = eps * MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + SMALL*SMALL) ;
	dlx3 = eps / (fabs(Kcon[3]) + SMALL*SMALL) ;

        idlx1 = 1./(fabs(dlx1) + SMALL*SMALL) ;
        idlx2 = 1./(fabs(dlx2) + SMALL*SMALL) ;
        idlx3 = 1./(fabs(dlx3) + SMALL*SMALL) ;

        dl = 1. / (idlx1 + idlx2 + idlx3) ;

        return (dl);
}

