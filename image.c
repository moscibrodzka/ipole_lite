#include "decs.h"
#include <ctype.h>

int compare_doubles(const void *a, const void *b)
{
	const double *da = (const double *) a;
	const double *db = (const double *) b;

	return (*da > *db) - (*da < *db);
}

#define TOP_FRAC	(0.005)	/* Brightest and dimmest pixels cut out of color scale */
#define BOT_FRAC	(0.005)

void make_ppm(double p[NX][NY], double freq, char filename[])
{

	int i, j, k, npixels;
	double *q, min, max;
	FILE *fp;
	double *alloc_double_array( int ndata )  ;

	npixels = NX * NY;

	q = alloc_double_array(npixels);
	k = 0 ;
        for (i = 0; i < NX; i++)
	for (j = 0; j < NY; j++) {
		q[k] = p[i][j];
		k++ ;
	}

	qsort(q, npixels, sizeof(double), compare_doubles);
	if (q[0] < 0) {		/* must be log scaling */
		min = q[(int) (npixels * BOT_FRAC)];
	} else {
		i = 0;
		while (i < npixels - 1 && q[i] == 0.)
			i++;
		min = q[i];
	}
	max = q[npixels - (int) (npixels * TOP_FRAC)];

	if ((fp = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "Failed to open ppm file %s\n", filename);
		exit(125);
	}

	/* write out header information */
	fprintf(fp,
		"P6\n#  min=%g  , max=%g \n#  frequency=%g \n%d %d\n%d\n",
		min, max, freq, NX, NY, 255);
	fflush(fp);
	
	int red,green,blue;
	for (j = NY-1; j >= 0; j--) 
	for (i = 0; i < NX; i++) 
	{
	    if(RAINBOW){
		rainbow_palette(p[i][j],min,max,&red,&green,&blue) ;
	    }
	    if(BW){
		monika_palette(p[i][j],min,max,&red,&green,&blue) ;
	    }
	    if(AFMHOT){
		afmhot_palette(p[i][j],min,max,&red,&green,&blue) ;
	    }
	    fputc((char) red, fp);
	    fputc((char) green, fp);
	    fputc((char) blue, fp);
	}

	fclose(fp);

	return;
}

#undef TOP_FRAC
#undef BOT_FRAC

/* 
        rainbow color palette, based on john.pal

        input: integer 0-255 
        output: red, green, blue integers 

        author: Bryan M. Johnson
*/

void rainbow_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue)
{
  double aa, b, c, d, e, f;
  double x, y;
  double max_min = max - min ;

  if(max_min > 0.0) { // trust no one
    x = (data - min)/(max_min) ;
    
    /* ========== Red ============ */
    aa = 4.0*x - 1.52549019607844;
    b = 4.52941176470589 - 4.0*x;
    y = aa < b ? aa : b;
    *pRed = (int)(255.0*y);
    *pRed = *pRed >   0 ? *pRed :   0;
    *pRed = *pRed < 255 ? *pRed : 255;

    /* ========== Green ========== */
    aa = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = aa < b ? aa : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    aa = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = aa < b ? aa : b;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  }
  else {
    *pRed = *pGreen = *pBlue = (data > max ? 255: 0) ;
  }

  return;
}



//black & white
void monika_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue)
{
  double aa, b, c, d, e, f;
  double x, y;
  double max_min = max - min ;

  if(max_min > 0.0) { // trust no one
    x = (data - min)/(max_min) ;
    
    /* ========== Red ============ */
     aa = 4.0*x - 1.52549019607844;
     b = 4.52941176470589 - 4.0*x;
     y = aa < b ? aa : b;
     y = x;
     *pRed = (int)(255.0*y);
     *pRed = *pRed >   0 ? *pRed :   0;
     *pRed = *pRed < 255 ? *pRed : 255;
     
    /* ========== Green ========== */
    aa = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = aa < b ? aa : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    y = x;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    aa = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = aa < b ? aa : b;
    y = x;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  }
  else {
    *pRed = *pGreen = *pBlue = (data > max ? 255: 0) ;
  }

  return;
}

// AFMHOT
void afmhot_palette(double data, double min, double max, int *pRed, int *pGreen, int *pBlue)
{
  double aa, b, c, d, e, f;
  double x, y;
//  max=max*0.1;
  double max_min = max - min ;

  if(max_min > 0.0) { // trust no one
    x = (data - min)/(max_min) ;

    /* ========== Red ============ */
     aa = 4.0*x - 1.52549019607844;
     b = 4.52941176470589 - 4.0*x;
     y = aa < b ? aa : b;
     y = 2*x;
     *pRed = (int)(255.0*y);
     *pRed = *pRed >   0 ? *pRed :   0;
     *pRed = *pRed < 255 ? *pRed : 255;
     
    /* ========== Green ========== */
    aa = 4.0*x - 0.521568627450979;
    b = 2.52549019607844 - 4.0*x;
    c = aa < b ? aa : b;
    d = 4.0*x - 1.53725490196073;
    e = 3.52941176470581 - 4.0*x;
    f = d < e ? d : e;
    y = c > f ? c : f;
    y = 2*x-0.5;
    *pGreen = (int)(255.0*y);
    *pGreen = *pGreen >   0 ? *pGreen :   0;
    *pGreen = *pGreen < 255 ? *pGreen : 255;

    /* ========== Blue =========== */
    aa = 4.0*x + 0.498039215686276;
    b = 2.50980392156862 - 4.0*x;
    y = aa < b ? aa : b;
    y = 2*x-1;
    *pBlue = (int)(255.0*y);
    *pBlue = *pBlue >   0 ? *pBlue :   0;
    *pBlue = *pBlue < 255 ? *pBlue : 255;

  }
  else {
    *pRed = *pGreen = *pBlue = (data > max ? 255: 0) ;
  }

  return;
}


/*********************************************************************************************
  alloc_double_array():
     -- returns with a pointer to an array of size "ndata"
 *********************************************************************************************/
double *alloc_double_array( int ndata ) 
{ 
  double *pa;

  pa = (double *) calloc(ndata,sizeof(double));
  if( pa == NULL ) { 
    fprintf(stderr,"Error allocating a double array of length = %d \n", ndata);
    fprintf(stderr,"....Exiting...\n");
    fflush(stderr);
    exit(1);
  }
  return(pa);
}


