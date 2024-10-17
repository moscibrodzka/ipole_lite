// or which source SGRA,M87,DABHB, if NT model these all 0
#define SOURCE_SGRA 1
#define SOURCE_M87  0
#define SOURCE_DABHB 0

// image resolution
#define NX 64
#define NY 64
// FOV in [M] 
#define DX 40.
#define DY 40.

//floor and ceiling of electron temperature, model dependent
#define THETAE_MAX 1000.
#define THETAE_MIN 0.01

#define NPRIM	10 // 10 variables in hdf5 in fmks
//#define NPRIM	8 // default

//chose integration scheme in ipolarray.c
#define INT_FULL  0 //full integration step
#define INT_SPLIT 1 //split: 1/2rot + ea + 1/2rot 

//colortheme for ppm files
#define RAINBOW 0
#define AFMHOT 1
#define BW 0

/********************** RATHER DO NOT CHANGE BELOW *************/

#define MAXNSTEP        60000
#define POLARIZATION_ON (1)







