/* == Definition of constants in CGS units == */

#define MPCL2                           (0.0015033              )

#define EE				(4.80320680e-10		) /* electron charge */
#define CL				(2.99792458e10		) /* speed of light */
#define ME				(9.1093826e-28		) /* electron mass */
#define MP				(1.67262171e-24		) /* proton mass */
#define MN				(1.67492728e-24		) /* neutron mass */
#define AMU				(1.66053886e-24		) /* atomic mass unit */
#define HPL				(6.6260693e-27		) /* Planck constant */
#define HPL_MECL2                       (8.09e-21               ) /* hpl/m_ec^2*/
#define HBAR			        (HPL/(2.*M_PI)		) /* Planck's consant / 2pi */
#define KBOL			        (1.3806505e-16		) /* Boltzmann constant */
#define GNEWT			        (6.6742e-8		) /* Gravitational constant */
#define SIG				(5.670400e-5		) /* Stefan-Boltzmann constant */
#define RGAS			        (8.3143e7		)	/* erg K^-1 mole^-1: ideal gas const */
#define EV				(1.60217653e-12		) /* electron volt in erg */
#define SIGMA_THOMSON	                (0.665245873e-24	) /* Thomson cross section in cm^2 */
#define JY				(1.e-23			) /* Jansky (flux/freq. unit) in cgs */

#define ALPHAF                          (7.29735e-3             ) /* (fine-structure constant)^-1 */
#define KEV                             (1.602e-9               ) /*kiloelectronovolt*/

#define PC				(3.085678e18		) /* parsec */
#define AU				(1.49597870691e13	) /* Astronomical Unit */

#define YEAR			        (31536000.		) /* No. of seconds in year */
#define DAY				(86400.			) /* No. of seconds in day  */
#define HOUR			        (3600.			) /* No. of seconds in hour */

#define MSUN			(1.989e33			) /* solar mass */
#define RSUN			(6.96e10			) /* Radius of Sun */
#define LSUN			(3.827e33				) /* Luminousity of Sun */
#define TSUN			(5.78e3				) /* Temperature of Sun's photosphere */

#define MEARTH			(5.976e27			) /* Earth's mass */
#define REARTH			(6.378e8			) /* Earth's radius */

#define TCBR			(2.726				) /* CBR temperature, from COBE */

#define DSGRA			(8.178e3 * PC			) /* Distance from Earth to Sgr A*  */
#define DM87			(16.9e6 * PC			) /* Distance from Earth to M87  */
#define DABHB			(1.06e3 * PC			) /* Distance from Earth to A0620-00  */

#define MSGRA			(4.154e6 	                ) /* Mass of Sgr A* in MSUN */
#define MM87			(6.5e9   			) /* Mass of M87  */
#define MABHB			(20 	         		) /* Mass of A0620-00  or cygx1*/

/* 
   abundances, from M & B, p. 99
*/

#define SOLX			(0.70				) /* H */
#define SOLY			(0.28				) /* He */
#define SOLZ			(0.02				) /* Metals */
