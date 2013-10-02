
#ifndef JH_STRUCT
  #define JH_STRUCT
/*	This file contains all the jh structures
 *
 *
*/
#include <stdlib.h>
#include <stddef.h>
#include <fftw3.h>


#define Pi 3.1415926535897932385
#define CONST_CLMB (167058.3436)
//#define CONST_CLMB (16539.46097)
#define CONST_JKM 8.316382
#define sqrt(a) pow( a, 0.5 )
#define IsZero(a) (fabs(a) < 1.0e20 * DBL_MIN )
#define DBL_MIN 1E-37



/*U_PAR*/
struct solute_param {	int num; char element[25]; 
				double ep; double sig; double charge;
				double x; double y; double z;
		       };
typedef struct solute_param U_PAR;

/*U_PAR*/
struct solute_param2 {	int num; char element[25]; int mol; 
				double ep12; double ep6; double sig; double charge;
				double x; double y; double z;
		       };
typedef struct solute_param2 U_PAR2;

/*ENV_PAR*/
struct environment_par { int nx; int ny; int nz;
			 int cx; int cy; int cz;
			double lx; double ly; double lz; 
			};
typedef struct environment_par ENV_PAR;

#endif

