#include "jh_struct.h"
#include <fftw3.h>

#ifndef GET
 #define GET

   	double get_dval( char in_file[], char val[] );
   	char * get_ext( char in_file[] );
   	char * get_sval( char in_file[], char val[] );

	double  * get_array_dval( char in_file[], char val[], int n );
	double ** get_array_sval( char in_file[], char val[], int n );

	int get_tag( char in_file[], char val[] );
	int get_num_atoms( char [] );

	double * get_dis( char in_file[], int col );

	double * get_jh3d( char in_name[], double temp, double pnd, ENV_PAR sys );
	fftw_complex * get_complex_jh3d( char in_name[], double temp, double pnd, ENV_PAR sys );
	double * get_sit( char in_name[], ENV_PAR sys );
	fftw_complex * get_complex_sit( char in_name[], ENV_PAR sys );

  	U_PAR  * get_par(  int *,  char * );
  	U_PAR2 * get_par2( int *,  char * );
	U_PAR2 * get_par2_par( char * );

	int get_file_stat( char [] );

#endif


