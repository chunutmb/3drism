
#ifndef JH_PRINT
  #define JH_PRINT

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>


void print_jh3d(char [], double *, ENV_PAR, double, double);     /*file_name, array, system val, temp, pnd*/
void print_sit(char [], double *, ENV_PAR);    /*file_name, array, system val*/

//void print_jh3d_mp( char [], double * , int, int, int, double ); /*file_name, array, nx, ny, nz, pnd*/

void print_cmplx_jh3d(char [], fftw_complex *, ENV_PAR, double, double);    /*file_name, array, system val, temp, pnd*/
void print_cmplx_sit(char [], fftw_complex *, ENV_PAR);    /*file_name, array, nx, ny, nz, pnd*/
#endif


