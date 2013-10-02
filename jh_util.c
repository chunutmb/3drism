#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

/****************************************************************************************/
/*											*/
/*				Utilities						*/		
/*											*/		
/****************************************************************************************/ 

double * add_arrays( double *v1, double *v2, int n )
{
	int i;	double *v3 = (double *) malloc( n *sizeof(double));
	for( i=0; i<=n-1; i++) 	v3[i] = v1[i] + v2[i];
	return v3;
}

double * assign_0( int n )
{	int i; 
	double *array; array = (double *) malloc( n* sizeof(double ));
	for( i=0; i<=n-1; i++) array[i] = 0.0;
	return array;
}

double * assign_1( int n )
{	int i; 
	double *array; array = (double *) malloc( n* sizeof(double ));
	for( i=0; i<=n-1; i++) array[i] = 1.0;
	return array;
}

void assign_array( double *out, double *in, int n )
{
	int i;
	for( i=0; i<=n-1; i++)	out[i] = in[i];
}

void check_array( double *in, int n)
{
	int i; 
	double tmp;
	if( in == NULL)
		printf("\nError in array\n"); fflush(stdout);
	for( i=0; i<=n-1; i++) tmp = in[i]; 
}

double ** matrix_malloc( int nn, int mm)
{
	double **tmp;
	int i;

	tmp = malloc( nn * sizeof(  double));
	if( tmp == NULL) printf("\nallocation failed\n");

	for( i=0; i<= nn-1; i++){
		 *(tmp + i) = malloc(mm * sizeof(  double));
		if( *(tmp + i) == NULL);
				};
	return tmp;
}

double ** matrix_malloc0( int nn, int mm)
{
	double **tmp;
	int i,j;

	tmp = malloc( nn * sizeof(  double));
	if( tmp == NULL) printf("\nallocation failed\n");

	for( i=0; i<= nn-1; i++){
		 *(tmp + i) = malloc(mm * sizeof(  double));
		if( *(tmp + i) == NULL);
				};
	for( i=0; i<=nn-1; i++)
	    for( j=0; j<=mm-1; j++)
		tmp[i][j] = 0.00;

	return tmp;
}

double * sub_arrays( double *v1, double *v2, int n )
{
	int i;
	double *v3 = (double *) malloc( n *sizeof(double));
	for( i=0; i<= n-1; i++)
		v3[i] = v1[i] - v2[i];
	return v3;
}

double svector_norm(  double *vn, int nq)
{
	int i;
	double temp=0;
	for( i=0; i<=nq-1; i++)	    temp += fabs( vn[i]*vn[i]/nq );
	temp = pow(temp, 0.5);
	return temp;
}

double vector_norm(  double *vn, int nn)
{
	int i;
	double temp=0;

	for( i=0; i<=nn-1; i++)
	    temp += fabs( vn[i]*vn[i] );

	return pow(temp, 0.5);
}


double vector_sum( double *v, int n)
{
	int i;
	double tmp=0;
	
	for( i=0; i<=n-1; i++)
	    tmp += fabs( v[i] );

	return tmp;
}
	

