#include "jh_struct.h"



void print_jh3d( char name[], double *v, ENV_PAR sys, double temp, double pnd )
{
	int x, y, z;
	int nx=sys.nx;
	int ny=sys.ny;
	int nz=sys.nz;

	FILE *out;
	if ( ( out = fopen( name, "w" )) == NULL)
		printf( "\nFile could not be opened\n");

	fprintf( out, "%d\n%d\n%d\n", sys.nx, sys.ny, sys.nz ); fflush(out);
	fprintf( out, "%.10f\n%.10f\n%.10f\n", sys.lx, sys.ly, sys.lz ); fflush(out);
	fprintf( out, "%.10f\n", temp ); fflush(out);
	fprintf( out, "%.10f\n", pnd ); fflush(out);

	for( x=0; x<=nx-1; x++){
	    for( y=0; y<=ny-1; y++){
	      	for( z=0; z<=nz-1; z++)
			fprintf( out, "%d\t%d\t%d\t%.15e\n", x, y, z, v[ nz*ny*x + nz*y + z ]);
	    	fprintf( out, "\n");}
	    fflush(out);}
	fclose( out );
}



void print_sit( char name[], double *v, ENV_PAR sys )
{
	int x, y, z;
	int nx=sys.nx;
	int ny=sys.ny;
	int nz=sys.nz;
	FILE *out;
	if ( ( out = fopen( name, "w" )) == NULL)
		printf( "\nFile could not be opened\n");

	double dx = sys.lx/(nx-1);
	double dy = sys.ly/(ny-1);
	double dz = sys.lz/(nz-1);

	double fx = ( -1.0 * sys.cx) * dx;
	double fy = ( -1.0 * sys.cy) * dy;
	double fz = ( -1.0 * sys.cz) * dz;

	fprintf( out, "%.15f\n", dx ); fflush(out);
	fprintf( out, "%.10f\t%.10f\t%.10f\n", fx, fy, fz ); fflush(out);
	fprintf( out, "%d\t%d\t%d\n", nx, ny, nz ); fflush(out);

	for( z=0; z<=nz-1; z++)
	    for( y=0; y<=ny-1; y++)
		for( x=0; x<=nx-1; x++)
			fprintf( out, "%.15e\n", v[ nz*ny*x + nz*y + z ]);

	fflush( out );
	fclose( out );
}



void print_cmplx_jh3d( char name[], fftw_complex *v, ENV_PAR sys, double temp, double pnd )
{
	int x, y, z;
	int nx=sys.nx;
	int ny=sys.ny;
	int nz=sys.nz;

	FILE *out;
	if ( ( out = fopen( name, "w" )) == NULL)
		printf( "\nFile could not be opened\n");

	fprintf( out, "%d\n%d\n%d\n", nx, ny, nz ); fflush(out);
	fprintf( out, "%.10f\n%.10f\n%.10f\n", sys.lx, sys.ly, sys.lz ); fflush(out);
	fprintf( out, "%.10f\n", temp ); fflush(out);
	fprintf( out, "%.10f\n", pnd ); fflush(out);

	for( x=0; x<=nx-1; x++){
	    for( y=0; y<=ny-1; y++){
	      	for( z=0; z<=nz-1; z++)
			fprintf( out, "%d\t%d\t%d\t%.14e\t%.14e\n", x, y, z, v[ nz*ny*x + nz*y + z ][0], v[ nz*ny*x + nz*y + z][1]);
	    	fprintf( out, "\n");}
	    fflush(out);}
	fclose( out );
}



void print_cmplx_sit( char name[], fftw_complex *v, ENV_PAR sys )
{
	int x, y, z;
	int nx=sys.nx;
	int ny=sys.ny;
	int nz=sys.nz;

	FILE *out;
	if ( ( out = fopen( name, "w" )) == NULL)
		printf( "\nFile could not be opened\n");

	double dx = sys.lx/(nx-1);
	double dy = sys.ly/(ny-1);
	double dz = sys.lz/(nz-1);

	double fx = ( -1.0 * sys.cx) * dx;
	double fy = ( -1.0 * sys.cy) * dy;
	double fz = ( -1.0 * sys.cz) * dz;

	fprintf( out, "%.15f\n", dx ); 				fflush(out);
	fprintf( out, "%.10f\t%.10f\t%.10f\n", fx, fy, fz ); 	fflush(out);
	fprintf( out, "%d\t%d\t%d\n", nx, ny, nz ); 		fflush(out);

  	for( z=0; z<=nz-1; z++)
	    for( y=0; y<=ny-1; y++)
		for( x=0; x<=nx-1; x++)
			fprintf( out, "%.14e\t%.14e\n", v[ nz*ny*x + nz*y + z ][0], v[ nz*ny*x + nz*y + z][1]);

	fflush( out );
	fclose( out );
}




/*
void print_par1d( U_PAR2 *u, int n_sites )
{
	int i;
	FILE *out;
		if( ( out = fopen("solute.dat", "w")) == NULL)
			fprintf( stdout, "Problem opening out file for lj parameters"); fflush(stdout);

	for( i=1; i<=n_sites; i++)
	{
		fprintf( out, "%d:%s\n", u[i].num, u[i].element); fflush(out);
		fprintf( out, "mol:%d\n", u[i].mol ); fflush(out);
		fprintf( out, "ep12:%f\n", u[i].ep12 ); fflush(out);
		fprintf( out, "ep6:%f\n", u[i].ep6 ); fflush(out);
		fprintf( out, "sig:%f\n", u[i].sig ); fflush(out);
		fprintf( out, "Coulomb Charge:%f\n", u[i].charge ); fflush( out );
		fprintf( out, "Cartesian Coordinates:\n"); fflush(out);
		fprintf( out, "%f\t%f\t%f\n\n", u[i].x, u[i].y, u[i].z); fflush(out);
	}
	fclose(out);
}
*/


#ifdef MPI
/***
void print_jh3d_mp( char name[], double *v, int nx, int ny, int nz, double pnd )
{
	int x, y, z;
	double l[3]={LX,LY,LZ};
	FILE *out;

	if( my_rank == 0 ){
		if ( ( out = fopen( name, "w" )) == NULL)
			printf( "\nFile could not be opened\n");
		fprintf( out, "%d\n%d\n%d\n", nx, ny, nz ); fflush(out);
		fprintf( out, "%.10f\n%.10f\n%.10f\n", l[0], l[1], l[2] ); fflush(out);
		fprintf( out, "%.10f\n", TEMP ); fflush(out);
		fprintf( out, "%.10f\n", pnd ); fflush(out);
		fclose( out );
	}

	for( i=0; i<=NP-1; i++){
	    if( RANK == i ){

		if ( ( out = fopen( name, "a" )) == NULL) printf( "\nFile could not be opened\n");

		for( x=0; x<=nx-1; x++){
	    	  for( y=0; y<=ny-1; y++){
	      	    for( z=0; z<=nz-1; z++)
			fprintf( out, "%d\t%d\t%d\t%.14e\n", x, y, z, v[ nz*ny*x + nz*y + z ]);
	    	    fprintf( out, "\n");}
	          fflush(out);}
		fclose( out );
	    }
	    MPI_Barrier( MPI_COMM_WORLD );
	}

}
*******/
#endif

