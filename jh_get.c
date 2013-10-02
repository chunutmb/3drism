#include <stdio.h>
#include <string.h>
#include <fftw3.h>
#include "jh_struct.h"



double get_dval( char in_file[], char val[] )
{
	FILE *in;
		if( ( in = fopen( in_file, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_file);

	int idx;
	char s1[100], c1;
	double tmp;
	int test1=0;

	while( !feof(in) ){
		while( (fgetc(in) != '#')  && !feof(in) );
		c1 = fgetc(in);
		if( c1 != '#' )	{
			idx=0;
			do{	s1[idx] = c1;	/*get full variable name*/
				idx++;	
				c1 = fgetc(in);  
			} while( (c1 != ' ') && (c1 != '\n') && (c1 != '\t') && (c1 != '\r'));
			s1[idx] = '\0';
				
			if( strncmp( val, s1, idx ) == 0 ){	
						fscanf( in, "%lf", &tmp );	/*set_temp*/
						test1=1;
			} 
		}
		if( test1 == 1 ) break;
	}
	return tmp;
	fclose( in );
}

	/* get string (val) from file (in_file) */
char * get_sval( char in_file[], char val[] )
{
	FILE *in;
		if( ( in = fopen( in_file, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_file);

	int  idx;
	char s1[100], c1;
	char *s2 = (char *) malloc( 100 * sizeof(char));
	int test1=0;

	while( !feof(in) ){
		while( (fgetc(in) != '#')  && !feof(in) );
		c1 = fgetc(in);
		if( c1 != '#' )	{
			idx=0;
			do{	s1[idx] = c1;	/*get full variable name*/
				idx++;	
				c1 = fgetc(in);  
			} while( (c1 != ' ') && (c1 != '\n') && (c1 != '\t') && (c1 != '\r'));
			s1[idx] = '\0';
				
			if( strncmp( val, s1, idx ) == 0 ){	
						fscanf( in, "%s", s2 );	/*set_temp*/
						test1=1;
			} 
		}
		if( test1 == 1 ) break;
	}

	return s2;
	fclose( in );
}



double * get_array_dval( char in_file[], char val[], int n )
{
	FILE *in;
		if( ( in = fopen( in_file, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_file);

	int i, idx;
	char s1[100], c1;
	double *tmp = (double *) malloc( n * sizeof( double ) );
	int test1=0;

	while( !feof(in) ){
		while( (fgetc(in) != '#')  && !feof(in) );
		c1 = fgetc(in);
		if( c1 != '#' )	{
			idx=0;
			do{	s1[idx] = c1;	/*get full variable name*/
				idx++;	
				c1 = fgetc(in);  
			} while( (c1 != ' ') && (c1 != '\n') && (c1 != '\t') && (c1 != '\r'));
			s1[idx] = '\0';
				
			if( strncmp( val, s1, idx ) == 0 ){	
				for( i=0; i<=n-1; i++)
					fscanf( in, "%lf", (tmp+i) );	/*set_temp*/
						test1=1;
			} 
		}
		if( test1 == 1 ) break;
	}

	fclose( in );

	return tmp;
}



char ** get_array_sval( char in_file[], char val[], int n )
{
	FILE *in;
		if( ( in = fopen( in_file, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_file);

	int  idx, i;
	char s1[100], c1;
	char ** s2 = (char **) malloc( n*10 * sizeof(char));
		for( i=0; i<=n-1; i++)
			*(s2+i) = (char *) malloc( 100 * sizeof(char));
	int test1=0;

	while( !feof(in) ){
		while( (fgetc(in) != '#')  && !feof(in) );
		c1 = fgetc(in);
		if( c1 != '#' )	{
			idx=0;
			do{	s1[idx] = c1;	/*get full variable name*/
				idx++;	
				c1 = fgetc(in);  
			} while( (c1 != ' ') && (c1 != '\n') && (c1 != '\t') && (c1 != '\r'));
			s1[idx] = '\0';
				
			if( strncmp( val, s1, idx ) == 0 ){	
				for( i=0; i<=n-1; i++)
					fscanf( in, "%s", *(s2+i) );	/*set_temp*/
				test1=1;
			} 
		}
		if( test1 == 1 ) break;
	}

	return s2;
	fclose( in );
}

int get_num_atoms( char infile[] )
{
	int n;
	FILE *in;
		if( ( in = fopen( infile, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", infile );

	fscanf( in,  "%d", &n );

	return n;
}


int get_tag( char in_file[], char val[] )
{
	FILE *in;
		if( ( in = fopen( in_file, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_file);

	int idx;
	char s1[100], c1;
	int test1=0;
	int slen = (int) strlen( val );


	while( !feof(in) ){
		while( (fgetc(in) != '#')  && !feof(in) );
		c1 = fgetc(in);
		if( c1 != '#' )	{
			idx=0;
			do{	s1[idx] = c1;	/*get full variable name*/
				idx++;	
				c1 = fgetc(in);  
			} while( (c1 != ' ') && (c1 != '\n') && (c1 != '\t') && (c1 != '\r'));
			s1[idx] = '\0';
				
			if( strncmp( val, s1, idx ) == 0 ){	
						test1=1;
			}else 
			if( (strncmp( "DIS", s1, idx ) == 0) || (strncmp( "END", s1, idx) == 0 )){
				break;
			}
		}
		if( test1 == 1 ) break;
	}


	return test1;

	fclose( in );
}

char * get_ext( char infile[] )
{

	int i;
	int slen = (int) strlen( infile );
	int idx = slen -1;
	int xlen;

		while( infile[idx] != '.' )  
			idx--;
		idx++;
		xlen = slen - idx;
		char *ext = (char *) malloc( 2*(xlen+1)*sizeof( char ));

		for( i=0; i<=xlen-1; i++)
			ext[i] = infile[ idx+i ];

		ext[xlen] = '\0';

	return ext;

}

double * get_dis( char in_file[], int col )
{

	FILE *in;
		if( ( in = fopen( in_file, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_file );

	int i, j, i2, idx;
	int n_pts, dis_num;
	double radius, tmp;
	double *hr ;
	char s1[100], c1;


	while( !feof(in) ){

		while( fgetc(in) != '#' );
		
		c1 = fgetc(in);
					/*ungetc(char, stream)*/
		if( c1 != '#' )	{
			idx=0;
			/*get full variable name*/
			do{  
				s1[idx] = c1;	
				idx++;	
				c1 = fgetc(in);  
			} while( (c1 != ' ') && (c1 != '\n') && (c1 != '\t') && (c1 != '\r'));
			s1[idx] = '\0';
							
			if( strncmp( "N_PTS", s1, idx ) == 0 ){		fscanf( in, "%d", &n_pts );
			} else
			if( strncmp( "RADIUS", s1, idx ) == 0 ){	fscanf( in, "%lf", &radius );
			} else
			if( strncmp( "DIS_NUM", s1, 7 ) == 0){		fscanf( in, "%d", &dis_num );
			} else
			if( strncmp( "DIS", s1, idx ) == 0 ){	

						hr = (double *) malloc( n_pts *sizeof(double));
		
						for( j=0; j<=n_pts-1; j++){
							fscanf( in,  "%d", &i2 );
							for( i=1; i<=dis_num; i++){
								fscanf( in, "%lf", &tmp );
								if( i == col )
									hr[i2] = tmp;
							}
						} break;
			} 
		}		
				/*fscanf( in, "%[^#], );*/
	}

	fclose( in );

	return hr;
}




/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/*												*/
/*					get_par							*/
/*												*/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/




U_PAR * get_par( int *nu_sites, char *in_name )
{
	int n_atoms=0, n;
	FILE *in;
		if( ( in = fopen( in_name, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_name );

	fscanf( in, "%d", nu_sites );

	U_PAR *u = (U_PAR *) malloc( ((*nu_sites)+1)*sizeof(U_PAR));

	while( !feof( in ))/*Read in Parameters*/
	{
		fscanf( in,  "%d", &n ); if( n == n_atoms) break;
		u[n].num = n;
		fscanf( in,  "%s", u[n].element);
		fscanf( in, "%lf", &u[n].ep );
		fscanf( in, "%lf", &u[n].sig );
		fscanf( in, "%lf", &u[n].charge );
		fscanf( in, "%lf%lf%lf", &u[n].x, &u[n].y, &u[n].z );
		n_atoms = n;
	}	
	
	if( n_atoms != *nu_sites ){	
		fprintf( stdout, "::Problem with number of sites specified and number read in, %s", in_name);
		exit(1);
	}

	fclose( in );

	return u;
}





U_PAR2 * get_par2( int *nu_sites, char *in_name )
{
	int n_atoms=0, n;
	FILE *in;
		if( ( in = fopen( in_name, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", in_name );

	fscanf( in, "%d", nu_sites );

	U_PAR2 *u = (U_PAR2 *) malloc( ((*nu_sites)+1) *sizeof(U_PAR2));

	while( !feof( in ))/*Read in Parameters*/
	{
		fscanf( in,  "%d", &n ); if( n == n_atoms) break;
		u[n].num = n;
		fscanf( in,  "%s", u[n].element);
		fscanf( in,  "%d", &u[n].mol);
		fscanf( in, "%lf", &u[n].ep12 );
		fscanf( in, "%lf", &u[n].ep6 );
		fscanf( in, "%lf", &u[n].sig );
		fscanf( in, "%lf", &u[n].charge );
		fscanf( in, "%lf%lf%lf", &u[n].x, &u[n].y, &u[n].z );
		n_atoms = n;
	}	
	
	if( n_atoms != *nu_sites ){	
		fprintf( stdout, "::Problem with number of sites specified and number read in, %s", in_name);
		exit(1);
				  }
	fclose( in );
	return u;
}




U_PAR2 * get_par2_par( char *infile )
{
	int n_atoms=0, n;
	int nu_sites;
	FILE *in;
		if( ( in = fopen( infile, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", infile );

	fscanf( in, "%d", &nu_sites );

	U_PAR2 *u = (U_PAR2 *) malloc( ((nu_sites)+1)*sizeof(U_PAR2));

	while( !feof( in ))/*Read in Parameters*/
	{
		fscanf( in,  "%d", &n ); if( n == n_atoms) break;
		u[n].num = n;
		fscanf( in,  "%s", u[n].element);
		u[n].mol = 1;
		fscanf( in, "%lf", &u[n].ep12 );
		u[n].ep6 = u[n].ep12;
		fscanf( in, "%lf", &u[n].sig );
		fscanf( in, "%lf", &u[n].charge );
		fscanf( in, "%lf%lf%lf", &u[n].x, &u[n].y, &u[n].z );
		n_atoms = n;
	}	
	
	if( n_atoms != nu_sites ){	
		fprintf( stdout, "::Problem with number of sites specified and number read in, %s", infile);
		exit(1);
	}

	fclose( in );

	return u;
}






double * get_jh3d( char in_name[], double temp, double pnd, ENV_PAR sys )
{
	/*GLOBAL*/
	double lx=sys.lx, ly=sys.ly, lz=sys.lz;
	int nx2=sys.nx, ny2=sys.ny, nz2=sys.nz;
	/*global*/

	int x, y, z, nx, ny, nz;
	int status=0;
	double lx2, ly2, lz2, temp2, pnd2;

	FILE *in;
		if( ( in = fopen( in_name, "r" )) == NULL )
			fprintf( stderr, "\nCan't open file:%s\n", in_name );

	/*Begin reading in jh3d.dat file*/
	fscanf( in, "%d%d%d"   , &nx, &ny, &nz );
	fscanf( in, "%lf%lf%lf", &lx2, &ly2, &lz2 );
	fscanf( in, "%lf", &temp2 );
	fscanf( in, "%lf", &pnd2 );

		if( (nx2!=nx) || (ny2!=ny) || (nz2!=nz) )	status++;
		if( (lx2!=lx) || (ly2!=ly) || (lz2!=lz) )	status++;
		if( temp2 != temp )	fprintf( stdout, "\n:: %s input temperature and runtime temperature are not the same\n", in_name);
		if( pnd2 != pnd )	fprintf( stdout, "\n:: %s input pnd and runtime pnd are not the same\n", in_name);
		if( status != 0 ) { fprintf( stdout, "\n:: %s input system values and runtime system values are not the same: exiting\n", in_name);
			    	    exit(1); }

	double *vec_3d = (double *) malloc( nx*ny*nz *sizeof(double));

	while( !feof( in )) /*Read in Parameters*/
	{
		fscanf( in,  "%d%d%d", &x, &y, &z );
		fscanf( in, "%lf", &vec_3d[ nz*ny*x + nz*y + z ] );
	}	
	return vec_3d;
}

fftw_complex * get_complex_jh3d( char in_name[], double temp, double pnd, ENV_PAR sys )
{
	/*GLOBAL*/
	double lx=sys.lx, ly=sys.ly, lz=sys.lz;
	int nx2=sys.nx, ny2=sys.ny, nz2=sys.nz;
	/*global*/

	int x, y, z, nx, ny, nz, id;
	int status=0;
	double lx2, ly2, lz2, temp2, pnd2;

	FILE *in;
		if( ( in = fopen( in_name, "r" )) == NULL )
			fprintf( stderr, "\nCan't open file:%s\n", in_name );

	/*Begin reading in jh3d.dat file*/
	fscanf( in, "%d%d%d"   , &nx, &ny, &nz );
	fscanf( in, "%lf%lf%lf", &lx2, &ly2, &lz2 );
	fscanf( in, "%lf", &temp2 );
	fscanf( in, "%lf", &pnd2 );

		if( (nx2!=nx) || (ny2!=ny) || (nz2!=nz) )	status++;
		if( (lx2!=lx) || (ly2!=ly) || (lz2!=lz) )	status++;
		if( temp2 != temp )	status++;
		if( pnd2 != pnd ) 	status++;
		if( status != 0 ) { fprintf( stdout, "\n::%s input values and runtime values are not the same\n", in_name);
				    exit(1);}

	fftw_complex *jh_3d = (fftw_complex *) malloc( nx*ny*nz *sizeof(fftw_complex));

	while( !feof( in ))  /*Read in Parameters*/
	{
		fscanf( in,  "%d%d%d", &x, &y, &z );
		id = nz*ny*x +nz*y + z;
		fscanf( in, "%lf%lf", &jh_3d[ id ][0], &jh_3d[ id ][1] );
	}	
	return jh_3d;
}



double * get_sit( char in_name[], ENV_PAR sys )
{

	int x, y, z, bx, by, bz, nx, ny, nz;
	int status=0, idx;
	double vox, pnd2;

	FILE *in;
		if( ( in = fopen( in_name, "r" )) == NULL )
			fprintf( stderr, "\nCan't open file:%s\n", in_name );

	/*Begin reading in sit header*/
	fscanf( in, "%lf", &vox );
	fscanf( in, "%lf%lf%lf", &bx, &by, &bz );
	fscanf( in, "%d%d%d"   , &nx, &ny, &nz );

		if( (sys.nx != nx) || (sys.ny != ny) || (sys.nz != nz) )	status++;
		if( status != 0 ) { fprintf( stdout, "\n:: %s input system values and runtime system values are not the same: exiting\n", in_name);
			    	    exit(1); }

	double *vec_3d = (double *) malloc( nx*ny*nz *sizeof(double));
	double *vec2_3d = (double *) malloc( nx*ny*nz *sizeof(double));

	idx=0;
	while( !feof( in ) ) /*Read in Parameters*/
	{
		if( fscanf( in, "%lf", &vec_3d[ idx ] ) == 1 )
		   idx++;
	}	

	if ( idx != (nx*ny*nz) )
	{
		fprintf( stdout, "\n:: %s, number of specified points and number of grid points don't match\n", in_name );
		exit(1);
	}

	for( x=0 ; x<=nx-1; x++ )
	    for( y=0; y<=ny-1; y++ )
		for( z=0; z<=nz-1; z++ )
			vec2_3d[ nz*ny*x + nz*y + z ] = vec_3d[ nx*ny*z + nx*y + x ];

	free( vec_3d );

	return vec2_3d;
}

fftw_complex * get_complex_sit( char in_name[], ENV_PAR sys )
{
	/*GLOBAL*/
	int x, y, z,bx, by, bz, nx, ny, nz, ix;
	int status=0, idx;
	double vox; 

	FILE *in;
		if( ( in = fopen( in_name, "r" )) == NULL )
			fprintf( stderr, "\nCan't open file:%s\n", in_name );

	/*Begin reading in jh3d.dat file*/
	fscanf( in, "%lf", &vox );
	fscanf( in, "%lf%lf%lf", &bx, &by, &bz );
	fscanf( in, "%d%d%d"   , &nx, &ny, &nz );

		if( (sys.nx != nx) || (sys.ny != ny) || (sys.nz != nz) )	status++;
		if( status != 0 ) { fprintf( stdout, "\n::%s input system values and runtime system values are not the same: exiting\n", in_name);
				    exit(1);}

	fftw_complex *jh_3d = (fftw_complex *) malloc( nx*ny*nz *sizeof(fftw_complex));
	fftw_complex *jh2_3d = (fftw_complex *) malloc( nx*ny*nz *sizeof(fftw_complex));

	idx=0;
	while( !feof( in ))  /*Read in Parameters*/
	{
		if( fscanf( in, "%lf%lf", &jh_3d[ idx ][0], &jh_3d[ idx ][1] ) == 2 )
		   idx++;
	}
	
	if ( idx != (nx*ny*nz) )
	{
		fprintf( stdout, "\n:: %s, number of specified points and number of grid points don't match\n", in_name );
		exit(1);
	}

	for( x=0 ; x<=nx-1; x++ )
	    for( y=0; y<=ny-1; y++ )
		for( z=0; z<=nz-1; z++ ){
			jh2_3d[ nz*ny*x + nz*y + z ][0] = jh_3d[ nx*ny*z + nx*y + x ][0];
			jh2_3d[ nz*ny*x + nz*y + z ][1] = jh_3d[ nx*ny*z + nx*y + x ][1];
		}

	free( jh_3d );

	return jh2_3d;

}







/*
void check_par( U_PAR *u, int n_sites )
{
	extern double LX, LY, LZ;
	extern double NX, NY, NZ;

	int i;
	double dx=LX/NX; 
	double dy=LY/NY;
	double dz=LZ/NZ;
	
	double x, y, z;

	FILE *out;
		if( ( out = fopen("solute_check.dat", "w")) == NULL)
			fprintf( stdout, "Problem opening out file for lj parameters"); fflush(stdout);

	for( i=1; i<=n_sites; i++)
	{
		fprintf( out, "%d:%s\n", u[i].num, u[i].element); fflush(out); 
		fprintf( out, "ep:%f\n", u[i].ep ); fflush(out);
		fprintf( out, "sig:%f\n", u[i].sig ); fflush(out);
		fprintf( out, "Coulomb Charge:%f\n", u[i].charge ); fflush( out );
		fprintf( out, "Cartesian Coordinates:\n"); fflush(out);
		fprintf( out, "%f\t%f\t%f\n\n", u[i].x, u[i].y, u[i].z); fflush(out);
	}

	for( i=1; i<=n_sites; i++)
	{
		x = u[i].x/dx; 
		y = u[i].y/dy;
		z = u[i].z/dz;
		
		if( x == 0.00 && y == 0.00 && z == 0.00 )
		{
		   printf("\n::ERROR - Element %d is on a grid point\n", i); fflush(stdout);  	
		   exit(1);
		}
	}

	fclose(out);
}
*/

int get_file_stat( char infile[] )
{
	int stat;	
	FILE *in;
		if( ( in = fopen( infile, "r")) == NULL )
			stat = 0;
		  else{	stat = 1;
			fclose( in );
		  }
	return stat;

}

