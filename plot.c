/****************************************************************************************/
/*This code does all the analysis of the input 3d gr file
 *
 * 1st arg:	Solute input file (.par)
 * 2nd arg: 	Solvent 3d gr (.3dgr)
 *
*/

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "jh_struct.h"
#include "jh_get.h"

#define A_ERF 1.08

#define ii(x,y,z) ( NZ*NY*(x) + NZ*(y) + (z) )
#define is(x,y,z) ( NX*NY*(z) + NX*(y) + (x) )

/*Global variables - Environment variables read in, a change of these must be announced */
int NX, NY, NZ, NNN;
int CX, CY, CZ;
int N_1D;
double LX, LY, LZ;
double DX, DY, DZ;
double FX, FY, FZ;
double P_D, TEMP;

double *gr_3d;

U_PAR2 *U;
int NU_SITES=0;
int PAR_STAT=0;

/************************************************************************/
/*				Subroutines				*/
/************************************************************************/
void get_solute_parameters(  char * );
void check_parameters( void );
void get_3d( char []);	/*name of file, n_pts, r_rad, column */

/*1*/ void print_2d_file( void );
/*2*/ void print_3d_to_1d_radial_average( void );
/*3*/ void print_3d_to_1d_cartesian_average( void );
/*4*/ void print_1d_unaveraged( void );
/*5*/ void print_box_jh3d( char [] );
/*6*/ void print_sit_file( void );
/*7*/ void print_peak_list( void );
/*8*/ void print_altered_data_set( void );
/*9*/ void print_2d_avg_sit( void );


void print_input_error( void );

double r0( int, int, int );			/*from origin  x, y, z*/
double rx( int, int, int, double, double, double );/* x, y, z, ui_x, ui_y, ui_z*/

double * add_arrays( double *, double *, int );

/*when executing program 1st arg solute file, 2nd arg is quick plotting of a jh3d file hr_vv*/
/****************************************************************************************/
/*											*/
/*					Main						*/
/*											*/
/****************************************************************************************/

int main(int argc, char *argv[])
{

	int i, j, k, l, x, y, z, cnt=0, *dum;
	int choice;
	char *parfile, jh3dfile;
	char *ext1;

	/*****************************GET SOLUTE PARAMETERS*****************************************/
							printf("\nReading from files...\n");fflush(stdout);
	if( *(argv+1) == NULL )
		print_input_error();

	ext1 = get_ext( *(argv+1) );

	if( (strncmp( "par", ext1, 3 ) == 0) && (strlen( ext1 ) == 3) )
	{
		U = get_par2_par( *(argv+1) );
		NU_SITES = get_num_atoms( *(argv+1) );

		if( *(argv+2) == NULL )
			print_input_error();
		else
			get_3d( *(argv+2) );

	} else
		if( (strncmp( "par2", ext1, 4 ) == 0) && (strlen( ext1 ) == 4) )
	 	{
			U = get_par2( dum, *(argv+1) );
			NU_SITES = get_num_atoms( *(argv+1) );

			if( *(argv+2) == NULL )
				print_input_error();
			else
				get_3d( *(argv+2) );

		} else {
			/* No par file in input */
			PAR_STAT = 1;
			get_3d( *(argv+1) );

	   	}

							printf("...done!!!\n\n");fflush(stdout);


	/*****************************KERNEL*********************************************/


	do{
		printf("\nSelect a choice:\n");
		printf("0: EXIT\n");
		printf("1: Print a 2D file (slice)(quick feat.)\n");
		printf("2: Print a 3D to 1D radial average file (quick feat.)\n");
		printf("3: Print a 3D to 1D cartesian average (.jh3d)\n");
		printf("4: Print a 1D unaveraged distribution\n");
		printf("5: Print a box.jh3d \n");
		printf("6: Print a sit file \n");
		printf("7: Print list of peaks in function \n" );
		printf("8: Print an altered data set \n" );
		printf("9: Print an averaged 2d sit file \n" );
		printf("\nENTER CHOICE:"); fflush(stdout);
		scanf("%d", &choice);
		printf("\nPerforming operation...\n"); fflush(stdout);

		switch( choice ){
				case 0:	break;
				case 1: print_2d_file(); 			break;
				case 2: print_3d_to_1d_radial_average();	break;
				case 3: print_3d_to_1d_cartesian_average();	break;
				case 4: print_1d_unaveraged();			break;
				case 5: print_box_jh3d(*(argv+1));		break;
				case 6: print_sit_file();			break;
				case 7: print_peak_list();			break;
				case 8: print_altered_data_set();		break;
				case 9: print_2d_avg_sit();			break;
				}

		printf("...DONE!!!\n"); fflush(stdout);

	} while( choice != 0 );

	return 0;
}

/****************************************************************************************/
/****************************************************************************************/
/**************************************END OF MAIN***************************************/
/****************************************************************************************/
/****************************************************************************************/


/****************************************************************************************/
/*											*/
/*				1:Print 2d slice					*/
/*											*/
/****************************************************************************************/

void print_2d_file( void )
{
	int i, x, y, z, axis, axval;
	int n[3]={NX,NY,NZ};
	char fname[50];

	printf("\nEnter name of output file:");
	scanf( "%s", fname );
	printf("\nEnter invariant axis:\n1: X-axis\n2: Y-axis\n3: Z-axis\nSelection:");fflush(stdout);
	scanf("%d", &axis);
	printf("\nEnter invariant axis value (%d->%d, Center->%d):", 0, n[axis-1], n[axis-1]/2); fflush( stdout );
	scanf("%d", &axval);

	printf("\nPrinting 2D file..."); fflush(stdout);

	printf( "\n%d\t%d\t%d\n", NX,NY,NZ ); fflush(stdout);

	FILE *out;
	if ( ( out = fopen( fname, "w" )) == NULL)   printf( "\nFile could not be opened\n");

	switch( axis ){
			case 1:
				for( y=0; y<=NY-1; y++){
	    			    for( z=0; z<=NZ-1; z++)
					fprintf( out, "%d\t%d\t%f\n", y, z, gr_3d[ ii(axval,y,z) ]);
	    			    fprintf( out, "\n");
				      			}
				break;
			case 2:
				for( x=0; x<=NX-1; x++){
	    			    for( z=0; z<=NZ-1; z++)
					fprintf( out, "%d\t%d\t%f\n", x, z, gr_3d[ ii(x,axval,z) ]);
	    			    fprintf( out, "\n");
				      			}
				break;
			case 3:
				for( x=0; x<=NX-1; x++){
	    			    for( y=0; y<=NY-1; y++)
					fprintf( out, "%d\t%d\t%f\n", x, y, gr_3d[ ii(x,y,axval) ]);
	    			    fprintf( out, "\n");
				      			}
				break;
			}
	fflush(out);
	fclose(out);

	printf("...done!!!\n"); fflush(stdout);

}

/****************************************************************************************/
/*											*/
/*				2: 3D to 1D radial averaging				*/
/*											*/
/****************************************************************************************/

void print_3d_to_1d_radial_average( void )
{

	int n[3]={NX,NY,NZ};
	double l[3]={LX,LY,LZ};

	int n_1d, i, id, x, y, z;
	int nu_sites = NU_SITES;
	double r, r1d, dr1d, lmax, x0, y0, z0;
	char fname[50];

	printf("\nEnter name of output file:");		scanf( "%s", fname );
	printf("\nEnter # of points to use on 1D grid:");	scanf( "%d", &n_1d);
	printf("\n\nSolute sites:");
	if( PAR_STAT == 0 )
	   for( i=1; i<=nu_sites; i++)
		printf("\n%s\tx:%f\ty:%f\tz:%f", U[i].element, U[i].x, U[i].y, U[i].z );
	else
	    printf( "\nNo par(2) file entered");
	printf("\nEnter X(0) coordinate:");		scanf("%lf", &x0 );
	printf("\nEnter Y(0) coordinate:");		scanf("%lf", &y0 );
	printf("\nEnter Z(0) coordinate:");		scanf("%lf", &z0 );
	printf("\nPrinting 1d radial average..."); fflush(stdout);

	int    *grnum = (int *) malloc( n_1d *sizeof(int));
	double *gr_1d = (double *) malloc( n_1d *sizeof(double));

	for( i=0; i<= n_1d-1; i++)  grnum[i] = 0;
	for( i=0; i<= n_1d-1; i++)  gr_1d[i]  = 0.00;

	lmax = ( l[0]*l[0] + l[1]*l[1] + l[2]*l[2] )/4;
	r1d = sqrt(lmax);
	dr1d = (double) r1d/( (double) (n_1d-1));

	for( x=0; x<=n[0]-1; x++)
	 for( y=0; y<=n[1]-1; y++)
	  for( z=0; z<=n[2]-1; z++)
		{
			r = rx( x, y, z, x0, y0, z0 );
			r += (dr1d/2.0);
			id = (int) (r/dr1d) ;
			gr_1d[id] += gr_3d[ ii(x,y,z) ];
			grnum[id]++;
		}

	for( i=0; i<= n_1d-1; i++)
		if( grnum[i] == 0 )  gr_1d[i] = 0.00;
		   else              gr_1d[i] = gr_1d[i] / ( (double) grnum[i]);

	FILE *out;
	  if ( ( out = fopen( fname, "w" )) == NULL)   printf( "\nFile could not be opened\n");

	for( i=0; i<=n_1d-1; i++)
		fprintf( out, "%f\t%f\n", dr1d *i, gr_1d[i]);

	fflush(out);

	printf("...done!!!\n"); fflush(stdout);
}




/****************************************************************************************/
/*											*/
/*				3: 3D to 1D cartesian averaging				*/
/*											*/
/****************************************************************************************/


void print_3d_to_1d_cartesian_average( void )
{
	int n[3]={NX,NY,NZ}, c[3]={CX,CY,CZ};
	double l[3]={LX,LY,LZ}, d[3]={DX,DY,DZ};

	int i, x, y, z, axis, chc, chc1;
	int imin[3] = { 0, 0, 0 };
	int imax[3] = { NX-1, NY-1, NZ-1 };
	double min[3];
			for( i=0; i<=2; i++)
				min[i] = (-c[i]*d[i] );
	double max[3];
			for( i=0; i<=2; i++)
				max[i] = (n[i]-1-c[i])*d[i];
	double r_x, r_y, r_z;
	char ax[]="xyz";

	double *gr_1d = (double *) malloc( NZ * sizeof(double));
	int *gr_num = ( int *) malloc( NZ * sizeof(int ));


		printf("\nEnter the axis to be plotted:  (1:x, 2:y, 3:z) : ");		scanf( "%d", &axis );
		axis--;

		printf("\n(1) use full data set or (2) specify reduced coordinate range : "); scanf( "%d", &chc );

		if( chc == 1 )
		{
				/*nothing to be done*/
		} else
		 if( chc == 2)
		 {
			printf("\nEnter values in Angstroms:\n\n");
			for( i=0; i<=2; i++ )
			    if( i != chc )
			    {
				printf("Enter %c-axis minimum (>%f): ", ax[i], min[i] );
				scanf("%lf", &min[i] );
			    }
			printf("\n");
			for( i=0; i<=2; i++ )
			    if( i != chc )
			    {
				printf("Enter %c-axis maximum (<%f): ", ax[i], max[i] );
				scanf("%lf", &max[i] );
			    }

			for( i=0; i<=2; i++)
			    if( i!= chc ){
				imin[i] = (int) (min[i]/d[i]) + c[i];
				imax[i] = (int) (max[i]/d[i]) + c[i];
			    }


		 } else
		  	printf("\nError reading user input.\n\n");


		printf("\nCalculating 1d distribution..."); 	fflush(stdout);

		if( axis == 0 )
		{
			for( x=0; x<=NX-1; x++)
			{
				gr_1d[x] = 0.00;
				gr_num[x] = 0;
				for( y=imin[1]; y<=imax[1]; y++ )
				 for( z=imin[2]; z<=imax[2]; z++ )
				 {
					gr_1d[x] += gr_3d[ ii(x,y,z) ];
					gr_num[x]++;
				 }
			}
		} else
		 if( axis == 1 )
		 {
			for( y=0; y<=NY-1; y++)
			{
				gr_1d[y] = 0.00;
				gr_num[y] = 0;
				for( x=imin[0]; x<=imax[0]; x++ )
				 for( z=imin[2]; z<=imax[2]; z++ )
				 {
					gr_1d[y] += gr_3d[ ii(x,y,z) ];
					gr_num[y]++;
				 }
			}


		 } else
		  if( axis == 2)
		  {
			for( z=0; z<=NZ-1; z++)
			{
				gr_1d[z] = 0.00;
				gr_num[z] = 0;
				for( x=imin[0]; x<=imax[0]; x++ )
				 for( y=imin[1]; y<=imax[1]; y++ )
				 {
					gr_1d[z] += gr_3d[ ii(x,y,z) ];
					gr_num[z]++;
				 }
			}

		  } else {
				printf( "\nError reading axis input\n\n" ); fflush(stdout);
		   	 }

		for( i=imin[axis]; i<=imax[axis]; i++)
			gr_1d[i] = gr_1d[i] / gr_num[i];

		printf("...Done!\n\n"); fflush(stdout);


		printf("Printing to gr1d_%c-avg.dat...", ax[axis] );

	char fname[50];
		sprintf( fname, "gr1d_%d-avg.dat", ax[axis] );

	FILE *out;
	  if ( ( out = fopen( fname, "w" )) == NULL)   printf( "\nFile could not be opened\n");

		for( i=imin[axis]; i<=imax[axis]; i++)
			fprintf( out, "%e\t%e\n", (i-c[axis]) *d[axis], gr_1d[i]);

		printf("...Done!!!\n\n");	fflush(stdout);

	fflush(out);
	fclose(out);

}


/****************************************************************************************/
/*											*/
/*				4: 1D unaveraged distribution 	 			*/
/*											*/
/****************************************************************************************/


void print_1d_unaveraged( void )
{

	int axis, x, y, z, i, j;
	char fname[50];
	FILE *fout;

	printf(" What axis is to be plotted (1->x, 2->y, 3->z): "); fflush( stdout );
	scanf("%d", &axis );

	if ( axis == 1 ){

		printf("\nEnter y value (0->%d): ", NY-1); fflush(stdout);
		scanf("%d", &y );
		printf("\nEnter z value (0->%d): ", NZ-1); fflush(stdout);
		scanf("%d", &z );

		sprintf( fname, "gr-_x-%dy-%dz.1d.dat", y, z );
		if( (fout=fopen( fname, "w")) == NULL )
			printf("\nProblem opening file\n"); fflush(stdout);

		for( i=0; i<=NX-1; i++)
		    fprintf( fout, "%d\t%lf\n", i, gr_3d[ ii(i,y,z)] );

	} else
	 if ( axis == 2 ){

		printf("\nEnter x value (0->%d): ", NX-1); fflush(stdout);
		scanf("%d", &x );
		printf("\nEnter z value (0->%d): ", NZ-1); fflush(stdout);
		scanf("%d", &z );

		sprintf( fname, "gr-%dx-_y-%dz.1d.dat", x, z );
		if( (fout=fopen( fname, "w")) == NULL )
			printf("\nProblem opening file\n"); fflush(stdout);

		for( i=0; i<=NY-1; i++)
		    fprintf( fout, "%d\t%lf\n", i, gr_3d[ ii(x,i,z)] );


	 } else
	  if ( axis == 3 ){

		printf("\nEnter x value (0->%d): ", NX-1); fflush(stdout);
		scanf("%d", &x );
		printf("\nEnter y value (0->%d): ", NY-1); fflush(stdout);
		scanf("%d", &y );

		sprintf( fname, "gr1d-%dx-%dy-_z.dat", x, y );
		if( (fout=fopen( fname, "w")) == NULL )
			printf("\nProblem opening file\n"); fflush(stdout);

		for( i=0; i<=NZ-1; i++)
		    fprintf( fout, "%d\t%lf\n", i, gr_3d[ ii(x,y,i)] );


	  } else { printf("\nERROR choosing axis\n\n" ); fflush(stdout); }


	fclose( fout );

}


/****************************************************************************************/
/*											*/
/*				5: print box_*.jh3d  	 			*/
/*											*/
/****************************************************************************************/


void print_box_jh3d( char fname[] )
{

	int x, y, z, stat;
	char s1[100];
	FILE *out;

	printf("\nWhat half of the box ( 1:1st, 2:2nd ) : ");		scanf( "%d", &stat );

	sprintf( s1, "box_%s", fname );

	if ( ( out = fopen( s1, "w" )) == NULL)
		printf( "\nFile could not be opened\n");

	fprintf( out, "%d\n%d\n%d\n", NX, NY, NZ/2 ); fflush(out);
	fprintf( out, "%.10f\n%.10f\n%.10f\n", LX, LY, LZ/2.0 ); fflush(out);
	fprintf( out, "%.10f\n", TEMP ); fflush(out);
	fprintf( out, "%.10f\n", P_D ); fflush(out);

	if ( stat == 1 ){

		for( x=0; x<=NX-1; x++){
		    for( y=0; y<=NY-1; y++){
	      		for( z=0; z<=CZ-1; z++)
				fprintf( out, "%d\t%d\t%d\t%.15e\n", x, y, z-CZ, gr_3d[ NZ*NY*x + NZ*y + z ]);
		    	fprintf( out, "\n");}
		    fflush(out);}

	} else
	 if ( stat == 2 ){

		for( x=0; x<=NX-1; x++){
		    for( y=0; y<=NZ-1; y++){
	      		for( z=CZ; z<=NZ-1; z++)
				fprintf( out, "%d\t%d\t%d\t%.15e\n", x, y, z-CZ, gr_3d[ NZ*NY*x + NZ*y + z ]);
		    	fprintf( out, "\n");}
		    fflush(out);}
	 }

	fclose( out );

}

/****************************************************************************************/
/*											*/
/*				6: print .sit file	  	 			*/
/*											*/
/****************************************************************************************/

void print_sit_file( void )
{

	int x, y, z, stat=0;
	int ix, iy, iz, lx, ly, lz;
	int nx, ny, nz;
	double fx, fy, fz;

	char s1[]="out.sit";

	FILE *fout;
		if( (fout = fopen( s1, "w" )) == NULL )
			printf( "\nfile could not be opened\n" );

		printf("\nPrint (1) full data set or (2) subset of data set w/ new dimensions\n");
		printf("Enter choice: ");
		scanf("%d", &stat );

		if( stat == 1 )
		{
			fprintf( fout, "%.10f\n", DX );
			fprintf( fout, "%.10f\t%.10f\t%.10f\n", FX, FY, FZ);
			fprintf( fout, "%d\t%d\t%d\n", NX, NY, NZ);
			for( z=0; z<=NZ-1; z++)
			 for( y=0; y<=NY-1; y++)
			  for( x=0; x<=NX-1; x++){
				fprintf( fout, "%.14e\n", gr_3d[ ii(x,y,z) ] );
			  }

		} else
		   if( stat == 2 )
		   {
			printf("\n");
			printf("Enter starting x grid point (%d-%d) : ", 0, NX-1 ); scanf( "%d", &ix );
			printf("Enter starting y grid point (%d-%d) : ", 0, NY-1 ); scanf( "%d", &iy );
			printf("Enter starting z grid point (%d-%d) : ", 0, NZ-1 ); scanf( "%d", &iz );
			printf("\n");
			printf("Enter final x grid point (%d-%d) : ", ix, NX-1); scanf( "%d", &lx );
			printf("Enter final y grid point (%d-%d) : ", iy, NY-1); scanf( "%d", &ly );
			printf("Enter final z grid point (%d-%d) : ", iz, NZ-1); scanf( "%d", &lz );

			fx = (ix-CX) * DX;
			fy = (iy-CY) * DY;
			fz = (iz-CZ) * DZ;

			nx = lx-ix+1;
			ny = ly-iy+1;
			nz = lz-iz+1;

			fprintf( fout, "%.10f\n", DX );
			fprintf( fout, "%.10f\t%.10f\t%.10f\n", fx, fy, fz );
			fprintf( fout, "%d\t%d\t%d\n", nx, ny, nz);
			for( z=iz; z<=lz; z++)
			 for( y=iy; y<=ly; y++)
			  for( x=ix; x<=lx; x++){
				fprintf( fout, "%.14e\n", gr_3d[ ii(x,y,z) ] );
			  }

		   } else {
				printf("\nError reading user input\n\n");
			  }

	fflush( fout );
	fclose( fout );

}



/****************************************************************************************/
/*											*/
/*				7: print list of peaks					*/
/*											*/
/****************************************************************************************/

void print_peak_list( void )
{

	int i, id=0, x, y, z;
	int chc;
	double r, rr, dr1d, x0=0.0, y0=0.0, z0=0.0;
	char fname[50];

	printf("\nEnter name of output file: ");					scanf( "%s", fname );
	FILE *out;
	  if ( ( out = fopen( fname, "w" )) == NULL)   printf( "\nFile could not be opened\n");

	printf("\nDefine a reference point other than box center (0-no, 1-yes): ");	scanf( "%d", &chc );

	if( chc == 1 ) {
		printf("\nEnter coordinates in angstroms");
		printf("\nEnter X(0) coordinate:");		scanf("%lf", &x0 );
		printf("\nEnter Y(0) coordinate:");		scanf("%lf", &y0 );
		printf("\nEnter Z(0) coordinate:");		scanf("%lf", &z0 );
	}

	fprintf(out,"%-10s%10s%10s%10s%10s%10s", "index", "g(r)", "r  ", "x  ", "y  ", "z  " );
	if( chc == 1 )
		fprintf( out, "%10s%10s%10s%10s", "r-ref", "x-ref", "y-ref", "z-ref"  );
	fprintf( out, "%10s%10s%10s", "x-index", "y-index", "z-index"  );

	fprintf( out, "\n" );

	for( x=1; x<=NX-2; x++)
	 for( y=1; y<=NY-2; y++)
	  for( z=1; z<=NZ-2; z++)
	  {
		if( gr_3d[ii(x,y,z)] > gr_3d[ii(x-1,y,z)] )
		if( gr_3d[ii(x,y,z)] > gr_3d[ii(x+1,y,z)] )
		if( gr_3d[ii(x,y,z)] > gr_3d[ii(x,y-1,z)] )
		if( gr_3d[ii(x,y,z)] > gr_3d[ii(x,y+1,z)] )
		if( gr_3d[ii(x,y,z)] > gr_3d[ii(x,y,z-1)] )
		if( gr_3d[ii(x,y,z)] > gr_3d[ii(x,y,z+1)] )
		{
			id++;
			fprintf( out, "%-10d", id );
			fprintf( out, "%10.3f", gr_3d[ ii(x,y,z) ] );
			fprintf( out, "%10.3f", r0(x,y,z ) );
			fprintf( out, "%10.3f%10.3f%10.3f", r0(x,CY,CZ), r0(CX,y,CZ), r0(CX,CY,z) );
			if( chc == 1 )
			{
				fprintf( out, "%10.3f", rx( x, y, z, x0, y0, z0 ) );
				fprintf( out, "%10.3f", ((x-CX)*DX - x0) );
				fprintf( out, "%10.3f", ((y-CY)*DY - y0) );
				fprintf( out, "%10.3f", ((z-CZ)*DZ - z0) );
			}
			fprintf( out, "%10d%10d%10d", x, y, z  );

			fprintf( out, "\n" );
		}
	   }

	fflush(out);

	printf("...done!!!\n"); fflush(stdout);
}



/****************************************************************************************/
/*											*/
/*				8: print altered data set  	 			*/
/*											*/
/****************************************************************************************/

void print_altered_data_set( void )
{

	int x, y, z, stat=0, schc, chc1;
	int ix, iy, iz, lx, ly, lz;
	int nx, ny, nz;
	double rix, riy, riz, rs;
	double xmin, ymin, zmin, xmax, ymax, zmax, rho, ngr, rad;

	char s1[]="altered_data_set.sit";

	FILE *fout;
		if( (fout = fopen( s1, "w" )) == NULL )
			printf( "\nfile could not be opened\n" );

		printf("\nEnter shape of region to change (0-box, 1-circle, 2-cylinder(z) ) : ");
		scanf("%d", &stat );
		printf("\n");
		printf("\nChange values (1) inside or (2) outside: " );
		scanf("%d", &chc1 );
		printf("\n");
		if( chc1 == 1 )
			printf("\nChange values inside of region to: ");
		if( chc1 == 2 )
			printf("\nChange values outside of region to: ");
		scanf("%lf", &ngr );
		printf("\n");

		if( stat == 0 )
		{
			printf("\nEnter coordinates of center of circle\n");
			printf("\nEnter x-min: "); 	scanf("%lf", &xmin );
			printf("\nEnter x-max: "); 	scanf("%lf", &xmax );
			printf("\nEnter y-min: "); 	scanf("%lf", &ymin );
			printf("\nEnter y-max: "); 	scanf("%lf", &ymax );
			printf("\nEnter z-min: "); 	scanf("%lf", &zmin );
			printf("\nEnter z-max: "); 	scanf("%lf", &zmax );

			for( z=0; z<=NZ-1; z++)
			 for( y=0; y<=NY-1; y++)
			  for( x=0; x<=NX-1; x++)
			  {
				rix = (x-CX)*DX;
				riy = (y-CY)*DY;
				riz = (z-CZ)*DZ;

				if( chc1 == 1)
				{
					if( rix > xmin && rix < xmax )
					if( riy > ymin && riy < ymax )
					if( riz > zmin && riz < zmax )
						gr_3d[ ii(x,y,z) ] = ngr;

				} else
				 if( chc1 == 2)
				 {
					if( rix < xmin || rix > xmax )
						gr_3d[ ii(x,y,z) ] = ngr;
					if( riy < ymin || riy > ymax )
						gr_3d[ ii(x,y,z) ] = ngr;
					if( riz < zmin || riz > zmax )
						gr_3d[ ii(x,y,z) ] = ngr;
				 }

			  }

		} else
		 if( stat == 1 )
		 {
			printf("\nEnter coordinates of center of circle\n");
			printf("\nEnter x: "); 		scanf("%lf", &rix );
			printf("\nEnter y: "); 		scanf("%lf", &riy );
			printf("\nEnter z: "); 		scanf("%lf", &riz );
			printf("\nEnter radius: ");	scanf("%lf", &rad );

			for( z=0; z<=NZ-1; z++)
			 for( y=0; y<=NY-1; y++)
			  for( x=0; x<=NX-1; x++)
			  {
				rs = rx(x,y,z,rix,riy,riz);

				if( chc1 == 1)
				{
					if( rs < rad )
						gr_3d[ ii(x,y,z) ] = ngr;
				}else
				 if( chc1 == 2)
				 {
					if( rs > rad )
						gr_3d[ ii(x,y,z) ] = ngr;
				 }
			  }
		 } else
		  if( stat == 2 )
		  {
			printf("\nEnter zmin: "); 	scanf("%lf", &zmin );
			printf("\nEnter zmax: "); 	scanf("%lf", &zmax );
			printf("\nEnter rho: ");	scanf("%lf", &rho );

			for( z=0; z<=NZ-1; z++)
			 for( y=0; y<=NY-1; y++)
			  for( x=0; x<=NX-1; x++)
			  {
				riz = (double) (z-CZ)*DZ;
				rs = r0(x,y,CZ);

				if( chc1 == 1)
				{
					if( riz > zmin && riz < zmax )
					 if( rs < rho )
						gr_3d[ ii(x,y,z) ] = ngr;
				} else
				 if( chc1 == 2 )
				 {
					if( riz < zmin || riz > zmax )
						gr_3d[ ii(x,y,z) ] = ngr;
					 if( rs > rho )
						gr_3d[ ii(x,y,z) ] = ngr;
				 }

			  }
		  }

			printf("\n\nPrinting to altered_data_set.sit\n\n"); fflush(stdout);

			fprintf( fout, "%.10f\n", DX );
			fprintf( fout, "%.10f\t%.10f\t%.10f\n", FX, FY, FZ);
			fprintf( fout, "%d\t%d\t%d\n", NX, NY, NZ);
			for( z=0; z<=NZ-1; z++)
			 for( y=0; y<=NY-1; y++)
			  for( x=0; x<=NX-1; x++){
				fprintf( fout, "%.14e\n", gr_3d[ ii(x,y,z) ] );
			  }


	fflush( fout );
	fclose( fout );

}


/****************************************************************************************/
/*											*/
/*				9:Print 2d avgeraged slice				*/
/*											*/
/****************************************************************************************/

void print_2d_avg_sit( void )
{
	int i, x, y, z, axis, axval, nslc;
	double div;
	int n[3]={NX,NY,NZ};
	double gval;
	char fname[50];

	printf("\nEnter name of output sit file:");
	scanf( "%s", fname );
	printf("\nEnter invariant axis:\n1: X-axis\n2: Y-axis\n3: Z-axis\nSelection:");fflush(stdout);
	scanf("%d", &axis);
	printf("\nEnter invariant axis value (%d->%d, Center->%d):", 0, n[axis-1], n[axis-1]/2); fflush( stdout );
	scanf("%d", &axval);
	printf("\nEnter number of slices +- to use: ");
	scanf("%d", &nslc);

	printf("\nPrinting 2D average file..."); fflush(stdout);

	printf( "\n%d\t%d\t%d\n", NX,NY,NZ ); fflush(stdout);

	FILE *out;
	if ( ( out = fopen( fname, "w" )) == NULL)   printf( "\nFile could not be opened\n");

	div = 2*nslc + 1;

	fprintf( out, "%f\n", DX );

	switch( axis ){
			case 1:
				fprintf( out, "%f\t%f\t%f\n", 0.00, -CY*DY, -CZ*DZ );
				fprintf( out, "%d\t%d\t%d\n", 1, NY, NZ);
	    			for( z=0; z<=NZ-1; z++)
				{
				  for( y=0; y<=NY-1; y++)
				  {
					gval = 0.00;
					for( x=axval-nslc; x<=axval+nslc; x++ )
						gval += gr_3d[ ii( x, y, z )] ;
					gval /= div;

					//fprintf( out, "%d\t%d\t%f\n", y, z, gval );
					fprintf( out, "%f\n", gval );
				  }
				}
				break;
			case 2:
				fprintf( out, "%f\t%f\t%f\n", -CX*DX, -CY*DY, 0.00 );
				fprintf( out, "%d\t%d\t%d\n", NX, 1, NZ);
	    			for( z=0; z<=NZ-1; z++)
				{
				  for( x=0; x<=NX-1; x++)
				  {
					gval = 0.00;
					for( y=axval-nslc; y<=axval+nslc; y++ )
						gval += gr_3d[ ii( x, y, z) ];
					gval /= div;

					//fprintf( out, "%d\t%d\t%f\n", x, z, gval );
					fprintf( out, "%f\n", gval );
				  }
				}
				break;
			case 3:
				fprintf( out, "%f\t%f\t%f\n", -CX*DX, -CY*DY, 0.00 );
				fprintf( out, "%d\t%d\t%d\n", NX, NY, 1);
	    			for( y=0; y<=NY-1; y++)
				{
				  for( x=0; x<=NX-1; x++)
				  {
					gval = 0.00;
					for( z=axval-nslc; z<=axval+nslc; z++ )
						gval += gr_3d[ ii( x, y, z) ];
					gval /= div;

					//fprintf( out, "%d\t%d\t%f\n", x, y, gval );
					fprintf( out, "%f\n", gval );
				  }
				}
				break;
			}
	fflush(out);
	fclose(out);

	printf("...done!!!\n"); fflush(stdout);

}


/****************************************************************************************/
/*											*/
/*				Get 3dgr for quick plots				*/
/*											*/
/****************************************************************************************/


void get_3d( char fname[]  )
{
	FILE *fin;
		if( ( fin = fopen( fname, "r")) == NULL )
			fprintf( stdout, "::Problem opening file containing solUte parameters:%s\n", fname );

	int x, y, z;
	int bx, by, bz, nx, ny, nz;
	int status=0, idx;
	double vox, pnd2;

	char *ext;
	ext = get_ext( fname );
	printf( "%s", ext);
	if ( (strncmp( "jh3d", ext, 4 ) == 0) )
	{
		fscanf( fin, "%d%d%d"   , &NX, &NY, &NZ );
		fscanf( fin, "%lf%lf%lf", &LX, &LY, &LZ );
		fscanf( fin, "%lf", &TEMP );
		fscanf( fin, "%lf", &P_D );

		gr_3d = (double *) malloc( NX*NY*NZ *sizeof(double));

		while( !feof( fin ))//Read in Parameters
		{
			fscanf( fin,  "%d", &x );
			fscanf( fin,  "%d", &y );
			fscanf( fin,  "%d", &z );
			fscanf( fin, "%lf", &gr_3d[ ii(x,y,z) ] );
		}

		DX = (LX/(NX-1));
		DY = (LY/(NY-1));
		DZ = (LZ/(NZ-1));

		if( (DX != DY) || (DY != DZ) || (DX != DZ) )
			printf("\nWarning!!! : sit files are not compatible with non-symmetric voxel dimensions\n");

		CX=(NX/2);
		CY=(NY/2);
		CZ=(NZ/2);

		FX = -1.0 *CX *DX;
		FY = -1.0 *CY *DY;
		FZ = -1.0 *CZ *DZ;

	} else
	 if( strncmp( "sit", ext, 3 ) == 0 )
	 {

		fscanf( fin, "%lf", &vox );
		fscanf( fin, "%lf%lf%lf", &FX, &FY, &FZ );
		fscanf( fin, "%d%d%d"   , &NX, &NY, &NZ );

		gr_3d = (double *) malloc( NX*NY*NZ *sizeof(double));
		double *gr2_3d = (double *) malloc( NX*NY*NZ *sizeof(double));

		idx=0;
		while( !feof( fin ) ) /*Read in Parameters*/
		{
			if( fscanf( fin, "%lf", &gr2_3d[ idx ] ) == 1 )
			   idx++;
		}

		if ( idx != (NX*NY*NZ) )
		{
			fprintf( stdout, "\n:: %s, number of specified points and number of grid points don't match\n", fname );
			exit(1);
		}


		for( x=0 ; x<=NX-1; x++ )
		    for( y=0; y<=NY-1; y++ )
			for( z=0; z<=NZ-1; z++ )
				gr_3d[ NZ*NY*x + NZ*y + z ] = gr2_3d[ NX*NY*z + NX*y + x ];

		free( gr2_3d );

		DX = vox; /*(-2.0*FX/(NX));	*/
		DY = vox; /*(-2.0*FY/(NY));	*/
		DZ = vox; /*(-2.0*FZ/(NZ));	*/

		LX = DX *(NX-1);
		LY = DY *(NY-1);
		LZ = DZ *(NZ-1);

		CX=(NX/2);
		CY=(NY/2);
		CZ=(NZ/2);

	 } else
		fprintf( stdout, "\n::File extension did not match compatible types\n"); fflush( stdout );

		NNN = NX*NY*NZ;

}



/****************************************************************************************/
/*					r(x,y,z)					*/
/*				Calculate the spatial distances				*/
/*											*/
/****************************************************************************************/

double rx( int x, int y, int z, double ux, double uy, double uz)
{
	double tmp;

	tmp  = pow(fabs( x*DX - (CX*DX + ux) ), 2) ;
	tmp += pow(fabs( y*DY - (CY*DY + uy) ), 2) ;
	tmp += pow(fabs( z*DZ - (CZ*DY + uz) ), 2) ;
	tmp = sqrt( tmp );

	return tmp;
}


double r0( int x, int y, int z)
{
	double tmp;

	tmp  = pow(fabs(x - CX),2) * pow(DX,2);
	tmp += pow(fabs(y - CY),2) * pow(DY,2);
	tmp += pow(fabs(z - CZ),2) * pow(DZ,2);
	tmp = sqrt( tmp );

	return tmp;
}


double * add_arrays( double *v1, double *v2, int n )
{
	int i;
	double *v3 = (double *) malloc( n *sizeof(double));
	for( i=0; i<=n-1; i++)
		v3[i] = v1[i] + v2[i];
	return v3;
}




/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				Unused subroutines					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/




void print_1d( char [], double *, int, double );
void print_2d( char [], double *, int );
void print_3d( char [], double * );
double k0( int, int, int );
double kx( int, int, int, double, double, double );
double ** matrix_malloc( int , int );
double * sub_arrays( double *, double *, int );
void assign_array( double *, double *, int); /*out, in, n*/
double * assign_1( int );
double * assign_0( int );
void check_array( double *, int );


double kx( int x, int y, int z, double ux, double uy, double uz)
{

	double tmp;
	double dkx = 2*Pi/LX;
	double dky = 2*Pi/LY;
	double dkz = 2*Pi/LZ;

	tmp  = pow(fabs( x - (CX + (ux/DX)) ),2) *pow(dkx,2);
	tmp += pow(fabs( y - (CY + (uy/DY)) ),2) *pow(dky,2);
	tmp += pow(fabs( z - (CZ + (uz/DZ)) ),2) *pow(dkz,2);
	tmp = sqrt( tmp );
	return tmp;
}

double k0( int x, int y, int z )
{
	double tmp;

	double dkx = 2*Pi/LX;
	double dky = 2*Pi/LY;
	double dkz = 2*Pi/LZ;

	tmp  = pow(fabs(x - CX),2) * pow(dkx,2);
	tmp += pow(fabs(y - CY),2) * pow(dky,2);
	tmp += pow(fabs(z - CZ),2) * pow(dkz,2);
	tmp = sqrt( tmp );

	return tmp;
}


/****************************************************************************************/
/*											*/
/*				Utilities						*/
/*											*/
/****************************************************************************************/

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

void assign_array( double *out, double *in, int n )
{
	int i;
	for( i=0; i<=n-1; i++)	out[i] = in[i];
}
double * assign_1( int n )
{	int i;
	double *array; array = (double *) malloc( n* sizeof(double ));
	for( i=0; i<=n-1; i++) array[i] = 1.0;
	return array;
}
double * assign_0( int n )
{	int i;
	double *array; array = (double *) malloc( n* sizeof(double ));
	for( i=0; i<=n-1; i++) array[i] = 0.0;
	return array;
}
double * sub_arrays( double *v1, double *v2, int n )
{
	int i;
	double *v3 = (double *) malloc( n *sizeof(double));
	for( i=0; i<= n-1; i++)
		v3[i] = v1[i] - v2[i];
	return v3;
}

void check_array( double *in, int n)
{
	int i;
	double tmp;
	if( in == NULL)
		printf("\nError in array\n"); fflush(stdout);
	for( i=0; i<=n-1; i++) tmp = in[i];
}


/****************************************************************************************/
/*											*/
/*				Printing Routines					*/
/*											*/
/****************************************************************************************/


void print_input_error( void )
{
	printf("\nMust enter one of the following combinations of files on command line:\n");
	printf(" 1) par(2) (*.jh3d/*.sit)\n");
	printf(" 2) (*.jh3d/sit)\n" );
	exit(1);
}


void print_1d( char nameq[], double * xq, int nq, double dq)
{
	int x;
	FILE *out;
	if ( ( out = fopen( nameq, "w" )) == NULL)   printf( "\nFile could not be opened\n");
	for( x=0; x<=nq-1; x++)   fprintf( out, "%f\t%f\n", (double) x*dq, xq[x]);
	fflush(out);
	fclose( out );
}

void print_2d( char nameq[], double *xq, int u_zq )
{
	int x, y, Nx=NX, Ny=NY, Nz=NZ;
	FILE *out;
	if ( ( out = fopen( nameq, "w" )) == NULL)   printf( "\nFile could not be opened\n");
	for( x=0; x<=Nx-1; x++){
	    for( y=0; y<=Ny-1; y++)   fprintf( out, "%d\t%d\t%f\n", x, y, xq[ii(x,y,u_zq)]);
	    fprintf( out, "\n");
	    fflush(out);    	}
	fclose( out );
}


void print_3d( char nameq[], double *xq )
{
	int x, y, z;
	int Nx=NX, Ny=NY, Nz=NZ;
	FILE *out;
	if ( ( out = fopen( nameq, "w" )) == NULL)
		printf( "\nFile could not be opened\n");
	for( x=0; x<=Nx-1; x++){
	    for( y=0; y<=Ny-1; y++){
	      	for( z=0; z<=Nz-1; z++)
			fprintf( out, "%d\t%d\t%d\t%f\n", x, y, z, xq[ ii(x,y,z)]);
	    	fprintf( out, "\n");}
	    fflush(out);}
	fclose( out );
}



