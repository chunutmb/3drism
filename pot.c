/********************************************************************************/
/*										*/
/*			Calculate potentials					*/
/*										*/
/********************************************************************************/
/*										*/
/*			Author:	Jesse Howard                                    */
/*			Date:	Jan 29, 2007					*/
/*										*/
/********************************************************************************/
/* This routine calculates the potential fields for the model input (U_par) in
 * the solvent model specified below
 *
 * This routine utilizes MPI routines
 * This routine is independent of all libraries
 *
 *
 *
*/

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <string.h>
#include "jh_get.h"
#include "jh_grid.h"
#include "jh_struct.h"
#include "jh_linalg.h"
#include "jh_util.h"
#include "jh_print.h"

/*#define MPI*/
#ifdef MPI
 #include <mpi.h>
#endif

/*#define WCA 0.0*/

/*env - Box parameters*/
double LX, LY, LZ;
int NX, NY, NZ;
int CX, CY, CZ;
char *EWALD_SUMS, *CLOSURE, *CONFIG_TYPE, *FILE_TYPE;
char *BRIDGE_FUNC0, *RBC_FUNC;
double GAUSS_SMEAR;
ENV_PAR SYS;


/*dis - solvent properties*/
int NSITES, NRSITES, TYPE, DIS_NUM;
double TEMP, A_ERF=1.08, *PND, *REDUN;
char **NAMES, **DIS_NAMES;
double *EP12, *EP6, *SIG, *CHARGE, *BOND_MAT;


/*par - solute parameters*/
int NU_SITES;
U_PAR2 *U;
int PAR_TYPE;

/*************************/

#define NNN NX*NY*NZ
#define NYZ NY*NZ

#define ii(x,y,z) ( NZ*NY*(x) + NZ*(y) + (z) )
#define ix(x) ( NZ*NY*(x) )

/*********************************************/

void set_env( char [] );
void set_dis( char [] );
void set_par( char [] );
void set_sys( void );

void check_dis( void );
void check_par( void );
void check_env( void );

void print_3d( char [], double *, int, int, int, double );
void print_cmplx_3d( char [], fftw_complex *, int, int, int, double); /*file_name, array, nx, ny, nz, pnd*/

/*_____________U(r) and U(k) subroutines_______________*/
void calc_and_print_ur_lj( U_PAR2 * , int, double , double, double, char [], int );	 /*solute parameters, n_u_sites, ep, sig*/
void calc_and_print_ur_lj12( U_PAR2 * , int, double , double, double, char [], int );	 /*solute parameters, n_u_sites, ep, sig*/
void calc_and_print_ur_clmb( U_PAR2 * , int, double, char [] , int );  	/* u parameters, nu_sites, z*/
void calc_and_print_ur_erf( U_PAR2 * , int, double, char [] , int );   			/* ur_clmb*/
void calc_and_print_uk_erf( U_PAR2 * , int, double, char [] , int );
void calc_and_print_ur_wca( U_PAR2 * , int, double , double, double, char [], int  );	 /*solute parameters, n_u_sites, ep, sig*/
void calc_and_print_ur_clmb_ewald( U_PAR2 * , int, double, char [] , int );  	/* u parameters, nu_sites, z*/
void calc_and_print_ur_clmb_ewald_rad( U_PAR2 * , int, double, char [] , int );  	/* u parameters, nu_sites, z*/
void calc_and_print_ewald_corrections( void );


/*______________FFTW routines________________*/
void fftw_3d( fftw_complex * , fftw_complex *);		/*in_r, out_k*/
void invfftw_3d( fftw_complex *, fftw_complex *);	/*in_k, out_r*/


/*_______Spatial distance_____*/
double r0( int, int, int );			/*from origin  x, y, z*/
double rx( int, int, int, double, double, double );/* x, y, z, ui_x, ui_y, ui_z*/
double rw( int, double );/* z, ui_z*/
double k0( int, int, int );
double k2( int, int, int );
double kk( int, int, int, double, double, double );


int RANK, NP;
char *ARGV1, *ARGV2, *ARGV3;

double DX, DY, DZ;
double DKX, DKY, DKZ;
double DKX2, DKY2, DKZ2;
/* 1st arg .dis, 2nd arg .env, 3rd arg .par */

int main( int argc, char *argv[] )
{
		set_env( *(argv +2) );
		set_dis( *(argv +1) );
		set_sys();

		ARGV1 = *(argv +1);
		ARGV2 = *(argv +2);
		ARGV3 = *(argv +3);


		int i;
		char s1[100];
		int my_rank=0, n_procs=1;

#ifdef MPI
		MPI_Init( &argc, &argv );
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &n_procs );
 #endif
		RANK = my_rank;
		NP = n_procs;

		if( my_rank == 0 ){
		   check_dis();
		   check_env();
		}
		#ifdef MPI
			MPI_Barrier( MPI_COMM_WORLD );
		 #endif

		/****Reading in Solute site parameters******************************/

		for( i=0; i<=n_procs-1; i++ ){
		    if( my_rank == i)
		    {
			printf("Processor-%d reading in solute parameters from: \"%s\" of type ", my_rank, ARGV3 );fflush(stdout);
			set_par( *(argv +3) );
			printf("Node %d, ...Done!\n", my_rank ); fflush(stdout);
		    }

		    #ifdef MPI
			MPI_Barrier( MPI_COMM_WORLD );
		     #endif
		}

		int nu_sites = NU_SITES;

		if( my_rank == 0 )
		   	printf("\n"); fflush(stdout);

		/*checking solute.par*/
		for( i=0; i<=n_procs-1; i++ ){
			if( my_rank == i )
			{
				printf("%d:# solute sites - %d ...", my_rank, NU_SITES );fflush(stdout);
				check_par( );
			}

			#ifdef MPI
				MPI_Barrier( MPI_COMM_WORLD );
			 #endif
		}

		/***********Calculating Energy Fields****************************************/
		if( my_rank == 0 )
		   	printf("\nCalculating potential fields with %d processors...\n", n_procs ); fflush(stdout);

		/********LENNARD-JONES term*****************/
		if( my_rank == 0 )
			printf("...u(r)_lj...");  fflush(stdout);

		for( i=0; i<=NRSITES-1; i++){
			sprintf( s1, "ur_%s_lj", NAMES[i] );

			#ifndef WCA
				calc_and_print_ur_lj( U, nu_sites, EP12[i], EP6[i], SIG[i], s1, i );
			#else
				calc_and_print_ur_wca( U, nu_sites, EP12[i], EP6[i], SIG[i], s1, i );
			#endif
		}

				#ifdef MPI
				 MPI_Barrier( MPI_COMM_WORLD );
 				  #endif

		/***** UR_lj-12 ******/
		if( (strncmp( "yes0", BRIDGE_FUNC0,4) == 0) || (strncmp( "yes", RBC_FUNC, 3) == 0)){

			if( my_rank == 0)
			   printf( "\n...u(r)_lj-12...");  fflush(stdout);

			for( i=0; i<=NRSITES-1; i++){
				sprintf( s1, "ur_%s_lj12", NAMES[i] );

				calc_and_print_ur_lj12( U, nu_sites, EP12[i], EP6[i], SIG[i], s1, i );
			}
		}

				#ifdef MPI
				 MPI_Barrier( MPI_COMM_WORLD );
				  #endif

		/***** NG's LONG RANGE METHOD ******/
		if( strncmp( "no", EWALD_SUMS, 2) == 0 ){

				/*****COULOMB term*****/
				if( my_rank == 0)
					printf("\n...u(r)_clmb...");  fflush(stdout);
				for( i=0; i<=NRSITES-1; i++){
					sprintf( s1, "ur_%s_clmb", NAMES[i] );
					calc_and_print_ur_clmb( U,  nu_sites, CHARGE[i], s1, i );
				}

				/*****long range (r) term*****/
				if( my_rank == 0)
					printf("\n...u(r)_l...");  fflush(stdout);
				for( i=0; i<=NRSITES-1; i++){
					sprintf( s1, "ur_%s_l", NAMES[i] );
					calc_and_print_ur_erf( U, nu_sites, CHARGE[i], s1 , i);
				}

				/*****long range (k) term*****/
				if( my_rank == 0)
					printf("\n...u(k)_l...");  fflush(stdout);
				for( i=0; i<=NRSITES-1; i++){
					sprintf( s1, "uk_%s_l", NAMES[i] );
					calc_and_print_uk_erf( U, nu_sites, CHARGE[i], s1 , i);
				}
		}else
		/***** EWALD LONG RANGE METHOD ******/
		 if( strncmp( "yes", EWALD_SUMS, 3) == 0) {

				/*Coulomb-Ewald term*/
				if( my_rank == 0)
					printf("\n...u(r)_clmb_ewald...");  fflush(stdout);
				for( i=0; i<=NRSITES-1; i++){
					sprintf( s1, "ur_%s_clmb", NAMES[i] );
					calc_and_print_ur_clmb_ewald( U,  nu_sites, CHARGE[i], s1, i );
				}

				if( my_rank == 0)
					printf("\n...u(r)_ewald_corrections...");  fflush(stdout);
					if( my_rank == 0 )
						calc_and_print_ewald_corrections();

		 /**** ERROR with Ewald ***/
		 }else
		  if( strncmp( "rad", EWALD_SUMS, 3) == 0 ) {

				/*Coulomb-Ewald term*/
				if( my_rank == 0)
					printf("\n...u(r)_clmb_ewald...");  fflush(stdout);
				for( i=0; i<=NRSITES-1; i++){
					sprintf( s1, "ur_%s_clmb", NAMES[i] );
					calc_and_print_ur_clmb_ewald_rad( U,  nu_sites, CHARGE[i], s1, i );
				}

				if( my_rank == 0)
					printf("\n...u(r)_ewald_corrections...");  fflush(stdout);
					if( my_rank == 0 )
						calc_and_print_ewald_corrections();

		  }else{
			printf("\nNo value for EWALD_SUMS was detected!!!\n");
			exit(1);
		      }
#ifdef MPI
		MPI_Barrier( MPI_COMM_WORLD );
 #endif


		if( my_rank == 0)
			printf("\n.....DONE!\n\n"); fflush(stdout);

#ifdef MPI
	MPI_Finalize();
 #endif

	return 0;
}




/****************************************************************************************/
/*					I/O						*/
/*				GET Parameter and solvent				*/
/*											*/
/****************************************************************************************/

void set_dis( char infile[] )
{
	TEMP 	= (double) get_dval( infile, "TEMP");
	TYPE	= (int) get_dval( infile, "TYPE" );
	NSITES 	= (int)  get_dval( infile, "NSITES");
	NRSITES = (int)  get_dval( infile, "NRSITES");
	DIS_NUM = (int ) get_dval( infile, "DIS_NUM" );
	REDUN  	= (double *) get_array_dval( infile, "REDUN", NRSITES );
	PND  	= (double *) get_array_dval( infile, "PND", NRSITES );
	EP12 	= (double *) get_array_dval( infile, "EP12", NRSITES );
	EP6 	= (double *) get_array_dval( infile, "EP6", NRSITES );
	SIG 	= (double *) get_array_dval( infile, "SIG", NRSITES );
	CHARGE 	= (double *) get_array_dval( infile, "CHARGE", NRSITES );
	BOND_MAT= (double *) get_array_dval( infile, "BOND_MAT", NRSITES );
	NAMES 	  = (char **) get_array_sval( infile, "NAMES", NRSITES);
	DIS_NAMES = (char **) get_array_sval( infile, "DIS_NAMES", NRSITES);


}

void check_dis( void )
{
	int i;
			printf("\n****************************************\n");
			printf("*		DIS		       *\n");
			printf("****************************************\n");
			printf("\nTEMP = %lf", TEMP ); fflush(stdout);
			printf("\nTYPE = %d", TYPE ); fflush(stdout);
			printf("\nNSITES = %d", NSITES ); fflush(stdout);
			printf("\nNRSITES = %d", NRSITES ); fflush(stdout);
			printf("\nDIS_NUM = %d", DIS_NUM ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nREDUN_%s = %lf", NAMES[i], REDUN[i] ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nPND_%s = %lf", NAMES[i], PND[i] ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nEP12_%s = %lf",NAMES[i], EP12[i] ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nEP6_%s = %lf",NAMES[i], EP6[i] ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nSIG_%s = %lf", NAMES[i], SIG[i] ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nCHARGE_%s = %lf",NAMES[i], CHARGE[i] ); fflush(stdout);
			for( i=0; i<=NRSITES-1; i++)
				printf("\nNAMES_%d = %s", i, NAMES[i]) ; fflush(stdout);
			printf("\n"); fflush(stdout);
			printf("****************************************\n\n"); fflush(stdout);

}


void set_env( char infile[] )
{
	NX = (int) get_dval( infile, "NX" );
	NY = (int) get_dval( infile, "NY" );
	NZ = (int) get_dval( infile, "NZ" );
	LX = (double) get_dval( infile, "LX" );
	LY = (double) get_dval( infile, "LY" );
	LZ = (double) get_dval( infile, "LZ" );
	CLOSURE = (char *) get_sval( infile, "CLOSURE");

	if( get_tag( infile, "CX" ) == 1)
		CX = (int) get_dval( infile, "CX" );
	   else CX = (int) NX / 2;

	if( get_tag( infile, "CY" ) == 1)
		CY = (int) get_dval( infile, "CY" );
	   else CY = (int) NY / 2;

	if( get_tag( infile, "CZ" ) == 1)
		CZ = (int) get_dval( infile, "CZ" );
	   else CZ = (int) NZ / 2;

	if( get_tag( infile, "A_ERF" ) == 1 )
		A_ERF = (double) get_dval( infile, "A_ERF");
	   else	A_ERF = 1.08;

	if( get_tag( infile, "GAUSS_SMEAR") == 1 )
		GAUSS_SMEAR = (double) get_dval( infile, "GAUSS_SMEAR");
	   else GAUSS_SMEAR = 0.65;

	if( get_tag( infile, "EWALD_SUMS") == 1 )
		EWALD_SUMS = (char *) get_sval( infile, "EWALD_SUMS" );
	   else EWALD_SUMS = "no";

	if( get_tag( infile, "CONFIG_TYPE") == 1 )
		CONFIG_TYPE = (char *) get_sval( infile, "CONFIG_TYPE" );
	   else CONFIG_TYPE = "none";

	if( get_tag( infile, "BRIDGE_FUNC0") == 1 )
		BRIDGE_FUNC0 = (char *) get_sval( infile, "BRIDGE_FUNC0");
	   else BRIDGE_FUNC0 = "no";

	if( get_tag( infile, "RBC_FUNC") == 1 )
		RBC_FUNC = (char *) get_sval( infile, "RBC_FUNC");
	   else RBC_FUNC = "no";

	if( get_tag( infile, "FILE_TYPE") == 1 )
		FILE_TYPE = (char *) get_sval( infile, "FILE_TYPE" );
	   else FILE_TYPE = "jh3d";

}


void check_env( void )
{
			printf("\n****************************************\n");
			printf("*		ENV		       *\n");
			printf("****************************************\n");
			printf("\nNX = %d\nNY = %d\nNZ = %d\n", NX, NY, NZ); 	fflush(stdout);
			printf("\nLX = %f\nLY = %f\nLZ = %f\n", LX, LY, LZ); 	fflush(stdout);
			printf("\nCX = %d\nCY = %d\nCZ = %d\n\n", CX, CY, CZ); 	fflush(stdout);

			printf("A_ERF = %f \n", A_ERF ); 			fflush(stdout);
			printf("EWALD_SUMS = %s \n", EWALD_SUMS  );		fflush(stdout);
			printf("CLOSURE = %s \n", CLOSURE  );			fflush(stdout);
			printf("GAUSS_SMEAR = %f \n", GAUSS_SMEAR );		fflush(stdout);
			printf("CONFIG_TYPE = %s \n", CONFIG_TYPE );		fflush(stdout);
			printf("BRIDGE_FUNC0 = %s \n", BRIDGE_FUNC0 );		fflush(stdout);
			printf("FILE_TYPE = %s \n", FILE_TYPE );		fflush(stdout);
			printf("RBC_FUNC = %s \n", RBC_FUNC );			fflush(stdout);
			printf("****************************************\n\n"); fflush(stdout);

}


void set_par( char infile[] )
{
		int *nu_s = (int *) malloc(sizeof(int));
		int i, idx=0, len;
		U_PAR *u1 = NULL;
		U_PAR2 *u2 = NULL;

		char vs;

		while( infile[idx] != '.' )  idx++;
		idx++;

		len = strlen( infile );
		vs = infile[len-1];

		if( vs == '2' ){
			printf( "par2" );
			PAR_TYPE = 2;
			u2 = (U_PAR2 *) get_par2( nu_s, infile );
		} else{
			printf( "par" );
			PAR_TYPE = 1;
			u1 = (U_PAR *) get_par( nu_s, infile );
		  }

		NU_SITES = *nu_s;

		if( PAR_TYPE == 1 ){

			U = (U_PAR2 *) malloc( (NU_SITES+1) * sizeof(U_PAR2));

			for( i=1; i<=NU_SITES; i++){
				U[i].num     = u1[i].num ;
			   	sprintf( U[i].element, "%s", u1[i].element ) ;
			   	U[i].mol     = 0 ;
				U[i].ep12    = u1[i].ep ;
			   	U[i].ep6     = u1[i].ep ;
			   	U[i].sig     = u1[i].sig ;
			   	U[i].charge  = u1[i].charge ;
		   		U[i].x 	     = u1[i].x ;
			   	U[i].y 	     = u1[i].y ;
			   	U[i].z 	     = u1[i].z ;
			}
		} else
 		 if( PAR_TYPE == 2)
			U = (U_PAR2 *) u2;
}


void check_par( void )
{
	int i;
	double dx=(double)LX/(NX-1);
	double dy=(double)LY/(NY-1);
	double dz=(double)LZ/(NZ-1);

	int ix, iy, iz, stat=0;
	double x, y, z;

		for( i=1; i<=NU_SITES; i++)
		{
			ix = U[i].x / dx;
			iy = U[i].y / dy;
			iz = U[i].z / dz;

			x = U[i].x - ix*dx;
			y = U[i].y - iy*dy;
			z = U[i].z - iz*dz;

			if( x == 0.00 && y == 0.00 && z == 0.00 )
			{
			   printf("\n::ERROR - Element %d is on a grid point\n", i); fflush(stdout);
			   /*exit(1);*/ stat = 1;
			}
		}

		if( stat != 1 )
			printf("...checked out OK!!!\n");fflush(stdout);

}

void set_sys( void )
{
		SYS.nx = NX;	SYS.ny = NY; 	SYS.nz = NZ;
		SYS.cx = CX;	SYS.cy = CY; 	SYS.cz = CZ;
		SYS.lx = LX;	SYS.ly = LY;	SYS.lz = LZ;

		DX=(double)LX/(NX-1);
		DY=(double)LY/(NY-1);
		DZ=(double)LZ/(NZ-1);

		DKX =(double) 2*Pi/LX;
		DKY =(double) 2*Pi/LY;
		DKZ =(double) 2*Pi/LZ;

		DKX2 = DKX;
		DKY2 = DKY;
//		if( strncmp( "wall", CONFIG_TYPE, 4) == 0 )
//			DKZ2 = 2.0 *DKZ;
//		 else
		DKZ2 = DKZ;

}


/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*					U(r) - LJ					*/
/*				calc LJ potential routines				*/
/*											*/
/****************************************************************************************/

void  calc_and_print_ur_lj( U_PAR2 * u, int nu_sites, double ep12, double ep6, double sig, char fname[], int ns )
{
	int i, x, y, z, indx;
	int nx=NX, ny=NY, nz=NZ;
	double r, t = TEMP;

	double ep12_ux[nu_sites +1], ep6_ux[nu_sites +1], sig_ux[nu_sites +1];
	for( i=1; i<=nu_sites; i++)	ep12_ux[i] = sqrt( u[i].ep12 * ep12 );
	for( i=1; i<=nu_sites; i++)	ep6_ux[i]  = sqrt( u[i].ep6 * ep6 );
	for( i=1; i<=nu_sites; i++)	sig_ux[i]  = (u[i].sig + sig ) /2;

	int xs,my_rank=0, np=1;
#ifdef MPI
        int root=0;
#endif

#ifdef MPI
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		double *ur = (double *) malloc( nx_red*NY*NZ *sizeof(double));
			for( i=0; i<= (nx_red*NY*NZ)-1; i++)	ur[i] = 0.00;

	for( x=my_rank; x<=nx-1; x += np )
	{
	    xs = (int) (x/np);
	    for( y=0; y<=ny-1; y++)
		for( z=0; z<=nz-1; z++)
		{
			indx = ii( xs, y, z );
			for( i=1; i<=nu_sites; i++){
				if( strncmp( "wall", u[i].element , 4 ) == 0 )
				 	r = rw(z, u[i].z);
				else
				 	r = rx(x,y,z, u[i].x, u[i].y, u[i].z);

				ur[ indx ] += 4.0 /t *( (ep12_ux[i] *pow( (sig_ux[i]/r) ,12) )- (ep6_ux[i] *pow( (sig_ux[i]/r) ,6)));
						   }
		}
	}

#ifdef MPI
	double *ur_ux = NULL;
	int nx_ext = nx_red * np;
	if( my_rank == 0 )	ur_ux = (double *) malloc( nx_ext *NYZ *sizeof(double));

	MPI_Barrier( MPI_COMM_WORLD );

	for( i=0; i<=nx_red-1; i++)
	    MPI_Gather( (ur + ix(i)), NYZ, MPI_DOUBLE, (ur_ux +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );

	if( my_rank == 0) print_3d( fname, ur_ux, nx, ny, nz, PND[ns]);

	MPI_Barrier( MPI_COMM_WORLD );
 #else
	print_3d(fname, ur, nx, ny, nz, PND[ns]);
  #endif

}

/****************************************************************************************/
/*					U(r)_12 - RBC					*/
/*			calc u(r)_term for repulsive bridge function			*/
/*											*/
/****************************************************************************************/

void  calc_and_print_ur_lj12( U_PAR2 * u, int nu_sites, double ep12, double ep6, double sig, char fname[], int ns )
{
	int i, x, y, z, indx;
	int nx=NX, ny=NY, nz=NZ;
	double r, t = TEMP;

	double ep12_ux[nu_sites +1], sig_ux[nu_sites +1];
	for( i=1; i<=nu_sites; i++)	ep12_ux[i] = sqrt( u[i].ep12 * ep12 );
	for( i=1; i<=nu_sites; i++)	sig_ux[i]  = (u[i].sig + sig ) /2;

	int xs,my_rank=0, np=1;

#ifdef MPI
        int root = 0;
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		double *ur = (double *) malloc( nx_red*NY*NZ *sizeof(double));
			for( i=0; i<= (nx_red*NY*NZ)-1; i++)	ur[i] = 0.00;

	for( x=my_rank; x<=nx-1; x += np )
	{
	    xs = (int) (x/np);
	    for( y=0; y<=ny-1; y++)
		for( z=0; z<=nz-1; z++)
		{
			indx = ii( xs, y, z );
			for( i=1; i<=nu_sites; i++){
				if( strncmp( "wall", u[i].element , 4 ) == 0 )
				 	r = rw(z, u[i].z);
				else
				 	r = rx(x,y,z, u[i].x, u[i].y, u[i].z);
				ur[ indx ] += 4 /t *( (ep12_ux[i] *pow( (sig_ux[i]/r) ,12) ) );
						   }
		}
	}

#ifdef MPI
	double *ur_ux = NULL;
	int nx_ext = nx_red * np;
	if( my_rank == 0 )	ur_ux = (double *) malloc( nx_ext *NYZ *sizeof(double));

	MPI_Barrier( MPI_COMM_WORLD );

	for( i=0; i<=nx_red-1; i++)
	    MPI_Gather( (ur + ix(i)), NYZ, MPI_DOUBLE, (ur_ux +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );

	if( my_rank == 0)	print_3d( fname, ur_ux, nx, ny, nz, PND[ns]);

	MPI_Barrier( MPI_COMM_WORLD );
 #else
	print_3d(fname, ur, nx, ny, nz, PND[ns]);
  #endif

}


/****************************************************************************************/
/*					U(r) - WCA					*/
/*				calc wca potential routines				*/
/*											*/
/****************************************************************************************/
#ifdef WCA

void  calc_and_print_ur_wca( U_PAR2 * u, int nu_sites, double ep12, double ep6, double sig, char fname[], int ns )
{
	int i, x, y, z, indx;
	int nx=NX, ny=NY, nz=NZ;
	double r, t = TEMP, u_lj, u_wca;

	double ep12_ux[nu_sites +1], ep6_ux[nu_sites +1], sig_ux[nu_sites +1], r_wca[ nu_sites +1];
	for( i=1; i<=nu_sites; i++)	ep12_ux[i] = sqrt( u[i].ep12 * ep12 );
	for( i=1; i<=nu_sites; i++)	ep6_ux[i]  = sqrt( u[i].ep6 * ep6 );
	for( i=1; i<=nu_sites; i++)	sig_ux[i]  = (u[i].sig + sig ) /2;
	for( i=1; i<=nu_sites; i++)	r_wca[i] = sig_ux[i]  * (1.122462048309372981);

	int xs,my_rank=0, np=1, root=0;

#ifdef MPI
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		int nx_ext = nx_red * np;
		double *ur = (double *) malloc( nx_red*NY*NZ *sizeof(double));
			for( i=0; i<= (nx_red*NY*NZ)-1; i++)	ur[i] = 0.00;

	for( x=my_rank; x<=nx-1; x += np )
	{
	    xs = (int) (x/np);
	    for( y=0; y<=ny-1; y++)
		for( z=0; z<=nz-1; z++)
		{
			indx = ii( xs, y, z );
			for( i=1; i<=nu_sites; i++){
				if( strncmp( "wall", u[i].element , 4 ) == 0 )
				 	r = rw(z, u[i].z);
				else
				 	r = rx(x,y,z, u[i].x, u[i].y, u[i].z);
				u_lj = 4 /t *( (ep12_ux[i] *pow( (sig_ux[i]/r) ,12) )- (ep6_ux[i] *pow( (sig_ux[i]/r) ,6)));

				if( r <= r_wca[i] )
					u_wca = u_lj + (ep6_ux[i]/t) ;
				else
					u_wca = 0.00;

				ur[ indx ] += ( ((1.0-WCA)*u_lj) + (WCA*u_wca) );

						   }
		}
	}

#ifdef MPI
	double *ur_ux;
	if( my_rank == 0 )	ur_ux = (double *) malloc( nx_ext *NYZ *sizeof(double));

	MPI_Barrier( MPI_COMM_WORLD );

	for( i=0; i<=nx_red-1; i++)
	    MPI_Gather( (ur + ix(i)), NYZ, MPI_DOUBLE, (ur_ux +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );

	if( my_rank == 0)	print_3d( fname, ur_ux, nx, ny, nz, PND[ns]);

	MPI_Barrier( MPI_COMM_WORLD );
 #else
	print_3d(fname, ur, nx, ny, nz, PND[ns]);
  #endif

}

#endif


/****************************************************************************************/
/*					U(r) - Coulomb					*/
/*				calc Coulomb potential routines				*/
/*											*/
/****************************************************************************************/
void calc_and_print_ur_clmb( U_PAR2 *u, int nu_sites, double z_x, char fname[], int ns )
/*	ur_clmb,		*/
{
	int i, x, y, z, nx = NX, ny = NY, nz = NZ;
	int indx;
	double r, t = TEMP, K = CONST_CLMB;

	int xs;
	int my_rank = 0, np = 1;

#ifdef MPI
        int root = 0;
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np ;
		if( nx_rmd != 0 )   nx_red++;
		double *ur = (double *) malloc( nx_red*NY*NZ *sizeof(double));
			for( i=0; i<= (nx_red*NY*NZ)-1; i++) 	ur[i] = 0.00;

	for( x=my_rank; x<=nx-1; x+=np )
	{
	    xs = (int) (x/np);
	    for( y=0; y<=ny-1; y++)
		for( z=0; z<=nz-1; z++)
		{
			indx = ii(xs,y,z);
			for( i=1; i<=nu_sites; i++)
			{
				if( u[i].charge == 0.0 ) continue;
				r = rx( x, y, z, u[i].x, u[i].y, u[i].z );
				ur[ indx ] += K /t *( u[i].charge * z_x /r);
			}
		}
	}

#ifdef MPI
	double *urx = NULL;
	int nx_ext = nx_red * np;
	if( my_rank == 0 )
		urx = (double *) malloc( nx_ext*NYZ *sizeof(double));

	MPI_Barrier( MPI_COMM_WORLD );

	for( i=0; i<=nx_red-1; i++)
	    MPI_Gather( (ur + ix(i)), NYZ, MPI_DOUBLE, (urx + ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );

	if( my_rank == root)	print_3d( fname, urx, nx, ny, nz, PND[ns]);

	MPI_Barrier( MPI_COMM_WORLD );
 #else
	print_3d(fname, ur, nx, ny, nz, PND[ns]);
  #endif


}

/****************************************************************************************/
/*											*/
/*				U(r)*erf(a*r) - Coulomb					*/
/*											*/
/****************************************************************************************/

void calc_and_print_ur_erf( U_PAR2 *u, int nu_sites, double z_x, char fname[], int ns )
{
	int i, x, y, z, indx;
	int nx = NX, ny = NY, nz = NZ;
	double t = TEMP, Kc = CONST_CLMB;
	double a = A_ERF, r;

	int xs;
	int my_rank = 0, np = 1;

#ifdef MPI
        int root = 0;
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		double *ur = (double *) malloc( nx_red*NY*NZ *sizeof(double));
			for( i=0; i<= (nx_red*NY*NZ)-1; i++)
		    		ur[i] = 0.00;

	for( x=my_rank; x<=nx-1; x+=np )
	{
	    xs = (int) (x/np);
	    for( y=0; y<=ny-1; y++)
		for( z=0; z<=nz-1; z++)
		{
			indx = ii(xs,y,z);
			for( i=1; i<=nu_sites; i++){
				if( u[i].charge == 0.0 ) continue;
				r = rx( x,y,z, u[i].x, u[i].y, u[i].z);
				ur[ indx ] += erf( a*r ) * Kc/t *(u[i].charge * z_x / r );
						   }
		}
	}

#ifdef MPI
	double *ur_erf = NULL;
	int nx_ext = nx_red * np;
	if( my_rank == root )
		ur_erf = (double *) malloc( nx_ext*NYZ *sizeof(double));
        MPI_Barrier( MPI_COMM_WORLD);

	for( i=0; i<=nx_red-1; i++)
	    MPI_Gather( (ur +ix(i)), NYZ, MPI_DOUBLE, (ur_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
			/*	ur_erf[ii(u_x,u_y,u_z)] = 0.00;*/

	if( my_rank == root)	print_3d( fname, ur_erf, nx, ny, nz, PND[ns] );
        MPI_Barrier( MPI_COMM_WORLD );
#else
	print_3d( fname, ur, nx, ny, nz, PND[ns] );
#endif



}

/****************************************************************************************/
/*											*/
/*			       U(k) -> U(r)*erf(a*r) in fourier space			*/
/*											*/
/****************************************************************************************/

void calc_and_print_uk_erf( U_PAR2 *u, int nu_sites, double z_x, char fname[], int ns )
{
	int i, x, y, z, indx, nx = NX, ny = NY, nz = NZ;
	int u_x=CX, u_y=CY, u_z=CZ;
	double al, Kc, t, k, kr, xf=1.0, yf=1.0, zf=1.0;
		al = A_ERF;
		Kc = CONST_CLMB;
		t = TEMP;

	int xs;
	int my_rank=0, np=1, root=0;

#ifdef MPI
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		int nx_ext = nx_red * np;
		double *ukr = (double *) malloc( nx_red*NY*NZ *sizeof(double));
		double *uki = (double *) malloc( nx_red*NY*NZ *sizeof(double));

		for( i=0; i<= (nx_red*NY*NZ)-1; i++){
		    ukr[i] = 0.00;
		    uki[i] = 0.00;
		}

	for( x=my_rank; x<=nx-1; x+=np )
	{
	     	xs = (int) (x/np);
		if( x < u_x ) xf=(-1.0); else xf=1.0;
	   for( y=0; y<=ny-1; y++)
	   {	 if( y < u_y ) yf=(-1.0); else yf=1.0;
	     for( z=0; z<=nz-1; z++)
	     {	  if( z < u_z ) zf=(-1.0); else zf=1.0;
			indx = ii(xs,y,z);
			k = k0( x,y,z);
			for( i=1; i<=nu_sites; i++){
				if( u[i].charge == 0.0 ) continue;
				kr = ( xf*k0(x,u_y,u_z)*u[i].x  + yf*k0(u_x, y, u_z)*u[i].y  + zf*k0(u_x,u_y,z)*u[i].z );
				ukr[indx] += (4*Pi/(k*k)) *exp( (-1)*k*k/ (4*al*al)) *(Kc/t)*(u[i].charge * z_x) *cos( kr );;
				uki[indx] += (4*Pi/(k*k)) *exp( (-1)*k*k/ (4*al*al)) *(Kc/t)*(u[i].charge * z_x) *(-1)*sin( kr );;
						   }
	  }
	 }
	}

	fftw_complex *uk_erf;

#ifdef MPI
        double *ukr_erf = NULL, *uki_erf = NULL;

	if( my_rank == root ){
		ukr_erf = (double *) malloc( nx_ext *NYZ *sizeof(double));
		uki_erf = (double *) malloc( nx_ext *NYZ *sizeof(double));
	}
	for( i=0; i<=nx_red-1; i++){
	    MPI_Gather( (ukr +ix(i)), NYZ, MPI_DOUBLE, (ukr_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
	    MPI_Gather( (uki +ix(i)), NYZ, MPI_DOUBLE, (uki_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
        }
	if( my_rank == 0 ){
		uk_erf = (fftw_complex *) malloc( nx_ext *NYZ *sizeof(fftw_complex));
		for( i=0; i<= (nx_ext*NYZ)-1; i++){
			uk_erf[i][0] = ukr_erf[i];
			uk_erf[i][1] = uki_erf[i];
		}
	}
 #else
	uk_erf = (fftw_complex *) malloc( nx_ext *NYZ *sizeof( fftw_complex ));
		for( i=0; i<= (nx_ext*NYZ) -1; i++){
				uk_erf[i][0] = ukr[i];
				uk_erf[i][1] = uki[i];
						}
  #endif

	if( my_rank == root )
	{
				indx = ii(u_x,u_y,u_z);
				uk_erf[indx][0]  = ( uk_erf[ ii(u_x+1, u_y, u_z) ][0] + uk_erf[ ii(u_x-1, u_y, u_z) ][0] )/6;
				uk_erf[indx][1]  = ( uk_erf[ ii(u_x+1, u_y, u_z) ][1] + uk_erf[ ii(u_x-1, u_y, u_z) ][1] )/6;
				uk_erf[indx][0] += ( uk_erf[ ii(u_x, u_y+1, u_z) ][0] + uk_erf[ ii(u_x, u_y-1, u_z) ][0] )/6;
				uk_erf[indx][1] += ( uk_erf[ ii(u_x, u_y+1, u_z) ][1] + uk_erf[ ii(u_x, u_y-1, u_z) ][1] )/6;
				uk_erf[indx][0] += ( uk_erf[ ii(u_x, u_y, u_z+1) ][0] + uk_erf[ ii(u_x, u_y, u_z-1) ][0] )/6;
				uk_erf[indx][1] += ( uk_erf[ ii(u_x, u_y, u_z+1) ][1] + uk_erf[ ii(u_x, u_y, u_z-1) ][1] )/6;
				uk_erf[indx][0] = 0.00;
				uk_erf[indx][1] = 0.00;
	}

	if( my_rank == root)
		print_cmplx_3d( fname, uk_erf, nx, ny, nz, PND[ns]);
}


/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*					U(r) - Coulomb-Ewald				*/
/*				Ewald Coulomb potential routines			*/
/*											*/
/****************************************************************************************/


void calc_and_print_ur_clmb_ewald( U_PAR2 *u, int nu_sites, double z_x, char fname[], int ns )
{
	int i, x, y, z, indx, nx = NX, ny = NY, nz = NZ;
	int cx=CX, cy=CY, cz=CZ;
	double  r, k, kr, kx, ky, kz;
	double	Kc = CONST_CLMB;
	double	t = TEMP;
	double gs = GAUSS_SMEAR;

	int xs;
	int my_rank = 0, np = 1, root = 0;

#ifdef MPI
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		int nx_ext = nx_red * np;

		double *ukr   = (double *) malloc( nx_red*NY*NZ *sizeof(double));
		double *uki   = (double *) malloc( nx_red*NY*NZ *sizeof(double));
		double *uerfc = (double *) malloc( nx_red*NY*NZ *sizeof(double));

		for( i=0; i<= (nx_red*NY*NZ)-1; i++){
		    ukr[i] = 0.00;
		    uki[i] = 0.00;
		    uerfc[i] = 0.00;
		}


		for( x=my_rank; x<=nx-1; x+=np )
		{
		     	xs = (int) (x/np);
			kx = (double) k0(x,cy,cz);
			if( x < cx ){ kx*=(-1.0); }

		   for( y=0; y<=ny-1; y++)
		   {
			ky = (double) k0(cx,y,cz);
			if( y < cy ) {ky*=(-1.0); }

		     for( z=0; z<=nz-1; z++)
		     {
			  	kz = (double) k0(cx,cy,z);
			  	if( z < cz ){ kz*=(-1.0);}

				indx = ii(xs,y,z);
				k = k2( x,y,z);
				for( i=1; i<=nu_sites; i++){
					if( u[i].charge == 0.0 ) continue;
					kr = ( kx * u[i].x  + ky * u[i].y  + kz * u[i].z );
					ukr[indx] += (u[i].charge * z_x) *( 1.0) *cos( kr );
					uki[indx] += (u[i].charge * z_x) *(-1.0) *sin( kr );

				 	r = rx(x,y,z, u[i].x, u[i].y, u[i].z);
					uerfc[indx] += (u[i].charge * z_x) *erfc( r/gs )/r;

				}

				ukr[indx] *= (4*Pi/(k*k)) *(Kc/t) *exp( (-1.0)*(k*k)*(gs*gs)/4.0 );
				uki[indx] *= (4*Pi/(k*k)) *(Kc/t) *exp( (-1.0)*(k*k)*(gs*gs)/4.0 );

		  }
		 }
		}

	/*collect from all processes*/
	fftw_complex *uk_erf = NULL;
        double *uerfc2 = NULL;

#ifdef MPI
	double *ukr_erf = NULL, *uki_erf = NULL;
	if( my_rank == root ){
		ukr_erf = (double *) malloc( nx_ext *NYZ *sizeof(double));
		uki_erf = (double *) malloc( nx_ext *NYZ *sizeof(double));
		uerfc2 = (double *) malloc( nx_ext *NYZ *sizeof(double));
			     }
	for( i=0; i<=nx_red-1; i++){
	    MPI_Gather( (ukr +ix(i)), NYZ, MPI_DOUBLE, (ukr_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
	    MPI_Gather( (uki +ix(i)), NYZ, MPI_DOUBLE, (uki_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
	    MPI_Gather( (uerfc +ix(i)), NYZ, MPI_DOUBLE, (uerfc2 +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
				   }
	if( my_rank == 0 ){
		uk_erf = (fftw_complex *) malloc( nx_ext *NYZ *sizeof(fftw_complex));
		for( i=0; i<= (nx_ext*NYZ)-1; i++){
				uk_erf[i][0] = ukr_erf[i];
				uk_erf[i][1] = uki_erf[i];
						  }
			  }
 #else
	uk_erf = (fftw_complex *) malloc( nx_ext *NYZ *sizeof( fftw_complex ));
		for( i=0; i<= (nx_ext*NYZ) -1; i++){
				uk_erf[i][0] = ukr[i];
				uk_erf[i][1] = uki[i];
						}
	uerfc2 = (double *) malloc( nx_ext *NYZ *sizeof(double));
		for( i=0; i<= NNN-1; i++)
			uerfc2[i] = uerfc[i];
  #endif

	fftw_complex *uk, *ur;
	double      *ur1 ;

	/*k=0 point*/
	if( my_rank == root )
	{
			indx = ii(cx,cy,cz);
			uk_erf[indx][0] = 0.00;
			uk_erf[indx][1] = 0.00;


		uk = (fftw_complex *) malloc( NNN *sizeof(fftw_complex) );
		ur = (fftw_complex *) malloc( NNN *sizeof(fftw_complex) );
		ur1 = (double *) malloc( NNN *sizeof(double) );

		for(i=0; i<=NNN-1; i++){
			uk[i][0] = uk_erf[i][0];
			uk[i][1] = uk_erf[i][1];
		}

		shift_origin_complex_inplace( uk, SYS );
		invfftw_3d( uk , ur );	/*in, out*/
		unshift_origin_complex_inplace( ur, SYS);

		for( i=0; i<=NNN-1; i++)
			ur1[i] = ur[i][0] + uerfc2[i];

	}



	/* Print to disk file*/
	if( my_rank == root)
		print_3d( fname, ur1, nx, ny, nz, PND[ns]);

	//fftw routines, maybe put these into a subroutine file

	if( my_rank == root){
		free(ukr);	free(uki);
		free(uk);	free(ur);
		free(ur1);
		free(uerfc);	free(uerfc2);
	}

}

/****************************************************************************************/
/*											*/
/*			Calculate EWALD 1D rad						*/
/*											*/
/****************************************************************************************/


void calc_and_print_ur_clmb_ewald_rad( U_PAR2 *u, int nu_sites, double z_x, char fname[], int ns )
{
	int i, x, y, z, indx, nx = NX, ny = NY, nz = NZ;
	int cx=CX, cy=CY, cz=CZ;
	double  r, k, kr, kx, ky, kz;
	double	Kc = CONST_CLMB;
	double	t = TEMP;
	double gs = GAUSS_SMEAR;

	int xs;
	int my_rank = 0, np = 1, root = 0;

#ifdef MPI
		MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
		MPI_Comm_size( MPI_COMM_WORLD, &np );
#endif
		int nx_red = (int) (nx/np) ;
		int nx_rmd = nx % np;
		if( nx_rmd != 0 )   nx_red++;
		int nx_ext = nx_red * np;

		double *ukr   = (double *) malloc( nx_red*NY*NZ *sizeof(double));
		double *uki   = (double *) malloc( nx_red*NY*NZ *sizeof(double));
		double *uerfc = (double *) malloc( nx_red*NY*NZ *sizeof(double));

		for( i=0; i<= (nx_red*NY*NZ)-1; i++){
		    ukr[i] = 0.00;
		    uki[i] = 0.00;
		    uerfc[i] = 0.00;
		}


		for( x=my_rank; x<=nx-1; x+=np )
		{
		     	xs = (int) (x/np);
			kx = (double) k0(x,cy,cz);
			if( x < cx ){ kx*=(-1.0); }

		   for( y=0; y<=ny-1; y++)
		   {
			ky = (double) k0(cx,y,cz);
			if( y < cy ) {ky*=(-1.0); }

		     for( z=0; z<=nz-1; z++)
		     {
			  	kz = (double) k0(cx,cy,z);
			  	if( z < cz ){ kz*=(-1.0);}

				indx = ii(xs,y,z);
				k = k2( x,y,z);
				for( i=1; i<=nu_sites; i++){
					if( u[i].charge == 0.0 ) continue;
					kr = ( kx * u[i].x  + ky * u[i].y  + kz * u[i].z );
					ukr[indx] += (u[i].charge * z_x) *( 1.0) *cos( kr );
					uki[indx] += (u[i].charge * z_x) *(-1.0) *sin( kr );

				 	r = rx(x,y,z, u[i].x, u[i].y, u[i].z);
					uerfc[indx] += (u[i].charge * z_x) *erfc( r/gs )/r;

				}

				ukr[indx] *= (4*Pi/(k*k)) *(Kc/t) *exp( (-1.0)*(k*k)*(gs*gs)/4.0 );
				uki[indx] *= (4*Pi/(k*k)) *(Kc/t) *exp( (-1.0)*(k*k)*(gs*gs)/4.0 );

		  }
		 }
		}

	/*collect from all processes*/
	fftw_complex *uk_erf = NULL;
	double *uerfc2 = NULL;

#ifdef MPI
        double *ukr_erf = NULL, *uki_erf = NULL;

	if( my_rank == root ){
		ukr_erf = (double *) malloc( nx_ext *NYZ *sizeof(double));
		uki_erf = (double *) malloc( nx_ext *NYZ *sizeof(double));
		uerfc2 = (double *) malloc( nx_ext *NYZ *sizeof(double));
                             }
	for( i=0; i<=nx_red-1; i++){
	    MPI_Gather( (ukr +ix(i)), NYZ, MPI_DOUBLE, (ukr_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
	    MPI_Gather( (uki +ix(i)), NYZ, MPI_DOUBLE, (uki_erf +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
	    MPI_Gather( (uerfc +ix(i)), NYZ, MPI_DOUBLE, (uerfc2 +ix(i*np)), NYZ, MPI_DOUBLE, root, MPI_COMM_WORLD );
				   }
	if( my_rank == 0 ){
		uk_erf = (fftw_complex *) malloc( nx_ext *NYZ *sizeof(fftw_complex));
		for( i=0; i<= (nx_ext*NYZ)-1; i++){
				uk_erf[i][0] = ukr_erf[i];
				uk_erf[i][1] = uki_erf[i];
						  }
			  }
 #else
	uk_erf = (fftw_complex *) malloc( nx_ext *NYZ *sizeof( fftw_complex ));
		for( i=0; i<= (nx_ext*NYZ) -1; i++){
				uk_erf[i][0] = ukr[i];
				uk_erf[i][1] = uki[i];
						}
	uerfc2 = (double *) malloc( nx_ext *NYZ *sizeof(double));
		for( i=0; i<= NNN-1; i++)
			uerfc2[i] = uerfc[i];
  #endif

	fftw_complex *uk, *ur;
	double      *ur1 ;

	/*k=0 point*/
	if( my_rank == root )
	{
			indx = ii(cx,cy,cz);
			uk_erf[indx][0] = 0.00;
			uk_erf[indx][1] = 0.00;


		uk = (fftw_complex *) malloc( NNN *sizeof(fftw_complex) );
		ur = (fftw_complex *) malloc( NNN *sizeof(fftw_complex) );
		ur1 = (double *) malloc( NNN *sizeof(double) );

		for(i=0; i<=NNN-1; i++){
			uk[i][0] = uk_erf[i][0];
			uk[i][1] = uk_erf[i][1];
		}

		shift_origin_complex_inplace( uk, SYS );
		invfftw_3d( uk , ur );	/*in, out*/
		unshift_origin_complex_inplace( ur, SYS);

		for( i=0; i<=NNN-1; i++)
			ur1[i] = ur[i][0] + uerfc2[i];

	}



	/* Print to disk file*/
	if( my_rank == root)
		print_3d( fname, ur1, nx, ny, nz, PND[ns]);

	//fftw routines, maybe put these into a subroutine file

	if( my_rank == root){
		free(ukr);	free(uki);
		free(uk);	free(ur);
		free(ur1);
		free(uerfc);	free(uerfc2);
	}

}

/****************************************************************************************/
/*											*/
/*			Calculate EWALD Correction Constant				*/
/*											*/
/****************************************************************************************/


double ** get_hr1vv_dis( void );
double ** calc_hk1vv_dis( double ** );
double ** calc_dkn_mat( int, double );
void print_ewald( double *, double , int );
void calc_hk_terms( double *, double *, double , int , int );  /* out, in, rad, N, #n of k */

void calc_and_print_ewald_corrections( void )
{

	int i, j, k, idx1=0;

	double Kc = CONST_CLMB;
	double t = TEMP;
	int n_pts  = (int) get_dval( ARGV1, "N_PTS" );
	double rad = (double) get_dval( ARGV1, "RADIUS");
	double dk = Pi/rad;
	double tc=0.0;
	int num_k = 10;

	double    *h0 = (double *) malloc( NRSITES *sizeof(double) );

	double **e_sum = (double **) matrix_malloc0( NRSITES, num_k );
	double **hr1vv = (double **) get_hr1vv_dis( );
#ifndef MPI
	double  *coef = (double *) malloc( num_k *sizeof(double) );
	double **a = (double **) calc_dkn_mat( num_k, dk );
#endif

	double **hk    = (double **) malloc( DIS_NUM *sizeof(double ));
			for( i=0; i<=DIS_NUM-1; i++){
				*(hk+i) = (double *) malloc( num_k *sizeof(double) );
				calc_hk_terms( *(hk+i), *(hr1vv +i) , rad, n_pts, num_k  );
			}

	double ***hk_mat = (double ***) malloc( NRSITES *sizeof(double) );
			for( i=0; i<=NRSITES-1; i++)
			    *(hk_mat +i) = (double **) malloc( NRSITES *sizeof(double));

			for( i=0; i<=NRSITES-1; i++)
			    for( j=i; j<=NRSITES-1; j++){
				*(*(hk_mat +i)+j) = *(hk + idx1);
				idx1++;
			    }

			for( i=1; i<=NRSITES-1; i++)
			    for( j=0; j<= i-1; j++)
				*(*(hk_mat +i) +j) = *(*(hk_mat +j) +i);


	double *wk_oh= (double *) malloc( num_k *sizeof(double) );
	double *wk_hh= (double *) malloc( num_k *sizeof(double) );
			if( TYPE == 1 || TYPE == 2 )
			   for( k=1; k<=num_k; k++){
				wk_oh[k-1] = sin( k*dk * BOND_MAT[0]) / (k*dk *BOND_MAT[0]);
				wk_hh[k-1] = sin( k*dk * BOND_MAT[1]) / (k*dk *BOND_MAT[1]);
			   }


		if( TYPE == 0 ){
				for( i=0; i<=NRSITES-1; i++)
					for( k=1; k<=num_k; k++){
					    e_sum[i][k-1] = CHARGE[i] / (k*k*dk*dk);
					    for( j=0; j<=NRSITES-1; j++)
						e_sum[i][k-1] += CHARGE[i]*PND[j]*hk_mat[j][i][k-1]*REDUN[j] /(k*k*dk*dk);
					}
		} else
		 if( TYPE == 1){
				for( k=1; k<=num_k; k++){
					e_sum[0][k-1] = CHARGE[0]/ (k*k*dk*dk);
					e_sum[1][k-1] = CHARGE[1]/ (k*k*dk*dk);
				}
				for( k=1; k<=num_k; k++){
					e_sum[0][k-1] += CHARGE[1] * wk_oh[k-1] * 2.00 /(k*k*dk*dk);
					e_sum[1][k-1] += (CHARGE[0] * wk_oh[k-1] + CHARGE[1] * wk_hh[k-1]) /(k*k*dk*dk);
				}

				for( i=0; i<=NRSITES-1; i++)
					for( k=1; k<=num_k; k++){
					    for( j=0; j<=NRSITES-1; j++)
						e_sum[i][k-1] += CHARGE[j]*PND[j]*hk_mat[j][i][k-1]*REDUN[j] /(k*k*dk*dk);
					}
		} else
		 if( TYPE == 2){
				for( k=1; k<=num_k; k++){
					e_sum[0][k-1] = CHARGE[1] * wk_oh[k-1] * 2.00 /(k*k*dk*dk);
					e_sum[1][k-1] = (CHARGE[0] * wk_oh[k-1] + CHARGE[1] * wk_hh[k-1]) /(k*k*dk*dk);
				}

				for( i=0; i<=NRSITES-1; i++)
					for( k=1; k<=num_k; k++){
					    e_sum[i][k-1] = CHARGE[i]/ (k*k*dk*dk);
					    for( j=0; j<=NRSITES-1; j++)
						e_sum[i][k-1] += CHARGE[j]*PND[j]*hk_mat[j][i][k-1]*REDUN[j] /(k*k*dk*dk);
					}


		} else
		 if( TYPE == 3){
					printf("\nHasn't been done yet.\n"); 	fflush(stdout);
		} else{
				printf("\nNO type specified - bye!\n"); 	fflush(stdout);
				exit(1);
		}
#ifndef MPI
		for( i=0; i<=NRSITES-1; i++){
			Newton_Solver( coef, a, *(e_sum +i), num_k );
			h0[i] = coef[0];
		}
#else
		printf("\newald correction can't be calculated with the mpi version\n");
#endif
			for( i=1; i<=NU_SITES; i++)
				tc += U[i].charge;

			for( i=0; i<=NRSITES-1; i++)
				h0[i] *= 4.0*Pi*tc*Kc/t/(LX*LY*LZ);

		print_ewald( h0, tc, NRSITES);
}

void print_ewald( double *h0, double tc, int nr)
{
	int i;
	FILE *out;
		if( (out = fopen("ewald.dat", "w")) == NULL)
			printf("\nCouldn't open ewald.dat\n");

	fprintf(out, "#TOTAL_SOLUTE_CHARGE	%f\n", tc);
	fprintf(out, "#EWALD_BGC\t");
	for( i=0; i<=nr-1; i++)
	    fprintf(out, "%f\t", h0[i]);
	fprintf(out, "\n");

}
double ** get_hr1vv_dis( void )
{
	int j;
	double **hr1vv = (double **) malloc( DIS_NUM *sizeof(double) );
			for( j=0; j<=DIS_NUM-1; j++)
				*(hr1vv +j) = (double *) get_dis( ARGV1, j+1);
	return hr1vv;
}

double ** calc_dkn_mat( int nk, double dk)
{
	int i, j;

	double **a = (double **) malloc( nk * sizeof(double) );

		for( i=0; i<=nk-1; i++)
			*(a+i) = (double *) malloc( nk *sizeof(double) );

		for( i=0; i<=nk-1; i++){
		    a[i][0] = 1.000;
		    for( j=1; j<=nk-1; j++)
			a[i][j] = pow( ((i+1)*dk), j);
		}

	return a;
}


void calc_hk_terms(  double * y, double * x, double R, int Nn, int Nk )   /* out, in, rad, N, #n of k */
{

	int j, l;
	double C2, dr;
	double Tmp=0;

	C2 = R*R/((Nn-1)*Pi);
        dr = (R/(Nn-1));

        /* End of initial point mess*/
	for( l=1; l<=Nk; l++){
		Tmp=0;
                for( j=1; j<=Nn-2; j++)
			Tmp += sin(Pi/(Nn-1)*j*l)*x[j]*C2*j/l;
             	y[l-1] = Tmp*4*Pi*dr;
        }
}


/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				FFTW_3D routines					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/



void fftw_3d( fftw_complex *in_r , fftw_complex *out_k )
{
	int    nx=NX, ny=NY, nz=NZ, nnn = NNN;		/*EXTERN*/
	double lx=LX, ly=LY, lz=LZ;			/*EXTERN*/

	int i;
	double dx = lx/(nx-1);
	double dy = ly/(ny-1);
	double dz = lz/(nz-1);
	double con = dx*dy*dz;

#ifdef FFTW_THREADS
  	fftw_init_threads();
	fftw_plan_with_nthreads(NUM_THREADS);
#endif

	fftw_plan rtok;
		rtok = fftw_plan_dft_3d( nx, ny, nz, in_r, out_k, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute( rtok );
	for( i=0; i<=nnn-1; i++){  out_k[i][0] *= con;  out_k[i][1] *= con;	}/*Normalizing*/
}



void invfftw_3d( fftw_complex *in_k, fftw_complex *out_r )
{
	int    nx=NX, ny=NY, nz=NZ, nnn = NNN;		/*EXTERN*/
	double lx=LX, ly=LY, lz=LZ;			/*EXTERN*/

	int i;
	double dx = lx/(nx-1),  dy = ly/(ny-1),  dz = lz/(nz-1);
	double con = dx*dy*dz*nnn;

#ifdef FFTW_THREADS
  	fftw_init_threads();
	fftw_plan_with_nthreads(NUM_THREADS);
#endif

	fftw_plan ktor;
		ktor = fftw_plan_dft_3d( nx, ny, nz, in_k, out_r, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute( ktor );
	for ( i=0; i<=nnn-1; i++){ out_r[i][0] *= (1/con);  out_r[i][1] *= (1/con); }/*Normalizing*/
}








/****************************************************************************************/
/*											*/
/*											*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

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
	tmp += pow(fabs( z*DZ - (CZ*DZ + uz) ), 2) ;
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

	if ( tmp == 0.0 )
		tmp = 1e-10 ;
	return tmp;
}


double rw( int z, double uz)
{
	double tmp;
	tmp = fabs( z*DZ - (CZ*DZ + uz) ) ;

	if ( tmp == 0.0 )
		tmp = 1e-10 ;
	return tmp;
}

double kk( int x, int y, int z, double ux, double uy, double uz)
{
	double tmp;

	tmp  = pow(fabs( x - (CX + (ux/DX)) ),2) *pow(DKX,2);
	tmp += pow(fabs( y - (CY + (uy/DY)) ),2) *pow(DKY,2);
	tmp += pow(fabs( z - (CZ + (uz/DZ)) ),2) *pow(DKZ,2);
	tmp = sqrt( tmp );
	if ( tmp == 0.0 )
		tmp = 1e-10 ;
	return tmp;
}

double k0( int x, int y, int z )
{
	double tmp;

	tmp  = pow(fabs(x - CX),2) * pow(DKX,2);
	tmp += pow(fabs(y - CY),2) * pow(DKY,2);
	tmp += pow(fabs(z - CZ),2) * pow(DKZ,2);
	if ( tmp == 0.0 )
		tmp = 1e-10 ;
	tmp = sqrt( tmp );

	return tmp;
}


double k2( int x, int y, int z )
{
	double tmp;

	tmp  = pow(fabs(x - CX),2) * pow(DKX2,2);
	tmp += pow(fabs(y - CY),2) * pow(DKY2,2);
	tmp += pow(fabs(z - CZ),2) * pow(DKZ2,2);
	tmp = sqrt( tmp );

	if ( tmp == 0.0 )
		tmp = 1e-10 ;
	return tmp;
}


void print_3d( char name[], double *v, int nx, int ny, int nz, double pnd )
{

	sprintf( name, "%s.%s", name, FILE_TYPE );

	if( strncmp( "sit", FILE_TYPE, 3 ) == 0 )
		print_sit( name, v, SYS );
	 else
	   if( strncmp( "jh3d", FILE_TYPE, 8 ) == 0 )
		print_jh3d( name, v, SYS, TEMP, pnd );
	   else
		printf("\nFile name not specified correctly\n" ); fflush(stdout);

}

void print_cmplx_3d( char name[], fftw_complex *v, int nx, int ny, int nz, double pnd )
{

	sprintf( name, "%s.%s", name, FILE_TYPE );

	if( strncmp( "sit", FILE_TYPE, 3 ) == 0 )
		print_cmplx_sit( name, v, SYS );
	 else
	   if( strncmp( "jh3d", FILE_TYPE, 8 ) == 0 )
		print_cmplx_jh3d( name, v, SYS, TEMP, pnd );
	   else
		printf("\nFile name not specified correctly\n" ); fflush(stdout);
}


