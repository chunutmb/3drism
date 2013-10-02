graveyard

void newton_gmres_ml_iter(  double *, double *, double *, double *, double *, double * );
void picard_iter( double **, double **, double ** );
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
/*				Picard Iterations					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void picard_iter( double **cr2_s,  double **cr_s, double **tr_s )
{
	int nnn=NNN;				/*EXTERN*/
	double *hk3_oo = HK3[0][0];
	double *hk3_oh = HK3[0][1];
	double *hk3_hh = HK3[1][1];
	double *wk_oh = WK_OH;
	double *wk_hh = WK_HH;
	double **ur_s = UR_S;
	fftw_complex **uk_l = UK_L ;

	int i, m, rc;
	double pnd=PND[0];
	static int cnt=0;

	fftw_complex **ck = (fftw_complex **) fftw_malloc( NRSITES *sizeof( fftw_complex )); /* [site][idx] */
	fftw_complex **tk = (fftw_complex **) fftw_malloc( NRSITES *sizeof( fftw_complex ));
	fftw_complex *xx_ux = (fftw_complex *) fftw_malloc( NNN *sizeof( fftw_complex ));

		for( m=0; m<=NRSITES-1; m++){	
			*( ck + m ) = (fftw_complex *) fftw_malloc( NNN *sizeof( fftw_complex ));
			*( tk + m ) = (fftw_complex *) fftw_malloc( NNN *sizeof( fftw_complex ));
		}

		/*ck_uh = cr_uo_s*/
		for( m=0; m<=NRSITES-1; m++)	
		    for( i=0; i<=nnn-1; i++){
			ck[m][i][0] = cr_s[m][i];	
			ck[m][i][1] = 0.00;
					}		

		/*FFT cr_uo_s and cr_uh_s*/
		for( m=0; m<=NRSITES-1; m++){	
			shift_origin_complex( *(ck+m), xx_ux, SYS );
			fftw_3d( xx_ux, *(ck+m) );
		}
		
		if( strncmp( "no", EWALD_SUMS,2) == 0)
		for( m=0; m<=NRSITES-1; m++)	
		    for( i=0; i<=nnn-1; i++){ 
			ck[m][i][0] = ck[m][i][0] - uk_l[m][i][0];
			ck[m][i][1] = ck[m][i][1] - uk_l[m][i][1];
		    }

	/****OZ EQUATION***********************************************************************************/

		if( TYPE == 1 ){ 
		
			for( i=0; i<=nnn-1; i++)
			    for( rc=0; rc<=1; rc++){
				tk[0][i][rc]  =  ck[0][i][rc] *(pnd*hk3_oo[i]) + 2*ck[1][i][rc] *(wk_oh[i] + pnd*hk3_oh[i]);	/*pch*/
		
				tk[1][i][rc]  =  ck[1][i][rc] *(wk_hh[i] + 2*pnd*hk3_hh[i] );		/*cw -c*/
				tk[1][i][rc] +=  ck[0][i][rc] *(wk_oh[i] +   pnd*hk3_oh[i] ); /*pch*/
						}
		}

	/****END OF OZ*************************************************************************************/

			/*Subtract long range pot. int k space, tk -> tk_s*/
		if( strncmp( "no", EWALD_SUMS,2) == 0)
		for( m=0; m<=NRSITES-1; m++)	
		    for( i=0; i<=nnn-1; i++){
	    		tk[m][i][0] = tk[m][i][0] - uk_l[m][i][0];
			tk[m][i][1] = tk[m][i][1] - uk_l[m][i][1];
		    }

		/*Inv_FFT and unshift*/
		/*tk_uo = tr_uo_s*/
		for( m=0; m<=NRSITES-1; m++){	
			invfftw_3d( *(tk+m), xx_ux);	
			unshift_origin_complex( xx_ux, *(tk+m), SYS);
		}
	
		for( m=0; m<=NRSITES-1; m++)	
		    for( i=0; i<=NNN-1; i++)  tr_s[m][i] = tk[m][i][0];
				
		/*HNC closure*/
		for( m=0; m<=NRSITES-1; m++)	
		    for( i=0; i<=NNN-1; i++)
			cr2_s[m][i] = exp( -ur_s[m][i] + tr_s[m][i] ) - tr_s[m][i] -1;

		for( m=0; m<=NRSITES-1; m++){	
			free( *(ck+m) );
			free( *(tk+m) );
		}

	free( ck );
	free( tk );
	free( xx_ux );

}

/************************************************************************************************/ 
/************************************************************************************************/ 
/************************************************************************************************/ 
/************************************************************************************************/ 
/************************************************************************************************/ 
/*												*/
/*					Newton_Gmres_ML_iter					*/
/*												*/
/************************************************************************************************/ 
/************************************************************************************************/
/************************************************************************************************/ 

/***************************************/
struct basis_range { int bx; int ex; 
		     int by; int ey; 
		     int bz; int ez; };
typedef struct basis_range BR;
int B_SIZE[3];
double *B_FUNC[3];
double *B3D_FUNC;
BR *B_RANGE;
double *B_NORM;
double *O_LAP;
/**************************************/
#define ibn(x,y,z) B_SIZE[2]*B_SIZE[1] *(x) + B_SIZE[2] * (y) + (z)

double *A_XR_OO, *A_XR_OH,  *A_XR_HH;
double *A_HR;

/**************************************/
void calc_XR_jac_constants( void );
void calc_basis_set_info( int [], int [], int [] );
void transfer_basis_to_fine( double *, double *, double * );
void reduce_to_A_and_q( double *, double *, double *);
double * calc_basis_coef( double * ); /*n_b, n, b_sp, cr*/
void gmres_iter(double *, double *, double  *, int ,  int  );
void print_basis_set_info( int [] , int [], int [] );


void newton_gmres_ml_iter( double *cr2_uo_s, double *cr2_uh_s, double *cr_uo_s, double *cr_uh_s, double *tr_uo_s, double *tr_uh_s )
{
	/*EXTERN*/
	int n[3]={NX, NY, NZ}, nx=NX, ny=NY, nz=NZ, n_b[3]={NX_B, NY_B, NZ_B};
	int n_sp[3]={NX_SP, NY_SP, NZ_SP}, u_x=CX, u_y=CY, u_z=CZ;
	double pres = PND[0], mp = PIC_MP;
	double l[3]={LX,LY,LZ};
	double *ur_uo_s = *(UR_S +0); 
	double *ur_uh_s = *(UR_S +1);
	/*extern*/

	int i, j, nnn= n[0]*n[1]*n[2], nnn_b=n_b[0]*n_b[1]*n_b[2];
	double tmp1=0.00;

	static int cnt = 0;
	/*basis set variables*/						printf("\nCalculating basis set information\n");fflush(stdout);
	if( cnt == 0 ){
			calc_basis_set_info( n_b, n_sp, n );
			print_basis_set_info(n_b, n_sp, n);
	 	      }

			print_2d( "3d_basis.2d.dat", B3D_FUNC, B_SIZE[0], B_SIZE[1], B_SIZE[2], B_SIZE[2]/2 );

	/*Calc c(r) components*/					printf("Calculating A and q for cr \n"); fflush(stdout);
	double *A_cr_uo = (double *) malloc( nnn_b *sizeof(double));
	double *A_cr_uh = (double *) malloc( nnn_b *sizeof(double));
	double *q_cr_uo = (double *) malloc( nnn *sizeof(double));
	double *q_cr_uh = (double *) malloc( nnn *sizeof(double));
		reduce_to_A_and_q( A_cr_uo, q_cr_uo, cr_uo_s );
		reduce_to_A_and_q( A_cr_uh, q_cr_uh, cr_uh_s );

		   	   print_2d(  "cr_uo_s.2d.dat",  cr_uo_s, CZ, nx, ny, nz );
			   print_2d( "cr2_uo_s_pic.2d.dat", cr2_uo_s, CZ, nx, ny, nz );
			   print_2d( "A_cr_uo.2d.dat", A_cr_uo, n_b[0], n_b[1], n_b[2], (n_b[2]/2));		
	   		   print_2d( "q_cr_uo.2d.dat", q_cr_uo, CZ, nx, ny, nz );

	/*Calc f(r) components*/					printf("Calculating A and q for fr\n"); fflush(stdout);
	double *fr_uo = (double *) malloc( nnn *sizeof(double));
	double *fr_uh = (double *) malloc( nnn *sizeof(double));
	double *q_fr_uo = (double *) malloc( nnn *sizeof(double));
	double *q_fr_uh = (double *) malloc( nnn *sizeof(double));
	double *A_fr_uo = (double *) malloc( nnn_b *sizeof(double));
	double *A_fr_uh = (double *) malloc( nnn_b *sizeof(double));
		for( i=0; i<=nnn-1; i++) fr_uo[i] = cr2_uo_s[i] - cr_uo_s[i];
		for( i=0; i<=nnn-1; i++) fr_uh[i] = cr2_uh_s[i] - cr_uh_s[i];
		reduce_to_A_and_q( A_fr_uo, q_fr_uo, fr_uo );
		reduce_to_A_and_q( A_fr_uh, q_fr_uh, fr_uh );
	double *A_fr = (double *) malloc( 2 *nnn_b *sizeof(double));
		for( i=0; i<=nnn_b-1; i++) 	A_fr[i] =  A_fr_uo[i];
		for( i=0; i<=nnn_b-1; i++) 	A_fr[nnn_b+i] =  A_fr_uh[i];

			   print_2d( "fr_uo.2d.dat", fr_uo, CZ, nx, ny, nz );
			   print_2d( "A_fr_uo.2d.dat", A_fr_uo, n_b[0], n_b[1], n_b[2], (n_b[2]/2) );
			   print_2d( "q_fr_uo.2d.dat", q_fr_uo, CZ, nx, ny, nz );

	/*Jacobian Elements __hr__*/						printf("Calculating A and q for hr\n"); fflush(stdout);
	double *hr_uo = (double *) malloc( nnn *sizeof(double));
	double *hr_uh = (double *) malloc( nnn *sizeof(double));
	double *q_hr_uo = (double *) malloc( nnn *sizeof(double));
	double *q_hr_uh = (double *) malloc( nnn *sizeof(double));
	double *A_hr_uo = (double *) malloc( nnn_b *sizeof(double));
	double *A_hr_uh = (double *) malloc( nnn_b *sizeof(double));
		for( i=0; i<= nnn-1; i++) hr_uo[i] = cr2_uo_s[i] + tr_uo_s[i];
		for( i=0; i<= nnn-1; i++) hr_uh[i] = cr2_uh_s[i] + tr_uh_s[i];
		reduce_to_A_and_q( A_hr_uo, q_hr_uo, hr_uo );
		reduce_to_A_and_q( A_hr_uh, q_hr_uh, hr_uh );
	double *A_hr = (double *)malloc( 2*nnn_b *sizeof( double ));
		for( i=0; i<=nnn_b-1; i++) A_hr[i] = A_hr_uo[i];
		for( i=0; i<=nnn_b-1; i++) A_hr[ nnn_b +i ] = A_hr_uh[i];


		/*adding the b_norm term to A_hr*/
		for( i=0; i<=(2*nnn_b)-1; i++) A_hr[i] *= (*B_NORM);
	A_HR = A_hr;
			   print_2d( "hr_uo.2d.dat", hr_uo, CZ, nx,ny,nz);
			   print_2d( "A_hr_uo.2d.dat", A_hr_uo, n_b[0], n_b[1], n_b[2], (n_b[2]/2));
			   print_2d( "q_hr_uo.2d.dat", q_hr_uo, CZ, nx, ny, nz);

							printf("Calculating overlap and XR values\n"); fflush(stdout);

	if( cnt == 0){	calc_XR_jac_constants();    }

	double *xA_cr = (double *) malloc( 2 *nnn_b *sizeof(double));
		for( i=0; i<= (2*nnn_b)-1; i++)	  xA_cr[ i ] = 0.00;

									printf("Calculating step length..."); fflush(stdout);
		double *Jac;
#ifdef JACOBIAN	
		Jac = calc_jacobian( tr_uo_s, tr_uh_s );
		gmres_iter( xA_cr, Jac, A_fr, Max_Gmres_Iter, 2*nnn_b );
 #else
		gmres_iter( xA_cr, Jac, A_fr, Max_Gmres_Iter, 2*nnn_b);
  #endif

		
		printf("...Done!!!\n"); fflush(stdout);
	double *xA_cr_uo = (double *) malloc( nnn_b *sizeof(double));
	double *xA_cr_uh = (double *) malloc( nnn_b *sizeof(double));
	for( i=0; i<=nnn_b-1; i++) {
		xA_cr_uo[i] = xA_cr[i];  
		xA_cr_uh[i] = xA_cr[i+nnn_b]; 
				   }
	print_2d( "A_xcr_uo.2d.dat", xA_cr_uo, n_b[0], n_b[1], n_b[2], (n_b[2]/2));
	print_2d( "A_xcr_uh.2d.dat", xA_cr_uh, n_b[0], n_b[1], n_b[2], (n_b[2]/2));

	double *A2_cr_uo = (double *) malloc( nnn_b *sizeof(double));
	double *A2_cr_uh = (double *) malloc( nnn_b *sizeof(double));

	for( i=0; i<=nnn_b-1; i++)	A2_cr_uo[i] = A_cr_uo[i] +  xA_cr_uo[i];
	for( i=0; i<=nnn_b-1; i++)	A2_cr_uh[i] = A_cr_uh[i] +  xA_cr_uh[i];

	print_2d( "A2_cr_uo.2d.dat", A2_cr_uo, n_b[0], n_b[1], n_b[2], (n_b[2]/2));

									printf("Transfering back to fine grid...");
	transfer_basis_to_fine( cr2_uo_s, A2_cr_uo, q_cr_uo );
	transfer_basis_to_fine( cr2_uh_s, A2_cr_uh, q_cr_uh );

/*	for( i=0; i<=nnn-1; i++)  cr2_uo_s[i] += q_cr_uo[i];*/
/*	for( i=0; i<=nnn-1; i++)  cr2_uh_s[i] += q_cr_uh[i];*/
									printf("...Done!!!\n"); fflush(stdout);
	print_2d( "cr2_uo_s.2d.dat", cr2_uo_s, CZ, nx, ny, nz );
			
	cnt++;

	free( A_cr_uo ); free( q_cr_uo ); 
	free( A_cr_uh ); free( q_cr_uh );

	free( fr_uo ); free( A_fr_uo ); free( q_fr_uo ); 
	free( fr_uh ); free( A_fr_uh ); free( q_fr_uh );

	free( hr_uo ); free( A_hr_uo ); free( q_hr_uo );
	free( hr_uh ); free( A_hr_uh ); free( q_hr_uh );

 	 free( xA_cr );
	free( A_fr ); free( A2_cr_uo ); free( xA_cr_uo ); 
	free( A_hr ); free( A2_cr_uh ); free( xA_cr_uh );
	
	
#ifdef JACOBIAN
	if( Jac != NULL )  free( Jac );
#endif



}

/****************************************************************************************/
/*				End of newton_gmres_ml routine				*/
/****************************************************************************************/ 


void calc_XR_jac_constants(void)
{
	/*EXTERN*/
	double pd=PND[0], l[3]={LX,LY,LZ};
	int n[3]={NX,NY,NZ}, n_b[3]={NX_B, NY_B, NZ_B}, nx=NX, ny=NY, nz=NZ;
	double *hk3_oo = HK3[0][0];
	double *hk3_oh = HK3[0][1];
	double *hk3_hh = HK3[1][1];
	double  *wk_oh = WK_OH,   *wk_hh =  WK_HH;
	double b_norm = *B_NORM;
	/*extern*/

	int i;
	int nnn = n[0]*n[1]*n[2];
	int nnn_b = n_b[0]*n_b[1]*n_b[2];

	double dl3=1.00;
		for( i=0; i>=2; i++)  dl3 *= (l[i]/(n[i] -1));
	
	double cnst1 = dl3;
	/**/
	double *hr3 = (double *) malloc( nnn *sizeof(double));
	fftw_complex *xx_ux = (fftw_complex *)malloc( nnn *sizeof(fftw_complex));
	fftw_complex *hc3   = (fftw_complex *)malloc( nnn *sizeof(fftw_complex));
	/**/

	/*calc A_X_oo*/
	for( i=0; i<=nnn-1; i++){  hc3[i][0] = pd *hk3_oo[i];	hc3[i][1] = 0.00;}
	  invfftw_3d(hc3, xx_ux);	
	    unshift_origin_complex( xx_ux, hc3, SYS );
	      for( i=0; i<=nnn-1; i++)  hr3[i] = hc3[i][0];		print_2d( "xr_oo.2d.dat", hr3, n[2]/2, nx,ny,nz );
	      double *A_oo = calc_basis_coef( hr3 );
		for( i=0; i<=nnn_b-1; i++)  A_oo[i] *= dl3;	
	A_XR_OO = A_oo;
									print_2d( "A_xr_oo.2d.dat", A_XR_OO, n_b[0], n_b[1], n_b[2], n_b[2]/2 );
	/*calc A_X_oh*/
	for( i=0; i<=nnn-1; i++){  hc3[i][0] = wk_oh[i] + pd *hk3_oh[i];   hc3[i][1] = 0.00; }
	  invfftw_3d(hc3, xx_ux);	
	    unshift_origin_complex( xx_ux, hc3, SYS );
	      for( i=0; i<=nnn-1; i++)   hr3[i] = hc3[i][0];		print_2d( "xr_oh.2d.dat", hr3, n[2]/2, nx, ny, nz );
	      double *A_oh = calc_basis_coef( hr3 );
		for( i=0; i<=nnn_b-1; i++)  A_oh[i] *= dl3;	
	A_XR_OH = A_oh;
									print_2d( "A_xr_oh.2d.dat", A_XR_OH, n_b[0], n_b[1], n_b[2], n_b[2]/2 );

	/*calc A_X_hh*/
	for( i=0; i<=nnn-1; i++){  hc3[i][0] = wk_hh[i]  + 2 *pd *hk3_hh[i];   hc3[i][1] = 0.00; }
	  invfftw_3d(hc3, xx_ux);	
	    unshift_origin_complex( xx_ux, hc3, SYS );
	      for( i=0; i<=nnn-1; i++)   hr3[i] = hc3[i][0];		print_2d( "xr_hh.2d.dat", hr3, n[2]/2, nx,ny,nz );
	      double *A_hh = calc_basis_coef( hr3 );
		for( i=0; i<=nnn_b-1; i++)  A_hh[i] *= dl3;	
	A_XR_HH = A_hh;
									print_2d( "A_xr_hh.2d.dat", A_XR_HH, n_b[0], n_b[1], n_b[2], n_b[2]/2 );

	/**/
	free( hr3 );
	free( xx_ux );
	free( hc3 );
	/**/
}

/****************************************************************************************/	
/*											*/
/*				Calc course grid values					*/
/*											*/
/****************************************************************************************/ 


void  reduce_to_A_and_q( double *A_cr, double *q_cr, double *cr )
{
	/*EXTERN*/
	int n[3]={NX,NY,NZ}, n_b[3]={NX_B,NY_B,NZ_B};
	BR *b_r = B_RANGE;
	int b_size[3];	 b_size[0]=B_SIZE[0]; b_size[1]=B_SIZE[1]; b_size[2]=B_SIZE[2];
	double *b3d_func = B3D_FUNC;
	double b_norm = *B_NORM;
	/*\extern*/

	int i, x, y, z, a, b, c, id1, id2, bx, by, bz;
	double norm=0, tmp1, tmp2;
	int nnn, nnn_b;  
		nnn = n[0]*n[1]*n[2];
		nnn_b = n_b[0] *n_b[1] *n_b[2];
	int b3d_s = b_size[0]*b_size[1]*b_size[2];

	for( i=0; i<=nnn-1; i++)	q_cr[i] = cr[i];

	for( a=0; a<=n_b[0]-1; a++)
	for( b=0; b<=n_b[1]-1; b++)
	for( c=0; c<=n_b[2]-1; c++)
	{
		id1 = ib(a,b,c);
		bx = b_r[ id1 ].bx;
		by = b_r[ id1 ].by;
		bz = b_r[ id1 ].bz;
		A_cr[ id1 ] = 0.00;
		tmp1 = 0.00;

		for( x=0; x<=b_size[0]-1; x++)
		for( y=0; y<=b_size[1]-1; y++)
		for( z=0; z<=b_size[2]-1; z++)
			tmp1 +=  cr[ ii( x+bx, y+by, z+bz ) ] *b3d_func[ ibn(x,y,z) ];
		
		A_cr[ id1 ] = tmp1/b_norm;
				
		for( x=0; x<=b_size[0]-1; x++)
		for( y=0; y<=b_size[1]-1; y++)
		for( z=0; z<=b_size[2]-1; z++)
			q_cr[ ii(x+bx,y+by,z+bz) ] -= (A_cr[ id1 ] * b3d_func[ ibn(x,y,z) ]);
	}
}


/****************************************************************************************/	
/*											*/
/*				Transfer basis coef to fine grid			*/
/*											*/
/****************************************************************************************/ 


void transfer_basis_to_fine( double *cr, double *A_cr, double *fd_cr )
{
	/*EXTERN*/
	int n_b[3]={NX_B,NY_B,NZ_B};
	int b_size[3] = {B_SIZE[0],B_SIZE[1],B_SIZE[2]};
	double *b3d_func = B3D_FUNC;	/*	b3d_func = calc_3d_basis_func( b_size, b_func ); */
	BR *b_r = B_RANGE;
	double *b_norm = B_NORM;
	/*extern*/
	
	int i, j, a, b, c, x, y, z, indx1,  bx, by, bz;
	int  nnn=NX*NY*NZ;

	for( i=0; i<=nnn-1; i++)  cr[i] = 0.00;

	for( a=0; a<=n_b[0]-1; a++)
	for( b=0; b<=n_b[1]-1; b++)
	for( c=0; c<=n_b[2]-1; c++)
	{
			indx1 = ib(a,b,c);
			bx = b_r[indx1].bx;
			by = b_r[indx1].by;
			bz = b_r[indx1].bz;
		for( x=0; x<=b_size[0]-1; x++)
		for( y=0; y<=b_size[1]-1; y++)
		for( z=0; z<=b_size[2]-1; z++)
			 	cr[ ii( x+bx, y+by, z+bz ) ] += b3d_func[ibn(x,y,z) ] * A_cr[ indx1 ];
	}

	for( i=0; i<=nnn-1; i++) cr[i] += fd_cr[i];
}


/****************************************************************************************/ 
/*											*/
/*			Routines to calculate the basis set functions 			*/
/*											*/
/****************************************************************************************/ 

double * calc_basis_coef_sub( int [], double * , BR * , int [], double *  );

double * calc_basis_coef( double * cr )
{
	int i;
	/*EXTERN*/
	int n[3]={NX,NY,NZ}, n_b[3]={NX_B,NY_B,NZ_B};
	int b_size[3]={B_SIZE[0],B_SIZE[1],B_SIZE[2]}; 
	double *b_func[3];	b_func[0] = B_FUNC[0]; b_func[1] = B_FUNC[1]; b_func[1] = B_FUNC[1];
	double *b3d_func = B3D_FUNC;
	BR *b_range = B_RANGE;
	/*\EXTERN*/

	double *b_coef = calc_basis_coef_sub( n_b, b3d_func, b_range, b_size, cr );
	return b_coef;
}

double * calc_basis_coef_sub( int n_b[], double * b3d_func, BR * b_r, int b_size[], double * cr )
{

	int i, x, y, z, a, b, c, indx1, indx2, bx, by, bz;
	double norm=0;
	int nnn_b = n_b[0] *n_b[1] *n_b[2];
	int b3d_s = b_size[0]*b_size[1]*b_size[2];
	double *b_coef = (double *) malloc( nnn_b *sizeof(double)); 

	for( i=0; i<=b3d_s-1; i++)	    norm += b3d_func[i];

	for( a=0; a<=n_b[0]-1; a++)
	    for( b=0; b<=n_b[1]-1; b++)
	        for( c=0; c<=n_b[2]-1; c++)
		{
			indx1 = ib(a,b,c);
			bx = b_r[indx1].bx;
			by = b_r[indx1].by;
			bz = b_r[indx1].bz;
			b_coef[ indx1 ] = 0.00;
			for( x=0; x<=b_size[0]-1; x++)
			    for( y=0; y<=b_size[1]-1; y++)
				for( z=0; z<=b_size[2]-1; z++)
					b_coef[ indx1 ] += cr[ ii( x+bx, y+by, z+bz ) ] *b3d_func[ibn(x,y,z) ];
			b_coef[indx1] /= norm;
		}
	return b_coef;
}




/****************************************************************************************/ 
/*											*/
/*				calc basis set informations				*/
/*											*/
/****************************************************************************************/ 
int calc_basis_size( int , int ); 
BR * calc_basis_range( int [], int [], int [], int [] );
double * calc_1d_basis_func( int );
double * calc_3d_basis_func( int [] , double * [] );

void calc_basis_set_info( int n_b[], int n_sp[] , int n[])
{
	int i;

	int b_size[3];
		b_size[0] = calc_basis_size( n_b[0], n_sp[0] ); 
		b_size[1] = calc_basis_size( n_b[1], n_sp[1] ); 
		b_size[2] = calc_basis_size( n_b[2], n_sp[2] );

	int nb_pts = (b_size[0]*b_size[1]*b_size[2]);

	double *b_func[3];
		b_func[0] = calc_1d_basis_func( b_size[0] );
		b_func[1] = calc_1d_basis_func( b_size[1] );
		b_func[2] = calc_1d_basis_func( b_size[2] );

	double *b3d_func = calc_3d_basis_func( b_size, b_func ); 

	BR *b_range = calc_basis_range( b_size, n_b, n_sp, n );

	double *b_norm = (double *) malloc( sizeof(double));
	double tmp=0.00;

	for( i=0; i<=nb_pts-1; i++)	tmp += b3d_func[i];
	(*b_norm) = tmp;

	/*Set Global variables*/
	B_SIZE[0] = b_size[0];  B_SIZE[1] = b_size[1];	B_SIZE[2] = b_size[2];
	B_FUNC[0] = b_func[0];	B_FUNC[1] = b_func[1]; 	B_FUNC[2] = b_func[2];
	B3D_FUNC = b3d_func;
	B_RANGE = b_range;
	B_NORM = b_norm;

}

#ifdef FLAT_TOP
int calc_basis_size( int n_b, int n_sp )
{
	int i, tmp, b_size;
	b_size  = (int) ( n_sp/n_b ); 		/* b_size + (n_b-1)*(b_size-1)/2 <= N */
	if( (b_size%2) == 0 )	b_size--;
	return b_size;
}

double * calc_1d_basis_func( int b_size )
{
	int i, mid;

	double inc;
	double *basis_func = (double *) malloc( b_size *sizeof(double));
	for( i=0; i<=b_size-1; i++)	basis_func[i] = 1.00;
	return basis_func;
}



double * calc_3d_basis_func( int b_size[] , double *b_func[] )
{
	int i, x, y, z, bt, indx;
		bt = b_size[0] *b_size[1] *b_size[2];
	double *b3d_func = (double *) malloc( bt *sizeof(double));

	for( z=0; z<= b_size[2] - 1; z++)
	    for( y=0; y<= b_size[1] - 1; y++)
		for( x=0; x<= b_size[0] -1; x++){
			indx = b_size[2]*b_size[1]*(x) + b_size[2]*(y) + (z);
			b3d_func[ indx ] = b_func[0][x] * b_func[1][y] * b_func[2][z];
						}
	return b3d_func;
}


BR * calc_basis_range( int b_size[], int n_b[], int n_sp[],  int n[] )
{
	int i, x, y, z, a, indx, n_bsp[3], n_rmd[3], ib_pnt[3], ie_pnt[3];
	int nnn_b = n_b[0]*n_b[1]*n_b[2];
	BR *b_r;
	b_r = (BR *) malloc( nnn_b *sizeof( BR )); 

	/*this calculates the actual number of full space pnts that are represented by the basis functions*/
	n_bsp[0] = b_size[0] * n_b[0] ;
	n_bsp[1] = b_size[1] * n_b[1] ;
	n_bsp[2] = b_size[2] * n_b[2] ;
	printf( "\n::Number of pts actual represented (x,y,z): (%d,%d,%d)\n", n_bsp[0], n_bsp[1], n_bsp[2] ); fflush(stdout);
	n_rmd[0] = n[0] - n_bsp[0];
	n_rmd[1] = n[1] - n_bsp[1];
	n_rmd[2] = n[2] - n_bsp[2];
	printf( "\n::Difference in the Number of pts actually represented (x,y,z): (%d,%d,%d)\n", n_rmd[0], n_rmd[1], n_rmd[2] ); fflush(stdout);
	ib_pnt[0] = (n_rmd[0] +1) / 2;
	ib_pnt[1] = (n_rmd[1] +1) / 2;
	ib_pnt[2] = (n_rmd[2] +1) / 2;
	printf( "\n::beginning pts for basis functions (x,y,z) : (%d, %d, %d)\n", ib_pnt[0], ib_pnt[1], ib_pnt[2] ); fflush(stdout);

	for( x=0; x<=n_b[0]-1; x++)
	    for( y=0; y<=n_b[1]-1; y++)
		for( z=0; z<=n_b[2]-1; z++)
		{
			indx = n_b[2]*n_b[1] *x + n_b[2] *y + z;
			b_r[indx].bx = ib_pnt[0] + x*b_size[0];
			b_r[indx].by = ib_pnt[1] + y*b_size[1];
			b_r[indx].bz = ib_pnt[2] + z*b_size[2];
			b_r[indx].ex = b_r[indx].bx + (b_size[0] -1);
			b_r[indx].ey = b_r[indx].by + (b_size[1] -1);
			b_r[indx].ez = b_r[indx].bz + (b_size[2] -1);
		}

	return b_r;
}

#endif

void print_basis_set_info( int n_b[], int n_sp[] , int n[] )
{
	int i, j, x, y, z;
	FILE *out;
	if( (out = fopen("basis_set_info.dat","w" )) == NULL )
		printf("Error opening file\n"); fflush( stdout );

	fprintf(out, "\n............BASIS SET INFO.............\n\n"); fflush(stdout);
	fprintf(out, "Number of basis functions\n");
		for( i=0; i<=2; i++)
			fprintf(out, "n_b[%d] : %d\n", i,n_b[i]);

	fprintf(out, "\nNumber of points spanned by basis function\n"); fflush(stdout);
		for( i=0; i<=2; i++)
			fprintf(out, "b_sp[%d] : %d\n", i,n_sp[i]);
		
	fprintf(out, "\nSize of basis function\n");
		for( i=0; i<=2; i++)
			fprintf(out, "b_size[%d] : %d\n", i,B_SIZE[i]);

	fprintf(out, "\nValues of 1D basis functions\n");
		for( i=0; i<=2; i++){
			fprintf(out, "\nbasis_func[%d]:\n", i );
			for( j=0; j<=B_SIZE[i]-1; j++)
				fprintf(out, "%.10f\n", B_FUNC[i][j]);
				    }

	fprintf(out, "\nNormalization constant: %lf\n", *B_NORM);

	fprintf(out, "\nRanges of basis functions\n"); fflush( stdout );
		for( x=0; x<=n_b[0]-1; x++)
	 	 for( y=0; y<=n_b[1]-1; y++)
		  for( z=0; z<=n_b[2]-1; z++)
			fprintf(out, "[%d,%d,%d]....\n\t%d->%d\n\t%d->%d\n\t%d->%d\n", x,y,z, 
					B_RANGE[ib(x,y,z)].bx, B_RANGE[ ib(x,y,z) ].ex,
					B_RANGE[ib(x,y,z)].by, B_RANGE[ ib(x,y,z) ].ey,
					B_RANGE[ib(x,y,z)].bz, B_RANGE[ ib(x,y,z) ].ez );

	fclose( out );
}


/****************************************************************************************/ 

/****************************************************************************************/
/*                     		          GMRESIter 					*/
/****************************************************************************************/
/*                                                                         		*/
/****************************************************************************************/

double * jh_Mul_MV( double *, double *, int, int );
double * jh_Mul_JV( double *, int, int);
		/*       M, V, n_rows, n_cols */

#define RTCEps 1e-8
#include "nrutil.h"


void gmres_iter(double *x, double *J, double  *b, int MaxIter,  int ni )
		/*vector x,mat_vec J,  vector b ,  Max iterations (total),  , size of A x  b*/
{

    	char vName[10];
    	int Iter, i, j, k;
    	double h1, h2, r;
    	double bNorm;
    	double **h, *y, *s, *c1, *c2;
    	double *vtmp1, *vtmp2;
    	int Dim;
    	double a1, a2;
    	double **v;

	Dim = ni;

    	/* allocation of matrix H and vectors y, s, c1 and c2 */
    	h = (double **)malloc((GMRESSteps + 1) * sizeof(double *));
        	for (i = 1; i <= GMRESSteps; i++) 
           		h[i] = (double *)malloc((GMRESSteps + 2) * sizeof(double));
    	y = (double *)malloc((GMRESSteps + 1) * sizeof(double));
    	s = (double *)malloc((GMRESSteps + 2) * sizeof(double));
    	c1 = (double *)malloc((GMRESSteps + 1) * sizeof(double));
    	c2 = (double *)malloc((GMRESSteps + 1) * sizeof(double));
	/*    vtmp1 = (double *)malloc( Dim * sizeof(double));*/
    	vtmp2 = (double *)malloc( Dim * sizeof(double));

    	/* ... and vectors u */   
    	v = malloc((GMRESSteps + 2) * sizeof(double));
       		for (i = 1; i <= GMRESSteps + 1; i++) 
	    		*(v+i) = malloc( Dim * sizeof(double));

        bNorm = vector_norm(b, Dim);
        /* loop for 'MaxIter' GMRES cycles */
        Iter = 0;
        /* v[1] = r = b - A * x(i) */
		a1 = vector_sum(x, Dim);
        if (!IsZero( a1 / Dim)) {
		#ifdef JACOBIAN
			vtmp1 = jh_Mul_MV(J, x, Dim, Dim );
		 #else
			vtmp1 = jh_Mul_JV( x, Dim, Dim );
		  #endif
	    	Sub_VV(vtmp2, b, vtmp1, Dim);
		free( vtmp1);
		vtmp1 = NULL;
            	Asgn_VV(*(v+1), vtmp2, Dim);
        } else {
            Asgn_VV(*(v+1), b, Dim);
	}
        while ( (a2=vector_norm(*(v+1), Dim)) > RTCEps * bNorm  && Iter < MaxIter) {/*out*/
	
            	s[1] = vector_norm(*(v+1), Dim);
            	MulAsgn_VS(*(v+1), 1.0 / s[1], Dim);

           /* GMRES iteration */
            	i = 0;
            	while (( fabs(s[i+1]) > RTCEps * bNorm) && i < GMRESSteps && Iter < MaxIter) 
		{
 	        	i++;
			Iter++;
                	/* w = v[i+1] */
			#ifdef JACOBIAN
				vtmp1 = jh_Mul_MV( J, *(v+i), Dim, Dim );
			 #else
				vtmp1 = jh_Mul_JV( *(v+i), Dim, Dim );
			  #endif
			Asgn_VV(*(v+i+1), vtmp1, Dim);
			free( vtmp1);
			vtmp1 = NULL;

                	/* modified Gram-Schmidt orthogonalization */
                	for (k = 1; k <= i; k++) {
                    		h[i][k] = (double ) Mul_VV(*(v+i+1), *(v+k), Dim);
				vtmp1 = (double *)Mul_SV(h[i][k], *(v+k), Dim);
                    		SubAsgn_VV(*(v+i+1), vtmp1, Dim);
				free(vtmp1);
				vtmp1 = NULL;
                				}

	                h[i][i+1] = l2Norm_V(*(v+i+1), Dim);
        	        MulAsgn_VS(*(v+i+1), 1.0 / h[i][i+1], Dim);

                	/* Q-R algorithm */
                	for (k = 1; k < i; k++) {
                    		h1 = c1[k] * h[i][k] + c2[k] * h[i][k+1];
                   		h2 = - c2[k] * h[i][k] + c1[k] * h[i][k+1];
                    		h[i][k] = h1;
                    		h[i][k+1] = h2;
                				}

                	r = sqrt(h[i][i] * h[i][i] + h[i][i+1] * h[i][i+1]);
                	c1[i] = h[i][i] / r;
                	c2[i] = h[i][i+1] / r;
                	h[i][i] = r;
                	h[i][i+1] = 0.0;
                	s[i+1] = - c2[i] * s[i];
               		s[i] = c1[i] * s[i];
            	}

            	/* Solving of the system of equations : H y = s */
            	for (j = i; j > 0; j--) {
                	y[j] = s[j] / h[j][j];
                	for (k = j - 1; k > 0; k--)
                    		s[k] -= h[j][k] * y[j];
            				}

            	/* updating solution */
            	for (j = i; j > 0; j--){
			vtmp1 = (double *) Mul_SV(y[j], *(v+j), Dim);
                	AddAsgn_VV(x, vtmp1, Dim);
			free( vtmp1);
			vtmp1 = NULL;
     					}

            	/* computing new residual */
		#ifdef JACOBIAN
			vtmp1 = jh_Mul_MV( J, x, Dim, Dim);
		 #else
			vtmp1 = jh_Mul_JV( x, Dim, Dim);
		  #endif			
	    	Sub_VV(vtmp2, b, vtmp1, Dim);
            	Asgn_VV(*(v+1), vtmp2, Dim);
		free( vtmp1 );
		vtmp1 = NULL;
        }/*out*/
    
    	/* release of vectors u, matrix H and vectors y, s, c1 and c2 */
    	if (v != NULL) {
        	for (i = 1; i <= GMRESSteps + 1; i++)
            		free(*(v+i));
        	free(v);
    			}

	if (h != NULL) {
        	for (i = 1; i <= GMRESSteps; i++)
	    		if (h[i] != NULL)
                		free(h[i]);
        	free(h);
    			}
 	if (y != NULL)  free(y);
    	if (s != NULL)  free(s);
    	if (c1 != NULL) free(c1);
    	if (c2 != NULL) free(c2);
  	if( vtmp2 != NULL) free( vtmp2 );

printf("# of iterations in gmres routine: %d\n", Iter); fflush(stdout);

}

/********************************************************************************/
/*										*/
/*			Gmres Subroutines					*/
/*										*/
/********************************************************************************/

double * jh_Mul_MV( double *M, double *V, int rows, int cols )
/*This is the same as a picard iteration, here J = Identity matrix*/
{
	int ir, ic;
	double tmp, cnst;
	double *MV = (double *) malloc( rows *sizeof(double));

#ifdef OMP
  #pragma omp parallel for num_threads(NUM_THREADS) private(ir,ic,tmp)
  
#endif
	for( ir=0; ir<=rows-1; ir++){
		tmp = 0.00;
	    	for( ic=0; ic<=cols-1; ic++)
			tmp+= M[ ir*cols +ic ] * V[ic];
		MV[ir] = tmp;
				    }
	return MV ;
}

double *jh_Mul_JV( double *V, int n_r, int n_c )
{

	/*Extern*/
	double den=PND[0];
	double b_norm= (*B_NORM);
	int n[3]={NX,NY,NZ}, n_b[3]={NX_B, NY_B, NZ_B};
	int nx_b = NX_B, ny_b = NY_B, nz_b = NZ_B;
	double *A_Xr_oo = A_XR_OO;	double *A_Xr_oh = A_XR_OH;	double *A_Xr_hh = A_XR_HH;
	double *o_lap = O_LAP;
	double *A_hr = A_HR;
	/*extern*/

	int ira, irb, ica, icb, i, id1;
	int dAx, dAy, dAz, id;
	int x,y,z, x1, y1, z1, a, b, di;
	int nnn_b=n_b[0]*n_b[1]*n_b[2];
	int b_0[3], sn_b=0.00;	
		for( i=0; i<=2; i++) b_0[i] = (n_b[i]-1)/2;
		for( i=0; i<=2; i++) sn_b += n_b[i];

	double *JV = (double *) malloc( n_r *sizeof( double ));
		for( i=0; i<=n_r-1; i++)
			JV[i] = 0.00;

	/*Calculate static constants*/
#ifdef OMP
  #pragma omp parallel for num_threads(NUM_THREADS) private(x,y,z,x1,y1,z1,ira,irb,ica,icb,dAx,dAy,dAz,id)
#endif
	for( x=0; x<=nx_b-1; x++)
	 for( y=0; y<=ny_b-1; y++)
	  for( z=0; z<=nz_b-1; z++)
	  {
		ira = ib(x,y,z); 
		irb = ira + nnn_b;

		for( x1=0; x1<=n_b[0]-1; x1++)
		 for( y1=0; y1<=n_b[1]-1; y1++)
		  for( z1=0; z1<=n_b[2]-1; z1++)
		  {
			ica = ib(x1,y1,z1);
			icb = ica + nnn_b;
			dAx = fabs(x-x1) + b_0[0];
			dAy = fabs(y-y1) + b_0[1]; 
			dAz = fabs(z-z1) + b_0[2];;
			
			/*new*/
			
			if( dAx<=n_b[0]-1 && dAy<=n_b[1]-1 && dAz<=n_b[2]-1 ) 
			{	
				id = ib(dAx,dAy,dAz);	

				JV[ ira ] += V[ica] *(   A_Xr_oo[id] );  /*   - 1.0*o_lap[di] );*/
				JV[ ira ] += V[icb] *( 2*A_Xr_oh[id] ); 
				JV[ irb ] += V[ica] *(   A_Xr_oh[id] ); 
				JV[ irb ] += V[icb] *(   A_Xr_hh[id]   ); /*- 1.0*o_lap[di] );*/
			}

		  }

		JV[ ira ] *= (1.0)*A_hr[ira];
		JV[ irb ] *= (1.0)*A_hr[irb];/*need a negative sign*/

		JV[ ira ] += (-1.0)*(b_norm)*V[ ira ];
		JV[ irb ] += (-1.0)*(b_norm)*V[ irb ];
	  }

	return JV;
}

/************************************************************************************************/	 
/*#ifdef OMP											*/
/*#pragma omp parallel for  num_threads(NUM_THREADS)   private(x,y,z,indx,i,k,xf,yf,zf,kr)	*/
/*#endif											*/
/************************************************************************************************/	 


double * XXXcalc_basis_coef( int n_b[], int n[], int b_sp[], double * cr )
{

	int i, j, k, x, y, z;
	int nnn = n[0]*n[1]*n[2];
	int nnn_b = n_b[0]*n_b[1]*n_b[2];

	int b_size[3];
		b_size[0] = calc_basis_size( n_b[0], b_sp[0] ); 
		b_size[1] = calc_basis_size( n_b[1], b_sp[1] ); 
		b_size[2] = calc_basis_size( n_b[2], b_sp[2] );

	double *b_func[3];
		b_func[0] = calc_1d_basis_func( b_size[0] );
		b_func[1] = calc_1d_basis_func( b_size[1] );
		b_func[2] = calc_1d_basis_func( b_size[2] );

	double *b3d_func;
		b3d_func = calc_3d_basis_func( b_size, b_func ); 
	BR *b_range;
		b_range = calc_basis_range( b_size, n_b, b_sp, n );
	double *b_coef;
		b_coef = calc_basis_coef_sub( n_b, b3d_func, b_range, b_size, cr );

	return b_coef;

}


/****************************************************************************************/ 
/*											*/
/*				Alternative basis functions				*/
/*											*/
/****************************************************************************************/ 

#ifdef ROOF_TOP

//int calc_basis_size( int n_b, int nn )
{
	int i, tmp, b_size;
	b_size  = 2 *(nn-1)/(n_b+1); 		/* b_size + (n_b-1)*(b_size-1)/2 <= N */
	if( (b_size%2) == 0 )	b_size++;
	return b_size;
}

//double * calc_1d_basis_func( int b_size )
{
	int i, mid;

	double inc;
	double *basis_func = (double *) malloc( b_size *sizeof(double));
	mid = (b_size+1) /2;
	inc = (double) 1.0/(mid-1);
	for( i=0; i<=mid-1; i++)	basis_func[i] = i*inc;
	for( i=0; i<=mid-2; i++)	basis_func[b_size-1-i] = basis_func[i];
	return basis_func;
}



//double * calc_3d_basis_func( int b_size[] , double *b_func[] )
{
	int i, x, y, z, bt, indx;
		bt = b_size[0] *b_size[1] *b_size[2];
	double *b3d_func = (double *) malloc( bt *sizeof(double));
	for( i=0; i<=bt-1; i++) b3d_func[i] = 0.00;

	for( z=0; z<= b_size[2] - 1; z++)
	    for( y=0; y<= b_size[1] - 1; y++)
		for( x=0; x<= b_size[0] -1; x++){
			indx = b_size[2]*b_size[1]*(x) + b_size[2]*(y) + (z);
			b3d_func[ indx ] = b_func[0][x] * b_func[1][y] * b_func[2][z];
						}
	return b3d_func;
}


//BR * calc_basis_range( int b_size[], int n_b[], int b_sp[],  int n[] )
{
	int i, x, y, z, a, indx, nb_pnts[3], n_rmd[3], ib_pnt[3], ie_pnt[3];
	int nnn_b = n_b[0]*n_b[1]*n_b[2];
	BR *b_r;
	b_r = (BR *) malloc( nnn_b *sizeof( BR )); 

	/*this calculates the actual number of full space pnts that are represented by the basis functions*/
	nb_pnts[0] = b_size[0] + (n_b[0] -1) * (b_size[0]-1)/2;
	nb_pnts[1] = b_size[1] + (n_b[1] -1) * (b_size[1]-1)/2;
	nb_pnts[2] = b_size[2] + (n_b[2] -1) * (b_size[2]-1)/2;
	printf( "\n::Number of pts actual represented (x,y,z): (%d,%d,%d)\n", nb_pnts[0], nb_pnts[1], nb_pnts[2] ); fflush(stdout);
	n_rmd[0] = n[0] - nb_pnts[0];
	n_rmd[1] = n[1] - nb_pnts[1];
	n_rmd[2] = n[2] - nb_pnts[2];
	printf( "\n::Difference in the Number of pts actually represented (x,y,z): (%d,%d,%d)\n", nb_pnts[0], nb_pnts[1], nb_pnts[2] ); fflush(stdout);
	ib_pnt[0] = (n_rmd[0] +1) / 2;
	ib_pnt[1] = (n_rmd[1] +1) / 2;
	ib_pnt[2] = (n_rmd[2] +1) / 2;
	printf( "\n::beginning pts for basis functions (x,y,z) : (%d, %d, %d)\n", ib_pnt[0], ib_pnt[1], ib_pnt[2] ); fflush(stdout);

	for( x=0; x<=n_b[0]-1; x++)
	    for( y=0; y<=n_b[1]-1; y++)
		for( z=0; z<=n_b[2]-1; z++)
		{
			indx = n_b[2]*n_b[1] *x + n_b[2] *y + z;
			b_r[indx].bx = ib_pnt[0] + x*(b_size[0]-1)/2;
			b_r[indx].by = ib_pnt[1] + y*(b_size[1]-1)/2;
			b_r[indx].bz = ib_pnt[2] + z*(b_size[2]-1)/2;
			b_r[indx].ex = b_r[indx].bx + (b_size[0] -1);
			b_r[indx].ey = b_r[indx].by + (b_size[1] -1);
			b_r[indx].ez = b_r[indx].bz + (b_size[2] -1);
		}

	return b_r;
}

#endif




