#include <fftw3.h>
#include "jh_struct.h"

#define ii(x,y,z) nz*ny*(x) + nz*(y) + (z)

/*This file is for code to manipulate a grid vector or what ever*/




/****************************************************************************************/ 
/*											*/
/*				Shifting routines					*/
/*											*/
/****************************************************************************************/ 

/*Shifts origin and wraps data around for FT routines*/

void shift_origin_complex( fftw_complex *in, fftw_complex *out,  ENV_PAR sys )
{
	int nx=sys.nx, ny=sys.ny, nz=sys.nz;		/*EXTERN*/
	int cx=sys.cx, cy=sys.cy, cz=sys.cz;
	int nnn = nx*ny*nz;

	int x, y, z, i;

	fftw_complex *tmp = (fftw_complex *) malloc( nnn * sizeof( fftw_complex ));

		/* Shifting x-coordinate*/
		if( cx != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( y=0; y<= ny-1; y++){
		     	for( x=0; x<=(cx-1); x++){
		       		out[ ii( (x+(nx-cx)),y,z) ][0] = in[ii(x,y,z)][0];
		       		out[ ii( (x+(nx-cx)),y,z) ][1] = in[ii(x,y,z)][1];
			}
		        for( x=cx; x<=nx-1; x++){
		    	    out[ ii( (x-cx),y,z) ][0] = in[ii(x,y,z)][0];
		    	    out[ ii( (x-cx),y,z) ][1] = in[ii(x,y,z)][1];
			}
		   }
		else
		    for( i=0; i<=nnn-1; i++){
			out[i][0] = in[i][0];
			out[i][1] = in[i][1];
		    }

		/* Shifting y-coordinate*/
		if( cy != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( x=0; x<= nx-1; x++){
		  	for( y=0; y<=(cy-1); y++){
		    	    tmp[ ii( x,(y+(ny-cy)),z) ][0] = out[ii(x,y,z)][0];
		    	    tmp[ ii( x,(y+(ny-cy)),z) ][1] = out[ii(x,y,z)][1];
			}
		  	for( y=cy; y<=ny-1; y++){
		    	    tmp[ ii( x,(y-cy),z) ][0] = out[ii(x,y,z)][0];
		    	    tmp[ ii( x,(y-cy),z) ][1] = out[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			tmp[i][0] = out[i][0];
			tmp[i][1] = out[i][1];
		    }

		/* Shifting z-coordinate*/
		if( cz != 0 )
		    for( x=0; x<= nx-1; x++)
		    for( y=0; y<= ny-1; y++){
		  	for( z=0; z<=(cz-1); z++){
		    	    out[ ii( x, y, (z+(nz-cz))) ][0] = tmp[ii(x,y,z)][0];
		    	    out[ ii( x, y, (z+(nz-cz))) ][1] = tmp[ii(x,y,z)][1];
			}
		  	for( z=cz; z<=nz-1; z++){
		    	    out[ ii( x, y, (z-cz)) ][0] = tmp[ii(x,y,z)][0];
		    	    out[ ii( x, y, (z-cz)) ][1] = tmp[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			out[i][0] = tmp[i][0];
			out[i][1] = tmp[i][1];
		    }

		free(tmp);
}


/**/
void unshift_origin_complex( fftw_complex *in, fftw_complex *out, ENV_PAR sys )
{
	int nx=sys.nx, ny=sys.ny, nz=sys.nz;		/*EXTERN*/
	int cx=sys.cx, cy=sys.cy, cz=sys.cz;
	int nnn = nx*ny*nz;

	int x, y, z, i;

	fftw_complex *tmp = (fftw_complex *) malloc( nnn * sizeof( fftw_complex ));

		/* unShifting z-coordinate*/
		if( cz != 0 )
		    for( x=0; x<= nx-1; x++)
		    for( y=0; y<= ny-1; y++){
		  	for( z=0; z<=(nz-1)-cz; z++){
		    	    out[ ii( x, y, (z+cz)) ][0] = in[ii(x,y,z)][0];
		    	    out[ ii( x, y, (z+cz)) ][1] = in[ii(x,y,z)][1];
			}
		  	for( z=nz-cz; z<=nz-1; z++){
		    	    out[ ii( x, y, (z-(nz-cz))) ][0] = in[ii(x,y,z)][0];
		    	    out[ ii( x, y, (z-(nz-cz))) ][1] = in[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			out[i][0] = in[i][0];
			out[i][1] = in[i][1];
		    }

		/* unShifting y-coordinate*/
		if( cy != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( x=0; x<= nx-1; x++){
		  	for( y=0; y<=(ny-1)-cy; y++){
		    	    tmp[ ii( x,(y+cy),z) ][0] = out[ii(x,y,z)][0];
		    	    tmp[ ii( x,(y+cy),z) ][1] = out[ii(x,y,z)][1];
			}
		  	for( y=ny-cy; y<=ny-1; y++){
		    	    tmp[ ii( x,(y-(ny-cy)),z) ][0] = out[ii(x,y,z)][0];
		    	    tmp[ ii( x,(y-(ny-cy)),z) ][1] = out[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			tmp[i][0] = out[i][0];
			tmp[i][1] = out[i][1];
		    }
	
		/* unShifting x-coordinate*/
		if( cx != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( y=0; y<= ny-1; y++){
		  	for( x=0; x<=(nx-1)-cx; x++){
		     	    out[ ii( (x+cx),y,z) ][0] = tmp[ii(x,y,z)][0];
		    	    out[ ii( (x+cx),y,z) ][1] = tmp[ii(x,y,z)][1];
			}
		  	for( x=nx-cx; x<=nx-1; x++){
		    	    out[ ii( (x-(nx-cx)),y,z) ][0] = tmp[ii(x,y,z)][0];
		    	    out[ ii( (x-(nx-cx)),y,z) ][1] = tmp[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			out[i][0] = tmp[i][0];
			out[i][1] = tmp[i][1];
		    }

	free(tmp);
}

void shift_origin_inplace( double *in, ENV_PAR sys )
{
	int nx=sys.nx, ny=sys.ny, nz=sys.nz;		/*EXTERN*/
	int cx=sys.cx, cy=sys.cy, cz=sys.cz;
	int nnn = nx*ny*nz;

	int x, y, z, i;

	double *tmp  = (double *) malloc( nnn *sizeof( double )); 
	double *tmp2 = (double *) malloc( nnn *sizeof( double )); 

		/* Shifting x-coordinate*/
		if( cx != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( y=0; y<= ny-1; y++){
		     	for( x=0; x<=(cx-1); x++)  tmp[ ii( (x+(nx-cx)),y,z) ] = in[ii(x,y,z)];
		     	for( x=cx;  x<=nx-1; x++)  tmp[ ii(      (x-cx),y,z) ] = in[ii(x,y,z)];
		    }
		else
		    for( i=0; i<=nnn-1; i++)
			tmp[i] = in[i];
		    
		/* Shifting y-coordinate*/
		if( cy != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( x=0; x<= nx-1; x++){
		    	for( y=0; y<=(cy-1); y++)  tmp2[ ii( x,(y+(ny-cy)),z) ] = tmp[ii(x,y,z)];
		    	for( y=cy;  y<=ny-1; y++)  tmp2[ ii( x,     (y-cy),z) ] = tmp[ii(x,y,z)];
		    }
		else
		    for( i=0; i<=nnn-1; i++)
			tmp2[i] = tmp[i];
		    
		/* Shifting z-coordinate*/
		if( cz != 0 )
		    for( x=0; x<= nx-1; x++)
		    for( y=0; y<= ny-1; y++){
		     	for( z=0; z<=(cz-1); z++)  in[ ii( x, y, (z+(nz-cz))) ] = tmp2[ii(x,y,z)];
		    	for( z=cz;  z<=nz-1; z++)  in[ ii( x, y,      (z-cz)) ] = tmp2[ii(x,y,z)];
		    }
		else
		    for( i=0; i<=nnn-1; i++)
			in[i] = tmp2[i];
		    
	free(tmp);
	free(tmp2);
}

void unshift_origin_inplace( double *in, ENV_PAR sys )
{
	int nx=sys.nx, ny=sys.ny, nz=sys.nz;		/*EXTERN*/
	int cx=sys.cx, cy=sys.cy, cz=sys.cz;
	int nnn = nx*ny*nz;

	int x, y, z, i;

	double *tmp  = (double *) malloc( nnn *sizeof( double )); 
	double *tmp2 = (double *) malloc( nnn *sizeof( double )); 

		/* unShifting z-coordinate*/
		if( cz != 0 )
		    for( x=0; x<= nx-1; x++)
		    for( y=0; y<= ny-1; y++){
			for( z=0; z<=(nz-1)-cz; z++)
		    	    tmp[ ii( x, y, (z+cz)) ] = in[ii(x,y,z)];
			for( z=nz-cz; z<=nz-1; z++)
			    tmp[ ii( x, y, (z-(nz-cz))) ] = in[ii(x,y,z)];
		    }
		else
		    for( i=0; i<=nnn-1; i++)
			tmp[i] = in[i];
		    
		/* unShifting y-coordinate*/
		if( cy != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( x=0; x<= nx-1; x++){
		    	for( y=0; y<=(ny-1)-cy; y++)
	    		    tmp2[ ii( x,(y+cy),z) ] = tmp[ii(x,y,z)];
	    	    	for( y=ny-cy; y<=ny-1; y++)
	    		    tmp2[ ii( x,(y-(ny-cy)),z) ] = tmp[ii(x,y,z)];
		    }
		else
		    for( i=0; i<=nnn-1; i++)
			tmp2[i] = tmp[i];
		    
		/* unShifting x-coordinate*/
		if( cx != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( y=0; y<= ny-1; y++){
		    	for( x=0; x<=(nx-1)-cx; x++)
	    		    in[ ii( (x+cx),y,z) ] = tmp2[ii(x,y,z)];
		    	for( x=nx-cx; x<=nx-1; x++)
		    	    in[ ii( (x-(nx-cx)),y,z) ] = tmp2[ii(x,y,z)];
		    }
		else
		    for( i=0; i<=nnn-1; i++)
			in[i] = tmp2[i];
		    
	free(tmp);
	free(tmp2);

}

void shift_origin_complex_inplace( fftw_complex *in, ENV_PAR sys)
{
	int nx=sys.nx, ny=sys.ny, nz=sys.nz;		/*EXTERN*/
	int cx=sys.cx, cy=sys.cy, cz=sys.cz;
	int nnn = nx*ny*nz;

	int x, y, z, i;

	fftw_complex *tmp  = (fftw_complex *) malloc( nnn *sizeof( fftw_complex )); 
	fftw_complex *tmp2 = (fftw_complex *) malloc( nnn *sizeof( fftw_complex )); 

		/* Shifting x-coordinate*/
		if( cx != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( y=0; y<= ny-1; y++){
		    	for( x=0; x<=(cx-1); x++){
	    		    tmp[ ii( (x+(nx-cx)),y,z) ][0] = in[ii(x,y,z)][0];
		    	    tmp[ ii( (x+(nx-cx)),y,z) ][1] = in[ii(x,y,z)][1];
			}
		        for( x=cx; x<=nx-1; x++){
		    	    tmp[ ii( (x-cx),y,z) ][0] = in[ii(x,y,z)][0];
	    		    tmp[ ii( (x-cx),y,z) ][1] = in[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			tmp[i][0] = in[i][0];
			tmp[i][1] = in[i][1];
		    }
		    
		/* Shifting y-coordinate*/
		if( cy != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( x=0; x<= nx-1; x++){
		    	for( y=0; y<=(cy-1); y++){
	        	    tmp2[ ii( x,(y+(ny-cy)),z) ][0] = tmp[ii(x,y,z)][0];
		            tmp2[ ii( x,(y+(ny-cy)),z) ][1] = tmp[ii(x,y,z)][1];
			}
		    	for( y=cy; y<=ny-1; y++){
		            tmp2[ ii( x,(y-cy),z) ][0] = tmp[ii(x,y,z)][0];
	        	    tmp2[ ii( x,(y-cy),z) ][1] = tmp[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			tmp2[i][0] = tmp[i][0];
			tmp2[i][1] = tmp[i][1];
		    }
		    
		/* Shifting z-coordinate*/
		if( cz != 0 )
		    for( x=0; x<= nx-1; x++)
		    for( y=0; y<= ny-1; y++){
		    	for( z=0; z<=(cz-1); z++){
	        	    in[ ii( x, y, (z+(nz-cz))) ][0] = tmp2[ii(x,y,z)][0];
		            in[ ii( x, y, (z+(nz-cz))) ][1] = tmp2[ii(x,y,z)][1];
			}
		    	for( z=cz; z<=nz-1; z++){
		            in[ ii( x, y, (z-cz)) ][0] = tmp2[ii(x,y,z)][0];
	        	    in[ ii( x, y, (z-cz)) ][1] = tmp2[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			in[i][0] = tmp2[i][0];
			in[i][1] = tmp2[i][1];
		    }
		    
		free(tmp);
		free(tmp2);
}

void unshift_origin_complex_inplace( fftw_complex *in, ENV_PAR sys)
{
	int nx=sys.nx, ny=sys.ny, nz=sys.nz;		/*EXTERN*/
	int cx=sys.cx, cy=sys.cy, cz=sys.cz;
	int nnn = nx*ny*nz;

	int x, y, z, i;

	fftw_complex *tmp  = (fftw_complex *) malloc( nnn *sizeof( fftw_complex )); 
	fftw_complex *tmp2 = (fftw_complex *) malloc( nnn *sizeof( fftw_complex )); 

		/* unShifting z-coordinate*/
		if( cz != 0 )
		    for( x=0; x<= nx-1; x++)
		    for( y=0; y<= ny-1; y++){
		    	for( z=0; z<=(nz-1)-cz; z++){	
			    tmp[ ii( x, y, (z+cz)) ][0] = in[ii(x,y,z)][0];
			    tmp[ ii( x, y, (z+cz)) ][1] = in[ii(x,y,z)][1]; 
		    	}
		    	for( z=nz-cz; z<=nz-1; z++){
			    tmp[ ii( x, y, (z-(nz-cz))) ][0] = in[ii(x,y,z)][0];
			    tmp[ ii( x, y, (z-(nz-cz))) ][1] = in[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			tmp[i][0] = in[i][0];
			tmp[i][1] = in[i][1];
		    }
		    
		/* unShifting y-coordinate*/
		if( cy != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( x=0; x<= nx-1; x++){
		    	for( y=0; y<=(ny-1)-cy; y++){
		    	    tmp2[ ii( x,(y+cy),z) ][0] = tmp[ii(x,y,z)][0];
		    	    tmp2[ ii( x,(y+cy),z) ][1] = tmp[ii(x,y,z)][1];
			}
		    	for( y=ny-cy; y<=ny-1; y++){
		    	    tmp2[ ii( x,(y-(ny-cy)),z) ][0] = tmp[ii(x,y,z)][0];
		    	    tmp2[ ii( x,(y-(ny-cy)),z) ][1] = tmp[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			tmp2[i][0] = tmp[i][0];
			tmp2[i][1] = tmp[i][1];
		    }
		    
		/* unShifting x-coordinate*/
		if( cz != 0 )
		    for( z=0; z<= nz-1; z++)
		    for( y=0; y<= ny-1; y++){
		    	for( x=0; x<=(nx-1)-cx; x++){
		    	    in[ ii( (x+cx),y,z) ][0] = tmp2[ii(x,y,z)][0];
		    	    in[ ii( (x+cx),y,z) ][1] = tmp2[ii(x,y,z)][1];
			}
		    	for( x=nx-cx; x<=nx-1; x++){
		    	    in[ ii( (x-(nx-cx)),y,z) ][0] = tmp2[ii(x,y,z)][0];
		    	    in[ ii( (x-(nx-cx)),y,z) ][1] = tmp2[ii(x,y,z)][1];
			}
		    }
		else
		    for( i=0; i<=nnn-1; i++){
			in[i][0] = tmp2[i][0];
			in[i][1] = tmp2[i][1];
		    }
		    
	free(tmp);
	free(tmp2);

}



