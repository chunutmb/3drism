
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "nrutil.h"
#include "jh_linalg.h"


void  AddAsgn_VV( double * A, double * B, int n)
{
	int i;
	for( i=0; i<=n-1; i++)
	    A[i] += B[i];

}


void Asgn_VV( double * A, double *B, int n)
{
	int i;

	for( i=0; i<=n-1; i++)
		A[i] = B[i];
}

double l2Norm_V( double * A, int dim )
{
	int i;
	double tmp=0;

	for( i=0; i<= dim-1; i++)
	    tmp += A[i] * A[i];

	tmp = pow( tmp, 0.5 );
	return tmp;

}


double Vec_Vec(  double * Aa,  double * Bb, int nn)
/*Multiply vec A x vec B and return value*/
{
	int i;
	double tmps=0;

	for( i=0; i<=nn-1; i++)	tmps += Aa[i]*Bb[i];
	return tmps;
}

void Set_VV( double * A, double *B, int n)
/*set vec A to vec B*/
{
	int i;

	for( i=0; i<=n-1; i++)
		A[i] = B[i];
}

void Sub_VV( double  * A, double * B, double * C, int n)
/*sub vec C from vec B and return new vec A*/
{
	int i;

	for( i=0; i<=n-1; i++)
		A[i] = B[i] - C[i];

}


void Matrix_Vec( double * Avj, double ** Aav, double * Bbv, int nr, int nc)
{
	int i, k;
	double  tmps;

	for( i=0; i <=nr-1; i++){ tmps = 0;
				  for( k=0; k<= nc-1; k++)  tmps += Aav[i][k]*Bbv[k];
				  Avj[i] = tmps;
				};
}

double * Mul_MV( double ** A, double *B, int nr, int nc)
{
	int i, k;
	double tmp;

	double *AB = (double *) malloc( nc * sizeof(double));

        for( i=0; i <=nr-1; i++){ tmp = 0;
                                  for( k=0; k<= nc-1; k++)  tmp += A[i][k]*B[k];
                                  AB[i] = tmp;
                                  };
	return AB ;
}

double Mul_VV( double * A, double * B, int dim)
{
	int i;
	double tmp=0;
	for( i=0; i<= dim-1; i++) tmp += A[i]*B[i];
	return tmp;
}

double * Mul_SV( double cnst , double * A, int dim )
{
	int i;
	double *B = (double *) malloc( dim *sizeof(double));
	for( i=0; i<=dim-1; i++)    B[i] = cnst * A[i];
	return B;
}

double * MulAsgn_VS( double * A, double cnst, int dim )
{
	int i;
	for( i=0; i<=dim-1; i++) A[i] *= cnst;
	return A;
}

double * SubAsgn_VV( double * A, double * B , int dim)
{
	int i;

	for( i=0; i<=dim-1; i++)
		A[i] = A[i] - B[i];

	return A;
}

void Matrix_Matrix( double ** Xl, double ** Al, double ** Bl, int il)
{

	int i, j, k;

                for( i=0; i<=il-1; i++)
                    for( j=0; j<=il-1; j++){
                        Xl[i][j] = 0;
                        for( k=0; k<=il-1; k++)
                            Xl[i][j] += Al[i][k]*Bl[k][j];   /*Multiply Matrices*/
						};
}


void JH_inp1_mat_mat( double ** A, double ** B, int n)
{

	int i, j, k;
	/*double **X = (double **) Matrix_Malloc( n, n);*/

	double **X = (double **) malloc( n *sizeof(double));
	for( i=0; i<=n-1; i++)
		*(X+i) = (double *) malloc( n *sizeof(double));

                for( i=0; i<=n-1; i++)
                    for( j=0; j<=n-1; j++){
                        X[i][j] = 0;
                        for( k=0; k<=n-1; k++)
                            X[i][j] += A[i][k]*B[k][j];   /*Multiply Matrices*/
						};
		for( i=0; i<=n-1; i++)
		    for( j=0; j<=n-1; j++)
			A[i][j] = X[i][j];

}


void JH_inp2_mat_mat( double ** A, double ** B, int n)
{

	int i, j, k;
	/*double **X = (double **) Matrix_Malloc( n, n);*/

	double **X = (double **) malloc( n *sizeof(double));
	for( i=0; i<=n-1; i++)
		*(X+i) = (double *) malloc( n *sizeof(double));

	        for( i=0; i<=n-1; i++)
                    for( j=0; j<=n-1; j++){
                        X[i][j] = 0;
                        for( k=0; k<=n-1; k++)
                            X[i][j] += A[i][k]*B[k][j];   /*Multiply Matrices*/
						};
		for( i=0; i<=n-1; i++)
		    for( j=0; j<=n-1; j++)
			B[i][j] = X[i][j];

}



void JMat_Inv_2x2( double **I, double **Mat)
{
	double tmp, M[2][2];

	M[0][0] = Mat[0][0];   M[0][1] = Mat[0][1];
	M[1][0] = Mat[1][0];   M[1][1] = Mat[1][1];

        I[0][0]=1; I[0][1]=0; I[1][0]=0; I[1][1]=1;

      	tmp = (M[1][0]/M[0][0]);
	M[1][0] -= (tmp * M[0][0]);	M[1][1] -= (tmp * M[0][1]);
	I[1][0] -= (tmp * I[0][0]);	I[1][1] -= (tmp * I[0][1]);

      	tmp = (M[0][1]/M[1][1]);
	M[0][0] -= (tmp * M[1][0]);	M[0][1] -= (tmp * M[1][1]);
	I[0][0] -= (tmp * I[1][0]);	I[0][1] -= (tmp * I[1][1]);

	tmp = M[0][0];
	M[0][0] = (M[0][0]/tmp);        M[0][1] = (M[0][1]/tmp);
        I[0][0] = (I[0][0]/tmp);        I[0][1] = (I[0][1]/tmp);

	tmp = M[1][1];
        M[1][0] = (M[1][0]/tmp);        M[1][1] = (M[1][1]/tmp);
        I[1][0] = (I[1][0]/tmp);        I[1][1] = (I[1][1]/tmp);

}

/*solve Ax=b -> Js=F with size vn */

void Newton_Solver(  double * x1,  double ** A,  double * b, int vn)
{
	int i, j;
	double ger;

	int *p = ivector( 1, vn);

	double  **LU = dmatrix( 1, vn, 1, vn);

		for( i=0; i<=vn-1; i++)
		    for( j=0; j<=vn-1; j++)
			LU[i+1][j+1] = A[i][j];

	double *X = dvector( 1, vn);
		for( i=0; i<=vn-1; i++)
		    X[i+1] = b[i];

 	ludcmp( LU, vn, p, &ger);

	lubksb( LU, vn, p, X );

	for( i=0; i<=vn-1; i++)
	    x1[i] = X[i+1];

}

void Matrix_Inverse( double **invMat, double **Mat, int nn)
{
        int i, j;
        double ger;

        int *p= ivector( 1, nn);


        double **LU = dmatrix( 1, nn, 1, nn);
                for( i=0; i<=nn-1; i++)
                    for( j=0; j<=nn-1; j++)
                        LU[i+1][j+1] = Mat[i][j];

        double *col = dvector( 1, nn);


        ludcmp( LU, nn, p, &ger);

	for( j=1; j<=nn; j++) {
	    for( i=1; i<=nn; i++) col[i] = 0.0;
 	    col[j] = 1.0;
	    lubksb( LU, nn, p, col );
	    for( i=1; i<=nn; i++) invMat[i-1][j-1] = (double) col[i];
			     };

	free_dvector( col, 1, nn);
	free_ivector( p, 1, nn);
	free_dmatrix( LU, 1, nn, 1, nn);

}

#define TINY 1.0e-20

void ludcmp(double **a, int n, int *index, double *d)
{
        int i,imax,j,k;
        double big,dum,sum,temp;
        double *vv;

        vv=dvector(1,n);
        *d=1.0;
        for (i=1;i<=n;i++) {
                big=0.0;
                for (j=1;j<=n;j++)
                        if ((temp=fabs(a[i][j])) > big) big=temp;
                if (big == 0.0) {nrerror("Singular matrix in routine ludcmp");
				 printf( "Singular Matrix in routine ludcmp");};
                vv[i]=1.0/big;
        }
        for (j=1;j<=n;j++) {
                for (i=1;i<j;i++) {
                        sum=a[i][j];
                        for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                }
                big=0.0;
                imax = j;
                for (i=j;i<=n;i++) {
                        sum=a[i][j];
                        for (k=1;k<j;k++)
                                sum -= a[i][k]*a[k][j];
                        a[i][j]=sum;
                        if ( (dum=vv[i]*fabs(sum)) >= big) {
                                big=dum;
                                imax=i;
                        }
                }
                if (j != imax) {
                        for (k=1;k<=n;k++) {
                                dum=a[imax][k];
                                a[imax][k]=a[j][k];
                                a[j][k]=dum;
                        }
                        *d = -(*d);
                        vv[imax]=vv[j];
                }
                index[j]=imax;
                if (a[j][j] == 0.0) a[j][j]=TINY;
                if (j != n) {
                        dum=1.0/(a[j][j]);
                        for (i=j+1;i<=n;i++) a[i][j] *= dum;
                }
        }
        free_dvector(vv,1,n);
}

#undef TINY

void lubksb( double **a, int n, int *index, double *b)
{
        int i,ii=0,ip,j;
        double sum;

        for (i=1;i<=n;i++) {
                ip=index[i];
                sum=b[ip];
                b[ip]=b[i];
                if (ii)
                        for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
                else if (sum) ii=i;
                b[i]=sum;
        }
        for (i=n;i>=1;i--) {
                sum=b[i];
                for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
                b[i]=sum/a[i][i];
        }
}


void Matrix_Inverse2( double **invMat, double **Mat, int nn)
{
        int i, j;

        double **LU = dmatrix( 1, nn, 1, nn);
                for( i=0; i<=nn-1; i++)
                    for( j=0; j<=nn-1; j++)
                        LU[i+1][j+1] = Mat[i][j];

        double **col = dmatrix( 1, nn, 1, 1);
	    for( i=1; i<=nn; i++) col[i][1] = 0.0;

	gaussj( LU, nn, col, 1);

	for( j=1; j<=nn; j++)
	    for( i=1; i<=nn; i++) invMat[i-1][j-1] = (double) LU[i][j];

	free_dmatrix( col, 1, nn, 1, 1);
	free_dmatrix( LU, 1, nn, 1, nn);

}

#define NRANSI
#define SWAP(a,b) {temp=(a); (a)=(b); (b)=temp;}
void gaussj(double **a, int n, double **b, int m)
{
        int *indxc,*indxr,*ipiv;
        int i,icol,irow,j,k,l,ll;
        double big,dum,pivinv,temp;

        indxc=ivector(1,n);
        indxr=ivector(1,n);
        ipiv=ivector(1,n);
        for (j=1;j<=n;j++) ipiv[j]=0;
        for (i=1;i<=n;i++) {
                big=0.0;
                irow = icol = 0;
                for (j=1;j<=n;j++)
                        if (ipiv[j] != 1)
                                for (k=1;k<=n;k++) {
                                        if (ipiv[k] == 0) {
                                                if (fabs(a[j][k]) >= big) {
                                                        big=fabs(a[j][k]);
                                                        irow=j;
                                                        icol=k;
                                                }
                                        } else if (ipiv[k] > 1) nrerror("gaussj: Singular Matrix-1");
                                }
                ++(ipiv[icol]);
                if (irow != icol) {
                        for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
                        for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
                }
                indxr[i]=irow;
                indxc[i]=icol;
                if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix-2");
                pivinv=1.0/a[icol][icol];
                a[icol][icol]=1.0;
                for (l=1;l<=n;l++) a[icol][l] *= pivinv;
                for (l=1;l<=m;l++) b[icol][l] *= pivinv;
                for (ll=1;ll<=n;ll++)
                        if (ll != icol) {
                                dum=a[ll][icol];
                                a[ll][icol]=0.0;
                                for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
                                for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
                        }
        }
        for (l=n;l>=1;l--) {
                if (indxr[l] != indxc[l])
                        for (k=1;k<=n;k++)
                                SWAP(a[k][indxr[l]],a[k][indxc[l]])
        }
        free_ivector(ipiv,1,n);
        free_ivector(indxr,1,n);
        free_ivector(indxc,1,n);
}
#undef SWAP
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software +11"?. */

