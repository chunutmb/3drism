/****************************************************************************************/
/*		TITLE: 3D-Rism code for a solute in H2O					*/
/****************************************************************************************/
/*											*/
/*		Author: Jesse Howard							*/
/*		Update: Oct 1, 2008							*/
/*		Version:  1.5								*/
/*											*/
/*											*/
/****************************************************************************************/


/********************NOTES*************************
*
* fully functioning mdiis routine is running
* dk factors = R / (lx/2)
*
* 1. This version reads in a dis2 or kdis2 type file
* 2. This version reads in a env type file
* 3. This version reads in a par type file
* 4. Ewald sums are implemented
* /
**************************************************/

/*********************OPTIONS****************************************
* define (OMP) -> include the OMP directives to utilize multiprocessor on a shared memory system
* define (
*
*
* ******************************************************************/

#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "jh_get.h"
#include "jh_struct.h"
#include "jh_grid.h"
#include "jh_util.h"
#include "jh_linalg.h"
#include "jh_print.h"

/*#define FFTW_THREADS*/
#ifdef FFTW_THREADS
        #include <pthread.h>
        #define NUM_THREADS 16
#endif

#ifdef MPI
        #include <mpi.h>
#endif

#define NNN NX * NY * NZ
#define ii(x,y,z) (NZ * NY * (x) + NZ * (y) + (z))

/*env - values*/
int NX, NY, NZ;
int CX, CY, CZ;
double LX, LY, LZ;
int MAX_ITER, N_DUMP, INIT_PIC_ITER, DIIS_SIZE;
double T_ERR, A_ERF, CHRG_PCT, PIC_MP, DIIS_MP, TEMP_FACTOR;
char *CONTINUE, *SOLVER, *CLOSURE, *EWALD_SUMS, *FILE_TYPE;
char *CONFIG_TYPE, *BRIDGE_FUNC1, *BRIDGE_FUNC0, *RBC_FUNC;
int RT_CHANGES;
double RT_ERR, *RT_TEMP_FACTOR;


/*dis2 - solvent properties*/
int NSITES, NRSITES, TYPE, DIS_NUM;
double TEMP, *PND, *REDUN;
double *EP12, *EP6, *SIG, *CHARGE, *BOND_MAT;
char **NAMES, **DIS_NAMES;


/*par properties*/
ENV_PAR SYS;
U_PAR2  *U, *U2;
U_PAR *U1;
int NU_SITES, PAR_TYPE;

double *EWALD_BGC;
int STAT = 1, CNT, PIC_CNT = 0, RT_CNT = 0, MDIIS_CNT = 0;
int UPB1 = 1, B1R2 = 0, UPB1_STAT = 0;

double ***HK3;          /* [row][col][idx] */
double *WK_OH;          /* [idx] */
double *WK_HH;
double **CR_S;          /* [site][idx] */
double **CR2_S;         /* [site][idx] */
double **TR_S;
double **EXP_BR;
double *HS, *HS2, **HSW;
double **B1R;

double **UR_S, **UR_L;  /* [site][idx] */
fftw_complex **UK_L;
double IERR = 9e9, FERR;

int UR_S_STAT = 0, UR_L_STAT = 0, UK_L_STAT = 0, HK3_STAT = 0;
int CR_S_STAT = 0, CR2_S_STAT = 0, TR_S_STAT = 0, HS_STAT = 0;
int B1R_STAT = 0, EXP_BR_STAT = 0, WK_STAT = 0;

void set_env(char []);
void set_dis(char []);
void set_par(char []);
void set_arrays(void);
void set_sys(void);

void check_env(void);
void check_dis(void);
void check_par(U_PAR2 *, int);

/*________________Printing Routines____________________*/
void print_1d(char [], double *, int, double);
void print_3d(char [], double *);
void print_jh3d_box(char [], double *, int, int, int);
void print_jh3d_box3(char [], double *, int, int, int);

double * get_3d(char [], double, double, ENV_PAR);
fftw_complex * get_complex_3d(char [], double, double, ENV_PAR);

int set_potential_fields(void);
int calc_intramolecular_functions(void);
int calc_hk_vv_correlations(char []);
int calc_rbc(void);
void calc_bridge_func1_fbond(double **);
void calc_bridge_func1(int);
void change_RT(void);

/*_____________Numerical Routines_______________*/
void full_picard_iter(int);   /*max_iter, r_err, t_err*/
void mdiis_iter(int);

/*______________FFTW routines________________*/
void fftw_3d(fftw_complex *, fftw_complex *);           /*in_r, out_k*/
void invfftw_3d(fftw_complex *, fftw_complex *);        /*in_k, out_r*/

/*_____________MPI variables_________________*/
int MY_RANK = 0;
int NP = 1;



/*Global variables that don't change, solvent variables are read in */
/*when executing program 1st arg solute file, 2nd arg solvent hr_vv*/
/* keep flow, printing, calc global variables */
/************************************************************************************************/
/*					--------						*/
/*					| MAIN |						*/
/*				        --------						*/
/*	1. 1st argument is the solvent file (.dis2) or (.kdis2)					*/
/*	2. 2nd argument is the env file (.env)                                                  */
/*	3. 3rd argument is the solute file (.par)						*/
/*												*/
/************************************************************************************************/

int main(int argc, char *argv[])
{
  char s1[100];
  int j, i;

  if (argc < 3) {
    fprintf(stderr, "need three parameters\n");
    exit(1);
  }

  set_env(*(argv + 2));
  check_env();
  set_dis(*(argv + 1));
  check_dis();
  set_sys();
  printf("\n%d:Reading in solute parameters from: \"%s\"... ", MY_RANK, *(argv + 3)); fflush(stdout);
  set_par(*(argv + 3));
  check_par(U, NU_SITES);            printf("%d:...done!!!\n", MY_RANK); fflush(stdout);

  /*U(r), U(k), W(k), H(k)_vv, exp(b(r))  **********************************************************/

  printf("\nCalculating Hk_vv, reading from \"%s\" ...", *(argv + 1)); fflush(stdout);
  calc_hk_vv_correlations(*(argv + 1));    printf("\n.......H(k)_vv...done!!!\n\n"); fflush(stdout);
  printf("Calculating Intramolecular function W(k)..."); fflush(stdout);
  if (TYPE == 1 || TYPE == 2)
    calc_intramolecular_functions(); printf("...done!!!\n\n"); fflush(stdout);

  printf("Calculating All Potential Energies U(r)..."); fflush(stdout);
  set_potential_fields();                 printf("...Done!!!\n\n"); fflush(stdout);

  if (TYPE == 1 || TYPE == 2)
    if (strncmp("yes", RBC_FUNC, 3) == 0) {
      printf("Calculating RBC function exp( b(r) )..."); fflush(stdout);
      calc_rbc();                 printf("...done!!!\n\n"); fflush(stdout);
    }
                                                                                                   #ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
                                                                                                   #endif
  /*Global Arrays - cr, cr2, tr************************************************/
  printf("%d:Allocating and setting arrays...", MY_RANK); fflush(stdout);
  set_arrays();                           printf("...done!!!\n\n"); fflush(stdout);
  printf("Initial picard iterations(%d)...\n", INIT_PIC_ITER); fflush(stdout);
  if (strncmp("no", CONTINUE, 2) == 0)
    full_picard_iter(INIT_PIC_ITER);
  printf("...done!!!\n"); fflush(stdout);

  /****************************************************************************************/
  /*                             Numerical routine                                        */
  /****************************************************************************************/

  if (strncmp("picard", SOLVER, 6) == 0) {
    printf("%d:starting picard iteration\n", MY_RANK); fflush(stdout);
    full_picard_iter(MAX_ITER);
  } else if (strncmp("mdiis", SOLVER, 5) == 0) {
    printf("%d:starting mdiis iteration\n", MY_RANK); fflush(stdout);
    mdiis_iter(MAX_ITER);
  } else if (strncmp("newton", SOLVER, 6) == 0) {
    printf("This part of code has been unkept: see ver.1.5");
    exit(1);
  }


  /****************************************************************************************/
  /******************************END Numerical Routine*************************************/
  /****************************************************************************************/


  double wr;
  double **gr = ( double **) malloc(NRSITES * sizeof(double));
  for (j = 0; j <= NRSITES - 1; j++)
    *(gr + j) = (double *) malloc(NNN * sizeof(double));


  if (STAT > 0 && FERR < 10.0) {
    printf("\n\n\nDone with iterations: Printing data..."); fflush(stdout);

    if (UR_S_STAT != 0) {
      for (j = 0; j <= NRSITES - 1; j++)
        unshift_origin_inplace(*(UR_S + j), SYS);
      UR_S_STAT = 0;
    }
    if (TR_S_STAT != 0) {
      for (j = 0; j <= NRSITES - 1; j++)
        unshift_origin_inplace(*(TR_S + j), SYS);
      TR_S_STAT = 0;
    }
    if (CR_S_STAT != 0) {
      for (j = 0; j <= NRSITES - 1; j++)
        unshift_origin_inplace(*(CR_S + j), SYS);
      CR_S_STAT = 0;
    }
    if (B1R2 == 1)
      if (B1R_STAT != 0) {
        for (j = 0; j <= NRSITES - 1; j++)
          unshift_origin_inplace(*(B1R + j), SYS);
        B1R_STAT = 0;
      }


    /*HNC*/ if ((strncmp("hnc", CLOSURE, 3) == 0) && (strlen(CLOSURE) == 3)) {
      if (B1R2 == 1) {
        for (j = 0; j <= NRSITES - 1; j++)
          for (i = 0; i <= NNN - 1; i++)
            gr[j][i] = exp(-UR_S[j][i] + TR_S[j][i] + B1R[j][i]);
      } else {
        for (j = 0; j <= NRSITES - 1; j++)
          for (i = 0; i <= NNN - 1; i++)
            gr[j][i] = exp(-UR_S[j][i] + TR_S[j][i]);
      }
    } else
    /*PY*/ if ((strncmp("py", CLOSURE, 2) == 0) && (strlen(CLOSURE) == 2)) {
      for (j = 0; j <= NRSITES - 1; j++)
        for (i = 0; i <= NNN - 1; i++)
          gr[j][i] = exp(-UR_S[j][i]) * (1 + TR_S[j][i]);
    } else
    /*KH*/ if ((strncmp("kh", CLOSURE, 2) == 0) && (strlen(CLOSURE) == 2)) {
      for (j = 0; j <= NRSITES - 1; j++)
        for (i = 0; i <= NNN - 1; i++) {
          wr = -UR_S[j][i] + TR_S[j][i];
          if (wr > 0.00)
            gr[j][i] = wr + 1.00;
          else
            gr[j][i] = exp(-UR_S[j][i] + TR_S[j][i]);
        }
    } else {
      printf("\n\nNo closure specified\n"); fflush(stdout);
      exit(1);
    }



    /*ADD RBC function*/
    if (strncmp("yes", RBC_FUNC, 3) == 0) {
      for (j = 0; j <= NRSITES - 1; j++)
        unshift_origin_inplace(*(EXP_BR + j), SYS);                                             /*(in, out), in*/
      for (j = 0; j <= NRSITES - 1; j++)
        for (i = 0; i <= NNN - 1; i++)
          gr[j][i] = gr[j][i] * EXP_BR[j][i];
    }


    /*CONFIG*/
    if ((strncmp("wall1", CONFIG_TYPE, 5) == 0) && (strlen(CONFIG_TYPE) == 5)) {
      if (HS_STAT != 0)
        unshift_origin_inplace(HS, SYS);                                /*in, out, in*/
      for (j = 0; j <= NRSITES - 1; j++)
        for (i = 0; i <= NNN - 1; i++) {
          if (HS[i] == 0.0)
            gr[j][i] = 0.00;
        }
    } else if ((strncmp("wall2", CONFIG_TYPE, 5) == 0) && (strlen(CONFIG_TYPE) == 5)) {
      if (HS_STAT != 0)
        for (j = 0; j <= NRSITES - 1; j++) {
          unshift_origin_inplace(*(HSW + j), SYS);                             /*in, out, in*/
          for (i = 0; i <= NNN - 1; i++) {
            if (HSW[j][i] == 0.0)
              gr[j][i] = 0.00;
          }
        }
    }


    /*___print gr___*/
    for (j = 0; j <= NRSITES - 1; j++) {
      sprintf(s1, "gr_%s", NAMES[j]);
      print_3d(s1, *(gr + j));
    }

    if (strncmp("wall", CONFIG_TYPE, 4) == 0)
      for (j = 0; j <= NRSITES - 1; j++) {
        sprintf(s1, "box_gr_%s.jh3d", NAMES[j]);
        if (strncmp("wall1", CONFIG_TYPE, 5) == 0)
          print_jh3d_box3(s1, *(gr + j), NX, NY, NZ);
        else
          print_jh3d_box(s1, *(gr + j), NX, NY, NZ);
      }

    /*___print cr_s___*/
    for (j = 0; j <= NRSITES - 1; j++) {
      sprintf(s1, "cr_%s_s", NAMES[j]);
      print_3d(s1, *(CR_S + j));
    }

    if (B1R2 == 1) {
      for (j = 0; j <= NRSITES - 1; j++) {
        sprintf(s1, "b1r_%s", NAMES[j]);
        print_3d(s1, *(B1R + j));
      }
    }
    printf("...done!"); fflush(stdout);
  }


  if (strncmp("picard", SOLVER, 6) == 0) {
    if (PIC_CNT >= MAX_ITER)
      printf("\n\nc[r] did NOT converge\n\n");
    else
      printf("\n\nc[r] converged in %d picard iterations\n\n", PIC_CNT);
  } else if (strncmp("mdiis", SOLVER, 5) == 0) {
    if (MDIIS_CNT >= MAX_ITER)
      printf("\n\nc[r] did NOT converge\n\n");
    else
      printf("\n\nc[r] converged in %d mdiis iterations\n\n", MDIIS_CNT);
  }

  printf("Program 3drism is done\n");

  return 0;
}



/**************************************END OF MAIN***************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
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
/*											*/
/*											*/
/*				SUBROUTINES						*/
/*											*/
/*											*/
/*											*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/



/****************************************************************************************/
/*					I/O						*/
/*				GET Parameter and solvent				*/
/*											*/
/****************************************************************************************/

void set_env(char infile[])
{
  NX = (int) get_dval(infile, "NX");
  NY = (int) get_dval(infile, "NY");
  NZ = (int) get_dval(infile, "NZ");
  LX = (double) get_dval(infile, "LX");
  LY = (double) get_dval(infile, "LY");
  LZ = (double) get_dval(infile, "LZ");

  T_ERR = (double) get_dval(infile, "T_ERR");
  CLOSURE = (char *)get_sval(infile, "CLOSURE");
  SOLVER = (char *)get_sval(infile, "SOLVER");
  PIC_MP = (double) get_dval(infile, "PIC_MP");
  if (strncmp("mdiis", SOLVER, 5) == 0) {
    DIIS_SIZE = (int) get_dval(infile,  "DIIS_SIZE");
    DIIS_MP = (double) get_dval(infile, "DIIS_MP");
  }


  if (get_tag(infile, "CX") == 1)
    CX = (int) get_dval(infile, "CX");
  else CX = (int) NX / 2;

  if (get_tag(infile, "CY") == 1)
    CY = (int) get_dval(infile, "CY");
  else CY = (int) NY / 2;

  if (get_tag(infile, "CZ") == 1)
    CZ = (int) get_dval(infile, "CZ");
  else CZ = (int) NZ / 2;

  if (get_tag(infile, "A_ERF") == 1)
    A_ERF = (double) get_dval(infile, "A_ERF");
  else A_ERF = 1.08;

  if (get_tag(infile, "BRIDGE_FUNC0") == 1)
    BRIDGE_FUNC0 = (char *) get_sval(infile, "BRIDGE_FUNC0");
  else BRIDGE_FUNC0 = "no";

  if (get_tag(infile, "BRIDGE_FUNC1") == 1)
    BRIDGE_FUNC1 = (char *) get_sval(infile, "BRIDGE_FUNC1");
  else BRIDGE_FUNC1 = "no";

  if (get_tag(infile, "CHRG_PCT") == 1)
    CHRG_PCT = (double) get_dval(infile, "CHRG_PCT");
  else CHRG_PCT = 1.00;

  if (get_tag(infile, "CONFIG_TYPE") == 1)
    CONFIG_TYPE = (char *) get_sval(infile, "CONFIG_TYPE");
  else CONFIG_TYPE = "none";

  if (get_tag(infile, "CONTINUE") == 1)
    CONTINUE = (char *)get_sval(infile, "CONTINUE");
  else CONTINUE = "no";

  if (get_tag(infile, "EWALD_SUMS") == 1)
    EWALD_SUMS = (char *) get_sval(infile, "EWALD_SUMS");
  else EWALD_SUMS = "no";

  if (get_tag(infile, "FILE_TYPE") == 1)
    FILE_TYPE = (char *) get_sval(infile, "FILE_TYPE");
  else FILE_TYPE = "jh3d";

  if (get_tag(infile, "INIT_PIC_ITER") == 1)
    INIT_PIC_ITER = (int) get_dval(infile, "INIT_PIC_ITER");
  else INIT_PIC_ITER = 1.0;

  if (get_tag(infile, "MAX_ITER") == 1)
    MAX_ITER = (int) get_dval(infile, "MAX_ITER");
  else MAX_ITER = (int) INT_MAX;

  if (get_tag(infile, "N_DUMP") == 1)
    N_DUMP = (int) get_dval(infile, "N_DUMP");
  else N_DUMP = (int) INT_MAX;

  if (get_tag(infile, "TEMP_FACTOR") == 1)
    TEMP_FACTOR = (double) get_dval(infile, "TEMP_FACTOR");
  else TEMP_FACTOR = 1.00;

  if (get_tag(infile, "RBC_FUNC") == 1)
    RBC_FUNC = (char *) get_sval(infile, "RBC_FUNC");
  else RBC_FUNC = "no";

  if (get_tag(infile, "RT_CHANGES") == 1)
    RT_CHANGES = (int) get_dval(infile, "RT_CHANGES");
  else RT_CHANGES = 0;

  if (RT_CHANGES > 0) {
    if (get_tag(infile, "RT_TEMP_FACTOR") == 1)
      RT_TEMP_FACTOR = (double *) get_array_dval(infile, "RT_TEMP_FACTOR", RT_CHANGES);
    else RT_TEMP_FACTOR = NULL;

    if (get_tag(infile, "RT_ERR") == 1)
      RT_ERR = (double) get_dval(infile, "RT_ERR");
    else RT_ERR = T_ERR;
  }
}



void set_dis(char infile[])
{
  TEMP = (double) get_dval(infile, "TEMP");
  TYPE = (int)  get_dval(infile, "TYPE");
  NSITES = (int)  get_dval(infile, "NSITES");
  NRSITES = (int)  get_dval(infile, "NRSITES");
  REDUN = (double *) get_array_dval(infile, "REDUN", NRSITES);
  PND = (double *) get_array_dval(infile, "PND", NRSITES);
  EP12 = (double *) get_array_dval(infile, "EP12", NRSITES);
  EP6 = (double *) get_array_dval(infile, "EP6", NRSITES);
  SIG = (double *) get_array_dval(infile, "SIG", NRSITES);
  CHARGE = (double *) get_array_dval(infile, "CHARGE", NRSITES);
  BOND_MAT = (double *) get_array_dval(infile, "BOND_MAT", NRSITES);
  NAMES = (char **) get_array_sval(infile, "NAMES", NRSITES);
  DIS_NUM = (int ) get_dval(infile, "DIS_NUM");
  DIS_NAMES = (char **) get_array_sval(infile, "DIS_NAMES", DIS_NUM);
}



void set_par(char infile[])
{
  int *nu_s = (int *) malloc(sizeof(int));
  int i, len;
  char vs;

  len = strlen(infile);
  vs = infile[len - 1];

  if (vs == '2') {
    printf("reading from a par2 file...");
    PAR_TYPE = 2;
    U2 = (U_PAR2 *) get_par2(nu_s, infile);
    printf("done \n"); fflush(stdout);
  } else {
    printf("reading from a par file...");
    PAR_TYPE = 1;
    U1 = (U_PAR *) get_par(nu_s, infile);
    printf("done \n"); fflush(stdout);
  }

  NU_SITES = *nu_s;

  if (PAR_TYPE == 1) {
    U = (U_PAR2 *) malloc((NU_SITES + 1) * sizeof(U_PAR2));

    for (i = 1; i <= NU_SITES; i++) {
      U[i].num = U1[i].num;
      sprintf(U[i].element, "%s", U1[i].element);
      U[i].mol = 0;
      U[i].ep12 = U1[i].ep;
      U[i].ep6 = U1[i].ep;
      U[i].sig = U1[i].sig;
      U[i].charge = U1[i].charge;
      U[i].x = U1[i].x;
      U[i].y = U1[i].y;
      U[i].z = U1[i].z;
    }
  } else if (PAR_TYPE == 2)
    U = (U_PAR2 *) U2;

  /*U = (U_PAR *) get_par( nu_s, infile );*/
}



int calc_upb1(void);

void set_sys(void)
{
  int i, j, x, y, z;
  int hz;
  double dz, hsig;

  SYS.nx = NX;    SYS.ny = NY;    SYS.nz = NZ;
  SYS.cx = CX;    SYS.cy = CY;    SYS.cz = CZ;
  SYS.lx = LX;    SYS.ly = LY;    SYS.lz = LZ;

  if (strncmp("wall2", CONFIG_TYPE, 5) == 0) {
    HSW = (double **) malloc(NRSITES * sizeof(double));
    for (j = 0; j <= NRSITES - 1; j++)
      *(HSW + j) = (double *) malloc(NNN * sizeof(double));
    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        HSW[j][i] = 1.00;

    dz = (double)LZ / (NZ - 1);
    for (j = 0; j <= NRSITES - 1; j++) {
      hsig = 0.5 * SIG[j];
      hz = (hsig / dz);
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = (NZ / 6) + 1; z <= CZ + hz; z++)
            HSW[j][ii(x,y,z)] = 0.00;
      shift_origin_inplace(*(HSW + j), SYS);                            /*in, out, in*/
    }
    HS_STAT = 1;
  } else if (strncmp("wall1", CONFIG_TYPE, 5) == 0) {
    HS = (double *) malloc(NNN * sizeof(double));
    for (i = 0; i <= NNN - 1; i++) HS[i] = 1.00;

    for (x = 0; x <= NX - 1; x++)
      for (y = 0; y <= NY - 1; y++) {
        for (z = (NZ / 6) + 1; z <= CZ; z++)
          HS[ii(x,y,z)] = 0.00;
      }
    shift_origin_inplace(HS, SYS);                           /*in, out, in*/
    HS_STAT = 1;
  } else if (strncmp("wall", CONFIG_TYPE, 4) == 0) {
    HS = (double *) malloc(NNN * sizeof(double));
    for (i = 0; i <= NNN - 1; i++) HS[i] = 0.00;

    for (x = 0; x <= NX - 1; x++)
      for (y = 0; y <= NY - 1; y++)
        for (z = CZ; z <= NZ - 1; z++)
          HS[ii(x,y,z)] = 1.00;

    shift_origin_inplace(HS, SYS);                        /*in, out, in*/
    HS_STAT = 1;
  }

  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    EWALD_BGC = (double *) get_array_dval("ewald.dat", "EWALD_BGC", NRSITES);
    for (i = 0; i <= NRSITES - 1; i++)
      EWALD_BGC[i] *= CHRG_PCT;
  } else {
    EWALD_BGC = (double *) malloc(NRSITES * sizeof(double));
    for (i = 0; i <= NRSITES - 1; i++)
      EWALD_BGC[i] = 0.00;
  }
   #ifdef MPI
  if (MY_RANK == 0)
   #endif
  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    for (i = 0; i <= NRSITES - 1; i++)
      printf("\nEWALD_BGC_%s = %lf",NAMES[i], EWALD_BGC[i]); fflush(stdout);
    printf("\n");
  }

  if ((strncmp("yes", BRIDGE_FUNC1, 3)) == 0)
    UPB1 = calc_upb1();
}



int calc_upb1(void)
{
  int i, m;
  int val = 32767, num = 0;
  char s1[10];

  for (m = 0; m <= 9; m++)
    for (i = 0; i <= 9; i++) {
      sprintf(s1, "yes%d%d", m, i);
      if (strncmp(s1, BRIDGE_FUNC1, 5) == 0) {
        val = num;
        break;
      } else
        num++;
    }

  printf("\nUpdate br1 after %d iterations\n", val); fflush(stdout);

  return val;
}



void set_arrays(void)
{
  int i, j;
  char s1[100];

  /*real-space arrays*/

  CR2_S = (double **) malloc(NRSITES * sizeof(double));
  for (j = 0; j <= NRSITES - 1; j++)
    *(CR2_S + j) = (double *) malloc(NNN * sizeof(double));

  TR_S = (double **) malloc(NRSITES * sizeof(double));
  for (j = 0; j <= NRSITES - 1; j++)
    *(TR_S + j) = (double *) malloc(NNN * sizeof(double));
  printf("...done!!!\n"); fflush(stdout);

  for (j = 0; j <= NRSITES - 1; j++)
    for (i = 0; i <= NNN - 1; i++) TR_S[j][i] = 0.00;

  /*INITIAL C(r) and Prep********************************************************************/

  CR_S = (double **) malloc(NRSITES * sizeof(double));

  if (strncmp("yes", CONTINUE, 3) == 0) {
    for (j = 0; j <= NRSITES - 1; j++) {
      sprintf(s1, "cr_%s_s", NAMES[j]);
      if (get_file_stat(s1) == 1)
        *(CR_S + j) = (double *) get_3d(s1, TEMP, PND[j], SYS);
      else {
        printf("\nFile %s.%s doesn't exist\n", s1, FILE_TYPE); fflush(stdout);
        *(CR_S + j) = (double *) malloc(NNN * sizeof(double));
        for (i = 0; i <= NNN - 1; i++) CR_S[j][i] = 0.00;
      }
    }
  } else if (strncmp("no", CONTINUE, 2) == 0) {
    for (j = 0; j <= NRSITES - 1; j++) {
      *(CR_S + j) = (double *) malloc(NNN * sizeof(double));

      for (i = 0; i <= NNN - 1; i++)
        CR_S[j][i] = 0.00;
      /*if( strncmp( "yes", EWALD_SUMS, 3 ) == 0)
          for( i=0; i<=NNN-1; i++)   CR_S[j][i] = exp( -UR_S[j][i]) -1.00;*/
    }
  } else {
    fprintf(stdout, "\nConfusion of what to use as initial input for cr vectors\n");
    exit(1);
  }



  if (CR_S_STAT == 0) {
    for (j = 0; j <= NRSITES - 1; j++)
      shift_origin_inplace(*(CR_S + j), SYS);
    CR_S_STAT = 1;
  }

  printf("\n"); fflush(stdout);
}



void check_env(void)
{
  int i;

  printf("\n****************************************\n"); fflush(stdout);
  printf("		ENV			\n");   fflush(stdout);
  printf("****************************************\n");   fflush(stdout);
  printf("\nNX = %d\nNY = %d\nNZ = %d\n", NX, NY, NZ);    fflush(stdout);
  printf("\nCX = %d\nCY = %d\nCZ = %d\n", CX, CY, CZ);    fflush(stdout);
  printf("\nLX = %f\nLY = %f\nLZ = %f\n", LX, LY, LZ);    fflush(stdout);
  printf("\nSOLVER = %s \n", SOLVER);                    fflush(stdout);
  printf("T_ERR = %.5e \n", T_ERR);
  printf("PIC_MP = %f \n", PIC_MP);                      fflush(stdout);
  printf("CHRG_PCT = %f \n", CHRG_PCT);                  fflush(stdout);
  printf("TEMP_FACTOR = %f \n", TEMP_FACTOR);            fflush(stdout);
  printf("A_ERF = %f \n", A_ERF);                        fflush(stdout);
  printf("MAX_ITER = %d \n", MAX_ITER);                  fflush(stdout);
  printf("N_DUMP = %d \n", N_DUMP);
  printf("INIT_PIC_ITER = %d \n", INIT_PIC_ITER);
  printf("CONTINUE = %s \n", CONTINUE);                 fflush(stdout);
  printf("FILE_TYPE = %s \n", FILE_TYPE);                fflush(stdout);
  printf("CLOSURE = %s \n", CLOSURE);                   fflush(stdout);
  printf("EWALD_SUMS = %s \n", EWALD_SUMS);             fflush(stdout);
  printf("CONFIG_TYPE = %s \n", CONFIG_TYPE);           fflush(stdout);
  printf("BRIDGE_FUNC0 = %s \n", BRIDGE_FUNC0);         fflush(stdout);
  printf("BRIDGE_FUNC1 = %s \n", BRIDGE_FUNC1);         fflush(stdout);
  printf("RBC_FUNC = %s \n", RBC_FUNC);                 fflush(stdout);
  if (strncmp("mdiis", SOLVER, 5) == 0) {
    printf("DIIS_SIZE = %d \n", DIIS_SIZE);
    printf("DIIS_MP = %lf \n", DIIS_MP);
  }
  printf("RT_CHANGES = %d \n", RT_CHANGES);
  if (RT_CHANGES > 0) {
    printf("RT_ERR = %.5e \n", RT_ERR);
    for (i = 0; i <= RT_CHANGES - 1; i++)
      printf("\nRT_TEMP_FACTOR_%d = %lf", i, RT_TEMP_FACTOR[i]); fflush(stdout);
  }

  printf("\n****************************************\n"); fflush(stdout);
}



void check_dis(void)
{
  int i;
  printf("\n****************************************\n"); fflush(stdout);
  printf("		DIS			\n"); fflush(stdout);
  printf("****************************************\n"); fflush(stdout);
  printf("\nTEMP = %lf", TEMP); fflush(stdout);
  printf("\nTYPE = %d", TYPE); fflush(stdout);
  printf("\nNSITES = %d", NSITES); fflush(stdout);
  printf("\nNRSITES = %d", NRSITES); fflush(stdout);
  printf("\nDIS_NUM = %d", DIS_NUM); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nNAMES_%d = %s", i, NAMES[i]); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nREDUN_%s = %lf", NAMES[i], REDUN[i]); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nPND_%s = %lf", NAMES[i], PND[i]); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nEP12_%s = %lf",NAMES[i], EP12[i]); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nEP6_%s = %lf",NAMES[i], EP6[i]); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nSIG_%s = %lf", NAMES[i], SIG[i]); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++) printf("\nCHARGE_%s = %lf",NAMES[i], CHARGE[i]); fflush(stdout);
  for (i = 0; i <= DIS_NUM - 1; i++) printf("\nDIS_NAMES_%d = %s", i, DIS_NAMES[i]); fflush(stdout);
  if (TYPE == 1)
    printf("\nBOND_MAT -> %lf \t %lf", BOND_MAT[0], BOND_MAT[1]);

  printf("\n****************************************\n"); fflush(stdout);
  printf("\n\n"); fflush(stdout);
}



void check_par(U_PAR2 *u, int n_sites)
{
  int i;
  double dx = LX / NX;
  double dy = LY / NY;
  double dz = LZ / NZ;

  double x, y, z;

  FILE *out;
  if ((out = fopen("solute_check.dat", "w")) == NULL)
    fprintf(stdout, "Problem opening out file for lj parameters"); fflush(stdout);

  for (i = 1; i <= n_sites; i++) {
    fprintf(out, "%d:%s(%d)\n", u[i].num, u[i].element, u[i].mol); fflush(out);
    fprintf(out, "ep12:%f  ep6:%f\n", u[i].ep12, u[i].ep6); fflush(out);
    fprintf(out, "sig:%f\n", u[i].sig); fflush(out);
    fprintf(out, "Coulomb Charge:%f\n", u[i].charge); fflush(out);
    fprintf(out, "Cartesian Coordinates:\n"); fflush(out);
    fprintf(out, "%f\t%f\t%f\n\n", u[i].x, u[i].y, u[i].z); fflush(out);
  }
  for (i = 1; i <= n_sites; i++) {
    x = u[i].x / dx;
    y = u[i].y / dy;
    z = u[i].z / dz;

    if (x == 0.00 && y == 0.00 && z == 0.00) {
      printf("\n::ERROR - Element %d is on a grid point\n", i); fflush(stdout);
      /* exit(1);*/
    }
  }

  fclose(out);
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				 UTILITIES						*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

/*_______Spatial distance_____*/
double r0(int, int, int);                       /*from origin  x, y, z*/
double rx(int, int, int, double, double, double);  /* x, y, z, ui_x, ui_y, ui_z*/
double k0(int, int, int);
double kx(int, int, int, double, double, double);


/*___________Utilities______________*/
/*check routines*/
void check_kr(U_PAR *);
double * calc_3d_to_1d_avg(double *, double, double, double);           /* gr, x, y, z */


/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				 POTENTIAL FIELDS					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/


void set_b1r(void);
void set_b0r(double **);

int set_potential_fields(void)
{
  int nnn = NNN;
  double temp = TEMP;

  printf("\n...Reading in Potential Functions...");
  int i, j;
  char s1[50];

  /*lj*/
  double **ur_lj = (double **) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++) {
    sprintf(s1, "ur_%s_lj", NAMES[i]);
    *(ur_lj + i) = (double *)  get_3d(s1, temp, PND[i], SYS);
  }
  /*clmb*/
  double **ur_clmb = (double **) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++) {
    sprintf(s1, "ur_%s_clmb", NAMES[i]);
    *(ur_clmb + i) = (double *)  get_3d(s1, temp, PND[i], SYS);
  }

  double **ur_l = (double **) malloc(NRSITES * sizeof(double));
  fftw_complex **uk_l = (fftw_complex **) malloc(NRSITES * sizeof(fftw_complex));

  if (strncmp("no", EWALD_SUMS, 2) == 0) {
    /*ur_l*/
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1, "ur_%s_l", NAMES[i]);
      *(ur_l + i) = (double *)  get_3d(s1, temp, PND[i], SYS);
    }
    /*uk_l*/
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1, "uk_%s_l", NAMES[i]);
      *(uk_l + i) = (fftw_complex *)  get_complex_3d(s1, temp, PND[i], SYS);
    }
  }


  /*scale charged arrays*/
  if (CHRG_PCT < 1.000) {
    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= nnn - 1; i++)
        ur_clmb[j][i] *= CHRG_PCT;

    if (strncmp("no", EWALD_SUMS, 2) == 0)
      for (j = 0; j <= NRSITES - 1; j++)
        for (i = 0; i <= nnn - 1; i++) {
          ur_l[j][i] *= CHRG_PCT;
          uk_l[j][i][0] *= CHRG_PCT;
          uk_l[j][i][1] *= CHRG_PCT;
        }
  }


  if (RT_CHANGES > 0) {
    if (RT_TEMP_FACTOR != NULL) {
      printf("\nStarting with a temp_factor = %f\n", RT_TEMP_FACTOR[0]);

      for (j = 0; j <= NRSITES - 1; j++) {
        for (i = 0; i <= nnn - 1; i++) ur_lj[j][i] /= RT_TEMP_FACTOR[0];
        for (i = 0; i <= nnn - 1; i++) ur_clmb[j][i] /= RT_TEMP_FACTOR[0];
      }
      if (strncmp("no", EWALD_SUMS, 2) == 0)
        for (j = 0; j <= NRSITES - 1; j++)
          for (i = 0; i <= nnn - 1; i++) {
            ur_l[j][i] /= RT_TEMP_FACTOR[0];
            uk_l[j][i][0] /= RT_TEMP_FACTOR[0];
            uk_l[j][i][1] /= RT_TEMP_FACTOR[0];
          }
    }
  } else {
    if (TEMP_FACTOR > 1.00) {
      for (j = 0; j <= NRSITES - 1; j++) {
        for (i = 0; i <= nnn - 1; i++) ur_lj[j][i] /= TEMP_FACTOR;
        for (i = 0; i <= nnn - 1; i++) ur_clmb[j][i] /= TEMP_FACTOR;
      }

      if (strncmp("no", EWALD_SUMS, 2) == 0)
        for (j = 0; j <= NRSITES - 1; j++)
          for (i = 0; i <= nnn - 1; i++) {
            ur_l[j][i] /= TEMP_FACTOR;
            uk_l[j][i][0] /= TEMP_FACTOR;
            uk_l[j][i][1] /= TEMP_FACTOR;
          }
    }
  }


  double **ur = (double **) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++)
    *(ur + i) = add_arrays(*(ur_lj + i), *(ur_clmb + i), nnn);

  double **ur_s = (double **) malloc(NRSITES * sizeof(double));

  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (i = 0; i <= NRSITES - 1; i++)
      *(ur_s + i) = sub_arrays(*(ur + i), *(ur_l + i), nnn);

  else if (strncmp("yes", EWALD_SUMS, 3) == 0)
    ur_s = ur;

  /*print ur*/
  /*for( i=0; i<=NRSITES-1; i++){
          sprintf( s1, "ur_%s", NAMES[i] );
          print_3d( s1, *(ur +i));
     }*/

  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (i = 0; i <= NRSITES - 1; i++)
      shift_origin_complex_inplace(*(uk_l + i), SYS);
  UK_L_STAT = 1;

  /*GLOBAL*/
  UR_S = ur_s;
  UR_L = ur_l;
  UK_L = uk_l;

  for (i = 0; i <= NRSITES - 1; i++)
    shift_origin_inplace(*(UR_S + i), SYS);
  UR_S_STAT = 1;

  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (i = 0; i <= NRSITES - 1; i++)
      shift_origin_inplace(*(UR_L + i), SYS);
  UR_L_STAT = 1;

  printf("...Done!!!\n");
  /*bridge_functions*/
  if (strncmp("yes", BRIDGE_FUNC1, 3) == 0) {
    set_b1r();
    B1R2 = 1;
  } else if (strncmp("yes", BRIDGE_FUNC0, 3) == 0) {
    set_b0r(ur_lj);
    B1R2 = 1;
  }
  B1R_STAT = 1;


  /*DECOMMISION*/
  for (i = 0; i <= NRSITES - 1; i++) {
    free(*(ur_lj + i));
    free(*(ur_clmb + i));
  }
  free(ur_lj);
  free(ur_clmb);

  if (strncmp("no", EWALD_SUMS, 2) == 0) {
    for (i = 0; i <= NRSITES - 1; i++)
      free(*(ur + i));
    free(ur);
  }

  return 0;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*											*/
/*		                bridge functions					*/
/*											*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void set_b1r(void)
{
  int i, j;
  char s1[100];
  int stat = 0;
  /*BRIDGE_FUNC1*/
  B1R = (double **) malloc(NRSITES * sizeof(double));
  if (strncmp("yes", CONTINUE, 3) == 0) {
    if (strncmp("yes00", BRIDGE_FUNC1, 5) == 0) {
      for (j = 0; j <= NRSITES - 1; j++) {
        sprintf(s1, "gr_%s", NAMES[j]);
        if (get_file_stat(s1) == 1) {
          *(B1R + j) = (double *) get_3d(s1, TEMP, PND[j], SYS);
          shift_origin_inplace(*(B1R + j), SYS);
        } else {
          printf("\nFile %s.%s doesn't exist:no bridge function\n", s1, FILE_TYPE); fflush(stdout);
          *(B1R + j) = (double *) malloc(NNN * sizeof(double));
          for (i = 0; i <= NNN - 1; i++) B1R[j][i] = 0.00;
          stat++;
        }
      }
      if (stat == 0) {
        printf("\nCalculating BR1 from gr files\n"); fflush(stdout);
        calc_bridge_func1(0);
      }
    } else if (strncmp("yes", BRIDGE_FUNC1, 3) == 0) {
      for (j = 0; j <= NRSITES - 1; j++) {
        sprintf(s1, "b1r_%s", NAMES[j]);
        if (get_file_stat(s1) == 1) {
          *(B1R + j) = (double *) get_3d(s1, TEMP, PND[j], SYS);
          shift_origin_inplace(*(B1R + j), SYS);
        } else {
          printf("\nFile %s.%s doesn't exist\n", s1, FILE_TYPE); fflush(stdout);
          *(B1R + j) = (double *) malloc(NNN * sizeof(double));
          for (i = 0; i <= NNN - 1; i++) B1R[j][i] = 0.00;
        }
      }
    }
  } else {
    for (j = 0; j <= NRSITES - 1; j++) {
      *(B1R + j) = (double *) malloc(NNN * sizeof(double));
      for (i = 0; i <= NNN - 1; i++)
        B1R[j][i] = 0.00;
    }
  }
}



void set_b0r(double ** ur_lj)
{
  int i, j;
  char s1[100];
  double **ur_lj12 = NULL;

  B1R = (double **) malloc(NRSITES * sizeof(double));

  for (j = 0; j <= NRSITES - 1; j++) {
    *(B1R + j) = (double *) malloc(NNN * sizeof(double));
    for (i = 0; i <= NNN - 1; i++)
      B1R[j][i] = 0.00;
  }

  if ((strncmp("yes0", BRIDGE_FUNC0, 4)) == 0) {
    printf("\nCalculating BR(ur_lj12) function\n"); fflush(stdout);
    ur_lj12 = (double **) malloc(NRSITES * sizeof(double));
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1, "ur_%s_lj12", NAMES[i]);
      *(ur_lj12 + i) = (double *)  get_3d(s1, TEMP, PND[i], SYS);
    }

    for (i = 0; i <= NRSITES - 1; i++)
      shift_origin_inplace(*(ur_lj12 + i), SYS);
    calc_bridge_func1_fbond(ur_lj12);
  } else if ((strncmp("yes1", BRIDGE_FUNC0, 4)) == 0) {
    printf("\nCalculating BR(ur_lj) function\n"); fflush(stdout);

    for (i = 0; i <= NRSITES - 1; i++)
      shift_origin_inplace(*(ur_lj + i), SYS);
    calc_bridge_func1_fbond(ur_lj);
  } else if ((strncmp("yes2", BRIDGE_FUNC0, 4)) == 0) {
    printf("\nCalculating BR(ur_s) function\n"); fflush(stdout);
    calc_bridge_func1_fbond(UR_S);
  } else {
    printf("\n<<< Error choosing bridge_func0 >>>\n");
  }


  if ((strncmp("yes0", BRIDGE_FUNC0, 4)) == 0) {
    for (i = 0; i <= NRSITES - 1; i++)
      free(*(ur_lj12 + i));
    free(ur_lj12);
  }
}



void calc_bridge_func1(int stat)
{
  int i, k, m, m1, m2, n;

  printf("\n\n<<< updating b1r >>>\n\n"); fflush(stdout);

  fftw_complex **hr = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  tk [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(hr + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex **hk = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  tk [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(hk + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex ***tkp = (fftw_complex ***) fftw_malloc(NRSITES * sizeof(fftw_complex));
  fftw_complex ***trp = (fftw_complex ***) fftw_malloc(NRSITES * sizeof(fftw_complex));
  for (m = 0; m <= NRSITES - 1; m++) {
    *(tkp + m) = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));
    *(trp + m) = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));
  }

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= NRSITES - 1; n++) {
      *(*(tkp + m) + n) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));
      *(*(trp + m) + n) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));
    }
  if (stat == 1) {
    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= NNN - 1; i++) {
        hr[m][i][0] = exp(-UR_S[m][i] + TR_S[m][i] + B1R[m][i]) - 1.00;
        hr[m][i][1] = 0.00;
      }
  } else if (stat == 0) {
    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= NNN - 1; i++) {
        hr[m][i][0] = B1R[m][i] - 1.00;
        hr[m][i][1] = 0.00;
      }
  }


  for (m = 0; m <= NRSITES - 1; m++)
    fftw_3d(*(hr + m),*(hk + m));

  if (TYPE == 1 || TYPE == 2) {
    /*O-O term [0][0]*/
    for (i = 0; i <= NNN - 1; i++)
      for (k = 0; k <= 1; k++)
        tkp[0][0][i][k] = hk[0][i][k] * PND[0] * HK3[0][0][i];

    /* O-H term [0][1] and [1][0] */
    for (i = 0; i <= NNN - 1; i++)
      for (k = 0; k <= 1; k++) {
        tkp[0][1][i][k] = hk[0][i][k] * PND[0] * (WK_OH[i] + HK3[0][1][i]);
        tkp[1][0][i][k] = hk[1][i][k] * PND[1] * (WK_OH[i] + HK3[1][0][i]);
      }

    /* H-H term [1][1] and [1][1] */
    for (i = 0; i <= NNN - 1; i++)
      for (k = 0; k <= 1; k++)
        tkp[1][1][i][k] = hk[1][i][k] * PND[1] * (WK_HH[i] + HK3[1][1][i]);

    /* All other terms */
    for (m = 0; m <= NRSITES - 1; m++)
      for (n = 0; n <= NRSITES - 1; n++) {
        if (m >= 2 || n >= 2)
          for (i = 0; i <= NNN - 1; i++)
            for (k = 0; k <= 1; k++) {
              tkp[m][n][i][k] = hk[m][i][k] * PND[m] * HK3[m][n][i];
            }
      }
  } else {
    for (m = 0; m <= NRSITES - 1; m++)
      for (n = 0; n <= NRSITES - 1; n++)
        for (i = 0; i <= NNN - 1; i++)
          for (k = 0; k <= 1; k++) {
            tkp[m][n][i][k] = hk[m][i][k] * PND[m] * HK3[m][n][i];
          }
  }

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= NRSITES - 1; n++)
      invfftw_3d(*(*(tkp + m) + n), *(*(trp + m) + n));


  for (n = 0; n <= NRSITES - 1; n++) {
    for (i = 0; i <= NNN - 1; i++)
      B1R[n][i] = 0.00;

    for (m1 = 0; m1 <= NRSITES - 1; m1++)
      for (m2 = 0; m2 <= NRSITES - 1; m2++)
        for (i = 0; i <= NNN - 1; i++)
          B1R[n][i] += (-1.00) * (0.5) * (REDUN[m1] * REDUN[m2]) * trp[m1][n][i][0] * trp[m2][n][i][0];
  }


  for (m1 = 0; m1 <= NRSITES - 1; m1++)
    for (m2 = 0; m2 <= NRSITES - 1; m2++) {
      free(*(*(trp + m1) + m2));
      free(*(*(tkp + m1) + m2));
    }

  for (m = 0; m <= NRSITES - 1; m++) {
    free(*(hr + m));
    free(*(hk + m));
    free(*(trp + m));
    free(*(tkp + m));
  }

  free(hr);
  free(hk);
  free(trp);
  free(tkp);

  B1R_STAT = 1;
}



void calc_bridge_func1_fbond(double **ur)
{
  int i, k, m, m1, m2, n;

  fftw_complex **fr = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  tk [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(fr + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex **fk = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  tk [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(fk + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex ***tkp = (fftw_complex ***) fftw_malloc(NRSITES * sizeof(fftw_complex));
  fftw_complex ***trp = (fftw_complex ***) fftw_malloc(NRSITES * sizeof(fftw_complex));
  for (m = 0; m <= NRSITES - 1; m++) {
    *(tkp + m) = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));
    *(trp + m) = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));
  }

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= NRSITES - 1; n++) {
      *(*(tkp + m) + n) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));
      *(*(trp + m) + n) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));
    }

  for (m = 0; m <= NRSITES - 1; m++)
    for (i = 0; i <= NNN - 1; i++) {
      fr[m][i][0] = exp(-1.0 * (ur[m][i])) - 1;
      fr[m][i][1] = 0.00;
    }

  for (m = 0; m <= NRSITES - 1; m++)
    fftw_3d(*(fr + m),*(fk + m));

  if (TYPE == 1 || TYPE == 2) {
    /*O-O term [0][0]*/
    for (i = 0; i <= NNN - 1; i++)
      for (k = 0; k <= 1; k++)
        tkp[0][0][i][k] = fk[0][i][k] * PND[0] * HK3[0][0][i];

    /* O-H term [0][1] and [1][0] */
    for (i = 0; i <= NNN - 1; i++)
      for (k = 0; k <= 1; k++) {
        tkp[0][1][i][k] = fk[0][i][k] * PND[0] * (WK_OH[i] + HK3[0][1][i]);
        tkp[1][0][i][k] = fk[1][i][k] * PND[1] * (WK_OH[i] + HK3[1][0][i]);
      }

    /* H-H term [1][1] and [1][1] */
    for (i = 0; i <= NNN - 1; i++)
      for (k = 0; k <= 1; k++)
        tkp[1][1][i][k] = fk[1][i][k] * PND[1] * (WK_HH[i] + HK3[1][1][i]);

    /* All other terms */
    for (m = 0; m <= NRSITES - 1; m++)
      for (n = 0; n <= NRSITES - 1; n++) {
        if (m >= 2 || n >= 2)
          for (i = 0; i <= NNN - 1; i++)
            for (k = 0; k <= 1; k++) {
              tkp[m][n][i][k] = fk[m][i][k] * PND[m] * HK3[m][n][i];
            }
      }
  } else {
    for (m = 0; m <= NRSITES - 1; m++)
      for (n = 0; n <= NRSITES - 1; n++)
        for (i = 0; i <= NNN - 1; i++)
          for (k = 0; k <= 1; k++) {
            tkp[m][n][i][k] = fk[m][i][k] * PND[m] * HK3[m][n][i];
          }
  }

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= NRSITES - 1; n++)
      invfftw_3d(*(*(tkp + m) + n), *(*(trp + m) + n));


  for (n = 0; n <= NRSITES - 1; n++) {
    for (i = 0; i <= NNN - 1; i++)
      B1R[n][i] = 0.00;

    for (m1 = 0; m1 <= NRSITES - 1; m1++)
      for (m2 = 0; m2 <= NRSITES - 1; m2++)
        for (i = 0; i <= NNN - 1; i++)
          B1R[n][i] += (-1.00) * (0.5) * (REDUN[m1] * REDUN[m2]) * trp[m1][n][i][0] * trp[m2][n][i][0];
  }

  for (m1 = 0; m1 <= NRSITES - 1; m1++)
    for (m2 = 0; m2 <= NRSITES - 1; m2++) {
      free(*(*(trp + m1) + m2));
      free(*(*(tkp + m1) + m2));
    }

  for (m = 0; m <= NRSITES - 1; m++) {
    free(*(fr + m));
    free(*(fk + m));
    free(*(trp + m));
    free(*(tkp + m));
  }

  free(fr);
  free(fk);
  free(trp);
  free(tkp);

  B1R_STAT = 1;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*											*/
/*			repulsive bridge correction					*/
/*											*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

int calc_rbc(void)
{
  /*EXTERN*/
  double temp = TEMP;
  /*EXTERN*/

  int i, m;
  char s1[50];

  /*lj12*/
  double **ur_lj12 = (double **) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++) {
    sprintf(s1, "ur_%s_lj12", NAMES[i]);
    *(ur_lj12 + i) = (double *)  get_3d(s1, temp, PND[i], SYS);
  }

  fftw_complex **exp_ur = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                /*  [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(exp_ur + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));
  fftw_complex **exp_uk = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                /*  [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(exp_uk + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex *exp_uk_w1 = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));           /*  [idx][r,c] */
  fftw_complex *exp_ur_w1 = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));           /*  [idx][r,c] */
  fftw_complex *exp_uk_w2 = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));           /*  [idx][r,c] */
  fftw_complex *exp_ur_w2 = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));           /*  [idx][r,c] */

  double **exp_br = (double **) fftw_malloc(NRSITES * sizeof(double));                  /*  [site][idx] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(exp_br + m) = (double  *) fftw_malloc(NNN * sizeof(double));

  for (m = 0; m <= NRSITES - 1; m++)
    for (i = 0; i <= NNN - 1; i++) {
      exp_ur[m][i][0] = exp(-1.0 * ur_lj12[m][i]);
      exp_ur[m][i][1] = 0.00;
    }

  /* shift exp_u(r) */
  for (m = 0; m <= NRSITES - 1; m++)
    shift_origin_complex_inplace(*(exp_ur + m), SYS);                       /*in, out, in*/

  for (m = 0; m <= NRSITES - 1; m++)
    fftw_3d(*(exp_ur + m), *(exp_uk + m));                       /*in, out*/

  /**** exp(b(r)) - oxygen ****/
  for (i = 0; i <= NNN - 1; i++) {
    exp_uk_w1[i][0] = exp_uk[1][i][0] * WK_OH[i];
    exp_uk_w1[i][1] = 0.0000;
  }

  invfftw_3d(exp_uk_w1, exp_ur_w1);

  for (i = 0; i <= NNN - 1; i++)
    exp_br[0][i] = pow(exp_ur_w1[i][0], 2);


  /**** exp(b(r)) - hydrogen ****/
  for (i = 0; i <= NNN - 1; i++) {
    exp_uk_w1[i][0] = exp_uk[1][i][0] * WK_HH[i];
    exp_uk_w1[i][1] = 0.000;
    exp_uk_w2[i][0] = exp_uk[0][i][0] * WK_OH[i];
    exp_uk_w2[i][1] = 0.000;
  }

  invfftw_3d(exp_uk_w1, exp_ur_w1);
  invfftw_3d(exp_uk_w2, exp_ur_w2);

  for (i = 0; i <= NNN - 1; i++)
    exp_br[1][i] = exp_ur_w1[i][0] * exp_ur_w2[i][0];

  for (m = 2; m <= NRSITES - 1; m++)
    for (i = 0; i <= NNN - 1; i++)
      exp_br[m][i] = 1.00000;

  EXP_BR = exp_br;

  for (m = 0; m <= NRSITES - 1; m++)
    unshift_origin_inplace(*(EXP_BR + m), SYS);                       /*in, out, in*/
   #ifdef MPI
  if (MY_RANK == 0)
   #endif
  for (m = 0; m <= NRSITES - 1; m++) {
    sprintf(s1, "br_%s", NAMES[m]);
    print_3d(s1, *(EXP_BR + m));
  }

  for (m = 0; m <= NRSITES - 1; m++)
    shift_origin_inplace(*(EXP_BR + m), SYS);                       /*in, out, in*/

  /*free*/
  for (m = 0; m <= NRSITES - 1; m++) {
    free(*(ur_lj12 + m));
    free(*(exp_ur + m));
    free(*(exp_uk + m));
  }

  free(ur_lj12);
  free(exp_ur);
  free(exp_uk);
  free(exp_ur_w1);
  free(exp_uk_w1);
  free(exp_ur_w2);
  free(exp_uk_w2);

  EXP_BR_STAT = 1;

  return 0;
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
/*			Intramolecular correlatin functions				*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/


double * calc_wk(double);

int calc_intramolecular_functions(void)
{
  double *wk_oh = calc_wk(BOND_MAT[0]);                 /*replaces L_OH*/
  double *wk_hh = calc_wk(BOND_MAT[1]);                 /*replaces L_HH */

  shift_origin_inplace(wk_oh, SYS);
  shift_origin_inplace(wk_hh, SYS);

  /*GLOBAL*/
  WK_OH = wk_oh;
  WK_HH = wk_hh;

  WK_STAT = 1;

  return 0;
}



/****************************************************************************************/
/*											*/
/*				Intramolecular ( WK ) function				*/
/*											*/
/****************************************************************************************/


double * calc_wk(double l_xx)
{
  /*EXTERN*/
  int nnn = NNN, nx = NX, ny = NY, nz = NZ;
  /*EXTERN*/

  int x, y, z;
  double k;
  double *wk = (double *) malloc(nnn * sizeof(double));

  for (x = 0; x <= nx - 1; x++)
    for (y = 0; y <= ny - 1; y++)
      for (z = 0; z <= nz - 1; z++) {
        k = k0(x, y, z);
        wk[ii(x,y,z)] = sin(k * l_xx) / (k * l_xx);
      }

  wk[ii(CX,CY,CZ)] = 1.00;

  return wk;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				HK_vv 3D correlations					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

double * calc_1d_to_3d_shift_fft(double *, int, double);
double * calc_1d_to_3d_shift(double *, int, double);

int calc_hk_vv_correlations(char hr_fname[])
{
  int i, j, k, slen;
  char s1[50], s2[50];

  char * ext = (char *) get_ext(hr_fname);
  int n_pts = (int) get_dval(hr_fname, "N_PTS");
  double rad = (double) get_dval(hr_fname, "RADIUS");

  double **hr1vv = (double **) malloc(DIS_NUM * sizeof(double));

  for (j = 0; j <= DIS_NUM - 1; j++)
    *(hr1vv + j) = get_dis(hr_fname, j + 1);

  printf("\n...done reading in solvent...\n"); fflush(stdout);

  if (MY_RANK == 0)
    for (j = 0; j <= DIS_NUM - 1; j++) {
      if (hr1vv[j] == NULL) {
        fprintf(stderr, "hr1vv[%d] is NULL\n", j);
        exit(1);
      }
      if (strncmp("dis2", ext, 3) == 0)
        sprintf(s1, "hrvv1d_%s.1d.dat", DIS_NAMES[j]);
      if (strncmp("kdis2", ext, 4) == 0)
        sprintf(s1, "hkvv1d_%s.1d.dat", DIS_NAMES[j]);
    }

  double **hk3 = (double **) malloc(DIS_NUM * sizeof(double));
  if (strncmp("dis2", ext, 3) == 0) {
    for (j = 0; j <= DIS_NUM - 1; j++)
      *(hk3 + j) = calc_1d_to_3d_shift_fft(*(hr1vv + j), n_pts, rad);
  } else if (strncmp("kdis2", ext, 4) == 0) {
    for (j = 0; j <= DIS_NUM - 1; j++)
      *(hk3 + j) = calc_1d_to_3d_shift(*(hr1vv + j), n_pts, rad);
  }

  double ***h3 = ( double ***) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++)
    *(h3 + i) = (double **) malloc(NRSITES * sizeof(double));

  for (i = 0; i <= NRSITES - 1; i++) {
    for (j = i; j <= NRSITES - 1; j++) {
      sprintf(s1, "%s-%s", NAMES[i], NAMES[j]);
      sprintf(s2, "%s-%s", NAMES[j], NAMES[i]);

      slen = strlen(s1);
      for (k = 0; k <= DIS_NUM - 1; k++)
        if ((strncmp(s1, DIS_NAMES[k], slen) == 0) || (strncmp(s2, DIS_NAMES[k], slen) == 0)) {
          h3[i][j] = *(hk3 + k);
          h3[j][i] = h3[i][j];
          printf("\n%d:Calculating h(k) for %s", MY_RANK, DIS_NAMES[k]); fflush(stdout);
          break;
        } else if (k == (DIS_NUM - 1))
          printf("\nh(r) distribution %s or %s not found\n", s1, s2);
    }
  }

  /* set global variable */
  HK3 = h3;

  HK3_STAT = 1;

  return 0;
}



/****************************************************************************************/
/****************************************************************************************/


double * calc_1d_to_3d_shift_fft(double *hr1d, int n_vv, double r_vv)
{
  int nx = NX, ny = NY, nz = NZ, nnn = NNN;             /*EXTERN*/

  int i, x, y, z;
  int indx;
  double r3, rmd;

  double dr = r_vv / (double) (n_vv - 1);
  double *hk3d = (double *) malloc(nnn * sizeof(double));
  fftw_complex *rk1 = (fftw_complex *) malloc(nnn * sizeof(fftw_complex));
  fftw_complex *rk2 = (fftw_complex *) malloc(nnn * sizeof(fftw_complex));

  for (x = 0; x <= nx - 1; x++)
    for (y = 0; y <= ny - 1; y++)
      for (z = 0; z <= nz - 1; z++) {
        r3 = r0(x, y, z);
        indx = r3 / dr;
        if (indx > (n_vv - 1)) {
          rk1[ii(x,y,z)][0] = 0.00; rk1[ii(x,y,z)][1] = 0.00; continue;
        }
        rmd = r3 - (indx * dr);
        rk1[ii(x,y,z)][0] = (dr - rmd) / dr * hr1d[indx] + (rmd / dr) * hr1d[indx + 1];
        rk1[ii(x,y,z)][1] = 0.00;
      }

  shift_origin_complex(rk1, rk2, SYS);
  fftw_3d(rk2, rk1);

  for (i = 0; i <= nnn - 1; i++) hk3d[i] = rk1[i][0];

  free(rk1);
  free(rk2);

  return hk3d;
}



double * calc_1d_to_3d_shift(double *hk1d, int n_vv, double r_vv)
{
  int nx = NX, ny = NY, nz = NZ, nnn = NNN;             /*EXTERN*/

  int x, y, z;
  int indx;
  double k3, rmd;

  double dk = (double) Pi / r_vv;
  double *hk3d = (double *) malloc(nnn * sizeof(double));

  for (x = 0; x <= nx - 1; x++)
    for (y = 0; y <= ny - 1; y++)
      for (z = 0; z <= nz - 1; z++) {
        k3 = k0(x, y, z);
        indx = k3 / dk;
        if (indx > (n_vv - 1)) {
          hk3d[ii(x,y,z)] = 0.00; continue;
        }
        rmd = k3 - (indx * dk);
        hk3d[ii(x,y,z)] = (dk - rmd) / dk * hk1d[indx] + (rmd / dk) * hk1d[indx + 1];
      }

  shift_origin_inplace(hk3d, SYS);

  return hk3d;
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
/*				FULL PICARD ITERATIONS					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void closure_cr(fftw_complex **, fftw_complex **);

void full_picard_iter(int max_iter)
{
  int nnn = NNN;

  double *hk3_oo = NULL;                /* [row][col][idx] */
  double *hk3_oh = NULL;
  double *hk3_hh = NULL;

  if (TYPE == 1 || TYPE == 2) {
    hk3_oo = HK3[0][0];                 /* [row][col][idx] */
    hk3_oh = HK3[0][1];
    hk3_hh = HK3[1][1];
  }

  fftw_complex **uk_l = UK_L;           /* [site][idx][r-c] */
  /*END EXTERN*/

  int i, k, m, n, counter = 1;
  double Test;
  char s1[100];

  double *t_vec = (double *) malloc(NNN * sizeof(double));

  fftw_complex **cr_s = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                  /*  cr_s [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(cr_s + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex **tr = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  tr [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(tr + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex **ck = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  ck [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(ck + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));

  fftw_complex **tk = (fftw_complex **) fftw_malloc(NRSITES * sizeof(fftw_complex));                    /*  tk [site][idx][r,c] */
  for (m = 0; m <= NRSITES - 1; m++)
    *(tk + m) = (fftw_complex *) fftw_malloc(NNN * sizeof(fftw_complex));


  /*  Initialize u(r)  */

  for (m = 0; m <= NRSITES - 1; m++) {
    for (i = 0; i <= nnn - 1; i++) {
      cr_s[m][i][0] = CR_S[m][i];
      cr_s[m][i][1] = 0.00;
    }                                                   /*cr = CR*/
  }

  /*#### picard iter ####*/
  /*#### picard iter ####*/
  /*#### picard iter ####*/
  /**********************************************************************************************/

  do {
    PIC_CNT++;
    if (MY_RANK == 0)
      printf("cr[0] = %f\n", cr_s[0][0][0]);

    /* set Im[ c(r) ] = 0 */
    for (m = 0; m <= NRSITES - 1; m++) {
      {
        for (i = 0; i <= nnn - 1; i++)
          cr_s[m][i][1] = 0.00;

        /*FFT cr_uo_s and cr_uh_s*/
        fftw_3d(*(cr_s + m), *(ck + m));                                 /*in, out*/

        /*Long Range*/
        if (strncmp("no", EWALD_SUMS,2) == 0)
          for (i = 0; i <= nnn - 1; i++) {
            ck[m][i][0] = ck[m][i][0] - uk_l[m][i][0];
            ck[m][i][1] = ck[m][i][1] - uk_l[m][i][1];
          }
      }
    }

    /*OZ EQUATION***********************************************************************************/
    /*OZ EQUATION***********************************************************************************/

    if (TYPE == 0) {
      for (k = 0; k <= 1; k++) {
        for (m = 0; m <= NRSITES - 1; m++) {
          {
            for (i = 0; i <= NNN - 1; i++)
              tk[m][i][k] = 0.00;

            for (n = 0; n <= NRSITES - 1; n++)
              for (i = 0; i <= NNN - 1; i++)
                tk[m][i][k] += ck[n][i][k] * PND[n] * HK3[n][m][i];
          }
        }
      }
    } else if (TYPE == 1) {
      for (i = 0; i <= nnn - 1; i++) {
        tk[0][i][0] = ck[0][i][0] * (PND[0] * hk3_oo[i]);                                                       /* O-O */
        tk[0][i][0] += 2 * ck[1][i][0] * (WK_OH[i] + PND[1] * hk3_oh[i]);                                       /* H-O, H2-O */
        tk[0][i][1] = ck[0][i][1] * (PND[0] * hk3_oo[i]);                                                       /* O-O */
        tk[0][i][1] += 2 * ck[1][i][1] * (WK_OH[i] + PND[1] * hk3_oh[i]);                               /* H-O, H2-O */
      }
      for (i = 0; i <= nnn - 1; i++) {
        tk[1][i][0] = ck[0][i][0] * (WK_OH[i] + PND[0] * hk3_oh[i]);                                           /* O-H */
        tk[1][i][0] += ck[1][i][0] * (WK_HH[i] + 2 * PND[1] * hk3_hh[i]);                                       /* H-H, H2-H */
        tk[1][i][1] = ck[1][i][1] * (WK_HH[i] + 2 * PND[0] * hk3_hh[i]);                                        /* O-H */
        tk[1][i][1] += ck[0][i][1] * (WK_OH[i] + PND[1] * hk3_oh[i]);                                           /* H-H, H2-H */
      }
    } else if (TYPE == 2) {
      /* k=0-real, k=1-imag */
      for (k = 0; k <= 1; k++) {
        /* H20 part */
        for (i = 0; i <= nnn - 1; i++) {
          tk[0][i][k] = ck[0][i][k] * (PND[0] * hk3_oo[i]);                                                             /* O-O */
          tk[0][i][k] += 2 * ck[1][i][k] * (WK_OH[i] + PND[1] * hk3_oh[i]);                                     /* H-O, H2-O */
        }
        for (i = 0; i <= nnn - 1; i++) {
          tk[1][i][k] = ck[0][i][k] * (WK_OH[i] + PND[0] * hk3_oh[i]);                                                  /* O-H */
          tk[1][i][k] += ck[1][i][k] * (WK_HH[i] + 2 * PND[1] * hk3_hh[i]);                                             /* H-H, H2-H */
        }

        /* other contributors to O and H part*/
        for (n = 2; n <= NRSITES - 1; n++)
          for (i = 0; i <= nnn - 1; i++) {
            tk[0][i][k] += REDUN[n] * ck[n][i][k] * PND[n] * HK3[n][0][i];
            tk[1][i][k] += REDUN[n] * ck[n][i][k] * PND[n] * HK3[n][1][i];
          }

        for (m = 2; m <= NRSITES - 1; m++)
          for (n = 0; n <= NRSITES - 1; n++)
            for (i = 0; i <= nnn - 1; i++)
              tk[m][i][k] += REDUN[n] * ck[n][i][k] * PND[n] * HK3[n][m][i];
      }
    }

    /*END OF OZ*************************************************************************************/
    /*END OF OZ*************************************************************************************/


    /*Subtract long range pot. int k space, tk -> tk_s*/
    if (strncmp("no", EWALD_SUMS,2) == 0)
      for (m = 0; m <= NRSITES - 1; m++)
        for (i = 0; i <= nnn - 1; i++) {
          tk[m][i][0] = tk[m][i][0] - uk_l[m][i][0];
          tk[m][i][1] = tk[m][i][1] - uk_l[m][i][1];
        }

    /*Inv_FFT */
    for (m = 0; m <= NRSITES - 1; m++)
      invfftw_3d(*(tk + m), *(tr + m));                                 /*in, out*/

    /****CLOSURE****/

    closure_cr(cr_s, tr);

    /****CLOSURE****/

    /*Residual vec*/
    Test = 0.00;

    for (m = 0; m <= NRSITES - 1; m++) {
      for (i = 0; i <= NNN - 1; i++)
        t_vec[i] = cr_s[m][i][1] - cr_s[m][i][0];
      Test += pow(svector_norm(t_vec, NNN), 2) / (double) NRSITES;
    }

    Test = sqrt(Test);

    if (MY_RANK == 0)
      printf("Picard Test[%d] = %.10e\t\t", PIC_CNT, Test);
    FERR = Test;

    /*new = mp*new + (1-mp)*old */
    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= NNN - 1; i++)
        cr_s[m][i][0] = (PIC_MP)*cr_s[m][i][1] + (1.0 - PIC_MP) * cr_s[m][i][0];


    if ((PIC_CNT % N_DUMP) == 0) {
      printf("\n\n<<< DUMPING cr_s with temp factor of: %f >>>\n", TEMP_FACTOR);
      /*cr are shifted therefore CR must be unshifted*/
      for (m = 0; m <= NRSITES - 1; m++)
        for (i = 0; i <= NNN - 1; i++) {
          CR_S[m][i] = cr_s[m][i][0];
        }
      if (CR_S_STAT != 0)
        for (m = 0; m <= NRSITES - 1; m++)
          unshift_origin_inplace(*(CR_S + m), SYS);                                               /*(in, out), in*/

      if (MY_RANK == 0)
        for (m = 0; m <= NRSITES - 1; m++) {
          sprintf(s1, "cr_%s_s", NAMES[m]);
          print_3d(s1, *(CR_S + m));
        }
      for (m = 0; m <= NRSITES - 1; m++)
        shift_origin_inplace(*(CR_S + m), SYS);                                                 /*(in, out), in*/
      CR_S_STAT = 1;
    }

    ++counter;

    /** RT_CHANGES and BRIDGE_FUNC*/
    if (strncmp("picard", SOLVER, 6) == 0) {
      if (Test < RT_ERR && RT_CNT < RT_CHANGES) {
        if (strncmp("yes00", BRIDGE_FUNC1, 5) == 0)
          calc_bridge_func1(1);
        change_RT();
        Test = T_ERR + 10.0;                                 /*needed if RT_ERR = T_ERR, premature exit*/
        UPB1_STAT = 0;                                       /* give some time before updating b1r*/
        IERR = 9e9;
      } else if (strncmp("yes", BRIDGE_FUNC1, 3) == 0 && UPB1_STAT >= UPB1)
        if (UPB1 >= 1)
          if (Test < IERR) {
            for (m = 0; m <= NRSITES - 1; m++)
              for (i = 0; i <= NNN - 1; i++)
                TR_S[m][i] = tr[m][i][0];
            calc_bridge_func1(1);
            IERR = Test;
            UPB1_STAT = 0;
          }
      UPB1_STAT++;
    }
  } while (Test >= T_ERR && counter <= max_iter);

  CNT = counter;

  if (Test >= T_ERR && counter <= max_iter) STAT = 0;
  if (Test >= T_ERR && counter >= max_iter) STAT = 1;
  if (Test <= T_ERR && counter <= max_iter) STAT = 2;

  for (m = 0; m <= NRSITES - 1; m++)
    for (i = 0; i <= NNN - 1; i++) {
      CR_S[m][i] = cr_s[m][i][0];
      CR2_S[m][i] = cr_s[m][i][1];
      TR_S[m][i] = tr[m][i][0];
    }
  TR_S_STAT = 1;
  CR_S_STAT = 1;

  for (m = 0; m <= NRSITES - 1; m++) {
    free(*(cr_s + m));
    free(*(tr + m));
    free(*(ck + m));
    free(*(tk + m));
  }

  free(cr_s);
  free(tr);
  free(ck);
  free(tk);
  free(t_vec);
}



void closure_cr(fftw_complex **cr_s, fftw_complex **tr)
{
  int m, i;
  double wr;

  if (UR_S_STAT != 1) {
    for (m = 0; m <= NRSITES - 1; m++)
      shift_origin_inplace(*(UR_S + m), SYS);
    UR_S_STAT = 1;
  }


  if ((strncmp("hnc", CLOSURE, 3) == 0) && (strlen(CLOSURE) == 3)) {
    if (strncmp("no", RBC_FUNC, 2) == 0) {
      if (B1R2 == 0) {
        for (m = 0; m <= NRSITES - 1; m++)
          for (i = 0; i <= NNN - 1; i++)
            cr_s[m][i][1] = exp(-UR_S[m][i] + tr[m][i][0]) - tr[m][i][0] - 1;
      } else if (B1R2 == 1) {
        for (m = 0; m <= NRSITES - 1; m++)
          for (i = 0; i <= NNN - 1; i++)
            cr_s[m][i][1] = exp(-UR_S[m][i] + tr[m][i][0] + B1R[m][i]) - tr[m][i][0] - 1;
      } else {
        printf("\nError on bridge_func1\n");
        exit(1);
      }
    } else if ((strncmp("yes", RBC_FUNC, 3) == 0)) {
      if (B1R2 == 0) {
        for (m = 0; m <= NRSITES - 1; m++)
          for (i = 0; i <= NNN - 1; i++)
            cr_s[m][i][1] = exp(-UR_S[m][i] + tr[m][i][0]) * EXP_BR[m][i] - tr[m][i][0] - 1;
      } else if (B1R2 == 1) {
        for (m = 0; m <= NRSITES - 1; m++)
          for (i = 0; i <= NNN - 1; i++)
            cr_s[m][i][1] = exp(-UR_S[m][i] + tr[m][i][0] + B1R[m][i]) * EXP_BR[m][i] - tr[m][i][0] - 1;
      } else {
        printf("\nError on bridge_func1\n");
        exit(1);
      }
    }
  } else
  /*PY*/
  if ((strncmp("py", CLOSURE, 2) == 0) && (strlen(CLOSURE) == 2)) {
    if (strncmp("no", RBC_FUNC, 2) == 0) {
      for (m = 0; m <= NRSITES - 1; m++)
        for (i = 0; i <= NNN - 1; i++)
          cr_s[m][i][1] = (exp(-UR_S[m][i]) - 1.00) * (1.00 + tr[m][i][0]);
    } else if ((strncmp("yes", RBC_FUNC, 3) == 0)) {
      for (m = 0; m <= NRSITES - 1; m++)
        for (i = 0; i <= NNN - 1; i++)
          cr_s[m][i][1] = exp(-UR_S[m][i]) * (1.00 + tr[m][i][0]) * EXP_BR[m][i] - (tr[m][i][0] + 1.00);
    }
  } else
  /*KH*/
  if ((strncmp("kh", CLOSURE, 2) == 0) && (strlen(CLOSURE) == 2)) {
    if (strncmp("no", RBC_FUNC, 2) == 0) {
      for (m = 0; m <= NRSITES - 1; m++)
        for (i = 0; i <= NNN - 1; i++) {
          wr = -UR_S[m][i] + tr[m][i][0];
          if (wr > 0.00)
            cr_s[m][i][1] = -UR_S[m][i];
          else
            cr_s[m][i][1] = exp(-UR_S[m][i] + tr[m][i][0]) - tr[m][i][0] - 1.0;
        }
    } else if ((strncmp("yes", RBC_FUNC, 3) == 0)) {
      for (m = 0; m <= NRSITES - 1; m++)
        for (i = 0; i <= NNN - 1; i++) {
          wr = exp(-UR_S[m][i] + tr[m][i][0]) * EXP_BR[m][i];
          if (wr > 1.00)
            cr_s[m][i][1] = -UR_S[m][i] + log(EXP_BR[m][i]);
          else
            cr_s[m][i][1] = wr - tr[m][i][0] - 1;
        }
    }
  } else {
    printf("\n\nNo closure specified\n"); fflush(stdout);
    exit(1);
  }

/*wall*/ if ((strncmp("wall1", CONFIG_TYPE, 5) == 0) && (strlen(CONFIG_TYPE) == 5))
    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= NNN - 1; i++) {
        if (HS[i] == 0.0)
          cr_s[m][i][1] = -1.0 - tr[m][i][0];
      }

/*wall*/ if ((strncmp("wall2", CONFIG_TYPE, 5) == 0) && (strlen(CONFIG_TYPE) == 5))
    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= NNN - 1; i++) {
        if (HSW[m][i] == 0.0)
          cr_s[m][i][1] = -1.0 - tr[m][i][0];
      }
}



/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/*												*/
/*					MDIIS							*/
/*												*/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/

double mul_vv(double *, double *, int);

int find_rmax(double *, int);

void mdiis_iter(int max_iter)
{
  /*EXTERN*/
  int nnn = NNN;
  double mp = DIIS_MP;
  int n_diis = DIIS_SIZE;
  /*END EXTERN*/

  int i, j, k, m, n, cnt, cnt2 = 0;
  int i_rmax = n_diis, io_rmax = n_diis;
  double any, test;

  double *r_mag = (double *) malloc((n_diis + 1) * sizeof(double));

  double ***cr = (double ***) malloc(NRSITES * sizeof(double));         /*[site][n_diis][i]*/
  double ***res = (double ***) malloc(NRSITES * sizeof(double));
  for (m = 0; m <= NRSITES - 1; m++) {
    *(cr + m) = (double **) malloc(n_diis * sizeof(double));
    *(res + m) = (double **) malloc((n_diis + 1) * sizeof(double));
  }

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= n_diis - 1; n++)
      cr[m][n] = (double *) malloc(nnn * sizeof(double));

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= n_diis; n++)
      res[m][n] = (double *) malloc(nnn * sizeof(double));

  double **s_mat, *c_vec, *b_vec;       /* sc-b=0 */
  c_vec = (double *) malloc((n_diis + 1) * sizeof(double));
  b_vec = (double *) malloc((n_diis + 1) * sizeof(double));
  s_mat = (double **) malloc((n_diis + 1) * sizeof(double));

  /*set s, c and b values*/
  for (i = 0; i <= n_diis; i++)
    *(s_mat + i) = (double *) malloc((n_diis + 1) * sizeof(double));

  for (i = 0; i <= n_diis; i++)
    b_vec[i] = 0.00;
  b_vec[n_diis] = -1.00;

  for (i = 0; i <= n_diis; i++) {
    s_mat[n_diis][i] = -1.00;
    s_mat[i][n_diis] = -1.00;
  }
  s_mat[n_diis][n_diis] = 0.00;


  do {
    if (i_rmax == io_rmax) {
      i_rmax = n_diis;
      io_rmax = n_diis;
      printf("\nInitial population or refill of cr and r vector \n"); fflush(stdout);
      /*initial population of cr and r*/
      for (j = 0; j <= n_diis - 1; j++) {
        for (m = 0; m <= NRSITES - 1; m++)
          for (i = 0; i <= nnn - 1; i++)
            cr[m][j][i] = CR_S[m][i];

        full_picard_iter(1);

        for (m = 0; m <= NRSITES - 1; m++)
          for (i = 0; i <= nnn - 1; i++)
            res[m][j][i] = CR2_S[m][i] - cr[m][j][i];

        r_mag[j] = svector_norm(res[0][j], nnn);

        for (m = 1; m <= NRSITES - 1; m++)
          r_mag[j] += svector_norm(res[m][j], nnn);
      }
      if (MY_RANK == 0)
        for (i = 0; i <= n_diis - 1; i++)
          printf("\nr_mag[%d] = %lf ", i, r_mag[i]); fflush(stdout);


      if (MY_RANK == 0)
        printf("\nCalculating initial s_mat \n"); fflush(stdout);
      /*s_mat***/
      for (j = 0; j <= n_diis - 1; j++)
        for (k = 0; k <= n_diis - 1; k++) {
          s_mat[j][k] = mul_vv(res[0][j], res[0][k], nnn);

          for (m = 1; m <= NRSITES - 1; m++)
            s_mat[j][k] += mul_vv(res[m][j], res[m][k], nnn);
        }
      if (MY_RANK == 0)
        printf("\n%d:Starting iterations\n", MY_RANK); fflush(stdout);
    }

    cnt = 0;
    MDIIS_CNT++;
    cnt++;
    cnt2++;
    CNT = cnt2;

    Newton_Solver(c_vec, s_mat, b_vec, n_diis + 1);

    if (MY_RANK == 0) {
      printf("\n\nc_vec -> "); fflush(stdout); any = 0;
      for (i = 0; i <= n_diis - 1; i++) {
        printf(" [%d]->%lf\t", i, c_vec[i]); any += c_vec[i];
      }
      ;
      printf("\ncoef total -> %lf\n\n", any); fflush(stdout);
    }

    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= nnn - 1; i++)
        res[m][n_diis][i] = 0.00;

    for (m = 0; m <= NRSITES - 1; m++)
      for (j = 0; j <= n_diis - 1; j++)
        for (i = 0; i <= nnn - 1; i++)
          res[m][n_diis][i] += (c_vec[j] * res[m][j][i]);

    r_mag[n_diis] = svector_norm(res[0][n_diis], nnn);
    for (m = 1; m <= NRSITES - 1; m++)
      r_mag[n_diis] += svector_norm(res[m][n_diis], nnn);

    if (MY_RANK == 0)
      printf("\nError[%d:%d] = %.10e\n", cnt2, cnt, r_mag[n_diis]);

    FERR = r_mag[n_diis];
    test = FERR;

    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= nnn - 1; i++)
        CR2_S[m][i] = 0.00;

    for (m = 0; m <= NRSITES - 1; m++)
      for (j = 0; j <= n_diis - 1; j++)
        for (i = 0; i <= nnn - 1; i++)
          CR2_S[m][i] += (c_vec[j] * cr[m][j][i]);

    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= nnn - 1; i++)
        CR_S[m][i] = CR2_S[m][i] + (mp * res[m][n_diis][i]);

    if (r_mag[n_diis] < T_ERR && RT_CNT >= RT_CHANGES) break;

    /*update diis vector with new c, r */
    io_rmax = i_rmax;
    i_rmax = find_rmax(r_mag, n_diis);

    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= nnn - 1; i++)
        cr[m][i_rmax][i] = CR_S[m][i];

    full_picard_iter(1);

    for (m = 0; m <= NRSITES - 1; m++)
      for (i = 0; i <= nnn - 1; i++)
        res[m][i_rmax][i] = CR2_S[m][i] - cr[m][i_rmax][i];

    r_mag[i_rmax] = svector_norm(res[0][i_rmax], nnn);
    for (m = 1; m <= NRSITES - 1; m++)
      r_mag[i_rmax] += svector_norm(res[m][i_rmax], nnn);

    /*update s_mat with new row and column i_rmax*/
    for (i = 0; i <= n_diis - 1; i++) {
      s_mat[i][i_rmax] = mul_vv(*(res[0] + i), *(res[0] + i_rmax), nnn);
      for (m = 1; m <= NRSITES - 1; m++)
        s_mat[i][i_rmax] += mul_vv(*(res[m] + i), *(res[m] + i_rmax), nnn);

      s_mat[i_rmax][i] = mul_vv(*(res[0] + i_rmax), *(res[0] + i), nnn);
      for (m = 1; m <= NRSITES - 1; m++)
        s_mat[i_rmax][i] += mul_vv(*(res[m] + i_rmax), *(res[m] + i), nnn);
    }

    /**RT_CHANGES and bridge functions **/
    if (test < RT_ERR && RT_CNT < RT_CHANGES) {
      if (strncmp("yes00", BRIDGE_FUNC1, 5) == 0)
        calc_bridge_func1(1);
      change_RT();
      i_rmax = io_rmax;                                 /*trick to refill diis vectors*/
      full_picard_iter(1);
      test = T_ERR + 10.0;                              /* ensure program doesn't prematurely exit */
      UPB1_STAT = 0;
      IERR = 9e9;
    } else if ((strncmp("yes", BRIDGE_FUNC1, 3)) == 0 && strncmp("mdiis", SOLVER, 5) == 0) {
      if (UPB1 >= 1)
        if (UPB1_STAT >= UPB1 && test < IERR) {
          calc_bridge_func1(1);
          IERR = test;
          UPB1_STAT = 0;
          i_rmax = io_rmax;                             /*trick to refill diis vectors*/
          full_picard_iter(1);
        }
    }
    UPB1_STAT++;
  } while ((test > T_ERR) && (cnt2 < max_iter));


  for (m = 0; m <= NRSITES - 1; m++)
    for (i = 0; i <= NNN - 1; i++)
      CR_S[m][i] = CR2_S[m][i];

  full_picard_iter(1);


  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= n_diis - 1; n++) free(cr[m][n]);

  for (m = 0; m <= NRSITES - 1; m++)
    for (n = 0; n <= n_diis; n++) free(res[m][n]);

  free(cr);
  free(res);
}



double mul_vv(double *v1, double *v2, int n)
{
  int i;
  double ret = 0;

  for (i = 0; i <= n - 1; i++) ret += v1[i] * v2[i];

  return ret;
}



int find_rmax(double *r_vec, int n)
{
  int i, k = 0;

  for (i = 1; i <= n - 1; i++) {
    if (r_vec[i] > r_vec[i - 1]) k = i;
  }

  return k;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				Utilities						*/
/*											*/
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


void fftw_3d(fftw_complex *in_r, fftw_complex *out_k)
{
  int nx = NX, ny = NY, nz = NZ, nnn = NNN;             /*EXTERN*/
  double lx = LX, ly = LY, lz = LZ;                     /*EXTERN*/

  int i;
  double dx = lx / (nx - 1);
  double dy = ly / (ny - 1);
  double dz = lz / (nz - 1);
  double con = dx * dy * dz;

#ifdef FFTW_THREADS
  fftw_init_threads();
  fftw_plan_with_nthreads(NUM_THREADS);
#endif

  fftw_plan rtok;
  rtok = fftw_plan_dft_3d(nx, ny, nz, in_r, out_k, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(rtok);
  for (i = 0; i <= nnn - 1; i++) {
    out_k[i][0] *= con;  out_k[i][1] *= con;
  }                                                                              /*Normalizing*/
}



void invfftw_3d(fftw_complex *in_k, fftw_complex *out_r)
{
  int nx = NX, ny = NY, nz = NZ, nnn = NNN;             /*EXTERN*/
  double lx = LX, ly = LY, lz = LZ;                     /*EXTERN*/

  int i;
  double dx = lx / (nx - 1),  dy = ly / (ny - 1),  dz = lz / (nz - 1);
  double con = dx * dy * dz * nnn;

#ifdef FFTW_THREADS
  fftw_init_threads();
  fftw_plan_with_nthreads(NUM_THREADS);
#endif

  fftw_plan ktor;
  ktor = fftw_plan_dft_3d(nx, ny, nz, in_k, out_r, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(ktor);
  for (i = 0; i <= nnn - 1; i++) {
    out_r[i][0] *= (1 / con);  out_r[i][1] *= (1 / con);
  }                                                                                  /*Normalizing*/
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				3D to 1D Averaging					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/


double * calc_3d_to_1d_avg(double *gr3d, double xx, double yy, double zz)
{
  int n_1d = 2 * NX, nx = NX, ny = NY, nz = NZ;         /*EXTERN*/
  double lx = LX, ly = LY, lz = LZ;                     /*EXTERN*/

  int i, x, y, z, id;
  double r, r1d, dr1d, lmax;

  int *grnum = (int *) malloc(n_1d * sizeof(int));
  double *gr1d = (double *) malloc(n_1d * sizeof(double));

  for (i = 0; i <= n_1d - 1; i++) gr1d[i] = 0.00;
  for (i = 0; i <= n_1d - 1; i++) grnum[i] = 0;

  lmax = lx;
  if (ly > lmax) lmax = ly;
  if (lz > lmax) lmax = lz;
  r1d = lmax;
  dr1d = (double) r1d / ((double) (n_1d - 1));

  for (x = 0; x <= nx - 1; x++)
    for (y = 0; y <= ny - 1; y++)
      for (z = 0; z < nz - 1; z++) {
        r = rx(x, y, z, xx, yy, zz);
        r += (dr1d / 2.0);
        id = (int) (r / dr1d);
        gr1d[id] += gr3d[ii(x,y,z)];
        grnum[id]++;
      }

  for (i = 0; i <= n_1d - 1; i++)
    if (grnum[i] == 0) gr1d[i] = 0.00;
    else gr1d[i] = gr1d[i] / ((double) grnum[i]);

  return gr1d;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*					r(x,y,z)					*/
/*				Calculate the spatial distances				*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

double rx(int x, int y, int z, double ux, double uy, double uz)
{
  double Nx = NX, Ny = NY, Nz = NZ, lx = LX, ly = LY, lz = LZ;  /*EXTERN*/
  int u_x = CX, u_y = CY, u_z = CZ;

  double tmp;
  double dx = lx / (Nx - 1);
  double dy = ly / (Ny - 1);
  double dz = lz / (Nz - 1);

  tmp = pow(fabs(x * dx - (u_x * dx + ux)), 2);
  tmp += pow(fabs(y * dy - (u_y * dy + uy)), 2);
  tmp += pow(fabs(z * dz - (u_z * dz + uz)), 2);
  tmp = sqrt(tmp);

  return tmp;
}



double r0(int x, int y, int z)
{
  double Nx = NX, Ny = NY, Nz = NZ, lx = LX, ly = LY, lz = LZ;  /*EXTERN*/
  int u_x = CX, u_y = CY, u_z = CZ;

  double tmp;
  double dx = lx / (Nx - 1);
  double dy = ly / (Ny - 1);
  double dz = lz / (Nz - 1);

  tmp = pow(fabs(x - u_x),2) * pow(dx,2);
  tmp += pow(fabs(y - u_y),2) * pow(dy,2);
  tmp += pow(fabs(z - u_z),2) * pow(dz,2);
  tmp = sqrt(tmp);

  return tmp;
}



double kx(int x, int y, int z, double ux, double uy, double uz)
{
  double nx = NX, ny = NY, nz = NZ, lx = LX, ly = LY, lz = LZ;  /*EXTERN*/
  int u_x = CX, u_y = CY, u_z = CZ;

  double tmp;
  double dkx = 2 * Pi / lx,  dky = 2 * Pi / ly,   dkz = 2 * Pi / lz;
  double dx = lx / (nx - 1), dy = ly / (ny - 1),  dz = lz / (nz - 1);

  tmp = pow(fabs(x - (u_x + (ux / dx))),2) * pow(dkx,2);
  tmp += pow(fabs(y - (u_y + (uy / dy))),2) * pow(dky,2);
  tmp += pow(fabs(z - (u_z + (uz / dz))),2) * pow(dkz,2);
  tmp = sqrt(tmp);
  return tmp;
}



double k0(int x, int y, int z)
{
  double lx = LX, ly = LY, lz = LZ;     /*EXTERN*/
  int u_x = CX, u_y = CY, u_z = CZ;

  double tmp;
  double dkx = 2 * Pi / lx;
  double dky = 2 * Pi / ly;
  double dkz = 2 * Pi / lz;

  tmp = pow((x - u_x) * dkx, 2);
  tmp += pow((y - u_y) * dky, 2);
  tmp += pow((z - u_z) * dkz, 2);
  tmp = sqrt(tmp);

  return tmp;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				Transformation routines					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void check_kr(U_PAR *u)
{
  int x,y;
  double kr, kx, ky, fx, fy;

  printf("\nCheck kr routines\n");
  printf("u[1].x=%f\tu[1].y=%f\tu[1].z=%f\n", u[1].x, u[1].y, u[1].z);

  for (x = 0; x <= 2; x++) {
    if (x == 0) fx = -1.0; else fx = 1.0;
    for (y = 0; y <= 2; y++) {
      if (y == 0) fy = -1.0;
      else fy = 1.0;
      kx = k0(x + (CX - 1), CY, CZ);
      ky = k0(CX, y + (CY - 1), CZ);
      kr = fx * kx * u[1].x + fy * ky * u[1].y + 0.00;

      printf("(kx,ky)-> %f\t%f\n", fx * kx, fy * ky);
      printf("\tkr->%f\n\n", kr); fflush(stdout);
    }
  }
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				Utilities						*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				get Routines						*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

double * get_3d(char fin[], double temp, double pnd, ENV_PAR sys)
{
  double *v3d = NULL;
  sprintf(fin, "%s.%s", fin, FILE_TYPE);

  if (strncmp("sit", FILE_TYPE, 3) == 0) {
    v3d = get_sit(fin, sys);
  } else if (strncmp("jh3d", FILE_TYPE, 4) == 0) {
    v3d = get_jh3d(fin, temp, pnd, sys);
  } else {
    fprintf(stderr, "\nError in get_3d with file type\n"); fflush(stdout);
  }

  return v3d;
}



fftw_complex * get_complex_3d(char fin[], double temp, double pnd, ENV_PAR sys)
{
  fftw_complex *v3d = NULL;
  sprintf(fin, "%s.%s", fin, FILE_TYPE);

  if (strncmp("sit", FILE_TYPE, 3) == 0) {
    v3d = get_complex_sit(fin, sys);
  } else if (strncmp("jh3d", FILE_TYPE, 4) == 0) {
    v3d = get_complex_jh3d(fin, temp, pnd, sys);
  } else {
    fprintf(stdout, "\nError in get_3d with file type\n"); fflush(stdout);
  }

  return v3d;
}



/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
/*											*/
/*				Printing Routines					*/
/*											*/
/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

void print_3d(char name[], double *v)
{
  sprintf(name, "%s.%s", name, FILE_TYPE);

  if (strncmp("sit", FILE_TYPE, 3) == 0)
    print_sit(name, v, SYS);
  else if (strncmp("jh3d", FILE_TYPE, 8) == 0)
    print_jh3d(name, v, SYS, TEMP, PND[0]);
  else
    printf("\nFile name not specified correctly\n"); fflush(stdout);
}



/**/ void print_jh3d_box(char name[], double *v, int nx, int ny, int nz)
{
  int x, y, z;
  double l[3] = {LX, LY, LZ};
  double PD = PND[0];
  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");


  fprintf(out, "%d\n%d\n%d\n", NX, NY, NZ / 2); fflush(out);
  fprintf(out, "%.10f\n%.10f\n%.10f\n", l[0], l[1], l[2] / 2.0); fflush(out);
  fprintf(out, "%.10f\n", TEMP); fflush(out);
  fprintf(out, "%.10f\n", PD); fflush(out);

  for (x = 0; x <= nx - 1; x++) {
    for (y = 0; y <= ny - 1; y++) {
      for (z = CZ; z <= NZ - 1; z++)
        fprintf(out, "%d\t%d\t%d\t%.15e\n", x, y, z - CZ, v[nz * ny * x + nz * y + z]);
      fprintf(out, "\n");
    }
    fflush(out);
  }
  fclose(out);
}



/**/ void print_jh3d_box3(char name[], double *v, int nx, int ny, int nz)
{
  int x, y, z;
  double l[3] = {LX, LY, LZ};
  double PD = PND[0];
  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");


  fprintf(out, "%d\n%d\n%d\n", NX, NY, NZ / 3); fflush(out);
  fprintf(out, "%.10f\n%.10f\n%.10f\n", l[0], l[1], l[2] / 3.0); fflush(out);
  fprintf(out, "%.10f\n", TEMP); fflush(out);
  fprintf(out, "%.10f\n", PD); fflush(out);

  for (x = 0; x <= nx - 1; x++) {
    for (y = 0; y <= ny - 1; y++) {
      for (z = CZ; z <= CZ + (NZ / 3) - 1; z++)
        fprintf(out, "%d\t%d\t%d\t%.15e\n", x, y, z - CZ, v[nz * ny * x + nz * y + z]);
      fprintf(out, "\n");
    }
    fflush(out);
  }
  fclose(out);
}



/****************************************************************************************/
/*											*/
/*				RT_CHANGES						*/
/*											*/
/****************************************************************************************/
void change_RT(void)
{
  int i, j;
  RT_CNT++;

  if (RT_CNT == RT_CHANGES) {
    printf("\n\n<<< Temp Factor Change %d of %d : %f -> %f >>> \n\n", RT_CNT, RT_CHANGES, RT_TEMP_FACTOR[RT_CNT - 1], TEMP_FACTOR);
    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++) {
        UR_S[j][i] *= RT_TEMP_FACTOR[RT_CNT - 1] / TEMP_FACTOR;
        if (strncmp("no", EWALD_SUMS, 2) == 0) {
          UR_L[j][i] *= RT_TEMP_FACTOR[RT_CNT - 1] / TEMP_FACTOR;
          UK_L[j][i][0] *= RT_TEMP_FACTOR[RT_CNT - 1] / TEMP_FACTOR;
          UK_L[j][i][1] *= RT_TEMP_FACTOR[RT_CNT - 1] / TEMP_FACTOR;
        }
      }
  } else {
    printf("\n\n<<< Temp Factor Change %d of %d : %f -> %f >>> \n\n", RT_CNT, RT_CHANGES, RT_TEMP_FACTOR[RT_CNT - 1], RT_TEMP_FACTOR[RT_CNT]);
    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++) {
        UR_S[j][i] *= RT_TEMP_FACTOR[RT_CNT - 1] / RT_TEMP_FACTOR[RT_CNT];
        if (strncmp("no", EWALD_SUMS, 2) == 0) {
          UR_L[j][i] *= RT_TEMP_FACTOR[RT_CNT - 1] / RT_TEMP_FACTOR[RT_CNT];
          UK_L[j][i][0] *= RT_TEMP_FACTOR[RT_CNT - 1] / RT_TEMP_FACTOR[RT_CNT];
          UK_L[j][i][1] *= RT_TEMP_FACTOR[RT_CNT - 1] / RT_TEMP_FACTOR[RT_CNT];
        }
      }
  }
}



/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************ OPTIONAL CODE *********************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/
/************************************************************************************************/


