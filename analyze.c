/****************************************************************************************/
/*											*/
/*			Version 0.1 by Jesse Howard					*/
/*			Date:   10-28-2008						*/
/*											*/
/****************************************************************************************/

/*This code calculates the thermodynamics of the input 3d gr file
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
#include "jh_util.h"
#include "jh_struct.h"
#include "jh_linalg.h"
#include "jh_print.h"

/*#define TEST*/

#define ii(x,y,z) (NZ * NY * (x) + NZ * (y) + (z))
#define NNN NX * NY * NZ

/*Global variables */
/*Environment variables read in, a change of these must be announced */
double DX, DY, DZ;

/*par properties*/
ENV_PAR SYS;
U_PAR2  *U, *U2;
U_PAR *U1;
int NU_SITES, PAR_TYPE;


/*env - values*/
int NX, NY, NZ;         /*set in get_jh3d*/
int CX, CY, CZ;
double LX, LY, LZ;              /*set in get_jh3d*/
double A_ERF, CHRG_PCT, TEMP_FACTOR;
char *EWALD_SUMS, *CLOSURE, *FILE_TYPE;
char *CONFIG_TYPE, *BRIDGE_FUNC1, *BRIDGE_FUNC0, *RBC_FUNC;

/*dis2 - solvent properties*/
int NSITES, NRSITES, TYPE;
double TEMP, *PND, DIS_NUM, *REDUN;
double *EP12, *EP6, *SIG, *CHARGE, *BOND_MAT;
char **NAMES, **DIS_NAMES;

/*file names*/
char **gr_fnames, **cr_s_fnames, **br_fnames, **b1r_fnames;
char **ur_lj_fnames, **ur_clmb_fnames, **ur_l_fnames, **ur_fnames, **uk_l_fnames;


/************************************************************************/
/*				Subroutines				*/
/************************************************************************/

/*1*/ void calc_internal_energy(void);
/*2*/ void calc_chem_potential(void);
/*3*/ void calc_coordination_number(void);
/*4*/ void calc_regional_internal_energy(void);
/*5*/ void calc_excluded_volume(void);
/*6*/ void calc_molar_volume(void);
/*7*/ void calc_solute_energy(void);
/*8*/ void calc_regional_chem_potential(void);
/*9*/ void calc_electrostatic_field(void);

void set_env(char []);
void set_dis(char []);
void set_par(char []);
void set_fnames(void);
void set_sys(void);

void get_solute_parameters(char *);
void check_parameters(void);
double * get_3d(char []);

double * get_ewald_bgc(void);

double r0(int, int, int);                       /*from origin  x, y, z*/
double rx(int, int, int, double, double, double);  /* x, y, z, ui_x, ui_y, ui_z*/
double rii(int, int, int, int, int, int);

double * calc_1d_to_3d_shift_fft(double *, int, double);
void print_1d(char [], double *, int, double);
void jh_print_2arrays(char [], double *, double *, int, double);
double k0(int, int, int);
double kx(int, int, int, double, double, double);

char *ARGV1, *ARGV2, *ARGV3;

/****************************************************************************************/
/*											*/
/*					Main						*/
/*											*/
/****************************************************************************************/

int main(int argc, char *argv[])
{
  int choice;

  ARGV1 = *(argv + 1);
  ARGV2 = *(argv + 2);
  ARGV3 = *(argv + 3);


  /*****************************GET SOLUTE PARAMETERS*****************************************/
  printf("\nReading in solute parameters...\n"); fflush(stdout);

  set_par(*(argv + 3));
  check_parameters();
  printf("...done!!!\n\n"); fflush(stdout);
  set_dis(*(argv + 1));
  set_env(*(argv + 2));
  set_sys();
  set_fnames();

  /*****************************KERNEL*********************************************/


  printf("\nSelect a choice:\n");
  printf("0: EXIT\n");
  printf("1: Calculate internal energy\n");
  printf("2: Calculate chemical potential\n");
  printf("3: Calculate coordination numbers N(r)\n");
  printf("4: Calculate regional internal energy\n");
  printf("5: Calculate excluded volume\n");
  printf("6: Calculate molar volume\n");
  printf("7: Calculate solute energy\n");
  printf("8: Calculate regional chemical potential\n");
  printf("9: Calculate electrostatic field (solute+solvent)\n");
  printf("\nENTER CHOICE:"); fflush(stdout);
  scanf("%d", &choice);
  printf("\nPerforming operation...\n"); fflush(stdout);

  switch (choice) {
  case 0: break;
  case 1: calc_internal_energy();                 break;
  case 2: calc_chem_potential();                  break;
  case 3: calc_coordination_number();             break;
  case 4: calc_regional_internal_energy();        break;
  case 5: calc_excluded_volume();                 break;
  case 6: calc_molar_volume();                    break;
  case 7: calc_solute_energy();                   break;
  case 8: calc_regional_chem_potential( );        break;
  case 9: calc_electrostatic_field( );            break;
  }

  printf("...DONE!!!\n\n"); fflush(stdout);


  return 0;
}



/****************************************************************************************/
/****************************************************************************************/
/**************************************END OF MAIN***************************************/
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
  CLOSURE = (char *)get_sval(infile, "CLOSURE");


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

  if (get_tag(infile, "EWALD_SUMS") == 1)
    EWALD_SUMS = (char *) get_sval(infile, "EWALD_SUMS");
  else EWALD_SUMS = "no";

  if (get_tag(infile, "CONFIG_TYPE") == 1)
    CONFIG_TYPE = (char *) get_sval(infile, "CONFIG_TYPE");
  else CONFIG_TYPE = "none";

  if (get_tag(infile, "FILE_TYPE") == 1)
    FILE_TYPE = (char *) get_sval(infile, "FILE_TYPE");
  else FILE_TYPE = "jh3d";

  if (get_tag(infile, "CHRG_PCT") == 1)
    CHRG_PCT = (double) get_dval(infile, "CHRG_PCT");
  else CHRG_PCT = 1.00;

  if (get_tag(infile, "TEMP_FACTOR") == 1)
    TEMP_FACTOR = (double) get_dval(infile, "TEMP_FACTOR");
  else TEMP_FACTOR = 1.00;

  if (get_tag(infile, "RBC_FUNC") == 1)
    RBC_FUNC = (char *) get_sval(infile, "RBC_FUNC");
  else RBC_FUNC = "no";

  if (get_tag(infile, "BRIDGE_FUNC0") == 1)
    BRIDGE_FUNC0 = (char *) get_sval(infile, "BRIDGE_FUNC0");
  else BRIDGE_FUNC0 = "no";

  if (get_tag(infile, "BRIDGE_FUNC1") == 1)
    BRIDGE_FUNC1 = (char *) get_sval(infile, "BRIDGE_FUNC1");
  else BRIDGE_FUNC1 = "no";
}



void set_dis(char infile[])
{
  TEMP = (double) get_dval(infile, "TEMP");
  TYPE = (int) get_dval(infile, "TYPE");
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
  DIS_NAMES = (char **) get_array_sval(infile, "DIS_NAMES", NRSITES);
}



void set_par(char infile[])
{
  int *nu_s = (int *) malloc(sizeof(int));
  int i, idx = 0, len;
  char vs;

  while (infile[idx] != '.') idx++;
  idx++;

  len = strlen(infile);
  vs = infile[len - 1];

  //if( strncmp( "par2", end , 4 ) == 0 ){
  if ('2' == vs) {
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
}



void set_sys(void)
{
  SYS.nx = NX;    SYS.ny = NY;    SYS.nz = NZ;
  SYS.cx = CX;    SYS.cy = CY;    SYS.cz = CZ;
  SYS.lx = LX;    SYS.ly = LY;    SYS.lz = LZ;

  DX = LX / (NX - 1);
  DY = LY / (NY - 1);
  DZ = LZ / (NZ - 1);
}



void set_fnames(void)
{
  int j;
  int fnlen = 100;

  cr_s_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  gr_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  ur_lj_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  ur_clmb_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  ur_l_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  ur_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  uk_l_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  br_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));
  b1r_fnames = (char **) malloc(NRSITES * 2 * sizeof(char));

  for (j = 0; j <= NRSITES - 1; j++) {
    *(cr_s_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(cr_s_fnames + j), "cr_%s_s.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(gr_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(gr_fnames + j), "gr_%s.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(ur_lj_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(ur_lj_fnames + j), "ur_%s_lj.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(ur_clmb_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(ur_clmb_fnames + j), "ur_%s_clmb.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(ur_l_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(ur_l_fnames + j), "ur_%s_l.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(ur_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(ur_fnames + j), "ur_%s.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(uk_l_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(uk_l_fnames + j), "uk_%s_l.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(br_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(br_fnames + j), "br_%s.%s", NAMES[j], FILE_TYPE);
  }

  for (j = 0; j <= NRSITES - 1; j++) {
    *(b1r_fnames + j) = (char *) malloc(fnlen * sizeof(char));
    sprintf(*(b1r_fnames + j), "b1r_%s.%s", NAMES[j], FILE_TYPE);
  }

#ifdef TEST
  for (i = 0; i <= NRSITES - 1; i++) {
    fprintf(stdout, "%s\n",    *(cr_s_fnames + i));
    fprintf(stdout, "%s\n",    *(gr_fnames + i));
    fprintf(stdout, "%s\n",    *(ur_lj_fnames + i));
    fprintf(stdout, "%s\n",    *(ur_clmb_fnames + i));
    fprintf(stdout, "%s\n",    *(ur_l_fnames + i));
    fprintf(stdout, "%s\n",    *(ur_fnames + i));
    fprintf(stdout, "%s\n",    *(uk_l_fnames + i));
    fprintf(stdout, "%s\n",    *(br_fnames + i));
    fprintf(stdout, "%s\n",    *(b1r_fnames + i));
    printf("\n");
  }
#endif
}



/****************************************************************************************/
/*											*/
/*				1: calc internal energy					*/
/*											*/
/****************************************************************************************/

void calc_internal_energy(void)
{
  int def, i, j, x, y, z;
  double jkm = CONST_JKM;
  double lj_energy = 0.00, clmb_energy = 0.00, tot_energy;
  double dr3d = DX * DY * DZ;
  char sdir[150], s1[150];

  double    **gr = (double **) malloc(NRSITES * sizeof(double));
  double **ur_lj = (double **) malloc(NRSITES * sizeof(double));
  double **ur_clmb = (double **) malloc(NRSITES * sizeof(double));
  double **ur_clmb2 = (double **) malloc(NRSITES * sizeof(double));
  double *int_lj = (double *) malloc(NNN * sizeof(double));
  double *int_clmb = (double *) malloc(NNN * sizeof(double));
  double *sum_lj = (double *) malloc(NRSITES * sizeof(double));
  double *sum_clmb = (double *) malloc(NRSITES * sizeof(double));
  double *gt_lj = (double *) malloc(NRSITES * sizeof(double));
  double *gt_clmb = (double *) malloc(NRSITES * sizeof(double));


  FILE *out;
  if ((out = fopen("internal_energy.dat", "w")) == NULL)
    printf("\nUnable to open internal energy file\n"); fflush(stdout);


  printf("\nReading gr files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(gr + i) = (double *) get_3d(*(gr_fnames + i));
  printf("\nReading ur_lj files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(ur_lj + i) = (double *) get_3d(*(ur_lj_fnames + i));
  printf("\nReading ur_clmb files\n"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(ur_clmb + i) = (double *) get_3d(*(ur_clmb_fnames + i));

  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    printf("\nUse ewald potential-0 or use real potential-1: ");
    scanf("%d", &def);
    printf("\nEnter directory to retrieve non-ewald potential: ");
    scanf("%s", sdir);
    printf("Reading new ur_clmb files"); fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1,"%s/%s", sdir, *(ur_clmb_fnames + i));
      *(ur_clmb2 + i) = (double *) get_3d(s1);
    }

    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        ur_clmb[j][i] = ur_clmb2[j][i];
  }


  for (j = 0; j <= NRSITES - 1; j++) {
    /*Filter gr_uo arrays in core region?????????*/
    for (i = 0; i <= NNN - 1; i++)
      if (gr[j][i] < 0.000001) gr[j][i] = 0.00;

    /*calculate integrand ur is in (KTs) */
    for (i = 0; i <= NNN - 1; i++) {
      int_lj[i] = ur_lj[j][i] * gr[j][i];
      int_clmb[i] = CHRG_PCT * ur_clmb[j][i] * gr[j][i];
    }

    sum_lj[j] = 0.00;
    sum_clmb[j] = 0.00;

    if (strncmp("wall1", CONFIG_TYPE, 5) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= (CZ + NZ / 3) - 1; z++) {
            sum_lj[j] += int_lj[ii(x,y,z)];
            sum_clmb[j] += int_clmb[ii(x,y,z)];
          }
    } else if (strncmp("wall", CONFIG_TYPE, 4) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= NZ - 1; z++) {
            sum_lj[j] += int_lj[ii(x,y,z)];
            sum_clmb[j] += int_clmb[ii(x,y,z)];
          }
    } else {
      for (i = 0; i <= NNN - 1; i++) {
        sum_lj[j] += int_lj[i];
        sum_clmb[j] += int_clmb[i];
      }
    }

    sum_lj[j] *= PND[j] * dr3d;
    sum_clmb[j] *= PND[j] * dr3d;

    gt_lj[j] = sum_lj[j] * jkm * TEMP * TEMP_FACTOR / 1000;
    gt_clmb[j] = sum_clmb[j] * jkm * TEMP * TEMP_FACTOR / 1000;

    lj_energy += REDUN[j] * gt_lj[j];
    clmb_energy += REDUN[j] * gt_clmb[j];
  }

  tot_energy = lj_energy + clmb_energy;


  printf("\nPrinting to internal_energy.dat\n"); fflush(stdout);
  fprintf(out, "*************************************************\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*             Internal Energy                  *\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*************************************************\n\n");

  fprintf(out, "Temperature (K) -> %f\n\n", TEMP * TEMP_FACTOR);

  fprintf(out, "____Lennard-Jones Energies____\n\n");

  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "%s(%lf)  =  %lf (kJ/mol)\n", NAMES[i], REDUN[i], gt_lj[i]);

  fprintf(out, "\n\tTotal LJ energy  =  %lf (kJ/mol)", lj_energy);
  fprintf(out, "\n\t                 =  %lf (kCal/mol)", lj_energy / 4.184);
  fprintf(out, "\n\t                 =  %lf (K)", lj_energy / 0.008315);


  fprintf(out, "\n\n____Coulomb Energies____\n\n");

  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "%s(%lf)  =  %lf (kJ/mol)\n", NAMES[i], REDUN[i], gt_clmb[i]);

  fprintf(out, "\n\tTotal coulomb energy = %lf (kJ/mol)", clmb_energy);
  fprintf(out, "\n\t                     = %lf (kCal/mol)", clmb_energy / 4.184);
  fprintf(out, "\n\t                     = %lf (K)\n\n", clmb_energy / 0.008315);


  fprintf(out, "\n____Total Energies____\n");
  fprintf(out, "\nTOTAL ENERGY = %lf (kJ/mol)", tot_energy);
  fprintf(out, "\n             = %lf (kCal/mol)", tot_energy / 4.184);
  fprintf(out, "\n             = %lf (K)\n\n", tot_energy / 0.008315);

  printf("\nTOTAL ENERGY = %lf (kJ/mol)", tot_energy);
  printf("\n             = %lf (kCal/mol)", tot_energy / 4.184);
  printf("\n             = %lf (K)\n\n", tot_energy / 0.008315);

  for (i = 0; i <= NRSITES - 1; i++) {
    free(*(gr + i));
    free(*(ur_lj + i));
    free(*(ur_clmb + i));
  }

  free(gr);
  free(ur_lj);   free(ur_clmb);
  free(int_lj);  free(int_clmb);
  free(sum_lj);  free(sum_clmb);
  free(gt_lj);   free(gt_clmb);
}



/****************************************************************************************/
/*											*/
/*				2: calc chemical potential				*/
/*											*/
/****************************************************************************************/

void calc_chem_potential(void)
{
  int i, j, x, y, z;

  int nnn = NX * NY * NZ;
  double jkm = CONST_JKM;
  double chem_pot = 0.00, b_chem_pot = 0.00;
  double dr3d = (LX / (NX - 1)) * (LY / (NY - 1)) * (LZ / (NZ - 1));


  /************************Read in g(r), c(r)_s, and u(r)_l terms**********************/
  char s1[150], sdir[150];

  double   **gr = (double **) malloc(NRSITES * sizeof(double));
  double **cr_s = (double **) malloc(NRSITES * sizeof(double));
  double **ur_l = (double **) malloc(NRSITES * sizeof(double));
  double   **br = (double **) malloc(NRSITES * sizeof(double));
  double   **b1r = (double **) malloc(NRSITES * sizeof(double));

  double **ur_clmb1 = (double **) malloc(NRSITES * sizeof(double));
  double **ur_clmb2 = (double **) malloc(NRSITES * sizeof(double));


  printf("\nReading gr files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(gr + i) = (double *) get_3d(*(gr_fnames + i));
  printf("\nReading cr_s  files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(cr_s + i) = (double *) get_3d(*(cr_s_fnames + i));
  printf("\nReading ur_l files"); fflush(stdout);
  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (i = 0; i <= NRSITES - 1; i++)
      *(ur_l + i) = (double *) get_3d(*(ur_l_fnames + i));

  if (strncmp("yes", RBC_FUNC, 3) == 0) {
    printf("\nReading rbc files");  fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++)
      *(br + i) = (double *) get_3d(*(br_fnames + i));
  }

  if ((strncmp("yes", BRIDGE_FUNC1, 3) == 0) || (strncmp("yes", BRIDGE_FUNC0, 3) == 0)) {
    printf("\nReading br1 files");  fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++)
      *(b1r + i) = (double *) get_3d(*(b1r_fnames + i));
  }

  printf("\n...Done!!!\n");       fflush(stdout);

  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        cr_s[j][i] = cr_s[j][i] - (CHRG_PCT * ur_l[j][i]);


  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    for (i = 0; i <= NRSITES - 1; i++)
      *(ur_clmb1 + i) = (double *) get_3d(*(ur_clmb_fnames + i));

    printf("\nEnter directory to retrieve non-ewald (real) potentials: ");
    scanf("%s", sdir);
    printf("Reading new ur_clmb files"); fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1,"%s/%s", sdir, *(ur_clmb_fnames + i));
      *(ur_clmb2 + i) = (double *) get_3d(s1);
    }

    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        cr_s[j][i] = cr_s[j][i] + CHRG_PCT * (ur_clmb1[j][i] - ur_clmb2[j][i]);
  }

  /***************************************************************************************/

  double *ux_gt = (double *) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++)
    ux_gt[i] = 0.00;

  double **ux_int = (double **) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++)
    *(ux_int + i) = (double *) malloc(NNN * sizeof(double));

  /*Filter gr_uo arrays in core region*/
  /*for( i=0; i<=nnn-1; i++)  if( gr_uo[i] < 0.00000001 )  gr_uo[i] = 0.00;*/
  /*for( i=0; i<=nnn-1; i++)  if( gr_uh[i] < 0.00000001 )  gr_uh[i] = 0.00;*/


  /*calculate integrand */
  for (j = 0; j <= NRSITES - 1; j++) {
    if ((strncmp("hnc", CLOSURE, 3) == 0) && (strlen(CLOSURE) == 3)) {
      for (i = 0; i <= nnn - 1; i++) {
        ux_int[j][i] = (0.5) * (gr[j][i] - 1.0) * (gr[j][i] - 1.0);                                                     /* h*h*/
        ux_int[j][i] += (-0.5) * (cr_s[j][i]) * (gr[j][i] - 1.0);                               /* c(r) * h(r) */
        ux_int[j][i] += (-1.0) * (cr_s[j][i]);
      }
      if (strncmp("yes", BRIDGE_FUNC1, 3) == 0 || strncmp("yes", BRIDGE_FUNC0,3) == 0) {
        for (i = 0; i <= nnn - 1; i++) {
          ux_int[j][i] += (2.0 / 3.0) * (b1r[j][i]) * (gr[j][i] - 1.0);                                 /* c(r) * h(r) */
          ux_int[j][i] += (1.0) * (b1r[j][i]);
        }
      }

      if ((strncmp("yes", RBC_FUNC, 3) == 0))
        for (i = 0; i <= nnn - 1; i++)
          ux_int[j][i] += gr[j][i] * (br[j][i] - 1.0);                                                  /* h*h*/
    } else if ((strncmp("kh", CLOSURE, 2) == 0) && (strlen(CLOSURE) == 2)) {
      for (i = 0; i <= nnn - 1; i++) {
        if ((-1.0 * (gr[j][i] - 1.0) >= 0))
          ux_int[j][i] = (0.5) * (gr[j][i] - 1.0) * (gr[j][i] - 1.0);                                                   /* h*h*/
        else
          ux_int[j][i] = 0.00;

        ux_int[j][i] += (-0.5) * (cr_s[j][i]) * (gr[j][i] - 1.0);                               /* c(r) * h(r) */
        ux_int[j][i] += (-1.0) * (cr_s[j][i]);
      }
    } else if ((strncmp("py", CLOSURE, 3) == 0) && (strlen(CLOSURE) == 2)) {
      for (i = 0; i <= nnn - 1; i++) {
        ux_int[j][i] = (cr_s[j][i]) / ((gr[j][i] - 1.0) - (cr_s[j][i]));                                 /* c(r) / (h(r)-c(r)) */
        ux_int[j][i] *= (-1.0) * log(gr[j][i] - (cr_s[j][i]));                                  /* -ln( 1+h(r)-c(r) )  */
      }
    } else
      printf("\nNO closure specified\n"); fflush(stdout);
  }
  /*calc integrals */

  for (j = 0; j <= NRSITES - 1; j++) {
    if (strncmp("wall1", CONFIG_TYPE, 5) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= (CZ + NZ / 3) - 1; z++) {
            ux_gt[j] += ux_int[j][ii(x,y,z)];
          }
    } else if (strncmp("wall", CONFIG_TYPE, 4) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= NZ - 1; z++) {
            ux_gt[j] += ux_int[j][ii(x,y,z)];
          }
    } else {
      for (i = 0; i <= nnn - 1; i++)
        ux_gt[j] += ux_int[j][i];
    }
  }

  /*calc chemical potential as beta chem_pot*/
  for (j = 0; j <= NRSITES - 1; j++)
    b_chem_pot += (REDUN[j] * PND[j] * ux_gt[j]);

  b_chem_pot *= dr3d;
  chem_pot = jkm * TEMP * TEMP_FACTOR * b_chem_pot / 1000;

  FILE *out;
  if ((out = fopen("chemical_potential.dat", "w")) == NULL) {
    printf("\nFile chemical_potential.dat could not be opened"); fflush(stdout);
  } else printf("\nWriting to chemical_potential.dat\n\n"); fflush(stdout);

  fprintf(out, "*************************************************\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*             Chemical Potential	       *\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*************************************************\n\n");


  fprintf(out, "Temperature (K) -> %f\n", TEMP);
  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "\nParticle density(%s) -> %lf", NAMES[i], PND[i]);

  fprintf(out, "\n\n");

  fprintf(out, "chemical potential = %.6f (kJ/mol)\n", chem_pot);
  fprintf(out, "                   = %.6f (kCal/mol)\n", chem_pot / 4.184);
  fprintf(out, "                   = %.6f (k)\n\n", chem_pot / 0.008315);

  fprintf(out, "*************************************************\n\n");

  printf("chemical potential = %.6f (kJ/mol)\n", chem_pot);
  printf("                   = %.6f (kCal/mol)\n", chem_pot / 4.184);
  printf("                   = %.6f (k)\n\n", chem_pot / 0.008315);

  for (i = 0; i <= NRSITES - 1; i++) {
    free(*(gr + i));
    free(*(cr_s + i));
    free(*(ur_l + i));
    free(*(ux_int + i));
  }

  free(gr);
  free(cr_s);
  free(ur_l);
  free(ux_int);
  free(ux_gt);
}



/****************************************************************************************/
/*											*/
/*				3: calculate coordination numbers			*/
/*											*/
/****************************************************************************************/

void calc_coordination_number(void)
{
  int il = 0, i, int_meth, sid, x = 0, y = 0, z = 0;
  double rad, r_xyz;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double xorg, yorg, zorg;
  double r_x, r_y, r_z;
  double coor_num = 0.0;
  int tmp_num = 0;
  double dr3d = DX * DY * DZ;

  printf("\n\nChoose number corresponding to interaction site:\n");
  for (i = 0; i <= NRSITES - 1; i++)
    printf("%d: %s\n", i + 1, *(NAMES + i));

  printf("\nENTER CHOICE: ");
  scanf("%d", &sid); printf("\n");
  sid--;
  printf("Reading file \"%s\" ...", *(gr_fnames + sid));         fflush(stdout);
  double *gr = (double *) get_3d(*(gr_fnames + sid));         printf("...Done!!!\n");

  printf("\nWhat type of Integral, (1: Cartesian)  (2: Spherical)\nEnter number #: ");    scanf("%d", &int_meth);

  if (int_meth == 1) {
    printf("\nUse \"%s\" for integral limits (1: yes) or (2: no) : ", ARGV2); fflush(stdout);
    scanf("%d", &il);
    printf("\n");

    if (il == 1) {
      xmin = -1.0 * CX * DX;      xmax = (NX - 1 - CX) * DX;
      ymin = -1.0 * CY * DY;      ymax = (NY - 1 - CY) * DY;
      zmin = -1.0 * CZ * DZ;      zmax = (NZ - 1 - CZ) * DZ;

      for (i = 0; i <= NNN - 1; i++) {
        coor_num += gr[ii(x,y,z)];
        tmp_num++;
      }
    } else if (il == 2) {
      printf("\n");
      printf("\nEnter units in Angstroms: \n");
      printf("Enter x min -> "); scanf("%lf", &xmin);
      printf("Enter x max -> "); scanf("%lf", &xmax);
      printf("Enter y min -> "); scanf("%lf", &ymin);
      printf("Enter y max -> "); scanf("%lf", &ymax);
      printf("Enter z min -> "); scanf("%lf", &zmin);
      printf("Enter z max -> "); scanf("%lf", &zmax);

      for (x = 0; x <= NX - 1; x++) {
        r_x = r0(x,CY,CZ);
        if (x < CX) r_x = -r_x;

        if (r_x >= xmin && r_x <= xmax)
          for (y = 0; y <= NY - 1; y++) {
            r_y = r0(CX, y, CZ);
            if (y < CY) r_y = -r_y;

            if (r_y >= ymin && r_y <= ymax)
              for (z = 0; z <= NZ - 1; z++) {
                r_z = r0(CX, CY, z);
                if (z < CZ) r_z = -r_z;

                if (r_z >= zmin && r_z <= zmax) {
                  coor_num += gr[ii(x,y,z)];
                  tmp_num++;
                }
              }
          }
      }
    } else {
      printf("\nERROR with integration limits\n\n");
    }
  } else if (int_meth == 2) {
    printf("\nEnter x origin: ");           scanf("%lf", &xorg);
    printf("Enter y origin: ");             scanf("%lf", &yorg);
    printf("Enter z origin: ");             scanf("%lf", &zorg);
    printf("Enter radius in Angstroms: ");  scanf("%lf", &rad);

    for (x = 0; x <= NX - 1; x++)
      for (y = 0; y <= NY - 1; y++)
        for (z = 0; z <= NZ - 1; z++) {
          r_xyz = rx(x,y,z, xorg, yorg, zorg);
          if (r_xyz <= rad) {
            coor_num += gr[ii(x,y,z)];
            tmp_num++;
          }
        }
  } else {
    printf("No integration method specified, exiting\n"); fflush(stdout);
  }

  coor_num *= (REDUN[sid] * PND[sid] * dr3d);

  FILE *out;
  char f_name[100];
  sprintf(f_name,"coordinated_%s_sites.dat", NAMES[sid]);
  if ((out = fopen(f_name, "w")) == NULL)
    printf("\nError opening file %s", f_name);

  printf("\nNumber of points -> %d \n", tmp_num);
  printf("\nprinting to \"%s\"\n", f_name);

  fprintf(out, "*************************************************\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*             Coordination number	       *\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*************************************************\n\n");

  fprintf(out, "Integral limits (Angstroms) : \n\n");
  if (int_meth == 1) {
    fprintf(out, "\tBox with the following limits:\n");
    fprintf(out, "\t\tx-min -> %f \t x-max -> %f \n", xmin, xmax);                        fflush(stdout);
    fprintf(out, "\t\ty-min -> %f \t y-max -> %f \n", ymin, ymax);                        fflush(stdout);
    fprintf(out, "\t\tz-min -> %f \t z-max -> %f \n\n", zmin, zmax);              fflush(stdout);
  } else if (int_meth == 2) {
    fprintf(out,"\tSphere center at:\n");
    fprintf(out,"\t\tx-origin -> %f\n", xorg);
    fprintf(out,"\t\ty-origin -> %f\n", yorg);
    fprintf(out,"\t\tz-origin -> %f\n", zorg);
    fprintf(out,"\t\twith radius -> %f\n\n", rad);
  }
  fprintf(out, "Number of points -> %d \n", tmp_num);
  fprintf(out, "Number density -> %f\n\n", PND[sid]);                                   fflush(stdout);
  fprintf(out, "Coordinated Number of \"%s\" sites -> %f\n", *(NAMES + sid), coor_num);        fflush(stdout);
  fprintf(out, "\n*************************************************\n\n");

  printf("\n\nCoordinated Number of \"%s\" sites -> %f\n\n", *(NAMES + sid), coor_num);        fflush(stdout);
}



/****************************************************************************************/
/*											*/
/*			4: calc regional internal energy				*/
/*											*/
/****************************************************************************************/


void calc_regional_internal_energy(void)
{
  int def, i, il, int_meth, j, tmp_num = 0, x, y, z;
  double jkm = CONST_JKM;
  double rad, r_xyz;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double xorg, yorg, zorg;
  double r_x, r_y, r_z;
  double lj_energy = 0.00, clmb_energy = 0.00, tot_energy;
  double dr3d = DX * DY * DZ;
  char sdir[150], s1[150];

  double    **gr = (double **) malloc(NRSITES * sizeof(double));
  double **ur_lj = (double **) malloc(NRSITES * sizeof(double));
  double **ur_clmb = (double **) malloc(NRSITES * sizeof(double));
  double **ur_clmb2 = (double **) malloc(NRSITES * sizeof(double));
  double *int_lj = (double *) malloc(NNN * sizeof(double));
  double *int_clmb = (double *) malloc(NNN * sizeof(double));
  double *sum_lj = (double *) malloc(NRSITES * sizeof(double));
  double *sum_clmb = (double *) malloc(NRSITES * sizeof(double));
  double *gt_lj = (double *) malloc(NRSITES * sizeof(double));
  double *gt_clmb = (double *) malloc(NRSITES * sizeof(double));
  double *yorn = (double *) malloc(NNN * sizeof(double));

  for (i = 0; i <= NNN - 1; i++)
    yorn[i] = 0.0;

  printf("\nReading gr files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(gr + i) = (double *) get_3d(*(gr_fnames + i));
  printf("\nReading ur_lj files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(ur_lj + i) = (double *) get_3d(*(ur_lj_fnames + i));
  printf("\nReading ur_clmb files\n"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(ur_clmb + i) = (double *) get_3d(*(ur_clmb_fnames + i));

  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    printf("\nUse ewald potential-0 or use real potential-1: ");
    scanf("%d", &def);
    printf("\nEnter directory to retrieve non-ewald potential: ");
    scanf("%s", sdir);
    printf("Reading new ur_clmb files"); fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1,"%s/%s", sdir, *(ur_clmb_fnames + i));
      *(ur_clmb2 + i) = (double *) get_3d(s1);
    }

    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        ur_clmb[j][i] = ur_clmb2[j][i];
  }


  printf("\nWhat type of Integral, (1: Cartesian)  (2: Spherical)\nEnter number #: ");    scanf("%d", &int_meth);


  if (int_meth == 1) {
    printf("\nUse \"%s\" for integral limits (1: yes) or (2: no) : ", ARGV2); fflush(stdout);
    scanf("%d", &il);
    printf("\n");

    if (il == 1) {
      xmin = -1.0 * CX * DX;      xmax = (NX - 1 - CX) * DX;
      ymin = -1.0 * CY * DY;      ymax = (NY - 1 - CY) * DY;
      zmin = -1.0 * CZ * DZ;      zmax = (NZ - 1 - CZ) * DZ;

      for (i = 0; i <= NNN - 1; i++) {
        yorn[i] = 1.0;
        tmp_num++;
      }
    } else if (il == 2) {
      printf("\n");
      printf("\nEnter units in Angstroms: \n");
      printf("Enter x min -> "); scanf("%lf", &xmin);
      printf("Enter x max -> "); scanf("%lf", &xmax);
      printf("Enter y min -> "); scanf("%lf", &ymin);
      printf("Enter y max -> "); scanf("%lf", &ymax);
      printf("Enter z min -> "); scanf("%lf", &zmin);
      printf("Enter z max -> "); scanf("%lf", &zmax);

      for (x = 0; x <= NX - 1; x++) {
        r_x = r0(x,CY,CZ);
        if (x < CX) r_x = -r_x;

        if (r_x >= xmin && r_x <= xmax)
          for (y = 0; y <= NY - 1; y++) {
            r_y = r0(CX, y, CZ);
            if (y < CY) r_y = -r_y;

            if (r_y >= ymin && r_y <= ymax)
              for (z = 0; z <= NZ - 1; z++) {
                r_z = r0(CX, CY, z);
                if (z < CZ) r_z = -r_z;

                if (r_z >= zmin && r_z <= zmax) {
                  yorn[ii(x,y,z)] = 1.0;
                  tmp_num++;
                }
              }
          }
      }
    }
  } else if (int_meth == 2) {
    printf("\nEnter x origin: ");           scanf("%lf", &xorg);
    printf("Enter y origin: ");             scanf("%lf", &yorg);
    printf("Enter z origin: ");             scanf("%lf", &zorg);
    printf("Enter radius in Angstroms: ");  scanf("%lf", &rad);

    for (x = 0; x <= NX - 1; x++)
      for (y = 0; y <= NY - 1; y++)
        for (z = 0; z <= NZ - 1; z++) {
          r_xyz = rx(x,y,z, xorg, yorg, zorg);
          if (r_xyz <= rad) {
            yorn[ii(x,y,z)] = 1.0;
            tmp_num++;
          }
        }
  } else {
    printf("No integration method specified, exiting\n"); fflush(stdout);
  }


  for (j = 0; j <= NRSITES - 1; j++) {
    /*Filter gr_uo arrays in core region?????????*/
    for (i = 0; i <= NNN - 1; i++)
      if (gr[j][i] < 0.000001) gr[j][i] = 0.00;

    /*calculate integrand ur is in (KTs) */
    for (i = 0; i <= NNN - 1; i++) {
      int_lj[i] = ur_lj[j][i] * gr[j][i] * yorn[i];
      int_clmb[i] = CHRG_PCT * ur_clmb[j][i] * gr[j][i] * yorn[i];
    }

    sum_lj[j] = 0.00;
    sum_clmb[j] = 0.00;

    if (strncmp("wall1", CONFIG_TYPE, 5) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= (CZ + NZ / 3) - 1; z++) {
            sum_lj[j] += int_lj[ii(x,y,z)];
            sum_clmb[j] += int_clmb[ii(x,y,z)];
          }
    } else if (strncmp("wall", CONFIG_TYPE, 4) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= NZ - 1; z++) {
            sum_lj[j] += int_lj[ii(x,y,z)];
            sum_clmb[j] += int_clmb[ii(x,y,z)];
          }
    } else {
      for (i = 0; i <= NNN - 1; i++) {
        sum_lj[j] += int_lj[i];
        sum_clmb[j] += int_clmb[i];
      }
    }

    sum_lj[j] *= PND[j] * dr3d;
    sum_clmb[j] *= PND[j] * dr3d;

    gt_lj[j] = sum_lj[j] * jkm * TEMP * TEMP_FACTOR / 1000;
    gt_clmb[j] = sum_clmb[j] * jkm * TEMP * TEMP_FACTOR / 1000;

    lj_energy += REDUN[j] * gt_lj[j];
    clmb_energy += REDUN[j] * gt_clmb[j];
  }

  tot_energy = lj_energy + clmb_energy;

  FILE *out;
  char f_name[100];
  sprintf(f_name,"regional_internal_energy.dat");
  if ((out = fopen(f_name, "w")) == NULL)
    printf("\nError opening file %s", f_name);


  printf("\nNumber of points -> %d \n", tmp_num);
  printf("\nPrinting to regional_internal_energy.dat\n"); fflush(stdout);
  fprintf(out, "*************************************************\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*            Regional Internal Energy	       *\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*************************************************\n\n");

  fprintf(out, "Temperature (K) -> %f\n\n", TEMP * TEMP_FACTOR);

  fprintf(out, "Integral limits (Angstroms) : \n\n");
  if (int_meth == 1) {
    fprintf(out, "\tBox with the following limits:\n");
    fprintf(out, "\t\tx-min -> %f \t x-max -> %f \n", xmin, xmax);                        fflush(stdout);
    fprintf(out, "\t\ty-min -> %f \t y-max -> %f \n", ymin, ymax);                        fflush(stdout);
    fprintf(out, "\t\tz-min -> %f \t z-max -> %f \n\n", zmin, zmax);              fflush(stdout);
  } else if (int_meth == 2) {
    fprintf(out,"\tSphere center at:\n");
    fprintf(out,"\t\tx-origin -> %f\n", xorg);
    fprintf(out,"\t\ty-origin -> %f\n", yorg);
    fprintf(out,"\t\tz-origin -> %f\n", zorg);
    fprintf(out,"\t\twith radius -> %f\n\n", rad);
  }
  fprintf(out, "\tNumber of points -> %d \n\n\n", tmp_num);

  fprintf(out, "____Lennard-Jones Energies____\n\n");

  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "%s(%lf)  =  %lf (kJ/mol)\n", NAMES[i], REDUN[i], gt_lj[i]);

  fprintf(out, "\n\tTotal LJ energy  =  %lf (kJ/mol)", lj_energy);
  fprintf(out, "\n\t                 =  %lf (kCal/mol)", lj_energy / 4.184);
  fprintf(out, "\n\t                 =  %lf (K)", lj_energy / 0.008315);


  fprintf(out, "\n\n____Coulomb Energies____\n\n");

  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "%s(%lf)  =  %lf (kJ/mol)\n", NAMES[i], REDUN[i], gt_clmb[i]);

  fprintf(out, "\n\tTotal coulomb energy = %lf (kJ/mol)", clmb_energy);
  fprintf(out, "\n\t                     = %lf (kCal/mol)", clmb_energy / 4.184);
  fprintf(out, "\n\t                     = %lf (K)\n\n", clmb_energy / 0.008315);


  fprintf(out, "\n____Total Energies____\n");
  fprintf(out, "\nTOTAL ENERGY = %lf (kJ/mol)", tot_energy);
  fprintf(out, "\n             = %lf (kCal/mol)", tot_energy / 4.184);
  fprintf(out, "\n             = %lf (K)\n\n", tot_energy / 0.008315);

  printf("\nTOTAL ENERGY = %lf (kJ/mol)", tot_energy);
  printf("\n             = %lf (kCal/mol)", tot_energy / 4.184);
  printf("\n             = %lf (K)\n\n", tot_energy / 0.008315);

  for (i = 0; i <= NRSITES - 1; i++) {
    free(*(gr + i));
    free(*(ur_lj + i));
    free(*(ur_clmb + i));
  }

  free(gr);
  free(ur_lj);   free(ur_clmb);
  free(int_lj);  free(int_clmb);
  free(sum_lj);  free(sum_clmb);
  free(gt_lj);   free(gt_clmb);
}



/****************************************************************************************/
/*											*/
/*			5: excluded volume						*/
/*											*/
/****************************************************************************************/


void calc_excluded_volume(void)
{
  int i, j;
  double crit;
  double dr3d = (LX / (NX - 1)) * (LY / (NY - 1)) * (LZ / (NZ - 1));
  double *sum = (double *) malloc(NRSITES * sizeof(double));
  double *ex_vol = (double *) malloc(NRSITES * sizeof(double));
  double **gr = (double **) malloc(NRSITES * sizeof(double));

  FILE *out;
  if ((out = fopen("excluded_volume.dat", "w")) == NULL)
    printf("\nUnable to open excluded volume file\n"); fflush(stdout);

  printf("\nReading gr files..."); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(gr + i) = (double *) get_3d(*(gr_fnames + i));
  printf("DONE!\n");

  printf("\nEnter g(r) max value to count as excluded volume (<=~0.001): "); scanf("%lf", &crit);

  for (i = 0; i <= NRSITES - 1; i++) {
    sum[i] = 0.0;
    ex_vol[i] = 0.0;
    for (j = 0; j <= NNN - 1; j++)
      if (gr[i][j] < crit)
        sum[i] += 1.0;

    ex_vol[i] = sum[i] * dr3d;
  }


  fprintf(out, "*********************************\n");
  fprintf(out, "*                               *\n");
  fprintf(out, "*      Excluded volume (A^3)    *\n");
  fprintf(out, "*                               *\n");
  fprintf(out, "*********************************\n");

  fprintf(out, "\nvoxel volume : %lf", dr3d);
  fprintf(out, "\ncritical g(r) value: %e\n\n", crit);

  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "Sum of empty %s voxels = %lf\n", NAMES[i], sum[i]);
  fprintf(out, "\n");

  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "Excluded volume of %s sites (A^3) -> %lf\n", NAMES[i], ex_vol[i]); fflush(stdout);
  fprintf(out, "\n");

  printf("\ncritical g(r) value: %e\n\n", crit);
  for (i = 0; i <= NRSITES - 1; i++)
    printf("Excluded volume of %s sites (A^3) -> %lf\n", NAMES[i], ex_vol[i]); fflush(stdout);
  printf("\n");
}



/****************************************************************************************/
/*											*/
/*			6: Molar volume							*/
/*											*/
/****************************************************************************************/

void calc_molar_volume(void)
{
  int i;
  double *gr_3d;
  char gr_3d_name[50];

  FILE *out;
  if ((out = fopen("molar_volume.dat", "w")) == NULL)
    printf("\nUnable to open molar volume file\n"); fflush(stdout);

  printf("\nEnter gr.jh3d file name: ");  scanf("%s", gr_3d_name); printf("\n");

  printf("\n\nReading in %s", gr_3d_name);
  gr_3d = (double *) get_3d(gr_3d_name);       printf("...DONE!!!\n");                fflush(stdout);

  double pnd = 1;       /* TODO: determine the pnd, missing from the original code */
  double lx = LX, ly = LY, lz = LZ;
  double nv = 0, nt, m_v;
  double vol = LX * LY * LZ;
  double dr3d = DX * DY * DZ;


  for (i = 0; i <= NNN - 1; i++)
    nv += gr_3d[i] * dr3d;

  nv = nv * pnd;
  nt = pnd * vol;

  m_v = (nt - nv) / pnd;

  fprintf(out, "*********************************\n");
  fprintf(out, "*	molar volume (L/mol)   *\n");
  fprintf(out, "*********************************\n");
  fprintf(out, "\nparameters -> (%lf, %lf, %lf, %lf)\n", pnd, lx, ly, lz); fflush(stdout);
  fprintf(out, "\nbox volume (A^3)           -> %lf\n", vol); fflush(stdout);
  fprintf(out, "\n# bulk water in volume     -> %lf", nt); fflush(stdout);
  fprintf(out, "\n# hydrated water in volume -> %lf\n", nv); fflush(stdout);
  fprintf(out, "\nchange in Volume (A^3)         -> %lf\n", m_v); fflush(stdout);

  m_v = m_v * 1e-30 * 6.02e23 * 1000;

  fprintf(out,    "\nPartial molar Volume (L/mol) -> %lf\n", m_v); fflush(out);
  fprintf(stdout, "\nPartial molar Volume (L/mol) -> %lf\n", m_v); fflush(stdout);
}



/****************************************************************************************/
/*											*/
/*			7: Solute Energy						*/
/*											*/
/****************************************************************************************/


void calc_solute_energy(void)
{
  int im, jm, is, js;
  int nmols = 0;
  double ep12, ep6, sig;
  double rx, ry, rz, r;
  double t = TEMP, K = CONST_CLMB, diel = 1.0;
  double jkm = CONST_JKM;
  double lj = 0.00, clmb = 0.00, tot = 0.00;
  FILE *out;

  if ((out = fopen("solute_energy.dat", "w")) == NULL)
    printf("\nError opening solute_energy.dat\n"); fflush(stdout);

  printf("\nEnter screening constant: ");         scanf("%lf", &diel); printf("\n");

  for (im = 1; im <= NU_SITES; im++)
    if (U[im].mol > nmols)
      nmols = U[im].mol;

  for (im = 1; im <= nmols - 1; im++)
    for (jm = im + 1; jm <= nmols; jm++) {
      for (is = 1; is <= NU_SITES; is++)
        for (js = 1; js <= NU_SITES; js++) {
          if ((U[is].mol == im) && (U[js].mol == jm)) {
            ep12 = sqrt(U[is].ep12 * U[js].ep12);
            ep6 = sqrt(U[is].ep6 * U[js].ep6);
            sig = (U[is].sig + U[js].sig) / 2.0;
            rx = fabs(U[is].x - U[js].x);
            ry = fabs(U[is].y - U[js].y);
            rz = fabs(U[is].z - U[js].z);
            if (strncmp("wall", U[is].element, 4) == 0 || strncmp("wall", U[js].element, 4) == 0)
              r = rz;
            else
              r = sqrt(rx * rx + ry * ry + rz * rz);

            lj += 4.0 / t * ((ep12 * pow((sig / r),12)) - (ep6 * pow((sig / r),6)));
            clmb += K / t * (U[is].charge * U[js].charge) / r;
          }
        }
    }

  lj = lj * jkm * TEMP * TEMP_FACTOR / 1000;
  clmb = (clmb / diel) * jkm * TEMP * TEMP_FACTOR / 1000;
  tot = lj + clmb;

  printf("\nprinting to \"solute_energy.dat\"\n");

  fprintf(out, "*************************************************\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*       Solute interaction Energy (kJ/mol)     *\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*************************************************\n\n");

  fprintf(out, "____Lennard-Jones Energies____\n");
  fprintf(out, "\n\tTotal LJ energy  =  %lf (kJ/mol)", lj);
  fprintf(out, "\n\t                 =  %lf (kCal/mol)", lj / 4.184);
  fprintf(out, "\n\t                 =  %lf (K)", lj / 0.008315);


  fprintf(out, "\n\n____Coulomb Energies____\n");
  fprintf(out, "\n\tTotal coulomb energy = %lf (kJ/mol)", clmb);
  fprintf(out, "\n\t                     = %lf (kCal/mol)", clmb / 4.184);
  fprintf(out, "\n\t                     = %lf (K)\n", clmb / 0.008315);


  fprintf(out, "\nTOTAL ENERGY = %lf (kJ/mol)", tot);
  fprintf(out, "\n             = %lf (kCal/mol)", tot / 4.184);
  fprintf(out, "\n             = %lf (K)\n\n", tot / 0.008315);

  printf("\nTOTAL ENERGY = %lf (kJ/mol)", tot);
  printf("\n             = %lf (kCal/mol)", tot / 4.184);
  printf("\n             = %lf (K)\n\n", tot / 0.008315);
}



/****************************************************************************************/
/*											*/
/*			8: Regional Chemical Potential					*/
/*											*/
/****************************************************************************************/


void calc_regional_chem_potential(void)
{
  int i, il, int_meth, j, x, y, z, tmp_num = 0;

  int nnn = NX * NY * NZ;
  double jkm = CONST_JKM;
  double chem_pot = 0.00, b_chem_pot = 0.00;
  double rad, r_xyz;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double xorg, yorg, zorg;
  double r_x, r_y, r_z;
  double dr3d = DX * DY * DZ;
  double *yorn = (double *) malloc(NNN * sizeof(double));
  for (i = 0; i <= NNN - 1; i++)
    yorn[i] = 0.0;


  /************************Read in g(r), c(r)_s, and u(r)_l terms**********************/
  char s1[150], sdir[150];

  double   **gr = (double **) malloc(NRSITES * sizeof(double));
  double **cr_s = (double **) malloc(NRSITES * sizeof(double));
  double **ur_l = (double **) malloc(NRSITES * sizeof(double));
  double   **br = (double **) malloc(NRSITES * sizeof(double));
  double   **b1r = (double **) malloc(NRSITES * sizeof(double));

  double **ur_clmb1 = (double **) malloc(NRSITES * sizeof(double));
  double **ur_clmb2 = (double **) malloc(NRSITES * sizeof(double));


  printf("\nWhat type of Integral, (1: Cartesian)  (2: Spherical)\nEnter number #: ");    scanf("%d", &int_meth);


  if (int_meth == 1) {
    printf("\nUse \"%s\" for integral limits (1: yes) or (2: no) : ", ARGV2); fflush(stdout);
    scanf("%d", &il);
    printf("\n");

    if (il == 1) {
      xmin = -1.0 * CX * DX;      xmax = (NX - 1 - CX) * DX;
      ymin = -1.0 * CY * DY;      ymax = (NY - 1 - CY) * DY;
      zmin = -1.0 * CZ * DZ;      zmax = (NZ - 1 - CZ) * DZ;

      for (i = 0; i <= NNN - 1; i++) {
        yorn[i] = 1.0;
        tmp_num++;
      }
    } else if (il == 2) {
      printf("\n");
      printf("\nEnter units in Angstroms: \n");
      printf("Enter x min -> "); scanf("%lf", &xmin);
      printf("Enter x max -> "); scanf("%lf", &xmax);
      printf("Enter y min -> "); scanf("%lf", &ymin);
      printf("Enter y max -> "); scanf("%lf", &ymax);
      printf("Enter z min -> "); scanf("%lf", &zmin);
      printf("Enter z max -> "); scanf("%lf", &zmax);

      for (x = 0; x <= NX - 1; x++) {
        r_x = r0(x,CY,CZ);
        if (x < CX) r_x = -r_x;

        if (r_x >= xmin && r_x <= xmax)
          for (y = 0; y <= NY - 1; y++) {
            r_y = r0(CX, y, CZ);
            if (y < CY) r_y = -r_y;

            if (r_y >= ymin && r_y <= ymax)
              for (z = 0; z <= NZ - 1; z++) {
                r_z = r0(CX, CY, z);
                if (z < CZ) r_z = -r_z;

                if (r_z >= zmin && r_z <= zmax) {
                  yorn[ii(x,y,z)] = 1.0;
                  tmp_num++;
                }
              }
          }
      }
    }
  } else if (int_meth == 2) {
    printf("\nEnter x origin: ");           scanf("%lf", &xorg);
    printf("Enter y origin: ");             scanf("%lf", &yorg);
    printf("Enter z origin: ");             scanf("%lf", &zorg);
    printf("Enter radius in Angstroms: ");  scanf("%lf", &rad);

    for (x = 0; x <= NX - 1; x++)
      for (y = 0; y <= NY - 1; y++)
        for (z = 0; z <= NZ - 1; z++) {
          r_xyz = rx(x,y,z, xorg, yorg, zorg);
          if (r_xyz <= rad) {
            yorn[ii(x,y,z)] = 1.0;
            tmp_num++;
          }
        }
  } else {
    printf("No integration method specified, exiting\n"); fflush(stdout);
  }


  printf("\nReading gr files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(gr + i) = (double *) get_3d(*(gr_fnames + i));
  printf("\nReading cr_s  files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(cr_s + i) = (double *) get_3d(*(cr_s_fnames + i));
  printf("\nReading ur_l files"); fflush(stdout);
  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (i = 0; i <= NRSITES - 1; i++)
      *(ur_l + i) = (double *) get_3d(*(ur_l_fnames + i));

  if (strncmp("yes", RBC_FUNC, 3) == 0) {
    printf("\nReading rbc files");  fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++)
      *(br + i) = (double *) get_3d(*(br_fnames + i));
  }

  if ((strncmp("yes", BRIDGE_FUNC1, 3) == 0) || (strncmp("yes", BRIDGE_FUNC0, 3) == 0)) {
    printf("\nReading br1 files");  fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++)
      *(b1r + i) = (double *) get_3d(*(b1r_fnames + i));
  }

  printf("\n...Done!!!\n");       fflush(stdout);

  if (strncmp("no", EWALD_SUMS, 2) == 0)
    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        cr_s[j][i] = cr_s[j][i] - (CHRG_PCT * ur_l[j][i]);


  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    for (i = 0; i <= NRSITES - 1; i++)
      *(ur_clmb1 + i) = (double *) get_3d(*(ur_clmb_fnames + i));

    printf("\nEnter directory to retrieve non-ewald (real) potentials: ");
    scanf("%s", sdir);
    printf("Reading new ur_clmb files"); fflush(stdout);
    for (i = 0; i <= NRSITES - 1; i++) {
      sprintf(s1,"%s/%s", sdir, *(ur_clmb_fnames + i));
      *(ur_clmb2 + i) = (double *) get_3d(s1);
    }

    for (j = 0; j <= NRSITES - 1; j++)
      for (i = 0; i <= NNN - 1; i++)
        cr_s[j][i] = cr_s[j][i] + CHRG_PCT * (ur_clmb1[j][i] - ur_clmb2[j][i]);
  }

  /***************************************************************************************/

  double *ux_gt = (double *) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++)
    ux_gt[i] = 0.00;

  double **ux_int = (double **) malloc(NRSITES * sizeof(double));
  for (i = 0; i <= NRSITES - 1; i++)
    *(ux_int + i) = (double *) malloc(NNN * sizeof(double));

  /*Filter gr_uo arrays in core region*/
  /*for( i=0; i<=nnn-1; i++)  if( gr_uo[i] < 0.00000001 )  gr_uo[i] = 0.00;*/
  /*for( i=0; i<=nnn-1; i++)  if( gr_uh[i] < 0.00000001 )  gr_uh[i] = 0.00;*/


  /*calculate integrand */
  for (j = 0; j <= NRSITES - 1; j++) {
    if ((strncmp("hnc", CLOSURE, 3) == 0) && (strlen(CLOSURE) == 3)) {
      for (i = 0; i <= nnn - 1; i++) {
        ux_int[j][i] = (0.5) * (gr[j][i] - 1.0) * (gr[j][i] - 1.0);                                                     /* h*h*/
        ux_int[j][i] += (-0.5) * (cr_s[j][i]) * (gr[j][i] - 1.0);                               /* c(r) * h(r) */
        ux_int[j][i] += (-1.0) * (cr_s[j][i]);
      }
      if (strncmp("yes", BRIDGE_FUNC1, 3) == 0 || strncmp("yes", BRIDGE_FUNC0,3) == 0) {
        for (i = 0; i <= nnn - 1; i++) {
          ux_int[j][i] += (2.0 / 3.0) * (b1r[j][i]) * (gr[j][i] - 1.0);                                 /* c(r) * h(r) */
          ux_int[j][i] += (1.0) * (b1r[j][i]);
        }
      }

      if ((strncmp("yes", RBC_FUNC, 3) == 0))
        for (i = 0; i <= nnn - 1; i++)
          ux_int[j][i] += gr[j][i] * (br[j][i] - 1.0);                                                  /* h*h*/
    } else if ((strncmp("kh", CLOSURE, 2) == 0) && (strlen(CLOSURE) == 2)) {
      for (i = 0; i <= nnn - 1; i++) {
        if ((-1.0 * (gr[j][i] - 1.0) >= 0))
          ux_int[j][i] = (0.5) * (gr[j][i] - 1.0) * (gr[j][i] - 1.0);                                                   /* h*h*/
        else
          ux_int[j][i] = 0.00;

        ux_int[j][i] += (-0.5) * (cr_s[j][i]) * (gr[j][i] - 1.0);                               /* c(r) * h(r) */
        ux_int[j][i] += (-1.0) * (cr_s[j][i]);
      }
    } else if ((strncmp("py", CLOSURE, 3) == 0) && (strlen(CLOSURE) == 2)) {
      for (i = 0; i <= nnn - 1; i++) {
        ux_int[j][i] = (cr_s[j][i]) / ((gr[j][i] - 1.0) - (cr_s[j][i]));                                 /* c(r) / (h(r)-c(r)) */
        ux_int[j][i] *= (-1.0) * log(gr[j][i] - (cr_s[j][i]));                                  /* -ln( 1+h(r)-c(r) )  */
      }
    } else
      printf("\nNO closure specified\n"); fflush(stdout);
  }
  /*calc integrals */

  for (j = 0; j <= NRSITES - 1; j++) {
    if (strncmp("wall1", CONFIG_TYPE, 5) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= (CZ + NZ / 3) - 1; z++) {
            ux_gt[j] += ux_int[j][ii(x,y,z)] * yorn[ii(x,y,z)];
          }
    } else if (strncmp("wall", CONFIG_TYPE, 4) == 0) {
      for (x = 0; x <= NX - 1; x++)
        for (y = 0; y <= NY - 1; y++)
          for (z = CZ; z <= NZ - 1; z++) {
            ux_gt[j] += ux_int[j][ii(x,y,z)] * yorn[ii(x,y,z)];
          }
    } else {
      for (i = 0; i <= nnn - 1; i++)
        ux_gt[j] += ux_int[j][i] * yorn[i];
    }
  }

  /*calc chemical potential as beta chem_pot*/
  for (j = 0; j <= NRSITES - 1; j++)
    b_chem_pot += (REDUN[j] * PND[j] * ux_gt[j]);

  b_chem_pot *= dr3d;
  chem_pot = jkm * TEMP * TEMP_FACTOR * b_chem_pot / 1000;

  FILE *out;
  if ((out = fopen("regional_chemical_potential.dat", "w")) == NULL) {
    printf("\nFile regional_chemical_potential.dat could not be opened"); fflush(stdout);
  } else printf("\nWriting to regional_chemical_potential.dat\n\n"); fflush(stdout);

  fprintf(out, "*************************************************\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*          Regional Chemical Potential	       *\n");
  fprintf(out, "*                                              *\n");
  fprintf(out, "*************************************************\n\n");


  fprintf(out, "Temperature (K) -> %f\n\n", TEMP);

  fprintf(out, "Integral limits (Angstroms) : \n\n");
  if (int_meth == 1) {
    fprintf(out, "\tBox with the following limits:\n");
    fprintf(out, "\t\tx-min -> %f \t x-max -> %f \n", xmin, xmax);                        fflush(stdout);
    fprintf(out, "\t\ty-min -> %f \t y-max -> %f \n", ymin, ymax);                        fflush(stdout);
    fprintf(out, "\t\tz-min -> %f \t z-max -> %f \n\n", zmin, zmax);              fflush(stdout);
  } else if (int_meth == 2) {
    fprintf(out,"\tSphere center at:\n");
    fprintf(out,"\t\tx-origin -> %f\n", xorg);
    fprintf(out,"\t\ty-origin -> %f\n", yorg);
    fprintf(out,"\t\tz-origin -> %f\n", zorg);
    fprintf(out,"\t\twith radius -> %f\n\n", rad);
  }
  fprintf(out, "\tNumber of points -> %d \n", tmp_num);
  for (i = 0; i <= NRSITES - 1; i++)
    fprintf(out, "\nParticle density(%s) -> %lf", NAMES[i], PND[i]);

  fprintf(out, "\n\n");

  fprintf(out, "chemical potential = %.6f (kJ/mol)\n", chem_pot);
  fprintf(out, "                   = %.6f (kCal/mol)\n", chem_pot / 4.184);
  fprintf(out, "                   = %.6f (k)\n\n", chem_pot / 0.008315);

  fprintf(out, "*************************************************\n\n");

  printf("chemical potential = %.6f (kJ/mol)\n", chem_pot);
  printf("                   = %.6f (kCal/mol)\n", chem_pot / 4.184);
  printf("                   = %.6f (k)\n\n", chem_pot / 0.008315);

  for (i = 0; i <= NRSITES - 1; i++) {
    free(*(gr + i));
    free(*(cr_s + i));
    free(*(ur_l + i));
    free(*(ux_int + i));
  }

  free(gr);
  free(cr_s);
  free(ur_l);
  free(ux_int);
  free(ux_gt);
}



/****************************************************************************************/
/*											*/
/*				9: calc electrostatic field				*/
/*											*/
/****************************************************************************************/

void calc_electrostatic_field(void)
{
  int def, i, x, y, z, x1, y1, z1;
  int idx, idx1;
  double r_ii, die;
  double fx = -1.0 * CX * DX;
  double fy = -1.0 * CY * DY;
  double fz = -1.0 * CZ * DZ;

  char sdir[150], s1[150];

  double    **gr = (double **) malloc(NRSITES * sizeof(double));

  double *ur_clmb_u;
  double *ur_clmb2_u;
  double *ur_clmb_v = (double *) malloc(NNN * sizeof(double));
  for (i = 0; i <= NNN - 1; i++)
    ur_clmb_v[i] = 0.00;

  printf("Enter dielectric constant of system: ");
  scanf("%lf", &die);
  printf("Enter units of ouput: (1) KT");
  printf("\nYou have chosen option-1 since it is the only choice\n");

  /* only need one */ printf("\nReading ur_clmb files\n"); fflush(stdout);
  ur_clmb_u = (double *) get_3d(*(ur_clmb_fnames + 0));
  printf("\nReading gr files"); fflush(stdout);
  for (i = 0; i <= NRSITES - 1; i++)
    *(gr + i) = (double *) get_3d(*(gr_fnames + i));

  if (strncmp("yes", EWALD_SUMS, 3) == 0) {
    printf("\nUse ewald potential-0 or use real potential-1: ");
    scanf("%d", &def);
    printf("\nEnter directory to retrieve non-ewald potential: ");
    scanf("%s", sdir);
    printf("Reading new ur_clmb files"); fflush(stdout);
    sprintf(s1,"%s/%s", sdir, *(ur_clmb_fnames + 0));
    ur_clmb2_u = (double *) get_3d(s1);

    for (i = 0; i <= NNN - 1; i++)
      ur_clmb_u[i] = ur_clmb2_u[i];
  }

  /* change potential energy to electrostatic energy */
  for (i = 0; i <= NNN - 1; i++)
    ur_clmb_u[i] /= (die * CHARGE[0] * TEMP_FACTOR);

  printf("\nCalculating electric field"); fflush(stdout);
  for (x = 0; x <= NX - 1; x++) {
    printf("."); fflush(stdout);
    for (y = 0; y <= NY - 1; y++)
      for (z = 0; z <= NZ - 1; z++) {
        idx = ii(x, y, z);

        for (x1 = 0; x1 <= NX - 1; x1++)
          for (y1 = 0; y1 <= NY - 1; y1++)
            for (z1 = 0; z1 <= NZ - 1; z1++) {
              idx1 = ii(x1, y1, z1);
              r_ii = rii(x,y,z,x1,y1,z1);

              for (i = 0; i <= NRSITES - 1; i++) {
                ur_clmb_v[idx] = (CHARGE[i] * REDUN[i] * PND[i] * gr[i][idx1]) / r_ii;
              }
            }
      }
  }
  printf("Done!\n"); fflush(stdout);

  for (i = 0; i <= NNN - 1; i++)
    ur_clmb_v[i] *= CONST_CLMB / (die * TEMP * TEMP_FACTOR);                            /*This make it units of KT */


  FILE *ufile, *vfile, *uvfile;
  if ((ufile = fopen("solute_electric_field.sit", "w")) == NULL)
    printf("\nUnable to open internal energy file\n"); fflush(stdout);
  if ((vfile = fopen("solvent_electric_field.sit", "w")) == NULL)
    printf("\nUnable to open internal energy file\n"); fflush(stdout);
  if ((uvfile = fopen("total_electric_field.sit", "w")) == NULL)
    printf("\nUnable to open internal energy file\n"); fflush(stdout);

  printf("\nPrinting solute electric field (KT) to \"solute_electric_field.sit\"\n");
  fprintf(ufile,"%f\n", DX);
  fprintf(ufile,"%f  %f  %f\n", fx, fy, fz);
  fprintf(ufile,"%d  %d  %d\n", NX, NY, NZ);
  for (z = 0; z <= NZ - 1; z++)
    for (y = 0; y <= NY - 1; y++)
      for (x = 0; x <= NX - 1; x++)
        fprintf(ufile, "%e\n", ur_clmb_u[ii(x,y,z)]);

  printf("\nPrinting solvent electric field (KT) to \"solvent_electric_field.sit\"\n");
  fprintf(vfile,"%f\n", DX);
  fprintf(vfile,"%f  %f  %f\n", fx, fy, fz);
  fprintf(vfile,"%d  %d  %d\n", NX, NY, NZ);
  for (z = 0; z <= NZ - 1; z++)
    for (y = 0; y <= NY - 1; y++)
      for (x = 0; x <= NX - 1; x++)
        fprintf(vfile, "%e\n", ur_clmb_v[ii(x,y,z)]);

  printf("\nPrinting total electric field (KT) to \"solvent_electric_field.sit\"\n");
  fprintf(uvfile,"%f\n", DX);
  fprintf(uvfile,"%f  %f  %f\n", fx, fy, fz);
  fprintf(uvfile,"%d  %d  %d\n", NX, NY, NZ);
  for (z = 0; z <= NZ - 1; z++)
    for (y = 0; y <= NY - 1; y++)
      for (x = 0; x <= NX - 1; x++)
        fprintf(uvfile, "%e\n", ur_clmb_u[ii(x,y,z)] + ur_clmb_v[ii(x,y,z)]);
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
/****************************************************************************************/


/****************************************************************************************/
/*											*/
/*			Calculate EWALD Correction Constant				*/
/*											*/
/****************************************************************************************/


double * get_ewald_bgc(void)
{
  int i;
  double *h0 = (double *) malloc(NRSITES * sizeof(double));

  //double *h0 = (double *) get_array_dval( "ewald.dat", "EWALD_BGC", NRSITES );

  for (i = 0; i <= NRSITES - 1; i++)
    h0[i] = 0.00;

  return h0;
}



void check_parameters(void)
{
  int i;
  FILE *out;
  if ((out = fopen("solute.analyze.check", "w")) == NULL)
    fprintf(stdout, "Problem opening out file for lj parameters"); fflush(stdout);
  U_PAR2 *u = U;
  int nu_sites = NU_SITES;

  for (i = 1; i <= nu_sites; i++) {
    fprintf(out, "\n%d:%s\n", u[i].num, u[i].element); fflush(out);
    fprintf(out, "ep12:%f\tep6:%f\tsig:%f\n", u[i].ep12, u[i].ep6, u[i].sig); fflush(out);
    fprintf(out, "Coulomb Charge:%f\n", u[i].charge); fflush(out);
    fprintf(out, "Cartesian Coordinates:\n"); fflush(out);
    fprintf(out, "%f\t%f\t%f\n\n", u[i].x, u[i].y, u[i].z); fflush(out);
  }
  fclose(out);
}



double * get_3d(char fname[])
{
  double *v3d = NULL;

  if (strncmp("sit", FILE_TYPE, 3) == 0) {
    v3d = get_sit(fname, SYS);
  } else if (strncmp("jh3d", FILE_TYPE, 4) == 0) {
    v3d = get_jh3d(fname, TEMP, 0.0, SYS);
  } else {
    fprintf(stdout, "\nError in get_3d with file type\n"); fflush(stdout);
  }

  return v3d;
}



/****************************************************************************************/
/*					r(x,y,z)					*/
/*				Calculate the spatial distances				*/
/*											*/
/****************************************************************************************/

double rx(int x, int y, int z, double ux, double uy, double uz)
{
  double tmp;

  tmp = pow(fabs(x * DX - (CX * DX + ux)), 2);
  tmp += pow(fabs(y * DY - (CY * DY + uy)), 2);
  tmp += pow(fabs(z * DZ - (CZ * DZ + uz)), 2);
  tmp = sqrt(tmp);

  return tmp;
}



double r0(int x, int y, int z)
{
  double tmp;

  tmp = pow(fabs(x - CX),2) * pow(DX,2);
  tmp += pow(fabs(y - CY),2) * pow(DY,2);
  tmp += pow(fabs(z - CZ),2) * pow(DZ,2);
  tmp = sqrt(tmp);

  return tmp;
}



double rii(int x, int y, int z, int x1, int y1, int z1)
{
  double tmp;
  tmp = pow(fabs(x - x1),2) * pow(DX,2);
  tmp += pow(fabs(y - y1),2) * pow(DY,2);
  tmp += pow(fabs(z - z1),2) * pow(DZ,2);
  tmp = sqrt(tmp);

  return tmp;
}



double kx(int x, int y, int z, double ux, double uy, double uz)
{
  double tmp;
  double dkx = 2 * Pi / LX;
  double dky = 2 * Pi / LY;
  double dkz = 2 * Pi / LZ;

  tmp = pow(fabs(x - (CX + (ux / DX))),2) * pow(dkx,2);
  tmp += pow(fabs(y - (CY + (uy / DY))),2) * pow(dky,2);
  tmp += pow(fabs(z - (CZ + (uz / DZ))),2) * pow(dkz,2);
  tmp = sqrt(tmp);
  return tmp;
}



double k0(int x, int y, int z)
{
  double tmp;
  double dkx = 2 * Pi / LX;
  double dky = 2 * Pi / LY;
  double dkz = 2 * Pi / LZ;

  tmp = pow(fabs(x - CX),2) * pow(dkx,2);
  tmp += pow(fabs(y - CY),2) * pow(dky,2);
  tmp += pow(fabs(z - CZ),2) * pow(dkz,2);
  tmp = sqrt(tmp);

  return tmp;
}



