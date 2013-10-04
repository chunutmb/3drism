#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "jh_get.h"
#include "jh_struct.h"

int NX, NY, NZ;
int CX, CY, CZ;
double LX, LY, LZ;
int MAX_ITER, N_DUMP, INIT_PIC_ITER, DIIS_SIZE;
double T_ERR, A_ERF, CHRG_PCT, PIC_MP, DIIS_MP;
double GAUSS_SMEAR;
char *CONTINUE, *SOLVER, *CLOSURE, *EWALD_SUMS;
char *CONFIG_TYPE;

void set_env(char []);
void print_env(char []);

int main(int argc, char *argv[])
{
  set_env(*(argv + 1));
  print_env(*(argv + 1));

  return 0;
}



void set_env(char infile[])
{
  NX = (int) get_dval(infile, "NX");
  NY = (int) get_dval(infile, "NY");
  NZ = (int) get_dval(infile, "NZ");
  LX = (double) get_dval(infile, "LX");
  LY = (double) get_dval(infile, "LY");
  LZ = (double) get_dval(infile, "LZ");
  CX = (int) get_dval(infile, "CX");
  CY = (int) get_dval(infile, "CY");
  CZ = (int) get_dval(infile, "CZ");

  A_ERF = (double) get_dval(infile, "A_ERF");
  CHRG_PCT = (double) get_dval(infile, "CHRG_PCT");
  CLOSURE = (char *)get_sval(infile, "CLOSURE");
  CONFIG_TYPE = (char *) get_sval(infile, "CONFIG_TYPE");
  CONTINUE = (char *)get_sval(infile, "CONTINUE");
  DIIS_SIZE = (int) get_dval(infile,  "DIIS_SIZE");
  DIIS_MP = (double) get_dval(infile, "DIIS_MP");
  EWALD_SUMS = (char *) get_sval(infile, "EWALD_SUMS");
  GAUSS_SMEAR = (double) get_dval(infile, "GAUSS_SMEAR");
  INIT_PIC_ITER = (int) get_dval(infile, "INIT_PIC_ITER");
  MAX_ITER = (int) get_dval(infile, "MAX_ITER");
  N_DUMP = (int) get_dval(infile, "N_DUMP");
  PIC_MP = (double) get_dval(infile, "PIC_MP");
  SOLVER = (char *)get_sval(infile, "SOLVER");
  T_ERR = (double) get_dval(infile, "T_ERR");
}



void print_env(char outfile[])
{
  int d;
  double f;
  char s1[50], s2[50];
  FILE *out;
  if ((out = fopen(outfile, "w")) == NULL)
    printf("Can't open %s to write\n", outfile);

  printf("\nEnter variable to change: ");
  scanf("%s", s1);
  printf("\nNew value: ");

  if (strncmp("NX", s1, 2) == 0) {
    scanf("%d", &d);
    fprintf(out, "#NX\t\t%d\n",  d);
  } else fprintf(out, "#NX\t\t%d\n", NX);

  if (strncmp("NY", s1, 2) == 0) {
    scanf("%d", &d);
    fprintf(out, "#NY\t\t%d\n",  d);
  } else fprintf(out, "#NY\t\t%d\n", NY);

  if (strncmp("NZ", s1, 2) == 0) {
    scanf("%d", &d);
    fprintf(out, "#NZ\t\t%d\n",  d);
  } else fprintf(out, "#NZ\t\t%d\n", NZ);

  if (strncmp("LX", s1, 2) == 0) {
    scanf("%lf", &f);
    fprintf(out, "#LX\t\t%f\n",  f);
  } else fprintf(out, "#LX\t\t%f\n", LX);

  if (strncmp("LY", s1, 2) == 0) {
    scanf("%lf", &f);
    fprintf(out, "#LY\t\t%f\n",  f);
  } else fprintf(out, "#LY\t\t%f\n", LY);

  if (strncmp("LZ", s1, 2) == 0) {
    scanf("%lf", &f);
    fprintf(out, "#LZ\t\t%f\n", f);
  } else fprintf(out, "#LZ\t\t%f\n", LZ);

  if (strncmp("CX", s1, 2) == 0) {
    scanf("%d", &d);
    fprintf(out, "#CX\t\t%d\n", d);
  } else fprintf(out, "#CX\t\t%d\n", CX);

  if (strncmp("CY", s1, 2) == 0) {
    scanf("%d", &d);    fprintf(out, "#CY\t\t%d\n", d);
  } else fprintf(out, "#CY\t\t%d\n", CY);

  if (strncmp("CZ", s1, 2) == 0) {
    scanf("%d", &d);    fprintf(out, "#CZ\t\t%d\n", d);
  } else fprintf(out, "#CZ\t\t%d\n", CZ);

  if (strncmp("A_ERF", s1, 5) == 0) {
    scanf("%lf", &f);  fprintf(out, "#A_ERF        \t%f\n",  f);
  } else fprintf(out, "#A_ERF        \t%f\n", A_ERF);

  if (strncmp("CHRG_PCT", s1, 8) == 0) {
    scanf("%lf", &f);  fprintf(out, "#CHRG_PCT     \t%f\n",  f);
  } else fprintf(out, "#CHRG_PCT     \t%f\n", CHRG_PCT);

  if (strncmp("CLOSURE", s1, 7) == 0) {
    scanf("%s", s2);  fprintf(out, "#CLOSURE      \t%s\n", s2);
  } else fprintf(out, "#CLOSURE      \t%s\n", CLOSURE);

  if (strncmp("CONFIG_TYPE", s1,11) == 0) {
    scanf("%s", s2);  fprintf(out, "#CONFIG_TYPE  \t%s\n", s2);
  } else fprintf(out, "#CONFIG_TYPE  \t%s\n", CONFIG_TYPE);

  if (strncmp("CONTINUE", s1, 8) == 0) {
    scanf("%s", s2);  fprintf(out, "#CONTINUE     \t%s\n", s2);
  } else fprintf(out, "#CONTINUE     \t%s\n", CONTINUE);

  if (strncmp("DIIS_SIZE", s1, 9) == 0) {
    scanf("%d", &d);  fprintf(out, "#DIIS_SIZE    \t%d\n", d);
  } else fprintf(out, "#DIIS_SIZE    \t%d\n", DIIS_SIZE);

  if (strncmp("DIIS_MP", s1, 7) == 0) {
    scanf("%d", &f);  fprintf(out, "#DIIS_MP      \t%f\n", f);
  } else fprintf(out, "#DIIS_MP      \t%f\n", DIIS_MP);

  if (strncmp("EWALD_SUMS", s1,10) == 0) {
    scanf("%s", s2);  fprintf(out, "#EWALD_SUMS   \t%s\n", s2);
  } else fprintf(out, "#EWALD_SUMS   \t%s\n", EWALD_SUMS);

  if (strncmp("GAUSS_SMEAR", s1,11) == 0) {
    scanf("%lf", &f);  fprintf(out, "#GAUSS_SMEAR  \t%f\n", f);
  } else fprintf(out, "#GAUSS_SMEAR  \t%f\n", GAUSS_SMEAR);

  if (strncmp("INIT_PIC_ITER", s1,13) == 0) {
    scanf("%d", &d);  fprintf(out, "#INIT_PIC_ITER\t%d\n", d);
  } else fprintf(out, "#INIT_PIC_ITER\t%d\n", INIT_PIC_ITER);

  if (strncmp("MAX_ITER", s1, 8) == 0) {
    scanf("%d", &d);  fprintf(out, "#MAX_ITER     \t%d\n", d);
  } else fprintf(out, "#MAX_ITER     \t%d\n", MAX_ITER);

  if (strncmp("N_DUMP", s1, 6) == 0) {
    scanf("%d", &d);  fprintf(out, "#N_DUMP       \t%d\n", d);
  } else fprintf(out, "#N_DUMP       \t%d\n", N_DUMP);

  if (strncmp("PIC_MP", s1, 6) == 0) {
    scanf("%lf", &f);  fprintf(out, "#PIC_MP       \t%f\n", f);;
  } else fprintf(out, "#PIC_MP       \t%f\n", PIC_MP);

  if (strncmp("SOLVER", s1, 6) == 0) {
    scanf("%s", s2);  fprintf(out, "#SOLVER       \t%s\n", s2);
  } else fprintf(out, "#SOLVER       \t%s\n", SOLVER);

  if (strncmp("T_ERR", s1, 5) == 0) {
    scanf("%lf", &f);  fprintf(out, "#T_ERR        \t%f\n", f);
  } else fprintf(out, "#T_ERR        \t%f\n", T_ERR);

  /*if( strncmp( "", s1 , ) == 0 ) scanf( "%", &);   = ;*/
}



int select_var(void)
{
  int i = 0;
  int sel;

  printf("\n%s-%d", "NX", ++i);
  printf("\n%s-%d", "NY", ++i);
  printf("\n%s-%d", "NZ", ++i);
  printf("\n%s-%d", "LX", ++i);
  printf("\n%s-%d", "LY", ++i);               /*5*/
  printf("\n%s-%d", "LZ", ++i);
  printf("\n%s-%d", "CX", ++i);
  printf("\n%s-%d", "CY", ++i);
  printf("\n%s-%d", "CZ", ++i);                /*9*/
  printf("\n%s-%d", "A_ERF", ++i);               /*10*/
  printf("\n%s-%d", "CHRG_PCT", ++i);
  printf("\n%s-%d", "CLOSURE", ++i);
  printf("\n%s-%d", "CONFIG_TYPE", ++i);
  printf("\n%s-%d", "CONTINUE", ++i);
  printf("\n%s-%d", "DIIS_SIZE", ++i);               /*15*/
  printf("\n%s-%d", "DIIS_MP", ++i);
  printf("\n%s-%d", "EWALD_SUMS", ++i);
  printf("\n%s-%d", "GAUSS_SMEAR", ++i);
  printf("\n%s-%d", "INIT_PIC_ITER", ++i);
  printf("\n%s-%d", "MAX_ITER", ++i);               /*20*/
  printf("\n%s-%d", "N_DUMP", ++i);
  printf("\n%s-%d", "PIC_MP", ++i);
  printf("\n%s-%d", "SOLVER", ++i);
  printf("\n%s-%d", "T_ERR", ++i);               /*24*/

  printf("\nEnter choice: ");
  scanf("%d", &sel);
}



