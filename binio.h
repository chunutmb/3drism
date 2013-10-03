#ifndef BINIO_H__
#define BINIO_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
/* definition of ENV_PAR */
#include "jh_struct.h"


/* binary input and output */


typedef double dblpair[2];

/* convenience wrappers */
#define writebin3dcomplex(arr, fn, sys, temp, den) \
  writebin3d((const double *) arr, fn, sys, 1, temp, den)
#define writebin3dreal(arr, fn, sys, temp, den) \
  writebin3d(arr, fn, sys, 0, temp, den)
#define readbin3dcomplex(fn, sys) \
  ((dblpair *) readbin3d(fn, sys, 1, 0, NULL, NULL, NULL))
#define readbin3dreal(fn, sys) \
  readbin3d(fn, sys, 0, 0, NULL, NULL, NULL)


/* write a 3D array from a binary file
 *  cpx = 1 for complex file, 0 for regular file
 * format:
 *  set nx ny nz
 *  lx ly lz temp density
 *  array
 * `set' is 1 for real or 2 for complex files */
__inline static int writebin3d(const double *arr, const char *fn,
   const ENV_PAR *sys, int cpx, double temp, double den)
{
  FILE *fp;
  int ival[4] = {cpx + 1, sys->nx, sys->ny, sys->nz};
  size_t cnt;
  double dval[5] = {sys->lx, sys->ly, sys->lz, temp, den};

  if ((fp = fopen(fn, "wb")) == NULL) {
    fprintf(stderr, "cannot write %s\n", fn);
    return -1;
  }
  fwrite(ival, sizeof(int), 4, fp);
  fwrite(dval, sizeof(double), 5, fp);
  cnt = sys->nx * sys->ny * sys->nz;
  if (cpx) cnt *= 2;
  fwrite(arr, sizeof(double), cnt, fp);
  return 0;
}




/* read a 3D array from a binary file
 *  cpx = 1 for complex file, 0 for regular file
 *  loadsys = 1 to load information into `sys'
 * format:
 *  version nx ny nz
 *  lx ly lz temp density
 *  array */
__inline static double *readbin3d(const char *fn, ENV_PAR *sys, int cpx,
    int loadsys, int *iscpx, double *temp, double *den)
{
  FILE *fp;
  int ival[4] = {0};
  size_t cnt;
  double dval[5] = {0}, *arr = NULL;

  if ((fp = fopen(fn, "rb")) == NULL) {
    fprintf(stderr, "cannot read %s\n", fn);
    return NULL;
  }
  if (fread(ival, sizeof(int), 4, fp) != 4
      || (!loadsys && ival[0] != cpx + 1)
      || (!loadsys && 
         (ival[1] != sys->nx
       || ival[2] != sys->ny
       || ival[3] != sys->nz) ) ) {
    fprintf(stderr, "%s: error reading int parameters %d vs. %d, (%d, %d, %d) vs (%d, %d, %d)\n",
        fn, ival[0], cpx + 1, ival[1], ival[2], ival[3],
        sys->nx, sys->ny, sys->nz);
    return NULL;
  }
  if (cpx < 0) cpx = ival[0] - 1;

  if (fread(dval, sizeof(double), 5, fp) != 5 
      || (!loadsys && 
         (fabs(dval[0] - sys->lx) > 0.001
       || fabs(dval[1] - sys->ly) > 0.001
       || fabs(dval[2] - sys->lz) > 0.001) ) ) {
    fprintf(stderr, "%s: error reading double parameters (%g, %g, %g) vs (%g, %g, %g)\n",
        fn, dval[0], dval[1], dval[2], sys->lx, sys->ly, sys->lz);
    return NULL;
  }
  
  if (loadsys) {
    sys->nx = ival[1];
    sys->ny = ival[2];
    sys->nz = ival[3];
    sys->lx = dval[0];
    sys->ly = dval[1];
    sys->lz = dval[2];
  }
  if (iscpx != NULL) *iscpx = cpx;
  if (temp != NULL) *temp = dval[3];
  if (den != NULL) *den = dval[4];

  cnt = ival[1] * ival[2] * ival[3];
  if (cpx) cnt *= 2;
  if ((arr = calloc(cnt, sizeof(double))) == NULL) {
    fprintf(stderr, "%s: cannot allocate memory", fn);
    return NULL;
  }
  if (fread(arr, sizeof(double), cnt, fp) != cnt) {
    fprintf(stderr, "%s: cannot read array", fn);
    return NULL;
  }
  return arr;
}



__inline int writejh3d(const double *arr, const char *fn, ENV_PAR *sys,
   int cpx, double temp, double den)
{
  int x, y, z, id;
  FILE *fp;

  if ((fp = fopen(fn, "w")) == NULL) {
    printf("cannot write %s\n", fn);
    return -1;
  }

  fprintf(fp, "%d\n%d\n%d\n", sys->nx, sys->ny, sys->nz);
  fprintf(fp, "%.10f\n%.10f\n%.10f\n", sys->lx, sys->ly, sys->lz);
  fprintf(fp, "%.10f\n", temp);
  fprintf(fp, "%.10f\n", den);

  for (x = 0; x < sys->nx; x++) {
    for (y = 0; y < sys->ny; y++) {
      for (z = 0; z < sys->nz; z++) {
        id = sys->nz * sys->ny * x + sys->nz * y + z;
        if (cpx) {
          fprintf(fp, "%d\t%d\t%d\t%.14e\t%.14e\n", x, y, z, arr[id*2], arr[id*2+1]);
        } else {
          fprintf(fp, "%d\t%d\t%d\t%.15e\n", x, y, z, arr[id]);
        }
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);
  return 0;
}



__inline double *readjh3d(const char *fn, ENV_PAR *sys,
    int *cpx, double *temp, double *den)
{
  int x, y, z, nx, ny, nz, id, cc;
  double lx, ly, lz, temp2, den2;
  double val = 0, val2 = 0;
  double *arr;
  FILE *fp;
  char s[128];

  if ((fp = fopen(fn, "r")) == NULL) {
    fprintf(stderr, "%s: can't open file\n", fn);
    return NULL;
  }

#define GETINT(a) { fgets(s, sizeof s, fp); sscanf(s, "%d", &a); }
#define GETDBL(a) { fgets(s, sizeof s, fp); sscanf(s, "%lf", &a); }

  GETINT(nx); GETINT(ny); GETINT(nz);
  GETDBL(lx); GETDBL(ly); GETDBL(lz);
  GETDBL(temp2); GETDBL(den2);

  sys->nx = nx;
  sys->ny = ny;
  sys->nz = nz;
  sys->lx = lx;
  sys->ly = ly;
  sys->lz = lz;
  if (temp != NULL) *temp = temp2;
  if (den != NULL) *den = den2;

  if ((arr = calloc(2 * nx * ny * nz, sizeof(double))) == NULL) {
    fprintf(stderr, "%s: no memory\n", fn);
    return NULL;
  }
  while (fgets(s, sizeof s, fp)) {
    if ( isspace(s[0]) ) continue; /* a blank line */
    cc = sscanf(s, "%d%d%d%lf%lf", &x, &y, &z, &val, &val2);
    id = nz * ny * x + nz * y + z;
    if (cc == 4) {
      arr[id] = val;
    } else if (cc == 5) {
      arr[id*2] = val;
      arr[id*2 + 1] = val2;
    } else { /* error */
      printf("cc %d:%s, id %d (%d, %d, %d)\n", cc, s, id, x, y, z);
      return arr;
    }
    if (id == 0) *cpx = cc - 4;
  }
  return arr;
}




#endif

