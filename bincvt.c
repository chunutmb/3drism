#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "binio.h"
#include "jh_get.c"
#include "jh_print.c"

char *fnin;
char *fnout;
int cpx = 0;  /* complex */



/* check if a file's extension is ext */
int endswith(const char *fn, const char *ext)
{
  const char *p;

  if ((p = strchr(fn, '.')) == NULL) return 0;
  return (strcmp(p + 1, ext) == 0);
}



int main(int argc, char **argv)
{
  double *arr;
  ENV_PAR sys;
  double temp = 298.15, den = 0.001;
  int type = 0; /* 1: binary, 0: text */
  clock_t t0;

  if (argc < 3) {
    fprintf(stderr, "%s from_file to_file\n", argv[0]);
    return -1;
  }
  if (argc >= 4) cpx = atoi(argv[3]);
  fnin = argv[1];
  fnout = argv[2];

  t0 = clock();
  if (endswith(fnin, "bin3d") || endswith(fnin, "bin")) { /* binary input */
    arr = readbin3d(fnin, &sys, cpx, 1, &cpx, &temp, &den);
    type = 1;
  } else { /* text input */
    arr = readjh3d(fnin, &sys, &cpx, &temp, &den);
    type = 0;
  }
  fprintf(stderr, "loaded %s file %s, %s, time %gs\n", 
      type ? "binary" : "text", fnin, cpx ? "complex" : "real",
      1.*(clock() - t0)/CLOCKS_PER_SEC);

  t0 = clock();
  if (endswith(fnout, "bin3d") || endswith(fnout, "bin")) { /* binary output */
    writebin3d(arr, fnout, &sys, cpx, temp, den);
    type = 1;
  } else { /* text output */
    writejh3d(arr, fnout, &sys, cpx, temp, den);
    type = 0;
  }
  fprintf(stderr, "saved %s file %s, %s, time %gs\n", 
      type ? "binary" : "text", fnout, cpx ? "complex" : "real",
      1.*(clock() - t0)/CLOCKS_PER_SEC);

  return 0; 
}
