/* conversion between .jh3d and .bin3d files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "binio.h"



/* check if a file's extension is ext */
int endswith(const char *fn, const char *ext)
{
  const char *p;

  if ((p = strchr(fn, '.')) == NULL) return 0;
  return (strncmp(p + 1, ext, strlen(ext)) == 0);
}



/* convert fnin to fnout, both files can be either text or binary */
static void bincvt(const char *fnin, const char *fnout)
{
  /* static variables are used to retain values from previous calls */
  static ENV_PAR sys;
  static double temp = 298.15, den = 0.001;
  static int cpx = 0;
  double *arr;
  int type = 0; /* 1: binary, 0: text */
  clock_t t0;

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
  free(arr);
}



/* deduce the output file name from the input file name */
static char *mkdefout(const char *fnin)
{
  char *fnout, *p;

  if ((fnout = malloc(strlen(fnin) + 8)) == NULL) return NULL;
  strcpy(fnout, fnin);
  p = strchr(fnout, '.');
  if (p != NULL) {
    p++; /* move behind the dot */
  } else {
    p = fnout + strlen(fnout);  /* move to the end end of the string */
  }
  if (strncmp(p, "bin", 3) == 0) { /* binary to text */
    strcpy(p, "jh3d");
  } else { /* text to binary */
    strcpy(p, "bin3d");
  }
  printf("deducing the default output for %s is %s\n", fnin, fnout);
  return fnout;
}



/* convert every file in the list */
static void dolist(const char *fnls)
{
  FILE *fp;
  char fn[FILENAME_MAX], *p, *fnout;
  int i = 0;

  if ((fp = fopen(fnls, "r")) == NULL) {
    fprintf(stderr, "cannot open the list %s\n", fnls);
    return;
  }
  for (i = 1; fgets(fn, sizeof fn, fp); i++) {
    p = fn + strlen(fn) - 1;
    while (isspace(*p)) *p-- = '\0'; /* remove trailing spaces */
    fprintf(stderr, "processing file %d: %s\n", i, fn);
    fnout = mkdefout(fn);
    bincvt(fn, fnout);
    free(fnout);
  }
  fclose(fp);
}



int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "%s from_file [to_file]\nor\n", argv[0]);
    fprintf(stderr, "%s -[l|L] [list_file]\n", argv[0]);
    fprintf(stderr, "-l: text to binary, -L: binary to text\n");
    return -1;
  }

  if (argv[1][0] == '-') { /* handle a list */
    char *fnls = argv[1] + 2;
    if (*fnls == '\0') {
      if (argc >= 3) {
        fnls = argv[2];
      } else if (argv[1][1] == 'l') { /* -l means text to binary */
#define TMPFN "tmp.ls"
        fnls = TMPFN;
        /* using the ls command to list all files */
        system("ls --color=no *.jh3d > " TMPFN);
      } else if (argv[1][1] == 'L') {
        fnls = TMPFN;
        system("ls --color=no *.bin* > " TMPFN);
      }
    }
    dolist(fnls);
    if (strcmp(fnls, TMPFN) == 0)
      remove(fnls);
  } else { /* a pair of input and output */
    char *fnin, *fnout;

    fnin = argv[1];
    if (argc >= 3) {
      fnout = argv[2];
    } else { /* deduce the output file */
      fnout = mkdefout(fnin);
    }
    bincvt(fnin, fnout);
  }
  return 0;
}
