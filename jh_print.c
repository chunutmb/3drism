#include "jh_struct.h"



void print_jh3d(char name[], double *v, ENV_PAR sys, double temp, double pnd)
{
  int x, y, z;
  int nx = sys.nx;
  int ny = sys.ny;
  int nz = sys.nz;

  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");

  fprintf(out, "%d\n%d\n%d\n", sys.nx, sys.ny, sys.nz);
  fprintf(out, "%.10f\n%.10f\n%.10f\n", sys.lx, sys.ly, sys.lz);
  fprintf(out, "%.10f\n", temp);
  fprintf(out, "%.10f\n", pnd);

  for (x = 0; x <= nx - 1; x++) {
    for (y = 0; y <= ny - 1; y++) {
      for (z = 0; z <= nz - 1; z++)
        fprintf(out, "%d\t%d\t%d\t%.15e\n", x, y, z, v[nz * ny * x + nz * y + z]);
      fprintf(out, "\n");
    }
  }
  fclose(out);
}



void print_sit(char name[], double *v, ENV_PAR sys)
{
  int x, y, z;
  int nx = sys.nx;
  int ny = sys.ny;
  int nz = sys.nz;
  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");

  double dx = sys.lx / (nx - 1);
  double dy = sys.ly / (ny - 1);
  double dz = sys.lz / (nz - 1);

  double fx = (-1.0 * sys.cx) * dx;
  double fy = (-1.0 * sys.cy) * dy;
  double fz = (-1.0 * sys.cz) * dz;

  fprintf(out, "%.15f\n", dx);
  fprintf(out, "%.10f\t%.10f\t%.10f\n", fx, fy, fz);
  fprintf(out, "%d\t%d\t%d\n", nx, ny, nz);

  for (z = 0; z <= nz - 1; z++)
    for (y = 0; y <= ny - 1; y++)
      for (x = 0; x <= nx - 1; x++)
        fprintf(out, "%.15e\n", v[nz * ny * x + nz * y + z]);

  fclose(out);
}



void print_cmplx_jh3d(char name[], fftw_complex *v, ENV_PAR sys, double temp, double pnd)
{
  int x, y, z;
  int nx = sys.nx;
  int ny = sys.ny;
  int nz = sys.nz;

  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");

  fprintf(out, "%d\n%d\n%d\n", nx, ny, nz);
  fprintf(out, "%.10f\n%.10f\n%.10f\n", sys.lx, sys.ly, sys.lz);
  fprintf(out, "%.10f\n", temp);
  fprintf(out, "%.10f\n", pnd);

  for (x = 0; x <= nx - 1; x++) {
    for (y = 0; y <= ny - 1; y++) {
      for (z = 0; z <= nz - 1; z++)
        fprintf(out, "%d\t%d\t%d\t%.14e\t%.14e\n", x, y, z, v[nz * ny * x + nz * y + z][0], v[nz * ny * x + nz * y + z][1]);
      fprintf(out, "\n");
    }
  }
  fclose(out);
}



void print_cmplx_sit(char name[], fftw_complex *v, ENV_PAR sys)
{
  int x, y, z;
  int nx = sys.nx;
  int ny = sys.ny;
  int nz = sys.nz;

  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");

  double dx = sys.lx / (nx - 1);
  double dy = sys.ly / (ny - 1);
  double dz = sys.lz / (nz - 1);

  double fx = (-1.0 * sys.cx) * dx;
  double fy = (-1.0 * sys.cy) * dy;
  double fz = (-1.0 * sys.cz) * dz;

  fprintf(out, "%.15f\n", dx);
  fprintf(out, "%.10f\t%.10f\t%.10f\n", fx, fy, fz);
  fprintf(out, "%d\t%d\t%d\n", nx, ny, nz);

  for (z = 0; z <= nz - 1; z++)
    for (y = 0; y <= ny - 1; y++)
      for (x = 0; x <= nx - 1; x++)
        fprintf(out, "%.14e\t%.14e\n", v[nz * ny * x + nz * y + z][0], v[nz * ny * x + nz * y + z][1]);
  fclose(out);
}

