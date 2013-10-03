#include <stdio.h>
#include <math.h>
#include "jh_get.h"
#include "jh_struct.h"

#define NNN NX * NY * NZ
#define NNN2 NX2 * NY2 * NZ2
#define ii(x,y,z) (NZ * NY * (x) + NZ * (y) + (z))
#define ii2(x,y,z) (NZ2 * NY2 * (x) + NZ2 * (y) + (z))

int NX, NY, NZ;
int NX2, NY2, NZ2;
int CX = 0.0, CY = 0.0, CZ = 0.0;
double LX, LY, LZ;
double TEMP, PND;

double *Vec;
double *Vec2;

ENV_PAR SYS;

void set_sys(char []);
void convert_grid(void);
void print_jh3d(char []);

int main(int argc, char *argv[])
{
  char *infile = *(argv + 1);

  set_sys(infile);

  Vec = (double *) get_jh3d_dat(infile, TEMP, PND, SYS);

  printf("\nEnter new NX: "); fscanf(stdin, "%d", &NX2);
  printf("\nEnter new NY: "); fscanf(stdin, "%d", &NY2);
  printf("\nEnter new NZ: "); fscanf(stdin, "%d", &NZ2);
  printf("\n");

  Vec2 = (double *)malloc(NNN2 * sizeof(double));

  convert_grid();

  print_jh3d(infile);

  return 0;
}



void set_sys(char infile[])
{
  FILE *in;

  if ((in = fopen(infile, "r")) == NULL)
    fprintf(stderr, "\nCan't open file:%s\n", infile);

  /*Begin reading in jh3d.dat file*/
  fscanf(in, "%d%d%d", &NX, &NY, &NZ);
  fscanf(in, "%lf%lf%lf", &LX, &LY, &LZ);
  fscanf(in, "%lf", &TEMP);
  fscanf(in, "%lf", &PND);

  SYS.nx = NX;    SYS.ny = NY;    SYS.nz = NZ;
  SYS.cx = CX;    SYS.cy = CY;    SYS.cz = CZ;
  SYS.lx = LX;    SYS.ly = LY;    SYS.lz = LZ;

  fclose(in);
}



void convert_grid(void)
{
  int x, y, z, id;
  int i, j;
  int ix, iy, iz;
  double rx, ry, rz, drx, dry, drz;
  double ddx, ddy, ddz;
  double tmp;
  double dx = ((double) 1 / (NX - 1)), dy = ((double) 1 / (NY - 1)), dz = ((double) 1 / (NZ - 1));
  double dx2 = ((double) 1 / (NX2 - 1)), dy2 = ((double) 1 / (NY2 - 1)), dz2 = ((double) 1 / (NZ2 - 1));

  for (x = 0; x <= NX2 - 1; x++)
    for (y = 0; y <= NY2 - 1; y++)
      for (z = 0; z <= NZ2 - 1; z++) {
        rx = x * dx2;    ix = (int) rx / dx;       drx = rx - (ix * dx);
        ry = y * dy2;    iy = (int) ry / dy;       dry = ry - (iy * dy);
        rz = z * dz2;    iz = (int) rz / dz;       drz = rz - (iz * dz);
        ddx = dx - drx;
        ddy = dy - dry;
        ddz = dz - drz;

        id = ii2(x,y,z);

        Vec2[id] = ddx * ddy * ddz * Vec[ii(ix,iy,iz)];

        if (iz + 1 <= NZ - 1)
          Vec2[id] += ddx * ddy * drz * Vec[ii(ix,iy,iz + 1)];

        if (iy + 1 <= NY - 1)
          Vec2[id] += ddx * dry * ddz * Vec[ii(ix,iy + 1,iz)];

        if (iz + 1 <= NZ - 1)
          if (iy + 1 <= NY - 1)
            Vec2[id] += ddx * dry * drz * Vec[ii(ix,iy + 1,iz + 1)];

        if (ix + 1 <= NX - 1) {
          Vec2[id] += drx * ddy * ddz * Vec[ii(ix + 1,iy,iz)];

          if (iz + 1 <= NZ - 1)
            Vec2[id] += drx * ddy * drz * Vec[ii(ix + 1,iy,iz + 1)];

          if (iy + 1 <= NY - 1)
            Vec2[id] += drx * dry * ddz * Vec[ii(ix + 1,iy + 1,iz)];

          if (iz + 1 <= NZ - 1)
            if (iy + 1 <= NY - 1)
              Vec2[id] += drx * dry * drz * Vec[ii(ix + 1,iy + 1,iz + 1)];
        }

        Vec2[id] *= 1 / (dx * dy * dz);
      }
}



/**/ void print_jh3d(char name[])
{
  int x, y, z;
  FILE *out;
  if ((out = fopen(name, "w")) == NULL)
    printf("\nFile could not be opened\n");

  fprintf(out, "%d\n%d\n%d\n", NX2, NY2, NZ2); fflush(out);
  fprintf(out, "%.10f\n%.10f\n%.10f\n", LX, LY, LZ); fflush(out);
  fprintf(out, "%.10f\n", TEMP); fflush(out);
  fprintf(out, "%.10f\n", PND); fflush(out);

  for (x = 0; x <= NX2 - 1; x++) {
    for (y = 0; y <= NY2 - 1; y++) {
      for (z = 0; z <= NZ2 - 1; z++)
        fprintf(out, "%d\t%d\t%d\t%.15e\n", x, y, z, Vec2[NZ2 * NY2 * x + NZ2 * y + z]);
      fprintf(out, "\n");
    }
    fflush(out);
  }
  fclose(out);
}



