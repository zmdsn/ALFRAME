/*
  Test function suite (last update:
  Mallipeddi Rammohan (email: mallipeddi.ram@gmail.com

  for linux:
  
  gcc -c fcnsuite.c
  ld -o fcnsuite.so -shared fcnsuite.o

  for windows:

  gcc -c fcnsuite.c -DWINDOWS
  dllwrap -o fcnsuite.dll fcnsuite.o
*/

#ifdef WINDOWS
#define DLLIMPORT __declspec (dllexport)
#else
#define DLLIMPORT
#endif

#include <math.h>


DLLIMPORT void
C01 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C02 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C03 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C04 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);



DLLIMPORT void
C05 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C06 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C07 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C08 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C09 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C10 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C11 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C12 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);

DLLIMPORT void
C13 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C14 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C15 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);

DLLIMPORT void
C16 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);


DLLIMPORT void
C17 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);

DLLIMPORT void
C18 (double *x, double *f, double *g, double *h, int nx, int nf, int ng, int nh);

