#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <omp.h>

#include "common.h"

#include <unistd.h>

typedef double Real;

/* function prototypes */
Real *createRealArray (int n);
Real **createReal2DArray (int m, int n);
void transpose (Real **A, int m, int n, int size, int bb, int bre);
void fillA (Real **A, Real *V, int *re, int *rd, int m, int n, int size);
void printArray (Real **a, int n, int m);
void fst_(Real *v, int *n, Real *w, int *nn);
void fstinv_(Real *v, int *n, Real *w, int *nn);

int main(int argc, char **argv )
{
  Real **A;
  int n, m, b, re, bb, bre, rank, size;

  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  n = atoi(argv[1]);
  b = floor(n/size);
  re = n - (size-1)*b;
  m = b;
  bb = b*b;
  bre = b*re;

  if(rank+1 == size) {
    m = re;
    bb = bre;
    bre = re*re;
  }

  A = createReal2DArray (m,n);

  int testC = n*b*rank;
  for (int j=0; j < m; j++) {
    for (int i=0; i < n; i++) {
      A[j][i] = testC; //(j+1)+(rank*b);
      testC++;
    }
  }

  transpose(A, m, n, size, bb, bre);

  sleep(1*rank);
  printf("rank: %i\n", rank);
  printArray(A, m, n);

  MPI_Finalize();
  return 0;
}

void transpose (Real **A, int m, int n, int size, int bb, int bre)
{
  int se[size], sd[size], re[size], rd[size];
  Real *V = createRealArray (n*m);

  for (int i = 0; i < size; ++i) {
    se[i] = bb;
    sd[i] = bb*i;
    re[i] = bb;
    rd[i] = bb*i;
  }
  se[size-1] = bre;
  re[size-1] = bre;

  for(int i = 0; i < n; i++) {
    for(int j = 0; j < m; j++) {
      V[j + i*m] = A[j][i];
    }
  }

  MPI_Alltoallv(V, se, sd, MPI_DOUBLE, V, re, rd, MPI_DOUBLE, MPI_COMM_WORLD);
  fillA(A, V, re, rd, m, n, size);
}

void fillA (Real **A, Real *V, int *re, int *rd, int m, int n, int size)
{
  int rem[size];

  for (int i = 0; i < size; ++i) {
    rem[i] = re[i]/m;
  }

  for (int i = 0; i < m; i++) {
    int r, k, k1;

    r = 0;
    k = rd[r] + rem[r]*i;
    k1 = k + rem[r] - 1;

    for (int j = 0; j < n; j++) {
      A[i][j] = V[k];

      if(k == k1) {
        r++;
        k = rd[r] + rem[r]*i;
        k1 = k + rem[r] - 1;
      } else {
        k++;
      }
    }
  }
}

Real *createRealArray (int n)
{
  Real *a;
  int i;
  a = (Real *)malloc(n*sizeof(Real));
  
  #pragma omp parallel for schedule(static)
  for (i=0; i < n; i++) {
    a[i] = 0.0;
  }
  return (a);
}

Real **createReal2DArray (int n1, int n2)
{
  int i, n;
  Real **a;
  a    = (Real **)malloc(n1   *sizeof(Real *));
  a[0] = (Real  *)malloc(n1*n2*sizeof(Real));
  
  for (i=1; i < n1; i++) {
    a[i] = a[i-1] + n2;
  }
  n = n1*n2;
  memset(a[0],0,n*sizeof(Real));
  return (a);
}

void printArray (Real **a, int m, int n)
{
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < n; j++) {
      printf("%lf ", a[i][j]);
    }
    printf("\n");
  } 
}