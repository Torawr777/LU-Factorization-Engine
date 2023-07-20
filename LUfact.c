/* Eric Le
 * LU Factorization Engine
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "LUfact.h"

double **createMatrix(int N) {
  double **M = (double **) malloc(N*sizeof(double*));
  for (int i = 0; i < N; i++)
    M[i] = (double*) malloc(N*sizeof(double));
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      M[i][j] = (i == j) ? 1.0 : 0.0;
  return M;
}

void destroyMatrix(int N, double **M) {
  for (int i = 0; i < N; i++)
    free(M[i]);
  free(M);
}

// Performs factorization and allocates, fills and returns a LUfact object
LUfact *LUfactor(int N, const double **A) {
  LUfact *LU = (LUfact*) malloc(sizeof(LUfact));
  LU->N = N;
  LU->LU = createMatrix(N);
  LU->mutate = (short *) malloc(N*sizeof(short));


  // Clone A into LU
  double **A_ = LU->LU;
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N; j++)
      A_[i][j] = A[i][j];

  for (int i = 0; i < N; i++)
    LU->mutate[i] = (short) i;

  // actual factorizing goes here
  double temp = 0;
  int i = 0;

  // Find pivot position
  for(int pos = 0; pos < N-1; pos++) {
    int size = pos;

    // Change mutation to be relevant to position 
    for(; i < N; i++) {
      double c = A[LU->mutate[i]][pos];

      if(fabs(c) > temp) {
        temp = c;
        size = i;
      }
      //If singular matrix
      if(c == 0) {
        LUdestroy(LU);
        return NULL;
      }
    }
    // Switch permutations according to position
    int shift = LU->mutate[pos];
    LU->mutate[pos] = LU->mutate[size];
    LU->mutate[size] = shift;

    // Find value of current position
    int n = pos+1;
    for(; n < N; n++) {
      double k = A_[LU->mutate[n]][pos] / A_[LU->mutate[pos]][pos];
      A_[LU->mutate[n]][pos] = k;

      int r = pos+1;
      for(; r < N; r++) {
        //Increasing values to avoid more bit lost during subtraction
        A_[LU->mutate[n]][r] = A_[LU->mutate[n]][r]*10;
        A_[LU->mutate[n]][r] -= (k*A_[LU->mutate[pos]][r])*10;
        A_[LU->mutate[n]][r] = (A_[LU->mutate[n]][r])/10;
      }
    }
  }
  return LU;
}

// Deallocates the data allocated in LUfactor()
void LUdestroy(LUfact *fact) {
  destroyMatrix(fact->N, fact->LU);
  free(fact->mutate);
  free(fact);
}

/* Solves the system Ax = b for x
 * First for y: Ly = b
 * Then for x: Ux = y
 */ 
void LUsolve(LUfact *fact, const double *b, double *x) {
  int index = fact->N;
  double y[index];
  int i;

  // Solve Ly = b for y (Forward subsitution)
  int pos = 0;
  for(; pos < index; pos++) {
    y[pos] = b[fact->mutate[pos]];
  }

  int pos2 = 0;
  for(; pos2 < index; pos2++) {
    i = pos2+1;

    for(; i < index; i++) {
      y[i] -= y[pos2] * fact->LU[fact->mutate[i]][pos2];
    }
  }
  
  // Solve Ux = y for x (Backward subsitution)
  pos = 0;
  for(; pos < index; pos++) {
    x[pos] = y[pos];
  }

  i = index-1;
  for(; i >= 0; i--) {
    x[i] = x[i] / fact->LU[fact->mutate[i]][i];

    pos2 = i-1;
    for(; pos2 >= 0; pos2--) {
      x[pos2] -= x[i] * fact->LU[fact->mutate[pos2]][i];
    }
  }
}
