#include <iostream>
#include <cmath>
#include "fdtd3d.h"

double**** memory_allocate4d(int m, int n, int o, int p)
{

  double ****array;
  array = new double***[m];
  for(int i = 0; i < m; i++){
    array[i] = new double**[n];
    for(int j = 0; j < n; j++){
      array[i][j] = new double*[o];
      for(int k = 0; k < o; k++){
        array[i][j][k] = new double[p];
        for(int l = 0; l < p; l++){
          array[i][j][k][l] = 0.0;
        }
      }
    }
  }

  return array;
  
}
