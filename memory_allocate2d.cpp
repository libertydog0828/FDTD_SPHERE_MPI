#include <iostream>
#include <cmath>

double** memory_allocate2d(int Num)
{
  double** array = new double*[Num];

  for(int i = 0; i < Num; i++){
      array[i] = new double[Num];
      for(int j = 0; j < Num; j++){
        array[i][j] = 0.0;
      }
  }

  return array;
}
