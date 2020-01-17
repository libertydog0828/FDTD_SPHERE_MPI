#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <omp.h>
#include "fdtd3d.h"

void D_update(double*** newD_r, double*** newD_theta, double*** newD_phi, 
              double*** H_r, double*** H_theta, double*** H_phi, 
              double*** oldD_r, double*** oldD_theta, double*** oldD_phi)
{
  double val_1, val_2;
  double ri_1, ri_2, ri_3;
  double sin_th1, sin_th2, sin_th3;

  std::cout << "num of threads  " << omp_get_max_threads() << std::endl;
  //D update (outside PML)//
#pragma omp parallel for private(i, j, k)
  for(int i = 0; i <= Nr - 1; i++){
    ri_2 = dist(i + 0.5);
    for(int j = L + 1; j <= Ntheta - L - 1; j++){
      sin_th1 = std::sin(th(j));
      sin_th2 = std::sin(th(j + 0.5));
      sin_th3 = std::sin(th(j - 0.5));
      for(int k = L + 1; k <= Nphi - L - 1; k++){
        val_1 = Dt/ri_2/sin_th1/delta_theta;
        val_2 = Dt/ri_2/sin_th1/delta_phi;

        newD_r[i][j][k] = oldD_r[i][j][k] + val_1*(sin_th2*H_phi[i][j][k] - sin_th3*H_phi[i][j-1][k])
          - val_2*(H_theta[i][j][k] - H_theta[i][j][k-1]);
      }
    }
  }
  
  for(int i = 1; i <= Nr - 1; i++){
    ri_1 = dist(i);
    ri_2 = dist(i + 0.5);
    ri_3 = dist(i - 0.5);
    for(int j = L; j <= Ntheta - L - 1; j++){
      sin_th2 = std::sin(th(j + 0.5));
      for(int k = L + 1; k <= Nphi - L - 1; k++){
        val_1 = Dt/ri_1/sin_th2/delta_phi;
        val_2 = Dt/ri_1/delta_r;

        newD_theta[i][j][k] = oldD_theta[i][j][k] + val_1*(H_r[i][j][k] - H_r[i][j][k-1])
          - val_2*(ri_2*H_phi[i][j][k] - ri_3*H_phi[i-1][j][k]);
      }
    }
  }
  
  for(int i = 1; i <= Nr - 1; i++){
    ri_1 = dist(i);
    ri_2 = dist(i + 0.5);
    ri_3 = dist(i - 0.5);
    for(int j = L + 1; j <= Ntheta - L - 1; j++){
      for(int k = L; k <= Nphi - L - 1; k++){
        val_1 = Dt/ri_1/delta_r;
        val_2 = Dt/ri_1/delta_theta;

        newD_phi[i][j][k] = oldD_phi[i][j][k] + val_1*(ri_2*H_theta[i][j][k] - ri_3*H_theta[i-1][j][k])
          - val_2*(H_r[i][j][k] - H_r[i][j-1][k]); 
      }
    }
  }
  
}
