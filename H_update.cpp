#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void H_update(double*** E_r, double*** E_theta, double*** E_phi, 
              double*** H_r, double*** H_theta, double*** H_phi)
{
  double val_1, val_2;
  double ri_1, ri_2, ri_3;
  double sin_th1, sin_th2, sin_th3;
  
  //update H outside PML area//
  for(int i = 0; i <= Nr; i++){
    ri_1 = dist(i);
    for(int j = L; j <= Ntheta - L - 1; j++){
      sin_th1 = std::sin(th(j));
      sin_th2 = std::sin(th(j+0.5));
      sin_th3 = std::sin(th(j+1.0));
      for(int k = L; k <= Nphi - L - 1; k++){
        val_1 = Dt/MU0/ri_1/sin_th2/delta_theta;
        val_2 = Dt/MU0/ri_1/sin_th2/delta_phi;
        
        H_r[i][j][k] = H_r[i][j][k] - val_1*(sin_th3*E_phi[i][j+1][k] - sin_th1*E_phi[i][j][k])
          + val_2*(E_theta[i][j][k+1] - E_theta[i][j][k]);
        
      }
    }
  }
  
  for(int i = 0; i <= Nr - 1; i++){
    ri_1 = dist(i);
    ri_2 = dist(i+0.5);
    ri_3 = dist(i+1.0);
    for(int j = L + 1; j <= Ntheta - L - 1; j++){
      sin_th1 = std::sin(th(j));
      for(int k = L; k <= Nphi - L - 1; k++){
        val_1 = Dt/MU0/ri_2/sin_th1/delta_phi;
        val_2 = Dt/MU0/ri_2/delta_r;
        
        H_theta[i][j][k] = H_theta[i][j][k] - val_1*(E_r[i][j][k+1] - E_r[i][j][k])
          + val_2*(ri_3*E_phi[i+1][j][k] - ri_1*E_phi[i][j][k]);
        
      }
    }
  }
  
  for(int i = 0; i <= Nr - 1; i++){
    ri_1 = dist(i);
    ri_2 = dist(i+0.5);
    ri_3 = dist(i+1.0);
    for(int j = L; j <= Ntheta - L - 1; j++){
      for(int k = L + 1; k <= Nphi - L - 1; k++){
        val_1 = Dt/MU0/ri_2/delta_r;
        val_2 = Dt/MU0/ri_2/delta_theta;
        
        H_phi[i][j][k] = H_phi[i][j][k] - val_1*(ri_3*E_theta[i+1][j][k] - ri_1*E_theta[i][j][k])
          + val_2*(E_r[i][j+1][k] - E_r[i][j][k]);
        
      }
    }
  }
  
}






