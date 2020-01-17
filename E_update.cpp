#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void E_update(double*** E_r, double*** E_theta, double*** E_phi, 
double*** newD_r, double*** newD_theta, double*** newD_phi, 
double*** oldD_r, double*** oldD_theta, double*** oldD_phi)
{
  double maxE;
  int I, J, K;
  
  for(int i = 0; i < Nr; i++){
    for(int j = 0; j < Ntheta + 1; j++){
      for(int k = 0; k < Nphi + 1; k++){
        E_r[i][j][k] = E_r[i][j][k] + (newD_r[i][j][k] - oldD_r[i][j][k])/EPS0;
        if(std::abs(E_phi[i][j][k])>std::abs(maxE)){
          maxE = E_phi[i][j][k];
          I = i;
          J = j;
          K = k;
        }
      }
    }
  }
  
  for(int i = 0; i < Nr + 1; i++){
    for(int j = 0; j < Ntheta; j++){
      for(int k = 0; k < Nphi + 1; k++){
        E_theta[i][j][k] = E_theta[i][j][k] + (newD_theta[i][j][k] - oldD_theta[i][j][k])/EPS0;
      }
    }
  }
  
  for(int i = 0; i < Nr + 1; i++){
    for(int j = 0; j < Ntheta + 1; j++){
      for(int k = 0; k < Nphi; k++){
        E_phi[i][j][k] = E_phi[i][j][k] + (newD_phi[i][j][k] - oldD_phi[i][j][k])/EPS0;
      }
    }
  }

  std::cout << I << " " << J << " " << K << " maxE = " << maxE << std::endl;
  
}










