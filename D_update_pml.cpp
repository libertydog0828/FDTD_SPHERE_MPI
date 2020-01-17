#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <omp.h>
#include "fdtd3d.h"

void D_update_pml(double*** newD_r, double*** newD_theta, double*** newD_phi,
                  double*** H_r, double*** H_theta, double*** H_phi, 
                  double**** Dr_theta1, double**** Dr_theta2, double**** Dr_phi, 
                  double**** Dtheta_phi, double**** Dtheta_r, 
                  double**** Dphi_r, double**** Dphi_theta,
                  double* sigma_theta, double* sigma_phi)
{
  double ri_1, ri_2, ri_3;
  double thetaj_1, thetaj_2;
  
  pml* idx_r = new pml[4];
  pml* idx_theta = new pml[4];
  pml* idx_phi = new pml[4];
  
  idx_r[0].set_point_1(1, 1);
  idx_r[0].set_point_2(L, Nphi - 1);
  idx_r[1].set_point_1(Ntheta - L, 1);
  idx_r[1].set_point_2(Ntheta - 1, Nphi - 1);
  idx_r[2].set_point_1(L + 1, 1);
  idx_r[2].set_point_2(Ntheta - L - 1, L);
  idx_r[3].set_point_1(L + 1, Nphi - L);
  idx_r[3].set_point_2(Ntheta - L - 1, Nphi - 1);

  idx_theta[0].set_point_1(0, 1);
  idx_theta[0].set_point_2(L - 1, Nphi - 1);
  idx_theta[1].set_point_1(Ntheta - L, 1);
  idx_theta[1].set_point_2(Ntheta - 1, Nphi - 1);
  idx_theta[2].set_point_1(L, 1);
  idx_theta[2].set_point_2(Ntheta - L - 1, L);
  idx_theta[3].set_point_1(L, Nphi - L);
  idx_theta[3].set_point_2(Ntheta - L - 1, Nphi - 1);
  
  idx_phi[0].set_point_1(1, 0);
  idx_phi[0].set_point_2(L, Nphi - 1);
  idx_phi[1].set_point_1(Ntheta - L, 0);
  idx_phi[1].set_point_2(Ntheta - 1, Nphi - 1);
  idx_phi[2].set_point_1(L + 1, 0);
  idx_phi[2].set_point_2(Ntheta - L - 1, L - 1);
  idx_phi[3].set_point_1(L + 1, Nphi - L);
  idx_phi[3].set_point_2(Ntheta - L - 1, Nphi - 1);
  
  //Update Dr using Dr_theta1, Dr_theta2, Dr_phi//
  for(int area = 0; area < 4; area++){
    for(int i = 0; i <= Nr - 1; i++){
      ri_2 = dist(i + 0.5);
      for(int j = idx_r[area].j1; j <= idx_r[area].j2; j++){
        int j_area = j - idx_r[area].j1;
        thetaj_1 = th(j);
        for(int k = idx_r[area].k1; k <= idx_r[area].k2; k++){
          int k_area = k - idx_r[area].k1;
          Dr_theta1[area][i][j_area][k_area] = C_1(sigma_theta[j])*Dr_theta1[area][i][j_area][k_area]
            + C_2(ri_2, sigma_theta[j])*(H_phi[i][j][k] - H_phi[i][j-1][k]);

          Dr_theta2[area][i][j_area][k_area] = Dr_theta2[area][i][j_area][k_area]
            + C_3(ri_2, thetaj_1)*(H_phi[i][j][k] + H_phi[i][j-1][k]);

          Dr_phi[area][i][j_area][k_area] = C_1(sigma_phi[k])*Dr_phi[area][i][j_area][k_area]
            - C_4(ri_2, thetaj_1, sigma_phi[k])*(H_theta[i][j][k] - H_theta[i][j][k-1]);

          newD_r[i][j][k] = Dr_theta1[area][i][j_area][k_area]
            + Dr_theta2[area][i][j_area][k_area] + Dr_phi[area][i][j_area][k_area];
        }
      }
    }
  }
  
  //Update Dtheta using Dtheta_phi, Dtheta_r//
  for(int area = 0; area < 4; area++){
    for(int i = 1; i <= Nr - 1; i++){

      ri_1 = dist(i - 0.5);
      ri_2 = dist(i);
      ri_3 = dist(i + 0.5);
      for(int j = idx_theta[area].j1; j <= idx_theta[area].j2; j++){
        int j_area = j - idx_theta[area].j1;
        thetaj_2 = th(j + 0.5);
        for(int k = idx_theta[area].k1; k <= idx_theta[area].k2; k++){
          int k_area = k - idx_theta[area].k1;
          Dtheta_phi[area][i][j_area][k_area] = C_1(sigma_phi[k])*Dtheta_phi[area][i][j_area][k_area]
          + C_4(ri_2, thetaj_2, sigma_phi[k])*(H_r[i][j][k] - H_r[i][j][k-1]);

          Dtheta_r[area][i][j_area][k_area] = Dtheta_r[area][i][j_area][k_area]
          - C_5(ri_2)*(ri_3*H_phi[i][j][k] - ri_1*H_phi[i-1][j][k]);

          newD_theta[i][j][k] = Dtheta_phi[area][i][j_area][k_area] + Dtheta_r[area][i][j_area][k_area];
        }
      }
    }
  }
  
  //Update Dphi using Dphi_r, Dphi_theta//
  for(int area = 0; area < 4; area++){
    for(int i = 1; i <= Nr - 1; i++){
      ri_1 = dist(i - 0.5);
      ri_2 = dist(i);
      ri_3 = dist(i + 0.5);
      for(int j = idx_phi[area].j1; j <= idx_phi[area].j2; j++){
        int j_area = j - idx_phi[area].j1;
        for(int k = idx_phi[area].k1; k <= idx_phi[area].k2; k++){
          int k_area = k - idx_phi[area].k1;
          Dphi_r[area][i][j_area][k_area] = Dphi_r[area][i][j_area][k_area]
            + C_5(ri_2)*(ri_3*H_theta[i][j][k] - ri_1*H_theta[i-1][j][k]);
          
          Dphi_theta[area][i][j_area][k_area] = C_1(sigma_theta[j])*Dphi_theta[area][i][j_area][k_area]
            - C_6(ri_2, sigma_theta[j])*(H_r[i][j][k] - H_r[i][j-1][k]);
          
          newD_phi[i][j][k] = Dphi_r[area][i][j_area][k_area] + Dphi_theta[area][i][j_area][k_area];
        }
      }
    }
  }

}





