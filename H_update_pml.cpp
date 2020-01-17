#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void H_update_pml(double*** E_r, double*** E_theta, double*** E_phi,
		  double*** H_r, double*** H_theta, double*** H_phi,
		  double**** Hr_theta1, double**** Hr_theta2, double**** Hr_phi,
		  double**** Htheta_phi, double**** Htheta_r,
		  double**** Hphi_r, double**** Hphi_theta,
		  double* sigma_theta_h, double* sigma_phi_h)
{
  double ri_1, ri_2, ri_3;
  double thetaj_1, thetaj_2;
  
  pml* idx_r = new pml[4];
  pml* idx_theta = new pml[4];
  pml* idx_phi = new pml[4];

  idx_r[0].set_point_1(0, 0);
  idx_r[0].set_point_2(L - 1, Nphi - 1);
  idx_r[1].set_point_1(Ntheta - L, 0);
  idx_r[1].set_point_2(Ntheta - 1, Nphi - 1);
  idx_r[2].set_point_1(L, 0);
  idx_r[2].set_point_2(Ntheta - L - 1, L - 1);
  idx_r[3].set_point_1(L, Nphi - L);
  idx_r[3].set_point_2(Ntheta - L - 1, Nphi - 1);

  idx_theta[0].set_point_1(0, 0);
  idx_theta[0].set_point_2(L, Nphi - 1);
  idx_theta[1].set_point_1(Ntheta - L, 0);
  idx_theta[1].set_point_2(Ntheta, Nphi - 1);
  idx_theta[2].set_point_1(L + 1, 0);
  idx_theta[2].set_point_2(Ntheta - L - 1, L - 1);
  idx_theta[3].set_point_1(L + 1, Nphi - L);
  idx_theta[3].set_point_2(Ntheta - L - 1, Nphi - 1);

  idx_phi[0].set_point_1(0, 0);
  idx_phi[0].set_point_2(L - 1, Nphi);
  idx_phi[1].set_point_1(Ntheta - L, 0);
  idx_phi[1].set_point_2(Ntheta - 1, Nphi);
  idx_phi[2].set_point_1(L, 0);
  idx_phi[2].set_point_2(Ntheta - L - 1, L);
  idx_phi[3].set_point_1(L, Nphi - L);
  idx_phi[3].set_point_2(Ntheta - L - 1, Nphi);
  
  //Update Hr using Hr_theta1, Hr_theta2, Hr_phi//
  for(int area = 0; area < 4; area++){
    for(int i = 0; i <= Nr; i++){
      ri_1 = dist(i);
      for(int j = idx_r[area].j1; j <= idx_r[area].j2; j++){
        int j_area = j - idx_r[area].j1;
        thetaj_2 = th(j + 0.5);
        for(int k = idx_r[area].k1; k <= idx_r[area].k2; k++){
          int k_area = k - idx_r[area].k1;
          Hr_theta1[area][i][j_area][k_area] = C_1(sigma_theta_h[j])*Hr_theta1[area][i][j_area][k_area]
            - C_2(ri_1, sigma_theta_h[j])/MU0*(E_phi[i][j+1][k] - E_phi[i][j][k]);
          
          Hr_theta2[area][i][j_area][k_area] = Hr_theta2[area][i][j_area][k_area]
            - C_3(ri_1, thetaj_2)/MU0*(E_phi[i][j+1][k] + E_phi[i][j][k]);
          
          Hr_phi[area][i][j_area][k_area] = C_1(sigma_phi_h[k])*Hr_phi[area][i][j_area][k_area]
            + C_4(ri_1, thetaj_2, sigma_phi_h[k])/MU0*(E_theta[i][j][k+1] - E_theta[i][j][k]);

          H_r[i][j][k] = Hr_theta1[area][i][j_area][k_area]
            + Hr_theta2[area][i][j_area][k_area] + Hr_phi[area][i][j_area][k_area];
        }
      }
    }
  }

  //Update Htheta using Htheta_phi, Htheta_r//
  for(int area = 0; area < 4; area++){
    for(int i = 0; i <= Nr - 1; i++){
      ri_1 = dist(i);
      ri_2 = dist(i + 0.5);
      ri_3 = dist(i + 1.0);
      for(int j = idx_theta[area].j1; j <= idx_theta[area].j2; j++){
        int j_area = j - idx_theta[area].j1;
        thetaj_1 = th(j);
        for(int k = idx_theta[area].k1; k <= idx_theta[area].k2; k++){
          int k_area = k - idx_theta[area].k1;
          
          Htheta_phi[area][i][j_area][k_area] = C_1(sigma_phi_h[k])*Htheta_phi[area][i][j_area][k_area]
          - C_4(ri_2, thetaj_1, sigma_phi_h[k])/MU0*(E_r[i][j][k+1] - E_r[i][j][k]);

          Htheta_r[area][i][j_area][k_area] = Htheta_r[area][i][j_area][k_area]
          + C_5(ri_2)/MU0*(ri_3*E_phi[i+1][j][k] - ri_1*E_phi[i][j][k]);

          H_theta[i][j][k] = Htheta_phi[area][i][j_area][k_area] + Htheta_r[area][i][j_area][k_area];
        }
      }
    }
  }
  
  //Update Hphi using Hphi_r, Hphi_theta//
  for(int area = 0; area < 4; area++){
    for(int i = 0; i <= Nr - 1; i++){
      ri_1 = dist(i);
      ri_2 = dist(i + 0.5);
      ri_3 = dist(i + 1.0);
      for(int j = idx_phi[area].j1; j <= idx_phi[area].j2; j++){
        int j_area = j - idx_phi[area].j1;
        for(int k = idx_phi[area].k1; k <= idx_phi[area].k2; k++){
          int k_area = k - idx_phi[area].k1;
          
          Hphi_r[area][i][j_area][k_area] = Hphi_r[area][i][j_area][k_area]
          - C_5(ri_2)/MU0*(ri_3*E_theta[i+1][j][k] - ri_1*E_theta[i][j][k]);

          Hphi_theta[area][i][j_area][k_area] = C_1(sigma_theta_h[j])*Hphi_theta[area][i][j_area][k_area]
          + C_6(ri_2, sigma_theta_h[j])/MU0*(E_r[i][j+1][k] - E_r[i][j][k]);

          H_phi[i][j][k] = Hphi_r[area][i][j_area][k_area] + Hphi_theta[area][i][j_area][k_area];
        }
      }
    }
  }
  
}











