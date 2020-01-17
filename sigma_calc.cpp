#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include "fdtd3d.h"

void sigma_calc(double* sigma_theta, double* sigma_phi, double* sigma_theta_h, double* sigma_phi_h)
{
  for(int i = 0; i <= Ntheta; i++){
    
    double theta = i*delta_theta;
    double theta_h = (i+0.5)*delta_theta;

    if(i <= L){
      sigma_theta[i] = sigma_th_max*std::pow((L*delta_theta - theta)/double(L)/delta_theta, M);
      sigma_theta_h[i] = sigma_th_max*std::pow((L*delta_theta - theta_h)/double(L)/delta_theta, M);
    } 
    else if(i >= Ntheta - L){
      sigma_theta[i] = sigma_th_max*std::pow((theta - (Ntheta - L)*delta_theta)/double(L)/delta_theta, M);
      sigma_theta_h[i] = sigma_th_max*std::pow((theta_h - (Ntheta - L)*delta_theta)/double(L)/delta_theta, M);
    } 
    else {
      sigma_theta[i] = 0.0;
      sigma_theta_h[i] = 0.0;
    }
  }

  for(int i = 0; i <= Nphi; i++){
    
    double phi = i*delta_phi;
    double phi_h = (i+0.5)*delta_phi;

    if(i <= L){
      sigma_phi[i] = sigma_phi_max*std::pow((L*delta_phi - phi)/double(L)/delta_phi, M);
      sigma_phi_h[i] = sigma_phi_max*std::pow((L*delta_phi - phi_h)/double(L)/delta_phi, M);
    } 
    else if(i >= Nphi - L) {
      sigma_phi[i] = sigma_phi_max*std::pow((phi - (Nphi - L)*delta_phi)/double(L)/delta_phi, M);
      sigma_phi_h[i] = sigma_phi_max*std::pow((phi_h - (Nphi - L)*delta_phi)/double(L)/delta_phi, M);
    }
    else {
      sigma_phi[i] = 0.0;
      sigma_phi_h[i] = 0.0;
    }
  }

  sigma_theta_h[L] = 0.0;
  sigma_phi_h[L] = 0.0;
}


