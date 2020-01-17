#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>

#include <mpi.h>
#include "fdtd3d.h"

//The number of R, Theta, Phi element//
const int Nr = 100;
const int Ntheta = 100;
const int Nphi = 200;

//Minute R, Theta, Phi, Time//
const double delta_r = 100.0e3/double(Nr);
const double delta_theta = 1.0e3/double(R0);
const double delta_phi = 1.0e3/double(R0);
const double Dt = 0.99/C0/std::sqrt(std::pow(1.0/delta_r, 2.0)
 + std::pow(1.0/R0/delta_theta, 2.0) + std::pow(1.0/R0/std::sin(THETA0)/delta_phi, 2.0));
const double inv_Dt = 1.0/Dt;

//PML info//
const int L = 10;
const double M = 3.5;
const double R = 1.0e-6;
const double sigma_th_max = -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_theta/R0;
const double sigma_phi_max = -(M + 1.0)*C0*std::log(R)/2.0/double(L)/delta_phi/R0;

int main(int argc, char** argv)
{
  int time_step = 300;
  double sigma_t = 7.0*Dt;
  double t0 = 6.0*sigma_t;
  double t;
  double I;
  double time_1, time_2, total_time;
  int NEW;
  int OLD;
  
  //center cordinate//
  int i0 = Nr/2;
  int j0 = Ntheta/2;
  int k0 = Nphi/2;
  
  double*** Hr, ***Htheta, ***Hphi;
  Hr = memory_allocate3d(Nr + 1, Ntheta, Nphi);
  Htheta = memory_allocate3d(Nr, Ntheta + 1, Nphi);
  Hphi = memory_allocate3d(Nr, Ntheta, Nphi + 1);
  
  double*** Er, *** Etheta, *** Ephi;
  Er = memory_allocate3d(Nr, Ntheta + 1, Nphi + 1);
  Etheta = memory_allocate3d(Nr + 1, Ntheta, Nphi + 1);
  Ephi = memory_allocate3d(Nr + 1, Ntheta + 1, Nphi);
  
  double**** Dr, **** Dtheta, **** Dphi;
  Dr = memory_allocate4d(2, Nr, Ntheta + 1, Nphi + 1);
  Dtheta = memory_allocate4d(2, Nr + 1, Ntheta, Nphi);
  Dphi = memory_allocate4d(2, Nr + 1, Ntheta + 1, Nphi);
  
  double**** Dr_theta1, **** Dr_theta2, **** Dr_phi;
  double**** Dtheta_phi, **** Dtheta_r; 
  double**** Dphi_r, **** Dphi_theta;
  
  Dr_theta1 = new double***[4];
  Dr_theta2 = new double***[4];
  Dr_phi = new double***[4];
  Dtheta_phi = new double***[4];
  Dtheta_r = new double***[4];
  Dphi_r = new double***[4];
  Dphi_theta = new double***[4];
  
  double**** Hr_theta1, ****Hr_theta2, ****Hr_phi;
  double**** Htheta_phi, ****Htheta_r;
  double**** Hphi_r, ****Hphi_theta;

  Hr_theta1 = new double***[4];
  Hr_theta2 = new double***[4];
  Hr_phi = new double***[4];
  Htheta_phi = new double***[4];
  Htheta_r = new double***[4];
  Hphi_r = new double***[4];
  Hphi_theta = new double***[4];

  //PML region (Theta direction)//
  for(int i = 0; i <= 1; i++){
    //D components in PML(Theta direction)//
    Dr_theta1[i] = memory_allocate3d(Nr, L, Nphi - 1);
    Dr_theta2[i] = memory_allocate3d(Nr, L, Nphi - 1);
    Dr_phi[i] = memory_allocate3d(Nr, L, Nphi - 1);
    Dtheta_phi[i] = memory_allocate3d(Nr, L, Nphi - 1);
    Dtheta_r[i] = memory_allocate3d(Nr, L, Nphi - 1);
    Dphi_r[i] = memory_allocate3d(Nr, L, Nphi);
    Dphi_theta[i] = memory_allocate3d(Nr, L, Nphi);

    //H compornents in PML(Theta direction)//
    Hr_theta1[i] = memory_allocate3d(Nr + 1, L, Nphi);
    Hr_theta2[i] = memory_allocate3d(Nr + 1, L, Nphi);
    Hr_phi[i] = memory_allocate3d(Nr + 1, L, Nphi);
    Htheta_phi[i] = memory_allocate3d(Nr, L + 1, Nphi);
    Htheta_r[i] = memory_allocate3d(Nr, L + 1, Nphi);
    Hphi_r[i] = memory_allocate3d(Nr, L, Nphi + 1);
    Hphi_theta[i] = memory_allocate3d(Nr, L, Nphi + 1);
  }

  //PML region (Phi direction)//
  for(int i = 2; i <= 3; i++){
    //D components in PML(Phi direction)//
    Dr_theta1[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);
    Dr_theta2[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);
    Dr_phi[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);
    Dtheta_phi[i] = memory_allocate3d(Nr, Ntheta - 2*L, Nphi - 1);
    Dtheta_r[i] = memory_allocate3d(Nr, Ntheta - 2*L, Nphi - 1);
    Dphi_r[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);
    Dphi_theta[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);

    //H components in PML(Phi direction)//
    Hr_theta1[i] = memory_allocate3d(Nr + 1, Ntheta - 2*L, L);
    Hr_theta2[i] = memory_allocate3d(Nr + 1, Ntheta - 2*L, L);
    Hr_phi[i] = memory_allocate3d(Nr + 1, Ntheta - 2*L, L);
    Htheta_phi[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);
    Htheta_r[i] = memory_allocate3d(Nr, Ntheta - 2*L - 1, L);
    Hphi_r[i] = memory_allocate3d(Nr, Ntheta - 2*L, L + 1);
    Hphi_theta[i] = memory_allocate3d(Nr, Ntheta - 2*L, L + 1);
  }
 
  double *sigma_theta, *sigma_phi, *sigma_theta_h, *sigma_phi_h;
  sigma_theta = new double[Ntheta + 1];
  sigma_phi = new double[Nphi + 1];
  sigma_theta_h = new double[Ntheta + 1];
  sigma_phi_h = new double[Nphi + 1];

  sigma_calc(sigma_theta, sigma_phi, sigma_theta_h, sigma_phi_h);
  
  std::ofstream ofs("E0.dat");
  std::ofstream ofs_2("graph.gpt");
  
  for(int i = 0; i < Nr; i++){
    double r = dist(i);
    for(int j = 0; j < Ntheta + 1; j++){
      double theta = th(j);
      ofs << r << " " << theta << " " << Er[i][j][k0];
    }
    ofs << std::endl;
  } 
  
  ofs.close();
  
  ofs_2 << 0 << " " << Er[i0][25][50] << std::endl;;

  MPI::Init(argc, argv);

  int myid = MPI::COMM_WORLD.Get_rank();
  int nprocs = MPI::COMM_WORLD.Get_size();

  if(myid == 0)std::cout << nprocs << " processers. " << std::endl;

  int band = Nphi/nprocs;
  int surplus = Nphi%nprocs;

  //各プロセス計算領域(E, H)の端番号//
  int peDphi_l = myid * band + myid;
  int 
  int pe_r = peE_l + bamd;
  
  ////////計測開始////////
  std::chrono..system_clock::time_point start
    = std:::chrono::system_clock::now();

  //FDTD_update//
  for(int n = 1; n < time_step + 1; n++){
    
    NEW = n%2;
    OLD = (n+1)%2;
    
    t = (double(n) - 0.5)*Dt;
    
    //Forced current//
    I = -((t - t0)/sigma_t/sigma_t)*std::exp(-(t - t0)*(t - t0)/2.0/sigma_t/sigma_t);
    std::cout << " I = " << I << std::endl;
    
    Er[i0][j0][k0] = Er[i0][j0][k0] - (Dt / EPS0)*I/delta_r/delta_theta/delta_phi;
    
    std::cout << "Er[" << i0 << "][" << j0 << "][" << k0 << "] = " << Er[i0][j0][k0] << std::endl;
    
    //update D (outside PML)//
    D_update(Dr[NEW], Dtheta[NEW], Dphi[NEW], Hr, Htheta, Hphi, Dr[OLD], Dtheta[OLD], Dphi[OLD]);
    
    //update D (inside PML)//
    D_update_pml(Dr[NEW], Dtheta[NEW], Dphi[NEW], Hr, Htheta, Hphi, 
    Dr_theta1, Dr_theta2, Dr_phi, Dtheta_phi, Dtheta_r, Dphi_r, Dphi_theta, 
    sigma_theta, sigma_phi);
    
    //update E//
    E_update(Er, Etheta, Ephi, Dr[NEW], Dtheta[NEW], Dphi[NEW], Dr[OLD], Dtheta[OLD], Dphi[OLD]);
    
    //update H (outside PML)//
    H_update(Er, Etheta, Ephi, Hr, Htheta, Hphi);
    
    //updata H (inside PML)//
    H_update_pml(Er, Etheta, Ephi, Hr, Htheta, Hphi, 
    Hr_theta1, Hr_theta2, Hr_phi, Htheta_phi, Htheta_r, Hphi_r, Hphi_theta, 
    sigma_theta_h, sigma_phi_h);
    
    std::string fn = "E" + std::to_string(n) + ".dat";
    std::ofstream ofs(fn.c_str());
    
    for(int i = 0; i < Nr; i++){
      double r = dist(i);
      for(int j = 0; j < Ntheta + 1; j++){
        double theta = th(j);
        ofs << r << " " << theta << " " << Er[i][j][k0] << std::endl;
      }
      ofs << std::endl;
    }
    
    ofs.close();

    ofs_2 << n << " " << Er[i0][25][50] << std::endl;
    
    std::cout << n << " / " << time_step << std::endl << std::endl;
    
  }

  std::chrono::system_clock::time_point end
    = std::chrono::system_clock::now();
  ///////計測終了///////

  total_time = std::chrono::duration_cast <std::chrono::milliseconds>
  
  std::cout << "elapsed_time = " << total_time*1e-3 << " [sec]"<< std::endl;

  ofs_2.close();
  
  return 0;
  
}



