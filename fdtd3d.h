//Physical quantity//
#define C0 (3.0e8)
#define MU0 (4.0*M_PI*1.0e-7)
#define EPS0 (1.0/MU0/C0/C0)
#define R0 (6370e3)
#define THETA0 (M_PI*0.5 - std::atan(50e3/R0))

//The number of R, Theta, Phi element//
extern const int Nr;
extern const int Ntheta;
extern const int Nphi;

//Minute R, Theta, Phi, Time//
extern const double delta_r;
extern const double delta_theta;
extern const double delta_phi;
extern const double Dt;
extern const double inv_Dt;

//PML information//
extern const int L;
extern const double M;
extern const double R;
extern const double sigma_th_max;
extern const double sigma_phi_max;

class pml{
 public:
  int j1, j2, k1, k2;
  void set_point_1(int, int);
  void set_point_2(int, int);
};

//function in main//
double*** memory_allocate3d(int, int, int);
double**** memory_allocate4d(int, int, int, int);

void sigma_calc(double* sigma_theta, double* sigma_phi, 
              double* sigma_theta_h, double* sigma_phi_h);

void D_update(double*** Dr_at_nDt, double*** Dtheta_at_nDt, double*** Dphi_at_nDt,
              double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_n_minus_halfDt,
              double*** Dr_at_n_minus_oneDt, double*** Dtheta_at_n_minus_halfDt, double*** Dphi_at_n_minus_oneDt);

void D_update_pml(double*** Dr_at_nDt, double*** Dtheta_at_nDt, double*** Dphi_at_nDt,
		  double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_n_minus_halfDt,
		  double**** Dr_theta1_n_minus_oneDt, double**** Dr_theta2_n_minus_halfDt, double**** Dr_phi_n_minus_oneDt,
		  double**** Dtheta_phi_at_n_minus_oneDt, double**** Dtheta_r_at_n_minus_oneDt,
		  double**** Dphi_r_at_n_minus_oneDt, double**** Dphi_theta_n_minus_oneDt,
                  double* sigma_theta, double* sigma_phi);

void E_update(double*** Er_at_n_minus_oneDt, double*** Etheta_at_n_minus_oneDt, double*** Ephi_at_n_minus_oneDt,
              double*** Dr_at_nDt, double*** Dtheta_at_n_Dt, double*** Dphi_at_nDt,
              double*** Dr_at_n_minus_oneDt, double*** Dtheta_at_nDt, double*** Dphi_at_n_minus_oneDt);

void H_update(double*** Er_at_nDt, double*** Etheta_at_nDt, double*** Ephi_at_nDt,
              double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_at_n_minus_halfDt);

void H_update_pml(double*** Er_at_nDt, double*** Etheta_at_nDt, double*** Ephi_at_n_Dt,
		  double*** Hr_at_n_minus_halfDt, double*** Htheta_at_n_minus_halfDt, double*** Hphi_at_n_minus_halfDt,
		  double**** Hr_theta1_at_n_minus_halfDt, double**** Hr_theta2_at_n_minus_halfDt, double**** Hr_phi_at_n_minus_halfDt,
		  double**** Htheta_phi_at_n_minus_halfDt, double**** Htheta_r_at_n_minus_halfDt,
		  double**** Hphi_r_at_n_minus_halfDt, double**** Hphi_theta_at_n_minus_halfDt,
		  double* sigma_theta_half, double* sigma_phi_half);

//inline function//
inline double dist(double i){return R0 + i*delta_r;};
inline double th(double j){return THETA0 + j*delta_theta;};

inline double C_1(double sig){return ((inv_Dt - sig/2.0)/(inv_Dt + sig/2.0));};
inline double C_2(double r, double sig){return 1.0/r/delta_theta/(inv_Dt + sig/2.0);};
inline double C_3(double r, double theta){return Dt*std::cos(theta)/std::sin(theta)/2.0/r;};
inline double C_4(double r, double theta, double sig){return 1.0/r/std::sin(theta)/delta_phi/(inv_Dt + sig/2.0);};
inline double C_5(double r){return Dt/r/delta_r;};
inline double C_6(double r, double sig){return 1.0/(inv_Dt + sig/2.0)/r/delta_theta;};



