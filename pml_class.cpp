#include <iostream>
#include <cmath>
#include "fdtd3d.h"

//PML class//
void pml::set_point_1(int py, int pz){
  j1 = py;
  k1 = pz;
}

void pml::set_point_2(int py, int pz){
  j2 = py;
  k2 = pz;
}
