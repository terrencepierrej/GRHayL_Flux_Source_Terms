#include "./NRPy_basic_defines.h"
/*
 * Set Cparameters to default values specified within NRPy+.
 */
void set_Cparameters_to_default(paramstruct *restrict params) {

  params->Nxx0 = 64;  // grid::Nxx0
  params->Nxx1 = 32;  // grid::Nxx1
  params->Nxx2 = 64;  // grid::Nxx2
  params->Nxx_plus_2NGHOSTS0 = 70;  // grid::Nxx_plus_2NGHOSTS0
  params->Nxx_plus_2NGHOSTS1 = 38;  // grid::Nxx_plus_2NGHOSTS1
  params->Nxx_plus_2NGHOSTS2 = 70;  // grid::Nxx_plus_2NGHOSTS2
  params->xx0 = 1e+300;  // grid::xx0
  params->xx1 = 1e+300;  // grid::xx1
  params->xx2 = 1e+300;  // grid::xx2
  params->dxx0 = 0.1;  // grid::dxx0
  params->dxx1 = 0.1;  // grid::dxx1
  params->dxx2 = 0.1;  // grid::dxx2
  params->invdx0 = 1.0;  // grid::invdx0
  params->invdx1 = 1.0;  // grid::invdx1
  params->invdx2 = 1.0;  // grid::invdx2
  params->xmin = -10.0;  // reference_metric::xmin
  params->xmax = 10.0;  // reference_metric::xmax
  params->ymin = -10.0;  // reference_metric::ymin
  params->ymax = 10.0;  // reference_metric::ymax
  params->zmin = -10.0;  // reference_metric::zmin
  params->zmax = 10.0;  // reference_metric::zmax
  params->GAMMA_SPEED_LIMIT = 10.0;  // GRMHD_equations_new_version::GAMMA_SPEED_LIMIT
  params->sqrt4pi = sqrt(4.0*M_PI);  // GRMHD_equations_new_version::sqrt4pi
}
