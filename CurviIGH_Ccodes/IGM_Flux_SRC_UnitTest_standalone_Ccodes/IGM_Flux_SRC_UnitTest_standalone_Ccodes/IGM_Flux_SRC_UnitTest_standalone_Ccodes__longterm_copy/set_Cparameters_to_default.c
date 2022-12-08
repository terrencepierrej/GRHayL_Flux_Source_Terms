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
  params->xxmin0 = 0.1;  // grid::xxmin0
  params->xxmin1 = 0.1;  // grid::xxmin1
  params->xxmin2 = 0.1;  // grid::xxmin2
  params->xxmax0 = 0.1;  // grid::xxmax0
  params->xxmax1 = 0.1;  // grid::xxmax1
  params->xxmax2 = 0.1;  // grid::xxmax2
  params->invdx0 = 1.0;  // grid::invdx0
  params->invdx1 = 1.0;  // grid::invdx1
  params->invdx2 = 1.0;  // grid::invdx2
  params->Cart_originx = 0.0;  // grid::Cart_originx
  params->Cart_originy = 0.0;  // grid::Cart_originy
  params->Cart_originz = 0.0;  // grid::Cart_originz
  params->Cart_CoM_offsetx = 0.0;  // grid::Cart_CoM_offsetx
  params->Cart_CoM_offsety = 0.0;  // grid::Cart_CoM_offsety
  params->Cart_CoM_offsetz = 0.0;  // grid::Cart_CoM_offsetz
  params->f0_of_xx0 = 1e+300;  // reference_metric::f0_of_xx0
  params->f1_of_xx1 = 1e+300;  // reference_metric::f1_of_xx1
  params->f2_of_xx1 = 1e+300;  // reference_metric::f2_of_xx1
  params->f2_of_xx0_xx1 = 1e+300;  // reference_metric::f2_of_xx0_xx1
  params->f3_of_xx0 = 1e+300;  // reference_metric::f3_of_xx0
  params->f4_of_xx2 = 1e+300;  // reference_metric::f4_of_xx2
  params->f0_of_xx0__D0 = 1e+300;  // reference_metric::f0_of_xx0__D0
  params->f0_of_xx0__DD00 = 1e+300;  // reference_metric::f0_of_xx0__DD00
  params->f0_of_xx0__DDD000 = 1e+300;  // reference_metric::f0_of_xx0__DDD000
  params->f1_of_xx1__D1 = 1e+300;  // reference_metric::f1_of_xx1__D1
  params->f1_of_xx1__DD11 = 1e+300;  // reference_metric::f1_of_xx1__DD11
  params->f1_of_xx1__DDD111 = 1e+300;  // reference_metric::f1_of_xx1__DDD111
  params->f2_of_xx1__D1 = 1e+300;  // reference_metric::f2_of_xx1__D1
  params->f2_of_xx1__DD11 = 1e+300;  // reference_metric::f2_of_xx1__DD11
  params->f2_of_xx1__DDD111 = 1e+300;  // reference_metric::f2_of_xx1__DDD111
  params->f2_of_xx0_xx1__D0 = 1e+300;  // reference_metric::f2_of_xx0_xx1__D0
  params->f2_of_xx0_xx1__D1 = 1e+300;  // reference_metric::f2_of_xx0_xx1__D1
  params->f2_of_xx0_xx1__DD00 = 1e+300;  // reference_metric::f2_of_xx0_xx1__DD00
  params->f2_of_xx0_xx1__DD11 = 1e+300;  // reference_metric::f2_of_xx0_xx1__DD11
  params->f3_of_xx0__D0 = 1e+300;  // reference_metric::f3_of_xx0__D0
  params->f3_of_xx0__DD00 = 1e+300;  // reference_metric::f3_of_xx0__DD00
  params->f4_of_xx2__D2 = 1e+300;  // reference_metric::f4_of_xx2__D2
  params->f4_of_xx2__DD22 = 1e+300;  // reference_metric::f4_of_xx2__DD22
  params->xmin = -10.0;  // reference_metric::xmin
  params->xmax = 10.0;  // reference_metric::xmax
  params->ymin = -10.0;  // reference_metric::ymin
  params->ymax = 10.0;  // reference_metric::ymax
  params->zmin = -10.0;  // reference_metric::zmin
  params->zmax = 10.0;  // reference_metric::zmax
  params->GAMMA_SPEED_LIMIT = 10.0;  // GRMHD_equations_new_version::GAMMA_SPEED_LIMIT
  params->sqrt4pi = sqrt(4.0*M_PI);  // GRMHD_equations_new_version::sqrt4pi
  params->has_outer_boundary = 0;  // CurviBoundaryConditions.CurviBoundaryConditions_new_way::has_outer_boundary
}
