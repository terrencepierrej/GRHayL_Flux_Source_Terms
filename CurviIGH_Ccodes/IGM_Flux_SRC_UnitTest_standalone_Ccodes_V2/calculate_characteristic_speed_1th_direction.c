#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Compute the characteristic speeds in1th direction
 */
void calculate_characteristic_speed_1th_direction(const reconstructed_prims_struct *restrict reconstructed_prims_r, const reconstructed_prims_struct *restrict reconstructed_prims_l,const metric_face_quantities_struct *restrict metric_face_quantities, conservative_fluxes_struct *restrict conservative_fluxes) {

{
const double u4_rU0 = reconstructed_prims_r->u4U0;
const double u4_rU1 = reconstructed_prims_r->u4U1;
const double u4_rU2 = reconstructed_prims_r->u4U2;
const double u4_rU3 = reconstructed_prims_r->u4U3;
const double B_rU0 = reconstructed_prims_r->BU0;
const double B_rU1 = reconstructed_prims_r->BU1;
const double B_rU2 = reconstructed_prims_r->BU2;
const double u4_lU0 = reconstructed_prims_l->u4U0;
const double u4_lU1 = reconstructed_prims_l->u4U1;
const double u4_lU2 = reconstructed_prims_l->u4U2;
const double u4_lU3 = reconstructed_prims_l->u4U3;
const double B_lU0 = reconstructed_prims_l->BU0;
const double B_lU1 = reconstructed_prims_l->BU1;
const double B_lU2 = reconstructed_prims_l->BU2;
const double P_r = reconstructed_prims_r->P;
const double P_l = reconstructed_prims_l->P;
const double h_r = reconstructed_prims_r->h;
const double h_l = reconstructed_prims_l->h;
const double rhob_r = reconstructed_prims_r->rhob;
const double rhob_l = reconstructed_prims_l->rhob;
const double Gamma_th_r = reconstructed_prims_r->Gamma_th;
const double Gamma_th_l = reconstructed_prims_l->Gamma_th;
const double epsilon_th_r = reconstructed_prims_r->epsilon_th;
const double epsilon_th_l = reconstructed_prims_l->epsilon_th;
const double dPcold_drhob_r = reconstructed_prims_r->dPcold_drhob;
const double dPcold_drhob_l = reconstructed_prims_l->dPcold_drhob;
const double alpha_face = metric_face_quantities->alpha_face;
const double beta_faceU0 = metric_face_quantities->beta_faceU0;
const double beta_faceU1 = metric_face_quantities->beta_faceU1;
const double beta_faceU2 = metric_face_quantities->beta_faceU2;
const double gamma_faceDD00 = metric_face_quantities->gamma_faceDD00;
const double gamma_faceDD01 = metric_face_quantities->gamma_faceDD01;
const double gamma_faceDD02 = metric_face_quantities->gamma_faceDD02;
const double gamma_faceDD11 = metric_face_quantities->gamma_faceDD11;
const double gamma_faceDD12 = metric_face_quantities->gamma_faceDD12;
const double gamma_faceDD22 = metric_face_quantities->gamma_faceDD22;
  const double tmp_1 = (1.0/((alpha_face)*(alpha_face)));
  const double tmp_2 = beta_faceU0*gamma_faceDD00 + beta_faceU1*gamma_faceDD01 + beta_faceU2*gamma_faceDD02;
  const double tmp_3 = beta_faceU0*gamma_faceDD01 + beta_faceU1*gamma_faceDD11 + beta_faceU2*gamma_faceDD12;
  const double tmp_4 = beta_faceU0*gamma_faceDD02 + beta_faceU1*gamma_faceDD12 + beta_faceU2*gamma_faceDD22;
  const double tmp_5 = B_lU0*(gamma_faceDD00*u4_lU1 + gamma_faceDD01*u4_lU2 + gamma_faceDD02*u4_lU3 + tmp_2*u4_lU0) + B_lU1*(gamma_faceDD01*u4_lU1 + gamma_faceDD11*u4_lU2 + gamma_faceDD12*u4_lU3 + tmp_3*u4_lU0) + B_lU2*(gamma_faceDD02*u4_lU1 + gamma_faceDD12*u4_lU2 + gamma_faceDD22*u4_lU3 + tmp_4*u4_lU0);
  const double tmp_6 = B_lU0 + tmp_5*u4_lU1;
  const double tmp_8 = tmp_1/((sqrt4pi)*(sqrt4pi));
  const double tmp_9 = tmp_8/((u4_lU0)*(u4_lU0));
  const double tmp_10 = B_lU1 + tmp_5*u4_lU2;
  const double tmp_11 = B_lU2 + tmp_5*u4_lU3;
  const double tmp_12 = tmp_8*(-((alpha_face)*(alpha_face)) + beta_faceU0*tmp_2 + beta_faceU1*tmp_3 + beta_faceU2*tmp_4);
  const double tmp_16 = 2*tmp_8;
  const double tmp_18 = tmp_5/u4_lU0;
  const double tmp_20 = gamma_faceDD00*((tmp_6)*(tmp_6))*tmp_9 + 2*gamma_faceDD01*tmp_10*tmp_6*tmp_9 + 2*gamma_faceDD02*tmp_11*tmp_6*tmp_9 + gamma_faceDD11*((tmp_10)*(tmp_10))*tmp_9 + 2*gamma_faceDD12*tmp_10*tmp_11*tmp_9 + gamma_faceDD22*((tmp_11)*(tmp_11))*tmp_9 + tmp_10*tmp_16*tmp_18*tmp_3 + tmp_11*tmp_16*tmp_18*tmp_4 + tmp_12*((tmp_5)*(tmp_5)) + tmp_16*tmp_18*tmp_2*tmp_6;
  const double tmp_21 = tmp_20/(h_l*rhob_l + tmp_20);
  const double tmp_23 = -(tmp_21 - 1)*(Gamma_th_l*epsilon_th_l*(Gamma_th_l - 1) + dPcold_drhob_l)/h_l;
  const double tmp_25 = tmp_1*(tmp_21 + tmp_23);
  const double tmp_26 = -tmp_21 - tmp_23 + 1;
  const double tmp_27 = tmp_26*((u4_lU0)*(u4_lU0));
  const double tmp_28 = (1.0/(tmp_25 + tmp_27));
  const double tmp_29 = -beta_faceU1*tmp_1*(2*tmp_21 + 2*tmp_23) + 2*tmp_26*u4_lU0*u4_lU2;
  const double tmp_30 = B_rU0*(gamma_faceDD00*u4_rU1 + gamma_faceDD01*u4_rU2 + gamma_faceDD02*u4_rU3 + tmp_2*u4_rU0) + B_rU1*(gamma_faceDD01*u4_rU1 + gamma_faceDD11*u4_rU2 + gamma_faceDD12*u4_rU3 + tmp_3*u4_rU0) + B_rU2*(gamma_faceDD02*u4_rU1 + gamma_faceDD12*u4_rU2 + gamma_faceDD22*u4_rU3 + tmp_4*u4_rU0);
  const double tmp_31 = B_rU0 + tmp_30*u4_rU1;
  const double tmp_33 = tmp_8/((u4_rU0)*(u4_rU0));
  const double tmp_34 = B_rU1 + tmp_30*u4_rU2;
  const double tmp_35 = B_rU2 + tmp_30*u4_rU3;
  const double tmp_37 = 2*tmp_35;
  const double tmp_38 = tmp_30/u4_rU0;
  const double tmp_39 = gamma_faceDD00*((tmp_31)*(tmp_31))*tmp_33 + 2*gamma_faceDD01*tmp_31*tmp_33*tmp_34 + gamma_faceDD02*tmp_31*tmp_33*tmp_37 + gamma_faceDD11*tmp_33*((tmp_34)*(tmp_34)) + gamma_faceDD12*tmp_33*tmp_34*tmp_37 + gamma_faceDD22*tmp_33*((tmp_35)*(tmp_35)) + tmp_12*((tmp_30)*(tmp_30)) + tmp_16*tmp_2*tmp_31*tmp_38 + tmp_16*tmp_3*tmp_34*tmp_38 + tmp_37*tmp_38*tmp_4*tmp_8;
  const double tmp_40 = tmp_39/(h_r*rhob_r + tmp_39);
  const double tmp_42 = -(tmp_40 - 1)*(Gamma_th_r*epsilon_th_r*(Gamma_th_r - 1) + dPcold_drhob_r)/h_r;
  const double tmp_44 = tmp_1*(tmp_40 + tmp_42);
  const double tmp_45 = -tmp_40 - tmp_42 + 1;
  const double tmp_46 = tmp_45*((u4_rU0)*(u4_rU0));
  const double tmp_47 = (1.0/(tmp_44 + tmp_46));
  const double tmp_48 = -beta_faceU1*tmp_1*(2*tmp_40 + 2*tmp_42) + 2*tmp_45*u4_rU0*u4_rU2;
  const double tmp_52 = -1.0/4.0*tmp_28*tmp_29;
  const double tmp_54 = (1.0/4.0)*tmp_47*tmp_48;
  const double tmp_55 = (1.0/2.0)*tmp_28;
  const double tmp_58 = -((beta_faceU1)*(beta_faceU1))*tmp_1 + (gamma_faceDD00*gamma_faceDD22 - ((gamma_faceDD02)*(gamma_faceDD02)))/(gamma_faceDD00*gamma_faceDD11*gamma_faceDD22 - gamma_faceDD00*((gamma_faceDD12)*(gamma_faceDD12)) - ((gamma_faceDD01)*(gamma_faceDD01))*gamma_faceDD22 + 2*gamma_faceDD01*gamma_faceDD02*gamma_faceDD12 - ((gamma_faceDD02)*(gamma_faceDD02))*gamma_faceDD11);
  const double tmp_59 = (4*tmp_25 + 4*tmp_27)*(tmp_26*((u4_lU2)*(u4_lU2)) - tmp_58*(tmp_21 + tmp_23));
  const double tmp_60 = fabs(tmp_28*sqrt((1.0/2.0)*((tmp_29)*(tmp_29)) - 1.0/2.0*tmp_59 + (1.0/2.0)*fabs(((tmp_29)*(tmp_29)) - tmp_59)));
  const double tmp_61 = (1.0/2.0)*tmp_47;
  const double tmp_63 = (4*tmp_44 + 4*tmp_46)*(tmp_45*((u4_rU2)*(u4_rU2)) - tmp_58*(tmp_40 + tmp_42));
  const double tmp_64 = fabs(tmp_47*sqrt((1.0/2.0)*((tmp_48)*(tmp_48)) - 1.0/2.0*tmp_63 + (1.0/2.0)*fabs(((tmp_48)*(tmp_48)) - tmp_63)));
  const double tmp_65 = (1.0/2.0)*tmp_60 - 1.0/2.0*tmp_64;
  const double tmp_66 = fabs((1.0/4.0)*tmp_28*tmp_29 - 1.0/4.0*tmp_47*tmp_48 - tmp_52 - tmp_54 - tmp_65);
  const double tmp_68 = -1.0/8.0*tmp_28*tmp_29;
  const double tmp_70 = -1.0/8.0*tmp_47*tmp_48;
  const double tmp_71 = (1.0/4.0)*tmp_60 + (1.0/4.0)*tmp_64;
  const double tmp_73 = (1.0/16.0)*tmp_28*tmp_29;
  const double tmp_75 = (1.0/16.0)*tmp_47*tmp_48;
  const double tmp_76 = -1.0/16.0*tmp_28*tmp_29;
  const double tmp_77 = -1.0/16.0*tmp_47*tmp_48;
  const double tmp_78 = (1.0/8.0)*tmp_60 + (1.0/8.0)*tmp_64;
  const double tmp_79 = fabs((1.0/4.0)*tmp_28*tmp_29 - 1.0/4.0*tmp_47*tmp_48 - tmp_52 - tmp_54 + tmp_65);
  conservative_fluxes->cmin_dirn1 = (1.0/4.0)*tmp_66 - tmp_73 - tmp_75 + tmp_76 + tmp_77 + tmp_78 + (1.0/2.0)*fabs((1.0/8.0)*tmp_28*tmp_29 + (1.0/8.0)*tmp_47*tmp_48 - 1.0/2.0*tmp_66 - tmp_68 - tmp_70 - tmp_71);
  conservative_fluxes->cmax_dirn1 = tmp_73 + tmp_75 - tmp_76 - tmp_77 + tmp_78 + (1.0/4.0)*tmp_79 + (1.0/2.0)*fabs((1.0/8.0)*tmp_28*tmp_29 + (1.0/8.0)*tmp_47*tmp_48 - tmp_68 - tmp_70 + tmp_71 + (1.0/2.0)*tmp_79);
}
}
