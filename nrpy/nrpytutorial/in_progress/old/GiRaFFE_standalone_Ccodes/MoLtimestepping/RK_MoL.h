// C code implementation of RK4 Method of Lines timestepping.
// ***k1 substep:***

calc_u0(Nxx_plus_2NGHOSTS,aux_gfs);
quantities_to_FD_for_rhs_eval(Nxx_plus_2NGHOSTS,dxx,xx,y_n_gfs,aux_gfs);
rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, y_n_gfs, aux_gfs, k_odd_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = k_odd_gfs[i]*dt*(1.0/6.0);
  k_odd_gfs[i] = y_n_gfs[i] + k_odd_gfs[i]*dt*(1.0/2.0);
}

GiRaFFE_HO_conserv_to_prims_FFE(Nxx, Nxx_plus_2NGHOSTS, dxx,xx, k_odd_gfs, aux_gfs);
apply_bcs(Nxx, Nxx_plus_2NGHOSTS, k_odd_gfs, aux_gfs);
driver_A_to_B(Nxx, Nxx_plus_2NGHOSTS, dxx, k_odd_gfs, aux_gfs);
//apply_bcs_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,k_odd_gfs,aux_gfs);


// ***k2 substep:***

calc_u0(Nxx_plus_2NGHOSTS,aux_gfs);
quantities_to_FD_for_rhs_eval(Nxx_plus_2NGHOSTS,dxx,xx,k_odd_gfs,aux_gfs);
rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, k_odd_gfs, aux_gfs, k_even_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = y_nplus1_running_total_gfs[i] + k_even_gfs[i]*dt*(1.0/3.0);
  k_even_gfs[i] = y_n_gfs[i] + k_even_gfs[i]*dt*(1.0/2.0);
}

GiRaFFE_HO_conserv_to_prims_FFE(Nxx, Nxx_plus_2NGHOSTS, dxx,xx, k_even_gfs, aux_gfs);
apply_bcs(Nxx, Nxx_plus_2NGHOSTS, k_even_gfs, aux_gfs);
driver_A_to_B(Nxx, Nxx_plus_2NGHOSTS, dxx, k_even_gfs, aux_gfs);
//apply_bcs_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,k_even_gfs,aux_gfs);


// ***k3 substep:***

calc_u0(Nxx_plus_2NGHOSTS,aux_gfs);
quantities_to_FD_for_rhs_eval(Nxx_plus_2NGHOSTS,dxx,xx,k_even_gfs,aux_gfs);
rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, k_even_gfs, aux_gfs, k_odd_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_nplus1_running_total_gfs[i] = y_nplus1_running_total_gfs[i] + k_odd_gfs[i]*dt*(1.0/3.0);
  k_odd_gfs[i] = y_n_gfs[i] + k_odd_gfs[i]*dt;
}

GiRaFFE_HO_conserv_to_prims_FFE(Nxx, Nxx_plus_2NGHOSTS, dxx,xx, k_odd_gfs, aux_gfs);
apply_bcs(Nxx, Nxx_plus_2NGHOSTS, k_odd_gfs, aux_gfs);
driver_A_to_B(Nxx, Nxx_plus_2NGHOSTS, dxx, k_odd_gfs, aux_gfs);
//apply_bcs_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,k_odd_gfs,aux_gfs);


// ***k4 substep:***

calc_u0(Nxx_plus_2NGHOSTS,aux_gfs);
quantities_to_FD_for_rhs_eval(Nxx_plus_2NGHOSTS,dxx,xx,k_odd_gfs,aux_gfs);
rhs_eval(Nxx,Nxx_plus_2NGHOSTS,dxx, xx, k_odd_gfs, aux_gfs, k_even_gfs);
LOOP_ALL_GFS_GPS(i) {
  y_n_gfs[i] = y_n_gfs[i] + y_nplus1_running_total_gfs[i] + k_even_gfs[i]*dt*(1.0/6.0);
}

GiRaFFE_HO_conserv_to_prims_FFE(Nxx, Nxx_plus_2NGHOSTS, dxx,xx, y_n_gfs, aux_gfs);
apply_bcs(Nxx, Nxx_plus_2NGHOSTS, y_n_gfs, aux_gfs);
driver_A_to_B(Nxx, Nxx_plus_2NGHOSTS, dxx, y_n_gfs, aux_gfs);
//apply_bcs_EXACT(Nxx,Nxx_plus_2NGHOSTS,xx,n,dt,y_n_gfs,aux_gfs);


