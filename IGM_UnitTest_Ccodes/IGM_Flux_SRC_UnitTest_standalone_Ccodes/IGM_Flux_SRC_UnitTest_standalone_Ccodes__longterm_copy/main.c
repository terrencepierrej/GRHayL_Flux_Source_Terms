#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include "inlined_functions.C"
// #include "reconstruct_set_of_prims_PPM.C"
// #include "loop_defines_reconstruction.h"
#include "add_fluxes_and_source_terms_to_hydro_rhss.C"
#include "compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu.C"
#include "mhdflux.C"
/*
 * // main() function:
 * // Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
 * // Step 1: Write test data to gridfunctions
 * // Step 2: Overwrite all data in ghost zones with NaNs
 * // Step 3: Apply curvilinear boundary conditions
 * // Step 4: Print gridfunction data after curvilinear boundary conditions have been applied
 * // Step 5: Free all allocated memory
 */
int main(int argc, const char *argv[]) {

  griddata_struct griddata;
  set_Cparameters_to_default(&griddata.params);

  // Step 0a: Read command-line input, error out if nonconformant
  /* if(argc != 5 || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < NGHOSTS) {
    printf("Error: Expected one command-line argument: ./ScalarWaveCurvilinear_Playground Nx0 Nx1 Nx2,\n");
    printf("where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
    printf("Nx[] MUST BE larger than NGHOSTS (= %d)\n",NGHOSTS);
    exit(1);
  } */
  
  // Step 0b: Set up numerical grid structure, first in space...
  const int Nxx[3] = { atoi(argv[1]), atoi(argv[2]), atoi(argv[3]) };
  if(Nxx[0]%2 != 0 || Nxx[1]%2 != 0 || Nxx[2]%2 != 0) {
    printf("Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
    printf("       For example, in case of angular directions, proper symmetry zones will not exist.\n");
    exit(1);
  }

  // Step 0c: Set free parameters, overwriting Cparameters defaults
  //          BU1 hand or with command-line input, as desired.
#include "free_parameters.h"

  {
    int EigenCoord = 0;
    // Step 0d.ii: Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    //             params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    //             chosen Eigen-CoordSystem.
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &griddata.params, griddata.xx);
    // Step 0e: Find ghostzone mappings; set up bcstruct
    // Step 0e.i: Free allocated space for xx[][] array
  }
  
  // Step 0.m: Allocate memory for non y_n_gfs. We do this here to free up
  //         memory for setting up initial data (for cases in which initial
  //         data setup is memory intensive.)
  const int Nxx_plus_2NGHOSTS0 = griddata.params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata.params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata.params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_tot = griddata.params.Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  double *restrict evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  double *restrict auxevol_gfs = (double *restrict)malloc(sizeof(double) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);


  // Now we need to properly populate our structs for the unit test.

  reconstructed_prims_struct reconstructed_prims;
  reconstructed_prims_struct reconstructed_prims_r;
  reconstructed_prims_struct reconstructed_prims_l;

  conservative_fluxes_struct conservative_fluxes;
  conservative_sources_struct conservative_sources;

  metric_quantities_struct metric_quantities;
  metric_face_quantities_struct metric_face_quantities;
  metric_quantities_derivatives_struct metric_quantities_derivatives;

  // We begin BU1 generating random inital data.

  for(int which_gf=0; which_gf<NUM_AUXEVOL_GFS; which_gf++) for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
    auxevol_gfs[IDX4ptS(which_gf, ijk)] = ((double)rand() - (double)rand())/(double)(RAND_MAX);
        }
        
  for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
      auxevol_gfs[IDX4ptS(GAMMADD00GF, ijk)] = 1.0;
      auxevol_gfs[IDX4ptS(GAMMADD01GF, ijk)] = 0.0;
      auxevol_gfs[IDX4ptS(GAMMADD02GF, ijk)] = 0.0;
      auxevol_gfs[IDX4ptS(GAMMADD11GF, ijk)] = 1.0;
      auxevol_gfs[IDX4ptS(GAMMADD12GF, ijk)] = 0.0;
      auxevol_gfs[IDX4ptS(GAMMADD22GF, ijk)] = 1.0;

      auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, ijk)] = 1.0;
      auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, ijk)] = 0.0;
      auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, ijk)] = 0.0;
      auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, ijk)] = 1.0;
      auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, ijk)] = 0.0;
      auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, ijk)] = 1.0;

    }

  for(int which_gf=0; which_gf<NUM_EVOL_GFS; which_gf++) for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
    evol_gfs[IDX4ptS(which_gf, ijk)] = 0.0/0.0;
        }

//   // Now we fill our structs and call the relevant functions
//   for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
  
//       initialize_structs(ijk, &griddata.params, auxevol_gfs, &reconstructed_prims,
//                          &reconstructed_prims_r,
//                          &reconstructed_prims_l,
//                          &metric_quantities,
//                          &metric_face_quantities,
//                          &metric_quantities_derivatives);
                         
         
//      calculate_all_source_terms(&reconstructed_prims,
//                                 &metric_quantities,
//                                 &metric_quantities_derivatives,
//                                 &conservative_sources);
                                
//     evol_gfs[IDX4ptS(STILDED0GF, ijk)]  = conservative_sources.StildeD0_src;
//     evol_gfs[IDX4ptS(STILDED1GF, ijk)]  = conservative_sources.StildeD1_src;
//     evol_gfs[IDX4ptS(STILDED2GF, ijk)]  = conservative_sources.StildeD2_src;
//     evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] = conservative_sources.tau_tilde_src;
//     evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  = 0.0;

//     calculate_HLLE_fluxes0(&reconstructed_prims_r, 
//                            &reconstructed_prims_l,
//                            &metric_face_quantities,
//                            &conservative_fluxes);
                           
//     evol_gfs[IDX4ptS(STILDED0GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD0;
//     evol_gfs[IDX4ptS(STILDED1GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD1;
//     evol_gfs[IDX4ptS(STILDED2GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD2;
//     evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] += conservative_fluxes.HLLE_flux_tau_tilde;
//     evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  += conservative_fluxes.HLLE_flux_rho_star;

//     calculate_HLLE_fluxes1(&reconstructed_prims_r, 
//                            &reconstructed_prims_l,
//                            &metric_face_quantities,
//                            &conservative_fluxes);

//     evol_gfs[IDX4ptS(STILDED0GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD0;
//     evol_gfs[IDX4ptS(STILDED1GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD1;
//     evol_gfs[IDX4ptS(STILDED2GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD2;
//     evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] += conservative_fluxes.HLLE_flux_tau_tilde;
//     evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  += conservative_fluxes.HLLE_flux_rho_star;

//     calculate_HLLE_fluxes2(&reconstructed_prims_r, 
//                            &reconstructed_prims_l,
//                            &metric_face_quantities,
//                            &conservative_fluxes);

//     evol_gfs[IDX4ptS(STILDED0GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD0;
//     evol_gfs[IDX4ptS(STILDED1GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD1;
//     evol_gfs[IDX4ptS(STILDED2GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD2;
//     evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] += conservative_fluxes.HLLE_flux_tau_tilde;
//     evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  += conservative_fluxes.HLLE_flux_rho_star;                         
//  }
 
 
//  for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
//     printf("%.6e %.6e %.6e %.6e %.6e\n", evol_gfs[IDX4ptS(STILDED0GF, ijk)], 
//                                               evol_gfs[IDX4ptS(STILDED1GF, ijk)], 
//                                               evol_gfs[IDX4ptS(STILDED2GF, ijk)], 
//                                               evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)],
//                                               evol_gfs[IDX4ptS(RHO_STARGF, ijk)]);
    
// }


        // in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAVU0GF;
        //   out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF;
        //   out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF;
        // ww++;
        // in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAVU1GF;
        //   out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF;
        //   out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF;

  // EOS parameters
  const double neos = 1;
  const double K_poly = 2.43;
  const double rho_tab = 2.43;
  const double P_tab = 2.43;
  const double gamma_th = 2.43;
  const double eps_tab = 2.43;
  const double k_tab = 2.43;
  const double gamma_tab = 2.43;

  eos_struct eos;
  eos.neos=neos;
  eos.K_poly=K_poly;
  eos.rho_tab[0]=rho_tab;
  eos.P_tab[0]=P_tab;
  eos.gamma_th=gamma_th;
  eos.eps_tab[0]=eps_tab;
  eos.k_tab[0]=k_tab;
  eos.gamma_tab[0]=gamma_tab;

  // in_prims,out_prims_r, and out_prims_l are arrays of pointers to the actual gridfunctions.
  gf_and_gz_struct in_prims[MAXNUMVARS],out_prims_r[MAXNUMVARS],out_prims_l[MAXNUMVARS];
  int which_prims_to_reconstruct[MAXNUMVARS],num_prims_to_reconstruct;

  // mhdflux(i,j,k,flux_dirn,Ul  ,Ur  ,FACEVAL  ,FACEVAL_LAPSE_PSI4  ,eos, cmax[index],cmin[index],
  //   rho_star_flux[index],tau_flux[index],st_x_flux[index],st_y_flux[index],st_z_flux[index]);




  int ww=0;
  in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOBGF; 
  out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOB_RGF; 
  out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOB_LGF; 
  ww++;

  in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*PGF;     
  out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*P_RGF;
  out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*P_LGF;
  ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v0GF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V0_RGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V0_LGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v1GF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V1_RGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V1_LGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v2GF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V2_RGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V2_LGF;    
  ww++;

  in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU0GF;    
  out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU0GF;    
  out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU0GF;    
  ww++;

  in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU1GF;    
  out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU1GF;    
  out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU1GF;    
  ww++;

  in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU2GF;    
  out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU2GF;    
  out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU2GF;    
  ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU0_STAGGERGF; 
  // out_prims_r[ww].gf=BU0_staggerr; 
  // out_prims_l[ww].gf=BU0_staggerl; 
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU1_STAGGERGF; 
  // out_prims_r[ww].gf=BU1_staggerr; 
  // out_prims_l[ww].gf=BU1_staggerl; 
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU2_STAGGERGF; 
  // out_prims_r[ww].gf=BU2_staggerr;
  // out_prims_l[ww].gf=BU2_staggerl;
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v0RGF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V0_RRGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V0_RLGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v1RGF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V1_RRGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V1_RLGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v2RGF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V2_RRGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V2_RLGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v0LGF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V0_LRGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V0_LLGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v1LGF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V1_LRGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V1_LLGF;    
  // ww++;

  // in_prims[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*v2LGF;    
  // out_prims_r[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V2_LRGF;    
  // out_prims_l[ww].gf=auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V2_LLGF;    
  // ww++;


  // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
  for(int i=0;i<MAXNUMVARS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
  // Left/right variables are not yet defined, yet we set the # of gz's to zero BU1 default:
  for(int i=0;i<MAXNUMVARS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
  for(int i=0;i<MAXNUMVARS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }

  ww=0;
  which_prims_to_reconstruct[ww]=RHOB;      ww++;
  which_prims_to_reconstruct[ww]=PRESSURE;  ww++;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  //which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER;ww++;
  num_prims_to_reconstruct=ww;



  double *restrict ftilde_gf = (double *restrict)malloc(sizeof(double) * Nxx_plus_2NGHOSTS_tot);
  double *restrict temporary = (double *restrict)malloc(sizeof(double) * Nxx_plus_2NGHOSTS_tot);

  int flux_dirn;
  flux_dirn=1;

  cGH cctkGH;

  // First compute ftilde, which is used for flattening left and right face values
  // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  // ftilde_gf_compute(&griddata.params, cctkGH,cctk_lsh,flux_dirn, in_prims, ftilde_gf);


  // // This function is housed in the file: "reconstruct_set_of_prims_PPM.C"
  reconstruct_set_of_prims_PPM(&griddata.params, &cctkGH,cctk_lsh,flux_dirn,num_prims_to_reconstruct,
                               which_prims_to_reconstruct,
                               &eos,in_prims,out_prims_r,out_prims_l,ftilde_gf,temporary);
  //


  // Step 4: Free all allocated memory
  free(evol_gfs);
  free(auxevol_gfs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
  return 0;
}
