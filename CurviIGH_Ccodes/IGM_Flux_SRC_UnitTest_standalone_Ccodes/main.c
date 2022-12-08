#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include "inlined_functions.C"
#include "reconstruct_set_of_prims_PPM.C"
#include "loop_defines_reconstruction.h"
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
  //          by hand or with command-line input, as desired.
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

  // We begin by generating random inital data.

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

  // Now we fill our structs and call the relevant functions
  for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
  
      initialize_structs(ijk, &griddata.params, auxevol_gfs, &reconstructed_prims,
                         &reconstructed_prims_r,
                         &reconstructed_prims_l,
                         &metric_quantities,
                         &metric_face_quantities,
                         &metric_quantities_derivatives);
                         
         
     calculate_all_source_terms(&reconstructed_prims,
                                &metric_quantities,
                                &metric_quantities_derivatives,
                                &conservative_sources);
                                
    evol_gfs[IDX4ptS(STILDED0GF, ijk)]  = conservative_sources.StildeD0_src;
    evol_gfs[IDX4ptS(STILDED1GF, ijk)]  = conservative_sources.StildeD1_src;
    evol_gfs[IDX4ptS(STILDED2GF, ijk)]  = conservative_sources.StildeD2_src;
    evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] = conservative_sources.tau_tilde_src;
    evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  = 0.0;

    calculate_HLLE_fluxes0(&reconstructed_prims_r, 
                           &reconstructed_prims_l,
                           &metric_face_quantities,
                           &conservative_fluxes);
                           
    evol_gfs[IDX4ptS(STILDED0GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD0;
    evol_gfs[IDX4ptS(STILDED1GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD1;
    evol_gfs[IDX4ptS(STILDED2GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD2;
    evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] += conservative_fluxes.HLLE_flux_tau_tilde;
    evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  += conservative_fluxes.HLLE_flux_rho_star;

    calculate_HLLE_fluxes1(&reconstructed_prims_r, 
                           &reconstructed_prims_l,
                           &metric_face_quantities,
                           &conservative_fluxes);

    evol_gfs[IDX4ptS(STILDED0GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD0;
    evol_gfs[IDX4ptS(STILDED1GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD1;
    evol_gfs[IDX4ptS(STILDED2GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD2;
    evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] += conservative_fluxes.HLLE_flux_tau_tilde;
    evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  += conservative_fluxes.HLLE_flux_rho_star;

    calculate_HLLE_fluxes2(&reconstructed_prims_r, 
                           &reconstructed_prims_l,
                           &metric_face_quantities,
                           &conservative_fluxes);

    evol_gfs[IDX4ptS(STILDED0GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD0;
    evol_gfs[IDX4ptS(STILDED1GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD1;
    evol_gfs[IDX4ptS(STILDED2GF, ijk)]  += conservative_fluxes.HLLE_flux_StildeD2;
    evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)] += conservative_fluxes.HLLE_flux_tau_tilde;
    evol_gfs[IDX4ptS(RHO_STARGF, ijk)]  += conservative_fluxes.HLLE_flux_rho_star;                         
 }
 
 
 for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
    printf("%.6e %.6e %.6e %.6e %.6e\n", evol_gfs[IDX4ptS(STILDED0GF, ijk)], 
                                              evol_gfs[IDX4ptS(STILDED1GF, ijk)], 
                                              evol_gfs[IDX4ptS(STILDED2GF, ijk)], 
                                              evol_gfs[IDX4ptS(TAU_TILDEGF, ijk)],
                                              evol_gfs[IDX4ptS(RHO_STARGF, ijk)]);
    
}
    

  FILE *infile = fopen("input_data3.txt", "r");
  int i, j, k;
  double x, y, z;
  double rho_b;
  double rho_bl;
  double rho_br;
  double P;
  double Pl;
  double Pr;
  double vx;
  double vxl;
  double vxr;
  double vy;
  double vyl;
  double vyr;
  double vz;
  double vzl;
  double vzr;
  double Bx;
  double Bxl;
  double Bxr;
  double By;
  double Byl;
  double Byr;
  double Bz;
  double Bzl;
  double Bzr;
  double gxx;
  double gxy;
  double gxz;
  double gyy;
  double gyz;
  double gzz;
  double lapm1;
  double betax;
  double betay;
  double betaz;
  double kxx;
  double kxy;
  double kxz;
  double kyy;
  double kyz;
  double kzz;
  while (fscanf(infile, "%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
    &i, &j, &k, &x, &y, &z, 
    &rho_b,
    &rho_bl,
    &rho_br,
    &P,
    &Pl,
    &Pr,
    &vx,
    &vxl,
    &vxr,
    &vy,
    &vyl,
    &vyr,
    &vz,
    &vzl,
    &vzr,
    &Bx,
    &Bxl,
    &Bxr,
    &By,
    &Byl,
    &Byr,
    &Bz,
    &Bzl,
    &Bzr,
    &gxx,
    &gxy,
    &gxz,
    &gyy,
    &gyz,
    &gzz,
    &lapm1,
    &betax,
    &betay,
    &betaz,
    &kxx,
    &kxy,
    &kxz,
    &kyy,
    &kyz,
    &kzz) == 1)
      printf("%d\n", k);

  fclose (infile);
  // Step 4: Free all allocated memory
  free(evol_gfs);
  free(auxevol_gfs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
  return 0;
}
