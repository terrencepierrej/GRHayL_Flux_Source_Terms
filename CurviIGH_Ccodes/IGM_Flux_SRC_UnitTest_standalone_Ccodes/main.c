#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
#include "inlined_functions.C"
#include "reconstruct_set_of_prims_PPM.C"
#include "loop_defines_reconstruction.h"
#include "add_fluxes_and_source_terms_to_hydro_rhss.C"
#include "compute_tau_rhs_extrinsic_curvature_terms_and_TUPmunu.C"
#include "mhdflux.C"

#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625
#define COMPUTE_FCVAL(METRICm2,METRICm1,METRIC,METRICp1) (AM2*(METRICm2) + AM1*(METRICm1) + A0*(METRIC) + A1*(METRICp1))

static inline void calculate_metric_face_values(const paramstruct *restrict params, 
                                  const int flux_dirn, const int metric_faces_gfs[],
                                  const int metric_gfs[],
                                  double *restrict auxevol_gfs) {
  #include "./set_Cparameters.h"

  LOOP_REGION(2, Nxx_plus_2NGHOSTS0 - 1,
              2, Nxx_plus_2NGHOSTS1 - 1,
              2, Nxx_plus_2NGHOSTS2 - 1) {

    int idxm2 = IDX3S(i0 - 2*kronecker_delta[flux_dirn+1][0], 
                      i1 - 2*kronecker_delta[flux_dirn+1][1], 
                      i2 - 2*kronecker_delta[flux_dirn+1][2]);
    
    int idxm1 = IDX3S(i0 - 1*kronecker_delta[flux_dirn+1][0], 
                      i1 - 1*kronecker_delta[flux_dirn+1][1], 
                      i2 - 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    for(int which_gf=0; which_gf<10; which_gf++){
      
      int gf_face = metric_faces_gfs[which_gf];
      int gf = metric_gfs[which_gf];

      auxevol_gfs[IDX4ptS(gf_face, idx)] = COMPUTE_FCVAL(auxevol_gfs[IDX4ptS(gf, idxm2)],
                                                         auxevol_gfs[IDX4ptS(gf, idxm1)],
                                                         auxevol_gfs[IDX4ptS(gf, idx  )],
                                                         auxevol_gfs[IDX4ptS(gf, idxp1)]);
      
      // printf("gf = %.15e, gf_face = %.15e\n", auxevol_gfs[IDX4ptS(gf     , idx)],
      //                                         auxevol_gfs[IDX4ptS(gf_face, idx)]);
    }
  }
}


static inline void read_from_binary_file_all(const char *restrict binary_file, const paramstruct *params, double *restrict auxevol_gfs){
  #include "./set_Cparameters.h"

  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  FILE *infile = fopen(binary_file, "rb");
  // double *Bz = (double  *)malloc(sizeof(Bx)*Nxx_plus_2NGHOSTS_tot);
  // fread(By, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);

  double correct_magic_number = 1.130814081305130e-9;
  double magic_number1, magic_number2, magic_number3, magic_number4, magic_number5;
  // fread(x, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  // fread(y, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  // fread(z, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOBGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*PGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOB_RGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*P_RGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOB_LGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*P_LGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number1, sizeof(double), 1, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_RU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_RU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_RU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_LU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_LU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_LU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number2, sizeof(double), 1, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number3, sizeof(double), 1, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BETAU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BETAU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BETAU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*GAMMADD00GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*GAMMADD01GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*GAMMADD02GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*GAMMADD11GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*GAMMADD12GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*GAMMADD22GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*ALPHAGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number4, sizeof(double), 1, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*KDD00GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*KDD01GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*KDD02GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*KDD11GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*KDD12GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*KDD22GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number5, sizeof(double), 1, infile);
  fclose(infile);
  if(magic_number1!=correct_magic_number){ printf("ERROR: magic_number1 does not match"); exit(1);}
  if(magic_number2!=correct_magic_number){ printf("ERROR: magic_number2 does not match"); exit(1);}
  if(magic_number3!=correct_magic_number){ printf("ERROR: magic_number3 does not match"); exit(1);}
  if(magic_number4!=correct_magic_number){ printf("ERROR: magic_number4 does not match"); exit(1);}
  if(magic_number5!=correct_magic_number){ printf("ERROR: magic_number5 does not match"); exit(1);}
}


static inline void read_from_binary_file_recons(const char *restrict binary_file, const paramstruct *params, double *restrict auxevol_gfs){
  #include "./set_Cparameters.h"

  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  FILE *infile = fopen(binary_file, "rb");
  // double *Bz = (double  *)malloc(sizeof(Bx)*Nxx_plus_2NGHOSTS_tot);
  // fread(By, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);

  double correct_magic_number = 1.130814081305130e-9;
  double magic_number1, magic_number2, magic_number3;
  // fread(x, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  // fread(y, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  // fread(z, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOB_RGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*P_RGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*RHOB_LGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*P_LGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number1, sizeof(double), 1, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_RU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_RU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_RU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_LU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_LU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*V_LU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number2, sizeof(double), 1, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
  fread(&magic_number3, sizeof(double), 1, infile);
  fclose(infile);
  if(magic_number1!=correct_magic_number){ printf("ERROR: magic_number1 does not match"); exit(1);}
  if(magic_number2!=correct_magic_number){ printf("ERROR: magic_number2 does not match"); exit(1);}
  if(magic_number3!=correct_magic_number){ printf("ERROR: magic_number3 does not match"); exit(1);}
}


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
  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  double *restrict     evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
  double *restrict etk_evol_gfs    = (double *restrict)malloc(sizeof(double) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);
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
    auxevol_gfs[IDX4ptS(which_gf, ijk)] = 0.0 / 0.0; // ((double)rand() - (double)rand())/(double)(RAND_MAX);
        }

  for(int which_gf=0; which_gf<NUM_EVOL_GFS; which_gf++) for(int ijk=0; ijk<Nxx_plus_2NGHOSTS_tot; ijk++) {
    evol_gfs[IDX4ptS(which_gf, ijk)] = 0.0/0.0;
        }

  /*
  After reading in the same random data used in the original IGM code, now we
  calculate face cell-face values and derivatives at cell centers
  */

  int metric_faces_gfs[10] = {ALPHA_FACEGF, BETA_FACEU0GF, BETA_FACEU1GF, BETA_FACEU2GF,
                          GAMMA_FACEDD00GF, GAMMA_FACEDD01GF, GAMMA_FACEDD02GF,
                          GAMMA_FACEDD11GF, GAMMA_FACEDD12GF, GAMMA_FACEDD22GF};

  int metric_gfs[10] = {ALPHAGF, BETAU0GF, BETAU1GF, BETAU2GF,
                          GAMMADD00GF, GAMMADD01GF, GAMMADD02GF,
                          GAMMADD11GF, GAMMADD12GF, GAMMADD22GF};

  double dxxi[4] = {0, griddata.params.dxx0, griddata.params.dxx1, griddata.params.dxx2};

  int flux_dirn = 0;
  read_from_binary_file_all("input_data1.txt", &griddata.params, auxevol_gfs);

  calculate_metric_face_values(&griddata.params, 
                                flux_dirn, 
                                metric_faces_gfs,
                                metric_gfs,
                                auxevol_gfs);

  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - (NGHOSTS),
              NGHOSTS, Nxx_plus_2NGHOSTS1 - (NGHOSTS),
              NGHOSTS, Nxx_plus_2NGHOSTS2 - (NGHOSTS)) {
    
    int ijk  = IDX3S(i0, i1, i2);

    auxevol_gfs[IDX4ptS(U4U0GF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(U4_RU0GF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(U4_LU0GF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(HGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(H_RGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(H_LGF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(GAMMA_TH_RGF, ijk)] = 4.567;
    auxevol_gfs[IDX4ptS(GAMMA_TH_LGF, ijk)] = 4.567;

    auxevol_gfs[IDX4ptS(EPSILON_TH_RGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(EPSILON_TH_LGF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(DPCOLD_DRHOB_RGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(DPCOLD_DRHOB_LGF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(U4U1GF, ijk)] = auxevol_gfs[IDX4ptS(VU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4U2GF, ijk)] = auxevol_gfs[IDX4ptS(VU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4U3GF, ijk)] = auxevol_gfs[IDX4ptS(VU2GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU1GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU2GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU3GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU2GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU1GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU2GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU3GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU2GF, ijk)];

    initialize_structs(ijk, &griddata.params, auxevol_gfs, 
                            &reconstructed_prims,
                            &reconstructed_prims_r,
                            &reconstructed_prims_l,
                            &metric_quantities,
                            &metric_face_quantities,
                            &metric_quantities_derivatives);
    // printf("indices = %d %d %d\n", i0, i1, i2);
    calculate_HLLE_fluxes0(&reconstructed_prims_r, 
                           &reconstructed_prims_l,
                           &metric_face_quantities, 
                           &conservative_fluxes);

    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD0;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD1;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD2;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, ijk)]  = conservative_fluxes.HLLE_flux_rho_star;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,ijk)]  = conservative_fluxes.HLLE_flux_tau_tilde;

    // printf("%d %d %d %.15e %.15e %.15e %.15e %.15e\n", i0, i1, i2,
    //                                                    conservative_fluxes.HLLE_flux_StildeD0,
    //                                                    conservative_fluxes.HLLE_flux_StildeD1,
    //                                                    conservative_fluxes.HLLE_flux_StildeD2,
    //                                                    conservative_fluxes.HLLE_flux_rho_star,
    //                                                    conservative_fluxes.HLLE_flux_tau_tilde);
  }

  double inv_dxx = 1.0 / dxxi[flux_dirn+1];
  printf("%.15e\n", inv_dxx);


  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {
    
    int idx  = IDX3S(i0, i1, i2);

    int idxp1 = IDX3S(i0 + kronecker_delta[flux_dirn+1][0], 
                      i1 + kronecker_delta[flux_dirn+1][1], 
                      i2 + kronecker_delta[flux_dirn+1][2]);

    evol_gfs[IDX4ptS(STILDED0GF, idx)]  = inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED1GF, idx)]  = inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED2GF, idx)]  = inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idxp1)]);
    evol_gfs[IDX4ptS(RHO_STARGF, idx)]  = inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idxp1)]);
    evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]  = inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idxp1)]);

    // printf("%d %d %d %.15e %.15e %.15e %.15e %.15e\n", i0, i1, i2, 
    //                                  evol_gfs[IDX4ptS(STILDED0GF, idx)],
    //                                  evol_gfs[IDX4ptS(STILDED1GF, idx)],
    //                                  evol_gfs[IDX4ptS(STILDED2GF, idx)],
    //                                  evol_gfs[IDX4ptS(RHO_STARGF, idx)],
    //                                  evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]);
  }

  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

      auxevol_gfs[IDX4ptS(ALPHA_DD0GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idxp1)]
                                                      - auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD00GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD10GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD20GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx  )]);


      auxevol_gfs[IDX4ptS(GAMMADD_DD000GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD010GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD020GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD110GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD120GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD220GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx  )]);
      
      // printf("gf = %.15e, gf_face = %.15e\n", auxevol_gfs[IDX4ptS(gf     , idx)],
      //                                         auxevol_gfs[IDX4ptS(gf_face, idx)]);
      // printf("%d %d %d %.15e\n", i0, i1, i2, auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx)] - 1);
  }


  flux_dirn = 1;
  read_from_binary_file_recons("input_data2.txt", &griddata.params, auxevol_gfs);

  calculate_metric_face_values(&griddata.params, 
                                flux_dirn, 
                                metric_faces_gfs,
                                metric_gfs,
                                auxevol_gfs);

  
  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int ijk  = IDX3S(i0, i1, i2);

    auxevol_gfs[IDX4ptS(U4U0GF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(U4_RU0GF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(U4_LU0GF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(HGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(H_RGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(H_LGF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(U4U1GF, ijk)] = auxevol_gfs[IDX4ptS(VU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4U2GF, ijk)] = auxevol_gfs[IDX4ptS(VU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4U3GF, ijk)] = auxevol_gfs[IDX4ptS(VU2GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU1GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU2GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU3GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU2GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU1GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU2GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU3GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU2GF, ijk)];

    initialize_structs(ijk, &griddata.params, auxevol_gfs, 
                            &reconstructed_prims,
                            &reconstructed_prims_r,
                            &reconstructed_prims_l,
                            &metric_quantities,
                            &metric_face_quantities,
                            &metric_quantities_derivatives);

    calculate_HLLE_fluxes1(&reconstructed_prims_r, 
                           &reconstructed_prims_l,
                           &metric_face_quantities, 
                           &conservative_fluxes);

    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD0;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD1;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD2;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, ijk)]  = conservative_fluxes.HLLE_flux_rho_star;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,ijk)]  = conservative_fluxes.HLLE_flux_tau_tilde;
  }

  inv_dxx = 1.0 / dxxi[flux_dirn+1];
  printf("%.15e\n", inv_dxx);


  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {
    
    int idx  = IDX3S(i0, i1, i2);

    int idxp1 = IDX3S(i0 + kronecker_delta[flux_dirn+1][0], 
                      i1 + kronecker_delta[flux_dirn+1][1], 
                      i2 + kronecker_delta[flux_dirn+1][2]);

    evol_gfs[IDX4ptS(STILDED0GF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED1GF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED2GF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idxp1)]);
    evol_gfs[IDX4ptS(RHO_STARGF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idxp1)]);
    evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idxp1)]);
  }

  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

      auxevol_gfs[IDX4ptS(ALPHA_DD1GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idxp1)]
                                                      - auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD01GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD11GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD21GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx  )]);


      auxevol_gfs[IDX4ptS(GAMMADD_DD001GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD011GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD021GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD111GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD121GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD221GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx  )]);

      
      // printf("gf = %.15e, gf_face = %.15e\n", auxevol_gfs[IDX4ptS(gf     , idx)],
      //                                         auxevol_gfs[IDX4ptS(gf_face, idx)]);

       // printf("%d %d %d %.15e\n", i0, i1, i2, auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx)] - 1);

  }

  flux_dirn = 2;
  read_from_binary_file_recons("input_data3.txt", &griddata.params, auxevol_gfs);

  calculate_metric_face_values(&griddata.params, 
                                flux_dirn, 
                                metric_faces_gfs,
                                metric_gfs,
                                auxevol_gfs);

  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {

    int ijk  = IDX3S(i0, i1, i2);

    auxevol_gfs[IDX4ptS(U4U0GF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(U4_RU0GF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(U4_LU0GF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(HGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(H_RGF, ijk)] = 1.0;
    auxevol_gfs[IDX4ptS(H_LGF, ijk)] = 1.0;

    auxevol_gfs[IDX4ptS(U4U1GF, ijk)] = auxevol_gfs[IDX4ptS(VU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4U2GF, ijk)] = auxevol_gfs[IDX4ptS(VU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4U3GF, ijk)] = auxevol_gfs[IDX4ptS(VU2GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU1GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU2GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_RU3GF, ijk)] = auxevol_gfs[IDX4ptS(V_RU2GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU1GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU0GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU2GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU1GF, ijk)];
    auxevol_gfs[IDX4ptS(U4_LU3GF, ijk)] = auxevol_gfs[IDX4ptS(V_LU2GF, ijk)];

    initialize_structs(ijk, &griddata.params, auxevol_gfs, 
                            &reconstructed_prims,
                            &reconstructed_prims_r,
                            &reconstructed_prims_l,
                            &metric_quantities,
                            &metric_face_quantities,
                            &metric_quantities_derivatives);

    calculate_HLLE_fluxes2(&reconstructed_prims_r, 
                           &reconstructed_prims_l,
                           &metric_face_quantities, 
                           &conservative_fluxes);

    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD0;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD1;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, ijk)]  = conservative_fluxes.HLLE_flux_StildeD2;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, ijk)]  = conservative_fluxes.HLLE_flux_rho_star;
    auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,ijk)]  = conservative_fluxes.HLLE_flux_tau_tilde;

    }

  inv_dxx = 1.0 / dxxi[flux_dirn+1];
  printf("%.15e\n", inv_dxx);


  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {
    
    int idx  = IDX3S(i0, i1, i2);

    int idxp1 = IDX3S(i0 + kronecker_delta[flux_dirn+1][0], 
                      i1 + kronecker_delta[flux_dirn+1][1], 
                      i2 + kronecker_delta[flux_dirn+1][2]);

    evol_gfs[IDX4ptS(STILDED0GF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED0GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED1GF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED1GF, idxp1)]);
    evol_gfs[IDX4ptS(STILDED2GF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_STILDED2GF, idxp1)]);
    evol_gfs[IDX4ptS(RHO_STARGF, idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_RHO_STARGF, idxp1)]);
    evol_gfs[IDX4ptS(TAU_TILDEGF,idx)]  += inv_dxx*(auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idx)] - auxevol_gfs[IDX4ptS(HLLE_FLUX_TAU_TILDEGF,idxp1)]);
  }

  LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
              NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
    
    int idxp1 = IDX3S(i0 + 1*kronecker_delta[flux_dirn+1][0], 
                      i1 + 1*kronecker_delta[flux_dirn+1][1], 
                      i2 + 1*kronecker_delta[flux_dirn+1][2]);

    int idx  = IDX3S(i0, i1, i2);

      auxevol_gfs[IDX4ptS(ALPHA_DD2GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idxp1)]
                                                      - auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD02GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD12GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx  )]);

      auxevol_gfs[IDX4ptS(BETAU_DD22GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idxp1)]
                                                       - auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx  )]);


      auxevol_gfs[IDX4ptS(GAMMADD_DD002GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD012GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD022GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD112GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD122GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx  )]);

      auxevol_gfs[IDX4ptS(GAMMADD_DD222GF, idx)] = inv_dxx*(auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idxp1)]
                                                          - auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx  )]);

      
      // printf("gf = %.15e, gf_face = %.15e\n", auxevol_gfs[IDX4ptS(gf     , idx)],
      //                                         auxevol_gfs[IDX4ptS(gf_face, idx)]);
      // printf("%d %d %d %.15e\n", i0, i1, i2, auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx)] - 1);
  }

    // Now we fill our structs and call the relevant functions
    LOOP_REGION(NGHOSTS - 1, Nxx_plus_2NGHOSTS0 - 2,
                NGHOSTS - 1, Nxx_plus_2NGHOSTS1 - 2,
                NGHOSTS - 1, Nxx_plus_2NGHOSTS2 - 2) {
        int idx  = IDX3S(i0, i1, i2);
    
      initialize_structs(idx, &griddata.params, auxevol_gfs, &reconstructed_prims,
                           &reconstructed_prims_r,
                           &reconstructed_prims_l,
                           &metric_quantities,
                           &metric_face_quantities,
                           &metric_quantities_derivatives);
                           
           
      calculate_all_source_terms(&reconstructed_prims,
                                 &metric_quantities,
                                 &metric_quantities_derivatives,
                                 &conservative_sources);
                                  
      evol_gfs[IDX4ptS(STILDED0GF, idx)]  += conservative_sources.StildeD0_src;
      evol_gfs[IDX4ptS(STILDED1GF, idx)]  += conservative_sources.StildeD1_src;
      evol_gfs[IDX4ptS(STILDED2GF, idx)]  += conservative_sources.StildeD2_src;
      evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] += conservative_sources.tau_tilde_src;
      // evol_gfs[IDX4ptS(RHO_STARGF, idx)]  = 0.0;

      // printf("%d %d %d :::: %.15e %.15e %.15e %.15e\n", i0, i1, i2, 
      //                                             evol_gfs[IDX4ptS(TAU_TILDEGF, idx)],
      //                                             evol_gfs[IDX4ptS(STILDED0GF, idx)],
      //                                             evol_gfs[IDX4ptS(STILDED1GF, idx)],
      //                                             evol_gfs[IDX4ptS(STILDED2GF, idx)]);
    }

    FILE *infile = fopen("output_rhs_data.txt", "rb");
    double rhs_correct_magic_number = 9.524300707856655e-3;
    double rhs_magic_number1, rhs_magic_number2, rhs_magic_number3;
    fread(&rhs_magic_number1, sizeof(double), 1, infile);
    fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*RHO_STARGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*TAU_TILDEGF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*STILDED0GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    fread(&rhs_magic_number2, sizeof(double), 1, infile);
    fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*STILDED1GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    fread(etk_evol_gfs + Nxx_plus_2NGHOSTS_tot*STILDED2GF, sizeof(double), Nxx_plus_2NGHOSTS_tot, infile);
    fread(&rhs_magic_number3, sizeof(double), 1, infile);
    fclose(infile);
    if(rhs_magic_number1!=rhs_correct_magic_number){ printf("ERROR: magic_number1 does not match"); exit(1);}
    if(rhs_magic_number2!=rhs_correct_magic_number){ printf("ERROR: magic_number2 does not match"); exit(1);}
    if(rhs_magic_number3!=rhs_correct_magic_number){ printf("ERROR: magic_number3 does not match"); exit(1);}

    LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS-1,
                NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS-1,
                NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS-1) {

      int idx  = IDX3S(i0, i1, i2);
      printf("%d %d %d ::: %.3f %.3f %.3f %.3f %.3f\n", 
      i0, i1, i2,
      log10(fabs(0.5*(evol_gfs[IDX4ptS(STILDED0GF, idx)] - etk_evol_gfs[IDX4ptS(STILDED0GF, idx)]) / (evol_gfs[IDX4ptS(STILDED0GF, idx)] + etk_evol_gfs[IDX4ptS(STILDED0GF, idx)]))),
      log10(fabs(0.5*(evol_gfs[IDX4ptS(STILDED1GF, idx)] - etk_evol_gfs[IDX4ptS(STILDED1GF, idx)]) / (evol_gfs[IDX4ptS(STILDED1GF, idx)] + etk_evol_gfs[IDX4ptS(STILDED1GF, idx)]))),
      log10(fabs(0.5*(evol_gfs[IDX4ptS(STILDED2GF, idx)] - etk_evol_gfs[IDX4ptS(STILDED2GF, idx)]) / (evol_gfs[IDX4ptS(STILDED2GF, idx)] + etk_evol_gfs[IDX4ptS(STILDED2GF, idx)]))),
      log10(fabs(0.5*(evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] - etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)]) / (evol_gfs[IDX4ptS(TAU_TILDEGF, idx)] + etk_evol_gfs[IDX4ptS(TAU_TILDEGF, idx)]))),
      log10(fabs(0.5*(evol_gfs[IDX4ptS(RHO_STARGF, idx)] - etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)]) / (evol_gfs[IDX4ptS(RHO_STARGF, idx)] + etk_evol_gfs[IDX4ptS(RHO_STARGF, idx)]))));
  }

 
  // Step 4: Free all allocated memory
  free(evol_gfs);
  free(auxevol_gfs);
  for(int i=0;i<3;i++) free(griddata.xx[i]);
  return 0;
}
