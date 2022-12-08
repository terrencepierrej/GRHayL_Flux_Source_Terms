# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
nrpy_dir_path = os.path.join("..")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step 1: The StildeD RHS *source* term
# Step P1: Import needed NRPy+ core modules:
from outputC import outputC, outCfunction, lhrh, add_to_Cfunction_dict, outC_function_dict # NRPy+: Core C code output module
import finite_difference as fin       # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par        # NRPy+: Parameter interface
import grid as gri                    # NRPy+: Functions having to do with numerical grids
import reference_metric as rfm        # NRPy+: Reference metric support
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import cmdline_helper as cmd          # NRPy+: Multi-platform Python command-line interface
import shutil, os, sys                # Standard Python modules for multiplatform OS-level functions
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends


# # Add struct for PPM algorithm
# with open(os.path.join(Ccodesrootdir,"NRPy_basic_defines.h"), "a") as file:
#     file.write(r"""
    
# typedef struct __gf_and_gz_struct__ {
#   REAL *gf;
#   int gz_lo[4],gz_hi[4];
# } gf_and_gz_struct;

# // Some additional constants needed for PPM:
# static const int VX=0,VY=1,VZ=2, BX_CENTER=3,BY_CENTER=4,BZ_CENTER=5,BX_STAGGER=6,BY_STAGGER=7,BZ_STAGGER=8,VXR=9,VYR=10,VZR=11,VXL=12,VYL=13,VZL=14;  //<-- Be _sure_ to define MAXNUMVARS appropriately!
# const int NUM_RECONSTRUCT_GFS = 15;

# """)
    
# import CurviGiRaFFE_NRPy.CurviGiRaFFE_C2P_P2C as CP
# # Add struct for current sheet algorithm
# with open(os.path.join(Ccodesrootdir,"NRPy_basic_defines.h"), "a") as file:
#     file.write(CP.data_struct)

# with open(os.path.join(Ccodesrootdir,"set_Cparameters.h"), "a") as file:
#     file.write(r"""

# const int kronecker_delta[4][3] = { { 0,0,0 },
#                                     { 1,0,0 },
#                                     { 0,1,0 },
#                                     { 0,0,1 } };
# """)


def set_CurviGiRaFFE_GFs():
    # We will pass values of the gridfunction on the cell faces into the function. This requires us
    # to declare them as C parameters in NRPy+. We will denote this with the _face infix/suffix.
    alpha_face = gri.register_gridfunctions("AUXEVOL","alpha_face")
    cf_face = gri.register_gridfunctions("AUXEVOL","cf_face")
    h_faceDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","h_faceDD","sym01",DIM=3)
    vet_faceU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","vet_faceU", DIM=3)

    alpha = gri.register_gridfunctions("AUXEVOL","alpha")
    cf = gri.register_gridfunctions("AUXEVOL","cf")
    hDD = ixp.register_gridfunctions_for_single_rank2("AUXEVOL","hDD","sym01",DIM=3)
    vetU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","vetU", DIM=3)

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    valenciav_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciav_rU",DIM=3)
    b_rU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","b_rU",DIM=3)
    valenciav_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciav_lU",DIM=3)
    b_lU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","b_lU",DIM=3)


    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciav_rrU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciav_rlU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciav_lrU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciav_llU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","bstagger_rU",DIM=3)
    ixp.register_gridfunctions_for_single_rank1("AUXEVOL","bstagger_lU",DIM=3)

    aD = ixp.register_gridfunctions_for_single_rank1("EVOL","aD")
    valenciavU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","valenciavU",DIM=3)
    bU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","bU",DIM=3)
    bstaggerU = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","bstaggerU",DIM=3)

    stildeD = ixp.register_gridfunctions_for_single_rank1("EVOL","stildeD")
    Stilde_flux_HLLED = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Stilde_flux_HLLED")

    psi6Phi = gri.register_gridfunctions("EVOL","psi6Phi")
    gri.register_gridfunctions("AUXEVOL","psi6_temp")
    gri.register_gridfunctions("AUXEVOL","psi6center")

    cmax_x = gri.register_gridfunctions("AUXEVOL","cmax_x")
    cmin_x = gri.register_gridfunctions("AUXEVOL","cmin_x")
    cmax_y = gri.register_gridfunctions("AUXEVOL","cmax_y")
    cmin_y = gri.register_gridfunctions("AUXEVOL","cmin_y")
    cmax_z = gri.register_gridfunctions("AUXEVOL","cmax_z")
    cmin_z = gri.register_gridfunctions("AUXEVOL","cmin_z")

    thismodule = __name__

    xi_damping = par.Cparameters("REAL",thismodule,"xi_damping",0.1)
    GAMMA_SPEED_LIMIT = par.Cparameters("REAL", thismodule, "GAMMA_SPEED_LIMIT", 10.0)  # Default 
    TINYDOUBLE = par.Cparameters("REAL",thismodule,"TINYDOUBLE",1e-100)
    M_PI  = par.Cparameters("#define",thismodule,["M_PI"], "")
    sqrt4pi = par.Cparameters("REAL",thismodule,"sqrt4pi","sqrt(4.0*M_PI)")

    # # And the function to which we'll write the output data:
    # Stilde_fluxD = ixp.register_gridfunctions_for_single_rank1("AUXEVOL","Stilde_fluxD",DIM=3)
    
def add_to_Cfunction_dict__CurviGiRaFFE_NRPy_all_functions(Ccodesdir):
    
    import CurviGiRaFFE_NRPy.CurviGiRaFFE_C2P_P2C as CP
#     import CurviGiRaFFE_NRPy.CurviGiRaFFE_NRPy_staggered_A2B as A2B
    import CurviGiRaFFE_NRPy.CurviGiRaFFE_NRPy_staggered_A2B_CC as A2B

    import CurviGiRaFFE_NRPy.CurviGiRaFFE_NRPy_Stilde_flux as SF

    import CurviGiRaFFE_NRPy.CurviGiRaFFE_Metric_Face_Values_Curvilinear_BSSN as FCVAL
    import CurviGiRaFFE_NRPy.GiRaFFE_NRPy_PPM as PPM

#     import CurviGiRaFFE_NRPy.GiRaFFE_NRPy_staggered_Afield_flux as AF
    import CurviGiRaFFE_NRPy.CurviGiRaFFE_NRPy_staggered_Afield_flux as AF

    import CurviGiRaFFE_NRPy.CurviGiRaFFE_Stilde_Source_Terms as source
    import CurviGiRaFFE_NRPy.CurviGiRaFFE_NRPy_staggered_Source_Terms as stgsrc

    FCVAL.add_to_Cfunction_dict__GiRaFFE_NRPy_FCVAL()
    PPM.add_to_Cfunction_dict__GiRaFFE_NRPy_PPM(Ccodesdir)
    
    SF.add_to_Cfunction_dict__Stilde_flux()
#     AF.add_to_Cfunction_dict__GiRaFFE_NRPy_staggered_Afield_flux()
    AF.add_to_Cfunction_dict__CurviGiRaFFE_NRPy_staggered_A_i_rhs_no_gauge_terms()
    source.add_to_Cfunction_dict__Stilde_SourceTerms()
    stgsrc.add_to_Cfunction_dict__CurviGiRaFFE_NRPy_staggered_Source_Terms()

    CP.add_to_Cfunction_dict__CurviGiRaFFE_NRPy_cons_to_prims()
    CP.add_to_Cfunction_dict__CurviGiRaFFE_NRPy_prims_to_cons()
    CP.add_to_Cfunction_dict__calculate_Cartesian_z_coordinate()
    CP.add_to_Cfunction_dict__count_pts_for_current_sheet_prescription()
    CP.add_to_Cfunction_dict__find_current_sheet_pts()
    CP.add_to_Cfunction_dict_apply_current_sheet_prescription()
    CP.add_to_Cfunction_dict_apply_obcs_valenciavU()
    CP.add_to_Cfunction_dict_apply_bcs_potentials()

    A2B.add_to_Cfunction_dict__GiRaFFE_NRPy_staggered_A2B()
    
    pre_body = r"""

#define WORKAROUND_ENABLED


void workaround_Valencia_to_Drift_velocity(const paramstruct *params, REAL *vU0, const REAL *alpha, const REAL *betaU0, const REAL flux_dirn) {
#include "set_Cparameters.h"
    // Converts Valencia 3-velocities to Drift 3-velocities for testing. The variable argument
    // vu0 is any Valencia 3-velocity component or reconstruction thereof.
#pragma omp parallel for
    for (int i2 = 2*(flux_dirn==3);i2 < Nxx_plus_2NGHOSTS2-1*(flux_dirn==3);i2++) for (int i1 = 2*(flux_dirn==2);i1 < Nxx_plus_2NGHOSTS1-1*(flux_dirn==2);i1++) for (int i0 = 2*(flux_dirn==1);i0 < Nxx_plus_2NGHOSTS0-1*(flux_dirn==1);i0++) {
        int ii = IDX3S(i0,i1,i2);
        // Here, we modify the velocity in place.
        vU0[ii] = alpha[ii]*vU0[ii]-betaU0[ii];
    }
}

void workaround_Drift_to_Valencia_velocity(const paramstruct *params, REAL *vU0, const REAL *alpha, const REAL *betaU0, const REAL flux_dirn) {
#include "set_Cparameters.h"
    // Converts Drift 3-velocities to Valencia 3-velocities for testing. The variable argument
    // vu0 is any drift (i.e. IllinoisGRMHD's definition) 3-velocity component or reconstruction thereof.
#pragma omp parallel for
    for (int i2 = 2*(flux_dirn==3);i2 < Nxx_plus_2NGHOSTS2-1*(flux_dirn==3);i2++) for (int i1 = 2*(flux_dirn==2);i1 < Nxx_plus_2NGHOSTS1-1*(flux_dirn==2);i1++) for (int i0 = 2*(flux_dirn==1);i0 < Nxx_plus_2NGHOSTS0-1*(flux_dirn==1);i0++) {
        int ii = IDX3S(i0,i1,i2);
        // Here, we modify the velocity in place.
        vU0[ii] = (vU0[ii]+betaU0[ii])/alpha[ii];
    }
}

"""
    
    body = r"""

    // First thing's first: initialize the RHSs to zero!
#pragma omp parallel for
    for(int ii=0;ii<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;ii++) {
        rhs_gfs[ii] = 0.0;
    }

    // Now, we set up a bunch of structs of pointers to properly guide the PPM algorithm.
    // They also count the number of ghostzones available.
    gf_and_gz_struct in_prims[NUM_RECONSTRUCT_GFS], out_prims_r[NUM_RECONSTRUCT_GFS], out_prims_l[NUM_RECONSTRUCT_GFS];
    int which_prims_to_reconstruct[NUM_RECONSTRUCT_GFS],num_prims_to_reconstruct;
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

    REAL *temporary = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*PSI6_TEMPGF; // Using dedicated temporary variables for the staggered grid
    REAL *psi6center = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*PSI6CENTERGF; // Because the prescription requires more acrobatics.
    // This sets pointers to the portion of auxevol_gfs containing the relevant gridfunction.
    int8_t ww=0;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAVU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAVU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAVU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*B_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGERU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGER_RU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGER_LU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGERU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGER_RU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGER_LU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGERU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGER_RU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*BSTAGGER_LU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RRU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RLU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RRU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RLU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RRU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_RLU2GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LRU0GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LLU0GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LRU1GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LLU1GF;
    ww++;
    in_prims[ww].gf      = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF;
      out_prims_r[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LRU2GF;
      out_prims_l[ww].gf = auxevol_gfs + Nxx_plus_2NGHOSTS_tot*VALENCIAV_LLU2GF;
    ww++;

    // Prims are defined AT ALL GRIDPOINTS, so we set the # of ghostzones to zero:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { in_prims[i].gz_lo[j]=0; in_prims[i].gz_hi[j]=0; }
    // Left/right variables are not yet defined, yet we set the # of gz's to zero by default:
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_r[i].gz_lo[j]=0; out_prims_r[i].gz_hi[j]=0; }
    for(int i=0;i<NUM_RECONSTRUCT_GFS;i++) for(int j=1;j<=3;j++) { out_prims_l[i].gz_lo[j]=0; out_prims_l[i].gz_hi[j]=0; }
"""
    
    
    body += r"""

  int flux_dirn;
  flux_dirn=0;
  interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
  // ftilde = 0 in GRFFE, since P=rho=0.

  /* There are two stories going on here:
   * 1) Computation of \partial_x on RHS of \partial_t {mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i-1/2,j,k), so that
   *    \partial_y F = [ F(i+1/2,j,k) - F(i-1/2,j,k) ] / dx
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) vx and vy are at (i,j,k), and we reconstruct them to (i-1/2,j,k) below. After
   *      this, we'll reconstruct again in the y-dir'n to get {vx,vy} at (i-1/2,j-1/2,k)
   * 2Ab) By_stagger is at (i,j+1/2,k), and we reconstruct below to (i-1/2,j+1/2,k). */
  ww=0;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  //which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^x face values to be consistent with BX_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexim1=IDX3S(i-1+(i==0),j,k); /* indexim1=0 when i=0 */
        out_prims_r[BX_CENTER].gf[index]=out_prims_l[BX_CENTER].gf[index]=in_prims[BX_STAGGER].gf[indexim1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  
  // TPJ SAYS: Set source terms first
  calculate_stildeDi_source_term(params,xx, auxevol_gfs, rhs_gfs);
  
  //calculate_StildeD0_source_term(params,auxevol_gfs,rhs_gfs);
  //calculate_Stilde_flux_D0(params,auxevol_gfs,rhs_gfs);
  //calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);
  
  calculate_Stilde_flux_D0(params, xx, auxevol_gfs,rhs_gfs);
  calculate_Stilde_rhsD(flux_dirn+1,params, xx, auxevol_gfs,rhs_gfs);

  // Note that we have already reconstructed vx and vy along the x-direction,
  //   at (i-1/2,j,k). That result is stored in v{x,y}{r,l}.  Bx_stagger data
  //   are defined at (i+1/2,j,k).
  
"""
    
    
    body += r"""

  // Next goal: reconstruct Bx, vx and vy at (i+1/2,j+1/2,k).
  flux_dirn=1;
  // ftilde = 0 in GRFFE, since P=rho=0.

  // in_prims[{VXR,VXL,VYR,VYL}].gz_{lo,hi} ghostzones are set to all zeros, which
  //    is incorrect. We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];

  /* There are two stories going on here:
   * 1) Computation of \partial_y on RHS of \partial_t {mhd_st_{x,y,z}},
   *    via PPM reconstruction onto (i,j-1/2,k), so that
   *    \partial_y F = [ F(i,j+1/2,k) - F(i,j-1/2,k) ] / dy
   * 2) Computation of \partial_t A_i, where A_i are *staggered* gridfunctions,
   *    where A_x is defined at (i,j+1/2,k+1/2), A_y at (i+1/2,j,k+1/2), etc.
   *    Ai_rhs = \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Az_rhs is defined at (i+1/2,j+1/2,k), and it depends on {Bx,By,vx,vy},
   *     so the trick is to reconstruct {Bx,By,vx,vy} cleverly to get to these
   *     staggered points. For example:
   * 2Aa) VXR = [right-face of vx reconstructed along x-direction above] is at (i-1/2,j,k),
   *      and we reconstruct it to (i-1/2,j-1/2,k) below. Similarly for {VXL,VYR,VYL}
   * 2Ab) Bx_stagger is at (i+1/2,j,k), and we reconstruct to (i+1/2,j-1/2,k) below
   * 2Ac) By_stagger is at (i-1/2,j+1/2,k) already for Az_rhs, from the previous step.
   * 2B) Ax_rhs is defined at (i,j+1/2,k+1/2), and it depends on {By,Bz,vy,vz}.
   *     Again the trick is to reconstruct these onto these staggered points.
   * 2Ba) Bz_stagger is at (i,j,k+1/2), and we reconstruct to (i,j-1/2,k+1/2) below */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VYR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VYL;       ww++;
  num_prims_to_reconstruct=ww;
#ifdef WORKAROUND_ENABLED
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
#ifdef WORKAROUND_ENABLED
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
  ww=0;
  // Reconstruct other primitives last!
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  //which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BX_STAGGER;ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^y face values to be consistent with BY_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexjm1=IDX3S(i,j-1+(j==0),k); /* indexjm1=0 when j=0 */
        out_prims_r[BY_CENTER].gf[index]=out_prims_l[BY_CENTER].gf[index]=in_prims[BY_STAGGER].gf[indexjm1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  
  //calculate_StildeD1_source_term(params,auxevol_gfs,rhs_gfs);
  //calculate_Stilde_flux_D1(params,auxevol_gfs,rhs_gfs);
  //calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);
  
  calculate_Stilde_flux_D1(params, xx, auxevol_gfs,rhs_gfs);
  calculate_Stilde_rhsD(flux_dirn+1,params, xx, auxevol_gfs,rhs_gfs);

"""
    
    import BSSN.BSSN_quantities as Bq

    Bq.phi_and_derivs()
    e6phi = (Bq.exp_m4phi**(0.5))**(-3.)

    sqrtgamma = e6phi*sp.sqrt(rfm.detgammahat)
    xxi_str = r"""

    REAL xx0 = xx[0][i];
    REAL xx1 = xx[1][j];
    REAL xx2 = xx[2][k];

    """
   #note: change name from psi6center to sqrtgamma 
    sqrtgamma_str = outputC(sqrtgamma, 'psi6center[index]',filename='returnstring',
                          prestring=xxi_str, params="includebraces=False").replace("cf", "CF[index]")

    
    body += r"""
  /*****************************************
   * COMPUTING RHS OF A_z, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_z - [gauge terms] = \psi^{6} (v^x B^y - v^y B^x).
   * A_z is defined at (i+1/2,j+1/2,k).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j-1/2,k)| {vxrr,vxrl,vxlr,vxll,vyrr,vyrl,vylr,vyll}
   * (i+1/2,j-1/2,k)| {Bx_stagger_r,Bx_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j+1/2,k)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Interpolates to i+1/2
#define IPH(METRICm1,METRICp0,METRICp1,METRICp2) (-0.0625*((METRICm1) + (METRICp2)) + 0.5625*((METRICp0) + (METRICp1)))
  // Next compute sqrt(gamma)=psi^6 at (i+1/2,j+1/2,k):
  // To do so, we first compute the sqrt of the metric determinant at all points:
  
  //TPJ SAYS: Could probably reuse interp_vars array
  //WARNING: divided by detgammahat above, but because of below interpolations may not work!
  REAL *CF = auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CFGF;
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k);
        /*
        const REAL gxx = auxevol_gfs[IDX4ptS(GAMMADD00GF,index)];
        const REAL gxy = auxevol_gfs[IDX4ptS(GAMMADD01GF,index)];
        const REAL gxz = auxevol_gfs[IDX4ptS(GAMMADD02GF,index)];
        const REAL gyy = auxevol_gfs[IDX4ptS(GAMMADD11GF,index)];
        const REAL gyz = auxevol_gfs[IDX4ptS(GAMMADD12GF,index)];
        const REAL gzz = auxevol_gfs[IDX4ptS(GAMMADD22GF,index)];
        psi6center[index] = sqrt( gxx*gyy*gzz
                               -  gxx*gyz*gyz
                               +2*gxy*gxz*gyz
                               -  gyy*gxz*gxz
                               -  gzz*gxy*gxy );
       */
       
            """+sqrtgamma_str+"""       
       
      }
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=1;j<Nxx_plus_2NGHOSTS1-2;j++) for(int i=1;i<Nxx_plus_2NGHOSTS0-2;i++) {
        temporary[IDX3S(i,j,k)]=
          IPH(IPH(psi6center[IDX3S(i-1,j-1,k)],psi6center[IDX3S(i,j-1,k)],psi6center[IDX3S(i+1,j-1,k)],psi6center[IDX3S(i+2,j-1,k)]),
              IPH(psi6center[IDX3S(i-1,j  ,k)],psi6center[IDX3S(i,j  ,k)],psi6center[IDX3S(i+1,j  ,k)],psi6center[IDX3S(i+2,j  ,k)]),
              IPH(psi6center[IDX3S(i-1,j+1,k)],psi6center[IDX3S(i,j+1,k)],psi6center[IDX3S(i+1,j+1,k)],psi6center[IDX3S(i+2,j+1,k)]),
              IPH(psi6center[IDX3S(i-1,j+2,k)],psi6center[IDX3S(i,j+2,k)],psi6center[IDX3S(i+1,j+2,k)],psi6center[IDX3S(i+2,j+2,k)]));
      }

  int A_directionz=3;
  A_i_rhs_no_gauge_terms(A_directionz,params,xx,out_prims_r,out_prims_l,temporary,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMAX_XGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMIN_XGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMAX_YGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMIN_YGF,
                         rhs_gfs+Nxx_plus_2NGHOSTS_tot*AD2GF);

  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not correct, so we fix
  //    this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VYR]=out_prims_r[VY];
  in_prims[VYL]=out_prims_l[VY];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VZL]=out_prims_l[VZ];

"""
    
    
    body += r"""

  flux_dirn=2;
  // ftilde = 0 in GRFFE, since P=rho=0.

  /* There are two stories going on here:
   * 1) Single reconstruction to (i,j,k-1/2) for {vx,vy,vz,Bx,By,Bz} to compute
   *    z-dir'n advection terms in \partial_t {mhd_st_{x,y,z}} at (i,j,k)
   * 2) Multiple reconstructions for *staggered* gridfunctions A_i:
   *    \partial_t A_i = \epsilon_{ijk} \psi^{6} (v^j B^k - v^j B^k),
   *    where \epsilon_{ijk} is the flat-space antisymmetric operator.
   * 2A) Ax_rhs is defined at (i,j+1/2,k+1/2), depends on v{y,z} and B{y,z}
   * 2Aa) v{y,z}{r,l} are at (i,j-1/2,k), so we reconstruct here to (i,j-1/2,k-1/2)
   * 2Ab) Bz_stagger{r,l} are at (i,j-1/2,k+1/2) already.
   * 2Ac) By_stagger is at (i,j+1/2,k), and below we reconstruct its value at (i,j+1/2,k-1/2)
   * 2B) Ay_rhs is defined at (i+1/2,j,k+1/2), depends on v{z,x} and B{z,x}.
   * 2Ba) v{x,z} are reconstructed to (i,j,k-1/2). Later we'll reconstruct again to (i-1/2,j,k-1/2).
   * 2Bb) Bz_stagger is at (i,j,k+1/2). Later we will reconstruct to (i-1/2,j,k+1/2).
   * 2Bc) Bx_stagger is at (i+1/2,j,k), and below we reconstruct its value at (i+1/2,j,k-1/2)
   */
  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VYR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VYL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  num_prims_to_reconstruct=ww;
#ifdef WORKAROUND_ENABLED
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
#ifdef WORKAROUND_ENABLED
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU1GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  interpolate_metric_gfs_to_cell_faces(params,auxevol_gfs,flux_dirn+1);
  // Reconstruct other primitives last!
  ww=0;
  which_prims_to_reconstruct[ww]=VX;        ww++;
  which_prims_to_reconstruct[ww]=VY;        ww++;
  which_prims_to_reconstruct[ww]=VZ;        ww++;
  which_prims_to_reconstruct[ww]=BX_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BY_CENTER; ww++;
  //which_prims_to_reconstruct[ww]=BZ_CENTER; ww++;
  which_prims_to_reconstruct[ww]=BX_STAGGER; ww++;
  which_prims_to_reconstruct[ww]=BY_STAGGER; ww++;
  num_prims_to_reconstruct=ww;
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
  //Right and left face values of BI_CENTER are used in GRFFE__S_i__flux computation (first to compute b^a).
  //   Instead of reconstructing, we simply set B^z face values to be consistent with BZ_STAGGER.
#pragma omp parallel for
  for(int k=0;k<Nxx_plus_2NGHOSTS2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        const int index=IDX3S(i,j,k), indexkm1=IDX3S(i,j,k-1+(k==0)); /* indexkm1=0 when k=0 */
        out_prims_r[BZ_CENTER].gf[index]=out_prims_l[BZ_CENTER].gf[index]=in_prims[BZ_STAGGER].gf[indexkm1]; }
  // Then add fluxes to RHS for hydro variables {vx,vy,vz}:
  // This function is housed in the file: "add_fluxes_and_source_terms_to_hydro_rhss.C"
  
  //calculate_StildeD2_source_term(params,auxevol_gfs,rhs_gfs);
  //calculate_Stilde_flux_D2(params,auxevol_gfs,rhs_gfs);
  //calculate_Stilde_rhsD(flux_dirn+1,params,auxevol_gfs,rhs_gfs);
  
  calculate_Stilde_flux_D2(params, xx, auxevol_gfs,rhs_gfs);
  calculate_Stilde_rhsD(flux_dirn+1,params, xx, auxevol_gfs,rhs_gfs);


  // in_prims[{VYR,VYL,VZR,VZL}].gz_{lo,hi} ghostzones are not set correcty.
  //    We fix this below.
  // [Note that this is a cheap operation, copying only 8 integers and a pointer.]
  in_prims[VXR]=out_prims_r[VX];
  in_prims[VZR]=out_prims_r[VZ];
  in_prims[VXL]=out_prims_l[VX];
  in_prims[VZL]=out_prims_l[VZ];
  // FIXME: lines above seem to be inconsistent with lines below.... Possible bug, not major enough to affect evolutions though.
  in_prims[VZR].gz_lo[1]=in_prims[VZR].gz_hi[1]=0;
  in_prims[VXR].gz_lo[1]=in_prims[VXR].gz_hi[1]=0;
  in_prims[VZL].gz_lo[1]=in_prims[VZL].gz_hi[1]=0;
  in_prims[VXL].gz_lo[1]=in_prims[VXL].gz_hi[1]=0;
  /*****************************************
   * COMPUTING RHS OF A_x, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_x - [gauge terms] = \psi^{6} (v^y B^z - v^z B^y).
   * A_x is defined at (i,j+1/2,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i,j-1/2,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i,j+1/2,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j-1/2,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i,j+1/2,k+1/2):
#pragma omp parallel for
  for(int k=1;k<Nxx_plus_2NGHOSTS2-2;k++) for(int j=1;j<Nxx_plus_2NGHOSTS1-2;j++) for(int i=0;i<Nxx_plus_2NGHOSTS0;i++) {
        temporary[IDX3S(i,j,k)]=
          IPH(IPH(psi6center[IDX3S(i,j-1,k-1)],psi6center[IDX3S(i,j,k-1)],psi6center[IDX3S(i,j+1,k-1)],psi6center[IDX3S(i,j+2,k-1)]),
              IPH(psi6center[IDX3S(i,j-1,k  )],psi6center[IDX3S(i,j,k  )],psi6center[IDX3S(i,j+1,k  )],psi6center[IDX3S(i,j+2,k  )]),
              IPH(psi6center[IDX3S(i,j-1,k+1)],psi6center[IDX3S(i,j,k+1)],psi6center[IDX3S(i,j+1,k+1)],psi6center[IDX3S(i,j+2,k+1)]),
              IPH(psi6center[IDX3S(i,j-1,k+2)],psi6center[IDX3S(i,j,k+2)],psi6center[IDX3S(i,j+1,k+2)],psi6center[IDX3S(i,j+2,k+2)]));
      }

  int A_directionx=1;
  A_i_rhs_no_gauge_terms(A_directionx,params,xx,out_prims_r,out_prims_l,temporary,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMAX_YGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMIN_YGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMAX_ZGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMIN_ZGF,
                         rhs_gfs+Nxx_plus_2NGHOSTS_tot*AD0GF);

"""
    
    
    body += r"""

  // We reprise flux_dirn=0 to finish up computations of Ai_rhs's!
  flux_dirn=0;
  // ftilde = 0 in GRFFE, since P=rho=0.

  ww=0;
  // NOTE! The order of variable reconstruction is important here,
  //   as we don't want to overwrite {vxr,vxl,vyr,vyl}!
  which_prims_to_reconstruct[ww]=VXR;       ww++;
  which_prims_to_reconstruct[ww]=VZR;       ww++;
  which_prims_to_reconstruct[ww]=VXL;       ww++;
  which_prims_to_reconstruct[ww]=VZL;       ww++;
  which_prims_to_reconstruct[ww]=BZ_STAGGER;ww++;
  num_prims_to_reconstruct=ww;
#ifdef WORKAROUND_ENABLED
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Valencia_to_Drift_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/
  // This function is housed in the file: "reconstruct_set_of_prims_PPM_GRFFE.C"
  reconstruct_set_of_prims_PPM_GRFFE_NRPy(params, auxevol_gfs, flux_dirn+1, num_prims_to_reconstruct,
                                          which_prims_to_reconstruct, in_prims, out_prims_r, out_prims_l, temporary);
#ifdef WORKAROUND_ENABLED
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU0GF,flux_dirn+1);
  workaround_Drift_to_Valencia_velocity(params,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHA_FACEGF,auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VET_FACEU2GF,flux_dirn+1);
#endif /*WORKAROUND_ENABLED*/

  /*****************************************
   * COMPUTING RHS OF A_y, BOOKKEEPING NOTE:
   * We want to compute
   * \partial_t A_y - [gauge terms] = \psi^{6} (v^z B^x - v^x B^z).
   * A_y is defined at (i+1/2,j,k+1/2).
   * ==========================
   * Where defined  | Variables
   * (i-1/2,j,k-1/2)| {vyrr,vyrl,vylr,vyll,vzrr,vzrl,vzlr,vzll}
   * (i+1/2,j,k-1/2)| {By_stagger_r,By_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i-1/2,j,k+1/2)| {Bz_stagger_r,Bz_stagger_l} (see Table 1 in arXiv:1007.2848)
   * (i,j,k)        | {phi}
   * ==========================
   ******************************************/
  // Next compute phi at (i+1/2,j,k+1/2):
#pragma omp parallel for
  for(int k=1;k<Nxx_plus_2NGHOSTS2-2;k++) for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) for(int i=1;i<Nxx_plus_2NGHOSTS0-2;i++) {
        temporary[IDX3S(i,j,k)]=
          IPH(IPH(psi6center[IDX3S(i-1,j,k-1)],psi6center[IDX3S(i,j,k-1)],psi6center[IDX3S(i+1,j,k-1)],psi6center[IDX3S(i+2,j,k-1)]),
              IPH(psi6center[IDX3S(i-1,j,k  )],psi6center[IDX3S(i,j,k  )],psi6center[IDX3S(i+1,j,k  )],psi6center[IDX3S(i+2,j,k  )]),
              IPH(psi6center[IDX3S(i-1,j,k+1)],psi6center[IDX3S(i,j,k+1)],psi6center[IDX3S(i+1,j,k+1)],psi6center[IDX3S(i+2,j,k+1)]),
              IPH(psi6center[IDX3S(i-1,j,k+2)],psi6center[IDX3S(i,j,k+2)],psi6center[IDX3S(i+1,j,k+2)],psi6center[IDX3S(i+2,j,k+2)]));
      }

  int A_directiony=2;
  A_i_rhs_no_gauge_terms(A_directiony,params,xx,out_prims_r,out_prims_l,temporary,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMAX_ZGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMIN_ZGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMAX_XGF,
                         auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CMIN_XGF,
                         rhs_gfs+Nxx_plus_2NGHOSTS_tot*AD1GF);


  // Next compute psi6phi_rhs, and add gauge terms to A_i_rhs terms!
  //   Note that in the following function, we don't bother with reconstruction, instead interpolating.
  // We need A^i, but only have A_i. So we add gtupij to the list of input variables.
  
  /*
  REAL *interp_vars[MAXNUMINTERP];
  ww=0;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BETAU0GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BETAU1GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BETAU2GF;   ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*GAMMADD00GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*GAMMADD01GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*GAMMADD02GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*GAMMADD11GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*GAMMADD12GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*GAMMADD22GF;  ww++;
  interp_vars[ww]=temporary;ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHAGF;   ww++;
  interp_vars[ww]=evol_gfs+Nxx_plus_2NGHOSTS_tot*AD0GF;      ww++;
  interp_vars[ww]=evol_gfs+Nxx_plus_2NGHOSTS_tot*AD1GF;      ww++;
  interp_vars[ww]=evol_gfs+Nxx_plus_2NGHOSTS_tot*AD2GF;      ww++;
  const int max_num_interp_variables=ww;
  
*/

  REAL *interp_vars[MAXNUMINTERP];
  ww=0;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VETU0GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VETU1GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VETU2GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*HDD00GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*HDD01GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*HDD02GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*HDD11GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*HDD12GF;  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*HDD22GF;  ww++;
  interp_vars[ww]=temporary;                                  ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*ALPHAGF;  ww++;
  interp_vars[ww]=evol_gfs+Nxx_plus_2NGHOSTS_tot*AD0GF;       ww++;
  interp_vars[ww]=evol_gfs+Nxx_plus_2NGHOSTS_tot*AD1GF;       ww++;
  interp_vars[ww]=evol_gfs+Nxx_plus_2NGHOSTS_tot*AD2GF;       ww++;
  interp_vars[ww]=auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CFGF;     ww++;
  const int max_num_interp_variables=ww;
  
  //   if(max_num_interp_variables>MAXNUMINTERP) {CCTK_VError(VERR_DEF_PARAMS,"Error: Didn't allocate enough space for interp_vars[]."); }
  // We are FINISHED with v{x,y,z}{r,l} and P{r,l} so we use these 8 gridfunctions' worth of space as temp storage.
  Lorenz_psi6phi_rhs__add_gauge_terms_to_A_i_rhs(params,xx, interp_vars,
                                                 evol_gfs+Nxx_plus_2NGHOSTS_tot*PSI6PHIGF,
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU0GF,  // .
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU1GF,  // .
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RU2GF,  // .
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU0GF,  // WARNING:
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU1GF,  // ALL VARIABLES
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_LU2GF,  // ON THESE LINES
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RRU0GF, // ARE OVERWRITTEN
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RLU0GF, // FOR TEMP STORAGE
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RRU1GF, // .
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RLU1GF, // .
                                                 auxevol_gfs+Nxx_plus_2NGHOSTS_tot*VALENCIAV_RRU2GF, // .
                                                 rhs_gfs+Nxx_plus_2NGHOSTS_tot*PSI6PHIGF,
                                                 rhs_gfs+Nxx_plus_2NGHOSTS_tot*AD0GF,
                                                 rhs_gfs+Nxx_plus_2NGHOSTS_tot*AD1GF,
                                                 rhs_gfs+Nxx_plus_2NGHOSTS_tot*AD2GF);

"""
    
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    c_type = "void"
    name = "CurviGiRaFFE_NRPy_RHSs"
    desc = "Compute CurviGiRaFFE RHSs"
    params = "const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs,REAL *restrict evol_gfs,REAL *restrict rhs_gfs"

    add_to_Cfunction_dict(
    includes=includes,
    desc=desc,
    c_type=c_type,
    name=name,
    params=params,
    prefunc=pre_body,
    body=body,
    rel_path_to_Cparams=os.path.join("."))
    
    
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    c_type = "void"
    name = "CurviGiRaFFE_NRPy_post_step"
    desc = "CurviGiRaFFE RHSs post time step"
    params = "paramstruct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict xx[3], REAL *restrict auxevol_gfs,REAL *restrict evol_gfs"

    body = r"""

    // First, apply BCs to AD and psi6Phi. Then calculate BU from AD
    //apply_bcs_potential(params,evol_gfs);
    apply_bcs_potentials(params, bcstruct, evol_gf_parity, evol_gfs);

    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
    
      
    GiRaFFE_compute_B_and_Bstagger_from_A(params, bcstruct, xx, auxevol_gf_parity,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*CFGF,
                                          evol_gfs+Nxx_plus_2NGHOSTS_tot*AD0GF,
                                          evol_gfs+Nxx_plus_2NGHOSTS_tot*AD1GF,
                                          evol_gfs+Nxx_plus_2NGHOSTS_tot*AD2GF,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BU0GF,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BU1GF,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BU2GF,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BSTAGGERU0GF,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BSTAGGERU1GF,
                                          auxevol_gfs+Nxx_plus_2NGHOSTS_tot*BSTAGGERU2GF);

    
    // Apply fixes to StildeD, then recompute the velocity at the new timestep.
    // Apply the current sheet prescription to the velocities
    //GiRaFFE_NRPy_cons_to_prims(params,xx,auxevol_gfs,evol_gfs);
    
    CurviGiRaFFE_NRPy_cons_to_prims(params, xx, evol_gfs, auxevol_gfs);

    /*Don't apply prescription
    if(params->num_grid_points_for_prescription==0) count_pts_for_current_sheet_prescription(params, xx);
    
    cs_struct csstruct;
    csstruct.cs = (cs_pt_struct *)malloc(sizeof(cs_pt_struct *)*(params->num_grid_points_for_prescription));
    
    if(params->num_grid_points_for_prescription==0) find_current_sheet_pts(params, &csstruct, xx);
    
    apply_current_sheet_prescription(params, &csstruct, xx, auxevol_gfs);

    free(csstruct.cs);
    //free(csstruct);
    */
    
    // Then, recompute StildeD to be consistent with the new velocities
    //GiRaFFE_NRPy_prims_to_cons(params,auxevol_gfs,evol_gfs);
    CurviGiRaFFE_NRPy_prims_to_cons(params, xx, auxevol_gfs, evol_gfs);
    // Finally, apply outflow boundary conditions to the velocities.
    //apply_bcs_velocity(params,auxevol_gfs);
    apply_obcs_valenciavU(params, bcstruct, xx, auxevol_gfs);

"""

    add_to_Cfunction_dict(
    includes=includes,
    desc=desc,
    c_type=c_type,
    name=name,
    params=params,
    body=body,
    rel_path_to_Cparams=os.path.join("."))

