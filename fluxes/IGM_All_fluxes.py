# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
GRMHD_dir_path = os.path.join("../GRMHD_formulation/")
sys.path.append(GRMHD_dir_path)
import GRMHD_equations_new_version as GRMHD    # NRPy+: Generate general relativistic magnetohydrodynamics equations

nrpy_dir_path = os.path.join("../nrpy/nrpytutorial/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)
    
from outputC import outputC, add_to_Cfunction_dict # NRPy+: Core C code output module
import sympy as sp               # SymPy: The Python computer algebra package upon which NRPy+ depends
import finite_difference as fin  # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import grid as gri               # NRPy+: Functions having to do with numerical grids
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import reference_metric as rfm   # NRPy+: Reference metric support
import cmdline_helper as cmd     # NRPy+: Multi-platform Python command-line interface


#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()


# We'll rewrite this assuming that we've passed the entire reconstructed
# gridfunctions. You could also do this with only one point, but then you'd
# need to declare everything as a Cparam in NRPy+
    
def calculate_GRMHD_Tmunu_and_contractions(formalism, flux_dirn, mom_comp, 
                                          gammaDD,betaU,alpha,
                                          rho_b, P, h, u4U, BU):
    GRMHD.set_up_base_vars(formalism=formalism)
    
    GRMHD.compute_vU_from_u4U__no_speed_limit(u4U)
    
    # First compute stress-energy tensor T4UU and T4UD:
    GRMHD.compute_sqrtgammaDET(gammaDD)
    GRMHD.compute_smallb4U(gammaDD, betaU, alpha, u4U, BU, GRMHD.sqrt4pi)
    GRMHD.compute_smallbsquared(gammaDD, betaU, alpha, GRMHD.smallb4U)

    # First compute stress-energy tensor T4UU and T4UD:
    GRMHD.compute_T4UU(gammaDD,betaU,alpha, rho_b, P, h, u4U, GRMHD.smallb4U, GRMHD.smallbsquared)
    GRMHD.compute_T4UD(gammaDD,betaU,alpha, GRMHD.T4UU)
    
    # Compute conservative variables in terms of primitive variables
    GRMHD.compute_rho_star(alpha, GRMHD.sqrtgammaDET, rho_b,u4U)
    GRMHD.compute_tau_tilde(alpha, GRMHD.sqrtgammaDET, GRMHD.T4UU,GRMHD.rho_star)
    GRMHD.compute_S_tildeD(alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    # Next compute fluxes of conservative variables
    GRMHD.compute_rho_star_fluxU(GRMHD.VU,GRMHD.rho_star)
    GRMHD.compute_tau_tilde_fluxU(alpha, GRMHD.sqrtgammaDET, GRMHD.VU, GRMHD.T4UU, GRMHD.rho_star)
    GRMHD.compute_S_tilde_fluxUD (alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    global U_rho_star, F_rho_star
    U_rho_star = GRMHD.rho_star
    F_rho_star = GRMHD.rho_star_fluxU[flux_dirn]

    global U_tau_tilde, F_tau_tilde
    U_tau_tilde = GRMHD.tau_tilde
    F_tau_tilde = GRMHD.tau_tilde_fluxU[flux_dirn]

    global U_S_tilde, F_S_tilde
    U_S_tilde = GRMHD.S_tildeD[mom_comp]
    F_S_tilde = GRMHD.S_tilde_fluxUD[flux_dirn][mom_comp]

def HLLE_solver(cmax, cmin, Fr, Fl, Ur, Ul):
    # This solves the Riemann problem for the mom_comp component of the momentum
    # flux StildeD in the flux_dirn direction.

    # st_j_flux = (c_\min f_R + c_\max f_L - c_\min c_\max ( st_j_r - st_j_l )) / (c_\min + c_\max)
    return (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin)




def calculate_HLLE_fluxes(formalism, flux_dirn, alpha_face, gamma_faceDD, beta_faceU,
                          u4_rU, u4_lU, B_rU, B_lU,
                          rho_b_r, rho_b_l,
                          P_r, P_l,
                          h_r, h_l,
                          cmin, cmax):

    global Stilde_flux_HLLED, rho_star_HLLE_flux, tau_tilde_HLLE_flux
    Stilde_flux_HLLED = ixp.zerorank1()
#     rescaled_Stilde_fluxD = ixp.zerorank1()
    for mom_comp in range(3):
        calculate_GRMHD_Tmunu_and_contractions(formalism, flux_dirn, mom_comp, 
                                               gamma_faceDD,beta_faceU,
                                               alpha_face,
                                               rho_b_r, P_r, h_r, u4_rU, 
                                               B_rU)
        
        F_S_tilde_r = F_S_tilde
        U_S_tilde_r = U_S_tilde
        
        if mom_comp==0:
            U_rho_star_r = U_rho_star
            F_rho_star_r = F_rho_star

            U_tau_tilde_r = U_tau_tilde
            F_tau_tilde_r = F_tau_tilde            
        
        calculate_GRMHD_Tmunu_and_contractions(formalism, flux_dirn, mom_comp, 
                                               gamma_faceDD,beta_faceU,
                                               alpha_face,
                                               rho_b_l, P_l, h_l, u4_lU, 
                                               B_lU)
        
        F_S_tilde_l = F_S_tilde
        U_S_tilde_l = U_S_tilde
        
        if mom_comp==0:
            U_rho_star_l = U_rho_star
            F_rho_star_l = F_rho_star

            U_tau_tilde_l = U_tau_tilde
            F_tau_tilde_l = F_tau_tilde
            
            # now calculate HLLE derived fluxes, and rescale all of them
            rho_star_HLLE_flux = HLLE_solver(cmax, cmin, 
                                      F_rho_star_r, F_rho_star_l, 
                                      U_rho_star_r, U_rho_star_l)
            
            tau_tilde_HLLE_flux = HLLE_solver(cmax, cmin, 
                                      F_tau_tilde_r, F_tau_tilde_l, 
                                      U_tau_tilde_r, U_tau_tilde_l)
        
        # Rescale the flux term, to be FD
        Stilde_flux_HLLED[mom_comp] = HLLE_solver(cmax, cmin, 
                                      F_S_tilde_r, F_S_tilde_l, 
                                      U_S_tilde_r, U_S_tilde_l)
        

def add_to_Cfunction_dict__GRMHD_fluxes(formalism="ADM", outCparams = "outCverbose=False,CSE_sorting=True"):
    
    sqrt4pi = sp.symbols("sqrt4pi")
    alpha_face = sp.symbols("alpha_face")
    if formalism=="BSSN":
        cf_face = sp.symbols("cf_face")
        h_faceDD = ixp.declarerank2("h_faceDD","sym01",DIM=3)
        vet_faceU = ixp.declarerank1("vet_faceU", DIM=3)
        beta_faceU = vet_faceU

        import BSSN.ADM_in_terms_of_BSSN as AitoB
        AitoB.ADM_in_terms_of_BSSN()

        import BSSN.BSSN_quantities as Bq

        gamma_faceDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gamma_faceDD[i][j] = AitoB.gammaDD[i][j].subs(Bq.hDD[i][j], h_faceDD[i][j]).subs(Bq.cf, cf_face)

    else:
        beta_faceU = ixp.declarerank1("beta_faceU") 
        gamma_faceDD = ixp.declarerank2("gamma_faceDD","sym01")

    # We'll need some more gridfunctions, now, to represent the reconstructions of BU and ValenciavU
    # on the right and left faces
    u4_rU = ixp.declarerank1("u4_rU", DIM=4)
    u4_lU = ixp.declarerank1("u4_lU", DIM=4)

    B_rU = ixp.declarerank1("B_rU", DIM=4)
    B_lU = ixp.declarerank1("B_lU", DIM=4)

    cmin_dirn0 = sp.symbols("cmin_dirn0")
    cmin_dirn1 = sp.symbols("cmin_dirn1")
    cmin_dirn2 = sp.symbols("cmin_dirn2")

    cmax_dirn0 = sp.symbols("cmax_dirn0")
    cmax_dirn1 = sp.symbols("cmax_dirn1")
    cmax_dirn2 = sp.symbols("cmax_dirn2")

    cmins = [cmin_dirn0, cmin_dirn1, cmin_dirn2]
    cmaxs = [cmax_dirn0, cmax_dirn1, cmax_dirn2]

    h_r = sp.symbols("h_r")
    h_l = sp.symbols("h_l")

    P_r = sp.symbols("P_r")
    P_l = sp.symbols("P_l")

    rho_b_r = sp.symbols("rhob_r")
    rho_b_l = sp.symbols("rhob_l")

    prims_velocities = ["u4_rU0", "u4_rU1", "u4_rU2", "u4_rU3",
                        "u4_lU0", "u4_lU1", "u4_lU2", "u4_lU3"]

    prims_mag_field = ["B_rU0", "B_rU1", "B_rU2",
                       "B_lU0", "B_lU1", "B_lU2"]


    other_prims = ["P", "h", "rhob"]

    characteristic_speeds = ["cmax", "cmin"]

    params = ["GAMMA_SPEED_LIMIT", "TINYDOUBLE", "sqrt4pi"]

    prestring = ""

    for var in prims_velocities + prims_mag_field:
        prestring += "const double "+str(var)+" = reconstructed_prims->"+str(var)+";\n"

    for var in other_prims:
        prestring += "const double "+str(var)+"_r = reconstructed_prims->"+str(var)+"_r;\n"
        prestring += "const double "+str(var)+"_l = reconstructed_prims->"+str(var)+"_l;\n"
    for var in params:
        prestring += "const double "+str(var)+" = rhss_params->"+str(var)+";\n"

    prestring += "const double "+str(alpha_face)+" = metric_face_quantities->"+str(alpha_face)+";\n"

    if formalism=="BSSN":
        prestring += "const double "+str(cf_face)+" = metric_face_quantities->"+str(cf_face)+";\n"

        for i in range(3):
            vetU_var = vet_faceU[i]
            prestring += "const double "+str(vetU_var)+" = metric_face_quantities->"+str(vetU_var)+";\n"

        for i in range(3):
            for j in range(3):
                hDD_var = h_faceDD[i][j]
                prestring += "const double "+str(hDD_var)+" = metric_face_quantities->"+str(hDD_var)+";\n"

    else:
        for i in range(3):
                betaU_var = beta_faceU[i]
                prestring += "const double "+str(betaU_var)+" = metric_face_quantities->"+str(betaU_var)+";\n"

        for i in range(3):
            for j in range(3):
                gammaDD_var = gamma_faceDD[i][j]
                prestring += "const double "+str(gammaDD_var)+" = metric_face_quantities->"+str(gammaDD_var)+";\n"
                
                
    vars_to_write = ["conservative_fluxes->StildeD0_rhs", "conservative_fluxes->StildeD1_rhs", "conservative_fluxes->StildeD2_rhs", 
                     "conservative_fluxes->rho_star_rhs", "conservative_fluxes->tau_tilde_rhs"]

    c_type = "void"
    params   = "const rhss_paramstruct *restrict rhss_params, "
    params  += "const prims_struct *restrict reconstructed_prims, "
    params  += "const metric_quantities_struct *restrict metric_face_quantities, "
    params  += "conservative_fluxes_struct *restrict conservative_fluxes"

    calc_char_speeds_params_str = "(rhss_params, reconstructed_prims, metric_face_quantities, conservative_fluxes)"
    
    for flux_dirn in range(3):
        var = cmins[flux_dirn]
        cmin_str = "const double "+str(var)+" = conservative_fluxes->"+str(var)+";\n"

        var = cmaxs[flux_dirn]
        cmax_str = "const double "+str(var)+" = conservative_fluxes->"+str(var)+";\n"
        
        calc_char_speeds_func_str = "calculate_characteristic_speed_"+str(flux_dirn)+"th_direction"

        calculate_HLLE_fluxes(formalism, flux_dirn, alpha_face, gamma_faceDD, beta_faceU,
                              u4_rU, u4_lU, B_rU, B_lU,
                              rho_b_r, rho_b_l,
                              P_r, P_l,
                              h_r, h_l,
                              cmins[flux_dirn], cmaxs[flux_dirn])

        vars_rhs = [Stilde_flux_HLLED[0], 
                    Stilde_flux_HLLED[1], 
                    Stilde_flux_HLLED[2], 
                    rho_star_HLLE_flux,
                    tau_tilde_HLLE_flux]

        body = outputC(vars_rhs, vars_to_write, params=outCparams, 
                   filename="returnstring", prestring=(calc_char_speeds_func_str+
                                                       calc_char_speeds_params_str+";\n\n"+
                                                       cmin_str+cmax_str+prestring))

        desc = "Compute the HLLE-derived fluxes on the left face in the " + str(flux_dirn) + "direction for all components."
        name = "calculate_HLLE_fluxes" + str(flux_dirn)

        add_to_Cfunction_dict(
#             includes=includes,
            desc=desc,
            name=name,
            params=params,
            body= body, 
            enableCparameters=False)