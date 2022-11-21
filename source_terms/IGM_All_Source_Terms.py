# Step 0: Add NRPy's directory to the path
# https://stackoverflow.com/questions/16780014/import-file-from-parent-directory
import os,sys
GRMHD_dir_path = os.path.join("../GRMHD_formulation/")
sys.path.append(GRMHD_dir_path)
import GRMHD_equations_new_version as GRMHD    # NRPy+: Generate general relativistic magnetohydrodynamics equations

nrpy_dir_path = os.path.join("../nrpy/nrpytutorial/")
if nrpy_dir_path not in sys.path:
    sys.path.append(nrpy_dir_path)

# Step P1: Import needed NRPy+ core modules:
from outputC import outputC, add_to_Cfunction_dict  # NRPy+: Core C code output module
import finite_difference as fin       # NRPy+: Finite difference C code generation module
import NRPy_param_funcs as par        # NRPy+: Parameter interface
import grid as gri                    # NRPy+: Functions having to do with numerical grids
import reference_metric as rfm        # NRPy+: Reference metric support
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import cmdline_helper as cmd          # NRPy+: Multi-platform Python command-line interface
import shutil, os, sys                # Standard Python modules for multiplatform OS-level functions
import sympy as sp                # SymPy: The Python computer algebra package upon which NRPy+ depends

thismodule = __name__

par.initialize_param(par.glb_param(type="bool", module=thismodule, 
       parname="using_Valencia_velocity", defaultval=True))

#Step 0: Set the spatial dimension parameter to 3.
par.set_parval_from_str("grid::DIM", 3)
DIM = par.parval_from_str("grid::DIM")

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()

def add_to_Cfunction_dict__Stilde_SourceTerms(formalism="ADM"):    
    # Generate SymPy symbolic expressions
    GRMHD.set_up_base_vars(formalism=formalism)

    GRMHD.compute_vU_from_u4U__no_speed_limit(GRMHD.u4U)   

    GRMHD.compute_sqrtgammaDET(GRMHD.gammaDD)
    GRMHD.compute_smallb4U(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.u4U, GRMHD.BU, GRMHD.sqrt4pi)
    GRMHD.compute_smallbsquared(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.smallb4U)

    # First compute stress-energy tensor T4UU and T4UD:
    GRMHD.compute_T4UU(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.rho_b,GRMHD.P,GRMHD.h,GRMHD.u4U, GRMHD.smallb4U, GRMHD.smallbsquared)
    GRMHD.compute_T4UD(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.T4UU)

    # Compute conservative variables in terms of primitive variables
    GRMHD.compute_rho_star(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.rho_b,GRMHD.u4U)
    GRMHD.compute_tau_tilde(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.T4UU,GRMHD.rho_star)
    GRMHD.compute_S_tildeD(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    # Next compute fluxes of conservative variables
    GRMHD.compute_rho_star_fluxU(GRMHD.VU,GRMHD.rho_star)
    GRMHD.compute_tau_tilde_fluxU(GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.VU, GRMHD.T4UU, GRMHD.rho_star)
    GRMHD.compute_S_tilde_fluxUD (GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.T4UD)

    # Then declare derivatives & compute g4DD_zerotimederiv_dD
    GRMHD.compute_g4DD_zerotimederiv_dD(GRMHD.gammaDD,GRMHD.betaU,GRMHD.alpha, GRMHD.gammaDD_dD,GRMHD.betaU_dD,GRMHD.alpha_dD)

    # Then compute source terms on tau_tilde and S_tilde equations
    GRMHD.compute_tau_tilde_source_term(GRMHD.KDD,GRMHD.betaU,GRMHD.alpha, GRMHD.sqrtgammaDET, GRMHD.alpha_dD, GRMHD.T4UU)
    GRMHD.compute_S_tilde_source_termD(GRMHD.alpha,GRMHD.sqrtgammaDET,GRMHD.g4DD_zerotimederiv_dD, GRMHD.T4UU)
    
    
    tau_tilde_source_term_free_symbols = GRMHD.tau_tilde_source_term.free_symbols
    S_tilde_source_termD0_free_symbols = GRMHD.S_tilde_source_termD[0].free_symbols
    S_tilde_source_termD1_free_symbols = GRMHD.S_tilde_source_termD[1].free_symbols
    S_tilde_source_termD2_free_symbols = GRMHD.S_tilde_source_termD[2].free_symbols

    all_free_sysmbols = tau_tilde_source_term_free_symbols.union(\
                                 S_tilde_source_termD0_free_symbols, 
                                 S_tilde_source_termD1_free_symbols,
                                 S_tilde_source_termD2_free_symbols)

    prims_velocities = ["u4U0", "u4U1", "u4U2", "u4U3"]
    
    prims = prims_velocities + ["BU0", "BU1", "BU2", "P", "h", "rhob"]
    params = ["GAMMA_SPEED_LIMIT", "TINYDOUBLE", "sqrt4pi"]


    prestring = ""

    for var in all_free_sysmbols:
        if str(var) in prims:
            prestring += "const double "+str(var)+" = prims->"+str(var)+";\n"
        if str(var) in params:
            prestring += "const double "+str(var)+" = rhss_params->"+str(var)+";\n"

    prestring += "const double "+str(GRMHD.alpha)+" = metric_quantities->"+str(GRMHD.alpha)+";\n"

    if formalism=="BSSN":
        prestring += "const double "+str(GRMHD.Bq.trK)+" = metric_quantities->"+str(GRMHD.Bq.trK)+";\n"
        prestring += "const double "+str(GRMHD.Bq.cf)+" = metric_quantities->"+str(GRMHD.Bq.cf)+";\n"

        for i in range(3):
            vetU_var = GRMHD.Bq.vetU[i]
            prestring += "const double "+str(vetU_var)+" = metric_quantities->"+str(vetU_var)+";\n"

        for i in range(3):
            for j in range(3):
                aDD_var = GRMHD.Bq.aDD[i][j]
                prestring += "const double "+str(aDD_var)+" = metric_quantities->"+str(aDD_var)+";\n"

        for i in range(3):
            for j in range(3):
                hDD_var = GRMHD.Bq.hDD[i][j]
                prestring += "const double "+str(hDD_var)+" = metric_quantities->"+str(hDD_var)+";\n"

        for var in all_free_sysmbols:
            if "_dD" in str(var):
                prestring += "const double "+str(var)+" = metric_quantities_derivatives->"+str(var)+";\n"

    else:
        for i in range(3):
                betaU_var = GRMHD.betaU[i]
                prestring += "const double "+str(betaU_var)+" = metric_quantities->"+str(betaU_var)+";\n"

        for i in range(3):
            for j in range(3):
                KDD_var = GRMHD.KDD[i][j]
                prestring += "const double "+str(KDD_var)+" = metric_quantities->"+str(KDD_var)+";\n"

        for i in range(3):
            for j in range(3):
                gammaDD_var = GRMHD.gammaDD[i][j]
                prestring += "const double "+str(gammaDD_var)+" = metric_quantities->"+str(gammaDD_var)+";\n"

        for var in all_free_sysmbols:
            if "_dD" in str(var):
                prestring += "const double "+str(var)+" = metric_quantities_derivatives->"+str(var)+";\n"
    
       
    outCparams = "outCverbose=False,CSE_sorting=canonical,CSE_enable=True"
    desc = "Adds source term and connection terms to Stilde, rho_star and tau_tilde"
#     includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    name = "calculate_all_source_terms"
    vars_to_write = ["conservative_sources->StildeD0_rhs", "conservative_sources->StildeD1_rhs", "conservative_sources->StildeD2_rhs", "conservative_rhs->tau_tilde_rhs"]

    vars_rhs = [GRMHD.S_tilde_source_termD[0], 
                GRMHD.S_tilde_source_termD[1], 
                GRMHD.S_tilde_source_termD[2], 
                GRMHD.tau_tilde_source_term]

    body = outputC(vars_rhs, vars_to_write, params=outCparams, 
               filename="returnstring", prestring=prestring)

    c_type = "void"
    params   = "const rhss_paramstruct *restrict rhss_params, "
    params  += "const prims_struct *restrict prims, "
    params  += "const metric_quantities_struct *restrict metric_quantities, "
    params  += "const metric_quantities_derivatives_struct *restrict metric_quantities_derivatives, "
    params  += "conservative_sources_struct *restrict conservative_sources"
    
    add_to_Cfunction_dict(
#         includes=includes,
        desc=desc,
        c_type=c_type, name=name, params=params,
        body=body,
        enableCparameters=False,)

