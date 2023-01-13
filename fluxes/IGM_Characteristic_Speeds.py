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
import indexedexp as ixp         # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import NRPy_param_funcs as par   # NRPy+: Parameter interface
import reference_metric as rfm   # NRPy+: Reference metric support

thismodule = __name__

CoordSystem = "Cartesian"

# Set coordinate system to dst_basis
par.set_parval_from_str("reference_metric::CoordSystem", CoordSystem)
rfm.reference_metric()

# We'll write this as a function so that we can calculate the expressions on-demand for any choice of i
def find_cp_cm(flux_dirn, g4UU, 
               smallbsquared, u4U, Gamma_th, 
               epsilon_th, dPcold_drhob, h, rho_b):
    # Outputs: cplus,cminus
    vA2 = smallbsquared / (rho_b*h + smallbsquared)
    cs2 = (dPcold_drhob + Gamma_th*(Gamma_th-1)*epsilon_th) / h
    v02 = vA2 + (cs2)*(1 - vA2)
    
    a = (1 - v02)*(u4U[0]**2) - v02*g4UU[0][0]
    b = 2*v02*g4UU[flux_dirn+1][0] - 2*u4U[flux_dirn+1]*u4U[0]*(1 - v02)
    c = (1 - v02)*(u4U[flux_dirn+1]**2) - v02*g4UU[flux_dirn+1][flux_dirn+1]

    # Now, we are free to solve the quadratic equation as usual. We take care to avoid passing a
    # negative value to the sqrt function.
    detm = b*b - sp.sympify(4)*a*c

    import Min_Max_and_Piecewise_Expressions as noif
    detm = sp.sqrt(noif.max_noif(sp.sympify(0),detm))
    global cplus,cminus
    # note that these correspond to a single interface, left or right
    cplus_tmp  = sp.Rational(1,2)* (-b/a + detm/a)
    cminus_tmp = sp.Rational(1,2)*-( b/a + detm/a)
    
    cminus = noif.min_noif(cplus_tmp, cminus_tmp)
    cplus  = noif.max_noif(cplus_tmp, cminus_tmp)

    # the above in C code
    # if (cplus < cminus) {
    # CCTK_REAL cp = cminus;
    # cminus = cplus;
    # cplus = cp;

# We'll write this as a function, and call it within HLLE_solver, below.
def find_cmax_cmin(flux_dirn, gamma_faceDD, beta_faceU, alpha_face,
                   smallbsquared_r, smallbsquared_l, u4U_r, u4U_l, 
                   Gamma_th_r, Gamma_th_l, 
                   epsilon_th_r, epsilon_th_l, dPcold_drhob_r, dPcold_drhob_l, 
                   h_r, h_l, rho_b_r, rho_b_l):
    # Inputs:  flux direction flux_dirn, Inverse metric g4_faceUU, shift beta_faceU,
    #          lapse alpha_face, metric determinant gammadet_face
    # Outputs: maximum and minimum characteristic speeds cmax and cmin
    # First, we need to find the characteristic speeds on each face
    import BSSN.ADMBSSN_tofrom_4metric as AB4m
    AB4m.g4UU_ito_BSSN_or_ADM("ADM",gamma_faceDD, beta_faceU, alpha_face)

    # Original needed for GRMHD
    find_cp_cm(flux_dirn, AB4m.g4UU, smallbsquared_r, u4U_r, Gamma_th_r, epsilon_th_r, dPcold_drhob_r, h_r, rho_b_r)
    cpr = cplus
    cmr = cminus
    
    find_cp_cm(flux_dirn, AB4m.g4UU, smallbsquared_l, u4U_l, Gamma_th_l, epsilon_th_l, dPcold_drhob_l, h_l, rho_b_l)
    cpl = cplus
    cml = cminus

    # The following algorithms have been verified with random floats:

    global cmax,cmin    
#   // Then compute cmax, cmin. This is required for the HLL flux.
#   original C code
#   CCTK_REAL cmaxL =  MAX(0.0,MAX(cplusl,cplusr));
#   CCTK_REAL cminL = -MIN(0.0,MIN(cminusl,cminusr));
    
    # Now, we need to set cmax to the larger of cpr,cpl, and 0
    import Min_Max_and_Piecewise_Expressions as noif
    cmax = noif.max_noif(0.0, noif.max_noif(cpl, cpr))
    # And then, set cmin to the smaller of cmr,cml, and 0
    cmin = -noif.min_noif(0.0, noif.min_noif(cml, cmr))
    
    
def add_to_Cfunction_dict__GRMHD_characteristic_speeds(formalism="ADM", outCparams = "outCverbose=False,CSE_sorting=True"):
    
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

    B_rU = ixp.declarerank1("B_rU")
    B_lU = ixp.declarerank1("B_lU")

    cmins = ["conservative_fluxes->cmin_dirn0", 
             "conservative_fluxes->cmin_dirn1", 
             "conservative_fluxes->cmin_dirn2"]

    cmaxs = ["conservative_fluxes->cmax_dirn0", 
             "conservative_fluxes->cmax_dirn1", 
             "conservative_fluxes->cmax_dirn2"]

    h_r = sp.symbols("h_r")
    h_l = sp.symbols("h_l")

    rho_b_r = sp.symbols("rhob_r")
    rho_b_l = sp.symbols("rhob_l")

    Gamma_th_r = sp.symbols("Gamma_th_r")
    Gamma_th_l = sp.symbols("Gamma_th_l")

    epsilon_th_r = sp.symbols("epsilon_th_r")
    epsilon_th_l = sp.symbols("epsilon_th_l")

    dPcold_drhob_r = sp.symbols("dPcold_drhob_r")
    dPcold_drhob_l = sp.symbols("dPcold_drhob_l")

    GRMHD.compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4_rU, B_rU, sqrt4pi)
    GRMHD.compute_smallbsquared(gamma_faceDD, beta_faceU, alpha_face, GRMHD.smallb4U)
    smallbsquared_r = GRMHD.smallbsquared

    GRMHD.compute_smallb4U(gamma_faceDD, beta_faceU, alpha_face, u4_lU, B_lU, sqrt4pi)
    GRMHD.compute_smallbsquared(gamma_faceDD, beta_faceU, alpha_face, GRMHD.smallb4U)
    smallbsquared_l = GRMHD.smallbsquared

    cmins_rhs = []
    cmaxs_rhs = []

    for flux_dirn in range(3):
        find_cmax_cmin(flux_dirn, gamma_faceDD, beta_faceU, alpha_face,
                       smallbsquared_r, smallbsquared_l, u4_rU, u4_lU, 
                       Gamma_th_r, Gamma_th_l, 
                       epsilon_th_r, epsilon_th_l, dPcold_drhob_r, dPcold_drhob_l, 
                       h_r, h_l, rho_b_r, rho_b_l)

        cmins_rhs.append(cmin)
        cmaxs_rhs.append(cmax)


    rho_b_r = sp.symbols("rhob_r")
    rho_b_l = sp.symbols("rhob_l")

    prims_velocities_r = ["u4_rU0", "u4_rU1", "u4_rU2", "u4_rU3"]
    prims_velocities_l = ["u4_lU0", "u4_lU1", "u4_lU2", "u4_lU3"]

    prims_mag_field_r = ["B_rU0", "B_rU1", "B_rU2"]
    prims_mag_field_l = ["B_lU0", "B_lU1", "B_lU2"]


    other_prims = ["P", "h", "rhob","Gamma_th", "epsilon_th", "dPcold_drhob"]

    characteristic_speeds = ["cmax", "cmin"]

    prestring = ""

    for var in prims_velocities_r + prims_mag_field_r:
        prestring += "const double "+var+" = reconstructed_prims_r->"+var.replace("_r", "")+";\n"
        
    for var in prims_velocities_l + prims_mag_field_l:
        prestring += "const double "+var+" = reconstructed_prims_l->"+var.replace("_l", "")+";\n"

    for var in other_prims:
        prestring += "const double "+var+"_r = reconstructed_prims_r->"+var+";\n"
        prestring += "const double "+var+"_l = reconstructed_prims_l->"+var+";\n"

    prestring += "const double "+str(alpha_face)+" = metric_face_quantities->"+str(alpha_face)+";\n"
    
    checker = []
    
    if formalism=="BSSN":
        prestring += "const double "+str(cf_face)+" = metric_face_quantities->"+str(cf_face)+";\n"

        for i in range(3):
            vetU_var = vet_faceU[i]
            prestring += "const double "+str(vetU_var)+" = metric_face_quantities->"+str(vetU_var)+";\n"

        for i in range(3):
            for j in range(3):
                hDD_var = h_faceDD[i][j]
                if hDD_var in checker:
                    continue
                prestring += "const double "+str(hDD_var)+" = metric_face_quantities->"+str(hDD_var)+";\n"
                checker.append(hDD_var)

    else:
        for i in range(3):
                betaU_var = beta_faceU[i]
                prestring += "const double "+str(betaU_var)+" = metric_face_quantities->"+str(betaU_var)+";\n"

        for i in range(3):
            for j in range(3):
                gammaDD_var = gamma_faceDD[i][j]
                if gammaDD_var in checker:
                    continue
                prestring += "const double "+str(gammaDD_var)+" = metric_face_quantities->"+str(gammaDD_var)+";\n"
                checker.append(gammaDD_var)

    c_type = "void"
    
    params  = "const reconstructed_prims_struct *restrict reconstructed_prims_r, const reconstructed_prims_struct *restrict reconstructed_prims_l,"
    params  += "const metric_face_quantities_struct *restrict metric_face_quantities, "
    params  += "conservative_fluxes_struct *restrict conservative_fluxes"

    for flux_dirn in range(3):
        write_speeds_str = [cmins[flux_dirn], cmaxs[flux_dirn]]
        write_speeds_rhs_str = [cmins_rhs[flux_dirn], cmaxs_rhs[flux_dirn]]

        body = outputC(write_speeds_rhs_str, write_speeds_str, params=outCparams, 
                       filename="returnstring", prestring=prestring)

        desc = "Compute the characteristic speeds in" + str(flux_dirn) +"th direction"
        name = "calculate_characteristic_speed_" + str(flux_dirn) +"th_direction"
        includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]

        add_to_Cfunction_dict(
                includes=includes,
                desc=desc,
                name=name,
                params=params,
                body= body, 
                enableCparameters=False)