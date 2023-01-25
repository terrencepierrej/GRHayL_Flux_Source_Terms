#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Initialize structs used for flux and source calculations
 */
void initialize_structs(const int ijk, const paramstruct *params, const double *restrict auxevol_gfs, reconstructed_prims_struct *restrict reconstructed_prims, reconstructed_prims_struct *restrict reconstructed_prims_r, reconstructed_prims_struct *restrict reconstructed_prims_l, metric_quantities_struct *restrict metric_quantities, metric_face_quantities_struct *restrict metric_face_quantities, metric_quantities_derivatives_struct *restrict metric_quantities_derivatives) {
#include "./set_Cparameters.h"


    
    reconstructed_prims->u4U0 = auxevol_gfs[IDX4ptS(U4U0GF, ijk)];
    reconstructed_prims_r->u4U0 = auxevol_gfs[IDX4ptS(U4_RU0GF, ijk)];
    reconstructed_prims_l->u4U0 = auxevol_gfs[IDX4ptS(U4_LU0GF, ijk)];
    
    reconstructed_prims->u4U1 = auxevol_gfs[IDX4ptS(U4U1GF, ijk)];
    reconstructed_prims_r->u4U1 = auxevol_gfs[IDX4ptS(U4_RU1GF, ijk)];
    reconstructed_prims_l->u4U1 = auxevol_gfs[IDX4ptS(U4_LU1GF, ijk)];
    
    reconstructed_prims->u4U2 = auxevol_gfs[IDX4ptS(U4U2GF, ijk)];
    reconstructed_prims_r->u4U2 = auxevol_gfs[IDX4ptS(U4_RU2GF, ijk)];
    reconstructed_prims_l->u4U2 = auxevol_gfs[IDX4ptS(U4_LU2GF, ijk)];
    
    reconstructed_prims->u4U3 = auxevol_gfs[IDX4ptS(U4U3GF, ijk)];
    reconstructed_prims_r->u4U3 = auxevol_gfs[IDX4ptS(U4_RU3GF, ijk)];
    reconstructed_prims_l->u4U3 = auxevol_gfs[IDX4ptS(U4_LU3GF, ijk)];
    
    reconstructed_prims->BU0 = auxevol_gfs[IDX4ptS(BU0GF, ijk)];
    reconstructed_prims_r->BU0 = auxevol_gfs[IDX4ptS(B_RU0GF, ijk)];
    reconstructed_prims_l->BU0 = auxevol_gfs[IDX4ptS(B_LU0GF, ijk)];
    
    reconstructed_prims->BU1 = auxevol_gfs[IDX4ptS(BU1GF, ijk)];
    reconstructed_prims_r->BU1 = auxevol_gfs[IDX4ptS(B_RU1GF, ijk)];
    reconstructed_prims_l->BU1 = auxevol_gfs[IDX4ptS(B_LU1GF, ijk)];
    
    reconstructed_prims->BU2 = auxevol_gfs[IDX4ptS(BU2GF, ijk)];
    reconstructed_prims_r->BU2 = auxevol_gfs[IDX4ptS(B_RU2GF, ijk)];
    reconstructed_prims_l->BU2 = auxevol_gfs[IDX4ptS(B_LU2GF, ijk)];

    reconstructed_prims_r->P = auxevol_gfs[IDX4ptS(P_RGF, ijk)];
    reconstructed_prims_l->P = auxevol_gfs[IDX4ptS(P_LGF, ijk)];

    reconstructed_prims->P = auxevol_gfs[IDX4ptS(PGF, ijk)];

    reconstructed_prims_r->h = auxevol_gfs[IDX4ptS(H_RGF, ijk)];
    reconstructed_prims_l->h = auxevol_gfs[IDX4ptS(H_LGF, ijk)];

    reconstructed_prims->h = auxevol_gfs[IDX4ptS(HGF, ijk)];

    reconstructed_prims_r->rhob = auxevol_gfs[IDX4ptS(RHOB_RGF, ijk)];
    reconstructed_prims_l->rhob = auxevol_gfs[IDX4ptS(RHOB_LGF, ijk)];

    reconstructed_prims->rhob = auxevol_gfs[IDX4ptS(RHOBGF, ijk)];

    reconstructed_prims_r->Gamma_th = auxevol_gfs[IDX4ptS(GAMMA_TH_RGF, ijk)];
    reconstructed_prims_l->Gamma_th = auxevol_gfs[IDX4ptS(GAMMA_TH_LGF, ijk)];

    reconstructed_prims_r->epsilon_th = auxevol_gfs[IDX4ptS(EPSILON_TH_RGF, ijk)];
    reconstructed_prims_l->epsilon_th = auxevol_gfs[IDX4ptS(EPSILON_TH_LGF, ijk)];

    reconstructed_prims_r->dPcold_drhob = auxevol_gfs[IDX4ptS(DPCOLD_DRHOB_RGF, ijk)];
    reconstructed_prims_l->dPcold_drhob = auxevol_gfs[IDX4ptS(DPCOLD_DRHOB_LGF, ijk)];
    
    metric_quantities->alpha = auxevol_gfs[IDX4ptS(ALPHAGF, ijk)];
 
    metric_face_quantities->alpha_face = auxevol_gfs[IDX4ptS(ALPHA_FACEGF, ijk)];
    
    metric_quantities->betaU0 = auxevol_gfs[IDX4ptS(BETAU0GF, ijk)];
 
    metric_face_quantities->beta_faceU0 = auxevol_gfs[IDX4ptS(BETA_FACEU0GF, ijk)];
    
    metric_quantities->betaU1 = auxevol_gfs[IDX4ptS(BETAU1GF, ijk)];
 
    metric_face_quantities->beta_faceU1 = auxevol_gfs[IDX4ptS(BETA_FACEU1GF, ijk)];
    
    metric_quantities->betaU2 = auxevol_gfs[IDX4ptS(BETAU2GF, ijk)];
 
    metric_face_quantities->beta_faceU2 = auxevol_gfs[IDX4ptS(BETA_FACEU2GF, ijk)];
    
    metric_quantities->gammaDD01 = auxevol_gfs[IDX4ptS(GAMMADD01GF, ijk)];
 
    metric_face_quantities->gamma_faceDD01 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, ijk)];
    
    metric_quantities->gammaDD02 = auxevol_gfs[IDX4ptS(GAMMADD02GF, ijk)];
 
    metric_face_quantities->gamma_faceDD02 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, ijk)];
    
    metric_quantities->gammaDD12 = auxevol_gfs[IDX4ptS(GAMMADD12GF, ijk)];
 
    metric_face_quantities->gamma_faceDD12 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, ijk)];
    
    metric_quantities->gammaDD00 = auxevol_gfs[IDX4ptS(GAMMADD00GF, ijk)];
 
    metric_face_quantities->gamma_faceDD00 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, ijk)];
    
    metric_quantities->gammaDD11 = auxevol_gfs[IDX4ptS(GAMMADD11GF, ijk)];
 
    metric_face_quantities->gamma_faceDD11 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, ijk)];
    
    metric_quantities->gammaDD22 = auxevol_gfs[IDX4ptS(GAMMADD22GF, ijk)];
 
    metric_face_quantities->gamma_faceDD22 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, ijk)];
    
    metric_quantities->KDD01 = auxevol_gfs[IDX4ptS(KDD01GF, ijk)];
    
    metric_quantities->KDD02 = auxevol_gfs[IDX4ptS(KDD02GF, ijk)];
    
    metric_quantities->KDD12 = auxevol_gfs[IDX4ptS(KDD12GF, ijk)];
    
    metric_quantities->KDD00 = auxevol_gfs[IDX4ptS(KDD00GF, ijk)];
    
    metric_quantities->KDD11 = auxevol_gfs[IDX4ptS(KDD11GF, ijk)];
    
    metric_quantities->KDD22 = auxevol_gfs[IDX4ptS(KDD22GF, ijk)];

    metric_quantities_derivatives->alpha_dD2 = auxevol_gfs[IDX4ptS(ALPHA_DD2GF, ijk)];

    metric_quantities_derivatives->betaU_dD00 = auxevol_gfs[IDX4ptS(BETAU_DD00GF, ijk)];

    metric_quantities_derivatives->betaU_dD10 = auxevol_gfs[IDX4ptS(BETAU_DD10GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD111 = auxevol_gfs[IDX4ptS(GAMMADD_DD111GF, ijk)];

    metric_quantities_derivatives->betaU_dD01 = auxevol_gfs[IDX4ptS(BETAU_DD01GF, ijk)];

    metric_quantities_derivatives->betaU_dD20 = auxevol_gfs[IDX4ptS(BETAU_DD20GF, ijk)];

    metric_quantities_derivatives->alpha_dD1 = auxevol_gfs[IDX4ptS(ALPHA_DD1GF, ijk)];

    metric_quantities_derivatives->betaU_dD02 = auxevol_gfs[IDX4ptS(BETAU_DD02GF, ijk)];

    metric_quantities_derivatives->alpha_dD0 = auxevol_gfs[IDX4ptS(ALPHA_DD0GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD110 = auxevol_gfs[IDX4ptS(GAMMADD_DD110GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD120 = auxevol_gfs[IDX4ptS(GAMMADD_DD120GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD221 = auxevol_gfs[IDX4ptS(GAMMADD_DD221GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD000 = auxevol_gfs[IDX4ptS(GAMMADD_DD000GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD122 = auxevol_gfs[IDX4ptS(GAMMADD_DD122GF, ijk)];

    metric_quantities_derivatives->betaU_dD11 = auxevol_gfs[IDX4ptS(BETAU_DD11GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD022 = auxevol_gfs[IDX4ptS(GAMMADD_DD022GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD222 = auxevol_gfs[IDX4ptS(GAMMADD_DD222GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD112 = auxevol_gfs[IDX4ptS(GAMMADD_DD112GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD012 = auxevol_gfs[IDX4ptS(GAMMADD_DD012GF, ijk)];

    metric_quantities_derivatives->betaU_dD22 = auxevol_gfs[IDX4ptS(BETAU_DD22GF, ijk)];

    metric_quantities_derivatives->betaU_dD21 = auxevol_gfs[IDX4ptS(BETAU_DD21GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD020 = auxevol_gfs[IDX4ptS(GAMMADD_DD020GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD121 = auxevol_gfs[IDX4ptS(GAMMADD_DD121GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD220 = auxevol_gfs[IDX4ptS(GAMMADD_DD220GF, ijk)];

    metric_quantities_derivatives->betaU_dD12 = auxevol_gfs[IDX4ptS(BETAU_DD12GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD011 = auxevol_gfs[IDX4ptS(GAMMADD_DD011GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD021 = auxevol_gfs[IDX4ptS(GAMMADD_DD021GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD010 = auxevol_gfs[IDX4ptS(GAMMADD_DD010GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD002 = auxevol_gfs[IDX4ptS(GAMMADD_DD002GF, ijk)];

    metric_quantities_derivatives->gammaDD_dD001 = auxevol_gfs[IDX4ptS(GAMMADD_DD001GF, ijk)];
}
