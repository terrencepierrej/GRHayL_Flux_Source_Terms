#include "./NRPy_basic_defines.h"
/*
 * Initialize structs used for flux and source calculations
 */
void initialize_structs(const int idx, const paramstruct *params, const double *restrict auxevol_gfs, reconstructed_prims_struct *restrict reconstructed_prims, reconstructed_prims_struct *restrict reconstructed_prims_r, reconstructed_prims_struct *restrict reconstructed_prims_l, metric_quantities_struct *restrict metric_quantities, metric_face_quantities_struct *restrict metric_face_quantities, metric_quantities_derivatives_struct *restrict metric_quantities_derivatives) {
#include "./set_Cparameters.h"


    
    reconstructed_prims->u4U0 = auxevol_gfs[IDX4ptS(U4U0GF, idx)];
    reconstructed_prims_r->u4U0 = auxevol_gfs[IDX4ptS(U4_RU0GF, idx)];
    reconstructed_prims_l->u4U0 = auxevol_gfs[IDX4ptS(U4_LU0GF, idx)];
    
    reconstructed_prims->u4U1 = auxevol_gfs[IDX4ptS(U4U1GF, idx)];
    reconstructed_prims_r->u4U1 = auxevol_gfs[IDX4ptS(U4_RU1GF, idx)];
    reconstructed_prims_l->u4U1 = auxevol_gfs[IDX4ptS(U4_LU1GF, idx)];
    
    reconstructed_prims->u4U2 = auxevol_gfs[IDX4ptS(U4U2GF, idx)];
    reconstructed_prims_r->u4U2 = auxevol_gfs[IDX4ptS(U4_RU2GF, idx)];
    reconstructed_prims_l->u4U2 = auxevol_gfs[IDX4ptS(U4_LU2GF, idx)];
    
    reconstructed_prims->u4U3 = auxevol_gfs[IDX4ptS(U4U3GF, idx)];
    reconstructed_prims_r->u4U3 = auxevol_gfs[IDX4ptS(U4_RU3GF, idx)];
    reconstructed_prims_l->u4U3 = auxevol_gfs[IDX4ptS(U4_LU3GF, idx)];
    
    reconstructed_prims->BU0 = auxevol_gfs[IDX4ptS(BU0GF, idx)];
    reconstructed_prims_r->BU0 = auxevol_gfs[IDX4ptS(B_RU0GF, idx)];
    reconstructed_prims_l->BU0 = auxevol_gfs[IDX4ptS(B_LU0GF, idx)];
    
    reconstructed_prims->BU1 = auxevol_gfs[IDX4ptS(BU1GF, idx)];
    reconstructed_prims_r->BU1 = auxevol_gfs[IDX4ptS(B_RU1GF, idx)];
    reconstructed_prims_l->BU1 = auxevol_gfs[IDX4ptS(B_LU1GF, idx)];
    
    reconstructed_prims->BU2 = auxevol_gfs[IDX4ptS(BU2GF, idx)];
    reconstructed_prims_r->BU2 = auxevol_gfs[IDX4ptS(B_RU2GF, idx)];
    reconstructed_prims_l->BU2 = auxevol_gfs[IDX4ptS(B_LU2GF, idx)];

    reconstructed_prims_r->P = auxevol_gfs[IDX4ptS(P_RGF, idx)];
    reconstructed_prims_l->P = auxevol_gfs[IDX4ptS(P_LGF, idx)];

    reconstructed_prims->P = auxevol_gfs[IDX4ptS(PGF, idx)];

    reconstructed_prims_r->h = auxevol_gfs[IDX4ptS(H_RGF, idx)];
    reconstructed_prims_l->h = auxevol_gfs[IDX4ptS(H_LGF, idx)];

    reconstructed_prims->h = auxevol_gfs[IDX4ptS(HGF, idx)];

    reconstructed_prims_r->rhob = auxevol_gfs[IDX4ptS(RHOB_RGF, idx)];
    reconstructed_prims_l->rhob = auxevol_gfs[IDX4ptS(RHOB_LGF, idx)];

    reconstructed_prims->rhob = auxevol_gfs[IDX4ptS(RHOBGF, idx)];

    reconstructed_prims_r->Gamma_th = auxevol_gfs[IDX4ptS(GAMMA_TH_RGF, idx)];
    reconstructed_prims_l->Gamma_th = auxevol_gfs[IDX4ptS(GAMMA_TH_LGF, idx)];

    reconstructed_prims_r->epsilon_th = auxevol_gfs[IDX4ptS(EPSILON_TH_RGF, idx)];
    reconstructed_prims_l->epsilon_th = auxevol_gfs[IDX4ptS(EPSILON_TH_LGF, idx)];

    reconstructed_prims_r->dPcold_drhob = auxevol_gfs[IDX4ptS(DPCOLD_DRHOB_RGF, idx)];
    reconstructed_prims_l->dPcold_drhob = auxevol_gfs[IDX4ptS(DPCOLD_DRHOB_LGF, idx)];
    
    metric_quantities->alpha = auxevol_gfs[IDX4ptS(ALPHAGF, idx)];
 
    metric_face_quantities->alpha_face = auxevol_gfs[IDX4ptS(ALPHA_FACEGF, idx)];
    
    metric_quantities->betaU0 = auxevol_gfs[IDX4ptS(BETAU0GF, idx)];
 
    metric_face_quantities->beta_faceU0 = auxevol_gfs[IDX4ptS(BETA_FACEU0GF, idx)];
    
    metric_quantities->betaU1 = auxevol_gfs[IDX4ptS(BETAU1GF, idx)];
 
    metric_face_quantities->beta_faceU1 = auxevol_gfs[IDX4ptS(BETA_FACEU1GF, idx)];
    
    metric_quantities->betaU2 = auxevol_gfs[IDX4ptS(BETAU2GF, idx)];
 
    metric_face_quantities->beta_faceU2 = auxevol_gfs[IDX4ptS(BETA_FACEU2GF, idx)];
    
    metric_quantities->gammaDD01 = auxevol_gfs[IDX4ptS(GAMMADD01GF, idx)];
 
    metric_face_quantities->gamma_faceDD01 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD01GF, idx)];
    
    metric_quantities->gammaDD02 = auxevol_gfs[IDX4ptS(GAMMADD02GF, idx)];
 
    metric_face_quantities->gamma_faceDD02 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD02GF, idx)];
    
    metric_quantities->gammaDD12 = auxevol_gfs[IDX4ptS(GAMMADD12GF, idx)];
 
    metric_face_quantities->gamma_faceDD12 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD12GF, idx)];
    
    metric_quantities->gammaDD00 = auxevol_gfs[IDX4ptS(GAMMADD00GF, idx)];
 
    metric_face_quantities->gamma_faceDD00 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD00GF, idx)];
    
    metric_quantities->gammaDD11 = auxevol_gfs[IDX4ptS(GAMMADD11GF, idx)];
 
    metric_face_quantities->gamma_faceDD11 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD11GF, idx)];
    
    metric_quantities->gammaDD22 = auxevol_gfs[IDX4ptS(GAMMADD22GF, idx)];
 
    metric_face_quantities->gamma_faceDD22 = auxevol_gfs[IDX4ptS(GAMMA_FACEDD22GF, idx)];
    
    metric_quantities->KDD01 = auxevol_gfs[IDX4ptS(KDD01GF, idx)];
    
    metric_quantities->KDD02 = auxevol_gfs[IDX4ptS(KDD02GF, idx)];
    
    metric_quantities->KDD12 = auxevol_gfs[IDX4ptS(KDD12GF, idx)];
    
    metric_quantities->KDD00 = auxevol_gfs[IDX4ptS(KDD00GF, idx)];
    
    metric_quantities->KDD11 = auxevol_gfs[IDX4ptS(KDD11GF, idx)];
    
    metric_quantities->KDD22 = auxevol_gfs[IDX4ptS(KDD22GF, idx)];

    metric_quantities_derivatives->alpha_dD1 = auxevol_gfs[IDX4ptS(ALPHA_DD1GF, idx)];

    metric_quantities_derivatives->gammaDD_dD111 = auxevol_gfs[IDX4ptS(GAMMADD_DD111GF, idx)];

    metric_quantities_derivatives->alpha_dD2 = auxevol_gfs[IDX4ptS(ALPHA_DD2GF, idx)];

    metric_quantities_derivatives->betaU_dD12 = auxevol_gfs[IDX4ptS(BETAU_DD12GF, idx)];

    metric_quantities_derivatives->gammaDD_dD020 = auxevol_gfs[IDX4ptS(GAMMADD_DD020GF, idx)];

    metric_quantities_derivatives->gammaDD_dD122 = auxevol_gfs[IDX4ptS(GAMMADD_DD122GF, idx)];

    metric_quantities_derivatives->gammaDD_dD001 = auxevol_gfs[IDX4ptS(GAMMADD_DD001GF, idx)];

    metric_quantities_derivatives->gammaDD_dD022 = auxevol_gfs[IDX4ptS(GAMMADD_DD022GF, idx)];

    metric_quantities_derivatives->gammaDD_dD121 = auxevol_gfs[IDX4ptS(GAMMADD_DD121GF, idx)];

    metric_quantities_derivatives->alpha_dD0 = auxevol_gfs[IDX4ptS(ALPHA_DD0GF, idx)];

    metric_quantities_derivatives->betaU_dD10 = auxevol_gfs[IDX4ptS(BETAU_DD10GF, idx)];

    metric_quantities_derivatives->betaU_dD01 = auxevol_gfs[IDX4ptS(BETAU_DD01GF, idx)];

    metric_quantities_derivatives->gammaDD_dD110 = auxevol_gfs[IDX4ptS(GAMMADD_DD110GF, idx)];

    metric_quantities_derivatives->betaU_dD00 = auxevol_gfs[IDX4ptS(BETAU_DD00GF, idx)];

    metric_quantities_derivatives->betaU_dD22 = auxevol_gfs[IDX4ptS(BETAU_DD22GF, idx)];

    metric_quantities_derivatives->gammaDD_dD221 = auxevol_gfs[IDX4ptS(GAMMADD_DD221GF, idx)];

    metric_quantities_derivatives->gammaDD_dD112 = auxevol_gfs[IDX4ptS(GAMMADD_DD112GF, idx)];

    metric_quantities_derivatives->gammaDD_dD021 = auxevol_gfs[IDX4ptS(GAMMADD_DD021GF, idx)];

    metric_quantities_derivatives->betaU_dD21 = auxevol_gfs[IDX4ptS(BETAU_DD21GF, idx)];

    metric_quantities_derivatives->betaU_dD20 = auxevol_gfs[IDX4ptS(BETAU_DD20GF, idx)];

    metric_quantities_derivatives->gammaDD_dD012 = auxevol_gfs[IDX4ptS(GAMMADD_DD012GF, idx)];

    metric_quantities_derivatives->betaU_dD02 = auxevol_gfs[IDX4ptS(BETAU_DD02GF, idx)];

    metric_quantities_derivatives->gammaDD_dD222 = auxevol_gfs[IDX4ptS(GAMMADD_DD222GF, idx)];

    metric_quantities_derivatives->gammaDD_dD220 = auxevol_gfs[IDX4ptS(GAMMADD_DD220GF, idx)];

    metric_quantities_derivatives->gammaDD_dD010 = auxevol_gfs[IDX4ptS(GAMMADD_DD010GF, idx)];

    metric_quantities_derivatives->betaU_dD11 = auxevol_gfs[IDX4ptS(BETAU_DD11GF, idx)];

    metric_quantities_derivatives->gammaDD_dD002 = auxevol_gfs[IDX4ptS(GAMMADD_DD002GF, idx)];

    metric_quantities_derivatives->gammaDD_dD000 = auxevol_gfs[IDX4ptS(GAMMADD_DD000GF, idx)];

    metric_quantities_derivatives->gammaDD_dD120 = auxevol_gfs[IDX4ptS(GAMMADD_DD120GF, idx)];

    metric_quantities_derivatives->gammaDD_dD011 = auxevol_gfs[IDX4ptS(GAMMADD_DD011GF, idx)];
}
