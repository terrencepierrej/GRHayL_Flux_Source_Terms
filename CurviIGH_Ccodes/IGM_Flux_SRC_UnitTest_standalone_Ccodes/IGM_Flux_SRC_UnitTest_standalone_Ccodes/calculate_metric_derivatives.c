#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Calculate metric derivatives
 */
void calculate_metric_derivatives(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict auxevol_gfs) {
#include "./set_Cparameters.h"

  #pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
        const REAL xx0 = xx[0][i0];
        {
          /*
           * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
           */
          const double alpha_i0_i1_i2m1 = auxevol_gfs[IDX4S(ALPHAGF, i0,i1,i2-1)];
          const double alpha_i0_i1m1_i2 = auxevol_gfs[IDX4S(ALPHAGF, i0,i1-1,i2)];
          const double alpha_i0m1_i1_i2 = auxevol_gfs[IDX4S(ALPHAGF, i0-1,i1,i2)];
          const double alpha_i0p1_i1_i2 = auxevol_gfs[IDX4S(ALPHAGF, i0+1,i1,i2)];
          const double alpha_i0_i1p1_i2 = auxevol_gfs[IDX4S(ALPHAGF, i0,i1+1,i2)];
          const double alpha_i0_i1_i2p1 = auxevol_gfs[IDX4S(ALPHAGF, i0,i1,i2+1)];
          const double gammaDD00_i0_i1_i2m1 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1,i2-1)];
          const double gammaDD00_i0_i1m1_i2 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1-1,i2)];
          const double gammaDD00_i0m1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD00GF, i0-1,i1,i2)];
          const double gammaDD00_i0p1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD00GF, i0+1,i1,i2)];
          const double gammaDD00_i0_i1p1_i2 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1+1,i2)];
          const double gammaDD00_i0_i1_i2p1 = auxevol_gfs[IDX4S(GAMMADD00GF, i0,i1,i2+1)];
          const double gammaDD01_i0_i1_i2m1 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1,i2-1)];
          const double gammaDD01_i0_i1m1_i2 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1-1,i2)];
          const double gammaDD01_i0m1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD01GF, i0-1,i1,i2)];
          const double gammaDD01_i0p1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD01GF, i0+1,i1,i2)];
          const double gammaDD01_i0_i1p1_i2 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1+1,i2)];
          const double gammaDD01_i0_i1_i2p1 = auxevol_gfs[IDX4S(GAMMADD01GF, i0,i1,i2+1)];
          const double gammaDD02_i0_i1_i2m1 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1,i2-1)];
          const double gammaDD02_i0_i1m1_i2 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1-1,i2)];
          const double gammaDD02_i0m1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD02GF, i0-1,i1,i2)];
          const double gammaDD02_i0p1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD02GF, i0+1,i1,i2)];
          const double gammaDD02_i0_i1p1_i2 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1+1,i2)];
          const double gammaDD02_i0_i1_i2p1 = auxevol_gfs[IDX4S(GAMMADD02GF, i0,i1,i2+1)];
          const double gammaDD11_i0_i1_i2m1 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1,i2-1)];
          const double gammaDD11_i0_i1m1_i2 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1-1,i2)];
          const double gammaDD11_i0m1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD11GF, i0-1,i1,i2)];
          const double gammaDD11_i0p1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD11GF, i0+1,i1,i2)];
          const double gammaDD11_i0_i1p1_i2 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1+1,i2)];
          const double gammaDD11_i0_i1_i2p1 = auxevol_gfs[IDX4S(GAMMADD11GF, i0,i1,i2+1)];
          const double gammaDD12_i0_i1_i2m1 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1,i2-1)];
          const double gammaDD12_i0_i1m1_i2 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1-1,i2)];
          const double gammaDD12_i0m1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD12GF, i0-1,i1,i2)];
          const double gammaDD12_i0p1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD12GF, i0+1,i1,i2)];
          const double gammaDD12_i0_i1p1_i2 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1+1,i2)];
          const double gammaDD12_i0_i1_i2p1 = auxevol_gfs[IDX4S(GAMMADD12GF, i0,i1,i2+1)];
          const double gammaDD22_i0_i1_i2m1 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1,i2-1)];
          const double gammaDD22_i0_i1m1_i2 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1-1,i2)];
          const double gammaDD22_i0m1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD22GF, i0-1,i1,i2)];
          const double gammaDD22_i0p1_i1_i2 = auxevol_gfs[IDX4S(GAMMADD22GF, i0+1,i1,i2)];
          const double gammaDD22_i0_i1p1_i2 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1+1,i2)];
          const double gammaDD22_i0_i1_i2p1 = auxevol_gfs[IDX4S(GAMMADD22GF, i0,i1,i2+1)];
          const double betaU0_i0_i1_i2m1 = auxevol_gfs[IDX4S(BETAU0GF, i0,i1,i2-1)];
          const double betaU0_i0_i1m1_i2 = auxevol_gfs[IDX4S(BETAU0GF, i0,i1-1,i2)];
          const double betaU0_i0m1_i1_i2 = auxevol_gfs[IDX4S(BETAU0GF, i0-1,i1,i2)];
          const double betaU0_i0p1_i1_i2 = auxevol_gfs[IDX4S(BETAU0GF, i0+1,i1,i2)];
          const double betaU0_i0_i1p1_i2 = auxevol_gfs[IDX4S(BETAU0GF, i0,i1+1,i2)];
          const double betaU0_i0_i1_i2p1 = auxevol_gfs[IDX4S(BETAU0GF, i0,i1,i2+1)];
          const double betaU1_i0_i1_i2m1 = auxevol_gfs[IDX4S(BETAU1GF, i0,i1,i2-1)];
          const double betaU1_i0_i1m1_i2 = auxevol_gfs[IDX4S(BETAU1GF, i0,i1-1,i2)];
          const double betaU1_i0m1_i1_i2 = auxevol_gfs[IDX4S(BETAU1GF, i0-1,i1,i2)];
          const double betaU1_i0p1_i1_i2 = auxevol_gfs[IDX4S(BETAU1GF, i0+1,i1,i2)];
          const double betaU1_i0_i1p1_i2 = auxevol_gfs[IDX4S(BETAU1GF, i0,i1+1,i2)];
          const double betaU1_i0_i1_i2p1 = auxevol_gfs[IDX4S(BETAU1GF, i0,i1,i2+1)];
          const double betaU2_i0_i1_i2m1 = auxevol_gfs[IDX4S(BETAU2GF, i0,i1,i2-1)];
          const double betaU2_i0_i1m1_i2 = auxevol_gfs[IDX4S(BETAU2GF, i0,i1-1,i2)];
          const double betaU2_i0m1_i1_i2 = auxevol_gfs[IDX4S(BETAU2GF, i0-1,i1,i2)];
          const double betaU2_i0p1_i1_i2 = auxevol_gfs[IDX4S(BETAU2GF, i0+1,i1,i2)];
          const double betaU2_i0_i1p1_i2 = auxevol_gfs[IDX4S(BETAU2GF, i0,i1+1,i2)];
          const double betaU2_i0_i1_i2p1 = auxevol_gfs[IDX4S(BETAU2GF, i0,i1,i2+1)];
          const double FDPart1_Rational_1_2 = 1.0/2.0;
          const double FDPart1_0 = FDPart1_Rational_1_2*invdx0;
          const double FDPart1_1 = FDPart1_Rational_1_2*invdx1;
          const double FDPart1_2 = FDPart1_Rational_1_2*invdx2;
          const double alpha_dD0 = FDPart1_0*(-alpha_i0m1_i1_i2 + alpha_i0p1_i1_i2);
          const double alpha_dD1 = FDPart1_1*(-alpha_i0_i1m1_i2 + alpha_i0_i1p1_i2);
          const double alpha_dD2 = FDPart1_2*(-alpha_i0_i1_i2m1 + alpha_i0_i1_i2p1);
          const double betaU_dD00 = FDPart1_0*(-betaU0_i0m1_i1_i2 + betaU0_i0p1_i1_i2);
          const double betaU_dD01 = FDPart1_1*(-betaU0_i0_i1m1_i2 + betaU0_i0_i1p1_i2);
          const double betaU_dD02 = FDPart1_2*(-betaU0_i0_i1_i2m1 + betaU0_i0_i1_i2p1);
          const double betaU_dD10 = FDPart1_0*(-betaU1_i0m1_i1_i2 + betaU1_i0p1_i1_i2);
          const double betaU_dD11 = FDPart1_1*(-betaU1_i0_i1m1_i2 + betaU1_i0_i1p1_i2);
          const double betaU_dD12 = FDPart1_2*(-betaU1_i0_i1_i2m1 + betaU1_i0_i1_i2p1);
          const double betaU_dD20 = FDPart1_0*(-betaU2_i0m1_i1_i2 + betaU2_i0p1_i1_i2);
          const double betaU_dD21 = FDPart1_1*(-betaU2_i0_i1m1_i2 + betaU2_i0_i1p1_i2);
          const double betaU_dD22 = FDPart1_2*(-betaU2_i0_i1_i2m1 + betaU2_i0_i1_i2p1);
          const double gammaDD_dD000 = FDPart1_0*(-gammaDD00_i0m1_i1_i2 + gammaDD00_i0p1_i1_i2);
          const double gammaDD_dD001 = FDPart1_1*(-gammaDD00_i0_i1m1_i2 + gammaDD00_i0_i1p1_i2);
          const double gammaDD_dD002 = FDPart1_2*(-gammaDD00_i0_i1_i2m1 + gammaDD00_i0_i1_i2p1);
          const double gammaDD_dD010 = FDPart1_0*(-gammaDD01_i0m1_i1_i2 + gammaDD01_i0p1_i1_i2);
          const double gammaDD_dD011 = FDPart1_1*(-gammaDD01_i0_i1m1_i2 + gammaDD01_i0_i1p1_i2);
          const double gammaDD_dD012 = FDPart1_2*(-gammaDD01_i0_i1_i2m1 + gammaDD01_i0_i1_i2p1);
          const double gammaDD_dD020 = FDPart1_0*(-gammaDD02_i0m1_i1_i2 + gammaDD02_i0p1_i1_i2);
          const double gammaDD_dD021 = FDPart1_1*(-gammaDD02_i0_i1m1_i2 + gammaDD02_i0_i1p1_i2);
          const double gammaDD_dD022 = FDPart1_2*(-gammaDD02_i0_i1_i2m1 + gammaDD02_i0_i1_i2p1);
          const double gammaDD_dD110 = FDPart1_0*(-gammaDD11_i0m1_i1_i2 + gammaDD11_i0p1_i1_i2);
          const double gammaDD_dD111 = FDPart1_1*(-gammaDD11_i0_i1m1_i2 + gammaDD11_i0_i1p1_i2);
          const double gammaDD_dD112 = FDPart1_2*(-gammaDD11_i0_i1_i2m1 + gammaDD11_i0_i1_i2p1);
          const double gammaDD_dD120 = FDPart1_0*(-gammaDD12_i0m1_i1_i2 + gammaDD12_i0p1_i1_i2);
          const double gammaDD_dD121 = FDPart1_1*(-gammaDD12_i0_i1m1_i2 + gammaDD12_i0_i1p1_i2);
          const double gammaDD_dD122 = FDPart1_2*(-gammaDD12_i0_i1_i2m1 + gammaDD12_i0_i1_i2p1);
          const double gammaDD_dD220 = FDPart1_0*(-gammaDD22_i0m1_i1_i2 + gammaDD22_i0p1_i1_i2);
          const double gammaDD_dD221 = FDPart1_1*(-gammaDD22_i0_i1m1_i2 + gammaDD22_i0_i1p1_i2);
          const double gammaDD_dD222 = FDPart1_2*(-gammaDD22_i0_i1_i2m1 + gammaDD22_i0_i1_i2p1);
          /*
           * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
           */
          auxevol_gfs[IDX4S(ALPHA_DD0GF, i0, i1, i2)] = alpha_dD0;
          auxevol_gfs[IDX4S(BETAU_DD00GF, i0, i1, i2)] = betaU_dD00;
          auxevol_gfs[IDX4S(BETAU_DD01GF, i0, i1, i2)] = betaU_dD01;
          auxevol_gfs[IDX4S(BETAU_DD02GF, i0, i1, i2)] = betaU_dD02;
          auxevol_gfs[IDX4S(ALPHA_DD1GF, i0, i1, i2)] = alpha_dD1;
          auxevol_gfs[IDX4S(BETAU_DD10GF, i0, i1, i2)] = betaU_dD10;
          auxevol_gfs[IDX4S(BETAU_DD11GF, i0, i1, i2)] = betaU_dD11;
          auxevol_gfs[IDX4S(BETAU_DD12GF, i0, i1, i2)] = betaU_dD12;
          auxevol_gfs[IDX4S(ALPHA_DD2GF, i0, i1, i2)] = alpha_dD2;
          auxevol_gfs[IDX4S(BETAU_DD20GF, i0, i1, i2)] = betaU_dD20;
          auxevol_gfs[IDX4S(BETAU_DD21GF, i0, i1, i2)] = betaU_dD21;
          auxevol_gfs[IDX4S(BETAU_DD22GF, i0, i1, i2)] = betaU_dD22;
          auxevol_gfs[IDX4S(GAMMADD_DD000GF, i0, i1, i2)] = gammaDD_dD000;
          auxevol_gfs[IDX4S(GAMMADD_DD001GF, i0, i1, i2)] = gammaDD_dD001;
          auxevol_gfs[IDX4S(GAMMADD_DD002GF, i0, i1, i2)] = gammaDD_dD002;
          auxevol_gfs[IDX4S(GAMMADD_DD010GF, i0, i1, i2)] = gammaDD_dD010;
          auxevol_gfs[IDX4S(GAMMADD_DD011GF, i0, i1, i2)] = gammaDD_dD011;
          auxevol_gfs[IDX4S(GAMMADD_DD012GF, i0, i1, i2)] = gammaDD_dD012;
          auxevol_gfs[IDX4S(GAMMADD_DD020GF, i0, i1, i2)] = gammaDD_dD020;
          auxevol_gfs[IDX4S(GAMMADD_DD021GF, i0, i1, i2)] = gammaDD_dD021;
          auxevol_gfs[IDX4S(GAMMADD_DD022GF, i0, i1, i2)] = gammaDD_dD022;
          auxevol_gfs[IDX4S(GAMMADD_DD110GF, i0, i1, i2)] = gammaDD_dD110;
          auxevol_gfs[IDX4S(GAMMADD_DD111GF, i0, i1, i2)] = gammaDD_dD111;
          auxevol_gfs[IDX4S(GAMMADD_DD112GF, i0, i1, i2)] = gammaDD_dD112;
          auxevol_gfs[IDX4S(GAMMADD_DD120GF, i0, i1, i2)] = gammaDD_dD120;
          auxevol_gfs[IDX4S(GAMMADD_DD121GF, i0, i1, i2)] = gammaDD_dD121;
          auxevol_gfs[IDX4S(GAMMADD_DD122GF, i0, i1, i2)] = gammaDD_dD122;
          auxevol_gfs[IDX4S(GAMMADD_DD220GF, i0, i1, i2)] = gammaDD_dD220;
          auxevol_gfs[IDX4S(GAMMADD_DD221GF, i0, i1, i2)] = gammaDD_dD221;
          auxevol_gfs[IDX4S(GAMMADD_DD222GF, i0, i1, i2)] = gammaDD_dD222;
        }
      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
