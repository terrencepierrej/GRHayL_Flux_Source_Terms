#include "./NRPy_basic_defines.h"
#include "./NRPy_function_prototypes.h"
/*
 * Adds source term and connection terms to Stilde, rho_star and tau_tilde
 */
void calculate_all_source_terms(const reconstructed_prims_struct *restrict reconstructed_prims, const metric_quantities_struct *restrict metric_quantities, const metric_quantities_derivatives_struct *restrict metric_quantities_derivatives, conservative_sources_struct *restrict conservative_sources) {

{
const double u4U0 = reconstructed_prims->u4U0;
const double h = reconstructed_prims->h;
const double BU0 = reconstructed_prims->BU0;
const double P = reconstructed_prims->P;
const double BU2 = reconstructed_prims->BU2;
const double rhob = reconstructed_prims->rhob;
const double u4U1 = reconstructed_prims->u4U1;
const double u4U3 = reconstructed_prims->u4U3;
const double BU1 = reconstructed_prims->BU1;
const double u4U2 = reconstructed_prims->u4U2;
const double alpha = metric_quantities->alpha;
const double betaU0 = metric_quantities->betaU0;
const double betaU1 = metric_quantities->betaU1;
const double betaU2 = metric_quantities->betaU2;
const double KDD00 = metric_quantities->KDD00;
const double KDD01 = metric_quantities->KDD01;
const double KDD02 = metric_quantities->KDD02;
const double KDD11 = metric_quantities->KDD11;
const double KDD12 = metric_quantities->KDD12;
const double KDD22 = metric_quantities->KDD22;
const double gammaDD00 = metric_quantities->gammaDD00;
const double gammaDD01 = metric_quantities->gammaDD01;
const double gammaDD02 = metric_quantities->gammaDD02;
const double gammaDD11 = metric_quantities->gammaDD11;
const double gammaDD12 = metric_quantities->gammaDD12;
const double gammaDD22 = metric_quantities->gammaDD22;
const double alpha_dD2 = metric_quantities_derivatives->alpha_dD2;
const double betaU_dD00 = metric_quantities_derivatives->betaU_dD00;
const double betaU_dD10 = metric_quantities_derivatives->betaU_dD10;
const double gammaDD_dD111 = metric_quantities_derivatives->gammaDD_dD111;
const double betaU_dD01 = metric_quantities_derivatives->betaU_dD01;
const double betaU_dD20 = metric_quantities_derivatives->betaU_dD20;
const double alpha_dD1 = metric_quantities_derivatives->alpha_dD1;
const double betaU_dD02 = metric_quantities_derivatives->betaU_dD02;
const double alpha_dD0 = metric_quantities_derivatives->alpha_dD0;
const double gammaDD_dD110 = metric_quantities_derivatives->gammaDD_dD110;
const double gammaDD_dD120 = metric_quantities_derivatives->gammaDD_dD120;
const double gammaDD_dD221 = metric_quantities_derivatives->gammaDD_dD221;
const double gammaDD_dD000 = metric_quantities_derivatives->gammaDD_dD000;
const double gammaDD_dD122 = metric_quantities_derivatives->gammaDD_dD122;
const double betaU_dD11 = metric_quantities_derivatives->betaU_dD11;
const double gammaDD_dD022 = metric_quantities_derivatives->gammaDD_dD022;
const double gammaDD_dD222 = metric_quantities_derivatives->gammaDD_dD222;
const double gammaDD_dD112 = metric_quantities_derivatives->gammaDD_dD112;
const double gammaDD_dD012 = metric_quantities_derivatives->gammaDD_dD012;
const double betaU_dD22 = metric_quantities_derivatives->betaU_dD22;
const double betaU_dD21 = metric_quantities_derivatives->betaU_dD21;
const double gammaDD_dD020 = metric_quantities_derivatives->gammaDD_dD020;
const double gammaDD_dD121 = metric_quantities_derivatives->gammaDD_dD121;
const double gammaDD_dD220 = metric_quantities_derivatives->gammaDD_dD220;
const double betaU_dD12 = metric_quantities_derivatives->betaU_dD12;
const double gammaDD_dD011 = metric_quantities_derivatives->gammaDD_dD011;
const double gammaDD_dD021 = metric_quantities_derivatives->gammaDD_dD021;
const double gammaDD_dD010 = metric_quantities_derivatives->gammaDD_dD010;
const double gammaDD_dD002 = metric_quantities_derivatives->gammaDD_dD002;
const double gammaDD_dD001 = metric_quantities_derivatives->gammaDD_dD001;
  const double tmp_0 = betaU0*gammaDD00 + betaU1*gammaDD01 + betaU2*gammaDD02;
  const double tmp_1 = betaU0*gammaDD01 + betaU1*gammaDD11 + betaU2*gammaDD12;
  const double tmp_2 = betaU0*gammaDD02 + betaU1*gammaDD12 + betaU2*gammaDD22;
  const double tmp_3 = BU0*(gammaDD00*u4U1 + gammaDD01*u4U2 + gammaDD02*u4U3 + tmp_0*u4U0) + BU1*(gammaDD01*u4U1 + gammaDD11*u4U2 + gammaDD12*u4U3 + tmp_1*u4U0) + BU2*(gammaDD02*u4U1 + gammaDD12*u4U2 + gammaDD22*u4U3 + tmp_2*u4U0);
  const double tmp_4 = BU0 + tmp_3*u4U1;
  const double tmp_7 = (1.0/((alpha)*(alpha)));
  const double tmp_8 = tmp_7/((sqrt4pi)*(sqrt4pi));
  const double tmp_9 = tmp_8/((u4U0)*(u4U0));
  const double tmp_10 = ((tmp_4)*(tmp_4))*tmp_9;
  const double tmp_12 = BU1 + tmp_3*u4U2;
  const double tmp_13 = ((tmp_12)*(tmp_12))*tmp_9;
  const double tmp_15 = BU2 + tmp_3*u4U3;
  const double tmp_16 = ((tmp_15)*(tmp_15))*tmp_9;
  const double tmp_18 = ((tmp_3)*(tmp_3))*tmp_8;
  const double tmp_19 = tmp_18*(-((alpha)*(alpha)) + betaU0*tmp_0 + betaU1*tmp_1 + betaU2*tmp_2);
  const double tmp_21 = tmp_12*tmp_4*tmp_9;
  const double tmp_22 = tmp_15*tmp_4*tmp_9;
  const double tmp_24 = tmp_12*tmp_15*tmp_9;
  const double tmp_26 = tmp_3*tmp_8/u4U0;
  const double tmp_27 = tmp_26*tmp_4;
  const double tmp_29 = tmp_12*tmp_26;
  const double tmp_31 = tmp_15*tmp_26;
  const double tmp_33 = gammaDD00*tmp_10 + 2*gammaDD01*tmp_21 + 2*gammaDD02*tmp_22 + gammaDD11*tmp_13 + 2*gammaDD12*tmp_24 + gammaDD22*tmp_16 + h*rhob + 2*tmp_0*tmp_27 + 2*tmp_1*tmp_29 + tmp_19 + 2*tmp_2*tmp_31;
  const double tmp_38 = gammaDD00*gammaDD11*gammaDD22 - gammaDD00*((gammaDD12)*(gammaDD12)) - ((gammaDD01)*(gammaDD01))*gammaDD22 + 2*gammaDD01*gammaDD02*gammaDD12 - ((gammaDD02)*(gammaDD02))*gammaDD11;
  const double tmp_39 = (1.0/(tmp_38));
  const double tmp_40 = P + (1.0/2.0)*gammaDD00*tmp_10 + gammaDD01*tmp_21 + gammaDD02*tmp_22 + (1.0/2.0)*gammaDD11*tmp_13 + gammaDD12*tmp_24 + (1.0/2.0)*gammaDD22*tmp_16 + tmp_0*tmp_27 + tmp_1*tmp_29 + (1.0/2.0)*tmp_19 + tmp_2*tmp_31;
  const double tmp_41 = -tmp_10 + tmp_33*((u4U1)*(u4U1)) + tmp_40*(-((betaU0)*(betaU0))*tmp_7 + tmp_39*(gammaDD11*gammaDD22 - ((gammaDD12)*(gammaDD12))));
  const double tmp_42 = alpha*sqrt(tmp_38);
  const double tmp_43 = (1.0/2.0)*tmp_42;
  const double tmp_44 = tmp_41*tmp_43;
  const double tmp_46 = -tmp_13 + tmp_33*((u4U2)*(u4U2)) + tmp_40*(-((betaU1)*(betaU1))*tmp_7 + tmp_39*(gammaDD00*gammaDD22 - ((gammaDD02)*(gammaDD02))));
  const double tmp_47 = tmp_43*tmp_46;
  const double tmp_49 = -tmp_16 + tmp_33*((u4U3)*(u4U3)) + tmp_40*(-((betaU2)*(betaU2))*tmp_7 + tmp_39*(gammaDD00*gammaDD11 - ((gammaDD01)*(gammaDD01))));
  const double tmp_50 = tmp_43*tmp_49;
  const double tmp_51 = betaU0*gammaDD_dD000 + betaU1*gammaDD_dD010 + betaU2*gammaDD_dD020 + betaU_dD00*gammaDD00 + betaU_dD10*gammaDD01 + betaU_dD20*gammaDD02;
  const double tmp_52 = tmp_40*tmp_7;
  const double tmp_54 = tmp_33*u4U1;
  const double tmp_56 = betaU0*tmp_52 - tmp_27 + tmp_54*u4U0;
  const double tmp_57 = tmp_42*tmp_56;
  const double tmp_58 = betaU0*gammaDD_dD010 + betaU1*gammaDD_dD110 + betaU2*gammaDD_dD120 + betaU_dD00*gammaDD01 + betaU_dD10*gammaDD11 + betaU_dD20*gammaDD12;
  const double tmp_61 = tmp_33*u4U0*u4U2;
  const double tmp_62 = betaU1*tmp_52 - tmp_29 + tmp_61;
  const double tmp_63 = tmp_42*tmp_62;
  const double tmp_64 = betaU0*gammaDD_dD020 + betaU1*gammaDD_dD120 + betaU2*gammaDD_dD220 + betaU_dD00*gammaDD02 + betaU_dD10*gammaDD12 + betaU_dD20*gammaDD22;
  const double tmp_66 = tmp_33*u4U0*u4U3;
  const double tmp_67 = betaU2*tmp_52 - tmp_31 + tmp_66;
  const double tmp_68 = tmp_42*tmp_67;
  const double tmp_69 = 2*alpha;
  const double tmp_70 = -tmp_18 + tmp_33*((u4U0)*(u4U0)) - tmp_52;
  const double tmp_71 = tmp_43*tmp_70;
  const double tmp_73 = -tmp_21 + tmp_40*(-betaU0*betaU1*tmp_7 + tmp_39*(-gammaDD01*gammaDD22 + gammaDD02*gammaDD12)) + tmp_54*u4U2;
  const double tmp_74 = tmp_42*tmp_73;
  const double tmp_75 = -tmp_22 + tmp_40*(-betaU0*betaU2*tmp_7 + tmp_39*(gammaDD01*gammaDD12 - gammaDD02*gammaDD11)) + tmp_54*u4U3;
  const double tmp_76 = tmp_42*tmp_75;
  const double tmp_77 = -tmp_24 + tmp_33*u4U2*u4U3 + tmp_40*(-betaU1*betaU2*tmp_7 + tmp_39*(-gammaDD00*gammaDD12 + gammaDD01*gammaDD02));
  const double tmp_78 = tmp_42*tmp_77;
  const double tmp_79 = betaU0*gammaDD_dD001 + betaU1*gammaDD_dD011 + betaU2*gammaDD_dD021 + betaU_dD01*gammaDD00 + betaU_dD11*gammaDD01 + betaU_dD21*gammaDD02;
  const double tmp_80 = betaU0*gammaDD_dD011 + betaU1*gammaDD_dD111 + betaU2*gammaDD_dD121 + betaU_dD01*gammaDD01 + betaU_dD11*gammaDD11 + betaU_dD21*gammaDD12;
  const double tmp_81 = betaU0*gammaDD_dD021 + betaU1*gammaDD_dD121 + betaU2*gammaDD_dD221 + betaU_dD01*gammaDD02 + betaU_dD11*gammaDD12 + betaU_dD21*gammaDD22;
  const double tmp_82 = betaU0*gammaDD_dD002 + betaU1*gammaDD_dD012 + betaU2*gammaDD_dD022 + betaU_dD02*gammaDD00 + betaU_dD12*gammaDD01 + betaU_dD22*gammaDD02;
  const double tmp_83 = betaU0*gammaDD_dD012 + betaU1*gammaDD_dD112 + betaU2*gammaDD_dD122 + betaU_dD02*gammaDD01 + betaU_dD12*gammaDD11 + betaU_dD22*gammaDD12;
  const double tmp_84 = betaU0*gammaDD_dD022 + betaU1*gammaDD_dD122 + betaU2*gammaDD_dD222 + betaU_dD02*gammaDD02 + betaU_dD12*gammaDD12 + betaU_dD22*gammaDD22;
  const double tmp_85 = betaU0*tmp_70;
  const double tmp_87 = 2*betaU0*tmp_52 - 2*tmp_27 + 2*tmp_54*u4U0;
  const double tmp_88 = 2*betaU1*tmp_52 - 2*tmp_29 + 2*tmp_61;
  const double tmp_89 = 2*betaU2*tmp_52 - 2*tmp_31 + 2*tmp_66;
  const double tmp_90 = betaU1*tmp_85 + tmp_73;
  const double tmp_91 = betaU2*tmp_85 + tmp_75;
  const double tmp_92 = betaU1*betaU2*tmp_70 + tmp_77;
  conservative_sources->StildeD0_src = gammaDD_dD000*tmp_44 + gammaDD_dD010*tmp_74 + gammaDD_dD020*tmp_76 + gammaDD_dD110*tmp_47 + gammaDD_dD120*tmp_78 + gammaDD_dD220*tmp_50 + tmp_51*tmp_57 + tmp_58*tmp_63 + tmp_64*tmp_68 + tmp_71*(-alpha_dD0*tmp_69 + betaU0*tmp_51 + betaU1*tmp_58 + betaU2*tmp_64 + betaU_dD00*tmp_0 + betaU_dD10*tmp_1 + betaU_dD20*tmp_2);
  conservative_sources->StildeD1_src = gammaDD_dD001*tmp_44 + gammaDD_dD011*tmp_74 + gammaDD_dD021*tmp_76 + gammaDD_dD111*tmp_47 + gammaDD_dD121*tmp_78 + gammaDD_dD221*tmp_50 + tmp_57*tmp_79 + tmp_63*tmp_80 + tmp_68*tmp_81 + tmp_71*(-alpha_dD1*tmp_69 + betaU0*tmp_79 + betaU1*tmp_80 + betaU2*tmp_81 + betaU_dD01*tmp_0 + betaU_dD11*tmp_1 + betaU_dD21*tmp_2);
  conservative_sources->StildeD2_src = gammaDD_dD002*tmp_44 + gammaDD_dD012*tmp_74 + gammaDD_dD022*tmp_76 + gammaDD_dD112*tmp_47 + gammaDD_dD122*tmp_78 + gammaDD_dD222*tmp_50 + tmp_57*tmp_82 + tmp_63*tmp_83 + tmp_68*tmp_84 + tmp_71*(-alpha_dD2*tmp_69 + betaU0*tmp_82 + betaU1*tmp_83 + betaU2*tmp_84 + betaU_dD02*tmp_0 + betaU_dD12*tmp_1 + betaU_dD22*tmp_2);
  conservative_sources->tau_tilde_src = tmp_42*(KDD00*(((betaU0)*(betaU0))*tmp_70 + betaU0*tmp_87 + tmp_41) + KDD01*(betaU0*tmp_88 + tmp_90) + KDD01*(betaU1*tmp_87 + tmp_90) + KDD02*(betaU0*tmp_89 + tmp_91) + KDD02*(betaU2*tmp_87 + tmp_91) + KDD11*(((betaU1)*(betaU1))*tmp_70 + betaU1*tmp_88 + tmp_46) + KDD12*(betaU1*tmp_89 + tmp_92) + KDD12*(betaU2*tmp_88 + tmp_92) + KDD22*(((betaU2)*(betaU2))*tmp_70 + betaU2*tmp_89 + tmp_49) + alpha_dD0*(-tmp_56 - tmp_85) + alpha_dD1*(-betaU1*tmp_70 - tmp_62) + alpha_dD2*(-betaU2*tmp_70 - tmp_67));
}
}
