#include "./NRPy_basic_defines.h"
/*
 * Read in binary data for metric quantities and primitives
 */
void read_from_binary_file_all(const char *restrict binary_file, const paramstruct *params, double *restrict auxevol_gfs) {
#include "./set_Cparameters.h"


  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  FILE *infile = fopen(binary_file, "rb");

  double correct_magic_number = 1.130814081305130e-9;
  double magic_number1, magic_number2, magic_number3, magic_number4, magic_number5;
  
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
