#include "./NRPy_basic_defines.h"
/*
 * Read in binary data for primitives
 */
void read_from_binary_file_recons(const char *restrict binary_file, const paramstruct *params, double *restrict auxevol_gfs) {
#include "./set_Cparameters.h"


  const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;

  FILE *infile = fopen(binary_file, "rb");
  
  double correct_magic_number = 1.130814081305130e-9;
  double magic_number1, magic_number2, magic_number3;
  
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
