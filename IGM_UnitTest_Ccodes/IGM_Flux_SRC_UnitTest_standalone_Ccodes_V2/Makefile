CC     = gcc
CFLAGS =  -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas -std=gnu99
#CFLAGS =  -O2 -march=native -g -fopenmp -Wall -Wno-unused-variable -std=gnu99
#CFLAGS =  -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas -std=gnu99
#CFLAGS =   -O2 -march=native -fopenmp -Wall -Wno-unused-variable -std=gnu99
INCLUDEDIRS = 
all: IGM_standalone  calculate_metric_derivatives.o calculate_all_source_terms.o calculate_characteristic_speed_0th_direction.o calculate_characteristic_speed_1th_direction.o calculate_characteristic_speed_2th_direction.o calculate_HLLE_fluxes0.o calculate_HLLE_fluxes1.o calculate_HLLE_fluxes2.o initialize_structs.o read_from_binary_file_all.o read_from_binary_file_recons.o set_Nxx_dxx_invdx_params__and__xx.o set_Cparameters_to_default.o main.o
calculate_metric_derivatives.o: calculate_metric_derivatives.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_metric_derivatives.c -o calculate_metric_derivatives.o

calculate_all_source_terms.o: calculate_all_source_terms.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_all_source_terms.c -o calculate_all_source_terms.o

calculate_characteristic_speed_0th_direction.o: calculate_characteristic_speed_0th_direction.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_characteristic_speed_0th_direction.c -o calculate_characteristic_speed_0th_direction.o

calculate_characteristic_speed_1th_direction.o: calculate_characteristic_speed_1th_direction.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_characteristic_speed_1th_direction.c -o calculate_characteristic_speed_1th_direction.o

calculate_characteristic_speed_2th_direction.o: calculate_characteristic_speed_2th_direction.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_characteristic_speed_2th_direction.c -o calculate_characteristic_speed_2th_direction.o

calculate_HLLE_fluxes0.o: calculate_HLLE_fluxes0.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_HLLE_fluxes0.c -o calculate_HLLE_fluxes0.o

calculate_HLLE_fluxes1.o: calculate_HLLE_fluxes1.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_HLLE_fluxes1.c -o calculate_HLLE_fluxes1.o

calculate_HLLE_fluxes2.o: calculate_HLLE_fluxes2.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c calculate_HLLE_fluxes2.c -o calculate_HLLE_fluxes2.o

initialize_structs.o: initialize_structs.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c initialize_structs.c -o initialize_structs.o

read_from_binary_file_all.o: read_from_binary_file_all.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c read_from_binary_file_all.c -o read_from_binary_file_all.o

read_from_binary_file_recons.o: read_from_binary_file_recons.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c read_from_binary_file_recons.c -o read_from_binary_file_recons.o

set_Nxx_dxx_invdx_params__and__xx.o: set_Nxx_dxx_invdx_params__and__xx.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c set_Nxx_dxx_invdx_params__and__xx.c -o set_Nxx_dxx_invdx_params__and__xx.o

set_Cparameters_to_default.o: set_Cparameters_to_default.c
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c set_Cparameters_to_default.c -o set_Cparameters_to_default.o

main.o: main.c free_parameters.h
	$(CC) $(CFLAGS) $(INCLUDEDIRS)  -c main.c -o main.o

IGM_standalone:   calculate_metric_derivatives.o calculate_all_source_terms.o calculate_characteristic_speed_0th_direction.o calculate_characteristic_speed_1th_direction.o calculate_characteristic_speed_2th_direction.o calculate_HLLE_fluxes0.o calculate_HLLE_fluxes1.o calculate_HLLE_fluxes2.o initialize_structs.o read_from_binary_file_all.o read_from_binary_file_recons.o set_Nxx_dxx_invdx_params__and__xx.o set_Cparameters_to_default.o main.o
	$(CC)   calculate_metric_derivatives.o calculate_all_source_terms.o calculate_characteristic_speed_0th_direction.o calculate_characteristic_speed_1th_direction.o calculate_characteristic_speed_2th_direction.o calculate_HLLE_fluxes0.o calculate_HLLE_fluxes1.o calculate_HLLE_fluxes2.o initialize_structs.o read_from_binary_file_all.o read_from_binary_file_recons.o  set_Nxx_dxx_invdx_params__and__xx.o set_Cparameters_to_default.o main.o -o IGM_standalone -lm

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png IGM_standalone
