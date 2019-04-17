MODE=SEQ

GGI_LIB_DIR=$(GGI_TLM_DIRECTORY)/GGI_TLM_LIBS

MOD_DIR=$(GGI_TLM_DIRECTORY)/GGI_TLM_MODULE_FILES

include $(GGI_TLM_DIRECTORY)/make_configuration_data

LIBS=$(GGI_LIB_DIR)/lib_GGI_optimisation.a 
NGSPICE_MODE=NO_NGSPICE

optimisation_process: optimisation_process.F90
	$(FC) $(FLAGS) -o optimisation_process optimisation_process.F90 $(LIBS)


