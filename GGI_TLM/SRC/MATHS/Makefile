
MATHS_OBJS=  \
dmatrix.o \
cmatrix.o \
Legendre.o \
bessel.o \
set_random_seed.o

default: $(MATHS_OBJS) GGI_TLM_generate_random_numbers

%.o: %.F90
	$(FC) $(FLAGS) -c $<
	
GGI_TLM_generate_random_numbers: GGI_TLM_generate_random_numbers.F90 
	$(FC) $(FLAGS) -o GGI_TLM_generate_random_numbers set_random_seed.o GGI_TLM_generate_random_numbers.F90
	mv GGI_TLM_generate_random_numbers $(EXECUTABLE_DIR)/GGI_TLM_generate_random_numbers

clean:
	( rm -f *.o ; rm -f *.mod ; rm -f GGI_TLM_generate_random_numbers )
