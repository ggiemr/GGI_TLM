
FILTER_FIT_OBJS=  \
read_filter_fit_information.o \
read_data.o \
filter_fit_memory_stuff.o \
get_initial_filter_coefficients.o \
get_testing_frequencies.o \
optimise_filter.o \
filter_param_conversion.o \
stabilise_FF_input_data.o \
stabilise_filter.o \
stabilise_material_filter.o \
stabilise_thin_layer_filter.o \
stabilise_impedance_filter.o \
test_stability.o \
write_FF_input_data.o \
write_filter_frequency_response.o \
write_filter.o \
calculate_mean_square_error.o

default: $(FILTER_FIT_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
