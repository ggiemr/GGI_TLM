
POST_PROCESSOR_OBJS=  \
allocate_post_data.o  \
find_output_files.o  \
read_write_time_domain_data.o  \
read_write_frequency_domain_data.o  \
extract_time_domain_data.o  \
plot_time_domain_data.o  \
fourier_transform.o  \
plot_frequency_domain_data.o \
create_surface_animation.o \
create_frequency_domain_animation.o \
create_vector_surface_frequency_domain_animation.o \
create_vector_volume_frequency_domain_animation.o \
extract_mode.o \
combine_frequency_domain_data.o \
sum_frequency_domain_data.o \
frequency_average.o \
oneport.o \
twoport.o \
transmission_cross_section.o \
IELF.o \
S_to_Z.o \
correlation.o \
antenna_factor.o \
fourier_transform_warp.o \
generate_random_time_sequence.o \
frequency_domain_filter.o \
PDF_CDF.o \
time_domain_filter.o \
create_time_domain_vector_animation.o \
generate_far_field_plot.o

default: $(POST_PROCESSOR_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
