
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
fourier_transform_warp.o 

default: $(POST_PROCESSOR_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
