
OUTPUT_OBJS=  \
output.o \
output_points.o \
output_surfaces.o \
output_volumes.o \
output_volume_averages.o \
output_modes.o \
frequency_output_surfaces.o \
frequency_output_volumes.o \
frequency_domain_power_surfaces.o \
far_field_surfaces.o \
RCS_surfaces.o \
SAR_volumes.o \
cable_outputs.o

default: $(OUTPUT_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o

