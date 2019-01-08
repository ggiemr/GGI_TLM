GERBER_CONVERSION_OBJS=  \
flash.o  \
set_aperture.o \
line.o \
set_aperture_pattern.o \
read_coordinate_and_operation.o   \
set_region.o \
read_format_scale_dcode.o  \
triangulate_surface.o \
read_numbers.o

default: $(GERBER_CONVERSION_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
