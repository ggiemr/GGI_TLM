
F90_GENERAL_OBJS=  \
convert_case.o \
strip_filename.o \
read_write_subroutines.o \
error_subroutines.o \
add_integer_to_filename.o \
write_time.o

default: $(F90_GENERAL_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
