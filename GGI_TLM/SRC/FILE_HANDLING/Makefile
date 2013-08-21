
FILE_HANDLING_OBJS=  \
file_stuff.o \
read_problem_name.o \
read_write_header_data.o 

default: $(FILE_HANDLING_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
