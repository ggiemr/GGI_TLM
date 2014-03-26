
OPT_OBJS= \
opt_module.o \
optimise.o \
opt_calc_error.o \
opt_bisect.o \
opt_simplex.o \
map.o

default: $(OPT_OBJS)

%.o: %.F90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o
