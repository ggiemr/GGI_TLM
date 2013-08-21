
EISPACK_OBJS=  \
eispack.o 

#
default: $(EISPACK_OBJS)
#
%.o: %.f90
	$(FC) $(FLAGS) -c $<

clean:
	rm -f *.o


