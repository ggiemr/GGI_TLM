
default: Mlayer

MLAYER_OBJS=  Mlayer_module.o Mlayer.o ABCD_calc.o   calc_Sparams.o \
       file_stuff.o memory_stuff.o

# make at top level
Mlayer: $(MLAYER_OBJS)
	$(FC) $(FLAGS) -o Mlayer $(MLAYER_OBJS) $(GGI_LIBS) $(LIBS)
	mv Mlayer $(EXECUTABLE_DIR)/GGI_Mlayer
# 
Mlayer_module.o: 	Mlayer_module.F90
	$(FC) $(FLAGS) -c Mlayer_module.F90 
# 
Mlayer.o: 	Mlayer.F90  Mlayer_module.o
	$(FC) $(FLAGS) -c Mlayer.F90 
# 
ABCD_calc.o: 	ABCD_calc.F90  Mlayer_module.o
	$(FC) $(FLAGS) -c ABCD_calc.F90 
# 
calc_Sparams.o: 	calc_Sparams.F90  Mlayer_module.o
	$(FC) $(FLAGS) -c calc_Sparams.F90 
# 
file_stuff.o: 	file_stuff.F90  Mlayer_module.o
	$(FC) $(FLAGS) -c file_stuff.F90 
# 
memory_stuff.o: 	memory_stuff.F90  Mlayer_module.o
	$(FC) $(FLAGS) -c memory_stuff.F90 
#		
clean:
	rm -f *.o
	rm -f *.mod
	rm -f Mlayer
# 
