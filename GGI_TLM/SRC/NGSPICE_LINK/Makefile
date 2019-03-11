
NGSPICE_LINK_OBJS=  \
ngspice_wrapper.o \
run_ngspice_TLM.o

default: $(NGSPICE_LINK_OBJS)

run_ngspice_TLM.o: run_ngspice_TLM.F90
	$(FC) $(FLAGS) -c run_ngspice_TLM.F90 

ngspice_wrapper.o: ngspice_wrapper.c
	gcc -c ngspice_wrapper.c

clean:
	rm -f *.o
