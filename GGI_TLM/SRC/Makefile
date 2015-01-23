SHELL = /bin/sh

#
# Makefile NOTES:
#
# GGI_TLM can be made for both sequenctial or parallel (mpi) operation
#
# This makefile controls the make process for the sequential and parallel modes, ensuring that
# the make process operates correctly.
#
# The sequential and parallel versions can be made explicitly with the commands:
# 
# for parallel (mpi) solver only
# 'make PARALLEL_MODE=MPI'
#
# for sequential solver and all supporting codes
# 'make PARALLEL_MODE=SEQ'
#
# 'make clean' cleans everything except executables
# 'make clean_all' cleans everything including executables
#
# note that you must 'make clean' before switching from making MPI to SEQ or vice-versa
#
# In addition the single command:
#
# make
#
# Has the following action:
# If the code has been newly installed or cleaned then the sequenctial version is made by default
# If there is a currently made version then this version is re-made
#
# The command:
#
# make all
# makes both the sequential and the mpi versions of the code.
#
# You may need to edit the file: ../make_configuration_data
# in order to set the compiler name, flags etc.
#
ifeq ($(wildcard ./most_recently_made_version), )
# file dose not exist so set LAST_MAKE to empty variable
  LAST_MAKE=
else
# File exists
  include ./most_recently_made_version
endif

#
ifeq ($(PARALLEL_MODE),MPI)
#
# Explicitly make the MPI version
#
  ifeq ($(LAST_MAKE),MPI)
action1:
	@echo "Explicitly make MPI, last make MPI"
	@echo "LAST_MAKE=$(LAST_MAKE)"
        make -f Makefile_GGI_TLM PARALLEL_MODE=MPI
  else
action2:
	@echo "Explicitly make MPI, last make *not* MPI"
	@echo "LAST_MAKE=$(LAST_MAKE)"
	make -f Makefile_GGI_TLM clean
	make -f Makefile_GGI_TLM PARALLEL_MODE=MPI
  endif
  
else ifeq ($(PARALLEL_MODE),SEQ)
#
# Explicitly make the SEQ version
#
  ifeq ($(LAST_MAKE),SEQ)
action3:
	@echo "Explicitly make SEQ, last make SEQ"
	@echo "LAST_MAKE=$(LAST_MAKE)"
        make -f Makefile_GGI_TLM PARALLEL_MODE=SEQ
  else
action4:
	@echo "Explicitly make SEQ, last make *not* SEQ"
	@echo "LAST_MAKE=$(LAST_MAKE)"
	@echo "Issuing command: make -f Makefile_GGI_TLM clean"
	make -f Makefile_GGI_TLM clean
	make -f Makefile_GGI_TLM PARALLEL_MODE=SEQ
  endif
else
#
# Re-make the most recently made version or the SEQ version if the most recent version is not specified
#
  ifeq ($(LAST_MAKE),MPI)
action5:
	@echo "re-make MPI"
	@echo "LAST_MAKE=$(LAST_MAKE)"
        make -f Makefile_GGI_TLM PARALLEL_MODE=MPI
  else ifeq ($(LAST_MAKE),SEQ)
action6:
	@echo "re-make SEQ"
	@echo "LAST_MAKE=$(LAST_MAKE)"
	make -f Makefile_GGI_TLM PARALLEL_MODE=SEQ
  else 
action7:
	@echo "make SEQ as default"
	@echo "LAST_MAKE=${LAST_MAKE}"
	make -f Makefile_GGI_TLM clean
	make -f Makefile_GGI_TLM PARALLEL_MODE=SEQ
  endif

endif

all:
	make -f Makefile_GGI_TLM clean
	make -f Makefile_GGI_TLM PARALLEL_MODE=MPI
	make -f Makefile_GGI_TLM clean
	make -f Makefile_GGI_TLM PARALLEL_MODE=SEQ

clean:
	make -f Makefile_GGI_TLM clean

clean_all:
	make -f Makefile_GGI_TLM clean_all

