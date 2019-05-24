#CF=ifort
CF=gfortran
#CC=icc
CC=gcc
FFLAGS=-g -O3 -march=native -mtune=native
BIN=../bin

all: test_orig test_maqao

wallclock.o: wallclock.c
	$(CC) -c $^ -o $@

test_orig: opernlb_ylm_pp.f90 test_choice1_cmplex2_paw3_nofunction.F90 wallclock.o
	$(CF) $(FFLAGS) $^ -o $(BIN)/$@

test_maqao: maqao_opernlb_ylm_pp_ref.f90 test_choice1_cmplex2_paw3_nofunction.F90 wallclock.o
	$(CF) $(FFLAGS) $^ -o $(BIN)/$@

test_maqao2:maqao_opernlb_ylm_pp.f90 test_choice1_cmplex2_paw3_nofunction.F90 wallclock.o
        $(CF) $(FFLAGS) $^ -o $(BIN)/$@


clean:
	rm *.o
