CC = gcc
INC =
LIB =
CFLAGS = -O3 -march=native -Wall $(INC) $(LIB)
PORF = -pg
LM = -lm -lfftw3

ICC = icc
IINC =
ILIB =
ICFLAGS = -O3 -march=native -ipo $(INC) $(LIB)
ILM = -lfftw3

TICC = icc
TIINC =
TILIB =
TICFLAGS = -DFFTW_THREADS -xhost -O3 -Wall -Wextra $(TIINC) $(TILIB)
TILM = -lm -lfftw3 -lfftw3_threads -lpthread

OMPICC = icc
INC =
LIB =
OMPCFLAGS = -openmp -xhost -O3 -Wall -Wextra $(INC) $(LIB)
##OMPCFLAGS = -DOMP -xhost -O3 -openmp -Wall -Wextra $(INC) $(LIB)
OMPLM = -lfftw3 -lfftw3_omp

IPROF = -profile-functions -profile-loops=all -profile-loops-report=2

# common code
sources = jh_get.c jh_grid.c jh_linalg.c jh_util.c jh_print.c nrutil.c

targets_omp_prof = 3drism_omp_prof
targets_omp = 3drism_omp 
targets = 3drism analyze pot plot
targets_thr = 3drism_thr
targets_icc = 3drism_icc
targets_icc_prof = 3drism_icc_prof
targets_gcc_prof = 3drism_gcc_prof

all: $(targets_omp)
##$(targets) $(targets_mpi) $(targets_thr) $(targets_icc) $(targets_icc_prof) $(targets_gcc_prof)

$(targets_omp) : %_omp : %.c $(sources)
	$(OMPICC) $(OMPCFLAGS) -o $@ $^ $(OMPLM)

$(targets) : % : %.c $(sources)
	$(CC) $(CFLAGS) -o $@ $^ $(LM)

$(targets_thr) : %_thr : %.c $(sources)
	$(TICC) $(TICFLAGS) -o $@ $^ $(TILM)

$(targets_icc) : %_icc : %.c $(sources)
	$(ICC) $(ICFLAGS) -o $@ $^ $(ILM)

$(targets_icc_prof) : %_icc_prof : %.c $(sources)
	$(ICC) $(ICFLAGS) $(IPROF) -o $@ $^ $(ILM)

$(targets_gcc_prof) : %_gcc_prof : %.c $(sources)
	$(CC) $(CFLAGS) $(PROF) -o $@ $^ $(LM)

pack:
	tar -cvf code_pac.tar *.c *.h Makefile test test-normal

clean:
	$(RM) *.o a.out *~ $(targets) $(targets_thr) $(targets_omp)

chmod:
	chmod a-x *.[ch] Makefile *_3drism *.env README

copy:
	cp 3drism_omp test 	
##cp 3drism analyze pot plot pot_mpi 3drism_thr 3drism_icc 3drism_icc_prof 3drism_gcc_prof test
