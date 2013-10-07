CC = gcc
INC =
LIB =
CFLAGS = -O3 -march=native -Wall -Wextra $(INC) $(LIB)
LM = -lm -lfftw3

MPICC = mpicc
MPIINC =
MPILIB =
MPICFLAGS = -DMPI -O3 -Wall $(MPIINC) $(MPILIB)
MPILM = -lm -lfftw3

TICC = gcc
TIINC =
TILIB =
TICFLAGS = -DFFTW_THREADS -O3 -Wall -Wextra $(TIINC) $(TILIB)
TILM = -lm -lfftw3 -lfftw3_threads -lpthread


PICC = gcc
PIINC =
PILIB =
PICFLAGS = -DFFTW_THREADS -O3 -Wall -Wextra $(PIINC) $(PILIB)
PILM = -lm -lfftw3 -lfftw3_threads -lpthread -pg



# common code
sources = jh_get.c jh_grid.c jh_linalg.c jh_util.c jh_print.c nrutil.c

targets = 3drism analyze pot plot
targets_mpi = pot_mpi
targets_thr = 3drism_thr
targets_gprof = 3drism_gprof
indeps = bincvt

all: $(targets) $(targets_mpi) $(targets_thr) $(indeps) $(targets_gprof)

$(targets) : % : %.c $(sources)
	$(CC) $(CFLAGS) -o $@ $^ $(LM)

$(targets_mpi) : %_mpi : %.c $(sources)
	$(MPICC) $(MPICFLAGS) -o $@ $^ $(MPILM)

$(targets_thr) : %_thr : %.c $(sources)
	$(TICC) $(TICFLAGS) -o $@ $^ $(TILM)

$(targets_gprof) : %_gprof : %.c $(sources)
	$(PICC) $(PICFLAGS) -o $@ $^ $(PILM)

bincvt : bincvt.c binio.h
	$(CC) $(CFLAGS) -o $@ $<

pack:
	tar -cvf code_pac.tar *.c *.h Makefile test

clean:
	$(RM) *.o a.out *~ $(targets) $(targets_thr) $(targets_mpi) $(indeps) $(targets_gprof)

chmod:
	chmod a-x *.[ch] Makefile *_3drism *.env README

copy:
	cp 3drism analyze pot plot pot_mpi bincvt 3drism_thr test
