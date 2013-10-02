

objects = nrutil.o JH_util.o JH_linalg.o JH_fourier.o \

VER=2.1

CC = gcc
MPCC = mpicc
CFLAGS =  -O3
CINC = 
CLIB = 
SYM = x

#CFLAGS = -fpic -Wall -O3

sysname := $(shell uname -s)
nodename := $(shell uname -n)


#export MACHINE_NAME
ifeq ($(nodename),atlantis)
         CC = gcc
     CFLAGS = -O3
       CINC = -I$(HOME)/include
       CLIB = -L$(HOME)/lib
        SYM = atl 
      MACHINE_NAME = Linux
endif


ifeq ($(nodename),F.Chem.UH.EDU)
         CC = gcc
     CFLAGS = -O3
       CINC = -I$(HOME)/include
       CLIB = -L$(HOME)/lib
        SYM = f 
      MACHINE_NAME = Linux
endif


ifeq ($(nodename),flexo.tlc2.uh.edu)
         CC = gcc
       #MPCC = icc
       MPCC = mpicc
     CFLAGS = -O3
       CINC = -I$(HOME)/include -I/opt/openmpi/intel/1.3/include
       CLIB = -L$(HOME)/lib -L/opt/openmpi/intel/1.3/lib
      #MPLIB = -lmpi 
      MPLIB =  
        SYM = flexo 
      MACHINE_NAME = Linux
endif


ifeq ($(nodename),nibbler.tlc2.uh.edu)
         CC = gcc
       #MPCC = icc
       MPCC = mpicc
     CFLAGS = -O3
       CINC = -I$(HOME)/include
       CLIB = -L$(HOME)/lib
      #MPLIB = -lmpi 
      MPLIB =  
        SYM = nib
      MACHINE_NAME = redhat-lam
endif


ifeq ($(nodename),Np-eth.Chem.UH.EDU)
         CC = gcc
       MPCC = mpicc	
     CFLAGS = -O3
    #   CINC = -I$(HOME)/np/include
     #  CLIB = -L$(HOME)/np/lib 
      MPLIB =  
        SYM = np
      MACHINE_NAME = redhat-lam
endif


ifeq ($(nodename),howard-home)
         CC = gcc
     CFLAGS = -O3
       CINC = -I$(HOME)/include
       CLIB = -L$(HOME)/lib 
        SYM = hh
      MACHINE_NAME = cygwin
endif


ifeq ($(nodename),Tesla2)
         CC = gcc
     CFLAGS = -O3
       CINC = -I$(HOME)/include
       CLIB = -L$(HOME)/lib
        SYM = t 
      MACHINE_NAME = Cygwin
endif


ifeq ($(MACHINE_NAME),redhat-gnu)
          CC = gcc
          FC = g77
      FFLAGS = -O3 -ffast-math -malign-double
      CFLAGS = -O3 -ffast-math
    CPPFLAGS =
      NO_MPI = 1
endif


#################################################
#3drism: $(objects)

jh_objects = jh_get.o jh_grid.o jh_linalg.o jh_util.o jh_print.o 


3drism: 3drism.o $(jh_objects) nrutil.o 
	$(CC) $(CFLAGS) $(CINC) $(CLIB) -o 3drism_$(SYM) 3drism.o $(jh_objects) nrutil.o  -lm -lfftw3 
	rm *.o

3drism_xx: 3drism_xx.o $(jh_objects) nrutil.o 
	$(CC) $(CFLAGS) $(CINC) $(CLIB) -o 3drism_xx_$(SYM) 3drism_xx.o $(jh_objects) nrutil.o  -lm -lfftw3 
	rm *.o

pot: pot.o $(jh_objects) nrutil.o 
	$(CC) $(CFLAGS) $(CINC) $(CLIB) -o pot_$(SYM) pot.o $(jh_objects) nrutil.o -lm -lfftw3
	rm *.o

analyze: analyze.o $(jh_objects) nrutil.o 
	$(CC) $(CFLAGS) $(CINC) $(CLIB) -o analyze_$(SYM) analyze.o $(jh_objects) nrutil.o -lm -lfftw3 
	rm *.o

convert_grid: convert_grid.o jh_get.o
	$(CC) $(CFLAGS) $(CINC) -o convert_grid convert_grid.o jh_get.o -lm
	rm *.o

edit_env: edit_env.o jh_get.o
	$(CC) $(CFLAGS) $(CINC) -o edit_env_$(SYM) edit_env.o jh_get.o -lm
	rm *.o

potmp:
	cat potmp.h pot.c > potmp.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_get.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_linalg.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_grid.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_util.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_print.c
	$(MPCC) $(CFLAGS) $(CINC) -c nrutil.c
	$(MPCC) $(CFLAGS) $(CINC) -c potmp.c
	$(MPCC) $(CFLAGS) $(CINC) $(CLIB) -o potmp_$(SYM) potmp.o jh_get.o jh_grid.o jh_linalg.o jh_util.o jh_print.o nrutil.o $(MPLIB) -lm -lfftw3
	rm *.o
	rm potmp.c

3drismmp:
	cat potmp.h 3drism.c > 3drismmp.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_get.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_linalg.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_grid.c
	$(MPCC) $(CFLAGS) $(CINC) -c jh_util.c
	$(MPCC) $(CFLAGS) $(CINC) -c nrutil.c
	$(MPCC) $(CFLAGS) $(CINC) -c 3drismmp.c
	$(MPCC) $(CFLAGS) $(CINC) $(CLIB) -o 3drismmp_$(SYM) 3drismmp.o jh_get.o jh_grid.o jh_linalg.o jh_util.o nrutil.o $(MPLIB) -lm -lfftw3
	rm *.o
	rm 3drismmp.c


plot: jh_get.o
	$(CC) $(CFLAGS) $(CINC) -o plot_$(SYM) plot.c jh_get.o -lm
	rm *.o

##############################################################3

analyze.o: analyze.c
	$(CC) $(CFLAGS) $(CINC) -c analyze.c

3drism.o: 3drism.c
	$(CC) $(CFLAGS) $(CINC) -c 3drism.c

3drism_xx.o: 3drism_xx.c
	$(CC) $(CFLAGS) $(CINC) -c 3drism_xx.c

pot.o: pot.c
	$(CC) $(CFLAGS) $(CINC) -c pot.c

jh_get.o: jh_get.c
	$(CC) $(CFLAGS) $(CINC) -c jh_get.c

jh_grid.o: jh_grid.c
	$(CC) $(CFLAGS) $(CINC) -c jh_grid.c

jh_util.o: jh_util.c
	$(CC) $(CFLAGS) $(CINC) -c jh_util.c

jh_linalg.o: jh_linalg.c
	$(CC) $(CFLAGS) $(CINC) -c jh_linalg.c

jh_print.o: jh_print.c
	$(CC) $(CFLAGS) $(CINC) -c jh_print.c

nrutil.o: nrutil.c
	$(CC) $(CFLAGS) $(CINC) -c nrutil.c

convert_grid.o: convert_grid.c
	$(CC) $(CFLAGS) $(CINC) -c convert_grid.c

##############################################################3
##############################################################3

code_pac:
	tar -cvf code_pac.tar *.c *.h makefile test

$(VER):  
	cp template.env $(HOME)/code/template_$(VER).env
	cp pot_$(SYM) $(HOME)/code/pot_$(VER)_$(SYM)
	cp potmp_$(SYM) $(HOME)/code/potmp_$(VER)_$(SYM)
	cp 3drism_$(SYM) $(HOME)/code/3drism_$(VER)_$(SYM)
	cp analyze_$(SYM) $(HOME)/code/analyze_$(VER)_$(SYM)
	cp plot_$(SYM) $(HOME)/code/plot_$(VER)_$(SYM)

$(VER)-sync:
	make pot
	make potmp
	make 3drism
	make analyze
	make plot  
	cp template.env $(HOME)/code/template_$(VER).env
	cp pot_$(SYM) $(HOME)/code/pot_$(VER)_$(SYM)
	cp potmp_$(SYM) $(HOME)/code/potmp_$(VER)_$(SYM)
	cp 3drism_$(SYM) $(HOME)/code/3drism_$(VER)_$(SYM)
	cp analyze_$(SYM) $(HOME)/code/analyze_$(VER)_$(SYM)
	cp plot_$(SYM) $(HOME)/code/plot_$(VER)_$(SYM)

all:
	make pot
	make potmp
	make 3drism
	make analyze
	make plot  

