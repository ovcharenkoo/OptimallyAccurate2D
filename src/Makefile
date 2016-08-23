#BIN = ../../bin
BIN = ../bin/
program = $(BIN)/multiOPT.x 
program2 = $(BIN)/modelgen.x
program3 = $(BIN)/frechetKernel.x
include ../etc/config.h
objs =   modules.o others.o vectors.o disconConfig.o seismicSources.o videos.o   cales.o cales_discon.o cales_free.o sub.inverse.o calmat.o OPT2D.o
objs2 = modelgenerator.o
objs3 = frechetKernel.o modules.o others.o vectors.o videos.o 

all : $(program) $(program2) $(program3)

.SUFFIXES: .o .f90

$(program3): $(objs3)
	$(FC) -o $@ $(objs3) $(FFLAGS)

$(program2): $(objs2)
	$(FC)  -o $@ $(objs2) $(FFLAGS)

$(program): $(objs)
	$(FC)  -o $@ $(objs) $(FFLAGS)
.f90.o:
	$(FC) -o $@ -c $< $(FFLAGS)

.PHONY: clean
clean:
	rm $(program) $(program2) $(objs) $(objs2) $(objs3) *.optrpt  

