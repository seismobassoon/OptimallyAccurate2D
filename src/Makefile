#BIN = ../../bin
BIN = ../bin/
program0 = $(BIN)/multiOPT_nodiscon.x
program = $(BIN)/multiOPT.x 
program2 = $(BIN)/modelgen.x
program3 = $(BIN)/frechetKernel.x
program4 = $(BIN)/FWI2D.x

include ../etc/config.h

objs0 =   modules.o others.o vectors.o  seismicSources.o videos.o   cales.o  OPT2D_nodiscon.o
objs =   modules.o others.o vectors.o disconConfig.o seismicSources.o videos.o   cales.o cales_discon.o cales_free.o sub.inverse.o calmat.o OPT2D.o
objs2 = modelgenerator4Anais.o
#objs2 = modelgenerator.o

objs3 = frechetKernel.o modules.o others.o vectors.o videos.o 
objs4 = forwardmodelling.o backpropagation.o gradientCalculation.o modules.o others.o vectors.o disconConfig.o seismicSources.o videos.o   cales.o cales_discon.o cales_free.o sub.inverse.o calmat.o FWI2D.o
all : $(program) $(program2) $(program3) $(program4)

.SUFFIXES: .o .f90

$(program4): $(objs4)
	$(FC) -o $@ $(objs4) $(FFLAGS)
$(program0): $(objs0)
	$(FC) -o $@ $(objs0) $(FFLAGS)

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
	rm $(program) $(program2) $(program3) $(program0)  $(objs0)  $(objs) $(objs2) $(objs3) *.optrpt  

