# Makefile for Car-Parrinello Program

COMPILER        = gfortran
LDOPTS		= -O
COMPOPTS	= -c -O   
LIBS		= -L/usr/local/lib -lfftw3 -lblas -llapack
OBJECTS		= fft.o globals.o utilities.o  \
	          excorr.o pseudo.o KS.o main.o

%.o : %.f90
	$(COMPILER) $(COMPOPTS) $<

%.o : %.f
	$(COMPILER) $(COMPOPTS) -I/usr/local/include $<

scpseud: $(OBJECTS)
	$(COMPILER) $(OBJECTS) $(LDOPTS) -o scpseud -L/usr/local/lib/ -lblas -llapack $(LIBS) 

clean: 
	rm -f *.o *.d *.mod *~

cleanall:
	rm -f *.o *.d *.mod scpseud last.dat work*.* *~
