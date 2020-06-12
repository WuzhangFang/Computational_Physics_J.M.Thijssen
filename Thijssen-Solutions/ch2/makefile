# COMPHY: path of COMPHY directory
COMPHY	=	$(HOME)/comphy
FC	=	gfortran
LFLAGS 	= 	-O -L$(COMPHY)/lib -lnumer
FFLAGS 	= 	-c -O

scatter: scatter.o special.o 
	$(FC) $(LFLAGS) scatter.o special.o -o scatter $(LFLAGS)

scatter.o: scatter.f globscatter
	$(FC) $(FFLAGS) scatter.f 

special.o: special.f 
	$(FC) $(FFLAGS) special.f 

numerov.o: numerov.f globscatter
	$(FC) $(FFLAGS) numerov.f 

clean:
	rm -f *.o scatter
