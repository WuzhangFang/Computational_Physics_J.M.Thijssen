FC=gfortran
FFLAGS = -o 
COMPHY=.
LIBS = -L -lnumer -llapack

apw: apw.o special.o apwglob
	$(FC) $(FFLAGS)  -o apw apw.o special.o $(LIBS)

#logapw: logapw.o special.o Pot.glob apwglob
#	$(FC) $(FFLAGS)  -o logapw logapw.o special.o $(LIBS)


%.o: %f
	$(FC) $(FFLAGS) $< 

clean:
	rm -f *.o apw logapw
