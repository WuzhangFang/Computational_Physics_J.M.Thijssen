# COMPHY directory path:
COMPHY	=	.
FC 	= 	gfortran
FFLAGS 	= 	-O
LDFLAGS = 	-O
GEN_OBJ	=	numinhom.o interpol.o
LIBS = -L$(COMPHY)/lib -lnumer

#hydro: hydro.o $(GEN_OBJ) globDFT
#	$(FC) $(FFLAGS) -o hydro hydro.o $(GEN_OBJ) $(LIBS)
test: test.o $(GEN_OBJ) globDFT
	$(FC) $(FFLAGS) -o test test.o $(GEN_OBJ) $(LIBS)

hartree: hartree.o $(GEN_OBJ) globDFT
	$(FC) $(FFLAGS) -o hartree hartree.o $(GEN_OBJ) $(LIBS)

exch: exch.o $(GEN_OBJ) globDFT
	$(FC) $(FFLAGS) -o exch exch.o $(GEN_OBJ) $(LIBS)

excorr: excorr.o $(GEN_OBJ) globDFT
	$(FC) $(FFLAGS) -o excorr excorr.o $(GEN_OBJ) $(LIBS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<
clean:
	rm -f hydro hartree excorr exch *.o
