# COMPHY directory path:
COMPHY	=	.
FC 	= 	gfortran
FOPTS	=	-O -c
LDOPTS	=	-O
LIBS	=	./lib/geneig.a -llapack

pseudo: pseudo.o 
	$(FC) $(LDOPTS) -o pseudo pseudo.o $(LIBS)

pseudo.o: pseudo.f pseudo.glb
	$(FC) $(FOPTS) pseudo.f
clean:
	rm -f pseudo *.o
