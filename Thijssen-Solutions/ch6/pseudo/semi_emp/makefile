# COMPHY directory path:
#COMPHY	=	$(HOME)/comphy
FC 	= 	gfortran
FOPTS	=	-O -c
LDOPTS	=	-O
LIBS	=	./geneig.a ./liblapack.a

pseudo: pseudo.o 
	$(FC) $(LDOPTS) -o pseudo pseudo.o $(LIBS)

pseudo.o: pseudo.f pseudo.glb
	$(FC) $(FOPTS) pseudo.f
clean:
	rm -f pseudo *.o
