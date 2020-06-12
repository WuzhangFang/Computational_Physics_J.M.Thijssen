FC	= 	gfortran
FOPTS	= 	-c -O
CC	= 	gcc
COPTS	= 	-c -O

libs: src/matrix/geneig.f
	cd src/matrix; $(FC) $(FOPTS) geneig.f; ar r ../../geneig.a geneig.o; cd ../..; 
	cd src/numerov; $(FC) $(FOPTS) numerov.f; ar r ../../libnumer.a numerov.o; cd ../..; 
	cd src/rangen; $(FC) $(FOPTS) rangen.f;	ar r ../../libran.a rangen.o;  cd ../..
	
clean:
	rm *.a; cd src/matrix; rm *.o; cd ../..;
	        cd src/numerov; rm *.o; cd ../..;
	        cd src/rangen; rm *.o; cd ../..;
