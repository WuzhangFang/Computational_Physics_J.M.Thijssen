# example makefile for Sun-Os, Linux

#
# If no environment variable XPSDIR is specified, 
# do it here, e.g.
#    XPSDIR = /home/jost/XPS
INCFLAG = -I/usr/include/X11 -I$(XPSDIR)
XLIBS =  -L$(XPSDIR) -L/usr/lib -lxps -lX11 -lm

.f:  
	f77 $< $(XLIBS) -o $@

.c: 
	cc $< $(XLIBS) -o $@

