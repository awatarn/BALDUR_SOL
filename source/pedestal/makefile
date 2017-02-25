# makefile for the BALDUR code
# Alexei Pankin, 2002
# pankin@fusion.physics.lehigh.edu
#
include ../../Make.flags
include ../../Make.local
#
ifdef DEBUG
	FFLAGS = $(DFFLAGS)
	CFLAGS = $(DCFLAGS)
	BINDIR = $(DEBDIR)
	LIBDIR = $(DBLDIR)
endif
#
L=$(LIBDIR)/pedestal.a

MODULES =  \
$(L)(bdboots.o) \
$(L)(bdcoll.o) \
$(L)(bdhden.o) \
$(L)(bdheight1.o) \
$(L)(bdheight2.o) \
$(L)(bdheight1.o) \
$(L)(bdheight2.o) \
$(L)(bdheight3.o) \
$(L)(bdheight4.o) \
$(L)(bdheight5.o) \
$(L)(bdheight6.o) \
$(L)(bdheight11.o) \
$(L)(bdheight12.o) \
$(L)(bdhtemp.o) \
$(L)(bdlhmode.o) \
$(L)(bdshear1.o) \
$(L)(bdshear2.o) \
$(L)(bdshear3.o) \
$(L)(lhtemode.o)

all:	$(MODULES)


.SUFFIXES:
.SUFFIXES: .o .f90 .a .F90 .for .f .F

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o

.F.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o

.f90.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o

.f90.o:
	$(FC) $(FFLAGS) -c $<
