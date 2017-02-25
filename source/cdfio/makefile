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
L=$(LIBDIR)/cdfio.a

MODULES =  \
$(LIBDIR)/cdfio.a(netcdf_int_mod.o) \
$(LIBDIR)/cdfio.a(ncprofile_mod.o)

all:	ezcdf $(MODULES)  


ezcdf: 
	cd ../ezcdf; \
	make

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
