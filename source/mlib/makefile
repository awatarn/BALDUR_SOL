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
L=$(LIBDIR)/mlib.a

MODULES =  \
$(L)(aitken.o) \
$(L)(amaxaf.o) \
$(L)(cubint.o) \
$(L)(dintlib.o) \
$(L)(dode.o) \
$(L)(fbes.o) \
$(L)(int1d.o) \
$(L)(ismax.o) \
$(L)(fft1.o) \
$(L)(lug.o) \
$(L)(matrx1.o) \
$(L)(matrx2.o) \
$(L)(random.o) \
$(L)(ranz.o) \
$(L)(realft.o) \
$(L)(scopy.o) \
$(L)(sdot.o) \
$(L)(sevalb.o) \
$(L)(smooth2.o) \
$(L)(sort1.o) \
$(L)(spntgl.o) \
$(L)(sscal.o) \
$(L)(ssum.o) \
$(L)(tomsqz.o) \
$(L)(u_erf.o) \
$(L)(u_lu_backsub.o) \
$(L)(u_lu_decomp.o) \
$(L)(zeroin.o) \
$(L)(zgeev.o) \
$(L)(zort.o)

all:	$(MODULES)  
	ld -o $(LMLIB) -whole-archive -shared $(L)

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
