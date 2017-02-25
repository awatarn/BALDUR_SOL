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
X         = $(BINDIR)/xbaldur

all:	LBALD LPORTLIB LMLIB LSAW LNEO LTRANSP LPEDESTAL
	$(FC) main.f $(FFLAGS) -o $(X) $(LBALD) $(LPORTLIB) $(LSAW) $(LPEDESTAL) $(LNEO) $(LTRANSP) $(LMLIB) $(LIB_TRAILER)
	@ echo 'new $(X) constructed'

LSAW:
	cd ../sawtooth; \
	make 

LPEDESTAL:
	cd ../pedestal; \
	make 

LNEO:
	cd ../neoclass; \
	make

LTRANSP:
	cd ../tmodels; \
	make

LMLIB:
	cd ../mlib; \
	make

LBALD:
	cd ../bald; \
	make 

LPORTLIB: 
	cd ../portlib; \
	make

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

