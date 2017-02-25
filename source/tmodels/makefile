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
L=$(LIBDIR)/tmodels.a

MODULES =  \
$(L)(callglf2d.o) \
$(L)(callglf2db.o) \
$(L)(etaw14a.o) \
$(L)(etaw14diff.o) \
$(L)(etaw14flux.o) \
$(L)(etaw17a.o) \
$(L)(etaw17diff.o) \
$(L)(etaw17flux.o) \
$(L)(etaw18.o) \
$(L)(glf2d.o) \
$(L)(kbmodels.o) \
$(L)(mixed_Bohm_gyro_Bohm.o) \
$(L)(mixed_merba.o) \
$(L)(mixerba.o) \
$(L)(mmm95.o) \
$(L)(mmm95a.o) \
$(L)(mmm98.o) \
$(L)(mmm98b.o) \
$(L)(mmm98c.o) \
$(L)(mmm98d.o) \
$(L)(mmm99.o) \
$(L)(mmm2001.o) \
$(L)(ohe.o) \
$(L)(ohe_model.o) \
$(L)(neo_pol.o) \
$(L)(ptheory.o) \
$(L)(sda01dif.o) \
$(L)(sda02flx.o) \
$(L)(sda04dif.o) \
$(L)(sda04flx.o) \
$(L)(sda05dif.o) \
$(L)(sda05flx.o) \
$(L)(theory.o) \
$(L)(vpolprof.o) \
$(L)(wexbint.o) \
$(L)(wexbprof.o) \
$(L)(weiland18.o) \
$(L)(weiland18diff.o) \
$(L)(weiland18disp.o) \
$(L)(weiland18eqns.o) \
$(L)(weiland18flux.o) \
$(L)(weiland18init.o) \
$(L)(weiland14.o) 

# $(L)(mixed_model.o) \

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
