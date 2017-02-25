# makefile for the BALDUR code, Glenn Bateman, Lehigh University
# bateman@plasma.physics.lehigh.edu
# Changes by Alexei Pankin, 31 Jul 2002 pankin@fusion.physics.lehigh.edu

include ../../Make.flags
include ../../Make.local
#

ifdef DEBUG
	FFLAGS = $(DFFLAGS)
	CFLAGS = $(DCFLAGS)
	BINDIR = $(DEBDIR)
	LIBDIR = $(DBLDIR)
endif

# modules 
MODU    = \
$(L)(abortb.o) \
$(L)(clear.o) \
$(L)(daytim.o) \
$(L)(dprtutil.o) \
$(L)(saxpy.o) \
$(L)(second.o) \
$(L)(stripx.o)

# modules derived from *.f files
MODF    = \
$(L)(d3hefus.o) \
$(L)(dadpak.o) \
$(L)(dauxheat.o) \
$(L)(dbeams.o) \
$(L)(dempirc.o) \
$(L)(deqbald.o) \
$(L)(deqmom2.o) \
$(L)(deqvmom.o) \
$(L)(dfusion.o) \
$(L)(dio.o) \
$(L)(dislbal.o) \
$(L)(divertor.o) \
$(L)(dneugas.o) \
$(L)(dolymp.o) \
$(L)(doutput.o) \
$(L)(dpellet.o) \
$(L)(dripple.o) \
$(L)(dscrunch.o) \
$(L)(dtransp.o) \
$(L)(labrun.o) \
$(L)(inital.o) \
$(L)(auxval.o) \
$(L)(errchk.o) \
$(L)(preset.o) \
$(L)(units.o) \
$(L)(oclear.o) \
$(L)(otable.o) \
$(L)(ofit.o) \
$(L)(cotrol.o) \
$(L)(stepon.o) \
$(L)(coef.o) \
$(L)(convrt.o) \
$(L)(cnvcof.o) \
$(L)(solveb.o) \
$(L)(noncor.o) \
$(L)(bounds.o) \
$(L)(reduce.o) \
$(L)(solve.o) \
$(L)(resolv.o) \
$(L)(getchi.o) \
$(L)(field2d.o) \
$(L)(ftrapfn.o) \
$(L)(ftrap_hughes.o) \
$(L)(boots.o) \
$(L)(imprad.o) \
$(L)(ncinfl.o) \
$(L)(ncinit.o) \
$(L)(ncdata.o) \
$(L)(ncsorc.o) \
$(L)(ncrats.o) \
$(L)(ncdifu.o) \
$(L)(nccoef.o) \
$(L)(ncsolv.o) \
$(L)(nccons.o) \
$(L)(ncprof0.o) \
$(L)(neudep.o) \
$(L)(neuset.o) \
$(L)(syedit.o) \
$(L)(synch.o) \
$(L)(syopac.o) \
$(L)(syndrv.o) \
$(L)(synch2.o) \
$(L)(syopac2.o) \
$(L)(syedit2.o) \
$(L)(hmode.o)\
$(L)(edge.o)\
$(L)(height1.o)\
$(L)(collision.o)

L       = $(LIBDIR)/libbald.a

MODULES = $(MODF) $(MODU)  # all modules

$(L)	:	$(MODULES)

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

