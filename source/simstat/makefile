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
L = $(LIBDIR)/simstat.a
X = $(BINDIR)/simstat
#
objects = \
$(L)(simstat.o) \
$(L)(maxlike.o) \
$(L)(maxlikfn.o) \
$(L)(chisqr.o) \
$(L)(expdata.o) \
$(L)(baldata.o) \
$(L)(rdata0.o) \
$(L)(prtext.o) \
$(L)(prtrve.o) \
$(L)(prtrvf.o) \
$(L)(int1d.o)
#
$(X) :	$(L)
	$(FC) simstat.f $(FFLAGS) -o $(X) $(L)
	rm -f $*.o

$(L) :	$(objects)

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o

#output : simstat input
#	-rm output
#	./simstat < input
