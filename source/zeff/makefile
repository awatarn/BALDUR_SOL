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
X = $(BINDIR)/zeff

$(X) :	
	$(FC) zeff.f $(FFLAGS) -o $(X) 

.f.a:
	$(FC) -c $(FFLAGS) -static -TENV:X=4 -g $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o


