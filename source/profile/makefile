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

L = $(LIBDIR)/getprofile.a
X = $(BINDIR)/getprofile

objects = \
$(L)(bpread.o) \
$(L)(br2read.o) \
$(L)(r2read.o) \
$(L)(int1d.o) \
$(L)(strim.o)

#output:	$(X) input
#	$(X) < input

$(X):	$(L) 
	$(FC) $(FFLAGS) getprofile.f -o $(X) $(L)
	rm -f $*.o

$(L):	$(objects)

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $(L) $*.o 
	rm -f $*.o
