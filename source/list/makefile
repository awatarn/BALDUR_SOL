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
L = $(LIBDIR)/list.a
X = $(BINDIR)/list

objects = \
$(L)(b2read.o) \
$(L)(r1read.o) \
$(L)(r2read.o) \
$(L)(int1d.o) \
$(L)(strim.o)

# output: input $(X)
#	$(X) < input

$(X) :	$(L)
	$(FC) list.f $(FFLAGS) -o $(X) $(L)
	rm -f $*.o

$(L) :	$(objects)

.f.a:
	$(FC) -c $(FFLAGS) -static -TENV:X=4 -g $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o


