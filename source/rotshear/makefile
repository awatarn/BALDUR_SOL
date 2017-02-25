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
L = $(LIBDIR)/rotshear.a
X = $(BINDIR)/rotshear

objects = \
$(L)(move1d.o) \
$(L)(wexbprof2.o) \
$(L)(ps_read.o) \
$(L)(wexbprof.o)

#$(L)(ps_read.o) \

# output: input $(X)
#	$(X) < input

.SUFFIXES:
.SUFFIXES: .o .f90 .a .F90 .for .f .F

$(X) :	$(L)
	$(FC) elf.f $(FFLAGS) -o $(X) $(L)

$(L) :	$(objects)

.f.a:
	$(FC) $(FFLAGS) -c $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o

.f90.a:
	$(FC) $(FFLAGS) -c $<
	$(AR) $(ARFLAGS) $(L) $*.o
	rm -f $*.o

