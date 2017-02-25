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
LL = $(LIBDIR)/list.a 
LE = $(LIBDIR)/ezcdf.a
X = $(BINDIR)/listcdf

$(X) :	$(LL) $(LE) source
	$(FC) listcdf.o $(FFLAGS) -o $(X) $(LL) $(LE) -lnetcdf
	rm -f listcdf.o

$(LL) :	$(objects)
	cd ../list; make

$(LE) :
	cd ../ezcdf; make

source:
	$(FC) -c $(FFLAGS) -static -TENV:X=4 -g listcdf.f


