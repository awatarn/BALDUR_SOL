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

LL = $(LIBDIR)/getprofile.a
#LE = $(LIBDIR)/ezcdf.a
X  = $(BINDIR)/getprofcdf

$(X):	$(LL) $(LE) src
	$(FC) $(FFLAGS) getprofcdf.o -o $(X) $(LL) /usr/local/lib/libnetcdf.so
	rm -f getprofcdf.o

$(LL):	
	cd ../getprofile;  make

$(LE):	
	cd ../ezcdf;  make

src:
	$(FC) -c $(FFLAGS) getprofcdf.f
