#GNU Makefile to build portlib.a 
#        and test driver porttest
#
#05/05/99 C. Ludescher
#

# define system dependent flags, etc.
-include ../../Make.local
include ../../Make.flags

# portlib special flags not covered by Make.depends:
FFLAGS := $(FFLAGS) $(STATIC)

# destination for processed source
# Check if user specified LOCAL or OBJ
# if run as a sub_make
ifneq ($(MAKELEVEL),0)
# check for LOCAL in case invoked via ufiles
ifeq ("${LOCAL}",".")
	OBJ = ..
endif
# check for LOCAL in case invoked via sglib
ifeq ("${LOCAL}","..")
	OBJ = ..
endif
# if envoced via anything else
ifeq ("$(OBJ)",".")
	override OBJ = ..
endif
endif

ifneq ("$(wildcard $(LIBROOT)/include/readline*)","")
	INCFLAGS := $(INCFLAGS) $(IOPT)$(LIBROOT)/include
endif

ifdef DEBUG
	FFLAGS   = $(DFFLAGS) $(STATIC)
	CFLAGS   = $(DCFLAGS) $(STATIC)
	BINDIR   = $(DEBDIR)
	LIBDIR   = $(DBLDIR)
endif

OBJDIR = $(OBJ)/obj/portlib
ARC = $(LIBDIR)/libportlib.a

ifdef EDITLIBS
	CDEFS := $(CDEFS) -D__GETLINE_EDITOR
endif

# don't rebuild library if using installed one in $PREFIX 
ifeq ($(MAKELEVEL),0)
	THISLIB=$(ARC)
endif

# .f9f and .F sources
F9 = $(subst .F90,.o, $(wildcard *.F90)) $(subst .f90,.o, $(wildcard *.f90))
F =  $(subst .F,.o, $(wildcard *.F)) $(subst .f,.o, $(wildcard *.f))
C =  $(subst .c,.o,$(wildcard *.c))
ALLM = $(C) $(F9) $(F)
MEMS = $(filter-out porttest.o,$(ALLM)) 
MEM = $(foreach file,$(MEMS),$(ARC)($(file)))



.PHONY: libs ckdirs clean realclean

libs:  ckdirs $(ARC)
	@if test ! -f $(OBJ)/lib/portlib.a; then \
	(cd $(OBJ)/lib; ln -s libportlib.a portlib.a); fi

$(LPORTLIB): $(ARC)
	ld -o $(LPORTLIB) -whole-archive -shared $(ARC)

$(ARC): $(MEM)
	@$(RANLIB) $@ > /dev/null

$(OLDLIB): timestamp.inf
	@echo "--- DETECTED portlib.a source update"
	@echo "--- Re-Making $(ARC)"
	@$(MAKE) libs

# pre-process / compile free form       
$(ARC)(%.o): %.F90
	$(F90_proc)
	$(make_ar)

# pre-pocess / compile fixed form
$(ARC)(%.o): %.F
	$(F_proc)
	$(make_ar)

# compile fixed
$(ARC)(%.o): %.f
	$(FC) -c $(FFLAGS) $(OUT) $(OBJDIR)/$*.o $<
	$(make_ar)

# compile c
$(ARC)(%.o): %.c
	$(c_proc)
	$(make_ar)

# pre-pocess / compile fixed
$(OBJDIR)/%.o: %.F
	$(F_proc)

# create destination directories if needed
ckdirs:
	@test -d $(OBJ)/lib || mkdir -p $(OBJ)/lib
	@test -d $(OBJ)/obj/portlib || mkdir -p $(OBJ)/obj/portlib


all:  libs exec
	@echo done  

exec: chkdir2 $(OBJ)/test/porttest

chkdir2:
	@test -d $(OBJ)/test || mkdir -p $(OBJ)/test
	@test -d $(OBJ)/obj/portlib || mkdir -p $(OBJ)/obj/portlib

$(OBJ)/test/porttest: $(OBJDIR)/porttest.o $(THISLIB) 
	$(FC) $(LDFLAGS) -o $@ $< -L$(OBJ)/lib $(LLOC) -lportlib $(XFLIB)	


install: 
	@test -d $(LIBDIR) || mkdir -p $(LIBDIR)
	@if test -f $(ARC); then \
	  umask 133; cp $(ARC) $(LIBDIR)/; \
	else test -f $(LIBDIR)/libportlib.a || echo ">>> can not find $(ARC)";fi;
	@test -d $(MANDIR)/man3 || mkdir -p $(MANDIR)/man3
	@umask 133; cp portlib.3 $(MANDIR)/man3/

uninstall:
	rm $(LIBDIR)/libportlib.a
	rm $(MANDIR)/man3/portlib.3

clean: 
	@rm -f $(OBJDIR)/*.o 

realclean: clean
	@rm -f $(ARC)
	@rm -f $(OBJ)/test/porttest

