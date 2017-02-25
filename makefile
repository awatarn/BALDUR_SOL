#include Make.flags
include Make.local
SRCDIR=./source
BINDIR=./bin
LIBDIR=./lib
DOCDIR=./doc

all: portlib ezcdf cdfio bald baldur rotshear list listcdf getprofile getprofcdf zeff simstat pedestal

baldur:
	cd ./$(SRCDIR)/baldur; \
	echo $(HDF); \
	make 

bald:	cdfio
	cd ./$(SRCDIR)/bald; \
	make 

cdfio:
	cd ./$(SRCDIR)/cdfio; \
	make 

ezcdf: 
	cd ./$(SRCDIR)/ezcdf; \
	make

portlib: 
	cd ./$(SRCDIR)/portlib; \
	make

list:
	cd ./$(SRCDIR)/list; \
	make 

listcdf:
	cd ./$(SRCDIR)/listcdf; \
	make 

getprofile:
	cd ./$(SRCDIR)/profile; \
	make 

getprofcdf:
	cd ./$(SRCDIR)/profilecdf; \
	make 

rotshear:
	cd ./$(SRCDIR)/rotshear; \
	make 

zeff:
	cd ./$(SRCDIR)/zeff; \
	make 

simstat:
	cd ./$(SRCDIR)/simstat; \
	make 

pedestal:
	cd./$(SRCDIR)/pedestal; \
        make

clean:
	rm -f bin/* bin/debug/*

test: baldur
	cd test; \
	../csh/run.sh test a01 a02 a03 a04
