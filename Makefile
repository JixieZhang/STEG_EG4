#
# Makefile basis, modified by Jixie @20180306
# compile every *.f in $(SRCDIR)
# identify program name by itself

# Linux compiler
COMPILER = gfortran

SRCDIR := src
INCDIR := src

###############################################################
ifndef PROG
  NFILES   = $(shell ls -1 $(SRCDIR)/*.f | wc -l) 

#kp: Gets the program name (e.g. in case of generate_map dir. it's 'steg' and in case of read_map/work/ dir. it's 'read_map'
  ifeq ($(NFILES),1)
    PROG     = $(shell grep -i program $(SRCDIR)/*.f | awk '$$1=="program"&&$$3==""{print $$2}' )
  else
    PROG     = $(shell grep -i program $(SRCDIR)/*.f | awk '$$2=="program"&&$$4==""{print $$3}' )
  endif
endif
###############################################################

SOURCES  = $(wildcard $(SRCDIR)/*.[fF])
INC      = $(wildcard $(SRCDIR)/*.inc)
OBJ     := $(addsuffix .o, $(basename $(SOURCES)))
OBJ     := $(patsubst  $(SRCDIR)/%.o,obj/%.o,$(OBJ))

###############################################################
# SYSTEM Libraries
ifndef SYSLIBS
  SYSLIBS  = -lnsl
endif

#CERN LIBS
ifndef CERNLIBS
  CERNLIBS = -lpdflib804 -lmathlib -lphtools -lgeant321 \
 -lpacklib -lkernlib -lpawlib
endif

#CLAS LIBS
ifndef CLASLIBS
  CLASLIBS = -lbankdefs -lclasutil -lbosio -lc_bos_io  \
  -lrecutl -lcc -lfputil -lbos -lfpack -lmapmanager 
 endif

#stitch all libs together
LIBS = $(SYSLIBS) -L$(CLAS_LIB) $(CLASLIBS) -L$(CERN_ROOT)/lib $(CERNLIBS)

###############################################################
#flags	
ifndef FFLAGS
  FFLAGS += -m64 -Wall 
  FFLAFS += -fno-f2c -DLinux -fno-automatic -finit-local-zero
  FFLAGS += -ffixed-line-length-none -fno-second-underscore -funroll-loops
  FFLAGS += -fomit-frame-pointer -fno-range-check -falign-commons
  FFLAGS += -I$(INCDIR)
endif	
#
# Note: The -fno-f2c option is very important for mixed FORTRAN/C environment
# where there are C functions that return a float to some FORTRAN code. The "f2c" 
# calling convention wants all routines to use double only. CLAS code is not
# consistend in this.
#	 

ifdef DEBUG
	FFLAGS +=  -g -ffloat-store
	ADD_DEBUG = _debug
#	EXE   = $(PROG)_$(OS_NAME)_${ADD_DUBUG}	
	EXE   = $(PROG)_rhel$(OS_VERSION)_${ADD_DUBUG}	
else 
	ADD_DEBUG = 
#	EXE   = $(PROG)_$(OS_NAME)
	 EXE   = $(PROG)_rhel${OS_VERSION}
endif

###############################################################
#compile rules
obj/%.o: $(SRCDIR)/%.f $(INC)
	$(COMPILER) $(FFLAGS) -c $< -o $@

%.o:    %.f $(INC)
	$(COMPILER) $(FFLAGS) -c $< -o $@

###############################################################
all:    $(EXE)

$(EXE): dir $(OBJ) $(INC)
	$(COMPILER) $(FFLAGS) -o $(EXE) $(OBJ) $(LIBS)

dir:
	@mkdir -p obj

clean:
	rm -f  obj/*.o core

distclean: clean
	rm -f $(EXE)

# 
# 8/14/17:
# https://stackoverflow.com/questions/16467718/how-to-print-out-a-variable-in-makefile
# https://www.cmcrossroads.com/article/printing-value-makefile-variable
# Command to use to check the value of X:
#        gmake -f Makefile -f helper.mak print-X 
#
# kp: Remember that we have a tab before @echo on the third line below
#
print-%:
	@echo $* = $($*)

