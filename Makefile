

FNMA_DIR ?= $(shell pwd)

# Default to debug flags
FNMA_FFLAGS ?= -O2 -g

FNMA_SRCDIR = $(FNMA_DIR)/src/
FNMA_EXASRCDIR = $(FNMA_DIR)/examples/


# Out-of-source build directories
FNMA_INCDIR = $(FNMA_DIR)/build/include/
FNMA_LIBDIR = $(FNMA_DIR)/build/lib/
FNMA_BINDIR = $(FNMA_DIR)/build/bin/
FNMA_EXADIR = $(FNMA_DIR)/build/examples/

vpath %.f90 $(FNMA_SRCDIR)

FNMA_F90_SRCS = fnma_Precision \
	        fnma_Mesh \
		fnma_Metadata \
		fnma_DataTypes

#FNMA_EXAMPLES = NeumannSquare \
#	       IrregularGeometry \
#	       IrregularGeometryNeumann \
#	       CircularGeometry \
#	       CircularGeometryNeumann \
#	       MITgcmGeometry

FNMA_LIBS = fnma

FNMA_OBJS = $(addprefix $(FNMA_SRCDIR), $(addsuffix .f.o, $(FNMA_F90_SRCS)))
FNMA_LIB_OBJS = $(addprefix $(FNMA_LIBDIR)lib, $(addsuffix .a, $(FNMA_LIBS)))
FNMA_BUILDDIRS = $(FNMA_INCDIR) $(FNMA_LIBDIR) $(FNMA_BINDIR) $(FNMA_EXADIR) $(FNMA_TESTDIR)

# Example Programs
#FNMA_EXAS = $(addprefix $(FNMA_EXADIR), $(FNMA_EXAMPLES))

# Add preprocessing to flags (Assumes gnu compilers)
FNMA_FFLAGS += -cpp


fnma: $(FNMA_LIB_OBJS)

#examples: $(FNMA_EXAS)

#$(FNMA_EXAS): $(FNMA_DIR)/build/%: %.f90 $(FNMA_OBJS)
#	$(FC) -c $(FNMA_FFLAGS) -I$(FNMA_INCDIR) $< -o $<.o
#	$(FC) $(FNMA_FFLAGS) -I$(FNMA_INCDIR) $(FNMA_SRCDIR)*.o  $<.o $(FNMA_FLIBS) -o $@
#rm $<.o

$(FNMA_LIBDIR)libfnma.a: $(FNMA_OBJS)
	rm -f $@
	$(AR) -cq $@ $^

$(FNMA_SRCDIR)%.f.o: %.f90
	$(FC) -J$(FNMA_INCDIR) $(FNMA_FFLAGS) $(FNMA_FLIBS) -c $< -o $@

clean:
	rm -f $(FNMA_BINDIR)*
	rm -f $(FNMA_LIBDIR)*.a
	rm -f $(FNMA_MODDIR)*.mod
	rm -f $(FNMA_SRCDIR)*.o
	rm -f $(FNMA_EXADIR)*.o
	rm -f $(FNMA_EXASRCDIR)*.o

# Dependency on build tree existence
$(FNMA_OBJS): | $(FNMA_BUILDDIRS)

$(FNMA_BUILDDIRS):
	mkdir -p $@

.PHONY: clean
