

SLS_DIR ?= $(shell pwd)

# Default to debug flags
SLS_FFLAGS ?= -O0 -g -Wall

SLS_LINK ?= -L$(LAPACK)/lib -lopenblas 
SLS_INCLUDE ?= -I$(LAPACK)/include 

SLS_FLIBS += $(SLS_LINK) $(SLS_INCLUDE)

SLS_SRCDIR = $(SLS_DIR)/src/


# Out-of-source build directories
SLS_INCDIR = $(SLS_DIR)/build/include/
SLS_LIBDIR = $(SLS_DIR)/build/lib/
SLS_BINDIR = $(SLS_DIR)/build/bin/
SLS_EXADIR = $(SLS_DIR)/build/examples/

vpath %.f90 $(SLS_SRCDIR)

SLS_F90_SRCS = SLSpectra_Precision \
	       SLSpectra_SupportRoutines \
	       SLSpectra_Stencil \
	       SLSpectra_Mesh \
   	       SLSpectra_AdjacencyGraph \
	       SLSpectra_Generators 

SLS_EXAMPLES = DirichletSquare \
	       NeumannSquare \
	       IrregularGeometry \
	       IrregularGeometryNeumann \
	       CircularGeometry \
	       CircularGeometryNeumann \
	       MITgcmGeometry

SLS_LIBS = slspectra

SLS_OBJS = $(addprefix $(SLS_SRCDIR), $(addsuffix .f.o, $(SLS_F90_SRCS)))
SLS_LIB_OBJS = $(addprefix $(SLS_LIBDIR)lib, $(addsuffix .a, $(SLS_LIBS)))
SLS_BUILDDIRS = $(SLS_INCDIR) $(SLS_LIBDIR) $(SLS_BINDIR) $(SLS_EXADIR) $(SLS_TESTDIR)

# Example Programs
SLS_EXAS = $(addprefix $(SLS_EXADIR), $(SLS_EXAMPLES))

# Add preprocessing to flags (Assumes gnu compilers)
SLS_FFLAGS += -cpp


slspectra: $(SLS_LIB_OBJS)

examples: $(SLS_EXAS)

$(SLS_EXAS): $(SLS_DIR)/build/%: %.f90 $(SLS_OBJS)
	$(FC) -c $(SLS_FFLAGS) -I$(SLS_INCDIR) $< -o $<.o
	$(FC) $(SLS_FFLAGS) -I$(SLS_INCDIR) $(SLS_SRCDIR)*.o  $<.o $(SLS_FLIBS) -o $@
	rm $<.o

$(SLS_LIBDIR)libslspectra.a: $(SLS_OBJS)
	rm -f $@
	$(AR) -cq $@ $^

$(SLS_SRCDIR)%.f.o: %.f90
	$(FC) -J$(SLS_INCDIR) $(SLS_FFLAGS) $(SLS_FLIBS) -c $< -o $@

clean:
	rm -f $(SLS_BINDIR)*
	rm -f $(SLS_LIBDIR)*.a
	rm -f $(SLS_MODDIR)*.mod
	rm -f $(SLS_SRCDIR)*.o
	rm -f $(SLS_EXADIR)*.o

# Dependency on build tree existence
$(SLS_OBJS): | $(SLS_BUILDDIRS)

$(SLS_BUILDDIRS):
	mkdir -p $@

.PHONY: clean
