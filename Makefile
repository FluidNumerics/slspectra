

XEMODES_DIR ?= $(shell pwd)

XEMODES_SRCDIR = $(XEMODES_DIR)/src/

# Out-of-source build directories
XEMODES_INCDIR = $(XEMODES_DIR)/build/include/
XEMODES_LIBDIR = $(XEMODES_DIR)/build/lib/
XEMODES_BINDIR = $(XEMODES_DIR)/build/bin/
XEMODES_EXADIR = $(XEMODES_DIR)/build/examples/

vpath %.f90 $(XEMODES_SRCDIR)

XEMODES_F90_SRCS = XEModes_Stencil 

XEMODES_LIBS = xemodes

XEMODES_OBJS = $(addprefix $(XEMODES_SRCDIR), $(addsuffix .f.o, $(XEMODES_F90_SRCS)))
XEMODES_LIB_OBJS = $(addprefix $(XEMODES_LIBDIR)lib, $(addsuffix .a, $(XEMODES_LIBS)))
XEMODES_BUILDDIRS = $(XEMODES_INCDIR) $(XEMODES_LIBDIR) $(XEMODES_BINDIR) $(XEMODES_EXADIR) $(XEMODES_TESTDIR)


xemodes: $(XEMODES_LIB_OBJS)


$(XEMODES_LIBDIR)libxemodes.a: $(XEMODES_OBJS)
	rm -f $@
	$(AR) -cq $@ $^

$(XEMODES_SRCDIR)%.f.o: %.f90
	$(FC) -J$(XEMODES_INCDIR) $(XEMODES_FFLAGS) $(XEMODES_FLIBS) -c $< -o $@

clean:
	rm -f $(XEMODES_BINDIR)*
	rm -f $(XEMODES_LIBDIR)*.a
	rm -f $(XEMODES_MODDIR)*.mod
	rm -f $(XEMODES_SRCDIR)*.o
	rm -f $(XEMODES_EXASRCDIR)*.o

# Dependency on build tree existence
$(XEMODES_OBJS): | $(XEMODES_BUILDDIRS)

$(XEMODES_BUILDDIRS):
	mkdir -p $@

.PHONY: clean
